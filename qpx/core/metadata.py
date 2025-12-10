"""
Metadata Schema Generator for QPX format.

This module generates metadata configuration files that document column mappings,
ontology terms, and computation methods for different QPX data models (PSM, Feature, PG).
"""

import csv
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Set
import pyarrow as pa

logger = logging.getLogger(__name__)


def _get_qpx_standard_fields() -> Set[str]:
    """Extract QPX standardized field names from PyArrow schemas in format.py."""
    from qpx.core.format import PSM_SCHEMA, FEATURE_SCHEMA, PG_SCHEMA

    qpx_fields = set()

    for schema in [PSM_SCHEMA, FEATURE_SCHEMA, PG_SCHEMA]:
        for field in schema:
            qpx_fields.add(field.name)

    return qpx_fields


QPX_STANDARD_FIELDS: Set[str] = _get_qpx_standard_fields()


class MetadataSchemaGenerator:
    """
    Generates metadata configuration files documenting QPX format columns.

    Format: model_view | column_name | ontology_accession | full_name | settings

    Example:
        psm | charge | MS:1000041 | charge state | computed as provided by MaxQuant
    """

    def __init__(self):
        self.metadata_records: List[Dict[str, str]] = []
        self.ontology_mapping = self._load_ontology_mapping()
        self.available_columns = None

    def add_column_metadata(
        self,
        model_view: str,
        column_name: str,
        ontology_accession: Optional[str] = None,
        full_name: Optional[str] = None,
        settings: Optional[str] = None,
        file_name: Optional[str] = None,
    ):
        """
        Add metadata for a single column.

        Args:
            model_view: Data model type (psm, feature, pg, ibaq)
            column_name: Name of the column in QPX format
            ontology_accession: Ontology term (e.g., MS:1000041, PRIDE:0000123)
            full_name: Full descriptive name of the column
            settings: Description of how the value is computed or obtained
            file_name: Source data file name (e.g., msms.txt, evidence.txt)
        """
        record = {
            "model_view": model_view,
            "file_name": file_name or "",
            "column_name": column_name,
            "ontology_accession": ontology_accession or "",
            "full_name": full_name or "",
            "settings": settings or "",
        }
        self.metadata_records.append(record)

    def add_from_schema(
        self,
        model_view: str,
        schema: pa.Schema,
        source_mapping: Optional[Dict] = None,
        settings_mapping: Optional[Dict] = None,
    ):
        """
        Extract metadata from PyArrow schema and add to records.

        Args:
            model_view: Data model type (psm, feature, pg, ibaq)
            schema: PyArrow schema to extract metadata from
            source_mapping: Dict mapping source column names to QPX column names
            settings_mapping: Dict mapping QPX column names to computation descriptions
        """
        source_mapping = source_mapping or {}
        settings_mapping = settings_mapping or {}

        for field in schema:
            field_metadata = field.metadata or {}
            description = field_metadata.get(b"description", b"").decode("utf-8")

            ontology_info = self.ontology_mapping.get(field.name, {})
            ontology_accession = ontology_info.get("ontology_accession", "")
            ontology_description = ontology_info.get("description", "")

            if not ontology_accession:
                ontology_accession = self._extract_ontology_term(description) or ""

            final_description = (
                ontology_description if ontology_description else description
            )

            reverse_source_map = {v: k for k, v in source_mapping.items()}
            source_col = reverse_source_map.get(field.name, "")
            settings = settings_mapping.get(field.name, "")

            # Skip fields that have no source column and no computation description
            # when available_columns filtering is active
            if self.available_columns is not None:
                if not source_col and not settings:
                    # This is a standard field with no mapping - skip it
                    continue

            if field.name in QPX_STANDARD_FIELDS:
                column_name = field.name
            elif source_col:
                column_name = source_col
            else:
                column_name = field.name

            if settings:
                if source_col:
                    full_settings = f"computed from '{source_col}' as {settings}"
                else:
                    full_settings = f"computed as {settings}"
            elif source_col:
                full_settings = f"provided by '{source_col}'"
            else:
                if field.name in QPX_STANDARD_FIELDS:
                    full_settings = "as defined in QPX specification"
                else:
                    full_settings = ""

            self.add_column_metadata(
                model_view=model_view,
                column_name=column_name,
                ontology_accession=ontology_accession,
                full_name=final_description,
                settings=full_settings,
                file_name=(
                    self._get_default_filename(model_view)
                    if hasattr(self, "_get_default_filename")
                    else ""
                ),
            )

    def _extract_ontology_term(self, text: str) -> Optional[str]:
        """
        Extract ontology term from description text.

        Looks for patterns like MS:XXXXX, PRIDE:XXXXX, etc.
        """
        import re

        if not text:
            return None

        patterns = [
            r"(MS:\d+)",
            r"(PRIDE:\d+)",
            r"(UNIMOD:\d+)",
            r"(PSI-MS:\d+)",
        ]

        for pattern in patterns:
            match = re.search(pattern, text)
            if match:
                return match.group(1)

        return None

    def _load_ontology_mapping(self) -> Dict[str, Dict[str, str]]:
        """
        Load ontology mapping from docs/field_msg.tsv file.

        Returns:
            Dict mapping field names to ontology information
        """
        try:
            current_file = Path(__file__)
            project_root = current_file.parent.parent.parent
            ontology_file = project_root / "docs" / "field_msg.tsv"

            if not ontology_file.exists():
                return {}

            df = pd.read_csv(ontology_file, sep="\t")

            df.columns = df.columns.str.strip()

            mapping = {}
            for _, row in df.iterrows():
                field_name = row.get("field", "")
                if pd.notna(field_name) and field_name:
                    field_name = str(field_name).strip()
                    mapping[field_name] = {
                        "ontology_accession": (
                            row.get("ontology_accession", "")
                            if pd.notna(row.get("ontology_accession"))
                            else ""
                        ),
                        "ontology_name": (
                            row.get("ontology_name", "")
                            if pd.notna(row.get("ontology_name"))
                            else ""
                        ),
                        "description": (
                            row.get("desc", "") if pd.notna(row.get("desc")) else ""
                        ),
                    }

            return mapping

        except Exception as e:
            logger.warning(f"Could not load ontology mapping: {e}")
            return {}

    def generate_file(self, output_path: str):
        """
        Generate metadata configuration file in CSV format.

        Args:
            output_path: Path to output CSV file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", newline="", encoding="utf-8") as f:
            fieldnames = [
                "model_view",
                "file_name",
                "column_name",
                "ontology_accession",
                "full_name",
                "settings",
            ]

            writer = csv.DictWriter(
                f,
                fieldnames=fieldnames,
                delimiter=",",
                quoting=csv.QUOTE_MINIMAL,
            )

            writer.writeheader()
            writer.writerows(self.metadata_records)

    def clear(self):
        """Clear all metadata records."""
        self.metadata_records = []

    def get_records(self) -> List[Dict[str, str]]:
        """Get all metadata records."""
        return self.metadata_records.copy()


class WorkflowMetadataGenerator(MetadataSchemaGenerator):
    """
    Workflow-specific metadata generator for QPX formats.

    Supports workflows: MaxQuant, DIANN, quantms-LFQ, quantms-TMT, quantms-PSM.
    Generates metadata configuration files documenting column mappings, ontology terms,
    and computation methods for each supported proteomics workflow.
    """

    def __init__(self, workflow: str, available_columns: Optional[Set[str]] = None):
        """
        Initialize workflow-specific metadata generator.

        Args:
            workflow: Workflow name (maxquant, diann, quantms-lfq, quantms-tmt, quantms-psm)
            available_columns: Set of column names available in input file (optional).
                             If provided, only metadata for these columns will be generated.
        """
        super().__init__()

        self.workflow = workflow.lower()
        self.available_columns = available_columns
        self.reset_mappings()
        self.init_workflow_mappings()

        if self.available_columns is not None:
            self._filter_mappings_by_available_columns()

    def reset_mappings(self):
        """Reset all mappings to empty state."""
        self.psm_source_map = {}
        self.feature_source_map = {}
        self.pg_source_map = {}
        self.psm_computed = {}
        self.feature_computed = {}
        self.pg_computed = {}

    def _filter_mappings_by_available_columns(self):
        """Filter source mappings to only include columns available in input file."""
        if self.available_columns is None:
            return

        self.psm_source_map = {
            k: v for k, v in self.psm_source_map.items() if k in self.available_columns
        }

        self.feature_source_map = {
            k: v
            for k, v in self.feature_source_map.items()
            if k in self.available_columns
        }

        self.pg_source_map = {
            k: v for k, v in self.pg_source_map.items() if k in self.available_columns
        }

        # Filter computed fields that depend on source columns
        # ion_mobility is mapped from '1/K0' column, only include if source exists
        if "1/K0" not in self.available_columns:
            self.psm_computed.pop("ion_mobility", None)
            self.feature_computed.pop("ion_mobility", None)

    def init_workflow_mappings(self):
        """Initialize workflow-specific column mappings."""
        if self.workflow == "maxquant":
            self._init_maxquant_mappings()
        elif self.workflow in ("diann", "quantms-lfq", "quantms-tmt", "quantms-psm"):
            pass
        else:
            raise ValueError(f"Unsupported workflow: {self.workflow}")

    def _init_maxquant_mappings(self):
        """Initialize MaxQuant-specific column mappings."""
        from qpx.core.common import (
            MAXQUANT_PSM_MAP,
            MAXQUANT_FEATURE_MAP,
            MAXQUANT_PG_MAP,
            MAXQUANT_PSM_COMPUTED,
            MAXQUANT_FEATURE_COMPUTED,
            MAXQUANT_PG_COMPUTED,
        )

        # Create copies to avoid modifying global dictionaries
        self.psm_source_map = MAXQUANT_PSM_MAP.copy()
        self.feature_source_map = MAXQUANT_FEATURE_MAP.copy()
        self.pg_source_map = MAXQUANT_PG_MAP.copy()
        self.psm_computed = MAXQUANT_PSM_COMPUTED.copy()
        self.feature_computed = MAXQUANT_FEATURE_COMPUTED.copy()
        self.pg_computed = MAXQUANT_PG_COMPUTED.copy()

    def generate_psm_metadata(self):
        """Generate metadata for PSM model."""
        from qpx.core.format import PSM_SCHEMA

        self.add_from_schema(
            model_view="psm",
            schema=PSM_SCHEMA,
            source_mapping=self.psm_source_map,
            settings_mapping=self.psm_computed,
        )

    def generate_feature_metadata(self):
        """Generate metadata for Feature model."""
        from qpx.core.format import FEATURE_SCHEMA

        self.add_from_schema(
            model_view="feature",
            schema=FEATURE_SCHEMA,
            source_mapping=self.feature_source_map,
            settings_mapping=self.feature_computed,
        )

    def generate_pg_metadata(self):
        """Generate metadata for PG model."""
        from qpx.core.format import PG_SCHEMA

        self.add_from_schema(
            model_view="pg",
            schema=PG_SCHEMA,
            source_mapping=self.pg_source_map,
            settings_mapping=self.pg_computed,
        )

    def generate_all_metadata(self):
        """Generate metadata for all data models (PSM, Feature, PG)."""
        self.generate_psm_metadata()
        self.generate_feature_metadata()
        self.generate_pg_metadata()

    def _get_default_filename(self, model_view: str) -> str:
        """Get default source filename for workflow and model view."""
        defaults = {
            "maxquant": {
                "psm": "msms.txt",
                "feature": "evidence.txt",
                "pg": "proteinGroups.txt",
            },
            "diann": {
                "psm": "report.tsv",
                "feature": "report.tsv",
                "pg": "report.pg_matrix.tsv",
            },
            "quantms-lfq": {
                "psm": "mztab (PSM section)",
                "feature": "mztab (PEP section)",
                "pg": "mztab (PRT section)",
            },
            "quantms-tmt": {
                "psm": "mztab (PSM section)",
                "feature": "mztab (PEP section)",
                "pg": "mztab (PRT section)",
            },
        }

        return defaults.get(self.workflow, {}).get(model_view, "unknown")
