"""Orchestrate QuantMS MSstats and SDRF conversion to QPX."""

from __future__ import annotations

import json
import logging
from pathlib import Path

from qpx._version import __version__
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.converters.quantms_msstats.feature_adapter import (
    FEATURE_IDENTITY_COMPOSITE,
    QuantmsMsstatsFeatureAdapter,
)
from qpx.converters.sdrf import SdrfConverter
from qpx.core.constants import FEATURE, ONTOLOGY, RUN, SAMPLE
from qpx.core.files import sha256_file
from qpx.core.scores import field_ontology_entries

logger = logging.getLogger(__name__)

# Feature rows per written batch; also the Parquet row-group size.
_BATCH_SIZE = 500_000


class QuantmsMsstatsConverter(BaseOrchestrator):
    """Convert QuantMS legacy MSstats plus SDRF input to QPX Parquet."""

    def __init__(
        self,
        max_memory: str | None = None,
        max_cpus: int | None = None,
    ):
        self._max_memory = max_memory
        self._max_cpus = max_cpus

    def convert(
        self,
        msstats_file: str | Path,
        sdrf_file: str | Path,
        output_folder: str | Path,
        output_prefix: str | None = None,
        project_accession: str | None = None,
    ) -> str:
        """Run Feature, SDRF, and metadata conversion."""
        msstats_file = Path(msstats_file)
        sdrf_file = Path(sdrf_file)
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        prefix = output_prefix or msstats_file.stem.removesuffix(".sdrf_openms_design_msstats_in")
        if not prefix:
            prefix = "quantms_msstats"

        feature_path = output_folder / f"{prefix}.feature.parquet"
        resolved = self._write_features(msstats_file, sdrf_file, feature_path)

        with SdrfConverter(
            duckdb_memory=self._max_memory,
            duckdb_threads=self._max_cpus,
        ) as sdrf_converter:
            sdrf_converter.convert(
                sdrf_path=str(sdrf_file),
                sample_output=str(output_folder / f"{prefix}.sample.parquet"),
                run_output=str(output_folder / f"{prefix}.run.parquet"),
            )
            ontology_entries = sdrf_converter.run_ontology_entries()

        ontology_entries.extend(
            field_ontology_entries(
                view=FEATURE,
                resolved_mappings=resolved,
                tool_name="QuantMS/MSstats",
            )
        )
        self._write_ontology(output_folder, prefix, ontology_entries)
        provenance_records = self.provenance_records(msstats_file, sdrf_file, resolved)
        self._write_provenance(output_folder, prefix, provenance_records)
        self._write_dataset(
            output_folder,
            prefix,
            project_accession,
            software_name="QuantMS/MSstats",
            software_version=None,
            provenance_records=provenance_records,
        )
        return prefix

    def _write_features(self, msstats_file: Path, sdrf_file: Path, feature_path: Path) -> dict[str, str]:
        """Write feature batches and return the resolved source columns."""
        with QuantmsMsstatsFeatureAdapter(
            duckdb_memory=self._max_memory,
            duckdb_threads=self._max_cpus,
        ) as adapter:
            adapter.convert(
                msstats_path=str(msstats_file),
                sdrf_path=str(sdrf_file),
                output_path=str(feature_path),
                chunksize=_BATCH_SIZE,
            )
            return adapter.get_resolved_columns()

    @staticmethod
    def provenance_records(
        msstats_file: Path,
        sdrf_file: Path,
        resolved: dict[str, str],
    ) -> list[dict]:
        """Build source and conversion provenance records."""
        output_views = [FEATURE, SAMPLE, RUN, ONTOLOGY]
        config = json.dumps(
            {
                "source_columns": resolved,
                "identity_composite": list(FEATURE_IDENTITY_COMPOSITE),
                "missing_values": "typed null; scan uses an empty list when unavailable",
                "protein_aggregation": "not performed",
                "unique_peptide": ("true only when a sequence maps to one ProteinName value containing one accession"),
            },
            sort_keys=True,
        )
        return [
            {
                "step_order": 1,
                "step_category": "quantification",
                "step_name": "quantms_msstats_feature_quantification",
                "tool_name": "QuantMS/MSstats",
                "tool_version": None,
                "tool_uri": "https://github.com/bigbio/quantms",
                "parameters": None,
                "config": None,
                "output_views": [FEATURE],
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "quantms_msstats_to_qpx",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": "https://github.com/bigbio/qpx",
                "parameters": [
                    {"key": "msstats_file", "value": msstats_file.name},
                    {"key": "msstats_sha256", "value": sha256_file(msstats_file)},
                    {"key": "sdrf_file", "value": sdrf_file.name},
                    {"key": "sdrf_sha256", "value": sha256_file(sdrf_file)},
                ],
                "config": config,
                "output_views": output_views,
            },
        ]
