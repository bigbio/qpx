"""BaseOrchestrator — shared ontology, provenance, and dataset writing for converters.

Orchestrators (MaxQuantConverter, OpenMSConverter, etc.) compose adapters and
write ontology.parquet, provenance.parquet, and dataset.parquet. This base class
consolidates the common writing logic.
"""

from __future__ import annotations

import logging
import uuid
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)


def build_dataset_record(
    *,
    project_accession: str | None = None,
    software_name: str = "unknown",
    software_version: str | None = None,
    provenance_records: list[dict] | None = None,
) -> dict:
    """Build a dataset.parquet record with consistent field structure.

    Reduces dict key drift across converters. If provenance_records is given,
    software_name and software_version are taken from the first step when not
    explicitly provided.
    """
    if provenance_records:
        first = provenance_records[0]
        if software_name == "unknown":
            software_name = first.get("tool_name") or software_name
        if software_version is None:
            software_version = first.get("tool_version")
    return {
        "project_accession": project_accession or "unknown",
        "project_title": None,
        "project_description": None,
        "pubmed_id": None,
        "software_name": software_name,
        "software_version": software_version,
        "analysis_id": str(uuid.uuid4()),
        "creation_date": datetime.now().isoformat(),
        "file_checksums": None,
        "file_row_counts": None,
        "file_sizes_bytes": None,
        "total_structures": None,
        "packaged_at": None,
    }


def build_provenance_step(
    *,
    step_category: str,
    step_name: str,
    tool_name: str,
    tool_version: str | None = None,
    tool_uri: str | None = None,
    parameters: list[dict] | None = None,
    config: str | None = None,
    output_views: list[str] | None = None,
    step_order: int = 0,
) -> dict:
    """Build a provenance.parquet record with consistent field structure (``tool_version`` '' becomes None)."""
    return {
        "step_order": step_order,
        "step_category": step_category,
        "step_name": step_name,
        "tool_name": tool_name,
        "tool_version": tool_version or None,
        "tool_uri": tool_uri,
        "parameters": parameters,
        "config": config,
        "output_views": list(output_views) if output_views is not None else None,
    }


class BaseOrchestrator:
    """Base class for converter orchestrators.

    Provides shared methods for writing ontology, provenance, and dataset
    Parquet files. Subclasses implement tool-specific conversion logic.
    """

    # Default compression; subclasses may override via __init__.
    _compression: str = "zstd"

    def _write_ontology(
        self,
        output_folder: Path,
        prefix: str,
        ontology_entries: list[dict],
    ) -> Path | None:
        """Write ontology.parquet for the views already written. Returns path if written, else None."""
        from qpx.converters.ontology_util import dedupe_ontology_entries, reconcile_ontology_entries
        from qpx.writers.ontology import OntologyWriter

        deduped = dedupe_ontology_entries(ontology_entries)
        dropped = len(ontology_entries) - len(deduped)
        entries = reconcile_ontology_entries(deduped, output_folder, prefix)
        if not entries:
            return None
        onto_path = output_folder / f"{prefix}.ontology.parquet"
        with OntologyWriter(onto_path, creator="qpx", compression=self._compression) as writer:
            writer.write_batch(entries)
        if dropped:
            logger.info("Collapsed %d duplicate (field_name, view) ontology entries", dropped)
        logger.info("Wrote %d ontology entries to %s", len(entries), onto_path)
        return onto_path

    def _write_provenance(
        self,
        output_folder: Path,
        prefix: str,
        records: list[dict],
    ) -> Path | None:
        """Write provenance.parquet. Returns path if written, else None."""
        if not records:
            return None
        from qpx.writers.provenance import ProvenanceWriter

        prov_path = output_folder / f"{prefix}.provenance.parquet"
        with ProvenanceWriter(prov_path, creator="qpx", compression=self._compression) as writer:
            writer.write_batch(records)
        logger.info("Wrote %d provenance steps to %s", len(records), prov_path)
        return prov_path

    def _write_dataset(
        self,
        output_folder: Path,
        prefix: str,
        project_accession: str | None,
        *,
        software_name: str = "unknown",
        software_version: str | None = None,
        provenance_records: list[dict] | None = None,
    ) -> Path:
        """Write dataset.parquet with project-level metadata."""
        from qpx.writers.dataset import DatasetWriter

        record = build_dataset_record(
            project_accession=project_accession,
            software_name=software_name,
            software_version=software_version,
            provenance_records=provenance_records,
        )
        ds_path = output_folder / f"{prefix}.dataset.parquet"
        with DatasetWriter(ds_path, creator="qpx", compression=self._compression) as writer:
            writer.write_batch([record])
        logger.info("Wrote dataset metadata to %s", ds_path)
        return ds_path

    def _write_mudata(self, output_folder: Path, prefix: str) -> Path | None:
        """
        Assemble and write the MuData (.h5mu) view of the converted dataset.

        Brings QuantMS/OpenMS QPX to parity with the DIA-NN (quantmsdiann) path,
        which emits muData. Must run after the core + metadata parquet are
        written. Best-effort: expected dependency, database, validation, and I/O
        failures are logged and skipped rather than failing the run.

        The implementation lives in ``qpx.mudata`` so the transform commands
        that rewrite those parquet can refresh the view the same way.
        """
        from qpx.mudata import write_dataset_mudata

        return write_dataset_mudata(output_folder, prefix)
