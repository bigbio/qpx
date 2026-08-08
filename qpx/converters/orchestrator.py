"""BaseOrchestrator — shared ontology, provenance, and dataset writing for converters.

Orchestrators (MaxQuantConverter, OpenMSConverter, etc.) compose adapters and
write ontology.parquet, provenance.parquet, and dataset.parquet. This base class
consolidates the common writing logic.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

import duckdb

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
        "creation_date": datetime.now().isoformat(),
        "file_checksums": None,
        "file_row_counts": None,
        "file_sizes_bytes": None,
        "total_structures": None,
        "packaged_at": None,
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
        """Write combined ontology.parquet. Returns path if written, else None."""
        if not ontology_entries:
            return None
        from qpx.converters.ontology_util import dedupe_ontology_entries
        from qpx.writers.ontology import OntologyWriter

        entries = dedupe_ontology_entries(ontology_entries)
        dropped = len(ontology_entries) - len(entries)
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
        """
        h5mu_path = output_folder / f"{prefix}.h5mu"
        h5mu_tmp = h5mu_path.with_name(f"{h5mu_path.stem}.tmp{h5mu_path.suffix}")
        try:
            # Core parquet files have already been refreshed, so an older MuData
            # view would be stale if this build fails.
            h5mu_path.unlink(missing_ok=True)

            from qpx.dataset import Dataset
            from qpx.mudata import build_mudata

            dataset = Dataset(str(output_folder), file_prefix=prefix)
            try:
                required_modalities = {
                    name
                    for name, structure in (
                        ("precursors", dataset.feature),
                        ("proteins", dataset.pg),
                    )
                    if structure is not None
                }
                if not required_modalities:
                    logger.info("Skipping muData for %s: no feature or pg quantification data", prefix)
                    return None

                modalities = required_modalities | {"expression", "differential"}
                mdata = build_mudata(
                    dataset,
                    modalities=sorted(modalities),
                    all_intensity_labels=True,
                )
                missing = required_modalities - set(mdata.mod)
                if missing:
                    raise ValueError(f"Missing required quantification modalities: {', '.join(sorted(missing))}")
                mdata.write(str(h5mu_tmp))
                h5mu_tmp.replace(h5mu_path)
            finally:
                dataset.close()
        except (ImportError, OSError, RuntimeError, TypeError, ValueError, duckdb.Error) as exc:
            logger.warning("Could not build muData for %s: %s", prefix, exc)
            return None
        finally:
            try:
                h5mu_tmp.unlink(missing_ok=True)
            except OSError as exc:
                logger.warning("Could not remove temporary muData file %s: %s", h5mu_tmp, exc)
        logger.info("Wrote muData to %s", h5mu_path)
        return h5mu_path
