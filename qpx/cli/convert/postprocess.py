"""Convert CLI postprocess implementation."""

from __future__ import annotations

import logging
from http.client import HTTPException
from pathlib import Path
from typing import Optional

import click

from .options import _DatasetMetadata

logger = logging.getLogger("qpx.cli.convert")


def _warn(message: str) -> None:
    """Print a warning to stderr with the shared ``WARNING:`` prefix."""
    click.echo(f"WARNING: {message}", err=True)


def _log_summary(output_folder) -> None:
    """Log a conversion summary for a just-written output folder.

    Shared by every convert subcommand. Guarded so a summary failure NEVER fails
    an otherwise-successful conversion (a warning is logged instead).
    """
    from qpx.converters.summary import log_conversion_summary

    log_conversion_summary(output_folder, logger=logger)


def _write_mudata(output_folder: Path, prefix: str, enabled: bool) -> None:
    """Refresh or remove the dataset's MuData view after a conversion.

    Uses qpx's own writer, which refuses a MuData missing a required
    quantification modality rather than writing a partial view, and is
    best-effort: the parquet views are the dataset's source of truth, so a view
    that cannot be built is reported and the conversion still succeeds.
    """
    if not enabled:
        try:
            (Path(output_folder) / f"{prefix}.h5mu").unlink(missing_ok=True)
        except OSError as exc:
            logger.warning("Could not remove stale MuData for %s: %s", prefix, exc)
        return
    from qpx.mudata import write_dataset_mudata

    written = write_dataset_mudata(Path(output_folder), prefix)
    if written is None:
        _warn("no MuData view was written (see the log); the Parquet views are complete")
    else:
        click.echo(f"MuData view: {written.name}")


def _annotate_protein_properties(output_folder: Path, fasta: Optional[Path], prefix: Optional[str] = None) -> None:
    """Fill null protein properties from an optional FASTA after a conversion.

    The conversion has already succeeded, so a dataset this step cannot handle
    (e.g. partitioned views) is reported and left as converted rather than failed.
    """
    if fasta is None:
        return
    from qpx.cli.transform import annotate_dataset_protein_properties

    click.echo(f"Filling protein properties from {fasta.name}")
    try:
        annotate_dataset_protein_properties(Path(output_folder), fasta, in_place=True, prefix=prefix)
    except click.ClickException as exc:
        _warn(f"protein properties from FASTA skipped: {exc.format_message()}")


def _normalize_accessions(output_folder: Path, prefix: str) -> None:
    """Store every protein as its bare accession (``sp|P04114|APOB_HUMAN`` -> ``P04114``); entry
    names stay in ``pg_names``."""
    from qpx.transforms.accession_normalizer import normalize_dataset_accessions

    rewritten = normalize_dataset_accessions(Path(output_folder), prefix)
    if rewritten:
        click.echo(f"Protein accessions stored as bare accessions in {rewritten} view file(s)")


def _maybe_enrich_pride(output_folder, project_accession: str | None, enrich: bool, prefix: str) -> None:
    """Optionally enrich a converted dataset with PRIDE metadata."""
    if not enrich:
        return
    if not project_accession:
        _warn("--enrich-pride requires --project-accession; skipping enrichment.")
        return
    try:
        from qpx.dataset import Dataset

        with Dataset(str(output_folder), file_prefix=prefix) as ds:
            ds.enrich_from_pride(project_accession)
        click.echo(f"Enriched dataset with PRIDE metadata for {project_accession}")
    except (OSError, HTTPException, ValueError) as exc:
        _warn(f"Could not enrich from PRIDE: {exc}")


def _apply_software(
    output_folder,
    prefix: str,
    software_name: str | None,
    software_version: str | None,
    config_json: str | None,
) -> None:
    """Record the workflow that produced the dataset in dataset and provenance parquet.

    Sets the software fields of dataset.parquet and prepends a 'workflow' step to
    provenance.parquet. Skipped with a warning when either metadata file is missing.
    """
    if not software_name:
        return
    import pandas as pd

    from qpx.converters.orchestrator import build_provenance_step
    from qpx.dataset import Dataset

    with Dataset(str(output_folder), file_prefix=prefix) as ds:
        if ds.dataset_meta is None or ds.provenance is None:
            _warn("no dataset/provenance metadata found; --software-name skipped")
            return
        record = ds.dataset_meta.to_df().iloc[0].to_dict()
        record["software_name"] = software_name
        record["software_version"] = software_version
        ds.save_structure([record], "dataset", prefix=prefix)

        workflow = build_provenance_step(
            step_category="workflow",
            step_name="pipeline_execution",
            tool_name=software_name,
            tool_version=software_version,
            config=config_json,
        )
        frame = pd.concat([pd.DataFrame([workflow]), ds.provenance.to_df()], ignore_index=True)
        frame["step_order"] = range(1, len(frame) + 1)
        ds.save_structure(frame, "provenance", prefix=prefix)


def _maybe_compute_integrity(output_folder, prefix: str, enabled: bool) -> None:
    """Record checksums, row counts and sizes of the output folder in dataset parquet.

    Call this after the last file write of a conversion, so the record covers every Parquet and h5ad file.
    """
    if not enabled:
        return
    from qpx.dataset import Dataset

    with Dataset(str(output_folder), file_prefix=prefix) as ds:
        if ds.dataset_meta is None:
            _warn("no dataset.parquet found; --compute-integrity skipped")
            return
        integrity = ds.compute_integrity()
        record = ds.dataset_meta.to_df().iloc[0].to_dict()
        record.update(integrity)
        ds.save_structure([record], "dataset", prefix=prefix)
    click.echo(f"Integrity recorded for {integrity['total_structures']} files")


def _finish_conversion(
    output_folder: Path,
    prefix: str,
    project_accession: str | None,
    metadata: _DatasetMetadata,
    *,
    fasta: Path | None = None,
    mudata: bool | None = None,
    written: bool = True,
) -> None:
    """Run the post-conversion steps shared by the convert commands, in order.

    Dataset-changing steps run only when ``written`` is true. The MuData view is
    refreshed when ``mudata`` is not None. The integrity record runs last, so it covers every
    Parquet and h5ad file (not the ``.h5mu`` view).
    """
    if written:
        _annotate_protein_properties(output_folder, fasta, prefix)
        _normalize_accessions(output_folder, prefix)
        _apply_software(output_folder, prefix, metadata.software_name, metadata.software_version, metadata.config_json)
        _maybe_enrich_pride(output_folder, project_accession, metadata.enrich_pride, prefix)
    if mudata is not None:
        _write_mudata(output_folder, prefix, mudata and written)
    if written:
        _maybe_compute_integrity(output_folder, prefix, metadata.compute_integrity)
    _log_summary(output_folder)
