"""Convert CLI postprocess implementation."""

from __future__ import annotations

import hashlib
import logging
import shutil
from http.client import HTTPException
from pathlib import Path
from typing import Optional

import click

from .options import _DatasetMetadata

logger = logging.getLogger("qpx.cli.convert")

# Options that say where or how fast a conversion ran, or that other fields already record.
_UNRECORDED_OPTIONS = frozenset(
    {
        "output_folder",
        "output_prefix",
        "verbose",
        "max_memory",
        "max_cpus",
        "temp_dir",
        "compression",
        "project_accession",
        "software_name",
        "software_version",
        "software_config",
        "enrich_pride",
        "compute_integrity",
    }
)
# Options whose unset value is itself a decision about the data, so it is recorded too.
_RECORDED_WHEN_UNSET = frozenset({"qvalue_threshold"})
_SDRF_OPTIONS = ("sdrf_file", "sdrf_path")
_SHA256_CHUNK_BYTES = 1024 * 1024


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


def _file_sha256(path: Path) -> str:
    """Return the SHA-256 hex digest of a file, read in chunks."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(_SHA256_CHUNK_BYTES):
            digest.update(chunk)
    return digest.hexdigest()


def _path_option_parameters(name: str, values: list) -> list[dict]:
    """Describe the input paths of one option by file name and, for files, SHA-256."""
    paths = [Path(value) for value in values]
    parameters = [{"key": name, "value": "; ".join(path.name for path in paths)}]
    digests = [_file_sha256(path) for path in paths if path.is_file()]
    if digests:
        parameters.append({"key": f"{name}_sha256", "value": "; ".join(digests)})
    return parameters


def _effective_parameters(ctx: click.Context) -> list[dict]:
    """Return the provenance parameters of the options a command ran with, defaults included."""
    parameters: list[dict] = []
    for option in ctx.command.params:
        name = option.name
        value = ctx.params.get(name)
        if name in _UNRECORDED_OPTIONS or (value in (None, ()) and name not in _RECORDED_WHEN_UNSET):
            continue
        values = list(value) if isinstance(value, (list, tuple)) else [value]
        if isinstance(option.type, click.Path):
            parameters.extend(_path_option_parameters(name, values))
        else:
            rendered = "none" if value is None else "; ".join(str(item) for item in values)
            parameters.append({"key": name, "value": rendered})
    return parameters


def _record_conversion_options(output_folder: Path, prefix: str) -> None:
    """Add the effective options and input checksums of this command to its ``format_conversion`` step."""
    from qpx._version import __version__
    from qpx.converters.orchestrator import build_provenance_step
    from qpx.dataset import Dataset

    ctx = click.get_current_context()
    with Dataset(str(output_folder), file_prefix=prefix) as ds:
        steps = ds.provenance.to_arrow().to_pylist() if ds.provenance is not None else []
        conversion = next((step for step in reversed(steps) if step["step_category"] == "format_conversion"), None)
        if conversion is None:
            conversion = build_provenance_step(
                step_category="format_conversion",
                step_name="qpx_conversion",
                tool_name="qpx",
                tool_version=__version__,
                step_order=len(steps) + 1,
            )
            steps.append(conversion)
        recorded = {parameter["key"]: parameter for parameter in conversion["parameters"] or []}
        recorded.update({parameter["key"]: parameter for parameter in _effective_parameters(ctx)})
        conversion["parameters"] = list(recorded.values())
        ds.save_structure(steps, "provenance", prefix=prefix)


def _preserve_sdrf(output_folder: Path, prefix: str) -> None:
    """Keep the SDRF this command was given next to the views as ``{prefix}.sdrf.tsv``."""
    params = click.get_current_context().params
    source = next((params[name] for name in _SDRF_OPTIONS if params.get(name) is not None), None)
    if source is None:
        return
    target = Path(output_folder) / f"{prefix}.sdrf.tsv"
    if not (target.exists() and target.samefile(source)):
        shutil.copyfile(source, target)


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
        _preserve_sdrf(output_folder, prefix)
        _record_conversion_options(output_folder, prefix)
        _annotate_protein_properties(output_folder, fasta, prefix)
        _apply_software(output_folder, prefix, metadata.software_name, metadata.software_version, metadata.config_json)
        _maybe_enrich_pride(output_folder, project_accession, metadata.enrich_pride, prefix)
    if mudata is not None:
        _write_mudata(output_folder, prefix, mudata and written)
    if written:
        _maybe_compute_integrity(output_folder, prefix, metadata.compute_integrity)
    _log_summary(output_folder)
