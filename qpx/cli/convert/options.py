"""Convert CLI options implementation."""

from __future__ import annotations

import functools
import json
from dataclasses import dataclass
from pathlib import Path

import click


@dataclass(frozen=True)
class _DatasetMetadata:
    """Validated dataset-metadata options of one convert command.

    ``config_json`` holds the workflow parameters as compact JSON, or None when no config file is given.
    """

    software_name: str | None
    software_version: str | None
    config_json: str | None
    enrich_pride: bool
    compute_integrity: bool


def _verbose_option(func):
    """Add the shared verbose flag."""
    return click.option("--verbose", help="Enable verbose logging", is_flag=True)(func)


def _output_folder_option(func):
    """Add the shared QPX output directory option."""
    return click.option(
        "--output-folder",
        help="Output directory for generated QPX files",
        required=True,
        type=click.Path(file_okay=False, path_type=Path),
    )(func)


def _output_prefix_option(func):
    """Add the shared optional output prefix."""
    return click.option("--output-prefix", help="Prefix for output file names")(func)


def _compression_option(func):
    """Add the shared Parquet compression option."""
    return click.option(
        "--compression",
        type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
        default="zstd",
        show_default=True,
        help="Parquet compression codec.",
    )(func)


def _sdrf_file_option(func):
    """Add the shared required SDRF file option."""
    return click.option(
        "--sdrf-file",
        help="SDRF metadata file path",
        required=True,
        type=click.Path(exists=True, dir_okay=False, path_type=Path),
    )(func)


def _max_memory_option(func):
    """Add the shared memory limit option."""
    return click.option(
        "--max-memory",
        help="DuckDB memory limit (e.g. '16GB'); the converter's Python row batches come on top",
    )(func)


def _max_cpus_option(func):
    """Add the shared thread limit option."""
    return click.option(
        "--max-cpus",
        help="Maximum number of threads",
        type=int,
    )(func)


def _project_accession_option(func):
    """Add the shared project accession option."""
    return click.option(
        "--project-accession",
        help="PRIDE / ProteomeXchange accession (e.g. PXD020192)",
    )(func)


def _openms_project_accession_option(func):
    """Add the project accession option shared by OpenMS commands."""
    return click.option(
        "--project-accession",
        help="PRIDE / ProteomeXchange accession (e.g. PXD001819)",
    )(func)


def _fasta_option(func):
    """Add the shared optional protein FASTA option."""
    return click.option(
        "--fasta",
        help=(
            "Optional protein FASTA used for the search. When given, fills null pg.sequence_coverage, "
            "pg.molecular_weight and feature.pg_positions for target rows after conversion; proteins "
            "absent from it (e.g. DIA-NN internal decoys) stay null. See 'qpxc transform protein-properties'."
        ),
        default=None,
        type=click.Path(exists=True, dir_okay=False, path_type=Path),
    )(func)


def _mudata_option(func):
    """Add the shared MuData output flag."""
    return click.option(
        "--mudata/--no-mudata",
        default=True,
        show_default=True,
        help=(
            "Write the dataset's MuData (.h5mu) view after conversion. Best-effort: a view that "
            "cannot be built is reported and the conversion still succeeds. Use --no-mudata to skip "
            "the build on very large datasets."
        ),
    )(func)


def _dataset_metadata_options(func):
    """Add the dataset-metadata options shared by the convert commands.

    The options set the workflow name and version, enrich from PRIDE and record
    file integrity. They are listed in the order declared here. The command receives
    them as one ``metadata`` argument.
    """
    options = (
        click.option(
            "--software-name",
            default=None,
            help=(
                "Workflow that produced the input (e.g. nf-core/mhcquant). "
                "Sets dataset.software_name/software_version and adds a 'workflow' provenance step."
            ),
        ),
        click.option("--software-version", default=None, help="Version of --software-name."),
        click.option(
            "--software-config",
            default=None,
            type=click.Path(exists=True, dir_okay=False, path_type=Path),
            help=(
                "JSON file with the workflow parameters (e.g. Nextflow params.json), "
                "stored in the workflow provenance step's config."
            ),
        ),
        click.option(
            "--enrich-pride",
            is_flag=True,
            default=False,
            help=(
                "Fetch project title, description and PubMed ID from the PRIDE API after "
                "conversion (needs --project-accession; best-effort)."
            ),
        ),
        click.option(
            "--compute-integrity",
            is_flag=True,
            default=False,
            help="Write file checksums, row counts and sizes for the output folder into dataset.parquet, as the last step.",
        ),
    )

    @functools.wraps(func)
    def wrapper(*, software_name, software_version, software_config, enrich_pride, compute_integrity, **kwargs):
        """Validate the workflow options, then call the command with them as ``metadata``."""
        metadata = _DatasetMetadata(
            software_name=software_name,
            software_version=software_version,
            config_json=_software_config_json(software_name, software_version, software_config),
            enrich_pride=enrich_pride,
            compute_integrity=compute_integrity,
        )
        return func(metadata=metadata, **kwargs)

    for option in reversed(options):
        wrapper = option(wrapper)
    return wrapper


def _software_config_json(software_name, software_version, software_config) -> str | None:
    """Validate the workflow options and return the parameters as compact, sorted JSON.

    Runs before the conversion, so an invalid option fails without writing output.
    """
    if (software_version is not None or software_config is not None) and not software_name:
        raise click.UsageError("--software-version and --software-config require --software-name.")
    if software_config is None:
        return None
    try:
        with open(software_config, encoding="utf-8") as handle:
            config = json.load(handle)
    except (OSError, ValueError) as exc:
        raise click.UsageError(f"--software-config {software_config} is not a readable JSON file: {exc}") from exc
    return json.dumps(config, separators=(",", ":"), sort_keys=True)
