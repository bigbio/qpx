"""Convert CLI maxquant implementation."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from .options import (
    _compression_option,
    _dataset_metadata_options,
    _max_cpus_option,
    _max_memory_option,
    _output_folder_option,
    _output_prefix_option,
    _project_accession_option,
    _verbose_option,
)
from .postprocess import _finish_conversion


def _maxquant_structures(structures, msms_file, evidence_file, protein_groups_file) -> set[str]:
    """Resolve requested MaxQuant views from options or available inputs."""
    # Determine which structures to produce
    if structures:
        requested = {s.strip() for s in structures.split(",")}
    else:
        # Auto-detect from provided files
        requested = set()
        if msms_file:
            requested.add("psm")
        if evidence_file:
            requested.add("feature")
        if protein_groups_file:
            requested.add("pg")
    if not requested:
        raise click.ClickException(
            "No input files provided. Supply at least one of --msms-file, --evidence-file, or --protein-groups-file."
        )

    return requested


@click.command("maxquant")
@click.option(
    "--msms-file",
    help="MaxQuant msms.txt file (for PSM conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--evidence-file",
    help="MaxQuant evidence.txt file (for feature conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--protein-groups-file",
    help="MaxQuant proteinGroups.txt file (for PG conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file (required for feature and PG; its fixed modifications are added to PSM and feature peptidoforms)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@_output_folder_option
@_output_prefix_option
@click.option(
    "--structures",
    help="Comma-separated list of structures to produce (psm, feature, pg). Default: all available.",
    default=None,
)
@_max_cpus_option
@_max_memory_option
@click.option(
    "--spectral-data",
    help="Include spectral data fields in PSM output",
    is_flag=True,
)
@_project_accession_option
@_compression_option
@_verbose_option
@_dataset_metadata_options
def convert_maxquant_cmd(**kwargs):
    """Convert MaxQuant output to QPX format.

    Reads MaxQuant result files (msms.txt, evidence.txt, proteinGroups.txt)
    and writes corresponding QPX Parquet files.

    \b
    Examples:
        # Convert everything
        qpxc convert maxquant \\
            --msms-file msms.txt \\
            --evidence-file evidence.txt \\
            --protein-groups-file proteinGroups.txt \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./qpx_output

        # Convert PSMs only
        qpxc convert maxquant \\
            --msms-file msms.txt \\
            --output-folder ./qpx_output \\
            --structures psm
    """
    msms_file = kwargs["msms_file"]
    evidence_file = kwargs["evidence_file"]
    protein_groups_file = kwargs["protein_groups_file"]
    sdrf_file = kwargs["sdrf_file"]
    output_folder = kwargs["output_folder"]
    output_prefix = kwargs["output_prefix"]
    structures = kwargs["structures"]
    max_cpus = kwargs["max_cpus"]
    max_memory = kwargs["max_memory"]
    spectral_data = kwargs["spectral_data"]
    project_accession = kwargs["project_accession"]
    compression = kwargs["compression"]
    verbose = kwargs["verbose"]
    metadata = kwargs["metadata"]

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder.mkdir(parents=True, exist_ok=True)

    requested = _maxquant_structures(structures, msms_file, evidence_file, protein_groups_file)

    from qpx.converters.maxquant import MaxQuantConverter

    converter = MaxQuantConverter(max_memory=max_memory, max_cpus=max_cpus, compression=compression)
    converter.convert(
        msms_file=msms_file,
        evidence_file=evidence_file,
        protein_groups_file=protein_groups_file,
        sdrf_file=sdrf_file,
        output_folder=output_folder,
        output_prefix=output_prefix,
        structures=list(requested),
        spectral_data=spectral_data,
        project_accession=project_accession,
    )

    prefix = output_prefix or "maxquant"
    _finish_conversion(output_folder, prefix, project_accession, metadata)

    click.echo(f"MaxQuant conversion complete. Output: {output_folder}")
