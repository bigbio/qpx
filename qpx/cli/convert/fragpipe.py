"""Convert CLI fragpipe implementation."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from qpx.cli.resources import convert_resource_options

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


@click.command("fragpipe")
@click.option(
    "--psm-file",
    help="FragPipe psm.tsv file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--ion-file",
    help="FragPipe combined_ion.tsv file, or a plex's ion.tsv with --annotation-file (for feature conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--pg-file",
    help="FragPipe combined_protein.tsv file, or a plex's protein.tsv with --annotation-file (for PG conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file (for sample/run conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--experiment-annotation-file",
    help=(
        "FragPipe experiment_annotation.tsv mapping experiments to member raw files; "
        "--ion-file and --pg-file need it unless the SDRF source names are the experiments"
    ),
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--annotation-file",
    help=(
        "annotation.txt of one TMT/iTRAQ plex, whose lines map a channel to its sample column; "
        "reads the reporter quantities of --ion-file and --pg-file, and needs --psm-file for the plex's raw runs"
    ),
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--tmt-report",
    help="TMT-Integrator abundance_protein_*.tsv, added to the pg view as the log2 abundance of each channel",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@_output_folder_option
@_output_prefix_option
@_max_memory_option
@_max_cpus_option
@convert_resource_options
@_project_accession_option
@_compression_option
@_verbose_option
@_dataset_metadata_options
def convert_fragpipe_cmd(**kwargs):
    """Convert FragPipe output to QPX format.

    Reads FragPipe result files and converts them into QPX Parquet format.
    Supports psm.tsv, combined_ion.tsv, and combined_protein.tsv.

    \b
    Examples:
        # Convert PSMs only
        qpxc convert fragpipe \\
            --psm-file psm.tsv \\
            --output-folder ./qpx_output

    \b
        # Convert features + protein groups
        qpxc convert fragpipe \\
            --ion-file combined_ion.tsv \\
            --pg-file combined_protein.tsv \\
            --sdrf-file metadata.sdrf.tsv \\
            --experiment-annotation-file experiment_annotation.tsv \\
            --output-folder ./qpx_output

    \b
        # One TMT plex: convert each plex with its own prefix
        qpxc convert fragpipe \\
            --psm-file plex1/psm.tsv \\
            --ion-file plex1/ion.tsv \\
            --pg-file plex1/protein.tsv \\
            --annotation-file plex1/annotation.txt \\
            --sdrf-file study.sdrf.tsv \\
            --output-folder ./qpx_output \\
            --output-prefix plex1
    """
    psm_file = kwargs["psm_file"]
    ion_file = kwargs["ion_file"]
    pg_file = kwargs["pg_file"]
    sdrf_file = kwargs["sdrf_file"]
    experiment_annotation_file = kwargs["experiment_annotation_file"]
    output_folder = kwargs["output_folder"]
    output_prefix = kwargs["output_prefix"]
    max_memory = kwargs["max_memory"]
    max_cpus = kwargs["max_cpus"]
    project_accession = kwargs["project_accession"]
    compression = kwargs["compression"]
    verbose = kwargs["verbose"]
    metadata = kwargs["metadata"]

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if not any([psm_file, ion_file, pg_file]):
        raise click.ClickException("No input files provided. Supply at least one of --psm-file, --ion-file, or --pg-file.")

    output_folder.mkdir(parents=True, exist_ok=True)

    from qpx.converters.fragpipe import FragPipeConverter

    converter = FragPipeConverter(
        output_directory=output_folder, max_memory=max_memory, max_cpus=max_cpus, compression=compression
    )
    if kwargs["annotation_file"]:
        converter.load_plex(kwargs["annotation_file"], kwargs["tmt_report"])
    converter.convert(
        psm_file=psm_file,
        ion_file=ion_file,
        pg_file=pg_file,
        sdrf_file=sdrf_file,
        experiment_annotation_file=experiment_annotation_file,
        output_prefix=output_prefix,
        project_accession=project_accession,
    )

    prefix = output_prefix or "fragpipe"
    _finish_conversion(output_folder, prefix, project_accession, metadata)

    click.echo(f"FragPipe conversion complete. Output: {output_folder}")
