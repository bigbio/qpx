"""Convert CLI spectronaut implementation."""

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


@click.command("spectronaut")
@click.option(
    "--report-path",
    help="Spectronaut report file path (TSV or Parquet)",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--qvalue-threshold",
    help=(
        "Optional q-value threshold. Unset (the default) converts the report as "
        "reported — no filtering. When given, the feature view filters on the "
        "precursor q-value."
    ),
    default=None,
    type=float,
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
def convert_spectronaut_cmd(**kwargs):
    r"""Convert Spectronaut report to QPX format.

    Reads a Spectronaut report TSV file and converts feature-level and
    protein-group quantification data into QPX Parquet format.

    \b
    Examples:
        # Feature + protein groups
        qpxc convert spectronaut \
            --report-path report.tsv \
            --output-folder ./qpx_output

        # With SDRF metadata
        qpxc convert spectronaut \
            --report-path report.tsv \
            --sdrf-file data.sdrf.tsv \
            --output-folder ./qpx_output
    """
    report_path = kwargs["report_path"]
    sdrf_file = kwargs["sdrf_file"]
    qvalue_threshold = kwargs["qvalue_threshold"]
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

    output_folder.mkdir(parents=True, exist_ok=True)
    prefix = output_prefix or "spectronaut"

    from qpx.converters.spectronaut import SpectronautConverter

    converter = SpectronautConverter(
        report_path=report_path,
        sdrf_path=sdrf_file,
        max_memory=max_memory,
        max_cpus=max_cpus,
        compression=compression,
    )
    converter.convert_features(
        qvalue_threshold=qvalue_threshold,
        output_folder=output_folder,
        output_prefix=prefix,
    )
    converter.convert_pg(
        output_folder=output_folder,
        output_prefix=prefix,
    )

    if sdrf_file:
        converter.convert_sdrf(output_folder=output_folder, prefix=prefix)
    converter.write_ontology(output_folder, prefix=prefix)
    converter.write_provenance(output_folder, prefix=prefix)
    converter.write_dataset(output_folder, prefix=prefix, project_accession=project_accession)

    _finish_conversion(output_folder, prefix, project_accession, metadata)

    click.echo(f"Spectronaut conversion complete. Output: {output_folder}")
