"""Convert CLI quantms msstats implementation."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from qpx.cli.resources import convert_resource_options

from .options import (
    _dataset_metadata_options,
    _output_folder_option,
    _output_prefix_option,
    _sdrf_file_option,
    _verbose_option,
)
from .postprocess import _finish_conversion


@click.command("quantms-msstats")
@click.option(
    "--msstats-file",
    help="QuantMS-generated *_msstats_in.csv file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@_sdrf_file_option
@_output_folder_option
@_output_prefix_option
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD007683)",
)
@click.option(
    "--max-memory",
    help="Maximum DuckDB memory limit (e.g. '16GB')",
)
@click.option(
    "--max-cpus",
    help="Maximum number of DuckDB threads",
    type=int,
)
@convert_resource_options
@_verbose_option
@_dataset_metadata_options
def convert_quantms_msstats_cmd(**kwargs):
    """Convert QuantMS MSstats and SDRF inputs to QPX.

    Produces Feature, Sample, Run, Dataset, Provenance, and Ontology Parquet
    files. MSstats does not contain the evidence required to produce PSM or
    protein-group views, so those views are intentionally not generated.

    \b
    Example:
        qpxc convert quantms-msstats \
            --msstats-file PXD007683.sdrf_openms_design_msstats_in.csv \
            --sdrf-file PXD007683.sdrf.tsv \
            --output-folder ./qpx_output
    """
    msstats_file = kwargs["msstats_file"]
    sdrf_file = kwargs["sdrf_file"]
    output_folder = kwargs["output_folder"]
    output_prefix = kwargs["output_prefix"]
    project_accession = kwargs["project_accession"]
    max_memory = kwargs["max_memory"]
    max_cpus = kwargs["max_cpus"]
    verbose = kwargs["verbose"]
    metadata = kwargs["metadata"]

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    from qpx.converters.quantms_msstats import QuantmsMsstatsConverter

    converter = QuantmsMsstatsConverter(
        max_memory=max_memory,
        max_cpus=max_cpus,
    )
    prefix = converter.convert(
        msstats_file=msstats_file,
        sdrf_file=sdrf_file,
        output_folder=output_folder,
        output_prefix=output_prefix,
        project_accession=project_accession,
    )

    _finish_conversion(output_folder, prefix, project_accession, metadata)

    click.echo(f"QuantMS MSstats conversion complete. Output prefix: {prefix}")
