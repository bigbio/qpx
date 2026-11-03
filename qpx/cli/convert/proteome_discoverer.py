"""Convert CLI proteome-discoverer implementation."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from qpx.cli.resources import convert_resource_options

from .options import (
    _dataset_metadata_options,
    _max_cpus_option,
    _max_memory_option,
    _output_folder_option,
    _output_prefix_option,
    _project_accession_option,
    _verbose_option,
)
from .postprocess import _finish_conversion


@click.command("proteome-discoverer")
@click.option(
    "--psm-file",
    help="Proteome Discoverer PSMs.txt export",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--proteins-file",
    help="Proteome Discoverer Proteins.txt export (for PG conversion and protein group members)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--input-files",
    help="Proteome Discoverer InputFiles.txt export; needed when the PSM export has no Spectrum File column",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file (for sample/run conversion); its comment[label] values must match the reporter channels",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option("--rank1-only", is_flag=True, help="Keep only the first-ranked hit of each spectrum in the PSM view")
@_output_folder_option
@_output_prefix_option
@_max_memory_option
@_max_cpus_option
@convert_resource_options
@_project_accession_option
@_verbose_option
@_dataset_metadata_options
def convert_proteome_discoverer_cmd(**kwargs):
    """Convert Proteome Discoverer text exports to QPX format.

    Reads the exported PSMs.txt and Proteins.txt tables, not the .pdResult file.

    Examples:

    \b
        # PSMs, features and protein groups of a TMT study
        qpxc convert proteome-discoverer \\
            --psm-file study_PSMs.txt \\
            --proteins-file study_Proteins.txt \\
            --input-files study_InputFiles.txt \\
            --sdrf-file study.sdrf.tsv \\
            --output-folder ./qpx_output
    """
    if kwargs["verbose"]:
        logging.getLogger().setLevel(logging.DEBUG)
    output_folder = kwargs["output_folder"]
    output_prefix = kwargs["output_prefix"]
    project_accession = kwargs["project_accession"]
    output_folder.mkdir(parents=True, exist_ok=True)

    from qpx.converters.proteome_discoverer import ProteomeDiscovererConverter

    converter = ProteomeDiscovererConverter(output_folder, max_memory=kwargs["max_memory"], max_cpus=kwargs["max_cpus"])
    converter.convert(
        kwargs["psm_file"],
        kwargs["proteins_file"],
        input_files=kwargs["input_files"],
        sdrf_file=kwargs["sdrf_file"],
        output_prefix=output_prefix,
        rank1_only=kwargs["rank1_only"],
        project_accession=project_accession,
    )
    prefix = output_prefix or "proteome_discoverer"
    _finish_conversion(output_folder, prefix, project_accession, kwargs["metadata"])
    click.echo(f"Proteome Discoverer conversion complete. Output: {output_folder}")
