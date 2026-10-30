"""Convert CLI cdap implementation."""

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
    _verbose_option,
)
from .postprocess import _finish_conversion


@click.command("cdap")
@click.option(
    "--psm-dir",
    help="Directory containing CDAP *.psm files for one study",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@_output_folder_option
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
    default="cdap",
    show_default=True,
)
@click.option(
    "--structures",
    help="Comma-separated list of structures to produce (psm, feature, pg). Default: all.",
    default=None,
)
@click.option(
    "--project-accession",
    help="PDC / ProteomeXchange accession (e.g. PDC000227)",
)
@_compression_option
@_max_memory_option
@_max_cpus_option
@convert_resource_options
@_verbose_option
@_dataset_metadata_options
def convert_cdap_cmd(**kwargs):
    """Convert CPTAC CDAP .psm files to QPX format.

    Reads all *.psm tab-separated files in a study directory (e.g. PDC000440/)
    and writes psm.parquet, feature.parquet, and pg.parquet.

    \b
    Examples:
        # Convert a single CPTAC study
        qpxc convert cdap \\
            --psm-dir /data/CPTAC/PDC000440 \\
            --output-folder ./qpx_output/PDC000440 \\
            --project-accession PDC000440

        # Only produce PSM and Feature (skip PG)
        qpxc convert cdap \\
            --psm-dir /data/CPTAC/PDC000440 \\
            --output-folder ./qpx_output \\
            --structures psm,feature
    """
    psm_dir = kwargs["psm_dir"]
    output_folder = kwargs["output_folder"]
    output_prefix = kwargs["output_prefix"]
    structures = kwargs["structures"]
    project_accession = kwargs["project_accession"]
    compression = kwargs["compression"]
    max_memory = kwargs["max_memory"]
    max_cpus = kwargs["max_cpus"]
    verbose = kwargs["verbose"]
    metadata = kwargs["metadata"]

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder.mkdir(parents=True, exist_ok=True)

    requested = None
    if structures:
        requested = [s.strip() for s in structures.split(",")]

    from qpx.converters.cdap import CdapConverter

    converter = CdapConverter(
        max_memory=max_memory,
        max_cpus=max_cpus,
        compression=compression,
    )
    converter.convert(
        psm_dir=psm_dir,
        output_folder=output_folder,
        output_prefix=output_prefix,
        structures=requested,
        project_accession=project_accession,
    )

    _finish_conversion(output_folder, output_prefix, project_accession, metadata)

    click.echo(f"CDAP conversion complete. Output: {output_folder}")
