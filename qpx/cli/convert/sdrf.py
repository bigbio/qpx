"""Convert CLI sdrf implementation."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from .options import _output_folder_option, _sdrf_file_option, _verbose_option
from .postprocess import _preserve_sdrf


@click.command("sdrf")
@_sdrf_file_option
@_output_folder_option
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
    default="sdrf",
)
@_verbose_option
def convert_sdrf_cmd(
    sdrf_file: Path,
    output_folder: Path,
    output_prefix: str,
    verbose: bool,
):
    """Convert SDRF metadata to QPX sample.parquet and run.parquet.

    Reads a Sample and Data Relationship Format (SDRF) file and produces
    the QPX sample and run data structures as Parquet files.

    \b
    Example:
        qpxc convert sdrf \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./qpx_output
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder.mkdir(parents=True, exist_ok=True)

    from qpx.converters.sdrf import SdrfConverter

    with SdrfConverter() as converter:
        converter.convert(
            sdrf_path=str(sdrf_file),
            sample_output=str(output_folder / f"{output_prefix}.sample.parquet"),
            run_output=str(output_folder / f"{output_prefix}.run.parquet"),
        )

    _preserve_sdrf(output_folder, output_prefix)

    click.echo(f"SDRF conversion complete. Output: {output_folder}")
