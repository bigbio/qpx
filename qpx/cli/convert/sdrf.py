"""Convert CLI sdrf implementation."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from .options import _output_folder_option, _sdrf_file_option, _verbose_option
from .postprocess import _preserve_sdrf


@click.command("sdrf")
@_sdrf_file_option
@click.option(
    "--mzml-dir",
    help="Directory containing mzML / mzML.gz files for acquisition starts missing from the SDRF",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
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
    mzml_dir: Path | None,
):
    """Convert SDRF metadata to QPX sample.parquet and run.parquet.

    Reads a Sample and Data Relationship Format (SDRF) file and produces
    the QPX sample and run data structures as Parquet files.
    Acquisition starts come from comment[acquisition date], or mzML headers
    with --mzml-dir, and retain the source text without timezone conversion.

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
            mzml_dir=mzml_dir,
        )

    _preserve_sdrf(output_folder, output_prefix)

    click.echo(f"SDRF conversion complete. Output: {output_folder}")
