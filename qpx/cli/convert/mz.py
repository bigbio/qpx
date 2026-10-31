"""Convert CLI mz implementation."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import click

from .options import _verbose_option


@click.command("mz")
@click.option(
    "--mzml-dir",
    help="Directory containing mzML / mzML.gz spectra files (one study)",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--output",
    "-o",
    help="Output .mz.parquet file path",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--ms-levels",
    default=None,
    help="Comma-separated MS levels to include (e.g. '2' or '1,2'). Default: all levels.",
)
@_verbose_option
def convert_mz_cmd(
    mzml_dir: Path,
    output: Path,
    ms_levels: Optional[str],
    verbose: bool,
):
    """Convert a directory of mzML spectra to a QPX mz.parquet (full spectra).

    Reads every .mzML / .mzML.gz file in --mzml-dir and writes scan-level
    spectral data (m/z and intensity arrays, precursor info) to a single
    mz.parquet. Each spectrum carries run_file_name + scan so it can be linked
    back to PSM / feature records.

    \b
    Examples:
        # All MS levels (full spectra; required for quantms precursor LFQ)
        qpxc convert mz \\
            --mzml-dir /data/CPTAC/PDC000109/mzml \\
            --output ./qpx_output/PDC000109.mz.parquet

    \b
        # Only MS2 spectra
        qpxc convert mz \\
            --mzml-dir /data/CPTAC/PDC000109/mzml \\
            --output ./PDC000109.mz.parquet \\
            --ms-levels 2
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    levels = [int(x.strip()) for x in ms_levels.split(",")] if ms_levels else None

    output.parent.mkdir(parents=True, exist_ok=True)

    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    mapper = SpectraMappingTransform(mzml_directory=mzml_dir)
    mapper.write_mz_parquet_from_dir(output, ms_levels=levels)

    click.echo(f"mz conversion complete. Output: {output}")
