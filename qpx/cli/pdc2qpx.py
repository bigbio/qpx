"""``qpxc pdc2qpx`` — download a PDC/CPTAC study and convert it to QPX."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from qpx.cli.resources import convert_resource_options


@click.command("pdc2qpx")
@click.option(
    "-a",
    "--accession",
    required=True,
    help="PDC study accession(s): a single ID, comma-separated IDs, or a CSV with a pdc_study_id/pdc_id column",
)
@click.option(
    "--download-dir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="Directory for downloaded PDC files (written under <download-dir>/<study>/)",
)
@click.option(
    "--output-folder",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="QPX output dir. One study writes here; multiple studies write to <output-folder>/<study>/",
)
@click.option(
    "--include-spectra",
    is_flag=True,
    default=False,
    help="Also download mzML and produce a full-spectra <study>.mz.parquet",
)
@click.option(
    "--no-metadata",
    is_flag=True,
    default=False,
    help="Skip building sample/run views from PDC metadata (built by default)",
)
@click.option(
    "--ms-levels",
    default=None,
    help="MS levels for the mz view (e.g. '2' or '1,2'). Default: all levels.",
)
@click.option("--max-cpus", type=int, help="Threads for the CDAP conversion")
@click.option("--max-memory", help="Memory limit for the CDAP conversion")
@convert_resource_options
@click.option(
    "--skip-download",
    is_flag=True,
    default=False,
    help="Use files already present under <download-dir>/<study>/ (skip downloading)",
)
@click.option(
    "--stop-on-error",
    is_flag=True,
    default=False,
    help="With multiple studies, abort on the first failure (default: continue and report)",
)
@click.option("--verbose", is_flag=True, help="Enable verbose logging")
def pdc2qpx_cmd(**kwargs):
    """Download PDC/CPTAC studies and convert them to QPX datasets.

    ``--accession`` accepts a single ID, comma-separated IDs, or a CSV with a
    ``pdc_study_id``/``pdc_id`` column (same input forms as ``pridepy``). Each
    study's CDAP ``.psm`` (and, with ``--include-spectra``, mzML) is converted to
    psm/feature/pg/sample/run/... views. With several studies the batch
    continues past a failed study (e.g. one with no ``.psm``) unless
    ``--stop-on-error`` is set. Requires the ``pdc`` extra (``pip install
    qpx[pdc]``) for the pridepy download/parse layer.

    \b
    Examples:
        # One study (full spectra for quantms reanalysis)
        qpxc pdc2qpx -a PDC000109 \\
            --download-dir ./downloads --output-folder ./qpx/PDC000109 \\
            --include-spectra

    \b
        # Many studies from a CSV (each -> ./qpx/<study>/)
        qpxc pdc2qpx -a cptac_lfq.csv \\
            --download-dir ./downloads --output-folder ./qpx --max-cpus 24
    """
    accession = kwargs["accession"]
    download_dir = kwargs["download_dir"]
    output_folder = kwargs["output_folder"]
    include_spectra = kwargs["include_spectra"]
    no_metadata = kwargs["no_metadata"]
    ms_levels = kwargs["ms_levels"]
    max_cpus = kwargs["max_cpus"]
    max_memory = kwargs["max_memory"]
    skip_download = kwargs["skip_download"]
    stop_on_error = kwargs["stop_on_error"]
    verbose = kwargs["verbose"]

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    levels = [int(x.strip()) for x in ms_levels.split(",")] if ms_levels else None

    from qpx.pipeline.pdc2qpx import parse_accessions, run_pdc2qpx, run_pdc2qpx_batch

    studies = parse_accessions(accession)

    common = {
        "include_spectra": include_spectra,
        "include_metadata": not no_metadata,
        "ms_levels": levels,
        "max_cpus": max_cpus,
        "max_memory": max_memory,
        "skip_download": skip_download,
    }

    if len(studies) == 1:
        outputs = run_pdc2qpx(studies[0], download_dir, output_folder, **common)
        click.echo(f"pdc2qpx complete for {studies[0]}. Output: {output_folder}")
        if "mz" in outputs:
            click.echo(f"Full spectra: {outputs['mz']}")
        return

    results = run_pdc2qpx_batch(studies, download_dir, output_folder, continue_on_error=not stop_on_error, **common)
    n_ok = sum(1 for record in results.values() if record["status"] == "ok")
    click.echo(f"\npdc2qpx: {n_ok}/{len(results)} studies ok (output under {output_folder}/<study>/)")
    for study, record in results.items():
        suffix = f" ({record['error']})" if record["error"] else ""
        click.echo(f"  {record['status']:6s} {study}{suffix}")
