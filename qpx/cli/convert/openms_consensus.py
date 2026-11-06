"""Convert CLI openms consensus implementation."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from qpx.cli.resources import convert_resource_options

from .options import (
    _dataset_metadata_options,
    _fasta_option,
    _mudata_option,
    _openms_project_accession_option,
)
from .postprocess import _finish_conversion


def _consensusxml_list(ctx, param, value) -> tuple[Path, ...]:
    """Split a comma-separated list of consensusXML paths and check that each file exists."""
    paths = tuple(Path(item.strip()) for item in value.split(",") if item.strip())
    if not paths:
        raise click.BadParameter("no consensusXML file given", ctx=ctx, param=param)
    for path in paths:
        if not path.is_file():
            raise click.BadParameter(f"file '{path}' does not exist", ctx=ctx, param=param)
    return paths


@click.command("openms-consensus")
@click.option(
    "--consensusxml",
    "consensusxml_paths",
    required=True,
    callback=_consensusxml_list,
    help=(
        "OpenMS .consensusXML file (peptide feature intensities + IDs + protein inference), or a "
        "comma-separated list of files (e.g. one per sample group) to convert into one dataset; "
        "no run may appear in more than one file."
    ),
)
@click.option(
    "--sdrf-file",
    "sdrf_path",
    required=False,
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="SDRF metadata (run/sample views + grouped_runs fraction grouping). Optional.",
)
@click.option(
    "--output-folder",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="Output directory for the QPX views.",
)
@click.option("--output-prefix", default="openms", show_default=True)
@click.option(
    "--structures",
    default="feature,psm,pg,run,sample",
    show_default=True,
    help="Comma-separated views to write (run/sample require --sdrf-file).",
)
@click.option(
    "--streaming/--no-streaming",
    default=None,
    help=(
        "Force the low-memory streaming reader on/off. Default: auto (streaming for "
        "consensusXML files > 1 GB). The in-memory reader peaks at ~7x the file "
        "size in RAM; streaming peaks at ~2.3x and takes ~1.7x as long."
    ),
)
@click.option("--verbose", is_flag=True, help="Enable verbose logging.")
@_openms_project_accession_option
@click.option(
    "--include-unassigned-psms/--no-include-unassigned-psms",
    default=True,
    show_default=True,
    help=(
        "Keep unassigned peptide identifications by default to preserve identification "
        "evidence. Use --no-include-unassigned-psms to exclude them. Feature links "
        "are only written when both feature and PSM views are exported; a null "
        "feature_id means no recorded link. Protein inference uses every "
        "identification either way."
    ),
)
@convert_resource_options
@_fasta_option
@_mudata_option
@_dataset_metadata_options
def convert_openms_consensus_cmd(**kwargs):
    """Convert one or more OpenMS consensusXML files (+ SDRF) to one QPX dataset.

    Interim quantms path while OpenMS -out_qpx is pre-1.1. The feature view carries
    the per-run/channel peptide intensities from the consensusXML; the pg view
    carries an interim unnormalized unique-peptide-sum protein intensity (stamped
    with a quantification_method cv_param) until OpenMS provides the authoritative
    protein quant. Several comma-separated --consensusxml inputs are written into the same views.
    """
    consensusxml_paths = kwargs["consensusxml_paths"]
    sdrf_path = kwargs["sdrf_path"]
    output_folder = kwargs["output_folder"]
    output_prefix = kwargs["output_prefix"]
    structures = kwargs["structures"]
    streaming = kwargs["streaming"]
    verbose = kwargs["verbose"]
    project_accession = kwargs["project_accession"]
    include_unassigned_psms = kwargs["include_unassigned_psms"]
    fasta = kwargs["fasta"]
    mudata = kwargs["mudata"]
    metadata = kwargs["metadata"]

    if verbose:
        logging.basicConfig(level=logging.INFO)
    from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter

    structs = tuple(s.strip() for s in structures.split(",") if s.strip())
    try:
        written = OpenMSConsensusConverter().convert(
            consensusxml_path=[str(path) for path in consensusxml_paths],
            output_folder=str(output_folder),
            output_prefix=output_prefix,
            include_unassigned_psms=include_unassigned_psms,
            sdrf_path=str(sdrf_path) if sdrf_path else None,
            structures=structs,
            streaming=streaming,
            project_accession=project_accession,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc
    _finish_conversion(
        output_folder,
        output_prefix,
        project_accession,
        metadata,
        fasta=fasta,
        mudata=mudata,
        written=bool(written),
    )

    click.echo(f"consensusXML conversion complete. Wrote: {sorted(written)}")
