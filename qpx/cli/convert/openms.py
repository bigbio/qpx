"""Convert CLI openms implementation."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from qpx.converters.openms import OpenMSConverter

from .options import _compression_option, _openms_project_accession_option, _verbose_option
from .postprocess import _maybe_enrich_pride, _preserve_sdrf, _record_conversion_options


@click.command("openms", deprecated=True)
@click.option(
    "--qpx-dir",
    help="Directory containing OpenMS -out_qpx parquet files (*.psm.parquet, *.feature.parquet, *.pg.parquet)",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file path (for sample/run generation)",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--consensusxml",
    "consensusxml_path",
    help="OpenMS .consensusXML (the -out_cxml companion of -out_qpx); its "
    "ColumnHeaders give the authoritative channel count/order for relabeling. "
    "Peptide identifications, when present, recover PSM runs using exact spectrum matches.",
    required=False,
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for the full QPX dataset",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
    default="openms",
)
@_openms_project_accession_option
@click.option(
    "--enrich-pride",
    help="Fetch project metadata from PRIDE API after conversion",
    is_flag=True,
    default=False,
)
@_compression_option
@_verbose_option
def convert_openms_cmd(**kwargs):
    """[DEPRECATED] Enrich OpenMS ProteomicsLFQ -out_qpx output into a full QPX dataset.

    DEPRECATED: Some OpenMS exporters assign PSMs to the first run or emit
    conflicting identities (OpenMS#9872, OpenMS#9871). Prefer
    ``qpxc convert openms-consensus`` to read the original consensusXML directly.

    With --consensusxml, uniquely matching spectrum evidence restores PSM runs.
    Missing or ambiguous spectrum matches fail before existing core outputs are
    replaced. Duplicate identities are retained with warnings; use
    ``qpxc validate`` to report them as errors. Rows are not discarded
    or assigned arbitrary IDs to resolve conflicts.

    Validates the existing psm/feature/pg parquet files, copies them to the
    output folder, and generates the missing metadata tables (run, sample,
    ontology, provenance, dataset) from the SDRF file.

    \b
    Examples:
        # Enrich OpenMS QPX output
        qpxc convert openms \\
            --qpx-dir ./openms_qpx_output \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./qpx_full

    \b
        # With project accession
        qpxc convert openms \\
            --qpx-dir ./openms_qpx_output \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./qpx_full \\
            --project-accession PXD001819
    """
    qpx_dir = kwargs["qpx_dir"]
    sdrf_file = kwargs["sdrf_file"]
    output_folder = kwargs["output_folder"]
    output_prefix = kwargs["output_prefix"]
    project_accession = kwargs["project_accession"]
    enrich_pride = kwargs["enrich_pride"]
    compression = kwargs["compression"]
    verbose = kwargs["verbose"]

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    converter = OpenMSConverter(
        qpx_dir=qpx_dir,
        sdrf_path=sdrf_file,
        consensusxml_path=kwargs.get("consensusxml_path"),
        compression=compression,
    )
    try:
        converter.convert(
            output_folder=output_folder,
            output_prefix=output_prefix,
            project_accession=project_accession,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    _preserve_sdrf(output_folder, output_prefix)
    _record_conversion_options(output_folder, output_prefix)
    _maybe_enrich_pride(output_folder, project_accession, enrich_pride, output_prefix)

    click.echo(f"OpenMS QPX enrichment complete. Output: {output_folder}")
