"""Convert CLI diann implementation."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from qpx.cli.resources import convert_resource_options

from .options import (
    _compression_option,
    _dataset_metadata_options,
    _fasta_option,
    _max_cpus_option,
    _max_memory_option,
    _mudata_option,
    _output_folder_option,
    _output_prefix_option,
    _project_accession_option,
    _sdrf_file_option,
    _verbose_option,
)
from .postprocess import _finish_conversion


def _write_diann_views(converter, options: dict, prefix: str) -> None:
    """Write DIA-NN data views and their metadata in conversion order."""
    output_folder = options["output_folder"]
    output_prefix = options["output_prefix"]
    mzml_info_folder = options["mzml_info_folder"]
    qvalue_threshold = options["qvalue_threshold"]
    pg_matrix_path = options["pg_matrix_path"]
    project_accession = options["project_accession"]
    converter.convert_features(
        mzml_info_folder=mzml_info_folder,
        qvalue_threshold=qvalue_threshold,
        output_folder=output_folder,
        output_prefix=output_prefix,
        fragments=options["fragments"],
    )
    converter.convert_pg(
        pg_matrix_path=pg_matrix_path,
        output_folder=output_folder,
        output_prefix=output_prefix,
        qvalue_threshold=qvalue_threshold,
    )

    converter.convert_mz(
        mzml_statistics_folder=mzml_info_folder,
        output_folder=output_folder,
        output_prefix=output_prefix,
        feature_parquet=output_folder / f"{prefix}.feature.parquet",
    )
    converter.convert_sdrf(output_folder=output_folder, prefix=prefix)
    converter.write_ontology(output_folder, prefix=prefix)
    converter.write_provenance(output_folder, prefix=prefix)
    converter.write_dataset(output_folder, prefix=prefix, project_accession=project_accession)


@click.command("diann")
@click.option(
    "--report-path",
    help="DIA-NN report file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@_sdrf_file_option
@click.option(
    "--mzml-info-folder",
    help=(
        "Folder containing mzML info files (optional; scan/mz fields left empty if omitted). "
        "Its *_ms2_info.parquet files, when present, also produce the mz view of the feature-matched MS2 spectra."
    ),
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--qvalue-threshold",
    help=(
        "Optional q-value threshold. Unset (the default) converts the DIA-NN "
        "report as reported — no filtering — since DIA-NN already FDR-filters its "
        "main report and every per-row q-value column is carried through for "
        "downstream filtering. When given, the feature and pg views keep only precursors "
        "with Q.Value below it; protein-group q-values are carried, not filtered."
    ),
    default=None,
    type=float,
)
@_output_folder_option
@click.option("--fragments", is_flag=True, help="Also write Fragment.Quant.* as the fragment view")
@_output_prefix_option
@click.option(
    "--pg-matrix-path",
    help=(
        "DIA-NN protein-group matrix (optional). The pg view comes from the report; the matrix fills "
        "a label-free group's MaxLFQ when the report lacks it."
    ),
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@_max_memory_option
@_max_cpus_option
@convert_resource_options
@_project_accession_option
@_compression_option
@click.option(
    "--diann-log",
    help="DIA-NN summary log file (version auto-detected from first line)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@_fasta_option
@_mudata_option
@_verbose_option
@_dataset_metadata_options
def convert_diann_cmd(**kwargs):
    """Convert DIA-NN report to QPX format.

    Reads a DIA-NN report.tsv file and converts its precursor and protein-group
    quantification data into QPX Parquet format.

    \b
    Examples:
        # Feature and protein-group conversion
        qpxc convert diann \\
            --report-path report.tsv \\
            --sdrf-file data.sdrf.tsv \\
            --mzml-info-folder ./mzml_info \\
            --output-folder ./qpx_output

        # With the protein-group matrix
        qpxc convert diann \\
            --report-path report.tsv \\
            --sdrf-file data.sdrf.tsv \\
            --mzml-info-folder ./mzml_info \\
            --pg-matrix-path report.pg_matrix.tsv \\
            --output-folder ./qpx_output
    """
    report_path = kwargs["report_path"]
    sdrf_file = kwargs["sdrf_file"]
    output_folder = kwargs["output_folder"]
    output_prefix = kwargs["output_prefix"]
    max_memory = kwargs["max_memory"]
    max_cpus = kwargs["max_cpus"]
    project_accession = kwargs["project_accession"]
    compression = kwargs["compression"]
    diann_log = kwargs["diann_log"]
    fasta = kwargs["fasta"]
    mudata = kwargs["mudata"]
    verbose = kwargs["verbose"]
    metadata = kwargs["metadata"]

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder.mkdir(parents=True, exist_ok=True)
    prefix = output_prefix or "diann"

    from qpx.converters.diann import DiaNNConverter

    converter = DiaNNConverter(
        report_path=report_path,
        sdrf_path=sdrf_file,
        max_memory=max_memory,
        max_cpus=max_cpus,
        compression=compression,
        diann_log=str(diann_log) if diann_log else None,
    )
    _write_diann_views(converter, kwargs, prefix)

    _finish_conversion(output_folder, prefix, project_accession, metadata, fasta=fasta, mudata=mudata)

    click.echo(f"DIA-NN conversion complete. Output: {output_folder}")
