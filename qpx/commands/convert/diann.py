"""
Command to convert DIA-NN report files to QPX format.
"""

import logging
import os
from pathlib import Path
from typing import Optional

import click

from qpx.core.diann.diann import DiaNNConvert
from qpx.core.project import create_uuid_filename
from qpx.utils.logger import get_logger


def convert_diann(
    report_path: Path,
    qvalue_threshold: float,
    mzml_info_folder: Path,
    sdrf_path: Path,
    output_folder: Path,
    protein_file: Optional[Path] = None,
    output_prefix: Optional[str] = None,
    partitions: Optional[str] = None,
    duckdb_max_memory: Optional[str] = None,
    duckdb_threads: Optional[int] = None,
    batch_size: int = 100,
    verbose: bool = False,
) -> None:
    """
    Convert DIA-NN report to QPX parquet format.

    Args:
        report_path: DIA-NN report file path
        qvalue_threshold: Q-value threshold for filtering
        mzml_info_folder: mzML info file folder
        sdrf_path: SDRF file path for metadata
        output_folder: Output directory for generated files
        protein_file: Optional protein file with specific requirements
        output_prefix: Optional prefix for output files
        partitions: Optional field(s) for splitting files (comma-separated)
        duckdb_max_memory: Optional maximum memory for DuckDB (e.g., "4GB")
        duckdb_threads: Optional number of threads for DuckDB
        file_num: Number of files to process simultaneously (default: 100)
        verbose: Enable verbose logging
    """
    logger = get_logger("qpx.commands.diann")
    if verbose:
        logger.setLevel(logging.DEBUG)

    try:
        if not all([report_path, mzml_info_folder, output_folder, sdrf_path]):
            raise click.UsageError("Please provide all required parameters")

        # Create output directory if it doesn't exist
        os.makedirs(output_folder, exist_ok=True)

        # Set default prefix if not provided
        prefix = output_prefix or "feature"
        filename = create_uuid_filename(prefix, ".feature.parquet")
        feature_output_path = output_folder / filename

        dia_nn = DiaNNConvert(
            diann_report=report_path,
            pg_matrix_path=None,
            sdrf_path=sdrf_path,
            duckdb_max_memory=duckdb_max_memory,
            duckdb_threads=duckdb_threads,
        )

        if not partitions:
            dia_nn.write_feature_to_file(
                qvalue_threshold=qvalue_threshold,
                mzml_info_folder=mzml_info_folder,
                output_path=str(feature_output_path),
                file_num=batch_size,
                protein_file=protein_file,
            )
        else:
            partition_list = partitions.split(",")
            dia_nn.write_features_to_file(
                qvalue_threshold=qvalue_threshold,
                mzml_info_folder=mzml_info_folder,
                output_folder=str(output_folder),
                filename=filename,
                partitions=partition_list,
                file_num=batch_size,
                protein_file=protein_file,
            )

    except Exception as e:
        logger.exception(f"Error in DIA-NN conversion: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


def convert_diann_pg(
    report_path: Path,
    pg_matrix_path: Path,
    sdrf_path: Path,
    output_folder: Path,
    output_prefix: Optional[str] = None,
    duckdb_max_memory: Optional[str] = None,
    duckdb_threads: Optional[int] = None,
    batch_size: int = 100,
    standardized_intensities: bool = False,
    verbose: bool = False,
) -> None:
    """
    Convert DIA-NN report to protein group format.

    Args:
        report_path: DIA-NN report file path
        pg_matrix_path: DIA-NN protein quantities table file path
        sdrf_path: SDRF file path for metadata
        output_folder: Output directory for generated files
        output_prefix: Optional prefix for output files
        duckdb_max_memory: Optional maximum memory for DuckDB (e.g., "4GB")
        duckdb_threads: Optional number of threads for DuckDB
        file_num: Number of files to process simultaneously (default: 100)
        standardized_intensities: Calculate total_all_peptides_intensity and
            top3_intensity (default: False)
        verbose: Enable verbose logging
    """
    logger = get_logger("qpx.commands.diann")
    if verbose:
        logger.setLevel(logging.DEBUG)

    try:
        if not all([report_path, pg_matrix_path, sdrf_path, output_folder]):
            raise click.UsageError("Please provide all required parameters")

        # Create output directory if it doesn't exist
        os.makedirs(output_folder, exist_ok=True)

        # Set default prefix if not provided
        prefix = output_prefix or "pg"
        filename = create_uuid_filename(prefix, ".pg.parquet")
        pg_output_path = output_folder / filename

        if standardized_intensities:
            logger.info(
                "Standardized intensities enabled: will calculate "
                "total_all_peptides_intensity and top3_intensity"
            )

        dia_nn = DiaNNConvert(
            diann_report=report_path,
            pg_matrix_path=pg_matrix_path,
            sdrf_path=sdrf_path,
            duckdb_max_memory=duckdb_max_memory,
            duckdb_threads=duckdb_threads,
        )

        dia_nn.write_pg_matrix_to_file(
            output_path=str(pg_output_path),
            file_num=batch_size,
            calculate_standardized_intensities=standardized_intensities,
        )

    except Exception as e:
        logger.exception(f"Error in DIA-NN PG conversion: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


@click.command(
    "diann",
    short_help="Convert DIA-NN report to QPX format",
)
@click.option(
    "--report-path",
    help="DIA-NN report file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--qvalue-threshold",
    help="Q-value threshold for filtering",
    required=True,
    default=0.05,
    type=float,
)
@click.option(
    "--mzml-info-folder",
    help="mzML info file folder",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-path",
    help="SDRF file path for metadata",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--protein-file",
    help="Protein file with specific requirements",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option("--output-prefix", help="Prefix for output files")
@click.option("--partitions", help="Field(s) for splitting files (comma-separated)")
@click.option("--duckdb-max-memory", help="Maximum memory for DuckDB (e.g., '4GB')")
@click.option("--duckdb-threads", help="Number of threads for DuckDB", type=int)
@click.option(
    "--batch-size",
    help="Number of files to process simultaneously",
    default=100,
    type=int,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_diann_cmd(**kwargs):
    """
    Convert DIA-NN report to QPX format.
    
    This command takes a DIA-NN report file and converts it to the QPX
    parquet format. The conversion includes feature data and can optionally
    split the output into multiple files based on specified fields.
    
    Example:
        qpxc convert diann \\
            --report-path report.tsv \\
            --qvalue-threshold 0.05 \\
            --mzml-info-folder ./mzml_info \\
            --sdrf-path data.sdrf.tsv \\
            --output-folder ./output
    """
    convert_diann(**kwargs)


@click.command(
    "diann-pg",
    short_help="Convert DIA-NN report to protein group format",
)
@click.option(
    "--report-path",
    help="DIA-NN report file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--pg-matrix-path",
    help="DIA-NN protein quantities table file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-path",
    help="SDRF file path for metadata",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option("--output-prefix", help="Prefix for output files")
@click.option("--duckdb-max-memory", help="Maximum memory for DuckDB (e.g., '4GB')")
@click.option("--duckdb-threads", help="Number of threads for DuckDB", type=int)
@click.option(
    "--batch-size",
    help="Number of files to process simultaneously",
    default=100,
    type=int,
)
@click.option(
    "--standardized-intensities",
    help="Calculate standardized intensity metrics (total_all_peptides_intensity and top3_intensity)",
    is_flag=True,
    default=False,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_diann_pg_cmd(**kwargs):
    """
    Convert DIA-NN report to protein group format.
    
    This command takes a DIA-NN report file and converts it to the QPX
    protein group format in parquet format.
    
    Example:
        qpxc convert diann-pg \\
            --report-path report.tsv \\
            --pg-matrix-path report.pg_matrix.tsv \\
            --sdrf-path data.sdrf.tsv \\
            --output-folder ./output \\
            --standardized-intensities
    """
    convert_diann_pg(**kwargs)
