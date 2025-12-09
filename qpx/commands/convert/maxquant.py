"""
MaxQuant-specific converters for QPX formats.
"""

import logging
from pathlib import Path
from typing import Optional

import click

from qpx.core.maxquant.maxquant import MaxQuant
from qpx.core.project import create_uuid_filename
from qpx.utils.logger import get_logger


@click.group()
def convert():
    """Convert MaxQuant formats to QPX format."""
    pass


@convert.command(
    "maxquant-psm",
    short_help="Convert PSM data from MaxQuant msms.txt to parquet format",
)
@click.option(
    "--msms-file",
    help="MaxQuant msms.txt file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output folder",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--batch-size",
    help="Batch size",
    default=100000,
    type=int,
)
@click.option(
    "--output-prefix",
    help="Output file prefix",
)
@click.option(
    "--spectral-data",
    help="Include spectral data fields",
    is_flag=True,
)
@click.option(
    "--n-workers",
    help="Number of parallel workers",
    default=None,
    type=int,
)
@click.option(
    "--memory-limit",
    help="Memory limit in GB",
    default=None,
    type=float,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_maxquant_psm_cmd(
    msms_file: Path,
    output_folder: Path,
    batch_size: int,
    output_prefix: Optional[str],
    spectral_data: bool = False,
    n_workers: Optional[int] = None,
    memory_limit: Optional[float] = None,
    verbose: bool = False,
):
    """
    Convert MaxQuant PSM data from msms.txt to parquet format.
    
    This command takes a MaxQuant msms.txt file and converts it to the QPX
    parquet format for PSM data.
    
    Example:
        qpxc convert maxquant-psm \\
            --msms-file msms.txt \\
            --output-folder ./output \\
            --n-workers 8 \\
            --memory-limit 16
    """
    logger = get_logger("qpx.commands.maxquant")
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    try:
        if not all([msms_file, output_folder]):
            raise click.UsageError("ERROR: Please provide all required parameters")

        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using output directory: {output_folder}")

        prefix = output_prefix or "psm"
        filename = create_uuid_filename(prefix, ".psm.parquet")
        output_path = output_folder / filename
        logger.info(f"Will save PSM file as: {filename}")

        logger.info("Initializing MaxQuant PSM converter...")
        processor = MaxQuant(spectral_data, memory_limit_gb=memory_limit)

        if memory_limit:
            logger.info(f"Memory limit set to {memory_limit:.1f} GB")
        else:
            logger.info(
                f"Memory limit set to {processor.memory_limit_gb:.1f} GB (all available memory)"
            )

        workers_msg = f"{n_workers}" if n_workers else "CPU cores + 1"
        logger.info(
            f"Starting PSM conversion with parallel processing (batch size: {batch_size:,}, workers: {workers_msg})..."
        )

        processor.process_psm_file(
            msms_path=str(msms_file),
            output_path=str(output_path),
            chunksize=batch_size,
            n_workers=n_workers,
        )
        logger.info(f"PSM file successfully saved to: {output_path}")

    except Exception as e:
        logger.error(f"Error in MaxQuant PSM conversion: {str(e)}", exc_info=True)
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


@convert.command(
    "maxquant-feature",
    short_help="Convert feature data from MaxQuant evidence.txt to parquet format",
)
@click.option(
    "--evidence-file",
    help="MaxQuant evidence.txt file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output folder",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--protein-file",
    help="Protein list file for filtering",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--protein-groups-file",
    help="MaxQuant proteinGroups.txt file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--batch-size",
    help="Batch size",
    default=100000,
    type=int,
)
@click.option(
    "--output-prefix",
    help="Output file prefix",
)
@click.option(
    "--n-workers",
    help="Number of parallel workers",
    default=None,
    type=int,
)
@click.option(
    "--memory-limit",
    help="Memory limit in GB",
    default=None,
    type=float,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_maxquant_feature_cmd(
    evidence_file: Path,
    sdrf_file: Path,
    output_folder: Path,
    protein_file: Optional[Path],
    protein_groups_file: Path,
    batch_size: int,
    output_prefix: Optional[str],
    n_workers: Optional[int] = None,
    memory_limit: Optional[float] = None,
    verbose: bool = False,
):
    """
    Convert MaxQuant feature data from evidence.txt to parquet format.

    This command takes a MaxQuant evidence.txt file and converts it to the QPX
    parquet format for feature data, using metadata from an SDRF file.

    Example:
        qpxc convert maxquant-feature \\
            --evidence-file evidence.txt \\
            --sdrf-file data.sdrf.tsv \\
            --protein-groups-file proteinGroups.txt \\
            --output-folder ./output \\
            --n-workers 8 \\
            --memory-limit 16
    """
    logger = get_logger("qpx.commands.maxquant")
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    try:
        if not all([evidence_file, sdrf_file, output_folder]):
            raise click.UsageError("ERROR: Please provide all required parameters")

        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using output directory: {output_folder}")

        prefix = output_prefix or "feature"
        filename = create_uuid_filename(prefix, ".feature.parquet")
        output_path = output_folder / filename
        logger.info(f"Will save feature file as: {filename}")

        logger.info("Initializing MaxQuant feature converter...")
        processor = MaxQuant(memory_limit_gb=memory_limit)

        if memory_limit:
            logger.info(f"Memory limit set to {memory_limit:.1f} GB")
        else:
            logger.info(f"Memory limit set to {processor.memory_limit_gb:.1f} GB")

        workers_msg = f"{n_workers}" if n_workers else "CPU cores + 1"
        logger.info(
            f"Starting feature conversion (batch size: {batch_size:,}, workers: {workers_msg})..."
        )

        logger.info(
            f"Using proteinGroups file for Q-value mapping: {protein_groups_file}"
        )
        processor._init_protein_group_qvalue_mapping(str(protein_groups_file))

        processor.process_feature_file(
            evidence_path=str(evidence_file),
            output_path=str(output_path),
            sdrf_path=str(sdrf_file),
            protein_file=str(protein_file) if protein_file else None,
            chunksize=batch_size,
            n_workers=n_workers,
        )
        logger.info(f"Feature file successfully saved to: {output_path}")

    except Exception as e:
        logger.error(f"Error in MaxQuant feature conversion: {str(e)}", exc_info=True)
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


@convert.command(
    "maxquant-pg",
    short_help="Convert MaxQuant proteinGroups.txt to QPX protein group format",
)
@click.option(
    "--protein-groups-file",
    help="MaxQuant proteinGroups.txt file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--evidence-file",
    help="MaxQuant evidence.txt file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output folder",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--batch-size",
    help="Batch size",
    default=10000,
    type=int,
)
@click.option(
    "--output-prefix",
    help="Output file prefix",
)
@click.option(
    "--n-workers",
    help="Number of parallel workers",
    default=None,
    type=int,
)
@click.option(
    "--memory-limit",
    help="Memory limit in GB",
    default=None,
    type=float,
)
@click.option(
    "--standardized-intensities",
    help="Calculate standardized intensity metrics (total_all_peptides_intensity and top3_intensity)",
    is_flag=True,
    default=False,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_maxquant_pg_cmd(
    protein_groups_file: Path,
    sdrf_file: Path,
    evidence_file: Path,
    output_folder: Path,
    batch_size: int,
    output_prefix: Optional[str],
    n_workers: Optional[int] = None,
    memory_limit: Optional[float] = None,
    standardized_intensities: bool = False,
    verbose: bool = False,
):
    """
    Convert MaxQuant proteinGroups.txt to QPX parquet format.

    Example:
        qpxc convert maxquant-pg \\
            --protein-groups-file proteinGroups.txt \\
            --sdrf-file data.sdrf.tsv \\
            --evidence-file evidence.txt \\
            --output-folder ./output \\
            --n-workers 8 \\
            --memory-limit 16 \\
            --standardized-intensities
    """
    logger = get_logger("qpx.commands.maxquant")
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    try:
        if not all([protein_groups_file, sdrf_file, output_folder]):
            raise click.UsageError("ERROR: Please provide all required parameters")

        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using output directory: {output_folder}")

        prefix = output_prefix or "pg"
        filename = create_uuid_filename(prefix, ".pg.parquet")
        output_path = output_folder / filename
        logger.info(f"Will save protein groups file as: {filename}")

        logger.info("Initializing MaxQuant converter...")
        processor = MaxQuant(memory_limit_gb=memory_limit)

        if memory_limit:
            logger.info(f"Memory limit set to {memory_limit:.1f} GB")
        else:
            logger.info(f"Memory limit set to {processor.memory_limit_gb:.1f} GB")

        workers_msg = f"{n_workers}" if n_workers else "CPU cores + 1"
        logger.info(
            f"Starting conversion (batch size: {batch_size:,}, workers: {workers_msg})..."
        )

        if standardized_intensities:
            logger.info(
                "Standardized intensities enabled: will calculate "
                "total_all_peptides_intensity and top3_intensity"
            )

        processor.process_pg_file(
            protein_groups_path=str(protein_groups_file),
            output_path=str(output_path),
            sdrf_path=str(sdrf_file),
            evidence_path=str(evidence_file),
            chunksize=batch_size,
            n_workers=n_workers,
            calculate_standardized_intensities=standardized_intensities,
        )
        logger.info(f"Protein groups file successfully saved to: {output_path}")

    except Exception as e:
        logger.error(
            f"Error in MaxQuant protein groups conversion: {str(e)}", exc_info=True
        )
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


if __name__ == "__main__":
    convert()
