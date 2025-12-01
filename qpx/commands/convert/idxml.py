import click
from pathlib import Path
from typing import Optional, Union, List, Dict, Tuple
import glob
import logging
import tempfile
import shutil

from qpx.core.project import create_uuid_filename
from qpx.core.idxml import IdXML, merge_idxml_parquet_files


@click.command(
    "convert-idxml",
    short_help="Convert IdXML to PSM parquet file in QPX",
)
@click.option(
    "--idxml-file",
    help="the IdXML file containing identifications",
    required=True,
)
@click.option(
    "--output-folder",
    help="Folder where the parquet file will be generated",
    required=True,
)
@click.option(
    "--mzml-file",
    help="Optional mzML to attach spectra by scan",
    required=False,
)
@click.option(
    "--output-prefix-file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
@click.option(
    "--spectral-data",
    help="Spectral data fields (optional)",
    is_flag=True,
)
def convert_idxml_file(
    idxml_file: Union[Path, str],
    output_folder: str,
    mzml_file: Optional[Union[Path, str]],
    output_prefix_file: Optional[str],
    spectral_data: bool = False,
) -> None:
    """Convert a single OpenMS idXML file to QPX PSM format.
    
    Converts PSM data from OpenMS idXML format to the QPX standardized parquet format. 
    Can optionally attach spectral information from corresponding mzML files.
    
    Example:
        qpxc convert idxml \\
            --idxml-file /path/to/data.idXML \\
            --output-folder ./output
    """
    if idxml_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = "psm"

    parser = IdXML(
        idxml_path=idxml_file, mzml_path=mzml_file, spectral_data=spectral_data
    )
    output_path = (
        f"{output_folder}/{create_uuid_filename(output_prefix_file, '.psm.parquet')}"
    )
    parser.to_parquet(output_path)


# Helper functions for batch conversion


def _validate_batch_parameters(
    idxml_folder: Optional[str],
    idxml_files: Optional[str],
    mzml_folder: Optional[str],
    mzml_files: Optional[str],
) -> None:
    """Validate input parameters for batch conversion."""
    if not idxml_folder and not idxml_files:
        raise click.UsageError("Please provide either --idxml-folder or --idxml-files")

    if idxml_folder and idxml_files:
        raise click.UsageError(
            "Please provide only one of --idxml-folder or --idxml-files"
        )

    if mzml_folder and mzml_files:
        raise click.UsageError(
            "Please provide only one of --mzml-folder or --mzml-files"
        )


def _get_idxml_files(
    idxml_folder: Optional[str], idxml_files: Optional[str], logger: logging.Logger
) -> List[str]:
    """Get list of idXML files to process."""
    if idxml_folder:
        idxml_pattern = str(Path(idxml_folder) / "*.idXML")
        idxml_file_paths = glob.glob(idxml_pattern)
        if not idxml_file_paths:
            raise click.UsageError(f"No .idXML files found in folder: {idxml_folder}")
        logger.info(f"Found {len(idxml_file_paths)} idXML files in {idxml_folder}")
        return idxml_file_paths
    else:
        idxml_file_paths = [f.strip() for f in idxml_files.split(",")]
        for idxml_file in idxml_file_paths:
            if not Path(idxml_file).exists():
                raise click.UsageError(f"IdXML file not found: {idxml_file}")
        logger.info(f"Processing {len(idxml_file_paths)} specified idXML files")
        return idxml_file_paths


def _setup_mzml_mapping(
    mzml_folder: Optional[str],
    mzml_files: Optional[str],
    idxml_files: Optional[str],
    logger: logging.Logger,
) -> Tuple[List[str], Optional[Dict[str, str]], bool]:
    """
    Setup mzML file mapping strategy.

    Returns:
        Tuple of (mzml_file_paths, mzml_map, use_index_mapping)
    """
    mzml_file_paths = []
    mzml_map = None
    use_index_mapping = False

    if mzml_files:
        mzml_file_paths = [f.strip() for f in mzml_files.split(",")]
        for mzml_file in mzml_file_paths:
            if not Path(mzml_file).exists():
                raise click.UsageError(f"mzML file not found: {mzml_file}")
        logger.info(f"Using {len(mzml_file_paths)} specified mzML files")

        if idxml_files:  # Both are file lists
            use_index_mapping = True
            logger.info(
                "Using index-based mapping: "
                "mzML files will be matched by position in the list"
            )
        else:  # idxml is folder, mzml is list - use basename matching
            mzml_map = {Path(f).stem: f for f in mzml_file_paths}
            logger.info("Using basename matching for mzML files")

    return mzml_file_paths, mzml_map, use_index_mapping


def _find_matching_mzml(
    idxml_file: str,
    index: int,
    use_index_mapping: bool,
    mzml_folder: Optional[str],
    mzml_file_paths: List[str],
    mzml_map: Optional[Dict[str, str]],
    logger: logging.Logger,
) -> Optional[str]:
    """Find matching mzML file for an idXML file."""
    if use_index_mapping:
        if index < len(mzml_file_paths):
            mzml_file = mzml_file_paths[index]
            logger.info(f"Matched by index [{index}]: {Path(mzml_file).name}")
            return mzml_file
        logger.warning(f"No mzML file at index {index} for {idxml_file}")
        return None

    if mzml_folder:
        idxml_basename = Path(idxml_file).stem
        mzml_pattern = str(Path(mzml_folder) / f"{idxml_basename}.mzML")
        mzml_matches = glob.glob(mzml_pattern)
        if mzml_matches:
            mzml_file = mzml_matches[0]
            logger.info(f"Matched by basename: {Path(mzml_file).name}")
            return mzml_file
        logger.warning(
            f"No matching mzML file found for {idxml_file} (basename: {idxml_basename})"
        )
        return None

    if mzml_map:
        idxml_basename = Path(idxml_file).stem
        if idxml_basename in mzml_map:
            mzml_file = mzml_map[idxml_basename]
            logger.info(f"Matched by basename: {Path(mzml_file).name}")
            return mzml_file
        logger.warning(
            f"No matching mzML file found for {idxml_file} (basename: {idxml_basename})"
        )
        return None

    return None


def _convert_single_idxml(
    idxml_file: str,
    mzml_file: Optional[str],
    temp_dir: str,
    index: int,
    logger: logging.Logger,
) -> str:
    """Convert a single idXML file to parquet."""
    spectral_data = mzml_file is not None

    parser = IdXML(
        idxml_path=idxml_file,
        mzml_path=mzml_file,
        spectral_data=spectral_data,
    )

    temp_parquet_file = (
        Path(temp_dir) / f"temp_{index}_{create_uuid_filename('psm', '.psm.parquet')}"
    )
    parser.to_parquet(str(temp_parquet_file))

    logger.info(
        f"Converted {idxml_file} -> {temp_parquet_file} ({parser.get_psm_count()} PSMs)"
    )
    return str(temp_parquet_file)


def _convert_batch_files(
    idxml_file_paths: List[str],
    temp_dir: str,
    use_index_mapping: bool,
    mzml_folder: Optional[str],
    mzml_file_paths: List[str],
    mzml_map: Optional[Dict[str, str]],
    logger: logging.Logger,
    verbose: bool,
) -> List[str]:
    """Convert all idXML files to parquet."""
    parquet_files = []

    for i, idxml_file in enumerate(idxml_file_paths):
        logger.info(f"Processing {i+1}/{len(idxml_file_paths)}: {idxml_file}")

        mzml_file = _find_matching_mzml(
            idxml_file,
            i,
            use_index_mapping,
            mzml_folder,
            mzml_file_paths,
            mzml_map,
            logger,
        )

        try:
            parquet_file = _convert_single_idxml(
                idxml_file, mzml_file, temp_dir, i, logger
            )
            parquet_files.append(parquet_file)
        except Exception as e:
            logger.error(f"Failed to convert {idxml_file}: {e}")
            if verbose:
                raise

    return parquet_files


@click.command(
    "convert-idxml-batch",
    short_help="Convert multiple IdXML files to a single merged PSM parquet file",
)
@click.option(
    "--idxml-folder",
    help="Folder containing IdXML files to convert",
    required=False,
    type=click.Path(exists=True, file_okay=False),
)
@click.option(
    "--idxml-files",
    help="Comma-separated list of IdXML file paths",
    required=False,
)
@click.option(
    "--output-folder",
    help="Folder where the merged parquet file will be generated",
    required=True,
)
@click.option(
    "--output-prefix-file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
    default="merged-psm",
)
@click.option(
    "--mzml-folder",
    help="Optional folder containing mzML files to attach spectra by scan",
    required=False,
    type=click.Path(exists=True, file_okay=False),
)
@click.option(
    "--mzml-files",
    help="Comma-separated list of mzML file paths",
    required=False,
)
@click.option(
    "--verbose",
    help="Enable verbose logging",
    is_flag=True,
)
def convert_idxml_batch(
    idxml_folder: Optional[str],
    idxml_files: Optional[str],
    output_folder: str,
    output_prefix_file: str,
    mzml_folder: Optional[str],
    mzml_files: Optional[str],
    verbose: bool,
) -> None:
    """Convert multiple OpenMS idXML files to a single merged PSM parquet file.
    
    Batch converts multiple idXML files and merges them into a single QPX PSM parquet file. 
    Supports both folder-based and file-list-based input, with flexible mzML matching strategies.
    
    Example:
        qpxc convert idxml-batch \\
            --idxml-folder ./idxml_files \\
            --output-folder ./output \\
            --output-prefix-file batch_psm
    """
    # Setup logging
    logger = logging.getLogger(__name__)
    log_level = logging.DEBUG if verbose else logging.INFO
    log_format = (
        "%(asctime)s - %(levelname)s - %(message)s"
        if verbose
        else "%(asctime)s - %(message)s"
    )
    logging.basicConfig(level=log_level, format=log_format)

    # Validate parameters
    _validate_batch_parameters(idxml_folder, idxml_files, mzml_folder, mzml_files)

    # Get list of idXML files
    idxml_file_paths = _get_idxml_files(idxml_folder, idxml_files, logger)

    # Setup mzML file mapping
    mzml_file_paths, mzml_map, use_index_mapping = _setup_mzml_mapping(
        mzml_folder, mzml_files, idxml_files, logger
    )

    # Create output folder
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    # Create temporary directory
    temp_dir = tempfile.mkdtemp(prefix="idxml_batch_")
    logger.info(f"Using temporary directory: {temp_dir}")

    try:
        # Convert all files
        parquet_files = _convert_batch_files(
            idxml_file_paths,
            temp_dir,
            use_index_mapping,
            mzml_folder,
            mzml_file_paths,
            mzml_map,
            logger,
            verbose,
        )

        # Validate results
        if not parquet_files:
            raise click.ClickException("No parquet files were successfully generated")

        # Merge all parquet files
        logger.info(f"Merging {len(parquet_files)} parquet files...")
        output_path = Path(output_folder) / create_uuid_filename(
            output_prefix_file, ".psm.parquet"
        )
        merge_idxml_parquet_files(parquet_files, str(output_path))

        logger.info(f"Successfully created merged parquet file: {output_path}")

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir)
        logger.info(f"Cleaned up temporary directory: {temp_dir}")
