"""
Project-level conversion for QPX formats.
"""

import os
import sys
from pathlib import Path
from typing import Optional

import click

from qpx.commands.convert.quantms import (
    convert_quantms_feature_cmd,
    convert_quantms_psm_cmd,
)
from qpx.commands.utils.attach import attach_file_to_json_cmd
from qpx.core.project import check_directory, create_uuid_filename
from qpx.operate.tools import write_ibaq_feature


def find_file(directory: str, pattern: str) -> Optional[Path]:
    """Find first file matching pattern in directory."""
    path = Path(directory)
    files = list(path.rglob(pattern))
    return files[0] if files else None


def get_project_prefix(sdrf_file: Path) -> str:
    """Extract project prefix from SDRF filename (e.g. 'PXD000865' from 'PXD000865.sdrf.tsv')."""
    filename = sdrf_file.name
    # Remove .sdrf.tsv and any variations like _openms_design.sdrf.tsv
    prefix = filename.split(".sdrf")[0].split("_openms")[0]
    return prefix


def check_dir(folder_path: str) -> None:
    """Create directory if it doesn't exist."""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


def _setup_project_paths(base_folder: str) -> tuple[Path, Path, Path, Path]:
    """Setup and validate project paths."""
    print("\nSetting up input paths...")
    quant_tables = Path(base_folder) / "quant_tables"
    sdrf_dir = Path(base_folder) / "sdrf"
    spectra_dir = Path(base_folder) / "spectra"

    # Find required files
    print("Searching for required files...")
    mztab_file = find_file(quant_tables, "*.mzTab")
    msstats_file = find_file(quant_tables, "*msstats_in.csv")
    sdrf_file = find_file(sdrf_dir, "*.sdrf.tsv")
    mzml_stats = spectra_dir / "mzml_statistics"

    return mztab_file, msstats_file, sdrf_file, mzml_stats


def _validate_required_files(mztab_file, msstats_file, sdrf_file, mzml_stats):
    """Validate that all required files exist."""
    if not all([mztab_file, msstats_file, sdrf_file, mzml_stats.exists()]):
        missing = []
        if not mztab_file:
            missing.append("mzTab file")
        if not msstats_file:
            missing.append("MSstats input file")
        if not sdrf_file:
            missing.append("SDRF file")
        if not mzml_stats.exists():
            missing.append("mzML statistics")
        raise click.UsageError(f"ERROR: Missing required files: {', '.join(missing)}")

    print("\nFound input files:")
    print(f"   - mzTab file: {mztab_file}")
    print(f"   - MSstats file: {msstats_file}")
    print(f"   - SDRF file: {sdrf_file}")
    print(f"   - mzML statistics: {mzml_stats}")


def _initialize_project(
    output_folder_path: Path,
    project_accession: str,
    sdrf_file: Path,
    qpx_version: str,
    quantms_version: str,
):
    """Initialize the project with metadata."""
    print("\n=== Initializing Project ===")
    try:
        project_handler = check_directory(str(output_folder_path), project_accession)
        project_handler.populate_from_pride_archive()
        project_handler.populate_from_sdrf(str(sdrf_file))
        project_handler.add_qpx_version(qpx_version=qpx_version)
        project_handler.add_software_provider(
            sortware_name="quantms", sortware_version=quantms_version
        )
        # Save initial project file
        project_json = str(output_folder_path / f"{project_accession}.project.json")
        project_handler.save_project_info(
            output_prefix_file=project_accession,
            output_folder=str(output_folder_path),
            delete_existing=True,
        )
        print("Project initialization completed successfully")
        return project_handler
    except Exception as e:
        print(f"ERROR: Project initialization failed: {str(e)}", file=sys.stderr)
        raise


def _convert_features(
    mztab_file: Path,
    sdrf_file: Path,
    output_folder_path: Path,
    project_accession: str,
    generate_ibaq_view: bool,
) -> list:
    """Convert features and optionally generate IBAQ view."""
    created_files = []

    print("\n=== Starting Feature Conversion ===")
    feature_file = output_folder_path / create_uuid_filename(
        project_accession, ".feature.parquet"
    )
    convert_quantms_feature_cmd.callback(
        input_file=mztab_file,
        output_file=feature_file,
        sdrf_file=sdrf_file,
        verbose=True,
    )

    if feature_file and feature_file.exists():
        created_files.append(("feature-file", str(feature_file)))
        print("Feature conversion completed successfully")

        # Generate IBAQ view if requested
        if generate_ibaq_view:
            _generate_ibaq_view(
                sdrf_file, feature_file, project_accession, output_folder_path
            )
    else:
        print("ERROR: Feature conversion failed: No output file was generated")

    return created_files


def _generate_ibaq_view(
    sdrf_file: Path,
    feature_file: Path,
    project_accession: str,
    output_folder_path: Path,
):
    """Generate IBAQ view from feature data."""
    print("\n=== Generating IBAQ View ===")
    try:
        ibaq_file = create_uuid_filename(project_accession, ".ibaq.parquet")
        ibaq_path = output_folder_path / ibaq_file
        write_ibaq_feature(str(sdrf_file), str(feature_file), str(ibaq_path))
        print("IBAQ view generation completed successfully")
    except Exception as e:
        print(f"ERROR: IBAQ view generation failed: {str(e)}", file=sys.stderr)


def _convert_psms(
    mztab_file: Path, sdrf_file: Path, output_folder_path: Path, project_accession: str
) -> list:
    """Convert PSMs."""
    created_files = []

    print("\n=== Starting PSM Conversion ===")
    psm_file = output_folder_path / create_uuid_filename(
        project_accession, ".psm.parquet"
    )
    convert_quantms_psm_cmd.callback(
        input_file=mztab_file,
        output_file=psm_file,
        sdrf_file=sdrf_file,
        verbose=True,
    )

    if psm_file and psm_file.exists():
        created_files.append(("psm-file", str(psm_file)))
        print("PSM conversion completed successfully")

    return created_files


def _register_files_in_project(
    created_files: list, output_folder_path: Path, project_accession: str
):
    """Register all created files in the project."""
    print("\n=== Registering Files in Project ===")
    project_json = str(output_folder_path / f"{project_accession}.project.json")

    for file_category, file_path in created_files:
        try:
            attach_file_to_json_cmd.callback(
                project_file=project_json,
                attach_file=file_path,
                category=file_category,
                is_folder=False,
                partitions=None,
                replace_existing=True,
            )
            print(f"Registered {file_category}: {file_path}")
        except Exception as e:
            print(
                f"ERROR: Failed to register {file_category}: {str(e)}", file=sys.stderr
            )


def qpx_workflow(
    base_folder: str,
    output_folder: str,
    project_accession: str,
    quantms_version: Optional[str] = None,
    qpx_version: Optional[str] = None,
    generate_ibaq_view: bool = False,
) -> None:
    """Convert quantms output to QPX format.

    Expected structure:
    base_folder/
        quant_tables/
            *.mzTab
            *msstats_in.csv
        sdrf/
            *.sdrf.tsv
        spectra/
            mzml_statistics/
    """
    print("\n=== Starting QPX Conversion Workflow ===")

    # Setup and validate paths
    mztab_file, msstats_file, sdrf_file, mzml_stats = _setup_project_paths(base_folder)
    _validate_required_files(mztab_file, msstats_file, sdrf_file, mzml_stats)
    print(f"\nUsing project accession: {project_accession}")

    # Create output directory
    output_folder_path = Path(output_folder).resolve()
    check_dir(str(output_folder_path))
    print(f"\nOutput directory: {output_folder_path}")

    # Initialize project
    project_handler = _initialize_project(
        output_folder_path,
        project_accession,
        sdrf_file,
        qpx_version,
        quantms_version,
    )

    created_files = []

    try:
        # Convert features
        feature_files = _convert_features(
            mztab_file,
            sdrf_file,
            output_folder_path,
            project_accession,
            generate_ibaq_view,
        )
        created_files.extend(feature_files)

        # Convert PSMs
        psm_files = _convert_psms(
            mztab_file, sdrf_file, output_folder_path, project_accession
        )
        created_files.extend(psm_files)

        # Register all created files in the project
        _register_files_in_project(created_files, output_folder_path, project_accession)

    except Exception as e:
        print(f"ERROR: Conversion failed: {str(e)}", file=sys.stderr)

    print("\n=== Conversion Workflow Complete ===\n")


@click.command(
    "quantms",
    short_help="Convert quantms project output to QPX format",
)
@click.option(
    "--quantms-dir",
    help="The quantms project directory containing quant_tables, sdrf, and spectra subdirectories",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--output-dir",
    help="Output directory for QPX files (defaults to 'QPX' in parent directory)",
    required=False,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--project-accession",
    help="PRIDE project accession (e.g. 'PXD000865')",
    required=True,
    type=str,
)
@click.option(
    "--quantms-version",
    help="Version of quantms used to generate the data",
    required=True,
    type=str,
)
@click.option(
    "--qpx-version",
    help="Version of QPX used for conversion",
    required=True,
    type=str,
)
@click.option(
    "--generate-ibaq-view",
    help="Generate IBAQ view from feature data",
    is_flag=True,
    default=False,
)
def convert_quantms_project_cmd(
    quantms_dir: Path,
    output_dir: Optional[Path] = None,
    project_accession: str = None,
    quantms_version: str = None,
    qpx_version: str = None,
    generate_ibaq_view: bool = False,
) -> None:
    """Convert a quantms project output to QPX format.

    The script expects a quantms output directory with:
    - quant_tables/ containing mzTab and MSstats files
    - sdrf/ containing SDRF files
    - spectra/ containing mzML statistics
    """
    # Default output to sibling QPX directory
    if not output_dir:
        output_dir = str(quantms_dir.parent / "QPX")

    qpx_workflow(
        str(quantms_dir),
        output_dir,
        project_accession,
        quantms_version,
        qpx_version,
        generate_ibaq_view,
    )
