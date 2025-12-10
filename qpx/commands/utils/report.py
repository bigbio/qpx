import logging
import os
from pathlib import Path
from typing import Dict, Optional, Set

import click
import pandas as pd

from qpx.operate.report import ProjectReportGenerator
from qpx.core.metadata import WorkflowMetadataGenerator
from qpx.utils.logger import get_logger

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


def _detect_workflow_files(
    input_path: Path, workflow: str, logger: logging.Logger
) -> Dict[str, Path]:
    """
    Detect workflow-specific files in a directory.

    Args:
        input_path: Path to file or directory
        workflow: Workflow type (maxquant, diann, etc.)
        logger: Logger instance

    Returns:
        Dictionary mapping model types to file paths
    """
    file_map = {}

    if input_path.is_file():
        model_type = _detect_file_model_type(input_path, workflow)
        if model_type:
            file_map[model_type] = input_path
        return file_map

    if workflow == "maxquant":
        for file_path in input_path.iterdir():
            if not file_path.is_file():
                continue
            filename = file_path.name.lower()
            if "msms" in filename and filename.endswith(".txt"):
                file_map["psm"] = file_path
                logger.debug(f"Found PSM file: {file_path.name}")
            elif "evidence" in filename and filename.endswith(".txt"):
                file_map["feature"] = file_path
                logger.debug(f"Found Feature file: {file_path.name}")
            elif "proteingroups" in filename and filename.endswith(".txt"):
                file_map["pg"] = file_path
                logger.debug(f"Found PG file: {file_path.name}")

    return file_map


def _detect_file_model_type(input_file: Path, workflow: str) -> Optional[str]:
    """
    Detect the model type based on input file name.

    Args:
        input_file: Path to input file
        workflow: Workflow type (maxquant, diann, etc.)

    Returns:
        Model type ('psm', 'feature', or 'pg'), or None if cannot detect
    """
    filename = input_file.name.lower()

    if workflow == "maxquant":
        if "msms" in filename:
            return "psm"
        elif "evidence" in filename:
            return "feature"
        elif "proteingroups" in filename or "protein_groups" in filename:
            return "pg"
    elif workflow == "diann":
        if "report" in filename:
            return "feature"

    return None


def _detect_available_columns(input_file: Path, logger: logging.Logger) -> Set[str]:
    """
    Detect available columns from input file.

    Args:
        input_file: Path to input file
        logger: Logger instance

    Returns:
        Set of column names found in the file
    """
    try:
        df_header = pd.read_csv(input_file, sep="\t", nrows=0)
        columns = set(df_header.columns)

        logger.debug(f"Successfully read header from {input_file.name}")
        return columns

    except Exception as e:
        logger.warning(f"Could not read columns from {input_file}: {e}")
        logger.warning("Falling back to generating all possible columns")
        return set()


def _validate_and_filter_file_map(
    file_map: Dict[str, Path],
    model: Optional[str],
    input_path: Path,
    logger: logging.Logger,
) -> Dict[str, Path]:
    """
    Validate file map and filter by model if specified.

    Args:
        file_map: Dictionary mapping model types to file paths
        model: Specific model to filter (optional)
        input_path: Original input path for error messages
        logger: Logger instance

    Returns:
        Filtered file map

    Raises:
        click.UsageError: If no files found or specified model not available
    """
    if not file_map:
        raise click.UsageError(
            f"No workflow files found in: {input_path}. "
            "For MaxQuant, expected files: msms.txt, evidence.txt, proteinGroups.txt"
        )

    logger.info(
        f"Found {len(file_map)} model file(s): {', '.join(file_map.keys()).upper()}"
    )

    if model:
        if model not in file_map:
            raise click.UsageError(
                f"Specified model '{model}' file not found in input. "
                f"Available models: {', '.join(file_map.keys())}"
            )
        file_map = {model: file_map[model]}
        logger.info(f"Filtering to specified model: {model.upper()}")

    return file_map


def _process_model_metadata(
    model_type: str,
    file_path: Path,
    workflow: str,
    logger: logging.Logger,
) -> list:
    """
    Process metadata for a single model type.

    Args:
        model_type: Type of model (psm, feature, pg)
        file_path: Path to the input file
        workflow: Workflow type
        logger: Logger instance

    Returns:
        List of metadata records
    """
    logger.info(f"Processing {model_type.upper()} from: {file_path.name}")

    available_columns = _detect_available_columns(file_path, logger)
    if available_columns is None:
        available_columns = set()
    logger.info(f"  Detected {len(available_columns)} columns")

    model_generator = WorkflowMetadataGenerator(
        workflow=workflow, available_columns=available_columns
    )

    if model_type == "psm":
        model_generator.generate_psm_metadata()
    elif model_type == "feature":
        model_generator.generate_feature_metadata()
    elif model_type == "pg":
        model_generator.generate_pg_metadata()

    return model_generator.get_records()


def _log_model_breakdown(records: list, logger: logging.Logger) -> None:
    """Log breakdown of records by model type."""
    model_counts: Dict[str, int] = {}
    for record in records:
        model = record["model_view"]
        model_counts[model] = model_counts.get(model, 0) + 1

    logger.debug(
        f"Model breakdown: {', '.join([f'{k.upper()}({v})' for k, v in model_counts.items()])}"
    )


@click.group(name="report", context_settings=CONTEXT_SETTINGS)
def report_cmd() -> None:
    """Report generation commands for QPX data"""


@report_cmd.command(
    "generate",
    short_help="Generate project report from ibaqpy pipeline results",
)
@click.option(
    "--ibaq-file",
    required=True,
    help="Path to IBAQ TSV or parquet file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--feature-file",
    help="Path to feature parquet file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--psm-file",
    help="Path to PSM parquet file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    required=True,
    help="Output folder for generated reports",
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--min-unique-peptides",
    default=2,
    type=int,
    help="Minimum unique peptides per protein (set 0 to disable)",
)
@click.option(
    "--memory-limit",
    default=None,
    type=float,
    help="Memory limit in GB for batch processing",
)
@click.option(
    "--n-workers",
    default=8,
    type=int,
    help="Number of parallel workers",
)
def generate_report_cmd(
    ibaq_file: Path,
    feature_file: Path,
    psm_file: Path,
    output_folder: Path,
    min_unique_peptides: int,
    memory_limit: float,
    n_workers: int,
) -> None:
    """Generate comprehensive statistical report for a QPX project.
    
    Analyzes IBAQ parquet files and generates detailed statistics about samples, proteins, 
    and intensity distributions. Optionally includes feature and PSM statistics for enhanced analysis.
    Supports dimensionality reduction visualizations (PCA, t-SNE, UMAP).
    
    Generates three report formats: .txt, .json, and .html in the specified folder.
    
    Example:
        qpxc report generate \\
            --ibaq-file ./path/to/*-ibaq.parquet \\
            --feature-file ./path/to/*-feature.parquet \\
            --psm-file ./path/to/*-psm.parquet \\
            --output-folder ./reports
    """
    max_workers = os.cpu_count()
    if n_workers > max_workers:
        click.echo(
            f"Error: n-workers ({n_workers}) exceeds available logical processors ({max_workers})",
            err=True,
        )
        raise click.Abort()

    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    base_name = (
        ibaq_file.stem.split("-")[0] if "-" in ibaq_file.stem else ibaq_file.stem
    )
    report_name = f"{base_name}-report"
    output_txt = output_folder / f"{report_name}.txt"
    output_json = output_folder / f"{report_name}.json"
    output_html = output_folder / f"{report_name}.html"

    generator = ProjectReportGenerator(
        ibaq_file=ibaq_file,
        feature_file=feature_file,
        psm_file=psm_file,
        min_unique_peptides=min_unique_peptides,
        memory_limit_gb=memory_limit,
        n_workers=n_workers,
    )

    if memory_limit:
        click.echo(f"Memory limit set to {memory_limit:.1f} GB")
    else:
        click.echo(f"Memory limit set to {generator.memory_limit_gb:.1f} GB")

    click.echo(f"Using {generator.n_workers} parallel workers")

    generator.compute_statistics()
    generator.generate_text_report(output_file=output_txt)
    generator.generate_json_report(output_file=output_json)
    generator.generate_html_report(output_file=output_html)
    click.echo(f"HTML report: {output_html}")

    click.echo("Report generation completed successfully")


@report_cmd.command(
    "metadata",
    short_help="Generate metadata configuration file for a workflow",
)
@click.option(
    "--workflow",
    required=True,
    type=click.Choice(
        ["maxquant", "diann", "quantms-lfq", "quantms-tmt", "quantms-psm"],
        case_sensitive=False,
    ),
    help="Workflow name",
)
@click.option(
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input file or folder path (folder mode auto-detects MaxQuant output files)",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output metadata configuration file path",
)
@click.option(
    "--model",
    type=click.Choice(["psm", "feature", "pg"], case_sensitive=False),
    help="Specific model to generate (optional, auto-detects if not specified)",
)
@click.option(
    "--verbose",
    is_flag=True,
    help="Enable verbose logging",
)
def generate_metadata_cmd(
    workflow: str,
    input_path: Path,
    output: Path,
    model: Optional[str] = None,
    verbose: bool = False,
) -> None:
    """Generate metadata configuration file documenting column mappings.
    
    This command generates a CSV file documenting column mappings between source 
    proteomics formats (MaxQuant, DIANN, quantms) and QPX format. Only columns 
    present in the input files will be documented. Supports both single file and 
    folder input modes for automatic detection of workflow-specific files.
    
    Example:
        qpxc report metadata \\
            --workflow maxquant \\
            --input ./maxquant_output \\
            --output metadata.csv
    """
    logger = get_logger("qpx.commands.report")
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    try:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using output file: {output}")

        input_path = Path(input_path)
        file_map = _detect_workflow_files(input_path, workflow, logger)
        file_map = _validate_and_filter_file_map(file_map, model, input_path, logger)

        logger.info(f"Initializing metadata generator for workflow: {workflow}")
        generator = WorkflowMetadataGenerator(workflow=workflow)

        for model_type, file_path in file_map.items():
            records = _process_model_metadata(model_type, file_path, workflow, logger)
            generator.metadata_records.extend(records)

        generator.generate_file(str(output))

        record_count = len(generator.get_records())
        logger.info(
            f"Successfully generated metadata configuration with {record_count:,} entries"
        )
        logger.info(f"Metadata configuration file saved to: {output}")

        if verbose and record_count > 0:
            _log_model_breakdown(generator.get_records(), logger)

    except Exception as e:
        logger.error(
            f"Error in metadata configuration generation: {str(e)}", exc_info=True
        )
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")
