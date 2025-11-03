import os
from pathlib import Path

import click

from quantmsio.operate.report import ProjectReportGenerator

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(name="report", context_settings=CONTEXT_SETTINGS)
def report_cmd() -> None:
    """Report generation commands for quantms.io data"""


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
    """Generate comprehensive statistical report for a quantms.io project.
    
    Analyzes IBAQ parquet files and generates detailed statistics about samples, proteins, 
    and intensity distributions. Optionally includes feature and PSM statistics for enhanced analysis.
    Supports dimensionality reduction visualizations (PCA, t-SNE, UMAP).
    
    Generates three report formats: .txt, .json, and .html in the specified folder.
    
    Example:
        quantmsioc report generate \\
            --ibaq-file ./path/to/*-ibaq.parquet \\
            --feature-file ./path/to/*-feature.parquet \\
            --psm-file ./path/to/*-psm.parquet \\
            --output-folder ./reports
    """
    # Validate n_workers against CPU count
    max_workers = os.cpu_count()
    if n_workers > max_workers:
        click.echo(
            f"Error: n-workers ({n_workers}) exceeds available logical processors ({max_workers})",
            err=True,
        )
        raise click.Abort()

    # Create output folder if it doesn't exist
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    # Generate report name from IBAQ file
    base_name = (
        ibaq_file.stem.split("-")[0] if "-" in ibaq_file.stem else ibaq_file.stem
    )
    report_name = f"{base_name}-report"

    # Generate output paths
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
