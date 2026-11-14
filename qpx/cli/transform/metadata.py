"""Transform CLI metadata implementation."""

import logging
from pathlib import Path
from typing import Optional

import click


@click.command("update-metadata")
@click.option(
    "--dataset",
    help="Path to a QPX dataset directory",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--sdrf",
    help="Path to the updated SDRF TSV file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--old-sdrf",
    help="Path to the original SDRF (for safety checks). If omitted, protected-field validation is skipped.",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
)
@click.option(
    "--force",
    help="Apply changes even if protected fields (instrument, enzymes, "
    "modifications, fractions, labels) have changed. Use with caution.",
    is_flag=True,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def transform_update_metadata_cmd(
    dataset: Path,
    sdrf: Path,
    old_sdrf: Optional[Path],
    force: bool,
    verbose: bool,
):
    """Update sample and run metadata from a revised SDRF file.

    Re-generates sample.parquet and run.parquet from the new SDRF.
    Original files are backed up as *.parquet.bak before overwriting.

    Safe to update (metadata-only, no impact on data):
        disease, organism_part, cell_type, cell_line, sex, age,
        treatment, individual, ancestry, developmental_stage,
        and any additional characteristics[*] columns.

    Protected fields (will BLOCK unless --force is used):
        instrument, enzymes, modification_parameters, fraction,
        dissociation_method, label/channel mapping, data file references.

    If --old-sdrf is provided, the tool compares old vs new SDRF and
    blocks if any protected fields changed. Without --old-sdrf, the
    safety check is skipped (useful for first-time metadata enrichment).

    \b
    Examples:
        # Update with safety check
        qpxc transform update-metadata \\
            --dataset ./my_project \\
            --sdrf ./updated_sdrf.tsv \\
            --old-sdrf ./original_sdrf.tsv

    \b
        # Update without safety check (first-time enrichment)
        qpxc transform update-metadata \\
            --dataset ./my_project \\
            --sdrf ./enriched_sdrf.tsv

    \b
        # Force update even if protected fields changed
        qpxc transform update-metadata \\
            --dataset ./my_project \\
            --sdrf ./new_sdrf.tsv \\
            --old-sdrf ./old_sdrf.tsv --force
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    from qpx.transforms.sdrf_updater import update_metadata

    try:
        result = update_metadata(
            dataset_path=dataset,
            new_sdrf_path=sdrf,
            old_sdrf_path=old_sdrf,
            force=force,
        )
    except ValueError as e:
        raise click.ClickException(str(e)) from e

    for w in result["warnings"]:
        if w.startswith("BLOCKED"):
            click.echo(f"  ⚠ {w}", err=True)
        elif w.startswith("WARNING"):
            click.echo(f"  ⚠ {w}", err=True)

    click.echo(f"\nDone — updated {result['samples_updated']} samples, {result['runs_updated']} runs in {dataset}")
