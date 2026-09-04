"""
Validate subcommand -- validate QPX datasets and individual structures
against their canonical schemas.

Usage:
    qpxc validate --dataset-path ./PXD014414
    qpxc validate --dataset-path ./PXD014414 --structure feature
    qpxc validate --file ./data.feature.parquet
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

import click

from qpx.dataset import Dataset

logger = logging.getLogger("qpx.cli.validate")

# Derived from the canonical registry so a new structure cannot be supported
# by the library while the CLI silently rejects it (bigbio/qpx#289).
_VALID_STRUCTURES = Dataset.structure_names()


@click.command("validate")
@click.option(
    "--dataset-path",
    help="Path to a QPX dataset directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--file",
    "file_path",
    help="Path to a single QPX Parquet file to validate",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--structure",
    "structures",
    help="Structure(s) to validate (repeatable). Default: all.",
    type=click.Choice(_VALID_STRUCTURES, case_sensitive=False),
    multiple=True,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def validate_cmd(
    dataset_path: Optional[Path],
    file_path: Optional[Path],
    structures: tuple[str, ...],
    verbose: bool,
):
    """Validate a QPX dataset or structure against the canonical schema.

    Checks column presence, type matching, null values in required columns,
    and primary key uniqueness.

    \b
    Examples:
        # Validate an entire dataset
        qpxc validate --dataset-path ./PXD014414

        # Validate specific structures
        qpxc validate --dataset-path ./PXD014414 --structure feature
        qpxc validate --dataset-path ./PXD014414 --structure feature --structure pg

        # Validate a single Parquet file
        qpxc validate --file ./data.feature.parquet
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if file_path:
        _validate_single_file(file_path)
        return

    if dataset_path:
        _validate_dataset(dataset_path, list(structures) if structures else None)
        return

    raise click.ClickException("Provide --dataset-path or --file.")


def _validate_single_file(file_path: Path):
    """Validate a single Parquet file by inferring its structure type."""
    from qpx.dataset import Dataset

    for _, (cls, suffix) in Dataset._STRUCTURE_REGISTRY.items():
        if file_path.name.endswith(suffix):
            struct = cls.from_file(file_path)
            result = struct.validate(strict=True)
            _print_result(result)
            sys.exit(0 if result.is_valid else 1)

    raise click.ClickException(
        f"Cannot infer structure type from filename '{file_path.name}'. "
        "Expected a filename ending with one of: " + ", ".join(suffix for _, suffix in Dataset._STRUCTURE_REGISTRY.values())
    )


def _validate_dataset(dataset_path: Path, structures: list[str] | None):
    """Validate a dataset directory."""
    import qpx

    with qpx.open_dataset(dataset_path) as ds:
        results = ds.validate(structures=structures, strict=True)

    all_valid = True
    for _, result in results.items():
        _print_result(result)
        if not result.is_valid:
            all_valid = False

    click.echo()
    if all_valid:
        click.echo("Overall: VALID")
    else:
        click.echo("Overall: INVALID")
    sys.exit(0 if all_valid else 1)


def _print_result(result):
    """Print a single ValidationResult to stdout."""
    click.echo(f"\n{result.summary}")
    for issue in result.issues:
        marker = "ERROR" if issue.severity == "error" else "WARN "
        col_info = f" [{issue.column}]" if issue.column else ""
        click.echo(f"  {marker}{col_info}: {issue.message}")
