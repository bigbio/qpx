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
import duckdb

from qpx.cli.resources import duckdb_resource_options
from qpx.core.data.validation import ValidationIssue, ValidationResult
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
@click.option(
    "--semantic",
    help="Also check PSM and feature peptidoforms (ProForma syntax, sequence, Unimod accessions), "
    "modification positions and calculated/observed m/z",
    is_flag=True,
)
@click.option(
    "--check-identity",
    help="Also re-derive each feature_id, psm_id and pg_id from the composite its file declares",
    is_flag=True,
)
@click.option(
    "--conformance",
    help="Run every check (--semantic and --check-identity included) and fail on warnings too",
    is_flag=True,
)
@duckdb_resource_options
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def validate_cmd(
    dataset_path: Optional[Path],
    file_path: Optional[Path],
    structures: tuple[str, ...],
    semantic: bool,
    check_identity: bool,
    conformance: bool,
    verbose: bool,
):
    """Validate a QPX dataset or structure against the canonical schema.

    Checks the columns, values and footer of each structure, and the
    references between structures.

    \b
    Examples:
        # Validate an entire dataset
        qpxc validate --dataset-path ./PXD014414

    \b
        # Validate specific structures
        qpxc validate --dataset-path ./PXD014414 --structure feature
        qpxc validate --dataset-path ./PXD014414 --structure feature --structure pg

    \b
        # Validate a single Parquet file
        qpxc validate --file ./data.feature.parquet

    \b
        # Also check peptidoforms, modification positions and m/z
        qpxc validate --dataset-path ./PXD014414 --semantic

    \b
        # Also re-derive the ids from their composites
        qpxc validate --dataset-path ./PXD014414 --check-identity

    \b
        # The gate for a producer's CI: every check, warnings fail too
        qpxc validate --dataset-path ./PXD014414 --conformance
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    semantic = semantic or conformance
    check_identity = check_identity or conformance
    if file_path:
        _validate_single_file(file_path, semantic, check_identity, conformance)
        return

    if dataset_path:
        _validate_dataset(dataset_path, list(structures) if structures else None, semantic, check_identity, conformance)
        return

    raise click.ClickException("Provide --dataset-path or --file.")


def _passes(result, conformance: bool) -> bool:
    return result.is_valid and not (conformance and result.issues)


def _validate_single_file(file_path: Path, semantic: bool, check_identity: bool, conformance: bool):
    """Validate a single Parquet file by inferring its structure type."""
    from qpx.dataset import Dataset

    for name, (cls, suffix) in Dataset._STRUCTURE_REGISTRY.items():
        if file_path.name.endswith(suffix):
            try:
                result = cls.from_file(file_path).validate(strict=True, semantic=semantic, check_identity=check_identity)
            except (OSError, ValueError, duckdb.Error) as exc:
                result = ValidationResult(structure=name)
                message = f"{file_path.name} could not be read: {exc}"
                result.issues.append(ValidationIssue(name, "unreadable_structure", "error", None, message))
            _print_result(result)
            sys.exit(0 if _passes(result, conformance) else 1)

    raise click.ClickException(
        f"Cannot infer structure type from filename '{file_path.name}'. "
        "Expected a filename ending with one of: " + ", ".join(suffix for _, suffix in Dataset._STRUCTURE_REGISTRY.values())
    )


def _validate_dataset(dataset_path: Path, structures: list[str] | None, semantic: bool, check_identity: bool, conformance: bool):
    """Validate a dataset directory."""
    import qpx

    with qpx.open_dataset(dataset_path) as ds:
        results = ds.validate(structures=structures, strict=True, semantic=semantic, check_identity=check_identity)

    all_valid = True
    for _, result in results.items():
        _print_result(result)
        if not _passes(result, conformance):
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
