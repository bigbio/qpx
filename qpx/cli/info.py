"""
Info subcommands — show dataset summary, schemas, and Parquet metadata.

These commands are fully functional because they only depend on
``qpx.Dataset``, ``qpx.open_dataset()``, and ``qpx.core.parquet_io`` which
are part of the core API.

Subcommands:
    qpxc info              — Show dataset summary (structures, row counts)
    qpxc info schema       — Show schema for a data structure
    qpxc info metadata     — Show Parquet footer metadata
"""

import logging
from pathlib import Path
from typing import Optional

import click

from qpx.dataset import Dataset

logger = logging.getLogger("qpx.cli.info")

# Derived from the canonical registry so a new structure cannot be supported
# by the library while the CLI silently rejects it (bigbio/qpx#289).
_VALID_STRUCTURES = Dataset.structure_names()


# ---------------------------------------------------------------------------
# Info group
# ---------------------------------------------------------------------------


@click.group(invoke_without_command=True)
@click.option(
    "--dataset-path",
    help="Path to a QPX dataset directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
@click.pass_context
def info(ctx, dataset_path: Optional[Path], verbose: bool):
    """Show information about a QPX dataset.

    When invoked without a subcommand, displays a summary of the
    dataset including available structures and row counts.

    \b
    Examples:
        # Show dataset summary
        qpxc info --dataset-path ./PXD014414

    \b
        # Show schema for a specific structure
        qpxc info schema --dataset-path ./PXD014414 --structure feature

    \b
        # Show Parquet footer metadata
        qpxc info metadata --file ./PXD014414/data.feature.parquet
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Store options for subcommands
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["dataset_path"] = dataset_path

    # If no subcommand given, show the dataset summary
    if ctx.invoked_subcommand is None:
        if dataset_path is None:
            raise click.ClickException("Please provide --dataset-path or use a subcommand (schema, metadata).")
        _show_dataset_summary(dataset_path, verbose)


# ---------------------------------------------------------------------------
# Default: dataset summary
# ---------------------------------------------------------------------------


def _show_dataset_summary(dataset_path: Path, verbose: bool):
    """Print a human-readable summary of a QPX dataset."""
    import qpx

    click.echo(f"QPX Dataset: {dataset_path.resolve()}")
    click.echo("=" * 60)

    with qpx.open_dataset(dataset_path) as ds:
        structures = ds.available_structures
        if not structures:
            click.echo("  (no QPX data structures found)")
            return

        click.echo(f"  Structures found: {len(structures)}")
        click.echo()

        # Table header
        click.echo(f"  {'Structure':<15} {'Rows':>12} {'File'}")
        click.echo(f"  {'-' * 15} {'-' * 12} {'-' * 40}")

        total_rows = 0
        for name in sorted(structures):
            struct = ds._structures[name]
            # Counts every shard or partition of the view, not just the first file.
            rows = struct.count()

            total_rows += rows
            rel_path = struct.file_paths[0].name
            click.echo(f"  {name:<15} {rows:>12,} {rel_path}")

        click.echo(f"  {'-' * 15} {'-' * 12}")
        click.echo(f"  {'TOTAL':<15} {total_rows:>12,}")
        click.echo()

        # Show brief file metadata if available
        if verbose:
            for name in sorted(structures):
                struct = ds._structures[name]
                meta = struct.file_metadata
                if meta:
                    click.echo(f"  [{name}] Parquet metadata:")
                    for k, v in sorted(meta.items()):
                        # Skip pandas metadata (very verbose)
                        if k == "pandas":
                            continue
                        display_val = v if len(v) < 80 else v[:77] + "..."
                        click.echo(f"    {k}: {display_val}")
                    click.echo()


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------


@info.command("schema")
@click.option(
    "--dataset-path",
    help="Path to a QPX dataset directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--file",
    "file_path",
    help="Path to a single Parquet file (alternative to --dataset-path)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--structure",
    help="QPX data structure name (when using --dataset-path)",
    type=click.Choice(_VALID_STRUCTURES, case_sensitive=False),
)
@click.option(
    "--canonical",
    help="Show the canonical QPX schema (from code) instead of the file schema",
    is_flag=True,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
@click.pass_context
def info_schema_cmd(
    ctx,
    dataset_path: Optional[Path],
    file_path: Optional[Path],
    structure: Optional[str],
    canonical: bool,
    verbose: bool,
):
    """Show the Arrow schema for a QPX data structure.

    Can read the schema from an on-disk Parquet file or display the
    canonical QPX schema as defined in code.

    \b
    Examples:
        # Schema from a dataset structure
        qpxc info schema --dataset-path ./PXD014414 --structure feature

    \b
        # Schema from a standalone Parquet file
        qpxc info schema --file ./data.feature.parquet

    \b
        # Canonical QPX schema (from code)
        qpxc info schema --structure feature --canonical
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Resolve dataset_path from parent context if not given
    if dataset_path is None:
        dataset_path = ctx.obj.get("dataset_path")

    if canonical and structure:
        _show_canonical_schema(structure)
        return

    if file_path:
        _show_file_schema(file_path)
        return

    if dataset_path and structure:
        import qpx

        with qpx.open_dataset(dataset_path, structures=[structure]) as ds:
            struct = getattr(ds, structure, None)
            if struct is None:
                raise click.ClickException(
                    f"Structure '{structure}' not found in dataset at {dataset_path}. Available: {ds.available_structures}"
                )
            _show_file_schema(struct._file_path)
        return

    if dataset_path and not structure:
        # Show schemas for all structures in the dataset
        import qpx

        with qpx.open_dataset(dataset_path) as ds:
            for name in sorted(ds.available_structures):
                struct = ds._structures[name]
                click.echo(f"=== {name} ===")
                _show_file_schema(struct._file_path)
                click.echo()
        return

    raise click.ClickException("Provide --dataset-path (+ optional --structure), --file, or --structure --canonical.")


def _show_file_schema(file_path: Path):
    """Read and display the Arrow schema from a Parquet file."""
    from qpx.core.parquet_io import read_parquet_schema

    schema = read_parquet_schema(file_path)
    click.echo(f"File: {file_path}")
    click.echo(f"Fields: {len(schema)}")
    click.echo()

    click.echo(f"  {'#':<4} {'Name':<35} {'Type'}")
    click.echo(f"  {'-' * 4} {'-' * 35} {'-' * 30}")
    for i, field in enumerate(schema):
        click.echo(f"  {i:<4} {field.name:<35} {field.type}")


def _show_canonical_schema(structure: str):
    """Display the canonical QPX schema from its YAML definition."""
    from qpx.core.data.loader import load_schema

    if structure not in _VALID_STRUCTURES:
        raise click.ClickException(f"Unknown structure: {structure}")

    try:
        schema = load_schema(structure).get_arrow_schema()
    except (AttributeError, TypeError) as e:
        raise click.ClickException(f"Cannot retrieve canonical schema for '{structure}': {e}")

    click.echo(f"Canonical QPX schema: {structure}")
    click.echo(f"Fields: {len(schema)}")
    click.echo()

    click.echo(f"  {'#':<4} {'Name':<35} {'Type':<25} {'Nullable'}")
    click.echo(f"  {'-' * 4} {'-' * 35} {'-' * 25} {'-' * 8}")
    for i, field in enumerate(schema):
        nullable = "yes" if field.nullable else "no"
        click.echo(f"  {i:<4} {field.name:<35} {str(field.type):<25} {nullable}")


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------


@info.command("metadata")
@click.option(
    "--file",
    "file_path",
    help="Path to a Parquet file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--dataset-path",
    help="Path to a QPX dataset directory (shows metadata for all files)",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
@click.pass_context
def info_metadata_cmd(
    ctx,
    file_path: Optional[Path],
    dataset_path: Optional[Path],
    verbose: bool,
):
    """Show Parquet footer metadata for QPX files.

    Displays the key-value metadata stored in the Parquet file footer,
    as well as basic file statistics (row groups, total rows, file size).

    \b
    Examples:
        # Single file
        qpxc info metadata --file ./data.feature.parquet

    \b
        # All files in a dataset
        qpxc info metadata --dataset-path ./PXD014414
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Resolve dataset_path from parent context if not given
    if dataset_path is None:
        dataset_path = ctx.obj.get("dataset_path")

    if file_path:
        _show_parquet_metadata(file_path)
        return

    if dataset_path:
        import qpx

        with qpx.open_dataset(dataset_path) as ds:
            for name in sorted(ds.available_structures):
                struct = ds._structures[name]
                click.echo(f"=== {name} ===")
                _show_parquet_metadata(struct._file_path)
                click.echo()
        return

    raise click.ClickException("Provide --file or --dataset-path.")


def _show_parquet_metadata(file_path: Path):
    """Display Parquet file footer metadata."""
    import pyarrow.parquet as pq

    try:
        pf = pq.ParquetFile(file_path)
    except Exception as e:
        raise click.ClickException(f"Cannot open Parquet file {file_path}: {e}")

    metadata = pf.metadata
    file_size = file_path.stat().st_size

    click.echo(f"File:        {file_path}")
    click.echo(f"Size:        {_format_size(file_size)}")
    click.echo(f"Rows:        {metadata.num_rows:,}")
    click.echo(f"Row groups:  {metadata.num_row_groups}")
    click.echo(f"Columns:     {metadata.num_columns}")
    click.echo(f"Format:      {metadata.format_version}")
    click.echo(f"Created by:  {metadata.created_by}")
    click.echo()

    # Key-value metadata
    kv_meta = pf.schema_arrow.metadata or {}
    if kv_meta:
        click.echo("Key-value metadata:")
        for k, v in sorted(kv_meta.items()):
            key = k.decode() if isinstance(k, bytes) else k
            val = v.decode() if isinstance(v, bytes) else v
            # Skip pandas metadata (very verbose JSON)
            if key == "pandas":
                click.echo(f"  {key}: <pandas schema, {len(val)} chars>")
                continue
            display_val = val if len(val) < 100 else val[:97] + "..."
            click.echo(f"  {key}: {display_val}")
    else:
        click.echo("Key-value metadata: (none)")

    # Row group summary
    if metadata.num_row_groups > 1:
        click.echo()
        click.echo("Row groups:")
        click.echo(f"  {'#':<5} {'Rows':>12} {'Compressed':>15} {'Uncompressed':>15}")
        click.echo(f"  {'-' * 5} {'-' * 12} {'-' * 15} {'-' * 15}")
        for i in range(metadata.num_row_groups):
            rg = metadata.row_group(i)
            compressed = sum(rg.column(j).total_compressed_size for j in range(rg.num_columns))
            uncompressed = sum(rg.column(j).total_uncompressed_size for j in range(rg.num_columns))
            click.echo(f"  {i:<5} {rg.num_rows:>12,} {_format_size(compressed):>15} {_format_size(uncompressed):>15}")


def _format_size(size_bytes: int) -> str:
    """Format a byte count as a human-readable string."""
    if size_bytes < 1024:
        return f"{size_bytes} B"
    elif size_bytes < 1024 * 1024:
        return f"{size_bytes / 1024:.1f} KB"
    elif size_bytes < 1024 * 1024 * 1024:
        return f"{size_bytes / (1024 * 1024):.1f} MB"
    else:
        return f"{size_bytes / (1024 * 1024 * 1024):.2f} GB"
