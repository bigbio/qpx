"""
Query subcommands — run SQL, filter, and inspect QPX datasets.

These commands are fully functional because they only depend on
``qpx.Dataset`` and ``qpx.open_dataset()`` which are part of the core API.

Subcommands:
    qpxc query sql     — Run arbitrary SQL against a dataset
    qpxc query filter  — Filter a data structure by condition
    qpxc query head    — Show first N rows of a structure
"""

import logging
import sys
from pathlib import Path
from typing import Optional

import click

from qpx.dataset import Dataset

logger = logging.getLogger("qpx.cli.query")

# Shared option for the dataset path
_dataset_path_option = click.option(
    "--dataset-path",
    help="Path to a QPX dataset directory",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)

# Valid QPX data structure names
# Derived from the canonical registry so a new structure cannot be supported
# by the library while the CLI silently rejects it (bigbio/qpx#289).
_VALID_STRUCTURES = Dataset.structure_names()


# ---------------------------------------------------------------------------
# Query group
# ---------------------------------------------------------------------------


@click.group("query")
def query_cmd():
    """Query and inspect QPX datasets."""
    pass


# ---------------------------------------------------------------------------
# SQL
# ---------------------------------------------------------------------------


@query_cmd.command("sql")
@_dataset_path_option
@click.option(
    "--sql",
    "sql_query",
    help="SQL query to execute against the dataset",
    required=True,
    type=str,
)
@click.option(
    "--output",
    help="Output file path (CSV). If not specified, results go to stdout.",
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--output-format",
    help="Output format for results",
    type=click.Choice(["csv", "tsv", "parquet", "json"]),
    default="csv",
)
@click.option(
    "--structures",
    help="Comma-separated list of structures to load (default: all found)",
)
@click.option(
    "--duckdb-memory",
    help="DuckDB memory limit (e.g., '16GB')",
    default="16GB",
)
@click.option(
    "--duckdb-threads",
    help="DuckDB thread count",
    default=4,
    type=int,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def query_sql_cmd(
    dataset_path: Path,
    sql_query: str,
    output: Optional[Path],
    output_format: str,
    structures: Optional[str],
    duckdb_memory: str,
    duckdb_threads: int,
    verbose: bool,
):
    """Run an arbitrary SQL query against a QPX dataset.

    Opens the dataset, registers all (or selected) data structures as
    DuckDB tables, then executes the provided SQL. Table names match
    the QPX structure names: psm, feature, pg, mz, sample, run, etc.

    \b
    Examples:
        # Count features per run
        qpxc query sql \\
            --dataset-path ./PXD014414 \\
            --sql "SELECT run_file_name, COUNT(*) AS n FROM feature GROUP BY 1"

        # Join feature with run, export to CSV
        qpxc query sql \\
            --dataset-path ./PXD014414 \\
            --sql "SELECT f.sequence, r.run_file_name FROM feature f JOIN run r USING (run_file_name)" \\
            --output results.csv

        # Query specific structures only
        qpxc query sql \\
            --dataset-path ./PXD014414 \\
            --sql "SELECT COUNT(*) FROM psm" \\
            --structures psm
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    import qpx

    structure_list = [s.strip() for s in structures.split(",")] if structures else None

    with qpx.open_dataset(
        dataset_path,
        structures=structure_list,
        duckdb_memory=duckdb_memory,
        duckdb_threads=duckdb_threads,
    ) as ds:
        logger.info(f"Opened dataset: {ds}")
        logger.info(f"Executing SQL: {sql_query}")

        result = ds.sql(sql_query)
        df = result.to_df()

        if output:
            _write_output(df, output, output_format)
            click.echo(f"Query results written to: {output} ({len(df)} rows)")
        else:
            _print_dataframe(df, output_format)


# ---------------------------------------------------------------------------
# Filter
# ---------------------------------------------------------------------------


@query_cmd.command("filter")
@_dataset_path_option
@click.option(
    "--structure",
    help="QPX data structure to filter (e.g., feature, psm, pg)",
    required=True,
    type=click.Choice(_VALID_STRUCTURES, case_sensitive=False),
)
@click.option(
    "--condition",
    help="SQL WHERE condition (e.g., 'charge > 2 AND global_qvalue < 0.01')",
    required=True,
    type=str,
)
@click.option(
    "--columns",
    help="Comma-separated list of columns to include in output (default: all)",
)
@click.option(
    "--limit",
    help="Maximum number of rows to return",
    type=int,
)
@click.option(
    "--output",
    help="Output file path. If not specified, results go to stdout.",
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--output-format",
    help="Output format for results",
    type=click.Choice(["csv", "tsv", "parquet", "json"]),
    default="csv",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def query_filter_cmd(
    dataset_path: Path,
    structure: str,
    condition: str,
    columns: Optional[str],
    limit: Optional[int],
    output: Optional[Path],
    output_format: str,
    verbose: bool,
):
    """Filter a QPX data structure by a SQL condition.

    Opens the dataset, applies the filter condition to the specified
    structure, and returns the matching rows.

    \b
    Examples:
        # Filter features by charge and q-value
        qpxc query filter \\
            --dataset-path ./PXD014414 \\
            --structure feature \\
            --condition "charge > 2 AND global_qvalue < 0.01"

        # Filter PSMs, select specific columns
        qpxc query filter \\
            --dataset-path ./PXD014414 \\
            --structure psm \\
            --condition "is_decoy = false" \\
            --columns "sequence,charge,score" \\
            --limit 100

        # Export filtered protein groups to parquet
        qpxc query filter \\
            --dataset-path ./PXD014414 \\
            --structure pg \\
            --condition "global_qvalue < 0.01" \\
            --output filtered_pg.parquet \\
            --output-format parquet
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    import qpx

    with qpx.open_dataset(dataset_path, structures=[structure]) as ds:
        structure_attr = "dataset_meta" if structure == "dataset" else structure
        struct = getattr(ds, structure_attr, None)
        if struct is None:
            raise click.ClickException(
                f"Structure '{structure}' not found in dataset at {dataset_path}. Available: {ds.available_structures}"
            )

        logger.info(f"Filtering {structure}: {condition}")
        query = struct.filter(condition)

        if columns:
            col_list = [c.strip() for c in columns.split(",")]
            query = query.select(*col_list)

        if limit:
            query = query.limit(limit)

        df = query.to_df()

        if output:
            _write_output(df, output, output_format)
            click.echo(f"Filtered {structure}: {len(df)} rows written to {output}")
        else:
            _print_dataframe(df, output_format)
            click.echo(f"\n--- {len(df)} rows ---", err=True)


# ---------------------------------------------------------------------------
# Head
# ---------------------------------------------------------------------------


@query_cmd.command("head")
@_dataset_path_option
@click.option(
    "--structure",
    help="QPX data structure to inspect (e.g., feature, psm, pg)",
    required=True,
    type=click.Choice(_VALID_STRUCTURES, case_sensitive=False),
)
@click.option(
    "-n",
    "--rows",
    help="Number of rows to display",
    default=10,
    type=int,
)
@click.option(
    "--columns",
    help="Comma-separated list of columns to include (default: all)",
)
@click.option(
    "--output-format",
    help="Output format for results",
    type=click.Choice(["csv", "tsv", "table"]),
    default="table",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def query_head_cmd(
    dataset_path: Path,
    structure: str,
    rows: int,
    columns: Optional[str],
    output_format: str,
    verbose: bool,
):
    """Show the first N rows of a QPX data structure.

    A quick way to peek at the contents of any QPX data structure
    within a dataset.

    \b
    Examples:
        # Show first 10 rows of features
        qpxc query head \\
            --dataset-path ./PXD014414 \\
            --structure feature

        # Show 5 PSM rows with specific columns
        qpxc query head \\
            --dataset-path ./PXD014414 \\
            --structure psm \\
            -n 5 \\
            --columns "sequence,charge,score"
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    import qpx

    with qpx.open_dataset(dataset_path, structures=[structure]) as ds:
        structure_attr = "dataset_meta" if structure == "dataset" else structure
        struct = getattr(ds, structure_attr, None)
        if struct is None:
            raise click.ClickException(
                f"Structure '{structure}' not found in dataset at {dataset_path}. Available: {ds.available_structures}"
            )

        query = struct.limit(rows)
        if columns:
            col_list = [c.strip() for c in columns.split(",")]
            query = query.select(*col_list)

        df = query.to_df()

        if output_format == "table":
            click.echo(df.to_string(index=False))
        elif output_format == "csv":
            click.echo(df.to_csv(index=False))
        elif output_format == "tsv":
            click.echo(df.to_csv(index=False, sep="\t"))

        total = struct.count()
        click.echo(
            f"\n--- Showing {len(df)} of {total} rows in {structure} ---",
            err=True,
        )


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------


def _write_output(df, path: Path, fmt: str):
    """Write a DataFrame to a file in the requested format."""
    path.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "csv":
        df.to_csv(path, index=False)
    elif fmt == "tsv":
        df.to_csv(path, index=False, sep="\t")
    elif fmt == "parquet":
        import pyarrow as pa
        import pyarrow.parquet as pq

        table = pa.Table.from_pandas(df)
        pq.write_table(table, str(path), compression="zstd")
    elif fmt == "json":
        df.to_json(path, orient="records", indent=2)
    else:
        df.to_csv(path, index=False)


def _print_dataframe(df, fmt: str):
    """Print a DataFrame to stdout in the requested format."""
    if fmt == "csv":
        sys.stdout.write(df.to_csv(index=False))
    elif fmt == "tsv":
        sys.stdout.write(df.to_csv(index=False, sep="\t"))
    elif fmt == "json":
        sys.stdout.write(df.to_json(orient="records", indent=2))
    else:
        sys.stdout.write(df.to_csv(index=False))
