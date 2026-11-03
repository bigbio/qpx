"""Transform CLI relayout implementation."""

from __future__ import annotations

import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional

import click

from qpx.cli.resources import duckdb_resource_options

from .dataset import _copy_dataset_files, _gene_map_dataset_prefix, _stage_dataset_integrity


def _sort_keys(view: str, sort_by: Optional[str]) -> tuple[str, ...]:
    """The columns a view is sorted by: ``--sort-by``, or the view's default."""
    from qpx.transforms.relayout import DEFAULT_SORT_KEYS

    return tuple(column.strip() for column in sort_by.split(",")) if sort_by else DEFAULT_SORT_KEYS[view]


@click.command("relayout")
@click.option(
    "--dataset",
    required=True,
    help="Path to a QPX dataset directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option("--in-place", is_flag=True, help="Replace the dataset's files instead of writing to --output-folder")
@click.option(
    "--output-folder",
    default=None,
    help="Write a sorted copy here (must be empty and outside the dataset)",
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option("--view", "views", multiple=True, help="View to sort; repeatable (default: every one of feature, pg, psm, mz)")
@click.option("--sort-by", default=None, help="Comma-separated sort columns; needs exactly one --view")
@click.option("--row-group-size", default=None, type=click.IntRange(min=1), help="Rows per row group (default: 131072)")
@click.option("--compression-level", default=None, type=int, help="Codec level (default: 9 for zstd and gzip)")
@duckdb_resource_options
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def transform_relayout_cmd(**options):
    """Sort a dataset's views for selective reads, for publishing.

    Rows are sorted within each view, written in row groups with a page index,
    and the order is recorded in the footer (``sorting_columns`` and
    ``qpx.sort_order``). By default feature is sorted by protein, peptidoform,
    charge and run, pg by protein, its runs and label, and psm by peptidoform,
    charge, run and scan, which suits protein and peptide lookups; mz is sorted
    by run then scan. Values and ids do not change.

    \b
    Example:
        qpxc transform relayout --dataset ./qpx_output --output-folder ./qpx_sorted
    """
    from qpx.dataset import Dataset
    from qpx.transforms.relayout import DEFAULT_SORT_KEYS, relayout_view
    from qpx.writers.encoding import ROW_GROUP_ROWS

    if options["verbose"]:
        logging.getLogger().setLevel(logging.DEBUG)
    dataset, output_folder = options["dataset"].resolve(), options["output_folder"]
    if options["in_place"] == (output_folder is not None):
        raise click.UsageError("Specify exactly one of --in-place and --output-folder")
    if options["sort_by"] and len(options["views"]) != 1:
        raise click.UsageError("--sort-by needs exactly one --view")
    prefix = _gene_map_dataset_prefix(dataset, command="relayout")
    out_dir = dataset if options["in_place"] else output_folder.resolve()
    if not options["in_place"]:
        if dataset == out_dir or dataset in out_dir.parents:
            raise click.UsageError("--output-folder must be outside the dataset")
        _copy_dataset_files(dataset, out_dir)

    views = options["views"] or [view for view in DEFAULT_SORT_KEYS if (dataset / f"{prefix}.{view}.parquet").is_file()]
    with TemporaryDirectory(prefix=".relayout_tmp-", dir=out_dir) as temporary:
        staging = Path(temporary)
        moves = []
        for view in views:
            name = f"{prefix}.{view}.parquet"
            if not (dataset / name).is_file():
                raise click.ClickException(f"{dataset} has no {name}")
            try:
                relayout_view(
                    dataset / name,
                    staging / name,
                    _sort_keys(view, options["sort_by"]),
                    row_group_size=options["row_group_size"] or ROW_GROUP_ROWS,
                    compression_level=options["compression_level"],
                )
            except ValueError as exc:
                raise click.ClickException(str(exc)) from exc
            moves.append((staging / name, out_dir / name))
            click.echo(f"  {name}: sorted by {', '.join(_sort_keys(view, options['sort_by']))}")
        with Dataset(str(dataset), file_prefix=prefix) as qpx_dataset:
            metadata = _stage_dataset_integrity(qpx_dataset, prefix, staging, out_dir)
        if metadata is not None:
            moves.append(metadata)
        for source, target in moves:
            source.replace(target)
    click.echo(f"Relayout complete. Output: {out_dir}")
