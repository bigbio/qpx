"""Transform CLI gene map implementation."""

import logging
from pathlib import Path
from typing import Optional

import click

from qpx.cli.resources import duckdb_resource_options

from .dataset import _copy_dataset_files, _gene_map_dataset_prefix, _refresh_dataset_mudata, _stage_dataset_integrity


def _validate_gene_map_inputs(
    parquet_path: Optional[Path],
    dataset: Optional[Path],
    in_place: bool,
    output_folder: Optional[Path],
) -> None:
    """Reject input/destination combinations gene-map cannot act on."""
    if (parquet_path is None) == (dataset is None):
        raise click.UsageError("Specify exactly one of --parquet-path or --dataset")
    if dataset is None and output_folder is None:
        raise click.UsageError("--output-folder is required with --parquet-path")
    if dataset is not None and not in_place and output_folder is None:
        raise click.UsageError("Specify --in-place or --output-folder with --dataset")
    if dataset is not None and in_place and output_folder is not None:
        raise click.UsageError("Specify either --in-place or --output-folder with --dataset, not both")
    if dataset is not None and output_folder is not None:
        source, destination = dataset.resolve(), output_folder.resolve()
        if destination == source:
            raise click.UsageError("--output-folder is the dataset itself; use --in-place")
        if source in destination.parents:
            raise click.UsageError("--output-folder must not be inside the dataset directory")


@click.command("gene-map")
@click.option(
    "--parquet-path",
    help="QPX PSM, feature or pg parquet file path (single-file mode)",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--dataset",
    help="Path to a QPX dataset directory (annotates pg + feature, refreshes the h5mu)",
    default=None,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--in-place",
    help="Dataset mode: overwrite the dataset's files instead of writing to --output-folder",
    is_flag=True,
)
@click.option(
    "--fasta",
    help="FASTA database file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    default=None,
    type=click.Path(file_okay=False, path_type=Path),
)
@duckdb_resource_options
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def transform_gene_map_cmd(
    parquet_path: Optional[Path],
    dataset: Optional[Path],
    in_place: bool,
    fasta: Path,
    output_folder: Optional[Path],
    verbose: bool,
):
    """Map gene names from a FASTA file to QPX parquet data.

    Enriches protein identifications with gene names read from the ``GN=`` field
    of the FASTA headers. The FASTA is the only source: nothing is fetched over
    the network.

    Given ``--dataset``, every quantification view that carries protein
    accessions (pg and feature) is annotated, and the dataset's MuData view is
    rebuilt when one is present, so the ``.h5mu`` never describes stale genes.
    Dataset mode requires flat Parquet views sharing one file prefix;
    partitioned pg/feature views are rejected before any files are changed.
    An output folder must be empty and outside the source dataset.

    \b
    Examples:
        # A single parquet file
        qpxc transform gene-map \\
            --parquet-path ./output/maxquant.psm.parquet \\
            --fasta proteins.fasta \\
            --output-folder ./output

    \b
        # A whole QPX dataset, refreshing its h5mu in place
        qpxc transform gene-map \\
            --dataset ./qpx_output \\
            --fasta proteins.fasta \\
            --in-place
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    _validate_gene_map_inputs(parquet_path, dataset, in_place, output_folder)

    from qpx.transforms.gene_mapping import GeneMappingTransform

    mapping = GeneMappingTransform(fasta_path=str(fasta))

    if dataset is not None:
        _gene_map_dataset(mapping, dataset, in_place, output_folder)
    else:
        _gene_map_parquet(mapping, parquet_path, Path(output_folder))


def _gene_map_parquet(mapping, parquet_path: Path, output_folder: Path) -> None:
    """Annotate one file while keeping its physical types and Parquet metadata."""
    import datetime
    import uuid

    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq

    from qpx.writers.encoding import SampledParquetWriter, row_group_rows

    output_folder.mkdir(parents=True, exist_ok=True)

    with pq.ParquetFile(parquet_path) as source:
        table = source.read()
        # The table is rebuilt from a new DataFrame, so the source's pandas metadata no longer applies.
        footer = {
            key: value for key, value in (source.metadata.metadata or {}).items() if key not in (b"ARROW:schema", b"pandas")
        }
    if b"uuid" in footer:
        # Like the other QPX rewrites, the annotated file gets its own uuid and date.
        footer |= {b"uuid": str(uuid.uuid4()).encode(), b"creation_date": datetime.datetime.now().isoformat().encode()}
    schema = table.schema
    df = table.to_pandas(ignore_metadata=True, types_mapper=pd.ArrowDtype)
    protein_col = "pg_accessions" if "pg_accessions" in df.columns else "protein_accessions"
    annotated = mapping.annotate_dataframe(df, protein_col=protein_col)
    table = pa.Table.from_pandas(annotated, schema=schema, preserve_index=False)
    for column in annotated.columns.difference(schema.names, sort=False):
        table = table.append_column(column, pa.array(annotated[column]))

    compression = footer.get(b"compression_format", b"snappy").decode()
    output_path = output_folder / parquet_path.name
    table = table.replace_schema_metadata(footer)
    with SampledParquetWriter(output_path, table.schema, compression) as writer:
        writer.write_table(table, row_group_size=row_group_rows(table))

    click.echo(f"Gene mapping complete. Output: {output_path}")


def _annotate_dataset_views(
    mapping,
    qpx_dataset,
    prefix: str,
    staging: Path,
    out_dir: Path,
) -> list[tuple[Path, Path]]:
    """Write gene-annotated pg + feature views into ``staging``.

    Returns the ``(staged, destination)`` pairs to move into place.
    """
    written: list[tuple[Path, Path]] = []
    views = (
        ("pg", qpx_dataset.pg, mapping.write_annotated_pg),
        ("feature", qpx_dataset.feature, mapping.write_annotated_features),
    )
    for view, structure, write_view in views:
        if structure is None:
            continue
        name = f"{prefix}.{view}.parquet"
        write_view(qpx_dataset, staging / name)
        share = mapping.last_mapped_share
        if share is not None:
            click.echo(f"  {name}: {share:.1f}% of rows matched a gene in the FASTA")
            if share == 0.0:
                click.echo(f"  WARNING: the FASTA resolved no genes for {view}; existing gene names were kept unchanged")
        written.append((staging / name, out_dir / name))
    return written


def _gene_map_dataset(
    mapping,
    dataset: Path,
    in_place: bool,
    output_folder: Optional[Path],
) -> None:
    """Annotate a QPX dataset's pg + feature views and refresh its MuData view.

    Views are written to a temporary directory first, so an in-place run never
    truncates a parquet that is still being read.
    """
    from tempfile import TemporaryDirectory

    from qpx.dataset import Dataset

    dataset = dataset.resolve()
    prefix = _gene_map_dataset_prefix(dataset)
    out_dir = dataset if in_place else Path(output_folder).resolve()
    if not in_place:
        _copy_dataset_files(dataset, out_dir)

    with TemporaryDirectory(prefix=".gene_map_tmp-", dir=out_dir) as temporary:
        staging = Path(temporary)
        with Dataset(str(dataset), file_prefix=prefix) as qpx_dataset:
            written = _annotate_dataset_views(mapping, qpx_dataset, prefix, staging, out_dir)
            if not written:
                raise click.ClickException(f"No pg or feature Parquet file found in {dataset}")
            metadata = _stage_dataset_integrity(qpx_dataset, prefix, staging, out_dir)
            if metadata is not None:
                written.append(metadata)
        for source, target in written:
            source.replace(target)
            click.echo(f"  {target.name}: updated")

    _refresh_dataset_mudata(out_dir, prefix)

    click.echo(f"\nGene mapping complete. Output: {out_dir}")
