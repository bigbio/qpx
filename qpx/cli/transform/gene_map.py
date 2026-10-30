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
            --parquet-path ./output/psm.parquet \\
            --fasta proteins.fasta \\
            --output-folder ./output

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
        return

    import pandas as pd

    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(str(parquet_path))
    protein_col = "pg_accessions" if "pg_accessions" in df.columns else "protein_accessions"
    annotated = mapping.annotate_dataframe(df, protein_col=protein_col)
    output_path = output_folder / parquet_path.name
    annotated.to_parquet(str(output_path), index=False)

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
