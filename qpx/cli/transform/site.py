"""Transform CLI site quantification implementation."""

import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional

import click

from qpx.cli.resources import duckdb_resource_options

from .dataset import _copy_dataset_files, _gene_map_dataset_prefix, _stage_dataset_integrity
from .protein_properties import _validate_destination


def _provenance_step(modifications, steps_path: Path, destination: Path) -> None:
    import pyarrow.parquet as pq

    from qpx._version import __version__
    from qpx.writers.provenance import ProvenanceWriter

    steps = pq.read_table(steps_path).to_pylist() if steps_path.is_file() else []
    step = {
        "step_order": max((s.get("step_order") or 0 for s in steps), default=0) + 1,
        "step_category": "annotation",
        "step_name": "site_quantification_from_features",
        "tool_name": "qpx",
        "tool_version": __version__,
        "parameters": [{"key": "modifications", "value": ",".join(modifications)}],
        "output_views": ["site"],
    }
    with ProvenanceWriter(destination, creator="qpx") as writer:
        writer.write_batch([*steps, step])


@click.command("site")
@click.option(
    "--dataset",
    required=True,
    help="Path to a QPX dataset directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--modification", "modifications", required=True, multiple=True, help="Unimod name or accession; repeat for several"
)
@click.option("--in-place", is_flag=True, help="Write into the dataset instead of --output-folder")
@click.option(
    "--output-folder",
    default=None,
    help="Write a copy with the site view here (must be empty and outside the dataset)",
    type=click.Path(file_okay=False, path_type=Path),
)
@duckdb_resource_options
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def transform_site_cmd(dataset: Path, modifications: tuple, in_place: bool, output_folder: Optional[Path], verbose: bool):
    """Quantify modification sites from the feature view into site.parquet.

    One row per site, run and label: the summed intensity of the peptides that place
    the modification there, and the apparent stoichiometry against all peptides covering
    the residue. Needs feature.pg_positions, which 'qpxc transform protein-properties
    --fasta' fills. The anchor protein is the one a shared peptide's site is placed on.

    \b
    Example:
        qpxc transform site --dataset ./qpx_output --modification Phospho --in-place
    """
    from qpx.dataset import Dataset
    from qpx.transforms.site_quant import derive_sites

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    _validate_destination(dataset, in_place, output_folder)
    dataset = dataset.resolve()
    prefix = _gene_map_dataset_prefix(dataset, command="site", prefix=None)
    feature = dataset / f"{prefix}.feature.parquet"
    if not feature.is_file():
        raise click.ClickException(f"No feature Parquet file found in {dataset}")
    out_dir = dataset if in_place else Path(output_folder).resolve()
    if not in_place:
        _copy_dataset_files(dataset, out_dir)
    with TemporaryDirectory(prefix=".site_tmp-", dir=out_dir) as temporary:
        staging = Path(temporary)
        try:
            rows = derive_sites(feature, list(modifications), staging / f"{prefix}.site.parquet")
        except ValueError as error:
            raise click.ClickException(str(error)) from error
        provenance = f"{prefix}.provenance.parquet"
        _provenance_step(modifications, dataset / provenance, staging / provenance)
        moves = [(staging / name, out_dir / name) for name in (f"{prefix}.site.parquet", provenance)]
        with Dataset(str(dataset), file_prefix=prefix) as qpx_dataset:
            metadata = _stage_dataset_integrity(qpx_dataset, prefix, staging, out_dir)
        if metadata is not None:
            moves.append(metadata)
        for source, target in moves:
            source.replace(target)
    click.echo(f"{rows} site rows written to {out_dir}")
