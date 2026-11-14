"""Transform CLI protein properties implementation."""

import logging
from pathlib import Path
from typing import Optional

import click

from qpx.cli.resources import duckdb_resource_options

from .dataset import _copy_dataset_files, _gene_map_dataset_prefix, _stage_dataset_integrity


def _validate_destination(dataset: Path, in_place: bool, output_folder: Optional[Path]) -> None:
    if in_place == (output_folder is not None):
        raise click.UsageError("Specify exactly one of --in-place or --output-folder")
    if output_folder is not None:
        source, destination = dataset.resolve(), output_folder.resolve()
        if destination == source:
            raise click.UsageError("--output-folder is the dataset itself; use --in-place")
        if source in destination.parents:
            raise click.UsageError("--output-folder must not be inside the dataset directory")


def annotate_dataset_protein_properties(
    dataset: Path,
    fasta: Path,
    *,
    in_place: bool,
    output_folder: Optional[Path] = None,
    prefix: Optional[str] = None,
    fasta_uri: Optional[str] = None,
):
    """Fill null protein properties of a QPX dataset from a FASTA; return the report.

    Shared by ``qpxc transform protein-properties`` and the converters' optional
    ``--fasta``. Views are staged and moved into place, so an in-place run never
    truncates a parquet that is still being read. The MuData view carries none of
    these fields, so it is left as is.
    """
    from tempfile import TemporaryDirectory

    from qpx.dataset import Dataset
    from qpx.transforms.protein_properties import annotate_protein_properties

    dataset = dataset.resolve()
    prefix = _gene_map_dataset_prefix(dataset, command="protein-properties", prefix=prefix)
    out_dir = dataset if in_place else Path(output_folder).resolve()
    if not in_place:
        _copy_dataset_files(dataset, out_dir)

    with TemporaryDirectory(prefix=".protein_properties_tmp-", dir=out_dir) as temporary:
        staging = Path(temporary)
        report, names = annotate_protein_properties(dataset, prefix, fasta, staging, fasta_uri)
        if not any(name.endswith((".pg.parquet", ".feature.parquet")) for name in names):
            raise click.ClickException(f"No pg or feature Parquet file found in {dataset}")
        moves = [(staging / name, out_dir / name) for name in names]
        with Dataset(str(dataset), file_prefix=prefix) as qpx_dataset:
            metadata = _stage_dataset_integrity(qpx_dataset, prefix, staging, out_dir)
        if metadata is not None:
            moves.append(metadata)
        for source, target in moves:
            source.replace(target)
    _echo_protein_properties_report(report)
    return report


def _echo_protein_properties_report(report) -> None:
    click.echo(f"  FASTA: {report.fasta_entries} entries ({report.fasta_decoy_entries} decoys skipped)")
    evidence = ", ".join(report.coverage_evidence) or "none"
    click.echo(f"  coverage evidence from: {evidence}")
    if "psm" not in report.coverage_evidence:
        click.echo("  note: no PSM view, so coverage counts only peptides assigned to each protein group")
    click.echo(
        f"  pg: {report.pg_coverage_filled} sequence_coverage and {report.pg_molecular_weight_filled} molecular_weight "
        f"filled of {report.pg_rows_eligible} target rows needing them; {report.pg_anchors_not_in_fasta} anchors not in the FASTA"
    )
    click.echo(
        f"  feature: {report.feature_positions_filled} pg_positions filled of {report.feature_rows_eligible} target rows "
        f"needing them; {report.feature_rows_without_match} with no occurrence in their group's FASTA sequences"
    )
    if report.unmatched_examples:
        click.echo(f"  not in the FASTA, e.g.: {', '.join(report.unmatched_examples)}")
    match_rate = report.pg_match_rate
    if match_rate is not None and match_rate < 0.5:
        click.echo(
            f"  WARNING: only {100 * match_rate:.1f}% of target protein groups were found in the FASTA; "
            "check that it is the database used for the search"
        )


@click.command("protein-properties")
@click.option(
    "--dataset",
    required=True,
    help="Path to a QPX dataset directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--fasta",
    required=True,
    help="Protein FASTA used for the search (optionally .gz)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option("--fasta-uri", default=None, help="Where others can get this FASTA, recorded in provenance beside its checksum")
@click.option("--in-place", is_flag=True, help="Overwrite the dataset's files instead of writing to --output-folder")
@click.option(
    "--output-folder",
    default=None,
    help="Write an annotated copy here (must be empty and outside the dataset)",
    type=click.Path(file_okay=False, path_type=Path),
)
@duckdb_resource_options
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def transform_protein_properties_cmd(
    dataset: Path, fasta: Path, fasta_uri: Optional[str], in_place: bool, output_folder: Optional[Path], verbose: bool
):
    """Fill protein properties the producer did not record, from a FASTA.

    For target rows only, and only where the value is null, fills
    ``pg.sequence_coverage`` (the anchor covered by the dataset's target peptides
    whose evidence names it: PSM protein accessions, feature groups, recorded positions),
    ``pg.molecular_weight`` (anchor average mass, kDa) and ``feature.pg_positions``
    (every one-based occurrence in each group member), and writes the ``protein``
    view of the accessions the FASTA has. Values a producer recorded
    are never overwritten. Proteins absent from the FASTA, such as DIA-NN's internal
    decoys, stay null and are counted in the report. The step is recorded in the
    provenance view with the FASTA's checksum.

    \b
    Example:
        qpxc transform protein-properties \\
            --dataset ./qpx_output \\
            --fasta search_database.fasta \\
            --in-place
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    _validate_destination(dataset, in_place, output_folder)
    out_dir = annotate_dataset_protein_properties(
        dataset, fasta, in_place=in_place, output_folder=output_folder, fasta_uri=fasta_uri
    )
    click.echo(f"\nProtein properties complete. Output: {dataset if in_place else output_folder}")
    return out_dir
