"""Transform CLI accessions implementation."""

import logging
from pathlib import Path
from typing import Optional

import click


def _accession_target_files(dataset: Path) -> list[str]:
    """Return the accession-bearing view files of a QPX dataset, relative to it."""
    from qpx.transforms.accession_normalizer import accession_view_files
    from qpx.transforms.utils import discover_qpx_file_prefix

    try:
        file_prefix = discover_qpx_file_prefix(dataset)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    if not any((dataset / f"{file_prefix}.{view}.parquet").is_file() for view in ("feature", "pg")):
        raise click.ClickException(f"No feature or protein-group Parquet file found in {dataset}")
    return [str(path.relative_to(dataset)) for path in accession_view_files(dataset, file_prefix)]


def _copy_unchanged_accession_files(dataset: Path, out_dir: Path, target_files: list[str]) -> None:
    """Copy files that accession normalization does not rewrite."""
    import shutil

    for f in dataset.iterdir():
        if f.name not in target_files and f.is_file():
            shutil.copy2(f, out_dir / f.name)


def _refresh_derived_files(out_dir: Path) -> None:
    """Rebuild the MuData view and the integrity record, when the dataset has them, from the rewritten views."""
    from qpx.cli.convert.postprocess import _maybe_compute_integrity
    from qpx.dataset import Dataset
    from qpx.transforms.utils import discover_qpx_file_prefix

    from .dataset import _refresh_dataset_mudata

    prefix = discover_qpx_file_prefix(out_dir)
    _refresh_dataset_mudata(out_dir, prefix)
    with Dataset(str(out_dir), file_prefix=prefix) as ds:
        meta = ds.dataset_meta.to_df() if ds.dataset_meta is not None else None
    record = meta.iloc[0].to_dict() if meta is not None and not meta.empty else {}
    if any(isinstance(record.get(field), dict) for field in ("file_checksums", "file_row_counts", "file_sizes_bytes")):
        _maybe_compute_integrity(out_dir, prefix, True)


@click.command("normalize-accessions")
@click.option(
    "--dataset",
    help="Path to a QPX dataset directory",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--direction",
    help="'forward' (sp|ACC|NAME → ACC) or 'reverse' (ACC → sp|ACC|NAME)",
    type=click.Choice(["forward", "reverse"]),
    default="forward",
)
@click.option(
    "--fasta",
    help="FASTA database file (required for --direction reverse)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
)
@click.option(
    "--in-place",
    help="Overwrite original files instead of writing to --output",
    is_flag=True,
)
@click.option(
    "--output",
    help="Output directory for the normalized views (required unless --in-place)",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def transform_normalize_accessions_cmd(
    dataset: Path,
    direction: str,
    fasta: Optional[Path],
    in_place: bool,
    output: Optional[Path],
    verbose: bool,
):
    """Normalize protein accessions in a QPX dataset.

    Forward (default): converts full UniProt identifiers to bare accessions.
        sp|P04114|APOB_HUMAN  →  P04114
        CONTAM_sp|CONTAM_P02768|...  →  CONTAM_P02768

    Reverse: converts bare accessions back to full UniProt format.
    Requires a FASTA database to look up the full identifiers.

    Normalizes protein accessions in the feature, pg, psm, pepmap, site and
    proteoform views. Every convert command already writes bare accessions;
    use this for datasets converted by an earlier qpx, or to restore full
    identifiers. Ids stay as they are, so links survive: a file whose identity
    composite holds an accession (pg_id) is then marked identity_source=provided.

    \b
    Examples:
        # Forward: strip sp|...|... to bare accessions
        qpxc transform normalize-accessions \\
            --dataset ./my_dataset --direction forward --in-place

    \b
        # Reverse: restore full UniProt identifiers from FASTA
        qpxc transform normalize-accessions \\
            --dataset ./my_dataset --direction reverse \\
            --fasta proteins.fasta --in-place

    \b
        # Forward to a new directory (non-destructive)
        qpxc transform normalize-accessions \\
            --dataset ./my_dataset --direction forward \\
            --output ./my_dataset_normalized
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if direction == "reverse" and fasta is None:
        raise click.UsageError("--fasta is required for --direction reverse")

    if not in_place and output is None:
        raise click.UsageError("Specify --in-place or --output")

    from qpx.transforms.accession_normalizer import (
        normalize_parquet,
        parse_fasta_headers,
    )

    # Build FASTA lookup for reverse
    fasta_lookup = None
    if direction == "reverse":
        fasta_lookup = parse_fasta_headers(fasta)
        click.echo(f"Loaded {len(fasta_lookup)} protein headers from FASTA")

    # Determine output directory
    out_dir = dataset if in_place else output
    if out_dir != dataset:
        out_dir.mkdir(parents=True, exist_ok=True)

    target_files = _accession_target_files(dataset)
    total_changed = 0
    files_changed = 0

    # Copy non-target files to output dir (once, outside the loop)
    if not in_place and out_dir != dataset:
        _copy_unchanged_accession_files(dataset, out_dir, target_files)

    for fname in target_files:
        src = dataset / fname
        dst = out_dir / fname
        dst.parent.mkdir(parents=True, exist_ok=True)
        result = normalize_parquet(
            parquet_path=src,
            output_path=dst,
            direction=direction,
            fasta_lookup=fasta_lookup,
        )
        total_changed += result["accessions_changed"]
        files_changed += bool(result["distinct_changed"])
        click.echo(f"  {fname}: {result['rows']:,} rows, {result['accessions_changed']:,} accessions normalized")

    if files_changed:
        _refresh_derived_files(out_dir)

    label = "stripped to bare" if direction == "forward" else "restored to full"
    click.echo(f"\nDone — {total_changed:,} accessions {label} in {out_dir}")
