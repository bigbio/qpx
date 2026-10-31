"""Transform CLI accessions implementation."""

import logging
from pathlib import Path
from typing import Optional

import click


def _accession_target_files(dataset: Path) -> list[str]:
    """Return the accession-bearing view files for a QPX dataset."""
    from qpx.transforms.utils import discover_qpx_file_prefix

    try:
        file_prefix = discover_qpx_file_prefix(dataset)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    target_files = [
        name for name in (f"{file_prefix}.feature.parquet", f"{file_prefix}.pg.parquet") if (dataset / name).is_file()
    ]
    if not target_files:
        raise click.ClickException(f"No feature or protein-group Parquet file found in {dataset}")
    target_files.extend(
        name
        for name in (f"{file_prefix}.{view}.parquet" for view in ("site", "proteoform", "psm", "pepmap"))
        if (dataset / name).is_file()
    )
    return target_files


def _copy_unchanged_accession_files(dataset: Path, out_dir: Path, target_files: list[str]) -> None:
    """Copy files that accession normalization does not rewrite."""
    import shutil

    for f in dataset.iterdir():
        if f.name not in target_files and f.is_file():
            shutil.copy2(f, out_dir / f.name)


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

    Normalizes protein accessions in feature.parquet, pg.parquet, site.parquet,
    proteoform.parquet, psm.parquet and pepmap.parquet files. Ids stay as they
    are, so links survive: a file whose identity composite holds an accession
    (pg_id) is then marked identity_source=provided.

    \b
    Examples:
        # Forward: strip sp|...|... to bare accessions
        qpxc transform normalize-accessions \\
            --dataset ./my_dataset --direction forward --in-place

        # Reverse: restore full UniProt identifiers from FASTA
        qpxc transform normalize-accessions \\
            --dataset ./my_dataset --direction reverse \\
            --fasta proteins.fasta --in-place

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

    # Copy non-target files to output dir (once, outside the loop)
    if not in_place and out_dir != dataset:
        _copy_unchanged_accession_files(dataset, out_dir, target_files)

    for fname in target_files:
        src = dataset / fname
        dst = out_dir / fname
        result = normalize_parquet(
            parquet_path=src,
            output_path=dst,
            direction=direction,
            fasta_lookup=fasta_lookup,
        )
        total_changed += result["accessions_changed"]
        click.echo(f"  {fname}: {result['rows']:,} rows, {result['accessions_changed']:,} accessions normalized")

    label = "stripped to bare" if direction == "forward" else "restored to full"
    click.echo(f"\nDone — {total_changed:,} accessions {label} in {out_dir}")
