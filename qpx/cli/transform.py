"""
Transform subcommands — compute derived data from QPX datasets.

Subcommands:
    qpxc transform gene-map   — Map genes from FASTA
"""

import logging
from pathlib import Path
from typing import Optional

import click

logger = logging.getLogger("qpx.cli.transform")


# ---------------------------------------------------------------------------
# Transform group
# ---------------------------------------------------------------------------


@click.group()
def transform():
    """Transform QPX data into derived representations."""
    pass


# ---------------------------------------------------------------------------
# Gene Mapping
# ---------------------------------------------------------------------------


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


@transform.command("gene-map")
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


def _copy_dataset_files(dataset: Path, out_dir: Path) -> None:
    """Copy the complete dataset, including directory-backed views."""
    import shutil

    if out_dir.exists() and any(out_dir.iterdir()):
        raise click.ClickException("--output-folder must be empty")
    shutil.copytree(dataset, out_dir, dirs_exist_ok=True)


def _gene_map_dataset_prefix(dataset: Path, command: str = "gene-map", prefix: Optional[str] = None) -> str:
    """Reject unsupported quantification layouts before copying or annotation.

    ``prefix`` skips discovery when the caller already knows it (a converter that
    just wrote the files), so a second dataset in the same folder cannot make the
    prefix ambiguous.
    """
    from qpx.transforms.utils import discover_qpx_file_prefix

    for view in ("pg", "feature"):
        if any((dataset / view).rglob("*.parquet")):
            raise click.ClickException(f"Partitioned {view} data is not supported by {command}; use flat Parquet views")
    if prefix:
        return prefix
    try:
        return discover_qpx_file_prefix(dataset)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc


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


def _stage_dataset_integrity(qpx_dataset, prefix: str, staging: Path, out_dir: Path) -> tuple[Path, Path] | None:
    """Update existing integrity records for the staged quantification files."""
    from qpx.dataset import Dataset
    from qpx.transforms.utils import count_staged_structures
    from qpx.writers import DatasetWriter

    if qpx_dataset.dataset_meta is None:
        return None
    metadata = qpx_dataset.dataset_meta.to_df()
    if metadata.empty:
        return None
    record = metadata.iloc[0].to_dict()
    fields = ("file_checksums", "file_row_counts", "file_sizes_bytes")
    if not any(isinstance(record.get(field), dict) for field in fields):
        return None

    with Dataset(str(staging), file_prefix=prefix) as staged_dataset:
        integrity = staged_dataset.compute_integrity()
    for field in fields:
        if isinstance(record.get(field), dict):
            record[field].update(integrity[field])
    record["total_structures"] = count_staged_structures(out_dir, staging, integrity["file_checksums"])
    record["packaged_at"] = integrity["packaged_at"]
    name = f"{prefix}.dataset.parquet"
    with DatasetWriter(staging / name) as writer:
        writer.write_batch([record])
    return staging / name, out_dir / name


def _refresh_dataset_mudata(out_dir: Path, prefix: str) -> None:
    """Rebuild the dataset's MuData view, so a stale h5mu never outlives the parquet."""
    from qpx.mudata import write_dataset_mudata

    if not (out_dir / f"{prefix}.h5mu").is_file():
        return

    refreshed = write_dataset_mudata(out_dir, prefix)
    if refreshed is None:
        click.echo("  h5mu: could not be rebuilt (see log); parquet views are annotated")
    else:
        click.echo(f"  {refreshed.name}: MuData view rebuilt")


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


# ---------------------------------------------------------------------------
# Protein properties from FASTA
# ---------------------------------------------------------------------------


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
        report, names = annotate_protein_properties(dataset, prefix, fasta, staging)
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


@transform.command("protein-properties")
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
@click.option("--in-place", is_flag=True, help="Overwrite the dataset's files instead of writing to --output-folder")
@click.option(
    "--output-folder",
    default=None,
    help="Write an annotated copy here (must be empty and outside the dataset)",
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def transform_protein_properties_cmd(dataset: Path, fasta: Path, in_place: bool, output_folder: Optional[Path], verbose: bool):
    """Fill protein properties the producer did not record, from a FASTA.

    For target rows only, and only where the value is null, fills
    ``pg.sequence_coverage`` (the anchor covered by the dataset's target peptides
    whose evidence names it: PSM protein accessions, feature groups, recorded positions),
    ``pg.molecular_weight`` (anchor average mass, kDa) and ``feature.pg_positions``
    (every one-based occurrence in each group member). Values a producer recorded
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
    out_dir = annotate_dataset_protein_properties(dataset, fasta, in_place=in_place, output_folder=output_folder)
    click.echo(f"\nProtein properties complete. Output: {dataset if in_place else output_folder}")
    return out_dir


# ---------------------------------------------------------------------------
# Accession Normalization
# ---------------------------------------------------------------------------


def _accession_target_files(dataset: Path) -> list[str]:
    """Return the Feature and protein-group files for a QPX dataset."""
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
    return target_files


@transform.command("normalize-accessions")
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
    help="Output directory (default: overwrites in place)",
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

    Normalizes anchor_protein and pg_accessions in both
    feature.parquet and pg.parquet files.

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
        import shutil

        for f in dataset.iterdir():
            if f.name not in target_files and f.is_file():
                shutil.copy2(f, out_dir / f.name)

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


# ---------------------------------------------------------------------------
# SDRF Metadata Update
# ---------------------------------------------------------------------------


@transform.command("update-metadata")
@click.option(
    "--dataset",
    help="Path to a QPX dataset directory",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--sdrf",
    help="Path to the updated SDRF TSV file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--old-sdrf",
    help="Path to the original SDRF (for safety checks). If omitted, protected-field validation is skipped.",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
)
@click.option(
    "--force",
    help="Apply changes even if protected fields (instrument, enzymes, "
    "modifications, fractions, labels) have changed. Use with caution.",
    is_flag=True,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def transform_update_metadata_cmd(
    dataset: Path,
    sdrf: Path,
    old_sdrf: Optional[Path],
    force: bool,
    verbose: bool,
):
    """Update sample and run metadata from a revised SDRF file.

    Re-generates sample.parquet and run.parquet from the new SDRF.
    Original files are backed up as *.parquet.bak before overwriting.

    Safe to update (metadata-only, no impact on data):
        disease, organism_part, cell_type, cell_line, sex, age,
        treatment, individual, ancestry, developmental_stage,
        and any additional characteristics[*] columns.

    Protected fields (will BLOCK unless --force is used):
        instrument, enzymes, modification_parameters, fraction,
        dissociation_method, label/channel mapping, data file references.

    If --old-sdrf is provided, the tool compares old vs new SDRF and
    blocks if any protected fields changed. Without --old-sdrf, the
    safety check is skipped (useful for first-time metadata enrichment).

    \b
    Examples:
        # Update with safety check
        qpxc transform update-metadata \\
            --dataset ./my_project \\
            --sdrf ./updated_sdrf.tsv \\
            --old-sdrf ./original_sdrf.tsv

        # Update without safety check (first-time enrichment)
        qpxc transform update-metadata \\
            --dataset ./my_project \\
            --sdrf ./enriched_sdrf.tsv

        # Force update even if protected fields changed
        qpxc transform update-metadata \\
            --dataset ./my_project \\
            --sdrf ./new_sdrf.tsv \\
            --old-sdrf ./old_sdrf.tsv --force
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    from qpx.transforms.sdrf_updater import update_metadata

    try:
        result = update_metadata(
            dataset_path=dataset,
            new_sdrf_path=sdrf,
            old_sdrf_path=old_sdrf,
            force=force,
        )
    except ValueError as e:
        raise click.ClickException(str(e)) from e

    for w in result["warnings"]:
        if w.startswith("BLOCKED"):
            click.echo(f"  ⚠ {w}", err=True)
        elif w.startswith("WARNING"):
            click.echo(f"  ⚠ {w}", err=True)

    click.echo(f"\nDone — updated {result['samples_updated']} samples, {result['runs_updated']} runs in {dataset}")


# ---------------------------------------------------------------------------
# Protein Quantification via mokume
# ---------------------------------------------------------------------------

QUANT_METHODS = ["directlfq", "maxlfq", "topn", "top3", "ibaq", "sum"]
_SUPPORTED_SUFFIXES = {".parquet", ".tsv", ".csv"}


def _validate_quantify_inputs(method_lower: str, fasta: Optional[Path]) -> None:
    """Validate CLI inputs before running quantification."""
    import importlib

    if method_lower == "ibaq" and not fasta:
        raise click.UsageError("The --fasta option is required for the ibaq method")
    try:
        importlib.import_module("mokume.quantification")
    except ImportError as exc:
        raise click.UsageError('mokume is not installed. Install with: pip install "qpx[quantify]"') from exc
    if method_lower == "directlfq":
        try:
            from mokume.quantification import is_directlfq_available
        except ImportError as exc:
            raise click.UsageError('DirectLFQ is not installed. Install with: pip install "qpx[quantify]"') from exc
        if not is_directlfq_available():
            raise click.UsageError('DirectLFQ is not installed. Install with: pip install "qpx[quantify]"')


def _run_ibaq(peptide_df, fasta, enzyme, normalize, min_aa, max_aa, ploidy, cpc, organism, output, verbose) -> None:
    """Run iBAQ quantification via mokume."""
    import tempfile

    from mokume.quantification import peptides_to_protein

    click.echo("Running iBAQ quantification ...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = str(Path(tmpdir) / "peptides.parquet")
        peptide_df.to_parquet(tmp_path, index=False)
        peptides_to_protein(
            fasta=str(fasta),
            peptides=tmp_path,
            enzyme=enzyme,
            normalize=normalize,
            min_aa=min_aa,
            max_aa=max_aa,
            tpa=False,
            ruler=False,
            ploidy=ploidy,
            cpc=cpc,
            organism=organism,
            output=str(output),
            verbose=verbose,
            qc_report=str(output.with_suffix(".qc.pdf")),
        )
    click.echo(f"iBAQ results saved to {output}")


def _run_generic_quantify(peptide_df, method_lower, method, topn_n, threads, normalize):
    """Run non-iBAQ quantification and return result DataFrame."""
    from mokume.quantification import get_quantification_method

    click.echo(f"Using {method} quantification method")

    method_kwargs = {"topn": {"n": topn_n}, "top3": {"n": 3}, "maxlfq": {"threads": threads}}
    resolved = "topn" if method_lower == "top3" else method_lower
    quant_method = get_quantification_method(resolved, **method_kwargs.get(method_lower, {}))

    click.echo(
        f"Quantifying {peptide_df['ProteinName'].nunique()} proteins across {peptide_df['SampleID'].nunique()} samples ..."
    )
    result_df = quant_method.quantify(
        peptide_df,
        protein_column="ProteinName",
        peptide_column="PeptideCanonical",
        intensity_column="NormIntensity",
        sample_column="SampleID",
    )

    if normalize:
        intensity_cols = [c for c in result_df.columns if "Intensity" in c]
        if intensity_cols:
            int_col = intensity_cols[-1]
            result_df[f"{int_col}Norm"] = result_df.groupby("SampleID")[int_col].transform(lambda x: x / x.sum())
    return result_df


def _save_result(result_df, output: Path) -> None:
    """Save quantification results to the specified output format."""
    suffix = output.suffix.lower()
    if suffix not in _SUPPORTED_SUFFIXES:
        raise click.UsageError(f"Unsupported output format '{suffix}'. Use .parquet, .tsv, or .csv")
    output.parent.mkdir(parents=True, exist_ok=True)
    if suffix == ".parquet":
        result_df.to_parquet(str(output), index=False)
    elif suffix == ".tsv":
        result_df.to_csv(str(output), sep="\t", index=False)
    else:
        result_df.to_csv(str(output), index=False)


def _quantify_sdrf_mapping(sdrf_path: Path):
    """Build an unambiguous QPX run/channel to sample/condition mapping."""
    import pandas as pd

    from qpx.converters.channel_labels import normalize_label
    from qpx.core.files import run_file_stem
    from qpx.core.sdrf import validate_sdrf_data_files

    sdrf = pd.read_csv(sdrf_path, sep="\t", header=0)
    sdrf.columns = sdrf.columns.str.lower()
    validate_sdrf_data_files(sdrf)

    source_col = "source name"
    label_col = "comment[label]"
    missing = [column for column in (source_col, label_col) if column not in sdrf.columns]
    if missing:
        raise click.UsageError(f"SDRF is missing required column(s): {', '.join(missing)}")

    sample_ids = sdrf[source_col].astype("string").str.strip()
    if sample_ids.isna().any() or sample_ids.eq("").any():
        raise click.UsageError("SDRF column 'source name' contains missing or blank values")

    factor_columns = [column for column in sdrf.columns if column.startswith("factor value[")]
    conditions = sdrf[factor_columns[0]].astype("string").str.strip() if factor_columns else sample_ids
    conditions = conditions.mask(conditions.isna() | conditions.eq(""), sample_ids)

    mapping = pd.DataFrame(
        {
            "_run_key": sdrf["comment[data file]"].map(run_file_stem),
            "_label_key": sdrf[label_col].map(normalize_label),
            "SampleID": sample_ids,
            "Condition": conditions,
        }
    )
    duplicate = mapping.duplicated(["_run_key", "_label_key"], keep=False)
    if duplicate.any():
        pairs = mapping.loc[duplicate, ["_run_key", "_label_key"]].drop_duplicates()
        rendered = ", ".join(f"{run_key}::{label_key}" for run_key, label_key in pairs.itertuples(index=False, name=None))
        raise click.UsageError(f"SDRF contains duplicate run/label mappings: {rendered}")
    return mapping


def _qpx_feature_to_peptide_df(feature_path: Path, sdrf_path: Optional[Path] = None):
    """Read QPX feature.parquet and flatten to mokume peptide DataFrame.

    Mokume expects columns: ``ProteinName``, ``PeptideCanonical``,
    ``NormIntensity``, ``SampleID`` and, for normalized iBAQ, ``Condition``.

    QPX stores intensities as a list of ``{label, intensity}`` structs.
    For label-free data each feature typically has a single intensity entry;
    for TMT/iTRAQ data the list contains one entry per channel.  This
    function explodes **all** intensity entries into separate rows so that
    multiplexed channels are preserved.  When an SDRF is provided,
    ``(run_file_name, label)`` is mapped to its real sample and condition.
    Without one, the existing ``run::label`` sample key is retained and
    condition is explicitly marked ``not available``.
    """
    import pandas as pd

    df = pd.read_parquet(str(feature_path))
    logger.info("Loaded %d QPX feature rows from %s", len(df), feature_path)

    # Filter decoys early
    if "is_decoy" in df.columns:
        df = df[~df["is_decoy"].astype(bool)]

    # Map protein column
    protein_col = "anchor_protein" if "anchor_protein" in df.columns else "pg_accessions"
    df["ProteinName"] = df[protein_col].apply(
        lambda x: x if isinstance(x, str) else (";".join(x) if isinstance(x, list) else str(x))
    )
    df["PeptideCanonical"] = df.get("sequence", df.get("peptidoform", ""))

    # Drop rows without intensities, then explode vectorised
    df = df[df["intensities"].apply(lambda x: x is not None and len(x) > 0)]
    if df.empty:
        return pd.DataFrame(columns=["ProteinName", "PeptideCanonical", "NormIntensity", "SampleID", "Condition"])

    exploded = df[["ProteinName", "PeptideCanonical", "run_file_name", "intensities"]].explode("intensities")
    # Extract struct fields
    exploded["NormIntensity"] = exploded["intensities"].apply(lambda e: e.get("intensity", 0.0) if isinstance(e, dict) else 0.0)
    exploded["_label"] = exploded["intensities"].apply(lambda e: e.get("label", "") if isinstance(e, dict) else "")
    if sdrf_path is not None:
        from qpx.converters.channel_labels import normalize_label
        from qpx.core.files import run_file_stem

        exploded["_run_key"] = exploded["run_file_name"].map(run_file_stem)
        exploded["_label_key"] = exploded["_label"].map(normalize_label)
        exploded = exploded.merge(
            _quantify_sdrf_mapping(sdrf_path),
            on=["_run_key", "_label_key"],
            how="left",
            validate="many_to_one",
        )
        unmapped = exploded["SampleID"].isna()
        if unmapped.any():
            pairs = exploded.loc[unmapped, ["_run_key", "_label_key"]].drop_duplicates()
            rendered = ", ".join(f"{run_key}::{label_key}" for run_key, label_key in pairs.itertuples(index=False, name=None))
            raise click.UsageError(f"QPX intensities have no matching SDRF run/label row: {rendered}")
    else:
        has_label = exploded["_label"].astype(bool)
        exploded["SampleID"] = exploded["run_file_name"]
        exploded.loc[has_label, "SampleID"] = exploded.loc[has_label, "run_file_name"] + "::" + exploded.loc[has_label, "_label"]
        exploded["Condition"] = "not available"
    # Filter invalid intensities
    result = exploded.loc[
        exploded["NormIntensity"].notna() & (exploded["NormIntensity"] > 0),
        ["ProteinName", "PeptideCanonical", "NormIntensity", "SampleID", "Condition"],
    ].copy()
    result = result[result["ProteinName"] != ""]

    logger.info(
        "Prepared %d peptide rows: %d proteins, %d samples",
        len(result),
        result["ProteinName"].nunique(),
        result["SampleID"].nunique(),
    )
    return result


@transform.command("quantify")
@click.option(
    "--feature-path",
    help="QPX feature.parquet file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--method",
    help="Quantification method (directlfq, maxlfq, topn, top3, ibaq, sum)",
    type=click.Choice(QUANT_METHODS, case_sensitive=False),
    default="directlfq",
)
@click.option(
    "--fasta",
    help="FASTA database (required for ibaq method)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
)
@click.option(
    "--sdrf",
    help="SDRF used to map QPX run/channel intensities to samples and conditions",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
)
@click.option(
    "--enzyme",
    help="Enzyme for iBAQ digestion (default: Trypsin)",
    default="Trypsin",
)
@click.option(
    "--topn-n",
    help="N for TopN method (default: 3)",
    default=3,
    type=int,
)
@click.option(
    "--threads",
    help="Parallel threads for MaxLFQ (-1 = all cores)",
    default=-1,
    type=int,
)
@click.option(
    "--output",
    "-o",
    help="Output file path (.parquet, .tsv, or .csv)",
    required=True,
    type=click.Path(path_type=Path),
)
@click.option("--normalize", help="Normalize quantification values", is_flag=True)
@click.option("--organism", help="Organism for iBAQ (default: human)", default="human")
@click.option("--ploidy", help="Ploidy for iBAQ ruler (default: 2)", default=2, type=int)
@click.option("--cpc", help="Cell copies per cell for iBAQ ruler (default: 200)", default=200, type=float)
@click.option("--min-aa", help="Min peptide length for iBAQ (default: 7)", default=7, type=int)
@click.option("--max-aa", help="Max peptide length for iBAQ (default: 30)", default=30, type=int)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def transform_quantify_cmd(
    feature_path: Path,
    method: str,
    fasta: Optional[Path],
    sdrf: Optional[Path],
    enzyme: str,
    topn_n: int,
    threads: int,
    output: Path,
    normalize: bool,
    organism: str,
    ploidy: int,
    cpc: float,
    min_aa: int,
    max_aa: int,
    verbose: bool,
):
    """Compute protein quantification from QPX feature data using mokume.

    Reads a QPX feature.parquet file, extracts peptide-level intensities,
    and computes protein-level quantification using the selected method.

    \b
    Supported methods:
      directlfq  — DirectLFQ intensity traces (default)
      maxlfq     — MaxLFQ delayed normalization
      topn       — Average of N most intense peptides
      top3       — Average of 3 most intense peptides
      ibaq       — Intensity-Based Absolute Quantification (requires --fasta)
      sum        — Sum of all peptide intensities

    \b
    Examples:
        # DirectLFQ from QPX feature
        qpxc transform quantify \\
            --feature-path ./qpx_output/feature.parquet \\
            --method directlfq \\
            -o proteins_directlfq.parquet

        # iBAQ (requires FASTA)
        qpxc transform quantify \\
            --feature-path ./qpx_output/feature.parquet \\
            --method ibaq --fasta proteome.fasta --sdrf experiment.sdrf.tsv \\
            -o proteins_ibaq.tsv

        # MaxLFQ with 8 threads
        qpxc transform quantify \\
            --feature-path ./qpx_output/feature.parquet \\
            --method maxlfq --threads 8 \\
            -o proteins_maxlfq.parquet
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    method_lower = method.lower()
    _validate_quantify_inputs(method_lower, fasta)

    # Step 1: Read QPX feature and flatten to mokume format
    click.echo(f"Reading QPX feature data from {feature_path}")
    peptide_df = _qpx_feature_to_peptide_df(feature_path, sdrf)
    if peptide_df.empty:
        raise click.UsageError("No valid peptide intensities found in the feature file")

    output.parent.mkdir(parents=True, exist_ok=True)

    # Step 2: Run quantification
    if method_lower == "ibaq":
        _run_ibaq(peptide_df, fasta, enzyme, normalize, min_aa, max_aa, ploidy, cpc, organism, output, verbose)
        return

    result_df = _run_generic_quantify(peptide_df, method_lower, method, topn_n, threads, normalize)

    # Step 3: Save output
    _save_result(result_df, output)
    click.echo(f"Quantification complete: {result_df['ProteinName'].nunique()} proteins -> {output}")
