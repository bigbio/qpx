"""Transform CLI quantify implementation."""

import logging
from pathlib import Path
from typing import Optional

import click

logger = logging.getLogger("qpx.cli.transform")


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


def _run_ibaq(peptide_df, options: dict) -> None:
    """Run iBAQ quantification via mokume."""
    fasta = options["fasta"]
    output = options["output"]

    import tempfile

    from mokume.quantification import peptides_to_protein

    click.echo("Running iBAQ quantification ...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = str(Path(tmpdir) / "peptides.parquet")
        peptide_df.to_parquet(tmp_path, index=False)
        peptides_to_protein(
            fasta=str(fasta),
            peptides=tmp_path,
            enzyme=options["enzyme"],
            normalize=options["normalize"],
            min_aa=options["min_aa"],
            max_aa=options["max_aa"],
            tpa=False,
            ruler=False,
            ploidy=options["ploidy"],
            cpc=options["cpc"],
            organism=options["organism"],
            output=str(output),
            verbose=options["verbose"],
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

    # Read as text: plexDIA labels such as 0/4/8 would otherwise become numbers.
    sdrf = pd.read_csv(sdrf_path, sep="\t", header=0, dtype=str)
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


def _struct_field(values, name: str):
    """Return one field of a struct array; null where the struct itself is null."""
    import pyarrow.compute as pc

    return pc.call_function("struct_field", [values], pc.StructFieldOptions(name))


def _flatten_intensities(batch):
    """Return one row per intensity entry with its feature's protein, sequence and run."""
    import pyarrow as pa
    import pyarrow.compute as pc

    intensities = batch.column("intensities")
    values = pc.call_function("list_flatten", [intensities])
    parents = pc.call_function("list_parent_indices", [intensities])
    features = pa.Table.from_batches([batch]).select(["anchor_protein", "sequence", "run_file_name"]).take(parents)
    # A null entry carries no label.
    label = pc.if_else(pc.is_valid(values), _struct_field(values, "label"), "")
    return pa.table(
        {
            "ProteinName": features["anchor_protein"],
            "PeptideCanonical": features["sequence"],
            "run_file_name": features["run_file_name"],
            "label": label,
            "NormIntensity": pc.cast(_struct_field(values, "intensity"), pa.float64()),
        }
    )


def _run_label_pairs(runs, labels):
    """Return the distinct ``(run, label)`` pairs in first-appearance order and each row's pair index."""
    import numpy as np
    import pyarrow.compute as pc

    options = pc.DictionaryEncodeOptions(null_encoding="encode")
    run_codes = pc.call_function("dictionary_encode", [runs], options).combine_chunks()
    label_codes = pc.call_function("dictionary_encode", [labels], options).combine_chunks()
    width = len(label_codes.dictionary)
    codes = run_codes.indices.to_numpy().astype(np.int64) * width + label_codes.indices.to_numpy()
    unique, first, inverse = np.unique(codes, return_index=True, return_inverse=True)
    order = np.argsort(first)
    position = np.empty(len(order), dtype=np.int64)
    position[order] = np.arange(len(order))
    pairs = unique[order]
    return (
        run_codes.dictionary.take(pairs // width).to_pylist(),
        label_codes.dictionary.take(pairs % width).to_pylist(),
        position[inverse],
    )


def _sdrf_samples(runs: list, labels: list, sdrf_path: Path):
    """Resolve each distinct ``(run, label)`` pair to its SDRF sample and condition."""
    import pandas as pd

    from qpx.converters.channel_labels import normalize_label
    from qpx.core.files import run_file_stem

    pairs = pd.DataFrame(
        {
            "_run_key": [run_file_stem(run) for run in runs],
            "_label_key": [normalize_label(label) for label in labels],
        }
    ).merge(
        _quantify_sdrf_mapping(sdrf_path),
        on=["_run_key", "_label_key"],
        how="left",
        validate="many_to_one",
    )
    unmapped = pairs["SampleID"].isna()
    if unmapped.any():
        missing = pairs.loc[unmapped, ["_run_key", "_label_key"]].drop_duplicates()
        rendered = ", ".join(f"{run_key}::{label_key}" for run_key, label_key in missing.itertuples(index=False, name=None))
        raise click.UsageError(f"QPX intensities have no matching SDRF run/label row: {rendered}")
    return pairs["SampleID"].array, pairs["Condition"].array


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

    Only the needed columns are read and the intensity lists are flattened in
    Arrow; samples are resolved once per distinct ``(run, label)`` pair.
    """
    import pandas as pd
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.dataset as ds
    import pyarrow.parquet as pq

    from qpx.version import check_feature_schema_compatible

    # A run name in a partition path is a string, whatever it looks like.
    partitioning = (
        ds.HivePartitioning.discover(schema=pa.schema([("run_file_name", pa.string())])) if feature_path.is_dir() else None
    )
    schema = pq.ParquetDataset(str(feature_path), partitioning=partitioning).schema
    check_feature_schema_compatible(schema, feature_path)
    has_decoy = "is_decoy" in schema.names
    columns = ["anchor_protein", "sequence", "run_file_name", "intensities"] + (["is_decoy"] if has_decoy else [])
    table = pq.read_table(str(feature_path), columns=columns, partitioning=partitioning)
    logger.info("Loaded %d QPX feature rows from %s", table.num_rows, feature_path)

    if has_decoy:
        table = table.filter(pc.invert(pc.fill_null(table["is_decoy"], False)))

    flattened = [_flatten_intensities(batch) for batch in table.to_batches()]
    entries = pa.concat_tables(flattened) if flattened else None
    if entries is None or entries.num_rows == 0:
        return pd.DataFrame(columns=["ProteinName", "PeptideCanonical", "NormIntensity", "SampleID", "Condition"])

    pair_runs, pair_labels, pair_index = _run_label_pairs(entries["run_file_name"], entries["label"])
    if sdrf_path is not None:
        sample_ids, conditions = _sdrf_samples(pair_runs, pair_labels, sdrf_path)
    else:
        sample_ids = pd.array(
            [run if not label else f"{run}::{label}" for run, label in zip(pair_runs, pair_labels)], dtype=object
        )
        conditions = None

    intensity = entries["NormIntensity"]
    protein = entries["ProteinName"]
    # A feature without an anchor protein (a peptide shared by several groups)
    # has no protein to be quantified under.
    keep = pc.and_kleene(
        pc.and_kleene(pc.is_valid(intensity), pc.greater(intensity, 0)),
        pc.not_equal(pc.fill_null(protein, ""), ""),
    )
    keep_rows = keep.to_numpy(zero_copy_only=False)
    result = (
        pa.table({"ProteinName": protein, "PeptideCanonical": entries["PeptideCanonical"], "NormIntensity": intensity})
        .filter(keep)
        .to_pandas()
    )
    rows = pair_index[keep_rows]
    result["SampleID"] = sample_ids.take(rows)
    result["Condition"] = conditions.take(rows) if conditions is not None else "not available"

    logger.info(
        "Prepared %d peptide rows: %d proteins, %d samples",
        len(result),
        result["ProteinName"].nunique(),
        result["SampleID"].nunique(),
    )
    return result


@click.command("quantify")
@click.option(
    "--feature-path",
    help="QPX feature.parquet file, or a directory partitioned by run_file_name",
    required=True,
    type=click.Path(exists=True, path_type=Path),
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
def transform_quantify_cmd(**kwargs):
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
    feature_path = kwargs["feature_path"]
    method = kwargs["method"]
    fasta = kwargs["fasta"]
    sdrf = kwargs["sdrf"]
    topn_n = kwargs["topn_n"]
    threads = kwargs["threads"]
    output = kwargs["output"]
    normalize = kwargs["normalize"]
    verbose = kwargs["verbose"]

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
        _run_ibaq(peptide_df, kwargs)
        return

    result_df = _run_generic_quantify(peptide_df, method_lower, method, topn_n, threads, normalize)

    # Step 3: Save output
    _save_result(result_df, output)
    click.echo(f"Quantification complete: {result_df['ProteinName'].nunique()} proteins -> {output}")
