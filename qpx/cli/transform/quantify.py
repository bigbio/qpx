"""Transform CLI quantify implementation."""

import json
import logging
from pathlib import Path
from typing import Optional

import click

logger = logging.getLogger("qpx.cli.transform")


QUANT_METHODS = ["directlfq", "maxlfq", "topn", "top3", "ibaq", "sum"]


_SUPPORTED_SUFFIXES = {".parquet", ".tsv", ".csv", ".h5ad"}


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
    """Run iBAQ quantification via mokume and save it in the format ``--output`` names."""
    import tempfile

    import pandas as pd
    from mokume.quantification import peptides_to_protein

    output = options["output"]
    click.echo("Running iBAQ quantification ...")
    with tempfile.TemporaryDirectory() as tmpdir:
        peptides = Path(tmpdir) / "peptides.parquet"
        table = Path(tmpdir) / "ibaq.tsv"
        peptide_df.to_parquet(peptides, index=False)
        # mokume writes tab-separated text, whatever the path is called.
        peptides_to_protein(
            fasta=str(options["fasta"]),
            peptides=str(peptides),
            enzyme=options["enzyme"],
            normalize=options["normalize"],
            min_aa=options["min_aa"],
            max_aa=options["max_aa"],
            tpa=False,
            ruler=False,
            ploidy=options["ploidy"],
            cpc=options["cpc"],
            organism=options["organism"],
            output=str(table),
            verbose=options["verbose"],
            qc_report=str(output.with_suffix(".qc.pdf")),
        )
        result = pd.read_csv(table, sep="\t")
    dropped = peptide_df["ProteinName"].nunique() - result["ProteinName"].nunique()
    if dropped:
        click.echo(f"WARNING: {dropped} protein(s) were not in the FASTA or had no iBAQ value and were dropped", err=True)
    recorded = ("enzyme", "normalize", "min_aa", "max_aa", "ploidy", "cpc", "organism")
    _save_result(result, output, {"method": "ibaq", **{name: options[name] for name in recorded}})
    click.echo(f"iBAQ results saved to {output}")


def _run_generic_quantify(peptide_df, method_lower, method, topn_n, threads, normalize):
    """Run non-iBAQ quantification and return result DataFrame."""
    from mokume.quantification import MaxLFQQuantification, get_quantification_method

    click.echo(f"Using {method} quantification method")

    if method_lower == "maxlfq":
        # Without force_builtin, mokume runs DirectLFQ for this method whenever it is installed.
        quant_method = MaxLFQQuantification(threads=threads, force_builtin=True)
    elif method_lower == "directlfq":
        quant_method = get_quantification_method("directlfq", num_cores=None if threads == -1 else threads)
    else:
        resolved = "topn" if method_lower == "top3" else method_lower
        quant_method = get_quantification_method(resolved, n=3 if method_lower == "top3" else topn_n)

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


def _save_ae_anndata(result_df, output: Path, parameters: dict) -> None:
    """Write iBAQ results as the absolute-expression AnnData view: ``X`` is log2 iBAQ, with raw values in a layer."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    from qpx._version import __version__
    from qpx.version import QPX_SPEC_VERSION

    if "Ibaq" not in result_df.columns:
        raise click.UsageError(".h5ad output is the absolute-expression view, so it needs --method ibaq")
    pivot = {
        column: result_df.pivot_table(index="SampleID", columns="ProteinName", values=column, aggfunc="first")
        for column in ("Ibaq", "IbaqPpb", "CopyNumber", "Concentration[nM]")
        if column in result_df.columns
    }
    raw = pivot["Ibaq"]
    layer_names = {
        "Ibaq": "ibaq_raw",
        "IbaqPpb": "ibaq_ppb",
        "CopyNumber": "copies_per_cell",
        "Concentration[nM]": "concentration_nm",
    }
    adata = ad.AnnData(
        X=np.log2(raw.where(raw > 0).to_numpy(dtype="float64")),
        obs=pd.DataFrame(index=raw.index.astype(str)),
        var=pd.DataFrame(index=raw.columns.astype(str)),
        layers={layer_names[column]: frame.reindex_like(raw).to_numpy(dtype="float64") for column, frame in pivot.items()},
    )
    adata.uns["qpx_version"] = QPX_SPEC_VERSION
    adata.uns["writer_version"] = __version__
    adata.uns["file_type"] = "absolute_expression"
    adata.uns["log_base"] = 2
    # quantify sums a peptide over the runs of a sample and treats all groups led by one protein as that protein.
    adata.uns["var_key"] = "anchor_protein"
    adata.uns["aggregation"] = "mokume_peptides_to_protein_per_sample"
    adata.uns["quantification_parameters"] = json.dumps(parameters, sort_keys=True)
    adata.write_h5ad(output)


def _save_result(result_df, output: Path, parameters: dict | None = None) -> None:
    """Save quantification results to the specified output format."""
    suffix = output.suffix.lower()
    if suffix not in _SUPPORTED_SUFFIXES:
        raise click.UsageError(f"Unsupported output format '{suffix}'. Use .parquet, .tsv, .csv, or .h5ad for iBAQ")
    output.parent.mkdir(parents=True, exist_ok=True)
    if suffix == ".h5ad":
        _save_ae_anndata(result_df, output, parameters or {})
    elif suffix == ".parquet":
        result_df.to_parquet(str(output), index=False)
    elif suffix == ".tsv":
        result_df.to_csv(str(output), sep="\t", index=False)
    else:
        result_df.to_csv(str(output), index=False)


def _condition_column(columns, factor: Optional[str]) -> Optional[str]:
    """Return the SDRF factor column to use as the condition; several without ``--factor`` warn on the first."""
    factor_columns = [column for column in columns if column.startswith("factor value[")]
    if factor is not None:
        wanted = factor if factor.startswith("factor value[") else f"factor value[{factor.lower()}]"
        if wanted not in factor_columns:
            raise click.UsageError(f"--factor {factor!r} is not in the SDRF; it has: {', '.join(factor_columns) or 'none'}")
        return wanted
    if len(factor_columns) > 1:
        click.echo(
            f"WARNING: the SDRF has {len(factor_columns)} factor columns; using {factor_columns[0]!r} as the "
            "condition. Choose one with --factor.",
            err=True,
        )
    return factor_columns[0] if factor_columns else None


def _quantify_sdrf_mapping(sdrf_path: Path, factor: Optional[str] = None):
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

    condition_column = _condition_column(sdrf.columns, factor)
    conditions = sdrf[condition_column].astype("string").str.strip() if condition_column else sample_ids
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


def _sdrf_samples(pairs, sdrf_path: Path, factor: Optional[str] = None):
    """Resolve each distinct ``(run, label)`` pair to its SDRF sample and condition."""
    import pandas as pd

    from qpx.converters.channel_labels import normalize_label
    from qpx.core.files import run_file_stem

    resolved = pd.DataFrame(
        {
            "_run_key": pairs["run_file_name"].map(run_file_stem),
            "_label_key": pairs["label"].map(normalize_label),
        }
    ).merge(
        _quantify_sdrf_mapping(sdrf_path, factor),
        on=["_run_key", "_label_key"],
        how="left",
        validate="many_to_one",
    )
    unmapped = resolved["SampleID"].isna()
    if unmapped.any():
        missing = resolved.loc[unmapped, ["_run_key", "_label_key"]].drop_duplicates()
        rendered = ", ".join(f"{run_key}::{label_key}" for run_key, label_key in missing.itertuples(index=False, name=None))
        raise click.UsageError(f"QPX intensities have no matching SDRF run/label row: {rendered}")
    return resolved["SampleID"].array, resolved["Condition"].array


def _feature_files(feature_path: Path) -> tuple[list[Path], bool]:
    """Return the Parquet files of a feature input, and whether its run names are partition values."""
    if feature_path.is_file():
        return [feature_path], False
    flat = sorted(feature_path.glob("*.feature.parquet"))
    return (flat, False) if flat else (sorted(feature_path.rglob("*.parquet")), True)


_PEPTIDE_ENTRIES_SQL = """
    SELECT anchor_protein AS "ProteinName", sequence AS "PeptideCanonical", run_file_name,
           coalesce(entry.label, '') AS label, SUM(CAST(entry.intensity AS DOUBLE)) AS "NormIntensity"
    FROM (
        SELECT anchor_protein, sequence, run_file_name, unnest(intensities) AS entry
        FROM read_parquet($files, union_by_name=true $hive)
        WHERE anchor_protein <> '' $target_only
    )
    WHERE entry.intensity > 0
    GROUP BY ALL
    ORDER BY ALL
"""


def _qpx_feature_to_peptide_df(feature_path: Path, sdrf_path: Optional[Path] = None, factor: Optional[str] = None):
    """Read QPX feature data into the peptide DataFrame mokume expects.

    The columns are ``ProteinName``, ``PeptideCanonical``, ``NormIntensity``,
    ``SampleID`` and ``Condition``. *feature_path* is a feature file, a dataset
    directory, or a directory partitioned by ``run_file_name``.

    DuckDB flattens the intensity entries and sums the charge states and
    peptidoforms of a peptide in each ``(run, label)``, which mokume's methods
    would otherwise do on a much larger frame. With an SDRF each pair maps to its
    real sample and condition. Without one the ``run::label`` sample key is kept
    and the condition is marked ``not available``.
    """
    import pandas as pd
    import pyarrow.parquet as pq

    from qpx.core.engine import create_connection
    from qpx.core.query import _escape_sql_string
    from qpx.core.sql import sql_build
    from qpx.version import check_feature_schema_compatible

    paths, partitioned = _feature_files(feature_path)
    if not paths:
        raise click.UsageError(f"No feature Parquet files found in {feature_path}")
    schema = pq.read_schema(paths[0])
    check_feature_schema_compatible(schema, feature_path)
    query = sql_build(
        _PEPTIDE_ENTRIES_SQL,
        files="[" + ", ".join(f"'{_escape_sql_string(str(path))}'" for path in paths) + "]",
        hive=", hive_partitioning=true, hive_types={'run_file_name': 'VARCHAR'}" if partitioned else "",
        target_only="AND NOT coalesce(is_decoy, false)" if "is_decoy" in schema.names else "",
    )
    with create_connection() as con:
        table = con.sql(query).to_arrow_table()
    logger.info("Loaded %d peptide intensities from %s", table.num_rows, feature_path)
    if table.num_rows == 0:
        return pd.DataFrame(columns=["ProteinName", "PeptideCanonical", "NormIntensity", "SampleID", "Condition"])

    frame = table.to_pandas()
    pairs = frame[["run_file_name", "label"]].drop_duplicates(ignore_index=True)
    if sdrf_path is not None:
        sample_ids, conditions = _sdrf_samples(pairs, sdrf_path, factor)
    else:
        sample_ids = pd.array(
            [run if not label else f"{run}::{label}" for run, label in pairs.itertuples(index=False, name=None)], dtype=object
        )
        conditions = None
    rows = frame.groupby(["run_file_name", "label"], sort=False).ngroup().to_numpy()
    frame["SampleID"] = sample_ids.take(rows)
    frame["Condition"] = conditions.take(rows) if conditions is not None else "not available"
    result = frame.drop(columns=["run_file_name", "label"])

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
    help="QPX feature.parquet file, a dataset directory, or a directory partitioned by run_file_name",
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
    help="Parallel threads for MaxLFQ and DirectLFQ (-1 = all cores)",
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
@click.option(
    "--factor",
    help="SDRF factor used as the condition, e.g. 'disease' for 'factor value[disease]' (default: the first factor column)",
    default=None,
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
            --feature-path ./qpx_output/diann.feature.parquet \\
            --method directlfq \\
            -o proteins_directlfq.parquet

    \b
        # iBAQ (requires FASTA)
        qpxc transform quantify \\
            --feature-path ./qpx_output/diann.feature.parquet \\
            --method ibaq --fasta proteome.fasta --sdrf experiment.sdrf.tsv \\
            -o proteins_ibaq.tsv

    \b
        # MaxLFQ with 8 threads
        qpxc transform quantify \\
            --feature-path ./qpx_output/diann.feature.parquet \\
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
    peptide_df = _qpx_feature_to_peptide_df(feature_path, sdrf, kwargs["factor"])
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
