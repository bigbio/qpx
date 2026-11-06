"""MuData export — assemble QPX dataset into a MuData object.

Each modality (precursors, proteins, expression, differential) becomes
an AnnData layer inside the MuData container.  Feature-to-protein
mapping is stored as a boolean sparse adjacency matrix in
``mdata.varp["feature_mapping"]``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import duckdb
import numpy as np
import pandas as pd
from scipy import sparse

from qpx._version import __version__
from qpx.core.mudata_obs import UNIT_SEPARATOR, _has_view, _load_run_obs, _load_unit_obs
from qpx.core.query import _escape_sql_string
from qpx.core.sql import sql_build
from qpx.version import QPX_SPEC_VERSION

if TYPE_CHECKING:
    from qpx.core.engine import DuckDBEngine
    from qpx.dataset import Dataset

logger = logging.getLogger(__name__)

_VALID_MODALITIES = {"precursors", "proteins", "expression", "differential"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _pivot_to_sparse(
    df: pd.DataFrame,
    row_col: str,
    col_col: str,
    value_col: str,
    row_index: pd.Index,
    col_index: pd.Index,
) -> sparse.csr_matrix:
    """Pivot a long-form DataFrame into a sparse CSR matrix.

    Parameters
    ----------
    df : DataFrame
        Long-form data with at least *row_col*, *col_col*, *value_col*.
    row_col, col_col, value_col : str
        Column names for rows, columns, and values.
    row_index, col_index : pd.Index
        Pre-built indices that define the matrix shape.

    Returns
    -------
    scipy.sparse.csr_matrix
    """
    row_pos = row_index.get_indexer(df[row_col])
    col_pos = col_index.get_indexer(df[col_col])

    # Drop entries that did not match (-1)
    mask = (row_pos >= 0) & (col_pos >= 0)
    data = df[value_col].values[mask].astype(np.float32)
    row_pos = row_pos[mask]
    col_pos = col_pos[mask]

    mat = sparse.coo_matrix(
        (data, (row_pos, col_pos)),
        shape=(len(row_index), len(col_index)),
    )
    return mat.tocsr()


def _dictionary_positions(column, index: pd.Index) -> np.ndarray:
    """Map an Arrow string column onto positions in *index*, via its dictionary.

    ``pd.Index.get_indexer`` hashes every element, so on a long table it hashes
    one Python string per row. Dictionary-encoding first means the hashing runs
    once per *distinct* value — for a feature table with hundreds of millions of
    rows over tens of thousands of precursors that is a different order of
    magnitude of work, and it never materialises the row-wise Python strings.

    Unmatched values keep ``get_indexer``'s -1, so callers filter exactly as before.
    """
    import pyarrow.compute as pc

    encoded = pc.dictionary_encode(_combine_chunks(column))
    dictionary_positions = index.get_indexer(encoded.dictionary.to_pylist())
    codes = encoded.indices.to_numpy(zero_copy_only=False)
    if not encoded.indices.null_count:
        return dictionary_positions[codes]

    valid = ~pd.isna(codes)
    safe_codes = np.where(valid, codes, 0).astype(np.intp, copy=False)
    positions = np.full(len(codes), -1, dtype=np.intp)
    positions[valid] = dictionary_positions[safe_codes[valid]]
    return positions


def _pivot_arrow_to_sparse(
    table,
    row_col: str,
    col_col: str,
    value_col: str,
    row_index: pd.Index,
    col_index: pd.Index,
) -> sparse.csr_matrix:
    """Arrow-native equivalent of :func:`_pivot_to_sparse`.

    Same result, without converting the long table to pandas object columns first.
    """
    matrix, _ = _pivot_arrow_combining_duplicates(table, row_col, col_col, value_col, (row_index, col_index), "sum")
    return matrix


def _pivot_arrow_combining_duplicates(
    table,
    row_col: str,
    col_col: str,
    value_col: str,
    axes: tuple[pd.Index, pd.Index],
    duplicates: str,
) -> tuple[sparse.csr_matrix, int]:
    """Pivot like :func:`_pivot_arrow_to_sparse` and count the cells measured twice.

    Several long rows can land in one matrix cell, e.g. a MaxQuant precursor
    quantified by an MS/MS and a match-between-runs feature in one run. They are
    combined with *duplicates* (``"sum"`` or ``"mean"``) and the number of such
    cells is returned, so the caller can report it instead of summing silently.
    """
    row_index, col_index = axes
    row_pos = _dictionary_positions(table.column(row_col), row_index)
    col_pos = _dictionary_positions(table.column(col_col), col_index)

    mask = (row_pos >= 0) & (col_pos >= 0)
    values = table.column(value_col).to_numpy(zero_copy_only=False)
    data = values[mask].astype(np.float32)
    positions = (row_pos[mask], col_pos[mask])
    shape = (len(row_index), len(col_index))

    # tocsr() sums duplicate coordinates; fewer stored cells than rows means some
    # cells got more than one value.
    matrix = sparse.coo_matrix((data, positions), shape=shape).tocsr()
    if matrix.nnz == len(data):
        return matrix, 0
    # Same coordinates, so the counts line up with matrix.data cell for cell.
    counts = sparse.coo_matrix((np.ones(len(data), dtype=np.float32), positions), shape=shape).tocsr()
    if duplicates == "mean":
        matrix.data /= counts.data
    return matrix, int(np.count_nonzero(counts.data > 1))


def _combine_chunks(column):
    """Concatenate an Arrow ChunkedArray's chunks, widening strings first.

    ``combine_chunks`` on a ``string``/``binary`` column builds one contiguous
    array with **int32** offsets, so it raises

        offset overflow while concatenating arrays, consider casting input from
        `string` to `large_string` first

    once the column's total bytes exceed 2 GiB. That is reachable on real data:
    a TMT feature table is one row per precursor *per channel*, so the
    run_file_name column repeats a long filename ten-plus times per precursor.
    MSV000085836 (552 runs x 10 channels) tripped it and MuData silently came
    back with no ``precursors`` modality (bigbio/qpx#316).

    ``large_string`` uses int64 offsets and is otherwise equivalent, so widen
    before combining rather than after — after is too late, the concat is what
    overflows.
    """
    import pyarrow as pa
    import pyarrow.compute as pc

    if isinstance(column.type, pa.DataType):
        if pa.types.is_string(column.type):
            column = pc.cast(column, pa.large_string())
        elif pa.types.is_binary(column.type):
            column = pc.cast(column, pa.large_binary())
    return column.combine_chunks()


def _sorted_unique(table, column: str) -> pd.Index:
    """Sorted distinct values of an Arrow column, named after the column."""
    import pyarrow.compute as pc

    values = pc.dictionary_encode(_combine_chunks(table.column(column))).dictionary.to_pylist()
    return pd.Index(sorted(v for v in values if v is not None), name=column)


# pg is flattened since 1.1 (scalar ``label``), so its label queries read the
# column directly; feature keeps the intensities list<struct>.
_LABEL_QUERIES: dict[str, str] = {
    "feature": "SELECT i.label FROM (SELECT UNNEST(intensities) AS i FROM feature) LIMIT 1",
    "pg": "SELECT label FROM pg WHERE label IS NOT NULL LIMIT 1",
}

_MULTIPLEXED_QUERIES: dict[str, str] = {
    "feature": "SELECT COUNT(DISTINCT i.label) FROM (SELECT UNNEST(intensities) AS i FROM feature)",
    "pg": "SELECT COUNT(DISTINCT label) FROM pg WHERE label IS NOT NULL",
}


_PRECURSOR_QUERY = """
    SELECT run_file_name, precursor_id, i.intensity
    FROM (
        SELECT f.run_file_name,
               f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
               UNNEST(f.intensities) AS i
        FROM feature f$decoy_filter
    )
    WHERE i.label = '$label'
    """

_PRECURSOR_ALL_QUERY = """
    SELECT run_file_name,
           CAST(i.label AS VARCHAR) AS intensity_label,
           precursor_id, i.intensity
    FROM (
        SELECT f.run_file_name,
               f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
               UNNEST(f.intensities) AS i
        FROM feature f$decoy_filter
    )
    WHERE i.label IS NOT NULL
    """

# A protein group is quantified once per quantification unit -- the grouped_runs
# set of fractions aggregated together -- and label (docs/spec/data-model.md), so
# the unit is the observation, named by its distinct member runs (sorted, joined
# by UNIT_SEPARATOR; a single-run unit keeps its run_file_name). Exploding
# grouped_runs instead gave every fraction an observation holding a copy of the
# unit's value, which multiplied any sum over observations by the fraction count.
# Since QPX 1.1 the pg view is flattened: scalar ``label`` / ``intensity`` columns,
# one row per label. ID-only rows (null label) are excluded by the label predicate,
# and rows without grouped_runs (empty or null) belong to no unit.
_PROTEIN_QUERY = """
    SELECT array_to_string(list_sort(list_distinct(pg.grouped_runs)), '$separator') AS grouped_runs,
           pg.anchor_protein,
           pg.intensity AS intensity
    FROM pg
    WHERE pg.label = '$label'
      AND len(pg.grouped_runs) > 0$decoy_filter
    """

_PROTEIN_ALL_QUERY = """
    SELECT array_to_string(list_sort(list_distinct(pg.grouped_runs)), '$separator') AS grouped_runs,
           CAST(pg.label AS VARCHAR) AS intensity_label,
           pg.anchor_protein,
           pg.intensity AS intensity
    FROM pg
    WHERE pg.label IS NOT NULL
      AND len(pg.grouped_runs) > 0$decoy_filter
    """

# Decoys are left out unless asked for: pg.md requires consumers to filter them,
# and a MuData matrix has no row-level flag a consumer could filter on.
_DECOY_FILTERS = {
    "feature": "\n        WHERE f.is_decoy IS NOT TRUE",
    "pg": "\n      AND pg.is_decoy IS NOT TRUE",
}

_DECOY_VAR_QUERIES = {
    "feature": """
        SELECT f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS name, bool_or(f.is_decoy) AS is_decoy
        FROM feature f GROUP BY 1""",
    "pg": "SELECT anchor_protein AS name, bool_or(is_decoy) AS is_decoy FROM pg GROUP BY 1",
}

_DUPLICATE_AGGREGATIONS = ("sum", "mean")


def _detect_intensity_label(engine: DuckDBEngine, table: str = "feature") -> str | None:
    """Return the first intensity label of the given table, or None when it has none."""
    rows = engine.execute(_LABEL_QUERIES[table]).fetchall()
    return rows[0][0] if rows else None


def _prepare_observations_arrow(
    table,
    intensity_label: str | None,
    key_col: str = "run_file_name",
):
    """Create matrix row IDs and their *key_col* (run or unit) / label lookup keys.

    Returns ``(table, row_col, obs_index, observation_keys)``. When every label is
    requested the ``observation_id`` column is appended to the Arrow table rather
    than built row-by-row in pandas: the old path copied the whole frame and then
    concatenated two object columns for every row, which on a large experiment is
    hundreds of millions of throwaway Python strings.
    """
    import pyarrow as pa
    import pyarrow.compute as pc

    if intensity_label is None:
        runs = _combine_chunks(table.column(key_col))
        labels = _combine_chunks(table.column("intensity_label"))
        observation_id = pc.binary_join_element_wise(
            pc.cast(runs, pa.large_string()),
            pc.cast(labels, pa.large_string()),
            pa.scalar("|", pa.large_string()),
        )
        table = table.append_column("observation_id", observation_id)
        obs_index = _sorted_unique(table, "observation_id")
        # Reduce to one row per observation in Arrow *before* touching pandas -
        # the distinct set is the number of run/label pairs, not the row count.
        distinct = (
            table.select(["observation_id", key_col, "intensity_label"])
            .group_by(["observation_id", key_col, "intensity_label"])
            .aggregate([])
        )
        pairs = (
            distinct.to_pandas()
            .drop_duplicates(subset="observation_id", keep="first")
            .set_index("observation_id")
            .reindex(obs_index)
        )
        return table, "observation_id", obs_index, pairs

    obs_index = _sorted_unique(table, key_col)
    keys = pd.DataFrame(
        {
            key_col: obs_index,
            "intensity_label": str(intensity_label),
        },
        index=obs_index,
    )
    return table, key_col, obs_index, keys


def _attach_uns_metadata(engine: DuckDBEngine, mdata) -> None:
    """Attach project metadata from dataset.parquet to MuData.uns."""
    try:
        df = engine.execute("SELECT * FROM dataset LIMIT 1").fetchdf()
    except duckdb.Error:
        logger.debug("No dataset table found; skipping dataset metadata.")
    else:
        if not df.empty:
            row = df.iloc[0]
            for col in df.columns:
                val = row[col]
                if val is None:
                    continue
                # pd.isna() catches np.nan, pd.NA, pd.NaT on scalars; guarded by is_scalar
                # so that dicts/lists/arrays do not raise or get misclassified.
                if pd.api.types.is_scalar(val) and pd.isna(val):
                    continue
                mdata.uns[col] = val

    mdata.uns["qpx_version"] = QPX_SPEC_VERSION
    mdata.uns["writer_version"] = __version__


# ---------------------------------------------------------------------------
# Modality builders
# ---------------------------------------------------------------------------

_LONG_QUERIES = {
    "feature": (_PRECURSOR_QUERY, _PRECURSOR_ALL_QUERY),
    "pg": (_PROTEIN_QUERY, _PROTEIN_ALL_QUERY),
}


def _long_table(engine: DuckDBEngine, table: str, intensity_label: str | None, include_decoys: bool):
    """Read the long precursor (``feature``) or protein (``pg``) table for one label or all."""
    one_label, all_labels = _LONG_QUERIES[table]
    decoy_filter = "" if include_decoys else _DECOY_FILTERS[table]
    if intensity_label is None:
        sql = sql_build(all_labels, decoy_filter=decoy_filter, separator=UNIT_SEPARATOR)
    else:
        label = _escape_sql_string(intensity_label)
        sql = sql_build(one_label, label=label, decoy_filter=decoy_filter, separator=UNIT_SEPARATOR)
    return engine.connection.sql(sql).to_arrow_table()


def _pivot_modality(name: str, table, row_col: str, var_col: str, axes, duplicates: str):
    """Pivot a modality's intensities and report cells that received several values."""
    matrix, n_cells = _pivot_arrow_combining_duplicates(table, row_col, var_col, "intensity", axes, duplicates)
    if n_cells:
        logger.warning(
            "MuData %s: %d cells received more than one intensity (several features or groups for one "
            "observation); combined with %r. Pass duplicates= to choose.",
            name,
            n_cells,
            duplicates,
        )
    return matrix, {"aggregation": duplicates, "cells": n_cells}


def _decoy_flags(engine: DuckDBEngine, table: str, names: pd.Index) -> np.ndarray:
    """Return ``is_decoy`` for each var name (precursor id or anchor protein)."""
    flags = engine.execute(_DECOY_VAR_QUERIES[table]).fetchdf().set_index("name")["is_decoy"]
    return flags.reindex(names).fillna(False).astype(bool).to_numpy()


def _build_precursor_adata(
    engine: DuckDBEngine,
    intensity_label: str | None,
    include_decoys: bool = False,
    duplicates: str = "sum",
):
    """Build AnnData from feature.parquet.

    obs = runs (one label) or run/label pairs (all labels),
    var = peptidoform|charge (with is_decoy when decoys are included),
    X = sparse intensity matrix.
    """
    import anndata as ad

    table = _long_table(engine, "feature", intensity_label, include_decoys)
    if table.num_rows == 0:
        return ad.AnnData()

    table, row_col, observations, observation_keys = _prepare_observations_arrow(table, intensity_label)
    precursors = _sorted_unique(table, "precursor_id")

    intensity_matrix, duplicate_info = _pivot_modality(
        "precursors", table, row_col, "precursor_id", (observations, precursors), duplicates
    )

    obs = _load_run_obs(engine, observation_keys, all_labels=intensity_label is None)
    var = pd.DataFrame(index=precursors)
    if include_decoys:
        var["is_decoy"] = _decoy_flags(engine, "feature", precursors)

    return ad.AnnData(X=intensity_matrix, obs=obs, var=var, uns={"qpx_duplicates": duplicate_info})


def _build_protein_adata(
    engine: DuckDBEngine,
    intensity_label: str | None,
    include_decoys: bool = False,
    duplicates: str = "sum",
):
    """Build AnnData from pg.parquet.

    obs = quantification units (one label) or unit/label pairs (all labels),
    var = anchor_protein, X = sparse intensity matrix.
    var includes gene_name column (and is_decoy when decoys are included).
    """
    import anndata as ad

    table = _long_table(engine, "pg", intensity_label, include_decoys)
    if table.num_rows == 0:
        return ad.AnnData()

    table, row_col, observations, observation_keys = _prepare_observations_arrow(table, intensity_label, key_col="grouped_runs")
    proteins = _sorted_unique(table, "anchor_protein")

    intensity_matrix, duplicate_info = _pivot_modality(
        "proteins", table, row_col, "anchor_protein", (observations, proteins), duplicates
    )

    obs = _load_unit_obs(engine, observation_keys, all_labels=intensity_label is None)

    # Gene names
    gene_sql = """
    SELECT DISTINCT pg.anchor_protein,
           pg.gg_names[1] AS gene_name
    FROM pg
    """
    gene_df = engine.execute(gene_sql).fetchdf()
    gene_df = gene_df.drop_duplicates(subset="anchor_protein", keep="first")
    gene_df = gene_df.set_index("anchor_protein").reindex(proteins)

    # A missing gene must still serialise as a string. When every protein lacks a
    # gene the column comes back float64 (all NaN) and h5py refuses to write var
    # ("Can't implicitly convert non-string objects to strings"). That is the norm
    # for TMT datasets whose consensusXML carries no GN= descriptions, so gg_names
    # is NULL on every row. Same reason _load_run_obs fills its string columns.
    gene_names = gene_df["gene_name"].fillna("").astype(str)

    var = pd.DataFrame({"gene_name": gene_names.values}, index=proteins)
    if include_decoys:
        var["is_decoy"] = _decoy_flags(engine, "pg", proteins)

    return ad.AnnData(X=intensity_matrix, obs=obs, var=var, uns={"qpx_duplicates": duplicate_info})


def _build_feature_mapping(engine: DuckDBEngine, mdata) -> sparse.csr_matrix:
    """Build boolean sparse adjacency linking precursor IDs to protein IDs.

    Returns a symmetric CSR matrix of shape (N, N) where N is ``mdata.n_vars``
    — the total number of var_names across **all** modalities present, in their
    global-axis order — so the matrix is assignable to ``mdata.varp`` even when
    expression / differential modalities are also built (bigbio/qpx#252). The
    precursor and protein blocks are positioned by their per-modality offset on
    the global var axis (not by a global name lookup), so a protein accession
    that also appears as an expression var is never mislocated. Non-zero entries
    indicate a precursor-protein relationship.
    """
    if "precursors" not in mdata.mod or "proteins" not in mdata.mod:
        raise ValueError("Both 'precursors' and 'proteins' modalities are required for feature mapping.")

    # Per-modality offsets on the global var axis (modalities are concatenated in
    # mdata.mod insertion order). Sizing/positioning against the full axis is what
    # keeps the mapping valid — and present — with a 3rd/4th modality around.
    offsets: dict[str, int] = {}
    running = 0
    for name, adata in mdata.mod.items():
        offsets[name] = running
        running += adata.n_vars
    n = running  # == mdata.n_vars

    precursor_names = pd.Index(mdata.mod["precursors"].var_names)
    protein_names = pd.Index(mdata.mod["proteins"].var_names)

    # Query mapping from feature table
    sql = """
    SELECT DISTINCT
           f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
           f.anchor_protein
    FROM feature f
    """
    df = engine.execute(sql).fetchdf()

    if df.empty:
        return sparse.csr_matrix((n, n), dtype=bool)

    precursor_local = precursor_names.get_indexer(df["precursor_id"])
    protein_local = protein_names.get_indexer(df["anchor_protein"])

    mask = (precursor_local >= 0) & (protein_local >= 0)
    precursor_pos = precursor_local[mask] + offsets["precursors"]
    protein_pos = protein_local[mask] + offsets["proteins"]

    data = np.ones(len(precursor_pos) * 2, dtype=bool)
    rows = np.concatenate([precursor_pos, protein_pos])
    cols = np.concatenate([protein_pos, precursor_pos])

    mat = sparse.coo_matrix((data, (rows, cols)), shape=(n, n))
    return mat.tocsr()


def _build_expression_adata(ds_path: Path, file_prefix: str | None = None):
    """Auto-detect and load ``.pe.h5ad`` or ``.pe.zarr``.

    Returns None if not found.
    """
    import anndata as ad

    # Try h5ad first
    h5ad_pattern = f"{file_prefix}.pe.h5ad" if file_prefix else "*.pe.h5ad"
    candidates = list(ds_path.glob(h5ad_pattern))
    if candidates:
        return ad.read_h5ad(candidates[0])

    # Try zarr
    zarr_pattern = f"{file_prefix}.pe.zarr" if file_prefix else "*.pe.zarr"
    candidates = list(ds_path.glob(zarr_pattern))
    if candidates:
        return ad.read_zarr(candidates[0])

    return None


def _build_differential_adata(ds_path: Path, file_prefix: str | None = None):
    """Auto-detect and load ``.de.h5ad`` or ``.de.zarr``.

    Restructures ``uns["de_results"]`` into matrix form:
    obs = contrasts, var = proteins, X = log2FC,
    layers = {pvals_adj, scores, pvals, se}.

    Returns None if not found.
    """
    import anndata as ad

    # Try h5ad first
    h5ad_pattern = f"{file_prefix}.de.h5ad" if file_prefix else "*.de.h5ad"
    candidates = list(ds_path.glob(h5ad_pattern))
    if not candidates:
        zarr_pattern = f"{file_prefix}.de.zarr" if file_prefix else "*.de.zarr"
        candidates = list(ds_path.glob(zarr_pattern))
        if not candidates:
            return None
        adata = ad.read_zarr(candidates[0])
    else:
        adata = ad.read_h5ad(candidates[0])

    # If it has de_results in uns, restructure into matrix form
    if "de_results" in adata.uns:
        de = adata.uns["de_results"]
        if isinstance(de, pd.DataFrame) and not de.empty:
            return _reshape_de_results(de)

    # Already in matrix form or no de_results — return as-is
    return adata


def _reshape_de_results(de: pd.DataFrame):
    """Reshape long-form DE results into AnnData matrix.

    Expects columns: contrast, protein, log2FC and optionally
    pvals_adj, scores, pvals, se.
    """
    import anndata as ad

    # Identify contrast and protein columns
    contrast_col = "contrast" if "contrast" in de.columns else de.columns[0]
    protein_col = "protein" if "protein" in de.columns else de.columns[1]

    contrasts = pd.Index(sorted(de[contrast_col].unique()), name="contrast")
    proteins = pd.Index(sorted(de[protein_col].unique()), name="protein")

    if "log2FC" not in de.columns:
        raise ValueError(f"DE results missing required 'log2FC' column. Found columns: {list(de.columns)}")

    X = _pivot_to_sparse(de, contrast_col, protein_col, "log2FC", contrasts, proteins)

    layers = {}
    for col in ("pvals_adj", "scores", "pvals", "se"):
        if col in de.columns:
            layers[col] = _pivot_to_sparse(de, contrast_col, protein_col, col, contrasts, proteins)

    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=contrasts),
        var=pd.DataFrame(index=proteins),
        layers=layers,
    )


# ---------------------------------------------------------------------------
# Top-level assembler
# ---------------------------------------------------------------------------


def _available_quant_modalities(dataset: Dataset) -> set[str]:
    """Return quantification modalities with registered QPX source structures."""
    sources = {"precursors": dataset.feature, "proteins": dataset.pg}
    return {name for name, structure in sources.items() if structure is not None}


def _try_build_modality(name: str, builder, mod: dict, failures: dict[str, str] | None = None) -> None:
    """Call *builder*, add the result to *mod* if non-empty.

    A modality that raises is logged and skipped so one bad view cannot cost the
    caller every other modality. The reason is also recorded in *failures*: a
    warning in a log is not something a caller can branch on, and a MuData that
    is quietly missing a modality is indistinguishable from one that never had
    it. MSV000085836 shipped an h5mu with proteins and no precursors that way
    (bigbio/qpx#316).
    """
    try:
        adata = builder()
        if adata is not None and (not hasattr(adata, "n_obs") or adata.n_obs > 0):
            mod[name] = adata
    except Exception as exc:
        logger.warning("Failed to build %s modality: %s", name, exc)
        if failures is not None:
            failures[name] = f"{type(exc).__name__}: {exc}"


def _is_multiplexed(engine: DuckDBEngine, table: str) -> bool:
    """Return True when *table* carries more than one distinct intensity label.

    Multiplexed data (isobaric TMT/iTRAQ, plexDIA) has many channel labels per
    run; label-free has a single label ("LFQ"). Detection drives whether obs is
    expanded over run x channel or kept at the run level.
    """
    return engine.execute(_MULTIPLEXED_QUERIES[table]).fetchall()[0][0] > 1


def _select_intensity_labels(
    engine: DuckDBEngine,
    intensity_label: str | None,
    all_intensity_labels: bool,
) -> tuple[str | None, str | None]:
    """Select feature/PG labels while preserving their table-specific defaults.

    ``None`` means "expand over every channel" (run x label observations); a
    concrete label means the single-label, run-level contract. Multiplexed
    tables auto-expand so isobaric/plexDIA data never collapses to one channel.
    """
    if all_intensity_labels:
        if intensity_label is not None:
            raise ValueError("intensity_label cannot be combined with all_intensity_labels=True")
        return None, None

    if intensity_label is not None:
        return intensity_label, intensity_label

    # Without a feature label (e.g. a proteins-only build) the precursors expand
    # over every label rather than aborting the whole MuData assembly.
    feature_label = None
    if _has_view(engine, "feature") and not _is_multiplexed(engine, "feature"):
        feature_label = _detect_intensity_label(engine, "feature")
    protein_label = feature_label
    if _has_view(engine, "pg"):
        if _is_multiplexed(engine, "pg"):
            protein_label = None
        else:
            detected = _detect_intensity_label(engine, "pg")
            protein_label = feature_label if detected is None else detected
    return feature_label, protein_label


def _requested_modalities(dataset: Dataset, modalities: list[str] | None) -> set[str]:
    """Validate *modalities*, defaulting to every available one."""
    if modalities is None:
        return _available_quant_modalities(dataset) | {"expression", "differential"}
    unknown = set(modalities) - _VALID_MODALITIES
    if unknown:
        raise ValueError(f"Unknown modalities: {unknown}. Valid: {sorted(_VALID_MODALITIES)}")
    return set(modalities)


def _fill_string_nans(mdata) -> None:
    """Replace NaN in object columns, which h5py cannot write in vlen-string arrays."""
    for adata in mdata.mod.values():
        for frame in (adata.obs, adata.var):
            str_cols = frame.select_dtypes(include=["object"]).columns
            if len(str_cols):
                frame[str_cols] = frame[str_cols].fillna("")


def build_mudata(
    dataset: Dataset,
    intensity_label: str | None = None,
    modalities: list[str] | None = None,
    all_intensity_labels: bool = False,
    include_decoys: bool = False,
    duplicates: str = "sum",
):
    """Assemble QPX dataset into a MuData object.

    Parameters
    ----------
    dataset : Dataset
        An open QPX Dataset instance.
    intensity_label : str, optional
        Intensity label to extract (e.g. "TMT126").  When omitted, label-free
        data auto-detects its single label, while multiplexed (isobaric/plexDIA)
        data expands over every channel as ``run|label`` (precursors) and
        ``grouped_runs|label`` (proteins) observations.
    modalities : list of str, optional
        Which modalities to include.  Valid values:
        ``"precursors"``, ``"proteins"``, ``"expression"``, ``"differential"``.
        Defaults to all available.
    all_intensity_labels : bool
        Force every intensity label into a separate ``run|label`` (precursors)
        or ``grouped_runs|label`` (proteins) observation; obs then carries the
        ``run_file_name`` or ``grouped_runs`` and ``intensity_label`` columns.
        Defaults to False: label-free data then keeps one observation per run
        (per unit for proteins), while multiplexed data still auto-expands over
        its channels (so isobaric/plexDIA never collapses to a single channel).
    include_decoys : bool
        Keep decoy precursors and protein groups, flagged in ``var["is_decoy"]``.
        Defaults to False: decoys are left out of both matrices.
    duplicates : {"sum", "mean"}
        How to combine several intensities that land in one cell, such as two
        features of one precursor in one run. The count of such cells is
        stored in each modality's ``uns["qpx_duplicates"]`` and logged.

    Returns
    -------
    mudata.MuData
        ``proteins`` observations are quantification units (one ``grouped_runs``
        set, per label when expanded), named by their member runs joined by
        ``;``; a single-run unit keeps its run name.
    """
    import mudata as mu

    if duplicates not in _DUPLICATE_AGGREGATIONS:
        raise ValueError(f"duplicates must be one of {_DUPLICATE_AGGREGATIONS}, got {duplicates!r}")
    engine = dataset._engine
    ds_path = Path(dataset.path)
    file_prefix = dataset.file_prefix
    requested = _requested_modalities(dataset, modalities)
    feat_intensity, pg_intensity = _select_intensity_labels(engine, intensity_label, all_intensity_labels)

    mod: dict = {}
    builders = {
        "precursors": lambda: _build_precursor_adata(engine, feat_intensity, include_decoys, duplicates),
        "proteins": lambda: _build_protein_adata(engine, pg_intensity, include_decoys, duplicates),
        "expression": lambda: _build_expression_adata(ds_path, file_prefix),
        "differential": lambda: _build_differential_adata(ds_path, file_prefix),
    }
    failures: dict[str, str] = {}
    for name, builder in builders.items():
        if name in requested:
            _try_build_modality(name, builder, mod, failures)

    mdata = mu.MuData(mod)

    # Surface build failures on the object itself, so a caller can detect an
    # incomplete view without scraping the log. Always set the key, so its
    # absence means "built by an older qpx", not "nothing failed".
    mdata.uns["qpx_failed_modalities"] = dict(failures)
    mdata.uns["qpx_include_decoys"] = bool(include_decoys)

    # Attach feature mapping if both precursors and proteins are present
    if "precursors" in mod and "proteins" in mod:
        try:
            mapping = _build_feature_mapping(engine, mdata)
            mdata.varp["feature_mapping"] = mapping
        except Exception as exc:
            logger.warning("Failed to build feature mapping: %s", exc)

    # Attach dataset-level metadata
    _attach_uns_metadata(engine, mdata)
    _fill_string_nans(mdata)
    return mdata


def write_dataset_mudata(
    output_folder: Path | str,
    prefix: str,
    *,
    all_intensity_labels: bool = True,
) -> Path | None:
    """Assemble and write the MuData (.h5mu) view of a QPX dataset on disk.

    Shared by the converters (which refresh the view right after writing the
    parquet) and by ``qpxc transform`` commands that rewrite those parquet, so a
    stale h5mu never outlives the data it describes.

    Best-effort: expected dependency, database, validation, and I/O failures are
    logged and skipped rather than raised, because the parquet views — not the
    h5mu — are the dataset's source of truth.
    """
    output_folder = Path(output_folder)
    h5mu_path = output_folder / f"{prefix}.h5mu"
    h5mu_tmp = h5mu_path.with_name(f"{h5mu_path.stem}.tmp{h5mu_path.suffix}")
    try:
        # Core parquet files have already been refreshed, so an older MuData
        # view would be stale if this build fails.
        h5mu_path.unlink(missing_ok=True)

        from qpx.dataset import Dataset

        dataset = Dataset(str(output_folder), file_prefix=prefix)
        try:
            required_modalities = _available_quant_modalities(dataset)
            if not required_modalities:
                logger.info("Skipping muData for %s: no feature or pg quantification data", prefix)
                return None

            modalities = required_modalities | {"expression", "differential"}
            mdata = build_mudata(
                dataset,
                modalities=sorted(modalities),
                all_intensity_labels=all_intensity_labels,
            )
            missing = required_modalities - set(mdata.mod)
            if missing:
                # Say why, not just what: the build failure is the actionable part.
                reasons = getattr(mdata, "uns", None) or {}
                reasons = reasons.get("qpx_failed_modalities") or {}
                detail = "; ".join(f"{name} ({reasons[name]})" if name in reasons else name for name in sorted(missing))
                raise ValueError(f"Missing required quantification modalities: {detail}")
            mdata.write(str(h5mu_tmp))
            h5mu_tmp.replace(h5mu_path)
        finally:
            dataset.close()
    except (ImportError, OSError, RuntimeError, TypeError, ValueError, duckdb.Error) as exc:
        logger.warning("Could not build muData for %s: %s", prefix, exc)
        return None
    finally:
        try:
            h5mu_tmp.unlink(missing_ok=True)
        except OSError as exc:
            logger.warning("Could not remove temporary muData file %s: %s", h5mu_tmp, exc)
    logger.info("Wrote muData to %s", h5mu_path)
    return h5mu_path
