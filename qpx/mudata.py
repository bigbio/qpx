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
from qpx.core.mudata_obs import _has_view, _load_run_obs
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
    """Pivot long rows into a sparse matrix, working on Arrow columns.

    Same result, without converting the long table to pandas object columns first.
    """
    row_pos = _dictionary_positions(table.column(row_col), row_index)
    col_pos = _dictionary_positions(table.column(col_col), col_index)

    values = table.column(value_col).to_numpy(zero_copy_only=False).astype(np.float32)
    # A cell without a stored value is not measured, so a null quantity stores nothing
    mask = (row_pos >= 0) & (col_pos >= 0) & ~np.isnan(values)
    data = values[mask]

    mat = sparse.coo_matrix(
        (data, (row_pos[mask], col_pos[mask])),
        shape=(len(row_index), len(col_index)),
    )
    return mat.tocsr()


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


# Long rows of the precursors modality. The pivot sums the rows of one cell: the
# features of one precursor in one run and label, as MaxQuant sums the evidence
# rows of a peptide in a run. Decoys and entries without a quantity are left out,
# so an absent cell is not quantified. $label selects one label or every label.
# UNNEST stays in the select list: a filter on a lateral CROSS JOIN UNNEST runs ~30x slower.
_PRECURSOR_QUERY = """
    SELECT run_file_name, CAST(i.label AS VARCHAR) AS intensity_label, precursor_id, i.intensity
    FROM (
        SELECT f.run_file_name,
               f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
               UNNEST(f.intensities) AS i
        FROM feature f
        WHERE NOT f.is_decoy
    )
    WHERE i.intensity IS NOT NULL AND $label
    """

# One row per additional intensity (normalized, LFQ, iBAQ), which become layers.
_PRECURSOR_LAYER_QUERY = """
    SELECT run_file_name, CAST(a.label AS VARCHAR) AS intensity_label, precursor_id,
           x.intensity_name, x.intensity_value AS intensity
    FROM (
        SELECT run_file_name, precursor_id, a, UNNEST(a.intensities) AS x
        FROM (
            SELECT f.run_file_name,
                   f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
                   UNNEST(f.additional_intensities) AS a
            FROM feature f
            WHERE NOT f.is_decoy
        )
    )
    WHERE x.intensity_value IS NOT NULL AND $label
    """

_PRECURSOR_IMPUTED_QUERY = """
    SELECT run_file_name, intensity_label, precursor_id, 1.0::FLOAT AS intensity
    FROM (
        SELECT f.run_file_name,
               f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
               UNNEST(f.imputed_labels) AS intensity_label
        FROM feature f
        WHERE NOT f.is_decoy
    )
    WHERE $label
    """

# A protein group is quantified once for its unit, the runs in grouped_runs (the
# fractions of one sample or plex), so an observation of the proteins modality is
# a unit, named after its first run. Variables are groups keyed by their sorted
# membership, since two groups can be led by the same protein. The pg view is
# flattened: one row per label, scalar label and intensity.
_UNIT_SQL = "list_sort(list_distinct(p.grouped_runs))[1]"
_GROUP_SQL = "array_to_string(list_sort(p.pg_accessions), ';')"

_PROTEIN_QUERY = f"""
    SELECT {_UNIT_SQL} AS run_file_name,
           CAST(p.label AS VARCHAR) AS intensity_label,
           {_GROUP_SQL} AS protein_group,
           p.intensity
    FROM pg p
    WHERE NOT p.is_decoy AND p.intensity IS NOT NULL AND $label
    """

_PROTEIN_IMPUTED_QUERY = f"""
    SELECT {_UNIT_SQL} AS run_file_name,
           CAST(p.label AS VARCHAR) AS intensity_label,
           {_GROUP_SQL} AS protein_group,
           1.0::FLOAT AS intensity
    FROM pg p
    WHERE NOT p.is_decoy AND p.is_imputed AND $label
    """

_PROTEIN_LAYER_QUERY = f"""
    SELECT run_file_name, intensity_label, protein_group, x.intensity_name, x.intensity_value AS intensity
    FROM (
        SELECT run_file_name, intensity_label, protein_group, UNNEST(a.intensities) AS x
        FROM (
            SELECT {_UNIT_SQL} AS run_file_name,
                   CAST(p.label AS VARCHAR) AS intensity_label,
                   {_GROUP_SQL} AS protein_group,
                   UNNEST(p.additional_intensities) AS a
            FROM pg p
            WHERE NOT p.is_decoy AND $label
        )
        WHERE a.label = intensity_label
    )
    WHERE x.intensity_value IS NOT NULL
    """

_PROTEIN_VAR_QUERY = f"""
    SELECT {_GROUP_SQL} AS protein_group,
           min(p.anchor_protein) AS anchor_protein,
           min(p.gg_names[1]) AS gene_name
    FROM pg p
    WHERE NOT p.is_decoy
    GROUP BY ALL
    """

_UNIT_RUNS_QUERY = f"""
    SELECT {_UNIT_SQL} AS unit,
           array_to_string(list_sort(list_distinct(flatten(list(p.grouped_runs)))), ';') AS runs
    FROM pg p
    WHERE NOT p.is_decoy
    GROUP BY ALL
    """

# Columns of a run, not of the unit the proteins modality observes.
_RUN_ONLY_OBS_COLUMNS = ["run_accession", "fraction"]


def _label_condition(column: str, intensity_label: str | None) -> str:
    """SQL predicate selecting *intensity_label* in *column*, or every label when it is None."""
    if intensity_label is None:
        return f"{column} IS NOT NULL"
    return f"{column} = '{_escape_sql_string(intensity_label)}'"


def _query_arrow(engine: DuckDBEngine, query: str, label_column: str, intensity_label: str | None):
    """Run one of the modality queries for the selected labels as an Arrow table."""
    sql = sql_build(query, label=_label_condition(label_column, intensity_label))
    return engine.connection.sql(sql).to_arrow_table()


def _detect_intensity_label(engine: DuckDBEngine, table: str = "feature") -> str | None:
    """Return the first intensity label of the given table, or None when it has none."""
    rows = engine.execute(_LABEL_QUERIES[table]).fetchall()
    return rows[0][0] if rows else None


def _with_observation_id(table):
    """Append ``observation_id`` (``run|label``) to an Arrow table of every label."""
    import pyarrow as pa
    import pyarrow.compute as pc

    runs = _combine_chunks(table.column("run_file_name"))
    labels = _combine_chunks(table.column("intensity_label"))
    observation_id = pc.binary_join_element_wise(
        pc.cast(runs, pa.large_string()),
        pc.cast(labels, pa.large_string()),
        pa.scalar("|", pa.large_string()),
    )
    return table.append_column("observation_id", observation_id)


def _prepare_observations_arrow(
    table,
    intensity_label: str | None,
):
    """Create matrix row IDs and their run/label lookup keys.

    Returns ``(table, row_col, obs_index, observation_keys)``. When every label is
    requested the ``observation_id`` column is appended to the Arrow table rather
    than built row-by-row in pandas: the old path copied the whole frame and then
    concatenated two object columns for every row, which on a large experiment is
    hundreds of millions of throwaway Python strings.
    """
    if intensity_label is None:
        table = _with_observation_id(table)
        obs_index = _sorted_unique(table, "observation_id")
        # Reduce to one row per observation in Arrow *before* touching pandas -
        # the distinct set is the number of run/label pairs, not the row count.
        distinct = (
            table.select(["observation_id", "run_file_name", "intensity_label"])
            .group_by(["observation_id", "run_file_name", "intensity_label"])
            .aggregate([])
        )
        pairs = (
            distinct.to_pandas()
            .drop_duplicates(subset="observation_id", keep="first")
            .set_index("observation_id")
            .reindex(obs_index)
        )
        return table, "observation_id", obs_index, pairs

    obs_index = _sorted_unique(table, "run_file_name")
    keys = pd.DataFrame(
        {
            "run_file_name": obs_index,
            "intensity_label": str(intensity_label),
        },
        index=obs_index,
    )
    return table, "run_file_name", obs_index, keys


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


def _additional_layers(
    engine: DuckDBEngine,
    query: str,
    label_column: str,
    intensity_label: str | None,
    var_col: str,
    observations: pd.Index,
    variables: pd.Index,
) -> dict[str, sparse.csr_matrix]:
    """One layer per additional intensity name, on the observations and variables of X."""
    import pyarrow.compute as pc

    table = _query_arrow(engine, query, label_column, intensity_label)
    row_col = "run_file_name"
    if intensity_label is None:
        table = _with_observation_id(table)
        row_col = "observation_id"
    layers = {}
    for name in sorted(table.column("intensity_name").unique().to_pylist()):
        rows = table.filter(pc.field("intensity_name") == name)
        layers[name] = _pivot_arrow_to_sparse(rows, row_col, var_col, "intensity", observations, variables)
    return layers


def _imputed_mask_layer(
    engine: DuckDBEngine,
    query: str,
    label_column: str,
    intensity_label: str | None,
    var_col: str,
    observations: pd.Index,
    variables: pd.Index,
) -> dict[str, sparse.csr_matrix]:
    """The ``imputed_mask`` layer, 1 in the cells the tool imputed; empty when none are."""
    table = _query_arrow(engine, query, label_column, intensity_label)
    if table.num_rows == 0:
        return {}
    row_col = "run_file_name"
    if intensity_label is None:
        table = _with_observation_id(table)
        row_col = "observation_id"
    # Summed cells stack the flags of the rows behind them.
    return {"imputed_mask": _pivot_arrow_to_sparse(table, row_col, var_col, "intensity", observations, variables).sign()}


def _build_precursor_adata(
    engine: DuckDBEngine,
    intensity_label: str | None,
):
    """Build AnnData from feature.parquet.

    obs = runs (one label) or run/label pairs (all labels),
    var = peptidoform|charge, X = sparse intensity matrix,
    layers = the additional intensities by name.
    """
    import anndata as ad

    table = _query_arrow(engine, _PRECURSOR_QUERY, "i.label", intensity_label)
    if table.num_rows == 0:
        return ad.AnnData()

    table, row_col, observations, observation_keys = _prepare_observations_arrow(table, intensity_label)
    precursors = _sorted_unique(table, "precursor_id")

    intensity_matrix = _pivot_arrow_to_sparse(table, row_col, "precursor_id", "intensity", observations, precursors)
    layers = _additional_layers(
        engine, _PRECURSOR_LAYER_QUERY, "a.label", intensity_label, "precursor_id", observations, precursors
    )
    if "imputed_labels" in engine.view_columns("feature"):
        layers |= _imputed_mask_layer(
            engine, _PRECURSOR_IMPUTED_QUERY, "intensity_label", intensity_label, "precursor_id", observations, precursors
        )

    obs = _load_run_obs(engine, observation_keys, all_labels=intensity_label is None)
    var = pd.DataFrame(index=precursors)

    return ad.AnnData(X=intensity_matrix, obs=obs, var=var, layers=layers)


def _unit_obs(engine: DuckDBEngine, obs: pd.DataFrame, all_labels: bool) -> pd.DataFrame:
    """Turn the run metadata of each unit's first run into the unit's metadata."""
    obs = obs.drop(columns=[column for column in _RUN_ONLY_OBS_COLUMNS if column in obs.columns])
    units = obs["run_file_name"] if all_labels else obs.index.to_series()
    runs = engine.execute(_UNIT_RUNS_QUERY).fetchdf().set_index("unit")["runs"]
    obs.insert(1 if all_labels else 0, "runs", units.map(runs).fillna("").values)
    if all_labels:
        obs = obs.rename(columns={"run_file_name": "unit"})
    return obs


def _build_protein_adata(
    engine: DuckDBEngine,
    intensity_label: str | None,
):
    """Build AnnData from pg.parquet.

    obs = quantification units (one label) or unit/label pairs (all labels), with
    their member runs in ``runs``; var = protein groups (sorted ``pg_accessions``
    joined by ``;``) with ``anchor_protein`` and ``gene_name``; X = sparse
    intensity matrix; layers = the additional intensities by name.
    """
    import anndata as ad

    table = _query_arrow(engine, _PROTEIN_QUERY, "p.label", intensity_label)
    if table.num_rows == 0:
        return ad.AnnData()

    table, row_col, observations, observation_keys = _prepare_observations_arrow(table, intensity_label)
    groups = _sorted_unique(table, "protein_group")

    intensity_matrix = _pivot_arrow_to_sparse(table, row_col, "protein_group", "intensity", observations, groups)
    layers = _additional_layers(engine, _PROTEIN_LAYER_QUERY, "p.label", intensity_label, "protein_group", observations, groups)
    if "is_imputed" in engine.view_columns("pg"):
        layers |= _imputed_mask_layer(
            engine, _PROTEIN_IMPUTED_QUERY, "p.label", intensity_label, "protein_group", observations, groups
        )

    obs = _load_run_obs(engine, observation_keys, all_labels=intensity_label is None)
    obs = _unit_obs(engine, obs, all_labels=intensity_label is None)

    # A missing gene must still serialise as a string. When every protein lacks a
    # gene the column comes back float64 (all NaN) and h5py refuses to write var
    # ("Can't implicitly convert non-string objects to strings"). That is the norm
    # for TMT datasets whose consensusXML carries no GN= descriptions, so gg_names
    # is NULL on every row. Same reason _load_run_obs fills its string columns.
    var = engine.execute(_PROTEIN_VAR_QUERY).fetchdf().set_index("protein_group").reindex(groups)
    var = var.fillna("").astype(str)

    return ad.AnnData(X=intensity_matrix, obs=obs, var=var, layers=layers)


def _build_feature_mapping(engine: DuckDBEngine, mdata) -> sparse.csr_matrix:
    """Build boolean sparse adjacency linking precursor IDs to protein groups.

    Returns a symmetric CSR matrix of shape (N, N) where N is ``mdata.n_vars``
    — the total number of var_names across **all** modalities present, in their
    global-axis order — so the matrix is assignable to ``mdata.varp`` even when
    expression / differential modalities are also built (bigbio/qpx#252). The
    precursor and protein blocks are positioned by their per-modality offset on
    the global var axis (not by a global name lookup), so a protein accession
    that also appears as an expression var is never mislocated. Non-zero entries
    indicate that a feature of the precursor belongs to the group
    (``feature.pg_accessions``).
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
           array_to_string(list_sort(list_transform(f.pg_accessions, lambda x: x.accession)), ';') AS protein_group
    FROM feature f
    WHERE NOT f.is_decoy
    """
    df = engine.execute(sql).fetchdf()

    if df.empty:
        return sparse.csr_matrix((n, n), dtype=bool)

    precursor_local = precursor_names.get_indexer(df["precursor_id"])
    protein_local = protein_names.get_indexer(df["protein_group"])

    mask = (precursor_local >= 0) & (protein_local >= 0)
    precursor_pos = precursor_local[mask] + offsets["precursors"]
    protein_pos = protein_local[mask] + offsets["proteins"]

    data = np.ones(len(precursor_pos) * 2, dtype=bool)
    rows = np.concatenate([precursor_pos, protein_pos])
    cols = np.concatenate([protein_pos, precursor_pos])

    mat = sparse.coo_matrix((data, (rows, cols)), shape=(n, n))
    return mat.tocsr()


def _build_expression_adata(ds_path: Path, file_prefix: str | None = None):
    """Auto-detect and load the absolute-expression AnnData: ``.ae``, or the older ``.pe``, as h5ad or zarr.

    Returns None if not found.
    """
    import anndata as ad

    stem = file_prefix or "*"
    for view in ("ae", "pe"):
        candidates = list(ds_path.glob(f"{stem}.{view}.h5ad"))
        if candidates:
            return ad.read_h5ad(candidates[0])
        candidates = list(ds_path.glob(f"{stem}.{view}.zarr"))
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
        if isinstance(de, dict):
            de = _long_de_results(de)
        if not de.empty:
            return _reshape_de_results(de)

    # Already in matrix form or no de_results — return as-is
    return adata


def _long_de_results(contrasts: dict) -> pd.DataFrame:
    """Flatten the specification's ``{contrast: {names, logfoldchanges, ...}}`` layout into long form."""
    frames = []
    for contrast, fields in contrasts.items():
        frame = pd.DataFrame(
            {"contrast": contrast, "protein": list(fields["names"]), "log2FC": fields.get("logfoldchanges", np.nan)}
        )
        for column in ("pvals_adj", "scores", "pvals", "se"):
            if column in fields:
                frame[column] = fields[column]
        frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _reshape_de_results(de: pd.DataFrame):
    """Reshape long-form DE results into an AnnData matrix.

    Expects columns: contrast, protein, log2FC and optionally pvals_adj, scores,
    pvals, se. A protein a contrast did not test is NaN, because 0 would read as
    no change and a p-value of 0.
    """
    import anndata as ad

    # Identify contrast and protein columns
    contrast_col = "contrast" if "contrast" in de.columns else de.columns[0]
    protein_col = "protein" if "protein" in de.columns else de.columns[1]

    if "log2FC" not in de.columns:
        raise ValueError(f"DE results missing required 'log2FC' column. Found columns: {list(de.columns)}")

    contrasts = pd.Index(sorted(de[contrast_col].unique()), name="contrast")
    proteins = pd.Index(sorted(de[protein_col].unique()), name="protein")

    def matrix(column: str) -> np.ndarray:
        return (
            de.pivot_table(index=contrast_col, columns=protein_col, values=column, aggfunc="first", dropna=False)
            .reindex(index=contrasts, columns=proteins)
            .to_numpy(dtype="float64")
        )

    return ad.AnnData(
        X=matrix("log2FC"),
        obs=pd.DataFrame(index=contrasts),
        var=pd.DataFrame(index=proteins),
        layers={column: matrix(column) for column in ("pvals_adj", "scores", "pvals", "se") if column in de.columns},
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


def build_mudata(
    dataset: Dataset,
    intensity_label: str | None = None,
    modalities: list[str] | None = None,
    all_intensity_labels: bool = False,
):
    """Assemble QPX dataset into a MuData object.

    Parameters
    ----------
    dataset : Dataset
        An open QPX Dataset instance.
    intensity_label : str, optional
        Intensity label to extract (e.g. "TMT126").  When omitted, label-free
        data auto-detects its single label, while multiplexed (isobaric/plexDIA)
        data expands over every channel as ``run|label`` observations.
    modalities : list of str, optional
        Which modalities to include.  Valid values:
        ``"precursors"``, ``"proteins"``, ``"expression"``, ``"differential"``.
        Defaults to all available.
    all_intensity_labels : bool
        Force every intensity label into a separate ``run|label`` observation.
        Defaults to False: label-free data then stays run-level, while
        multiplexed data still auto-expands over its channels (so isobaric/
        plexDIA never collapses to a single channel).

    Returns
    -------
    mudata.MuData
    """
    import mudata as mu

    engine = dataset._engine
    ds_path = Path(dataset.path)
    file_prefix = dataset.file_prefix

    if modalities is not None:
        unknown = set(modalities) - _VALID_MODALITIES
        if unknown:
            raise ValueError(f"Unknown modalities: {unknown}. Valid: {sorted(_VALID_MODALITIES)}")
        requested = set(modalities)
    else:
        requested = _available_quant_modalities(dataset) | {"expression", "differential"}

    feat_intensity, pg_intensity = _select_intensity_labels(engine, intensity_label, all_intensity_labels)

    mod: dict = {}
    builders = {
        "precursors": lambda: _build_precursor_adata(engine, feat_intensity),
        "proteins": lambda: _build_protein_adata(engine, pg_intensity),
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

    # Attach feature mapping if both precursors and proteins are present
    if "precursors" in mod and "proteins" in mod:
        try:
            mapping = _build_feature_mapping(engine, mdata)
            mdata.varp["feature_mapping"] = mapping
        except Exception as exc:
            logger.warning("Failed to build feature mapping: %s", exc)

    # Attach dataset-level metadata
    _attach_uns_metadata(engine, mdata)

    # Sanitize NaN in object columns across all modalities to prevent
    # h5py serialization errors (can't write NaN in vlen-string arrays).
    for adata in mdata.mod.values():
        for frame in (adata.obs, adata.var):
            str_cols = frame.select_dtypes(include=["object"]).columns
            if len(str_cols):
                frame[str_cols] = frame[str_cols].fillna("")

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
