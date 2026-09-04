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
from qpx.core.files import run_file_stem
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

    encoded = pc.dictionary_encode(column.combine_chunks())
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
    row_pos = _dictionary_positions(table.column(row_col), row_index)
    col_pos = _dictionary_positions(table.column(col_col), col_index)

    mask = (row_pos >= 0) & (col_pos >= 0)
    values = table.column(value_col).to_numpy(zero_copy_only=False)
    data = values[mask].astype(np.float32)

    mat = sparse.coo_matrix(
        (data, (row_pos[mask], col_pos[mask])),
        shape=(len(row_index), len(col_index)),
    )
    return mat.tocsr()


def _sorted_unique(table, column: str) -> pd.Index:
    """Sorted distinct values of an Arrow column, named after the column."""
    import pyarrow.compute as pc

    values = pc.dictionary_encode(table.column(column).combine_chunks()).dictionary.to_pylist()
    return pd.Index(sorted(v for v in values if v is not None), name=column)


# pg is flattened since 1.1 (scalar ``label``), so its label queries read the
# column directly; feature keeps the intensities list<struct>.
_LABEL_FIELD_QUERIES: dict[tuple[str, str], str] = {
    ("channel", "feature"): "SELECT i.channel FROM feature, UNNEST(intensities) AS _t(i) LIMIT 1",
    ("channel", "pg"): "SELECT label FROM pg WHERE label IS NOT NULL LIMIT 1",
    ("label", "feature"): "SELECT i.label FROM feature, UNNEST(intensities) AS _t(i) LIMIT 1",
    ("label", "pg"): "SELECT label FROM pg WHERE label IS NOT NULL LIMIT 1",
}

_MULTIPLEXED_QUERIES: dict[tuple[str, str], str] = {
    ("channel", "feature"): "SELECT COUNT(DISTINCT i.channel) FROM feature, UNNEST(intensities) AS _t(i)",
    ("channel", "pg"): "SELECT COUNT(DISTINCT label) FROM pg WHERE label IS NOT NULL",
    ("label", "feature"): "SELECT COUNT(DISTINCT i.label) FROM feature, UNNEST(intensities) AS _t(i)",
    ("label", "pg"): "SELECT COUNT(DISTINCT label) FROM pg WHERE label IS NOT NULL",
}


_PRECURSOR_QUERIES: dict[str, str] = {
    "channel": """
    SELECT f.run_file_name,
           f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
           i.intensity
    FROM feature f, UNNEST(f.intensities) AS _t(i)
    WHERE i.channel = $1
    """,
    "label": """
    SELECT f.run_file_name,
           f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
           i.intensity
    FROM feature f, UNNEST(f.intensities) AS _t(i)
    WHERE i.label = $1
    """,
}

_PRECURSOR_ALL_QUERIES: dict[str, str] = {
    "channel": """
    SELECT f.run_file_name,
           CAST(i.channel AS VARCHAR) AS intensity_label,
           f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
           i.intensity
    FROM feature f, UNNEST(f.intensities) AS _t(i)
    WHERE i.channel IS NOT NULL
    """,
    "label": """
    SELECT f.run_file_name,
           CAST(i.label AS VARCHAR) AS intensity_label,
           f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
           i.intensity
    FROM feature f, UNNEST(f.intensities) AS _t(i)
    WHERE i.label IS NOT NULL
    """,
}

# A protein group's quantification is reported once for a group of runs. Explode
# grouped_runs so EVERY member run receives the group's intensity, rather than
# only the first element (the old ``grouped_runs[1]``, which silently dropped all
# but one run of a genuine multi-run group).
# Since QPX 1.1 the pg view is flattened: scalar ``label`` / ``intensity`` columns,
# one row per label, so there is no intensities list to UNNEST. grouped_runs is
# exploded through ``list_distinct`` so a member run counts once even if a producer
# repeated it. ID-only rows (null label) are excluded by the label predicate. The
# ``channel`` and ``label`` keys are identical now (pg carries only ``label``); both
# are kept so the shared label-field dispatch does not need a pg special-case.
_PROTEIN_QUERIES: dict[str, str] = {
    "channel": """
    SELECT gr AS run_file_name,
           pg.anchor_protein,
           pg.intensity AS intensity
    FROM pg,
         UNNEST(list_distinct(pg.grouped_runs)) AS _g(gr)
    WHERE pg.label = $1
    """,
    "label": """
    SELECT gr AS run_file_name,
           pg.anchor_protein,
           pg.intensity AS intensity
    FROM pg,
         UNNEST(list_distinct(pg.grouped_runs)) AS _g(gr)
    WHERE pg.label = $1
    """,
}

_PROTEIN_ALL_QUERIES: dict[str, str] = {
    "channel": """
    SELECT gr AS run_file_name,
           CAST(pg.label AS VARCHAR) AS intensity_label,
           pg.anchor_protein,
           pg.intensity AS intensity
    FROM pg,
         UNNEST(list_distinct(pg.grouped_runs)) AS _g(gr)
    WHERE pg.label IS NOT NULL
    """,
    "label": """
    SELECT gr AS run_file_name,
           CAST(pg.label AS VARCHAR) AS intensity_label,
           pg.anchor_protein,
           pg.intensity AS intensity
    FROM pg,
         UNNEST(list_distinct(pg.grouped_runs)) AS _g(gr)
    WHERE pg.label IS NOT NULL
    """,
}


_TYPEOF_QUERIES: dict[str, str] = {
    "feature": "SELECT typeof(intensities) FROM feature LIMIT 1",
    "pg": "SELECT typeof(label) FROM pg LIMIT 1",
}


def _detect_label_field(engine: DuckDBEngine, table: str = "feature") -> str:
    """Return ``"channel"`` or ``"label"`` by inspecting the intensities struct type."""
    try:
        row = engine.execute(_TYPEOF_QUERIES[table]).fetchone()
        if row:
            type_str = row[0].lower()
            if "channel" in type_str and "label" not in type_str:
                return "channel"
        return "label"
    except Exception:
        logger.debug("Could not detect label field from %s; defaulting to 'label'", table)
        return "label"


def _detect_intensity_label(engine: DuckDBEngine, table: str = "feature") -> str:
    """Auto-detect the first intensity label from the given table."""
    label_field = _detect_label_field(engine, table)

    query = _LABEL_FIELD_QUERIES.get((label_field, table))
    if query is None:
        raise ValueError(f"Invalid label_field={label_field!r} or table={table!r}")
    try:
        row = engine.execute(query).fetchone()
    except Exception as exc:
        raise ValueError(f"Failed to read intensity labels from {table}: {exc}") from exc
    if row is None:
        raise ValueError(f"No intensity labels found in {table}")
    return row[0]


def _prepare_observations(
    df: pd.DataFrame,
    intensity_label: str | None,
) -> tuple[pd.DataFrame, str, pd.Index, pd.DataFrame]:
    """Create matrix row IDs and their run/label lookup keys."""
    if intensity_label is None:
        df = df.copy()
        df["observation_id"] = df["run_file_name"].astype(str) + "|" + df["intensity_label"].astype(str)
        obs_index = pd.Index(sorted(df["observation_id"].unique()), name="observation_id")
        keys = (
            df[["observation_id", "run_file_name", "intensity_label"]]
            .drop_duplicates(subset="observation_id", keep="first")
            .set_index("observation_id")
            .reindex(obs_index)
        )
        return df, "observation_id", obs_index, keys

    obs_index = pd.Index(sorted(df["run_file_name"].unique()), name="run_file_name")
    keys = pd.DataFrame(
        {
            "run_file_name": obs_index,
            "intensity_label": str(intensity_label),
        },
        index=obs_index,
    )
    return df, "run_file_name", obs_index, keys


def _prepare_observations_arrow(
    table,
    intensity_label: str | None,
):
    """Arrow-native :func:`_prepare_observations`.

    Returns ``(table, row_col, obs_index, observation_keys)``. When every label is
    requested the ``observation_id`` column is appended to the Arrow table rather
    than built row-by-row in pandas: the old path copied the whole frame and then
    concatenated two object columns for every row, which on a large experiment is
    hundreds of millions of throwaway Python strings.
    """
    import pyarrow as pa
    import pyarrow.compute as pc

    if intensity_label is None:
        runs = table.column("run_file_name").combine_chunks()
        labels = table.column("intensity_label").combine_chunks()
        observation_id = pc.binary_join_element_wise(
            pc.cast(runs, pa.string()),
            pc.cast(labels, pa.string()),
            "|",
        )
        table = table.append_column("observation_id", observation_id)
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


def _load_run_obs(
    engine: DuckDBEngine,
    observation_keys: pd.DataFrame,
    all_labels: bool,
) -> pd.DataFrame:
    """Load run/sample metadata for matrix observations."""
    try:
        df = engine.execute(
            """
            SELECT r.run_file_name,
                   r.run_accession,
                   rs.sample_accession,
                   rs.label,
                   rs.biological_replicate,
                   rs.technical_replicate,
                   r.fraction,
                   r.instrument
            FROM run r, UNNEST(r.samples) AS _t(rs)
            """
        ).fetchdf()
    except Exception:
        if all_labels:
            return observation_keys.copy()
        return pd.DataFrame(index=observation_keys.index)

    if df.empty:
        if all_labels:
            return observation_keys.copy()
        return pd.DataFrame(index=observation_keys.index)

    exact: dict[tuple[str, str], dict] = {}
    stem_exact: dict[tuple[str, str], dict] = {}
    first_by_run: dict[str, dict] = {}
    first_by_stem: dict[str, dict] = {}
    labels_by_run: dict[str, set[str]] = {}
    labels_by_stem: dict[str, set[str]] = {}
    for record in df.to_dict("records"):
        run_name = str(record["run_file_name"])
        stem = run_file_stem(run_name)
        label = str(record["label"])
        exact.setdefault((run_name, label), record)
        stem_exact.setdefault((stem, label), record)
        first_by_run.setdefault(run_name, record)
        first_by_stem.setdefault(stem, record)
        labels_by_run.setdefault(run_name, set()).add(label)
        labels_by_stem.setdefault(stem, set()).add(label)

    rows: list[dict] = []
    for _, key in observation_keys.iterrows():
        run_name = str(key["run_file_name"])
        stem = run_file_stem(run_name)
        label = str(key["intensity_label"])
        record = exact.get((run_name, label))
        if record is None:
            record = stem_exact.get((stem, label))
        # Fall back to the run's sole sample when the intensity label does not
        # match any SDRF sample label. This is the label-free case (intensity
        # label "LFQ" vs sample label "label free sample"), detected by the run
        # carrying at most one distinct SDRF label. For multiplexed runs (more
        # than one distinct label / channel) we never guess, so each channel
        # keeps its own real sample_accession.
        if record is None and len(labels_by_run.get(run_name, ())) <= 1:
            record = first_by_run.get(run_name)
        if record is None and len(labels_by_stem.get(stem, ())) <= 1:
            record = first_by_stem.get(stem)

        values = dict(record) if record is not None else {}
        values.pop("run_file_name", None)
        if all_labels:
            values = {
                "run_file_name": run_name,
                "intensity_label": label,
                **values,
            }
        rows.append(values)

    obs = pd.DataFrame(rows, index=observation_keys.index)

    # Fill remaining NaN in string columns to avoid h5py serialization errors
    str_cols = obs.select_dtypes(include=["object"]).columns
    obs[str_cols] = obs[str_cols].fillna("")
    return obs


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


def _build_precursor_adata(
    engine: DuckDBEngine,
    intensity_label: str | None,
    label_field: str = "label",
):
    """Build AnnData from feature.parquet.

    obs = runs (one label) or run/label pairs (all labels),
    var = peptidoform|charge, X = sparse intensity matrix.
    """
    import anndata as ad

    if intensity_label is None:
        table = engine.execute(_PRECURSOR_ALL_QUERIES[label_field]).fetch_arrow_table()
    else:
        table = engine.execute(_PRECURSOR_QUERIES[label_field], [intensity_label]).fetch_arrow_table()

    if table.num_rows == 0:
        return ad.AnnData()

    table, row_col, observations, observation_keys = _prepare_observations_arrow(table, intensity_label)
    precursors = _sorted_unique(table, "precursor_id")

    intensity_matrix = _pivot_arrow_to_sparse(table, row_col, "precursor_id", "intensity", observations, precursors)

    obs = _load_run_obs(engine, observation_keys, all_labels=intensity_label is None)
    var = pd.DataFrame(index=precursors)

    return ad.AnnData(X=intensity_matrix, obs=obs, var=var)


def _build_protein_adata(
    engine: DuckDBEngine,
    intensity_label: str | None,
    label_field: str = "label",
):
    """Build AnnData from pg.parquet.

    obs = runs (one label) or run/label pairs (all labels),
    var = anchor_protein, X = sparse intensity matrix.
    var includes gene_name column.
    """
    import anndata as ad

    if intensity_label is None:
        table = engine.execute(_PROTEIN_ALL_QUERIES[label_field]).fetch_arrow_table()
    else:
        table = engine.execute(_PROTEIN_QUERIES[label_field], [intensity_label]).fetch_arrow_table()

    if table.num_rows == 0:
        return ad.AnnData()

    table, row_col, observations, observation_keys = _prepare_observations_arrow(table, intensity_label)
    proteins = _sorted_unique(table, "anchor_protein")

    intensity_matrix = _pivot_arrow_to_sparse(table, row_col, "anchor_protein", "intensity", observations, proteins)

    obs = _load_run_obs(engine, observation_keys, all_labels=intensity_label is None)

    # Gene names
    gene_sql = """
    SELECT DISTINCT pg.anchor_protein,
           pg.gg_names[1] AS gene_name
    FROM pg
    """
    gene_df = engine.execute(gene_sql).fetchdf()
    gene_df = gene_df.drop_duplicates(subset="anchor_protein", keep="first")
    gene_df = gene_df.set_index("anchor_protein").reindex(proteins)

    var = pd.DataFrame({"gene_name": gene_df["gene_name"].values}, index=proteins)

    return ad.AnnData(X=intensity_matrix, obs=obs, var=var)


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


def _try_build_modality(name: str, builder, mod: dict) -> None:
    """Call *builder*, add the result to *mod* if non-empty."""
    try:
        adata = builder()
        if adata is not None and (not hasattr(adata, "n_obs") or adata.n_obs > 0):
            mod[name] = adata
    except Exception as exc:
        logger.warning("Failed to build %s modality: %s", name, exc)


def _is_multiplexed(engine: DuckDBEngine, table: str, label_field: str) -> bool:
    """Return True when *table* carries more than one distinct intensity label.

    Multiplexed data (isobaric TMT/iTRAQ, plexDIA) has many channel labels per
    run; label-free has a single label ("LFQ"). Detection drives whether obs is
    expanded over run x channel or kept at the run level.
    """
    query = _MULTIPLEXED_QUERIES.get((label_field, table))
    if query is None:
        return False
    try:
        row = engine.execute(query).fetchone()
    except Exception:
        return False
    return bool(row) and row[0] is not None and row[0] > 1


def _select_intensity_labels(
    engine: DuckDBEngine,
    intensity_label: str | None,
    all_intensity_labels: bool,
    feat_label_field: str,
    pg_label_field: str,
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

    if _is_multiplexed(engine, "feature", feat_label_field):
        feature_label = None
    else:
        try:
            feature_label = _detect_intensity_label(engine, "feature")
        except ValueError:
            # No detectable feature label (e.g. a proteins-only build): fall back
            # to expansion rather than aborting the whole MuData assembly.
            feature_label = None
    if _is_multiplexed(engine, "pg", pg_label_field):
        protein_label = None
    else:
        try:
            protein_label = _detect_intensity_label(engine, "pg")
        except ValueError:
            protein_label = feature_label
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
    file_prefix = getattr(dataset, "_file_prefix", None)

    if modalities is not None:
        unknown = set(modalities) - _VALID_MODALITIES
        if unknown:
            raise ValueError(f"Unknown modalities: {unknown}. Valid: {sorted(_VALID_MODALITIES)}")
        requested = set(modalities)
    else:
        requested = _VALID_MODALITIES.copy()

    feat_label_field = _detect_label_field(engine, "feature")
    pg_label_field = _detect_label_field(engine, "pg")
    feat_intensity, pg_intensity = _select_intensity_labels(
        engine, intensity_label, all_intensity_labels, feat_label_field, pg_label_field
    )

    mod: dict = {}
    builders = {
        "precursors": lambda: _build_precursor_adata(engine, feat_intensity, feat_label_field),
        "proteins": lambda: _build_protein_adata(engine, pg_intensity, pg_label_field),
        "expression": lambda: _build_expression_adata(ds_path, file_prefix),
        "differential": lambda: _build_differential_adata(ds_path, file_prefix),
    }
    for name, builder in builders.items():
        if name in requested:
            _try_build_modality(name, builder, mod)

    mdata = mu.MuData(mod)

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
