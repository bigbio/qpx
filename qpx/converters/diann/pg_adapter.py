"""DIA-NN PG adapter -- report.tsv (+ optional pg_matrix.tsv) to pg.parquet.

Loads a DIA-NN ``report.tsv`` and optionally the ``pg_matrix.tsv`` into
DuckDB, aggregates protein-group quantities and scores per group, run and
channel, and writes through ``PgWriter``.

Key schema changes:
    - ``reference_file_name`` -> ``run_file_name``
    - ``is_decoy`` (int, 0/1) -> ``is_decoy`` (bool)
    - Intensities struct: ``{sample_accession, channel, intensity}`` ->
      ``{label, intensity}``
"""

from __future__ import annotations

import logging
import math
import re
from typing import Optional

import numpy as np
import pandas as pd

from qpx.converters.base import resolve_columns
from qpx.converters.diann.base_adapter import (
    BATCH_ROWS,
    DiaNNBaseAdapter,
    decoy_protein_group_sql,
    run_name_sql,
)
from qpx.converters.mappings import get_field_mappings
from qpx.converters.utils import is_contaminant_accession
from qpx.core.files import RUN_FILE_SUFFIX_PATTERN
from qpx.core.sql import sql_build, validate_identifier
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)

_BATCH_SQL = """
            SELECT
                $select_clause
            FROM report
            $channel_join
            WHERE $run_col IN ($placeholders)
              AND $pg_col IS NOT NULL
              $qvalue_filter
            """

_EXT_RE = re.compile(RUN_FILE_SUFFIX_PATTERN)

# Extra columns needed for PG aggregation but not in the field mappings
_PG_EXTRA_COLS = [
    ('"Proteotypic"', "proteotypic"),
    ('"Stripped.Sequence"', "stripped_sequence"),
    ('"Precursor.Id"', "precursor_id"),
    # Raw protein-group quantity (label-free primary intensity). Kept distinct
    # from PG.MaxLFQ (additional_intensities) so the primary value is not the
    # MaxLFQ matrix number duplicated under the "LFQ" label.
    ('"PG.Quantity"', "pg_quantity_raw"),
    ('"PG.PEP"', "pg_pep"),
    ('"PG.Normalised"', "pg_normalised"),
    ('"PG.TopN"', "pg_topn"),
    ('"Genes.MaxLFQ"', "genes_maxlfq"),
    ('"Genes.MaxLFQ.Unique"', "genes_maxlfq_unique"),
    ('"Genes.TopN"', "genes_topn"),
    ('"PG.MaxLFQ.Quality"', "pg_maxlfq_quality"),
    ('"Genes.MaxLFQ.Quality"', "genes_maxlfq_quality"),
    ('"Lib.PG.Q.Value"', "lib_pg_qvalue"),
    ('"Protein.Q.Value"', "protein_qvalue"),
]

# Precursor q-value column filtered by ``qvalue_threshold``, as in the feature view.
_PRECURSOR_QVALUE_COL = "Q.Value"

# Present when DIA-NN ran with --report-decoys (1 = decoy precursor).
_DECOY_COL = "Decoy"

# DIA-NN decoy protein-accession prefixes (mirrors feature_adapter is_decoy logic).
_DECOY_PREFIXES = ("DECOY_", "decoy_", "rev_", "REV_")

# Proteotypic flags as they appear after ``astype(str)`` for int, float and bool columns.
_PROTEOTYPIC_VALUES = ("1", "1.0", "True")

# Protein-group scores DIA-NN repeats on every precursor row of a group and
# channel, with whether higher is better. A plexDIA report gives each channel its
# own. Older reports may omit any of them.
_CHANNEL_SCORES = {
    "global_qvalue": False,
    "pg_qvalue": False,
    "gg_qvalue": False,
    "pg_pep": False,
    "lib_pg_qvalue": False,
    "protein_qvalue": False,
    "pg_maxlfq_quality": True,
    "genes_maxlfq_quality": True,
}

# Protein-group scores carried in additional_scores (pg.md).
_ADDITIONAL_SCORES = ("pg_qvalue", "pg_pep", "lib_pg_qvalue", "protein_qvalue", "pg_maxlfq_quality", "genes_maxlfq_quality")

# Additional protein-group quantities and their intensity names (pg.md). PG.MaxLFQ
# ("lfq") is the primary intensity instead when the report has no PG.Quantity.
_ADDITIONAL_QUANTITIES = (
    ("lfq", "maxlfq"),
    ("pg_normalised", "normalize_intensity"),
    ("pg_topn", "topn"),
    ("genes_maxlfq", "genes_maxlfq"),
    ("genes_maxlfq_unique", "genes_maxlfq_unique"),
    ("genes_topn", "genes_topn"),
)

# Protein-group quantities, aggregated per (group, run, channel).
_QUANTITY_COLUMNS = ("pg_quantity_raw", *(alias for alias, _ in _ADDITIONAL_QUANTITIES))

_GROUP_KEYS = ["pg_membership_key", "run_file_name"]


def _pg_members(value) -> list[str]:
    """Return the deduplicated, sorted protein-group membership."""
    if pd.isna(value):
        return []
    return sorted({part.strip() for part in str(value).split(";") if part.strip()})


def _canonical_pg_membership(value) -> str:
    """Return one order-independent key for a DIA-NN protein group."""
    return ";".join(_pg_members(value))


def _numeric(series: pd.Series) -> pd.Series:
    """Vectorised ``safe_float``: float64 values, NaN where missing or not numeric."""
    return pd.to_numeric(series, errors="coerce").astype("float64")


def _measured(series: pd.Series) -> pd.Series:
    """Quantities as float64, NaN where missing or not positive: DIA-NN writes 0 for a quantity it has not estimated."""
    values = _numeric(series)
    return values.where(values > 0)


def _optional(value: float) -> float | None:
    """Map a NaN from a float64 column to None."""
    return None if math.isnan(value) else value


def _split_annotation(value) -> list[str] | None:
    """Split a ``;``-separated name/gene annotation, None when absent."""
    return str(value).split(";") if pd.notna(value) and value else None


def _aggregate_groups(frame: pd.DataFrame) -> pd.DataFrame:
    """Aggregate one value per (protein group, run) for the fields shared by its channels.

    The leader is taken from the smallest producer representation so it does not
    depend on row order. Names and genes take the first non-empty value, so
    inconsistent per-precursor annotations do not split a group into
    duplicate-identity records. Counts are distinct stripped sequences and
    precursors, overall and among proteotypic rows.
    """
    proteotypic = frame["proteotypic"].astype(str).isin(_PROTEOTYPIC_VALUES)
    feature_col = "precursor_id" if "precursor_id" in frame.columns else "stripped_sequence"
    # Rank the distinct producer representations once, so the smallest per group
    # is an integer minimum instead of a string comparison per group.
    leader_codes, leader_sources = pd.factorize(frame["pg_accessions"], sort=True)
    groups = (
        frame.assign(
            leader_code=leader_codes,
            pg_names=frame["pg_names"].where(frame["pg_names"] != ""),
            gg_accessions=frame["gg_accessions"].where(frame["gg_accessions"] != ""),
            unique_sequence=frame["stripped_sequence"].where(proteotypic),
            unique_feature=frame[feature_col].where(proteotypic),
        )
        .groupby(_GROUP_KEYS, sort=True, dropna=False)
        .agg(
            leader_code=("leader_code", "min"),
            pg_names=("pg_names", "first"),
            gg_accessions=("gg_accessions", "first"),
            total_sequences=("stripped_sequence", "nunique"),
            unique_sequences=("unique_sequence", "nunique"),
            total_features=(feature_col, "nunique"),
            unique_features=("unique_feature", "nunique"),
        )
    )
    groups["leader_source"] = np.asarray(leader_sources)[groups.pop("leader_code").to_numpy()]
    return groups


def _matrix_fallback(groups: pd.DataFrame, pg_matrix_indexed: pd.DataFrame) -> list[float | None]:
    """Return each group's PG matrix quantity, None when missing, NaN or zero."""
    runs = groups.index.get_level_values("run_file_name")
    matrix_runs = [run for run in runs.unique() if run in pg_matrix_indexed.columns]
    values = np.full(len(groups), np.nan)
    if matrix_runs:
        rows = pg_matrix_indexed.index.get_indexer(groups.index.get_level_values("pg_membership_key"))
        cols = pd.Index(matrix_runs).get_indexer(runs)
        found = (rows >= 0) & (cols >= 0)
        values[found] = pg_matrix_indexed[matrix_runs].to_numpy(dtype="float64")[rows[found], cols[found]]
    return [None if math.isnan(value) or value <= 0 else value for value in values.tolist()]


def _channel_quantities(label: str, values: dict[str, float | None]) -> tuple[dict, dict | None]:
    """Build one channel's intensity and additional-intensity entries.

    Raw PG.Quantity is preferred for the primary value. Reports that omit it use
    MaxLFQ as a compatibility fallback; the experimental channel remains the
    label, while "maxlfq" is recorded as the algorithm name in
    additional_intensities. A channel with no quantity keeps its label with a
    null intensity, since the group was identified in the run.
    """
    raw_quantity, maxlfq = values["pg_quantity_raw"], values["lfq"]
    primary = raw_quantity if raw_quantity is not None else maxlfq
    intensity = {"label": label, "intensity": primary, "cv_params": None}
    if primary is not None:
        method = "PG.Quantity" if raw_quantity is not None else "PG.MaxLFQ"
        intensity["cv_params"] = [{"cv_name": "quantification_method", "cv_value": method}]
    entries = [
        {"intensity_name": name, "intensity_value": values[alias]}
        for alias, name in _ADDITIONAL_QUANTITIES
        if values[alias] is not None and not (alias == "lfq" and raw_quantity is None)
    ]
    return intensity, {"label": label, "intensities": entries} if entries else None


def _membership_fields(key: str) -> tuple[list[str], bool, bool]:
    """Return the members, decoy flag and contaminant flag of a canonical membership key.

    A PG made up ENTIRELY of decoy proteins is a decoy group. DIA-NN reports a
    decoy precursor under its target's protein group and flags it only in the
    ``Decoy`` column, so ``_process_batch`` selects those accessions with a
    ``DECOY_`` prefix; the prefix also covers decoy entries of a target-decoy FASTA.
    """
    members = _pg_members(key)
    is_decoy = all(acc.startswith(_DECOY_PREFIXES) for acc in members)
    return members, is_decoy, any(is_contaminant_accession(acc) for acc in members)


def _pg_record(
    row,
    channel,
    membership: tuple[list[str], bool, bool],
    quantities: tuple[dict, dict | None],
    additional_scores: list[dict] | None,
) -> dict:
    """Build one natural pg record from an aggregated (protein group, run) row and one of its channels."""
    pg_accessions, is_decoy, contaminant = membership
    intensity, additional = quantities
    gg_accessions = _split_annotation(row.gg_accessions)
    peptide_count = max(row.total_sequences, 1)
    return {
        "pg_accessions": pg_accessions,
        "pg_names": _split_annotation(row.pg_names),
        "gg_accessions": gg_accessions,
        "gg_names": gg_accessions,  # Gene symbols serve as both accession and name
        "gg_qvalue": _optional(channel.gg_qvalue),
        "anchor_protein": row.leader_source.split(";", 1)[0].strip() or pg_accessions[0],
        "grouped_runs": [row.run_file_name],
        "global_qvalue": _optional(channel.global_qvalue),
        "pg_qvalue": _optional(channel.pg_qvalue),
        "intensities": [intensity],
        "additional_intensities": [additional] if additional else None,
        "is_decoy": is_decoy,
        "contaminant": contaminant,
        "peptides": [{"protein_name": acc, "peptide_count": peptide_count} for acc in pg_accessions],
        "peptide_counts": {
            "unique_sequences": row.unique_sequences,
            "total_sequences": row.total_sequences,
        },
        "feature_counts": {
            "unique_features": row.unique_features,
            "total_features": row.total_features,
        },
        "sequence_coverage": None,
        "molecular_weight": None,
        "additional_scores": additional_scores,
        "cv_params": None,
    }


class DiannPgAdapter(DiaNNBaseAdapter):
    """Convert DIA-NN report + PG matrix to ``pg.parquet``.

    Usage::

        with DiannPgAdapter() as adapter:
            adapter.convert(
                diann_report="report.tsv",
                pg_matrix_path="pg_matrix.tsv",
                output_path="pg.parquet",
                sdrf_path="sdrf.tsv",
            )
    """

    def convert(
        self,
        diann_report: str,
        pg_matrix_path: Optional[str],
        output_path: str,
        sdrf_path: Optional[str] = None,
        batch_rows: int = BATCH_ROWS,
        creator: str = "diann",
        qvalue_threshold: Optional[float] = None,
    ) -> None:
        """Run the DIA-NN report+matrix -> pg.parquet conversion.

        Args:
            diann_report: Path to the DIA-NN ``report.tsv``.
            pg_matrix_path: Optional DIA-NN ``pg_matrix.tsv``. A label-free group
                without PG.MaxLFQ in the report takes the matrix value; the
                report alone gives the pg view, and DIA-NN writes no matrix for
                a channel (plexDIA) report.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF file for sample mapping.
            batch_rows: Report rows to aggregate per batch of runs.
            creator: Creator tag in Parquet metadata.
            qvalue_threshold: Optional precursor ``Q.Value`` cutoff. Filtering is
                opt-in (bigbio/qpx#241): with the default ``None`` every protein
                group in the report is emitted. When a threshold is provided,
                only precursor rows with ``Q.Value`` below it are aggregated, the
                rows the feature view keeps. PG-level q-values and ``PG.PEP`` are
                carried for every group and never filtered (bigbio/qpx#342).
        """
        # Step 1: Load report into DuckDB
        self._load_diann_report(diann_report)

        # Step 2: Load PG matrix and pre-index for O(1) lookups
        pg_matrix_indexed = self._load_pg_matrix(pg_matrix_path).set_index("pg_accessions") if pg_matrix_path else None

        # Step 3: Cache report columns and resolve mappings (once, not per batch)
        report_cols = {
            c[0]
            for c in self._conn.execute("SELECT column_name FROM information_schema.columns WHERE table_name='report'").fetchall()
        }
        self._resolved_pg = resolve_columns(get_field_mappings("diann", "pg"), report_cols)
        channel_col = self._register_channel_labels(report_cols, sdrf_path)

        # Step 4: Get unique runs and process in batches
        runs = self._get_unique_runs()

        quantified = 0
        with PgWriter(output_path, creator=creator, compression=self._compression) as writer:
            done = 0
            for batch_runs in self._run_batches(validate_identifier(self._resolved_pg["run_file_name"]), runs, batch_rows):
                self.logger.info(f"Processing PG runs {done + 1}-{done + len(batch_runs)} of {len(runs)}")
                done += len(batch_runs)
                records = self._process_batch(
                    batch_runs,
                    pg_matrix_indexed,
                    report_cols,
                    channel_col,
                    qvalue_threshold=qvalue_threshold,
                )
                if records:
                    quantified += sum(record["intensities"][0]["intensity"] is not None for record in records)
                    writer.write_batch(records)

        if runs and not quantified:
            self.logger.warning(
                "No DIA-NN protein group has a quantity: the report carries no positive PG.Quantity or PG.MaxLFQ, "
                "so every pg intensity is null"
            )
        self.logger.info(f"DIA-NN PG conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_pg_matrix(self, path: str) -> pd.DataFrame:
        """Load the DIA-NN PG matrix TSV."""
        # Resolve against TSV header columns (may differ from report)
        header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        header_set = set(header)
        pg_matrix_resolved = resolve_columns(get_field_mappings("diann", "pg"), header_set)
        pg_col = pg_matrix_resolved["pg_accessions"]
        names_col = pg_matrix_resolved["pg_names"]
        genes_col = pg_matrix_resolved["gg_accessions"]

        meta_cols = {pg_col, names_col, genes_col, "Protein.Ids", "First.Protein.Description"}
        run_cols = [c for c in header if c not in meta_cols]
        usecols = [pg_col, names_col, genes_col] + run_cols
        df = pd.read_csv(path, sep="\t", usecols=usecols)
        df.rename(
            columns={
                pg_col: "pg_accessions",
                names_col: "pg_names",
                genes_col: "gg_accessions",
            },
            inplace=True,
        )
        # Strip common file extensions from run column names
        df.columns = [_EXT_RE.sub("", c) for c in df.columns]
        df["pg_accessions"] = df["pg_accessions"].map(_canonical_pg_membership)
        duplicate_keys = sorted(df.loc[df["pg_accessions"].duplicated(keep=False), "pg_accessions"].unique())
        if duplicate_keys:
            raise ValueError(
                "DIA-NN PG matrix contains duplicate canonical protein-group memberships: " + ", ".join(duplicate_keys[:10])
            )
        return df

    def _get_unique_runs(self) -> list[str]:
        """Get sorted list of unique Run values from the report."""
        run_col = self._resolved_pg["run_file_name"]
        qcol = validate_identifier(run_col)
        rows = self._conn.execute(
            sql_build(
                "SELECT DISTINCT $col FROM report ORDER BY $col",
                col=qcol,
            )
        ).fetchall()
        return [r[0] for r in rows]

    # ------------------------------------------------------------------
    # Batch processing
    # ------------------------------------------------------------------

    def _batch_select_parts(self, actual_report_cols: set[str], channel_col: str | None) -> list[str]:
        """Build the source and aggregation columns for a PG batch."""
        # Build SQL SELECT clause from resolved field mappings
        r = self._resolved_pg

        # Decoy rows form their own group with DECOY_ accessions.
        pg_group = f'"{r["pg_accessions"]}"'
        if _DECOY_COL in actual_report_cols:
            pg_group = decoy_protein_group_sql(pg_group, f'"{_DECOY_COL}"')
        select_parts = [f"{pg_group} AS pg_accessions"]
        select_parts += [f'"{col}" AS {qpx_field}' for qpx_field, col in r.items() if qpx_field != "pg_accessions"]
        # Add extra columns needed for aggregation (guarded by presence)
        for src_col, alias in _PG_EXTRA_COLS:
            raw_name = src_col.strip('"')
            if raw_name in actual_report_cols:
                select_parts.append(f'"{raw_name}" AS {alias}')
        if channel_col:
            select_parts.append("cl.canonical_label AS channel_label")

        # Push run_file_name extension stripping into SQL
        run_col = r["run_file_name"]
        select_parts.append(f"{run_name_sql(validate_identifier(run_col))} AS run_file_name_clean")

        return select_parts

    def _build_batch_sql(
        self,
        runs: list[str],
        actual_report_cols: set[str],
        channel_col: str | None,
        qvalue_threshold: Optional[float],
    ) -> str:
        """Build a PG batch query with run, channel, and q-value constraints."""
        r = self._resolved_pg
        select_parts = self._batch_select_parts(actual_report_cols, channel_col)

        select_clause = ",\n                ".join(select_parts)

        # Use resolved column names for filtering
        pg_col = r["pg_accessions"]

        # The q-value filter is opt-in (bigbio/qpx#241) and precursor-level, so
        # pg aggregates exactly the rows the feature view keeps. PG-level
        # q-values are carried for every group, never used as a cutoff: a
        # protein FDR is left to downstream tools (bigbio/qpx#342).
        qvalue_filter = ""
        if qvalue_threshold is not None:
            if _PRECURSOR_QVALUE_COL not in actual_report_cols:
                raise ValueError(f"qvalue_threshold requires the DIA-NN report column {_PRECURSOR_QVALUE_COL}")
            qvalue_filter = sql_build(
                "AND report.$qv_col < $threshold",
                qv_col=validate_identifier(_PRECURSOR_QVALUE_COL),
                threshold=str(float(qvalue_threshold)),
            )

        placeholders = ", ".join(["?" for _ in runs])
        channel_join = ""
        if channel_col:
            channel_join = sql_build(
                """
                JOIN diann_channel_labels cl
                  ON CAST(report.$channel_col AS VARCHAR) = cl.channel_value
                """,
                channel_col=validate_identifier(channel_col),
            )
        stmt = sql_build(
            _BATCH_SQL,
            select_clause=select_clause,
            channel_join=channel_join,
            run_col=validate_identifier(r["run_file_name"]),
            placeholders=placeholders,
            pg_col=validate_identifier(pg_col),
            qvalue_filter=qvalue_filter,
        )
        return stmt

    def _process_batch(
        self,
        runs: list[str],
        pg_matrix_indexed: pd.DataFrame | None,
        actual_report_cols: set[str],
        channel_col: str | None = None,
        qvalue_threshold: Optional[float] = None,
    ) -> list[dict]:
        """Process a batch of runs for PG quantification."""
        stmt = self._build_batch_sql(runs, actual_report_cols, channel_col, qvalue_threshold)
        report_df = self._conn.execute(stmt, runs).df()

        if report_df.empty:
            return []

        # Use SQL-computed clean run names
        report_df["run_file_name"] = report_df["run_file_name_clean"]
        report_df.drop(columns=["run_file_name_clean"], inplace=True)
        # Canonicalise each distinct producer representation once, not once per row.
        raw_groups = report_df["pg_accessions"]
        report_df["pg_membership_key"] = raw_groups.map({raw: _canonical_pg_membership(raw) for raw in raw_groups.unique()})
        report_df = report_df[report_df["pg_membership_key"] != ""]

        # Aggregate per protein group per run — one record per (pg_accessions, run).
        # Grouping additionally on pg_names/gg_accessions would split a single
        # protein group into multiple records when DIA-NN emits inconsistent
        # name/gene annotations across a group's precursor rows; those records
        # share the same (pg_accessions, grouped_runs, label) identity and would
        # collide on pg_id. Names/genes are instead derived within the group
        # (first non-null), so annotation noise no longer splits a group.
        groups = _aggregate_groups(report_df)
        channels = self._aggregate_channels(report_df, groups.index)
        fallback = None
        if pg_matrix_indexed is not None and "channel_label" not in report_df.columns:
            fallback = _matrix_fallback(groups, pg_matrix_indexed)
        return self._build_pg_records(groups, channels, fallback)

    def _aggregate_channels(self, frame: pd.DataFrame, group_index: pd.MultiIndex) -> pd.DataFrame:
        """Collapse each (group, run, channel) to one quantity and one value per score.

        PG.Quantity / PG.MaxLFQ are protein-group-level values that DIA-NN repeats
        identically on every precursor row of a group within a run and channel,
        so a group's rows should agree. When annotation noise merges rows that
        disagree on the quantity, picking the "first" row makes the emitted
        intensity depend on SQL row order, which is undefined. The max of the
        measured values is deterministic, order-independent, and equal to the
        shared value in the normal (consistent) case; a warning keeps a real
        disagreement visible. A score takes its best value. Rows are indexed by
        the group's position in *group_index* and the label.
        """
        present = [column for column in _QUANTITY_COLUMNS if column in frame.columns]
        values = pd.DataFrame({column: _measured(frame[column]) for column in present}, index=frame.index)
        for column in _CHANNEL_SCORES:
            values[column] = _numeric(frame[column]) if column in frame.columns else float("nan")
        values["group"] = frame.groupby(_GROUP_KEYS, sort=True, dropna=False).ngroup()
        values["label"] = frame["channel_label"] if "channel_label" in frame.columns else "LFQ"
        grouped = values.groupby(["group", "label"], sort=True)
        channels = pd.DataFrame(index=grouped.size().index)
        conflicting = np.zeros(len(channels), dtype=bool)
        for column in _QUANTITY_COLUMNS:
            channels[column] = grouped[column].max() if column in present else float("nan")
            if column in present:
                conflicting |= (grouped[column].nunique() > 1).to_numpy()
        for column, higher_better in _CHANNEL_SCORES.items():
            channels[column] = grouped[column].max() if higher_better else grouped[column].min()
        if conflicting.any():
            self._warn_conflicting_quantities(values, channels.index[conflicting], present, group_index)
        return channels

    def _warn_conflicting_quantities(
        self,
        values: pd.DataFrame,
        conflicts: pd.MultiIndex,
        present: list[str],
        group_index: pd.MultiIndex,
    ) -> None:
        """Warn for each (group, run, channel) whose finite quantities disagree."""
        rows = values[values["group"].isin(conflicts.get_level_values("group"))]
        for (group, _label), channel_rows in rows.groupby(["group", "label"], sort=True):
            pg_acc, run_file_name = group_index[group]
            for column in present:
                finite = sorted(set(channel_rows[column].dropna().tolist()))
                if len(finite) > 1:
                    self.logger.warning(
                        "Protein group %s in run %s has conflicting %s values %s; using max (%s) deterministically.",
                        pg_acc,
                        run_file_name,
                        column,
                        finite,
                        finite[-1],
                    )

    def _build_pg_records(
        self,
        groups: pd.DataFrame,
        channels: pd.DataFrame,
        fallback: list[float | None] | None,
    ) -> list[dict]:
        """Build one natural pg record per channel of each (protein group, run) row of *groups*.

        *fallback* holds the label-free PG matrix quantity of each group, used for
        a missing MaxLFQ.
        """
        group_rows = list(groups.reset_index().itertuples(index=False))
        membership = {key: _membership_fields(key) for key in groups.index.unique("pg_membership_key")}
        records: list[dict] = []
        for channel in channels.reset_index().itertuples(index=False):
            row = group_rows[channel.group]
            values = {column: _optional(getattr(channel, column)) for column in _QUANTITY_COLUMNS}
            if values["lfq"] is None and fallback is not None:
                values["lfq"] = fallback[channel.group]
            quantities = _channel_quantities(str(channel.label), values)
            records.append(
                _pg_record(row, channel, membership[row.pg_membership_key], quantities, self._additional_scores(channel))
            )
        return records

    def _additional_scores(self, channel) -> list[dict] | None:
        """Return the PG scores without a dedicated pg column, tracking their names."""
        scores = []
        for name in _ADDITIONAL_SCORES:
            value = _optional(getattr(channel, name))
            if value is not None:
                self._discovered_scores.add(name)
                scores.append({"score_name": name, "score_value": value, "higher_better": _CHANNEL_SCORES[name]})
        return scores or None
