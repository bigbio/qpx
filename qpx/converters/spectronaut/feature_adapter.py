"""Spectronaut Feature adapter -- report to feature.parquet.

Performance strategy:
    - DuckDB VIEW over TSV/Parquet for memory efficiency
    - GROUP BY to deduplicate fragment-level rows to precursor level
    - Precursor lookup (peptidoform, calculated_mz, missed_cleavages) via DuckDB temp table JOIN
    - All scalar and nested fields (scores, intensities, cv_params) constructed in SQL
    - Only modifications + pg_accessions built in Python (requires PTM parsing)
    - Arrow tables written directly via write_table()
"""

from __future__ import annotations

import logging
from typing import Optional

import pyarrow as pa

from qpx.converters.base import resolve_columns
from qpx.converters.channel_labels import read_sdrf_labels
from qpx.converters.fragments import write_fragment_view
from qpx.converters.localization import LOCALIZATION_SCORE, localize_modifications
from qpx.converters.mappings import get_field_mappings
from qpx.converters.spectronaut.base_adapter import SpectronautBaseAdapter
from qpx.converters.spectronaut.constants import localization_sites, to_modifications, to_proforma
from qpx.converters.spectronaut.fragments import fragment_sql
from qpx.converters.spectronaut.report_layout import channel_columns, channel_labels
from qpx.converters.utils import pg_accessions_array
from qpx.core.constants import FEATURE
from qpx.core.data import FeatureSchema
from qpx.core.query import _escape_sql_string
from qpx.core.sql import sql_build, validate_identifier
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

# Score mappings: (Spectronaut column, output score name, higher_better)
_SCORE_MAPPINGS = [
    ("EG.Qvalue", "qvalue", False),
    ("PG.Qvalue", "pg_qvalue", False),
    ("PG.QValue (Run-Wise)", "pg_qvalue_runwise", False),
    ("EG.Cscore", "spectronaut_cscore", True),
    ("PG.Cscore", "spectronaut_pg_cscore", True),
    ("PG.Cscore (Run-Wise)", "spectronaut_pg_cscore_runwise", True),
    ("EG.PTMAssayProbability", "spectronaut_ptm_assay_probability", True),
]

# CV param mappings: (Spectronaut column, output cv_name)
_CV_MAPPINGS = [
    ("EG.iRTPredicted", "predicted_irt"),
    ("EG.iRTEmpirical", "empirical_irt"),
]

_RT_SECONDS_FACTOR = 60.0

# Producer-specific feature identity composite (bigbio/qpx#229). Spectronaut reports
# one feature per precursor per run, so the measured key is just peptidoform +
# charge + run (all in feature.yaml). Passed to FeatureWriter so feature_id hashes
# exactly these columns instead of the schema default.
_FEATURE_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name")


def _safe_float_sql(col: str) -> str:
    """Safe float via FIRST() aggregate — NULL/NaN → NULL."""
    return (
        f'CASE WHEN FIRST(r."{col}") IS NOT NULL '
        f'AND NOT isnan(CAST(FIRST(r."{col}") AS DOUBLE)) '
        f'THEN CAST(FIRST(r."{col}") AS FLOAT) END'
    )


def _measured_float_sql(col: str) -> str:
    """Positive float via FIRST() aggregate — 0, NULL and NaN (no measurement) → NULL."""
    value = f'CAST(FIRST(r."{col}") AS DOUBLE)'
    return f"CASE WHEN {value} > 0 AND NOT isnan({value}) THEN CAST({value} AS FLOAT) END"


def _safe_double_sql(col: str) -> str:
    """Safe double via FIRST() aggregate — NULL/NaN → NULL."""
    return (
        f'CASE WHEN FIRST(r."{col}") IS NOT NULL '
        f'AND NOT isnan(CAST(FIRST(r."{col}") AS DOUBLE)) '
        f'THEN CAST(FIRST(r."{col}") AS DOUBLE) END'
    )


class SpectronautFeatureAdapter(SpectronautBaseAdapter):
    """Convert Spectronaut report to ``feature.parquet``.

    Usage::

        with SpectronautFeatureAdapter() as adapter:
            adapter.convert(
                spectronaut_report="report.tsv",
                output_path="feature.parquet",
                sdrf_path="sdrf.tsv",
            )
    """

    def __init__(self, **kwargs):
        """Initialize the feature adapter."""
        super().__init__(**kwargs)
        self._resolved: dict | None = None
        # Column of each label, when the report quantifies channels.
        self._label_columns: dict[str, str] = {}

    def convert(
        self,
        spectronaut_report: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        qvalue_threshold: Optional[float] = None,
        file_num: int = 100,
        creator: str = "spectronaut",
    ) -> None:
        """Run the Spectronaut report -> feature.parquet conversion.

        Filtering is opt-in: with the default ``qvalue_threshold=None`` the
        report is converted as reported (no q-value filter). When a threshold is
        provided, precursors above it are dropped.
        """
        # 1. Create VIEW over report (lazy — no data loaded into memory)
        self._require_view(spectronaut_report, FEATURE)
        self._load_spectronaut_report(spectronaut_report)

        # 2. SDRF for enzyme info
        enzyme_name: str | None = None
        if sdrf_path:
            enzyme_name = self._load_sdrf_enzyme(sdrf_path)

        # 3. Detect report columns (once, reused for all batches)
        report_cols = self._get_report_columns()
        self._resolved = resolve_columns(get_field_mappings("spectronaut", "feature"), report_cols)
        self._load_channel_labels(report_cols, sdrf_path)

        # 4. Pre-compute precursors: DuckDB temp table + Python modifications dict
        mods_lookup = self._register_precursor_lookup(self._resolved, enzyme_name, to_proforma, to_modifications)

        # 5. Get target Arrow schema
        target_schema = FeatureSchema.get_arrow_schema()
        mods_type = target_schema.field("modifications").type
        pg_acc_type = target_schema.field("pg_accessions").type

        # 6. Track discovered scores upfront
        for sn_col, score_name, _ in _SCORE_MAPPINGS:
            if sn_col in report_cols:
                self._discovered_scores.add(score_name)

        # 7. Build SQL template (once, reused for all batches)
        sql_template = self._build_batch_sql(report_cols, qvalue_threshold)

        # 8. Discover runs
        run_names = self._discover_runs()

        # 9. Process in batches → Arrow tables → write directly
        with FeatureWriter(
            output_path,
            creator=creator,
            compression=self._compression,
            identity_composite=_FEATURE_IDENTITY_COMPOSITE,
        ) as writer:
            for i in range(0, len(run_names), file_num):
                batch_runs = run_names[i : i + file_num]
                self.logger.info("Processing runs %d-%d of %d", i + 1, min(i + file_num, len(run_names)), len(run_names))
                table = self._process_batch_arrow(
                    batch_runs,
                    sql_template,
                    mods_lookup,
                    mods_type,
                    pg_acc_type,
                    target_schema,
                )
                if table is not None:
                    writer.write_table(table)

        self.logger.info("Spectronaut feature conversion complete -> %s", output_path)

    def write_fragments(self, output_path: str) -> bool:
        """Write the report's fragment rows as the fragment view; False when it has none."""
        report_cols = self._get_report_columns()
        sql = fragment_sql(self._resolved, report_cols, labelled=bool(self._label_columns))
        if sql is None:
            self.logger.warning("The report has no F.FrgIon or F.PeakArea column: no fragment view")
            return False
        rows = write_fragment_view(
            self._conn, sql, output_path, _FEATURE_IDENTITY_COMPOSITE, creator="spectronaut", compression=self._compression
        )
        self.logger.info("Wrote %d fragment rows to %s", rows, output_path)
        return True

    # ------------------------------------------------------------------
    # Data loading helpers
    # ------------------------------------------------------------------

    def _load_channel_labels(self, report_cols: set[str], sdrf_path: Optional[str]) -> None:
        """Assign the report's quantity channels their labels, lightest first."""
        channels = channel_columns(report_cols)
        if not channels:
            return
        labels = channel_labels(channels, read_sdrf_labels(sdrf_path))
        self._label_columns = {labels[channel]: column for channel, column in channels.items()}
        if len(channels) > 1:
            self.logger.warning(
                "Spectronaut does not name its channels: assuming channels %s are %s, lightest first",
                list(channels),
                list(labels.values()),
            )

    def _discover_runs(self) -> list[str]:
        """Discover raw run values from the report for SQL filtering."""
        raw_runs = self._query_distinct_runs(self._resolved["run_file_name"])
        self.logger.info("Discovered %d runs from report", len(raw_runs))
        return raw_runs

    # ------------------------------------------------------------------
    # SQL template builder
    # ------------------------------------------------------------------

    @staticmethod
    def _build_core_select_parts(
        r: dict,
        report_cols: set[str],
    ) -> list[str]:
        """Build SQL SELECT parts for core scalar columns."""
        _has = report_cols.__contains__
        parts = []

        parts.append(f'FIRST(r."{r["sequence"]}") AS sequence')
        parts.append("FIRST(lk.peptidoform) AS peptidoform")
        parts.append(f'CAST(r."{r["charge"]}" AS SMALLINT) AS charge')

        # PEP
        pep_col = r.get("posterior_error_probability")
        if pep_col and _has(pep_col):
            parts.append(f"{_safe_double_sql(pep_col)} AS posterior_error_probability")
        else:
            parts.append("NULL::DOUBLE AS posterior_error_probability")

        # is_decoy
        is_decoy_col = r.get("is_decoy")
        if is_decoy_col and _has(is_decoy_col):
            parts.append(f'CAST(FIRST(r."{is_decoy_col}") AS BOOLEAN) AS is_decoy')
        else:
            parts.append("false AS is_decoy")

        # calculated_mz / observed_mz (FG.PrecMz is theoretical; the calibrated m/z is the run's)
        parts.append("CAST(FIRST(lk.calculated_mz) AS FLOAT) AS calculated_mz")
        omz_col = r.get("observed_mz")
        if omz_col and _has(omz_col):
            parts.append(f"{_safe_float_sql(omz_col)} AS observed_mz")
        else:
            parts.append("NULL::FLOAT AS observed_mz")

        # predicted_rt / rt (minutes → seconds)
        prt_col = r.get("predicted_rt")
        if prt_col and _has(prt_col):
            parts.append(f"CAST({_safe_float_sql(prt_col)} * {_RT_SECONDS_FACTOR} AS FLOAT) AS predicted_rt")
        else:
            parts.append("NULL::FLOAT AS predicted_rt")

        run_col = r["run_file_name"]
        # Match the pg adapter's case-insensitive, extended extension strip so
        # feature / pg / run keys agree on .wiff / .htrms / .RAW (bigbio/qpx#252).
        parts.append(f"regexp_replace(FIRST(r.\"{run_col}\"), '\\.(mzML|raw|d|wiff|htrms)$', '', 'i') AS run_file_name")
        parts.append("[]::INTEGER[] AS scan")

        rt_col = r.get("rt")
        if rt_col and _has(rt_col):
            parts.append(f"CAST({_safe_float_sql(rt_col)} * {_RT_SECONDS_FACTOR} AS FLOAT) AS rt")
        else:
            parts.append("NULL::FLOAT AS rt")

        for field in ("rt_start", "rt_stop"):
            column = r.get(field)
            parts.append(
                f"CAST({_safe_float_sql(column)} * {_RT_SECONDS_FACTOR} AS FLOAT) AS {field}"
                if column and _has(column)
                else f"NULL::FLOAT AS {field}"
            )

        im_col = r.get("ion_mobility")
        parts.append(
            f"{_measured_float_sql(im_col)} AS ion_mobility" if im_col and _has(im_col) else "NULL::FLOAT AS ion_mobility"
        )

        parts.append("CAST(FIRST(lk.missed_cleavages) AS SMALLINT) AS missed_cleavages")
        return parts

    @staticmethod
    def _build_meta_select_parts(
        r: dict,
        report_cols: set[str],
    ) -> list[str]:
        """Build SQL SELECT parts for metadata columns."""
        _has = report_cols.__contains__
        parts = []
        pg_col = r["pg_accessions"]

        parts.append(f"SPLIT_PART(FIRST(r.\"{pg_col}\"), ';', 1) AS anchor_protein")

        # unique (proteotypic)
        proteotypic_col = r.get("is_proteotypic")
        if proteotypic_col and _has(proteotypic_col):
            parts.append(
                f'NOT contains(FIRST(r."{pg_col}"), \';\') AND COALESCE(CAST(FIRST(r."{proteotypic_col}") AS BOOLEAN), true) AS "unique"'
            )
        else:
            parts.append(f'NOT contains(FIRST(r."{pg_col}"), \';\') AS "unique"')

        # pg_global_qvalue
        pgqv_col = r.get("pg_qvalue")
        if pgqv_col and _has(pgqv_col):
            parts.append(f"{_safe_double_sql(pgqv_col)} AS pg_global_qvalue")
        else:
            parts.append("NULL::DOUBLE AS pg_global_qvalue")

        # gg_accessions / gg_names
        genes_col = r.get("gg_names")
        if genes_col and _has(genes_col):
            ge = (
                f'CASE WHEN FIRST(r."{genes_col}") IS NOT NULL '
                f"AND FIRST(r.\"{genes_col}\") != '' "
                f"THEN STRING_SPLIT(FIRST(r.\"{genes_col}\"), ';') "
                f"ELSE NULL END"
            )
        else:
            ge = "NULL"
        parts.append(f"{ge} AS gg_accessions")
        parts.append(f"{ge} AS gg_names")

        return parts

    def _build_identification_select_parts(self, r: dict, report_cols: set[str]) -> list[str]:
        """Build SQL for the identifying run, and how the tool arrived at it when ``EG.IsImputed`` says."""
        run_col = r["run_file_name"]
        run_name = f"regexp_replace(FIRST(r.\"{run_col}\"), '\\.(mzML|raw|d|wiff|htrms)$', '', 'i')"
        imputed_col = r.get("is_imputed")
        if not imputed_col or imputed_col not in report_cols:
            return [f"{run_name} AS id_run_file_name"]
        # An imputed peak was not identified in this run: its identification came from another one.
        imputed = f'TRY_CAST(FIRST(r."{imputed_col}") AS BOOLEAN)'
        labels = ", ".join(f"'{label}'" for label in self._label_columns or ("LFQ",))
        return [
            f"CASE WHEN {imputed} THEN NULL ELSE {run_name} END AS id_run_file_name",
            f"CASE {imputed} WHEN TRUE THEN 'transferred' WHEN FALSE THEN 'direct' END AS identification_type",
            f"CASE WHEN {imputed} THEN [{labels}] END AS imputed_labels",
        ]

    def _build_batch_sql(
        self,
        report_cols: set[str],
        qvalue_threshold: Optional[float],
    ) -> str:
        """Build the SQL query template for batch processing.

        Spectronaut report has one row per fragment ion. We GROUP BY the
        precursor key (run_file_name, modified_sequence, sequence, charge)
        to get one row per precursor, using FIRST() for non-aggregated columns.
        """
        r = self._resolved
        parts = self._build_core_select_parts(r, report_cols)
        parts.extend(self._build_meta_select_parts(r, report_cols))
        parts.extend(self._build_identification_select_parts(r, report_cols))

        # --- Nested columns ---
        parts.append(self._build_intensities_sql(r))
        parts.append(self._build_additional_intensities_sql(r, report_cols))
        parts.append(self._build_additional_scores_sql(report_cols))
        parts.append(self._build_cv_params_sql(report_cols))

        # --- Helper columns for Python post-processing ---
        mod_col = r["modified_sequence"]
        pg_col = r["pg_accessions"]
        parts.append(f'FIRST(r."{mod_col}") AS _modified_sequence')
        parts.append(f'FIRST(r."{pg_col}") AS _pg_group')
        parts.append(self._ptm_localization_sql(report_cols))

        select_clause = ",\n        ".join(parts)

        # GROUP BY precursor key to deduplicate fragment rows
        run_id = validate_identifier(r["run_file_name"])
        seq_id = validate_identifier(r["sequence"])
        chg_id = validate_identifier(r["charge"])
        mod_id = validate_identifier(mod_col)
        pg_id = validate_identifier(pg_col)

        qv_col = r.get("qvalue")
        having_clause = "1=1"
        if qvalue_threshold is not None and qv_col and qv_col in report_cols:
            qv_id = validate_identifier(qv_col)
            having_clause = f"1=1 AND CAST(FIRST(r.{qv_id}) AS DOUBLE) < {float(qvalue_threshold)}"
        elif qvalue_threshold is not None:
            self.logger.warning("The report has no EG.Qvalue column: --qvalue-threshold is not applied")

        sql = sql_build(
            "SELECT $select FROM report r"
            " JOIN precursor_lookup lk"
            " ON r.$mod = lk.modified_sequence"
            " AND r.$seq = lk.sequence"
            " AND r.$chg = lk.charge"
            " WHERE r.$run IN ({run_placeholders})"
            " AND r.$pg IS NOT NULL"
            " GROUP BY r.$run, r.$mod, r.$seq, r.$chg"
            " HAVING $having"
            " ORDER BY r.$mod, r.$chg, r.$run",
            select=select_clause,
            mod=mod_id,
            seq=seq_id,
            chg=chg_id,
            run=run_id,
            pg=pg_id,
            having=having_clause,
        )

        return sql

    def _build_intensities_sql(self, r: dict) -> str:
        """Build SQL for the intensities column: one entry per channel, or the label-free quantity."""
        columns = self._label_columns or {"LFQ": r["intensity"]}
        entries = [
            f"STRUCT_PACK(label := '{label}', intensity := {_safe_float_sql(column)})" for label, column in columns.items()
        ]
        return f"[{', '.join(entries)}] AS intensities"

    def _build_additional_intensities_sql(self, r: dict, report_cols: set[str]) -> str:
        """Build SQL for additional_intensities: ``PG.Quantity``, which belongs to no single channel of a labelled report."""
        column = r.get("pg_quantity")
        if self._label_columns or not column or column not in report_cols:
            return "NULL AS additional_intensities"
        value = _safe_float_sql(column)
        return (
            f"CASE WHEN {value} IS NOT NULL THEN [STRUCT_PACK(label := 'LFQ', "
            f"intensities := [STRUCT_PACK(intensity_name := 'pg_quantity', intensity_value := {value})])] "
            f"END AS additional_intensities"
        )

    @staticmethod
    def _ptm_localization_sql(report_cols: set[str]) -> str:
        """Build SQL for the helper column of ``EG.PTMLocalizationProbabilities``."""
        column = "EG.PTMLocalizationProbabilities"
        value = f'FIRST(r."{column}")' if column in report_cols else "NULL"
        return f"CAST({value} AS VARCHAR) AS _ptm_localization"

    @staticmethod
    def _build_additional_scores_sql(report_cols: set[str]) -> str:
        """Build SQL for the additional_scores nested column."""
        entries, checks = [], []
        for sn_col, score_name, higher_better in _SCORE_MAPPINGS:
            if sn_col in report_cols:
                hb = "true" if higher_better else "false"
                chk = f'(FIRST(r."{sn_col}") IS NOT NULL AND NOT isnan(CAST(FIRST(r."{sn_col}") AS DOUBLE)))'
                checks.append(chk)
                entries.append(
                    f"CASE WHEN {chk} THEN STRUCT_PACK("
                    f"score_name := '{score_name}', "
                    f'score_value := CAST(FIRST(r."{sn_col}") AS DOUBLE), '
                    f"higher_better := {hb}) END"
                )
        if not entries:
            return "NULL AS additional_scores"
        sl = ", ".join(entries)
        sc = " OR ".join(checks)
        return f"CASE WHEN {sc} THEN LIST_FILTER([{sl}], lambda x: x IS NOT NULL) ELSE NULL END AS additional_scores"

    @staticmethod
    def _build_cv_params_sql(report_cols: set[str]) -> str:
        """Build SQL for the cv_params nested column."""
        entries, checks = [], []
        for sn_col, cv_name in _CV_MAPPINGS:
            if sn_col in report_cols:
                chk = f'(FIRST(r."{sn_col}") IS NOT NULL AND CAST(FIRST(r."{sn_col}") AS VARCHAR) != \'\')'
                checks.append(chk)
                entries.append(
                    f"CASE WHEN {chk} THEN STRUCT_PACK("
                    f"cv_name := '{cv_name}', "
                    f'cv_value := CAST(FIRST(r."{sn_col}") AS VARCHAR)) END'
                )
        if not entries:
            return "NULL AS cv_params"
        cl = ", ".join(entries)
        cc = " OR ".join(checks)
        return f"CASE WHEN {cc} THEN LIST_FILTER([{cl}], lambda x: x IS NOT NULL) ELSE NULL END AS cv_params"

    # ------------------------------------------------------------------
    # Batch processing
    # ------------------------------------------------------------------

    @staticmethod
    def _build_modifications_array(
        table: pa.Table,
        mods_lookup: dict,
        mods_type: pa.DataType,
    ) -> pa.Array:
        """Build modifications column from Python lookup."""
        mod_seqs = table.column("_modified_sequence").to_pylist()
        seqs = table.column("sequence").to_pylist()
        charges = table.column("charge").to_pylist()
        mods_list = [mods_lookup.get((str(ms), str(s), int(c))) for ms, s, c in zip(mod_seqs, seqs, charges)]
        return pa.array(mods_list, type=mods_type)

    def _localized(self, table: pa.Table, mods_array: pa.Array, target_schema: pa.Schema) -> tuple[pa.Array, pa.Array]:
        """Modifications with their site scores, and the candidates, from each row's localization probabilities."""
        localized = [
            localize_modifications(mods, sequence, localization_sites(localization, sequence))
            if localization and "%" in localization
            else (mods, None)
            for mods, sequence, localization in zip(
                mods_array.to_pylist(), table.column("sequence").to_pylist(), table.column("_ptm_localization").to_pylist()
            )
        ]
        candidates = [row_candidates for _, row_candidates in localized]
        if any(candidates):
            self._discovered_scores.add(LOCALIZATION_SCORE)
        return (
            pa.array([mods for mods, _ in localized], type=mods_array.type),
            pa.array(candidates, type=target_schema.field("modification_candidates").type),
        )

    def _process_batch_arrow(
        self,
        runs: list[str],
        sql_template: str,
        mods_lookup: dict,
        mods_type: pa.DataType,
        pg_acc_type: pa.DataType,
        target_schema: pa.Schema,
    ) -> Optional[pa.Table]:
        """Execute SQL batch, add Python-only columns, cast to schema."""
        placeholders = ", ".join(f"'{_escape_sql_string(r)}'" for r in runs)
        sql = sql_template.replace("{run_placeholders}", placeholders)

        table = self._conn.sql(sql).to_arrow_table()
        if len(table) == 0:
            return None

        mods_array = self._build_modifications_array(table, mods_lookup, mods_type)
        pg_array = pg_accessions_array(table.column("_pg_group"), pg_acc_type)
        mods_array, candidates = self._localized(table, mods_array, target_schema)

        table = table.drop("_modified_sequence").drop("_pg_group").drop("_ptm_localization")
        table = table.append_column("modification_candidates", candidates)

        return self._assemble_feature_table(table, target_schema, mods_array, pg_array)
