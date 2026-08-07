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

import pandas as pd
import pyarrow as pa

from qpx.converters.base import resolve_columns
from qpx.converters.mappings import get_field_mappings
from qpx.converters.ptm import compute_precursor_mz
from qpx.converters.spectronaut.base_adapter import SpectronautBaseAdapter
from qpx.converters.spectronaut.constants import to_modifications, to_proforma
from qpx.core.cleavage import count_missed_cleavages
from qpx.core.data import FeatureSchema
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
]

# CV param mappings: (Spectronaut column, output cv_name)
_CV_MAPPINGS = [
    ("EG.iRTPredicted", "predicted_irt"),
    ("EG.iRTEmpirical", "empirical_irt"),
    ("EG.IsImputed", "is_imputed"),
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
        """Initialize feature adapter with empty caches."""
        super().__init__(**kwargs)
        self._mz_cache: dict[tuple[str, int], float | None] = {}
        self._resolved: dict | None = None

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
        self._load_spectronaut_report(spectronaut_report)

        # 2. SDRF for enzyme info
        enzyme_name: str | None = None
        if sdrf_path:
            enzyme_name = self._load_sdrf_enzyme(sdrf_path)

        # 3. Detect report columns (once, reused for all batches)
        report_cols = self._get_report_columns()
        self._resolved = resolve_columns(get_field_mappings("spectronaut", "feature"), report_cols)

        # 4. Pre-compute precursors: DuckDB temp table + Python modifications dict
        mods_lookup = self._register_precursor_lookup(enzyme_name)

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

    # ------------------------------------------------------------------
    # Data loading helpers
    # ------------------------------------------------------------------

    def _load_sdrf_enzyme(self, sdrf_path: str) -> str | None:
        """Load the first enzyme name from SDRF."""
        try:
            from qpx.core.sdrf import SDRFHandler

            handler = SDRFHandler(sdrf_path)
            enzymes = handler.get_enzymes()
            if enzymes:
                return str(enzymes[0])
        except (ValueError, KeyError, OSError, TypeError):
            self.logger.debug("Could not load enzyme from SDRF")
        return None

    def _discover_runs(self) -> list[str]:
        """Discover raw run values from the report for SQL filtering."""
        raw_runs = self._query_distinct_runs(self._resolved["run_file_name"])
        self.logger.info("Discovered %d runs from report", len(raw_runs))
        return raw_runs

    # ------------------------------------------------------------------
    # Precursor lookup
    # ------------------------------------------------------------------

    def _register_precursor_lookup(self, enzyme_name: str | None) -> dict:
        """Pre-compute per-precursor values and register in DuckDB.

        Spectronaut report is fragment-level; we need distinct precursors.
        The lookup table provides peptidoform, calculated_mz, and
        missed_cleavages for each unique (modified_sequence, sequence, charge).

        Returns:
            Dict mapping (modified_seq, sequence, charge) -> modifications list.
        """
        r = self._resolved
        mod_col = r["modified_sequence"]
        seq_col = r["sequence"]
        chg_col = r["charge"]

        precursors = self._conn.execute(
            sql_build(
                "SELECT DISTINCT $mc, $sc, $cc FROM report",
                mc=validate_identifier(mod_col),
                sc=validate_identifier(seq_col),
                cc=validate_identifier(chg_col),
            )
        ).fetchall()

        self.logger.info("Pre-computing %d unique precursors", len(precursors))

        mods_lookup: dict[tuple[str, str, int], list | None] = {}
        rows = []

        for modified_seq, sequence, charge in precursors:
            modified_seq = str(modified_seq)
            sequence = str(sequence)
            charge = int(charge)

            peptidoform = to_proforma(modified_seq)
            _, modifications = to_modifications(modified_seq, sequence)

            cache_key = (peptidoform, charge)
            if cache_key not in self._mz_cache:
                self._mz_cache[cache_key] = compute_precursor_mz(peptidoform, charge)
            calculated_mz = self._mz_cache[cache_key] or 0.0

            missed = count_missed_cleavages(sequence, enzyme_name) if enzyme_name else None

            mods_lookup[(modified_seq, sequence, charge)] = modifications
            rows.append(
                (
                    modified_seq,
                    sequence,
                    charge,
                    peptidoform,
                    calculated_mz,
                    missed,
                )
            )

        self._conn.execute("DROP TABLE IF EXISTS precursor_lookup")
        self._conn.from_df(
            pd.DataFrame(
                rows,
                columns=[
                    "modified_sequence",
                    "sequence",
                    "charge",
                    "peptidoform",
                    "calculated_mz",
                    "missed_cleavages",
                ],
            )
        ).create("precursor_lookup")
        self.logger.info("Registered precursor_lookup temp table in DuckDB")

        return mods_lookup

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

        # calculated_mz / observed_mz
        parts.append("COALESCE(CAST(FIRST(lk.calculated_mz) AS FLOAT), 0.0::FLOAT) AS calculated_mz")
        omz_col = r.get("observed_mz")
        if omz_col and _has(omz_col):
            parts.append(f"COALESCE({_safe_float_sql(omz_col)}, 0.0::FLOAT) AS observed_mz")
        else:
            parts.append("0.0::FLOAT AS observed_mz")

        # predicted_rt / rt (minutes → seconds)
        prt_col = r.get("predicted_rt")
        if prt_col and _has(prt_col):
            parts.append(f"CAST({_safe_float_sql(prt_col)} * {_RT_SECONDS_FACTOR} AS FLOAT) AS predicted_rt")
        else:
            parts.append("NULL::FLOAT AS predicted_rt")

        run_col = r["run_file_name"]
        parts.append(f"regexp_replace(FIRST(r.\"{run_col}\"), '\\.(mzML|raw|d)$', '') AS run_file_name")
        parts.append("[]::INTEGER[] AS scan")

        rt_col = r.get("rt")
        if rt_col and _has(rt_col):
            parts.append(f"CAST({_safe_float_sql(rt_col)} * {_RT_SECONDS_FACTOR} AS FLOAT) AS rt")
        else:
            parts.append("NULL::FLOAT AS rt")

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
            parts.append(f'CAST(FIRST(r."{proteotypic_col}") AS BOOLEAN) AS "unique"')
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

        # id_run_file_name
        run_col = r["run_file_name"]
        parts.append(f"regexp_replace(FIRST(r.\"{run_col}\"), '\\.(mzML|raw|d)$', '') AS id_run_file_name")
        return parts

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

        # --- Nested columns ---
        int_col = r["intensity"]
        parts.append(
            f"[STRUCT_PACK(label := 'LFQ', intensity := COALESCE({_safe_float_sql(int_col)}, 0.0::FLOAT))] AS intensities"
        )
        parts.append(self._build_additional_intensities_sql(r, report_cols))
        parts.append(self._build_additional_scores_sql(report_cols))
        parts.append(self._build_cv_params_sql(report_cols))

        # --- Helper columns for Python post-processing ---
        mod_col = r["modified_sequence"]
        pg_col = r["pg_accessions"]
        parts.append(f'FIRST(r."{mod_col}") AS _modified_sequence')
        parts.append(f'FIRST(r."{pg_col}") AS _pg_group')

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

        sql = sql_build(
            "SELECT $select FROM report r"
            " JOIN precursor_lookup lk"
            " ON r.$mod = lk.modified_sequence"
            " AND r.$seq = lk.sequence"
            " AND r.$chg = lk.charge"
            " WHERE r.$run IN ({run_placeholders})"
            " AND r.$pg IS NOT NULL"
            " GROUP BY r.$run, r.$mod, r.$seq, r.$chg"
            " HAVING $having",
            select=select_clause,
            mod=mod_id,
            seq=seq_id,
            chg=chg_id,
            run=run_id,
            pg=pg_id,
            having=having_clause,
        )

        return sql

    @staticmethod
    def _build_additional_intensities_sql(r: dict, report_cols: set[str]) -> str:
        """Build SQL for the additional_intensities nested column."""
        candidates = []
        pg_qty_col = r.get("pg_quantity")
        if pg_qty_col and pg_qty_col in report_cols:
            candidates.append((pg_qty_col, "pg_quantity"))
        omz_cal = r.get("observed_mz_calibrated")
        if omz_cal and omz_cal in report_cols:
            candidates.append((omz_cal, "precmz_calibrated"))

        if not candidates:
            return "NULL AS additional_intensities"

        entries, checks = [], []
        for col, name in candidates:
            chk = f'(FIRST(r."{col}") IS NOT NULL AND NOT isnan(CAST(FIRST(r."{col}") AS DOUBLE)))'
            checks.append(chk)
            entries.append(
                f"CASE WHEN {chk} THEN STRUCT_PACK("
                f"intensity_name := '{name}', "
                f'intensity_value := CAST(FIRST(r."{col}") AS FLOAT)) END'
            )
        ai_list = ", ".join(entries)
        ai_check = " OR ".join(checks)
        return (
            f"CASE WHEN {ai_check} "
            f"THEN [STRUCT_PACK(label := 'LFQ', "
            f"intensities := LIST_FILTER("
            f"[{ai_list}], x -> x IS NOT NULL))] "
            f"ELSE NULL END AS additional_intensities"
        )

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
        return f"CASE WHEN {sc} THEN LIST_FILTER([{sl}], x -> x IS NOT NULL) ELSE NULL END AS additional_scores"

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
        return f"CASE WHEN {cc} THEN LIST_FILTER([{cl}], x -> x IS NOT NULL) ELSE NULL END AS cv_params"

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

    @staticmethod
    def _build_pg_accessions_array(
        table: pa.Table,
        pg_acc_type: pa.DataType,
    ) -> pa.Array:
        """Build pg_accessions column from protein group strings."""
        pg_groups = table.column("_pg_group").to_pylist()
        pg_list = [
            [
                {"accession": acc, "start": None, "end": None, "pre": None, "post": None}
                for acc in (str(pg).split(";") if pg else [])
            ]
            or None
            for pg in pg_groups
        ]
        return pa.array(pg_list, type=pg_acc_type)

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
        placeholders = ", ".join(["?" for _ in runs])
        sql = sql_template.replace("{run_placeholders}", placeholders)

        table = self._conn.execute(sql, runs).fetch_arrow_table()
        if len(table) == 0:
            return None

        n = len(table)
        mods_array = self._build_modifications_array(table, mods_lookup, mods_type)
        pg_array = self._build_pg_accessions_array(table, pg_acc_type)

        table = table.drop("_modified_sequence").drop("_pg_group")

        # Assemble final table matching target schema
        columns = {}
        for field in target_schema:
            if field.name == "modifications":
                columns[field.name] = mods_array
            elif field.name == "pg_accessions":
                columns[field.name] = pg_array
            elif field.name in table.schema.names:
                col = table.column(field.name)
                if col.type != field.type:
                    try:
                        col = col.cast(field.type, safe=False)
                    except (pa.ArrowInvalid, pa.ArrowNotImplementedError):
                        col = pa.array(col.to_pylist(), type=field.type)
                columns[field.name] = col
            else:
                columns[field.name] = pa.nulls(n, type=field.type)

        return pa.table(columns, schema=target_schema)
