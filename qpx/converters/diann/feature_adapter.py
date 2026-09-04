"""DIA-NN Feature adapter -- report to feature.parquet.

Performance strategy:
    - DuckDB VIEW over parquet for memory efficiency (no in-memory table)
    - Precursor lookup (peptidoform, mz, missed_cleavages) via DuckDB temp table JOIN
    - All scalar and nested fields (scores, intensities, cv_params) constructed in SQL
    - Only modifications + pg_accessions built in Python (requires PTM parsing)
    - Arrow tables written directly via write_table() — no Python dict intermediaries
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pandas as pd
import pyarrow as pa

from qpx.converters.base import resolve_columns
from qpx.converters.diann.base_adapter import DiaNNBaseAdapter
from qpx.converters.diann.constants import to_modifications, to_proforma
from qpx.converters.mappings import get_field_mappings
from qpx.converters.ptm import compute_precursor_mz
from qpx.core.cleavage import count_missed_cleavages
from qpx.core.sql import sql_build, validate_identifier
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

# Score mappings: (DIA-NN column, output score name, higher_better)
_SCORE_MAPPINGS = [
    ("Q.Value", "precursor_qvalue", False),
    ("PG.Q.Value", "pg_qvalue", False),
    ("Global.Q.Value", "global_qvalue", False),
    ("CScore", "diann_cscore", True),
    ("Evidence", "diann_evidence", True),
    ("Spectrum.Similarity", "diann_spectrum_similarity", True),
    ("Ms1.Profile.Corr", "diann_ms1_profile_corr", True),
    ("Mass.Evidence", "diann_mass_evidence", True),
    ("Averagine", "diann_averagine", True),
    ("Quantity.Quality", "diann_quantity_quality", True),
    ("Lib.Q.Value", "diann_lib_qvalue", False),
    ("Lib.PG.Q.Value", "diann_lib_pg_qvalue", False),
    ("Protein.Q.Value", "diann_protein_qvalue", False),
    ("Decoy.CScore", "diann_decoy_cscore", True),
    ("Decoy.Evidence", "diann_decoy_evidence", True),
    ("Translated.Quality", "diann_translated_quality", True),
    ("Translated.Q.Value", "diann_translated_qvalue", False),
]

# CV param mappings: (DIA-NN column, output cv_name)
_CV_MAPPINGS = [
    ("Proteotypic", "proteotypic"),
    ("iRT", "irt"),
    ("Predicted.iRT", "predicted_irt"),
    ("Predicted.IM", "predicted_im"),
    ("Predicted.iIM", "predicted_iim"),
    ("iIM", "iim"),
]


_RT_SECONDS_FACTOR = 60.0

# Producer-specific feature identity composite (bigbio/qpx#229). DIA-NN reports one
# feature per precursor per run, so the measured key is just peptidoform + charge +
# run (all in feature.yaml). Passed to FeatureWriter so feature_id hashes exactly
# these columns instead of the schema default.
_FEATURE_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name")


def _safe_float_sql(column: str) -> str:
    """Build a SQL expression that converts a finite value to FLOAT."""
    return f'CASE WHEN r."{column}" IS NOT NULL AND NOT isnan(CAST(r."{column}" AS DOUBLE)) THEN CAST(r."{column}" AS FLOAT) END'


def _safe_double_sql(column: str) -> str:
    """Build a SQL expression that converts a finite value to DOUBLE."""
    return f'CASE WHEN r."{column}" IS NOT NULL AND NOT isnan(CAST(r."{column}" AS DOUBLE)) THEN CAST(r."{column}" AS DOUBLE) END'


class DiannFeatureAdapter(DiaNNBaseAdapter):
    """Convert DIA-NN report to ``feature.parquet``.

    Usage::

        with DiannFeatureAdapter() as adapter:
            adapter.convert(
                diann_report="report.tsv",
                output_path="feature.parquet",
                mzml_info_folder="ms_info/",
                sdrf_path="sdrf.tsv",
            )
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._mz_cache: dict[tuple[str, int], float | None] = {}

    def convert(
        self,
        diann_report: str,
        output_path: str,
        mzml_info_folder: Optional[str] = None,
        sdrf_path: Optional[str] = None,
        qvalue_threshold: Optional[float] = None,
        file_num: int = 100,  # VIEW-based lazy IO reduces per-batch memory
        creator: str = "diann",
    ) -> None:
        """Run the DIA-NN report -> feature.parquet conversion.

        Filtering is opt-in (bigbio/qpx#241). ``qvalue_threshold`` defaults to
        ``None``, which emits every feature the report contains — DIA-NN already
        FDR-filters its main report and all per-row q-value columns (``precursor_qvalue``,
        ``global_qvalue``, ``diann_lib_qvalue``, ``diann_translated_qvalue`` …)
        are carried through unconditionally for downstream filtering. When a
        threshold is provided, the precursor ``Q.Value`` filter is applied.
        """
        # 1. Create VIEW over report (lazy — no data loaded into memory)
        self._load_diann_report(diann_report)

        # 2. SDRF for sample mapping and enzyme info
        enzyme_name: str | None = None
        if sdrf_path:
            enzyme_name = self._load_sdrf_enzyme(sdrf_path)

        # 3. Detect report columns (once, reused for all batches)
        report_cols = self._get_report_columns()
        self._resolved = resolve_columns(get_field_mappings("diann", "feature"), report_cols)
        decoy_col = self._resolved.get("decoy")
        has_decoy_col = decoy_col is not None
        channel_col = self._register_channel_labels(report_cols, sdrf_path)

        # 4. Pre-compute precursors: DuckDB temp table + Python modifications dict
        mods_lookup = self._register_precursor_lookup(enzyme_name)

        # 5. Get target Arrow schema
        target_schema = FeatureWriter._schema_class.get_arrow_schema()
        mods_type = target_schema.field("modifications").type
        pg_acc_type = target_schema.field("pg_accessions").type

        # 6. Track discovered scores upfront (all score columns present in report)
        for diann_col, score_name, _ in _SCORE_MAPPINGS:
            if diann_col in report_cols:
                self._discovered_scores.add(score_name)

        # 7. Build SQL template (once, reused for all batches)
        sql_template = self._build_batch_sql(
            report_cols,
            has_decoy_col,
            qvalue_threshold,
            channel_col,
            target_schema,
        )

        # 8. Discover runs
        run_names = self._discover_runs(mzml_info_folder)

        # 9. Process in batches → Arrow tables → write directly
        with FeatureWriter(
            output_path,
            creator=creator,
            compression=self._compression,
            identity_composite=_FEATURE_IDENTITY_COMPOSITE,
        ) as writer:
            for i in range(0, len(run_names), file_num):
                batch_runs = run_names[i : i + file_num]
                self.logger.info(f"Processing runs {i + 1}-{min(i + file_num, len(run_names))} of {len(run_names)}")
                table = self._process_batch_arrow(
                    batch_runs,
                    sql_template,
                    mods_lookup,
                    mods_type,
                    pg_acc_type,
                    target_schema,
                    mzml_info_folder,
                )
                if table is not None:
                    writer.write_table(table)

        self.logger.info(f"DIA-NN feature conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading helpers
    # ------------------------------------------------------------------

    def _get_report_columns(self) -> set[str]:
        """Return the set of column names in the DuckDB report table/view."""
        return {
            c[0]
            for c in self._conn.execute("SELECT column_name FROM information_schema.columns WHERE table_name='report'").fetchall()
        }

    def _load_sdrf_enzyme(self, sdrf_path: str) -> str | None:
        """Load the first enzyme name from SDRF for missed-cleavage computation."""
        try:
            from qpx.core.sdrf import SDRFHandler

            handler = SDRFHandler(sdrf_path)
            enzymes = handler.get_enzymes()
            if enzymes:
                return str(enzymes[0])
        except (OSError, KeyError, TypeError, ValueError):
            self.logger.debug("Could not load enzyme from SDRF")
        return None

    def _discover_runs(self, mzml_info_folder: Optional[str]) -> list[str]:
        """Discover run names from ms_info files or from the report."""
        if mzml_info_folder:
            info_files = list(Path(mzml_info_folder).glob("*_ms_info.parquet"))
        else:
            info_files = []

        if info_files:
            run_names = [f.stem.replace("_ms_info", "") for f in info_files]
            self.logger.info(f"Found {len(run_names)} MS info files")
            # The feature SQL filters `run IN (run_names)`, so any report Run
            # without a matching *_ms_info.parquet is silently dropped from
            # feature.parquet. Surface the discrepancy (bigbio/qpx#244) so the
            # loss is visible; the drop behaviour itself is unchanged.
            self._warn_on_dropped_runs(set(run_names))
        else:
            run_col = self._resolved["run_file_name"]
            qcol = validate_identifier(run_col)
            rows = self._conn.execute(
                sql_build(
                    """
                    SELECT DISTINCT regexp_replace(
                        CAST($col AS VARCHAR),
                        '(?i)\\.(mzML|raw|d|wiff|htrms)$',
                        ''
                    ) AS run_file_name
                    FROM report
                    ORDER BY run_file_name
                    """,
                    col=qcol,
                )
            ).fetchall()
            run_names = [row[0] for row in rows]
            self.logger.info(f"Discovered {len(run_names)} runs from report")

        return run_names

    def _report_run_names(self) -> set[str]:
        """Return the distinct extension-stripped Run values in the report."""
        run_col = self._resolved["run_file_name"]
        rows = self._conn.execute(
            sql_build(
                """
                SELECT DISTINCT regexp_replace(
                    CAST($col AS VARCHAR),
                    '(?i)\\.(mzML|raw|d|wiff|htrms)$',
                    ''
                ) AS run_file_name
                FROM report
                """,
                col=validate_identifier(run_col),
            )
        ).fetchall()
        return {row[0] for row in rows}

    def _warn_on_dropped_runs(self, kept_runs: set[str]) -> None:
        """Warn about report runs with no matching ms_info file.

        When runs are discovered from ``*_ms_info.parquet`` stems, the feature
        SQL keeps only those runs. Any report ``Run`` absent from *kept_runs* is
        dropped from ``feature.parquet`` without reconciliation; log the dropped
        run names so the loss is visible (bigbio/qpx#244).
        """
        dropped = sorted(self._report_run_names() - kept_runs)
        if dropped:
            self.logger.warning(
                "%d report run(s) have no matching *_ms_info.parquet and will be dropped from feature.parquet: %s",
                len(dropped),
                ", ".join(dropped),
            )

    # ------------------------------------------------------------------
    # Precursor lookup (DuckDB temp table + Python modifications dict)
    # ------------------------------------------------------------------

    def _register_precursor_lookup(self, enzyme_name: str | None) -> dict:
        """Pre-compute per-precursor values and register flat fields in DuckDB.

        Returns:
            Dict mapping (modified_seq, sequence, charge) -> modifications list.
        """
        mod_col = self._resolved["modified_sequence"]
        seq_col = self._resolved["sequence"]
        chg_col = self._resolved["charge"]

        precursors = self._conn.execute(
            sql_build(
                "SELECT DISTINCT $mc, $sc, $cc FROM report",
                mc=validate_identifier(mod_col),
                sc=validate_identifier(seq_col),
                cc=validate_identifier(chg_col),
            )
        ).fetchall()

        self.logger.info(f"Pre-computing {len(precursors):,} unique precursors")

        mods_lookup: dict[tuple[str, str, int], list | None] = {}
        rows = []

        for modified_seq, sequence, charge in precursors:
            modified_seq = str(modified_seq)
            sequence = str(sequence)
            charge = int(charge)

            peptidoform = to_proforma(modified_seq)
            _, modifications = to_modifications(modified_seq, sequence)

            cache_key = (modified_seq, charge)
            if cache_key not in self._mz_cache:
                self._mz_cache[cache_key] = compute_precursor_mz(modified_seq, charge)
            calculated_mz = self._mz_cache[cache_key] or 0.0

            missed = count_missed_cleavages(sequence, enzyme_name) if enzyme_name else None

            mods_lookup[(modified_seq, sequence, charge)] = modifications
            rows.append((modified_seq, sequence, charge, peptidoform, calculated_mz, missed))

        # Register flat lookup as DuckDB temp table for SQL JOIN
        self._conn.execute("DROP TABLE IF EXISTS precursor_lookup")
        self._conn.from_df(
            pd.DataFrame(
                rows,
                columns=["modified_sequence", "sequence", "charge", "peptidoform", "calculated_mz", "missed_cleavages"],
            )
        ).create("precursor_lookup")
        self.logger.info("Registered precursor_lookup temp table in DuckDB")

        return mods_lookup

    # ------------------------------------------------------------------
    # SQL template builder
    # ------------------------------------------------------------------

    @staticmethod
    def _build_core_select_parts(
        resolved: dict,
        report_cols: set[str],
        has_decoy_col: bool,
    ) -> list[str]:
        """Build SQL SELECT parts for core feature columns."""
        has_column = report_cols.__contains__
        parts = [
            f'r."{resolved["sequence"]}" AS sequence',
            "lk.peptidoform",
            f'CAST(r."{resolved["charge"]}" AS SMALLINT) AS charge',
        ]

        pep_col = resolved.get("posterior_error_probability")
        parts.append(
            f"{_safe_double_sql(pep_col)} AS posterior_error_probability"
            if has_column(pep_col)
            else "NULL::DOUBLE AS posterior_error_probability"
        )

        # DIA-NN Q.Value is precursor-level and remains in additional_scores;
        # DIA-NN does not report a peptide-level q-value.
        parts.append("NULL::DOUBLE AS peptide_qvalue")

        pg_col = resolved["pg_accessions"]
        if has_decoy_col:
            parts.append(f'CAST(r."{resolved["decoy"]}" AS BOOLEAN) AS is_decoy')
        else:
            parts.append(
                f"(starts_with(r.\"{pg_col}\", 'DECOY_') OR starts_with(r.\"{pg_col}\", 'decoy_') "
                f"OR starts_with(r.\"{pg_col}\", 'rev_') OR starts_with(r.\"{pg_col}\", 'REV_')) AS is_decoy"
            )

        parts.append("COALESCE(CAST(lk.calculated_mz AS FLOAT), 0.0::FLOAT) AS calculated_mz")

        # observed_mz: emit a truthful NULL when the report carries no Precursor.Mz
        # (bigbio/qpx#244). A sentinel 0.0 looks like a real measurement and
        # corrupts downstream mass-error / PPM math; NULL propagates as "unknown".
        # The per-row _safe_float_sql already yields NULL for null/NaN cells, so no
        # COALESCE-to-zero is applied even when the column is present.
        observed_mz_col = resolved.get("observed_mz")
        parts.append(
            f"{_safe_float_sql(observed_mz_col)} AS observed_mz"
            if observed_mz_col and has_column(observed_mz_col)
            else "NULL::FLOAT AS observed_mz"
        )

        # mass_error_ppm: prefer a reported column, otherwise DERIVE it. DIA-NN
        # emits no mass-error column, so this was always null even though both
        # inputs are present — on PXD030304 22% of rows carry an observed m/z that
        # genuinely differs from the theoretical one, and that error was discarded
        # (bigbio/qpx#298). Guarded on both values being positive: calculated_mz is
        # COALESCEd to 0.0 above, and observed_mz is null when the report has no
        # Precursor.Mz, so an unmeasured row stays null rather than dividing by zero
        # or reporting a fabricated error. The formula matches docs/spec/feature.md.
        mass_error_col = resolved.get("mass_error_ppm")
        if mass_error_col and has_column(mass_error_col):
            parts.append(f"{_safe_float_sql(mass_error_col)} AS mass_error_ppm")
        elif observed_mz_col and has_column(observed_mz_col):
            observed_expr = _safe_float_sql(observed_mz_col)
            parts.append(
                "CASE WHEN lk.calculated_mz > 0 AND "
                f"({observed_expr}) > 0 "
                f"THEN CAST((({observed_expr}) - lk.calculated_mz) / lk.calculated_mz * 1e6 AS FLOAT) "
                "END AS mass_error_ppm"
            )
        else:
            parts.append("NULL::FLOAT AS mass_error_ppm")

        predicted_rt_col = resolved.get("predicted_rt")
        parts.append(
            f"CAST({_safe_float_sql(predicted_rt_col)} * {_RT_SECONDS_FACTOR} AS FLOAT) AS predicted_rt"
            if predicted_rt_col and has_column(predicted_rt_col)
            else "NULL::FLOAT AS predicted_rt"
        )

        run_col = resolved["run_file_name"]
        parts.append(f"regexp_replace(r.\"{run_col}\", '(?i)\\.(mzML|raw|d|wiff|htrms)$', '') AS run_file_name")

        ms2_col = resolved.get("ms2_scan")
        if ms2_col and has_column(ms2_col):
            parts.append(
                f'CASE WHEN r."{ms2_col}" IS NOT NULL THEN [CAST(r."{ms2_col}" AS INTEGER)] ELSE []::INTEGER[] END AS scan'
            )
        else:
            parts.append("[]::INTEGER[] AS scan")

        rt_col = resolved.get("rt")
        parts.append(
            f"CAST({_safe_float_sql(rt_col)} * {_RT_SECONDS_FACTOR} AS FLOAT) AS rt"
            if rt_col and has_column(rt_col)
            else "NULL::FLOAT AS rt"
        )

        ion_mobility_col = resolved.get("ion_mobility")
        parts.append(
            f"{_safe_float_sql(ion_mobility_col)} AS ion_mobility"
            if ion_mobility_col and has_column(ion_mobility_col)
            else "NULL::FLOAT AS ion_mobility"
        )
        parts.extend(
            [
                "CAST(lk.missed_cleavages AS SMALLINT) AS missed_cleavages",
                # Blank/whitespace Protein.Group → NULL (e.g. DIA-NN unmapped rows).
                f"NULLIF(TRIM(SPLIT_PART(r.\"{pg_col}\", ';', 1)), '') AS anchor_protein",
            ]
        )
        return parts

    @staticmethod
    def _build_metadata_select_parts(
        resolved: dict,
        report_cols: set[str],
    ) -> list[str]:
        """Build SQL SELECT parts for feature metadata columns."""
        has_column = report_cols.__contains__
        parts: list[str] = []
        pg_col = resolved["pg_accessions"]

        proteotypic_col = resolved.get("proteotypic")
        if proteotypic_col and has_column(proteotypic_col):
            parts.append(f'CAST(r."{proteotypic_col}" AS INTEGER) = 1 AS "unique"')
        else:
            parts.append(f'NOT contains(r."{pg_col}", \';\') AS "unique"')

        pg_qvalue_col = resolved.get("pg_global_qvalue")
        parts.append(
            f"{_safe_double_sql(pg_qvalue_col)} AS pg_global_qvalue"
            if pg_qvalue_col and has_column(pg_qvalue_col)
            else "NULL::DOUBLE AS pg_global_qvalue"
        )

        for field in ("ion_mobility_start", "ion_mobility_stop"):
            column = resolved.get(field)
            parts.append(f"{_safe_float_sql(column)} AS {field}" if column and has_column(column) else f"NULL::FLOAT AS {field}")

        genes_col = resolved.get("gg_names")
        if genes_col and has_column(genes_col):
            genes_expr = (
                f'CASE WHEN r."{genes_col}" IS NOT NULL '
                f"AND r.\"{genes_col}\" != '' "
                f"THEN STRING_SPLIT(r.\"{genes_col}\", ';') "
                "ELSE NULL END"
            )
        else:
            genes_expr = "NULL"
        parts.append(f"{genes_expr} AS gg_accessions")
        parts.append(f"{genes_expr} AS gg_names")

        run_col = resolved["run_file_name"]
        parts.append(f"regexp_replace(r.\"{run_col}\", '(?i)\\.(mzML|raw|d|wiff|htrms)$', '') AS id_run_file_name")

        for field in ("rt_start", "rt_stop"):
            column = resolved.get(field)
            parts.append(
                f"CAST({_safe_float_sql(column)} * {_RT_SECONDS_FACTOR} AS FLOAT) AS {field}"
                if column and has_column(column)
                else f"NULL::FLOAT AS {field}"
            )
        return parts

    def _build_batch_sql(
        self,
        report_cols: set[str],
        has_decoy_col: bool,
        qvalue_threshold: Optional[float],
        channel_col: str | None,
        target_schema: pa.Schema,
    ) -> str:
        """Build the SQL query template for batch processing.

        Constructs all output columns (including nested structs) in DuckDB SQL.
        Only modifications and pg_accessions are handled in Python post-processing.
        """
        r = self._resolved
        run_col = r["run_file_name"]
        qv_col = r.get("precursor_qvalue")
        if qvalue_threshold is not None and qv_col is None:
            raise ValueError("qvalue_threshold requires the DIA-NN report column Q.Value")
        pg_col = r["pg_accessions"]

        has_column = report_cols.__contains__
        parts = self._build_core_select_parts(r, report_cols, has_decoy_col)
        parts.extend(self._build_metadata_select_parts(r, report_cols))

        # --- Nested: intensities ---
        label_expr, channel_helper_parts, channel_join = self._build_channel_sql_parts(channel_col, qv_col)
        int_col = r["intensity"]
        parts.append(
            f"[STRUCT_PACK(label := {label_expr}, intensity := COALESCE({_safe_float_sql(int_col)}, 0.0::FLOAT))] AS intensities"
        )

        # --- Nested columns built by helpers ---
        parts.append(self._build_additional_intensities_sql(r, has_column, label_expr))
        parts.append(self._build_additional_scores_sql(has_column))
        parts.append(self._build_cv_params_sql(r, has_column))

        # --- Helper columns (for Python post-processing, dropped later) ---
        mod_col = r["modified_sequence"]
        seq_col = r["sequence"]
        chg_col = r["charge"]
        parts.append(f'r."{mod_col}" AS _modified_sequence')
        parts.append(f'r."{pg_col}" AS _pg_group')
        parts.extend(channel_helper_parts)

        # --- Build full query ---
        select_clause = ",\n        ".join(parts)

        # Precursor Q.Value filter is opt-in (bigbio/qpx#241). With no threshold
        # the report is converted as reported (DIA-NN already FDR-filters it);
        # only apply the cutoff when a threshold is explicitly provided.
        qvalue_filter = ""
        if qvalue_threshold is not None:
            qvalue_filter = sql_build(
                "AND r.$qv_col < $qv_threshold",
                qv_col=validate_identifier(qv_col),
                qv_threshold=str(float(qvalue_threshold)),
            )

        row_sql = sql_build(
            """
        SELECT
            $select_clause
        FROM report r
        JOIN precursor_lookup lk
            ON r.$mod_col = lk.modified_sequence
            AND r.$seq_col = lk.sequence
            AND r.$chg_col = lk.charge
        $channel_join
        WHERE regexp_replace(
                  CAST(r.$run_col AS VARCHAR),
                  '(?i)\\.(mzML|raw|d|wiff|htrms)$',
                  ''
              ) IN ({run_placeholders})
          $qvalue_filter
          AND r.$pg_col IS NOT NULL
        """,
            select_clause=select_clause,
            mod_col=validate_identifier(mod_col),
            seq_col=validate_identifier(seq_col),
            chg_col=validate_identifier(chg_col),
            channel_join=channel_join,
            run_col=validate_identifier(run_col),
            qvalue_filter=qvalue_filter,
            pg_col=validate_identifier(pg_col),
        )

        if channel_col:
            return self._aggregate_channel_rows_sql(row_sql, target_schema)
        return row_sql

    @staticmethod
    def _build_channel_sql_parts(channel_col: str | None, qvalue_col: str | None) -> tuple[str, list[str], str]:
        """Return the label, helper columns, and join used for channel-aware rows."""
        if channel_col is None:
            return "'LFQ'", [], ""

        qvalue_expr = f"{_safe_double_sql(qvalue_col)} AS _qvalue" if qvalue_col else "NULL::DOUBLE AS _qvalue"
        channel_join = sql_build(
            """
            JOIN diann_channel_labels cl
              ON CAST(r.$channel_col AS VARCHAR) = cl.channel_value
            """,
            channel_col=validate_identifier(channel_col),
        )
        return "cl.canonical_label", [qvalue_expr], channel_join

    @staticmethod
    def _aggregate_channel_rows_sql(row_sql: str, target_schema: pa.Schema) -> str:
        """Collapse DIA-NN's one-row-per-channel output to one QPX feature row."""
        group_fields = {"sequence", "charge", "run_file_name"}
        python_fields = {
            "modifications",
            "pg_accessions",
            "pg_positions",
            # Derived identity + optional cross-refs are not produced by the
            # channel row SQL: feature_id is stamped by the writer, psm_ids/pg_ids
            # stay null. Skip them here; they are filled as NULL in final assembly.
            "feature_id",
            "psm_ids",
            "pg_ids",
            # Producer-specific optional fields not emitted by DIA-NN.
            "quantification_unit_id",
            "compensation_voltage",
            "consensus_rt",
        }
        select_parts: list[str] = []
        for field in target_schema:
            name = field.name
            quoted = validate_identifier(name)
            if name in python_fields:
                continue
            if name in group_fields:
                select_parts.append(quoted)
            elif name == "intensities":
                select_parts.append("LIST(intensities[1] ORDER BY intensities[1].label) AS intensities")
            elif name == "additional_intensities":
                select_parts.append(
                    "LIST(additional_intensities[1] "
                    "ORDER BY intensities[1].label) "
                    "FILTER (WHERE additional_intensities IS NOT NULL) "
                    "AS additional_intensities"
                )
            else:
                select_parts.append(f"FIRST({quoted} ORDER BY _qvalue NULLS LAST) AS {quoted}")

        select_parts.extend(["_modified_sequence", "_pg_group"])
        select_clause = ",\n            ".join(select_parts)
        return sql_build(
            """
            WITH channel_rows AS (
                $row_sql
            ),
            best_channel_rows AS (
                SELECT *
                FROM channel_rows
                QUALIFY ROW_NUMBER() OVER (
                    PARTITION BY sequence,
                                 charge,
                                 run_file_name,
                                 _modified_sequence,
                                 _pg_group,
                                 intensities[1].label
                    ORDER BY _qvalue NULLS LAST
                ) = 1
            )
            SELECT
                $select_clause
            FROM best_channel_rows
            GROUP BY sequence,
                     charge,
                     run_file_name,
                     _modified_sequence,
                     _pg_group
            """,
            row_sql=row_sql,
            select_clause=select_clause,
        )

    @staticmethod
    def _build_additional_intensities_sql(r: dict, _has, label_expr: str = "'LFQ'") -> str:
        """Build SQL for the additional_intensities nested column."""
        candidates = [
            ("PG.MaxLFQ", "maxlfq"),
        ]
        if "normalize_intensity" in r:
            candidates.append((r["normalize_intensity"], "precursor_normalised"))
        if "ms1_area" in r:
            candidates.append((r["ms1_area"], "ms1_area"))
        entries, checks = [], []
        for col, name in candidates:
            if _has(col):
                chk = f'(r."{col}" IS NOT NULL AND NOT isnan(CAST(r."{col}" AS DOUBLE)))'
                checks.append(chk)
                entries.append(
                    f"CASE WHEN {chk} THEN STRUCT_PACK("
                    f"intensity_name := '{name}', "
                    f'intensity_value := CAST(r."{col}" AS FLOAT)) END'
                )
        if not entries:
            return "NULL AS additional_intensities"
        ai_list = ", ".join(entries)
        ai_check = " OR ".join(checks)
        return (
            f"CASE WHEN {ai_check} "
            f"THEN [STRUCT_PACK(label := {label_expr}, "
            f"intensities := LIST_FILTER([{ai_list}], x -> x IS NOT NULL))] "
            f"ELSE NULL END AS additional_intensities"
        )

    @staticmethod
    def _build_additional_scores_sql(_has) -> str:
        """Build SQL for the additional_scores nested column."""
        entries, checks = [], []
        for diann_col, score_name, higher_better in _SCORE_MAPPINGS:
            if _has(diann_col):
                hb = "true" if higher_better else "false"
                chk = f'(r."{diann_col}" IS NOT NULL AND NOT isnan(CAST(r."{diann_col}" AS DOUBLE)))'
                checks.append(chk)
                entries.append(
                    f"CASE WHEN {chk} THEN STRUCT_PACK("
                    f"score_name := '{score_name}', "
                    f'score_value := CAST(r."{diann_col}" AS DOUBLE), '
                    f"higher_better := {hb}) END"
                )
        if not entries:
            return "NULL AS additional_scores"
        sl = ", ".join(entries)
        sc = " OR ".join(checks)
        return f"CASE WHEN {sc} THEN LIST_FILTER([{sl}], x -> x IS NOT NULL) ELSE NULL END AS additional_scores"

    @staticmethod
    def _build_cv_params_sql(r: dict, _has) -> str:
        """Build SQL for the cv_params nested column."""
        entries, checks = [], []
        pqs_col = r.get("precursor_quantification_score")
        if pqs_col and _has(pqs_col):
            chk = f'(r."{pqs_col}" IS NOT NULL)'
            checks.append(chk)
            entries.append(
                f"CASE WHEN {chk} THEN STRUCT_PACK("
                f"cv_name := 'precursor_quantification_score', "
                f'cv_value := CAST(r."{pqs_col}" AS VARCHAR)) END'
            )
        for diann_col, cv_name in _CV_MAPPINGS:
            if _has(diann_col):
                chk = f'(r."{diann_col}" IS NOT NULL AND CAST(r."{diann_col}" AS VARCHAR) != \'\')'
                checks.append(chk)
                entries.append(
                    f"CASE WHEN {chk} THEN STRUCT_PACK("
                    f"cv_name := '{cv_name}', "
                    f'cv_value := CAST(r."{diann_col}" AS VARCHAR)) END'
                )
        if not entries:
            return "NULL AS cv_params"
        cl = ", ".join(entries)
        cc = " OR ".join(checks)
        return f"CASE WHEN {cc} THEN LIST_FILTER([{cl}], x -> x IS NOT NULL) ELSE NULL END AS cv_params"

    # ------------------------------------------------------------------
    # Batch processing (SQL → Arrow → write)
    # ------------------------------------------------------------------

    def _process_batch_arrow(
        self,
        runs: list[str],
        sql_template: str,
        mods_lookup: dict,
        mods_type: pa.DataType,
        pg_acc_type: pa.DataType,
        target_schema: pa.Schema,
        mzml_info_folder: Optional[str],
    ) -> Optional[pa.Table]:
        """Execute SQL batch, add Python-only columns, cast to target schema."""
        placeholders = ", ".join(["?" for _ in runs])
        sql = sql_template.replace("{run_placeholders}", placeholders)

        # Execute SQL — DuckDB pushes predicates to parquet, only reads needed data
        table = self._conn.execute(sql, runs).fetch_arrow_table()

        if len(table) == 0:
            return None

        n = len(table)

        # --- Merge scan info from mzml files if available ---
        if mzml_info_folder:
            table = self._merge_scan_info_arrow(table, mzml_info_folder)

        # --- Build modifications column from Python lookup ---
        mod_seqs = table.column("_modified_sequence").to_pylist()
        seqs = table.column("sequence").to_pylist()
        charges = table.column("charge").to_pylist()

        mods_list = [mods_lookup.get((str(ms), str(s), int(c))) for ms, s, c in zip(mod_seqs, seqs, charges)]
        mods_array = pa.array(mods_list, type=mods_type)

        # --- Build pg_accessions column ---
        pg_groups = table.column("_pg_group").to_pylist()
        pg_list = [
            [
                {"accession": acc, "start": None, "end": None, "pre": None, "post": None}
                for acc in (str(pg).split(";") if pg else [])
            ]
            or None
            for pg in pg_groups
        ]
        pg_array = pa.array(pg_list, type=pg_acc_type)

        # feature.pg_ids is left null: the feature->pg association is computed on
        # read via the Dataset softlink (bigbio/qpx#269), not materialized here.

        # --- Drop helper columns ---
        table = table.drop("_modified_sequence")
        table = table.drop("_pg_group")

        # --- Assemble final table matching target schema ---
        columns = {}
        for field in target_schema:
            if field.name == "modifications":
                columns[field.name] = mods_array
            elif field.name == "pg_accessions":
                columns[field.name] = pg_array
            elif field.name in table.schema.names:
                col = table.column(field.name)
                # Cast to target type if needed (e.g. large_string → string)
                if col.type != field.type:
                    try:
                        col = col.cast(field.type, safe=False)
                    except (pa.ArrowInvalid, pa.ArrowNotImplementedError):
                        # For complex nested types, rebuild via Python round-trip
                        col = pa.array(col.to_pylist(), type=field.type)
                columns[field.name] = col
            else:
                # Missing column → NULL array
                columns[field.name] = pa.nulls(n, type=field.type)

        return pa.table(columns, schema=target_schema)

    def _merge_scan_info_arrow(self, table: pa.Table, mzml_info_folder: str) -> pa.Table:
        """Merge MS scan info from mzml files into the Arrow table."""
        df = table.to_pandas()
        merged_parts = []

        for run_name, group_df in df.groupby("run_file_name", sort=False):
            ms_info_path = next(
                Path(mzml_info_folder).glob(f"*{run_name}_ms_info.parquet"),
                None,
            )
            if ms_info_path:
                group_df = self._merge_scan_info(group_df, ms_info_path)
            merged_parts.append(group_df)

        merged_df = pd.concat(merged_parts, ignore_index=True)

        # Recompute the mass error for rows whose observed_mz was only just
        # backfilled from the mzML. The SQL derives it before this merge runs, so
        # without this those rows keep a null error while carrying both m/z
        # (bigbio/qpx#298). Only null errors are filled; a reported value wins.
        if {"calculated_mz", "observed_mz", "mass_error_ppm"} <= set(merged_df.columns):
            calc = pd.to_numeric(merged_df["calculated_mz"], errors="coerce")
            obs = pd.to_numeric(merged_df["observed_mz"], errors="coerce")
            derivable = merged_df["mass_error_ppm"].isna() & (calc > 0) & (obs > 0)
            if derivable.any():
                merged_df.loc[derivable, "mass_error_ppm"] = ((obs[derivable] - calc[derivable]) / calc[derivable] * 1e6).astype(
                    "float32"
                )
        # Rebuild Arrow table preserving the original schema for non-scan columns
        return pa.Table.from_pandas(merged_df, schema=table.schema, preserve_index=False)

    def _merge_scan_info(self, run_data: pd.DataFrame, ms_info_path: Path) -> pd.DataFrame:
        """Fill missing scan metadata from one run's MS-info table."""
        target = pd.read_parquet(ms_info_path, columns=["rt", "scan", "precursor_mz"])
        run_data = run_data.copy()
        run_data["_merge_row"] = range(len(run_data))
        run_data["rt"] = run_data["rt"].astype("float64")
        run_data["observed_mz"] = pd.to_numeric(run_data["observed_mz"], errors="coerce").astype("float64")
        target["scan"] = pd.to_numeric(target["scan"], errors="coerce")
        target["precursor_mz"] = pd.to_numeric(target["precursor_mz"], errors="coerce")
        target["rt"] = target["rt"].astype("float64")

        def _first_scan(value) -> int | None:
            if hasattr(value, "tolist"):
                value = value.tolist()
            if isinstance(value, (list, tuple)):
                value = value[0] if value else None
            try:
                return int(value) if value is not None else None
            except (TypeError, ValueError):
                return None

        report_scans = run_data["scan"].map(_first_scan)
        scan_to_mz = target.dropna(subset=["scan"]).drop_duplicates("scan").set_index("scan")["precursor_mz"]
        exact_mz = report_scans.map(scan_to_mz)
        observed_mz = run_data["observed_mz"]
        fill_exact_mz = (observed_mz.isna() | (observed_mz <= 0)) & exact_mz.notna()
        run_data.loc[fill_exact_mz, "observed_mz"] = exact_mz[fill_exact_mz]

        # Older DIA-NN reports may not carry MS2.Scan.  For only those rows,
        # select the nearest MS2 entry by RT, and never extend past the measured
        # acquisition range.  Existing report scans remain authoritative.
        missing_scan = report_scans.isna() & run_data["rt"].notna()
        rt_target = target.dropna(subset=["rt", "scan", "precursor_mz"]).sort_values("rt")
        if missing_scan.any() and not rt_target.empty:
            left = run_data.loc[missing_scan, ["_merge_row", "rt"]].sort_values("rt")
            right = rt_target[["rt", "scan", "precursor_mz"]].rename(
                columns={"rt": "_ms_rt", "scan": "_matched_scan", "precursor_mz": "_matched_mz"}
            )
            nearest = pd.merge_asof(left, right, left_on="rt", right_on="_ms_rt", direction="nearest")
            outside_acquisition = (nearest["rt"] < right["_ms_rt"].min()) | (nearest["rt"] > right["_ms_rt"].max())
            nearest.loc[outside_acquisition, ["_matched_scan", "_matched_mz"]] = None

            scan_matches = nearest.set_index("_merge_row")["_matched_scan"]
            mz_matches = nearest.set_index("_merge_row")["_matched_mz"]
            matched_scan = run_data["_merge_row"].map(scan_matches)
            matched_mz = run_data["_merge_row"].map(mz_matches)
            has_match = matched_scan.notna()
            run_data.loc[has_match, "scan"] = matched_scan[has_match].map(lambda value: [int(value)])

            observed_mz = run_data["observed_mz"]
            fill_nearest_mz = has_match & (observed_mz.isna() | (observed_mz <= 0)) & matched_mz.notna()
            run_data.loc[fill_nearest_mz, "observed_mz"] = matched_mz[fill_nearest_mz]

        return run_data.drop(columns="_merge_row")
