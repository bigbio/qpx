"""SQL expressions and batch queries for DIA-NN features."""

from __future__ import annotations

from typing import Optional

import pyarrow as pa

from qpx.converters.diann.base_adapter import decoy_protein_group_sql, run_name_sql
from qpx.core.sql import sql_build, validate_identifier

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

# Ion mobility fields and cv columns. For data acquired without ion mobility
# separation, DIA-NN reports IM (IM.Start / IM.Stop), Predicted.IM and
# Predicted.iIM as 0; DIA-NN 1.8 also writes iIM = 0, while DIA-NN 2.x reports a
# library iIM. A 1/K0 or drift time is never 0, so each value is checked on its
# own and a 0 is written as null (or left out of cv_params).
_ION_MOBILITY_FIELDS = ("ion_mobility", "ion_mobility_start", "ion_mobility_stop")
_ION_MOBILITY_CV_COLUMNS = frozenset({"Predicted.IM", "Predicted.iIM", "iIM"})


_RT_SECONDS_FACTOR = 60.0

_BATCH_SQL = """
    SELECT
        $select_clause
    FROM report r
    JOIN precursor_lookup lk
        ON r.$mod_col = lk.modified_sequence
        AND r.$seq_col = lk.sequence
        AND r.$chg_col = lk.charge
    $channel_join
    WHERE $run_name IN ({run_placeholders})
      $qvalue_filter
      AND r.$pg_col IS NOT NULL
    """


# Per-channel nested columns collected into one feature row, and the channel
# rows each one keeps: an unquantified channel has a null intensity and no entry.
# Every channel row keeps its labelled entry so that the deduplication below can
# partition on intensities[1].label. With a plain label column there, DuckDB
# rewrites ROW_NUMBER() = 1 into an arg_min aggregate over the whole wide row,
# which made the plexDIA feature step about 20 % slower.
_CHANNEL_LIST_FILTERS = {
    "intensities": "intensities[1].intensity IS NOT NULL",
    "additional_intensities": "additional_intensities IS NOT NULL",
}

_CHANNEL_ROWS_SQL = """
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
        """


def _safe_float_sql(column: str) -> str:
    """Build a SQL expression that converts a finite value to FLOAT."""
    return f'CASE WHEN r."{column}" IS NOT NULL AND NOT isnan(CAST(r."{column}" AS DOUBLE)) THEN CAST(r."{column}" AS FLOAT) END'


def _safe_double_sql(column: str) -> str:
    """Build a SQL expression that converts a finite value to DOUBLE."""
    return f'CASE WHEN r."{column}" IS NOT NULL AND NOT isnan(CAST(r."{column}" AS DOUBLE)) THEN CAST(r."{column}" AS DOUBLE) END'


def _quantity_sql(column: str) -> str:
    """Build a SQL expression that keeps a positive quantity as FLOAT, NULL otherwise.

    DIA-NN writes 0 for a quantity it did not measure (an unquantified channel's
    ``Precursor.Quantity``, an ``Ms1.Area`` without MS1 signal, QuantUMS's
    ``PG.MaxLFQ``); a tool's zero sentinel is not a measured zero. DuckDB orders
    NaN above every number, so NaN is excluded explicitly.
    """
    value = f'CAST(r."{column}" AS DOUBLE)'
    return f'CASE WHEN {value} > 0 AND NOT isnan({value}) THEN CAST(r."{column}" AS FLOAT) END'


def _build_mass_select_parts(resolved: dict, report_cols: set[str]) -> list[str]:
    """Build observed m/z and mass error expressions."""
    has_column = report_cols.__contains__
    parts: list[str] = []
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
    # null when the peptidoform's mass cannot be computed, and observed_mz is null
    # when the report has no Precursor.Mz, so an unmeasured row stays null rather
    # than reporting a fabricated error. The formula matches docs/spec/feature.md.
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

    return parts


def _build_acquisition_select_parts(resolved: dict, report_cols: set[str]) -> list[str]:
    """Build run, scan, retention time, and ion mobility expressions."""
    has_column = report_cols.__contains__
    parts: list[str] = []
    predicted_rt_col = resolved.get("predicted_rt")
    parts.append(
        f"CAST({_safe_float_sql(predicted_rt_col)} * {_RT_SECONDS_FACTOR} AS FLOAT) AS predicted_rt"
        if predicted_rt_col and has_column(predicted_rt_col)
        else "NULL::FLOAT AS predicted_rt"
    )

    run_col = "r." + validate_identifier(resolved["run_file_name"])
    parts.append(f"{run_name_sql(run_col)} AS run_file_name")

    # MS2.Scan (DIA-NN 1.8) is a 0-based position among the run's MS2 spectra,
    # not a native scan number: it becomes a scan only through the MS-info table.
    parts.append("[]::INTEGER[] AS scan")
    ms2_col = resolved.get("ms2_scan")
    parts.append(
        f'CAST(r."{ms2_col}" AS INTEGER) AS _ms2_index' if ms2_col and has_column(ms2_col) else "NULL::INTEGER AS _ms2_index"
    )

    rt_col = resolved.get("rt")
    parts.append(
        f"CAST({_safe_float_sql(rt_col)} * {_RT_SECONDS_FACTOR} AS FLOAT) AS rt"
        if rt_col and has_column(rt_col)
        else "NULL::FLOAT AS rt"
    )

    parts.extend(_build_ion_mobility_select_parts(resolved, report_cols))
    return parts


def _build_ion_mobility_select_parts(resolved: dict, report_cols: set[str]) -> list[str]:
    """Build ion mobility expressions: NULL when the column is absent or holds DIA-NN's 0 (no ion mobility data)."""
    parts: list[str] = []
    for field in _ION_MOBILITY_FIELDS:
        column = resolved.get(field)
        parts.append(
            f"NULLIF({_safe_float_sql(column)}, 0) AS {field}" if column and column in report_cols else f"NULL::FLOAT AS {field}"
        )
    return parts


def _qvalue_filter_sql(qvalue_threshold: Optional[float], qv_col: str | None) -> str:
    """Build the opt-in precursor q-value filter."""
    qvalue_filter = ""
    if qvalue_threshold is not None:
        qvalue_filter = sql_build(
            "AND r.$qv_col < $qv_threshold",
            qv_col=validate_identifier(qv_col),
            qv_threshold=str(float(qvalue_threshold)),
        )

    return qvalue_filter


def _build_core_select_parts(
    resolved: dict,
    report_cols: set[str],
    has_decoy_col: bool,
    pg_group: str,
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

    parts.append("CAST(lk.calculated_mz AS FLOAT) AS calculated_mz")

    parts.extend(_build_mass_select_parts(resolved, report_cols))

    parts.extend(_build_acquisition_select_parts(resolved, report_cols))
    parts.extend(
        [
            "CAST(lk.missed_cleavages AS SMALLINT) AS missed_cleavages",
            # Blank/whitespace Protein.Group → NULL (e.g. DIA-NN unmapped rows).
            f"NULLIF(TRIM(SPLIT_PART({pg_group}, ';', 1)), '') AS anchor_protein",
        ]
    )
    return parts


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

    run_col = "r." + validate_identifier(resolved["run_file_name"])
    parts.append(f"{run_name_sql(run_col)} AS id_run_file_name")

    for field in ("rt_start", "rt_stop"):
        column = resolved.get(field)
        parts.append(
            f"CAST({_safe_float_sql(column)} * {_RT_SECONDS_FACTOR} AS FLOAT) AS {field}"
            if column and has_column(column)
            else f"NULL::FLOAT AS {field}"
        )
    return parts


def _build_batch_sql(
    resolved: dict,
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
    r = resolved
    run_col = r["run_file_name"]
    qv_col = r.get("precursor_qvalue")
    if qvalue_threshold is not None and qv_col is None:
        raise ValueError("qvalue_threshold requires the DIA-NN report column Q.Value")
    pg_col = r["pg_accessions"]
    # Decoy rows get the DECOY_ accessions of their pg decoy group.
    pg_group = f'r."{pg_col}"'
    if has_decoy_col:
        pg_group = decoy_protein_group_sql(pg_group, f'r."{r["decoy"]}"')

    has_column = report_cols.__contains__
    parts = _build_core_select_parts(r, report_cols, has_decoy_col, pg_group)
    parts.extend(_build_metadata_select_parts(r, report_cols))

    # --- Nested: intensities ---
    label_expr, channel_helper_parts, channel_join = _build_channel_sql_parts(channel_col, qv_col)
    parts.append(_build_intensities_sql(r["intensity"], label_expr, bool(channel_col)))

    # --- Nested columns built by helpers ---
    parts.append(_build_additional_intensities_sql(r, has_column, label_expr))
    parts.append(_build_additional_scores_sql(has_column))
    parts.append(_build_cv_params_sql(r, has_column))

    # --- Helper columns (for Python post-processing, dropped later) ---
    mod_col = r["modified_sequence"]
    seq_col = r["sequence"]
    chg_col = r["charge"]
    parts.append(f'r."{mod_col}" AS _modified_sequence')
    parts.append(f"{pg_group} AS _pg_group")
    parts.extend(channel_helper_parts)

    # --- Build full query ---
    # Precursor Q.Value filter is opt-in (bigbio/qpx#241). With no threshold
    # the report is converted as reported (DIA-NN already FDR-filters it);
    # only apply the cutoff when a threshold is explicitly provided.
    qvalue_filter = _qvalue_filter_sql(qvalue_threshold, qv_col)

    row_sql = sql_build(
        _BATCH_SQL,
        select_clause=",\n        ".join(parts),
        mod_col=validate_identifier(mod_col),
        seq_col=validate_identifier(seq_col),
        chg_col=validate_identifier(chg_col),
        channel_join=channel_join,
        run_name=run_name_sql("r." + validate_identifier(run_col)),
        qvalue_filter=qvalue_filter,
        pg_col=validate_identifier(pg_col),
    )

    if channel_col:
        return _aggregate_channel_rows_sql(row_sql, target_schema)
    return row_sql


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
        "modification_candidates",
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
        elif name in _CHANNEL_LIST_FILTERS:
            select_parts.append(
                f"LIST({quoted}[1] ORDER BY intensities[1].label) FILTER (WHERE {_CHANNEL_LIST_FILTERS[name]}) AS {quoted}"
            )
        else:
            select_parts.append(f"FIRST({quoted} ORDER BY _qvalue NULLS LAST) AS {quoted}")

    select_parts.extend(["_modified_sequence", "_pg_group"])
    select_clause = ",\n            ".join(select_parts)
    return sql_build(
        _CHANNEL_ROWS_SQL,
        row_sql=row_sql,
        select_clause=select_clause,
    )


def _build_intensities_sql(intensity_col: str, label_expr: str, channel_rows: bool) -> str:
    """Build SQL for the intensities nested column: no entry for an unquantified precursor.

    The intensity struct cannot hold a null, so a feature DIA-NN identified
    without a quantity keeps its row with ``intensities = null``. A channel row
    keeps its labelled entry, because the per-channel deduplication partitions on
    ``intensities[1].label``; ``_aggregate_channel_rows_sql`` then drops the
    entries without a quantity.
    """
    quantity = _quantity_sql(intensity_col)
    entry = f"[STRUCT_PACK(label := {label_expr}, intensity := {quantity})]"
    if channel_rows:
        return f"{entry} AS intensities"
    return f"CASE WHEN {quantity} IS NOT NULL THEN {entry} END AS intensities"


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
            quantity = _quantity_sql(col)
            chk = f"({quantity} IS NOT NULL)"
            checks.append(chk)
            entries.append(f"CASE WHEN {chk} THEN STRUCT_PACK(intensity_name := '{name}', intensity_value := {quantity}) END")
    if not entries:
        return "NULL AS additional_intensities"
    ai_list = ", ".join(entries)
    ai_check = " OR ".join(checks)
    return (
        f"CASE WHEN {ai_check} "
        f"THEN [STRUCT_PACK(label := {label_expr}, "
        f"intensities := LIST_FILTER([{ai_list}], lambda x: x IS NOT NULL))] "
        f"ELSE NULL END AS additional_intensities"
    )


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
    return f"CASE WHEN {sc} THEN LIST_FILTER([{sl}], lambda x: x IS NOT NULL) ELSE NULL END AS additional_scores"


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
            if diann_col in _ION_MOBILITY_CV_COLUMNS:
                chk = f'({chk} AND TRY_CAST(r."{diann_col}" AS DOUBLE) IS DISTINCT FROM 0)'
            checks.append(chk)
            entries.append(
                f"CASE WHEN {chk} THEN STRUCT_PACK(cv_name := '{cv_name}', cv_value := CAST(r.\"{diann_col}\" AS VARCHAR)) END"
            )
    if not entries:
        return "NULL AS cv_params"
    cl = ", ".join(entries)
    cc = " OR ".join(checks)
    return f"CASE WHEN {cc} THEN LIST_FILTER([{cl}], lambda x: x IS NOT NULL) ELSE NULL END AS cv_params"
