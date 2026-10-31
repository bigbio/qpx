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
    ("PTM.Site.Confidence", "diann_ptm_site_confidence", True),
    ("Lib.PTM.Site.Confidence", "diann_lib_ptm_site_confidence", True),
    ("Peptidoform.Q.Value", "diann_peptidoform_qvalue", False),
    ("Global.Peptidoform.Q.Value", "diann_global_peptidoform_qvalue", False),
    ("Lib.Peptidoform.Q.Value", "diann_lib_peptidoform_qvalue", False),
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

# Ion mobility columns, 0 in a report of data acquired without ion mobility.
_ION_MOBILITY_CV = {"Predicted.IM", "Predicted.iIM", "iIM"}


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
    WHERE r.$run_col IN ({run_placeholders})
      $qvalue_filter
      AND r.$pg_col IS NOT NULL
    """


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
                ORDER BY _qvalue NULLS LAST, rt NULLS LAST, intensities[1].intensity NULLS LAST
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


def _positive_sql(column: str) -> str:
    """Whether a value is finite and above zero: DIA-NN writes 0 for a quantity or ion mobility it has not measured."""
    return f'(isfinite(CAST(r."{column}" AS DOUBLE)) AND CAST(r."{column}" AS DOUBLE) > 0)'


def _measured_float_sql(column: str) -> str:
    """Build a SQL expression that converts a measured (positive) value to FLOAT, null for DIA-NN's 0."""
    return f'CASE WHEN {_positive_sql(column)} THEN CAST(r."{column}" AS FLOAT) END'


def _build_mass_select_parts(resolved: dict, report_cols: set[str]) -> list[str]:
    """Build observed m/z and mass error expressions."""
    has_column = report_cols.__contains__
    parts: list[str] = []
    # observed_mz: DIA-NN's Precursor.Mz is the theoretical m/z of the precursor;
    # the measured one is Precursor.Mz plus Ms1.Apex.Mz.Delta, the offset of the
    # MS1 apex. A zero offset marks a precursor without a measured MS1 apex (all
    # of them where Ms1.Apex.Area is 0), and a report without the offset (DIA-NN
    # 1.8, 1.9) has no measured m/z, so both stay NULL rather than repeating the
    # theoretical value (bigbio/qpx#244).
    precursor_mz_col = resolved.get("precursor_mz")
    delta_col = resolved.get("ms1_apex_mz_delta")
    observed_expr = None
    if precursor_mz_col and delta_col and has_column(precursor_mz_col) and has_column(delta_col):
        delta = f'CAST(r."{delta_col}" AS DOUBLE)'
        observed_expr = (
            f'CASE WHEN isfinite({delta}) AND {delta} <> 0 AND isfinite(CAST(r."{precursor_mz_col}" AS DOUBLE)) '
            f'THEN CAST(r."{precursor_mz_col}" AS DOUBLE) + {delta} END'
        )
    parts.append(f"CAST({observed_expr} AS FLOAT) AS observed_mz" if observed_expr else "NULL::FLOAT AS observed_mz")

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
    elif observed_expr:
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

    # DIA-NN writes ion mobility 0 for data acquired without it.
    ion_mobility_col = resolved.get("ion_mobility")
    parts.append(
        f"{_measured_float_sql(ion_mobility_col)} AS ion_mobility"
        if ion_mobility_col and has_column(ion_mobility_col)
        else "NULL::FLOAT AS ion_mobility"
    )
    return parts


# The proteins each target peptidoform maps to. A DIA-NN decoy precursor has
# its target's sequence, so decoy rows are left out.
_PEPMAP_SQL = """
    SELECT lk.peptidoform,
           any_value(lk.sequence) AS sequence,
           list_sort(list_distinct(list_filter(
               list_transform(flatten(list(string_split(r.$ids_col, ';'))), lambda x: trim(x)),
               lambda x: x <> ''
           ))) AS accessions
    FROM report r
    JOIN precursor_lookup lk
        ON r.$mod_col = lk.modified_sequence
        AND r.$seq_col = lk.sequence
        AND r.$chg_col = lk.charge
    WHERE r.$run_col IN ({run_placeholders})
      $qvalue_filter
      $target_filter
      AND r.$pg_col IS NOT NULL
    GROUP BY lk.peptidoform
    ORDER BY lk.peptidoform
    """


def build_pepmap_sql(resolved: dict, report_cols: set[str], qvalue_threshold: Optional[float]) -> str | None:
    """The pepmap query over the feature rows: each peptidoform with its ``Protein.Ids``; None without that column."""
    ids_col = resolved.get("mp_accessions")
    if ids_col not in report_cols:
        return None
    decoy_col = resolved.get("decoy")
    return sql_build(
        _PEPMAP_SQL,
        ids_col=validate_identifier(ids_col),
        mod_col=validate_identifier(resolved["modified_sequence"]),
        seq_col=validate_identifier(resolved["sequence"]),
        chg_col=validate_identifier(resolved["charge"]),
        run_col=validate_identifier(resolved["run_file_name"]),
        qvalue_filter=_qvalue_filter_sql(qvalue_threshold, resolved.get("precursor_qvalue")),
        target_filter=f'AND CAST(r."{decoy_col}" AS BOOLEAN) IS NOT TRUE' if decoy_col in report_cols else "",
        pg_col=validate_identifier(resolved["pg_accessions"]),
    )


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

    # unique: a one-protein group of a peptide that maps to no other protein.
    # Protein.Ids lists every protein the peptide maps to; DIA-NN's Proteotypic
    # flag depends on the protein-grouping level (genes by default), so it stays
    # a cv_param.
    members = f"list_filter(string_split(r.\"{pg_col}\", ';'), lambda x: trim(x) <> '')"
    single = f"len({members}) = 1"
    ids_col = resolved.get("mp_accessions")
    if ids_col and has_column(ids_col):
        single += f" AND len(list_filter(string_split(r.\"{ids_col}\", ';'), lambda x: trim(x) <> '')) = 1"
    parts.append(f'CASE WHEN len({members}) > 0 THEN {single} END AS "unique"')

    pg_qvalue_col = resolved.get("pg_global_qvalue")
    parts.append(
        f"{_safe_double_sql(pg_qvalue_col)} AS pg_global_qvalue"
        if pg_qvalue_col and has_column(pg_qvalue_col)
        else "NULL::DOUBLE AS pg_global_qvalue"
    )

    for field in ("ion_mobility_start", "ion_mobility_stop"):
        column = resolved.get(field)
        parts.append(f"{_measured_float_sql(column)} AS {field}" if column and has_column(column) else f"NULL::FLOAT AS {field}")

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
    int_col = r["intensity"]
    # A precursor DIA-NN reports without a quantity (0) keeps its label with a
    # null intensity.
    parts.append(f"[STRUCT_PACK(label := {label_expr}, intensity := {_measured_float_sql(int_col)})] AS intensities")

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
    occupancy = 'r."Site.Occupancy.Probabilities"' if has_column("Site.Occupancy.Probabilities") else "NULL"
    parts.append(f"CAST({occupancy} AS VARCHAR) AS _site_occupancy")
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
        run_col=validate_identifier(run_col),
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
    helpers = [qvalue_expr, "cl.channel_mass AS _channel_mass"]
    channel_join = sql_build(
        """
        JOIN diann_channel_labels cl
          ON CAST(r.$channel_col AS VARCHAR) = cl.channel_value
        """,
        channel_col=validate_identifier(channel_col),
    )
    return "cl.canonical_label", helpers, channel_join


# Fields describing the multiplet's peptidoform, which is the lightest partner's,
# and so taken from the lightest channel's row.
_LIGHTEST_CHANNEL_FIELDS = {"observed_mz", "mass_error_ppm"}


def _aggregate_channel_rows_sql(row_sql: str, target_schema: pa.Schema) -> str:
    """Collapse DIA-NN's one-row-per-channel output to one QPX feature row.

    The other scalar fields come from the channel with the best precursor
    q-value, ties broken by label, so a conversion is reproducible.
    """
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
        "modification_candidates",
        # Derived by the writer from additional_scores and pg_accessions.
        "precursor_qvalue",
        "global_qvalue",
        "pg_group_id",
        # Producer-specific optional fields not emitted by DIA-NN.
        "identified_labels",
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
        elif name in _LIGHTEST_CHANNEL_FIELDS:
            select_parts.append(f"FIRST({quoted} ORDER BY _channel_mass NULLS LAST, intensities[1].label) AS {quoted}")
        else:
            select_parts.append(f"FIRST({quoted} ORDER BY _qvalue NULLS LAST, intensities[1].label) AS {quoted}")

    select_parts.extend(
        ["_modified_sequence", "_pg_group", "FIRST(_site_occupancy ORDER BY _qvalue NULLS LAST) AS _site_occupancy"]
    )
    select_clause = ",\n            ".join(select_parts)
    return sql_build(
        _CHANNEL_ROWS_SQL,
        row_sql=row_sql,
        select_clause=select_clause,
    )


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
            chk = _positive_sql(col)
            checks.append(chk)
            entries.append(
                f"CASE WHEN {chk} THEN STRUCT_PACK(intensity_name := '{name}', intensity_value := CAST(r.\"{col}\" AS FLOAT)) END"
            )
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
            if diann_col in _ION_MOBILITY_CV:
                chk = _positive_sql(diann_col)
            checks.append(chk)
            entries.append(
                f"CASE WHEN {chk} THEN STRUCT_PACK(cv_name := '{cv_name}', cv_value := CAST(r.\"{diann_col}\" AS VARCHAR)) END"
            )
    if not entries:
        return "NULL AS cv_params"
    cl = ", ".join(entries)
    cc = " OR ".join(checks)
    return f"CASE WHEN {cc} THEN LIST_FILTER([{cl}], lambda x: x IS NOT NULL) ELSE NULL END AS cv_params"
