"""SQL of the fragment view from a DIA-NN report's ``Fragment.Quant.*`` lists."""

from __future__ import annotations

from qpx.converters.diann.base_adapter import run_name_sql
from qpx.converters.diann.feature_sql import _build_channel_sql_parts, _qvalue_filter_sql
from qpx.core.sql import sql_build, validate_identifier

_FRAGMENT_SQL = """
    SELECT run_file_name, peptidoform, charge, label, fragment_index,
           regexp_extract(info, '^([a-z]+[0-9]+)', 1) AS ion_type,
           TRY_CAST(regexp_extract(info, '\\^([0-9]+)', 1) AS INTEGER) AS fragment_charge,
           NULLIF(NULLIF(regexp_extract(info, '^[a-z]+[0-9]+-([^^]+)', 1), ''), 'unknown') AS neutral_loss,
           TRY_CAST(regexp_extract(info, '/([0-9.]+)', 1) AS DOUBLE) AS mz,
           NULLIF(TRY_CAST(raw AS FLOAT), 0) AS intensity,
           NULLIF(TRY_CAST(corrected AS FLOAT), 0) AS corrected_intensity,
           TRY_CAST(correlations AS FLOAT) AS correlation
    FROM (
        SELECT $run AS run_file_name, lk.peptidoform, CAST(r.$charge AS SMALLINT) AS charge, $label AS label,
               generate_subscripts($raw_list, 1) AS fragment_index,
               unnest($raw_list) AS raw, unnest($corrected_list) AS corrected,
               unnest($correlation_list) AS correlations, unnest($info_list) AS info
        FROM report r
        JOIN precursor_lookup lk
          ON r.$modified = lk.modified_sequence AND r.$sequence = lk.sequence AND r.$charge = lk.charge
        $channel_join
        WHERE r.$pg_col IS NOT NULL $qvalue_filter
    )
    WHERE NULLIF(TRY_CAST(raw AS FLOAT), 0) IS NOT NULL OR NULLIF(TRY_CAST(corrected AS FLOAT), 0) IS NOT NULL
"""


def _fragment_list_sql(column: str, report_cols: set[str]) -> str:
    """A ``;``-terminated report cell as a list, the same length as the raw list when absent."""
    if column in report_cols:
        return f"string_split(rtrim(CAST(r.{validate_identifier(column)} AS VARCHAR), ';'), ';')"
    return "list_transform(string_split(rtrim(CAST(r.\"Fragment.Quant.Raw\" AS VARCHAR), ';'), ';'), lambda x: NULL::VARCHAR)"


def build_fragment_sql(
    resolved: dict, report_cols: set[str], channel_col: str | None, qvalue_threshold: float | None
) -> str | None:
    """SQL of the fragment view's rows, or None for a report without fragment quantities."""
    if "Fragment.Quant.Raw" not in report_cols:
        return None
    qvalue_col = resolved.get("qvalue") if resolved.get("qvalue") in report_cols else None
    label, _, channel_join = _build_channel_sql_parts(channel_col, qvalue_col, report_cols.__contains__)
    return sql_build(
        _FRAGMENT_SQL,
        run=run_name_sql("r." + validate_identifier(resolved["run_file_name"])),
        charge=validate_identifier(resolved["charge"]),
        modified=validate_identifier(resolved["modified_sequence"]),
        sequence=validate_identifier(resolved["sequence"]),
        label=label,
        channel_join=channel_join,
        raw_list=_fragment_list_sql("Fragment.Quant.Raw", report_cols),
        corrected_list=_fragment_list_sql("Fragment.Quant.Corrected", report_cols),
        correlation_list=_fragment_list_sql("Fragment.Correlations", report_cols),
        info_list=_fragment_list_sql("Fragment.Info", report_cols),
        pg_col=validate_identifier(resolved["pg_accessions"]),
        qvalue_filter=_qvalue_filter_sql(qvalue_threshold, qvalue_col),
    )
