"""SQL of the fragment view from a DIA-NN report.

DIA-NN 1.x lists every library fragment in the ``Fragment.Quant.*`` cells; 2.x, with
``--export-quant``, writes one ``Fr.N.Index/Id/Quantity/Score`` group of columns per fragment.
"""

from __future__ import annotations

import re

from qpx.converters.diann.base_adapter import run_name_sql
from qpx.converters.diann.feature_sql import _build_channel_sql_parts, _qvalue_filter_sql, decoy_flag_sql
from qpx.core.sql import sql_build, validate_identifier

# The report rows of the precursors a fragment view keeps.
_PRECURSORS_SQL = """
        FROM report r
        JOIN precursor_lookup lk
          ON r.$modified = lk.modified_sequence AND r.$sequence = lk.sequence AND r.$charge = lk.charge
        $channel_join
        WHERE r.$pg_col IS NOT NULL $qvalue_filter"""

_LIST_SQL = """
    SELECT run_file_name, peptidoform, charge, is_decoy, is_library_decoy, label, fragment_index, $annotation,
           NULLIF(TRY_CAST(raw AS FLOAT), 0) AS intensity,
           NULLIF(TRY_CAST(corrected AS FLOAT), 0) AS corrected_intensity,
           TRY_CAST(correlations AS FLOAT) AS correlation
    FROM (
        SELECT $run AS run_file_name, lk.peptidoform, CAST(r.$charge AS SMALLINT) AS charge,
               $decoy AS is_decoy, $library_decoy AS is_library_decoy, $label AS label,
               generate_subscripts($raw_list, 1) AS fragment_index,
               unnest($raw_list) AS raw, unnest($corrected_list) AS corrected,
               unnest($correlation_list) AS correlations, unnest($info_list) AS info
        $precursors
    )
    WHERE NULLIF(TRY_CAST(raw AS FLOAT), 0) IS NOT NULL OR NULLIF(TRY_CAST(corrected AS FLOAT), 0) IS NOT NULL
"""

# Fr.N.Index is the 0-based position in the library's fragment list, -1 for an unused slot.
# Fr.N.Quantity is the fragment's peak height, not normalised, and Fr.N.Score the quality of
# its extracted elution profile (DIA-NN README, --export-quant).
_SLOT_SQL = """
    SELECT run_file_name, peptidoform, charge, is_decoy, is_library_decoy, label, fr.fr_index + 1 AS fragment_index,
           $annotation,
           fr.height AS peak_height,
           [STRUCT_PACK(score_name := 'diann_fragment_score', score_value := fr.score, higher_better := true)]
               AS additional_scores
    FROM (
        SELECT $run AS run_file_name, lk.peptidoform, CAST(r.$charge AS SMALLINT) AS charge,
               $decoy AS is_decoy, $library_decoy AS is_library_decoy, $label AS label,
               unnest([$slots]) AS fr
        $precursors
    )
    WHERE fr.height > 0
"""

_SLOT_ID = re.compile(r"Fr\.(\d+)\.Id")


def _annotation_sql(info: str) -> str:
    """Ion, charge, neutral loss and theoretical m/z of an annotation such as ``y7-H2O^1/826.4``."""
    return sql_build(
        """regexp_extract($info, '^([a-z]+[0-9]+)', 1) AS ion_type,
           TRY_CAST(regexp_extract($info, '\\^([0-9]+)', 1) AS INTEGER) AS fragment_charge,
           NULLIF(NULLIF(regexp_extract($info, '^[a-z]+[0-9]+-([^^]+)', 1), ''), 'unknown') AS neutral_loss,
           TRY_CAST(regexp_extract($info, '/([0-9.]+)', 1) AS DOUBLE) AS mz""",
        info=info,
    )


def _fragment_list_sql(column: str, report_cols: set[str]) -> str:
    """A ``;``-terminated report cell as a list, the same length as the raw list when absent."""
    if column in report_cols:
        return f"string_split(rtrim(CAST(r.{validate_identifier(column)} AS VARCHAR), ';'), ';')"
    return "list_transform(string_split(rtrim(CAST(r.\"Fragment.Quant.Raw\" AS VARCHAR), ';'), ';'), lambda x: NULL::VARCHAR)"


def _slots_sql(report_cols: set[str]) -> str:
    """One struct per ``Fr.N`` group of columns of a 2.x report, empty when it has none."""
    slots = sorted(int(match.group(1)) for column in report_cols if (match := _SLOT_ID.fullmatch(column)))
    return ", ".join(
        f"STRUCT_PACK(fr_index := CAST(r.{validate_identifier(f'Fr.{n}.Index')} AS INTEGER), "
        f"id := CAST(r.{validate_identifier(f'Fr.{n}.Id')} AS VARCHAR), "
        f"height := CAST(r.{validate_identifier(f'Fr.{n}.Quantity')} AS FLOAT), "
        f"score := CAST(r.{validate_identifier(f'Fr.{n}.Score')} AS DOUBLE))"
        for n in slots
    )


def build_fragment_sql(
    resolved: dict, report_cols: set[str], channel_col: str | None, qvalue_threshold: float | None
) -> str | None:
    """SQL of the fragment view's rows, or None for a report without fragment quantities."""
    slots = _slots_sql(report_cols)
    if "Fragment.Quant.Raw" not in report_cols and not slots:
        return None
    qvalue_col = resolved.get("precursor_qvalue") if resolved.get("precursor_qvalue") in report_cols else None
    label, _, channel_join = _build_channel_sql_parts(channel_col, qvalue_col, report_cols.__contains__)
    charge = validate_identifier(resolved["charge"])
    modified = validate_identifier(resolved["modified_sequence"])
    precursors = sql_build(
        _PRECURSORS_SQL,
        charge=charge,
        modified=modified,
        sequence=validate_identifier(resolved["sequence"]),
        channel_join=channel_join,
        pg_col=validate_identifier(resolved["pg_accessions"]),
        qvalue_filter=_qvalue_filter_sql(qvalue_threshold, qvalue_col),
    )
    common = {
        "run": run_name_sql("r." + validate_identifier(resolved["run_file_name"])),
        "charge": charge,
        "label": label,
        # Both flags are part of the feature identity (feature_sql).
        "decoy": decoy_flag_sql(resolved, resolved.get("decoy") in report_cols),
        "library_decoy": f"starts_with(r.{modified}, '^')",
        "precursors": precursors,
    }
    if "Fragment.Quant.Raw" not in report_cols:
        return sql_build(_SLOT_SQL, annotation=_annotation_sql("fr.id"), slots=slots, **common)
    return sql_build(
        _LIST_SQL,
        annotation=_annotation_sql("info"),
        raw_list=_fragment_list_sql("Fragment.Quant.Raw", report_cols),
        corrected_list=_fragment_list_sql("Fragment.Quant.Corrected", report_cols),
        correlation_list=_fragment_list_sql("Fragment.Correlations", report_cols),
        info_list=_fragment_list_sql("Fragment.Info", report_cols),
        **common,
    )
