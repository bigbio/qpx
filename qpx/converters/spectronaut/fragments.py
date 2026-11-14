"""The fragment view of a Spectronaut report, which already has one row per fragment."""

from __future__ import annotations

from qpx.core.sql import validate_identifier

_REQUIRED = {"F.FrgIon", "F.PeakArea"}
_COLUMNS = (
    ("F.TheoreticalMz", "mz", "DOUBLE"),
    ("F.MeasuredMz", "observed_mz", "DOUBLE"),
    ("F.PeakArea", "intensity", "FLOAT"),
    ("F.NormalizedPeakArea", "corrected_intensity", "FLOAT"),
    ("F.PeakHeight", "peak_height", "FLOAT"),
)
_SCORES = (
    ("F.InterferenceScore", "spectronaut_interference_score", False),
    ("F.Log10SignalToNoise", "spectronaut_log10_signal_to_noise", True),
    ("F.MeasuredRelativeIntensity", "spectronaut_measured_relative_intensity", True),
    ("F.PredictedRelativeIntensity", "spectronaut_predicted_relative_intensity", True),
    ("F.CalibratedMassAccuracy (PPM)", "spectronaut_calibrated_mass_accuracy_ppm", False),
)


def _measured_sql(column: str, kind: str) -> str:
    """A finite value above zero, else null: NaN stands for not measured."""
    value = f"TRY_CAST(r.{validate_identifier(column)} AS {kind})"
    return f"CASE WHEN isfinite({value}) AND {value} > 0 THEN {value} END"


def _optional_cast(column: str, kind: str, report_cols: set[str]) -> str:
    return f"TRY_CAST(r.{validate_identifier(column)} AS {kind})" if column in report_cols else f"NULL::{kind}"


def _scores_sql(report_cols: set[str]) -> str:
    entries = []
    for column, name, higher in _SCORES:
        if column in report_cols:
            value = f'TRY_CAST(r."{column}" AS DOUBLE)'
            entries.append(
                f"CASE WHEN isfinite({value}) THEN "
                f"STRUCT_PACK(score_name := '{name}', score_value := {value}, higher_better := {str(higher).lower()}) END"
            )
    return f"list_filter([{', '.join(entries)}], lambda x: x IS NOT NULL)" if entries else "NULL"


def fragment_sql(resolved: dict, report_cols: set[str], labelled: bool) -> str | None:
    """SQL of the fragment view's rows, or None for a report without fragments.

    A labelled report has its channels in columns of the precursor, so its fragments have no label.
    """
    if not _REQUIRED <= report_cols:
        return None
    run = validate_identifier(resolved["run_file_name"])
    charge = validate_identifier(resolved["charge"])
    modified = validate_identifier(resolved["modified_sequence"])
    sequence = validate_identifier(resolved["sequence"])
    label = "NULL::VARCHAR" if labelled else "'LFQ'"
    parts = [
        f"regexp_replace(r.{run}, '\\.(mzML|raw|d|wiff|htrms)$', '', 'i') AS run_file_name",
        "lk.peptidoform AS peptidoform",
        f"CAST(r.{charge} AS SMALLINT) AS charge",
        f"{label} AS label",
        'r."F.FrgIon" AS ion_type',
        f"{_optional_cast('F.Charge', 'INTEGER', report_cols)} AS fragment_charge",
        "NULLIF(CAST(r.\"F.FrgLossType\" AS VARCHAR), 'noloss') AS neutral_loss"
        if "F.FrgLossType" in report_cols
        else "NULL::VARCHAR AS neutral_loss",
        *(
            f"{_measured_sql(column, kind)} AS {name}" if column in report_cols else f"NULL::{kind} AS {name}"
            for column, name, kind in _COLUMNS
        ),
        f"{_optional_cast('F.ExcludedFromQuantification', 'BOOLEAN', report_cols)} AS excluded",
        f"{_scores_sql(report_cols)} AS additional_scores",
    ]
    return (
        f"SELECT {', '.join(parts)} FROM report r JOIN precursor_lookup lk "
        f"ON r.{modified} = lk.modified_sequence AND r.{sequence} = lk.sequence AND r.{charge} = lk.charge "
        'WHERE r."F.FrgIon" IS NOT NULL'
    )
