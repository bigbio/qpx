"""Row repairs of ``qpxc migrate`` that apply to the files of one producer."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from qpx.converters.ptm import _add_nterm_dash
from qpx.version import parse_spec_version

# The cv_params entries in which DIA-NN and Spectronaut files kept the indexed retention time before its column.
_IRT_PARAMS = ("irt", "empirical_irt")
# The MSNet psm columns that scripts/convert_msnet.py renames, in files it did not convert.
_RAW_MSNET_COLUMNS = {
    "precursor_charge": "charge",
    "exp_mass_to_charge": "observed_mz",
    "cal_mass_to_charge": "calculated_mz",
    "retention_time": "rt",
    "reference_file_name": "run_file_name",
}
# DIA-NN feature intensity names before specification 1.3.
_DIANN_INTENSITY_NAMES = {"maxlfq": "pg_maxlfq", "precursor_normalised": "normalize_intensity"}
# DIA-NN reports 0 for these without ion mobility data; a conversion writes null, or leaves the cv entry out.
_DIANN_ION_MOBILITY_COLUMNS = ("ion_mobility", "ion_mobility_start", "ion_mobility_stop")
_DIANN_ION_MOBILITY_PARAMS = ("iim", "predicted_im", "predicted_iim")
# The first specification whose migration rescales the coverage of the quantms converter (removed in 1.1).
_PERCENT_QUANTMS_COVERAGE = (1, 3)


def _msnet_row(row: dict) -> dict:
    """MSNet's consensus_support, is_decoy and cv_params as scripts/convert_msnet.py writes them.

    consensus_support becomes a score of additional_scores. MSNet holds targets only; writer 1.2.1.dev5
    left the is_decoy of some files null and kept cv_params without a value, which convert_msnet leaves out.
    """
    value = row.pop("consensus_support")
    if value is not None:
        score = {"score_name": "consensus_support", "score_value": float(value), "higher_better": None}
        row["additional_scores"] = [*(row.get("additional_scores") or []), score]
    if row.get("is_decoy") is None:
        row["is_decoy"] = False
    if row.get("cv_params"):
        row["cv_params"] = [entry for entry in row["cv_params"] if entry["cv_value"] is not None]
    return row


def _raw_msnet_row(record: dict) -> dict:
    """An MSNet psm row that scripts/convert_msnet.py did not convert, in the columns it writes.

    mass_error_ppm is computed from the float32 m/z values the row keeps. The global q-value
    becomes the first score; _msnet_row adds consensus_support after it.
    """
    row = {_RAW_MSNET_COLUMNS.get(key, key): value for key, value in record.items()}
    observed, calculated = float(np.float32(row["observed_mz"])), float(np.float32(row["calculated_mz"]))
    row["mass_error_ppm"] = 1e6 * (observed - calculated) / calculated
    qvalue = row.pop("global_qvalue")
    row["additional_scores"] = [{"score_name": "global_qvalue", "score_value": float(qvalue), "higher_better": None}]
    row["cv_params"] = [{"cv_name": key, "cv_value": str(value)} for key, value in row["cv_params"].items() if value is not None]
    row["scan"] = [row["scan"]]
    row["modifications"] = row["modifications"] and [
        {
            "name": entry["name"],
            "accession": None,
            "positions": [{"position": int(site["position"]), "amino_acid": None, "scores": None} for site in entry["positions"]],
        }
        for entry in row["modifications"]
    ]
    return row


def _diann_intensity_names(entries: list[dict]) -> list[dict]:
    return [
        {
            **entry,
            "intensities": [
                {**item, "intensity_name": _DIANN_INTENSITY_NAMES.get(item["intensity_name"], item["intensity_name"])}
                for item in entry["intensities"]
            ],
        }
        for entry in entries
    ]


def _diann_feature_row(record: dict, measured_mz: bool) -> dict:
    """Recover the library-decoy prefix separately from the report's is_decoy flag, and the current intensity names.

    Without *measured_mz*, observed_mz held the theoretical Precursor.Mz and mass_error_ppm a float32 rounding
    step; both become null, as a conversion writes them where DIA-NN measured no MS1 apex.
    """
    peptidoform = record["peptidoform"]
    library_prefix = peptidoform.startswith("^")
    row = {**record, "is_library_decoy": record.get("is_library_decoy") or library_prefix}
    if library_prefix:
        row["peptidoform"] = _add_nterm_dash(peptidoform.removeprefix("^"))
    if not measured_mz:
        row["observed_mz"] = row["mass_error_ppm"] = None
    if row.get("additional_intensities"):
        row["additional_intensities"] = _diann_intensity_names(row["additional_intensities"])
    for column in _DIANN_ION_MOBILITY_COLUMNS:
        if row.get(column) == 0:
            row[column] = None
    if row.get("cv_params"):
        row["cv_params"] = _diann_cv_params(row["cv_params"], row.get("additional_scores"))
    return row


def _diann_cv_params(params: list[dict], scores: list[dict] | None) -> list[dict] | None:
    """The cv entries a conversion writes: none for a zero ion mobility, and Quantity.Quality only as its score.

    Writer 1.1 also kept Quantity.Quality as precursor_quantification_score, a float32 string of the score
    diann_quantity_quality; an entry with another value is kept.
    """
    quality = next((score["score_value"] for score in scores or [] if score["score_name"] == "diann_quantity_quality"), None)

    def repeated(entry: dict) -> bool:
        name, value = entry["cv_name"], entry["cv_value"]
        if value is None:
            return False
        if name in _DIANN_ION_MOBILITY_PARAMS:
            return float(value) == 0
        return name == "precursor_quantification_score" and quality is not None and np.float32(value) == np.float32(quality)

    return [entry for entry in params if not repeated(entry)] or None


def _diann_pg_row(row: dict) -> dict:
    """A DIA-NN pg row without the score pg_qvalue that writer 1.1 repeated from the pg_qvalue column."""
    scores = row.get("additional_scores")
    if not scores:
        return row
    kept = [
        score for score in scores if not (score["score_name"] == "pg_qvalue" and score["score_value"] == row.get("pg_qvalue"))
    ]
    return {**row, "additional_scores": kept or None}


def _quantms_pg_row(row: dict) -> dict:
    """A quantms pg row with sequence_coverage in percent: QPX 1.0 wrote the mzTab fraction, and 0 for none."""
    coverage = row.get("sequence_coverage")
    return {**row, "sequence_coverage": coverage * 100 if coverage else None}


def _pg_row_fix(footer: dict[bytes, bytes]) -> Callable[[dict], dict] | None:
    """The producer fix of the pg rows of a file, or None."""
    creator = footer.get(b"creator")
    if creator == b"diann":
        return _diann_pg_row
    version = parse_spec_version(footer.get(b"qpx_version", b"").decode()) or (0, 0)
    if creator == b"quantms" and version < _PERCENT_QUANTMS_COVERAGE:
        return _quantms_pg_row
    return None


def _irt_row(row: dict) -> dict:
    """The irt of a feature row in its column, as a conversion writes it, and no longer among its cv_params."""
    params = row.get("cv_params") or []
    values = [entry["cv_value"] for entry in params if entry["cv_name"] in _IRT_PARAMS]
    if not values:
        return row
    irt = row.get("irt")
    if irt is None and values[0] is not None:
        irt = float(values[0])
    return {**row, "irt": irt, "cv_params": [entry for entry in params if entry["cv_name"] not in _IRT_PARAMS] or None}
