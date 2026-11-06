"""MaxQuant ``<modification>Sites.txt``: one site.parquet row per site, experiment and label."""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

from qpx.converters.maxquant.silac import multiplet_labels, ratio_letters
from qpx.converters.ptm import name_to_unimod, unimod_mass, unimod_title
from qpx.converters.utils import canonical_accession, safe_float
from qpx.writers import SiteWriter

logger = logging.getLogger(__name__)

_SITES_FILE = re.compile(r"^(.+)Sites\.txt(?:\.gz)?$")
_REPORTER = re.compile(r"Reporter intensity corrected (\d+)(?: (.+))?")
_SCORES = (
    ("PEP", "posterior_error_probability", False),
    ("Score", "andromeda_score", True),
    ("Score diff", "maxquant_score_diff", True),
    ("Delta score", "maxquant_delta_score", True),
    ("Score for localization", "maxquant_score_for_localization", True),
)
_TRACE_COLUMNS = {
    "Proteins": "proteins",
    "Positions within proteins": "positions_within_proteins",
    "Best localization raw file": "best_localization_raw_file",
    "Best localization scan number": "best_localization_scan_number",
}


def _modification(path: Path) -> dict:
    """The modification a file holds, from its name: ``Phospho (STY)Sites.txt``."""
    match = _SITES_FILE.match(path.name)
    uid = name_to_unimod(match.group(1)) if match else None
    if uid is None:
        raise ValueError(f"Cannot tell the modification of {path.name}: expected <Unimod name>Sites.txt")
    return {"modification_name": unimod_title(uid), "modification_accession": f"UNIMOD:{uid}", "mass_delta": unimod_mass(uid)}


def _cell(row: dict, *parts: str):
    return row.get(" ".join(filter(None, parts)))


def _experiment_fields(row: dict, experiment: str) -> dict:
    """Scores of the site in one experiment, else over all of them."""

    def value(column: str):
        specific = safe_float(_cell(row, column, experiment))
        return specific if specific is not None else safe_float(row.get(column))

    scores = [
        {"score_name": name, "score_value": score, "higher_better": higher}
        for column, name, higher in _SCORES
        if (score := value(column)) is not None
    ]
    return {"modified_probability": value("Localization prob"), "additional_scores": scores or None}


def _quantities(row: dict, experiment: str, labels: list[tuple[str, str, str]]) -> list[tuple]:
    """``(label, intensity, stoichiometry)`` of each measured label of an experiment."""
    quantities = []
    for label, intensity_column, occupancy_column in labels:
        intensity = safe_float(_cell(row, intensity_column, experiment))
        if intensity:
            quantities.append((label, intensity, safe_float(_cell(row, occupancy_column, experiment))))
    return quantities


def _site_rows(row: dict, modification: dict, experiment_runs: dict[str, list[str]], labels: list) -> list[dict]:
    positions = [cv for column, name in _TRACE_COLUMNS.items() if (cv := _trace(row, column, name))]
    if row.get("Potential contaminant") == "+":
        positions.append({"cv_name": "potential_contaminant", "cv_value": "true"})
    base = {
        **modification,
        "protein_accession": canonical_accession(row["Protein"]),
        "position": int(float(row["Position"])),
        "amino_acid": row.get("Amino acid"),
        "psm_count": len(str(row["MS/MS IDs"]).split(";")) if pd.notna(row.get("MS/MS IDs")) else None,
        "is_decoy": row.get("Reverse") == "+",
        "cv_params": positions or None,
    }
    rows = []
    for experiment, runs in experiment_runs.items():
        scores = _experiment_fields(row, experiment)
        rows.extend(
            {**base, **scores, "grouped_runs": runs, "label": label, "intensity": intensity, "stoichiometry": stoichiometry}
            for label, intensity, stoichiometry in _quantities(row, experiment, labels)
        )
    if rows:
        return rows
    all_runs = sorted({run for runs in experiment_runs.values() for run in runs})
    return [{**base, **_experiment_fields(row, ""), "grouped_runs": all_runs}]


def _trace(row: dict, column: str, name: str) -> dict | None:
    value = row.get(column)
    return {"cv_name": name, "cv_value": str(value)} if pd.notna(value) else None


def _labels(columns: list[str], channels: list[tuple[str, int]], run_labels, experiment_type: str) -> list[tuple[str, str, str]]:
    """``(label, intensity column, occupancy column)``, the columns without their experiment."""
    multiplet = multiplet_labels(ratio_letters(columns), run_labels, experiment_type)
    if multiplet:
        return [(label, f"Intensity {letter}", f"Occupancy {letter}") for letter, label in multiplet.items()]
    numbers = sorted({int(match.group(1)) for column in columns if (match := _REPORTER.fullmatch(column))})
    if channels and numbers:
        return [(label, f"Reporter intensity corrected {numbers[index]}", "") for label, index in channels]
    return [("LFQ", "Intensity", "Occupancy")]


def convert_sites_file(path: Path, experiment_runs: dict[str, list[str]], load_sdrf, writer: SiteWriter) -> int:
    """Write one ``Sites.txt``; *load_sdrf* gives the run labels, experiment type and channels from its columns."""
    table = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    columns = list(table.columns)
    run_labels, experiment_type, channels = load_sdrf(columns)
    labels = _labels(columns, channels, run_labels, experiment_type)
    modification = _modification(path)
    # One experiment has no suffix on its columns; "Intensity <experiment>" then does not exist.
    if not any(_cell({column: 1 for column in columns}, labels[0][1], experiment) for experiment in experiment_runs):
        experiment_runs = {"": sorted({run for runs in experiment_runs.values() for run in runs})}
    rows = [
        site
        for row in table.to_dict("records")
        if pd.notna(row.get("Protein")) and pd.notna(row.get("Position"))
        for site in _site_rows(row, modification, experiment_runs, labels)
    ]
    if rows:
        writer.write_batch(rows)
    logger.info("%s: %d site rows", path.name, len(rows))
    return len(rows)
