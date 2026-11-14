"""MaxQuant SILAC multiplets: channel labels, ratios and the partners identified by MS/MS."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import duckdb
import pandas as pd
from sdrf_pipelines.converters.channel_map import CHANNEL_LABELS

from qpx.converters.maxquant.constants import MQ_LABEL_LETTERS
from qpx.converters.utils import safe_float
from qpx.core.sql import sql_build

# MaxQuant's whole-multiplet ratio columns: "Ratio H/L".
_RATIO_COLUMN = re.compile(r"Ratio ([LMH])/([LMH])")

# SILAC labels by plex position, lightest first.
_SILAC_ORDER = list(CHANNEL_LABELS["silac3plex"].values())


def ratio_letters(columns) -> set[str]:
    """Multiplet letters of the ratio columns (``Ratio H/L``) MaxQuant writes for SILAC; empty otherwise."""
    return {letter for column in columns if (match := _RATIO_COLUMN.fullmatch(column)) for letter in match.groups()}


def multiplet_labels(letters, run_labels: dict[str, set[str]], experiment_type: str) -> dict[str, str]:
    """MaxQuant multiplet letter (``L``, ``M``, ``H``) -> canonical SILAC label, lightest first.

    The letters take the SDRF's SILAC labels in plex order, or ``SILAC
    light``/``medium``/``heavy`` without an SDRF. Empty for output without
    multiplets; an SDRF that disagrees with the output raises.
    """
    letters = sorted(set(letters), key=MQ_LABEL_LETTERS.index)
    if not letters:
        if experiment_type == "SILAC":
            raise ValueError("The SDRF declares SILAC labels, but the MaxQuant output has no SILAC multiplets")
        return {}
    if not run_labels:
        return dict(zip(letters, ["SILAC light", "SILAC heavy"] if len(letters) == 2 else _SILAC_ORDER))
    declared = set().union(*run_labels.values())
    labels = [label for label in _SILAC_ORDER if label in declared]
    if experiment_type != "SILAC" or len(labels) != len(letters):
        raise ValueError(
            f"The MaxQuant output has SILAC multiplets {'/'.join(letters)}, "
            f"but the SDRF declares {experiment_type} labels {sorted(declared)}"
        )
    return dict(zip(letters, labels))


def ratio_entries(row: dict, multiplet: dict[str, str], letter: str, experiment: str = "") -> list[dict]:
    """``ratio_<num>_<den>`` and ``ratio_<num>_<den>_normalized`` of *letter* over each lighter label.

    MaxQuant divides the heavier label by the lighter (``Ratio H/L``, ``Ratio H/M``).
    *experiment* names the proteinGroups.txt columns of one experiment (``Ratio H/L Mix1``).
    """
    entries = []
    letters = list(multiplet)
    for denominator in letters[: letters.index(letter)]:
        for column_suffix, name_suffix in (("", ""), (" normalized", "_normalized")):
            column = " ".join(filter(None, (f"Ratio {letter}/{denominator}{column_suffix}", experiment)))
            value = safe_float(row.get(column))
            if value is not None:
                name = f"ratio_{letter.lower()}_{denominator.lower()}{name_suffix}"
                entries.append({"intensity_name": name, "intensity_value": value})
    return entries


def multiplet_intensities(row: dict, multiplet: dict[str, str], run_labels: Optional[set[str]]) -> tuple:
    """One evidence intensity per SILAC label the SDRF assigns in the run, and each label's ratios as numerator.

    MaxQuant leaves the label intensities empty when the partners cannot be
    told apart (no labelled residue); the labels then have null intensities.
    """
    intensities: list[dict] = []
    additional_intensities: list[dict] = []
    for letter, label in multiplet.items():
        if run_labels is not None and label not in run_labels:
            continue
        value = safe_float(row.get(f"Intensity {letter}"))
        intensities.append({"label": label, "intensity": value if value is not None and value > 0 else None})
        entries = ratio_entries(row, multiplet, letter)
        if entries:
            additional_intensities.append({"label": label, "intensities": entries})
    return intensities, additional_intensities or None


def identified_states(conn: duckdb.DuckDBPyConnection, columns: set[str], evidence_path: str) -> dict[int, list[int]]:
    """evidence ``id`` -> the SILAC labeling states identified by MS/MS, from the loaded ``evidence`` table.

    evidence.txt leaves ``Labeling State`` empty on a multiplet with several
    sequenced partners (``MULTI-MSMS``), whose states are in the msms.txt
    MaxQuant writes beside it.
    """
    if "id" not in columns:
        return {}
    sources = ['SELECT TRY_CAST("id" AS BIGINT) AS id, TRY_CAST("Labeling State" AS INTEGER) AS state FROM evidence']
    params = []
    msms_path = Path(evidence_path).with_name("msms.txt")
    if msms_path.exists():
        sources.append(
            'SELECT TRY_CAST("Evidence ID" AS BIGINT), TRY_CAST("Labeling state" AS INTEGER)'
            " FROM read_csv_auto($1, delim='\t', header=true, auto_detect=true, null_padding=true)"
        )
        params.append(str(msms_path))
    rows = conn.execute(
        sql_build(
            "SELECT id, LIST(DISTINCT state) FROM ($sources) WHERE state >= 0 GROUP BY id",
            sources=" UNION ALL ".join(sources),
        ),
        params,
    ).fetchall()
    return dict(rows)


def identified_labels(evidence_ids, state, multiplet: dict[str, str], states_by_evidence: dict) -> Optional[list[str]]:
    """Labels of the SILAC partners identified by MS/MS; ``None`` when no labeling state is known.

    *evidence_ids* are the evidence rows behind a feature, looked up in
    *states_by_evidence*; without it, the row's own ``Labeling State`` *state*.
    """
    if states_by_evidence:
        states = {
            found
            for evidence_id in evidence_ids
            if pd.notna(evidence_id)
            for found in states_by_evidence.get(int(evidence_id), ())
        }
    else:
        states = {int(state)} if pd.notna(state) and state >= 0 else set()
    letters = list(multiplet)
    return [multiplet[letters[found]] for found in sorted(states)] or None
