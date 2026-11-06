"""Per-site localization probabilities: from a producer's values to scores and candidates."""

from __future__ import annotations

import re
from typing import Optional, Tuple

from qpx.converters.ptm import parse_unimod_delta, unimod_mass, unimod_title

LOCALIZATION_SCORE = "localization_probability"

_SITE_PROBABILITY_RE = re.compile(r"([A-Za-z])(?:\(([\d.eE+-]+)\))?")
# Candidates of one alternative group sum to at most 1, as the validator allows for rounding.
_ALTERNATIVE_GROUP_SUM = 1.01


def parse_site_probabilities(text: object, sequence: str) -> dict[int, float]:
    """Position -> probability from a sequence with probabilities in parentheses, ``S(0.06)RT(0.94)PPRR``.

    MaxQuant's ``<mod> Probabilities`` and PTMProphet's ``STY:79.9663`` columns have
    this form; PTMProphet writes the termini as ``n`` and ``c``. Empty for a missing
    value and for one that does not spell *sequence*.
    """
    if not isinstance(text, str):
        return {}
    residues = []
    probabilities = {}
    for residue, probability in _SITE_PROBABILITY_RE.findall(text.strip("_")):
        if residue == "n":
            position = 0
        elif residue == "c":
            position = len(sequence) + 1
        else:
            residues.append(residue)
            position = len(residues)
        if probability:
            probabilities[position] = float(probability)
    return probabilities if "".join(residues) == sequence else {}


def site_candidate(
    tag: str,
    position: int,
    sequence: str,
    probability: Optional[float],
    assigned: bool,
    alternative_group: Optional[int] = None,
) -> dict:
    """The ``modification_candidates`` entry of a canonical tag at one position."""
    name, accession, mass = tag, None, parse_unimod_delta(tag)
    if tag.startswith("UNIMOD:"):
        uid = int(tag.removeprefix("UNIMOD:"))
        name, accession, mass = unimod_title(uid) or tag, tag, unimod_mass(uid)
    return {
        "position": position,
        "amino_acid": sequence[position - 1] if 0 < position <= len(sequence) else None,
        "name": name,
        "accession": accession,
        "mass_delta": mass,
        "probability": probability,
        "assigned": assigned,
        "alternative_group": alternative_group,
        "scores": None,
    }


def _score_placed_sites(mod: dict, probabilities: dict[int, float]) -> dict:
    """A copy of *mod* with a ``localization_probability`` score on each site that has a probability."""
    positions = []
    for pos in mod["positions"]:
        if pos["position"] in probabilities:
            score = {"score_name": LOCALIZATION_SCORE, "score_value": probabilities[pos["position"]], "higher_better": True}
            pos = {**pos, "scores": [*(pos.get("scores") or []), score]}
        positions.append(pos)
    return {**mod, "positions": positions}


def _mod_key(mod: dict) -> str:
    return mod["accession"] or mod["name"]


def _modification_candidates(
    tag: str, probabilities: dict[int, float], occupied: set[int], sequence: str, group: Optional[int]
) -> list[dict]:
    """The candidates of one modification; a site with no chance and no placement is left out."""
    return [
        site_candidate(tag, position, sequence, probability, position in occupied, group)
        for position, probability in sorted(probabilities.items())
        if probability > 0 or position in occupied
    ]


def localize_modifications(
    modifications: Optional[list[dict]],
    sequence: str,
    sites: dict[str, dict[int, float]],
    first_group: int = 1,
) -> Tuple[Optional[list[dict]], Optional[list[dict]]]:
    """Attach per-site probabilities to ``modifications`` and list every scored site as a candidate.

    ``sites`` maps a canonical tag (``UNIMOD:21``, ``+79.9663``) to position -> probability.
    The sites a modification occupies get a ``localization_probability`` score; a
    candidate is ``assigned`` there. A modification placed once whose probabilities
    sum to about 1 is one alternative group, numbered from *first_group*. Several
    placements or a larger sum (k of n sites, or independent occupancies) stay ungrouped.
    The input list may come from the ``from_proforma`` cache and is not changed.
    """
    if not sites:
        return modifications, None
    mods = modifications or []
    placed: dict[str, set[int]] = {}
    for mod in mods:
        placed.setdefault(_mod_key(mod), set()).update(pos["position"] for pos in mod["positions"])
    candidates = []
    group = first_group
    for tag, probabilities in sites.items():
        occupied = placed.get(tag, set())
        grouped = len(occupied) == 1 and sum(probabilities.values()) <= _ALTERNATIVE_GROUP_SUM
        candidates += _modification_candidates(tag, probabilities, occupied, sequence, group if grouped else None)
        group += grouped
    localized = [_score_placed_sites(mod, sites.get(_mod_key(mod), {})) for mod in mods]
    return localized or modifications, candidates or None
