"""mzIdentML numeric and modification-localization scores."""

from __future__ import annotations

import re

from qpx.converters.localization import localize_modifications, site_candidate
from qpx.core.cv_terms import (
    CV_PEP,
    CV_PEP_GLOBAL,
    SITE_LOCALIZATION_ACCESSIONS,
    SITE_PROBABILITY_ACCESSIONS,
    SKIP_SCORE_ACCESSIONS,
)
from qpx.core.scores import is_higher_better, normalize_score_name

_SKIP_SCORE_ACCESSIONS = SKIP_SCORE_ACCESSIONS
_PEP_ACCESSIONS = frozenset({CV_PEP, CV_PEP_GLOBAL})
_INDEXED_SCORE_RE = re.compile(r"(\d+):([\d.]+(?:[Ee][+-]?\d+)?):(\d+(?:\|\d+)*):(?:true|false)")


def _parse_position_scores(value_str: str) -> dict[int, float]:
    """Parse residue-position localization scores from one cvParam value."""
    pos_scores: dict[int, float] = {}
    for token in re.split(r"[,;]\s*", value_str):
        m = re.match(r"([A-Za-z])(\d+)\s*:\s*([\d.eE+-]+)", token.strip())
        if m:
            pos = int(m.group(2))
            try:
                pos_scores[pos] = float(m.group(3))
            except ValueError:
                pass

    return pos_scores


def _parse_indexed_scores(value_str: str) -> list[tuple[int, float, list[int]]]:
    """Modification index, score and positions of each MS:1002505 entry in one cvParam value."""
    return [
        (int(index), float(score), [int(position) for position in positions.split("|")])
        for index, score, positions in _INDEXED_SCORE_RE.findall(value_str)
    ]


def _position_score(score_name: str, value: float) -> dict:
    return {"score_name": score_name, "score_value": value, "higher_better": is_higher_better(score_name)}


def _attach_position_scores(modifications, pos_scores, score_name, tag=None) -> None:
    """Attach localization values to the matching positions, of the modification *tag* when it is known."""
    for mod in modifications or []:
        if tag is not None and (mod["accession"] or mod["name"]) != tag:
            continue
        for pos_info in mod.get("positions", []):
            position = pos_info.get("position", 0)
            if position in pos_scores:
                pos_info["scores"] = [*(pos_info.get("scores") or []), _position_score(score_name, pos_scores[position])]


def _add_candidate_score(candidates: list[dict], modifications, sequence: str, tag: str, position: int, score: dict) -> None:
    """Record a non-probability score of a modification at a position on its candidate."""
    for candidate in candidates:
        if (candidate["accession"] or candidate["name"]) == tag and candidate["position"] == position:
            candidate["scores"] = [*(candidate["scores"] or []), score]
            return
    assigned = any(
        (mod["accession"] or mod["name"]) == tag and any(pos["position"] == position for pos in mod["positions"])
        for mod in modifications or []
    )
    candidates.append({**site_candidate(tag, position, sequence, None, assigned), "scores": [score]})


def _attach_site_localization_scores(
    modifications: list[dict] | None,
    sequence: str,
    cv_params: list[dict],
    modification_tags: dict[int, str],
) -> tuple[list[dict] | None, list[dict] | None]:
    """Modifications with their site scores, and the candidates, from the SII's localization cvParams.

    Scores follow MS:1002505, ``1:0.95:4:true``, where the first number is the
    ``modification index`` of a SearchModification (*modification_tags*); older
    files write ``S3: 95.3, T5: 2.1``. Probabilities and scores of other residues
    become candidates, which needs the index. A scalar value says nothing about a
    site and stays in ``additional_scores``.
    """
    probabilities: dict[str, dict[int, float]] = {}
    other_scores: list[tuple[str, int, dict]] = []
    for cv in (cv for cv in cv_params if cv["accession"] in SITE_LOCALIZATION_ACCESSIONS):
        score_name = normalize_score_name(cv["name"])
        value_str = cv.get("value", "")
        for index, score, positions in _parse_indexed_scores(value_str):
            tag = modification_tags.get(index)
            pos_scores = dict.fromkeys(positions, score)
            if tag is not None and cv["accession"] in SITE_PROBABILITY_ACCESSIONS:
                probabilities.setdefault(tag, {}).update(pos_scores)
                continue
            _attach_position_scores(modifications, pos_scores, score_name, tag)
            if tag is not None:
                other_scores.extend((tag, position, _position_score(score_name, score)) for position in positions)
        _attach_position_scores(modifications, _parse_position_scores(value_str), score_name)

    modifications, candidates = localize_modifications(modifications, sequence, probabilities)
    candidates = candidates or []
    for tag, position, score in other_scores:
        _add_candidate_score(candidates, modifications, sequence, tag, position, score)
    return modifications, candidates or None


def _extract_scores(cv_params: list[dict]) -> list[dict]:
    """Extract numeric scores from SII cvParams."""
    scores = []
    for cv in cv_params:
        acc = cv["accession"]
        if acc in _SKIP_SCORE_ACCESSIONS:
            continue
        value_str = cv.get("value", "")
        if not value_str:
            continue
        # Try to parse as float
        try:
            value = float(value_str)
        except (ValueError, TypeError):
            continue
        normalized = normalize_score_name(cv["name"])
        scores.append(
            {
                "score_name": normalized,
                "score_value": value,
                "higher_better": is_higher_better(normalized),
            }
        )
    return scores


def _extract_pep(cv_params: list[dict]) -> float | None:
    """Extract posterior error probability from SII cvParams.

    Looks for MS:1001493 (PEP) or MS:1002352 (global PEP).
    Returns float value or None.
    """
    for cv in cv_params:
        if cv["accession"] in _PEP_ACCESSIONS:
            try:
                return float(cv.get("value", ""))
            except (ValueError, TypeError):
                pass
    return None
