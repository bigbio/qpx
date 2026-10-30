"""MaxQuant converter constants — PTM parsing and label maps."""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Optional

from qpx.converters.ptm import NTERM, build_proforma, compute_precursor_mz, unimod_from_name

PROTON_MASS: float = 1.007_276_466_77

# ---------------------------------------------------------------------------
# MaxQuant modification name -> UNIMOD accession mapping
# ---------------------------------------------------------------------------

# Older MaxQuant versions abbreviate modifications in "Modified sequence"; full
# names such as "Oxidation (M)" resolve through :func:`unimod_from_name`.
MQ_SHORT_TO_UNIMOD: dict[str, int] = {
    "ox": 35,
    "ac": 1,
    "ph": 21,
    "de": 7,
    "gl": 28,
    "me": 34,
}


def _resolve_mq_mod(tag: str) -> str:
    """Resolve a MaxQuant modification name to ``UNIMOD:N`` or return as-is.

    Handles full names (``Oxidation (M)``), short forms (``ac``),
    and already-resolved ``UNIMOD:`` prefixes.
    """
    m = re.match(r"(?i)unimod:(\d+)", tag)
    if m:
        return f"UNIMOD:{m.group(1)}"

    uid = MQ_SHORT_TO_UNIMOD.get(tag.lower()) or unimod_from_name(tag)
    if uid is not None:
        return f"UNIMOD:{uid}"

    return tag


FixedMods = tuple[tuple[str, str], ...]


def _site_list(sites: str) -> list[str]:
    """Residues and ``N-term`` named by a MaxQuant or SDRF site field; protein N-term sites are skipped."""
    text = sites.strip()
    if text.lower() in ("n-term", "any n-term"):
        return [NTERM]
    if "term" in text.lower():
        return []
    return [aa for aa in text.replace(",", "") if aa.isalpha() and aa.isupper()]


def sdrf_fixed_mods(entries: Optional[list[dict]]) -> FixedMods:
    """Fixed modifications from the SDRF ``comment[modification parameters]`` entries.

    *entries* is :meth:`SdrfConverter.parse_modification_params` output.
    """
    fixed: list[tuple[str, str]] = []
    for entry in entries or []:
        if not entry["fixed"]:
            continue
        accession = re.match(r"(?i)unimod:(\d+)$", entry["accession"] or "")
        uid = int(accession.group(1)) if accession else unimod_from_name(entry["name"] or "")
        if uid is None:
            continue
        position = (entry["position"] or "").strip().lower()
        if position in ("any n-term", "n-term"):
            sites = [NTERM]
        elif "term" in position:
            sites = []
        else:
            sites = _site_list(entry["target_amino_acid"] or "")
        fixed.extend((site, f"UNIMOD:{uid}") for site in sites)
    return tuple(fixed)


def parameters_fixed_mods(parameters_path) -> FixedMods:
    """Fixed modifications from the ``Fixed modifications`` row of MaxQuant's ``parameters.txt``."""
    try:
        with open(parameters_path, encoding="utf-8") as handle:
            row = next((line for line in handle if line.startswith("Fixed modifications\t")), None)
    except FileNotFoundError:
        return ()
    if row is None:
        return ()
    fixed: list[tuple[str, str]] = []
    for name in row.split("\t", 1)[1].strip().split(";"):
        match = re.match(r"(.+?)\s+\(([^()]*)\)$", name.strip())
        uid = unimod_from_name(name) if match else None
        if uid is not None:
            fixed.extend((site, f"UNIMOD:{uid}") for site in _site_list(match.group(2)))
    return tuple(fixed)


def _with_fixed_mods(sequence: str, mods: list[tuple[int, str]], fixed_mods: FixedMods) -> list[tuple[int, str]]:
    """Add the fixed modifications MaxQuant leaves out of ``Modified sequence``.

    A fixed modification goes on every matching site that carries no other modification.
    """
    taken = {position for position, _ in mods}
    extra = [(0, tag) for site, tag in fixed_mods if site == NTERM and 0 not in taken]
    extra += [
        (position, tag)
        for position, aa in enumerate(sequence, start=1)
        if position not in taken
        for site, tag in fixed_mods
        if site == aa
    ]
    return mods + extra


@lru_cache(maxsize=8192)
def to_proforma(modified_sequence: Optional[str], fixed_mods: FixedMods = ()) -> str:
    """Convert a MaxQuant Modified sequence to ProForma notation.

    MaxQuant wraps sequences in underscores and uses parenthetical
    modification names (possibly nested)::

        _PEPTM(Oxidation (M))IDEK_         -> PEPTM[UNIMOD:35]IDEK
        _(Acetyl (Protein N-term))PEPTIDEK_ -> [UNIMOD:1]-PEPTIDEK
        _(ac)PEPTIDEK_                      -> [UNIMOD:1]-PEPTIDEK

    *fixed_mods* (``(site, tag)`` pairs) are written inline, since MaxQuant
    omits fixed modifications and isobaric labels from the sequence.
    """
    if not isinstance(modified_sequence, str) or not modified_sequence:
        return ""

    seq = modified_sequence.strip("_")
    if not seq:
        return ""

    mods: list[tuple[int, str]] = []
    plain_chars: list[str] = []
    i = 0
    n = len(seq)
    seq_pos = 0

    while i < n:
        if seq[i] == "(":
            # Walk forward tracking parenthesis depth to handle nested parens
            depth = 1
            j = i + 1
            while j < n and depth > 0:
                if seq[j] == "(":
                    depth += 1
                elif seq[j] == ")":
                    depth -= 1
                j += 1
            mod_name = seq[i + 1 : j - 1]
            tag = _resolve_mq_mod(mod_name)
            position = 0 if seq_pos == 0 else seq_pos
            mods.append((position, tag))
            i = j
        else:
            plain_chars.append(seq[i])
            seq_pos += 1
            i += 1

    plain_seq = "".join(plain_chars)
    return build_proforma(plain_seq, _with_fixed_mods(plain_seq, mods, fixed_mods))


def mq_calculated_mz(peptidoform: Optional[str], charge: int, neutral_mass: Optional[float]) -> Optional[float]:
    """Theoretical m/z of *peptidoform*; MaxQuant's ``Mass`` when a modification has no known mass."""
    calculated = compute_precursor_mz(peptidoform, charge) if peptidoform else None
    if calculated is None and neutral_mass and charge:
        calculated = (neutral_mass + charge * PROTON_MASS) / charge
    return calculated


def parse_phospho_probabilities(
    prob_string: str,
) -> dict[int, list[dict]] | None:
    """Parse MaxQuant ``Phospho (STY) Probabilities`` into per-position scores.

    MaxQuant encodes per-position phosphorylation probabilities using the same
    underscore-wrapped format as Modified sequence, but with floating-point
    probabilities in parentheses instead of modification names::

        _PEPTIDES(0.95)T(0.03)Y(0.02)K_

    Only positions with non-zero probability are included in the output.

    Args:
        prob_string: The ``Phospho (STY) Probabilities`` column value.

    Returns:
        Dict mapping 1-indexed position -> list of score dicts, or ``None``
        if no probabilities found.
    """
    if not isinstance(prob_string, str) or not prob_string:
        return None

    seq = prob_string.strip("_")
    if not seq:
        return None

    site_scores: dict[int, list[dict]] = {}
    seq_pos = 0
    i = 0
    n = len(seq)

    while i < n:
        if seq[i] == "(":
            end = seq.index(")", i)
            prob_str = seq[i + 1 : end]
            try:
                prob = float(prob_str)
            except ValueError:
                i = end + 1
                continue
            if prob > 0:
                site_scores[seq_pos] = [
                    {
                        "score_name": "phospho_sty_probability",
                        "score_value": prob,
                        "higher_better": True,
                    }
                ]
            i = end + 1
        else:
            seq_pos += 1
            i += 1

    return site_scores if site_scores else None


PHOSPHO_MASS = 79.966331


def phospho_candidates(
    site_scores: dict[int, list[dict]] | None, sequence: str, modifications: list[dict] | None
) -> list[dict] | None:
    """Every residue MaxQuant gives a phosphorylation probability, as modification candidates.

    ``site_scores`` is :func:`parse_phospho_probabilities`'s output. A candidate
    is ``assigned`` when the peptidoform places the phosphorylation there; the
    others are the residues it competes with, which ``modifications`` cannot hold.
    With one phosphate, the candidates are one alternative group (their
    probabilities sum to 1); with several, MaxQuant's probabilities are
    reported per site and sum to the number of phosphates, which no single
    group can express, so the candidates are left ungrouped. Both sums hold
    to within the rounding of MaxQuant's three decimals (PXD071171, PXD072739).
    """
    if not site_scores:
        return None
    placed = {
        pos["position"]
        for mod in modifications or []
        if mod.get("accession") == "UNIMOD:21" or mod.get("name") == "Phospho"
        for pos in mod.get("positions") or []
    }
    group = 1 if len(placed) == 1 else None
    candidates = []
    for position in sorted(site_scores):
        probability = site_scores[position][0]["score_value"]
        candidates.append(
            {
                "position": position,
                "amino_acid": sequence[position - 1] if 0 < position <= len(sequence) else None,
                "name": "Phospho",
                "accession": "UNIMOD:21",
                "mass_delta": PHOSPHO_MASS,
                "probability": probability,
                "assigned": position in placed,
                "alternative_group": group,
                "scores": None,
            }
        )
    return candidates


# MaxQuant orders reporter ions strictly by mass (N-isomer before C-isomer), so
# alphabetical sorting of channel labels (C < N) would map each label to the
# wrong column.  The i+1 fallback that was previously used to handle 0- vs
# 1-indexed MaxQuant columns caused adjacent channels to read from the same
# column, producing duplicate intensities (e.g. TMT126 ≡ TMT127C in 96 % of
# features for PXD016999).  Using this explicit map avoids both issues.
TMT_LABEL_TO_MQ_COL: dict[str, int] = {
    "TMT126": 0,
    "TMT127N": 1,
    "TMT127C": 2,
    "TMT128N": 3,
    "TMT128C": 4,
    "TMT129N": 5,
    "TMT129C": 6,
    "TMT130N": 7,
    "TMT130C": 8,
    "TMT131": 9,
    "TMT131N": 9,
    "TMT131C": 10,
    "TMT132N": 11,
    "TMT132C": 12,
    "TMT133N": 13,
    "TMT133C": 14,
    "TMT134N": 15,
    "TMT134C": 16,
    "TMT135N": 17,
    # iTRAQ 4-plex
    "ITRAQ4PLEX-114": 0,
    "ITRAQ4PLEX-115": 1,
    "ITRAQ4PLEX-116": 2,
    "ITRAQ4PLEX-117": 3,
    # iTRAQ 8-plex
    "ITRAQ8PLEX-113": 0,
    "ITRAQ8PLEX-114": 1,
    "ITRAQ8PLEX-115": 2,
    "ITRAQ8PLEX-116": 3,
    "ITRAQ8PLEX-117": 4,
    "ITRAQ8PLEX-118": 5,
    "ITRAQ8PLEX-119": 6,
    "ITRAQ8PLEX-121": 7,
}
