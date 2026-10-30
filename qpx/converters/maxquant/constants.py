"""MaxQuant converter constants — PTM parsing and fixed modifications."""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Optional

from qpx.converters.ptm import build_proforma

PROTON_MASS: float = 1.007_276_466_77

# ---------------------------------------------------------------------------
# MaxQuant modification name -> UNIMOD accession mapping
# ---------------------------------------------------------------------------

MQ_MOD_TO_UNIMOD: dict[str, int] = {
    "oxidation (m)": 35,
    "acetyl (protein n-term)": 1,
    "carbamidomethyl (c)": 4,
    "phospho (sty)": 21,
    "phospho (s)": 21,
    "phospho (t)": 21,
    "phospho (y)": 21,
    "deamidation (nq)": 7,
    "deamidation (n)": 7,
    "deamidation (q)": 7,
    "gln->pyro-glu": 28,
    "methyl (kr)": 34,
    "methyl (k)": 34,
    "methyl (r)": 34,
    "dimethyl (kr)": 36,
    "dimethyl (k)": 36,
    "dimethyl (r)": 36,
    "trimethyl (k)": 37,
    "gg (k)": 121,
    "formyl (protein n-term)": 122,
    "carbamyl (protein n-term)": 5,
    "tmt6plex (k)": 737,
    "tmt6plex (n-term)": 737,
    "tmt10plex (k)": 737,
    "tmt10plex (n-term)": 737,
    "tmt11plex (k)": 737,
    "tmt11plex (n-term)": 737,
    "tmt16plex (k)": 730,
    "tmt16plex (n-term)": 730,
    "tmtpro (k)": 730,
    "tmtpro (n-term)": 730,
    "itraq4plex (k)": 214,
    "itraq4plex (n-term)": 214,
    "itraq8plex (k)": 385,
    "itraq8plex (n-term)": 385,
    "pyro-glu from e": 214,
}

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

    uid = MQ_MOD_TO_UNIMOD.get(tag.lower())
    if uid is None:
        uid = MQ_SHORT_TO_UNIMOD.get(tag.lower())
    if uid is not None:
        return f"UNIMOD:{uid}"

    return tag


def _closing_paren_end(seq: str, start: int) -> int:
    """Index just past the parenthesis closing the one at ``start``, tracking nesting."""
    depth = 1
    j = start + 1
    while j < len(seq) and depth > 0:
        if seq[j] == "(":
            depth += 1
        elif seq[j] == ")":
            depth -= 1
        j += 1
    return j


def _split_modified_sequence(seq: str) -> tuple[str, list[tuple[int, str]]]:
    """Split an underscore-free Modified sequence into residues and ``(position, tag)`` pairs.

    A modification takes the position of the residue before it, 0 before the
    first residue (N-term).
    """
    mods: list[tuple[int, str]] = []
    plain_chars: list[str] = []
    i = 0
    while i < len(seq):
        if seq[i] == "(":
            j = _closing_paren_end(seq, i)
            mods.append((len(plain_chars), _resolve_mq_mod(seq[i + 1 : j - 1])))
            i = j
        else:
            plain_chars.append(seq[i])
            i += 1
    return "".join(plain_chars), mods


FixedModification = tuple[str, str, str]
"""``(accession, target residues, terminus)``; terminus is ``""``, ``"N"`` or ``"C"``."""


def _fixed_positions(sequence: str, residues: str, terminus: str) -> list[int]:
    """Where a fixed modification applies on ``sequence`` (0 = N-term, ``len + 1`` = C-term)."""
    if not terminus:
        return [position for position, residue in enumerate(sequence, start=1) if residue in residues]
    position, residue = (0, sequence[:1]) if terminus == "N" else (len(sequence) + 1, sequence[-1:])
    return [position] if not residues or residue in residues else []


def fixed_modification_sites(
    sequence: str, mods: list[tuple[int, str]], fixed_mods: tuple[FixedModification, ...]
) -> list[tuple[int, str]]:
    """The ``(position, accession)`` pairs the fixed modifications add to a peptide.

    MaxQuant leaves fixed modifications out of ``Modified sequence``. Each one
    is placed on every matching site that carries no reported modification,
    so a site never gets a second modification (an N-terminal acetyl keeps a
    fixed N-term label off, as in the search).
    """
    occupied = {position for position, _ in mods}
    sites: list[tuple[int, str]] = []
    for accession, residues, terminus in fixed_mods:
        for position in _fixed_positions(sequence, residues, terminus):
            if position not in occupied:
                occupied.add(position)
                sites.append((position, accession))
    return sites


@lru_cache(maxsize=8192)
def to_proforma(modified_sequence: Optional[str], fixed_mods: tuple[FixedModification, ...] = ()) -> str:
    """Convert a MaxQuant Modified sequence to ProForma notation.

    MaxQuant wraps sequences in underscores and uses parenthetical
    modification names (possibly nested)::

        _PEPTM(Oxidation (M))IDEK_         -> PEPTM[UNIMOD:35]IDEK
        _(Acetyl (Protein N-term))PEPTIDEK_ -> [UNIMOD:1]-PEPTIDEK
        _(ac)PEPTIDEK_                      -> [UNIMOD:1]-PEPTIDEK

    ``fixed_mods`` adds the fixed modifications MaxQuant leaves out (see
    :func:`fixed_modification_sites`)::

        _PEPCK_ with (("UNIMOD:4", "C", ""),) -> PEPC[UNIMOD:4]K
    """
    if not isinstance(modified_sequence, str) or not modified_sequence:
        return ""

    seq = modified_sequence.strip("_")
    if not seq:
        return ""

    plain_seq, mods = _split_modified_sequence(seq)
    return build_proforma(plain_seq, mods + fixed_modification_sites(plain_seq, mods, fixed_mods))


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
