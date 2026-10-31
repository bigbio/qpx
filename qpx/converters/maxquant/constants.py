"""MaxQuant converter constants — PTM parsing and fixed modifications."""

from __future__ import annotations

from functools import lru_cache
from typing import Optional

from qpx.converters.ptm import build_proforma, modification_tag, position_site, split_modified_sequence, vendor_position

PROTON_MASS: float = 1.007_276_466_77

# ---------------------------------------------------------------------------
# MaxQuant modification names
# ---------------------------------------------------------------------------

# Short forms of older MaxQuant versions. Full names (``Oxidation (M)``) are
# resolved through Unimod by :func:`qpx.converters.ptm.modification_tag`.
MQ_SHORT_TO_UNIMOD: dict[str, int] = {
    "ox": 35,
    "ac": 1,
    "ph": 21,
    "de": 7,
    "gl": 28,
    "me": 34,
}


def _resolve_mq_mod(tag: str, residue: Optional[str] = None, terminus: Optional[str] = None) -> str:
    """Resolve a MaxQuant modification name at *residue* / *terminus* to a ProForma tag.

    Handles short forms (``ac``), full names (``Oxidation (M)``, ``GlyGly (K)``),
    ``UniMod:N`` accessions and signed masses; an unknown name is kept as a name.
    """
    uid = MQ_SHORT_TO_UNIMOD.get(tag.lower())
    if uid is not None:
        return f"UNIMOD:{uid}"
    return modification_tag(tag, residue, terminus)


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
        _PEPK(GlyGly (K))R_                 -> PEPK[UNIMOD:121]R
        _PEPK(Amidated (Protein C-term))_   -> PEPK-[UNIMOD:2]

    ``fixed_mods`` adds the fixed modifications MaxQuant leaves out (see
    :func:`fixed_modification_sites`)::

        _PEPCK_ with (("UNIMOD:4", "C", ""),) -> PEPC[UNIMOD:4]K
    """
    if not isinstance(modified_sequence, str) or not modified_sequence.strip("_"):
        return ""
    sequence, mods = split_modified_sequence(modified_sequence.strip("_"))
    tags: list[tuple[int, str]] = []
    for position, name in mods:
        position = vendor_position(name, position, sequence)
        tags.append((position, _resolve_mq_mod(name, *position_site(sequence, position))))
    return build_proforma(sequence, tags + fixed_modification_sites(sequence, tags, fixed_mods))


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
