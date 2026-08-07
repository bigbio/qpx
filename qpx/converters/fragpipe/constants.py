"""FragPipe converter constants — PTM parsing."""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Optional

from qpx.converters.ptm import build_proforma, from_proforma, mass_to_unimod

# ---------------------------------------------------------------------------
# Decoy detection
# ---------------------------------------------------------------------------
# FragPipe/Philosopher tag decoy accessions with one of these prefixes. Keep the
# set identical across the psm/feature/pg adapters so a decoy is labelled the
# same everywhere. Detection must run on the RAW accession, before any
# ``sp|ACC|NAME`` parsing strips the prefix.
DECOY_PREFIXES = ("rev_", "REV_", "DECOY_", "decoy_")


def is_decoy_accession(protein_field: str) -> bool:
    """True if any accession in a (possibly comma-joined) protein field is a decoy.

    Operates on the raw string so ``rev_sp|P123|NAME`` is still recognised as a
    decoy (``parse_uniprot_id`` would otherwise drop the ``rev_`` prefix).
    """
    if not protein_field:
        return False
    return any(part.strip().startswith(DECOY_PREFIXES) for part in str(protein_field).split(","))


# ---------------------------------------------------------------------------
# PTM parsing: FragPipe "Assigned Modifications" -> ProForma
# ---------------------------------------------------------------------------

_FP_MOD_TOKEN_RE = re.compile(r"(?:(\d+)([A-Z])|(N-term))\(([^)]+)\)")


@lru_cache(maxsize=8192)
def to_proforma(assigned_mods: str, sequence: str) -> str:
    """Build ProForma from a FragPipe sequence and Assigned Modifications string.

    Args:
        assigned_mods: Assigned Modifications string, e.g. ``5M(15.9949)`` or
            ``N-term(42.0106), 5M(15.9949)``.
        sequence: Stripped peptide sequence.

    Returns:
        ProForma string, e.g. ``PEPTM[UNIMOD:35]IDEK``.
    """
    if not sequence:
        return ""
    if not assigned_mods or not assigned_mods.strip():
        return sequence

    pos_mods: dict[int, str] = {}
    nterm_tag: Optional[str] = None

    for token in assigned_mods.split(","):
        token = token.strip()
        if not token:
            continue

        m = _FP_MOD_TOKEN_RE.match(token)
        if not m:
            continue

        pos_str, _aa, nterm, mass_str = m.groups()
        try:
            mass = float(mass_str)
        except ValueError:
            continue

        unimod_id = mass_to_unimod(mass)
        if unimod_id is not None:
            tag = f"UNIMOD:{unimod_id}"
        else:
            tag = f"+{mass}" if mass >= 0 else str(mass)

        if nterm:
            nterm_tag = tag
        elif pos_str is not None:
            pos_mods[int(pos_str)] = tag

    mods: list[tuple[int, str]] = []
    if nterm_tag:
        mods.append((0, nterm_tag))
    for pos, tag in pos_mods.items():
        mods.append((pos, tag))

    return build_proforma(sequence, mods)


def to_modifications(assigned_mods: str, sequence: str) -> tuple[str, list[dict] | None]:
    """Parse modifications from FragPipe Assigned Modifications format.

    Args:
        assigned_mods: Assigned Modifications string, e.g. ``5M(15.9949)``.
        sequence: Stripped peptide sequence.

    Returns:
        Tuple of (peptidoform, modifications) where modifications is a list
        of modification dicts per QPX schema, or ``None`` if unmodified.
    """
    proforma = to_proforma(assigned_mods, sequence)
    return from_proforma(proforma, sequence, meta=None)
