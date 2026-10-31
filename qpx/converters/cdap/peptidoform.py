"""CDAP peptidoform parsing -- delta-mass annotations to ProForma.

CDAP encodes modifications as signed delta masses inline with the peptide
sequence::

    +229.163PEPTM+15.995IDEK+229.163
    ↓
    [UNIMOD:737]-PEPTM[UNIMOD:35]IDEK[UNIMOD:737]

The N-terminal label tag (e.g. TMT/iTRAQ on the alpha amine) appears before
the first residue.  Subsequent ``±mass`` tokens annotate the residue that
immediately precedes them.  Tokens are mapped to Unimod accessions at their
residue or N-terminus with :func:`qpx.converters.ptm.mass_to_unimod`; masses
without a Unimod entry there are passed through unchanged as ``CHEMMOD``
strings so that downstream tooling can still report them.
"""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Iterable

from qpx.converters.ptm import build_proforma, mass_to_unimod, position_site

# Regex capturing a single signed delta-mass token (e.g. ``+229.163``,
# ``-17.027``).  Anchored to start-of-string when scanning the prefix; the
# token must contain a decimal point so we don't grab residue numbering.
_MASS_TOKEN_RE = re.compile(r"^([+-]\d+\.\d+)")


def _resolve_mass_to_unimod(mass_token: str, sequence: str, position: int) -> str:
    """Map a CDAP delta-mass token (e.g. ``+229.163``) at *position* to a UNIMOD or CHEMMOD tag.

    Position 0 is the N-terminus.  A mass without a Unimod entry at that residue
    or terminus is preserved as a ``CHEMMOD:±mass`` string so the converter never
    silently drops a modification.
    """
    accession = mass_to_unimod(float(mass_token), *position_site(sequence, position))
    if accession is not None:
        return f"UNIMOD:{accession}"
    return f"CHEMMOD:{mass_token}"


@lru_cache(maxsize=65536)
def cdap_to_proforma(peptide_seq: str) -> str:
    """Convert a CDAP ``PeptideSequence`` value to ProForma.

    Examples
    --------
    ``+229.163DTHEDHD``                  -> ``[UNIMOD:737]-DTHEDHD``
    ``+229.163HNM+15.995VNGR``           -> ``[UNIMOD:737]-HNM[UNIMOD:35]VNGR``
    ``+229.163K+229.163PEPTIDE``         -> ``[UNIMOD:737]-K[UNIMOD:737]PEPTIDE``
    ``DTHEDHD``                          -> ``DTHEDHD`` (no modifications)
    """
    if not isinstance(peptide_seq, str) or not peptide_seq:
        return ""

    sequence_chars: list[str] = []
    mods: list[tuple[int, str]] = []  # (position, mass token)
    i = 0
    n = len(peptide_seq)
    last_residue_pos = 0  # 1-indexed position of the most recent amino acid

    while i < n:
        ch = peptide_seq[i]
        if ch in "+-":
            match = _MASS_TOKEN_RE.match(peptide_seq[i:])
            if not match:
                # Stray sign without a numeric payload -- preserve and advance
                # by one character so we never get stuck in an infinite loop.
                i += 1
                continue
            mass_token = match.group(1)
            # 0 means N-term in build_proforma(); resolved once the residues are known.
            mods.append((last_residue_pos, mass_token))
            i += len(mass_token)
        elif ch.isalpha():
            sequence_chars.append(ch.upper())
            last_residue_pos += 1
            i += 1
        else:
            # Drop unexpected characters (whitespace, punctuation) silently --
            # the survey saw none of these but stay defensive.
            i += 1

    sequence = "".join(sequence_chars)
    return build_proforma(sequence, [(pos, _resolve_mass_to_unimod(token, sequence, pos)) for pos, token in mods])


@lru_cache(maxsize=65536)
def strip_label_tags(peptide_seq: str) -> str:
    """Return the bare peptide sequence with all delta-mass tokens removed.

    Every ``±mass`` token (label tags and PTMs alike) is dropped, leaving a
    residue-only string for the QPX ``sequence`` field.
    """
    if not isinstance(peptide_seq, str) or not peptide_seq:
        return ""

    out: list[str] = []
    i = 0
    n = len(peptide_seq)
    while i < n:
        ch = peptide_seq[i]
        if ch in "+-":
            match = _MASS_TOKEN_RE.match(peptide_seq[i:])
            if match:
                i += len(match.group(1))
                continue
            i += 1
        elif ch.isalpha():
            out.append(ch.upper())
            i += 1
        else:
            i += 1
    return "".join(out)


def iter_modification_masses(peptide_seq: str) -> Iterable[str]:
    """Yield every delta-mass token in *peptide_seq* in left-to-right order.

    Useful for label auto-detection -- the caller can inspect the masses to
    decide whether to attempt a TMT/iTRAQ identification when channel column
    names are absent or ambiguous.
    """
    if not isinstance(peptide_seq, str) or not peptide_seq:
        return
    i = 0
    n = len(peptide_seq)
    while i < n:
        ch = peptide_seq[i]
        if ch in "+-":
            match = _MASS_TOKEN_RE.match(peptide_seq[i:])
            if match:
                yield match.group(1)
                i += len(match.group(1))
                continue
        i += 1


__all__ = [
    "cdap_to_proforma",
    "iter_modification_masses",
    "strip_label_tags",
]
