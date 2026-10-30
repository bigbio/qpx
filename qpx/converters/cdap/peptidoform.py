"""CDAP peptidoform parsing -- delta-mass annotations to ProForma.

CDAP encodes modifications as signed delta masses inline with the peptide
sequence::

    +229.163PEPTM+15.995IDEK+229.163
    ↓
    [UNIMOD:737]-PEPTM[UNIMOD:35]IDEK[UNIMOD:737]

The N-terminal label tag (e.g. TMT/iTRAQ on the alpha amine) appears before
the first residue.  Subsequent ``±mass`` tokens annotate the residue that
immediately precedes them.  Tokens are mapped to UNIMOD accessions via the
mass tables in :mod:`qpx.converters.cdap.constants`; unknown masses are kept
as ProForma mass shifts (``[+123.4560]``) so that downstream tooling can still
report them.
"""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Iterable

from qpx.converters.cdap.constants import ALL_MOD_MASSES
from qpx.converters.ptm import build_proforma, format_mass_shift

# Regex capturing a single signed delta-mass token (e.g. ``+229.163``,
# ``-17.027``).  Anchored to start-of-string when scanning the prefix; the
# token must contain a decimal point so we don't grab residue numbering.
_MASS_TOKEN_RE = re.compile(r"^([+-]\d+\.\d+)")


def _resolve_mass_to_unimod(mass_token: str) -> str:
    """Map a CDAP delta-mass token (e.g. ``+229.163``) to a UNIMOD tag or a mass shift.

    The :data:`ALL_MOD_MASSES` table is the authoritative source.  Anything not
    found there is kept as a signed mass shift so the converter never silently
    drops a modification.
    """
    accession = ALL_MOD_MASSES.get(mass_token)
    if accession is not None:
        return f"UNIMOD:{accession}"
    return format_mass_shift(float(mass_token))


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
    mods: list[tuple[int, str]] = []
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
            tag = _resolve_mass_to_unimod(mass_token)
            position = last_residue_pos  # 0 means N-term in build_proforma()
            mods.append((position, tag))
            i += len(mass_token)
        elif ch.isalpha():
            sequence_chars.append(ch.upper())
            last_residue_pos += 1
            i += 1
        else:
            # Drop unexpected characters (whitespace, punctuation) silently --
            # the survey saw none of these but stay defensive.
            i += 1

    return build_proforma("".join(sequence_chars), mods)


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
