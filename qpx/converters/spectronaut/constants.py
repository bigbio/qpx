"""Spectronaut converter constants — PTM parsing."""

from __future__ import annotations

import re
from functools import lru_cache

from qpx.converters.ptm import build_proforma, from_proforma, unimod_from_name

_NTERM_RE = re.compile(r"N-term|Protein N-term")

# ---------------------------------------------------------------------------
# PTM parsing: Spectronaut Modified.Sequence -> ProForma
# ---------------------------------------------------------------------------

_MOD_TOKEN_RE = re.compile(r"\[([^\]]+)\]")


@lru_cache(maxsize=8192)
def to_proforma(modified_sequence: str) -> str:
    """Convert Spectronaut ``EG.ModifiedSequence`` to ProForma notation.

    Spectronaut encodes modifications in brackets, wrapped in underscores,
    e.g. ``_PIGLC[Carbamidomethyl (C)]IAPVLAAK_``.

    Args:
        modified_sequence: The ``EG.ModifiedSequence`` value.

    Returns:
        ProForma string.
    """
    if not modified_sequence:
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
        if seq[i] == "[":
            end = seq.index("]", i)
            mod_name = seq[i + 1 : end]
            unimod_id = unimod_from_name(mod_name)
            if unimod_id is not None:
                tag = f"UNIMOD:{unimod_id}"
            else:
                tag = mod_name

            position = 0 if seq_pos == 0 else seq_pos
            if _NTERM_RE.search(mod_name) and seq_pos > 0:
                position = 0
            mods.append((position, tag))
            i = end + 1
        else:
            plain_chars.append(seq[i])
            seq_pos += 1
            i += 1

    plain_seq = "".join(plain_chars)
    return build_proforma(plain_seq, mods)


def to_modifications(modified_sequence: str, sequence: str) -> tuple[str, list[dict] | None]:
    """Parse modifications from a Spectronaut EG.ModifiedSequence string.

    Converts to ProForma first, then delegates to the shared ``from_proforma``
    parser to produce the QPX modification list.

    Args:
        modified_sequence: The ``EG.ModifiedSequence`` value from Spectronaut.
        sequence: Stripped peptide sequence (no modification annotations).

    Returns:
        Tuple of (peptidoform, modifications) where modifications is a list
        of modification dicts per QPX schema, or ``None`` if unmodified.
    """
    proforma = to_proforma(modified_sequence)
    return from_proforma(proforma, sequence, meta=None)
