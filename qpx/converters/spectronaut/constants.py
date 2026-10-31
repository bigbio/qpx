"""Spectronaut converter constants — PTM parsing."""

from __future__ import annotations

from functools import lru_cache

from qpx.converters.ptm import (
    build_proforma,
    from_proforma,
    modification_tag,
    position_site,
    split_modified_sequence,
    vendor_position,
)

# ---------------------------------------------------------------------------
# PTM parsing: Spectronaut Modified.Sequence -> ProForma
# ---------------------------------------------------------------------------


@lru_cache(maxsize=8192)
def to_proforma(modified_sequence: str) -> str:
    """Convert Spectronaut ``EG.ModifiedSequence`` to ProForma notation.

    Spectronaut encodes modifications in brackets, wrapped in underscores,
    e.g. ``_PIGLC[Carbamidomethyl (C)]IAPVLAAK_``. Names are resolved through
    Unimod by :func:`qpx.converters.ptm.modification_tag`; a name for the N- or
    C-terminus is placed on that terminus.

    Args:
        modified_sequence: The ``EG.ModifiedSequence`` value.

    Returns:
        ProForma string.
    """
    if not modified_sequence:
        return ""
    sequence, mods = split_modified_sequence(modified_sequence.strip("_"), "[", "]")
    tags: list[tuple[int, str]] = []
    for position, name in mods:
        position = vendor_position(name, position, sequence)
        tags.append((position, modification_tag(name, *position_site(sequence, position))))
    return build_proforma(sequence, tags)


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
