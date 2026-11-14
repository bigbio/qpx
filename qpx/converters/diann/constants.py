"""DIA-NN converter constants — PTM parsing."""

from __future__ import annotations

import re
from functools import lru_cache

from qpx.converters.ptm import (
    build_proforma,
    format_mass_shift,
    from_proforma,
    modification_tag,
    normalize_peptidoform,
    position_site,
)

# DIA-NN's plexDIA SILAC channels tag every labelled residue with a zero-mass
# (SILAC) modification; the feature peptidoform is the light partner's and has
# no SILAC label.
_LABEL_TAGS = frozenset({"silac"})

# ---------------------------------------------------------------------------
# PTM parsing: DIA-NN Modified.Sequence -> ProForma
# ---------------------------------------------------------------------------


@lru_cache(maxsize=8192)
def to_proforma(modified_sequence: str, mod_masses: tuple[tuple[str, float], ...] = ()) -> str:
    """Convert a DIA-NN Modified.Sequence to canonical ProForma notation.

    DIA-NN encodes modifications with parentheses::

        PEPTM(UniMod:35)IDE       -> PEPTM[UNIMOD:35]IDE
        (UniMod:1)PEPTIDE         -> [UNIMOD:1]-PEPTIDE
        (mTRAQ)PEPTIDEK(mTRAQ)    -> [UNIMOD:888]-PEPTIDEK[UNIMOD:888]

    A modification named in *mod_masses* (the ``--fixed-mod`` / ``--var-mod``
    declarations of the DIA-NN log, which ``--original-mods`` reports keep by
    name) is resolved from its declared mass; other names are resolved through
    Unimod. The (SILAC) channel tag is dropped. DIA-NN 2.x prefixes a decoy
    searched from the spectral library (a first-pass decoy kept in an empirical
    library) with ``^``, which is not part of the sequence and is dropped::

        ^(UniMod:1)MEPSK          -> [UNIMOD:1]-MEPSK

    Args:
        modified_sequence: The ``Modified.Sequence`` value from DIA-NN output.
        mod_masses: Declared ``(name, mass delta)`` pairs.

    Returns:
        ProForma string, e.g. ``PEPTM[UNIMOD:35]IDEK``.
    """
    if not modified_sequence:
        return ""

    modified_sequence = modified_sequence.removeprefix("^")
    masses = dict(mod_masses)
    mods: list[tuple[int, str]] = []
    plain_chars: list[str] = []
    i = 0
    n = len(modified_sequence)
    seq_pos = 0

    while i < n:
        if modified_sequence[i] == "(":
            end = modified_sequence.index(")", i)
            inner = modified_sequence[i + 1 : end]
            i = end + 1
            mass = masses.get(inner)
            if inner.casefold() in _LABEL_TAGS or mass == 0:
                continue
            tag = format_mass_shift(mass) if mass is not None else re.sub(r"(?i)unimod:", "UNIMOD:", inner)
            mods.append((seq_pos, tag))
        else:
            plain_chars.append(modified_sequence[i])
            seq_pos += 1
            i += 1

    plain_seq = "".join(plain_chars)
    return normalize_peptidoform(build_proforma(plain_seq, mods))


# Site.Occupancy.Probabilities: S(UniMod:21){0.92}T(UniMod:21){0.08}K
_OCCUPANCY_TOKEN_RE = re.compile(r"([A-Z])|\(([^()]*)\)(?:\{([\d.eE+-]+)\})?")


def occupancy_sites(occupancy: str, sequence: str, mod_masses: tuple[tuple[str, float], ...] = ()) -> dict[str, dict[int, float]]:
    """Per-modification site occupancy probabilities from ``Site.Occupancy.Probabilities``.

    Keys are the canonical tags of :func:`to_proforma`. Only variable modifications
    have a probability in curly brackets.
    """
    masses = dict(mod_masses)
    sites: dict[str, dict[int, float]] = {}
    position = 0
    for residue, name, probability in _OCCUPANCY_TOKEN_RE.findall(occupancy):
        if residue:
            position += 1
        elif probability:
            mass = masses.get(name)
            tag = format_mass_shift(mass) if mass is not None else modification_tag(name, *position_site(sequence, position))
            sites.setdefault(tag, {})[position] = float(probability)
    return sites


def to_modifications(
    modified_sequence: str, sequence: str, mod_masses: tuple[tuple[str, float], ...] = ()
) -> tuple[str, list[dict] | None]:
    """Parse modifications from a DIA-NN Modified.Sequence string.

    Converts to ProForma first, then delegates to the shared ``from_proforma``
    parser to produce the QPX modification list.

    Args:
        modified_sequence: The ``Modified.Sequence`` value from DIA-NN output.
        sequence: Stripped peptide sequence (no modification annotations).
        mod_masses: Declared ``(name, mass delta)`` pairs, as for :func:`to_proforma`.

    Returns:
        Tuple of (peptidoform, modifications) where modifications is a list
        of modification dicts per QPX schema, or ``None`` if unmodified.
    """
    proforma = to_proforma(modified_sequence, mod_masses)
    return from_proforma(proforma, sequence, meta=None)
