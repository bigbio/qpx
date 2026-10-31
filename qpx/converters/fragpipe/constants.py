"""FragPipe converter constants — PTM parsing."""

from __future__ import annotations

import re
from functools import lru_cache

from qpx.converters.ptm import (
    PROTON_MASS,
    build_proforma,
    compute_precursor_mz,
    format_mass_shift,
    from_proforma,
    mass_tag,
    parse_unimod_delta,
    position_site,
)
from qpx.converters.utils import parse_uniprot_id

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


def protein_accessions(*fields: object) -> list[str]:
    """Distinct accessions of comma-separated FragPipe protein fields, in order.

    ``sp|P12345|NAME`` gives ``P12345``, the form of FragPipe's ``Protein ID``; a
    decoy keeps its prefix (``rev_P12345``) so it never merges with its target.
    Missing fields are skipped.
    """
    accessions: dict[str, None] = {}
    for field in fields:
        if not isinstance(field, str):
            continue
        for part in field.split(","):
            part = part.strip()
            if part:
                prefix = next((p for p in DECOY_PREFIXES if part.startswith(p)), "")
                accessions[prefix + parse_uniprot_id(part[len(prefix) :])[0]] = None
    return list(accessions)


# ---------------------------------------------------------------------------
# PTM parsing: FragPipe "Assigned Modifications" -> ProForma
# ---------------------------------------------------------------------------

_FP_MOD_TOKEN_RE = re.compile(r"(?:(\d+)([A-Z])|([NC]-term))\(([^)]+)\)")


def _token_modification(token: str, sequence: str) -> tuple[int, str] | None:
    """Position and ProForma tag of one Assigned Modifications token, or ``None`` if unparseable.

    The mass is resolved to a Unimod accession at the token's residue or terminus
    (see :func:`qpx.converters.ptm.mass_to_unimod`); a mass without a Unimod
    entry there is kept as the signed mass.
    """
    m = _FP_MOD_TOKEN_RE.match(token.strip())
    if not m:
        return None
    pos_str, _aa, terminus, mass_str = m.groups()
    try:
        mass = float(mass_str)
    except ValueError:
        return None
    if terminus == "N-term":
        position = 0
    elif terminus == "C-term":
        position = len(sequence) + 1
    else:
        position = int(pos_str)
    return position, mass_tag(mass, *position_site(sequence, position))


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
    mods = [mod for token in assigned_mods.split(",") if (mod := _token_modification(token, sequence)) is not None]
    return build_proforma(sequence, mods)


# A precursor up to this many 13C peaks off the monoisotopic mass, within the
# window (Da or ppm of the peptide mass, the wider), is an isotope error of the
# precursor selection rather than a modification. MSFragger's isotope_error
# goes from -1 to 3 and its default precursor tolerance is 20 ppm.
_C13_SPACING = 1.0033548
_ISOTOPE_ERRORS = range(-1, 4)
_ISOTOPE_ERROR_DA = 0.1
_ISOTOPE_ERROR_PPM = 20.0


def unassigned_mass_shift(
    assigned_mods: str, sequence: str, calculated_mass: float | None, delta_mass: float | None
) -> float | None:
    """The measured mass that Assigned Modifications leave unexplained, as in an open or mass-offset search.

    The measured mass is ``Calculated Peptide Mass + Delta Mass``, which holds
    whether PTM-Shepherd moved a glycan out of ``Delta Mass`` or only added it to
    Assigned Modifications. ``None`` when there is no shift or it is an isotope error.
    """
    if calculated_mass is None or delta_mass is None:
        return None
    assigned_mz = compute_precursor_mz(to_proforma(assigned_mods, sequence), 1)
    if assigned_mz is None:
        return None
    assigned_mass = assigned_mz - PROTON_MASS
    shift = calculated_mass + delta_mass - assigned_mass
    window = max(_ISOTOPE_ERROR_DA, assigned_mass * _ISOTOPE_ERROR_PPM * 1e-6)
    isotopes = round(shift / _C13_SPACING)
    if isotopes in _ISOTOPE_ERRORS and abs(shift - isotopes * _C13_SPACING) <= window:
        return None
    return shift


def place_mass_shift(assigned_mods: str, sequence: str, shift: float, localization: object) -> tuple[str, list[dict] | None]:
    """Peptidoform with an unassigned mass shift, and its candidate sites when there are several.

    ``MSFragger Localization`` marks the best sites of the shift in lower case.
    One site takes the shift; otherwise its position is unknown (``[+129.0804]?PEPTIDE``),
    and several sites are listed as ``modification_candidates``.
    """
    sites = []
    if isinstance(localization, str) and localization.upper() == sequence:
        sites = [position for position, residue in enumerate(localization, start=1) if residue.islower()]
    if len(sites) == 1:
        site = f"{sites[0]}{sequence[sites[0] - 1]}({shift:.4f})"
        return to_proforma(f"{assigned_mods}, {site}" if assigned_mods else site, sequence), None
    tag = format_mass_shift(shift)
    candidates = [
        {
            "position": position,
            "amino_acid": sequence[position - 1],
            "name": tag,
            "accession": None,
            "mass_delta": shift,
            "probability": None,
            "assigned": False,
            "alternative_group": 1,
            "scores": None,
        }
        for position in sites
    ]
    return f"[{tag}]?{to_proforma(assigned_mods, sequence)}", candidates or None


# Total Glycan Composition of PTM-Shepherd: HexNAc(2)Hex(9) % 1864.6341
_GLYCAN_COMPOSITION_RE = re.compile(r"((?:[A-Za-z][\w-]*\(\d+\))+) % ([\d.]+)")
_GLYCAN_RESIDUE_RE = re.compile(r"([A-Za-z][\w-]*)\((\d+)\)")
# The assigned glycan's mass and its composition's come from one number, rounded.
_GLYCAN_MASS_TOL = 0.001


def name_glycan(modifications: list[dict] | None, composition: object) -> list[dict] | None:
    """Name the modification that holds a glycan's mass after its ProForma composition, ``Glycan:HexNAc2Hex9``.

    The modification is found by mass, whether the peptidoform holds it as a mass or
    a Unimod accession. The input list comes from the ``from_proforma`` cache and is
    not changed.
    """
    match = _GLYCAN_COMPOSITION_RE.fullmatch(composition) if isinstance(composition, str) else None
    if not match or not modifications:
        return modifications
    name = "Glycan:" + "".join(f"{residue}{count}" for residue, count in _GLYCAN_RESIDUE_RE.findall(match.group(1)))
    mass = float(match.group(2))
    return [
        {**mod, "name": name} if abs(parse_unimod_delta(mod["name"]) - mass) <= _GLYCAN_MASS_TOL else mod for mod in modifications
    ]


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
