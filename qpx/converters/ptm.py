"""Shared PTM utilities: UNIMOD tables, ProForma builder/parser, precursor m/z.

This flat module consolidates shared logic previously spread across
``qpx.converters.ptm.unimod``, ``qpx.converters.ptm.proforma``,
``qpx.converters.ptm.modifications``, and ``qpx.converters.ptm.__init__``.

Tool-specific constants (MaxQuant _MQ_MOD_TO_UNIMOD, _MQ_SHORT_TO_UNIMOD,
resolve_to_unimod) are **not** included here -- they belong in the respective
converter's ``constants.py``.
"""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Optional, Tuple

# ---------------------------------------------------------------------------
# UNIMOD mass registry
# ---------------------------------------------------------------------------

NTERM = "N-term"


def _sites(residues: str, *termini: str) -> frozenset[str]:
    return frozenset(residues) | frozenset(termini)


# Unimod accession -> (PSI-MS name, monoisotopic mass delta, sites). Sites are
# residues plus "N-term"; they only separate modifications of equal mass.
# tests/converters/test_ptm.py checks every name and mass against the Unimod
# copy bundled with pyopenms.
UNIMOD_MODS: dict[int, tuple[str, float, frozenset[str]]] = {
    1: ("Acetyl", 42.010565, _sites("KSTYCH", NTERM)),
    2: ("Amidated", -0.984016, _sites("")),
    4: ("Carbamidomethyl", 57.021464, _sites("CKHDE", NTERM)),
    5: ("Carbamyl", 43.005814, _sites("KRCM", NTERM)),
    6: ("Carboxymethyl", 58.005479, _sites("CKW", NTERM)),
    7: ("Deamidated", 0.984016, _sites("NQRF")),
    21: ("Phospho", 79.966331, _sites("STYHDCRK")),
    23: ("Dehydrated", -18.010565, _sites("STYD")),
    24: ("Propionamide", 71.037114, _sites("CK", NTERM)),
    26: ("Pyro-carbamidomethyl", 39.994915, _sites("C")),
    27: ("Glu->pyro-Glu", -18.010565, _sites("E")),
    28: ("Gln->pyro-Glu", -17.026549, _sites("Q")),
    34: ("Methyl", 14.01565, _sites("KRHEDCSTNQ", NTERM)),
    35: ("Oxidation", 15.994915, _sites("MWHCPKYF")),
    36: ("Dimethyl", 28.0313, _sites("KRNP", NTERM)),
    37: ("Trimethyl", 42.04695, _sites("KR")),
    39: ("Methylthio", 45.987721, _sites("C")),
    121: ("GG", 114.042927, _sites("KSTC")),
    122: ("Formyl", 27.994915, _sites("KST", NTERM)),
    188: ("Label:13C(6)", 6.020129, _sites("KRLI")),
    199: ("Dimethyl:2H(4)", 32.056407, _sites("K", NTERM)),
    214: ("iTRAQ4plex", 144.102063, _sites("KY", NTERM)),
    259: ("Label:13C(6)15N(2)", 8.014199, _sites("K")),
    267: ("Label:13C(6)15N(4)", 10.008269, _sites("R")),
    312: ("Cysteinyl", 119.004099, _sites("C")),
    385: ("Ammonia-loss", -17.026549, _sites("CN")),
    481: ("Label:2H(4)", 4.025107, _sites("K")),
    510: ("Dimethyl:2H(4)13C(2)", 34.063117, _sites("K", NTERM)),
    730: ("iTRAQ8plex", 304.20536, _sites("KY", NTERM)),
    737: ("TMT6plex", 229.162932, _sites("KSTH", NTERM)),
    738: ("TMT2plex", 225.155833, _sites("K", NTERM)),
    739: ("TMT", 224.152478, _sites("K", NTERM)),
    888: ("mTRAQ", 140.094963, _sites("KY", NTERM)),
    2016: ("TMTpro", 304.207146, _sites("KSTH", NTERM)),
    2017: ("TMTpro_zero", 295.189592, _sites("K", NTERM)),
}

UNIMOD_MASS: dict[int, float] = {uid: mass for uid, (_, mass, _) in UNIMOD_MODS.items()}

# Tool spellings that are not Unimod PSI-MS names, lower-cased.
_UNIMOD_SYNONYMS: dict[str, int] = {
    "acetylation": 1,
    "carbamidomethylation": 4,
    "deamidation": 7,
    "phosphorylation": 21,
    "pyro-glu from e": 27,
    "pyro-glu from q": 28,
    "methylation": 34,
    "dimethylation": 36,
    "trimethylation": 37,
    "glygly": 121,
    "formylation": 122,
    "tmt10plex": 737,
    "tmt11plex": 737,
    "tmt16plex": 2016,
    "tmt18plex": 2016,
}

_UNIMOD_BY_NAME: dict[str, int] = {name.lower(): uid for uid, (name, _, _) in UNIMOD_MODS.items()} | _UNIMOD_SYNONYMS

# MaxQuant and Spectronaut append the allowed sites: "Oxidation (M)".
_SITE_SUFFIX_RE = re.compile(r"\s+\([^()]*\)$")

_MASS_TOL = 0.01

# Candidates whose distances to the observed mass differ by less than this
# cannot be told apart.
_TIE_TOL = 0.0005


def mass_to_unimod(mass: float, residue: Optional[str] = None, nterm: bool = False) -> Optional[int]:
    """Find the UNIMOD accession whose delta is closest to *mass*.

    A *residue* (one-letter code) or *nterm* restricts the search to
    modifications allowed at that site, when any of them is. Returns ``None``
    when nothing lies within ``_MASS_TOL`` or the two closest candidates are
    equally close (Gln->pyro-Glu and Ammonia-loss without a site, for example),
    so the caller keeps the mass shift instead of guessing.
    """
    hits = sorted((abs(ref - mass), uid, sites) for uid, (_, ref, sites) in UNIMOD_MODS.items() if abs(ref - mass) <= _MASS_TOL)
    if residue or nterm:
        hits = [h for h in hits if residue in h[2] or (nterm and NTERM in h[2])] or hits
    if not hits or (len(hits) > 1 and hits[1][0] - hits[0][0] < _TIE_TOL):
        return None
    return hits[0][1]


def unimod_from_name(name: str) -> Optional[int]:
    """Resolve a modification name such as ``GlyGly (K)`` or ``Phospho`` to its UNIMOD accession."""
    return _UNIMOD_BY_NAME.get(_SITE_SUFFIX_RE.sub("", name.strip()).lower())


def format_mass_shift(mass: float) -> str:
    """Return the canonical ProForma tag body for an unknown mass shift, e.g. ``+79.9663``."""
    return f"{mass:+.4f}"


@lru_cache(maxsize=1)
def _pyopenms_unimod_masses() -> dict[int, float]:
    from pyopenms import ModificationsDB  # pylint: disable=no-name-in-module

    db = ModificationsDB()
    masses: dict[int, float] = {}
    for index in range(db.getNumberOfModifications()):
        mod = db.getModification(index)
        masses.setdefault(mod.getUniModRecordId(), mod.getDiffMonoMass())
    return masses


def unimod_mass(uid: int) -> Optional[float]:
    """Monoisotopic mass delta of a UNIMOD accession, from :data:`UNIMOD_MODS` or pyopenms' Unimod copy."""
    if uid in UNIMOD_MASS:
        return UNIMOD_MASS[uid]
    return _pyopenms_unimod_masses().get(uid)


def canonical_mod_tag(
    accession: Optional[str],
    name: Optional[str],
    mass: Optional[float],
    residue: Optional[str] = None,
    nterm: bool = False,
) -> str:
    """ProForma tag body for a modification reported by accession, name and/or mass.

    A UNIMOD accession is kept when its mass agrees with *mass*, and replaced
    by the mass when it does not. An accession QPX has no mass for keeps the
    mass alongside (``UNIMOD:1842|+54.0106``). Other accessions (PSI-MOD,
    PSI-MS "unknown modification") and names are mapped to UNIMOD by name,
    then by mass; anything still unmapped is written as the signed mass.
    """
    unimod = re.match(r"(?i)unimod:(\d+)$", accession or "")
    if unimod:
        return _checked_unimod_tag(int(unimod.group(1)), mass)
    uid = unimod_from_name(name or "")
    if uid is not None and mass is not None and abs(UNIMOD_MASS[uid] - mass) > _MASS_TOL:
        uid = None
    if uid is None and mass is not None:
        uid = mass_to_unimod(mass, residue, nterm)
    if uid is not None:
        return f"UNIMOD:{uid}"
    return format_mass_shift(mass) if mass is not None else accession or name or ""


def _checked_unimod_tag(uid: int, mass: Optional[float]) -> str:
    """``UNIMOD:uid`` if *mass* agrees with it, the mass if not, both if its mass is unknown."""
    ref = unimod_mass(uid)
    if mass is None or (ref is not None and abs(ref - mass) <= _MASS_TOL):
        return f"UNIMOD:{uid}"
    return format_mass_shift(mass) if ref is not None else f"UNIMOD:{uid}|{format_mass_shift(mass)}"


# Monoisotopic amino acid residue masses (Da)
AA_MASS: dict[str, float] = {
    "A": 71.03711,
    "R": 156.10111,
    "N": 114.04293,
    "D": 115.02694,
    "C": 103.00919,
    "E": 129.04259,
    "Q": 128.05858,
    "G": 57.02146,
    "H": 137.05891,
    "I": 113.08406,
    "L": 113.08406,
    "K": 128.09496,
    "M": 131.04049,
    "F": 147.06841,
    "P": 97.05276,
    "S": 87.03203,
    "T": 101.04768,
    "W": 186.07931,
    "Y": 163.06333,
    "V": 99.06841,
}

WATER_MASS = 18.010565
PROTON_MASS = 1.007276


def parse_unimod_delta(mod_str: str) -> Optional[float]:
    """Extract the mass delta from a modification token.

    Accepts:
        - ``UniMod:35`` / ``UNIMOD:35`` -> lookup in ``UNIMOD_MASS``
        - ``+15.9949`` / ``-17.0265`` -> parse directly
        - ``UNIMOD:1842|+54.0106`` -> the explicit mass of a pipe-separated tag
    """
    parts = mod_str.split("|")
    for part in parts:
        try:
            return float(part)
        except ValueError:
            continue
    m = re.match(r"(?:UniMod|UNIMOD):(\d+)", parts[0], re.IGNORECASE)
    if m:
        uid = int(m.group(1))
        return UNIMOD_MASS.get(uid)
    return None


# ---------------------------------------------------------------------------
# ProForma builder
# ---------------------------------------------------------------------------


def build_proforma(sequence: str, mods: list[tuple[int, str]]) -> str:
    """Build a ProForma string from sequence and modification list.

    Args:
        sequence: Stripped peptide sequence (no modification annotations).
        mods: List of (position, tag) where position is 0 for N-term,
            1..N for residue positions (1-indexed), and ``len(sequence)+1``
            for the C-term (the mzIdentML/FragPipe/CDAP convention). Tag is
            e.g. "UNIMOD:35" or "CHEMMOD:+15.99". Multiple mods may share a
            position; all are emitted (order preserved) rather than collapsed.

    Returns:
        ProForma string, e.g. ``PEPTM[UNIMOD:35]IDEK``, ``[UNIMOD:1]-PEPTIDEK``
        (N-term) or ``PEPTIDEK-[UNIMOD:2]`` (C-term).
    """
    if not sequence:
        return ""

    cterm_pos = len(sequence) + 1

    # Preserve every mod, including two at the same position; do NOT collapse
    # via ``dict(mods)`` (which would silently drop all but the last tag).
    pos_to_tags: dict[int, list[str]] = {}
    for position, tag in mods:
        pos_to_tags.setdefault(position, []).append(tag)

    parts: list[str] = []

    # N-terminus: ``[tag][tag]-`` before the first residue.
    if 0 in pos_to_tags:
        for tag in pos_to_tags[0]:
            parts.append(f"[{tag}]")
        parts.append("-")

    for i, aa in enumerate(sequence, start=1):
        parts.append(aa)
        for tag in pos_to_tags.get(i, []):
            parts.append(f"[{tag}]")

    # C-terminus: ``-[tag][tag]`` after the last residue. Any position beyond
    # the last residue is treated as C-terminal.
    cterm_tags = [tag for pos, tags in pos_to_tags.items() if pos >= cterm_pos for tag in tags]
    if cterm_tags:
        parts.append("-")
        for tag in cterm_tags:
            parts.append(f"[{tag}]")

    return "".join(parts)


# ---------------------------------------------------------------------------
# ProForma parser -> list<modification>
# ---------------------------------------------------------------------------


def _normalize_peptidoform(peptidoform: str) -> str:
    """Normalise mzTab-style ``(ModName)`` to ProForma ``[ModName]``.

    mzTab encodes modifications with parentheses (e.g.
    ``C(Carbamidomethyl)R``), whereas ProForma uses square brackets.
    This function converts parenthetical mod tags that follow an amino-acid
    residue (or appear at position 0 for N-term mods) to bracket notation
    so that :func:`_from_proforma_impl` can parse them. Parentheses inside a
    bracketed tag (``K[GlyGly (K)]``, ``K[Label:13C(6)15N(2)]``) are part of
    the tag and stay as they are.
    """
    if "(" not in peptidoform:
        return peptidoform
    out: list[str] = []
    peptidoform = peptidoform.removeprefix(".")
    n = len(peptidoform)
    i = 0
    while i < n:
        if peptidoform[i] == "[":
            end = peptidoform.find("]", i)
            end = n if end == -1 else end + 1
            out.append(peptidoform[i:end])
            i = end
        elif peptidoform[i] == "(":
            # Find matching closing parenthesis respecting nesting depth
            depth = 1
            j = i + 1
            while j < n and depth > 0:
                if peptidoform[j] == "(":
                    depth += 1
                elif peptidoform[j] == ")":
                    depth -= 1
                j += 1
            if depth != 0:
                out.append(peptidoform[i:])
                break
            # Convert outer delimiters; keep inner content as-is
            out.append("[")
            out.append(peptidoform[i + 1 : j - 1])
            out.append("]")
            i = j
        else:
            out.append(peptidoform[i])
            i += 1

    result = "".join(out)

    # For N-term
    if result.startswith("["):
        idx = result.find("]")
        if idx != -1 and idx + 1 < len(result) and result[idx + 1] not in "-?":
            result = result[: idx + 1] + "-" + result[idx + 1 :]

    return result


def _canonical_tag(tag: str, residue: Optional[str], nterm: bool) -> str:
    """Canonical form of one bracketed tag of :func:`normalize_peptidoform`."""
    head, _, rest = tag.partition("|")
    unimod = re.match(r"(?i)unimod:(\d+)$", head)
    if unimod:
        return f"UNIMOD:{unimod.group(1)}" + (f"|{rest}" if rest else "")
    mass = parse_unimod_delta(head.removeprefix("CHEMMOD:"))
    if mass is not None:
        return canonical_mod_tag(None, None, mass, residue, nterm)
    return canonical_mod_tag(None, head, None)


def normalize_peptidoform(peptidoform: str) -> str:
    """Rewrite a ProForma or mzTab-style peptidoform into the canonical form.

    Names and ``UniMod:`` accessions become ``UNIMOD:n``; ``CHEMMOD:`` and
    bare mass shifts become ``UNIMOD:n`` when a Unimod entry of that mass is
    allowed at the site, else ``[+x.xxxx]``. Fixed modifications a producer
    left out cannot be recovered from the string alone.
    """
    text = _normalize_peptidoform(peptidoform)
    out: list[str] = []
    residue: Optional[str] = None
    i = 0
    while i < len(text):
        if text[i] == "[":
            end = text.find("]", i)
            if end == -1:
                out.append(text[i:])
                break
            nterm = residue is None and text[end + 1 : end + 2] == "-"
            site = next((ch for ch in text[end:] if ch.isalpha()), None) if nterm else residue
            out.append(f"[{_canonical_tag(text[i + 1 : end], site, nterm)}]")
            i = end + 1
        else:
            if text[i].isalpha():
                residue = text[i]
            out.append(text[i])
            i += 1
    return "".join(out)


_UNLOCALIZED_TAGS = re.compile(r"(?:\[[^\[\]]+\])+\?")


def _modification_site(seq_pos: int, last_aa: str | None, cterm: bool, sequence_length: int) -> tuple[int, str | None]:
    """Return the QPX position and residue for a ProForma modification."""
    if seq_pos == 0:
        return 0, None
    if cterm:
        return sequence_length + 1, None
    return seq_pos, last_aa


def _modification_identity(mod_str: str, aa: str | None, meta: Optional[dict]) -> tuple[str, str | None]:
    """Resolve a modification tag to its name and accession."""
    mod_str = mod_str.split("|", 1)[0]
    name = mod_str
    accession = None
    if mod_str.startswith("UNIMOD:"):
        accession = mod_str
        if meta:
            entry = meta.get(accession)
            if entry:
                name = entry[0]
        else:
            name = accession
    elif mod_str.startswith("+") or mod_str.startswith("-"):
        if meta:
            for acc, (mod_name, sites, _) in meta.items():
                if aa in sites or "X" in sites:
                    name = mod_name
                    accession = acc
                    break
    return name, accession


def _from_proforma_impl(
    peptidoform: str,
    sequence: str,
    meta: Optional[dict] = None,
    site_scores: Optional[dict[int, list[dict]]] = None,
) -> Tuple[str, Optional[list[dict]]]:
    """Core implementation of ProForma modification parsing.

    See :func:`from_proforma` for full documentation.
    """
    # Normalise mzTab parenthetical notation to ProForma brackets
    peptidoform = _normalize_peptidoform(peptidoform)
    if not peptidoform or peptidoform == sequence:
        return peptidoform, None

    mods: dict[str, dict] = {}
    seq_pos = 0
    n = len(peptidoform)
    last_aa = None
    cterm = False

    # Tags before "?" have an unknown position, so they list no positions.
    i = 0
    unlocalized = _UNLOCALIZED_TAGS.match(peptidoform)
    if unlocalized:
        for mod_str in re.findall(r"\[([^\[\]]+)\]", unlocalized.group()):
            name, accession = _modification_identity(mod_str, None, meta)
            mods.setdefault(accession or name, {"name": name, "accession": accession, "positions": []})
        i = unlocalized.end()
    while i < n:
        if peptidoform[i] == "[":
            try:
                end = peptidoform.index("]", i)
            except ValueError:
                return peptidoform, None  # Malformed ProForma
            mod_str = peptidoform[i + 1 : end]

            position, aa = _modification_site(seq_pos, last_aa, cterm, len(sequence))

            name, accession = _modification_identity(mod_str, aa, meta)

            key = accession or name
            if key not in mods:
                mods[key] = {"name": name, "accession": accession, "positions": []}
            pos_scores = site_scores.get(position) if site_scores else None
            mods[key]["positions"].append({"position": position, "amino_acid": aa, "scores": pos_scores})

            i = end + 1
        elif peptidoform[i] == "-":
            cterm = seq_pos > 0
            i += 1
        else:
            cterm = False
            last_aa = peptidoform[i]
            seq_pos += 1
            i += 1

    mods = list(mods.values()) if mods else None

    return peptidoform, mods


@lru_cache(maxsize=8192)
def _from_proforma_cached(peptidoform: str, sequence: str) -> Tuple[str, Optional[list[dict]]]:
    """Cached fast path for from_proforma when no meta or site_scores."""
    return _from_proforma_impl(peptidoform, sequence)


def from_proforma(
    peptidoform: str,
    sequence: str,
    meta: Optional[dict] = None,
    site_scores: Optional[dict[int, list[dict]]] = None,
) -> Tuple[str, Optional[list[dict]]]:
    """Parse modifications from a ProForma-style peptidoform string.

    Handles: ``M[UNIMOD:35]PEPTIDEK``, ``M[+15.9949]PEPTIDEK``,
    ``[UNIMOD:1]-PEPTIDEK``, ``PEPTIDEK-[UNIMOD:2]``

    Uses an LRU cache for the common case (no *meta*, no *site_scores*).
    Callers must **not** mutate the returned list or its nested dicts, as
    the cached path returns the same object for repeated inputs.

    Args:
        peptidoform: ProForma string with [tag] annotations.
        sequence: Stripped (unmodified) sequence.
        meta: Optional mzTab modifications metadata:
            ``{accession: (name, sites, positions)}``.
            When provided, used to resolve UNIMOD accessions to human-readable
            names and to match mass shifts by site. When None, uses the tag as
            name.
        site_scores: Optional per-position site localization scores.
            Dict mapping position (1-indexed, 0 for N-term) to a list of
            score dicts ``[{score_name, score_value, higher_better}]``.
            Used for phospho site localization probabilities.

    Returns:
        Tuple of (peptidoform, modifications) where peptidoform is the
        normalised ProForma string and modifications is a list of dicts
        (``{name, accession, positions}``) per QPX schema, or ``None``
        if no modifications.
    """
    if meta is None and site_scores is None:
        return _from_proforma_cached(peptidoform, sequence)
    return _from_proforma_impl(peptidoform, sequence, meta, site_scores)


# ---------------------------------------------------------------------------
# Precursor m/z computation
# ---------------------------------------------------------------------------


def _sum_sequence_mass(modified_seq: str) -> Optional[float]:
    """Sum residue + modification masses for a modified sequence string.

    Supports ``(UniMod:N)`` and ``[UNIMOD:N]`` bracket formats,
    as well as the ProForma N-term separator ``-`` and the ``?`` that ends
    modifications of unknown position.
    """
    mass = 0.0
    i = 0
    n = len(modified_seq)

    while i < n:
        ch = modified_seq[i]
        if ch in ("(", "["):
            close = ")" if ch == "(" else "]"
            end = modified_seq.index(close, i)
            delta = parse_unimod_delta(modified_seq[i + 1 : end])
            if delta is None:
                return None
            mass += delta
            i = end + 1
        elif ch in "-?":
            i += 1
        else:
            aa_mass = AA_MASS.get(ch.upper())
            if aa_mass is None:
                return None
            mass += aa_mass
            i += 1
    return mass


@lru_cache(maxsize=16384)
def compute_precursor_mz(modified_seq: str, charge: int) -> Optional[float]:
    """Compute theoretical precursor m/z for a modified peptide.

    Supports ``(UniMod:N)``, ProForma ``[UNIMOD:N]``, and N-term ``-``.

    Args:
        modified_seq: Modified sequence, e.g. ``PEPTM(UniMod:35)IDE``
            or ``PEPTM[UNIMOD:35]IDE`` or ``[UNIMOD:1]-PEPTIDE``.
        charge: Precursor charge state.

    Returns:
        Theoretical m/z or ``None`` when unknown modification or charge <= 0.
    """
    if not modified_seq or charge <= 0:
        return None
    mass = _sum_sequence_mass(modified_seq)
    if mass is None:
        return None
    neutral_mass = mass + WATER_MASS
    return (neutral_mass + charge * PROTON_MASS) / charge


__all__ = [
    "AA_MASS",
    "PROTON_MASS",
    "UNIMOD_MASS",
    "WATER_MASS",
    "UNIMOD_MODS",
    "build_proforma",
    "canonical_mod_tag",
    "compute_precursor_mz",
    "format_mass_shift",
    "from_proforma",
    "mass_to_unimod",
    "normalize_peptidoform",
    "parse_unimod_delta",
    "unimod_from_name",
    "unimod_mass",
]
