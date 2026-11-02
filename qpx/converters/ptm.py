"""Shared PTM utilities: Unimod lookups, ProForma builder/parser, precursor m/z.

This flat module consolidates shared logic previously spread across
``qpx.converters.ptm.unimod``, ``qpx.converters.ptm.proforma``,
``qpx.converters.ptm.modifications``, and ``qpx.converters.ptm.__init__``.

Unimod accessions, names and masses come from pyopenms' ``ModificationsDB``
(the Unimod release pyopenms ships), not from hand-written tables. Tool-specific
parsing (MaxQuant short names such as ``ox``) stays in the converter's
``constants.py``.
"""

from __future__ import annotations

import logging
import math
import re
from bisect import bisect_left, bisect_right
from collections.abc import Iterator, Mapping
from functools import lru_cache
from typing import NamedTuple, Optional, Tuple

# ---------------------------------------------------------------------------
# Unimod lookups (pyopenms ModificationsDB)
# ---------------------------------------------------------------------------

logger = logging.getLogger(__name__)

_MASS_TOL = 0.01

# OpenMS adds its own entries to ModificationsDB from ``CHEMISTRY/custom_mods.xml``
# (CUSTOM0-9, PhosphoDecoy, ...) under record ids from 99900 up, which Unimod
# leaves free. They are not Unimod accessions.
_OPENMS_CUSTOM_RECORD_ID = 99900

# Equal masses at the same site (same composition) are ambiguous. A tie goes to
# the modification here when exactly one of the tied entries is listed; any other
# tie stays unresolved and the mass is written. Listed are modifications search
# engines are configured with, against isobaric amino-acid substitutions,
# artefacts and cross-linker remnants: Acetyl 1 (Ser->Glu), Amidated 2,
# Carbamidomethyl 4 (Gly 1263, Carbofuran), Carbamyl 5, Carboxymethyl 6,
# Deamidated 7 (Asn->Asp 621, Gln->Glu 632), Methyl 34, Oxidation 35 (Phe->Tyr),
# Dimethyl 36 (Ethyl 280, Delta:H(4)C(2) 255), Trimethyl 37 (Propyl 1305),
# Propionyl 58, GG 121 (Dicarbamidomethyl 1290), Formyl 122, iTRAQ4plex 214
# (mTRAQ:13C(3)15N(1) 889), mTRAQ 888 (Xlink:DMP[140]), Butyryl 1289
# (Crotonaldehyde 253), Gluratylation 1848 (Xlink:BS2G[114]), Lactylation 2114
# (Carboxyethyl 378) and Hexanoyl 2121 (MesitylOxide 1873).
_PREFERRED_ON_TIE = frozenset({1, 2, 4, 5, 6, 7, 34, 35, 36, 37, 58, 121, 122, 214, 888, 1289, 1848, 2114, 2121})

# Vendor modification names that are neither Unimod names nor ModificationsDB synonyms.
_NAME_ALIASES: dict[str, str] = {
    "glygly": "GG",
    "tmt10plex": "TMT6plex",
    "tmt11plex": "TMT6plex",
    "tmt16plex": "TMTpro",
    "tmt18plex": "TMTpro",
}

# MaxQuant/Spectronaut "Name (Sites)" suffix, e.g. " (STY)" or " (Protein N-term)".
_SITE_SUFFIX_RE = re.compile(r"\s+\([^()]*\)$")
_ACCESSION_RE = re.compile(r"(?i)unimod:(\d+)")
_SIGNED_MASS_RE = re.compile(r"[+-]\d+(?:\.\d+)?")


class _UnimodSite(NamedTuple):
    """One site of a Unimod modification: ``residue`` is ``"X"`` for any residue, ``terminus`` "N", "C" or None."""

    accession: int
    mass: float
    residue: str
    terminus: Optional[str]


def _modification_names(mod) -> set[str]:
    """Lower-cased names of a ``ResidueModification``: Unimod title, full id, full name and synonyms."""
    names = {mod.getId(), mod.getFullId(), mod.getFullName()}
    names.update(name.decode() if isinstance(name, bytes) else name for name in mod.getSynonyms())
    return {name.lower() for name in names if name}


class _UnimodIndex(NamedTuple):
    """The Unimod entries of pyopenms' ``ModificationsDB``, sorted by delta mass and keyed by name."""

    sites: list[_UnimodSite]
    masses: list[float]
    by_name: dict[str, list[_UnimodSite]]
    mass_by_accession: dict[int, float]


@lru_cache(maxsize=1)
def _unimod_index() -> _UnimodIndex:
    """Build the Unimod index once per process.

    ModificationsDB holds one entry per Unimod accession and site, read from the
    Unimod release that pyopenms ships. Entries without a Unimod accession (from
    PSI-MOD or XLMOD only) and OpenMS' own entries (record id 99900 and up) are
    left out.
    """
    import pyopenms as oms  # pylint: disable=import-outside-toplevel

    spec = oms.ResidueModification.TermSpecificity
    termini = {int(spec.N_TERM): "N", int(spec.PROTEIN_N_TERM): "N", int(spec.C_TERM): "C", int(spec.PROTEIN_C_TERM): "C"}
    db = oms.ModificationsDB()
    sites: list[_UnimodSite] = []
    by_name: dict[str, list[_UnimodSite]] = {}
    for index in range(db.getNumberOfModifications()):
        mod = db.getModification(index)
        if 0 < mod.getUniModRecordId() < _OPENMS_CUSTOM_RECORD_ID:
            terminus = termini.get(int(mod.getTermSpecificity()))
            sites.append(_UnimodSite(mod.getUniModRecordId(), mod.getDiffMonoMass(), mod.getOrigin(), terminus))
            for name in _modification_names(mod):
                by_name.setdefault(name, []).append(sites[-1])
    sites.sort(key=lambda site: site.mass)
    return _UnimodIndex(
        sites=sites,
        masses=[site.mass for site in sites],
        by_name=by_name,
        mass_by_accession={site.accession: site.mass for site in sites},
    )


class _UnimodMasses(Mapping):
    """Read-only ``{Unimod accession: monoisotopic delta mass}`` view of ``ModificationsDB``."""

    def __getitem__(self, accession: int) -> float:
        return _unimod_index().mass_by_accession[accession]

    def __iter__(self) -> Iterator[int]:
        return iter(_unimod_index().mass_by_accession)

    def __len__(self) -> int:
        return len(_unimod_index().mass_by_accession)


UNIMOD_MASS: Mapping[int, float] = _UnimodMasses()


def _site_rank(site: _UnimodSite, residue: Optional[str], terminus: Optional[str]) -> Optional[int]:
    """How well a Unimod site fits a modification at *residue* / *terminus*; ``None`` if it does not.

    0: same placement (residue, or the N-/C-terminus) and the site names the residue;
    1: same placement and the site allows any residue; 2: the site names the
    residue but has the other placement, as when a tool reports pyro-Glu, which
    Unimod defines on the N-terminal Q, on the Q itself. Without a residue and a
    terminus every site fits; without a residue, every site of the placement.
    """
    if residue is None:
        return 0 if terminus is None or site.terminus == terminus else None
    same_placement = site.terminus == terminus
    if site.residue == residue:
        return 0 if same_placement else 2
    return 1 if site.residue == "X" and same_placement else None


def _break_tie(accessions: set[int]) -> Optional[int]:
    """The accession of equally good candidates: the only one, or the only one in ``_PREFERRED_ON_TIE``."""
    if len(accessions) == 1:
        return next(iter(accessions))
    preferred = accessions & _PREFERRED_ON_TIE
    return next(iter(preferred)) if len(preferred) == 1 else None


@lru_cache(maxsize=8192)
def mass_to_unimod(mass: float, residue: Optional[str] = None, terminus: Optional[str] = None) -> Optional[int]:
    """Find the Unimod accession of a delta *mass* at a residue or terminus.

    Args:
        mass: Monoisotopic delta mass in Da.
        residue: One-letter code of the modified residue; for a terminal
            modification the terminal residue. ``None`` when unknown.
        terminus: ``"N"`` or ``"C"`` for a modification of the peptide
            terminus, ``None`` for a residue modification.

    The candidates are the Unimod sites within ``_MASS_TOL`` Da that fit the site
    (see :func:`_site_rank`). The closest mass wins, then the better site fit.
    Entries left equal have the same composition at the same site; they resolve
    only to a modification in ``_PREFERRED_ON_TIE`` (iTRAQ4plex 214, not
    mTRAQ:13C(3)15N(1) 889; Butyryl 1289, not Crotonaldehyde 253). Returns
    ``None`` if no site matches, the tie stays open (Isovaleryl, 2-methylbutyryl
    and Pentanoyl on K) or *mass* is not finite.
    """
    if not math.isfinite(mass):
        return None
    index = _unimod_index()
    low = bisect_left(index.masses, mass - _MASS_TOL)
    high = bisect_right(index.masses, mass + _MASS_TOL)
    keys: dict[int, tuple[float, int]] = {}
    for site in index.sites[low:high]:
        rank = _site_rank(site, residue, terminus)
        if rank is not None:
            key = (round(abs(site.mass - mass), 6), rank)
            keys[site.accession] = min(key, keys.get(site.accession, key))
    if not keys:
        return None
    best = min(keys.values())
    return _break_tie({accession for accession, key in keys.items() if key == best})


def format_mass_delta(mass: float) -> str:
    """Signed ProForma mass tag for a delta mass without a Unimod entry, e.g. ``+999.1234``."""
    return f"+{mass}" if mass >= 0 else str(mass)


def mass_tag(mass: float, residue: Optional[str] = None, terminus: Optional[str] = None) -> str:
    """ProForma tag for a delta mass: ``UNIMOD:n`` (see :func:`mass_to_unimod`), else the signed mass."""
    accession = mass_to_unimod(mass, residue, terminus)
    return f"UNIMOD:{accession}" if accession is not None else format_mass_delta(mass)


def _name_keys(name: str) -> list[str]:
    """Lookup keys for a modification name: as given, without a vendor site suffix, and its alias."""
    key = name.strip().lower()
    base = _SITE_SUFFIX_RE.sub("", key)
    keys = [key, base]
    if base in _NAME_ALIASES:
        keys.append(_NAME_ALIASES[base].lower())
    return keys


@lru_cache(maxsize=4096)
def name_to_unimod(name: str, residue: Optional[str] = None, terminus: Optional[str] = None) -> Optional[int]:
    """Find the Unimod accession of a modification name.

    Accepts the Unimod names and synonyms ModificationsDB knows (``Phospho``,
    ``Deamidation``, ``Pyro-glu from E``), the vendor form with a site suffix
    (MaxQuant and Spectronaut ``Phospho (STY)``, ``Acetyl (Protein N-term)``) and
    the vendor names in ``_NAME_ALIASES`` (``GlyGly``, ``TMT10plex``,
    ``TMT16plex``). When a name has several accessions, those with a site fitting
    *residue* / *terminus* (see :func:`mass_to_unimod`) win, then the lowest.
    Returns ``None`` for an unknown name.
    """
    by_name = _unimod_index().by_name
    for key in _name_keys(name):
        sites = by_name.get(key)
        if sites:
            ranked = [(_site_rank(site, residue, terminus), site.accession) for site in sites]
            fitting = [entry for entry in ranked if entry[0] is not None]
            return min(fitting)[1] if fitting else min(site.accession for site in sites)
    return None


@lru_cache(maxsize=1)
def _names_by_accession() -> dict[int, frozenset[str]]:
    names: dict[int, set[str]] = {}
    for name, sites in _unimod_index().by_name.items():
        for site in sites:
            names.setdefault(site.accession, set()).add(name)
    return {accession: frozenset(found) for accession, found in names.items()}


def unimod_names(accession: int) -> frozenset[str]:
    """Lower-cased Unimod title, full name and synonyms of a Unimod *accession*; empty when it is unknown."""
    return _names_by_accession().get(accession, frozenset())


@lru_cache(maxsize=4096)
def modification_tag(name: str, residue: Optional[str] = None, terminus: Optional[str] = None) -> str:
    """ProForma tag for a vendor modification name at *residue* / *terminus*.

    ``UNIMOD:n`` for a Unimod accession (``UniMod:35``) or name (see
    :func:`name_to_unimod`); a signed mass (``+15.9949``) is resolved like a
    FragPipe mass (see :func:`mass_tag`). Any other name is kept, with square
    brackets replaced by parentheses so it cannot end the ProForma tag early,
    and logged once.
    """
    match = _ACCESSION_RE.fullmatch(name.strip())
    if match:
        return f"UNIMOD:{int(match.group(1))}"
    accession = name_to_unimod(name, residue, terminus)
    if accession is not None:
        return f"UNIMOD:{accession}"
    if _SIGNED_MASS_RE.fullmatch(name.strip()):
        return mass_tag(float(name), residue, terminus)
    logger.warning("Modification %r has no Unimod entry; it is kept as a name in the peptidoform", name)
    return name.replace("[", "(").replace("]", ")")


def _group_end(text: str, start: int, open_char: str, close_char: str) -> Optional[int]:
    """Index after the delimiter closing the one at *start*, counting nested pairs; ``None`` if unclosed."""
    depth = 0
    for index in range(start, len(text)):
        if text[index] == open_char:
            depth += 1
        elif text[index] == close_char:
            depth -= 1
            if depth == 0:
                return index + 1
    return None


def split_modified_sequence(modified: str, open_char: str = "(", close_char: str = ")") -> tuple[str, list[tuple[int, str]]]:
    """Split a vendor modified sequence into its residues and ``(position, modification)`` pairs.

    A modification gets the position of the residue before it, or 0 before the
    first residue (the :func:`build_proforma` convention). Delimiters nest, so
    MaxQuant ``M(Oxidation (M))`` gives ``Oxidation (M)``; an unclosed
    delimiter takes the rest of the string.
    """
    residues: list[str] = []
    mods: list[tuple[int, str]] = []
    i = 0
    while i < len(modified):
        if modified[i] != open_char:
            residues.append(modified[i])
            i += 1
            continue
        end = _group_end(modified, i, open_char, close_char) or len(modified) + 1
        mods.append((len(residues), modified[i + 1 : end - 1]))
        i = end
    return "".join(residues), mods


def position_site(sequence: str, position: int) -> tuple[Optional[str], Optional[str]]:
    """Residue and terminus of a :func:`build_proforma` position in *sequence*.

    Position 0 is the N-terminus and ``len(sequence) + 1`` the C-terminus; both
    return the terminal residue with ``"N"`` / ``"C"``.
    """
    if position <= 0:
        return (sequence[:1] or None), "N"
    if position > len(sequence):
        return (sequence[-1:] or None), "C"
    return sequence[position - 1], None


def vendor_position(name: str, position: int, sequence: str) -> int:
    """:func:`build_proforma` position of a MaxQuant or Spectronaut modification called *name*.

    Both tools name the site (``Acetyl (Protein N-term)``, ``Amidated (C-term)``)
    and write a C-terminal modification after the last residue, where
    :func:`split_modified_sequence` gives it that residue's position. A name
    for the N-terminus gets position 0, one for the C-terminus
    ``len(sequence) + 1``; any other name keeps *position*.
    """
    if "N-term" in name:
        return 0
    if "C-term" in name:
        return len(sequence) + 1
    return position


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
        - ``UniMod:35`` / ``UNIMOD:35`` -> Unimod delta mass (``UNIMOD_MASS``)
        - ``+15.9949`` / ``-17.0265`` -> parse directly
    """
    m = re.match(r"(?:UniMod|UNIMOD):(\d+)", mod_str, re.IGNORECASE)
    if m:
        uid = int(m.group(1))
        return UNIMOD_MASS.get(uid)
    try:
        return float(mod_str)
    except ValueError:
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
    so that :func:`_from_proforma_impl` can parse them. Only parentheses
    outside square brackets are converted: a ProForma tag such as
    ``K[GlyGly (K)]`` is copied as it is instead of becoming the nested
    ``K[GlyGly [K]]``.
    """
    if "(" not in peptidoform:
        return peptidoform
    peptidoform = peptidoform.removeprefix(".")
    out: list[str] = []
    i = 0
    while i < len(peptidoform):
        char = peptidoform[i]
        if char not in "([":
            out.append(char)
            i += 1
            continue
        end = _group_end(peptidoform, i, char, ")" if char == "(" else "]")
        if end is None:
            out.append(peptidoform[i:])
            break
        # Convert the outer delimiters; keep the inner content as it is
        out.append(f"[{peptidoform[i + 1 : end - 1]}]")
        i = end
    return _add_nterm_dash("".join(out))


def _add_nterm_dash(peptidoform: str) -> str:
    """Insert the ProForma ``-`` after the leading N-terminal tags when it is missing."""
    end = 0
    while end < len(peptidoform) and peptidoform[end] == "[":
        end = _group_end(peptidoform, end, "[", "]") or len(peptidoform)
    if 0 < end < len(peptidoform) and peptidoform[end] != "-":
        return f"{peptidoform[:end]}-{peptidoform[end:]}"
    return peptidoform


def _modification_site(seq_pos: int, last_aa: str | None, cterm: bool, sequence_length: int) -> tuple[int, str | None]:
    """Return the QPX position and residue for a ProForma modification."""
    if seq_pos == 0:
        return 0, None
    if cterm:
        return sequence_length + 1, None
    return seq_pos, last_aa


def _modification_identity(mod_str: str, aa: str | None, meta: Optional[dict]) -> tuple[str, str | None]:
    """Resolve a modification tag to its name and accession."""
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
            else:
                name = f"CHEMMOD:{mod_str}"
        else:
            name = f"CHEMMOD:{mod_str}"
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

    i = 0
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
    as well as the ProForma N-term separator ``-``.
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
        elif ch == "-":
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
    "build_proforma",
    "compute_precursor_mz",
    "format_mass_delta",
    "from_proforma",
    "mass_tag",
    "mass_to_unimod",
    "modification_tag",
    "name_to_unimod",
    "parse_unimod_delta",
    "position_site",
    "split_modified_sequence",
    "vendor_position",
]
