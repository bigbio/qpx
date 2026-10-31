"""MaxQuant converter constants — PTM parsing and label maps."""

from __future__ import annotations

import csv
import re
from functools import lru_cache
from typing import Optional

from qpx.converters.localization import parse_site_probabilities
from qpx.converters.ptm import (
    NTERM,
    build_proforma,
    compute_precursor_mz,
    modification_tag,
    name_to_unimod,
    position_site,
    split_modified_sequence,
    vendor_position,
)

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


FixedMods = tuple[tuple[str, str], ...]


def _site_list(sites: str) -> list[str]:
    """Residues and ``N-term`` named by a MaxQuant or SDRF site field; protein N-term sites are skipped."""
    text = sites.strip()
    if text.lower() in ("n-term", "any n-term"):
        return [NTERM]
    if "term" in text.lower():
        return []
    return [aa for aa in text.replace(",", "") if aa.isalpha() and aa.isupper()]


def sdrf_fixed_mods(entries: Optional[list[dict]]) -> FixedMods:
    """Fixed modifications from the SDRF ``comment[modification parameters]`` entries.

    *entries* is :meth:`SdrfConverter.parse_modification_params` output.
    """
    fixed: list[tuple[str, str]] = []
    for entry in entries or []:
        if not entry["fixed"]:
            continue
        accession = re.match(r"(?i)unimod:(\d+)$", entry["accession"] or "")
        uid = int(accession.group(1)) if accession else name_to_unimod(entry["name"] or "")
        if uid is None:
            continue
        position = (entry["position"] or "").strip().lower()
        if position in ("any n-term", "n-term"):
            sites = [NTERM]
        elif "term" in position:
            sites = []
        else:
            sites = _site_list(entry["target_amino_acid"] or "")
        fixed.extend((site, f"UNIMOD:{uid}") for site in sites)
    return tuple(fixed)


def _named_fixed_mods(names: str) -> FixedMods:
    """``(site, "UNIMOD:n")`` pairs of a MaxQuant modification list such as ``Carbamidomethyl (C);Oxidation (M)``."""
    fixed: list[tuple[str, str]] = []
    for name in names.strip().split(";"):
        match = re.match(r"(.+?)\s+\(([^()]*)\)$", name.strip())
        uid = name_to_unimod(name) if match else None
        if uid is not None:
            fixed.extend((site, f"UNIMOD:{uid}") for site in _site_list(match.group(2)))
    return tuple(fixed)


def parameters_fixed_mods(parameters_path) -> FixedMods:
    """Fixed modifications from the ``Fixed modifications`` row of MaxQuant's ``parameters.txt``."""
    try:
        with open(parameters_path, encoding="utf-8") as handle:
            row = next((line for line in handle if line.startswith("Fixed modifications\t")), None)
    except FileNotFoundError:
        return ()
    return _named_fixed_mods(row.split("\t", 1)[1]) if row is not None else ()


# MaxQuant SILAC labels, as named in summary.txt ``Labels<n>``: (residue, UNIMOD accession).
MQ_SILAC_LABELS: dict[str, tuple[str, int]] = {
    "Arg6": ("R", 188),
    "Arg10": ("R", 267),
    "Lys4": ("K", 481),
    "Lys6": ("K", 188),
    "Lys8": ("K", 259),
}

# Multiplet letters of MaxQuant's column names (``Intensity H``), lightest first.
MQ_LABEL_LETTERS = "LMH"


def _summary_runs(summary_path) -> list[dict]:
    """The per-raw-file rows of MaxQuant's ``summary.txt``; empty when the file is absent."""
    try:
        with open(summary_path, encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except FileNotFoundError:
        return []
    # The closing "Total" row has no multiplicity.
    return [row for row in rows if (row.get("Multiplicity") or "").isdigit()]


def summary_fixed_mods(summary_path) -> FixedMods:
    """Fixed modifications from MaxQuant's ``summary.txt``; MaxQuant 2 leaves them out of ``parameters.txt``."""
    names = next((row["Fixed modifications"] for row in _summary_runs(summary_path) if row.get("Fixed modifications")), "")
    return _named_fixed_mods(names)


def summary_label_mods(summary_path) -> dict[str, tuple[FixedMods, ...]]:
    """Label modifications of each labeling state, per raw file, from MaxQuant's ``summary.txt``.

    ``Labels<n>`` names the labels of state *n* (``Arg6;Lys6``); a state without
    labels has none. Empty when ``summary.txt`` is absent.
    """
    states: dict[str, tuple[FixedMods, ...]] = {}
    for row in _summary_runs(summary_path):
        per_state = []
        for state in range(int(row["Multiplicity"])):
            names = [name for name in (row.get(f"Labels{state}") or "").split(";") if name]
            unknown = [name for name in names if name not in MQ_SILAC_LABELS]
            if unknown:
                raise ValueError(f"{summary_path}: unsupported MaxQuant label(s) {', '.join(unknown)}")
            per_state.append(tuple((MQ_SILAC_LABELS[name][0], f"UNIMOD:{MQ_SILAC_LABELS[name][1]}") for name in names))
        states[row["Raw file"]] = tuple(per_state)
    return states


def _with_fixed_mods(sequence: str, mods: list[tuple[int, str]], fixed_mods: FixedMods) -> list[tuple[int, str]]:
    """Add the fixed modifications MaxQuant leaves out of ``Modified sequence``.

    A fixed modification goes on every matching site that carries no other modification.
    """
    taken = {position for position, _ in mods}
    extra = [(0, tag) for site, tag in fixed_mods if site == NTERM and 0 not in taken]
    extra += [
        (position, tag)
        for position, aa in enumerate(sequence, start=1)
        if position not in taken
        for site, tag in fixed_mods
        if site == aa
    ]
    return mods + extra


@lru_cache(maxsize=8192)
def to_proforma(modified_sequence: Optional[str], fixed_mods: FixedMods = (), label_mods: FixedMods = ()) -> str:
    """Convert a MaxQuant Modified sequence to ProForma notation.

    MaxQuant wraps sequences in underscores and uses parenthetical
    modification names (possibly nested)::

        _PEPTM(Oxidation (M))IDEK_         -> PEPTM[UNIMOD:35]IDEK
        _(Acetyl (Protein N-term))PEPTIDEK_ -> [UNIMOD:1]-PEPTIDEK
        _(ac)PEPTIDEK_                      -> [UNIMOD:1]-PEPTIDEK
        _PEPK(GlyGly (K))R_                 -> PEPK[UNIMOD:121]R
        _PEPK(Amidated (Protein C-term))_   -> PEPK-[UNIMOD:2]

    *fixed_mods* (``(site, tag)`` pairs) are written inline, since MaxQuant
    omits fixed modifications and isobaric labels from the sequence. The SILAC
    *label_mods* of a labelled PSM go on every matching residue, beside any
    other modification there.
    """
    if not isinstance(modified_sequence, str) or not modified_sequence:
        return ""
    sequence, mods = split_modified_sequence(modified_sequence.strip("_"))
    tags: list[tuple[int, str]] = []
    for position, name in mods:
        position = vendor_position(name, position, sequence)
        tags.append((position, _resolve_mq_mod(name, *position_site(sequence, position))))
    labels = [(position, tag) for position, aa in enumerate(sequence, start=1) for site, tag in label_mods if site == aa]
    return build_proforma(sequence, labels + _with_fixed_mods(sequence, tags, fixed_mods))


def mq_calculated_mz(peptidoform: Optional[str], charge: int, neutral_mass: Optional[float]) -> Optional[float]:
    """Theoretical m/z of *peptidoform*; MaxQuant's ``Mass`` when a modification has no known mass."""
    calculated = compute_precursor_mz(peptidoform, charge) if peptidoform else None
    if calculated is None and neutral_mass and charge:
        calculated = (neutral_mass + charge * PROTON_MASS) / charge
    return calculated


def mq_observed_mz(mz: Optional[float], error_ppm: Optional[float]) -> Optional[float]:
    """Measured precursor m/z, or ``None`` when MaxQuant reports no mass error.

    MaxQuant's ``m/z`` column is the theoretical value of its peptide mass; its
    ``Mass error [ppm]`` is the recalibrated measurement's deviation from it.
    """
    if mz is None or error_ppm is None:
        return None
    return mz * (1 + error_ppm * 1e-6)


def localization_columns(columns) -> list[tuple[str, str]]:
    """The ``<mod> Probabilities`` columns of a MaxQuant table, with their modification names."""
    suffix = " Probabilities"
    return [(column, column.removesuffix(suffix)) for column in columns if column.endswith(suffix)]


def localization_sites(row: dict, columns: list[tuple[str, str]], sequence: str) -> dict[str, dict[int, float]]:
    """Per-modification site probabilities of a msms.txt or evidence.txt row.

    MaxQuant reports one column per variable modification, e.g.
    ``Phospho (STY) Probabilities`` = ``PEPS(0.95)T(0.03)Y(0.02)K``. With several
    phosphates the probabilities sum to their number (PXD071171, PXD072739).
    """
    sites: dict[str, dict[int, float]] = {}
    for column, name in columns:
        for position, probability in parse_site_probabilities(row.get(column), sequence).items():
            sites.setdefault(_resolve_mq_mod(name, *position_site(sequence, position)), {})[position] = probability
    return sites
