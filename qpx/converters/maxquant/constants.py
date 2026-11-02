"""MaxQuant converter constants — PTM parsing and fixed modifications."""

from __future__ import annotations

import csv
import logging
import re
from functools import lru_cache
from typing import Optional

from qpx.converters.localization import parse_site_probabilities
from qpx.converters.ptm import (
    build_proforma,
    compute_precursor_mz,
    modification_tag,
    name_to_unimod,
    position_site,
    split_modified_sequence,
    unimod_names,
    vendor_position,
)

logger = logging.getLogger(__name__)

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


FixedModification = tuple[str, str, str]
"""``(accession, target residues, terminus)``; terminus is ``""``, ``"N"`` or ``"C"``."""

FixedMods = tuple[FixedModification, ...]

LabelMods = tuple[tuple[str, str], ...]
"""``(residue, "UNIMOD:n")`` pairs of the SILAC labels of one labeling state."""

# Unimod names of quantification labels, which MaxQuant configures as labels, not as fixed
# modifications: isobaric tags (TMT, iTRAQ) and the isotopic ``Label:`` modifications of SILAC.
_LABEL_NAME = re.compile(r"(?i)tmt|itraq|^label:")

# SDRF ``PP`` values and MaxQuant site names a fixed modification can be placed for from the
# peptide alone; protein termini need the peptide's protein position, which MaxQuant does not give.
_TERMINUS = {"anywhere": "", "any n-term": "N", "n-term": "N", "any c-term": "C", "c-term": "C"}

_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTUVWY")


def _is_quantification_label(name: Optional[str], uid: Optional[int]) -> bool:
    """Whether an SDRF modification is a TMT, iTRAQ or SILAC label, by its ``NT`` or its Unimod names."""
    if name and _LABEL_NAME.search(name):
        return True
    return uid is not None and any(_LABEL_NAME.search(known) for known in unimod_names(uid))


def _sdrf_target(entry: dict) -> Optional[tuple[str, str]]:
    """``(residues, terminus)`` of an SDRF fixed modification, or ``None`` when it cannot be placed.

    ``TA`` lists the target residues (comma-separated); a terminal ``PP`` without
    a residue (``TA=N-term``, ``TA=X``) applies to any terminal residue.
    """
    terminus = _TERMINUS.get((entry["position"] or "Anywhere").strip().lower())
    tokens = (entry["target_amino_acid"] or "").upper().replace(" ", "").split(",")
    residues = "".join(token for token in tokens if token in _RESIDUES)
    if terminus is None or (terminus == "" and not residues):
        return None
    return residues, terminus


def sdrf_fixed_mods(entries: Optional[list[dict]]) -> FixedMods:
    """Fixed modifications from the SDRF ``comment[modification parameters]`` entries.

    *entries* is :func:`qpx.core.sdrf.parse_modification_parameters` output. Each one
    is taken by its accession (``AC=UNIMOD:n``; without one, the Unimod accession of
    its ``NT`` name), on its target residues or peptide terminus. Quantification
    labels are left out: MaxQuant carries a SILAC label on the labelled rows only,
    and the isobaric label is added from the SDRF channel labels.
    """
    fixed: list[FixedModification] = []
    for entry in entries or []:
        if not entry["fixed"]:
            continue
        accession = re.fullmatch(r"(?i)\s*unimod:(\d+)\s*", entry["accession"] or "")
        uid = int(accession.group(1)) if accession else name_to_unimod((entry["name"] or "").strip())
        if _is_quantification_label(entry["name"], uid):
            continue
        target = _sdrf_target(entry)
        if uid is None or target is None:
            logger.warning("Cannot place SDRF fixed modification %r on MaxQuant peptides; peptidoforms omit it", entry)
        else:
            fixed.append((f"UNIMOD:{uid}", *target))
    return tuple(fixed)


def _named_fixed_mods(names: str) -> FixedMods:
    """Fixed modifications of a MaxQuant modification list such as ``Carbamidomethyl (C);Oxidation (M)``."""
    fixed: list[FixedModification] = []
    for name in names.strip().split(";"):
        match = re.match(r"(.+?)\s+\(([^()]*)\)$", name.strip())
        uid = name_to_unimod(name) if match else None
        if uid is None:
            continue
        site = match.group(2).strip()
        terminus = _TERMINUS.get(site.lower())
        if terminus is not None:
            fixed.append((f"UNIMOD:{uid}", "", terminus))
        elif "term" not in site.lower():
            fixed.append((f"UNIMOD:{uid}", "".join(aa for aa in site if aa in _RESIDUES), ""))
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


def summary_label_mods(summary_path) -> dict[str, tuple[LabelMods, ...]]:
    """Label modifications of each labeling state, per raw file, from MaxQuant's ``summary.txt``.

    ``Labels<n>`` names the labels of state *n* (``Arg6;Lys6``); a state without
    labels has none. Empty when ``summary.txt`` is absent.
    """
    states: dict[str, tuple[LabelMods, ...]] = {}
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


def _fixed_positions(sequence: str, residues: str, terminus: str) -> list[int]:
    """Where a fixed modification applies on ``sequence`` (0 = N-term, ``len + 1`` = C-term)."""
    if not terminus:
        return [position for position, residue in enumerate(sequence, start=1) if residue in residues]
    position, residue = (0, sequence[:1]) if terminus == "N" else (len(sequence) + 1, sequence[-1:])
    return [position] if not residues or residue in residues else []


def fixed_modification_sites(sequence: str, mods: list[tuple[int, str]], fixed_mods: FixedMods) -> list[tuple[int, str]]:
    """The ``(position, accession)`` pairs the fixed modifications add to a peptide.

    MaxQuant leaves fixed modifications out of ``Modified sequence``. Each one
    is placed on every matching site that carries no reported modification,
    so a site never gets a second modification (an N-terminal acetyl keeps a
    fixed N-term label off, as in the search).
    """
    occupied = {position for position, _ in mods}
    sites: list[tuple[int, str]] = []
    for accession, residues, terminus in fixed_mods:
        for position in _fixed_positions(sequence, residues, terminus):
            if position not in occupied:
                occupied.add(position)
                sites.append((position, accession))
    return sites


@lru_cache(maxsize=8192)
def to_proforma(modified_sequence: Optional[str], fixed_mods: FixedMods = (), label_mods: LabelMods = ()) -> str:
    """Convert a MaxQuant Modified sequence to ProForma notation.

    MaxQuant wraps sequences in underscores and uses parenthetical
    modification names (possibly nested)::

        _PEPTM(Oxidation (M))IDEK_         -> PEPTM[UNIMOD:35]IDEK
        _(Acetyl (Protein N-term))PEPTIDEK_ -> [UNIMOD:1]-PEPTIDEK
        _(ac)PEPTIDEK_                      -> [UNIMOD:1]-PEPTIDEK
        _PEPK(GlyGly (K))R_                 -> PEPK[UNIMOD:121]R
        _PEPK(Amidated (Protein C-term))_   -> PEPK-[UNIMOD:2]

    *fixed_mods* are written inline (see :func:`fixed_modification_sites`), since
    MaxQuant omits fixed modifications and isobaric labels from the sequence. The
    SILAC *label_mods* of a labelled PSM go on every matching residue, beside any
    other modification there.
    """
    if not isinstance(modified_sequence, str) or not modified_sequence.strip("_"):
        return ""
    sequence, mods = split_modified_sequence(modified_sequence.strip("_"))
    tags: list[tuple[int, str]] = []
    for position, name in mods:
        position = vendor_position(name, position, sequence)
        tags.append((position, _resolve_mq_mod(name, *position_site(sequence, position))))
    labels = [(position, tag) for position, aa in enumerate(sequence, start=1) for site, tag in label_mods if site == aa]
    return build_proforma(sequence, labels + tags + fixed_modification_sites(sequence, tags, fixed_mods))


def mq_calculated_mz(peptidoform: Optional[str], charge: int, neutral_mass: Optional[float]) -> Optional[float]:
    """Theoretical m/z of *peptidoform*; MaxQuant's ``Mass`` when a modification has no known mass."""
    calculated = compute_precursor_mz(peptidoform, charge) if peptidoform else None
    if calculated is None and neutral_mass and charge:
        calculated = (neutral_mass + charge * PROTON_MASS) / charge
    return calculated


def mq_observed_mz(mz: Optional[float], error_ppm: Optional[float], calculated_mz: Optional[float] = None) -> Optional[float]:
    """Measured precursor m/z, or ``None`` when MaxQuant reports no mass error.

    MaxQuant's ``m/z`` column is the theoretical value of its peptide mass, which
    leaves isobaric labels out; its ``Mass error [ppm]`` is the recalibrated
    measurement's deviation from it. In evidence.txt that deviation includes the
    labels the measured ion carries, so ``m/z`` shifted by it is the measurement.
    msms.txt leaves the labels out of the deviation too; give *calculated_mz*,
    the theoretical m/z of the peptidoform that carries them, to shift instead.
    """
    if mz is None or error_ppm is None:
        return None
    return (mz if calculated_mz is None else calculated_mz) + mz * error_ppm * 1e-6


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
