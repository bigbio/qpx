"""Opt-in semantic checks of the PSM and feature views: peptidoforms, modification positions and masses."""

from __future__ import annotations

import math
import re

import pyarrow as pa

from qpx.converters.ptm import AA_MASS, PROTON_MASS, WATER_MASS, unimod_mass
from qpx.core.data.validation import ValidationIssue, _issue
from qpx.core.sql import sql_build

# calculated_mz may differ from the peptidoform's theoretical m/z by rounding
# (float32, four-decimal mass shifts), not by a modification.
_CALCULATED_MZ_PPM = 20.0
# A wider window for measured precursors.
_MASS_ERROR_PPM = 50.0
# Canonical ProForma of the specification (docs/spec/peptidoform.md): optional
# tags of unknown position and "?", optional N-terminal tags and a hyphen,
# residues each followed by tags, optional hyphen and C-terminal tags.
_TAGS = r"(?:\[[^\[\]]+\])"
_CANONICAL = re.compile(rf"(?:{_TAGS}+\?)?(?:{_TAGS}+-)?(?:[A-Z]{_TAGS}*)+(?:-{_TAGS}+)?")
_TAG = re.compile(r"\[([^\[\]]+)\]")
_MASS_SHIFT = re.compile(r"[+-]\d+(?:\.\d+)?")
_UNIMOD_TAG = re.compile(r"UNIMOD:(\d+)(?:\|([+-]\d+(?:\.\d+)?))?")
_SEMANTIC_TABLE = "_qpx_semantic_peptidoforms"

_POSITION_SQL = """
    SELECT count(*)
    FROM (
        SELECT sequence, UNNEST(m.positions) AS p
        FROM (SELECT sequence, UNNEST(modifications) AS m FROM ($query))
    )
    WHERE p.position < 0
       OR p.position > len(sequence) + 1
       OR (p.position BETWEEN 1 AND len(sequence)
           AND p.amino_acid IS NOT NULL
           AND p.amino_acid <> substr(sequence, p.position, 1))
"""

# Rows whose calculated_mz is off the theoretical m/z of their peptidoform;
# cross-linked PSMs carry their partner's mass and are skipped.
_CALCULATED_MZ_SQL = """
    WITH rows AS (
        SELECT v.peptidoform, v.charge, v.calculated_mz,
               (t.mass + v.charge * $proton) / v.charge AS theoretical_mz
        FROM ($query) v
        JOIN $table t USING (peptidoform)
        WHERE v.calculated_mz IS NOT NULL AND v.charge > 0 $no_cross_links
    )
    SELECT count(*), arg_max(peptidoform, abs(calculated_mz - theoretical_mz)),
           arg_max(charge, abs(calculated_mz - theoretical_mz)),
           arg_max(calculated_mz, abs(calculated_mz - theoretical_mz)),
           arg_max(theoretical_mz, abs(calculated_mz - theoretical_mz))
    FROM rows
    WHERE abs(calculated_mz - theoretical_mz) > theoretical_mz * $ppm * 1e-6
"""

_MASS_ERROR_SQL = """
    SELECT count(*)
    FROM ($query) v
    WHERE v.observed_mz IS NOT NULL AND v.calculated_mz > 0
      AND abs(v.observed_mz - v.calculated_mz) > v.calculated_mz * $ppm * 1e-6 $skip
"""


def _tag_mass(tag: str) -> tuple[bool, float | None]:
    """Whether a tag is canonical (a known UNIMOD accession, a mass shift or a name), and its mass when known."""
    if _MASS_SHIFT.fullmatch(tag):
        return True, float(tag)
    unimod = _UNIMOD_TAG.fullmatch(tag)
    if unimod:
        if unimod.group(2):
            return True, float(unimod.group(2))
        mass = unimod_mass(int(unimod.group(1)))
        return mass is not None, mass
    # Another accession namespace (PSI-MOD, PSI-MS, CHEMMOD, mixed-case UniMod)
    # is not canonical; a plain name is kept as written and has no mass.
    return ":" not in tag, None


def peptidoform_problem(peptidoform: str, sequence: str | None) -> tuple[str | None, float | None]:
    """The first problem of a peptidoform ("malformed", "sequence" or "accession"), and its neutral mass when known."""
    if not _CANONICAL.fullmatch(peptidoform):
        return "malformed", None
    if re.sub(rf"{_TAGS}|[-?]", "", peptidoform) != sequence:
        return "sequence", None
    mass = WATER_MASS + sum(AA_MASS.get(residue, math.nan) for residue in sequence)
    for tag in _TAG.findall(peptidoform):
        canonical, delta = _tag_mass(tag)
        if not canonical:
            return "accession", None
        mass += math.nan if delta is None else delta
    return None, None if math.isnan(mass) else mass


_MESSAGES = {
    "malformed": ("peptidoform_malformed", "are not canonical ProForma"),
    "sequence": ("peptidoform_sequence", "do not match their sequence once the modifications are removed"),
    "accession": ("modification_accession", "carry a tag that is neither a known UNIMOD accession, a mass shift nor a name"),
}


def _peptidoform_issues(engine, query: str, structure: str, severity: str) -> tuple[list[ValidationIssue], pa.Table]:
    """Issues of the distinct peptidoforms, and their neutral masses where every part has one."""
    pairs = engine.connection.sql(sql_build("SELECT DISTINCT peptidoform, sequence FROM ($query)", query=query)).fetchall()
    found: dict[str, list[str]] = {}
    masses: dict[str, float] = {}
    for peptidoform, sequence in pairs:
        if peptidoform is None:
            continue
        problem, mass = peptidoform_problem(peptidoform, sequence)
        if problem:
            found.setdefault(problem, []).append(peptidoform)
        elif mass is not None:
            masses[peptidoform] = mass
    issues = []
    for problem, peptidoforms in found.items():
        check, text = _MESSAGES[problem]
        message = f"{len(peptidoforms)} distinct peptidoform(s) {text}, e.g. {peptidoforms[0]!r}"
        issues.append(_issue(structure, check, severity, "peptidoform", message))
    table = pa.table({"peptidoform": list(masses), "mass": list(masses.values())})
    return issues, table


def _mass_issues(engine, query: str, columns: list[str], structure: str, severity: str) -> list[ValidationIssue]:
    no_cross_links = "AND v.cross_links IS NULL" if "cross_links" in columns else ""
    count, peptidoform, charge, calculated, theoretical = engine.connection.sql(
        sql_build(
            _CALCULATED_MZ_SQL,
            query=query,
            table=_SEMANTIC_TABLE,
            proton=repr(PROTON_MASS),
            ppm=repr(_CALCULATED_MZ_PPM),
            no_cross_links=no_cross_links,
        )
    ).fetchall()[0]
    issues = []
    if count:
        message = (
            f"{count} row(s) have a calculated_mz more than {_CALCULATED_MZ_PPM:g} ppm from the theoretical m/z of "
            f"their peptidoform, e.g. {peptidoform!r} at charge {charge}: {calculated:.4f} vs {theoretical:.4f}"
        )
        issues.append(_issue(structure, "calculated_mz_mismatch", severity, "calculated_mz", message))
    # A multiplet row's m/z values describe its lightest partner, whatever partner was measured.
    skip = "AND v.identified_labels IS NULL" if "identified_labels" in columns else ""
    sql = sql_build(_MASS_ERROR_SQL, query=query, ppm=repr(_MASS_ERROR_PPM), skip=f"{skip} {no_cross_links}")
    (off,) = engine.connection.sql(sql).fetchall()[0]
    if off:
        message = f"{off} row(s) have an observed_mz more than {_MASS_ERROR_PPM:g} ppm from calculated_mz"
        issues.append(_issue(structure, "mass_error", "warning", "observed_mz", message))
    return issues


def semantic_issues(engine, query: str, columns: list[str], structure: str, severity: str) -> list[ValidationIssue]:
    """Peptidoform syntax, sequence and accessions; modification positions; calculated and observed m/z."""
    if structure not in ("psm", "feature") or not {"peptidoform", "sequence"} <= set(columns):
        return []
    issues, masses = _peptidoform_issues(engine, query, structure, severity)
    if "modifications" in columns:
        (misplaced,) = engine.connection.sql(sql_build(_POSITION_SQL, query=query)).fetchall()[0]
        if misplaced:
            message = f"{misplaced} modification position(s) lie outside the peptide or name another residue"
            issues.append(_issue(structure, "modification_position", severity, "modifications", message))
    if {"charge", "calculated_mz", "observed_mz"} <= set(columns):
        engine.connection.register(_SEMANTIC_TABLE, masses)
        try:
            issues.extend(_mass_issues(engine, query, columns, structure, severity))
        finally:
            engine.connection.unregister(_SEMANTIC_TABLE)
    return issues
