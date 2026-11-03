"""Row-level value checks, site and proteoform link queries, and validation result types."""

from __future__ import annotations

from dataclasses import dataclass, field

import pyarrow as pa
import pyarrow.compute as pc

# ---------------------------------------------------------------------------
# Validation result types
# ---------------------------------------------------------------------------


@dataclass
class ValidationIssue:
    """A single validation issue found during schema validation."""

    structure: str
    check: str  # "missing_column" | "type_mismatch" | "null_values" | "duplicate_pk"
    severity: str  # "error" | "warning"
    column: str | None
    message: str


@dataclass
class ValidationResult:
    """Result of validating a structure against its schema."""

    structure: str
    issues: list[ValidationIssue] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        """True if no errors (warnings are acceptable)."""
        return not any(i.severity == "error" for i in self.issues)

    @property
    def errors(self) -> list[ValidationIssue]:
        """The issues of severity error."""
        return [i for i in self.issues if i.severity == "error"]

    @property
    def warnings(self) -> list[ValidationIssue]:
        """The issues of severity warning."""
        return [i for i in self.issues if i.severity == "warning"]

    @property
    def summary(self) -> str:
        """The structure name, VALID or INVALID, and the error and warning counts."""
        n_err = len(self.errors)
        n_warn = len(self.warnings)
        status = "VALID" if self.is_valid else "INVALID"
        parts = []
        if n_err:
            parts.append(f"{n_err} error{'s' if n_err != 1 else ''}")
        if n_warn:
            parts.append(f"{n_warn} warning{'s' if n_warn != 1 else ''}")
        detail = f" ({', '.join(parts)})" if parts else ""
        return f"{self.structure}: {status}{detail}"


# Unit-interval columns of the site and proteoform views: q-values,
# probabilities and evidence fractions.
_UNIT_INTERVAL_COLUMNS = {
    "site": ("site_qvalue", "modified_probability", "evidence_fraction"),
    "proteoform": ("proteoform_qvalue",),
}
# Tolerance for a residue's evidence fractions summing to 1.
_FRACTION_SUM_TOLERANCE = 1e-6
# Producers round candidate probabilities (MaxQuant to three decimals), so an
# alternative group can sum to a few thousandths above 1.
_CANDIDATE_SUM_TOLERANCE = 0.01
# Closed vocabularies of the site and proteoform views. ``prior_evidence`` and
# ``source`` are open: producers name their own evidence classes and sources.
_VOCABULARIES = {
    "proteoform": {"kind": ("canonical", "isoform", "variant", "chain"), "evidence": ("specific", "shared")},
    "xic": {"transition_type": ("ms1", "fragment")},
}


def _out_of_range_count(values: pa.ChunkedArray | pa.Array) -> int:
    """Non-null values that are not finite numbers in [0, 1]; 0 for non-float input."""
    if not pa.types.is_floating(values.type):
        return 0
    bad = pc.or_(
        pc.or_(pc.less(values, 0.0), pc.greater(values, 1.0)),
        pc.invert(pc.is_finite(values)),
    )
    return pc.sum(pc.cast(pc.fill_null(bad, False), pa.int64())).as_py() or 0


def _issue(structure: str, check: str, severity: str, column: str, message: str) -> ValidationIssue:
    return ValidationIssue(structure=structure, check=check, severity=severity, column=column, message=message)


def _unit_interval_issues(counts: dict[str, int], structure: str, severity: str) -> list[ValidationIssue]:
    """Flag q-values, probabilities and evidence fractions outside [0, 1] or not finite."""
    issues = []
    for name, bad in counts.items():
        if bad:
            message = f"{bad} value(s) of {name} are not finite numbers in [0, 1]"
            issues.append(_issue(structure, "value_out_of_range", severity, name, message))
    return issues


def _alternative_fractions(table: pa.Table) -> list[list | None] | None:
    """Each site row's alternative evidence fractions, or None when the column is absent or malformed."""
    if "alternatives" not in table.schema.names:
        return None
    alternatives = table.column("alternatives")
    if not (pa.types.is_list(alternatives.type) and pa.types.is_struct(alternatives.type.value_type)):
        return None
    index = alternatives.type.value_type.get_field_index("evidence_fraction")
    if index < 0 or not pa.types.is_floating(alternatives.type.value_type.field(index).type):
        return None
    return [None if alts is None else [a and a.get("evidence_fraction") for a in alts] for alts in alternatives.to_pylist()]


def _evidence_fraction_counts(table: pa.Table | pa.RecordBatch) -> tuple[int, int, int]:
    """Alternatives out of [0, 1], rows not summing to 1, and rows listed out of order.

    A residue's evidence is one pool: a row's fraction and its alternatives' sum to 1.
    Without alternatives, the row holds the whole pool, so its fraction must be 1
    (or null). Alternatives are listed largest first and are each in [0, 1].
    """
    fractions = _alternative_fractions(table)
    if fractions is None or "evidence_fraction" not in table.schema.names:
        return 0, 0, 0
    if not pa.types.is_floating(table.column("evidence_fraction").type):
        return 0, 0, 0
    own = table.column("evidence_fraction").to_pylist()
    out_of_range = unsummed = unordered = 0
    for row_own, alts in zip(own, fractions):
        alts = [a for a in (alts or []) if a is not None]
        out_of_range += sum(not (0.0 <= a <= 1.0) for a in alts)
        unordered += any(x < y for x, y in zip(alts, alts[1:]))
        if row_own is not None and abs(row_own + sum(alts) - 1.0) > _FRACTION_SUM_TOLERANCE:
            unsummed += 1
    return out_of_range, unsummed, unordered


def _evidence_fraction_issues(counts: tuple[int, int, int], structure: str, severity: str) -> list[ValidationIssue]:
    out_of_range, unsummed, unordered = counts
    issues = []
    if out_of_range:
        message = f"{out_of_range} alternative evidence fraction(s) are not finite numbers in [0, 1]"
        issues.append(_issue(structure, "value_out_of_range", severity, "alternatives", message))
    if unsummed:
        message = (
            f"{unsummed} row(s) have an evidence fraction and alternatives that do not sum to 1 "
            "(a row without alternatives holds the whole residue's evidence: fraction 1)"
        )
        issues.append(_issue(structure, "evidence_fraction_sum", severity, "evidence_fraction", message))
    if unordered:
        message = f"{unordered} row(s) list alternatives that are not in decreasing evidence fraction"
        issues.append(_issue(structure, "alternatives_order", severity, "alternatives", message))
    return issues


def _vocabulary_count(values: pa.Array | pa.ChunkedArray, allowed: tuple[str, ...]) -> int:
    """Non-null values outside a closed vocabulary; 0 for non-string input."""
    if not (pa.types.is_string(values.type) or pa.types.is_large_string(values.type)):
        return 0
    outside = pc.invert(pc.is_in(values, value_set=pa.array(allowed)))
    return pc.sum(pc.cast(pc.and_(outside, pc.is_valid(values)), pa.int64())).as_py() or 0


def _chromatogram_length_counts(table: pa.Table | pa.RecordBatch) -> int:
    """XIC rows whose parallel ``rt`` and ``intensity`` arrays differ in length.

    The two arrays are positionally aligned, so a mismatch means the trace cannot
    be plotted or integrated. Null lists are reported by the non-nullable check.
    """
    names = set(table.schema.names)
    if not {"rt", "intensity"} <= names:
        return 0
    rt = table.column("rt")
    intensity = table.column("intensity")
    if not (pa.types.is_list(rt.type) and pa.types.is_list(intensity.type)):
        return 0
    mismatch = pc.not_equal(pc.list_value_length(rt), pc.list_value_length(intensity))
    return pc.sum(pc.cast(pc.fill_null(mismatch, False), pa.int64())).as_py() or 0


def _chromatogram_length_issues(count: int, structure: str, severity: str) -> list[ValidationIssue]:
    if not count:
        return []
    message = f"{count} row(s) have rt and intensity arrays of different length"
    return [_issue(structure, "chromatogram_length_mismatch", severity, "rt", message)]


def _vocabulary_issues(counts: dict[str, int], structure: str, severity: str) -> list[ValidationIssue]:
    """Flag values outside a closed vocabulary (null is allowed)."""
    issues = []
    for name, bad in counts.items():
        if bad:
            message = f"{bad} value(s) of {name} are not one of {', '.join(_VOCABULARIES[structure][name])}"
            issues.append(_issue(structure, "value_not_in_vocabulary", severity, name, message))
    return issues


# Columns holding sequence edits, per view, and the keys one sequence is
# identified by within a row (pepmap: one per protein mapping).
_SEQUENCE_EDIT_COLUMNS = {"pepmap": ("sequence_edits", ("protein_accession",)), "proteoform": ("edits", ())}


def _is_float(t: pa.DataType) -> bool:
    return pa.types.is_floating(t)


def _is_int(t: pa.DataType) -> bool:
    return pa.types.is_integer(t)


def _is_str(t: pa.DataType) -> bool:
    return pa.types.is_string(t) or pa.types.is_large_string(t)


def _flattened_struct_column(table: pa.Table | pa.RecordBatch, column: str, fields: dict) -> pa.Table | None:
    """One row per list element of a list<struct> column, with its parent row index.

    ``fields`` maps each needed struct field to a type predicate. None when the
    column is absent or empty, not a list of structs, or a field is missing or of
    another type: the type check reports that, and semantic checks must not run
    arithmetic or string kernels on it.
    """
    if column not in table.schema.names or len(table) == 0:
        return None
    col = table.column(column)
    if not (pa.types.is_list(col.type) and pa.types.is_struct(col.type.value_type)):
        return None
    struct = col.type.value_type
    for name, ok in fields.items():
        index = struct.get_field_index(name)
        if index < 0 or not ok(struct.field(index).type):
            return None
    flat = pc.list_flatten(col)
    return pa.table({"row": pc.list_parent_indices(col), **{f: pc.struct_field(flat, f) for f in fields}})


def _modification_candidate_counts(table: pa.Table | pa.RecordBatch) -> tuple[int, int]:
    """Candidate probabilities out of [0, 1], and alternative groups whose probabilities sum above 1.

    Candidates without a group may co-occur (a label and an acetyl on one
    residue), so no sum is implied across them.
    """
    flat = _flattened_struct_column(
        table,
        "modification_candidates",
        {"position": _is_int, "probability": _is_float, "alternative_group": _is_int},
    )
    if flat is None or flat.num_rows == 0:
        return 0, 0
    bad = _out_of_range_count(flat.column("probability"))
    over = 0
    grouped = flat.filter(pc.is_valid(flat.column("alternative_group")))
    if grouped.num_rows:
        sums = grouped.group_by(["row", "alternative_group"]).aggregate([("probability", "sum")])
        over = pc.greater(sums.column("probability_sum"), 1.0 + _CANDIDATE_SUM_TOLERANCE)
        over = pc.sum(pc.cast(pc.fill_null(over, False), pa.int64())).as_py() or 0
    return bad, over


def _modification_candidate_issues(counts: tuple[int, int], structure: str, severity: str) -> list[ValidationIssue]:
    """Candidate probabilities are in [0, 1]; within an alternative group they sum to at most 1."""
    bad, over = counts
    issues = []
    if bad:
        message = f"{bad} modification candidate probabilities are not finite numbers in [0, 1]"
        issues.append(_issue(structure, "value_out_of_range", severity, "modification_candidates", message))
    if over:
        message = f"{over} alternative group(s) have candidate probabilities that sum to more than 1.01"
        issues.append(_issue(structure, "candidate_probability_sum", severity, "modification_candidates", message))
    return issues


def _malformed_edit_count(flat: pa.Table) -> int:
    """Edits that do not name a real reference span and residues, or change nothing."""
    start, end = flat.column("start"), flat.column("end")
    ref, alt = flat.column("reference"), flat.column("alternative")
    span = pc.add(pc.subtract(end, start), 1)
    letters = r"^[A-Z]*$"
    bad = pc.or_kleene(
        pc.or_kleene(pc.less(start, 1), pc.not_equal(pc.utf8_length(ref), span)),
        pc.or_kleene(
            pc.or_kleene(
                pc.invert(pc.match_substring_regex(ref, letters)),
                pc.invert(pc.match_substring_regex(alt, letters)),
            ),
            pc.equal(ref, alt),
        ),
    )
    return pc.sum(pc.cast(pc.fill_null(bad, True), pa.int64())).as_py() or 0


def _conflicting_edit_count(flat: pa.Table, keys: tuple[str, ...]) -> int:
    """Sequences whose edits are out of order or overlap, so they define no single sequence.

    An edit covers the reference interval ``[start, end]``; an insertion before
    residue ``s`` (``end = s - 1``) is the point ``s - 1/2``. Within one sequence
    (a row, or a row and protein mapping), each edit must start strictly after
    the previous one ends: two edits of one residue, overlapping spans and two
    insertions at one point all conflict.
    """
    previous_end: dict[tuple, float] = {}
    conflicting: set[tuple] = set()
    for record in flat.select(["row", *keys, "start", "end"]).to_pylist():
        start, end = record["start"], record["end"]
        if start is None or end is None:
            continue
        key = tuple(record[k] for k in ("row", *keys))
        low, high = (start - 0.5, start - 0.5) if end == start - 1 else (start, end)
        if key in previous_end and low <= previous_end[key]:
            conflicting.add(key)
        previous_end[key] = max(high, previous_end.get(key, high))
    return len(conflicting)


def _sequence_edit_counts(table: pa.Table | pa.RecordBatch, structure: str) -> tuple[int, int]:
    """Malformed sequence edits, and sequences whose edits overlap or are out of order.

    ``reference`` is the span ``start..end`` (1-based, inclusive) of the reference
    sequence, so its length is ``end - start + 1``; an insertion has an empty
    reference and ``end = start - 1``. Residues are one-letter codes, and an edit
    changes something.
    """
    column, keys = _SEQUENCE_EDIT_COLUMNS.get(structure, (None, ()))
    if column is None:
        return 0, 0
    fields = {"start": _is_int, "end": _is_int, "reference": _is_str, "alternative": _is_str}
    fields.update({k: _is_str for k in keys})
    flat = _flattened_struct_column(table, column, fields)
    if flat is None or flat.num_rows == 0:
        return 0, 0
    return _malformed_edit_count(flat), _conflicting_edit_count(flat, keys)


def _sequence_edit_issues(counts: tuple[int, int], structure: str, severity: str) -> list[ValidationIssue]:
    """Each sequence edit names a real reference span and residues, and the edits of
    one sequence are ordered and do not overlap.
    """
    malformed, conflicting = counts
    column = _SEQUENCE_EDIT_COLUMNS.get(structure, (None, ()))[0]
    issues = []
    if malformed:
        message = (
            f"{malformed} sequence edit(s) are malformed: start must be >= 1, reference must be the "
            "end - start + 1 residues of the span (empty with end = start - 1 for an insertion), "
            "residues must be one-letter codes, and an edit must change something"
        )
        issues.append(_issue(structure, "malformed_sequence_edit", severity, column, message))
    if conflicting:
        message = (
            f"{conflicting} sequence(s) have edits that overlap or are out of order, so they do not define "
            "one sequence: each edit must start after the previous one ends"
        )
        issues.append(_issue(structure, "conflicting_sequence_edits", severity, column, message))
    return issues


def _anchor_without_membership(table: pa.Table | pa.RecordBatch) -> int:
    """Count rows that carry an ``anchor_protein`` but no ``pg_accessions`` membership.

    The protein group (``pg_accessions``) is the analyte and the feature->pg join
    key (see docs/spec/pg.md -> Protein group semantics); ``anchor_protein`` is only
    a descriptive representative. A row that has been protein-mapped (non-null
    ``anchor_protein``) MUST therefore carry its full group membership, otherwise it
    cannot be attributed to a group. Applies to any view exposing both columns
    (``feature``, where both are nullable; ``pg`` already requires both).
    """
    names = table.schema.names
    if "anchor_protein" not in names or "pg_accessions" not in names or len(table) == 0:
        return 0
    anchor = table.column("anchor_protein")
    pg = table.column("pg_accessions")
    # A non-list pg_accessions is a type mismatch already reported by the type check
    # above; skip here rather than let list_value_length raise on malformed input.
    if not (pa.types.is_list(pg.type) or pa.types.is_large_list(pg.type)):
        return 0
    # Non-string anchor_protein is a type mismatch reported elsewhere.
    if not (pa.types.is_string(anchor.type) or pa.types.is_large_string(anchor.type)):
        return 0
    # list_value_length is null where pg_accessions is null; fill_null(0) folds a
    # null membership into "empty" so the boolean logic never propagates nulls.
    pg_len = pc.fill_null(pc.list_value_length(pg), 0)
    pg_missing = pc.equal(pg_len, 0)
    # Blank/whitespace anchors are unset (not a membership violation).
    anchor_set = pc.greater(pc.utf8_length(pc.utf8_trim_whitespace(pc.fill_null(anchor, ""))), 0)
    violation = pc.and_(anchor_set, pg_missing)
    return pc.sum(pc.cast(violation, pa.int64())).as_py() or 0


def _anchor_membership_issues(bad: int, structure: str, severity: str) -> list[ValidationIssue]:
    if not bad:
        return []
    message = (
        f"{bad} row(s) have anchor_protein set but no pg_accessions "
        "(protein-group membership). The group is the analyte and the "
        "feature->pg join key, so a row with a leading protein must carry "
        "its full pg_accessions membership."
    )
    return [_issue(structure, "anchor_without_membership", severity, "pg_accessions", message)]


# Link checks of the site and proteoform views, as constant queries.
_LINK_SQL = {
    "site": {
        "dangling_pg": (
            'SELECT count(*) FROM "site" v WHERE v.pg_id IS NOT NULL AND NOT EXISTS (SELECT 1 FROM pg p WHERE p.pg_id = v.pg_id)'
        ),
        "pg_context": (
            'SELECT count(*) FROM "site" v JOIN pg p ON p.pg_id = v.pg_id '
            "WHERE p.label IS NOT NULL "
            "OR list_sort(list_distinct(v.grouped_runs)) <> list_sort(list_distinct(p.grouped_runs))"
        ),
        "pg_membership": (
            'SELECT count(*) FROM "site" v JOIN pg p ON p.pg_id = v.pg_id '
            "WHERE NOT list_contains(p.pg_accessions, v.protein_accession)"
        ),
        "dangling_psm": (
            'SELECT count(*) FROM "site" v WHERE v.best_psm_id IS NOT NULL '
            "AND NOT EXISTS (SELECT 1 FROM psm s WHERE s.psm_id = v.best_psm_id)"
        ),
    },
    "proteoform": {
        "dangling_pg": (
            'SELECT count(*) FROM "proteoform" v WHERE v.pg_id IS NOT NULL '
            "AND NOT EXISTS (SELECT 1 FROM pg p WHERE p.pg_id = v.pg_id)"
        ),
        "pg_context": (
            'SELECT count(*) FROM "proteoform" v JOIN pg p ON p.pg_id = v.pg_id '
            "WHERE NOT (v.label IS NOT DISTINCT FROM p.label) "
            "OR list_sort(list_distinct(v.grouped_runs)) <> list_sort(list_distinct(p.grouped_runs))"
        ),
        "pg_membership": (
            'SELECT count(*) FROM "proteoform" v JOIN pg p ON p.pg_id = v.pg_id '
            "WHERE NOT list_contains(p.pg_accessions, v.family)"
        ),
        "dangling_psm": (
            'SELECT count(*) FROM "proteoform" v WHERE v.best_psm_id IS NOT NULL '
            "AND NOT EXISTS (SELECT 1 FROM psm s WHERE s.psm_id = v.best_psm_id)"
        ),
    },
}
_LINK_CONTEXT_MESSAGE = {
    "site": "link to a pg row that has a label or other grouped_runs (a site links to the null-label row of its runs)",
    "proteoform": "link to a pg row whose label or grouped_runs differ from the row's",
}
_LINK_MEMBERSHIP_MESSAGE = {
    "site": "link to a pg row whose pg_accessions do not hold the protein_accession",
    "proteoform": "link to a pg row whose pg_accessions do not hold the family",
}
