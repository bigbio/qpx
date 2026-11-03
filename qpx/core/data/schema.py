"""Schema base classes for QPX data structures.

ViewSchema instances are loaded from YAML specs at import time.
FieldDef is a simple container for a resolved Arrow field.
"""

from __future__ import annotations

import json

import pyarrow as pa
import pyarrow.compute as pc

from qpx.core.data.validation import (
    _SEQUENCE_EDIT_COLUMNS,
    _UNIT_INTERVAL_COLUMNS,
    _VOCABULARIES,
    ValidationIssue,
    ValidationResult,
    _anchor_membership_issues,
    _anchor_without_membership,
    _chromatogram_length_counts,
    _chromatogram_length_issues,
    _evidence_fraction_counts,
    _evidence_fraction_issues,
    _issue,
    _modification_candidate_counts,
    _modification_candidate_issues,
    _out_of_range_count,
    _sequence_edit_counts,
    _sequence_edit_issues,
    _unit_interval_issues,
    _vocabulary_count,
    _vocabulary_issues,
)
from qpx.core.sql import sql_build

# ---------------------------------------------------------------------------
# Field and schema definitions
# ---------------------------------------------------------------------------


class FieldDef:
    """A resolved schema field definition."""

    __slots__ = ("name", "arrow_type", "nullable", "optional", "doc")

    def __init__(
        self,
        name: str,
        arrow_type: pa.DataType,
        *,
        nullable: bool = True,
        optional: bool = False,
        doc: str = "",
    ):
        self.name = name
        self.arrow_type = arrow_type
        self.nullable = nullable
        self.optional = optional
        self.doc = doc

    def to_arrow_field(self) -> pa.Field:
        metadata = {"doc": self.doc} if self.doc else None
        return pa.field(self.name, self.arrow_type, nullable=self.nullable, metadata=metadata)


def _types_compatible(expected: pa.DataType, actual: pa.DataType) -> bool:
    """Check type compatibility, tolerating nullability differences in nested structs.

    Parquet does not preserve inner-struct field nullability, so we compare
    types structurally while ignoring the ``nullable`` flag on struct children.
    """
    if expected == actual:
        return True

    # list<X> — compare value types recursively
    if isinstance(expected, pa.ListType) and isinstance(actual, pa.ListType):
        return _types_compatible(expected.value_type, actual.value_type)

    # map<K, V>
    if isinstance(expected, pa.MapType) and isinstance(actual, pa.MapType):
        return _types_compatible(expected.key_type, actual.key_type) and _types_compatible(expected.item_type, actual.item_type)

    # struct — compare field names, types (ignoring nullable)
    if isinstance(expected, pa.StructType) and isinstance(actual, pa.StructType):
        if expected.num_fields != actual.num_fields:
            return False
        for i in range(expected.num_fields):
            ef = expected.field(i)
            af = actual.field(i)
            if ef.name != af.name:
                return False
            if not _types_compatible(ef.type, af.type):
                return False
        return True

    return False


def _canonicalize_primary_key_column(column: pa.ChunkedArray) -> pa.Array | pa.ChunkedArray:
    """Normalize list-valued keys to sorted JSON for stable, collision-safe comparison."""
    if not pa.types.is_list(column.type):
        return column

    values = [
        json.dumps(sorted(value, key=lambda item: (item is None, str(item)))) if value is not None else None
        for value in column.to_pylist()
    ]
    return pa.array(values, type=pa.string())


_SCAN_FORMATS = ("scan", "index", "nativeId")


def _declared_scan_format(schema: pa.Schema, structure: str) -> str | None:
    """The scan_format a PSM view declares in its footer; None when there is nothing to check."""
    if structure != "psm" or "scan" not in schema.names:
        return None
    declared = (schema.metadata or {}).get(b"scan_format")
    return declared.decode(errors="replace") if declared else None


def _scan_format_bounds(scan_format: str) -> tuple[int, int]:
    return (2, 4) if scan_format == "nativeId" else (1, 1)


def _scan_format_violations(scans: pa.Array | pa.ChunkedArray, scan_format: str) -> int:
    """Non-empty scan arrays whose component count the declared format does not allow."""
    if scan_format not in _SCAN_FORMATS or not (pa.types.is_list(scans.type) or pa.types.is_large_list(scans.type)):
        return 0  # A non-list scan is reported by the type check.
    lengths = pc.call_function("list_value_length", [scans])
    minimum, maximum = _scan_format_bounds(scan_format)
    outside = pc.call_function(
        "or",
        [pc.call_function("less", [lengths, pa.scalar(minimum)]), pc.call_function("greater", [lengths, pa.scalar(maximum)])],
    )
    invalid = pc.call_function("and", [pc.call_function("greater", [lengths, pa.scalar(0)]), outside])
    return pc.call_function("sum", [invalid]).as_py() or 0


def _scan_format_issues(scan_format: str | None, count: int, structure: str, severity: str) -> list[ValidationIssue]:
    """Check a declared PSM format without guessing legacy or Feature encodings."""
    if scan_format is None:
        return []
    if scan_format not in _SCAN_FORMATS:
        return [
            _issue(structure, "scan_format", "warning", "scan", f"Unknown scan_format {scan_format!r}; cardinality not checked")
        ]
    if not count:
        return []
    minimum, maximum = _scan_format_bounds(scan_format)
    message = f"{count} non-empty PSM scan array(s) violate scan_format={scan_format!r}: expected {minimum}–{maximum} components"
    return [_issue(structure, "scan_format", severity, "scan", message)]


_PG_RUN_DOUBLE_COUNT_SQL = """
        WITH exploded AS (
            -- Deduplicate each row's own grouped_runs first (a within-row
            -- duplicate is already reported by duplicate_grouped_run); this
            -- way run_double_count measures only cross-row disjointness.
            -- Key on the canonical (order-independent) group MEMBERSHIP, so
            -- distinct groups sharing a leader are not conflated.
            SELECT list_sort(list_distinct(pg_accessions)) AS members, label,
                   UNNEST(list_distinct(grouped_runs)) AS run
            FROM pg_validate
        ),
        repeated AS (
            SELECT COUNT(*) AS c
            FROM exploded
            GROUP BY members, label, run
            HAVING COUNT(*) > 1
        )
        SELECT COALESCE(SUM(c - 1), 0) FROM repeated
        """


def _query_pg_referential_issues(
    con,
    structure: str,
    severity: str,
) -> list[ValidationIssue]:
    """Query protein-group referential invariants from ``pg_validate``.

    Both are intrinsic to the pg table — no run→sample resolution:

    * **duplicate_grouped_run** — ``grouped_runs`` is a set of distinct raw files;
      it must contain no duplicate element.
    * **run_double_count** (run-disjointness) — within one
      ``(pg_accessions, label)`` — the protein-group *membership* that keys the
      pg_id, not merely the leader — the ``grouped_runs`` sets across rows must be
      disjoint, so every raw file contributes to at most one row and no
      measurement is counted twice. Keying on the full membership (canonicalized
      order-independently) matches the pg_id identity, so two genuinely distinct
      groups that share a leading protein are not conflated into a false
      double-count. This is the join-free form of the double-count check.
    """
    issues: list[ValidationIssue] = []
    dup_lists = con.execute(
        "SELECT COUNT(*) FROM pg_validate WHERE grouped_runs IS NOT NULL "
        "AND length(grouped_runs) <> length(list_distinct(grouped_runs))"
    ).fetchall()[0][0]
    if dup_lists:
        issues.append(
            ValidationIssue(
                structure=structure,
                check="duplicate_grouped_run",
                severity=severity,
                column="grouped_runs",
                message=f"{dup_lists} pg row(s) have duplicate raw files within grouped_runs (must be a set of distinct files)",
            )
        )
    double = con.execute(_PG_RUN_DOUBLE_COUNT_SQL).fetchall()[0][0]
    if double:
        issues.append(
            ValidationIssue(
                structure=structure,
                check="run_double_count",
                severity=severity,
                column=None,
                message=(
                    f"{double} run occurrence(s) repeat across pg rows sharing the same "
                    "(pg_accessions, label): grouped_runs sets must be disjoint or protein "
                    "intensity is double-counted"
                ),
            )
        )
    return issues


def _pg_referential_issues(
    table: pa.Table,
    structure: str,
    severity: str,
) -> list[ValidationIssue]:
    """Validate PG referential invariants for an in-memory Arrow table."""
    names = set(table.schema.names)
    if not {"pg_accessions", "grouped_runs", "label"} <= names or len(table) == 0:
        return []

    import duckdb

    con = duckdb.connect()
    try:
        con.register("pg_validate", table)
        return _query_pg_referential_issues(con, structure, severity)
    finally:
        con.close()


def pg_referential_issues_from_parquet(
    path: str,
    structure: str,
    severity: str,
) -> list[ValidationIssue]:
    """Validate PG referential invariants against a file-backed Parquet relation."""
    import duckdb

    con = duckdb.connect()
    try:
        relation = con.read_parquet(path)
        if not {"pg_accessions", "grouped_runs", "label"} <= set(relation.columns):
            return []
        relation.create_view("pg_validate")
        return _query_pg_referential_issues(con, structure, severity)
    finally:
        con.close()


def pg_referential_issues_from_query(con, query: str, structure: str, severity: str) -> list[ValidationIssue]:
    """Validate PG referential invariants over the rows of *query* on *con*."""
    con.execute(sql_build("CREATE OR REPLACE TEMP VIEW pg_validate AS $query", query=query))
    try:
        return _query_pg_referential_issues(con, structure, severity)
    finally:
        con.execute("DROP VIEW pg_validate")


def _primary_key_issues(
    keys: pa.Table | None,
    primary_key: tuple[str, ...],
    structure: str,
    severity: str,
) -> list[ValidationIssue]:
    """Return a duplicate-key issue for the primary-key columns of a whole view.

    *keys* is None when a primary-key column is missing.
    """
    if keys is None or keys.num_rows == 0:
        return []

    key_table = pa.table({column: _canonicalize_primary_key_column(keys.column(column)) for column in primary_key})
    total_count = len(key_table)
    unique_count = key_table.group_by(list(primary_key)).aggregate([]).num_rows
    if unique_count == total_count:
        return []

    duplicate_count = total_count - unique_count
    return [
        ValidationIssue(
            structure=structure,
            check="duplicate_pk",
            severity=severity,
            column=None,
            message=(
                f"Primary key ({', '.join(primary_key)}) has "
                f"{duplicate_count} duplicate row(s) "
                f"({unique_count} unique out of {total_count})"
            ),
        )
    ]


# Columns _RowCounts.add reads besides the required, unit-interval and vocabulary
# ones. Validation streams no other column, so a check that reads a new column
# must list it here.
_ROW_CHECK_COLUMNS = (
    "anchor_protein",
    "pg_accessions",
    "modification_candidates",
    "evidence_fraction",
    "alternatives",
    "scan",
    *(column for column, _ in _SEQUENCE_EDIT_COLUMNS.values()),
)


class _RowCounts:
    """Counts of the checks that read one row at a time, summed over a view's batches."""

    def __init__(self, structure: str, schema: pa.Schema, required: list[str]):
        names = schema.names
        self.structure = structure
        self.nulls = dict.fromkeys(required, 0)
        self.out_of_range = dict.fromkeys((name for name in _UNIT_INTERVAL_COLUMNS.get(structure, ()) if name in names), 0)
        self.vocabulary = dict.fromkeys((name for name in _VOCABULARIES.get(structure, {}) if name in names), 0)
        # Counts of the checks that are not kept per column.
        self.totals = {"anchor": 0, "fractions": (0, 0, 0), "candidates": (0, 0), "edits": (0, 0), "scan": 0, "chromatogram": 0}

    def add(self, batch: pa.Table | pa.RecordBatch, scan_format: str | None) -> None:
        """Add one batch's counts."""
        totals = self.totals
        for name in self.nulls:
            self.nulls[name] += batch.column(name).null_count
        totals["anchor"] += _anchor_without_membership(batch)
        for name in self.out_of_range:
            self.out_of_range[name] += _out_of_range_count(batch.column(name))
        if self.structure == "site":
            totals["fractions"] = tuple(a + b for a, b in zip(totals["fractions"], _evidence_fraction_counts(batch)))
        totals["candidates"] = tuple(a + b for a, b in zip(totals["candidates"], _modification_candidate_counts(batch)))
        totals["edits"] = tuple(a + b for a, b in zip(totals["edits"], _sequence_edit_counts(batch, self.structure)))
        totals["chromatogram"] += _chromatogram_length_counts(batch) if self.structure == "xic" else 0
        for name in self.vocabulary:
            self.vocabulary[name] += _vocabulary_count(batch.column(name), _VOCABULARIES[self.structure][name])
        if scan_format is not None:
            totals["scan"] += _scan_format_violations(batch.column("scan"), scan_format)

    def null_issues(self, rows: int, severity: str) -> list[ValidationIssue]:
        """Issues for the non-nullable columns that hold nulls."""
        return [
            _issue(
                self.structure,
                "null_values",
                severity,
                name,
                f"Column '{name}' is non-nullable but contains {count} null value(s) out of {rows} rows",
            )
            for name, count in self.nulls.items()
            if count
        ]

    def value_issues(self, severity: str) -> list[ValidationIssue]:
        """Issues from the value checks counted over the batches."""
        return [
            *_anchor_membership_issues(self.totals["anchor"], self.structure, severity),
            *_unit_interval_issues(self.out_of_range, self.structure, severity),
            *_evidence_fraction_issues(self.totals["fractions"], self.structure, severity),
            *_vocabulary_issues(self.vocabulary, self.structure, severity),
            *_modification_candidate_issues(self.totals["candidates"], self.structure, severity),
            *_sequence_edit_issues(self.totals["edits"], self.structure, severity),
            *_chromatogram_length_issues(self.totals["chromatogram"], self.structure, severity),
        ]


class ViewSchema:
    """Schema for a QPX view, loaded from YAML."""

    def __init__(
        self,
        view_name: str,
        file_type: str,
        primary_key: list[str],
        fields: dict[str, FieldDef],
        extra_columns: bool = False,
        identity_composite: list[str] | None = None,
    ):
        self._view_name = view_name
        self._file_type = file_type
        self._primary_key = tuple(primary_key)
        self._fields = fields
        self._extra_columns = extra_columns
        self._identity_composite = tuple(identity_composite) if identity_composite else None
        self._arrow_schema_cache: pa.Schema | None = None

    @property
    def primary_key(self) -> tuple[str, ...]:
        """The view's primary-key column(s)."""
        return self._primary_key

    @property
    def view_name(self) -> str:
        """The schema view name."""
        return self._view_name

    @property
    def identity_composite(self) -> tuple[str, ...] | None:
        """Composite of existing fields the single-column identity id is derived
        from, as declared by the schema's ``identity_composite`` key. ``None``
        for views without a derived identity (mz, run, sample, ...)."""
        return self._identity_composite

    @property
    def __name__(self) -> str:
        """For compatibility with parametrized tests and repr."""
        parts = self._view_name.split("_")
        return "".join(p.capitalize() for p in parts) + "Schema"

    def __repr__(self) -> str:
        return f"ViewSchema({self._view_name!r}, fields={len(self._fields)})"

    def get_arrow_schema(self) -> pa.Schema:
        """Generate PyArrow schema from field definitions."""
        if self._arrow_schema_cache is None:
            arrow_fields = [fdef.to_arrow_field() for fdef in self._fields.values()]
            self._arrow_schema_cache = pa.schema(arrow_fields)
        return self._arrow_schema_cache

    def get_extended_arrow_schema(self, extra_column_names: list[str]) -> pa.Schema:
        """Return the base schema extended with additional string columns.

        Only allowed when ``extra_columns=True`` in the YAML spec.
        Extra columns are appended as nullable ``pa.string()`` fields.
        """
        if not self._extra_columns:
            raise ValueError(f"Schema '{self._view_name}' does not allow extra columns")
        base = self.get_arrow_schema()
        known = set(base.names)
        extra_fields = [pa.field(name, pa.string(), nullable=True) for name in extra_column_names if name not in known]
        if not extra_fields:
            return base
        return pa.schema(list(base) + extra_fields)

    def _iter_fields(self):
        """Yield (name, FieldDef) pairs in definition order."""
        yield from self._fields.items()

    # --- Validation --------------------------------------------------------

    def validation_columns(self, schema: pa.Schema) -> list[str]:
        """The columns of *schema* that the row-level checks of validate_batches() read."""
        needed = {name for name, fdef in self._fields.items() if not fdef.nullable}
        needed.update(_ROW_CHECK_COLUMNS, _UNIT_INTERVAL_COLUMNS.get(self._view_name, ()), _VOCABULARIES.get(self._view_name, {}))
        return [name for name in schema.names if name in needed]

    def validate(self, table: pa.Table, *, strict: bool = False) -> list[str]:
        """Validate an Arrow table against this schema. Returns list of error strings.

        For structured results, use validate_full() instead.
        """
        result = self.validate_full(table, strict=strict)
        return [issue.message for issue in result.errors]

    def validate_full(self, table: pa.Table, *, strict: bool = False) -> ValidationResult:
        """Full schema validation returning structured results.

        Checks:
            1. Required columns are present
            2. Column types match the schema
            3. Non-nullable columns contain no null values
            4. Primary key is unique
            5. Non-empty PSM scans match a declared scan_format

        Args:
            table: The Arrow table to validate.
            strict: When True (used by ``qpxc validate`` and CI), nulls in
                non-nullable columns and duplicate primary keys are reported as
                ``error`` (so ``is_valid`` is False and the CLI exits non-zero).
                Defaults to False so that *writers* and *converters*, which
                persist source data as-produced, are not blocked by legitimate
                duplicate keys (e.g. a DIA-NN feature whose stripped-sequence PK
                repeats across modiforms) — strictness is an audit concern, not
                a write-time one. Missing columns and type mismatches are always
                errors regardless of ``strict``.
        """
        severity = "error" if strict else "warning"
        keys = table.select(list(self._primary_key)) if set(self._primary_key) <= set(table.schema.names) else None
        referential = _pg_referential_issues(table, self._view_name, severity) if self._view_name == "pg" else []
        return self.validate_batches(table.schema, [table], keys, referential, strict=strict)

    def validate_batches(
        self,
        schema: pa.Schema,
        batches,
        keys: pa.Table | None,
        referential: list[ValidationIssue],
        *,
        strict: bool = False,
    ) -> ValidationResult:
        """Validate a view read as Arrow batches sharing *schema*, with the checks of validate_full().

        The checks that read one row at a time are summed batch by batch, so the
        whole view is never held at once. The two that span rows are passed in for
        the whole view: *keys* holds its primary-key columns (None when one is
        missing) and *referential* its pg referential issues.
        """
        result = ValidationResult(structure=self._view_name)
        severity = "error" if strict else "warning"
        result.issues.extend(self._column_issues(schema))

        required = [name for name, fdef in self._fields.items() if name in schema.names and not fdef.nullable]
        counts = _RowCounts(self._view_name, schema, required)
        scan_format = _declared_scan_format(schema, self._view_name)
        rows = 0
        for batch in batches:
            rows += batch.num_rows
            counts.add(batch, scan_format)

        result.issues.extend(counts.null_issues(rows, severity))
        result.issues.extend(_primary_key_issues(keys, self._primary_key, self._view_name, severity))
        result.issues.extend(referential)
        result.issues.extend(counts.value_issues(severity))
        result.issues.extend(_scan_format_issues(scan_format, counts.totals["scan"], self._view_name, severity))
        return result

    def _column_issues(self, schema: pa.Schema) -> list[ValidationIssue]:
        """Required columns that are missing and columns whose type does not match."""
        issues = []
        for expected_field in self.get_arrow_schema():
            name = expected_field.name
            if name not in schema.names:
                if not self._is_optional(name):
                    issues.append(_issue(self._view_name, "missing_column", "error", name, f"Missing required column: {name}"))
                continue
            actual = schema.field(name)
            if not _types_compatible(expected_field.type, actual.type):
                message = f"Column '{name}': expected {expected_field.type}, got {actual.type}"
                issues.append(_issue(self._view_name, "type_mismatch", "error", name, message))
        return issues

    def _is_optional(self, field_name: str) -> bool:
        fdef = self._fields.get(field_name)
        return fdef is not None and fdef.optional
