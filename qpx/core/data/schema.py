"""Schema base classes for QPX data structures.

ViewSchema instances are loaded from YAML specs at import time.
FieldDef is a simple container for a resolved Arrow field.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import pyarrow as pa

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
        return [i for i in self.issues if i.severity == "error"]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [i for i in self.issues if i.severity == "warning"]

    @property
    def summary(self) -> str:
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


def _pg_referential_issues(
    table: pa.Table,
    structure: str,
    severity: str,
) -> list[ValidationIssue]:
    """Protein-group referential invariants, stated at the measurement level.

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
    names = set(table.schema.names)
    if not {"pg_accessions", "grouped_runs", "label"} <= names or len(table) == 0:
        return []

    import duckdb

    con = duckdb.connect()
    try:
        con.register("pg_validate", table)
        issues: list[ValidationIssue] = []
        dup_lists = con.execute(
            "SELECT COUNT(*) FROM pg_validate WHERE grouped_runs IS NOT NULL "
            "AND length(grouped_runs) <> length(list_distinct(grouped_runs))"
        ).fetchone()[0]
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
        double = con.execute(
            """
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
        ).fetchone()[0]
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
    finally:
        con.close()


def _primary_key_issues(
    table: pa.Table,
    primary_key: tuple[str, ...],
    structure: str,
    severity: str,
) -> list[ValidationIssue]:
    """Return a duplicate-key issue when every primary-key column is available."""
    columns = [column for column in primary_key if column in table.schema.names]
    if len(columns) != len(primary_key) or len(table) == 0:
        return []

    key_table = pa.table({column: _canonicalize_primary_key_column(table.column(column)) for column in columns})
    total_count = len(key_table)
    unique_count = key_table.group_by(columns).aggregate([]).num_rows
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


class ViewSchema:
    """Schema for a QPX view, loaded from YAML."""

    def __init__(
        self,
        view_name: str,
        file_type: str,
        primary_key: list[str],
        fields: dict[str, FieldDef],
        doc: str = "",
        extra_columns: bool = False,
        identity_composite: list[str] | None = None,
    ):
        self._view_name = view_name
        self._file_type = file_type
        self._primary_key = tuple(primary_key)
        self._fields = fields
        self._doc = doc
        self._extra_columns = extra_columns
        self._identity_composite = tuple(identity_composite) if identity_composite else None
        self._arrow_schema_cache: pa.Schema | None = None

    @property
    def primary_key(self) -> tuple[str, ...]:
        """The view's primary-key column(s)."""
        return self._primary_key

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
        result = ValidationResult(structure=self._view_name)
        gated_severity = "error" if strict else "warning"
        expected = self.get_arrow_schema()

        # Check 1 & 2: Column presence and type matching
        for expected_field in expected:
            if expected_field.name not in table.schema.names:
                if not self._is_optional(expected_field.name):
                    result.issues.append(
                        ValidationIssue(
                            structure=self._view_name,
                            check="missing_column",
                            severity="error",
                            column=expected_field.name,
                            message=f"Missing required column: {expected_field.name}",
                        )
                    )
                continue
            actual = table.schema.field(expected_field.name)
            if not _types_compatible(expected_field.type, actual.type):
                result.issues.append(
                    ValidationIssue(
                        structure=self._view_name,
                        check="type_mismatch",
                        severity="error",
                        column=expected_field.name,
                        message=(f"Column '{expected_field.name}': expected {expected_field.type}, got {actual.type}"),
                    )
                )

        # Check 3: Null values in non-nullable columns
        for field_name, fdef in self._fields.items():
            if field_name not in table.schema.names:
                continue
            if not fdef.nullable:
                col = table.column(field_name)
                null_count = col.null_count
                if null_count > 0:
                    result.issues.append(
                        ValidationIssue(
                            structure=self._view_name,
                            check="null_values",
                            severity=gated_severity,
                            column=field_name,
                            message=(
                                f"Column '{field_name}' is non-nullable but contains "
                                f"{null_count} null value(s) out of {len(col)} rows"
                            ),
                        )
                    )

        result.issues.extend(
            _primary_key_issues(
                table,
                self._primary_key,
                self._view_name,
                gated_severity,
            )
        )

        if self._view_name == "pg":
            result.issues.extend(_pg_referential_issues(table, self._view_name, gated_severity))

        return result

    def _is_optional(self, field_name: str) -> bool:
        fdef = self._fields.get(field_name)
        return fdef is not None and fdef.optional
