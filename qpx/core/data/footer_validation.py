"""Checks of the Parquet footer contract: specification version, file type and identity keys."""

from __future__ import annotations

from pathlib import Path

from qpx.core.data.identity import SET_VALUED_FIELDS, derive_composite_id
from qpx.core.data.validation import ValidationIssue, _issue
from qpx.core.sql import sql_build, validate_identifier
from qpx.version import QPX_SPEC_VERSION, parse_spec_version

_READER_VERSION = parse_spec_version(QPX_SPEC_VERSION)
_CHECK_IDENTITY_BATCH_ROWS = 100_000

_DUPLICATE_COMPOSITE_SQL = """
    SELECT count(*), sum(ids)
    FROM (
        SELECT count(DISTINCT $id_field) AS ids
        FROM ($query)
        GROUP BY $keys
        HAVING count(DISTINCT $id_field) > 1
    )
"""


def _version_problem(raw: bytes | None, severity: str) -> tuple[str, str, str] | None:
    """(check, severity, message) for a footer qpx_version this reader does not support."""
    if raw is None:
        return "footer_version", severity, "has no qpx_version in its footer"
    text = raw.decode(errors="replace")
    version = parse_spec_version(text)
    if version is None:
        return "footer_version", "error", f"has qpx_version {text!r}, which is not a <major>.<minor> specification version"
    if version[0] != _READER_VERSION[0]:
        return "footer_version", "error", f"has qpx_version {text}; this reader supports specification {_READER_VERSION[0]}.x"
    if version > _READER_VERSION:
        message = f"has qpx_version {text}, newer than this reader's {QPX_SPEC_VERSION}; columns added since are not checked"
        return "footer_version", "warning", message
    return None


def _file_type_problem(view, raw: bytes | None, severity: str) -> tuple[str, str, str] | None:
    if raw is None:
        return "footer_file_type", severity, "has no file_type in its footer"
    text = raw.decode(errors="replace")
    if text != view.file_type:
        return "footer_file_type", "error", f"has file_type {text!r}; a {view.view_name} file is {view.file_type!r}"
    return None


def _identity_problems(view, meta: dict, severity: str) -> list[tuple[str, str, str]]:
    """Missing identity keys, a composite naming fields outside the schema, or another primary key."""
    if view.identity_composite is None:
        return []
    problems = []
    for key in (b"identity_composite", b"primary_key"):
        if key not in meta:
            problems.append(
                ("footer_identity", severity, f"has no {key.decode()} in its footer, so its ids cannot be re-derived")
            )
    composite = meta.get(b"identity_composite")
    if composite is not None:
        known = set(view.get_arrow_schema().names)
        unknown = [name for name in composite.decode(errors="replace").split(",") if name not in known]
        if unknown:
            message = f"declares identity_composite fields outside the {view.view_name} schema: {unknown}"
            problems.append(("footer_identity", "error", message))
    primary_key = meta.get(b"primary_key")
    if primary_key is not None and primary_key.decode(errors="replace") != view.primary_key[0]:
        message = (
            f"declares primary_key {primary_key.decode(errors='replace')!r}; the {view.view_name} key is {view.primary_key[0]!r}"
        )
        problems.append(("footer_identity", "error", message))
    return problems


def footer_issues(view, footers: list[tuple[str, dict]], severity: str) -> list[ValidationIssue]:
    """Footer contract problems of the files backing a view, one issue per distinct problem."""
    files_by_problem: dict[tuple[str, str, str], list[str]] = {}
    for source, meta in footers:
        problems = [
            _version_problem(meta.get(b"qpx_version"), severity),
            _file_type_problem(view, meta.get(b"file_type"), severity),
            *_identity_problems(view, meta, severity),
        ]
        for problem in filter(None, problems):
            files_by_problem.setdefault(problem, []).append(Path(source).name)
    issues = []
    for (check, level, text), names in files_by_problem.items():
        where = names[0] if len(names) == 1 else f"{len(names)} files (e.g. {names[0]})"
        issues.append(_issue(view.view_name, check, level, None, f"{where} {text}"))
    return issues


def mixed_identity_issues(view, footers: list[tuple[str, dict]]) -> list[ValidationIssue]:
    """Files of one view whose ids were derived differently, so equal rows do not have equal ids."""
    issues = []
    for key in (b"identity_composite", b"identity_hash_version"):
        declared = sorted({meta[key].decode(errors="replace") for _, meta in footers if key in meta})
        if len(declared) > 1:
            message = f"The files declare different {key.decode()} values ({', '.join(declared)}): their ids are not comparable"
            issues.append(_issue(view.view_name, "mixed_identity", "warning", view.primary_key[0], message))
    return issues


def stale_identity_issues(engine, query: str, view, columns: list[str], footers, severity: str) -> list[ValidationIssue]:
    """Ids that no longer match the composite their footers declare them derived from."""
    composite = _declared_composite(view, footers, columns) if view.identity_composite else None
    if composite is None:
        return []
    id_field = view.primary_key[0]
    if any(meta.get(b"identity_source", b"derived") != b"derived" for _, meta in footers):
        message = f"Not all {id_field} values are derived (footer identity_source), so they were not re-derived"
        return [_issue(view.view_name, "identity_not_checked", "warning", id_field, message)]
    selected = ", ".join(validate_identifier(name) for name in (id_field, *composite))
    reader = engine.connection.sql(sql_build("SELECT $selected FROM ($query)", selected=selected, query=query))
    rows = stale = 0
    for batch in reader.to_arrow_reader(_CHECK_IDENTITY_BATCH_ROWS):
        for row in batch.to_pylist():
            rows += 1
            stale += row[id_field] != derive_composite_id(row, composite)
    if not stale:
        return []
    message = f"{stale} of {rows} {id_field} values do not match the declared composite ({', '.join(composite)})"
    return [_issue(view.view_name, "stale_identity", severity, id_field, message)]


def _declared_composite(view, footers: list[tuple[str, dict]], columns: list[str]) -> tuple[str, ...] | None:
    """The identity composite the files declare (the schema default when none does), if its columns are present."""
    declared = {meta[b"identity_composite"].decode(errors="replace") for _, meta in footers if b"identity_composite" in meta}
    if len(declared) > 1:
        return None
    composite = tuple(declared.pop().split(",")) if declared else view.identity_composite
    if not composite or not set(composite) <= set(columns) or view.primary_key[0] not in columns:
        return None
    return composite


def duplicate_composite_issues(engine, query: str, view, columns: list[str], footers, severity: str) -> list[ValidationIssue]:
    """Identity composites that carry several ids: one composite value must name one row identity."""
    composite = _declared_composite(view, footers, columns) if view.identity_composite else None
    if composite is None:
        return []
    keys = ", ".join(
        f"list_sort(list_distinct({validate_identifier(name)}))" if name in SET_VALUED_FIELDS else validate_identifier(name)
        for name in composite
    )
    id_field = view.primary_key[0]
    sql = sql_build(_DUPLICATE_COMPOSITE_SQL, id_field=validate_identifier(id_field), query=query, keys=keys)
    groups, ids = engine.connection.sql(sql).fetchall()[0]
    if not groups:
        return []
    message = (
        f"{groups} identity composite value(s) ({', '.join(composite)}) carry {ids} different {id_field} values; "
        f"one composite value must name one {view.view_name} row identity"
    )
    return [_issue(view.view_name, "duplicate_composite", severity, id_field, message)]
