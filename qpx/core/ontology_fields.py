"""Names a QPX view emits, the subjects of its ontology rows."""

from __future__ import annotations

from qpx.core.data.schema import ValidationIssue
from qpx.core.sql import sql_build

_NESTED_NAME_QUERIES = {
    "additional_scores": "SELECT DISTINCT unnest(additional_scores).score_name AS name FROM $relation",
    "cv_params": "SELECT DISTINCT unnest(cv_params).cv_name AS name FROM $relation",
    "additional_intensities": (
        "SELECT DISTINCT unnest(pairs).intensity_name AS name FROM "
        "(SELECT unnest(additional_intensities).intensities AS pairs FROM $relation)"
    ),
}


def emitted_field_names(con, relation: str) -> set[str]:
    """Return the top-level columns of a relation and the score, CV and intensity names its rows use.

    ``con`` is anything with a DuckDB-style ``execute`` method, and ``relation`` a
    table, view or ``read_parquet`` expression.
    """
    columns = {row[0] for row in con.execute(sql_build("DESCRIBE SELECT * FROM $relation", relation=relation)).fetchall()}
    names = set(columns)
    for column in sorted(columns & set(_NESTED_NAME_QUERIES)):
        rows = con.execute(sql_build(_NESTED_NAME_QUERIES[column], relation=relation)).fetchall()
        names.update(row[0] for row in rows if row[0] is not None)
    return names


def ontology_field_issues(engine, present: set[str]) -> list[ValidationIssue]:
    """Flag ontology rows naming a field that their view does not emit."""
    rows = engine.execute("SELECT field_name, view FROM ontology").fetchall()
    issues = []
    for view in sorted({row[1] for row in rows} & present - {"ontology"}):
        unknown = sorted({row[0] for row in rows if row[1] == view} - emitted_field_names(engine, view))
        if unknown:
            message = f"{len(unknown)} ontology row(s) of view {view!r} name no field it emits: {', '.join(unknown[:5])}"
            issues.append(
                ValidationIssue(
                    structure="ontology", check="ontology_unknown_field", severity="warning", column=None, message=message
                )
            )
    return issues
