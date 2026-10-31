"""Names a QPX view emits, the subjects of its ontology rows."""

from __future__ import annotations

_NESTED_NAME_QUERIES = {
    "additional_scores": "SELECT DISTINCT unnest(additional_scores).score_name AS name FROM {relation}",
    "cv_params": "SELECT DISTINCT unnest(cv_params).cv_name AS name FROM {relation}",
    "additional_intensities": (
        "SELECT DISTINCT unnest(pairs).intensity_name AS name FROM "
        "(SELECT unnest(additional_intensities).intensities AS pairs FROM {relation})"
    ),
}


def emitted_field_names(con, relation: str) -> set[str]:
    """Return the top-level columns of a relation and the score, CV and intensity names its rows use.

    ``con`` is anything with a DuckDB-style ``execute`` method, and ``relation`` a
    table, view or ``read_parquet`` expression.
    """
    columns = {row[0] for row in con.execute(f"DESCRIBE SELECT * FROM {relation}").fetchall()}
    names = set(columns)
    for column in sorted(columns & set(_NESTED_NAME_QUERIES)):
        rows = con.execute(_NESTED_NAME_QUERIES[column].format(relation=relation)).fetchall()
        names.update(row[0] for row in rows if row[0] is not None)
    return names
