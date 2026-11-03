"""Name-value lists stored as structs with one field per name.

A list such as ``additional_scores`` names each of its values. When every row of
a file holds the same names in the same order, as the scores of a DIA-NN feature
do, ``qpxc transform relayout`` stores the column as a struct with one field per
name, so that the values of one name share a Parquet column and its encoding.
The footer key ``qpx.struct_lists`` records each such column's names, and QPX
readers turn the struct back into the list.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import pyarrow as pa

from qpx.core.sql import _SAFE_IDENTIFIER, sql_build, validate_identifier

STRUCT_LISTS_FOOTER_KEY = b"qpx.struct_lists"

# The names of a row's entries, and for scores their direction, as one value per row.
_LAYOUT_SQL = {
    "additional_scores": "list_transform($column, lambda e: struct_pack(name := e.score_name, higher_better := e.higher_better))",
    "cv_params": "list_transform($column, lambda e: struct_pack(name := e.cv_name))",
    "additional_intensities": (
        "list_transform($column, lambda e: struct_pack(name := e.label, "
        "names := list_transform(e.intensities, lambda i: i.intensity_name)))"
    ),
}
STRUCT_LIST_COLUMNS = tuple(_LAYOUT_SQL)

# A value is stored as the first smaller type it reads back from unchanged, else as it is.
_SMALLER = {"VARCHAR": ("FLOAT", "DOUBLE"), "DOUBLE": ("FLOAT",), "FLOAT": ()}
_ARROW_TYPES = {"FLOAT": pa.float32(), "DOUBLE": pa.float64(), "VARCHAR": pa.string()}


@dataclass(frozen=True)
class StructList:
    """A list column with one layout in every row, and how relayout stores it."""

    column: str
    layout: list[dict]
    select: str
    arrow_type: pa.DataType


def _fields(names: list) -> bool:
    """Whether *names* can be the fields of one struct; DuckDB renames fields that differ only in case."""
    return (
        bool(names)
        and all(isinstance(name, str) and _SAFE_IDENTIFIER.match(name) for name in names)
        and len({name.lower() for name in names}) == len(names)
    )


def _representable(column: str, layout: list[dict]) -> bool:
    """Whether every name of *layout* can be a struct field."""
    if not _fields([entry["name"] for entry in layout]):
        return False
    return column != "additional_intensities" or all(_fields(entry["names"] or []) for entry in layout)


def _leaves(column: str, layout: list[dict]) -> list[tuple[tuple[str, ...], str, str]]:
    """The struct path, list-form expression and type of every value of *layout*."""
    name = validate_identifier(column)
    if column == "additional_intensities":
        return [
            ((entry["name"], inner), f"{name}[{i}].intensities[{k}].intensity_value", "FLOAT")
            for i, entry in enumerate(layout, start=1)
            for k, inner in enumerate(entry["names"], start=1)
        ]
    value, value_type = ("score_value", "DOUBLE") if column == "additional_scores" else ("cv_value", "VARCHAR")
    return [((entry["name"],), f"{name}[{i}].{value}", value_type) for i, entry in enumerate(layout, start=1)]


def _unchanged_sql(expression: str, value_type: str, smaller: str) -> str:
    """SQL true when every value of *expression* reads back unchanged from *smaller*; a null row passes."""
    stored = f"TRY_CAST({expression} AS {smaller})"
    if value_type == "VARCHAR":
        stored = f"CAST({stored} AS VARCHAR)"
    return f"bool_and({expression} IS NOT DISTINCT FROM {stored})"


def _stored_types(con, source: str, leaves) -> list[str]:
    """The type each leaf is stored as: the first smaller type that keeps every value."""
    checks = [(index, smaller) for index, (_, _, value_type) in enumerate(leaves) for smaller in _SMALLER[value_type]]
    if not checks:
        return [value_type for _, _, value_type in leaves]
    aggregates = ", ".join(_unchanged_sql(leaves[index][1], leaves[index][2], smaller) for index, smaller in checks)
    row = con.sql(sql_build("SELECT $aggregates FROM read_parquet($source)", aggregates=aggregates, source=source)).fetchall()[0]
    kept = {check for check, unchanged in zip(checks, row) if unchanged}
    return [
        next((smaller for smaller in _SMALLER[value_type] if (index, smaller) in kept), value_type)
        for index, (_, _, value_type) in enumerate(leaves)
    ]


def _struct(fields: list[tuple[tuple[str, ...], str, str]]) -> tuple[str, pa.DataType]:
    """The SQL and Arrow type of a struct of *fields*, nested along their paths."""
    groups: dict[str, list] = {}
    for path, expression, stored in fields:
        groups.setdefault(path[0], []).append((path[1:], expression, stored))
    sql, arrow = [], []
    for name, members in groups.items():
        if members[0][0]:
            member_sql, member_type = _struct(members)
        else:
            _, expression, stored = members[0]
            member_sql, member_type = f"CAST({expression} AS {stored})", _ARROW_TYPES[stored]
        sql.append(f"{validate_identifier(name)} := {member_sql}")
        arrow.append(pa.field(name, member_type))
    return f"struct_pack({', '.join(sql)})", pa.struct(arrow)


def struct_lists(con, source: str, schema: pa.Schema) -> list[StructList]:
    """The list columns of the Parquet *source* (a quoted path) with one layout in every row.

    A column qualifies when its non-null rows all list the same names in the same
    order, none twice, and scores keep one direction per name. Each value is
    stored as ``FLOAT``, then ``DOUBLE``, when every one of its values reads back
    unchanged, and as it is otherwise; a ``cv_value`` string must read back as the
    same string.
    """
    found = []
    for column in STRUCT_LIST_COLUMNS:
        if column not in schema.names or not pa.types.is_list(schema.field(column).type):
            continue
        name = validate_identifier(column)
        layouts = con.sql(
            sql_build(
                "SELECT $layout AS layout FROM read_parquet($source) WHERE $column IS NOT NULL GROUP BY layout LIMIT 2",
                layout=sql_build(_LAYOUT_SQL[column], column=name),
                source=source,
                column=name,
            )
        ).fetchall()
        if len(layouts) != 1 or not _representable(column, layouts[0][0]):
            continue
        layout = layouts[0][0]
        leaves = _leaves(column, layout)
        stored = _stored_types(con, source, leaves)
        struct_sql, arrow_type = _struct([(path, f"v.{expression}", kind) for (path, expression, _), kind in zip(leaves, stored)])
        select = f"CASE WHEN v.{name} IS NULL THEN NULL ELSE {struct_sql} END AS {name}"
        found.append(StructList(column, layout, select, arrow_type))
    return found


def footer_value(found: list[StructList], previous: bytes | None) -> bytes:
    """The ``qpx.struct_lists`` footer value: the layouts in *previous* and those of *found*."""
    layouts = json.loads(previous) if previous else {}
    layouts.update({struct_list.column: struct_list.layout for struct_list in found})
    return json.dumps(layouts, separators=(",", ":")).encode()


def _literal(value: str) -> str:
    """A SQL string literal of a field name, which never holds a quote."""
    return f"'{value}'"


def _field(expression: str, name: str) -> str:
    return f"struct_extract({expression}, {_literal(name)})"


def restore_sql(column: str, layout: list[dict]) -> str:
    """The select item that reads the struct *column* back as its list, entries in *layout* order."""
    name = validate_identifier(column)
    if column == "additional_scores":
        entries = [
            f"struct_pack(score_name := {_literal(entry['name'])}, "
            f"score_value := CAST({_field(name, entry['name'])} AS DOUBLE), "
            f"higher_better := CAST({json.dumps(entry['higher_better'])} AS BOOLEAN))"
            for entry in layout
        ]
    elif column == "cv_params":
        entries = [
            f"struct_pack(cv_name := {_literal(entry['name'])}, cv_value := CAST({_field(name, entry['name'])} AS VARCHAR))"
            for entry in layout
        ]
    else:
        entries = []
        for entry in layout:
            label = _field(name, entry["name"])
            values = ", ".join(
                f"struct_pack(intensity_name := {_literal(inner)}, intensity_value := CAST({_field(label, inner)} AS FLOAT))"
                for inner in entry["names"]
            )
            entries.append(f"struct_pack(label := {_literal(entry['name'])}, intensities := [{values}])")
    return f"CASE WHEN {name} IS NULL THEN NULL ELSE [{', '.join(entries)}] END AS {name}"


def file_layouts(metadata: dict[bytes, bytes] | None) -> dict[str, list[dict]]:
    """The struct-list layouts a file's footer *metadata* declares."""
    return json.loads((metadata or {}).get(STRUCT_LISTS_FOOTER_KEY, b"{}"))


def restore_table(con, table: pa.Table, layouts: dict[str, list[dict]]) -> pa.Table:
    """*table*, read from a file whose footer declares *layouts*, with those struct lists read back as lists."""
    columns = [column for column in STRUCT_LIST_COLUMNS if column in layouts]
    if not columns:
        return table
    con.register("struct_list_rows", table)
    try:
        items = ", ".join(restore_sql(column, layouts[column]) for column in columns)
        restored = con.sql(sql_build("SELECT $items FROM struct_list_rows", items=items)).to_arrow_table()
    finally:
        con.unregister("struct_list_rows")
    for column in columns:
        table = table.set_column(table.schema.get_field_index(column), column, restored.column(column))
    return table


def restore_items(con, files: str, column_types: dict[str, str]) -> list[str]:
    """Select items reading the struct lists that the footer of *files* (SQL) declares back as lists."""
    if not any(column_types.get(column, "").startswith("STRUCT") for column in STRUCT_LIST_COLUMNS):
        return []
    rows = con.execute(
        sql_build(
            "SELECT decode(value) FROM parquet_kv_metadata($files) WHERE key = '$key' LIMIT 1",
            files=files,
            key=STRUCT_LISTS_FOOTER_KEY.decode(),
        )
    ).fetchall()
    layouts = json.loads(rows[0][0]) if rows else {}
    return [
        restore_sql(column, layouts[column])
        for column in STRUCT_LIST_COLUMNS
        if column in layouts and column_types.get(column, "").startswith("STRUCT")
    ]
