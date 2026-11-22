"""Restore name-value lists from structs using the layout declared in the Parquet footer."""

from __future__ import annotations

import json

import pyarrow as pa

from qpx.core.sql import sql_build, validate_identifier

STRUCT_LISTS_FOOTER_KEY = b"qpx.struct_lists"

STRUCT_LIST_COLUMNS = ("additional_scores", "cv_params", "additional_intensities")


def _literal(value: str) -> str:
    """A SQL string literal of a field name, which never holds a quote."""
    return f"'{value}'"


def _field(expression: str, name: str) -> str:
    return f"struct_extract({expression}, {_literal(name)})"


def restore_sql(column: str, layout: list[dict]) -> str:
    """The select item that reads the struct *column* back as its list, entries in *layout* order.

    A null score or cv field is a name the row does not hold, so it has no entry.
    """
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
    restored = f"[{', '.join(entries)}]"
    if column in ("additional_scores", "cv_params"):
        present = ", ".join(f"{_field(name, entry['name'])} IS NOT NULL" for entry in layout)
        restored = f"list_where({restored}, [{present}])"
    return f"CASE WHEN {name} IS NULL THEN NULL ELSE {restored} END AS {name}"


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


def is_struct(column_type: str) -> bool:
    """Whether the DuckDB *column_type* is a struct; a list or array of structs ends with ``]``."""
    return column_type.startswith("STRUCT") and not column_type.endswith("]")


def footer_layouts(con, files: str) -> dict[str, dict[str, list[dict]]]:
    """The struct-list layouts that the footers of the Parquet *files* (SQL) declare, by file name."""
    rows = con.execute(
        sql_build(
            "SELECT file_name, decode(value) FROM parquet_kv_metadata($files) WHERE key = '$key'",
            files=files,
            key=STRUCT_LISTS_FOOTER_KEY.decode(),
        )
    ).fetchall()
    return {name: json.loads(value) for name, value in rows}


def layout_groups(con, files: str) -> list[tuple[list[str], dict[str, list[dict]]]]:
    """The Parquet *files* (SQL), in file order, in groups that one ``read_parquet`` reads as each file alone.

    The files of a group declare the same struct-list layouts in their footers
    (``{}`` for none), store the fields of those structs as the same types, and
    hold the same struct-list columns as structs without a layout. Read together,
    a ``FLOAT`` field would widen to another file's ``DOUBLE``, which changes the
    string of a ``cv_value`` (``52.822754`` reads ``52.82275390625``), and a struct
    next to a list of one column fails once the column is read.
    """
    names = [row[0] for row in con.execute(sql_build("SELECT file FROM glob($files)", files=files)).fetchall()]
    declared = footer_layouts(con, files)
    if not declared.keys() <= set(names):
        raise ValueError(f"Cannot match the footers of {', '.join(sorted(declared.keys() - set(names)))} to the files read")
    columns = ", ".join(_literal(column) for column in STRUCT_LIST_COLUMNS)
    # The physical type of each struct-list leaf, by its path: ('cv_params, irt', 'FLOAT').
    leaves: dict[str, list[tuple[str, str]]] = {}
    if declared:
        for name, path, kind in con.execute(
            sql_build(
                "SELECT DISTINCT file_name, path_in_schema, type FROM parquet_metadata($files) "
                "WHERE split_part(path_in_schema, ', ', 1) IN ($columns) ORDER BY ALL",
                files=files,
                columns=columns,
            )
        ).fetchall():
            leaves.setdefault(name, []).append((path, kind))
    # The struct-list columns each file holds as a struct: a group node that is no list.
    structs: dict[str, set[str]] = {}
    for name, column in con.execute(
        sql_build(
            "SELECT file_name, name FROM parquet_schema($files) WHERE name IN ($columns) "
            "AND num_children > 0 AND converted_type IS NULL AND repetition_type <> 'REPEATED'",
            files=files,
            columns=columns,
        )
    ).fetchall():
        structs.setdefault(name, set()).add(column)
    groups: dict[str, tuple[list[str], dict]] = {}
    for name in names:
        layouts = declared.get(name, {})
        stored = [leaf for leaf in leaves.get(name, []) if leaf[0].split(", ")[0] in layouts]
        bare = sorted(structs.get(name, set()) - layouts.keys())
        groups.setdefault(json.dumps([layouts, stored, bare], sort_keys=True), ([], layouts))[0].append(name)
    return list(groups.values())


def restore_items(layouts: dict[str, list[dict]], column_types: dict[str, str]) -> list[str]:
    """Select items reading the struct lists of *layouts* that *column_types* (DuckDB) holds as structs back as lists."""
    return [
        restore_sql(column, layouts[column])
        for column in STRUCT_LIST_COLUMNS
        if column in layouts and is_struct(column_types.get(column, ""))
    ]
