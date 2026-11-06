"""Rewrite a dataset's views sorted for selective reads."""

from __future__ import annotations

from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from qpx.core.engine import create_connection
from qpx.core.sql import sql_build, validate_identifier
from qpx.core.struct_lists import STRUCT_LISTS_FOOTER_KEY, footer_value, struct_lists
from qpx.transforms.mantissa import MANTISSA_BITS, MANTISSA_BITS_FOOTER_KEY, mantissa_footer, round_table, rounding_plan
from qpx.transforms.protein_properties import _stamped_schema
from qpx.writers.base import _cast_column
from qpx.writers.encoding import ROW_GROUP_ROWS, SampledParquetWriter, row_group_rows

# Protein and peptide lookups read feature, pg and psm, and run-by-run processing reads mz.
# Rows of one peptide sit together, which also makes feature and psm about a third
# smaller on TMT data; pg keeps a group's runs together for the same reason.
DEFAULT_SORT_KEYS = {
    "feature": ("anchor_protein", "peptidoform", "charge", "run_file_name"),
    "pg": ("anchor_protein", "grouped_runs", "label"),
    "psm": ("peptidoform", "charge", "run_file_name", "scan"),
    "mz": ("run_file_name", "scan"),
}
SORT_ORDER_FOOTER_KEY = b"qpx.sort_order"

# A sorted view's ids become 1, 2, 3, ... in its new row order. A hashed id is eight
# random bytes per row that no encoding compresses; ordinals delta-encode to almost
# nothing (MSV000085836 sorted pg 55.8 to 17.0 MB, psm 245.4 to 150.7 MB). Views are
# sorted in this order, so a view's references to an earlier one are rewritten as it
# is sorted.
ORDINAL_IDS = {"pg": "pg_id", "feature": "feature_id", "psm": "psm_id"}

# The columns of each view that hold the ids of another view. A dotted name is
# a field of the structs in a list column.
ID_REFERENCES = {
    "feature": {"pg_ids": "pg", "psm_ids": "psm"},
    "psm": {"feature_id": "feature"},
    "fragment": {"feature_id": "feature"},
    "site": {"pg_id": "pg", "best_psm_id": "psm"},
    "proteoform": {"pg_id": "pg", "best_psm_id": "psm"},
    "xic": {"feature_ids": "feature", "feature_scores.feature_id": "feature"},
}


def _sortable(field: pa.Field) -> bool:
    """Whether DuckDB can order rows by *field*: a flat column, or a list of flat values."""
    if pa.types.is_list(field.type) or pa.types.is_large_list(field.type):
        return not pa.types.is_nested(field.type.value_type)
    return not pa.types.is_nested(field.type)


def _quoted(path: Path) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def _holds_values(con, source: Path, column: str, is_list: bool) -> bool:
    """Whether some row of *source* has a value in *column* (a non-empty list, for a list column)."""
    template = "len($column) > 0" if is_list else "$column IS NOT NULL"
    sql = sql_build(
        "SELECT EXISTS (SELECT 1 FROM read_parquet($source) WHERE $condition)",
        source=_quoted(source),
        condition=sql_build(template, column=validate_identifier(column)),
    )
    return bool(con.sql(sql).fetchall()[0][0])


_LIST_REMAP_JOIN = (
    "LEFT JOIN (SELECT e.row_number, list($new ORDER BY e.position) AS ids FROM ("
    "SELECT file_row_number AS row_number, unnest($column) AS old, unnest(range(1, len($column) + 1)) AS position "
    "FROM read_parquet($source, file_row_number=true)"
    ") AS e LEFT JOIN read_parquet($ids) AS m ON $old = m.old GROUP BY e.row_number) AS $alias "
    "ON v.file_row_number = $alias.row_number"
)


def _list_remap_join(column: str, field: str, source: Path, id_map: Path, alias: str) -> str:
    """The join giving each row of *source* its list *column* with new ids, in its elements or their *field*."""
    if field:
        old = f"e.old.{validate_identifier(field)}"
        new = f"struct_update(e.old, {validate_identifier(field)} := coalesce(m.new, {old}))"
    else:
        old, new = "e.old", "coalesce(m.new, e.old)"
    return sql_build(
        _LIST_REMAP_JOIN,
        new=new,
        old=old,
        column=validate_identifier(column),
        source=_quoted(source),
        ids=_quoted(id_map),
        alias=alias,
    )


def _select_sql(source: Path, schema: pa.Schema, remaps: dict[str, Path], order_by: str, items: list[str]) -> str:
    """The rows of *source* in *order_by*, each reference of *remaps* holding the new ids its map gives.

    A reference that its map does not name keeps its value. *items* replace
    further columns.
    """
    joins, replace = [], list(items)
    for index, (reference, id_map) in enumerate(sorted(remaps.items())):
        column, _, field = reference.partition(".")
        name, alias = validate_identifier(column), f"m{index}"
        if pa.types.is_list(schema.field(column).type):
            joins.append(_list_remap_join(column, field, source, id_map, alias))
            replace.append(
                sql_build("CASE WHEN len(v.$column) > 0 THEN $alias.ids ELSE v.$column END AS $column", column=name, alias=alias)
            )
        else:
            joins.append(
                sql_build(
                    "LEFT JOIN read_parquet($ids) AS $alias ON v.$column = $alias.old",
                    ids=_quoted(id_map),
                    alias=alias,
                    column=name,
                )
            )
            replace.append(sql_build("coalesce($alias.new, v.$column) AS $column", alias=alias, column=name))
    columns = "v.* EXCLUDE (file_row_number)" + (f" REPLACE ({', '.join(replace)})" if replace else "")
    return sql_build(
        "SELECT $columns FROM read_parquet($source, file_row_number=true) AS v $joins ORDER BY $order_by",
        columns=columns,
        source=_quoted(source),
        joins=" ".join(joins),
        order_by=order_by,
    )


def _footer_schema(source: Path) -> tuple[pa.Schema, str]:
    """The schema of *source* with every key of its footer, re-stamped, and its compression.

    Keys a writer adds at close, such as ``identity_source``, are in the footer only.
    """
    schema, compression = _stamped_schema(pq.read_schema(source))
    footer = {key: value for key, value in (pq.read_metadata(source).metadata or {}).items() if key != b"ARROW:schema"}
    return schema.with_metadata({**footer, **schema.metadata}), compression


@dataclass(frozen=True)
class _Rewrite:
    """How a view is written again: its row order, footer keys, ids and rounding.

    *sorting* lists the flat keys of *order_by* that ``sorting_columns`` can state.
    *remaps* gives the id maps of the reference columns, and *renumber* names the
    primary key to number 1, 2, 3, ... in the new order and the file that receives
    each old id with its new one. *mantissa_bits* gives the bits measured values keep.
    """

    order_by: str
    sorting: list[tuple[str, str]]
    footer: dict[bytes, bytes]
    remaps: dict[str, Path]
    renumber: tuple[str, Path] | None
    mantissa_bits: dict[str, int]


def _sorted(source: Path, sort_by: tuple[str, ...], *, remaps, renumber, mantissa_bits) -> _Rewrite:
    """*source* sorted by *sort_by*, ties broken by the old id when its ids are renumbered.

    Parquet's ``sorting_columns`` cannot name a list column, so it states the keys
    before the first one; ``qpx.sort_order`` states them all.
    """
    schema = pq.read_schema(source)
    unsortable = [column for column in sort_by if column not in schema.names or not _sortable(schema.field(column))]
    if unsortable:
        raise ValueError(f"{source.name} has no sortable column named {', '.join(unsortable)}")
    flat = []
    for column in sort_by:
        if pa.types.is_nested(schema.field(column).type):
            break
        flat.append((column, "ascending"))
    keys = [*sort_by, renumber[0]] if renumber else list(sort_by)
    footer = {SORT_ORDER_FOOTER_KEY: ",".join(sort_by).encode()}
    if renumber:
        footer[b"identity_source"] = b"provided"
    order_by = ", ".join(f"v.{validate_identifier(column)}" for column in keys)
    return _Rewrite(order_by, flat, footer, remaps, renumber, mantissa_bits)


def _kept_order(source: Path, remaps: dict[str, Path]) -> _Rewrite:
    """*source* in its row order, with the reference columns of *remaps* rewritten."""
    metadata = pq.read_metadata(source)
    sorting = metadata.row_group(0).sorting_columns if metadata.num_row_groups else ()
    keys = list(pq.SortingColumn.to_ordering(metadata.schema.to_arrow_schema(), sorting)[0]) if sorting else []
    return _Rewrite("v.file_row_number", keys, {}, remaps, renumber=None, mantissa_bits={})


def _written_schema(con, source: Path, rewrite: _Rewrite) -> tuple[pa.Schema, str, list[str], dict[str, int]]:
    """The schema and compression *source* is written with, the select items of its struct lists and its rounding.

    Name-value lists with one layout in every row become structs
    (:mod:`qpx.core.struct_lists`), and measured values are rounded
    (:mod:`qpx.transforms.mantissa`).
    """
    schema, compression = _footer_schema(source)
    footer = dict(rewrite.footer)
    found = struct_lists(con, _quoted(source), schema)
    for struct_list in found:
        schema = schema.set(schema.get_field_index(struct_list.column), pa.field(struct_list.column, struct_list.arrow_type))
    if found:
        footer[STRUCT_LISTS_FOOTER_KEY] = footer_value(found, schema.metadata.get(STRUCT_LISTS_FOOTER_KEY))
    plan = rounding_plan(rewrite.mantissa_bits, schema)
    if plan:
        footer[MANTISSA_BITS_FOOTER_KEY] = mantissa_footer(plan)
    return schema.with_metadata({**schema.metadata, **footer}), compression, [item.select for item in found], plan


def _write_view(source: Path, destination: Path, rewrite: _Rewrite, *, row_group_size: int) -> None:
    """Write the rows of *source* to *destination* as *rewrite* says.

    Rows stream from a DuckDB sort that spills to disk, so the view may be larger
    than memory.
    """
    with create_connection() as con:
        schema, compression, items, plan = _written_schema(con, source, rewrite)
    sql = _select_sql(source, schema, rewrite.remaps, rewrite.order_by, items)
    # Struct lists change the leaf columns that sorting_columns counts.
    sorting_columns = pq.SortingColumn.from_ordering(schema, rewrite.sorting)
    renumber = rewrite.renumber
    id_index = schema.get_field_index(renumber[0]) if renumber else -1
    id_schema = pa.schema([("old", pa.int64()), ("new", pa.int64())])
    with (
        create_connection() as con,
        SampledParquetWriter(destination, schema, compression, sorting_columns=sorting_columns) as writer,
        pq.ParquetWriter(str(renumber[1]), id_schema) if renumber else nullcontext() as ids,
    ):
        next_id = 1
        for batch in con.sql(sql).to_arrow_reader(row_group_size):
            columns = [_cast_column(pa.chunked_array([batch.column(field.name)]), field.type) for field in schema]
            if renumber:
                new = pa.chunked_array([np.arange(next_id, next_id + batch.num_rows, dtype=np.int64)])
                ids.write_table(pa.Table.from_arrays([columns[id_index], new], schema=id_schema))
                columns[id_index] = new
                next_id += batch.num_rows
            table = round_table(pa.Table.from_arrays(columns, schema=schema), plan)
            writer.write_table(table, row_group_size=row_group_rows(table, row_group_size))


def relayout_dataset(
    dataset: Path,
    prefix: str,
    staging: Path,
    sort_keys: dict[str, tuple[str, ...]],
    *,
    row_group_size: int = ROW_GROUP_ROWS,
) -> dict[str, str]:
    """Write the views of *sort_keys* to *staging* sorted, with ordinal ids, and the references to those ids.

    The ids of a sorted feature, pg or psm view are renumbered in its new row
    order, and every column that holds them, in a sorted view or not, is
    rewritten. The measured values of a sorted view keep the mantissa bits of
    :data:`~qpx.transforms.mantissa.MANTISSA_BITS`. Returns what was done to each
    staged ``{prefix}.{view}.parquet``.
    """
    done: dict[str, str] = {}
    with TemporaryDirectory(prefix=".relayout_ids-", dir=staging.parent) as temporary, create_connection() as con:
        id_maps: dict[str, Path] = {}
        remapped: dict[str, set[str]] = {}

        def pending(view: str, source: Path) -> dict[str, Path]:
            schema = pq.read_schema(source)
            found = {}
            for reference, target in ID_REFERENCES.get(view, {}).items():
                column = reference.partition(".")[0]
                if target not in id_maps or column not in schema.names or reference in remapped.get(view, set()):
                    continue
                if _holds_values(con, source, column, pa.types.is_list(schema.field(column).type)):
                    found[reference] = id_maps[target]
            return found

        order = [view for view in ORDINAL_IDS if view in sort_keys] + [view for view in sort_keys if view not in ORDINAL_IDS]
        for view in order:
            name = f"{prefix}.{view}.parquet"
            source = dataset / name
            if not source.is_file():
                raise ValueError(f"{dataset} has no {name}")
            remaps = pending(view, source)
            renumber = (ORDINAL_IDS[view], Path(temporary) / f"{view}.parquet") if view in ORDINAL_IDS else None
            rewrite = _sorted(
                source, sort_keys[view], remaps=remaps, renumber=renumber, mantissa_bits=MANTISSA_BITS.get(view, {})
            )
            _write_view(source, staging / name, rewrite, row_group_size=row_group_size)
            remapped[view] = set(remaps)
            if renumber:
                id_maps[view] = renumber[1]
            done[name] = f"sorted by {', '.join(sort_keys[view])}" + ("; ids renumbered" if renumber else "")

        # References to a view sorted later, and in views that were not sorted.
        for view in ID_REFERENCES:
            name = f"{prefix}.{view}.parquet"
            source = staging / name if name in done else dataset / name
            if not source.is_file() or not (remaps := pending(view, source)):
                continue
            rewrite = _kept_order(source, remaps)
            _write_view(source, Path(temporary) / name, rewrite, row_group_size=row_group_size)
            (Path(temporary) / name).replace(staging / name)
            done[name] = done.get(name, "kept its order") + f"; {', '.join(sorted(remaps))} renumbered"
    return done
