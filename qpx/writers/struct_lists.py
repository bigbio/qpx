"""Pack logical name-value lists into the writer's fixed Parquet structs."""

import json

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc

from qpx.core.struct_lists import STRUCT_LISTS_FOOTER_KEY

_ENTRY_FIELDS = {
    "additional_scores": ("score_name", "score_value", pa.float64()),
    "cv_params": ("cv_name", "cv_value", pa.string()),
}


def struct_schema(schema: pa.Schema, layouts: dict[str, list[dict]]) -> pa.Schema:
    """Declare the nonempty fixed layouts and give their fields the fixed value types."""
    physical = {column: layout for column, layout in layouts.items() if layout}
    if not physical:
        return schema
    fields = [
        field.with_type(pa.struct([pa.field(entry["name"], _ENTRY_FIELDS[field.name][2]) for entry in physical[field.name]]))
        if field.name in physical
        else field
        for field in schema
    ]
    metadata = {**schema.metadata, STRUCT_LISTS_FOOTER_KEY: json.dumps(physical, separators=(",", ":")).encode()}
    return pa.schema(fields, metadata=metadata)


def _pack_column(column: pa.ListArray, layout: list[dict], field: pa.Field, source: str) -> pa.Array:
    """Scatter each named entry's value to its parent row in the corresponding struct field."""
    name_field, value_field, _ = _ENTRY_FIELDS[field.name]
    entries = column.flatten()
    names = entries.field(name_field)
    layout_names = pa.array([entry["name"] for entry in layout], type=pa.string())
    codes = pc.call_function("index_in", [names], options=pc.SetLookupOptions(layout_names))
    if codes.null_count:
        unknown = names.filter(codes.is_null())[0].as_py()
        raise ValueError(f"{source} {field.name}: name {unknown!r} is outside the fixed layout")
    values = entries.field(value_field)
    if values.null_count:
        # A null field means no entry, so a listed name must carry a value.
        raise ValueError(f"{source} {field.name}: name {names.filter(values.is_null())[0].as_py()!r} has a null value")
    if not layout:
        return column

    # Use the visible list lengths: flatten skips children hidden by a null list.
    lengths = pc.fill_null(pc.call_function("list_value_length", [column]), 0).to_numpy()
    parents = np.repeat(np.arange(len(column)), lengths)
    codes = codes.to_numpy()
    children = []
    for code, child in enumerate(field.type):
        positions = np.flatnonzero(codes == code)
        rows = parents[positions]
        if (np.diff(rows) == 0).any():
            raise ValueError(f"{source} {field.name}: duplicate name {child.name!r} in one row")
        indices = np.full(len(column), -1, dtype=np.int64)
        indices[rows] = positions
        children.append(pc.take(values, pa.array(indices, mask=indices < 0)))
    return pa.StructArray.from_arrays(children, fields=list(field.type), mask=column.is_null())


def pack_struct_lists(table: pa.Table, schema: pa.Schema, layouts: dict[str, list[dict]], source: str) -> pa.Table:
    """Apply the selected layouts while keeping every other column as supplied."""
    for column, layout in layouts.items():
        index = table.schema.get_field_index(column)
        field = schema.field(column)
        packed = _pack_column(table.column(index).combine_chunks(), layout, field, source)
        table = table.set_column(index, field, packed)
    return table
