"""Rewrite a view with its rows sorted for selective reads."""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from qpx.core.engine import create_connection
from qpx.core.sql import validate_identifier
from qpx.transforms.protein_properties import _stamped_schema
from qpx.writers.base import _cast_column
from qpx.writers.encoding import ROW_GROUP_ROWS, SampledParquetWriter, row_group_rows

# Protein lookups read feature and pg, and run-by-run processing reads psm and mz.
DEFAULT_SORT_KEYS = {
    "feature": ("anchor_protein", "run_file_name"),
    "pg": ("anchor_protein", "label"),
    "psm": ("run_file_name", "scan"),
    "mz": ("run_file_name", "scan"),
}
SORT_ORDER_FOOTER_KEY = b"qpx.sort_order"


def relayout_view(
    source: Path,
    destination: Path,
    sort_by: tuple[str, ...],
    *,
    row_group_size: int = ROW_GROUP_ROWS,
    compression_level: int | None = None,
) -> None:
    """Write *source* to *destination* sorted by *sort_by*, recording the order in the footer.

    Rows stream from a DuckDB sort that spills to disk, so the view may be larger
    than memory. Values, ids and the schema are unchanged.
    """
    schema, compression = _stamped_schema(pq.read_schema(source))
    unsortable = [column for column in sort_by if column not in schema.names or pa.types.is_nested(schema.field(column).type)]
    if unsortable:
        raise ValueError(f"{source.name} has no flat column named {', '.join(unsortable)}")
    schema = schema.with_metadata({**schema.metadata, SORT_ORDER_FOOTER_KEY: ",".join(sort_by).encode()})
    sorting = pq.SortingColumn.from_ordering(schema, [(column, "ascending") for column in sort_by])
    order_by = ", ".join(validate_identifier(column) for column in sort_by)
    quoted_source = "'" + str(source).replace("'", "''") + "'"
    with (
        create_connection() as con,
        SampledParquetWriter(destination, schema, compression, compression_level, sorting_columns=sorting) as writer,
    ):
        rows = con.sql(f"SELECT * FROM read_parquet({quoted_source}) ORDER BY {order_by}")  # nosec B608
        for batch in rows.to_arrow_reader(row_group_size):
            columns = [_cast_column(pa.chunked_array([batch.column(field.name)]), field.type) for field in schema]
            table = pa.Table.from_arrays(columns, schema=schema)
            writer.write_table(table, row_group_size=min(row_group_size, row_group_rows(table)))
