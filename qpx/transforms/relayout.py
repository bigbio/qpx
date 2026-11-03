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


def _sortable(field: pa.Field) -> bool:
    """Whether DuckDB can order rows by *field*: a flat column, or a list of flat values."""
    if pa.types.is_list(field.type) or pa.types.is_large_list(field.type):
        return not pa.types.is_nested(field.type.value_type)
    return not pa.types.is_nested(field.type)


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
    than memory. Values, ids and the schema are unchanged. Parquet's
    ``sorting_columns`` cannot name a list column, so it records the keys before
    the first one; ``qpx.sort_order`` records them all.
    """
    schema, compression = _stamped_schema(pq.read_schema(source))
    unsortable = [column for column in sort_by if column not in schema.names or not _sortable(schema.field(column))]
    if unsortable:
        raise ValueError(f"{source.name} has no sortable column named {', '.join(unsortable)}")
    schema = schema.with_metadata({**schema.metadata, SORT_ORDER_FOOTER_KEY: ",".join(sort_by).encode()})
    flat = []
    for column in sort_by:
        if pa.types.is_nested(schema.field(column).type):
            break
        flat.append((column, "ascending"))
    sorting = pq.SortingColumn.from_ordering(schema, flat)
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
