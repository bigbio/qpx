"""Parquet encoding options shared by the QPX writers."""

from __future__ import annotations

import pyarrow as pa
import pyarrow.parquet as pq

# ZSTD level 9 (vs the pyarrow default of 3) buys a few extra percent at modest
# write-cost; only codecs that accept a level get one.
_COMPRESSION_LEVELS: dict[str, int] = {"zstd": 9, "gzip": 9}

# Row groups stay small enough for selective reads: 128k rows with a page index made
# protein lookups 3-4x faster and 5-7x cheaper remotely than 1 Mi-row groups (783k rows).
# Spectra rows are large, so a byte cap bounds what one row group makes a reader decode.
ROW_GROUP_ROWS = 128 * 1024
ROW_GROUP_BYTES = 256 * 1024 * 1024

# The smallest encoding depends on the data, not only on the column: TMT reporter
# intensities of MSV000085836 took 253 MB as a dictionary, 225 MB plain and 200 MB
# byte-stream-split, and PEP is smallest as a dictionary in DIA-NN files but
# byte-stream-split in OpenMS ones. Each leaf therefore takes the encoding its first
# rows measure smallest. A dictionary is rebuilt in every row group, so it loses on
# high-entropy values once row groups are small.
_DICTIONARY = "dictionary"
_PLAIN = "PLAIN"
_FLOAT_ENCODINGS = (_DICTIONARY, _PLAIN, "BYTE_STREAM_SPLIT")
_INTEGER_ENCODINGS = (_DICTIONARY, _PLAIN, "DELTA_BINARY_PACKED")
_OTHER_ENCODINGS = (_DICTIONARY, _PLAIN)
# Strings keep the dictionary readers group and filter on unless PLAIN is more than 20% smaller.
_STRING_PLAIN_RATIO = 0.8


def _leaf_types(schema: pa.Schema) -> dict[str, pa.DataType]:
    """Every leaf column path, in pyarrow's dotted convention (``a.list.element``), with its type."""

    def _walk(dtype: pa.DataType, path: str) -> dict[str, pa.DataType]:
        if pa.types.is_struct(dtype):
            leaves: dict[str, pa.DataType] = {}
            for index in range(dtype.num_fields):
                child = dtype.field(index)
                leaves.update(_walk(child.type, f"{path}.{child.name}"))
            return leaves
        if pa.types.is_map(dtype):
            return {**_walk(dtype.key_type, f"{path}.key_value.key"), **_walk(dtype.item_type, f"{path}.key_value.value")}
        if pa.types.is_list(dtype) or pa.types.is_large_list(dtype):
            return _walk(dtype.value_type, f"{path}.list.element")
        return {path: dtype}

    leaves: dict[str, pa.DataType] = {}
    for field in schema:
        leaves.update(_walk(field.type, field.name))
    return leaves


def _candidates(dtype: pa.DataType) -> tuple[str, ...]:
    """The encodings tried for a leaf of type *dtype*."""
    if pa.types.is_float32(dtype) or pa.types.is_float64(dtype):
        return _FLOAT_ENCODINGS
    if pa.types.is_int32(dtype) or pa.types.is_int64(dtype):
        return _INTEGER_ENCODINGS
    return _OTHER_ENCODINGS


def _leaf_sizes(sample: pa.Table, encoding: str, leaves: list[str], options: dict) -> dict[str, int]:
    """Compressed bytes of each of *leaves* when *sample* is written as one row group with *encoding* on them."""
    trial = {**options, "write_page_index": False, "use_dictionary": encoding == _DICTIONARY}
    if encoding not in (_DICTIONARY, _PLAIN):
        trial["column_encoding"] = dict.fromkeys(leaves, encoding)
        # Only the top-level columns holding such a leaf need writing.
        sample = sample.select(sorted({leaf.split(".", 1)[0] for leaf in leaves}, key=sample.schema.get_field_index))
    sink = pa.BufferOutputStream()
    pq.write_table(sample, sink, row_group_size=sample.num_rows, **trial)
    row_group = pq.read_metadata(pa.BufferReader(sink.getvalue())).row_group(0)
    written = {
        row_group.column(i).path_in_schema: row_group.column(i).total_compressed_size for i in range(row_group.num_columns)
    }
    return {leaf: written[leaf] for leaf in leaves}


def _smallest_encodings(schema: pa.Schema, sample: pa.Table, options: dict) -> dict[str, str]:
    """The encoding of each leaf of *schema* that measures smallest on *sample*."""
    types = _leaf_types(schema)
    tried: dict[str, list[str]] = {}
    for leaf, dtype in types.items():
        for encoding in _candidates(dtype):
            tried.setdefault(encoding, []).append(leaf)
    sizes = {encoding: _leaf_sizes(sample, encoding, leaves, options) for encoding, leaves in tried.items()}
    picked: dict[str, str] = {}
    for leaf, dtype in types.items():
        measured = {encoding: sizes[encoding][leaf] for encoding in _candidates(dtype)}
        best = min(measured, key=measured.__getitem__)
        if best == _PLAIN and pa.types.is_string(dtype) and measured[_PLAIN] > _STRING_PLAIN_RATIO * measured[_DICTIONARY]:
            best = _DICTIONARY
        picked[leaf] = best
    return picked


def row_group_rows(data: pa.Table | pa.RecordBatch, rows: int | None = None) -> int:
    """Rows per row group for *data*: *rows*, by default ``ROW_GROUP_ROWS``, fewer when rows are large."""
    bytes_per_row = max(1, data.nbytes // max(1, data.num_rows))
    return max(1, min(rows or ROW_GROUP_ROWS, ROW_GROUP_BYTES // bytes_per_row))


def parquet_write_options(
    schema: pa.Schema, compression: str, compression_level: int | None = None, sample: pa.Table | None = None
) -> dict:
    """Build QPX ParquetWriter options for *schema* and *compression*; the level overrides the codec's default.

    *sample* holds the first rows to be written, with *schema*: each leaf takes the
    encoding that measures smallest on its first row group. Without rows, every leaf
    keeps pyarrow's default, a dictionary.
    """
    options: dict = {"compression": compression, "version": "2.6", "write_page_index": True}
    level = _COMPRESSION_LEVELS.get(str(compression).lower()) if compression_level is None else compression_level
    if level is not None:
        options["compression_level"] = level
    if sample is None or not sample.num_rows:
        return options

    picked = _smallest_encodings(schema, sample.slice(0, row_group_rows(sample)), options)
    options["use_dictionary"] = [leaf for leaf, encoding in picked.items() if encoding == _DICTIONARY]
    explicit = {leaf: encoding for leaf, encoding in picked.items() if encoding not in (_DICTIONARY, _PLAIN)}
    if explicit:
        options["column_encoding"] = explicit
    return options


class SampledParquetWriter:
    """A ``pq.ParquetWriter`` opened at the first table written, encoded as :func:`parquet_write_options` picks for it."""

    def __init__(self, where, schema: pa.Schema, compression: str, compression_level: int | None = None, **writer_options):
        self._where = str(where)
        self._schema = schema
        self._compression = compression
        self._compression_level = compression_level
        self._writer_options = writer_options
        self._writer: pq.ParquetWriter | None = None

    def _open(self, sample: pa.Table | None) -> pq.ParquetWriter:
        options = parquet_write_options(self._schema, self._compression, self._compression_level, sample)
        return pq.ParquetWriter(self._where, self._schema, **self._writer_options, **options)

    def write_table(self, table: pa.Table, row_group_size: int | None = None) -> None:
        """Write *table*, opening the file with the encodings *table* measures smallest when it is the first."""
        if self._writer is None:
            self._writer = self._open(table)
        self._writer.write_table(table, row_group_size=row_group_size)

    def close(self) -> None:
        """Close the file, writing an empty one with the schema when no table was written."""
        if self._writer is None:
            self._writer = self._open(None)
        self._writer.close()

    def __enter__(self) -> SampledParquetWriter:
        return self

    def __exit__(self, *exc_info) -> None:
        self.close()
