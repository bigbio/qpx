"""PgWriter — writes pg.parquet files.

Since QPX 1.1 the pg view is **flattened**: instead of one row per
``(pg_accessions, grouped_runs)`` carrying an ``intensities: list<{label,
intensity}>`` column, there is **one row per label** with scalar ``label`` and
``intensity`` columns. Since the mandatory-identity change the primary key is the
single derived ``pg_id`` (hashed from the ``identity_composite``
``(pg_accessions, grouped_runs, label)`` — the full group membership keys the id,
not just the leader). Converters
still build the natural ``intensities`` list per protein group; the writer is the
single place that materializes the flat physical layout, so no converter needs to
know about the on-disk shape. Identification-only groups (no intensity, e.g.
mzIdentML) become a single row with ``label``/``intensity`` null.
"""

from __future__ import annotations

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc

from qpx.core.data import PgSchema
from qpx.writers.base import BaseWriter


def _field(values: pa.Array, name: str) -> pa.Array:
    """One field of a struct array; null where the struct itself is null."""
    return pc.call_function("struct_field", [values], pc.StructFieldOptions(name))


def _overlay(base: pa.Array, values: pa.Array, rows: np.ndarray, picks: np.ndarray) -> pa.Array:
    """*base* with ``values[picks]`` in place of ``base[rows]``."""
    take = np.arange(len(base))
    take[rows] = len(base) + picks
    return pa.concat_arrays([base, values]).take(take)


def _additional_for_labels(additional: pa.ListArray, owners: np.ndarray, labels: pa.Array) -> pa.ListArray:
    """For each (record, label) pair, the record's additional intensities carrying that label.

    *owners* gives each pair's index in *additional*. Null entries are dropped, a
    null label matches a null label, entries keep their order, and a pair without
    any entry gets null.
    """
    values = pc.call_function("list_flatten", [additional])
    parents = pc.call_function("list_parent_indices", [additional]).to_numpy()
    valid = values.is_valid().to_numpy(zero_copy_only=False)
    values, parents = values.filter(valid), parents[valid]
    options = pc.DictionaryEncodeOptions(null_encoding="encode")
    encoded = pc.call_function("dictionary_encode", [pa.concat_arrays([_field(values, "label"), labels])], options)
    codes = encoded.indices.to_numpy().astype(np.int64)
    width = len(encoded.dictionary)
    value_keys = parents * width + codes[: len(values)]
    pair_keys = owners * width + codes[len(values) :]
    order = np.argsort(value_keys, kind="stable")
    first = np.searchsorted(value_keys[order], pair_keys, "left")
    counts = np.searchsorted(value_keys[order], pair_keys, "right") - first
    offsets = np.concatenate([[0], np.cumsum(counts)])
    positions = np.repeat(first - offsets[:-1], counts) + np.arange(offsets[-1])
    return pa.ListArray.from_arrays(
        pa.array(offsets, pa.int32()), values.take(order[positions]), type=additional.type, mask=pa.array(counts == 0)
    )


def _explode_pg_table(records: list[dict], schema: pa.Schema) -> pa.Table:
    """Expand each protein-group record into one row per label, as a table on *schema*.

    ``intensities`` (``list<{label, intensity}>``) is unrolled into scalar
    ``label`` + ``intensity`` columns; ``additional_intensities`` is filtered to
    the entries carrying that label (kept as ``list<additional_intensity>``).
    Non-null ``cv_params`` on an intensity entry override the record-level
    parameters for that label.
    A record with no primary intensity (identification-only) yields a single row
    with null ``label``/``intensity`` and its whole ``additional_intensities``;
    a record that is already flat keeps its own ``label``/``intensity``.
    The records are converted to Arrow once and the rows are taken from them, so
    a group is never copied into one Python dict per label.
    """
    entry = pa.struct([schema.field("label"), schema.field("intensity"), schema.field("cv_params")])
    table = pa.Table.from_pylist(records, schema=schema.append(pa.field("intensities", pa.list_(entry))))
    lists = table.column("intensities").combine_chunks()
    counts = pc.fill_null(pc.list_value_length(lists), 0).to_numpy()
    repeats = np.maximum(counts, 1)
    parents = np.repeat(np.arange(len(lists)), repeats)
    rows = np.flatnonzero(np.repeat(counts > 0, repeats))
    entries = pc.call_function("list_flatten", [lists])
    picks = np.arange(len(entries))
    out = table.select(schema.names).take(parents)

    labels = _field(entries, "label")
    cv_params = _field(entries, "cv_params")
    with_cv = cv_params.is_valid().to_numpy(zero_copy_only=False)
    additional = out.column("additional_intensities").combine_chunks()
    empty = np.flatnonzero(pc.fill_null(pc.equal(pc.list_value_length(additional), 0), True).to_numpy(zero_copy_only=False))
    additional = _overlay(additional, pa.nulls(1, additional.type), empty, np.zeros(len(empty), np.int64))
    by_label = _additional_for_labels(table.column("additional_intensities").combine_chunks(), parents[rows], labels)
    replaced = {
        "label": _overlay(out.column("label").combine_chunks(), labels, rows, picks),
        "intensity": _overlay(out.column("intensity").combine_chunks(), _field(entries, "intensity"), rows, picks),
        "cv_params": _overlay(out.column("cv_params").combine_chunks(), cv_params, rows[with_cv], picks[with_cv]),
        "additional_intensities": _overlay(additional, by_label, rows, picks),
    }
    for name, values in replaced.items():
        index = out.schema.get_field_index(name)
        out = out.set_column(index, out.schema.field(index), values)
    return out


class PgWriter(BaseWriter):
    _schema_class = PgSchema
    _record_only_keys = ("intensities",)

    def _write_batch(self, records: list[dict]):
        """Explode ``intensities`` into one row per label, then buffer/flush.

        The buffer holds exploded Arrow tables instead of records; full batches of
        ``batch_size`` rows are written as for any other view.
        """
        self._buffer.append(_explode_pg_table(records, self._relaxed_id_schema()))
        rows = sum(table.num_rows for table in self._buffer)
        if rows < self._batch_size:
            return
        table = pa.concat_tables(self._buffer)
        full = rows - rows % self._batch_size
        for start in range(0, full, self._batch_size):
            self._write_arrow_batch([table.slice(start, self._batch_size)])
        self._buffer = [table.slice(full)] if full < rows else []

    def _write_arrow_batch(self, records: list[pa.Table]):
        """Write buffered exploded tables."""
        self._write_identity_table(pa.concat_tables(records))

    def _write_dataframe(self, df):
        """Explode ``intensities`` before writing a flat pg DataFrame.

        The base ``write_dataframe`` builds the table against the flat pg schema
        (which has no ``intensities`` column), so a DataFrame that still carries
        the natural ``intensities`` list would have its quant silently dropped
        (bigbio/qpx#252). Route such frames through the same explode the
        ``write_batch`` / ``write_table`` paths use; frames already flattened
        (scalar ``label``/``intensity``) fall through unchanged.
        """
        if "intensities" in df.columns:
            self._write_batch(df.to_dict("records"))
            return
        super()._write_dataframe(df)

    def _write_table(self, table: pa.Table):
        """Explode a pg table that still carries the ``intensities`` list column."""
        if "intensities" in table.schema.names:
            # Built with a nullable id; the base implementation derives pg_id.
            super()._write_table(_explode_pg_table(table.to_pylist(), self._relaxed_id_schema()))
            return
        super()._write_table(table)
