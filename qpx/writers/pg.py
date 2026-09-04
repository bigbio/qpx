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

import pyarrow as pa

from qpx.core.data import PgSchema
from qpx.writers.base import BaseWriter


def _explode_pg_records(records: list[dict]) -> list[dict]:
    """Expand each protein-group record into one row per label.

    ``intensities`` (``list<{label, intensity}>``) is unrolled into scalar
    ``label`` + ``intensity`` columns; ``additional_intensities`` is filtered to
    the entries carrying that label (kept as ``list<additional_intensity>``).
    A record with no primary intensity (identification-only) yields a single row
    with null ``label``/``intensity``/``additional_intensities``.
    """
    exploded: list[dict] = []
    for record in records:
        intensities = record.get("intensities")
        base = {k: v for k, v in record.items() if k != "intensities"}
        additional = base.pop("additional_intensities", None) or []
        additional_by_label: dict[object, list] = {}
        for entry in additional:
            if isinstance(entry, dict):
                additional_by_label.setdefault(entry.get("label"), []).append(entry)

        if not intensities:
            # No intensities *list* to explode. The record may already be flat
            # (scalar label/intensity, e.g. the openms-consensus interim pg where
            # the label slot is populated but intensity is null) — pass those
            # through. A record with neither is identification-only (null label).
            row = dict(base)
            row.setdefault("label", None)
            row.setdefault("intensity", None)
            row["additional_intensities"] = additional or None
            exploded.append(row)
            continue

        for entry in intensities:
            label = entry.get("label") if isinstance(entry, dict) else None
            intensity = entry.get("intensity") if isinstance(entry, dict) else None
            row = dict(base)
            row["label"] = label
            row["intensity"] = intensity
            row["additional_intensities"] = additional_by_label.get(label) or None
            exploded.append(row)
    return exploded


class PgWriter(BaseWriter):
    _schema_class = PgSchema

    def _write_batch(self, records: list[dict]):
        """Explode ``intensities`` into one row per label, then buffer/flush."""
        super()._write_batch(_explode_pg_records(records))

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
            exploded = _explode_pg_records(table.to_pylist())
            # Build with a nullable id; the base implementation derives pg_id.
            relaxed = self._relaxed_id_schema()
            batch = pa.RecordBatch.from_pylist(exploded, schema=relaxed)
            super()._write_table(pa.Table.from_batches([batch], schema=relaxed))
            return
        super()._write_table(table)
