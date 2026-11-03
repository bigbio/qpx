"""BaseConverter -- abstract base for all tool-specific converters.

Each converter creates its own DuckDB in-memory connection, loads
tool-specific output, transforms it via SQL into the QPX schema, and pipes
rows into the appropriate Writer.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections.abc import Callable, Iterable

import duckdb
import pandas as pd
import pyarrow as pa

from qpx.converters.ptm import compute_precursor_mz
from qpx.core.cleavage import count_missed_cleavages
from qpx.core.engine import create_converter_connection
from qpx.core.sql import escape_path, sql_build, validate_identifier, validate_table
from qpx.writers.pepmap import PepMapWriter


def resolve_columns(
    field_mappings: dict[str, list[str]],
    available_columns: set[str],
) -> dict[str, str]:
    """Resolve QPX fields to actual tool column names.

    For each QPX field, tries candidates in order and returns the first match.
    Skips fields where no candidate matches the available columns.

    Args:
        field_mappings: QPX field name -> ordered list of candidate column names.
        available_columns: Set of column names present in the input data.

    Returns:
        Dict mapping QPX field name -> resolved tool column name.
    """
    resolved = {}
    for qpx_field, candidates in field_mappings.items():
        for candidate in candidates:
            if candidate in available_columns:
                resolved[qpx_field] = candidate
                break
    return resolved


class BaseConverter(ABC):
    """Abstract base class for converter adapters.

    Subclasses implement ``convert()`` which performs:
        1. Load tool output into DuckDB tables
        2. SQL transform into QPX schema columns
        3. Stream results into a Writer (FeatureWriter, PsmWriter, etc.)
    """

    def __init__(
        self,
        duckdb_memory: str | None = None,
        duckdb_threads: int | None = None,
        conn: duckdb.DuckDBPyConnection | None = None,
    ):
        self.logger = logging.getLogger(f"{self.__class__.__module__}.{self.__class__.__name__}")
        if conn is not None:
            self._conn = conn
            self._owns_conn = False
        else:
            self._conn = create_converter_connection(
                memory_limit=duckdb_memory,
                threads=duckdb_threads,
            )
            self._owns_conn = True
        # Accumulate discovered score names during conversion
        self._discovered_scores: set[str] = set()

    # ------------------------------------------------------------------
    # Public adapter API (used by orchestrators)
    # ------------------------------------------------------------------

    def get_discovered_scores(self) -> set[str]:
        """Return the set of score names discovered during conversion."""
        return self._discovered_scores

    def get_resolved_columns(self) -> dict[str, str]:
        """Return QPX field -> resolved tool column mappings.

        Subclasses that perform column resolution (e.g. _resolved, _resolved_pg)
        provide their mappings; default is empty dict.
        """
        return getattr(self, "_resolved_pg", None) or getattr(self, "_resolved", None) or {}

    def get_table_columns(self, table_name: str) -> set[str]:
        """Return column names for a DuckDB table (for column resolution).

        Used by orchestrators that need to resolve columns without accessing
        adapter internals.
        """
        rows = self._conn.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name=?",
            [table_name],
        ).fetchall()
        return {r[0] for r in rows}

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    @abstractmethod
    def convert(self, *args, **kwargs) -> None:
        """Run the full conversion pipeline and write output files.

        Subclasses define their own parameter signatures.
        """
        ...

    # ------------------------------------------------------------------
    # Helpers shared by all converters
    # ------------------------------------------------------------------

    def _write_pepmap(self, output_path: str, pepmap: dict[str, tuple[str, list[str]]], creator: str) -> bool:
        """Write ``{peptidoform: (sequence, proteins)}`` as the pepmap view; returns whether there was anything to write."""
        if not pepmap:
            return False
        records = [
            {
                "sequence": sequence,
                "peptidoform": peptidoform,
                "pg_accessions": [
                    {"accession": acc, "start": None, "end": None, "pre": None, "post": None} for acc in accessions
                ],
                "is_unique": len(accessions) == 1,
            }
            for peptidoform, (sequence, accessions) in pepmap.items()
        ]
        with PepMapWriter(output_path, creator=creator) as writer:
            writer.write_batch(records)
        self.logger.info("Wrote %d pepmap entries to %s", len(records), output_path)
        return True

    @staticmethod
    def _load_sdrf_enzyme(sdrf_path: str) -> str | None:
        """Load the first enzyme name from SDRF for missed-cleavage computation."""
        from qpx.core.sdrf import SDRFHandler

        enzymes = SDRFHandler(sdrf_path).get_enzymes()
        return str(enzymes[0]) if enzymes else None

    def _register_precursor_lookup(self, resolved: dict, enzyme_name: str | None, to_proforma, to_modifications) -> dict:
        """Register the distinct precursors of the ``report`` table as ``precursor_lookup``.

        ``resolved`` maps the field names to the report columns; ``to_proforma`` and
        ``to_modifications`` parse the tool's modified-sequence notation. Returns
        ``(modified_seq, sequence, charge) -> modifications``.
        """
        precursors = self._conn.execute(
            sql_build(
                "SELECT DISTINCT $mc, $sc, $cc FROM report",
                mc=validate_identifier(resolved["modified_sequence"]),
                sc=validate_identifier(resolved["sequence"]),
                cc=validate_identifier(resolved["charge"]),
            )
        ).fetchall()
        self.logger.info("Pre-computing %d unique precursors", len(precursors))

        mods_lookup: dict[tuple[str, str, int], list | None] = {}
        rows = []
        for modified_seq, sequence, charge in precursors:
            modified_seq = str(modified_seq)
            sequence = str(sequence)
            charge = int(charge)
            peptidoform = to_proforma(modified_seq)
            _, modifications = to_modifications(modified_seq, sequence)
            calculated_mz = compute_precursor_mz(peptidoform, charge)
            missed = count_missed_cleavages(sequence, enzyme_name) if enzyme_name else None
            mods_lookup[(modified_seq, sequence, charge)] = modifications
            rows.append((modified_seq, sequence, charge, peptidoform, calculated_mz, missed))

        self._conn.execute("DROP TABLE IF EXISTS precursor_lookup")
        self._conn.from_df(
            pd.DataFrame(
                rows,
                columns=["modified_sequence", "sequence", "charge", "peptidoform", "calculated_mz", "missed_cleavages"],
            )
        ).create("precursor_lookup")
        self.logger.info("Registered precursor_lookup temp table in DuckDB")
        return mods_lookup

    @staticmethod
    def _assemble_feature_table(table: pa.Table, target_schema: pa.Schema, mods_array: pa.Array, pg_array: pa.Array) -> pa.Table:
        """Assemble feature columns in schema order, casting or filling nulls."""
        columns = {}
        for field in target_schema:
            if field.name == "modifications":
                columns[field.name] = mods_array
            elif field.name == "pg_accessions":
                columns[field.name] = pg_array
            elif field.name in table.schema.names:
                col = table.column(field.name)
                if col.type != field.type:
                    try:
                        col = col.cast(field.type, safe=False)
                    except (pa.ArrowInvalid, pa.ArrowNotImplementedError):
                        col = pa.array(col.to_pylist(), type=field.type)
                columns[field.name] = col
            else:
                columns[field.name] = pa.nulls(len(table), type=field.type)
        return pa.table(columns, schema=target_schema)

    def _transform_rows_with_skip(
        self,
        rows: Iterable[dict],
        transform_fn: Callable[[dict], dict | list[dict] | None],
        debug_label: str,
        warning_label: str = "",
    ) -> list[dict]:
        """Transform rows, skipping and counting only the supported row errors."""
        records: list[dict] = []
        skipped = 0
        for row in rows:
            try:
                rec = transform_fn(row)
                if isinstance(rec, list):
                    records.extend(rec)
                elif rec is not None:
                    records.append(rec)
            except (ValueError, TypeError, KeyError, IndexError, ArithmeticError) as exc:
                skipped += 1
                self.logger.debug("Skipping %s row: %s", debug_label, exc)
        if skipped:
            total = skipped + len(records)
            self.logger.warning(
                "Skipped %d / %d%s rows (%.1f%%) in batch",
                skipped,
                total,
                warning_label,
                100 * skipped / total,
            )
        return records

    def _table_exists(self, name: str) -> bool:
        """Return True if *name* is registered as a DuckDB table or view."""
        rows = self._conn.execute(
            "SELECT 1 FROM information_schema.tables WHERE table_name = ? LIMIT 1",
            [name],
        ).fetchall()
        return len(rows) > 0

    def _load_parquet(self, table_name: str, path: str) -> None:
        """Load a Parquet file into a DuckDB table."""
        safe = escape_path(path)
        tbl = validate_table(table_name)
        self._conn.execute(
            sql_build(
                "CREATE OR REPLACE TABLE $tbl AS SELECT * FROM read_parquet('$path')",
                tbl=tbl,
                path=safe,
            )
        )
        count = self._conn.execute(sql_build("SELECT COUNT(*) FROM $tbl", tbl=tbl)).fetchall()[0][0]
        self.logger.info(f"Loaded {count:,} rows into '{table_name}' from {path}")

    def _query_batched(self, sql: str, batch_size: int = 500_000):
        """Execute *sql* and yield result batches as Arrow record batches.

        ``to_arrow_reader()`` returns a ``RecordBatchReader`` (iterable of
        ``RecordBatch``) that streams the query instead of materializing it.
        """
        reader = self._conn.sql(sql).to_arrow_reader(batch_size)
        for batch in reader:
            if batch.num_rows > 0:
                yield batch

    # ------------------------------------------------------------------
    # Score tracking
    # ------------------------------------------------------------------

    def _track_scores(self, records: list[dict]) -> None:
        """Collect score names from a batch of records into ``_discovered_scores``.

        Scans both ``additional_scores`` and per-position scores inside
        ``modifications`` so that site localization scores (e.g.
        localization_probability) are included in ontology.parquet.
        """
        for rec in records:
            scores = rec.get("additional_scores")
            if scores:
                for s in scores:
                    name = s.get("score_name")
                    if name:
                        self._discovered_scores.add(name)
            # Also collect scores from modification positions
            modifications = rec.get("modifications")
            if modifications:
                for mod in modifications:
                    positions = mod.get("positions")
                    if positions:
                        for pos in positions:
                            pos_scores = pos.get("scores")
                            if pos_scores:
                                for s in pos_scores:
                                    name = s.get("score_name")
                                    if name:
                                        self._discovered_scores.add(name)

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self):
        """Close the DuckDB connection (only if we created it)."""
        if self._conn and self._owns_conn:
            try:
                self._conn.close()
            except Exception:
                self.logger.debug("Failed to close DuckDB connection", exc_info=True)
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
