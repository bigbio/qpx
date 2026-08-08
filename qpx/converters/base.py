"""BaseConverter -- abstract base for all tool-specific converters.

Each converter creates its own DuckDB in-memory connection, loads
tool-specific output, transforms it via SQL into the QPX schema, and pipes
rows into the appropriate Writer.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path

import duckdb

from qpx.core.engine import create_converter_connection
from qpx.core.scores import score_ontology_entries
from qpx.core.sql import escape_path, sql_build, validate_table


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

    After conversion, call ``write_score_ontology()`` to emit ontology
    entries for all discovered score names.
    """

    def __init__(
        self,
        duckdb_memory: str = "16GB",
        duckdb_threads: int = 4,
        conn: duckdb.DuckDBPyConnection | None = None,
        compression: str = "zstd",
    ):
        self.logger = logging.getLogger(f"{self.__class__.__module__}.{self.__class__.__name__}")
        self._compression = compression
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

    @staticmethod
    def _escape_path(path: str) -> str:
        """Escape a filesystem path for safe embedding in DuckDB SQL.

        Doubles any single-quote characters so that paths like
        ``O'Brien/data.txt`` don't break or inject SQL.
        """
        return str(path).replace("'", "''")

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
        count = self._conn.execute(sql_build("SELECT COUNT(*) FROM $tbl", tbl=tbl)).fetchone()[0]
        self.logger.info(f"Loaded {count:,} rows into '{table_name}' from {path}")

    def _query_batched(self, sql: str, batch_size: int = 500_000):
        """Execute *sql* and yield result batches as Arrow record batches.

        DuckDB >= 1.0 ``fetch_record_batch()`` returns a
        ``RecordBatchReader`` (iterable of ``RecordBatch``).
        """
        result = self._conn.execute(sql)
        reader = result.fetch_record_batch(batch_size)
        for batch in reader:
            if batch.num_rows > 0:
                yield batch

    # ------------------------------------------------------------------
    # Score tracking and ontology writing
    # ------------------------------------------------------------------

    def _track_scores(self, records: list[dict]) -> None:
        """Collect score names from a batch of records into ``_discovered_scores``.

        Scans both ``additional_scores`` and per-position scores inside
        ``modifications`` so that phospho site localization scores (e.g.
        phospho_sty_probability) are included in ontology.parquet.
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

    def write_score_ontology(
        self,
        output_path: str | Path,
        view: str = "psm",
        extra_entries: list[dict] | None = None,
    ) -> Path | None:
        """Write ontology.parquet with CV term entries for discovered scores.

        Args:
            output_path: Path for the ontology Parquet file.
            view: QPX view name (``"psm"``, ``"feature"``, ``"pg"``).
            extra_entries: Additional ontology entries to include.

        Returns:
            Path to the written file, or ``None`` if no entries to write.
        """
        entries = score_ontology_entries(self._discovered_scores, view=view)
        if extra_entries:
            entries.extend(extra_entries)
        if not entries:
            return None

        from qpx.converters.ontology_util import dedupe_ontology_entries
        from qpx.writers.ontology import OntologyWriter

        entries = dedupe_ontology_entries(entries)
        output_path = Path(output_path)
        with OntologyWriter(output_path, creator="qpx", compression=self._compression) as writer:
            writer.write_batch(entries)
        self.logger.info(
            "Wrote %d ontology entries (%d scores) to %s",
            len(entries),
            len(self._discovered_scores),
            output_path,
        )
        return output_path

    def write_ontology(
        self,
        output_path: str | Path,
        view: str = "psm",
        resolved_mappings: dict[str, str] | None = None,
        tool_name: str | None = None,
        extra_entries: list[dict] | None = None,
    ) -> Path | None:
        """Write ontology.parquet with score + field provenance entries.

        Combines:
        1. Score ontology entries (from _discovered_scores)
        2. Field provenance entries (from resolved column mappings)
        3. Any extra entries passed by the converter

        Args:
            output_path: Path for the ontology Parquet file.
            view: QPX view name.
            resolved_mappings: QPX field -> resolved tool column name.
            tool_name: Tool name string (e.g. "DIA-NN").
            extra_entries: Additional ontology entries to include.

        Returns:
            Path to the written file, or None if no entries to write.
        """
        from qpx.core.scores import field_ontology_entries

        entries = score_ontology_entries(self._discovered_scores, view=view)
        if resolved_mappings:
            entries.extend(
                field_ontology_entries(
                    view=view,
                    resolved_mappings=resolved_mappings,
                    tool_name=tool_name,
                )
            )
        if extra_entries:
            entries.extend(extra_entries)
        if not entries:
            return None

        from qpx.converters.ontology_util import dedupe_ontology_entries
        from qpx.writers.ontology import OntologyWriter

        entries = dedupe_ontology_entries(entries)
        output_path = Path(output_path)
        with OntologyWriter(output_path, creator="qpx", compression=self._compression) as writer:
            writer.write_batch(entries)
        self.logger.info(
            "Wrote %d ontology entries (%d scores, %d field mappings) to %s",
            len(entries),
            len(self._discovered_scores),
            len(resolved_mappings or {}),
            output_path,
        )
        return output_path

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
