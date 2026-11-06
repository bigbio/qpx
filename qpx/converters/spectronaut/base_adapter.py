"""Shared base adapter for Spectronaut converters.

Provides common data-loading helpers used by both
:class:`SpectronautFeatureAdapter` and :class:`SpectronautPgAdapter`.
"""

from __future__ import annotations

from qpx.converters.base import BaseConverter
from qpx.converters.spectronaut.report_layout import (
    long_columns,
    melted_report_sql,
    pivot_columns,
    report_header,
    unavailable_views,
)
from qpx.core.sql import escape_path, sql_build, validate_identifier


def _detect_decimal_separator(path: str, sample_lines: int = 20) -> str:
    """Auto-detect decimal separator by inspecting the first data lines.

    European Spectronaut exports use ``,`` (comma) as decimal separator,
    while English exports use ``.`` (period).  We look at a known numeric
    column (``EG.Qvalue``) to decide.

    Returns:
        ``','`` or ``'.'``.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().strip().split("\t")
        try:
            qv_idx = header.index("EG.Qvalue")
        except ValueError:
            return "."
        for _ in range(sample_lines):
            line = fh.readline()
            if not line:
                break
            fields = line.strip().split("\t")
            if qv_idx < len(fields):
                val = fields[qv_idx].strip()
                if "," in val:
                    return ","
    return "."


class SpectronautBaseAdapter(BaseConverter):
    """Base adapter with Spectronaut-specific loading utilities.

    Subclasses inherit ``_load_spectronaut_report()`` so the logic is defined
    once rather than duplicated across the feature and PG adapters.
    """

    def _load_spectronaut_report(self, path: str) -> None:
        """Create a DuckDB view over a Spectronaut report file.

        Uses ``CREATE VIEW`` so DuckDB reads from the file lazily with
        column pruning and predicate pushdown.

        Supports both TSV and Parquet formats.  If the ``report`` view
        already exists it is skipped.  Decimal separator is auto-detected
        for TSV files.  A pivot report (one column per run) is melted to the
        long layout in a temporary table.

        Args:
            path: Filesystem path to the Spectronaut report file.
        """
        if self._table_exists("report"):
            self.logger.debug("report view already loaded -- skipping reload")
            return
        safe_path = escape_path(path)
        header = report_header(path)
        if pivot_columns(header):
            self._load_pivot_report(path, header)
        elif path.endswith(".parquet"):
            self._conn.execute(
                sql_build(
                    "CREATE VIEW report AS SELECT * FROM read_parquet('$path')",
                    path=safe_path,
                )
            )
        else:
            dec_sep = _detect_decimal_separator(path)
            self.logger.info("Detected decimal separator: '%s'", dec_sep)
            self._conn.execute(
                sql_build(
                    """CREATE VIEW report AS
                SELECT * FROM read_csv_auto('$path',
                    delim='\t', header=true, auto_detect=true,
                    null_padding=true, decimal_separator='$dec')""",
                    path=safe_path,
                    dec=dec_sep,
                )
            )
        count = self._conn.execute("SELECT COUNT(*) FROM report").fetchall()[0][0]
        self.logger.info("Spectronaut report view created (%d rows)", count)

    def _load_pivot_report(self, path: str, header: list[str]) -> None:
        """Melt a pivot report into the ``report`` table; a quantity column can hold ``Filtered``, so all are read as text."""
        safe_path = escape_path(path)
        reader = (
            f"read_parquet('{safe_path}')"
            if path.endswith(".parquet")
            else f"read_csv('{safe_path}', delim='\t', header=true, all_varchar=true, null_padding=true)"
        )
        self.logger.info("Melting Spectronaut pivot report: %d runs", len(pivot_columns(header)))
        self._conn.execute(f"CREATE TEMP TABLE report_pivot AS SELECT * FROM {reader}")
        self._conn.execute(f"CREATE TEMP TABLE report AS {melted_report_sql(header, 'report_pivot')}")
        self._conn.execute("DROP TABLE report_pivot")

    def _require_view(self, path: str, view: str) -> None:
        """Raise a clear error when the report at *path* cannot produce *view*."""
        reason = unavailable_views(long_columns(report_header(path))).get(view)
        if reason:
            raise ValueError(f"Cannot convert {path} to the {view} view: {reason}")

    def _get_report_columns(self) -> set[str]:
        """Return column names in the DuckDB report view."""
        return {
            c[0]
            for c in self._conn.execute("SELECT column_name FROM information_schema.columns WHERE table_name='report'").fetchall()
        }

    def _query_distinct_runs(self, run_col: str) -> list[str]:
        """Query distinct values of a run column from the report view."""
        qcol = validate_identifier(run_col)
        rows = self._conn.execute(
            sql_build(
                "SELECT DISTINCT $col FROM report ORDER BY $col",
                col=qcol,
            )
        ).fetchall()
        return [r[0] for r in rows]
