"""Xic data structure — extracted-ion chromatograms (XIC)."""

from __future__ import annotations

from qpx.core.convert import QueryResult
from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.query import _escape_sql_string
from qpx.core.sql import sql_build

XicSchema = load_schema("xic")  # pylint: disable=invalid-name


class Xic(BaseStructure):
    """Extracted-ion chromatograms, one row per traced precursor/fragment ion.

    A row stores a whole trace as aligned ``rt`` and ``intensity`` arrays; use
    :meth:`points` to expand it to the long (one point per row) form used by
    plotting and integration code.
    """

    _schema_class = XicSchema

    def by_run(self, run_file_name: str) -> "Xic":
        """Filter traces by raw file name."""
        return self.filter(f"run_file_name = '{_escape_sql_string(run_file_name)}'")

    def by_precursor(self, precursor_id: str) -> "Xic":
        """Filter traces by the producer's precursor identifier."""
        return self.filter(f"precursor_id = '{_escape_sql_string(precursor_id)}'")

    def by_peptidoform(self, peptidoform: str) -> "Xic":
        """Filter traces by ProForma peptidoform."""
        return self.filter(f"peptidoform = '{_escape_sql_string(peptidoform)}'")

    def ms1(self) -> "Xic":
        """Filter to precursor (MS1) traces."""
        return self.filter("transition_type = 'ms1'")

    def fragments(self) -> "Xic":
        """Filter to fragment-ion (MS2) traces."""
        return self.filter("transition_type = 'fragment'")

    def points(self) -> QueryResult:
        """Expand each trace to one row per retention-time point.

        Returns one row per retention-time point with the trace identity columns
        plus ``trace_rt`` and ``trace_intensity``. The two arrays are positionally
        aligned, so DuckDB's zipped ``UNNEST`` keeps each intensity with its
        retention time.
        """
        stmt = sql_build(
            """SELECT xic_id, run_file_name, label, precursor_id, peptidoform,
                      charge, feature_id, transition, transition_type,
                      UNNEST(rt) AS trace_rt,
                      UNNEST(intensity) AS trace_intensity
               FROM $src""",
            src=self._query.source,
        )
        return QueryResult(self._engine.connection, stmt)
