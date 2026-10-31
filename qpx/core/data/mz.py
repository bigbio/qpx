"""MzSpectra data structure — mass spectrometry spectral data."""

from qpx.core.convert import QueryResult
from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.sql import sql_build

MzSchema = load_schema("mz")


class MzSpectra(BaseStructure):
    """Mass spectrometry spectral data (scan-level)."""

    _schema_class = MzSchema

    def ms1(self) -> "MzSpectra":
        """Filter to MS1 scans only."""
        return self.filter("ms_level = 1")

    def ms2(self) -> "MzSpectra":
        """Filter to MS2 scans only."""
        return self.filter("ms_level = 2")

    def run_summary(self) -> QueryResult:
        """Per run: spectra per MS level, total ion current, injection time and the covered time (minutes)."""
        stmt = sql_build(
            """SELECT run_file_name,
                      COUNT(*) FILTER (WHERE ms_level = 1) AS ms1_spectra,
                      COUNT(*) FILTER (WHERE ms_level = 2) AS ms2_spectra,
                      SUM(total_ion_current) AS total_ion_current,
                      AVG(ion_injection_time) AS mean_ion_injection_time,
                      MIN(scan_start_time) AS first_scan_time,
                      MAX(scan_start_time) AS last_scan_time
               FROM $src GROUP BY run_file_name ORDER BY run_file_name""",
            src=self._query.source,
        )
        return QueryResult(self._engine.connection, stmt)

    def by_rt_range(self, rt_start: float, rt_end: float) -> "MzSpectra":
        """Filter scans by retention time range (minutes)."""
        return self.filter(f"scan_start_time >= {rt_start} AND scan_start_time <= {rt_end}")
