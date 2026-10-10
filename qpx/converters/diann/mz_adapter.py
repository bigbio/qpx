"""DIA-NN m/z adapter — generate the ``mz`` view from mzML-derived spectra."""

from __future__ import annotations

import logging
from pathlib import Path

import duckdb
import numpy as np

from qpx.core.spectra_codec import encode_intensity, encode_mz
from qpx.core.sql import escape_path, sql_build
from qpx.writers.mz import MzWriter

logger = logging.getLogger(__name__)


class DiannMzAdapter:
    """Generate the ``mz`` view from AMS_info (mzML-parsed) fragment spectra.

    Reads per-run ``*_AMS_info.parquet`` files (``run_name``, ``scan``,
    ``precursor_rt``, ``mz_array``, ``intensity_array``), deduplicates by
    ``(run_name, scan)``, losslessly encodes the peak arrays into the ``mz`` /
    ``intensity`` binary columns, and writes the ``mz`` view. Only the
    feature-matched MS2 scans are kept when ``matched_scans`` is provided.
    """

    def __init__(self, compression: str = "zstd", batch_size: int = 100_000):
        self._compression = compression
        self._batch_size = batch_size

    def convert(
        self,
        ams_info_folder: str,
        output_path: str,
        matched_scans: set[tuple[str, int]] | None = None,
    ) -> int:
        """Write ``{output_path}`` from a folder of ``*_AMS_info.parquet`` files.

        Returns the number of unique spectra written.
        """
        ams_files = sorted(Path(ams_info_folder).glob("*_AMS_info.parquet"))
        if not ams_files:
            raise FileNotFoundError(f"no *_AMS_info.parquet files in {ams_info_folder}")

        seen: set[tuple[str, int]] = set()
        pending: list[dict] = []
        written = 0

        con = duckdb.connect()
        writer = MzWriter(output_path, compression=self._compression)
        try:
            for ams_file in ams_files:
                safe_path = escape_path(str(ams_file))
                cursor = con.execute(
                    sql_build(
                        "SELECT run_name, scan, precursor_rt, mz_array, intensity_array FROM parquet_scan('$path')",
                        path=safe_path,
                    )
                )
                for run_name, scan, precursor_rt, mz, intensity in cursor.fetchall():
                    scan0 = int(scan[0]) if isinstance(scan, (list, tuple)) else int(scan)
                    key = (run_name, scan0)
                    if key in seen:
                        continue
                    if matched_scans is not None and key not in matched_scans:
                        continue
                    seen.add(key)
                    pending.append(self.build_record(run_name, scan0, precursor_rt, mz, intensity))
                    if len(pending) >= self._batch_size:
                        writer.write_batch(pending)
                        written += len(pending)
                        pending = []
            if pending:
                writer.write_batch(pending)
                written += len(pending)
            writer.close()
        finally:
            con.close()
        logger.info("wrote %s unique spectra to %s", written, output_path)
        return written

    @staticmethod
    def build_record(run_name: str, scan: int, precursor_rt, mz, intensity) -> dict:
        """Build one ``mz``-view record with losslessly-encoded peak arrays."""
        mz = np.asarray(mz, dtype=np.float32)
        intensity = np.asarray(intensity, dtype=np.float32)
        return {
            "id": f"{run_name}:scan={scan}",
            "run_file_name": run_name,
            "scan": scan,
            "ms_level": 2,
            "centroid": True,
            "scan_start_time": float(precursor_rt) / 60.0,
            "inverse_ion_mobility": None,
            "ion_injection_time": None,
            "total_ion_current": None,
            "precursors": [],
            "mz": encode_mz(mz),
            "intensity": encode_intensity(intensity),
            "cv_params": [],
        }
