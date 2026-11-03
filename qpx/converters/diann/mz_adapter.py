"""DIA-NN m/z adapter — generate the ``mz`` view from feature-matched MS2 spectra."""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc

from qpx.converters.diann.feature_adapter import _ms_info_scans
from qpx.core.engine import create_connection
from qpx.core.sql import escape_path, sql_build
from qpx.writers.mz import MzWriter

logger = logging.getLogger(__name__)

# Rows per Arrow batch fetched from DuckDB, so only a bounded slice of peak arrays is in memory.
_FETCH_ROWS = 500


def _scan_numbers(con: duckdb.DuckDBPyConnection, path: Path) -> pd.DataFrame:
    """Distinct ``scan`` values of a stats table (``raw``) with their scan numbers (``scan``).

    Scan numbers are resolved as the feature view resolves them; IDs without a
    single scan number are dropped.
    """
    raw = con.execute(sql_build("SELECT DISTINCT scan FROM parquet_scan('$path')", path=escape_path(str(path)))).df()["scan"]
    scans = pd.DataFrame({"raw": raw, "scan": _ms_info_scans(raw)}).dropna(subset=["scan"])
    return scans.astype({"scan": "int64"})


def _ms_info_lookup(con: duckdb.DuckDBPyConnection, ms_info: Path) -> dict[int, tuple]:
    """Map scan -> ``(rt, precursor_mz, precursor_charge)`` from one run's ``*_ms_info.parquet``.

    The charge is kept only when positive, as the mzML spectra mapping keeps
    it. A run without that table maps no scan.
    """
    if not ms_info.exists():
        return {}
    rows = con.execute(
        sql_build(
            "SELECT scan, rt, precursor_mz, "
            "CASE WHEN precursor_charge > 0 THEN precursor_charge END "
            "FROM parquet_scan('$path') WHERE rt IS NOT NULL",
            path=escape_path(str(ms_info)),
        )
    ).fetchall()
    scans = _ms_info_scans(pd.Series([row[0] for row in rows], dtype=object))
    return {int(scan): row[1:] for scan, row in zip(scans.tolist(), rows) if pd.notna(scan)}


def _precursors(precursor_mz: float | None, charge: int | None) -> list[dict]:
    """The precursor list of one MS2 spectrum, empty without a precursor m/z.

    Fields follow the mzML spectra mapping. MS-info tables carry no
    isolation-window offsets, so the window bounds stay null. The selected-ion
    intensity also stays null: a DIA window centre is not a selected ion, and
    MS-info replaces a zero mzML intensity with an MS1 purity estimate.
    """
    if precursor_mz is None:
        return []
    return [
        {
            "selected_ion_mz": precursor_mz,
            "selected_ion_charge": charge,
            "selected_ion_intensity": None,
            "isolation_window_target": precursor_mz,
            "isolation_window_lower": None,
            "isolation_window_upper": None,
            "spectrum_ref": None,
        }
    ]


def _float32_peaks(column) -> tuple[np.ndarray, np.ndarray]:
    """Values and offsets of a list column of peak values cast to float32."""
    peaks = pc.cast(column, pa.list_(pa.float32()))
    return peaks.values.to_numpy(zero_copy_only=False), peaks.offsets.to_numpy()


def _matched_spectra(con: duckdb.DuckDBPyConnection, ms2_file: Path, matched_scans: set[int] | None):
    """Yield ``(scan, ms_level, mz_array, intensity_array)`` once per matched scan of one run.

    DuckDB applies the scan filter and rows are fetched as Arrow batches, so only
    the matched spectra are read and peak arrays reach Python as float32 numpy
    arrays. ``matched_scans`` holds the run's scan numbers; ``None`` keeps every scan.
    """
    keep = _scan_numbers(con, ms2_file)
    if matched_scans is not None:
        keep = keep[keep["scan"].isin(matched_scans)]
    con.register("keep_scans", keep.drop_duplicates("scan"))
    try:
        reader = con.sql(
            sql_build(
                "SELECT k.scan, m.ms_level, m.mz_array, m.intensity_array "
                "FROM parquet_scan('$path') m JOIN keep_scans k ON m.scan = k.raw",
                path=escape_path(str(ms2_file)),
            )
        ).to_arrow_reader(_FETCH_ROWS)
        seen: set[int] = set()
        for batch in reader:
            mz, mz_offsets = _float32_peaks(batch.column("mz_array"))
            intensity, intensity_offsets = _float32_peaks(batch.column("intensity_array"))
            scans = batch.column("scan").to_pylist()
            ms_levels = batch.column("ms_level").to_pylist()
            for row, scan in enumerate(scans):
                if scan not in seen:
                    seen.add(scan)
                    yield (
                        scan,
                        ms_levels[row],
                        mz[mz_offsets[row] : mz_offsets[row + 1]],
                        intensity[intensity_offsets[row] : intensity_offsets[row + 1]],
                    )
    finally:
        con.unregister("keep_scans")


class DiannMzAdapter:
    """Generate the ``mz`` view from mzML-statistics MS2 spectra.

    Reads per-run ``*_ms2_info.parquet`` files (``scan``, ``ms_level``,
    ``mz_array``, ``intensity_array``) produced by the quantmsdiann workflow,
    keeping only the feature-matched MS2 scans (``matched_scans``) and writing
    their peak arrays to the ``mz`` / ``intensity`` float32 list columns.
    ``scan_start_time`` and ``precursors`` are filled from the per-run
    ``*_ms_info.parquet`` (``rt`` and the ``precursor_*`` columns).
    """

    def __init__(self, batch_size: int = 100_000):
        self._batch_size = batch_size

    def convert(
        self,
        mzml_statistics_folder: str,
        output_path: str,
        matched_scans: set[tuple[str, int]] | None = None,
    ) -> int:
        """Write ``{output_path}`` from ``*_ms2_info.parquet`` / ``*_ms_info.parquet``.

        Returns the number of unique feature-matched spectra written.
        """
        stats = Path(mzml_statistics_folder)
        ms2_files = sorted(stats.glob("*_ms2_info.parquet"))
        if not ms2_files:
            raise FileNotFoundError(f"no *_ms2_info.parquet files in {mzml_statistics_folder}")

        scans_by_run: dict[str, set[int]] | None = None
        if matched_scans is not None:
            scans_by_run = defaultdict(set)
            for run, scan in matched_scans:
                scans_by_run[run].add(scan)

        con = create_connection()

        pending: list[dict] = []
        written = 0
        writer = MzWriter(output_path, batch_size=self._batch_size)
        writer.add_footer_metadata(
            {b"spectrum_content": b"selected_ms_levels" if scans_by_run is None else b"identified_spectra"}
        )
        try:
            for ms2_file in ms2_files:
                run_name = ms2_file.name.removesuffix("_ms2_info.parquet")
                run_info = _ms_info_lookup(con, stats / f"{run_name}_ms_info.parquet")
                run_scans = None if scans_by_run is None else scans_by_run.get(run_name, set())
                for scan, ms_level, mz, intensity in _matched_spectra(con, ms2_file, run_scans):
                    rt, *precursor = run_info.get(scan, (None, None, None))
                    pending.append(self.build_record(run_name, scan, ms_level, mz, intensity, rt, _precursors(*precursor)))
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
        logger.info("wrote %s unique feature-matched spectra to %s", written, output_path)
        return written

    @staticmethod
    def build_record(
        run_name: str,
        scan: int,
        ms_level: int,
        mz,
        intensity,
        rt: float | None = None,
        precursors: list[dict] | None = None,
    ) -> dict:
        """Build one ``mz``-view record with float32 peak arrays."""
        mz = np.asarray(mz, dtype=np.float32)
        intensity = np.asarray(intensity, dtype=np.float32)
        return {
            "id": f"{run_name}:scan={scan}",
            "run_file_name": run_name,
            "scan": [scan],
            "ms_level": int(ms_level) if ms_level is not None else 2,
            "centroid": True,
            "scan_start_time": (float(rt) / 60.0) if rt is not None else 0.0,
            "inverse_ion_mobility": None,
            "ion_injection_time": None,
            "total_ion_current": float(intensity.sum()) if len(intensity) else 0.0,
            "precursors": precursors or [],
            "mz": mz,
            "intensity": intensity,
            "cv_params": [],
        }
