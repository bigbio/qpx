"""MzSpectra data structure — mass spectrometry spectral data."""

from __future__ import annotations

import numpy as np
import pandas as pd

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.spectra_codec import decode_intensity, decode_mz

MzSchema = load_schema("mz")


class MzSpectra(BaseStructure):
    """Mass spectrometry spectral data (scan-level).

    The ``mz`` and ``intensity`` columns hold losslessly-encoded blobs (see
    :mod:`qpx.core.spectra_codec`); decode them back to ``float32`` arrays with
    :meth:`decode_spectrum`.
    """

    _schema_class = MzSchema

    def ms1(self) -> "MzSpectra":
        """Filter to MS1 scans only."""
        return self.filter("ms_level = 1")

    def ms2(self) -> "MzSpectra":
        """Filter to MS2 scans only."""
        return self.filter("ms_level = 2")

    def by_rt_range(self, rt_start: float, rt_end: float) -> "MzSpectra":
        """Filter scans by retention time range (minutes)."""
        return self.filter(f"scan_start_time >= {rt_start} AND scan_start_time <= {rt_end}")

    def decode_spectrum(self, row: dict | pd.Series) -> tuple[np.ndarray, np.ndarray]:
        """Decode one row into ``(mz, intensity)`` ``float32`` arrays.

        ``mz`` holds the uint32-delta+varint blob; ``intensity`` holds the
        bitshuffled blob. Both decode bit-exactly back to the original arrays.
        """
        mz_blob = row.get("mz")
        intensity_blob = row.get("intensity")
        n_peaks = len(intensity_blob) // 4  # bitshuffle is 4 bytes per peak
        return (
            decode_mz(bytes(mz_blob), n_peaks),
            decode_intensity(bytes(intensity_blob), n_peaks),
        )
