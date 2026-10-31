"""Lightweight MGF (Mascot Generic Format) parser for spectrum lookup by scan number.

Parses MGF files (plain or gzip-compressed) and builds an in-memory index
mapping scan numbers to spectra data (mz_array, intensity_array).  Designed
for attaching identified spectra to PSM records during mzIdentML conversion.

The parser streams through the file line-by-line so it never loads the
entire file into memory at once -- only the accumulated index is kept.
"""

from __future__ import annotations

import gzip
import logging
import re
from array import array
from pathlib import Path
from typing import IO

logger = logging.getLogger(__name__)

# Regex patterns for extracting scan numbers
_SCAN_PATTERNS = [
    re.compile(r"(?:scan|index|spectrum)=(\d+)", re.IGNORECASE),
    re.compile(r"\bScan\s+(\d+)\b"),  # "Scan 154" (Distiller style)
]
_BARE_INT = re.compile(r"^(\d+)$")
_POSITION_TITLE = re.compile(r"index=\d+", re.IGNORECASE)


def _extract_scan(title: str, scans_value: str | None) -> int | None:
    """Extract a scan number from TITLE string and/or SCANS= header value.

    Priority:
      1. ``SCANS=`` header value (if present and numeric)
      2. ``scan=N``, ``index=N``, ``spectrum=N`` in TITLE
      3. ``Scan N`` (case-sensitive, Distiller style) in TITLE
      4. Bare integer TITLE
      5. None (no scan found)
    """
    # 1. Try SCANS= header value first (most reliable)
    if scans_value is not None:
        scans_value = scans_value.strip()
        if scans_value.isdigit():
            return int(scans_value)

    # 2-3. Try regex patterns on TITLE
    for pat in _SCAN_PATTERNS:
        m = pat.search(title)
        if m:
            return int(m.group(1))

    # 4. Bare integer TITLE
    m = _BARE_INT.match(title.strip())
    if m:
        return int(m.group(1))

    return None


class MgfSpectraIndex:
    """Index MGF file by scan number for O(1) spectrum lookup.

    Parameters
    ----------
    mgf_path : str or Path
        Path to an MGF file (plain-text ``.mgf`` or gzip-compressed ``.mgf.gz``).
    keep_peaks : bool
        Hold each spectrum's peaks; without them only its scan number and RT are kept.

    Examples
    --------
    >>> idx = MgfSpectraIndex("spectra.mgf.gz")
    >>> spec = idx.get_spectrum(1234)
    >>> spec["mz_array"]
    [200.123, 201.456, 202.789]
    """

    def __init__(self, mgf_path: str | Path, keep_peaks: bool = True) -> None:
        self.path = Path(mgf_path)
        self.keep_peaks = keep_peaks
        self._index: dict[int, dict] = {}
        self._by_position: dict[int, dict] = {}
        self._position_counter: int = 0
        self._parse()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_spectrum(self, scan: int) -> dict | None:
        """Return spectrum data for *scan*, or ``None`` if not indexed.

        Returns
        -------
        dict or None
            ``{"mz_array": list[float], "intensity_array": list[float]}``
        """
        return self._as_lists(self._index.get(scan))

    def get_spectrum_by_index(self, index: int) -> dict | None:
        """Return spectrum data by 0-based position in the MGF file.

        Useful as a fallback when scan-based lookup fails (e.g., when
        mzIdentML uses ``index=N`` spectrumIDs that correspond to the
        Nth spectrum in the MGF file).
        """
        return self._as_lists(self._by_position.get(index))

    @staticmethod
    def _as_lists(spectrum: dict | None) -> dict | None:
        """The spectrum with plain lists; compact arrays hold every spectrum of the file in memory."""
        if spectrum is None or spectrum["mz_array"] is None:
            return spectrum
        return {**spectrum, "mz_array": spectrum["mz_array"].tolist(), "intensity_array": spectrum["intensity_array"].tolist()}

    def __len__(self) -> int:
        """Number of indexed spectra."""
        return len(self._index)

    def __contains__(self, scan: int) -> bool:
        """Check if *scan* is in the index."""
        return scan in self._index

    @property
    def scans(self) -> list[int]:
        """Sorted list of all indexed scan numbers."""
        return sorted(self._index)

    # ------------------------------------------------------------------
    # Parsing internals
    # ------------------------------------------------------------------

    def _parse(self) -> None:
        """Stream-parse the MGF file and populate ``self._index``."""
        opener = gzip.open if self.path.suffix == ".gz" else open
        with opener(self.path, "rt", encoding="utf-8", errors="replace") as fh:
            self._parse_stream(fh)
        logger.info(
            "Indexed %d spectra from %s",
            len(self._index),
            self.path.name,
        )

    def _parse_stream(self, fh: IO[str]) -> None:
        """Parse an open text stream of MGF data."""
        in_spectrum = False
        title: str | None = None
        scans_value: str | None = None
        rt_value: float | None = None
        mz_list: list[float] = []
        int_list: list[float] = []

        for line in fh:
            line = line.rstrip("\n\r")

            if line.startswith("BEGIN IONS"):
                in_spectrum = True
                title = None
                scans_value = None
                rt_value = None
                mz_list = []
                int_list = []
                continue

            if line.startswith("END IONS"):
                if in_spectrum:
                    self._store_spectrum(title, scans_value, mz_list, int_list, rt_value)
                in_spectrum = False
                continue

            if not in_spectrum:
                continue

            # Header lines inside a spectrum block
            if line.startswith("TITLE="):
                title = line[6:]
                continue
            if line.startswith("SCANS="):
                scans_value = line[6:]
                continue
            if line.startswith("RTINSECONDS="):
                try:
                    rt_value = float(line[12:].strip())
                except (ValueError, TypeError):
                    pass
                continue
            # Skip other header lines (PEPMASS, CHARGE, etc.)
            if "=" in line and not line[0].isdigit():
                continue

            # Peak data: "mz intensity" (or possibly "mz\tintensity")
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mz = float(parts[0])
                    intensity = float(parts[1])
                    mz_list.append(mz)
                    int_list.append(intensity)
                except ValueError:
                    continue

    def _store_spectrum(
        self,
        title: str | None,
        scans_value: str | None,
        mz_list: list[float],
        int_list: list[float],
        rt_value: float | None = None,
    ) -> None:
        """Extract scan number and store spectrum in the index."""
        pos = self._position_counter
        self._position_counter += 1

        scan = _extract_scan(title or "", scans_value)
        spectrum = {
            "mz_array": array("d", mz_list) if self.keep_peaks else None,
            "intensity_array": array("d", int_list) if self.keep_peaks else None,
            "rt": rt_value,
            # A TITLE with index=N names a file position, not an acquisition scan.
            "scan": None if scans_value is None and _POSITION_TITLE.search(title or "") else scan,
        }

        # Always store by position
        self._by_position[pos] = spectrum

        # Also store by scan if parseable
        if scan is None:
            logger.debug("Spectrum at position %d has no parseable scan: TITLE=%s", pos, title)
            return

        if scan in self._index:
            logger.debug("Duplicate scan %d — keeping first occurrence", scan)
            return

        self._index[scan] = spectrum  # Same dict object, not a copy
