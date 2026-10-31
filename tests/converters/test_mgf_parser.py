"""Tests for the MGF spectra parser."""

import gzip
from pathlib import Path

import pytest

from qpx.converters.mzidentml.mgf_parser import MgfSpectraIndex, _extract_scan

# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

EXAMPLES_DIR = Path(__file__).parent.parent / "examples" / "mzidentml"
PXD054720_MGF = EXAMPLES_DIR / "PXD054720" / "F001234.mgf.gz"

SYNTHETIC_MGF = """\
BEGIN IONS
TITLE=F001234.mgf scan=1234
PEPMASS=456.789 1234.56
CHARGE=2+
200.123 5000.0
201.456 3000.0
202.789 1500.0
END IONS

BEGIN IONS
TITLE=F001234.mgf index=5678
PEPMASS=789.012 9876.54
CHARGE=3+
300.111 8000.0
301.222 6000.0
END IONS

BEGIN IONS
TITLE=spectrum=42
PEPMASS=500.0 1000.0
CHARGE=1+
400.0 100.0
END IONS

BEGIN IONS
TITLE=99
PEPMASS=600.0 2000.0
CHARGE=2+
500.0 200.0
END IONS
"""

DISTILLER_MGF = """\
BEGIN IONS
CHARGE=3+
PEPMASS=758.2215 36576.469 3+
RAWSCANS=sn154
RTINSECONDS=1137.0618
SCANS=154
TITLE=1: Scan 154 (rt=18.951) [G:\\MSData\\file.raw]
150.509360 26.36
299.060930 150.7
END IONS

BEGIN IONS
CHARGE=2+
PEPMASS=400.123 5000.0 2+
SCANS=353
TITLE=2: Scan 353 (rt=20.19) [G:\\MSData\\file.raw]
200.0 100.0
300.0 200.0
400.0 300.0
END IONS
"""


# ---------------------------------------------------------------------------
# Unit tests: _extract_scan helper
# ---------------------------------------------------------------------------


class TestExtractScan:
    def test_scan_equals(self):
        assert _extract_scan("F001234.mgf scan=1234", None) == 1234

    def test_index_equals(self):
        assert _extract_scan("some title index=5678", None) == 5678

    def test_spectrum_equals(self):
        assert _extract_scan("spectrum=42", None) == 42

    def test_scan_equals_case_insensitive(self):
        assert _extract_scan("SCAN=999", None) == 999

    def test_distiller_scan_format(self):
        assert _extract_scan("1: Scan 154 (rt=18.951) [file.raw]", None) == 154

    def test_bare_integer(self):
        assert _extract_scan("99", None) == 99

    def test_scans_header_priority(self):
        """SCANS= header value takes priority over TITLE patterns."""
        assert _extract_scan("some title scan=100", "200") == 200

    def test_scans_header_alone(self):
        assert _extract_scan("no scan here", "154") == 154

    def test_no_scan_found(self):
        assert _extract_scan("no_parseable_value", None) is None

    def test_empty_title(self):
        assert _extract_scan("", None) is None


# ---------------------------------------------------------------------------
# Unit tests: synthetic MGF string (plain text)
# ---------------------------------------------------------------------------


class TestMgfParserSynthetic:
    """Parse a small synthetic MGF string to validate core logic."""

    @pytest.fixture()
    def mgf_index(self, tmp_path: Path) -> MgfSpectraIndex:
        mgf_file = tmp_path / "test.mgf"
        mgf_file.write_text(SYNTHETIC_MGF)
        return MgfSpectraIndex(mgf_file)

    def test_len(self, mgf_index: MgfSpectraIndex):
        assert len(mgf_index) == 4

    def test_contains(self, mgf_index: MgfSpectraIndex):
        assert 1234 in mgf_index
        assert 5678 in mgf_index
        assert 42 in mgf_index
        assert 99 in mgf_index
        assert 9999 not in mgf_index

    def test_get_spectrum_scan_equals(self, mgf_index: MgfSpectraIndex):
        spec = mgf_index.get_spectrum(1234)
        assert spec is not None
        assert spec["mz_array"] == pytest.approx([200.123, 201.456, 202.789])
        assert spec["intensity_array"] == pytest.approx([5000.0, 3000.0, 1500.0])

    def test_get_spectrum_index_equals(self, mgf_index: MgfSpectraIndex):
        spec = mgf_index.get_spectrum(5678)
        assert spec is not None
        assert len(spec["mz_array"]) == 2
        assert spec["mz_array"] == pytest.approx([300.111, 301.222])
        assert spec["intensity_array"] == pytest.approx([8000.0, 6000.0])

    def test_get_spectrum_spectrum_equals(self, mgf_index: MgfSpectraIndex):
        spec = mgf_index.get_spectrum(42)
        assert spec is not None
        assert spec["mz_array"] == pytest.approx([400.0])

    def test_get_spectrum_bare_int(self, mgf_index: MgfSpectraIndex):
        spec = mgf_index.get_spectrum(99)
        assert spec is not None
        assert spec["mz_array"] == pytest.approx([500.0])

    def test_get_spectrum_missing(self, mgf_index: MgfSpectraIndex):
        assert mgf_index.get_spectrum(9999) is None

    def test_scans_property(self, mgf_index: MgfSpectraIndex):
        assert mgf_index.scans == [42, 99, 1234, 5678]


# ---------------------------------------------------------------------------
# Unit tests: Distiller-style TITLE + SCANS= header
# ---------------------------------------------------------------------------


class TestMgfParserDistillerFormat:
    """Parse Distiller-style MGF with SCANS= header."""

    @pytest.fixture()
    def mgf_index(self, tmp_path: Path) -> MgfSpectraIndex:
        mgf_file = tmp_path / "distiller.mgf"
        mgf_file.write_text(DISTILLER_MGF)
        return MgfSpectraIndex(mgf_file)

    def test_len(self, mgf_index: MgfSpectraIndex):
        assert len(mgf_index) == 2

    def test_scan_from_scans_header(self, mgf_index: MgfSpectraIndex):
        assert 154 in mgf_index
        assert 353 in mgf_index

    def test_spectrum_data(self, mgf_index: MgfSpectraIndex):
        spec = mgf_index.get_spectrum(154)
        assert spec is not None
        assert spec["mz_array"] == pytest.approx([150.509360, 299.060930])
        assert spec["intensity_array"] == pytest.approx([26.36, 150.7])


# ---------------------------------------------------------------------------
# Unit tests: gzip-compressed MGF
# ---------------------------------------------------------------------------


class TestMgfParserGzip:
    """Validate that gzip-compressed MGF files are handled correctly."""

    @pytest.fixture()
    def mgf_index(self, tmp_path: Path) -> MgfSpectraIndex:
        mgf_gz_file = tmp_path / "test.mgf.gz"
        with gzip.open(mgf_gz_file, "wt", encoding="utf-8") as fh:
            fh.write(SYNTHETIC_MGF)
        return MgfSpectraIndex(mgf_gz_file)

    def test_len(self, mgf_index: MgfSpectraIndex):
        assert len(mgf_index) == 4

    def test_spectra_match_plain(self, mgf_index: MgfSpectraIndex):
        """Gzip parsing should produce identical results to plain parsing."""
        spec = mgf_index.get_spectrum(1234)
        assert spec is not None
        assert spec["mz_array"] == pytest.approx([200.123, 201.456, 202.789])
        assert spec["intensity_array"] == pytest.approx([5000.0, 3000.0, 1500.0])

    def test_contains(self, mgf_index: MgfSpectraIndex):
        assert 1234 in mgf_index
        assert 5678 in mgf_index
        assert 42 in mgf_index
        assert 99 in mgf_index


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


class TestMgfParserEdgeCases:
    def test_empty_mgf(self, tmp_path: Path):
        mgf_file = tmp_path / "empty.mgf"
        mgf_file.write_text("")
        idx = MgfSpectraIndex(mgf_file)
        assert len(idx) == 0

    def test_spectrum_with_no_peaks(self, tmp_path: Path):
        mgf = "BEGIN IONS\nTITLE=scan=1\nPEPMASS=500.0\nCHARGE=2+\nEND IONS\n"
        mgf_file = tmp_path / "no_peaks.mgf"
        mgf_file.write_text(mgf)
        idx = MgfSpectraIndex(mgf_file)
        assert 1 in idx
        spec = idx.get_spectrum(1)
        assert spec["mz_array"] == []
        assert spec["intensity_array"] == []

    def test_spectrum_with_unparseable_title(self, tmp_path: Path):
        mgf = "BEGIN IONS\nTITLE=no_scan_here\nPEPMASS=500.0\n100.0 200.0\nEND IONS\n"
        mgf_file = tmp_path / "bad_title.mgf"
        mgf_file.write_text(mgf)
        idx = MgfSpectraIndex(mgf_file)
        assert len(idx) == 0  # spectrum skipped

    def test_duplicate_scan_keeps_first(self, tmp_path: Path):
        mgf = "BEGIN IONS\nTITLE=scan=1\n100.0 1000.0\nEND IONS\nBEGIN IONS\nTITLE=scan=1\n200.0 2000.0\nEND IONS\n"
        mgf_file = tmp_path / "dup.mgf"
        mgf_file.write_text(mgf)
        idx = MgfSpectraIndex(mgf_file)
        assert len(idx) == 1
        spec = idx.get_spectrum(1)
        assert spec["mz_array"] == pytest.approx([100.0])  # first occurrence


# ---------------------------------------------------------------------------
# Integration test: real PXD054720 MGF file
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestMgfParserPXD054720:
    """Integration test with the real PXD054720 MGF.gz file."""

    @pytest.fixture(scope="class")
    def mgf_index(self) -> MgfSpectraIndex:
        if not PXD054720_MGF.exists():
            pytest.skip(f"Test data not available: {PXD054720_MGF}")
        return MgfSpectraIndex(PXD054720_MGF)

    def test_spectra_count(self, mgf_index: MgfSpectraIndex):
        """The file should contain a substantial number of spectra."""
        assert len(mgf_index) > 40000

    def test_first_scan_exists(self, mgf_index: MgfSpectraIndex):
        """Scan 154 is the first spectrum in the file (from SCANS=154)."""
        assert 154 in mgf_index

    def test_spectrum_has_peaks(self, mgf_index: MgfSpectraIndex):
        """First spectrum should have peaks."""
        spec = mgf_index.get_spectrum(154)
        assert spec is not None
        assert len(spec["mz_array"]) > 0
        assert len(spec["intensity_array"]) > 0
        assert len(spec["mz_array"]) == len(spec["intensity_array"])

    def test_mz_values_positive(self, mgf_index: MgfSpectraIndex):
        """All m/z values should be positive."""
        spec = mgf_index.get_spectrum(154)
        assert all(mz > 0 for mz in spec["mz_array"])

    def test_intensity_values_nonnegative(self, mgf_index: MgfSpectraIndex):
        """All intensity values should be non-negative."""
        spec = mgf_index.get_spectrum(154)
        assert all(i >= 0 for i in spec["intensity_array"])

    def test_random_scan_lookup(self, mgf_index: MgfSpectraIndex):
        """Spot-check a scan from the middle of the file."""
        scans = mgf_index.scans
        mid = scans[len(scans) // 2]
        spec = mgf_index.get_spectrum(mid)
        assert spec is not None
        assert len(spec["mz_array"]) == len(spec["intensity_array"])


# ---------------------------------------------------------------------------
# Unit tests: RTINSECONDS parsing
# ---------------------------------------------------------------------------


class TestMgfParserRtInseconds:
    """Test RTINSECONDS parsing in MGF."""

    def test_rtinseconds_parsed(self, tmp_path):
        mgf = tmp_path / "test.mgf"
        mgf.write_text("BEGIN IONS\nTITLE=spectrum1\nSCANS=1\nRTINSECONDS=456.78\nPEPMASS=500.0\n100.0 1000\nEND IONS\n")
        idx = MgfSpectraIndex(str(mgf))
        spec = idx.get_spectrum(1)
        assert spec["rt"] == pytest.approx(456.78)

    def test_rtinseconds_absent(self, tmp_path):
        mgf = tmp_path / "test.mgf"
        mgf.write_text("BEGIN IONS\nTITLE=spectrum2\nSCANS=2\nPEPMASS=500.0\n100.0 1000\nEND IONS\n")
        idx = MgfSpectraIndex(str(mgf))
        spec = idx.get_spectrum(2)
        assert spec.get("rt") is None

    def test_rtinseconds_distiller_format(self, tmp_path):
        """Distiller-style MGF with RTINSECONDS should be parsed."""
        mgf = tmp_path / "distiller_rt.mgf"
        mgf.write_text(DISTILLER_MGF)
        idx = MgfSpectraIndex(str(mgf))
        spec = idx.get_spectrum(154)
        assert spec is not None
        assert spec["rt"] == pytest.approx(1137.0618)
        # Spectrum without RTINSECONDS should have rt=None
        spec2 = idx.get_spectrum(353)
        assert spec2 is not None
        assert spec2.get("rt") is None


# ---------------------------------------------------------------------------
# Unit tests: position-based (0-based index) spectrum lookup
# ---------------------------------------------------------------------------


class TestMgfParserPositionIndex:
    """Test position-based (0-based index) spectrum lookup."""

    def test_get_spectrum_by_index(self, tmp_path):
        """Spectra should be accessible by 0-based position."""
        mgf = tmp_path / "test.mgf"
        mgf.write_text(
            "BEGIN IONS\nTITLE=first\nSCANS=100\n200.0 1000\nEND IONS\n"
            "BEGIN IONS\nTITLE=second\nSCANS=200\n300.0 2000\nEND IONS\n"
            "BEGIN IONS\nTITLE=third\nSCANS=300\n400.0 3000\nEND IONS\n"
        )
        idx = MgfSpectraIndex(str(mgf))

        # Position-based access
        spec0 = idx.get_spectrum_by_index(0)
        assert spec0 is not None
        assert spec0["mz_array"] == [200.0]

        spec1 = idx.get_spectrum_by_index(1)
        assert spec1 is not None
        assert spec1["mz_array"] == [300.0]

        spec2 = idx.get_spectrum_by_index(2)
        assert spec2 is not None
        assert spec2["mz_array"] == [400.0]

        # Out of bounds returns None
        assert idx.get_spectrum_by_index(3) is None
        assert idx.get_spectrum_by_index(-1) is None

    def test_position_index_same_object_as_scan_index(self, tmp_path):
        """Position and scan indices should reference the same dict object."""
        mgf = tmp_path / "test.mgf"
        mgf.write_text("BEGIN IONS\nTITLE=spec\nSCANS=42\n100.0 500\nEND IONS\n")
        idx = MgfSpectraIndex(str(mgf))
        assert idx._index[42] is idx._by_position[0]  # Same dict object
        assert idx.get_spectrum(42) == idx.get_spectrum_by_index(0)

    def test_position_tracks_unparseable_scans(self, tmp_path):
        """Position counter should advance even for spectra with no parseable scan."""
        mgf = tmp_path / "test.mgf"
        mgf.write_text(
            "BEGIN IONS\nTITLE=has_scan\nSCANS=10\n100.0 500\nEND IONS\n"
            "BEGIN IONS\nTITLE=no_parseable_scan_here\n200.0 600\nEND IONS\n"
            "BEGIN IONS\nTITLE=has_scan_too\nSCANS=20\n300.0 700\nEND IONS\n"
        )
        idx = MgfSpectraIndex(str(mgf))

        # Position 0 = first spectrum (scan=10)
        assert idx.get_spectrum_by_index(0)["mz_array"] == [100.0]
        # Position 1 = second spectrum (no parseable scan, but still indexed by position)
        assert idx.get_spectrum_by_index(1)["mz_array"] == [200.0]
        # Position 2 = third spectrum (scan=20)
        assert idx.get_spectrum_by_index(2)["mz_array"] == [300.0]

        # scan-based: second spectrum not accessible by scan
        assert idx.get_spectrum(10) is not None
        assert idx.get_spectrum(20) is not None
