"""Tests for the DIA-NN m/z adapter (mzML-statistics MS2 -> mz view with float32 arrays)."""

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

try:
    from qpx.converters.diann.mz_adapter import DiannMzAdapter
    from qpx.core.data import MzSpectra
except ImportError:  # pragma: no cover - requires pyopenms/sdrf_pipelines
    DiannMzAdapter = None
    MzSpectra = None

pytestmark = pytest.mark.skipif(
    DiannMzAdapter is None,
    reason="qpx converter dependencies (pyopenms / sdrf_pipelines) not installed",
)


def _write_ms2_info(path, spectra):
    """Write an ``*_ms2_info.parquet`` file.

    ``spectra`` is a list of ``(scan, ms_level, mz, intensity)``.
    """
    table = pa.table(
        {
            "scan": pa.array([str(s[0]) for s in spectra], type=pa.string()),
            "ms_level": pa.array([s[1] for s in spectra], type=pa.int32()),
            "mz_array": pa.array([np.asarray(s[2], np.float64) for s in spectra], type=pa.list_(pa.float64())),
            "intensity_array": pa.array([np.asarray(s[3], np.float64) for s in spectra], type=pa.list_(pa.float64())),
        }
    )
    pq.write_table(table, str(path))


def _write_ms_info(path, rows):
    """Write an ``*_ms_info.parquet`` file. ``rows`` is a list of ``(scan, rt)``; precursor columns are null."""
    table = pa.table(
        {
            "scan": pa.array([str(r[0]) for r in rows], type=pa.string()),
            "rt": pa.array([r[1] for r in rows], type=pa.float64()),
            "precursor_charge": pa.nulls(len(rows), type=pa.int32()),
            "precursor_mz": pa.nulls(len(rows), type=pa.float64()),
        }
    )
    pq.write_table(table, str(path))


def _one_spectrum():
    return (np.array([100.0, 200.0, 300.0], np.float64), np.array([10.0, 20.0, 30.0], np.float64))


def test_adapter_dedups_by_run_and_scan(tmp_path):
    stats_dir = tmp_path / "mzml_statistics"
    stats_dir.mkdir()
    mz, it = _one_spectrum()
    spectra = [
        (1, 2, mz, it),
        (1, 2, mz, it),  # duplicate scan
        (2, 2, mz, it),
    ]
    _write_ms2_info(stats_dir / "run_a_ms2_info.parquet", spectra)

    out = tmp_path / "out.mz.parquet"
    written = DiannMzAdapter().convert(str(stats_dir), str(out))
    assert written == 2  # duplicate scan collapsed

    table = pq.read_table(str(out))
    assert table.num_rows == 2
    assert sorted(table.column("scan").to_pylist()) == [[1], [2]]
    assert table.column("run_file_name").to_pylist() == ["run_a", "run_a"]


def test_adapter_stores_float32_peak_arrays(tmp_path):
    stats_dir = tmp_path / "mzml_statistics"
    stats_dir.mkdir()
    mz, it = _one_spectrum()
    _write_ms2_info(stats_dir / "run_a_ms2_info.parquet", [(1, 2, mz, it)])

    out = tmp_path / "out.mz.parquet"
    DiannMzAdapter().convert(str(stats_dir), str(out))

    table = pq.read_table(str(out))
    assert pa.types.is_list(table.schema.field("mz").type)
    assert pa.types.is_floating(table.schema.field("mz").type.value_type)
    assert pa.types.is_list(table.schema.field("intensity").type)
    assert pa.types.is_floating(table.schema.field("intensity").type.value_type)
    assert all(len(v) > 0 for v in table.column("mz").to_pylist())


def test_adapter_roundtrip_preserves_peaks(tmp_path):
    stats_dir = tmp_path / "mzml_statistics"
    stats_dir.mkdir()
    mz = np.array([300.0, 100.0, 200.0, 400.0], np.float64)
    it = np.array([30.0, 10.0, 20.0, 40.0], np.float64)
    _write_ms2_info(stats_dir / "run_a_ms2_info.parquet", [(1, 2, mz, it)])

    out = tmp_path / "out.mz.parquet"
    DiannMzAdapter().convert(str(stats_dir), str(out))

    table = pq.read_table(str(out))
    mz_out = np.asarray(table.column("mz").to_pylist()[0])
    it_out = np.asarray(table.column("intensity").to_pylist()[0])
    np.testing.assert_array_equal(mz_out, mz.astype(np.float32))
    np.testing.assert_array_equal(it_out, it.astype(np.float32))


def test_adapter_matched_scans_filter(tmp_path):
    stats_dir = tmp_path / "mzml_statistics"
    stats_dir.mkdir()
    mz, it = _one_spectrum()
    spectra = [(1, 2, mz, it), (2, 2, mz, it), (3, 2, mz, it)]
    _write_ms2_info(stats_dir / "run_a_ms2_info.parquet", spectra)

    out = tmp_path / "out.mz.parquet"
    matched = {("run_a", 1), ("run_a", 3)}
    written = DiannMzAdapter().convert(str(stats_dir), str(out), matched_scans=matched)
    assert written == 2

    table = pq.read_table(str(out))
    assert sorted(table.column("scan").to_pylist()) == [[1], [3]]


def test_adapter_fills_scan_start_time_from_ms_info(tmp_path):
    stats_dir = tmp_path / "mzml_statistics"
    stats_dir.mkdir()
    mz, it = _one_spectrum()
    _write_ms2_info(stats_dir / "run_a_ms2_info.parquet", [(1, 2, mz, it)])
    _write_ms_info(stats_dir / "run_a_ms_info.parquet", [(1, 301.2)])  # 301.2 s -> 5.02 min

    out = tmp_path / "out.mz.parquet"
    DiannMzAdapter().convert(str(stats_dir), str(out))

    table = pq.read_table(str(out))
    rt = table.column("scan_start_time").to_pylist()[0]
    assert rt is not None
    assert abs(rt - 301.2 / 60.0) < 1e-4


def test_adapter_empty_folder_raises(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError):
        DiannMzAdapter().convert(str(empty), str(tmp_path / "out.mz.parquet"))
