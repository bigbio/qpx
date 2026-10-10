"""Tests for the DIA-NN m/z adapter (AMS_info -> mz view with lossless blobs)."""

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


def _write_ams_info(path, spectra):
    """Write an ``*_AMS_info.parquet`` file.

    ``spectra`` is a list of ``(run_name, scan, precursor_rt, mz, intensity)``.
    """
    table = pa.table(
        {
            "run_name": pa.array([s[0] for s in spectra], type=pa.string()),
            "scan": pa.array([[s[1]] for s in spectra], type=pa.list_(pa.int64())),
            "precursor_rt": pa.array([s[2] for s in spectra], type=pa.float64()),
            "mz_array": pa.array([np.asarray(s[3], np.float64) for s in spectra], type=pa.list_(pa.float64())),
            "intensity_array": pa.array([np.asarray(s[4], np.float32) for s in spectra], type=pa.list_(pa.float32())),
        }
    )
    pq.write_table(table, str(path))


def _one_spectrum():
    return (np.array([100.0, 200.0, 300.0], np.float64), np.array([10.0, 20.0, 30.0], np.float32))


def test_adapter_dedups_by_run_and_scan(tmp_path):
    ams_dir = tmp_path / "ams"
    ams_dir.mkdir()
    mz, it = _one_spectrum()
    spectra = [
        ("run_a", 1, 60.0, mz, it),
        ("run_a", 1, 60.0, mz, it),  # duplicate (same run+scan)
        ("run_a", 2, 61.0, mz, it),
    ]
    _write_ams_info(ams_dir / "run_a_AMS_info.parquet", spectra)

    out = tmp_path / "out.mz.parquet"
    written = DiannMzAdapter(compression="zstd").convert(str(ams_dir), str(out))
    assert written == 2  # the duplicate collapsed

    table = pq.read_table(str(out))
    assert table.num_rows == 2
    assert sorted(table.column("scan").to_pylist()) == [1, 2]


def test_adapter_stores_encoded_blobs_in_mz_intensity(tmp_path):
    ams_dir = tmp_path / "ams"
    ams_dir.mkdir()
    mz, it = _one_spectrum()
    _write_ams_info(ams_dir / "run_a_AMS_info.parquet", [("run_a", 1, 60.0, mz, it)])

    out = tmp_path / "out.mz.parquet"
    DiannMzAdapter(compression="zstd").convert(str(ams_dir), str(out))

    table = pq.read_table(str(out))
    # mz / intensity are binary (losslessly encoded), no separate blob columns
    assert pa.types.is_binary(table.schema.field("mz").type)
    assert pa.types.is_binary(table.schema.field("intensity").type)
    assert "mz_blob" not in table.column_names
    assert "intensity_blob" not in table.column_names
    assert all(b is not None and len(b) > 0 for b in table.column("mz").to_pylist())


def test_adapter_roundtrip_is_bit_exact(tmp_path):
    ams_dir = tmp_path / "ams"
    ams_dir.mkdir()
    mz = np.array([100.0, 200.0, 300.0, 400.0], np.float64)
    it = np.array([10.0, 20.0, 30.0, 40.0], np.float32)
    _write_ams_info(ams_dir / "run_a_AMS_info.parquet", [("run_a", 1, 60.0, mz, it)])

    out = tmp_path / "out.mz.parquet"
    DiannMzAdapter(compression="zstd").convert(str(ams_dir), str(out))

    spec = MzSpectra.from_file(out)
    df = spec.to_df()
    dmz, dit = spec.decode_spectrum(df.iloc[0])
    np.testing.assert_array_equal(dmz, mz.astype(np.float32))
    np.testing.assert_array_equal(dit, it)


def test_adapter_matched_scans_filter(tmp_path):
    ams_dir = tmp_path / "ams"
    ams_dir.mkdir()
    mz, it = _one_spectrum()
    spectra = [
        ("run_a", 1, 60.0, mz, it),
        ("run_a", 2, 61.0, mz, it),
        ("run_a", 3, 62.0, mz, it),
    ]
    _write_ams_info(ams_dir / "run_a_AMS_info.parquet", spectra)

    out = tmp_path / "out.mz.parquet"
    matched = {("run_a", 1), ("run_a", 3)}
    written = DiannMzAdapter(compression="zstd").convert(str(ams_dir), str(out), matched_scans=matched)
    assert written == 2

    table = pq.read_table(str(out))
    assert sorted(table.column("scan").to_pylist()) == [1, 3]


def test_adapter_empty_folder_raises(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError):
        DiannMzAdapter().convert(str(empty), str(tmp_path / "out.mz.parquet"))
