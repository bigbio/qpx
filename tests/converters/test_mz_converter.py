"""Tests for mzML -> mz.parquet conversion (full spectra).

Covers ``SpectraMappingTransform.write_mz_parquet_from_dir``, the directory
entry point used by ``qpxc convert mz``. Mini mzML files are generated with
pyopenms so the tests are self-contained and fast.
"""

import gzip
import shutil

import pyarrow.parquet as pq
import pytest

pytest.importorskip("pyopenms")


def _make_mini_mzml(path, specs):
    """Write a minimal mzML at *path*.

    *specs* is a list of ``(ms_level, scan, rt_seconds)`` tuples. MS2+ spectra
    get a precursor. Native IDs use the Thermo ``scan=N`` form so scan-number
    extraction can be exercised.
    """
    import pyopenms as oms

    exp = oms.MSExperiment()
    for ms_level, scan, rt in specs:
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(ms_level)
        spectrum.setRT(rt)
        spectrum.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
        if ms_level >= 2:
            precursor = oms.Precursor()
            precursor.setMZ(500.0 + scan)
            precursor.setCharge(2)
            spectrum.setPrecursors([precursor])
        spectrum.set_peaks(([100.0, 200.0, 300.0], [10.0, 20.0, 30.0]))
        exp.addSpectrum(spectrum)
    oms.MzMLFile().store(str(path), exp)


@pytest.fixture
def mzml_dir(tmp_path):
    """Directory with one mini mzML: 1 MS1 + 2 MS2 spectra."""
    d = tmp_path / "mzml"
    d.mkdir()
    _make_mini_mzml(d / "run_a.mzML", [(1, 1, 60.0), (2, 2, 60.5), (2, 3, 61.0)])
    return d


def test_from_dir_all_levels(mzml_dir, tmp_path):
    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    out = tmp_path / "out.mz.parquet"
    SpectraMappingTransform(mzml_directory=mzml_dir).write_mz_parquet_from_dir(out)

    table = pq.read_table(str(out))
    assert table.num_rows == 3  # 1 MS1 + 2 MS2
    assert {"run_file_name", "scan", "ms_level", "mz", "intensity"} <= set(table.column_names)
    assert set(table.column("run_file_name").to_pylist()) == {"run_a"}
    assert sorted(table.column("scan").to_pylist()) == [1, 2, 3]
    # every spectrum keeps its peaks (decoded from the lossless binary blobs)
    from qpx.core.spectra_codec import decode_mz

    mz_blobs = table.column("mz").to_pylist()
    int_blobs = table.column("intensity").to_pylist()
    assert all(len(decode_mz(mz_blobs[i], len(int_blobs[i]) // 4)) == 3 for i in range(table.num_rows))


def test_from_dir_ms2_only(mzml_dir, tmp_path):
    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    out = tmp_path / "out_ms2.mz.parquet"
    SpectraMappingTransform(mzml_directory=mzml_dir).write_mz_parquet_from_dir(out, ms_levels=[2])

    table = pq.read_table(str(out))
    assert table.num_rows == 2
    assert set(table.column("ms_level").to_pylist()) == {2}
    precursors = table.column("precursors").to_pylist()
    assert all(p and p[0]["selected_ion_mz"] > 0 for p in precursors)


def test_from_dir_reads_gzip(tmp_path):
    """pyopenms reads .mzML.gz directly; run name strips the .gz extension."""
    d = tmp_path / "mzml_gz"
    d.mkdir()
    plain = d / "run_b.mzML"
    _make_mini_mzml(plain, [(2, 10, 60.0)])
    with open(plain, "rb") as fin, gzip.open(str(d / "run_b.mzML.gz"), "wb") as fout:
        shutil.copyfileobj(fin, fout)
    plain.unlink()  # leave only the .gz

    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    out = tmp_path / "out_gz.mz.parquet"
    SpectraMappingTransform(mzml_directory=d).write_mz_parquet_from_dir(out)

    table = pq.read_table(str(out))
    assert table.num_rows == 1
    assert table.column("run_file_name").to_pylist() == ["run_b"]
    assert table.column("scan").to_pylist() == [10]


def test_from_dir_empty_raises(tmp_path):
    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError):
        SpectraMappingTransform(mzml_directory=empty).write_mz_parquet_from_dir(tmp_path / "x.mz.parquet")
