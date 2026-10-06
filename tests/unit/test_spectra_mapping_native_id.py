"""Spectral annotation resolves complete native IDs without guessing components."""

import pytest

from qpx import Dataset
from qpx.transforms.spectra_mapping import SpectraMappingTransform
from qpx.writers.psm import PsmWriter
from tests.conftest import make_psm_record


def _write_spectra(path, references):
    """Write distinct peaks for each supplied native ID."""
    oms = pytest.importorskip("pyopenms")
    experiment = oms.MSExperiment()
    for position, reference in enumerate(references):
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(2)
        spectrum.setRT(float(position))
        spectrum.setNativeID(reference)
        spectrum.set_peaks(([100.0 + position], [10.0 + position]))
        experiment.addSpectrum(spectrum)
    oms.MzMLFile().store(str(path), experiment)


def test_annotation_matches_whole_native_id_and_preserves_single_scan(tmp_path):
    """Equal scan ordinals in different frames receive their own spectral peaks."""
    _write_spectra(
        tmp_path / "run_01.mzML",
        ["frame=120 scan=475", "frame=121 scan=475", "controllerType=0 controllerNumber=1 scan=120"],
    )
    scans = [[120, 475], [121, 475], [120], [120, 999], []]
    with PsmWriter(tmp_path / "native.psm.parquet") as writer:
        writer.write_batch([make_psm_record(scan=scan) for scan in scans])

    with Dataset(tmp_path, structures=["psm"]) as dataset, SpectraMappingTransform(tmp_path) as transform:
        annotated = transform.annotate_dataset_psms(dataset)

    peaks = {tuple(row.scan): row.mz_array for row in annotated.itertuples()}
    assert peaks == {(120, 475): [100.0], (121, 475): [101.0], (120,): [102.0], (120, 999): [], (): []}


def test_ambiguous_components_do_not_attach_a_spectrum(tmp_path, caplog):
    """Different native IDs encoding the same tuple cannot be uniquely mapped."""
    _write_spectra(tmp_path / "run_01.mzML", ["frame=120 scan=475", "sample=120 cycle=475"])
    with SpectraMappingTransform(tmp_path) as transform:
        count, mz, intensity = transform.get_spectrum("run_01", (120, 475))

    assert count == 0
    assert mz.size == intensity.size == 0
    assert "No unique spectrum" in caplog.text
