"""Tests for the pdc2qpx orchestration (download + CDAP + mz).

Downloads are skipped (``skip_download=True``); the CDAP ``.psm`` fixture is
staged into a ``<download_dir>/<study>/`` layout and, for the spectra case, a
mini mzML is generated with pyopenms. This exercises the orchestration wiring
without hitting the network.
"""

import shutil
from pathlib import Path

import pyarrow.parquet as pq
import pytest

_FIXTURE_CDAP = Path(__file__).resolve().parent.parent / "examples" / "cdap"


def _stage_study(tmp_path, study, with_mzml=False):
    """Build a ``download_dir/<study>/`` with fixture .psm (+ optional mini mzML)."""
    download_dir = tmp_path / "downloads"
    study_dir = download_dir / study
    study_dir.mkdir(parents=True)
    for psm in _FIXTURE_CDAP.glob("*.psm"):
        shutil.copy(psm, study_dir / psm.name)

    if with_mzml:
        pytest.importorskip("pyopenms")
        import pyopenms as oms

        exp = oms.MSExperiment()
        for ms_level, scan in [(1, 1), (2, 2)]:
            spectrum = oms.MSSpectrum()
            spectrum.setMSLevel(ms_level)
            spectrum.setRT(60.0 + scan)
            spectrum.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
            if ms_level >= 2:
                precursor = oms.Precursor()
                precursor.setMZ(500.0)
                precursor.setCharge(2)
                spectrum.setPrecursors([precursor])
            spectrum.set_peaks(([100.0, 200.0], [10.0, 20.0]))
            exp.addSpectrum(spectrum)
        oms.MzMLFile().store(str(study_dir / "run_x.mzML"), exp)

    return download_dir


def test_pdc2qpx_base(tmp_path):
    from qpx.pipeline.pdc2qpx import run_pdc2qpx

    study = "PDC_TEST"
    download_dir = _stage_study(tmp_path, study)
    out = tmp_path / "qpx"

    result = run_pdc2qpx(study, download_dir, out, skip_download=True, include_metadata=False)

    assert "mz" not in result
    feature = out / f"{study}.feature.parquet"
    assert feature.exists()
    assert pq.read_table(str(feature)).num_rows > 0


def test_pdc2qpx_stores_bare_accessions(tmp_path):
    """The fixture lists sp|Q29443|TRFE_BOVIN, which qpxc convert cdap stores as Q29443."""
    from qpx.pipeline.pdc2qpx import run_pdc2qpx

    study = "PDC_TEST"
    out = tmp_path / "qpx"
    run_pdc2qpx(study, _stage_study(tmp_path, study), out, skip_download=True, include_metadata=False)

    rows = pq.read_table(str(out / f"{study}.psm.parquet"), columns=["protein_accessions"]).column(0).to_pylist()
    accessions = {accession for row in rows for accession in row or []}
    assert "Q29443" in accessions
    assert "sp|Q29443|TRFE_BOVIN" not in accessions


def test_pdc2qpx_with_spectra(tmp_path):
    from qpx.pipeline.pdc2qpx import run_pdc2qpx

    study = "PDC_TEST"
    download_dir = _stage_study(tmp_path, study, with_mzml=True)
    out = tmp_path / "qpx"

    result = run_pdc2qpx(study, download_dir, out, include_spectra=True, skip_download=True, include_metadata=False)

    assert result["mz"] == out / f"{study}.mz.parquet"
    mz = pq.read_table(str(result["mz"]))
    assert mz.num_rows == 2  # 1 MS1 + 1 MS2
    assert {"run_file_name", "scan"} <= set(mz.column_names)
    # The ontology view, written with the CDAP views, also describes the mz view.
    ontology = pq.read_table(out / f"{study}.ontology.parquet").to_pylist()
    terms = {row["field_name"]: row["ontology_accession"] for row in ontology if row["view"] == "mz"}
    assert terms["ms_level"] == "MS:1000511"
    assert {row["view"] for row in ontology} >= {"psm", "mz"}


def test_pdc2qpx_describes_the_metadata_views(tmp_path, monkeypatch):
    """The sample and run views built from PDC metadata come after the CDAP ontology view and get its rows."""
    from qpx.pipeline.pdc2qpx import run_pdc2qpx
    from qpx.writers import RunWriter
    from tests.conftest import make_run_record

    def build(_study, output_folder, prefix, **_options):
        path = output_folder / f"{prefix}.run.parquet"
        with RunWriter(path) as writer:
            writer.write_batch([make_run_record()])
        return {"run": path}

    monkeypatch.setattr("qpx.pipeline.pdc2qpx._fetch_run_to_plex", lambda study: None)
    monkeypatch.setattr("qpx.converters.cdap.pdc_sample_run.build_sample_run_from_pdc", build)
    study = "PDC_TEST"
    out = tmp_path / "qpx"
    run_pdc2qpx(study, _stage_study(tmp_path, study), out, skip_download=True)

    views = {row["view"] for row in pq.read_table(out / f"{study}.ontology.parquet").to_pylist()}
    assert {"psm", "run"} <= views


def test_pdc2qpx_missing_psm_raises(tmp_path):
    from qpx.pipeline.pdc2qpx import run_pdc2qpx

    (tmp_path / "downloads" / "EMPTY").mkdir(parents=True)
    with pytest.raises(FileNotFoundError):
        run_pdc2qpx("EMPTY", tmp_path / "downloads", tmp_path / "qpx", skip_download=True)


def test_pdc2qpx_non_cdap_psm_raises(tmp_path):
    """A .psm whose header lacks CDAP signature columns is rejected early."""
    from qpx.pipeline.pdc2qpx import run_pdc2qpx

    study = "PDC_NONCDAP"
    study_dir = tmp_path / "downloads" / study
    study_dir.mkdir(parents=True)
    # A non-CDAP .psm: header without FileName/ScanNum/PeptideSequence/Protein.
    (study_dir / "bad.psm").write_text("colA\tcolB\tcolC\n1\t2\t3\n", encoding="utf-8")

    with pytest.raises(ValueError, match="not a CDAP-format"):
        run_pdc2qpx(study, tmp_path / "downloads", tmp_path / "qpx", skip_download=True)


def test_assert_cdap_psm_accepts_crlf(tmp_path):
    """A valid CDAP .psm with Windows CRLF line endings passes pre-check."""
    from qpx.pipeline.pdc2qpx import _assert_cdap_psm

    study_dir = tmp_path / "study"
    study_dir.mkdir()
    header = "FileName\tScanNum\tPeptideSequence\tProtein\tCharge\r\n"
    row = "run1.raw\t100\tPEPTIDE\tNP_000001.1\t2\r\n"
    (study_dir / "valid_crlf.psm").write_bytes((header + row).encode("utf-8"))

    _assert_cdap_psm(study_dir)


def _stage_psm_only(download_dir, study):
    """Stage download_dir/<study>/ with the fixture .psm files."""
    study_dir = download_dir / study
    study_dir.mkdir(parents=True)
    for psm in _FIXTURE_CDAP.glob("*.psm"):
        shutil.copy(psm, study_dir / psm.name)


def test_run_pdc2qpx_batch_isolates_failures(tmp_path):
    """One bad study (no .psm) is recorded as failed; the batch still finishes."""
    from qpx.pipeline.pdc2qpx import run_pdc2qpx_batch

    download_dir = tmp_path / "downloads"
    _stage_psm_only(download_dir, "PDC_GOOD")
    (download_dir / "PDC_BAD").mkdir(parents=True)  # empty -> no .psm
    out_root = tmp_path / "qpx"

    results = run_pdc2qpx_batch(["PDC_GOOD", "PDC_BAD"], download_dir, out_root, skip_download=True, include_metadata=False)

    assert results["PDC_GOOD"]["status"] == "ok"
    assert results["PDC_BAD"]["status"] == "failed"
    assert results["PDC_BAD"]["error"]
    assert (out_root / "PDC_GOOD" / "PDC_GOOD.feature.parquet").exists()
    assert sum(1 for r in results.values() if r["status"] == "ok") == 1


def test_run_pdc2qpx_batch_stop_on_error(tmp_path):
    """With continue_on_error=False the first failure propagates."""
    from qpx.pipeline.pdc2qpx import run_pdc2qpx_batch

    download_dir = tmp_path / "downloads"
    (download_dir / "PDC_BAD").mkdir(parents=True)
    with pytest.raises(FileNotFoundError):
        run_pdc2qpx_batch(
            ["PDC_BAD"],
            download_dir,
            tmp_path / "qpx",
            skip_download=True,
            include_metadata=False,
            continue_on_error=False,
        )


def test_pdc2qpx_cli_comma_accessions_runs_batch(tmp_path):
    """`-a PDC_A,PDC_B` resolves to multiple studies and writes one folder each."""
    from click.testing import CliRunner

    from qpx.cli.pdc2qpx import pdc2qpx_cmd

    download_dir = tmp_path / "downloads"
    _stage_psm_only(download_dir, "PDC_A")
    _stage_psm_only(download_dir, "PDC_B")
    out = tmp_path / "qpx"

    result = CliRunner().invoke(
        pdc2qpx_cmd,
        [
            "-a",
            "PDC_A,PDC_B",
            "--download-dir",
            str(download_dir),
            "--output-folder",
            str(out),
            "--skip-download",
            "--no-metadata",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "2/2 studies ok" in result.output
    assert (out / "PDC_A" / "PDC_A.feature.parquet").exists()
    assert (out / "PDC_B" / "PDC_B.feature.parquet").exists()


def test_parse_accessions(tmp_path):
    from qpx.pipeline.pdc2qpx import parse_accessions

    # single / comma+whitespace / dedup / ignore comments
    assert parse_accessions("PDC000109") == ["PDC000109"]
    assert parse_accessions("PDC1, PDC2 PDC1 #c") == ["PDC1", "PDC2"]
    # CSV with a pdc_study_id column
    csvf = tmp_path / "list.csv"
    csvf.write_text("pdc_study_id,quant_method\nPDC5,TMT10\nPDC6,Label Free\n", encoding="utf-8")
    assert parse_accessions(str(csvf)) == ["PDC5", "PDC6"]


def test_pdc2qpx_cli_does_not_require_pridepy_when_skipping(tmp_path, monkeypatch):
    """The parse + skip-download + no-metadata path must work without the pridepy extra.

    Regression guard: pridepy is an optional extra (absent in CI), so accession
    parsing and the offline conversion path must not import it.
    """
    import builtins

    real_import = builtins.__import__

    def _no_pridepy(name, *args, **kwargs):
        if name == "pridepy" or name.startswith("pridepy."):
            raise ImportError("pridepy blocked for this test")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _no_pridepy)

    from click.testing import CliRunner

    from qpx.cli.pdc2qpx import pdc2qpx_cmd

    download_dir = tmp_path / "downloads"
    _stage_psm_only(download_dir, "PDC_X")
    out = tmp_path / "qpx"

    result = CliRunner().invoke(
        pdc2qpx_cmd,
        ["-a", "PDC_X", "--download-dir", str(download_dir), "--output-folder", str(out), "--skip-download", "--no-metadata"],
    )

    assert result.exit_code == 0, result.output
    assert (out / "PDC_X.feature.parquet").exists()
