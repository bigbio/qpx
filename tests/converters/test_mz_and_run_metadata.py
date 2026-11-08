"""mzML metadata reaches the mz view, and SDRF metadata reaches the run it belongs to."""

import gzip
import json

import duckdb
import pandas as pd
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.converters.sdrf import SdrfConverter
from qpx.core.data.mz import MzSpectra
from qpx.transforms.spectra_mapping import SpectraMappingTransform


@pytest.fixture
def mzml_dir(tmp_path):
    oms = pytest.importorskip("pyopenms")
    experiment = oms.MSExperiment()
    ms1 = oms.MSSpectrum()
    ms1.setMSLevel(1)
    ms1.setRT(30.0)
    ms1.setNativeID("scan=1")
    ms1.set_peaks(([100.0], [1.0]))
    ms2 = oms.MSSpectrum()
    ms2.setMSLevel(2)
    ms2.setRT(60.0)
    ms2.setNativeID("scan=2")
    precursor = oms.Precursor()
    precursor.setMZ(500.5)
    precursor.setMetaValue("isolation window target m/z", 500.0)
    precursor.setIsolationWindowLowerOffset(12.5)
    precursor.setIsolationWindowUpperOffset(0.0)
    precursor.setActivationEnergy(30.0)
    precursor.setActivationMethods({oms.Precursor.ActivationMethod.HCID})
    ms2.setPrecursors([precursor])
    ms2.setMetaValue("filter string", "FTMS + c NSI d Full ms2 500.00@hcd30.00 [100-2000]")
    ms2.set_peaks(([100.0, 200.0], [1.0, 2.0]))
    experiment.addSpectrum(ms1)
    experiment.addSpectrum(ms2)
    date = oms.DateTime()
    date.set("2024-01-02 03:04:05")
    experiment.setDateTime(date)
    folder = tmp_path / "mzml"
    folder.mkdir()
    oms.MzMLFile().store(str(folder / "run.mzML"), experiment)
    return folder


def _mz_rows(mzml_dir, output, **kwargs):
    SpectraMappingTransform(str(mzml_dir)).write_mz_parquet_from_dir(output, **kwargs)
    return {row["ms_level"]: row for row in pq.read_table(output).to_pylist()}


def test_isolation_window_keeps_its_target_and_a_zero_offset(mzml_dir, tmp_path):
    (precursor,) = _mz_rows(mzml_dir, tmp_path / "x.mz.parquet")[2]["precursors"]
    assert precursor["selected_ion_mz"] == pytest.approx(500.5)
    assert precursor["isolation_window_target"] == pytest.approx(500.0)
    assert precursor["isolation_window_lower"] == pytest.approx(12.5)
    assert precursor["isolation_window_upper"] == 0.0


def test_spectrum_records_activation_energy_filter_and_analyzer(mzml_dir, tmp_path):
    params = {p["cv_name"]: p["cv_value"] for p in _mz_rows(mzml_dir, tmp_path / "x.mz.parquet")[2]["cv_params"]}
    assert params["activation_method"] == "High-energy collision-induced dissociation"
    assert float(params["collision_energy"]) == 30.0
    assert params["analyzer"] == "FTMS"
    assert params["filter_string"].startswith("FTMS")


def test_footer_states_the_content_and_the_acquisition_time(mzml_dir, tmp_path):
    output = tmp_path / "x.mz.parquet"
    _mz_rows(mzml_dir, output)
    footer = pq.read_schema(output).metadata
    written = pq.read_metadata(output).metadata
    assert (footer | written)[b"spectrum_content"] == b"all_spectra"
    assert json.loads((footer | written)[b"acquisition_start"]) == {"run": "2024-01-02 03:04:05"}
    partial = tmp_path / "y.mz.parquet"
    _mz_rows(mzml_dir, partial, ms_levels=[2])
    assert pq.read_metadata(partial).metadata[b"spectrum_content"] == b"selected_ms_levels"


def test_run_summary_counts_spectra_per_level(mzml_dir, tmp_path):
    output = tmp_path / "x.mz.parquet"
    _mz_rows(mzml_dir, output)
    (row,) = MzSpectra.from_file(output).run_summary().to_df().to_dict("records")
    assert (row["ms1_spectra"], row["ms2_spectra"]) == (1, 1)
    assert row["last_scan_time"] == pytest.approx(1.0)


def _sdrf(path, rows):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def _base_row(source, data_file, label, **extra):
    return {
        "source name": source,
        "characteristics[organism]": "Homo sapiens",
        "characteristics[organism part]": "liver",
        "characteristics[biological replicate]": "1",
        "comment[data file]": data_file,
        "comment[label]": label,
        "comment[technical replicate]": "1",
        **extra,
    }


def test_each_run_keeps_its_own_modifications_and_label_definitions(tmp_path):
    heavy = "NT=Label:13C(6)15N(2);AC=UNIMOD:259;TA=K;MT=fixed"
    phospho = "NT=Phospho;AC=UNIMOD:21;TA=S;MT=variable"
    oxidation = "NT=Oxidation;AC=UNIMOD:35;TA=M;MT=variable"
    path = _sdrf(
        tmp_path / "in.sdrf.tsv",
        [
            _base_row("S1", "a.raw", "SILAC heavy", **{"comment[modification parameters]": heavy}),
            _base_row("S2", "a.raw", "SILAC light", **{"comment[modification parameters]": oxidation}),
            _base_row("S3", "b.raw", "SILAC light", **{"comment[modification parameters]": phospho}),
        ],
    )
    output = tmp_path / "x.run.parquet"
    SdrfConverter().convert(str(path), run_output=str(output))
    runs = {row["run_file_name"]: row for row in pq.read_table(output).to_pylist()}
    assert {mod["name"] for mod in runs["a"]["modification_parameters"]} == {"Label:13C(6)15N(2)", "Oxidation"}
    assert {mod["name"] for mod in runs["b"]["modification_parameters"]} == {"Phospho"}
    (definition,) = runs["a"]["label_definitions"]
    assert definition["label"] == "SILAC heavy"
    assert definition["modifications"][0]["target_amino_acid"] == "K"
    assert definition["modifications"][0]["mass_delta"] == pytest.approx(8.014199, abs=1e-3)
    assert runs["b"]["label_definitions"] is None


def test_channel_role_mixture_and_gradient_length(tmp_path):
    path = _sdrf(
        tmp_path / "in.sdrf.tsv",
        [
            _base_row(
                "bridge",
                "a.raw",
                "TMT126",
                **{"comment[reference channel]": "TMT126", "comment[mixture]": "plex1", "comment[gradient length]": "2 h"},
            ),
            _base_row(
                "S1",
                "a.raw",
                "TMT127N",
                **{"comment[reference channel]": "TMT126", "comment[mixture]": "plex1", "comment[gradient length]": "2 h"},
            ),
        ],
    )
    output = tmp_path / "x.run.parquet"
    SdrfConverter().convert(str(path), run_output=str(output))
    (run,) = pq.read_table(output).to_pylist()
    assert run["gradient_length"] == 7200.0
    roles = {sample["sample_accession"]: sample for sample in run["samples"]}
    assert roles["bridge"]["channel_role"] == "reference"
    assert roles["S1"]["channel_role"] == "sample"
    assert roles["S1"]["mixture"] == "plex1"
    assert duckdb.connect().execute(f"SELECT COUNT(*) FROM '{output}'").fetchall() == [(1,)]


@pytest.mark.parametrize(
    "values, expected",
    [
        (["2011-05-31"], "2011-05-31"),
        (["2014-07-23T00:14:06"], "2014-07-23T00:14:06"),
        (["not available", "2016-10-06T22:36:56.803000"], "2016-10-06T22:36:56.803000"),
        (["not applicable", "2015-02-19T20:08:48Z"], "2015-02-19T20:08:48Z"),
        (["2011-05-31", "2011-06-01"], "2011-05-31"),
        (["not available", "not applicable"], None),
    ],
)
def test_run_acquisition_start_keeps_first_valid_sdrf_text(tmp_path, values, expected):
    """SDRF acquisition dates retain their precision and ignore missing-value sentinels."""
    rows = [
        _base_row(f"S{i}", "a.raw", label, **{"comment[acquisition date]": value})
        for i, (label, value) in enumerate(zip(("TMT126", "TMT127N"), values))
    ]
    path = _sdrf(tmp_path / "in.sdrf.tsv", rows)
    output = tmp_path / "x.run.parquet"
    with SdrfConverter() as converter:
        converter.convert(str(path), run_output=str(output))
    assert pq.read_table(output)["acquisition_start"].to_pylist() == [expected]


def test_run_acquisition_start_without_sources_is_null(tmp_path):
    """A run without an SDRF date or mzML directory has no acquisition start."""
    path = _sdrf(tmp_path / "in.sdrf.tsv", [_base_row("S1", "a.raw", "LFQ")])
    output = tmp_path / "x.run.parquet"
    with SdrfConverter() as converter:
        converter.convert(str(path), run_output=str(output))
    assert pq.read_table(output)["acquisition_start"].to_pylist() == [None]


def _mzml_header(path, acquisition_start):
    """Write only an mzML header; reading beyond the run start would fail on the incomplete XML."""
    attribute = f' startTimeStamp="{acquisition_start}"' if acquisition_start is not None else ""
    header = f'<mzML xmlns="http://psi.hupo.org/ms/mzml"><run id="run"{attribute}>'
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as stream:
            stream.write(header)
    else:
        path.write_text(header)


@pytest.mark.parametrize(
    "extension, timestamp",
    [(".mzML", "2016-10-06T22:36:56.803000"), (".mzML.gz", "2015-02-19T20:08:48Z")],
)
def test_run_acquisition_start_reads_only_mzml_header(tmp_path, extension, timestamp):
    """Plain and gzipped mzML headers supply exact timestamps without reading spectra."""
    _mzml_header(tmp_path / f"a{extension}", timestamp)
    path = _sdrf(tmp_path / "in.sdrf.tsv", [_base_row("S1", "a.raw", "LFQ")])
    output = tmp_path / "x.run.parquet"
    with SdrfConverter() as converter:
        converter.convert(str(path), run_output=str(output), mzml_dir=tmp_path)
    assert pq.read_table(output)["acquisition_start"].to_pylist() == [timestamp]


@pytest.mark.parametrize("sdrf_date, expected", [("2011-05-31", "2011-05-31"), ("not available", "2015-02-19T20:08:48Z")])
def test_run_acquisition_start_prefers_sdrf_over_mzml(tmp_path, sdrf_date, expected):
    """A valid SDRF date takes precedence; a sentinel allows the mzML source."""
    _mzml_header(tmp_path / "a.mzML", "2015-02-19T20:08:48Z")
    path = _sdrf(tmp_path / "in.sdrf.tsv", [_base_row("S1", "a.raw", "LFQ", **{"comment[acquisition date]": sdrf_date})])
    output = tmp_path / "x.run.parquet"
    with SdrfConverter() as converter:
        converter.convert(str(path), run_output=str(output), mzml_dir=tmp_path)
    assert pq.read_table(output)["acquisition_start"].to_pylist() == [expected]


def test_run_acquisition_start_warns_once_for_missing_mzml_files(tmp_path, caplog):
    """Missing files get one warning; missing attributes and dates already in SDRF do not."""
    _mzml_header(tmp_path / "no_date.mzML", None)
    path = _sdrf(
        tmp_path / "in.sdrf.tsv",
        [
            _base_row("S1", "missing_a.raw", "LFQ"),
            _base_row("S2", "missing_b.raw", "LFQ"),
            _base_row("S3", "no_date.raw", "LFQ"),
            _base_row("S4", "sdrf_only.raw", "LFQ", **{"comment[acquisition date]": "2011-05-31"}),
        ],
    )
    output = tmp_path / "x.run.parquet"
    with caplog.at_level("WARNING"), SdrfConverter() as converter:
        converter.convert(str(path), run_output=str(output), mzml_dir=tmp_path)
    assert {row["run_file_name"]: row["acquisition_start"] for row in pq.read_table(output).to_pylist()} == {
        "missing_a": None,
        "missing_b": None,
        "no_date": None,
        "sdrf_only": "2011-05-31",
    }
    (warning,) = caplog.records
    assert warning.getMessage() == "No mzML file found for 2 runs; acquisition_start left null: missing_a, missing_b"


def test_run_acquisition_start_without_mzml_attribute_is_null(tmp_path, caplog):
    """An mzML run without startTimeStamp leaves a null value without warning."""
    _mzml_header(tmp_path / "a.mzML", None)
    path = _sdrf(tmp_path / "in.sdrf.tsv", [_base_row("S1", "a.raw", "LFQ")])
    output = tmp_path / "x.run.parquet"
    with caplog.at_level("WARNING"), SdrfConverter() as converter:
        converter.convert(str(path), run_output=str(output), mzml_dir=tmp_path)
    assert pq.read_table(output)["acquisition_start"].to_pylist() == [None]
    assert not caplog.records


def test_convert_sdrf_cli_reads_acquisition_start_from_mzml(tmp_path):
    """The SDRF CLI passes its optional mzML directory through to run metadata."""
    _mzml_header(tmp_path / "a.mzML.gz", "2015-02-19T20:08:48Z")
    path = _sdrf(tmp_path / "in.sdrf.tsv", [_base_row("S1", "a.raw", "LFQ")])
    output = tmp_path / "out"
    result = CliRunner().invoke(
        qpx_main,
        ["convert", "sdrf", "--sdrf-file", str(path), "--mzml-dir", str(tmp_path), "--output-folder", str(output)],
    )
    assert result.exit_code == 0, result.output
    assert pq.read_table(output / "sdrf.run.parquet")["acquisition_start"].to_pylist() == ["2015-02-19T20:08:48Z"]
