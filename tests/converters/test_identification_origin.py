"""How a feature was identified, whether its value was imputed, and what it was read from."""

import re
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from qpx.converters.cdap.psm_adapter import _parent_ion_fraction
from qpx.converters.maxquant.evidence_fields import identification_fields
from qpx.converters.openms_consensus.feature_adapter import _feature_identification_fields
from qpx.converters.openms_consensus.provenance import _step_parameters
from qpx.converters.openms_consensus.quant_scan import quant_scan_fields
from qpx.converters.openms_consensus.streaming import read_provenance_header
from tests.conftest import make_feature_record

REPORT = Path(__file__).parents[1] / "examples" / "spectronaut" / "small" / "spectronaut_report.tsv"


@pytest.mark.parametrize(
    ("evidence_type", "expected"),
    [
        ("MSMS", "direct"),
        ("MULTI-MSMS", "direct"),
        ("ISO-MSMS", "direct"),
        ("MULTI-SECPEP", "secondary"),
        ("MULTI-MATCH", "transferred"),
        ("MULTI-MATCH-MSMS", "transferred"),
        ("SOMETHING-NEW", None),
    ],
)
def test_maxquant_type_gives_the_identification_type(evidence_type, expected):
    assert identification_fields({}, evidence_type)["identification_type"] == expected


def test_maxquant_match_metadata_is_kept_in_spec_units():
    row = {"Match score": 12.5, "Match q-value": 0.01, "Match time difference": -0.5, "Match m/z difference": 0.002}
    fields = identification_fields(row, "MULTI-MATCH")
    assert fields["match_score"] == 12.5
    assert fields["match_qvalue"] == 0.01
    assert fields["match_rt_delta"] == -30.0
    assert fields["match_mz_delta"] == 0.002


@pytest.mark.parametrize(("pif", "expected"), [(0.8, 0.8), (0.0, None), (float("nan"), None)])
def test_maxquant_zero_pif_means_not_computed(pif, expected):
    assert identification_fields({"PIF": pif}, "MSMS")["parent_ion_fraction"] == expected


@pytest.mark.parametrize(
    ("cell", "expected"),
    [("95.5,88.0", 0.955), ("100.0", 1.0), ("", None), (None, None)],
)
def test_cdap_purity_is_the_first_percentage_as_a_fraction(cell, expected):
    assert _parent_ion_fraction(cell) == expected


def _evidence(consistent=True, has_spectrum=True):
    return {"consistent": consistent, "has_spectrum": has_spectrum, "positions": set()}


@pytest.mark.parametrize(
    ("evidence_by_run", "expected_type", "expected_run"),
    [
        ({"a": _evidence()}, "direct", "a"),
        ({"a": _evidence(has_spectrum=False)}, None, None),
        ({"a": _evidence(consistent=False)}, None, None),
        ({"b": _evidence()}, "transferred", None),
        ({"b": _evidence(), None: _evidence()}, None, None),
    ],
)
def test_openms_separates_a_transfer_from_an_unknown_origin(evidence_by_run, expected_type, expected_run):
    fields = _feature_identification_fields("a", {}, evidence_by_run, set())
    assert fields["identification_type"] == expected_type
    assert fields["id_run_file_name"] == expected_run


class _Meta:
    def __init__(self, **values):
        self.values = values

    def metaValueExists(self, key):  # noqa: N802
        return key in self.values

    def getMetaValue(self, key):  # noqa: N802
        return self.values[key]


def test_openms_reporter_scan_metadata():
    ms3 = _Meta(
        precursor_purity=0.5,
        quant_precursor_purity=0.9,
        scan_id="controllerType=0 controllerNumber=1 scan=20",
        id_scan_id="controllerType=0 controllerNumber=1 scan=19",
    )
    assert quant_scan_fields(ms3) == {"parent_ion_fraction": 0.9, "quant_scan": [20], "quant_ms_level": 3}

    same = _Meta(precursor_purity=0.5, scan_id="scan=20", id_scan_id="scan=20")
    assert quant_scan_fields(same) == {"parent_ion_fraction": 0.5, "quant_scan": [20], "quant_ms_level": 2}

    # Only the two ids together tell the level.
    assert quant_scan_fields(_Meta(scan_id="scan=20")) == {"quant_scan": [20]}
    assert quant_scan_fields(_Meta()) == {}


def test_openms_header_parameters_reach_the_step(tmp_path):
    path = tmp_path / "x.consensusXML"
    path.write_text(
        '<?xml version="1.0"?><consensusXML>'
        '<dataProcessing><software name="IsobaricAnalyzer" version="3.5"/><processingAction name="Quantitation"/>'
        '<UserParam type="string" name="parameter: in" value="a.mzML"/>'
        '<UserParam type="int" name="parameter: threads" value="8"/>'
        '<UserParam type="string" name="parameter: type" value="tmt10plex"/>'
        '<UserParam type="string" name="parameter: quantification:isotope_correction" value="false"/>'
        '<UserParam type="string" name="parameter: quantification:normalization" value="true"/>'
        '<UserParam type="string" name="unrelated" value="x"/>'
        "</dataProcessing><mapList/></consensusXML>"
    )
    ((name, _, _, parameters),) = read_provenance_header(str(path)).data_processing
    assert name == "IsobaricAnalyzer"
    assert _step_parameters(parameters) == [
        {"key": "isotope_corrected", "value": "false"},
        {"key": "normalization", "value": "true"},
        {"key": "type", "value": "tmt10plex"},
        {"key": "quantification:isotope_correction", "value": "false"},
        {"key": "quantification:normalization", "value": "true"},
    ]
    assert _step_parameters(()) is None


def test_fragpipe_match_type_gives_the_identification_type(tmp_path):
    from qpx.converters.fragpipe.feature_adapter import FragPipeFeatureAdapter

    ion = tmp_path / "combined_ion.tsv"
    ion.write_text(
        "Peptide Sequence\tProtein\tM/Z\tCharge\tAssigned Modifications\t"
        "expA Intensity\texpA Match Type\texpB Intensity\texpB Match Type\n"
        "PEPTIDEK\tsp|P1|A_HUMAN\t500.0\t2\t\t1000\tMS/MS\t2000\tMBR\n"
    )
    out = tmp_path / "feature.parquet"
    with FragPipeFeatureAdapter() as adapter:
        adapter.convert(feature_path=str(ion), output_path=str(out), experiment_to_runs={"expA": ["run_A"], "expB": ["run_B"]})
    rows = {row["quantification_unit_id"]: row["identification_type"] for row in pq.read_table(out).to_pylist()}
    assert rows == {"expA": "direct", "expB": "transferred"}


@pytest.fixture(name="imputed_output")
def _imputed_output(tmp_path):
    """A Spectronaut report whose first run's peaks are imputed."""
    from qpx.converters.spectronaut.converter import SpectronautConverter

    lines = REPORT.read_text().splitlines()
    header = lines[0].split("\t")
    imputed, run = header.index("EG.IsImputed"), header.index("R.FileName")
    first_run = lines[1].split("\t")[run]
    rows = []
    for line in lines[1:]:
        cells = line.split("\t")
        cells[imputed] = "True" if cells[run] == first_run else "False"
        rows.append("\t".join(cells))
    report = tmp_path / "report.tsv"
    report.write_text("\n".join([lines[0], *rows]) + "\n")
    out = tmp_path / "out"
    out.mkdir()
    converter = SpectronautConverter(report_path=str(report))
    converter.convert_features(qvalue_threshold=1.0, output_folder=str(out), output_prefix="sn")
    converter.convert_pg(output_folder=str(out), output_prefix="sn")
    return out, first_run.rsplit(".", 1)[0]


def test_spectronaut_imputation_flags(imputed_output):
    out, first_run = imputed_output
    features = pq.read_table(next(out.glob("*.feature.parquet"))).to_pylist()
    assert {row["identification_type"] for row in features} == {"direct", "transferred"}
    for row in features:
        is_first = row["run_file_name"] == first_run
        assert row["identification_type"] == ("transferred" if is_first else "direct")
        assert row["imputed_labels"] == (["LFQ"] if is_first else None)
        assert row["id_run_file_name"] == (None if is_first else row["run_file_name"])
        assert not any(param["cv_name"] == "is_imputed" for param in row["cv_params"] or ())

    pg = pq.read_table(next(out.glob("*.pg.parquet"))).to_pylist()
    assert {row["is_imputed"] for row in pg} == {True, False}
    assert all(row["is_imputed"] == (row["grouped_runs"] == [first_run]) for row in pg)


def test_imputed_cells_become_a_mudata_layer(imputed_output):
    anndata = pytest.importorskip("anndata")
    from qpx.dataset import Dataset
    from qpx.mudata import build_mudata

    assert anndata
    out, first_run = imputed_output
    with Dataset(out) as dataset:
        mdata = build_mudata(dataset)
    for modality in ("precursors", "proteins"):
        adata = mdata[modality]
        mask = adata.layers["imputed_mask"].toarray().astype(bool)
        assert mask.any()
        assert not mask[adata.obs_names != first_run].any()
        # Every imputed cell is an observation of the first run.
        assert mask[adata.obs_names == first_run].sum() == mask.sum()


def test_validate_rejects_an_unknown_identification_type(tmp_path):
    import qpx
    from qpx.writers import FeatureWriter

    with FeatureWriter(tmp_path / "exp.feature.parquet") as writer:
        writer.write_batch([{**make_feature_record(), "identification_type": "maybe"}])
    with qpx.open_dataset(str(tmp_path)) as dataset:
        result = dataset.validate(["feature"], strict=True)["feature"]
    assert any(issue.check == "value_not_in_vocabulary" for issue in result.errors)


def test_mudata_without_the_column_has_no_layer(tmp_path):
    pytest.importorskip("anndata")
    from qpx.dataset import Dataset
    from qpx.mudata import build_mudata
    from qpx.writers import FeatureWriter

    with FeatureWriter(tmp_path / "exp.feature.parquet") as writer:
        writer.write_batch([make_feature_record()])
    table = pq.read_table(tmp_path / "exp.feature.parquet").drop_columns(["imputed_labels"])
    pq.write_table(table, tmp_path / "exp.feature.parquet")
    with Dataset(tmp_path) as dataset:
        mdata = build_mudata(dataset)
    assert "imputed_mask" not in mdata["precursors"].layers


def _documented_names(path: str) -> set[str]:
    """The ``intensity_name`` column of the name tables of a spec page."""
    text = (Path(__file__).parents[2] / "docs" / "spec" / path).read_text(encoding="utf-8")
    names = re.findall(r"^\| `([a-z0-9_<>]+)` \|", text, flags=re.MULTILINE)
    return set(names)


def _names_in(path: Path) -> set[str]:
    rows = pq.read_table(path, columns=["additional_intensities"]).column(0).to_pylist()
    return {entry["intensity_name"] for row in rows for label in row or () for entry in label["intensities"]}


def test_feature_intensity_names_are_documented(tmp_path):
    """Names drift unless each one a converter writes is in the spec's table."""
    from qpx.converters.diann.converter import DiaNNConverter

    examples = Path(__file__).parents[1] / "examples" / "diann" / "small"
    out = tmp_path / "out"
    out.mkdir()
    converter = DiaNNConverter(report_path=str(examples / "diann_report.tsv"), sdrf_path=str(examples / "PXD019909-DIA.sdrf.tsv"))
    converter.convert_features(qvalue_threshold=1.0, output_folder=str(out), output_prefix="d")
    names = _names_in(next(out.glob("*.feature.parquet")))
    assert "pg_maxlfq" in names
    assert names <= _documented_names("intensities.md")
