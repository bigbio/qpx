"""Spectronaut pivot and labelled reports, and the columns the feature view maps from the long layout."""

from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import pytest

from qpx.converters.spectronaut.converter import SpectronautConverter
from qpx.converters.spectronaut.report_layout import channel_labels, long_columns, melted_report_sql, pivot_columns

_REPORT = Path(__file__).parents[1] / "examples" / "spectronaut" / "small" / "spectronaut_report.tsv"
_RUNS = ("runA.raw", "runB.raw")


def _pivot_report(path, metric_values):
    """A pivot report of one precursor, as Spectronaut exports it; *metric_values* maps a metric to its per-run values."""
    row = {
        "PG.ProteinAccessions": "P12345;P67890",
        "PG.Genes": "GENE1",
        "PEP.StrippedSequence": "PEPTIDEK",
        "EG.PrecursorId": "_PEPTIDEM[Oxidation (M)]K_.3",
        "EG.ModifiedSequence": "_PEPTIDEM[Oxidation (M)]K_",
    }
    for metric, values in metric_values.items():
        for number, (run, value) in enumerate(zip(_RUNS, values), start=1):
            row[f"[{number}] {run}.{metric}"] = value
    pd.DataFrame([row]).to_csv(path, sep="\t", index=False)
    return path


def _features(tmp_path, report, sdrf=None):
    converter = SpectronautConverter(report_path=str(report), sdrf_path=sdrf)
    converter.convert_features(output_folder=str(tmp_path), output_prefix="sn")
    return pq.read_table(tmp_path / "sn.feature.parquet").to_pylist()


def _sdrf(tmp_path, labels):
    path = tmp_path / "sdrf.tsv"
    pd.DataFrame({"comment[label]": labels}).to_csv(path, sep="\t", index=False)
    return str(path)


def test_pivot_channels_melt_to_one_feature_per_run_with_labels(tmp_path):
    report = _pivot_report(
        tmp_path / "pivot.tsv",
        {"EG.Channel1Quantity": ["1000.5", "Filtered"], "EG.Channel2Quantity": ["2000.5", "3000.5"]},
    )

    features = {row["run_file_name"]: row for row in _features(tmp_path, report)}

    assert set(features) == {"runA", "runB"}
    assert features["runA"]["charge"] == 3
    assert features["runA"]["peptidoform"] == "PEPTIDEM[UNIMOD:35]K"
    assert features["runA"]["intensities"] == [
        {"label": "SILAC light", "intensity": 1000.5},
        {"label": "SILAC heavy", "intensity": 2000.5},
    ]
    assert features["runB"]["intensities"] == [
        {"label": "SILAC light", "intensity": None},
        {"label": "SILAC heavy", "intensity": 3000.5},
    ]
    assert [p["accession"] for p in features["runA"]["pg_accessions"]] == ["P12345", "P67890"]


def test_pivot_drops_a_run_where_the_precursor_was_not_quantified(tmp_path):
    report = _pivot_report(
        tmp_path / "pivot.tsv",
        {"EG.Channel1Quantity": ["1000.5", "Filtered"], "EG.Channel2Quantity": ["2000.5", "Filtered"]},
    )

    assert [row["run_file_name"] for row in _features(tmp_path, report)] == ["runA"]


def test_label_free_pivot_reads_comma_decimals(tmp_path):
    report = _pivot_report(tmp_path / "pivot.tsv", {"FG.Quantity": ["1000,5", "NaN"]})

    features = _features(tmp_path, report)

    assert [(row["run_file_name"], row["intensities"]) for row in features] == [("runA", [{"label": "LFQ", "intensity": 1000.5}])]


def test_channels_take_the_sdrf_labels_lightest_first(tmp_path):
    report = _pivot_report(
        tmp_path / "pivot.tsv",
        {"EG.Channel1Quantity": ["1.5", "1.5"], "EG.Channel2Quantity": ["2.5", "2.5"], "EG.Channel3Quantity": ["3.5", "3.5"]},
    )
    sdrf = _sdrf(tmp_path, ["SILAC heavy", "SILAC light", "SILAC medium"])

    feature = _features(tmp_path, report, sdrf)[0]

    assert [entry["label"] for entry in feature["intensities"]] == ["SILAC light", "SILAC medium", "SILAC heavy"]


def test_an_sdrf_that_disagrees_with_the_channels_raises():
    with pytest.raises(ValueError, match="3 quantity channels"):
        channel_labels([1, 2, 3], {"SILAC light", "SILAC heavy"})
    with pytest.raises(ValueError, match="2 quantity channels"):
        channel_labels([1, 2], {"label free sample"})


def test_labels_without_an_sdrf():
    assert channel_labels([1], None) == {1: "LFQ"}
    assert channel_labels([1, 2], None) == {1: "SILAC light", 2: "SILAC heavy"}
    assert channel_labels([1, 2, 3], None) == {1: "SILAC light", 2: "SILAC medium", 3: "SILAC heavy"}
    with pytest.raises(ValueError, match="pass an SDRF"):
        channel_labels([1, 2, 3, 4], None)


def test_long_labelled_report_maps_channels_and_has_no_channel_free_quantities(tmp_path):
    report = pd.read_csv(_REPORT, sep="\t", nrows=1, low_memory=False, decimal=",")
    report["EG.Channel1Quantity"] = 11.5
    report["EG.Channel2Quantity"] = 22.5
    path = tmp_path / "report.tsv"
    report.to_csv(path, sep="\t", index=False)

    feature = _features(tmp_path, path)[0]

    assert feature["intensities"] == [
        {"label": "SILAC light", "intensity": 11.5},
        {"label": "SILAC heavy", "intensity": 22.5},
    ]
    assert feature["additional_intensities"] is None


def test_labelled_report_has_no_pg_view_and_protein_pivot_no_feature_view(tmp_path):
    labelled = _pivot_report(
        tmp_path / "labelled.tsv", {"EG.Channel1Quantity": ["1.5", "2.5"], "EG.Channel2Quantity": ["3.5", "4.5"]}
    )
    assert SpectronautConverter(report_path=str(labelled)).available_views() == ["feature"]

    protein = tmp_path / "protein.tsv"
    pd.DataFrame([{"PG.ProteinGroups": "P12345", "PG.Genes": "GENE1", "[1] runA.raw.PG.Quantity": "12.5"}]).to_csv(
        protein, sep="\t", index=False
    )
    converter = SpectronautConverter(report_path=str(protein))
    assert converter.available_views() == ["pg"]
    converter.convert_pg(output_folder=str(tmp_path), output_prefix="sn")
    row = pq.read_table(tmp_path / "sn.pg.parquet").to_pylist()[0]
    assert row["grouped_runs"] == ["runA"]
    assert (row["label"], row["intensity"]) == ("LFQ", 12.5)
    assert row["peptides"] == [] and row["peptide_counts"] is None


def test_report_without_required_columns_raises_a_clear_error(tmp_path):
    path = tmp_path / "report.tsv"
    pd.DataFrame([{"PEP.StrippedSequence": "PEPTIDEK", "R.FileName": "run"}]).to_csv(path, sep="\t", index=False)

    converter = SpectronautConverter(report_path=str(path))
    with pytest.raises(ValueError, match="charge.*pg_accessions"):
        converter.convert_features(output_folder=str(tmp_path), output_prefix="sn")
    with pytest.raises(ValueError, match="not a convertible Spectronaut report"):
        converter.available_views()


def test_feature_view_maps_measured_mz_rt_bounds_and_ion_mobility(tmp_path):
    report = pd.read_csv(_REPORT, sep="\t", nrows=1, low_memory=False, decimal=",")
    report["EG.StartRT"] = 10.0
    report["EG.EndRT"] = 11.0
    report["FG.ApexIonMobility"] = 0.85
    path = tmp_path / "report.tsv"
    report.to_csv(path, sep="\t", index=False)

    feature = _features(tmp_path, path)[0]

    assert feature["observed_mz"] == pytest.approx(float(report["FG.PrecMzCalibrated"].iloc[0]))
    assert feature["observed_mz"] != pytest.approx(feature["calculated_mz"], abs=1e-3)
    assert (feature["rt_start"], feature["rt_stop"]) == (600.0, 660.0)
    assert feature["ion_mobility"] == pytest.approx(0.85)
    names = [e["intensity_name"] for entry in feature["additional_intensities"] for e in entry["intensities"]]
    assert names == ["pg_quantity"]


def test_a_zero_ion_mobility_is_no_measurement(tmp_path):
    report = pd.read_csv(_REPORT, sep="\t", nrows=1, low_memory=False, decimal=",")
    report["FG.ApexIonMobility"] = 0.0
    path = tmp_path / "report.tsv"
    report.to_csv(path, sep="\t", index=False)

    assert _features(tmp_path, path)[0]["ion_mobility"] is None


def test_pivot_parsing_keeps_dots_in_run_names():
    columns = ["EG.PrecursorId", "[1] a.b.raw.EG.Channel1Quantity", "[2] c.raw.EG.TotalQuantity (Settings)"]

    assert pivot_columns(columns) == {
        "a.b.raw": {"EG.Channel1Quantity": columns[1]},
        "c.raw": {"EG.TotalQuantity (Settings)": columns[2]},
    }
    assert long_columns(columns) == {
        "EG.PrecursorId", "FG.Charge", "R.FileName", "EG.Channel1Quantity", "EG.TotalQuantity (Settings)",
    }  # fmt: skip
    assert "report_pivot" in melted_report_sql(columns, "report_pivot")
