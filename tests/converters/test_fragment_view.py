"""The fragment view of DIA-NN and Spectronaut reports links to their features."""

from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pytest

from qpx.converters.diann.feature_adapter import DiannFeatureAdapter
from qpx.converters.spectronaut.feature_adapter import SpectronautFeatureAdapter
from qpx.dataset import Dataset

_EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
_DIANN = _EXAMPLES / "diann" / "full" / "diann_report.tsv.gz"
_SPECTRONAUT = _EXAMPLES / "spectronaut" / "small" / "spectronaut_report.tsv"


def _diann(report, tmp_path, **options):
    with DiannFeatureAdapter() as adapter:
        adapter.convert(diann_report=str(report), output_path=str(tmp_path / "x.feature.parquet"), **options)
        return adapter.write_fragments(str(tmp_path / "x.fragment.parquet"))


def _assert_links(tmp_path):
    features = set(pq.read_table(tmp_path / "x.feature.parquet", columns=["feature_id"]).column(0).to_pylist())
    table = pq.read_table(tmp_path / "x.fragment.parquet")
    assert table.num_rows
    assert set(table.column("feature_id").to_pylist()) <= features
    assert table.column("fragment_id").combine_chunks().unique().to_pylist()
    assert len(set(table.column("fragment_id").to_pylist())) == table.num_rows
    return table


def test_diann_lists_become_rows_that_link_to_features(tmp_path):
    assert _diann(_DIANN, tmp_path)
    table = _assert_links(tmp_path)
    first = table.filter(pc.field("fragment_index") == 1).slice(0, 1).to_pylist()[0]
    assert first["ion_type"] is None
    assert first["label"] == "LFQ"
    assert first["intensity"] > 0
    assert 0 <= first["correlation"] <= 1


def test_diann_unmeasured_fragments_have_no_row_and_named_ones_are_parsed(tmp_path):
    frame = pd.read_csv(_DIANN, sep="\t", nrows=1)
    frame["Fragment.Quant.Raw"] = "10;0;30;"
    frame["Fragment.Quant.Corrected"] = "9;0;29;"
    frame["Fragment.Correlations"] = "0.9;0;0.8;"
    frame["Fragment.Info"] = "y7-unknown^1/826.4;b2-H2O^1/100.1;y3-H2O^2/300.5;"
    report = tmp_path / "report.tsv"
    frame.to_csv(report, sep="\t", index=False)

    assert _diann(report, tmp_path)

    rows = pq.read_table(tmp_path / "x.fragment.parquet").sort_by("fragment_index").to_pylist()
    assert [row["fragment_index"] for row in rows] == [1, 3]
    assert rows[0]["ion_type"] == "y7" and rows[0]["neutral_loss"] is None and rows[0]["fragment_charge"] == 1
    assert rows[0]["mz"] == pytest.approx(826.4)
    assert rows[1]["neutral_loss"] == "H2O" and rows[1]["fragment_charge"] == 2
    assert rows[1]["corrected_intensity"] == 29.0


def test_diann_2_export_quant_columns_become_rows(tmp_path):
    """DIA-NN 2.x --export-quant writes Fr.N.Index/Id/Quantity/Score, the quantity being a peak height."""
    frame = pd.read_csv(_DIANN, sep="\t", nrows=1).drop(
        columns=["Fragment.Quant.Raw", "Fragment.Quant.Corrected", "Fragment.Correlations"]
    )
    slots = [(0, "y7^1/826.4", 10.0, 0.9), (2, "b2-H2O^2/100.1", 30.0, -0.2), (3, "y3^1/300.5", 0.0, 0.5), (-1, "", 0.0, 0.0)]
    for n, (index, name, quantity, score) in enumerate(slots):
        frame[f"Fr.{n}.Index"], frame[f"Fr.{n}.Id"] = index, name
        frame[f"Fr.{n}.Quantity"], frame[f"Fr.{n}.Score"] = quantity, score
    report = tmp_path / "report.tsv"
    frame.to_csv(report, sep="\t", index=False)

    assert _diann(report, tmp_path)

    rows = pq.read_table(tmp_path / "x.fragment.parquet").sort_by("fragment_index").to_pylist()
    assert [row["fragment_index"] for row in rows] == [1, 3]
    assert rows[0]["ion_type"] == "y7" and rows[0]["fragment_charge"] == 1
    assert rows[0]["mz"] == pytest.approx(826.4)
    assert rows[0]["peak_height"] == 10.0 and rows[0]["intensity"] is None
    assert rows[1]["neutral_loss"] == "H2O" and rows[1]["fragment_charge"] == 2
    assert rows[1]["additional_scores"] == [{"score_name": "diann_fragment_score", "score_value": -0.2, "higher_better": True}]


@pytest.mark.parametrize("prefix", ["^", ""])
def test_diann_decoy_fragments_link_to_their_own_feature(tmp_path, prefix):
    """A decoy with a target's peptidoform, a DIA-NN 2.x library decoy (^) or a generated one, is its own feature."""
    frame = pd.read_csv(_DIANN, sep="\t", nrows=1)
    decoy = frame.assign(**{column: prefix + frame[column] for column in ("Modified.Sequence", "Precursor.Id")})
    report = tmp_path / "report.tsv"
    pd.concat([frame.assign(Decoy=0), decoy.assign(Decoy=1)]).to_csv(report, sep="\t", index=False)

    assert _diann(report, tmp_path)

    features = pq.read_table(tmp_path / "x.feature.parquet", columns=["feature_id", "is_decoy"]).to_pylist()
    assert sorted(feature["is_decoy"] for feature in features) == [False, True]
    assert len({feature["feature_id"] for feature in features}) == 2
    linked = set(pq.read_table(tmp_path / "x.fragment.parquet", columns=["feature_id"]).column(0).to_pylist())
    assert linked == {feature["feature_id"] for feature in features}


def test_a_report_without_fragments_writes_nothing(tmp_path):
    frame = pd.read_csv(_DIANN, sep="\t", nrows=2).drop(columns=["Fragment.Quant.Raw"])
    report = tmp_path / "report.tsv"
    frame.to_csv(report, sep="\t", index=False)

    assert not _diann(report, tmp_path)
    assert not (tmp_path / "x.fragment.parquet").exists()


def test_diann_q_value_threshold_keeps_the_fragments_of_kept_features(tmp_path):
    frame = pd.read_csv(_DIANN, sep="\t", nrows=2)
    frame["Q.Value"] = [0.02, 0.001]
    report = tmp_path / "report.tsv"
    frame.to_csv(report, sep="\t", index=False)

    with DiannFeatureAdapter() as adapter:
        adapter.convert(diann_report=str(report), output_path=str(tmp_path / "x.feature.parquet"), qvalue_threshold=0.01)
        assert adapter.write_fragments(str(tmp_path / "x.fragment.parquet"), qvalue_threshold=0.01)

    features = pq.read_table(tmp_path / "x.feature.parquet", columns=["feature_id"]).column(0).to_pylist()
    assert len(features) == 1
    linked = set(pq.read_table(tmp_path / "x.fragment.parquet", columns=["feature_id"]).column(0).to_pylist())
    assert linked == set(features)


def test_spectronaut_rows_keep_each_fragment(tmp_path):
    with SpectronautFeatureAdapter() as adapter:
        adapter.convert(spectronaut_report=str(_SPECTRONAUT), output_path=str(tmp_path / "x.feature.parquet"))
        assert adapter.write_fragments(str(tmp_path / "x.fragment.parquet"))

    table = _assert_links(tmp_path)
    report = pd.read_csv(_SPECTRONAUT, sep="\t", dtype=str)
    assert table.num_rows == len(report)
    row = table.slice(0, 1).to_pylist()[0]
    assert row["ion_type"] and row["fragment_charge"] and row["mz"] > 0
    names = {score["score_name"] for score in row["additional_scores"]}
    assert "spectronaut_interference_score" in names


def test_validate_reports_a_fragment_of_no_feature(tmp_path):
    assert _diann(_DIANN, tmp_path)
    fragment = pq.read_table(tmp_path / "x.fragment.parquet").slice(0, 3)
    table = fragment.set_column(fragment.schema.get_field_index("feature_id"), "feature_id", pa.array([1, 2, 3]))
    pq.write_table(table, tmp_path / "x.fragment.parquet")

    with Dataset(str(tmp_path), file_prefix="x") as dataset:
        result = dataset.validate(["feature", "fragment"], strict=True)["fragment"]

    assert any("3 fragment row(s)" in issue.message for issue in result.errors)
