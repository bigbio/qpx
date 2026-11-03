"""Tests for the XIC view: schema, validation, writer, and point expansion."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.core.data import Xic, XicSchema
from qpx.writers import XicWriter


def _trace(**overrides):
    row = {
        "xic_id": "run1|PEPTIDEK2|ms1",
        "run_file_name": "run1",
        "label": "LFQ",
        "feature_ids": [1],
        "compensation_voltage": None,
        "precursor_id": "PEPTIDEK2",
        "peptidoform": "PEPTIDEK",
        "charge": 2,
        "transition": "ms1",
        "transition_type": "ms1",
        "fragment_type": None,
        "fragment_number": None,
        "fragment_charge": None,
        "neutral_loss": None,
        "isotope": 0,
        "target_mz": None,
        "target_rt": None,
        "rt_start": None,
        "rt_stop": None,
        "target_ion_mobility": None,
        "ion_mobility_start": None,
        "ion_mobility_stop": None,
        "library_intensity": None,
        "used_for_quantification": None,
        "rt": [60.0, 60.1, 60.2],
        "intensity": [100.0, 200.0, 150.0],
        "mz": None,
        "is_decoy": False,
        "additional_scores": None,
        "cv_params": None,
    }
    row.update(overrides)
    return row


def _fragment(**overrides):
    row = _trace(
        xic_id="run1|PEPTIDEK2|y7^1",
        transition="y7^1",
        transition_type="fragment",
        fragment_type="y",
        fragment_number=7,
        fragment_charge=1,
        target_mz=842.4512,
    )
    row.update(overrides)
    return row


def test_xic_schema_types():
    schema = XicSchema.get_arrow_schema()
    assert XicSchema.primary_key == ("xic_id",)
    assert schema.field("xic_id").type == pa.string()
    assert schema.field("charge").type == pa.int16()
    assert schema.field("rt").type == pa.list_(pa.float32())
    assert schema.field("intensity").type == pa.list_(pa.float32())
    assert schema.field("mz").type == pa.list_(pa.float32())
    assert schema.field("feature_ids").type == pa.list_(pa.int64())
    assert schema.field("feature_ids").nullable is True
    assert schema.field("isotope").type == pa.int16()
    cv_type = schema.field("cv_params").type
    assert isinstance(cv_type, pa.ListType) and isinstance(cv_type.value_type, pa.StructType)
    assert cv_type.value_type.get_field_index("cv_name") >= 0
    assert cv_type.value_type.get_field_index("cv_value") >= 0


def test_valid_xic_passes_strict_validation():
    table = pa.Table.from_pylist([_trace(), _fragment()], schema=XicSchema.get_arrow_schema())
    result = XicSchema.validate_full(table, strict=True)
    assert not result.issues


def test_rt_and_intensity_length_mismatch_is_reported():
    table = pa.Table.from_pylist(
        [_trace(rt=[60.0, 60.1], intensity=[100.0, 200.0, 150.0])],
        schema=XicSchema.get_arrow_schema(),
    )
    result = XicSchema.validate_full(table, strict=True)
    assert any(issue.check == "chromatogram_length_mismatch" for issue in result.issues)


def test_mz_array_length_mismatch_is_reported():
    table = pa.Table.from_pylist(
        [_trace(mz=[500.0, 500.1])],
        schema=XicSchema.get_arrow_schema(),
    )
    result = XicSchema.validate_full(table, strict=True)
    assert any(issue.check == "chromatogram_length_mismatch" for issue in result.issues)


def test_rt_must_be_strictly_increasing():
    table = pa.Table.from_pylist(
        [_trace(rt=[60.0, 60.0, 60.2])],
        schema=XicSchema.get_arrow_schema(),
    )
    result = XicSchema.validate_full(table, strict=True)
    assert any(issue.check == "rt_not_increasing" for issue in result.issues)


def test_transition_type_vocabulary():
    table = pa.Table.from_pylist(
        [_fragment(transition_type="precursor")],
        schema=XicSchema.get_arrow_schema(),
    )
    result = XicSchema.validate_full(table, strict=True)
    assert any(issue.check == "value_not_in_vocabulary" for issue in result.issues)


def test_writer_and_structure_round_trip(tmp_path):
    path = tmp_path / "test.xic.parquet"
    table = pa.Table.from_pylist([_trace(), _fragment()], schema=XicSchema.get_arrow_schema())
    with XicWriter(path) as writer:
        writer.write_table(table)

    reloaded = pq.read_table(path)
    assert reloaded.num_rows == 2

    xic = Xic.from_file(path)
    assert xic.ms1().count() == 1
    assert xic.fragments().count() == 1
    assert xic.by_precursor("PEPTIDEK2").count() == 2

    points = xic.points().to_df()
    assert len(points) == 6
    ms1_points = points[points["xic_id"] == "run1|PEPTIDEK2|ms1"].sort_values("trace_rt")
    assert ms1_points["trace_rt"].tolist() == pytest.approx([60.0, 60.1, 60.2], rel=1e-6)
    assert ms1_points["trace_intensity"].tolist() == [100.0, 200.0, 150.0]
