"""Value ranges declared in the YAML schemas, and the single-row dataset view."""

import math

import pyarrow as pa
import pytest

from qpx.core.data import DatasetSchema, FeatureSchema, PgSchema, PsmSchema
from tests.conftest import make_dataset_record, make_feature_record, make_pg_record, make_psm_record


def _range_issues(schema_class, record, *, strict):
    table = pa.Table.from_pylist([record], schema=schema_class.get_arrow_schema())
    result = schema_class.validate_full(table, strict=strict)
    return [issue for issue in result.issues if issue.check == "value_out_of_range"]


def _with(record, **overrides):
    return {**record, **overrides}


def _feature_intensity(value):
    return _with(make_feature_record(), intensities=[{"label": "TMT126", "intensity": value}])


@pytest.mark.parametrize(
    ("schema_class", "record", "path"),
    [
        (FeatureSchema, _with(make_feature_record(), posterior_error_probability=1.5), "posterior_error_probability"),
        (FeatureSchema, _with(make_feature_record(), pg_global_qvalue=math.nan), "pg_global_qvalue"),
        (FeatureSchema, _feature_intensity(math.nan), "intensities.intensity"),
        (FeatureSchema, _feature_intensity(-1.0), "intensities.intensity"),
        (FeatureSchema, _feature_intensity(math.inf), "intensities.intensity"),
        (PsmSchema, _with(make_psm_record(), posterior_error_probability=-0.1), "posterior_error_probability"),
        (PgSchema, _with(make_pg_record(), intensity=-5.0), "intensity"),
        (PgSchema, _with(make_pg_record(), sequence_coverage=150.0), "sequence_coverage"),
        (PgSchema, _with(make_pg_record(), global_qvalue=math.nan), "global_qvalue"),
    ],
)
def test_value_outside_its_declared_range(schema_class, record, path):
    """NaN, infinities and values outside the declared range are warnings, errors when strict."""
    [issue] = _range_issues(schema_class, record, strict=False)
    assert issue.severity == "warning"
    assert f"of {path} are not finite numbers" in issue.message
    [strict_issue] = _range_issues(schema_class, record, strict=True)
    assert strict_issue.severity == "error"


@pytest.mark.parametrize(
    ("schema_class", "record"),
    [
        (FeatureSchema, _with(make_feature_record(), posterior_error_probability=1.0, pg_global_qvalue=0.0)),
        (FeatureSchema, _feature_intensity(0.0)),
        (FeatureSchema, _feature_intensity(None)),
        (PgSchema, _with(make_pg_record(), intensity=None, sequence_coverage=100.0, global_qvalue=0.0)),
    ],
)
def test_values_at_the_bounds_and_nulls_pass(schema_class, record):
    """The bounds are inclusive, and a null is a missing value, not an out-of-range one."""
    assert _range_issues(schema_class, record, strict=True) == []


@pytest.mark.parametrize("rows", [0, 2])
def test_dataset_view_must_have_one_row(rows):
    """dataset.parquet describes the dataset in exactly one row."""
    records = [make_dataset_record(project_accession=f"PXD00000{i}") for i in range(rows)]
    table = pa.Table.from_pylist(records, schema=DatasetSchema.get_arrow_schema())
    [issue] = [i for i in DatasetSchema.validate_full(table).issues if i.check == "row_count"]
    assert issue.severity == "error"
    assert f"{rows} rows" in issue.message
