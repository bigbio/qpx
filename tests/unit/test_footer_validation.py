"""Validation of the footer contract and of identity composites."""

import pyarrow.parquet as pq
import pytest

from qpx.core.data.feature import Feature
from qpx.writers import FeatureWriter
from tests.conftest import make_feature_record


def _write_feature(path, records=None, **footer):
    """Write features, then replace or drop (None) footer keys."""
    with FeatureWriter(path, override_provided_ids=False) as w:
        w.write_batch(records or [make_feature_record()])
    table = pq.read_table(path)
    meta = dict(table.schema.metadata)
    for key, value in footer.items():
        if value is None:
            meta.pop(key.encode())
        else:
            meta[key.encode()] = value.encode()
    pq.write_table(table.replace_schema_metadata(meta), path)
    return path


def _issues(path, check, *, strict=False):
    return [issue for issue in Feature.from_file(path).validate(strict=strict).issues if issue.check == check]


def test_writer_footer_passes(tmp_path):
    """A footer written by the qpx writers satisfies the contract."""
    path = _write_feature(tmp_path / "a.feature.parquet")
    assert Feature.from_file(path).validate(strict=True).issues == []


@pytest.mark.parametrize(
    ("footer", "check", "severity", "text"),
    [
        ({"qpx_version": "2.0"}, "footer_version", "error", "supports specification 1.x"),
        ({"qpx_version": "abc"}, "footer_version", "error", "not a <major>.<minor>"),
        ({"qpx_version": "1.9"}, "footer_version", "warning", "newer than this reader"),
        ({"qpx_version": None}, "footer_version", "warning", "no qpx_version"),
        ({"file_type": "psm_file"}, "footer_file_type", "error", "'psm_file'"),
        ({"file_type": None}, "footer_file_type", "warning", "no file_type"),
        ({"identity_composite": "peptidoform,charge,bogus"}, "footer_identity", "error", "['bogus']"),
        ({"identity_composite": None}, "footer_identity", "warning", "no identity_composite"),
        ({"primary_key": "psm_id"}, "footer_identity", "error", "'psm_id'"),
    ],
)
def test_footer_contract(tmp_path, footer, check, severity, text):
    """Unknown major or malformed versions, another file type and bogus identity keys are errors."""
    path = _write_feature(tmp_path / "a.feature.parquet", **footer)
    [issue] = _issues(path, check)
    assert issue.severity == severity
    assert text in issue.message


def test_missing_footer_key_fails_strict_validation(tmp_path):
    """A footer without qpx_version is an error under strict validation."""
    path = _write_feature(tmp_path / "a.feature.parquet", qpx_version=None)
    [issue] = _issues(path, "footer_version", strict=True)
    assert issue.severity == "error"


def test_one_composite_with_several_ids(tmp_path):
    """Producer ids must not split one identity composite value into several rows."""
    records = [make_feature_record(), make_feature_record()]
    records[0]["feature_id"], records[1]["feature_id"] = 1, 2
    path = _write_feature(tmp_path / "a.feature.parquet", records)
    [issue] = _issues(path, "duplicate_composite")
    assert issue.severity == "warning"
    assert "1 identity composite value(s) (peptidoform, charge, run_file_name, rt) carry 2" in issue.message
    [strict_issue] = _issues(path, "duplicate_composite", strict=True)
    assert strict_issue.severity == "error"
