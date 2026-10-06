"""Native duplicates persist with warnings while strict audit reports errors."""

import logging
from copy import deepcopy

import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.validate import validate_cmd
from qpx.converters.openms.converter import OpenMSConverter
from qpx.core.data.loader import load_schema
from qpx.writers import FeatureWriter, PgWriter, PsmWriter
from tests.conftest import make_feature_record, make_pg_record, make_psm_record

_VIEW_INPUTS = {
    "feature": (FeatureWriter, make_feature_record),
    "psm": (PsmWriter, make_psm_record),
    "pg": (PgWriter, make_pg_record),
}


def _write_source_bundle(folder, conflict_view, conflict):
    """Put the conflicting records in separate row groups of a full bundle."""
    for view, (writer_class, factory) in _VIEW_INPUTS.items():
        first = factory()
        records = [first]
        if view == conflict_view:
            if conflict == "unidentified":
                first.update(sequence="", peptidoform="")
            if conflict == "provided_id":
                first["feature_id"] = 101
            second = deepcopy(first)
            if conflict == "provided_id":
                second["rt"] += 1.0
            elif conflict == "different_match":
                second.update(rt=121.5, posterior_error_probability=0.2)
            elif conflict == "different_quantity":
                second["intensities"][0]["intensity"] += 1.0
            records.append(second)
        kwargs = {"override_provided_ids": False} if view == "feature" else {}
        with writer_class(folder / f"native.{view}.parquet", batch_size=1, **kwargs) as writer:
            for record in records:
                writer.write_batch([record])


@pytest.mark.parametrize(
    ("view", "conflict"),
    [
        ("feature", "identical"),
        ("feature", "unidentified"),
        ("feature", "provided_id"),
        ("psm", "different_match"),
        ("pg", "different_quantity"),
    ],
)
def test_native_duplicate_warning_preserves_source_records(tmp_path, caplog, view, conflict):
    """Preserve duplicate rows, identities and quantities for a separate strict audit."""
    source = tmp_path / "source"
    source.mkdir()
    _write_source_bundle(source, view, conflict)
    original_inputs = {path: path.read_bytes() for path in source.iterdir()}
    output = tmp_path / "output"
    output.mkdir()
    for core_view in _VIEW_INPUTS:
        path = output / f"openms.{core_view}.parquet"
        path.write_bytes(f"existing {core_view} output".encode())

    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="qpx.writers.base"):
        OpenMSConverter(source).convert(output)

    assert f"Primary key ({view}_id) has 1 duplicate row" in caplog.text
    for core_view in _VIEW_INPUTS:
        before = pq.read_table(source / f"native.{core_view}.parquet")
        after = pq.read_table(output / f"openms.{core_view}.parquet")
        assert after.equals(before, check_metadata=False)
    assert {path: path.read_bytes() for path in source.iterdir()} == original_inputs
    assert not list(output.rglob("*.tmp"))

    duplicate_path = output / f"openms.{view}.parquet"
    table = pq.read_table(duplicate_path)
    assert table.num_rows == 2
    for strict, severity in ((False, "warning"), (True, "error")):
        result = load_schema(view).validate_full(table, strict=strict)
        assert result.is_valid is not strict
        assert [issue.severity for issue in result.issues if issue.check == "duplicate_pk"] == [severity]
    if conflict == "identical":
        audit = CliRunner().invoke(validate_cmd, ["--file", str(duplicate_path)])
        assert audit.exit_code == 1
        assert "duplicate row" in audit.output
