"""Regression tests for converter-independent QPX file-prefix handling."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.transforms.sdrf_updater import _check_protected_changes, update_metadata
from qpx.transforms.utils import discover_qpx_file_prefix


@pytest.mark.parametrize("prefix", ["openms", "quantms"])
def test_discover_qpx_file_prefix(tmp_path, prefix):
    (tmp_path / f"{prefix}.feature.parquet").touch()
    (tmp_path / f"{prefix}.pg.parquet").touch()

    assert discover_qpx_file_prefix(tmp_path) == prefix


def test_discover_qpx_file_prefix_rejects_ambiguous_dataset(tmp_path):
    (tmp_path / "openms.feature.parquet").touch()
    (tmp_path / "legacy.pg.parquet").touch()

    with pytest.raises(ValueError, match="Multiple QPX file prefixes"):
        discover_qpx_file_prefix(tmp_path)


def test_update_metadata_preserves_discovered_prefix(tmp_path, monkeypatch):
    original = pa.table({"value": ["old"]})
    pq.write_table(original, tmp_path / "openms.sample.parquet")
    pq.write_table(original, tmp_path / "openms.run.parquet")
    new_sdrf = tmp_path / "metadata.sdrf.tsv"
    new_sdrf.write_text("source name\nS1\n", encoding="utf-8")

    class FakeSdrfConverter:
        def convert(self, *, sdrf_path, sample_output, run_output):
            del sdrf_path
            updated = pa.table({"value": ["new"]})
            pq.write_table(updated, sample_output)
            pq.write_table(updated, run_output)

    monkeypatch.setattr("qpx.converters.sdrf.SdrfConverter", FakeSdrfConverter)

    result = update_metadata(tmp_path, new_sdrf)

    assert result["samples_updated"] == 1
    assert result["runs_updated"] == 1
    assert (tmp_path / "openms.sample.parquet.bak").is_file()
    assert (tmp_path / "openms.run.parquet.bak").is_file()
    assert not (tmp_path / "quantms.sample.parquet").exists()
    assert not (tmp_path / "quantms.run.parquet").exists()


def test_protected_check_compares_sdrf_values_as_text(tmp_path):
    """Renaming source name 001 to 1 changes the sample-run-label mapping."""
    header = "source name\tcomment[data file]\tcomment[label]\n"
    old, new = tmp_path / "old.sdrf.tsv", tmp_path / "new.sdrf.tsv"
    old.write_text(header + "001\ta.raw\tlabel free sample\n")
    new.write_text(header + "1\ta.raw\tlabel free sample\n")

    assert any("sample-run-label mapping changed" in warning for warning in _check_protected_changes(old, new))
