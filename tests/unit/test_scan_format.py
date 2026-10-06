"""Declared PSM scan formats survive streaming and reach public validation."""

import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

import qpx
from qpx.cli.main import qpx_main
from qpx.core.data import FeatureSchema, PsmSchema
from qpx.core.scan import scan_format_from_native_id
from qpx.writers.psm import PsmWriter
from tests.conftest import _valid_arrays, make_psm_record


@pytest.mark.parametrize(
    ("native_id", "expected"),
    [
        ("scan=42", "scan"),
        ("index=0", "index"),
        ("controllerType=0 controllerNumber=1 scan=42", "scan"),
        ("controllerType=5 controllerNumber=1 scan=7", "nativeId"),
        ("frame=120 scan=475", "nativeId"),
        ("frame=120 scan=475 precursor=3", "nativeId"),
        ("frame=120 windowGroup=2 scan=475", "nativeId"),
        ("merged=0 frame=120 scanStart=4 scanEnd=8", "nativeId"),
        ("function=10 process=1 scan=345", "nativeId"),
        ("sample=1 period=1 cycle=2740 experiment=10", "nativeId"),
        ("SCAN=42", "scan"),
        ("42", None),
        ("spectrum=42", None),
        ("prefix scan=42", None),
        ("vendor:scan=42", None),
        ("uuid=opaque", None),
        ("", None),
    ],
)
def test_scan_format_uses_native_keys(native_id, expected):
    """Format classification must not guess from a stored integer's value."""
    assert scan_format_from_native_id(native_id) == expected


def _write_psm(path, scan, scan_format=None):
    """Write one complete PSM with optional legacy or canonical metadata."""
    with PsmWriter(path, scan_format=scan_format) as writer:
        writer.write_batch([make_psm_record(scan=scan)])


def test_scan_format_can_be_declared_after_multiple_batches(tmp_path):
    """Late metadata must be visible in both the footer and restored Arrow schema."""
    path = tmp_path / "test.psm.parquet"
    with PsmWriter(path, batch_size=1) as writer:
        writer.write_batch([make_psm_record(scan=[1])])
        writer.write_batch([make_psm_record(scan=[2])])
        writer.set_scan_format("index")
    parquet = pq.ParquetFile(path)
    assert parquet.metadata.num_row_groups == 2
    expected = b"index" if hasattr(pq.ParquetWriter, "add_key_value_metadata") else None
    assert parquet.metadata.metadata.get(b"scan_format") == expected
    assert parquet.schema_arrow.metadata.get(b"scan_format") == expected
    assert parquet.read().column("scan").to_pylist() == [[1], [2]]


def test_metadata_failure_keeps_previous_destination(tmp_path, monkeypatch):
    """A caught metadata update failure still prevents publishing a partial file."""
    path = tmp_path / "test.psm.parquet"
    _write_psm(path, [1])
    previous = path.read_bytes()

    def fail_metadata_update(*_args):
        raise RuntimeError("metadata update failed")

    monkeypatch.setattr(pq.ParquetWriter, "add_key_value_metadata", fail_metadata_update, raising=False)
    with PsmWriter(path, batch_size=1) as writer:
        writer.write_batch([make_psm_record(scan=[2])])
        with pytest.raises(RuntimeError, match="metadata update failed"):
            writer.set_scan_format("scan")
    assert path.read_bytes() == previous
    assert not list(tmp_path.glob(".*.tmp"))


@pytest.mark.parametrize("flushed", [False, True])
def test_older_arrow_preserves_output_without_late_metadata(tmp_path, monkeypatch, caplog, flushed):
    """Older Arrow still writes rows, and supports declarations before opening."""
    monkeypatch.delattr(pq.ParquetWriter, "add_key_value_metadata", raising=False)
    path = tmp_path / "test.psm.parquet"
    with PsmWriter(path, batch_size=1 if flushed else 10) as writer:
        writer.write_batch([make_psm_record(scan=[1])])
        writer.set_scan_format("scan")
        expected = None if flushed else b"scan"
        assert writer.arrow_schema.metadata.get(b"scan_format") == expected
    assert pq.read_schema(path).metadata.get(b"scan_format") == expected
    assert pq.read_table(path).column("scan").to_pylist() == [[1]]
    assert len([record for record in caplog.records if "Cannot add scan_format" in record.message]) == int(flushed)


@pytest.mark.parametrize(
    ("scan_format", "scan"), [("scan", [1, 2]), ("index", [1, 2]), ("nativeId", [1]), ("nativeId", [1, 2, 3, 4, 5])]
)
def test_declared_format_reaches_schema_dataset_and_cli(tmp_path, scan_format, scan):
    """Public file and dataset validation must not lose the footer in DuckDB."""
    path = tmp_path / "test.psm.parquet"
    _write_psm(path, scan, scan_format)
    table = pq.read_table(path)
    result = PsmSchema.validate_full(table, strict=True)
    assert any(issue.check == "scan_format" for issue in result.errors)
    assert PsmSchema.validate_full(table).is_valid
    with qpx.Dataset(tmp_path, duckdb_threads=6) as dataset:
        assert not dataset.validate(structures=["psm"], strict=True)["psm"].is_valid
        assert dataset.psm.to_arrow().schema.metadata is None
    for args in (["--file", str(path)], ["--dataset-path", str(tmp_path), "--structure", "psm"]):
        result = CliRunner().invoke(qpx_main, ["validate", *args])
        assert result.exit_code == 1
        assert "scan_format" in result.output


@pytest.mark.parametrize(
    ("scan_format", "scan"),
    [("scan", [1]), ("index", [0]), ("nativeId", [1, 2]), ("nativeId", [1, 2, 3]), ("nativeId", [1, 2, 3, 4])],
)
def test_valid_declared_cardinalities_pass(tmp_path, scan_format, scan):
    """Every documented length is accepted under strict validation."""
    path = tmp_path / "test.psm.parquet"
    _write_psm(path, scan, scan_format)
    assert PsmSchema.validate_full(pq.read_table(path), strict=True).is_valid


@pytest.mark.skipif(
    not hasattr(pq.ParquetWriter, "add_key_value_metadata"), reason="Creating a footer-only fixture requires PyArrow 17+"
)
def test_footer_only_declaration_is_validated(tmp_path):
    """External producers may update footer keys without updating ARROW:schema."""
    path = tmp_path / "test.psm.parquet"
    _write_psm(path, [1, 2])
    table = pq.read_table(path)
    with pq.ParquetWriter(path, table.schema) as writer:
        writer.write_table(table)
        writer.add_key_value_metadata({b"scan_format": b"scan"})
    assert b"scan_format" not in pq.read_schema(path).metadata
    with qpx.Dataset(tmp_path, duckdb_threads=6) as dataset:
        assert any(issue.check == "scan_format" for issue in dataset.psm.validate(strict=True).errors)


@pytest.mark.parametrize(("scan_format", "scan"), [(None, [1, 2, 3, 4, 5]), ("native", [1]), ("scan", [])])
def test_legacy_or_empty_scans_remain_compatible(tmp_path, scan_format, scan):
    """Absent/unknown metadata and empty arrays do not become strict errors."""
    path = tmp_path / "test.psm.parquet"
    _write_psm(path, scan, scan_format)
    with qpx.Dataset(tmp_path, duckdb_threads=6) as dataset:
        result = dataset.validate(structures=["psm"], strict=True)["psm"]
    assert result.is_valid
    assert bool([issue for issue in result.warnings if issue.check == "scan_format"]) == (scan_format == "native")


def test_shards_do_not_inherit_first_files_format(tmp_path):
    """Mixed shard declarations must not be checked against the first footer."""
    _write_psm(tmp_path / "part0.psm.parquet", [1], "scan")
    _write_psm(tmp_path / "part1.psm.parquet", [1, 2], "nativeId")
    with qpx.Dataset(tmp_path, duckdb_threads=6) as dataset:
        assert dataset.validate(structures=["psm"], strict=True)["psm"].is_valid


def test_feature_supporting_scans_are_not_single_spectrum_ids():
    """Feature scans may flatten several native spectra, so PSM rules do not apply."""
    import pyarrow as pa

    schema = FeatureSchema.get_arrow_schema()
    arrays = _valid_arrays(schema)
    arrays["scan"] = pa.array([[1, 2, 3, 4, 5]], type=schema.field("scan").type)
    table = pa.table(arrays, schema=schema).replace_schema_metadata({b"scan_format": b"scan"})
    assert not [issue for issue in FeatureSchema.validate_full(table, strict=True).issues if issue.check == "scan_format"]
