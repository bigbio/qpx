"""Tests for QPX spec-version stamping and the PG/feature layout load guards."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

import qpx
from qpx.version import (
    QPX_SPEC_VERSION,
    QpxVersionError,
    check_pg_columns_compatible,
    check_pg_file_compatible,
    parse_spec_version,
)
from qpx.writers import PgWriter
from tests.conftest import make_pg_record


def test_spec_version_is_the_format_version():
    """The spec version is the two-component format version, distinct from the build."""
    assert QPX_SPEC_VERSION == "1.3"
    assert parse_spec_version(QPX_SPEC_VERSION) == (1, 3)
    assert parse_spec_version("1.0") < parse_spec_version("1.1")


def test_written_pg_footer_carries_spec_version(tmp_path):
    """A freshly written 1.1 PG file passes the load-compatibility guard unchanged.

    (The qpx_version/writer_version footer stamps are asserted by
    test_writers.py::test_writer_footer_metadata.)
    """
    path = tmp_path / "test.pg.parquet"
    with PgWriter(path) as w:
        w.write_batch([make_pg_record()])

    # A well-formed 1.1 file passes the guard unchanged.
    check_pg_file_compatible(path)


def _write_old_pg_file(path, qpx_version="1.0"):
    """Synthesize a pre-1.1 PG file: scalar run_file_name, no grouped_runs."""
    table = pa.table({"anchor_protein": ["P12345"], "run_file_name": ["run_01"], "global_qvalue": [0.005]})
    table = table.replace_schema_metadata({b"qpx_version": qpx_version.encode()})
    pq.write_table(table, str(path))


def test_guard_raises_on_old_schema_pg_file(tmp_path):
    """A pre-1.1 PG file (scalar run_file_name) raises a clear version error."""
    path = tmp_path / "old.pg.parquet"
    _write_old_pg_file(path, qpx_version="1.0")

    with pytest.raises(QpxVersionError) as exc:
        check_pg_file_compatible(path)

    msg = str(exc.value)
    assert "grouped_runs" in msg
    assert "run_file_name" in msg
    assert "1.0" in msg  # names the offending file version
    assert "1.1" in msg  # names the required version


def test_read_pg_raises_on_old_schema_file(tmp_path):
    """The public read_pg() load path surfaces the clear error, not a raw DuckDB one."""
    path = tmp_path / "old.pg.parquet"
    _write_old_pg_file(path, qpx_version="1.0")

    with pytest.raises(QpxVersionError):
        qpx.read_pg(str(path))


def test_guard_noop_on_unstamped_complete_grouped_file(tmp_path):
    """An unstamped file with the complete flattened layout is allowed."""
    path = tmp_path / "grouped.pg.parquet"
    table = pa.table(
        {
            "anchor_protein": ["P1"],
            "grouped_runs": [["run_01"]],
            "label": [None],
            "intensity": [None],
        }
    )
    pq.write_table(table, str(path))
    check_pg_file_compatible(path)  # must not raise


def test_column_guard_rejects_pre_1_1_and_passes_1_1():
    """The column-based guard (used on the S3 path, where PyArrow cannot read the
    footer of a glob) rejects a scalar run_file_name layout and passes grouped_runs."""
    # pre-1.1: scalar run_file_name, no grouped_runs -> hard error
    with pytest.raises(QpxVersionError):
        check_pg_columns_compatible(
            ["anchor_protein", "run_file_name", "intensities"],
            source="s3://bucket/quantms/datasets/x/pg.parquet",
        )
    # explicit pre-1.1 stamp, still missing grouped_runs -> error
    with pytest.raises(QpxVersionError):
        check_pg_columns_compatible(["anchor_protein", "run_file_name"], source="s3://b/x", version_str="1.0")
    # pre-flatten 1.1: has grouped_runs but still the old intensities list and no
    # scalar label/intensity -> error (would otherwise pass the grouped_runs check
    # and then fail every pg query, silently)
    with pytest.raises(QpxVersionError, match="flattened|label"):
        check_pg_columns_compatible(["anchor_protein", "grouped_runs", "intensities"], source="s3://b/x")
    # incomplete 1.1: grouped_runs but neither the old list nor both flat columns
    with pytest.raises(QpxVersionError, match="label.*intensity|scalar"):
        check_pg_columns_compatible(["anchor_protein", "grouped_runs"], source="s3://b/x")
    # flat 1.1 layout -> no raise
    check_pg_columns_compatible(["anchor_protein", "grouped_runs", "label", "intensity"], source="s3://b/x")


def test_column_guard_rejects_old_s3_layout():
    """S3 discovery can reject old PG columns without a local PyArrow file."""
    with pytest.raises(QpxVersionError, match="s3://bucket/data"):
        check_pg_columns_compatible(
            ["anchor_protein", "run_file_name"],
            "s3://bucket/data/*.pg.parquet",
        )


def test_dataset_s3_discovery_fails_soft_on_version_error(monkeypatch, caplog):
    """The S3 registration path skips an incompatible PG layout with a warning.

    Dataset construction must not throw at __init__ because one structure file
    is an old/incompatible version; the pg structure is simply skipped.
    """
    from qpx.dataset import Dataset

    class FakeEngine:
        """Minimal DuckDBEngine replacement for the public S3 load path."""

        dropped = False

        def __init__(self, **_kwargs):
            pass

        @staticmethod
        def register_s3_parquet(name, path):
            """Validate the PG view registration request."""
            assert name == "pg"
            assert path == "s3://bucket/data/*.pg.parquet"

        @classmethod
        def execute(cls, sql):
            """Validate and return the simulated DESCRIBE result."""
            if sql == 'DESCRIBE "pg"':
                return cls()
            assert sql == "DROP VIEW IF EXISTS pg"
            cls.dropped = True
            return cls()

        @staticmethod
        def fetchall():
            """Return columns from an incompatible pre-1.1 PG layout."""
            return [
                ("anchor_protein", "VARCHAR", "YES", None, None, None),
                ("run_file_name", "VARCHAR", "YES", None, None, None),
            ]

    monkeypatch.setattr("qpx.dataset.DuckDBEngine", FakeEngine)

    with caplog.at_level("WARNING"):
        ds = Dataset("s3://bucket/data", structures=["pg"])

    assert ds.pg is None
    assert "pg" not in ds.available_structures
    assert FakeEngine.dropped
    assert any("pg" in rec.message for rec in caplog.records)


def test_partitioned_pg_checks_every_part_file(tmp_path, caplog):
    """A later old-schema partition is detected and the pg structure is skipped soft."""
    from qpx.dataset import Dataset

    part_dir = tmp_path / "pg"
    part_dir.mkdir()
    pq.write_table(
        pa.table(
            {
                "anchor_protein": ["P1"],
                "grouped_runs": [["run_01"]],
                "label": [None],
                "intensity": [None],
            }
        ),
        part_dir / "a.parquet",
    )
    _write_old_pg_file(part_dir / "b.parquet")

    with caplog.at_level("WARNING"):
        ds = Dataset(tmp_path, structures=["pg"])

    assert ds.pg is None
    assert any("pg" in rec.message for rec in caplog.records)


def test_parse_spec_version_is_graceful_on_garbage():
    """parse_spec_version never raises; garbage/partial inputs return None or a tuple."""
    assert parse_spec_version("1.1") == (1, 1)
    assert parse_spec_version("1") == (1, 0)
    assert parse_spec_version("1.1.3") == (1, 1)  # stray patch ignored
    for bad in ("v1.1", "", None, "abc", "1.x", ".", "  "):
        assert parse_spec_version(bad) is None


def test_dataset_with_one_bad_pg_file_still_opens_other_structures(tmp_path):
    """A Dataset with an old pg file still opens and exposes its other structures."""
    from qpx.dataset import Dataset

    prefix = tmp_path.name
    # A well-formed sample structure alongside a pre-1.1 pg file.
    pq.write_table(
        pa.table({"sample_accession": ["S1"], "condition": ["control"]}),
        tmp_path / f"{prefix}.sample.parquet",
    )
    _write_old_pg_file(tmp_path / f"{prefix}.pg.parquet")

    ds = Dataset(tmp_path)  # must not raise at construction

    assert ds.pg is None  # bad pg skipped
    assert ds.sample is not None  # other structure usable
    assert "sample" in ds.available_structures
    assert ds.sample.to_df().shape[0] == 1


def _write_channel_feature_file(path):
    """Synthesize a quantms.io-era feature file whose intensities are keyed by ``channel``."""
    intensity = pa.list_(pa.struct([("sample_accession", pa.string()), ("channel", pa.string()), ("intensity", pa.float32())]))
    pq.write_table(
        pa.table(
            {
                "sequence": ["PEPTIDE"],
                "run_file_name": ["run_01"],
                "intensities": pa.array([[{"sample_accession": "S1", "channel": "legacy", "intensity": 42.0}]], type=intensity),
            }
        ),
        str(path),
    )


def test_channel_keyed_feature_is_rejected(tmp_path, caplog):
    """Feature.from_file raises on the ``channel`` layout and a Dataset skips that view."""
    from qpx.core.data.feature import Feature
    from qpx.dataset import Dataset

    path = tmp_path / "old.feature.parquet"
    _write_channel_feature_file(path)

    with pytest.raises(QpxVersionError, match="'channel'"):
        Feature.from_file(path)

    with caplog.at_level("WARNING"):
        ds = Dataset(tmp_path, file_prefix="old")
    assert ds.feature is None
    assert any("feature" in rec.message for rec in caplog.records)


def test_feature_with_another_intensity_layout_reaches_validation(tmp_path):
    """A feature file whose intensities are not label structs opens, and validation reports the type."""
    from qpx.dataset import Dataset

    pq.write_table(
        pa.table(
            {
                "sequence": ["PEPTIDE"],
                "run_file_name": ["run_01"],
                "intensities": pa.array([[1.0]], type=pa.list_(pa.float64())),
            }
        ),
        str(tmp_path / "bad.feature.parquet"),
    )

    issues = Dataset(tmp_path, file_prefix="bad").validate()["feature"].issues
    assert any(issue.check == "type_mismatch" and issue.column == "intensities" for issue in issues)
