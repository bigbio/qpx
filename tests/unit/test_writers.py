"""Tests for QPX writer layer: BaseWriter subclasses and round-trip validation."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.core.data import FeatureSchema, PsmSchema
from qpx.core.parquet_io import parquet_row_count, read_parquet_metadata
from qpx.version import QPX_SPEC_VERSION
from qpx.writers import (
    DatasetWriter,
    FeatureWriter,
    OntologyWriter,
    PgWriter,
    ProvenanceWriter,
    PsmWriter,
    RunWriter,
    SampleWriter,
)
from tests.conftest import (
    make_dataset_record,
    make_feature_record,
    make_ontology_record,
    make_pg_record,
    make_provenance_record,
    make_psm_record,
    make_run_record,
    make_sample_record,
)


def test_writers_create_valid_parquet(tmp_path):
    """All writers produce readable Parquet files with correct row count."""
    configs = [
        (FeatureWriter, make_feature_record, "feature"),
        (PsmWriter, make_psm_record, "psm"),
        (PgWriter, make_pg_record, "pg"),
        (SampleWriter, make_sample_record, "sample"),
        (RunWriter, make_run_record, "run"),
        (DatasetWriter, make_dataset_record, "dataset"),
        (OntologyWriter, make_ontology_record, "ontology"),
        (ProvenanceWriter, make_provenance_record, "provenance"),
    ]
    for writer_cls, record_fn, name in configs:
        path = tmp_path / f"test.{name}.parquet"
        with writer_cls(path) as w:
            w.write_batch([record_fn()])
        assert path.exists()
        table = pq.read_table(path)
        assert table.num_rows == 1


def test_writer_footer_metadata(tmp_path):
    """Footer metadata contains file_type, version, creator, software_provider, scan_format."""
    path = tmp_path / "test.feature.parquet"
    with FeatureWriter(path, creator="test_suite", software_provider="my_tool", scan_format="native") as w:
        w.write_batch([make_feature_record()])

    meta = read_parquet_metadata(path)
    assert meta["file_type"] == "feature_file"
    assert meta["creator"] == "test_suite"
    # qpx_version is the on-disk SPEC version, not the package build version.
    assert meta["qpx_version"] == QPX_SPEC_VERSION == "1.1"
    # The package build version is stamped separately for provenance.
    assert "writer_version" in meta
    assert "uuid" in meta
    assert "creation_date" in meta
    assert meta["software_provider"] == "my_tool"
    assert meta["scan_format"] == "native"

    # Metadata should be in footer, not in data columns
    table = pq.read_table(path)
    data_columns = table.schema.names
    assert "qpx_version" not in data_columns
    assert "file_type" not in data_columns


def test_writer_footer_supports_composite_primary_key(tmp_path):
    """A multi-column primary key must not be treated as a derived id field."""
    path = tmp_path / "composite.ontology.parquet"
    with OntologyWriter(path, identity_composite=("field_name", "view")) as writer:
        writer.write_batch([make_ontology_record()])

    metadata = read_parquet_metadata(path)
    assert metadata["identity_composite"] == "field_name,view"
    assert "primary_key" not in metadata


def test_writer_schema_validation(tmp_path):
    """Writer rejects bad types and accepts valid tables."""
    # Reject bad types
    path = tmp_path / "bad.feature.parquet"
    writer = FeatureWriter(path)
    try:
        schema = FeatureSchema.get_arrow_schema()
        wrong_fields = [pa.field("charge", pa.string()) if f.name == "charge" else f for f in schema]
        wrong_schema = pa.schema(wrong_fields)
        arrays = {f.name: pa.nulls(1, type=f.type) for f in wrong_schema}
        arrays["peptidoform"] = pa.array(["PEPTIDE"], type=pa.string())
        arrays["charge"] = pa.array(["wrong"], type=pa.string())
        table = pa.table(arrays, schema=wrong_schema)
        with pytest.raises(ValueError, match="Schema validation failed"):
            writer.write_table(table)
    finally:
        writer.close()

    # Accept valid table
    path2 = tmp_path / "good.feature.parquet"
    writer = FeatureWriter(path2)
    schema = writer.arrow_schema
    batch = pa.RecordBatch.from_pylist([make_feature_record()], schema=schema)
    table = pa.Table.from_batches([batch])
    with writer:
        writer.write_table(table)
    assert parquet_row_count(path2) == 1


def test_writers_support_de_novo_records_without_database_fields(tmp_path):
    """PSM and feature writers accept records from a database-free workflow."""
    psm_record = make_psm_record(sequence="DENOVO", peptidoform="DENOVO")
    psm_record["protein_accessions"] = None

    psm_path = tmp_path / "denovo.psm.parquet"
    with PsmWriter(psm_path) as writer:
        writer.write_batch([psm_record])
    psm_table = pq.read_table(psm_path)
    assert psm_table.column("is_decoy").to_pylist() == [False]
    assert psm_table.column("protein_accessions").to_pylist() == [None]
    assert PsmSchema.validate(psm_table) == []

    feature_record = make_feature_record(sequence="DENOVO", peptidoform="DENOVO")
    feature_record["anchor_protein"] = None
    feature_record["pg_accessions"] = None
    feature_record["unique"] = None

    feature_path = tmp_path / "denovo.feature.parquet"
    with FeatureWriter(feature_path) as writer:
        writer.write_batch([feature_record])
    feature_table = pq.read_table(feature_path)
    assert feature_table.column("is_decoy").to_pylist() == [False]
    assert feature_table.column("anchor_protein").to_pylist() == [None]
    assert FeatureSchema.validate(feature_table) == []


def test_writer_batching(tmp_path):
    """Batch flush, remaining buffer flush on close, and multiple write_batch calls."""
    # Batch size flush
    path = tmp_path / "batch.feature.parquet"
    with FeatureWriter(path, batch_size=2) as w:
        records = [make_feature_record(sequence=f"SEQ{i}", peptidoform=f"SEQ{i}") for i in range(5)]
        w.write_batch(records)
    assert parquet_row_count(path) == 5

    # Remaining buffer flushed on close
    path2 = tmp_path / "buffer.feature.parquet"
    with FeatureWriter(path2, batch_size=100) as w:
        records = [make_feature_record(sequence=f"SEQ{i}", peptidoform=f"SEQ{i}") for i in range(3)]
        w.write_batch(records)
    assert parquet_row_count(path2) == 3

    # Multiple write_batch calls
    path3 = tmp_path / "multi.feature.parquet"
    with FeatureWriter(path3) as w:
        w.write_batch([make_feature_record(sequence="SEQ1", peptidoform="SEQ1")])
        w.write_batch([make_feature_record(sequence="SEQ2", peptidoform="SEQ2")])
        w.write_batch([make_feature_record(sequence="SEQ3", peptidoform="SEQ3")])
    assert parquet_row_count(path3) == 3


def test_writer_discards_buffer_when_context_body_fails(tmp_path):
    """A body error must not flush incomplete buffered records."""
    path = tmp_path / "failed.feature.parquet"

    with pytest.raises(RuntimeError, match="body failure"):
        with FeatureWriter(path, batch_size=2) as writer:
            writer.write_batch([make_feature_record()])
            raise RuntimeError("body failure")

    assert not path.exists()


def test_writer_compression(tmp_path):
    """Default compression is zstd; custom compression is recorded."""
    path = tmp_path / "default.feature.parquet"
    with FeatureWriter(path) as w:
        w.write_batch([make_feature_record()])
    assert read_parquet_metadata(path)["compression_format"] == "zstd"

    path2 = tmp_path / "snappy.feature.parquet"
    with FeatureWriter(path2, compression="snappy") as w:
        w.write_batch([make_feature_record()])
    assert read_parquet_metadata(path2)["compression_format"] == "snappy"


def test_pg_write_batch_warns_on_run_double_count(tmp_path, caplog):
    """bigbio/qpx#242: the streaming write_batch path must run the pg
    run-disjointness / referential check at close() (as a warning), not skip it.

    Two pg rows sharing the same (pg_accessions, label) with an overlapping
    grouped_runs double-count that run's intensity.
    """
    import logging

    path = tmp_path / "dup.pg.parquet"
    rec_a = make_pg_record(run_file_name="run_01")
    rec_b = make_pg_record(run_file_name="run_01")  # same members + label + run
    with caplog.at_level(logging.WARNING, logger="qpx.writers.base"):
        with PgWriter(path) as w:
            w.write_batch([rec_a, rec_b])
    assert path.exists()
    assert any("run_double_count" in rec.message for rec in caplog.records), (
        "streaming write_batch should warn on cross-row run double-count at close"
    )


def test_pg_write_batch_clean_has_no_referential_warning(tmp_path, caplog):
    """Disjoint grouped_runs across rows must NOT trigger the referential warning."""
    import logging

    path = tmp_path / "ok.pg.parquet"
    rec_a = make_pg_record(run_file_name="run_01")
    rec_b = make_pg_record(run_file_name="run_02")  # disjoint runs
    with caplog.at_level(logging.WARNING, logger="qpx.writers.base"):
        with PgWriter(path) as w:
            w.write_batch([rec_a, rec_b])
    assert not any("run_double_count" in rec.message for rec in caplog.records)


def test_pg_write_dataframe_explodes_intensities(tmp_path):
    """bigbio/qpx#252: PgWriter.write_dataframe must explode the natural
    ``intensities`` list into flat label/intensity rows instead of silently
    dropping quant (the base write_dataframe builds against the flat schema,
    which has no ``intensities`` column)."""
    import pandas as pd

    rec = make_pg_record(
        intensities=[
            {"label": "TMT126", "intensity": 5000.0},
            {"label": "TMT127N", "intensity": 6000.0},
        ]
    )
    df = pd.DataFrame([rec])
    path = tmp_path / "df.pg.parquet"
    with PgWriter(path) as w:
        w.write_dataframe(df)

    table = pq.read_table(path)
    assert table.num_rows == 2, "intensities list should explode into one row per label"
    assert set(table.column("label").to_pylist()) == {"TMT126", "TMT127N"}
    assert set(table.column("intensity").to_pylist()) == {5000.0, 6000.0}


class TestAtomicWrites:
    """A failed conversion must not destroy an existing valid file.

    The writer used to open the destination directly, so a converter that died
    mid-run left a truncated file where a good one had been (bigbio/qpx#288).
    This is not hypothetical: a 5,798-run DIA conversion was killed three times
    mid-write, twice by a full filesystem.
    """

    # batch_size defaults to 1,000,000, so a single buffered record never reaches
    # disk and nothing is at risk. Flush on every record so the writer really has
    # an open file when the failure lands — that is the situation being tested.
    BATCH = 1

    @classmethod
    def _write_one(cls, path):
        with FeatureWriter(str(path), batch_size=cls.BATCH) as w:
            w.write_batch([make_feature_record()])

    def test_existing_file_survives_a_failed_write(self, tmp_path):
        """A failed replacement leaves the existing destination unchanged."""
        path = tmp_path / "x.feature.parquet"
        self._write_one(path)
        good_bytes = path.read_bytes()

        with pytest.raises(RuntimeError):
            with FeatureWriter(str(path), batch_size=self.BATCH) as w:
                w.write_batch([make_feature_record()])
                raise RuntimeError("converter blew up mid-run")

        assert path.read_bytes() == good_bytes, "the previous valid file was modified"
        assert parquet_row_count(str(path)) == 1

    def test_no_staging_files_are_left_behind(self, tmp_path):
        """A failed write removes its temporary staging file."""
        path = tmp_path / "y.feature.parquet"
        self._write_one(path)
        with pytest.raises(RuntimeError):
            with FeatureWriter(str(path), batch_size=self.BATCH) as w:
                w.write_batch([make_feature_record()])
                raise RuntimeError("boom")

        leftovers = [p.name for p in tmp_path.iterdir() if p.name != path.name]
        assert leftovers == [], f"staging files left behind: {leftovers}"

    def test_successful_write_publishes_to_the_destination(self, tmp_path):
        """A successful close publishes only the requested destination."""
        path = tmp_path / "z.feature.parquet"
        self._write_one(path)
        assert path.exists()
        assert parquet_row_count(str(path)) == 1
        assert [p.name for p in tmp_path.iterdir()] == [path.name]

    def test_failed_write_leaves_no_file_when_there_was_none(self, tmp_path):
        """A failed first write does not publish a partial destination."""
        path = tmp_path / "fresh.feature.parquet"
        with pytest.raises(RuntimeError):
            with FeatureWriter(str(path), batch_size=self.BATCH) as w:
                w.write_batch([make_feature_record()])
                raise RuntimeError("boom")
        assert not path.exists(), "a partial file was published for a failed run"

    def test_close_after_a_caught_write_error_does_not_publish(self, tmp_path):
        """A converter that swallows its own write error must not publish.

        __exit__ discards when it sees an exception, but a caller that catches the
        error itself and then calls close() normally reached the promote path —
        publishing the truncated file the staging exists to prevent.
        """
        path = tmp_path / "caught.feature.parquet"
        self._write_one(path)
        good_bytes = path.read_bytes()

        writer = FeatureWriter(str(path), batch_size=self.BATCH)
        writer.write_batch([make_feature_record()])
        with pytest.raises(ValueError):
            writer.write_batch([{"not": "a valid record"}])
        writer.close()

        assert path.read_bytes() == good_bytes, "a partial file was published after a caught error"
        assert [p.name for p in tmp_path.iterdir()] == [path.name]

    def test_context_after_a_caught_write_error_does_not_publish(self, tmp_path):
        """A caught write error must still prevent context-manager promotion."""
        path = tmp_path / "caught-context.feature.parquet"
        self._write_one(path)
        good_bytes = path.read_bytes()

        with FeatureWriter(str(path), batch_size=self.BATCH) as writer:
            writer.write_batch([make_feature_record(sequence="PARTIAL", peptidoform="PARTIAL")])
            with pytest.raises(ValueError):
                writer.write_batch([{"not": "a valid record"}])

        assert path.read_bytes() == good_bytes, "context exit published a partial file after a caught error"
        assert [p.name for p in tmp_path.iterdir()] == [path.name]

    def test_dataframe_conversion_error_latches_failure(self, tmp_path):
        """A caught DataFrame conversion error must preserve the destination."""
        import pandas as pd

        path = tmp_path / "caught-dataframe.feature.parquet"
        self._write_one(path)
        good_bytes = path.read_bytes()

        writer = FeatureWriter(str(path), batch_size=self.BATCH)
        writer.write_batch([make_feature_record(sequence="PARTIAL", peptidoform="PARTIAL")])
        with pytest.raises(KeyError):
            writer.write_dataframe(pd.DataFrame([{"not": "a valid record"}]))
        writer.close()

        assert path.read_bytes() == good_bytes, "DataFrame failure allowed a partial file to be published"
        assert [p.name for p in tmp_path.iterdir()] == [path.name]

    def test_pg_preprocessing_error_latches_failure(self, tmp_path):
        """PgWriter preprocessing runs inside the same atomic-write guard."""
        path = tmp_path / "caught-preprocess.pg.parquet"
        with PgWriter(str(path), batch_size=self.BATCH) as writer:
            writer.write_batch([make_pg_record()])
        good_bytes = path.read_bytes()

        bad_record = make_pg_record()
        bad_record["additional_intensities"] = 1
        with PgWriter(str(path), batch_size=self.BATCH) as writer:
            writer.write_batch([make_pg_record(anchor_protein="PARTIAL")])
            with pytest.raises(TypeError):
                writer.write_batch([bad_record])

        assert path.read_bytes() == good_bytes, "PgWriter preprocessing failure published a partial file"
        assert [p.name for p in tmp_path.iterdir()] == [path.name]

    def test_final_buffer_failure_discards_staging(self, tmp_path):
        """A failure while flushing the final buffer leaves no partial file."""
        path = tmp_path / "final-buffer.feature.parquet"
        self._write_one(path)
        good_bytes = path.read_bytes()

        writer = FeatureWriter(str(path), batch_size=2)
        writer.write_batch(
            [
                make_feature_record(sequence="FIRST", peptidoform="FIRST"),
                make_feature_record(sequence="SECOND", peptidoform="SECOND"),
            ]
        )
        writer.write_batch([{"not": "a valid record"}])

        with pytest.raises(ValueError):
            writer.close()

        assert path.read_bytes() == good_bytes
        assert [p.name for p in tmp_path.iterdir()] == [path.name]

    @pytest.mark.parametrize("failure_target", ["_validate_identity_uniqueness", "_promote_staged"])
    def test_finalization_failure_discards_staging(self, tmp_path, monkeypatch, failure_target):
        """Validation and promotion failures both clean up staged output."""
        path = tmp_path / "finalization.feature.parquet"
        self._write_one(path)
        good_bytes = path.read_bytes()
        writer = FeatureWriter(str(path), batch_size=self.BATCH)
        writer.write_batch([make_feature_record(sequence="PARTIAL", peptidoform="PARTIAL")])

        def fail_finalization():
            """Simulate a failure after the staged Parquet file is closed."""
            raise RuntimeError("finalization failed")

        monkeypatch.setattr(writer, failure_target, fail_finalization)
        with pytest.raises(RuntimeError, match="finalization failed"):
            writer.close()

        assert path.read_bytes() == good_bytes
        assert [p.name for p in tmp_path.iterdir()] == [path.name]

    def test_two_writers_in_one_process_do_not_share_a_staging_file(self, tmp_path):
        """The staging name must be unique per writer, not per process."""
        path = tmp_path / "shared.feature.parquet"
        with FeatureWriter(str(path), batch_size=self.BATCH) as first, FeatureWriter(str(path), batch_size=self.BATCH) as second:
            first.write_batch([make_feature_record(sequence="FIRST", peptidoform="FIRST")])
            second.write_batch([make_feature_record(sequence="SECOND", peptidoform="SECOND")])
            staging_files = list(tmp_path.glob(f".{path.name}.*.tmp"))
            assert len(staging_files) == 2, "two writers share a staging file"
