"""Hive partitioning tests."""

import pyarrow as pa
import pyarrow.parquet as pq


class TestPartitionedRegistration:
    def test_engine_reads_partitioned_dir(self, tmp_path):
        from qpx.core.engine import DuckDBEngine

        # Create Hive-partitioned structure
        for run in ["run_01.raw", "run_02.raw"]:
            run_dir = tmp_path / "feature" / f"run_file_name={run}"
            run_dir.mkdir(parents=True)
            table = pa.table({"sequence": ["PEPTIDEK"], "charge": [2]})
            pq.write_table(table, run_dir / "part-0.parquet")

        engine = DuckDBEngine()
        engine.register_partitioned_parquet("feature", tmp_path / "feature")
        result = engine.execute("SELECT COUNT(*) AS cnt FROM feature").fetchone()
        assert result[0] == 2

        # Verify the partition column is accessible
        result = engine.execute("SELECT DISTINCT run_file_name FROM feature ORDER BY run_file_name").fetchall()
        assert len(result) == 2
        engine.close()

    def test_dataset_discovers_partitioned_dir(self, dataset_dir):
        """Move feature into partitioned layout and verify Dataset discovers it."""
        import qpx

        feature_file = list(dataset_dir.glob("*.feature.parquet"))[0]
        table = pq.read_table(feature_file)
        original_rows = table.num_rows
        feature_file.unlink()

        # Write as partitioned
        from qpx.writers.base import BaseWriter

        part_dir = dataset_dir / "feature"
        BaseWriter.write_partitioned(table, part_dir, ["run_file_name"])

        ds = qpx.open_dataset(str(dataset_dir))
        assert "feature" in ds.available_structures
        assert ds.feature.count() == original_rows
        ds.close()

    def test_single_file_preferred_over_partitioned(self, dataset_dir):
        """When both a single file and partitioned dir exist, single file wins."""
        import qpx

        # Create a partitioned dir alongside the existing single file
        part_dir = dataset_dir / "feature"
        part_dir.mkdir()
        run_dir = part_dir / "run_file_name=fake_run"
        run_dir.mkdir()
        table = pa.table(
            {
                "sequence": ["FAKEPEP"],
                "charge": pa.array([1], type=pa.int16()),
            }
        )
        pq.write_table(table, run_dir / "part-0.parquet")

        ds = qpx.open_dataset(str(dataset_dir))
        # Should use the single file (3 rows), not the partitioned dir (1 row)
        assert ds.feature.count() == 3
        ds.close()


class TestPartitionedWriting:
    def test_write_partitioned_creates_hive_dirs(self, tmp_path):
        from qpx.writers.base import BaseWriter

        table = pa.table(
            {
                "run_file_name": ["run_01", "run_01", "run_02"],
                "sequence": ["PEP1", "PEP2", "PEP3"],
            }
        )
        out = tmp_path / "partitioned"
        BaseWriter.write_partitioned(table, out, ["run_file_name"])
        dirs = sorted(d.name for d in out.iterdir() if d.is_dir())
        assert "run_file_name=run_01" in dirs
        assert "run_file_name=run_02" in dirs

    def test_write_partitioned_pg_raises_actionable_error(self, tmp_path):
        """A 1.1 pg table (grouped_runs list, no run_file_name) is not Hive-
        partitionable — the default path must raise a clear ValueError, not a
        raw KeyError (regression for bigbio/qpx#220 review)."""
        import pytest

        from qpx.writers.base import BaseWriter

        pg = pa.table(
            {
                "anchor_protein": ["P1"],
                "grouped_runs": [["r1", "r2"]],
                "intensities": [[{"label": "LFQ", "intensity": 1.0}]],
            }
        )
        # Default path (would default to the missing run_file_name).
        with pytest.raises(ValueError, match="not Hive-partitionable"):
            BaseWriter.write_partitioned(pg, tmp_path / "pg")
        # Explicitly partitioning on the list column is also rejected clearly.
        with pytest.raises(ValueError, match="non-scalar|not partitionable"):
            BaseWriter.write_partitioned(pg, tmp_path / "pg2", ["grouped_runs"])
        # A missing scalar column also gives an actionable message, not KeyError.
        with pytest.raises(ValueError, match="missing column"):
            BaseWriter.write_partitioned(pg, tmp_path / "pg3", ["nope"])

    def test_write_partitioned_data_readable(self, tmp_path):
        from qpx.core.engine import DuckDBEngine
        from qpx.writers.base import BaseWriter

        table = pa.table(
            {
                "run_file_name": ["run_01", "run_01", "run_02"],
                "sequence": ["PEP1", "PEP2", "PEP3"],
                "charge": [2, 3, 2],
            }
        )
        out = tmp_path / "partitioned"
        BaseWriter.write_partitioned(table, out, ["run_file_name"])

        engine = DuckDBEngine()
        engine.register_partitioned_parquet("test", out)
        result = engine.execute("SELECT COUNT(*) FROM test").fetchone()
        assert result[0] == 3

        # Verify data integrity
        result = engine.execute("SELECT sequence FROM test WHERE run_file_name = 'run_02'").fetchone()
        assert result[0] == "PEP3"
        engine.close()

    def test_write_partitioned_default_partition_col(self, tmp_path):
        from qpx.writers.base import BaseWriter

        table = pa.table(
            {
                "run_file_name": ["run_01", "run_02"],
                "data": ["a", "b"],
            }
        )
        out = tmp_path / "default_part"
        BaseWriter.write_partitioned(table, out)  # No partition_cols arg
        dirs = sorted(d.name for d in out.iterdir() if d.is_dir())
        assert any("run_file_name=" in d for d in dirs)

    def test_write_partitioned_via_subclass_derives_id_and_validates(self, tmp_path, caplog):
        """bigbio/qpx#252: called on a concrete writer subclass, write_partitioned
        must derive the identity id (when absent) and run schema validation,
        instead of bypassing both. Called on the bare BaseWriter it stays a raw
        dump (covered by the other tests)."""
        import logging

        import pyarrow.dataset as pds

        from qpx.writers import FeatureWriter
        from tests.conftest import make_feature_record

        # Build a proper feature file, then strip the derived feature_id so the
        # partitioned write has to re-derive it.
        src = tmp_path / "src.feature.parquet"
        with FeatureWriter(src) as w:
            w.write_batch(
                [
                    make_feature_record(peptidoform="PEPTIDEK", run_file_name="run_01"),
                    make_feature_record(peptidoform="ELVISK", run_file_name="run_02"),
                ]
            )
        table = pq.read_table(src)
        assert "feature_id" in table.schema.names
        stripped = table.drop_columns(["feature_id"])

        out = tmp_path / "feat_part"
        with caplog.at_level(logging.WARNING, logger="qpx.writers.base"):
            FeatureWriter.write_partitioned(stripped, out, ["run_file_name"])

        back = pds.dataset(str(out), format="parquet", partitioning="hive").to_table()
        assert "feature_id" in back.schema.names
        ids = back.column("feature_id").to_pylist()
        assert all(i is not None for i in ids), "feature_id must be derived, not left null"
        assert len(set(ids)) == 2
