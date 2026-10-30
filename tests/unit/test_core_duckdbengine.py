"""Tests for QPX core layer: DuckDBEngine, LazyQuery, QueryResult, Parquet I/O."""

import os

import pandas as pd
import pyarrow as pa
import pytest

from qpx.core.convert import QueryResult
from qpx.core.engine import DuckDBEngine, create_connection, create_engine, duckdb_defaults
from qpx.core.parquet_io import (
    parquet_row_count,
    read_parquet_metadata,
    read_parquet_schema,
)
from qpx.core.query import LazyQuery


def test_duckdb_engine(feature_parquet, tmp_path):
    """Engine: create, custom config, register parquet, context manager, params."""
    # Basic creation
    engine = DuckDBEngine()
    try:
        result = engine.execute("SELECT 42 AS answer")
        assert result.fetchone()[0] == 42
    finally:
        engine.close()

    # Custom memory and threads
    engine = DuckDBEngine(memory_limit="1GB", threads=2)
    try:
        result = engine.execute("SELECT current_setting('threads') AS t")
        assert int(result.fetchone()[0]) == 2
    finally:
        engine.close()

    # Register parquet creates view
    engine = DuckDBEngine()
    try:
        engine.register_parquet("feature", feature_parquet)
        count = engine.execute("SELECT COUNT(*) FROM feature").fetchone()[0]
        assert count == 3
    finally:
        engine.close()

    # Register replaces existing
    from qpx.writers import FeatureWriter
    from tests.conftest import make_feature_record

    path2 = tmp_path / "other.feature.parquet"
    with FeatureWriter(path2) as w:
        w.write_batch([make_feature_record()])

    engine = DuckDBEngine()
    try:
        engine.register_parquet("feature", feature_parquet)
        assert engine.execute("SELECT COUNT(*) FROM feature").fetchone()[0] == 3
        engine.register_parquet("feature", path2)
        assert engine.execute("SELECT COUNT(*) FROM feature").fetchone()[0] == 1
    finally:
        engine.close()

    # Context manager
    with DuckDBEngine() as engine:
        assert engine.execute("SELECT 1").fetchone()[0] == 1

    # Factory function
    engine = create_engine(memory_limit="2GB", threads=1)
    try:
        assert isinstance(engine, DuckDBEngine)
        assert engine.execute("SELECT 'ok'").fetchone()[0] == "ok"
    finally:
        engine.close()

    # Parameterized queries
    engine = DuckDBEngine()
    try:
        assert engine.execute("SELECT ? + ? AS sum", [10, 20]).fetchone()[0] == 30
    finally:
        engine.close()


@pytest.mark.skipif(not hasattr(os, "sched_getaffinity"), reason="no CPU affinity API")
def test_default_threads_follow_cpu_affinity(monkeypatch):
    """Without a thread count, DuckDB is kept to the CPUs the process may run on."""
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {0})
    engine = DuckDBEngine()
    try:
        assert engine.execute("SELECT current_setting('threads')").fetchall()[0][0] == 1
    finally:
        engine.close()


def test_command_resources_reach_connections_opened_in_the_block(tmp_path):
    """duckdb_defaults sets memory, threads and the spill folder of connections without their own, until the block ends."""
    spill = tmp_path / "spill dir"
    settings = "SELECT current_setting('memory_limit'), current_setting('threads'), current_setting('temp_directory')"
    with duckdb_defaults(memory_limit="1GiB", threads=1, temp_directory=spill):
        engine = DuckDBEngine()
        own = create_connection(threads=2)
    try:
        assert engine.execute(settings).fetchall()[0] == ("1.0 GiB", 1, str(spill))
        assert own.execute("SELECT current_setting('threads')").fetchall()[0][0] == 2
    finally:
        engine.close()
        own.close()
    after = create_connection()
    try:
        assert after.execute("SELECT current_setting('temp_directory')").fetchall()[0][0] == ".tmp"
    finally:
        after.close()


def test_lazy_query(feature_parquet):
    """LazyQuery: count, filter, execute, distinct_values."""
    engine = DuckDBEngine()
    engine.register_parquet("feature", feature_parquet)
    try:
        q = LazyQuery(engine, "feature")
        assert q.count() == 3
        assert q.filter("is_decoy = false").count() == 2

        rows = q.execute().fetchall()
        assert len(rows) == 3

        runs = q.distinct_values("run_file_name")
        assert set(runs) == {"run_01", "run_02"}
    finally:
        engine.close()


@pytest.mark.parametrize(
    ("list_column", "value_column"),
    [
        ("grouped_runs); DROP TABLE run; --", "run_file_name"),
        ("grouped_runs", "run_file_name); DROP TABLE pg; --"),
    ],
)
def test_join_list_membership_validates_column_names(list_column, value_column):
    """List-membership joins reject SQL passed as a column identifier."""
    with DuckDBEngine() as engine:
        left = LazyQuery(engine, "pg")
        right = LazyQuery(engine, "run")

        with pytest.raises(ValueError, match="Unsafe SQL identifier"):
            left.join_list_membership(
                right,
                list_column=list_column,
                value_column=value_column,
            )


def test_query_result(feature_parquet):
    """QueryResult: to_df, to_arrow, fetchone, fetchall, iter."""
    engine = DuckDBEngine()
    engine.register_parquet("feature", feature_parquet)
    try:
        # to_df
        result = QueryResult(engine.connection, "SELECT sequence, charge FROM feature")
        df = result.to_df()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert "sequence" in df.columns

        # to_arrow
        result = QueryResult(engine.connection, "SELECT sequence, charge FROM feature")
        table = result.to_arrow()
        assert isinstance(table, pa.Table)
        assert table.num_rows == 3
        assert result.to_arrow().num_rows == table.num_rows

        # fetchone
        assert QueryResult(engine.connection, "SELECT COUNT(*) FROM feature").fetchone()[0] == 3

        # fetchall
        rows = QueryResult(engine.connection, "SELECT sequence FROM feature ORDER BY sequence").fetchall()
        assert len(rows) == 3
        result = QueryResult(engine.connection, "SELECT sequence FROM feature ORDER BY sequence")
        assert result.fetchone() == rows[0]
        assert result.fetchall() == rows[1:]

        # iter
        assert len(list(QueryResult(engine.connection, "SELECT sequence FROM feature"))) == 3

        # statements run once
        QueryResult(engine.connection, "CREATE TABLE t AS SELECT 1 AS x").fetchall()
        QueryResult(engine.connection, "INSERT INTO t VALUES (2)").to_df()
        assert QueryResult(engine.connection, "SELECT COUNT(*) FROM t").fetchone()[0] == 2

        # a partly read result does not block later reads
        partial = QueryResult(engine.connection, "SELECT sequence FROM feature")
        partial.fetchone()
        assert QueryResult(engine.connection, "SELECT sequence FROM feature").to_arrow().num_rows == 3
        assert len(partial.fetchall()) == 2
    finally:
        engine.close()


def test_parquet_io(feature_parquet, dataset_parquet):
    """Parquet I/O: read_metadata, read_schema, row_count."""
    meta = read_parquet_metadata(feature_parquet)
    assert isinstance(meta, dict)
    assert meta["file_type"] == "feature_file"
    assert "qpx_version" in meta
    assert "creator" in meta
    assert "uuid" in meta
    assert len(meta["uuid"]) > 0

    schema = read_parquet_schema(feature_parquet)
    assert isinstance(schema, pa.Schema)
    assert "sequence" in schema.names

    assert parquet_row_count(feature_parquet) == 3
    assert parquet_row_count(dataset_parquet) == 1
