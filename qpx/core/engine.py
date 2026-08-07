"""DuckDB connection management, configuration, and resource settings."""

from __future__ import annotations

import logging
import re
from pathlib import Path

import duckdb

from qpx.core.sql import escape_path, sql_build, validate_table

logger = logging.getLogger(__name__)

_SAFE_SET_VALUE = re.compile(r"^[\w./:@\-]+$")


def _validate_set_value(value: str, name: str) -> str:
    """Validate a value used in a DuckDB SET statement."""
    if not _SAFE_SET_VALUE.match(value):
        raise ValueError(f"Invalid character in DuckDB config '{name}': {value!r}")
    return value


class DuckDBEngine:
    """
    Manages a DuckDB connection with configurable resources.

    Used by Dataset (shared connection) and standalone structures
    (own connection). All DuckDB tuning lives here.
    """

    def __init__(
        self,
        memory_limit: str = "16GB",
        threads: int = 4,
        temp_directory: str | None = None,
        s3_config: dict | None = None,
    ):
        self._conn = duckdb.connect(":memory:")
        self._conn.execute(
            sql_build(
                "SET memory_limit='$val'",
                val=_validate_set_value(memory_limit, "memory_limit"),
            )
        )
        self._conn.execute(sql_build("SET threads=$val", val=str(int(threads))))
        if temp_directory:
            self._conn.execute(
                sql_build(
                    "SET temp_directory='$val'",
                    val=_validate_set_value(temp_directory, "temp_directory"),
                )
            )
        self._conn.execute("SET enable_progress_bar=true")
        if s3_config is not None:
            self._setup_s3(s3_config)

    def _setup_s3(self, config: dict) -> None:
        """Configure DuckDB httpfs for S3 access."""
        self._conn.execute("INSTALL httpfs")
        self._conn.execute("LOAD httpfs")
        if "region" in config:
            self._conn.execute(
                sql_build(
                    "SET s3_region='$val'",
                    val=_validate_set_value(config["region"], "s3_region"),
                )
            )
        if "access_key_id" in config and "secret_access_key" in config:
            self._conn.execute(
                sql_build(
                    "SET s3_access_key_id='$val'",
                    val=_validate_set_value(config["access_key_id"], "s3_access_key_id"),
                )
            )
            self._conn.execute(
                sql_build(
                    "SET s3_secret_access_key='$val'",
                    val=_validate_set_value(config["secret_access_key"], "s3_secret_access_key"),
                )
            )
        if "endpoint" in config:
            self._conn.execute(
                sql_build(
                    "SET s3_endpoint='$val'",
                    val=_validate_set_value(config["endpoint"], "s3_endpoint"),
                )
            )
        if config.get("anonymous", False):
            self._conn.execute("SET s3_url_style='path'")

    @property
    def connection(self) -> duckdb.DuckDBPyConnection:
        return self._conn

    def register_parquet(self, name: str, file_path: str | Path) -> None:
        """Register a Parquet file as a lazy DuckDB view."""
        self._conn.execute(
            sql_build(
                "CREATE OR REPLACE VIEW $view AS SELECT * FROM read_parquet('$path')",
                view=validate_table(name),
                path=escape_path(str(file_path)),
            )
        )

    def register_parquet_files(self, name: str, file_paths) -> None:
        """Register several Parquet shards as one unioned DuckDB view.

        DuckDB reads the list of files positionally-unioned, matching the S3
        glob path's behaviour so a multi-shard structure yields the same rows
        locally and over S3 (bigbio/qpx#252).
        """
        paths = list(file_paths)
        if len(paths) == 1:
            self.register_parquet(name, paths[0])
            return
        path_list = ", ".join(f"'{escape_path(str(p))}'" for p in paths)
        self._conn.execute(
            sql_build(
                "CREATE OR REPLACE VIEW $view AS SELECT * FROM read_parquet([$paths])",
                view=validate_table(name),
                paths=path_list,
            )
        )

    def register_partitioned_parquet(self, name: str, directory: str | Path) -> None:
        """Register a Hive-partitioned Parquet directory as a DuckDB view."""
        glob_pattern = str(Path(directory) / "**" / "*.parquet")
        self._conn.execute(
            sql_build(
                "CREATE OR REPLACE VIEW $view AS SELECT * FROM read_parquet('$path', hive_partitioning=true)",
                view=validate_table(name),
                path=escape_path(glob_pattern),
            )
        )

    def register_s3_parquet(self, name: str, s3_path: str) -> None:
        """Register an S3 Parquet file as a DuckDB view."""
        self._conn.execute(
            sql_build(
                "CREATE OR REPLACE VIEW $view AS SELECT * FROM read_parquet('$path')",
                view=validate_table(name),
                path=escape_path(s3_path),
            )
        )

    def execute(self, sql: str, params=None):
        if params:
            return self._conn.execute(sql, params)
        return self._conn.execute(sql)

    def close(self):
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def create_engine(**kwargs) -> DuckDBEngine:
    """Factory for creating a DuckDB engine with default settings."""
    return DuckDBEngine(**kwargs)


def create_converter_connection(
    memory_limit: str = "16GB",
    threads: int = 4,
    temp_directory: str | None = None,
) -> duckdb.DuckDBPyConnection:
    """Create a validated DuckDB in-memory connection for converter use.

    Unifies DuckDB configuration across converters. Uses the same validation
    as DuckDBEngine to prevent injection via SET values.
    """
    conn = duckdb.connect(":memory:")
    conn.execute(
        sql_build(
            "SET memory_limit='$val'",
            val=_validate_set_value(memory_limit, "memory_limit"),
        )
    )
    conn.execute(sql_build("SET threads=$val", val=str(int(threads))))
    conn.execute("SET enable_progress_bar=true")
    conn.execute("SET preserve_insertion_order=false")
    if temp_directory:
        conn.execute(
            sql_build(
                "SET temp_directory='$val'",
                val=_validate_set_value(temp_directory, "temp_directory"),
            )
        )
    return conn
