"""DuckDB connection management, configuration, and resource settings."""

from __future__ import annotations

import contextlib
import logging
import os
import re
from pathlib import Path

import duckdb

from qpx.core.sql import escape_path, sql_build, validate_identifier, validate_table

logger = logging.getLogger(__name__)

_SAFE_SET_VALUE = re.compile(r"^[\w./:@\-]+$")

# Columns that became list<int32> after files were written. Older files are read
# in the current shape: mz.scan was a single int32 before it followed the
# PSM/feature scan convention.
_SCALAR_TO_LIST_COLUMNS: dict[str, tuple[str, ...]] = {"mz": ("scan",)}


def _validate_set_value(value: str, name: str) -> str:
    """Validate a value used in a DuckDB SET statement."""
    if not _SAFE_SET_VALUE.match(value):
        raise ValueError(f"Invalid character in DuckDB config '{name}': {value!r}")
    return value


# Resource settings of the running command (see duckdb_defaults), used by every
# connection that does not set its own.
_command_resources: dict[str, str | int | Path | None] = {}


@contextlib.contextmanager
def duckdb_defaults(
    memory_limit: str | None = None,
    threads: int | None = None,
    temp_directory: str | Path | None = None,
):
    """Apply these resource settings to the DuckDB connections opened inside the block that do not set their own.

    The settings are process-wide for the duration of the block, so they also reach
    connections that a command opens in helpers and post-processing steps.
    """
    previous = dict(_command_resources)
    _command_resources.update(memory_limit=memory_limit, threads=threads, temp_directory=temp_directory)
    try:
        yield
    finally:
        _command_resources.clear()
        _command_resources.update(previous)


def _set_resources(
    conn: duckdb.DuckDBPyConnection,
    memory_limit: str | None,
    threads: int | None,
    temp_directory: str | Path | None = None,
) -> None:
    """Apply the given memory limit, thread count and spill directory; ``None`` keeps the command's or DuckDB's default.

    DuckDB defaults to 80% of RAM and one thread per CPU, both capped by a cgroup
    limit such as a container's CPU or memory quota, and spills an in-memory
    database to ``.tmp`` in the working directory.
    """
    memory_limit = memory_limit or _command_resources.get("memory_limit")
    threads = threads or _command_resources.get("threads")
    temp_directory = temp_directory or _command_resources.get("temp_directory")
    if temp_directory:
        conn.execute(sql_build("SET temp_directory='$path'", path=escape_path(str(temp_directory))))
    if memory_limit is not None:
        conn.execute(
            sql_build(
                "SET memory_limit='$val'",
                val=_validate_set_value(memory_limit, "memory_limit"),
            )
        )
    if threads is None and hasattr(os, "sched_getaffinity"):
        # DuckDB counts every CPU of the host and ignores the affinity mask that
        # taskset or a SLURM/LSF cpuset sets, so it would oversubscribe those CPUs.
        default = conn.execute("SELECT current_setting('threads')").fetchall()[0][0]
        threads = min(default, len(os.sched_getaffinity(0)))
    if threads is not None:
        conn.execute(sql_build("SET threads=$val", val=str(int(threads))))


class DuckDBEngine:
    """
    Manages a DuckDB connection with configurable resources.

    Used by Dataset (shared connection) and standalone structures
    (own connection). All DuckDB tuning lives here.
    """

    def __init__(
        self,
        memory_limit: str | None = None,
        threads: int | None = None,
        temp_directory: str | None = None,
        s3_config: dict | None = None,
    ):
        self._conn = create_connection(memory_limit, threads, temp_directory)
        self._conn.execute("SET enable_progress_bar=true")
        if s3_config is not None:
            self._setup_s3(s3_config)

    def _setup_s3(self, config: dict) -> None:
        """Configure DuckDB httpfs for S3 access.

        The settings are GLOBAL so that cursors of the connection, which stream
        query rows, read S3 too; the in-memory database is this engine's own.
        """
        self._conn.execute("INSTALL httpfs")
        self._conn.execute("LOAD httpfs")
        if "region" in config:
            self._conn.execute(
                sql_build(
                    "SET GLOBAL s3_region='$val'",
                    val=_validate_set_value(config["region"], "s3_region"),
                )
            )
        if "access_key_id" in config and "secret_access_key" in config:
            self._conn.execute(
                sql_build(
                    "SET GLOBAL s3_access_key_id='$val'",
                    val=_validate_set_value(config["access_key_id"], "s3_access_key_id"),
                )
            )
            self._conn.execute(
                sql_build(
                    "SET GLOBAL s3_secret_access_key='$val'",
                    val=_validate_set_value(config["secret_access_key"], "s3_secret_access_key"),
                )
            )
        if "endpoint" in config:
            self._conn.execute(
                sql_build(
                    "SET GLOBAL s3_endpoint='$val'",
                    val=_validate_set_value(config["endpoint"], "s3_endpoint"),
                )
            )
        if config.get("anonymous", False):
            self._conn.execute("SET GLOBAL s3_url_style='path'")

    @property
    def connection(self) -> duckdb.DuckDBPyConnection:
        return self._conn

    def _create_view(self, name: str, source: str) -> None:
        """Create view *name* over a ``read_parquet`` *source*, reading legacy scalar columns as lists."""
        view = validate_table(name)
        self._conn.execute(sql_build("CREATE OR REPLACE VIEW $view AS SELECT * FROM $source", view=view, source=source))
        columns = _SCALAR_TO_LIST_COLUMNS.get(view, ())
        if not columns:
            return
        types = dict(self._conn.execute(sql_build("SELECT column_name, column_type FROM (DESCRIBE $view)", view=view)).fetchall())
        legacy = [column for column in columns if column in types and not types[column].endswith("]")]
        if legacy:
            replace = ", ".join(
                sql_build("CASE WHEN $column IS NULL THEN NULL ELSE [CAST($column AS INTEGER)] END AS $column", column=column)
                for column in map(validate_identifier, legacy)
            )
            self._conn.execute(
                sql_build(
                    "CREATE OR REPLACE VIEW $view AS SELECT * REPLACE ($replace) FROM $source",
                    view=view,
                    replace=replace,
                    source=source,
                )
            )

    def register_parquet(self, name: str, file_path: str | Path) -> None:
        """Register a Parquet file as a lazy DuckDB view."""
        self._create_view(name, sql_build("read_parquet('$path')", path=escape_path(str(file_path))))

    def register_parquet_files(self, name: str, file_paths) -> None:
        """Register several Parquet shards as one unioned DuckDB view.

        DuckDB unions the files by column name, matching the S3 glob path so a
        multi-shard structure yields the same rows locally and over S3
        (bigbio/qpx#252); a shard without a column reads as null for it.
        """
        paths = list(file_paths)
        if len(paths) == 1:
            self.register_parquet(name, paths[0])
            return
        path_list = ", ".join(f"'{escape_path(str(p))}'" for p in paths)
        self._create_view(name, sql_build("read_parquet([$paths], union_by_name=true)", paths=path_list))

    def register_partitioned_parquet(self, name: str, directory: str | Path) -> None:
        """Register a Hive-partitioned Parquet directory as a DuckDB view.

        A run name is a string however it looks: autocast, ``1/10`` would be a
        BIGINT and ``2024-01-15`` a DATE.
        """
        glob_pattern = str(Path(directory) / "**" / "*.parquet")
        hive_types = ", hive_types={'run_file_name': 'VARCHAR'}" if next(Path(directory).rglob("run_file_name=*"), None) else ""
        self._create_view(
            name,
            sql_build("read_parquet('$path', hive_partitioning=true$types)", path=escape_path(glob_pattern), types=hive_types),
        )

    def register_s3_parquet(self, name: str, s3_path: str) -> None:
        """Register an S3 Parquet file as a DuckDB view."""
        self._create_view(name, sql_build("read_parquet('$path', union_by_name=true)", path=escape_path(s3_path)))

    def register_view(self, name: str, sql: str) -> None:
        """Register ``sql`` as a named DuckDB view.

        ``name`` is validated as a table identifier before use. The view is created
        through ``execute`` rather than the relation API: DuckDB 2.0 refuses a
        relation while an earlier ``execute`` result on the connection is unread.
        """
        view = validate_table(name)
        self._conn.execute(sql_build("CREATE OR REPLACE VIEW $view AS $query", view=view, query=sql))

    def view_columns(self, view: str) -> set[str]:
        """Names of the columns of ``view``."""
        describe = sql_build("DESCRIBE $view", view=validate_table(view))
        return {row[0] for row in self._conn.execute(describe).fetchall()}

    def optional_columns(self, view: str, alias: str, names: tuple[str, ...]) -> str:
        """Select list for optional columns of ``view``: ``alias.name``, or ``NULL AS name`` when the view lacks it."""
        present = self.view_columns(view)
        return ", ".join(f"{alias}.{name}" if name in present else f"NULL AS {name}" for name in names)

    def execute(self, sql: str, params=None):
        """Run ``sql`` through the DB-API cursor.

        Read the result to the end: DuckDB 2.0 refuses relation work on the
        connection while an ``execute`` result is unread.
        """
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


def create_connection(
    memory_limit: str | None = None,
    threads: int | None = None,
    temp_directory: str | Path | None = None,
) -> duckdb.DuckDBPyConnection:
    """Create an in-memory DuckDB connection with validated resource settings."""
    conn = duckdb.connect(":memory:")
    _set_resources(conn, memory_limit, threads, temp_directory)
    return conn


def create_converter_connection(
    memory_limit: str | None = None,
    threads: int | None = None,
    temp_directory: str | None = None,
) -> duckdb.DuckDBPyConnection:
    """Create a validated DuckDB in-memory connection for converter use.

    Unifies DuckDB configuration across converters. Uses the same validation
    as DuckDBEngine to prevent injection via SET values.
    """
    conn = create_connection(memory_limit, threads, temp_directory)
    conn.execute("SET enable_progress_bar=true")
    conn.execute("SET preserve_insertion_order=false")
    return conn
