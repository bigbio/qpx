"""BaseStructure — lazy DuckDB-backed data structure over a Parquet file."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from qpx.core.convert import QueryResult
from qpx.core.data.footer_validation import duplicate_composite_issues, footer_issues
from qpx.core.data.schema import ValidationResult, pg_quantity_issues, pg_referential_issues_from_query
from qpx.core.engine import DuckDBEngine, create_engine
from qpx.core.parquet_io import read_parquet_metadata
from qpx.core.query import LazyQuery, _escape_sql_string
from qpx.core.sql import sql_build, validate_identifier

# Rows per Arrow batch while validating: bounds the memory a validation holds.
_VALIDATE_BATCH_ROWS = 100_000


class BaseStructure:
    """
    Lazy DuckDB-backed data structure over a Parquet file.

    Thin wrapper on core.query.LazyQuery — domain-specific methods only.
    All query building, materialization, and I/O lives in qpx.core.
    """

    _schema_class = None  # Override in subclasses (e.g., FeatureSchema)

    def __init__(
        self,
        engine: DuckDBEngine,
        table_name: str,
        file_path: Path,
        file_paths=None,
        query: LazyQuery | None = None,
    ):
        self._engine = engine
        self._table_name = table_name
        self._file_path = file_path
        # Every file backing this structure. A sharded structure has more than
        # one, and ``_file_path`` names only the first (kept for provenance), so
        # anything that re-registers the data from a path MUST use this list or
        # it silently reads a subset (bigbio/qpx#286). Defaults to the single
        # file so producers that pass only ``file_path`` stay correct.
        # Keep locators exactly as given. Path() would rewrite an "s3://bucket/x"
        # URI to "s3:/bucket/x" and it would no longer resolve.
        self._file_paths = list(file_paths) if file_paths else [file_path]
        self._query = query if query is not None else LazyQuery(engine, table_name)

    @property
    def file_paths(self) -> list:
        """Every file backing this structure, in registration order.

        A sharded structure has more than one; ``_file_path`` names only the
        first. Anything re-registering the data from a path must use this, or it
        silently reads a subset (bigbio/qpx#286).
        """
        return list(self._file_paths)

    @classmethod
    def from_file(cls, path: str | Path, **engine_kwargs) -> BaseStructure:
        """Open a standalone Parquet file (creates its own DuckDB engine)."""
        path = Path(path)
        engine = create_engine(**engine_kwargs)
        table_name = cls._schema_class.view_name
        engine.register_parquet(table_name, path)
        return cls(engine=engine, table_name=table_name, file_path=path)

    # --- Schema ---
    @classmethod
    def schema(cls) -> pa.Schema:
        return cls._schema_class.get_arrow_schema()

    # --- Lazy query building (delegates to core.query.LazyQuery) ---
    def select(self, *columns: str) -> BaseStructure:
        return self._with_query(self._query.select(*columns))

    def filter(self, condition: str) -> BaseStructure:
        return self._with_query(self._query.filter(condition))

    def limit(self, n: int) -> BaseStructure:
        return self._with_query(self._query.limit(n))

    def join(self, other: BaseStructure, on: str | None = None) -> BaseStructure:
        if on is None:
            on = self._auto_join_key(other)
        return self._with_query(self._query.join(other.query, on=on))

    def _join_list_membership(
        self,
        other: BaseStructure,
        list_column: str,
        value_column: str,
    ) -> BaseStructure:
        """Join by expanding one of this structure's list columns."""
        return self._with_query(
            self._query.join_list_membership(
                other.query,
                list_column=list_column,
                value_column=value_column,
            )
        )

    def order_by(self, *columns: str, desc: bool = False) -> BaseStructure:
        return self._with_query(self._query.order_by(*columns, desc=desc))

    # --- Materialization (delegates to core.convert.QueryResult) ---
    def to_df(self, dtype_backend: str | None = None) -> pd.DataFrame:
        """Execute query and return pandas DataFrame; ``dtype_backend="pyarrow"`` keeps nested columns in Arrow memory."""
        return QueryResult(self._engine.connection, self._query.build_sql()).to_df(dtype_backend)

    def to_arrow(self) -> pa.Table:
        """Execute query and return Arrow Table."""
        return self._engine.connection.sql(self._query.build_sql()).to_arrow_table()

    def to_polars(self):
        """Execute query and return Polars DataFrame."""
        return QueryResult(self._engine.connection, self._query.build_sql()).to_polars()

    def count(self) -> int:
        return self._query.count()

    @property
    def query(self) -> LazyQuery:
        """Return the immutable query backing this structure."""
        return self._query

    # --- Batch iteration (for sample-level processing) ---
    def iter_batches(self, partition_by: str, batch_size: int = 20) -> Iterator[tuple[list[str], pd.DataFrame]]:
        """
        Iterate over data in batches, partitioned by a column.

        Yields (partition_values, DataFrame) tuples. Useful for
        memory-bounded processing of large datasets.
        """
        all_values = self._query.distinct_values(partition_by)
        for i in range(0, len(all_values), batch_size):
            batch_values = all_values[i : i + batch_size]
            placeholders = ", ".join(f"'{_escape_sql_string(str(v))}'" for v in batch_values)
            filtered = self._query.filter(f"{partition_by} IN ({placeholders})")
            df = QueryResult(self._engine.connection, filtered.build_sql()).to_df()
            yield batch_values, df

    # --- Validation ---
    def validate(self, *, strict: bool = False, semantic: bool = False) -> ValidationResult:
        """Validate this structure's data against its schema.

        ``strict=True`` (used by ``qpxc validate``) promotes duplicate-PK and
        required-null issues to errors; the default is lenient. Rows are read in
        batches; only the primary-key columns are held for the whole view.
        ``semantic=True`` also checks PSM and feature peptidoforms, modification
        positions and m/z values.
        """
        view = self._schema_class
        severity = "error" if strict else "warning"
        query = self._query.build_sql()
        schema = self._engine.connection.sql(sql_build("SELECT * FROM ($query) LIMIT 0", query=query)).to_arrow_table().schema
        # DuckDB drops footer metadata, so the footers are read from the local files.
        footers = [(str(path), pq.read_metadata(path).metadata or {}) for path in self._file_paths if Path(path).is_file()]
        # A single PSM file has one format; never apply the first shard's
        # declaration to a union of several files.
        if view.view_name == "psm" and len(self._file_paths) == 1 and footers:
            scan_format = footers[0][1].get(b"scan_format")
            if scan_format:
                schema = schema.with_metadata({**(schema.metadata or {}), b"scan_format": scan_format})
        # Whole-view queries run before the batches stream: a new query on the
        # connection would end the stream.
        keys = None
        if set(view.primary_key) <= set(schema.names):
            columns = ", ".join(validate_identifier(name) for name in view.primary_key)
            keys = self._engine.connection.sql(
                sql_build("SELECT $columns FROM ($query)", columns=columns, query=query)
            ).to_arrow_table()
        referential = footer_issues(view, footers, severity)
        referential += duplicate_composite_issues(self._engine, query, view, schema.names, footers, severity)
        if view.view_name == "pg":
            referential += self._pg_view_issues(query, set(schema.names), severity)
        # Only the columns the row checks read are streamed; the constant keeps the
        # select valid for a view that has none of them.
        columns = ", ".join(["1 AS _validation_row", *map(validate_identifier, view.validation_columns(schema))])
        batches = self._engine.connection.sql(
            sql_build("SELECT $columns FROM ($query)", columns=columns, query=query)
        ).to_arrow_reader(_VALIDATE_BATCH_ROWS)
        result = view.validate_batches(schema, batches, keys, referential, strict=strict)
        if semantic:
            from qpx.core.semantic_validation import semantic_issues

            result.issues.extend(semantic_issues(self._engine, query, schema.names, view.view_name, severity))
        return result

    def _pg_view_issues(self, query: str, names: set[str], severity: str) -> list:
        """Whole-view pg checks: referential invariants and the presence of a quantity."""
        issues = []
        if {"pg_accessions", "grouped_runs", "label"} <= names:
            issues += pg_referential_issues_from_query(self._engine, query, "pg", severity)
        if {"label", "intensity"} <= names:
            issues += pg_quantity_issues(self._engine, query, "pg")
        return issues

    # --- Parquet metadata (delegates to core.parquet_io) ---
    @property
    def file_metadata(self) -> dict[str, str]:
        return read_parquet_metadata(self._file_path)

    # --- Internals ---
    def _with_query(self, new_query: LazyQuery) -> BaseStructure:
        """Return a clone with a different LazyQuery (immutable pattern)."""
        return self.__class__(
            engine=self._engine,
            table_name=self._table_name,
            file_path=self._file_path,
            file_paths=self.file_paths,
            query=new_query,
        )

    def _auto_join_key(self, other: BaseStructure) -> str:
        my_cols = set(self.schema().names)
        their_cols = set(other.schema().names)
        common = (my_cols & their_cols) - {"is_decoy"}
        if "run_file_name" in common:
            return "run_file_name"
        if "sample_accession" in common:
            return "sample_accession"
        if common:
            return next(iter(common))
        raise ValueError("No common column found for auto-join. Specify 'on' explicitly.")

    def __repr__(self):
        return f"{self.__class__.__name__}('{self._file_path}', rows={self.count()})"

    def __len__(self):
        return self.count()
