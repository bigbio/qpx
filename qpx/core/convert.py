"""QueryResult — thin wrapper over DuckDB query results for materialization."""

import pandas as pd
import pyarrow as pa

# Rows fetched per batch when a result is iterated.
_BATCH_ROWS = 10_000


class QueryResult:
    """Read pandas, Arrow and row output through one DuckDB relation."""

    def __init__(self, connection, sql: str):
        self._connection = connection
        self._sql = sql
        # sql() runs a statement that is not a query at once; DuckDB 1.5 then returns None.
        self._relation = connection.sql(sql)
        self._rows = None

    def filter(self, condition: str) -> "QueryResult":
        """Keep the rows matching a SQL *condition*; DuckDB pushes it into the scans of the query."""
        return QueryResult(self._connection, f"SELECT * FROM ({self._sql}) WHERE {condition}")

    def to_df(self, dtype_backend: str | None = None) -> pd.DataFrame:
        """Materialize as a pandas DataFrame.

        A nested column becomes Python objects by default, about 8 KB per TMT16
        row. ``dtype_backend="pyarrow"`` keeps columns in Arrow memory instead.
        """
        if self._relation is None:
            return pd.DataFrame()
        if dtype_backend == "pyarrow":
            return self._relation.to_arrow_table().to_pandas(types_mapper=pd.ArrowDtype)
        if dtype_backend is not None:
            raise ValueError(f"Unknown dtype_backend {dtype_backend!r}; use 'pyarrow'")
        return self._relation.df()

    def to_arrow(self) -> pa.Table:
        """Materialize as an Arrow Table."""
        if self._relation is None:
            return pa.table({})
        return self._relation.to_arrow_table()

    def to_polars(self):
        """Materialize as a Polars DataFrame."""
        import polars as pl

        return pl.from_arrow(self.to_arrow())

    def _row_iter(self):
        if self._rows is None:
            if self._relation is None:
                self._rows = iter([])
            elif self._relation.type == "QUERY_RELATION":
                self._rows = self._stream_rows()
            else:
                # A statement DuckDB has already run (DDL and DML on 2.0) holds its result.
                self._rows = iter(self._relation.fetchall())
        return self._rows

    def _stream_rows(self):
        # Rows stream in batches through a cursor of their own: reading part of them
        # neither loads the rest nor blocks the connection, on which DuckDB 2.0
        # refuses relation work while a result is unread.
        cursor = self._connection.cursor()
        try:
            cursor.execute(self._sql)
            while batch := cursor.fetchmany(_BATCH_ROWS):
                yield from batch
        finally:
            cursor.close()

    def fetchone(self):
        """Fetch a single row as a tuple."""
        return next(self._row_iter(), None)

    def fetchall(self):
        """Fetch all rows as a list of tuples."""
        return list(self._row_iter())

    def __iter__(self):
        """Iterate over rows."""
        return self._row_iter()
