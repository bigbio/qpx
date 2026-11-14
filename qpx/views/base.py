"""BaseView — common logic for computed views."""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
from typing import TYPE_CHECKING

import pyarrow as pa
import pyarrow.parquet as pq

if TYPE_CHECKING:
    from qpx.dataset import Dataset

# Maximum number of cached query results per view instance.
_MAX_CACHE_ENTRIES = 16


class CachedQueryResult:
    """A pre-materialized query result that supports repeated access.

    QueryResult reruns its query on each materialization; this class keeps the Arrow table in memory.
    """

    def __init__(self, arrow_table: pa.Table):
        self._table = arrow_table

    def to_df(self):
        return self._table.to_pandas()

    def to_arrow(self) -> pa.Table:
        return self._table

    def to_polars(self):
        import polars as pl

        return pl.from_arrow(self._table)


class BaseView:
    """Base class for computed views over QPX datasets.

    Provides:
        - Shared DuckDB engine reference from the parent Dataset
        - LRU-bounded result caching (materialized Arrow tables)
        - Convenience materialization helpers (to_parquet, to_df, to_arrow)
    """

    def __init__(self, dataset: Dataset):
        self._engine = dataset._engine
        self._structures = set(dataset.available_structures)
        self._cache: OrderedDict[str, CachedQueryResult] = OrderedDict()

    def _require(self, *structures: str) -> None:
        """Raise a ValueError naming the structures the dataset lacks for this query."""
        missing = [name for name in structures if name not in self._structures]
        if missing:
            raise ValueError(f"{type(self).__name__} requires the {' and '.join(missing)} structure(s)")

    def _execute_cached(self, key: str, sql: str) -> CachedQueryResult:
        """Execute SQL with LRU caching. Returns a reusable CachedQueryResult."""
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        table = self._engine.connection.sql(sql).to_arrow_table()
        self._cache[key] = CachedQueryResult(table)
        if len(self._cache) > _MAX_CACHE_ENTRIES:
            self._cache.popitem(last=False)
        return self._cache[key]

    def invalidate_cache(self):
        """Clear all cached results."""
        self._cache.clear()

    @staticmethod
    def to_parquet(result, path: str | Path) -> Path:
        """Write a query result to Parquet."""
        table = result.to_arrow()
        pq.write_table(table, str(path), compression="zstd")
        return Path(path)

    @staticmethod
    def to_df(result):
        """Materialize a query result as a pandas DataFrame."""
        return result.to_df()

    @staticmethod
    def to_arrow(result):
        """Materialize a query result as an Arrow Table."""
        return result.to_arrow()
