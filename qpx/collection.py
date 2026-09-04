"""DatasetCollection — virtual and physical multi-dataset operations."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from qpx.core.convert import QueryResult
from qpx.core.engine import DuckDBEngine

if TYPE_CHECKING:
    from qpx.dataset import Dataset


class DatasetCollection:
    """Analyze multiple QPX datasets together.

    Virtual mode: register all datasets in one DuckDB engine for cross-dataset SQL.
    Physical mode: merge structures from multiple datasets into a new output directory.

    Usage::

        coll = DatasetCollection([ds1, ds2])
        coll.sql("SELECT COUNT(*) FROM feature_0 UNION ALL SELECT COUNT(*) FROM feature_1")

        # Merge into a new dataset
        coll.merge("merged_output/", structures=["feature"])
    """

    def __init__(
        self,
        datasets: list[Dataset],
        duckdb_memory: str = "32GB",
        duckdb_threads: int = 4,
    ):
        self.datasets = datasets
        self._engine = DuckDBEngine(
            memory_limit=duckdb_memory,
            threads=duckdb_threads,
        )
        self._register_all()

    def _register_all(self):
        """Register all dataset structures with indexed names (e.g. feature_0, feature_1)."""
        for i, ds in enumerate(self.datasets):
            for name in ds.available_structures:
                struct = ds._structures[name]
                indexed_name = f"{name}_{i}"
                if hasattr(struct, "_file_path") and struct._file_path:
                    # Register EVERY file backing the structure, not just
                    # ``_file_path``. That attribute names only the first shard;
                    # a sharded structure registered from it returns a subset of
                    # the rows the dataset itself returns, silently
                    # (bigbio/qpx#286). ``_file_paths`` falls back to the single
                    # file for non-sharded structures.
                    file_paths = struct.file_paths
                    first = str(file_paths[0])
                    if first.startswith("s3://"):
                        # S3 locators are globs, not filesystem paths: hand them
                        # to the S3 registration so the collection reads what the
                        # dataset reads.
                        self._engine.register_s3_parquet(indexed_name, first)
                    elif Path(first).is_dir():
                        # Hive-partitioned structures are a single directory.
                        self._engine.register_partitioned_parquet(indexed_name, file_paths[0])
                    else:
                        self._engine.register_parquet_files(indexed_name, file_paths)

    def sql(self, query: str) -> QueryResult:
        """Execute SQL across all registered datasets."""
        return QueryResult(self._engine.execute(query))

    @property
    def structure_names(self) -> dict[int, list[str]]:
        """Return available structures per dataset index."""
        return {i: ds.available_structures for i, ds in enumerate(self.datasets)}

    def merge(
        self,
        output_dir: str | Path,
        structures: list[str] | None = None,
        prefix: str = "merged",
    ) -> Path:
        """Physically merge datasets into a new directory.

        Concatenates matching structures across all datasets.
        Adds a ``source_dataset`` column to distinguish origins.
        """
        import pyarrow as pa
        import pyarrow.parquet as pq

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        target_structures = structures or self._common_structures()

        for name in target_structures:
            tables = []
            for ds in self.datasets:
                struct = ds._structures.get(name)
                if struct is None:
                    continue
                tbl = struct.to_arrow()
                source_col = pa.array([str(ds.path)] * tbl.num_rows, type=pa.string())
                tbl = tbl.append_column("source_dataset", source_col)
                tables.append(tbl)

            if tables:
                merged = pa.concat_tables(tables, promote_options="permissive")
                _, suffix = self.datasets[0]._STRUCTURE_REGISTRY.get(name, (None, f".{name}.parquet"))
                out_path = output_dir / f"{prefix}{suffix}"
                pq.write_table(merged, str(out_path), compression="zstd")

        return output_dir

    def _common_structures(self) -> list[str]:
        """Return structure names present in ALL datasets."""
        if not self.datasets:
            return []
        common = set(self.datasets[0].available_structures)
        for ds in self.datasets[1:]:
            common &= set(ds.available_structures)
        return sorted(common)

    def close(self):
        self._engine.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
