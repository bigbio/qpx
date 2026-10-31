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

        # Merge into a new dataset (refuses overlapping datasets; see ``merge``)
        coll.merge("merged_output/", structures=["feature"])
    """

    def __init__(
        self,
        datasets: list[Dataset],
        duckdb_memory: str | None = None,
        duckdb_threads: int | None = None,
        duckdb_temp_dir: str | Path | None = None,
    ):
        self.datasets = datasets
        self._engine = DuckDBEngine(
            memory_limit=duckdb_memory,
            threads=duckdb_threads,
            temp_directory=duckdb_temp_dir,
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
        return QueryResult(self._engine.connection, query)

    @property
    def structure_names(self) -> dict[int, list[str]]:
        """Return available structures per dataset index."""
        return {i: ds.available_structures for i, ds in enumerate(self.datasets)}

    def merge(
        self,
        output_dir: str | Path,
        structures: list[str] | None = None,
        prefix: str = "merged",
        on_collision: str = "error",
    ) -> Path:
        """Physically merge datasets into a new directory.

        Rows are streamed through the QPX writers, so the output has QPX footers
        and row groups, whatever the memory limit of the collection allows. Data
        views get a ``source_dataset`` column, the project accession or directory
        name of the row's dataset. ``dataset`` becomes one row, ``provenance`` the
        steps of every source followed by the merge, ``ontology`` and ``pepmap``
        one row per key. Without ``structures``, only views present in every
        dataset are merged, and the others are named in a warning.

        Datasets with different identity composites give one file per composite,
        because a footer declares one.

        ``on_collision`` decides what happens when datasets share a run, a sample
        or a primary key, which would count those rows twice:

        - ``"error"`` (default): refuse, naming what they share, before writing.
        - ``"namespace"``: keep both. Run and sample names and other names that
          link rows become ``<source>:<name>``, and ids are rewritten per source.
        - ``"dedupe-identical"``: for copies of the same data. Rows with an
          already written key are dropped if identical; a difference raises.
        """
        from qpx.core.dataset_merge import merge_datasets

        return merge_datasets(self, output_dir, structures, prefix, on_collision)

    def close(self):
        self._engine.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
