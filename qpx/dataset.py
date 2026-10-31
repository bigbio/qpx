"""Dataset class — the central entry point for opening QPX datasets."""

from __future__ import annotations

import logging
from pathlib import Path

import duckdb

from qpx.core.data import (
    PG,
    PSM,
    BaseStructure,
    DatasetMeta,
    Feature,
    MzSpectra,
    Ontology,
    PepMap,
    Protein,
    Proteoform,
    Provenance,
    Run,
    Sample,
    Site,
)
from qpx.core.dataset_queries import _DatasetQueryMixin
from qpx.core.dataset_validation import _DatasetValidationMixin
from qpx.core.dataset_writes import _DatasetWriteMixin
from qpx.core.engine import DuckDBEngine
from qpx.core.sql import sql_build, validate_table
from qpx.version import (
    QpxVersionError,
    check_feature_file_compatible,
    check_feature_schema_compatible,
    check_pg_columns_compatible,
    check_pg_file_compatible,
)

_log = logging.getLogger(__name__)

# Per-file layout guards run on every local file of a structure before it is registered.
_FILE_LAYOUT_GUARDS = {"pg": check_pg_file_compatible, "feature": check_feature_file_compatible}


def _require_local_directory(path: str | Path) -> None:
    """Raise for a location a local dataset cannot be read from; discovery would find nothing in it and say so nowhere."""
    if "://" in str(path):
        raise ValueError(f"Cannot open {path}: only local directories and s3:// locations are supported")
    if not Path(path).is_dir():
        raise FileNotFoundError(f"QPX dataset directory not found: {path}")


class Dataset(_DatasetValidationMixin, _DatasetQueryMixin, _DatasetWriteMixin):
    """
    A QPX dataset — a directory of related Parquet/H5AD files.

    Usage:
        ds = qpx.Dataset("PXD014414/")
        ds.feature.filter("charge > 2").to_df()
        ds.pg.join(ds.run).to_df()

    Partial access:
        ds = qpx.Dataset("path/", structures=["feature", "sample", "run"])
    """

    # Data structure registry: name → (Class, file_suffix)
    _STRUCTURE_REGISTRY = {
        "psm": (PSM, ".psm.parquet"),
        "feature": (Feature, ".feature.parquet"),
        "pg": (PG, ".pg.parquet"),
        "mz": (MzSpectra, ".mz.parquet"),
        "sample": (Sample, ".sample.parquet"),
        "run": (Run, ".run.parquet"),
        "dataset": (DatasetMeta, ".dataset.parquet"),
        "ontology": (Ontology, ".ontology.parquet"),
        "provenance": (Provenance, ".provenance.parquet"),
        "pepmap": (PepMap, ".pepmap.parquet"),
        "site": (Site, ".site.parquet"),
        "proteoform": (Proteoform, ".proteoform.parquet"),
        "protein": (Protein, ".protein.parquet"),
    }

    @classmethod
    def structure_names(cls) -> list[str]:
        """Every structure name QPX knows about, in registry order.

        The single source of truth for anything that needs to enumerate
        structures — CLI choices in particular. Hand-maintained copies drifted:
        ``pepmap`` was registered here and exposed on the API while both CLIs
        rejected it (bigbio/qpx#289).
        """
        return list(cls._STRUCTURE_REGISTRY)

    def __init__(
        self,
        path: str | Path,
        structures: list[str] | None = None,
        duckdb_memory: str | None = None,
        duckdb_threads: int | None = None,
        s3_config: dict | None = None,
        file_prefix: str | None = None,
        duckdb_temp_dir: str | Path | None = None,
    ):
        self._is_s3 = isinstance(path, str) and path.startswith("s3://")
        if not self._is_s3:
            _require_local_directory(path)
        self.path = path if self._is_s3 else Path(path)
        self._file_prefix = file_prefix
        self._engine = DuckDBEngine(
            memory_limit=duckdb_memory,
            threads=duckdb_threads,
            temp_directory=duckdb_temp_dir,
            s3_config=s3_config if self._is_s3 else None,
        )

        self._structures: dict[str, BaseStructure] = {}
        # Structures whose files exist but could not be registered, with the reason.
        self._skipped: dict[str, str] = {}
        self._discover_and_register(structures)

    def _discover_and_register(self, requested: list[str] | None):
        """Scan directory, find QPX files, register as DuckDB tables.

        Checks for single Parquet files first, then falls back to
        Hive-partitioned directories (e.g., feature/ with run_file_name= subdirs).
        For S3 paths, attempts to register each structure via S3 glob.
        """
        if self._is_s3:
            self._discover_s3(requested)
        else:
            self._discover_local(requested)

    def _discover_s3(self, requested: list[str] | None):
        """Register structures from S3 path."""
        for name, (cls, suffix) in self._STRUCTURE_REGISTRY.items():
            if requested and name not in requested:
                continue
            file_name = f"{self._file_prefix}{suffix}" if self._file_prefix else f"*{suffix}"
            s3_glob = f"{self.path.rstrip('/')}/{file_name}"
            try:
                # PyArrow cannot open an S3 glob to read the footer, so register
                # first and run the PG guard against the DuckDB view's columns.
                self._engine.register_s3_parquet(name, s3_glob)
                if name == "pg":
                    cols = [row[0] for row in self._engine.execute('DESCRIBE "pg"').fetchall()]
                    check_pg_columns_compatible(cols, source=s3_glob)
                elif name == "feature":
                    schema = self._engine.connection.sql('SELECT * FROM "feature" LIMIT 0').to_arrow_table().schema
                    check_feature_schema_compatible(schema, s3_glob)
                self._structures[name] = cls(
                    engine=self._engine,
                    table_name=name,
                    file_path=f"{self.path}/{name}",
                    # The glob is what was actually registered; a consumer
                    # re-registering from the display path would read nothing.
                    file_paths=[s3_glob],
                )
            except QpxVersionError as exc:
                # One incompatible/old structure file must not abort the whole
                # dataset — skip it with a warning so the other structures stay
                # usable. (A direct read_pg()/PG.from_file() still raises.)
                self._engine.execute(
                    sql_build(
                        "DROP VIEW IF EXISTS $view",
                        view=validate_table(name),
                    )
                )
                _log.warning("Skipping incompatible structure '%s': %s", name, exc)
                self._skipped[name] = str(exc)
            except (FileNotFoundError, duckdb.IOException):
                pass  # Structure not present in S3
            except (OSError, ValueError, duckdb.Error) as exc:
                _log.warning("Failed to register '%s' from S3: %s", name, exc)
                self._skipped[name] = str(exc)

    def _discover_local(self, requested: list[str] | None):
        """Register structures from local filesystem.

        Each structure is registered independently: a bad/old/garbage file for
        one structure is skipped with a warning so the rest of the dataset stays
        usable, instead of aborting the whole Dataset construction.
        """
        for name, (cls, suffix) in self._STRUCTURE_REGISTRY.items():
            if requested and name not in requested:
                continue
            try:
                self._register_local_structure(name, cls, suffix)
            except QpxVersionError as exc:
                # One incompatible/old structure file must not abort the whole
                # dataset — skip it with a warning so the other structures stay
                # usable. (A direct read_pg()/PG.from_file() still raises.)
                _log.warning("Skipping incompatible structure '%s': %s", name, exc)
                self._skipped[name] = str(exc)
            except (OSError, ValueError, duckdb.Error) as exc:
                _log.warning("Failed to register structure '%s': %s", name, exc)
                self._skipped[name] = str(exc)

    def _register_local_structure(self, name: str, cls, suffix: str):
        """Register a single local structure (single file or Hive-partitioned dir)."""
        # Check for single file first
        pattern = f"{self._file_prefix}{suffix}" if self._file_prefix else f"*{suffix}"
        matches = sorted(self.path.glob(pattern))
        if matches:
            # Read ALL matching shards unioned, matching the S3 path (which globs
            # and reads every match). Previously only the first was taken, so a
            # sharded structure returned different rows locally vs over S3
            # (bigbio/qpx#252). The registered structure's file_path points at the
            # first shard for provenance; the view unions them all.
            guard = _FILE_LAYOUT_GUARDS.get(name)
            if guard is not None:
                for match in matches:
                    guard(match)
            if len(matches) > 1:
                _log.warning(
                    "%s is read from %d files as one view (%s, ...). If they are not shards of one dataset, "
                    "their rows are counted together: pass file_prefix to read one.",
                    name,
                    len(matches),
                    matches[0].name,
                )
            self._engine.register_parquet_files(name, matches)
            self._structures[name] = cls(
                engine=self._engine,
                table_name=name,
                file_path=matches[0],
                file_paths=matches,
            )
        elif self._file_prefix is None:
            # Check for Hive-partitioned directory
            part_dir = self.path / name
            part_files = list(part_dir.glob("**/*.parquet")) if part_dir.is_dir() else []
            if part_files:
                guard = _FILE_LAYOUT_GUARDS.get(name)
                if guard is not None:
                    for part_file in sorted(part_files):
                        guard(part_file)
                self._engine.register_partitioned_parquet(name, part_dir)
                self._structures[name] = cls(
                    engine=self._engine,
                    table_name=name,
                    file_path=part_dir,
                )

    # --- Data structure accessors (lazy, return None if not present) ---
    @property
    def psm(self) -> PSM | None:
        return self._structures.get("psm")

    @property
    def feature(self) -> Feature | None:
        return self._structures.get("feature")

    @property
    def pg(self) -> PG | None:
        return self._structures.get("pg")

    @property
    def mz(self) -> MzSpectra | None:
        return self._structures.get("mz")

    @property
    def sample(self) -> Sample | None:
        return self._structures.get("sample")

    @property
    def run(self) -> Run | None:
        return self._structures.get("run")

    @property
    def dataset_meta(self) -> DatasetMeta | None:
        return self._structures.get("dataset")

    @property
    def ontology(self) -> Ontology | None:
        return self._structures.get("ontology")

    @property
    def provenance(self) -> Provenance | None:
        return self._structures.get("provenance")

    @property
    def pepmap(self) -> PepMap | None:
        return self._structures.get("pepmap")

    @property
    def site(self) -> Site | None:
        """The site view, or None when the dataset has no site.parquet."""
        return self._structures.get("site")

    @property
    def protein(self) -> Protein | None:
        """The protein view, or None when the dataset has no protein.parquet."""
        return self._structures.get("protein")

    @property
    def proteoform(self) -> Proteoform | None:
        """The proteoform view, or None when the dataset has no proteoform.parquet."""
        return self._structures.get("proteoform")

    @property
    def available_structures(self) -> list[str]:
        return list(self._structures.keys())

    @property
    def skipped_structures(self) -> dict[str, str]:
        """Structures whose files are present but could not be read (an older layout, a corrupt file), with the reason."""
        return dict(self._skipped)

    @property
    def file_prefix(self) -> str | None:
        """Prefix of the QPX file names this dataset was opened with, if any."""
        return self._file_prefix

    # --- External registration (for downstream tool integration) ---
    def register_external(self, name: str, file_path: str | Path) -> None:
        """Register an external Parquet file in the DuckDB engine."""
        self._engine.register_parquet(name, file_path)

    # --- Lifecycle ---
    def refresh(self) -> None:
        """Re-scan the dataset directory for new or updated files.

        Clears all cached view instances so next access gets fresh data.
        Call this after writing new structures or AnnData files.
        """
        for attr in (
            "protein_view",
            "peptide_view",
            "identification_summary",
            "run_summary",
            "modification_view",
            "qc_view",
            "sample_summary",
            "ae_view",
        ):
            self.__dict__.pop(attr, None)

        self._structures.clear()
        self._skipped.clear()
        self._discover_and_register(None)

    def close(self):
        self._engine.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __repr__(self):
        structs = ", ".join(self.available_structures)
        return f"Dataset('{self.path}', structures=[{structs}])"
