"""Physical merge of QPX datasets: streamed through the QPX writers, one view at a time, with explicit collision rules."""

from __future__ import annotations

import datetime
import hashlib
import logging
from pathlib import Path
from typing import TYPE_CHECKING

import pyarrow as pa
import pyarrow.parquet as pq

from qpx._version import __version__
from qpx.core.data.schema import SOURCE_DATASET_COLUMN
from qpx.core.sql import validate_identifier

if TYPE_CHECKING:
    from qpx.collection import DatasetCollection

logger = logging.getLogger(__name__)

ON_COLLISION = ("error", "namespace", "dedupe-identical")
SOURCE_COLUMN = SOURCE_DATASET_COLUMN
_BATCH_ROWS = 100_000
_METADATA_VIEWS = ("dataset", "ontology", "provenance", "pepmap")
_RUN_VIEWS = {
    "feature": "run_file_name",
    "psm": "run_file_name",
    "mz": "run_file_name",
    "run": "run_file_name",
    "pg": "unnest(grouped_runs)",
    "site": "unnest(grouped_runs)",
    "proteoform": "unnest(grouped_runs)",
}
# Namespacing: every id and every name that links rows of one source to each other.
_NAMESPACED_IDS = ("feature_id", "psm_id", "pg_id", "site_id", "proteoform_id", "best_psm_id")
_NAMESPACED_NAMES = ("run_file_name", "id_run_file_name", "quantification_unit_id", "run_accession", "sample_accession", "id")
_ONTOLOGY_FIRST = ("ontology_name", "ontology_accession", "ontology_source", "ontology_version", "description")


def _source_labels(datasets) -> list[str]:
    """One distinct label per dataset: its project accession, else its directory name."""
    labels: list[str] = []
    for index, ds in enumerate(datasets):
        accession = None
        if ds.dataset_meta is not None:
            accession = ds.dataset_meta.to_df()["project_accession"].iloc[0]
        name = str(ds.path).rstrip("/").rsplit("/", 1)[-1]
        candidates = [accession, name, f"{name}#{index}"]
        labels.append(next(c for c in candidates if c and c != "unknown" and c not in labels))
    return labels


def _namespace_constant(label: str) -> int:
    """A signed 64-bit constant per source: XOR keeps ids unique within it and links intact."""
    return int.from_bytes(hashlib.blake2b(label.encode(), digest_size=8).digest(), "big", signed=True)


def _literal(text: str) -> str:
    """A SQL string literal."""
    return "'" + text.replace("'", "''") + "'"


def _quoted(name: str) -> str:
    """A quoted SQL identifier."""
    return validate_identifier(name)


class _Merge:
    """One merge of a collection into *output_dir*."""

    def __init__(self, collection: "DatasetCollection", output_dir: Path, prefix: str, on_collision: str):
        self._collection = collection
        self._conn = collection._engine.connection
        self._datasets = collection.datasets
        self._output_dir = output_dir
        self._prefix = prefix
        self._on_collision = on_collision
        self._labels = _source_labels(self._datasets)

    # --- Selection ---

    def _sources(self, name: str) -> list[int]:
        """Indexes of the datasets that have view *name*."""
        return [i for i, ds in enumerate(self._datasets) if name in ds.available_structures]

    def _columns(self, name: str, index: int) -> list[str]:
        """The columns of a source's view."""
        return [row[0] for row in self._conn.execute(f'DESCRIBE "{name}_{index}"').fetchall()]

    def _writer_class(self, name: str):
        return self._datasets[0].writer_for(name)[0]

    def _suffix(self, name: str) -> str:
        return self._datasets[0].writer_for(name)[1]

    # --- Collisions ---

    def _key_sql(self, name: str) -> str:
        """The view's primary key columns, as SQL."""
        view = self._writer_class(name).view_schema()
        return ", ".join(_quoted(column) for column in view.primary_key)

    def _shared(self, values_sql: dict[int, str]) -> tuple[int, list[str]]:
        """The values seen in more than one source, as (count, examples)."""
        if len(values_sql) < 2:
            return 0, []
        union = " UNION ALL ".join(f"SELECT DISTINCT v, {index} AS s FROM ({sql}) t(v)" for index, sql in values_sql.items())
        sql = f"SELECT v FROM ({union}) GROUP BY v HAVING COUNT(*) > 1"
        rows = self._conn.execute(f"SELECT COUNT(*), list(CAST(v AS VARCHAR))[1:3] FROM ({sql})").fetchall()[0]
        return rows[0], rows[1] or []

    def _run_values(self, views: list[str], index: int) -> str | None:
        """SQL of the run names a source has in *views*, if any."""
        parts = [
            f'SELECT {expression} FROM "{name}_{index}"'
            for name, expression in _RUN_VIEWS.items()
            if name in views and index in self._sources(name)
        ]
        return " UNION ".join(parts) or None

    def collisions(self, views: list[str]) -> list[str]:
        """What the sources share: their runs, samples and the primary keys of their data views."""
        found = []
        runs = {i: sql for i in range(len(self._datasets)) if (sql := self._run_values(views, i))}
        count, examples = self._shared(runs)
        if count:
            found.append(f"{count} run(s), e.g. {examples}")
        for name in views:
            if name in _METADATA_VIEWS:
                continue
            key = self._key_sql(name)
            if "," in key:
                continue
            count, examples = self._shared({i: f'SELECT {key} FROM "{name}_{i}"' for i in self._sources(name)})
            if count:
                found.append(f"{count} {key.strip(chr(34))} value(s) in {name}, e.g. {examples}")
        return found

    # --- Row sources ---

    def _namespaced_select(self, columns: list[str], index: int) -> str:
        """SELECT list of a source with its ids and names rewritten into its namespace."""
        label = self._labels[index]
        constant = _namespace_constant(label)
        prefix = _literal(f"{label}:")
        parts = []
        for column in columns:
            quoted = _quoted(column)
            if column in _NAMESPACED_IDS:
                parts.append(f"xor({quoted}, {constant}::BIGINT) AS {quoted}")
            elif column in _NAMESPACED_NAMES:
                parts.append(f"{prefix} || {quoted} AS {quoted}")
            elif column == "grouped_runs":
                parts.append(f"list_transform({quoted}, lambda value: {prefix} || value) AS {quoted}")
            elif column == "samples":
                parts.append(
                    f"list_transform({quoted}, lambda value: struct_pack("
                    f"sample_accession := {prefix} || value.sample_accession, label := value.label, "
                    "biological_replicate := value.biological_replicate, "
                    f"technical_replicate := value.technical_replicate)) AS {quoted}"
                )
            else:
                parts.append(quoted)
        return ", ".join(parts)

    def _rows_sql(self, name: str, index: int) -> str:
        """SQL of a source's rows as they are written."""
        columns = self._columns(name, index)
        select = self._namespaced_select(columns, index) if self._on_collision == "namespace" else "*"
        return f'SELECT {select}, {_literal(self._labels[index])} AS {SOURCE_COLUMN} FROM "{name}_{index}"'

    def _deduplicated_rows(self, name: str, index: int, first: bool) -> str:
        """Rows of a source not already written; a key seen with different values raises."""
        key = self._key_sql(name)
        seen = f"merge_seen_{name}"
        source = f'"{name}_{index}"'
        if first:
            self._conn.execute(f"CREATE OR REPLACE TEMP TABLE {seen} AS SELECT {key}, hash(t) AS row_hash FROM {source} t")
            return self._rows_sql(name, index)
        different = self._conn.execute(
            f"SELECT COUNT(*) FROM (SELECT {key}, hash(t) AS row_hash FROM {source} t) a "
            f"JOIN {seen} b USING ({key}) WHERE a.row_hash <> b.row_hash"
        ).fetchall()[0][0]
        if different:
            raise ValueError(
                f"{different} {name} row(s) of {self._labels[index]} share a key with an earlier dataset but differ from it; "
                "they are not identical, so on_collision='dedupe-identical' cannot merge them. Use on_collision='namespace'."
            )
        rows = f"SELECT r.* FROM ({self._rows_sql(name, index)}) r ANTI JOIN {seen} USING ({key})"
        run = _RUN_VIEWS.get(name)
        if run:
            # A run shared without identical rows is two analyses of it: writing both double-counts it.
            shared = self._conn.execute(
                f"SELECT COUNT(*) FROM (SELECT {run} AS v FROM ({rows})) WHERE v IN (SELECT run_name FROM merge_runs_{name})"
            ).fetchall()[0][0]
            if shared:
                raise ValueError(
                    f"{shared} {name} row(s) or run(s) of {self._labels[index]} are in a run an earlier dataset also has, "
                    "and are not identical to its rows. Use on_collision='namespace'."
                )
        return rows

    def _record_seen(self, name: str, index: int) -> None:
        """Remember the keys and row hashes of a written source."""
        key = self._key_sql(name)
        self._conn.execute(
            f'INSERT INTO merge_seen_{name} SELECT a.* FROM (SELECT {key}, hash(t) AS row_hash FROM "{name}_{index}" t) a '
            f"ANTI JOIN merge_seen_{name} USING ({key})"
        )

    # --- Data views ---

    def _composite_groups(self, name: str, sources: list[int]) -> dict[tuple[str, ...] | None, list[int]]:
        """Sources by the identity composite their footers declare; only one composite fits one file's footer."""
        groups: dict[tuple[str, ...] | None, list[int]] = {}
        for index in sources:
            struct = getattr(self._datasets[index], name)
            composite = None
            first = str(struct.file_paths[0]) if struct.file_paths else ""
            if first.endswith(".parquet") and not first.startswith("s3://") and "*" not in first:
                metadata = pq.read_schema(first).metadata or {}
                if b"identity_composite" in metadata:
                    composite = tuple(metadata[b"identity_composite"].decode().split(","))
            groups.setdefault(composite, []).append(index)
        return groups

    def _writer(self, name: str, path: Path, composite, extra: list[str]):
        """The view's writer with the merge's columns, keeping every provided id."""
        return self._writer_class(name)(path, identity_composite=composite, override_provided_ids=False, extra_columns=extra)

    def _write_batches(self, writer, sql: str) -> None:
        """Stream the rows of *sql* into *writer*."""
        reader = self._conn.sql(sql).to_arrow_reader(_BATCH_ROWS)
        for batch in reader:
            if batch.num_rows:
                writer.write_table(writer.align_table_to_schema(pa.Table.from_batches([batch])))

    def _write_data_view(self, name: str) -> None:
        """Stream a view's rows from every source, one file per identity composite."""
        sources = self._sources(name)
        extra = [SOURCE_COLUMN]
        if name == "sample":
            known = set(self._writer_class(name).view_schema().get_arrow_schema().names)
            extra = sorted({c for i in sources for c in self._columns(name, i)} - known) + extra
        dedupe = self._on_collision == "dedupe-identical"
        run = _RUN_VIEWS.get(name)
        if dedupe and run:
            self._conn.execute(f"CREATE OR REPLACE TEMP TABLE merge_runs_{name} (run_name VARCHAR)")
        groups = self._composite_groups(name, sources)
        suffix = self._suffix(name)
        first_source = True
        for number, (composite, members) in enumerate(groups.items()):
            file_prefix = self._prefix if len(groups) == 1 else f"{self._prefix}-{number + 1}"
            with self._writer(name, self._output_dir / f"{file_prefix}{suffix}", composite, extra) as writer:
                for index in members:
                    sql = self._deduplicated_rows(name, index, first_source) if dedupe else self._rows_sql(name, index)
                    self._write_batches(writer, sql)
                    if dedupe:
                        if not first_source:
                            self._record_seen(name, index)
                        if run:
                            self._conn.execute(f'INSERT INTO merge_runs_{name} SELECT DISTINCT {run} FROM "{name}_{index}"')
                    first_source = False

    # --- Metadata views ---

    def _write_records(self, name: str, table: pa.Table) -> None:
        """Write a metadata view's table."""
        writer_class = self._writer_class(name)
        with writer_class(self._output_dir / f"{self._prefix}{self._suffix(name)}") as writer:
            writer.write_table(writer.align_table_to_schema(table))

    def _write_dataset_view(self) -> None:
        """One row: the view describes the merged dataset, and each source stays in provenance and ``source_dataset``."""
        rows = [self._conn.sql(f'SELECT * FROM "dataset_{i}"').to_arrow_table().to_pylist()[0] for i in self._sources("dataset")]
        accessions = list(dict.fromkeys(row["project_accession"] for row in rows))
        names = list(dict.fromkeys(row.get("software_name") for row in rows if row.get("software_name")))
        record = {
            "project_accession": ";".join(accessions),
            "project_title": rows[0].get("project_title") if len(accessions) == 1 else f"Merge of {', '.join(accessions)}",
            "software_name": ";".join(names) or None,
            "creation_date": datetime.date.today().isoformat(),
        }
        self._write_records("dataset", pa.Table.from_pylist([record]))

    def _write_ontology_view(self) -> None:
        """One row per field and view; the tools and columns that report it are joined."""
        sources = self._sources("ontology")
        union = " UNION ALL BY NAME ".join(f'SELECT * FROM "ontology_{i}"' for i in sources)
        firsts = ", ".join(f"first({column}) AS {column}" for column in _ONTOLOGY_FIRST)
        conflicts = self._conn.execute(
            f'SELECT field_name, "view" FROM ({union}) GROUP BY field_name, "view" HAVING COUNT(DISTINCT ontology_accession) > 1'
        ).fetchall()
        if conflicts:
            logger.warning("%d ontology field(s) have different terms across the datasets; keeping the first: %s",
                           len(conflicts), conflicts[:5])  # fmt: skip
        table = self._conn.sql(
            f'SELECT field_name, "view", {firsts}, '
            "string_agg(DISTINCT source_column_name, ';') AS source_column_name, "
            "string_agg(DISTINCT source_tool, ';') AS source_tool "
            f'FROM ({union}) GROUP BY field_name, "view" ORDER BY "view", field_name'
        ).to_arrow_table()
        self._write_records("ontology", table)

    def _write_provenance_view(self, views: list[str]) -> None:
        """Each source's steps in order, numbered through, then the merge itself."""
        rows: list[dict] = []
        for index in self._sources("provenance"):
            steps = self._conn.sql(f'SELECT * FROM "provenance_{index}" ORDER BY step_order').to_arrow_table().to_pylist()
            rows.extend({**step, SOURCE_COLUMN: self._labels[index]} for step in steps)
        rows.append(
            {
                "step_category": "format_conversion",
                "step_name": "merge_datasets",
                "tool_name": "qpx",
                "tool_version": __version__,
                "parameters": [{"key": "on_collision", "value": self._on_collision}],
                "output_views": views,
                SOURCE_COLUMN: None,
            }
        )
        for number, row in enumerate(rows, start=1):
            row["step_order"] = number
        with self._writer(
            "provenance", self._output_dir / f"{self._prefix}{self._suffix('provenance')}", None, [SOURCE_COLUMN]
        ) as writer:
            schema = writer.arrow_schema
            writer.write_table(pa.Table.from_pylist([{f: row.get(f) for f in schema.names} for row in rows], schema=schema))

    def _write_pepmap_view(self) -> None:
        """One row per peptidoform, mapping to the proteins of every source."""
        union = " UNION ALL BY NAME ".join(f'SELECT * FROM "pepmap_{i}"' for i in self._sources("pepmap"))
        table = self._conn.sql(
            "SELECT first(sequence) AS sequence, peptidoform, list_distinct(flatten(list(pg_accessions))) AS pg_accessions, "
            "list_distinct(flatten(list(sequence_edits))) AS sequence_edits "
            f"FROM ({union}) GROUP BY peptidoform"
        )
        reader = table.to_arrow_reader(_BATCH_ROWS)
        with self._writer_class("pepmap")(self._output_dir / f"{self._prefix}{self._suffix('pepmap')}") as writer:
            for batch in reader:
                self._write_pepmap_batch(writer, batch)

    @staticmethod
    def _write_pepmap_batch(writer, batch: pa.RecordBatch) -> None:
        accessions = batch.column("pg_accessions").to_pylist()
        unique = pa.array([len({p["accession"] for p in row or []}) == 1 for row in accessions], type=pa.bool_())
        table = pa.Table.from_batches([batch]).append_column("is_unique", unique)
        writer.write_table(writer.align_table_to_schema(table))

    # --- Entry point ---

    def run(self, views: list[str]) -> None:
        """Check for collisions, then write every view."""
        if self._on_collision == "error":
            collisions = self.collisions(views)
            if collisions:
                raise ValueError(
                    "The datasets overlap: "
                    + "; ".join(collisions)
                    + ". Merging them would count those rows twice. Pass on_collision='namespace' to keep both, each under "
                    "its source's name, or on_collision='dedupe-identical' if they are copies of the same data."
                )
        self._output_dir.mkdir(parents=True, exist_ok=True)
        for name in views:
            if name == "dataset":
                self._write_dataset_view()
            elif name == "ontology":
                self._write_ontology_view()
            elif name == "provenance":
                self._write_provenance_view(views)
            elif name == "pepmap":
                self._write_pepmap_view()
            else:
                self._write_data_view(name)


def merge_views(collection: "DatasetCollection", requested: list[str] | None) -> list[str]:
    """The views to merge: warns about those that are left out, raises when none remain."""
    datasets = collection.datasets
    available = [set(ds.available_structures) for ds in datasets]
    if requested is None:
        views = sorted(set.intersection(*available)) if available else []
        dropped = sorted(set().union(*available) - set(views)) if available else []
        if dropped:
            logger.warning("Not merging %s: only views present in every dataset are merged by default", dropped)
    else:
        known = datasets[0].structure_names() if datasets else []
        unknown = sorted(set(requested) - set(known))
        if unknown:
            raise ValueError(f"Unknown structure(s) {unknown}. Valid: {sorted(known)}")
        for name in requested:
            missing = [str(ds.path) for ds, names in zip(datasets, available) if name not in names]
            if missing:
                logger.warning("%s is missing from %s", name, missing)
        views = [name for name in requested if any(name in names for names in available)]
    if not views:
        raise ValueError(f"Nothing to merge: the datasets have no {'common ' if requested is None else 'requested '}views")
    return views


def merge_datasets(
    collection: "DatasetCollection",
    output_dir: str | Path,
    structures: list[str] | None,
    prefix: str,
    on_collision: str,
) -> Path:
    """Merge the datasets of *collection* into *output_dir*; see :meth:`DatasetCollection.merge`."""
    if on_collision not in ON_COLLISION:
        raise ValueError(f"on_collision must be one of {ON_COLLISION}, got {on_collision!r}")
    views = merge_views(collection, structures)
    output_dir = Path(output_dir)
    _Merge(collection, output_dir, prefix, on_collision).run(views)
    return output_dir
