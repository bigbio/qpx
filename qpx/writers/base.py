"""BaseWriter — batched Parquet writer with schema validation."""

from __future__ import annotations

import base64
import datetime
import logging
import os
import uuid
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Sequence

import pyarrow as pa
import pyarrow.parquet as pq

from qpx._version import __version__
from qpx.core.data.identity import IDENTITY_HASH_VERSION, derive_composite_id
from qpx.core.data.schema import _scan_format_issues, _scan_format_violations, pg_referential_issues_from_parquet
from qpx.core.engine import create_connection
from qpx.core.sql import sql_build, validate_identifier
from qpx.version import QPX_SPEC_VERSION
from qpx.writers.encoding import parquet_write_options, row_group_rows

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)


def _cast_column(column: pa.ChunkedArray, target: pa.DataType) -> pa.ChunkedArray | pa.Array:
    """Cast *column*; a null struct read back from DuckDB has null children, which a required child field refuses."""
    try:
        return column.cast(target)
    except pa.ArrowInvalid:
        return pa.array(column.to_pylist(), type=target)


@dataclass
class _IdentityState:
    """A writer's identity composite override and where the ids written so far came from."""

    composite_override: tuple[str, ...] | None
    sources: set[str] = dataclass_field(default_factory=set)


class BaseWriter:
    """
    Batched Parquet writer with schema validation.

    Usage:
        writer = FeatureWriter("output.feature.parquet", creator="quantms")
        writer.write_batch(records)   # list[dict] or pd.DataFrame
        writer.write_table(table)     # pa.Table
        writer.close()
    """

    _schema_class = None  # Override in subclasses
    # Keys a record may carry that the writer turns into schema columns.
    _record_only_keys: tuple[str, ...] = ()

    def __init__(
        self,
        path: str | Path,
        creator: str = "qpx",
        software_provider: str = "qpx",
        compression: str = "zstd",
        batch_size: int = 100_000,
        scan_format: str | None = None,
        extra_columns: list[str] | None = None,
        override_provided_ids: bool = False,
        identity_composite: Optional[Sequence[str]] = None,
    ):
        self._path = Path(path)
        # A per-instance nonce, not just the pid: two writers in the SAME process
        # targeting one destination would otherwise share a staging file and
        # corrupt each other.
        # Bytes go to a sibling temp file and are moved onto _path only after a
        # clean close. Opening the destination directly meant a converter that
        # failed mid-run replaced an existing valid file with a truncated one and
        # there was no way back (bigbio/qpx#288). A sibling keeps the rename on
        # the same filesystem, so os.replace is atomic.
        self._staging_path = self._path.with_name(f".{self._path.name}.{os.getpid()}.{uuid.uuid4().hex[:8]}.tmp")
        self._compression = compression
        self._batch_size = batch_size
        self._buffer: list[dict] = []
        self._writer: pq.ParquetWriter | None = None
        # Set when a write raises. __exit__ sees the exception and discards, but a
        # converter that catches its own write error and then calls close()
        # normally would otherwise promote the partial file — publishing exactly
        # the truncated output the staging was added to prevent (bigbio/qpx#288).
        self._write_failed = False
        self._warned_checks: set[str] = set()
        self._extra_columns = extra_columns or []
        # Override mode (e.g. consensusXML): re-derive the id even when the
        # producer supplied one, stashing the original as a cv_param so it is
        # preserved/traceable.
        self._override_provided_ids = override_provided_ids
        # Producer-specific identity composite override. When set, the writer
        # hashes THIS composite into the derived id (feature_id / psm_id / pg_id)
        # instead of the schema default, letting each converter declare the
        # measured keys its producer emits (bigbio/qpx#229). Stored as a tuple so
        # the footer stamps the effective composite via the ``_identity_composite``
        # property below.
        self._identity = _IdentityState(tuple(identity_composite) if identity_composite else None)
        self._validate_identity_composite()

        # Build Parquet footer metadata. Two distinct versions are stamped:
        #   qpx_version    — the on-disk *specification* version (QPX_SPEC_VERSION),
        #                    what readers check for format compatibility.
        #   writer_version — the *package* build (__version__, hatch-vcs), for
        #                    provenance/debugging only.
        self._file_metadata = {
            b"qpx_version": QPX_SPEC_VERSION.encode(),
            b"writer_version": __version__.encode(),
            b"creator": creator.encode(),
            b"file_type": self._schema_class._file_type.encode(),
            b"creation_date": datetime.datetime.now().isoformat().encode(),
            b"uuid": str(uuid.uuid4()).encode(),
            b"software_provider": software_provider.encode(),
            b"compression_format": compression.encode(),
        }
        if scan_format:
            self._file_metadata[b"scan_format"] = scan_format.encode()

        # Stamp the effective composite in the footer. A single-column derived
        # identity also records its id field as the primary key.
        identity_composite = self._identity_composite
        if identity_composite:
            id_field = self._id_field
            if id_field is not None:
                self._file_metadata[b"primary_key"] = id_field.encode()
            self._file_metadata[b"identity_composite"] = ",".join(identity_composite).encode()
            self._file_metadata[b"identity_hash_version"] = IDENTITY_HASH_VERSION.encode()

    def _validate_identity_composite(self) -> None:
        """Reject composite fields absent from the effective writer schema."""
        effective_composite = self._identity_composite
        if effective_composite:
            identity_schema = (
                self._schema_class.get_extended_arrow_schema(self._extra_columns)
                if self._extra_columns
                else self._schema_class.get_arrow_schema()
            )
            available_fields = set(identity_schema.names)
            unknown_fields = [field for field in effective_composite if field not in available_fields]
            if unknown_fields:
                raise ValueError(
                    f"identity_composite contains fields outside the {self._schema_class.view_name} schema: {unknown_fields}"
                )

    # --- Derived identity --------------------------------------------------

    @property
    def _identity_composite(self) -> tuple[str, ...] | None:
        """The effective identity_composite the id is derived from.

        A converter-supplied override (producer-specific composite) wins over the
        schema class attribute; ``None`` for views without any derived identity.
        """
        if self._identity.composite_override is not None:
            return self._identity.composite_override
        return getattr(self._schema_class, "identity_composite", None)

    @property
    def _id_field(self) -> str | None:
        """The single-column identity primary key (the id field), else None."""
        if not self._identity_composite:
            return None
        pk = self._schema_class.primary_key
        return pk[0] if len(pk) == 1 else None

    def _relaxed_id_schema(self) -> pa.Schema:
        """``arrow_schema`` with the derived-id field made nullable, so a batch or
        table can be materialized before :meth:`_fill_identity_table` fills the id
        (casting a null id to the required non-nullable field would otherwise fail).
        """
        id_field = self._id_field
        if id_field is None:
            return self.arrow_schema
        return pa.schema([f.with_nullable(True) if f.name == id_field else f for f in self.arrow_schema])

    def _check_required_non_null(self, obj) -> None:
        """Raise if any non-nullable schema column has nulls (Table or RecordBatch)."""
        for i, field in enumerate(self.arrow_schema):
            if not field.nullable and obj.column(i).null_count > 0:
                raise ValueError(
                    f"Column '{field.name}' has {obj.column(i).null_count} null(s) but is marked as required in the schema"
                )

    def _persisted_composite_values(self, table: pa.Table, composite: Sequence[str]) -> dict[str, list]:
        """Return composite columns cast to their persisted schema types."""
        names = set(table.schema.names)
        schema = self.arrow_schema
        values: dict[str, list] = {}
        for field_name in composite:
            if field_name not in names:
                values[field_name] = [None] * table.num_rows
                continue
            column = table.column(field_name)
            field_type = schema.field(field_name).type
            if column.type != field_type:
                try:
                    column = column.cast(field_type)
                except (pa.ArrowInvalid, pa.ArrowTypeError, pa.ArrowNotImplementedError):
                    # Schema validation below reports incompatible source types.
                    pass
            values[field_name] = column.to_pylist()
        return values

    def _identity_values(
        self,
        existing: list,
        composite_values: dict[str, list],
        composite: Sequence[str],
        cv_lists: list | None,
        id_field: str,
    ) -> list[int]:
        """Preserve or derive row IDs, recording overridden producer IDs."""
        ids: list[int] = []
        for index, provided in enumerate(existing):
            if provided is not None and not self._should_override_provided_id(index, composite_values):
                self._identity.sources.add("provided")
                ids.append(provided)
                continue
            self._identity.sources.add("derived")
            if provided is not None:
                if cv_lists is not None:
                    params = list(cv_lists[index] or [])
                    provided_param = {"cv_name": f"provided_{id_field}", "cv_value": str(provided)}
                    if provided_param not in params:
                        params.append(provided_param)
                    cv_lists[index] = params
            ids.append(derive_composite_id({field: composite_values[field][index] for field in composite}, composite))
        return ids

    def _should_override_provided_id(self, _index: int, _composite_values: dict[str, list]) -> bool:
        """Return whether a producer id should be replaced for one row."""
        return self._override_provided_ids

    def _keeps_every_id(self, table: pa.Table, id_field: str) -> bool:
        """Whether every row keeps its provided id, so the composite need not be read."""
        return not self._override_provided_ids and id_field in table.column_names and table.column(id_field).null_count == 0

    def _derive_columns(self, table: pa.Table) -> pa.Table:
        """Fill the columns this view derives from the others; none by default."""
        return table

    def _fill_identity_table(self, table: pa.Table) -> pa.Table:
        """Return *table* with its derived-id column present and non-null.

        The table is assumed to already match ``arrow_schema`` in every other
        column (callers build via ``arrow_schema`` or align first). Rows with a
        non-null id keep it; the rest are derived from the composite columns.
        In ``override_provided_ids`` mode the id is re-derived even when present,
        stashing the original as a ``provided_<id>`` cv_param. This is the single
        derivation point for every write path (record, table, DataFrame).
        """
        table = self._derive_columns(table)
        composite = self._identity_composite
        id_field = self._id_field
        if not composite or id_field is None:
            return table
        names = set(table.schema.names)
        id_arrow_field = pa.field(id_field, pa.int64(), nullable=False)
        if self._keeps_every_id(table, id_field):
            self._identity.sources.add("provided")
            return table.set_column(
                table.schema.get_field_index(id_field), id_arrow_field, table.column(id_field).cast(pa.int64())
            )
        n = table.num_rows
        # Hash the PERSISTED (schema-cast) value of each composite field so the id
        # is identical regardless of which write path built the table — e.g. a
        # float32 ``rt`` hashes to the same value whether it arrived as a raw
        # Python float (record path) or an already-cast float32 column (table path).
        composite_values = self._persisted_composite_values(table, composite)
        existing = table.column(id_field).to_pylist() if id_field in names else [None] * n
        override = self._override_provided_ids
        any_provided = any(provided is not None for provided in existing)
        if override and "cv_params" not in names and any_provided:
            raise ValueError(f"override_provided_ids requires a 'cv_params' column to preserve provided {id_field} values")
        # Only round-trip cv_params when there is a provided id to stash.
        cv_lists = table.column("cv_params").to_pylist() if (override and any_provided and "cv_params" in names) else None
        ids = self._identity_values(existing, composite_values, composite, cv_lists, id_field)
        if cv_lists is not None:
            cv_type = table.schema.field("cv_params").type
            table = table.set_column(
                table.schema.get_field_index("cv_params"),
                table.schema.field("cv_params"),
                pa.array(cv_lists, type=cv_type),
            )
        id_array = pa.array(ids, type=pa.int64())
        if id_field in names:
            return table.set_column(table.schema.get_field_index(id_field), id_arrow_field, id_array)
        return table.add_column(0, id_arrow_field, id_array)

    def align_table_to_schema(self, table: pa.Table) -> pa.Table:
        """Project *table* onto the schema: keep matching columns, fill absent ones
        (e.g. the derived id or the optional cross-ref columns) with nulls, then
        cast. The id is left nullable here so a not-yet-derived (null) id is valid;
        :meth:`write_table` is the single place that derives it, so the id is never
        derived twice (which in override mode would double-stash the provided id).
        """
        relaxed = self._relaxed_id_schema()
        present = set(table.schema.names)
        n = table.num_rows
        columns = [_cast_column(table.column(f.name), f.type) if f.name in present else pa.nulls(n, type=f.type) for f in relaxed]
        return pa.Table.from_arrays(columns, schema=relaxed)

    @property
    def arrow_schema(self) -> pa.Schema:
        if self._extra_columns:
            base = self._schema_class.get_extended_arrow_schema(self._extra_columns)
        else:
            base = self._schema_class.get_arrow_schema()
        return base.with_metadata(self._file_metadata)

    def _guard(self, fn, *args, **kwargs):
        """Run a write, latching failure so close() can never promote the result."""
        try:
            return fn(*args, **kwargs)
        except BaseException:
            self._write_failed = True
            raise

    def write_batch(self, records: list[dict]):
        """Accumulate records and flush when batch_size is reached."""
        self._guard(self._check_record_keys, records)
        return self._guard(self._write_batch, records)

    def set_scan_format(self, scan_format: str) -> None:
        """Declare a confirmed format before close, including after streamed batches."""
        self._guard(self._set_scan_format, scan_format)

    def identity_source_metadata(self) -> dict[bytes, bytes]:
        """The ``identity_source`` footer key for the ids written so far; empty before any."""
        sources = self._identity.sources
        if not sources:
            return {}
        source = next(iter(sources)) if len(sources) == 1 else "mixed"
        return {b"identity_source": source.encode()}

    def add_footer_metadata(self, entries: dict[bytes, bytes]) -> None:
        """Add footer keys before close, including after streamed batches."""
        self._guard(self._add_footer_metadata, entries)

    def _add_footer_metadata(self, entries: dict[bytes, bytes]) -> None:
        if self._writer is not None and not hasattr(self._writer, "add_key_value_metadata"):
            logger.warning("Cannot add %s after batches were written: PyArrow 17 or newer is required", sorted(entries))
            return
        self._file_metadata.update(entries)
        if self._writer is not None:
            self._writer.add_key_value_metadata(entries)

    def _set_scan_format(self, scan_format: str) -> None:
        """Keep footer keys and the serialized Arrow schema in sync atomically."""
        if self._writer is not None and not hasattr(self._writer, "add_key_value_metadata"):
            logger.warning(
                "Cannot add scan_format after batches were written: PyArrow 17 or newer is required; "
                "keeping the existing footer declaration for %s",
                self._path,
            )
            return
        self._file_metadata[b"scan_format"] = scan_format.encode()
        if self._writer is not None:
            self._writer.add_key_value_metadata(
                {
                    b"scan_format": self._file_metadata[b"scan_format"],
                    b"ARROW:schema": base64.b64encode(self.arrow_schema.serialize().to_pybytes()),
                }
            )

    def _check_record_keys(self, records: list[dict]) -> None:
        """Refuse a mistyped key, which would otherwise be an all-null column.

        Keys starting with an underscore are a converter's transient hints.
        """
        allowed = {*self.arrow_schema.names, *self._record_only_keys}
        unknown = {key for key in set().union(*records) if not key.startswith("_")} - allowed
        if unknown:
            raise ValueError(f"{self._schema_class.view_name} records have keys outside the schema: {sorted(unknown)}")

    def _write_batch(self, records: list[dict]):
        self._buffer.extend(records)
        while len(self._buffer) >= self._batch_size:
            batch = self._buffer[: self._batch_size]
            self._buffer = self._buffer[self._batch_size :]
            self._write_arrow_batch(batch)

    def write_table(self, table: pa.Table):
        """Write a complete Arrow table. Validates schema.

        For views with a derived identity, the id column is stamped
        (derived from the identity_composite) before validation, so callers may
        pass tables whose id column is absent or null.
        """
        return self._guard(self._write_table, table)

    def _write_table(self, table: pa.Table):
        table = self._fill_identity_table(table)
        result = self._schema_class.validate_full(table, strict=False)
        if result.errors:
            raise ValueError("Schema validation failed:\n" + "\n".join(issue.message for issue in result.errors))
        self._log_validation_warnings(result.warnings)
        self._ensure_writer()
        self._writer.write_table(table, row_group_size=row_group_rows(table))

    def write_dataframe(self, df: "pd.DataFrame"):
        """Write a pandas DataFrame. Converts to Arrow and validates.

        Built against a schema whose derived-id field is **nullable**, so a frame
        without a precomputed id is accepted; :meth:`write_table` then derives the
        id before validating against the strict schema. Only columns the schema
        declares as nullable — plus the derived id — are backfilled when the frame
        omits them; a genuinely missing *required* column is left absent so
        ``from_pandas`` surfaces it as a clear error.
        """
        return self._guard(self._write_dataframe, df)

    def _write_dataframe(self, df: "pd.DataFrame"):
        id_field = self._id_field
        relaxed = self._relaxed_id_schema()
        missing = [f.name for f in relaxed if f.name not in df.columns and (f.nullable or f.name == id_field)]
        if missing:
            df = df.copy()
            for name in missing:
                df[name] = None
        table = pa.Table.from_pandas(df, schema=relaxed, preserve_index=False)
        self.write_table(table)

    def close(self):
        """Flush buffered rows, close the file, and validate its identity PK.

        A write that already failed is never promoted, however close() is reached.
        """
        self._close(validate=not self._write_failed)

    def _close(self, *, validate: bool) -> None:
        """Close the writer, discarding buffered rows after a body error."""
        completed = False
        writer = None
        try:
            if validate and self._buffer:
                self._write_arrow_batch(self._buffer)
            self._buffer = []
            if validate and self._identity.sources:
                self._add_footer_metadata(self.identity_source_metadata())
            wrote_file = self._writer is not None
            writer = self._writer
            self._writer = None
            if writer:
                writer.close()
                writer = None
            if not wrote_file or not validate:
                # The body raised: leave the destination untouched and discard
                # any partial staging file.
                self._discard_staged()
            else:
                # Validate before the atomic publish so invalid output never
                # reaches the destination.
                self._validate_identity_uniqueness()
                self._validate_referential()
                self._validate_scan_format()
                self._promote_staged()
            completed = True
        finally:
            if not completed:
                self._buffer = []
                self._writer = None
                writer = None
                self._discard_staged()

    @property
    def _validation_path(self) -> Path:
        """Where the bytes currently live: the staging file until it is promoted."""
        return self._staging_path if self._staging_path.exists() else self._path

    def _promote_staged(self) -> None:
        """Atomically move the staged file onto the destination."""
        if not self._staging_path.exists():
            return
        os.replace(self._staging_path, self._path)

    def _discard_staged(self) -> None:
        """Remove the staged file, leaving any existing destination intact."""
        try:
            self._staging_path.unlink(missing_ok=True)
        except OSError as exc:  # pragma: no cover - best effort cleanup
            logger.warning("Could not remove staging file %s: %s", self._staging_path, exc)

    def _validate_referential(self) -> None:
        """Warn on pg referential / run-disjointness problems over the full file.

        The streaming ``write_batch`` path validates each batch's schema but
        cannot see cross-row invariants (run-disjointness and duplicate
        ``grouped_runs`` span the whole table), so ``write_table``'s referential
        checks were effectively off for every streamed pg converter
        (bigbio/qpx#242). They are re-run here at ``close()`` over the complete
        written file, so the record/batch and table paths converge. Per the
        maintainer policy (1.1.2) these are **warnings** on the write/convert
        path, never a raise; ``qpxc validate`` (always strict) stays the strict gate.
        """
        if getattr(self._schema_class, "view_name", None) != "pg":
            return
        issues = pg_referential_issues_from_parquet(str(self._validation_path), self._schema_class.view_name, "warning")
        for issue in issues:
            logger.warning("pg referential check '%s': %s", issue.check, issue.message)

    def _validate_scan_format(self) -> None:
        """Warn about scans the declared format does not allow; a streamed file declares it only at the end."""
        declared = self._file_metadata.get(b"scan_format")
        if self._schema_class.view_name != "psm" or not declared:
            return
        declared = declared.decode()
        with pq.ParquetFile(self._validation_path) as file:
            count = sum(_scan_format_violations(batch.column(0), declared) for batch in file.iter_batches(columns=["scan"]))
        self._log_validation_warnings(_scan_format_issues(declared, count, "psm", "warning"))

    def _validate_identity_uniqueness(self) -> None:
        """Reject duplicate or null identity ids across the complete output file."""
        id_field = self._id_field
        if id_field is None:
            return
        quoted_id = validate_identifier(id_field)
        query = sql_build(
            """
            SELECT
                COUNT(*) AS total_count,
                COUNT(DISTINCT $id_field) AS unique_count,
                COUNT(*) FILTER (WHERE $id_field IS NULL) AS null_count
            FROM read_parquet(?)
            """,
            id_field=quoted_id,
        )
        connection = create_connection()
        try:
            total_count, unique_count, null_count = connection.execute(query, [str(self._validation_path)]).fetchall()[0]
        finally:
            connection.close()

        if null_count:
            raise ValueError(f"Primary key ({id_field}) contains {null_count} null row(s) out of {total_count}")
        if unique_count != total_count:
            duplicate_count = total_count - unique_count
            logger.warning(
                "Primary key (%s) has %d duplicate row(s) (%d unique out of %d). "
                "Duplicate ids are tolerated as a warning while the qpx format stabilises "
                "and usually indicate a converter identity bug.",
                id_field,
                duplicate_count,
                unique_count,
                total_count,
            )

    def _log_validation_warnings(self, issues) -> None:
        """Log each distinct non-fatal finding once per file, or a violation would be published unnoticed."""
        for issue in issues:
            if issue.check not in self._warned_checks:
                self._warned_checks.add(issue.check)
                logger.warning("%s: %s check '%s': %s", self._path.name, self._schema_class.view_name, issue.check, issue.message)

    def _records_problem(self, records: list[dict], schema: pa.Schema) -> str | None:
        """Name the first column and row of *records* that the schema cannot hold."""
        for field in schema:
            try:
                pa.array([record.get(field.name) for record in records], type=field.type)
            except (pa.ArrowInvalid, pa.ArrowTypeError, pa.ArrowNotImplementedError):
                for row, record in enumerate(records):
                    try:
                        pa.array([record.get(field.name)], type=field.type)
                    except (pa.ArrowInvalid, pa.ArrowTypeError, pa.ArrowNotImplementedError) as exc:
                        return f"column '{field.name}', row {row} of the batch, value {record.get(field.name)!r}: {exc}"
        return None

    def _build_table(self, records: list[dict], schema: pa.Schema) -> pa.Table:
        """Arrow table of *records*, naming the column and row of a value the schema cannot hold."""
        try:
            return pa.Table.from_pylist(records, schema=schema)
        except (pa.ArrowInvalid, pa.ArrowTypeError, pa.ArrowNotImplementedError) as exc:
            problem = self._records_problem(records, schema)
            raise ValueError(f"Cannot write {self._schema_class.view_name} records: {problem or exc}") from exc

    def _write_arrow_batch(self, records: list[dict]):
        if self._identity_composite and self._id_field is not None:
            # Build with a nullable id, then derive from the persisted (cast)
            # values via _fill_identity_table — the same path table/DataFrame
            # writes use, so a row's id is identical no matter how it was written.
            self._write_identity_table(self._build_table(records, self._relaxed_id_schema()))
            return
        batch = self._build_table(records, self.arrow_schema).combine_chunks().to_batches()[0]
        self._check_required_non_null(batch)
        self._ensure_writer()
        self._writer.write_batch(batch, row_group_size=row_group_rows(batch))

    def _write_identity_table(self, table: pa.Table) -> None:
        """Derive the ids of a table built on the nullable-id schema, check it and write it."""
        table = self._fill_identity_table(table)
        self._check_required_non_null(table)
        self._ensure_writer()
        self._writer.write_table(table, row_group_size=row_group_rows(table))

    def _ensure_writer(self):
        if self._writer is None:
            self._staging_path.parent.mkdir(parents=True, exist_ok=True)
            self._writer = pq.ParquetWriter(
                str(self._staging_path),
                schema=self.arrow_schema,
                **self._parquet_write_options(),
            )

    def _parquet_write_options(self) -> dict:
        """Build ParquetWriter kwargs: codec level + selective BYTE_STREAM_SPLIT.

        BYTE_STREAM_SPLIT is applied only to the high-entropy float leaves
        present in this writer's schema; dictionary encoding is kept on every
        other column (pyarrow forbids dictionary + explicit encoding on the same
        column). All changes are encoding-only and lossless — the schema and
        values are untouched, so existing readers (pyarrow, DuckDB) read the
        output unchanged.
        """
        return parquet_write_options(self.arrow_schema, self._compression)

    def identified_table(self, table: pa.Table) -> pa.Table:
        """*table* with its id column filled as this writer would write it."""
        return self._fill_identity_table(table)

    @property
    def footer_metadata(self) -> dict[bytes, bytes]:
        """The footer this writer stamps on its files."""
        return dict(self._file_metadata)

    @classmethod
    def view_schema(cls):
        """The schema of the view this writer writes."""
        return cls._schema_class

    @classmethod
    def derive_ids(cls, table: pa.Table, identity_composite: Optional[Sequence[str]] = None) -> pa.Table:
        """Return *table* with its id column derived from the identity composite for the rows that have none.

        Lets a producer know a row's id before writing, for instance to fill
        ``psm.feature_id``. Pass the composite the file will declare.
        """
        return cls("unused", identity_composite=identity_composite, override_provided_ids=False).identified_table(table)

    @classmethod
    def write_partitioned(
        cls,
        table: pa.Table,
        output_dir: str | Path,
        partition_cols: list[str] | None = None,
        compression: str = "zstd",
        file_metadata: dict | None = None,
        mode: str = "error",
        identity_composite: Optional[Sequence[str]] = None,
    ) -> Path:
        """Write Arrow table as Hive-partitioned Parquet.

        Called on a concrete writer subclass (e.g. ``FeatureWriter``) the table
        is routed through the **same schema validation** as the single-file
        writers before it is written, so a partitioned dump no longer bypasses
        those checks (bigbio/qpx#252): the derived id is stamped for rows without
        one, provided ids are kept, structural validation errors are rejected,
        and non-strict findings such as duplicate identities are logged as
        warnings. The footer declares the writer's identity. Called on the
        un-bound ``BaseWriter`` (no schema is known) it stays a raw partitioned
        dump, unchanged.

        Repeated calls with ``mode="append"`` stream a large table into one
        directory, as every call adds files with unique names.

        Encoding parity with the single-file writers: the same
        :func:`parquet_write_options` path is used, so partitioned output gets
        the tuned ZSTD level, selective BYTE_STREAM_SPLIT and dictionary
        encoding, Parquet v2.6, plus a sane row-group size (so compression is
        not defeated by ~80-row groups). The qpx footer KV metadata is stamped
        onto every file so partitioned datasets keep ``qpx_version`` /
        ``file_type`` / ``uuid`` like single files do.

        Parameters
        ----------
        table : pa.Table
            Data to write.
        output_dir : str | Path
            Root directory for partitioned output.
        partition_cols : list[str] | None
            Columns to partition by. Defaults to ["run_file_name"].
        compression : str
            Compression algorithm. Defaults to "zstd".
        file_metadata : dict | None
            Extra footer KV metadata (bytes->bytes) to stamp, e.g. a writer's
            ``_file_metadata`` carrying ``file_type`` / ``uuid``. Merged on top
            of any metadata already on the table schema and a minimal qpx footer.
        mode : str
            With data already in *output_dir*: ``"error"`` (default) raises,
            ``"append"`` adds to it, ``"replace_partitions"`` replaces the
            partitions this table writes and keeps the others.
        identity_composite : Sequence[str] | None
            The composite that ids are derived from and the footer declares;
            defaults to the one in the table's footer metadata, as a table read
            from a QPX file carries, else to the schema's.

        Returns
        -------
        Path
            The output directory.
        """
        from qpx.writers.partitioned import write_partitioned_table

        return write_partitioned_table(
            cls, table, output_dir, partition_cols, compression, file_metadata, mode, identity_composite
        )

    def __enter__(self):
        return self

    def __exit__(self, exc_type, _exc_value, _traceback):
        # Preserve an exception raised inside the context instead of masking it
        # with validation of a partial output file.
        self._close(validate=exc_type is None and not self._write_failed)
