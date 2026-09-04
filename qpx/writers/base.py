"""BaseWriter — batched Parquet writer with schema validation."""

from __future__ import annotations

import datetime
import logging
import os
import uuid
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Sequence

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from qpx._version import __version__
from qpx.core.data.identity import derive_id
from qpx.core.data.schema import pg_referential_issues_from_parquet
from qpx.core.sql import sql_build, validate_identifier
from qpx.version import QPX_SPEC_VERSION

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)


# High-entropy float leaves that compress poorly under the default
# dictionary/PLAIN encoding (≈1.1x with ZSTD on continuous values). Parquet's
# BYTE_STREAM_SPLIT encoding reorders the bytes of IEEE floats so that ZSTD can
# exploit the per-byte regularity, measuring 15-28% smaller per column on real
# data. Paths use pyarrow's dotted leaf convention (``col``, ``a.list.element``,
# nested struct fields joined by ``.``). Only leaves present in a given schema
# are touched, so this list can name columns from any QPX view.
#
# NOTE: this is deliberately NOT applied to low-cardinality / dictionary-friendly
# float columns (e.g. ``additional_scores.score_value``,
# ``posterior_error_probability``), where BYTE_STREAM_SPLIT measured *larger*
# than the default dictionary encoding. Those keep their default encoding.
#
# PEP was a candidate for float32 + BYTE_STREAM_SPLIT (audit projected ~52%),
# but empirical validation over the real qpx corpus found genuine PEP values as
# small as ~1e-300..1e-318 (DIA/quantms feature+psm files): these underflow to
# zero in float32 and would corrupt any downstream -log10(PEP)/q-value use. The
# float32 narrowing was therefore held (PEP stays float64), so it is not split.
#
# The intensity leaves (``intensities.list.element.intensity`` and the nested
# ``additional_intensities…intensity_value``) are intentionally NOT split.
# BYTE_STREAM_SPLIT helps continuous LFQ intensities (~-20%) but slightly *hurts*
# TMT reporter-ion columns (+4-7%), which dominate TMT feature files (~64% of
# the file). Rather than probe the data at write time, we simply omit them here:
# they fall back to dictionary encoding, which is safe everywhere (never worse
# than PLAIN) and avoids the TMT regression. RT/mz leaves (a proven win) stay.
_BYTE_STREAM_SPLIT_LEAVES: tuple[str, ...] = (
    "predicted_rt",
    "rt",
    "rt_start",
    "rt_stop",
    "calculated_mz",
    "observed_mz",
)

# ZSTD level 9 (vs the pyarrow default of 3) buys a few extra percent at modest
# write-cost; only codecs that accept a level get one.
_COMPRESSION_LEVELS: dict[str, int] = {"zstd": 9, "gzip": 9}

# Row-group sizing for partitioned (dataset) writes. Left unset, pyarrow's
# dataset writer flushes a group per input batch, which for streamed msnet
# conversion fragments into ~80-row groups and collapses ZSTD (1.12x vs 1.75x).
# A 200k-1M window matches the single-file writer's effective grouping.
_PARTITIONED_MIN_ROWS_PER_GROUP: int = 200_000
_PARTITIONED_MAX_ROWS_PER_GROUP: int = 1_000_000


def _stamp_footer_metadata(
    table: pa.Table,
    compression: str,
    file_metadata: dict | None = None,
) -> pa.Table:
    """Return *table* with qpx footer KV metadata on its schema.

    Layering: any metadata already on the table schema, then defaults that a
    caller may override (``creation_date``/``uuid``), then *file_metadata* (e.g.
    a writer's ``_file_metadata`` carrying ``file_type``/``uuid``), and finally
    the **canonical identity** (``qpx_version``/``writer_version``/
    ``compression_format``) which always wins so a stale caller value cannot
    mis-stamp the file. Used so partitioned datasets carry the same footer
    identity as single files.
    """
    merged: dict[bytes, bytes] = dict(table.schema.metadata or {})
    # Overridable defaults.
    merged.update(
        {
            b"creation_date": datetime.datetime.now().isoformat().encode(),
            b"uuid": str(uuid.uuid4()).encode(),
        }
    )
    if file_metadata:
        merged.update(
            {
                (k if isinstance(k, bytes) else str(k).encode()): (v if isinstance(v, bytes) else str(v).encode())
                for k, v in file_metadata.items()
            }
        )
    # Canonical identity — reserved, must not be overridden by caller metadata.
    merged.update(
        {
            b"qpx_version": QPX_SPEC_VERSION.encode(),
            b"writer_version": __version__.encode(),
            b"compression_format": str(compression).encode(),
        }
    )
    return table.replace_schema_metadata(merged)


def _validate_partitioned_table(schema_class, table: pa.Table) -> None:
    """Reject structural errors and log non-strict partitioned-write findings."""
    result = schema_class.validate_full(table, strict=False)
    if result.errors:
        raise ValueError("Schema validation failed:\n" + "\n".join(issue.message for issue in result.errors))
    for issue in result.warnings:
        logger.warning(
            "%s partitioned-write validation '%s': %s",
            schema_class.view_name,
            issue.check,
            issue.message,
        )


def _schema_leaf_paths(schema: pa.Schema) -> list[str]:
    """Return every leaf column path in pyarrow's dotted convention."""

    def _walk(field: pa.Field, prefix: str) -> list[str]:
        dtype = field.type
        if pa.types.is_struct(dtype):
            out: list[str] = []
            for i in range(dtype.num_fields):
                child = dtype.field(i)
                out.extend(_walk(child, f"{prefix}.{child.name}"))
            return out
        if pa.types.is_list(dtype) or pa.types.is_large_list(dtype):
            return _walk(dtype.value_field, f"{prefix}.list.element")
        return [prefix]

    leaves: list[str] = []
    for top in schema:
        leaves.extend(_walk(top, top.name))
    return leaves


def parquet_write_options(schema: pa.Schema, compression: str) -> dict:
    """Build QPX ParquetWriter options for *schema* and *compression*."""
    options: dict = {"compression": compression}
    level = _COMPRESSION_LEVELS.get(str(compression).lower())
    if level is not None:
        options["compression_level"] = level

    leaves = _schema_leaf_paths(schema)
    bss = {leaf: "BYTE_STREAM_SPLIT" for leaf in _BYTE_STREAM_SPLIT_LEAVES if leaf in leaves}
    # Always stamp the explicit v2.6 format and dictionary-encode non-BSS leaves so
    # every artifact (feature/psm and pg alike) uses the same encoding, rather than
    # falling back to pyarrow defaults when no byte-stream-split leaf is present.
    options["version"] = "2.6"
    options["use_dictionary"] = [leaf for leaf in leaves if leaf not in bss]
    if bss:
        options["column_encoding"] = bss
    return options


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

    def __init__(
        self,
        path: str | Path,
        creator: str = "qpx",
        software_provider: str = "qpx",
        compression: str = "zstd",
        batch_size: int = 1_000_000,
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
        self._extra_columns = extra_columns or []
        # Override mode (e.g. consensusXML): re-derive the id even when the
        # producer supplied one, stashing the original as a cv_param so it is
        # preserved/traceable. ``overridden_id_count`` lets the converter record
        # the override in provenance.
        self._override_provided_ids = override_provided_ids
        self.overridden_id_count = 0
        # Producer-specific identity composite override. When set, the writer
        # hashes THIS composite into the derived id (feature_id / psm_id / pg_id)
        # instead of the schema default, letting each converter declare the
        # measured keys its producer emits (bigbio/qpx#229). Stored as a tuple so
        # the footer stamps the effective composite via the ``_identity_composite``
        # property below.
        self._identity_composite_override = tuple(identity_composite) if identity_composite else None
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

    # --- Derived identity --------------------------------------------------

    @property
    def _identity_composite(self) -> tuple[str, ...] | None:
        """The effective identity_composite the id is derived from.

        A converter-supplied override (producer-specific composite) wins over the
        schema class attribute; ``None`` for views without any derived identity.
        """
        if self._identity_composite_override is not None:
            return self._identity_composite_override
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
        # Identity list fields whose element order is not meaningful: grouped_runs
        # is a set of raw files, and pg_accessions is a set of group members (the
        # anchor/leader is carried separately). Hashing them order-independently
        # keeps the derived id stable across producers that emit the same
        # membership in a different order.
        unordered_list_indices = tuple(
            index for index, field in enumerate(composite) if field in ("grouped_runs", "pg_accessions")
        )
        for index, provided in enumerate(existing):
            if provided is not None and not self._should_override_provided_id(index, composite_values):
                ids.append(provided)
                continue
            if provided is not None:
                if cv_lists is not None:
                    params = list(cv_lists[index] or [])
                    params.append({"cv_name": f"provided_{id_field}", "cv_value": str(provided)})
                    cv_lists[index] = params
                self.overridden_id_count += 1
            ids.append(
                derive_id(
                    [composite_values[field][index] for field in composite],
                    unordered_list_indices=unordered_list_indices,
                )
            )
        return ids

    def _should_override_provided_id(self, _index: int, _composite_values: dict[str, list]) -> bool:
        """Return whether a producer id should be replaced for one row."""
        return self._override_provided_ids

    def _fill_identity_table(self, table: pa.Table) -> pa.Table:
        """Return *table* with its derived-id column present and non-null.

        The table is assumed to already match ``arrow_schema`` in every other
        column (callers build via ``arrow_schema`` or align first). Rows with a
        non-null id keep it; the rest are derived from the composite columns.
        In ``override_provided_ids`` mode the id is re-derived even when present,
        stashing the original as a ``provided_<id>`` cv_param. This is the single
        derivation point for every write path (record, table, DataFrame).
        """
        composite = self._identity_composite
        id_field = self._id_field
        if not composite or id_field is None:
            return table
        names = set(table.schema.names)
        n = table.num_rows
        # Hash the PERSISTED (schema-cast) value of each composite field so the id
        # is identical regardless of which write path built the table — e.g. a
        # float32 ``rt`` hashes to the same value whether it arrived as a raw
        # Python float (record path) or an already-cast float32 column (table path).
        composite_values = self._persisted_composite_values(table, composite)
        existing = table.column(id_field).to_pylist() if id_field in names else [None] * n
        override = self._override_provided_ids
        if override and "cv_params" not in names and any(provided is not None for provided in existing):
            raise ValueError(f"override_provided_ids requires a 'cv_params' column to preserve provided {id_field} values")
        cv_lists = table.column("cv_params").to_pylist() if (override and "cv_params" in names) else None
        ids = self._identity_values(existing, composite_values, composite, cv_lists, id_field)
        if cv_lists is not None:
            cv_type = table.schema.field("cv_params").type
            table = table.set_column(
                table.schema.get_field_index("cv_params"),
                table.schema.field("cv_params"),
                pa.array(cv_lists, type=cv_type),
            )
        id_array = pa.array(ids, type=pa.int64())
        id_arrow_field = pa.field(id_field, pa.int64(), nullable=False)
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
        columns = [table.column(f.name) if f.name in present else pa.nulls(n, type=f.type) for f in relaxed]
        intermediate = pa.Table.from_arrays(
            columns,
            schema=pa.schema([pa.field(f.name, columns[i].type) for i, f in enumerate(relaxed)]),
        )
        return intermediate.cast(relaxed)

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
        return self._guard(self._write_batch, records)

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
        errors = self._schema_class.validate(table, strict=False)
        if errors:
            raise ValueError("Schema validation failed:\n" + "\n".join(errors))
        self._ensure_writer()
        self._writer.write_table(table)

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
        path, never a raise; ``qpxc validate --strict`` stays the strict gate.
        """
        if getattr(self._schema_class, "view_name", None) != "pg":
            return
        issues = pg_referential_issues_from_parquet(str(self._validation_path), self._schema_class.view_name, "warning")
        for issue in issues:
            logger.warning("pg referential check '%s': %s", issue.check, issue.message)

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
        connection = duckdb.connect()
        try:
            total_count, unique_count, null_count = connection.execute(query, [str(self._validation_path)]).fetchone()
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

    def _write_arrow_batch(self, records: list[dict]):
        if self._identity_composite and self._id_field is not None:
            # Build with a nullable id, then derive from the persisted (cast)
            # values via _fill_identity_table — the same path table/DataFrame
            # writes use, so a row's id is identical no matter how it was written.
            table = self._fill_identity_table(pa.Table.from_pylist(records, schema=self._relaxed_id_schema()))
            self._check_required_non_null(table)
            self._ensure_writer()
            self._writer.write_table(table)
            return
        batch = pa.RecordBatch.from_pylist(records, schema=self.arrow_schema)
        self._check_required_non_null(batch)
        self._ensure_writer()
        self._writer.write_batch(batch)

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

    @classmethod
    def write_partitioned(
        cls,
        table: pa.Table,
        output_dir: str | Path,
        partition_cols: list[str] | None = None,
        compression: str = "zstd",
        file_metadata: dict | None = None,
    ) -> Path:
        """Write Arrow table as Hive-partitioned Parquet.

        Called on a concrete writer subclass (e.g. ``FeatureWriter``) the table
        is routed through the **same id-derivation and schema validation** as
        the single-file writers before it is written, so a partitioned dump no
        longer bypasses those checks (bigbio/qpx#252): the derived id is stamped
        when absent, structural validation errors are rejected, and non-strict
        findings such as duplicate identities are logged as warnings. Called on
        the un-bound ``BaseWriter`` (no schema is known) it stays a raw
        partitioned dump, unchanged.

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

        Returns
        -------
        Path
            The output directory.
        """
        import pyarrow.dataset as ds

        output_dir = Path(output_dir)
        cols = partition_cols or ["run_file_name"]

        # When called on a concrete writer subclass, route through the same
        # id-derivation + validation the single-file writers use so partitioned
        # output is not a validation-bypass (bigbio/qpx#252). An instance is
        # built only to reuse _fill_identity_table / the schema; no file is
        # opened until write_dataset below. On the un-bound BaseWriter
        # (_schema_class is None) this is skipped and the raw dump is unchanged.
        schema_class = cls._schema_class
        if schema_class is not None:
            deriver = cls(output_dir, compression=compression)
            table = deriver._fill_identity_table(table)
            _validate_partitioned_table(schema_class, table)

        # Validate the partition columns up front so the default path (or an
        # unsupported key) fails with an actionable message instead of a raw
        # KeyError deep in pyarrow. Since QPX 1.1 the pg view is keyed on
        # grouped_runs (list<string>) with no scalar run_file_name, and Hive
        # partitioning encodes each value as a directory name — it cannot
        # represent a list — so the pg view is not Hive-partitionable.
        available = set(table.schema.names)
        missing = [c for c in cols if c not in available]
        if missing:
            if "grouped_runs" in available and "run_file_name" in missing:
                raise ValueError(
                    "The protein-group (pg) view is not Hive-partitionable: it is keyed on "
                    "'grouped_runs' (list<string>, since QPX 1.1) and has no scalar "
                    "'run_file_name' column. Hive partitioning encodes each value as a "
                    "directory name and cannot represent a list. Write pg as a single "
                    "Parquet file (write_table), or partition the psm/feature views, which "
                    "keep 'run_file_name'."
                )
            raise ValueError(f"Cannot partition on missing column(s) {missing}; available columns: {sorted(available)}")
        for col in cols:
            col_type = table.schema.field(col).type
            if pa.types.is_list(col_type) or pa.types.is_large_list(col_type) or pa.types.is_struct(col_type):
                raise ValueError(
                    f"Cannot Hive-partition on non-scalar column '{col}' (type {col_type}); partition "
                    "columns must be scalar. The pg view's 'grouped_runs' is a list, so pg is not "
                    "partitionable — write it as a single file."
                )

        table = _stamp_footer_metadata(table, compression, file_metadata)

        part_schema = pa.schema([table.schema.field(c) for c in cols])
        partitioning = ds.partitioning(part_schema, flavor="hive")

        opts = parquet_write_options(table.schema, compression)
        file_options = ds.ParquetFileFormat().make_write_options(**opts)
        ds.write_dataset(
            table,
            str(output_dir),
            format="parquet",
            partitioning=partitioning,
            existing_data_behavior="overwrite_or_ignore",
            file_options=file_options,
            min_rows_per_group=_PARTITIONED_MIN_ROWS_PER_GROUP,
            max_rows_per_group=_PARTITIONED_MAX_ROWS_PER_GROUP,
        )
        return output_dir

    def __enter__(self):
        return self

    def __exit__(self, exc_type, _exc_value, _traceback):
        # Preserve an exception raised inside the context instead of masking it
        # with validation of a partial output file.
        self._close(validate=exc_type is None and not self._write_failed)
