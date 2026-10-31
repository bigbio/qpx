"""Hive-partitioned writes of a QPX view."""

from __future__ import annotations

import datetime
import logging
import uuid
from pathlib import Path
from typing import Optional, Sequence

import pyarrow as pa
import pyarrow.dataset as ds

from qpx._version import __version__
from qpx.version import QPX_SPEC_VERSION
from qpx.writers.encoding import parquet_write_options

logger = logging.getLogger(__name__)

# Row-group sizing for partitioned (dataset) writes. Left unset, pyarrow's
# dataset writer flushes a group per input batch, which for streamed msnet
# conversion fragments into ~80-row groups and collapses ZSTD (1.12x vs 1.75x).
# A 200k-1M window matches the single-file writer's effective grouping.
_PARTITIONED_MIN_ROWS_PER_GROUP: int = 200_000
_PARTITIONED_MAX_ROWS_PER_GROUP: int = 1_000_000
# write_partitioned mode -> pyarrow behaviour; "error" is checked up front with a clear message.
_PARTITIONED_MODES = {"error": "overwrite_or_ignore", "append": "overwrite_or_ignore", "replace_partitions": "delete_matching"}


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


def _validate_partition_columns(table: pa.Table, cols: list[str]) -> None:
    """Require existing scalar columns for Hive partitioning."""
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


def write_partitioned_table(
    writer_class,
    table: pa.Table,
    output_dir: str | Path,
    partition_cols: list[str] | None,
    compression: str,
    file_metadata: dict | None,
    mode: str,
    identity_composite: Optional[Sequence[str]],
) -> Path:
    """Write *table* under *output_dir*; see :meth:`BaseWriter.write_partitioned`."""
    output_dir = Path(output_dir)
    cols = partition_cols or ["run_file_name"]
    if mode not in _PARTITIONED_MODES:
        raise ValueError(f"mode must be one of {sorted(_PARTITIONED_MODES)}, got {mode!r}")
    if mode == "error" and output_dir.is_dir() and next(output_dir.rglob("*.parquet"), None):
        raise FileExistsError(f"{output_dir} already holds Parquet files; pass mode='append' or mode='replace_partitions'")

    # When called on a concrete writer subclass, route through the same
    # id-derivation + validation the single-file writers use so partitioned
    # output is not a validation-bypass (bigbio/qpx#252). An instance is
    # built only to reuse its id derivation / the schema; no file is
    # opened until write_dataset below. On the un-bound BaseWriter
    # (no view schema) this is skipped and the raw dump is unchanged.
    schema_class = writer_class.view_schema()
    if schema_class is not None:
        declared = (table.schema.metadata or {}).get(b"identity_composite")
        composite = identity_composite or (declared.decode().split(",") if declared else None)
        deriver = writer_class(output_dir, compression=compression, override_provided_ids=False, identity_composite=composite)
        table = deriver.identified_table(table)
        _validate_partitioned_table(schema_class, table)
        footer = {k: v for k, v in deriver.footer_metadata.items() if k not in (b"creator", b"software_provider")}
        footer |= deriver.identity_source_metadata()
        file_metadata = {**footer, **(file_metadata or {})}

    _validate_partition_columns(table, cols)

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
        basename_template=f"part-{uuid.uuid4().hex[:8]}-{{i}}.parquet",
        existing_data_behavior=_PARTITIONED_MODES[mode],
        file_options=file_options,
        min_rows_per_group=_PARTITIONED_MIN_ROWS_PER_GROUP,
        max_rows_per_group=_PARTITIONED_MAX_ROWS_PER_GROUP,
    )
    return output_dir
