"""Bring a QPX dataset from an older specification to the current one, without its source files."""

from __future__ import annotations

import logging
import shutil
import uuid
from pathlib import Path

import pyarrow.parquet as pq

from qpx.converters.channel_labels import normalize_label
from qpx.converters.ontology_util import dedupe_ontology_entries, reconcile_ontology_entries
from qpx.converters.ptm import normalize_peptidoform
from qpx.dataset import Dataset
from qpx.version import QPX_SPEC_VERSION

logger = logging.getLogger(__name__)

_BATCH_ROWS = 50_000
# Derived again from the rewritten rows, and computed links that are not stored any more.
_DROPPED_ID_COLUMNS = ("feature_id", "psm_id", "pg_id", "psm_ids", "pg_ids")
# The feature composite each producer declares (docs/spec/identity.md), by footer `creator`.
_PRODUCER_FEATURE_COMPOSITES = {
    "diann": ("peptidoform", "charge", "run_file_name"),
    "spectronaut": ("peptidoform", "charge", "run_file_name"),
    "cdap": ("peptidoform", "charge", "run_file_name"),
    "maxquant": ("peptidoform", "charge", "run_file_name", "rt_start", "rt_stop", "scan"),
    "quantms-msstats": ("run_file_name", "peptidoform", "charge", "rt", "scan"),
    "fragpipe": ("quantification_unit_id", "peptidoform", "charge", "compensation_voltage"),
}
# Non-conformant label-free labels of 1.0.2: an erratum, not an older vocabulary.
_LABEL_ERRATA = {"raw": "LFQ", "precursor_area": "LFQ"}
_FILE_LEVEL_DATASET_COLUMNS = ("file_checksums", "file_row_counts", "file_sizes_bytes", "total_structures", "packaged_at")
_NOT_MIGRATED_SUFFIXES = (".h5mu", ".h5ad")


def _label(value):
    return value if value is None else normalize_label(_LABEL_ERRATA.get(value, value))


def _entries(values):
    return values and [{**entry, "label": _label(entry["label"])} for entry in values]


def _footer(path: Path) -> dict[bytes, bytes]:
    return pq.read_schema(path).metadata or {}


def _batches(path: Path):
    source = pq.ParquetFile(path)
    try:
        for batch in source.iter_batches(_BATCH_ROWS):
            yield batch.to_pylist()
    finally:
        source.close()


def _keep_schema_fields(writer, records: list[dict]) -> list[dict]:
    """Drop columns the current schema no longer has, naming them once."""
    names = set(writer.arrow_schema.names)
    dropped = {key for record in records for key in record} - names
    if dropped:
        logger.warning("Dropping columns outside the schema: %s", sorted(dropped))
    return [{key: value for key, value in record.items() if key in names} for record in records]


def _identity_row(record: dict) -> dict:
    """The rewritten peptidoform and labels of a psm or feature row, without its old id."""
    row = {key: value for key, value in record.items() if key not in _DROPPED_ID_COLUMNS}
    if record.get("feature_id") is not None and "psm_id" in record:
        row["feature_id"] = record["feature_id"]
    if row.get("peptidoform"):
        row["peptidoform"] = normalize_peptidoform(row["peptidoform"])
    for column in ("calculated_mz", "observed_mz"):
        # 0.0 was the old sentinel for an m/z the producer did not report
        if row.get(column) == 0.0:
            row[column] = None
    for column in ("intensities", "additional_intensities"):
        if column in row:
            row[column] = _entries(row[column])
    if row.get("label") is not None:
        row["label"] = _label(row["label"])
    return row


def _flat_pg_rows(record: dict) -> list[dict]:
    """One pg row per label of a pre-1.1 row, keyed by a one-run grouped_runs."""
    row = {key: value for key, value in record.items() if key not in _DROPPED_ID_COLUMNS}
    if "run_file_name" in row:
        row["grouped_runs"] = [row.pop("run_file_name")]
    if "intensities" not in row:
        if row.get("label") is not None:
            row["label"] = _label(row["label"])
        return [row]
    entries = row.pop("intensities")
    additional = row.pop("additional_intensities", None) or []
    if not entries:
        return [{**row, "label": None, "intensity": None, "additional_intensities": None}]
    return [
        {
            **row,
            "label": _label(entry["label"]),
            "intensity": entry["intensity"],
            "additional_intensities": [
                {**item, "label": _label(item["label"])} for item in additional if item["label"] == entry["label"]
            ]
            or None,
        }
        for entry in entries
    ]


def _run_row(record: dict) -> dict:
    return {**record, "samples": _entries(record.get("samples"))}


def _dataset_row(record: dict) -> dict:
    """File-level columns describe the old files, and this is a new analysis id."""
    row = {key: value for key, value in record.items() if key not in _FILE_LEVEL_DATASET_COLUMNS}
    return {**row, "analysis_id": str(uuid.uuid4())}


def _mz_row(record: dict) -> dict:
    scan = record.get("scan")
    return {**record, "scan": [scan] if isinstance(scan, int) else scan}


_ROW_MIGRATIONS = {
    "feature": lambda record: [_identity_row(record)],
    "psm": lambda record: [_identity_row(record)],
    "pg": _flat_pg_rows,
    "run": lambda record: [_run_row(record)],
    "dataset": lambda record: [_dataset_row(record)],
    "mz": lambda record: [_mz_row(record)],
}


def _feature_id_map(source: Path, target: Path) -> dict[int, int]:
    """Old id to new id, by row position: the writer keeps the order of its rows."""
    if "feature_id" not in pq.read_schema(source).names:
        return {}
    old = pq.read_table(source, columns=["feature_id"]).column(0).to_pylist()
    new = pq.read_table(target, columns=["feature_id"]).column(0).to_pylist()
    return dict(zip(old, new))


def _writer_options(view: str, footer: dict[bytes, bytes]) -> dict:
    creator = footer.get(b"creator", b"qpx").decode()
    options = {"creator": creator}
    if b"scan_format" in footer:
        options["scan_format"] = footer[b"scan_format"].decode()
    if view == "feature":
        declared = footer.get(b"identity_composite")
        composite = tuple(declared.decode().split(",")) if declared else _PRODUCER_FEATURE_COMPOSITES.get(creator)
        if composite is not None:
            options["identity_composite"] = composite
    return options


def _migrate_ontology(source: Path, target: Path) -> int:
    """One row per field and view, naming only fields the migrated views emit."""
    writer_class, suffix = Dataset.writer_for("ontology")
    prefix = target.name.removesuffix(suffix)
    records = pq.read_table(source).to_pylist()
    entries = reconcile_ontology_entries(dedupe_ontology_entries(records), target.parent, prefix)
    with writer_class(target, **_writer_options("ontology", _footer(source))) as writer:
        writer.write_batch(entries)
    return len(entries)


def _migrate_view(view: str, source: Path, target: Path, feature_ids: dict[int, int]) -> int:
    if view == "ontology":
        return _migrate_ontology(source, target)
    writer_class, _ = Dataset.writer_for(view)
    migrate_row = _ROW_MIGRATIONS.get(view, lambda record: [record])
    keep_ids = view == "psm" and _footer(source).get(b"identity_source") == b"provided"
    rows = 0
    with writer_class(target, **_writer_options(view, _footer(source))) as writer:
        for records in _batches(source):
            migrated = [row for record in records for row in migrate_row(record)]
            if keep_ids:
                for row, record in zip(migrated, records):
                    row["psm_id"] = record["psm_id"]
            if view == "psm":
                for row in migrated:
                    if row.get("feature_id") is not None:
                        row["feature_id"] = feature_ids.get(row["feature_id"])
            writer.write_batch(_keep_schema_fields(writer, migrated))
            rows += len(migrated)
    return rows


def _append_migration_step(target: Path, source_version: str) -> None:
    """Record the migration as the last provenance step."""
    from qpx._version import __version__
    from qpx.writers.provenance import ProvenanceWriter

    rows = pq.read_table(target).to_pylist()
    step = {
        "step_order": len(rows) + 1,
        "step_category": "format_conversion",
        "step_name": "qpxc migrate",
        "tool_name": "qpx",
        "tool_version": __version__,
        "tool_uri": None,
        "parameters": [{"key": "source_qpx_version", "value": source_version}],
        "config": None,
        "output_views": None,
    }
    with ProvenanceWriter(target) as writer:
        writer.write_batch([*rows, step])


def migrate_dataset(source: Path, output: Path) -> dict[str, int]:
    """Write *source* under the current specification to the new folder *output*.

    Returns the rows written per view. Views keep their file names. Ids are derived
    again, `psm.feature_id` follows the features, and labels and peptidoforms take
    their canonical spelling. A fixed modification the old file left out cannot be
    recovered: re-convert from the source files to get one.
    """
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    names = Dataset.structure_names()
    first, last = ("feature", "psm"), ("ontology",)
    order = [*first, *[name for name in names if name not in (*first, *last)], *last]
    written: dict[str, int] = {}
    feature_ids: dict[int, int] = {}
    migrated_files = set()
    source_version = "unknown"
    for view in order:
        for path in sorted(source.glob(f"*{Dataset.writer_for(view)[1]}")):
            source_version = _footer(path).get(b"qpx_version", source_version.encode()).decode()
            target = output / path.name
            written[view] = written.get(view, 0) + _migrate_view(view, path, target, feature_ids)
            migrated_files.add(path)
            if view == "feature":
                feature_ids.update(_feature_id_map(path, target))
    for path in sorted(source.iterdir()):
        if path in migrated_files or not path.is_file():
            continue
        if path.suffix in _NOT_MIGRATED_SUFFIXES:
            logger.warning("%s is not migrated: build it again from the migrated views", path.name)
        else:
            shutil.copy2(path, output / path.name)
    for path in sorted(output.glob(f"*{Dataset.writer_for('provenance')[1]}")):
        _append_migration_step(path, source_version)
    logger.info("Migrated %s from QPX %s to %s", source, source_version, QPX_SPEC_VERSION)
    return written
