"""Streaming and in-memory consensusXML conversion drivers."""

from __future__ import annotations

import logging
import os
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pyarrow as pa
import pyarrow.parquet as pq

from qpx.converters.openms_consensus.feature_adapter import (
    check_channels_vs_sdrf,
    column_runs,
    identification_context,
    load_consensus_map,
)
from qpx.converters.openms_consensus.feature_dedup import FEATURE_IDENTITY_COMPOSITE as _FEATURE_IDENTITY_COMPOSITE
from qpx.converters.openms_consensus.feature_dedup import FeatureDeduplicator as _FeatureDeduplicator
from qpx.converters.openms_consensus.pg_adapter import (
    ProteinGroupAccumulator,
    accumulate_cf_maps,
    accumulate_consensus_map,
    accumulate_unassigned_maps,
)
from qpx.converters.openms_consensus.psm_adapter import (
    collect_psm_scan_formats,
    declare_psm_scan_format,
    psm_records_for_pid,
    uniform_scan_format,
)
from qpx.converters.openms_consensus.record_linking import _cf_feature_psm_records, _PsmKeys, _RecordContext
from qpx.writers.base import parquet_write_options
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter

_log = logging.getLogger("qpx.converters.openms_consensus.converter")
_PSM_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name", "scan")
# Pending record dicts cost far more than their row groups: keep few of them.
_STREAM_BATCH = 20_000


@dataclass(frozen=True)
class _OutputConfig:
    """Destination and writer metadata shared by every exported view."""

    folder: Path
    prefix: str
    creator: str
    compression: str


@dataclass(frozen=True)
class _ConversionConfig:
    """Requested views and identification/quantification settings."""

    output: _OutputConfig
    structures: tuple[str, ...]
    sdrf_path: Optional[str]
    pg_top: int
    include_unassigned_psms: bool


def _patch_record_columns(table, patches, offset):
    """Apply sparse record patches to one row group without changing other columns."""
    local = {i - offset: rec for i, rec in patches.items() if offset <= i < offset + table.num_rows}
    for name in {name for rec in local.values() for name in rec}:
        col = table.schema.get_field_index(name)
        field = table.schema.field(col)
        values = table.column(col).to_pylist()
        for row, record in local.items():
            if name in record:
                values[row] = record[name]
        table = table.set_column(col, field, pa.array(values, type=field.type))
    return table


def _patch_feature_ids(table, feature_ids, remap_feature_id, offset):
    """Apply PSM link replacements, then resolve links to the final Feature IDs."""
    patched = {i: feature_ids[offset + i] for i in range(table.num_rows) if offset + i in feature_ids}
    if not patched and remap_feature_id is None:
        return table
    col = table.schema.get_field_index("feature_id")
    field = table.schema.field(col)
    values = table.column(col).to_pylist()
    for i, fid in patched.items():
        values[i] = fid
    if remap_feature_id is not None:
        values = [None if fid is None else remap_feature_id(fid) for fid in values]
    return table.set_column(col, field, pa.array(values, type=field.type))


def _rewrite_parquet_rows(
    path: Path,
    rows: set[int],
    compression: str,
    feature_ids: Optional[dict] = None,
    remap_feature_id=None,
    record_patches=None,
) -> None:
    """Rewrite a Parquet file without the rows at the given positions, one row group at a time.

    ``feature_ids`` (output position -> feature_id) replaces the ``feature_id`` of those rows;
    ``remap_feature_id`` (callable) then maps every remaining ``feature_id`` value.
    """
    feature_ids = feature_ids or {}
    tmp = path.with_name(f"{path.name}.pruned")
    source = pq.ParquetFile(str(path))
    try:
        schema = source.schema_arrow
        with pq.ParquetWriter(str(tmp), schema=schema, **parquet_write_options(schema, compression)) as writer:
            offset = 0
            for index in range(source.num_row_groups):
                table = source.read_row_group(index)
                n = table.num_rows
                table = _patch_record_columns(table, record_patches or {}, offset)
                table = _patch_feature_ids(table, feature_ids, remap_feature_id, offset)
                if rows:
                    table = table.filter(pa.array([offset + i not in rows for i in range(n)]))
                writer.write_table(table)
                offset += n
    finally:
        source.close()
    os.replace(tmp, path)


class _RowPruning:  # pylint: disable=too-few-public-methods
    """Writer mixin: drop ``drop_rows`` and apply ``feature_id_patch`` (output positions) before validation."""

    drop_rows: frozenset = frozenset()
    feature_id_patch: dict = {}
    remap_feature_id = None
    record_patches: dict = {}

    def _validate_identity_uniqueness(self) -> None:
        # pylint: disable=no-member  # mixed into FeatureWriter / PsmWriter, which define these
        if self.drop_rows or self.feature_id_patch or self.remap_feature_id is not None or self.record_patches:
            _rewrite_parquet_rows(
                self._validation_path,
                set(self.drop_rows),
                self._compression,
                self.feature_id_patch,
                self.remap_feature_id,
                self.record_patches,
            )
        super()._validate_identity_uniqueness()


def _prunable_writer(writer_cls):
    """Build a pruning writer from the currently configured view writer class.

    Creating the small mixin class at use time keeps test and downstream writer
    substitutions effective while retaining the row-pruning hook required by the
    streaming deduplication path.
    """
    return type(f"_Prunable{writer_cls.__name__}", (_RowPruning, writer_cls), {})


def _claim_runs(run_owner: dict[str, str], path: str, cm) -> None:
    """Record this input's runs, rejecting any run an earlier consensusXML input already contributed."""
    runs = set(column_runs(cm).values())
    for run in sorted(runs):
        if run in run_owner:
            raise ValueError(f"Run {run!r} appears in both {run_owner[run]} and {path}; each run must come from one consensusXML")
    run_owner.update(dict.fromkeys(runs, path))


def _warn_channel_mismatch(cm, sdrf_path) -> None:
    if sdrf_path:
        for msg in check_channels_vs_sdrf(cm, sdrf_path):
            _log.warning("consensusXML/SDRF channel mismatch: %s", msg)


def _write_view(writer_cls, path, records, *, creator, compression, identity_composite=None, scan_format=None):
    """Write records via a view writer; empty input does not create a file."""
    kwargs = {"creator": creator, "compression": compression, "scan_format": scan_format}
    if identity_composite is not None:
        kwargs["identity_composite"] = identity_composite
    with writer_cls(str(path), **kwargs) as w:
        if records:
            w.write_batch(records)
    return path


def _flush_records(writer, records: list[dict], min_size: int = 1) -> list[dict]:
    """Write a complete batch or a nonempty final remainder, returning the pending records."""
    if writer is not None and len(records) >= min_size:
        writer.write_batch(records)
        return []
    return records


def _finalize_in_memory_records(feat_recs, psm_recs, seen, dedup):
    """Consolidate evidence and resolve PSM links before deriving final Feature IDs."""
    for row, fid in seen.feature_id_patch.items():
        psm_recs[row]["feature_id"] = fid
    if dedup.redirect:
        for rec in psm_recs:
            if rec.get("feature_id") is not None:
                rec["feature_id"] = dedup.resolve(rec["feature_id"])
    for row, patch in dedup.row_patches.items():
        if feat_recs:
            feat_recs[row].update(patch)
            # The writer derives the final ID; it is not a producer-supplied ID.
            feat_recs[row].pop("feature_id")
    feat_recs = [rec for row, rec in enumerate(feat_recs) if row not in dedup.superseded]
    return feat_recs


def _write_pg(pg, config, written) -> None:
    """Write accumulated protein groups when they contain exportable rows."""
    if pg is not None:
        records = pg.build(config.sdrf_path, config.pg_top)
        if records:
            output = config.output
            written["pg"] = _write_view(
                PgWriter,
                output.folder / f"{output.prefix}.pg.parquet",
                records,
                creator=output.creator,
                compression=output.compression,
            )


def _stream_feature_psm(cm, writers, context, pg, config, scan_formats: set[str | None], *, batch=_STREAM_BATCH) -> int:
    """Write one ordered element/unassigned pass in batches and accumulate PG evidence."""
    fw, pw = writers
    maps = pg.maps if pg is not None else None
    map_run = column_runs(cm)
    seen = context.seen
    feat_buf: list[dict] = []
    psm_buf: list[dict] = []
    feature_count = 0
    track = isinstance(seen, _PsmKeys)
    removed_index, confidence = identification_context(cm)
    for kind, obj in cm.iter_all():
        if kind == "element":
            if track:
                seen.unassigned = seen.before_features = False
            cf_feats, cf_psms = _cf_feature_psm_records(
                obj,
                context,
                config.structures,
                removed_index=removed_index,
                confidence=confidence,
                include_unassigned_psms=config.include_unassigned_psms,
            )
            feat_buf.extend(cf_feats)
            feature_count += len(cf_feats)
            psm_buf.extend(cf_psms)
            collect_psm_scan_formats(obj.getPeptideIdentifications(), scan_formats, enabled=pw is not None)
            if maps is not None:
                accumulate_cf_maps(obj, map_run, context.resolve_run, maps)
        else:  # unassigned peptide identification
            if track:
                seen.unassigned = True
            if pw is not None and config.include_unassigned_psms:
                # Unassigned PSMs map to no feature -> feature_id stays null.
                psm_buf.extend(
                    psm_records_for_pid(
                        obj, context.resolve_run, seen, enzyme=context.enzyme, confidence=confidence, ms1=context.ms1
                    )
                )
                collect_psm_scan_formats([obj], scan_formats)
            if maps is not None:
                # Protein inference always sees every identification, whether or
                # not the PSM rows are emitted: dropping evidence would change the
                # protein groups, which is not what this option is for.
                accumulate_unassigned_maps(obj, context.resolve_run, maps)
        feat_buf = _flush_records(fw, feat_buf, batch)
        psm_buf = _flush_records(pw, psm_buf, batch)
    _flush_records(fw, feat_buf)
    _flush_records(pw, psm_buf)
    return feature_count


def _open_stream_writers(stack, config) -> tuple:
    """Open requested pruning writers in their original context-manager order."""
    fw = pw = None
    output = config.output
    if "feature" in config.structures:
        fw = stack.enter_context(
            _prunable_writer(FeatureWriter)(
                str(output.folder / f"{output.prefix}.feature.parquet"),
                creator=output.creator,
                identity_composite=_FEATURE_IDENTITY_COMPOSITE,
                compression=output.compression,
                batch_size=_STREAM_BATCH,
            )
        )
    if "psm" in config.structures:
        pw = stack.enter_context(
            _prunable_writer(PsmWriter)(
                str(output.folder / f"{output.prefix}.psm.parquet"),
                creator=output.creator,
                identity_composite=_PSM_IDENTITY_COMPOSITE,
                compression=output.compression,
                batch_size=_STREAM_BATCH,
            )
        )
    return fw, pw


def _finalize_stream_writers(writers, seen, dedup, output) -> dict:
    """Attach final row/link patches and register nonempty staged views."""
    fw, pw = writers
    written = {}
    if fw is not None:
        fw.drop_rows = frozenset(dedup.superseded)
        fw.record_patches = dedup.row_patches
        if dedup.rows > len(dedup.superseded):
            written["feature"] = output.folder / f"{output.prefix}.feature.parquet"
    if pw is not None:
        pw.drop_rows = frozenset(seen.dropped_rows)
        pw.feature_id_patch = seen.feature_id_patch
        if dedup.redirect:
            pw.remap_feature_id = dedup.resolve
        if seen.rows > len(seen.dropped_rows):
            written["psm"] = output.folder / f"{output.prefix}.psm.parquet"
    return written


def _convert_streaming(consensusxml_paths, config) -> dict:
    """Run the identification pre-pass and batched output without retaining a whole XML map."""
    from qpx.converters.openms_consensus.streaming import StreamingConsensusMap

    pg = ProteinGroupAccumulator() if "pg" in config.structures else None
    seen = _PsmKeys()
    dedup = _FeatureDeduplicator(pg.pep_intensity if pg is not None else None)
    run_owner: dict[str, str] = {}
    scan_formats: set[str | None] = set()
    with ExitStack() as stack:
        writers = _open_stream_writers(stack, config)
        for path in consensusxml_paths:
            cm = StreamingConsensusMap(path)
            _claim_runs(run_owner, path, cm)
            _warn_channel_mismatch(cm, config.sdrf_path)
            context = _RecordContext.from_map(cm, config.structures, seen, dedup, config.sdrf_path)
            if pg is not None:
                pg.add_source(cm, context.map_info, path)
            seen.start_file()
            _stream_feature_psm(cm, writers, context, pg, config, scan_formats)
        declare_psm_scan_format(writers[1], scan_formats)
        dedup.log()
        written = _finalize_stream_writers(writers, seen, dedup, config.output)
    _write_pg(pg, config, written)
    return written


def _memory_map_records(cm, path, seen, dedup, pg, config, scan_formats) -> tuple[list, list]:
    """Collect assigned records before unassigned PSMs from one loaded map."""
    context = _RecordContext.from_map(cm, config.structures, seen, dedup, config.sdrf_path)
    seen.start_file()
    seen.before_features = False
    feat_recs, psm_recs = [], []
    removed_index, confidence = identification_context(cm)
    for cf in cm:
        cf_feats, cf_psms = _cf_feature_psm_records(
            cf,
            context,
            config.structures,
            removed_index=removed_index,
            confidence=confidence,
            include_unassigned_psms=config.include_unassigned_psms,
        )
        feat_recs.extend(cf_feats)
        psm_recs.extend(cf_psms)
        collect_psm_scan_formats(cf.getPeptideIdentifications(), scan_formats, enabled="psm" in config.structures)
    seen.unassigned = True
    if "psm" in config.structures and config.include_unassigned_psms:
        for pid in cm.getUnassignedPeptideIdentifications():
            psm_recs.extend(
                psm_records_for_pid(pid, context.resolve_run, seen, enzyme=context.enzyme, confidence=confidence, ms1=context.ms1)
            )
            collect_psm_scan_formats([pid], scan_formats)
    if pg is not None:
        pg.add_source(cm, context.map_info, path)
        accumulate_consensus_map(cm, context.map_info, context.resolve_run, pg.maps)
    return feat_recs, psm_recs


def _write_memory_records(feat_recs, psm_recs, scan_formats, config) -> dict:
    """Write collected feature/PSM records with their canonical identities."""
    output = config.output
    written = {}
    if feat_recs:
        written["feature"] = _write_view(
            FeatureWriter,
            output.folder / f"{output.prefix}.feature.parquet",
            feat_recs,
            creator=output.creator,
            compression=output.compression,
            identity_composite=_FEATURE_IDENTITY_COMPOSITE,
        )
    if "psm" in config.structures and psm_recs:
        written["psm"] = _write_view(
            PsmWriter,
            output.folder / f"{output.prefix}.psm.parquet",
            psm_recs,
            creator=output.creator,
            compression=output.compression,
            identity_composite=_PSM_IDENTITY_COMPOSITE,
            scan_format=uniform_scan_format(scan_formats),
        )
    return written


def _convert_in_memory(consensusxml_paths, config) -> dict:
    """Load each pyopenms map once and collect its records in memory."""
    pg = ProteinGroupAccumulator() if "pg" in config.structures else None
    seen = _PsmKeys()
    dedup = _FeatureDeduplicator(pg.pep_intensity if pg is not None else None)
    run_owner: dict[str, str] = {}
    feat_recs: list[dict] = []
    psm_recs: list[dict] = []
    scan_formats: set[str | None] = set()
    for path in consensusxml_paths:
        cm = load_consensus_map(path)
        _claim_runs(run_owner, path, cm)
        _warn_channel_mismatch(cm, config.sdrf_path)
        features, psms = _memory_map_records(cm, path, seen, dedup, pg, config, scan_formats)
        feat_recs.extend(features)
        psm_recs.extend(psms)
    dedup.log()
    feat_recs = _finalize_in_memory_records(feat_recs, psm_recs, seen, dedup)
    written = _write_memory_records(feat_recs, psm_recs, scan_formats, config)
    _write_pg(pg, config, written)
    return written
