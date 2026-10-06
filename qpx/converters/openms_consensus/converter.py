"""Orchestrate consensusXML (+ SDRF) -> a QPX dataset (feature/psm/pg [+ run/sample]).

Interim quantms path while OpenMS ``-out_qpx`` is pre-1.1. pg carries an interim
unnormalized unique-peptide-sum intensity (until ``-out_qpx`` provides the real
quant). run/sample come from the SDRF (reusing :class:`SdrfConverter`) when an
SDRF is provided; when it is, the consensusXML channels are also checked against
the SDRF ``comment[label]`` and mismatches are logged as warnings.
"""

from __future__ import annotations

import logging
import math
import os
from collections.abc import Sequence
from pathlib import Path
from typing import Optional, Union

import pyarrow as pa
import pyarrow.parquet as pq

from qpx._version import __version__
from qpx.converters.openms_consensus.feature_adapter import (
    check_channels_vs_sdrf,
    column_runs,
    load_consensus_map,
)
from qpx.converters.openms_consensus.feature_dedup import (
    FEATURE_IDENTITY_COMPOSITE as _FEATURE_IDENTITY_COMPOSITE,
)
from qpx.converters.openms_consensus.feature_dedup import (
    FeatureDeduplicator as _FeatureDeduplicator,
)
from qpx.converters.openms_consensus.feature_dedup import (
    feature_ids as _feature_ids,
)
from qpx.converters.openms_consensus.pg_adapter import ProteinGroupAccumulator, protein_group_maps
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.core.constants import DATASET, FEATURE, ONTOLOGY, PG, PROVENANCE, PSM, RUN, SAMPLE
from qpx.core.scan import scan_format_from_native_id
from qpx.writers.base import parquet_write_options
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter

_log = logging.getLogger(__name__)

_STRUCTURE_ALL = ("feature", "psm", "pg", "run", "sample")

# Schema identity for PSM rows.
_PSM_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name", "scan")


_NO_LINK = (0,)


def _link_rank(peptidoform, charge, rt, mz, feature: Optional[dict], feature_id: Optional[int], quality: float) -> tuple:
    """Order-independent preference of one candidate feature for a PSM (larger is better).

    A feature in the PSM's own run beats no feature; then a feature quantifying the
    PSM's own peptidoform/charge; then the feature whose RT and m/z lie closest to
    the precursor; then consensus quality; the feature_id breaks any remaining tie.
    """
    if feature is None or feature_id is None:
        return _NO_LINK
    same = feature.get("peptidoform") == peptidoform and feature.get("charge") == charge
    drt = abs(feature["rt"] - rt) if rt is not None and feature.get("rt") is not None else math.inf
    dmz = abs((feature.get("observed_mz") or 0.0) - (mz or 0.0))
    return (1, same, -drt, -dmz, quality, -feature_id)


def _link_psm_feature(
    feature_records: list[dict],
    psm_records: list[dict],
    feature_ids: Optional[list[int]] = None,
    *,
    links: Optional["_PsmKeys"] = None,
    duplicates: Sequence[tuple] = (),
    quality: float = 0.0,
) -> None:
    """Stamp ``psm.feature_id`` in place for one consensus feature.

    ``psm.feature_id`` is the authoritative producer assignment — which consensus
    feature a PSM belongs to. The id is derived with the SAME composite the
    FeatureWriter uses, so this converter-populated foreign key matches the
    writer-derived ``feature_id`` byte-for-byte (``derive_id`` is deterministic).
    A consensus feature yields one feature record per run; each PSM links to the
    feature record of its own run. PSMs whose run has no feature record keep
    ``feature_id`` null (resolves #182).

    Merged identifications put copies of one spectrum into several consensus
    features. ``links`` remembers each emitted PSM's best candidate (see
    :func:`_link_rank`); ``duplicates`` are this feature's copies of PSMs emitted
    earlier, and a better one re-points the earlier row, so the link does not
    depend on which copy is read first.

    The inverse ``feature.psm_ids`` is NOT materialized here: it is the pure
    inverse of ``psm.feature_id`` (group PSMs by their ``feature_id``) and is
    computed on read via :meth:`qpx.dataset.Dataset.link_feature_psm`.
    """
    ids = feature_ids if feature_ids is not None else _feature_ids(feature_records)
    by_run: dict[str, tuple[dict, int]] = {rec["run_file_name"]: (rec, fid) for rec, fid in zip(feature_records, ids)}
    for prec in psm_records:
        rec, feat_id = by_run.get(prec.get("run_file_name"), (None, None))
        if feat_id is not None:
            prec["feature_id"] = feat_id
        if links is not None:
            key = (prec["peptidoform"], prec["charge"], prec["run_file_name"], tuple(prec["scan"]))
            rank = _link_rank(prec["peptidoform"], prec["charge"], prec.get("rt"), prec.get("observed_mz"), rec, feat_id, quality)
            links.offer(key, rank, feat_id)
    if links is not None:
        for key, rt, mz in duplicates:
            rec, feat_id = by_run.get(key[2], (None, None))
            links.offer(key, _link_rank(key[0], key[1], rt, mz, rec, feat_id, quality), feat_id)


class _PsmKeys:  # pylint: disable=too-many-instance-attributes
    """PSM identity keys already emitted (the dedup ``seen`` set shared by both paths).

    OpenMS writes unassigned identifications before the consensus features, while the
    in-memory path visits assigned ones first. A key first emitted by an unassigned
    identification ahead of the features stays claimable: when an assigned copy of the
    same spectrum is looked up, the membership test releases the key and records the
    unassigned row in ``dropped_rows``, so both paths keep the same (linked) PSM.

    Keys emitted from a consensus feature also keep their output row and link rank:
    :meth:`offer` re-points that row (``feature_id_patch``) when a later copy of the
    spectrum ranks better, so the final link does not depend on reading order.
    """

    def __init__(self):
        self._keys: set = set()
        self._claimable: dict = {}
        self._links: dict = {}
        self.rows = 0
        self.dropped_rows: set[int] = set()
        self.feature_id_patch: dict[int, Optional[int]] = {}
        self.unassigned = False
        self.before_features = True

    def __contains__(self, key) -> bool:
        if key not in self._keys:
            return False
        if not self.unassigned and key in self._claimable:
            self.dropped_rows.add(self._claimable.pop(key))
            return False
        return True

    def add(self, key) -> None:
        """Register an emitted key at the next output row."""
        if self.unassigned and self.before_features:
            self._claimable[key] = self.rows
        elif not self.unassigned:
            self._links[key] = (None, self.rows)
        self._keys.add(key)
        self.rows += 1

    def offer(self, key, rank: tuple, feature_id: Optional[int]) -> None:
        """Record a candidate link for an emitted key; keep it when it outranks the current one."""
        entry = self._links.get(key)
        if entry is None:
            return
        current, row = entry
        if current is None:
            self._links[key] = (rank, row)
        elif rank > current:
            self._links[key] = (rank, row)
            self.feature_id_patch[row] = feature_id

    def start_file(self) -> None:
        """Reset the per-file state before the next input."""
        # Runs never span inputs, so no later file can offer a link for these keys.
        self._claimable.clear()
        self._links.clear()
        self.unassigned = False
        self.before_features = True


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


def _as_paths(consensusxml_path) -> list[str]:
    """Normalise one path or a sequence of paths to a non-empty list of strings."""
    if isinstance(consensusxml_path, (str, os.PathLike)):
        return [str(consensusxml_path)]
    paths = [str(p) for p in consensusxml_path]
    if not paths:
        raise ValueError("At least one consensusXML input is required")
    return paths


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


# Above this consensusXML size, auto-select the streaming reader: the pyopenms
# in-memory load needs ~0.8x the file in RAM, so ~4 GB (≈3 GB RAM) is a safe point
# to switch to the low-memory path on a normal node.
_STREAM_THRESHOLD_BYTES = 4 * 1024**3


def _should_stream(path: str) -> bool:
    try:
        return Path(path).stat().st_size > _STREAM_THRESHOLD_BYTES
    except OSError:
        return False


def _cf_feature_psm_records(
    cf,
    map_info,
    group_map,
    resolve_run,
    seen,
    *,
    want_feature,
    want_psm,
    enzyme=None,
    group_meta=None,
    dedup=None,
    removed_index=None,
    confidence=None,
    include_unassigned_psms=True,
    want_pg=False,
):
    """Feature + PSM records for one consensus feature, cross-linked when both views
    are emitted. Shared by the streaming and in-memory paths so their output matches.
    ``dedup`` (a :class:`_FeatureDeduplicator`) filters duplicate feature rows.
    ``include_unassigned_psms`` controls whether recovered source-unassigned IDs
    are emitted as PSM rows; they remain available for Feature metadata recovery.
    """
    from qpx.converters.openms_consensus.feature_adapter import feature_records_for_cf, removed_same_peptide_ids
    from qpx.converters.openms_consensus.psm_adapter import _cf_element_runs, psm_records_for_pid

    removed = removed_same_peptide_ids(cf, removed_index)
    cf_feats = (
        feature_records_for_cf(
            cf,
            map_info,
            group_map,
            enzyme=enzyme,
            group_meta=group_meta,
            resolve_run=resolve_run,
            removed_pids=removed,
            confidence=confidence,
        )
        if want_feature or want_pg
        else []
    )
    cf_psms: list[dict] = []
    link = want_feature and want_psm
    duplicates: Optional[list] = [] if link and isinstance(seen, _PsmKeys) else None
    if want_psm:
        # Multi-run isobaric PIDs carry a local id_merge_index; the feature's
        # element runs disambiguate which run they belong to.
        cf_runs = _cf_element_runs(cf, map_info)
        pids = list(cf.getPeptideIdentifications())
        if include_unassigned_psms:
            pids.extend(removed)
        for pid in pids:
            cf_psms.extend(
                psm_records_for_pid(
                    pid, resolve_run, seen, enzyme=enzyme, cf_runs=cf_runs, duplicates=duplicates, confidence=confidence
                )
            )
    feature_ids = _feature_ids(cf_feats) if cf_feats and (want_psm or dedup is not None) else None
    # Stamp psm.feature_id only when both views are emitted (the FK references a
    # feature row written in this dataset). feature.psm_ids is the computed inverse.
    if link:
        _link_psm_feature(
            cf_feats,
            cf_psms,
            feature_ids,
            links=seen if duplicates is not None else None,
            duplicates=duplicates or (),
            quality=float(cf.getQuality()),
        )
    if dedup is not None and cf_feats:
        cf_feats = dedup.keep(cf_feats, feature_ids, float(cf.getQuality()))
    return cf_feats if want_feature else [], cf_psms


def _collect_psm_scan_formats(identifications, formats: set[str | None], *, enabled=True) -> None:
    """Accumulate explicit PID formats; unknown references prevent a file declaration."""
    if not enabled:
        return
    for pid in identifications:
        if not pid.getHits():
            continue
        reference = pid.getSpectrumReference() if hasattr(pid, "getSpectrumReference") else ""
        if not reference and pid.metaValueExists("spectrum_reference"):
            reference = pid.getMetaValue("spectrum_reference")
        if reference:
            formats.add(scan_format_from_native_id(str(reference)))


def _uniform_scan_format(formats: set[str | None]) -> str | None:
    """Return the declaration only after the full input confirms one known format."""
    return next(iter(formats)) if len(formats) == 1 else None


def _declare_psm_scan_format(writer, formats: set[str | None]) -> None:
    """Finalize a streamed PSM writer's confirmed declaration before it closes."""
    scan_format = _uniform_scan_format(formats)
    if writer is not None and scan_format is not None:
        writer.set_scan_format(scan_format)


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


def _stream_feature_psm(
    cm,
    fw,
    pw,
    *,
    map_info,
    group_map,
    resolve_run,
    maps,
    map_run,
    seen,
    batch,
    scan_formats: set[str | None],
    include_unassigned_psms=True,
    enzyme=None,
    group_meta=None,
    dedup=None,
) -> int:
    """One ordered element/unassigned pass: write feature/psm in batches and
    accumulate the pg maps and PSM scan formats in place. Return the number of
    emitted feature records."""
    from qpx.converters.openms_consensus.feature_adapter import (
        identification_context,
    )
    from qpx.converters.openms_consensus.pg_adapter import (
        accumulate_cf_maps,
        accumulate_unassigned_maps,
    )
    from qpx.converters.openms_consensus.psm_adapter import psm_records_for_pid

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
                map_info,
                group_map,
                resolve_run,
                seen,
                want_feature=fw is not None,
                want_pg=maps is not None,
                want_psm=pw is not None,
                enzyme=enzyme,
                group_meta=group_meta,
                dedup=dedup,
                removed_index=removed_index,
                confidence=confidence,
                include_unassigned_psms=include_unassigned_psms,
            )
            feat_buf.extend(cf_feats)
            feature_count += len(cf_feats)
            psm_buf.extend(cf_psms)
            _collect_psm_scan_formats(obj.getPeptideIdentifications(), scan_formats, enabled=pw is not None)
            if maps is not None:
                accumulate_cf_maps(obj, map_run, maps)
        else:  # unassigned peptide identification
            if track:
                seen.unassigned = True
            if pw is not None and include_unassigned_psms:
                # Unassigned PSMs map to no feature -> feature_id stays null.
                psm_buf.extend(psm_records_for_pid(obj, resolve_run, seen, enzyme=enzyme, confidence=confidence))
                _collect_psm_scan_formats([obj], scan_formats)
            if maps is not None:
                # Protein inference always sees every identification, whether or
                # not the PSM rows are emitted: dropping evidence would change the
                # protein groups, which is not what this option is for.
                accumulate_unassigned_maps(obj, resolve_run, maps)
        feat_buf = _flush_records(fw, feat_buf, batch)
        psm_buf = _flush_records(pw, psm_buf, batch)
    _flush_records(fw, feat_buf)
    _flush_records(pw, psm_buf)
    return feature_count


def _convert_streaming(
    consensusxml_paths,
    out,
    output_prefix,
    structures,
    sdrf_path,
    creator,
    pg_top,
    compression="zstd",
    include_unassigned_psms=True,
) -> dict:
    """Two-pass, low-memory feature/psm/pg from streamed consensusXML input(s).

    A pre-pass indexes recovered IDs and peptide confidence together. An output pass feeds
    the same per-element builders the in-memory adapters use (so output is identical),
    writing feature/psm in batches into one set of writers and accumulating the pg
    maps; pg records are built at the end. No input map is held whole, and each file
    has two complete XML traversals. Rows that a later record supersedes (a better duplicate feature,
    or an assigned copy of an unassigned PSM) are pruned from the staged files.
    """
    from contextlib import ExitStack

    from qpx.converters.openms_consensus.feature_adapter import feature_map_info, resolve_enzyme
    from qpx.converters.openms_consensus.psm_adapter import _run_resolver
    from qpx.converters.openms_consensus.streaming import StreamingConsensusMap

    want_feature, want_psm, want_pg = ("feature" in structures, "psm" in structures, "pg" in structures)
    pg = ProteinGroupAccumulator() if want_pg else None
    seen = _PsmKeys()
    dedup = _FeatureDeduplicator(pg.pep_intensity if pg is not None else None)
    run_owner: dict[str, str] = {}
    written: dict[str, Path] = {}
    # One PSM writer spans every input, so its format is declared once for all of them.
    scan_formats: set[str | None] = set()

    with ExitStack() as stack:
        fw = pw = None
        if want_feature:
            fw = stack.enter_context(
                _prunable_writer(FeatureWriter)(
                    str(out / f"{output_prefix}.feature.parquet"),
                    creator=creator,
                    identity_composite=_FEATURE_IDENTITY_COMPOSITE,
                    compression=compression,
                )
            )
        if want_psm:
            pw = stack.enter_context(
                _prunable_writer(PsmWriter)(
                    str(out / f"{output_prefix}.psm.parquet"),
                    creator=creator,
                    identity_composite=_PSM_IDENTITY_COMPOSITE,
                    compression=compression,
                )
            )
        for path in consensusxml_paths:
            cm = StreamingConsensusMap(path)
            _claim_runs(run_owner, path, cm)
            _warn_channel_mismatch(cm, sdrf_path)
            map_info = feature_map_info(cm)
            group_map, group_meta = protein_group_maps(cm) if want_feature or want_pg else (None, None)
            if pg is not None:
                pg.add_source(cm, map_info)
            seen.start_file()
            _stream_feature_psm(
                cm,
                fw,
                pw,
                map_info=map_info,
                group_map=group_map,
                resolve_run=_run_resolver(cm),
                maps=pg.maps if pg is not None else None,
                map_run=column_runs(cm),
                seen=seen,
                batch=100_000,
                scan_formats=scan_formats,
                include_unassigned_psms=include_unassigned_psms,
                enzyme=resolve_enzyme(cm, sdrf_path),
                group_meta=group_meta,
                dedup=dedup,
            )
        _declare_psm_scan_format(pw, scan_formats)
        dedup.log()
        if fw is not None:
            fw.drop_rows = frozenset(dedup.superseded)
            fw.record_patches = dedup.row_patches
            if dedup.rows > len(dedup.superseded):
                written["feature"] = out / f"{output_prefix}.feature.parquet"
        if pw is not None:
            pw.drop_rows = frozenset(seen.dropped_rows)
            pw.feature_id_patch = seen.feature_id_patch
            if dedup.redirect:
                pw.remap_feature_id = dedup.resolve
            if seen.rows > len(seen.dropped_rows):
                written["psm"] = out / f"{output_prefix}.psm.parquet"

    if pg is not None:
        records = pg.build(sdrf_path, pg_top)
        if records:
            written["pg"] = _write_view(
                PgWriter, out / f"{output_prefix}.pg.parquet", records, creator=creator, compression=compression
            )
    return written


def _validate_structures(structures: tuple[str, ...], sdrf_path: Optional[str]) -> None:
    """Reject unknown structure names and run/sample requested without an SDRF."""
    unknown = [s for s in structures if s not in _STRUCTURE_ALL]
    if unknown:
        raise ValueError(f"Unknown structure(s) {unknown}; valid values are {list(_STRUCTURE_ALL)}")
    metadata = set(structures).intersection({"run", "sample"})
    if metadata and not sdrf_path:
        raise ValueError(f"An SDRF is required to write {sorted(metadata)}")


def _collect_score_names(table_path: Path) -> set[str]:
    """Read ``additional_scores`` score names from a written parquet file."""
    table = pq.read_table(str(table_path), columns=["additional_scores"])
    names: set[str] = set()
    for row in table.column("additional_scores").to_pylist():
        if row:
            for score in row:
                if score and score.get("score_name"):
                    names.add(score["score_name"])
    return names


def _consensus_is_isobaric(consensusxml_path: str) -> bool:
    """Return whether a consensusXML's map labels indicate an isobaric (TMT/iTRAQ) run.

    Reads only the small column-map header via the streaming reader, so it never
    loads the whole ConsensusMap into memory — the full pyopenms load would defeat
    the streaming path and can exhaust RAM on multi-GB files.
    """
    from qpx.converters.openms_consensus.feature_adapter import consensus_channels
    from qpx.converters.openms_consensus.streaming import StreamingConsensusMap

    try:
        return bool(consensus_channels(StreamingConsensusMap(consensusxml_path)))
    except (OSError, ValueError, KeyError, AttributeError, RuntimeError) as exc:
        # Best-effort provenance hint: any read/parse failure -> assume label-free.
        _log.warning("could not determine isobaric labeling from %s: %s", consensusxml_path, exc)
        return False


def _write_sdrf_metadata(
    output_folder: Path,
    output_prefix: str,
    sdrf_path: str,
    requested: set[str],
    compression: str = "zstd",
) -> tuple[dict[str, Path], list[dict]]:
    """Write the requested SDRF-backed run/sample structures.

    Returns ``(paths, run_ontology_entries)`` — the latter collected from the
    same ``SdrfConverter`` instance (its parse state populates the run ontology),
    so callers need not re-parse the SDRF.
    """
    metadata = requested.intersection({"run", "sample"})
    if not metadata:
        return {}, []

    from qpx.converters.sdrf import SdrfConverter

    paths = {
        "sample": output_folder / f"{output_prefix}.sample.parquet",
        "run": output_folder / f"{output_prefix}.run.parquet",
    }
    with SdrfConverter(compression=compression) as sdrf_converter:
        sdrf_converter.convert(
            sdrf_path=sdrf_path,
            sample_output=str(paths["sample"]) if "sample" in metadata else None,
            run_output=str(paths["run"]) if "run" in metadata else None,
        )
        run_ontology = list(sdrf_converter.run_ontology_entries())
    return {name: paths[name] for name in metadata}, run_ontology


def _remove_orphaned_metadata(output_folder: Path, output_prefix: str) -> None:
    """Remove metadata after an empty rerun only if no same-prefix data remains."""
    if any((output_folder / f"{output_prefix}.{view}.parquet").is_file() for view in _STRUCTURE_ALL):
        return
    for view in (DATASET, ONTOLOGY, PROVENANCE):
        (output_folder / f"{output_prefix}.{view}.parquet").unlink(missing_ok=True)


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


class OpenMSConsensusConverter(BaseOrchestrator):  # pylint: disable=too-few-public-methods
    """consensusXML + SDRF -> QPX views.

    A single-entry orchestrator (``convert``) — the interim counterpart to the
    other converter classes, kept as a class for call-site symmetry with them.
    """

    def convert(
        self,
        consensusxml_path: Union[str, os.PathLike, Sequence[Union[str, os.PathLike]]],
        output_folder: str,
        output_prefix: str = "openms",
        sdrf_path: Optional[str] = None,
        structures: tuple[str, ...] = _STRUCTURE_ALL,
        creator: str = "openms-consensus",
        pg_top: int = 0,
        streaming: Optional[bool] = None,
        project_accession: Optional[str] = None,
        compression: str = "zstd",
        include_unassigned_psms: bool = True,
    ) -> dict[str, Path]:
        """Write the requested QPX views and return ``{structure: parquet path}``.

        ``consensusxml_path`` is one consensusXML or a sequence of them. Several
        inputs (e.g. one per sample group) become ONE dataset: each file keeps its
        own column -> run mapping and identification metadata, rows go to the same
        views, pg spans all inputs and provenance lists every input. Two inputs
        must not contain the same run (``ValueError``).

        ``include_unassigned_psms`` (default ``True``) controls whether
        unassigned PeptideIdentifications from the source reach ``psm.parquet``.
        They are retained to preserve identification evidence; pass ``False``
        to exclude them.

        Requested feature/PSM/PG views with no exportable records are skipped
        with a warning and excluded from output metadata. After successful core
        conversion, remove any existing file for each skipped view and this
        output prefix so Dataset cannot discover stale records. Other requested
        views are still exported; an entirely empty export returns an empty dict.
        When no same-prefix core or run/sample views remain, their orphaned
        ontology/provenance/dataset metadata is also removed.

        ``feature_id`` records a link in the exported dataset, not quantification
        status. It is only populated when both feature and PSM views are emitted;
        PSM-only output leaves it null even for assigned identifications.
        ``PSM.with_feature()`` and ``PSM.without_feature()`` filter these links.

        The feature view always omits consensus features without peptide hits.
        Protein inference always uses every identification, regardless of this
        option.

        ``structures`` selects which of feature/psm/pg/run/sample to emit. pg
        carries an interim unnormalized unique-peptide-sum intensity; ``pg_top``
        bounds the peptides used (0 = all; 3 mirrors ProteomicsLFQ/IsobaricWorkflow).

        ``streaming`` picks the consensusXML reader: ``None`` (default) auto-selects
        — the low-memory streaming reader when any input is above ``_STREAM_THRESHOLD_BYTES``
        (which pyopenms would otherwise load whole into ~0.8x-file RAM), else the
        faster in-memory pyopenms load. ``True``/``False`` forces the choice.

        ``project_accession`` (e.g. ``PXD001819``) is stamped into dataset.parquet.

        ``compression`` is the Parquet codec for every written view (default zstd).
        """
        self._compression = compression
        _validate_structures(structures, sdrf_path)
        consensusxml_paths = _as_paths(consensusxml_path)
        requested = set(structures)

        out = Path(output_folder)
        out.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}

        if {"feature", "psm", "pg"}.intersection(structures):
            use_stream = streaming if streaming is not None else any(_should_stream(p) for p in consensusxml_paths)
            if use_stream:
                # A shared identification pre-pass precedes batched output of all
                # views; the complete XML map is never retained in memory.
                written.update(
                    _convert_streaming(
                        consensusxml_paths,
                        out,
                        output_prefix,
                        structures,
                        sdrf_path,
                        creator,
                        pg_top,
                        compression,
                        include_unassigned_psms,
                    )
                )
            else:
                # In-memory path: pyopenms loads the map once (fast for smaller files);
                # the adapters iterate it cheaply. Output is identical either way.
                written.update(
                    self._convert_in_memory(
                        consensusxml_paths,
                        out,
                        output_prefix,
                        structures,
                        sdrf_path,
                        creator,
                        pg_top,
                        compression,
                        include_unassigned_psms,
                    )
                )

        for view in (FEATURE, PSM, PG):
            if view in requested and view not in written:
                (out / f"{output_prefix}.{view}.parquet").unlink(missing_ok=True)
                _log.warning("No exportable %s records; skipping %s Parquet output.", view.upper(), view.upper())

        sdrf_paths, run_ontology = _write_sdrf_metadata(out, output_prefix, sdrf_path, requested, compression)
        written.update(sdrf_paths)

        # Metadata tables (ontology / provenance / dataset) — written when the
        # core + run/sample structures they describe are present. Best-effort:
        # ontology entries only materialise when a PSI-MS term resolves.
        if written:
            self._write_metadata(out, output_prefix, written, consensusxml_paths, run_ontology, requested, project_accession)
        else:
            _remove_orphaned_metadata(out, output_prefix)

        return written

    def _write_metadata(
        self,
        out: Path,
        output_prefix: str,
        written: dict[str, Path],
        consensusxml_paths: list[str],
        run_ontology: list[dict],
        requested: set[str],
        project_accession: Optional[str] = None,
    ) -> None:
        """Write ontology/provenance/dataset metadata tables (best-effort).

        Mirrors the ``openms`` enrichment path: ontology collects the run-level
        terms from the SDRF plus the score names discovered in the written core
        tables; provenance records the consensusXML -> QPX conversion steps; the
        dataset table carries the project-level metadata.
        """
        ontology_entries = list(run_ontology)
        for view in (PSM, FEATURE, PG):
            path = written.get(view)
            if not path:
                continue
            try:
                names = _collect_score_names(path)
            except (KeyError, pa.ArrowInvalid):
                continue
            if names:
                from qpx.core.scores import score_ontology_entries

                ontology_entries.extend(score_ontology_entries(names, view=view))  # noqa: PERF401
        self._write_ontology(out, output_prefix, ontology_entries)

        structures = sorted(requested.intersection(written, {PSM, FEATURE, PG, RUN, SAMPLE}))
        provenance_records = self._build_provenance(structures, consensusxml_paths)
        self._write_provenance(out, output_prefix, provenance_records)
        if provenance_records:
            self._write_dataset(
                out,
                output_prefix,
                project_accession,
                software_name=provenance_records[0]["tool_name"],
                software_version=None,
                provenance_records=provenance_records,
            )

    @staticmethod
    def _build_provenance(structures: list[str], consensusxml_paths: list[str]) -> list[dict]:
        """Provenance records: consensusXML parsing + QPX conversion (listing every input)."""
        is_isobaric = any(_consensus_is_isobaric(path) for path in consensusxml_paths)
        step_name = "isobaric_quantification" if is_isobaric else "label_free_quantification"
        tool_name = "OpenMS/IsobaricWorkflow" if is_isobaric else "OpenMS/ProteomicsLFQ"
        return [
            {
                "step_order": 1,
                "step_category": "quantification",
                "step_name": step_name,
                "tool_name": tool_name,
                "tool_version": None,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": [v for v in (FEATURE, PSM, PG) if v in structures],
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "openms_consensus_qpx_conversion",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": [{"key": "consensusxml", "value": Path(path).name} for path in consensusxml_paths],
                "config": None,
                "output_views": [v for v in (SAMPLE, RUN) if v in structures] + [ONTOLOGY],
            },
        ]

    @staticmethod
    def _convert_in_memory(
        consensusxml_paths,
        out,
        output_prefix,
        structures,
        sdrf_path,
        creator,
        pg_top,
        compression="zstd",
        include_unassigned_psms=True,
    ) -> dict:
        """feature/psm/pg via in-memory pyopenms maps (each input loaded once, iterated cheaply)."""
        from qpx.converters.openms_consensus.feature_adapter import (
            feature_map_info,
            identification_context,
            resolve_enzyme,
        )
        from qpx.converters.openms_consensus.pg_adapter import accumulate_consensus_map
        from qpx.converters.openms_consensus.psm_adapter import _run_resolver, psm_records_for_pid

        written: dict[str, Path] = {}
        want_feature, want_psm, want_pg = "feature" in structures, "psm" in structures, "pg" in structures
        pg = ProteinGroupAccumulator() if want_pg else None
        seen = _PsmKeys()
        dedup = _FeatureDeduplicator(pg.pep_intensity if pg is not None else None)
        run_owner: dict[str, str] = {}
        feat_recs: list[dict] = []
        psm_recs: list[dict] = []
        scan_formats: set[str | None] = set()
        for path in consensusxml_paths:
            cm = load_consensus_map(path)
            _claim_runs(run_owner, path, cm)
            _warn_channel_mismatch(cm, sdrf_path)
            map_info = feature_map_info(cm)
            resolve_run = _run_resolver(cm)
            seen.start_file()
            seen.before_features = False
            if want_feature or want_psm or want_pg:
                # Build feature and psm per consensus feature (mirroring the streaming
                # path's element loop) so their cross-refs can be linked identically:
                # assigned PSMs first (cf order), then the unassigned PSMs.
                # Share the full protein-group membership so feature.anchor_protein and
                # feature.pg_accessions match pg (unambiguous even for shared leaders).
                group_map, group_meta = protein_group_maps(cm) if want_feature or want_pg else (None, None)
                enzyme = resolve_enzyme(cm, sdrf_path)
                removed_index, confidence = identification_context(cm)
                for cf in cm:
                    cf_feats, cf_psms = _cf_feature_psm_records(
                        cf,
                        map_info,
                        group_map,
                        resolve_run,
                        seen,
                        want_feature=want_feature,
                        want_pg=want_pg,
                        want_psm=want_psm,
                        enzyme=enzyme,
                        group_meta=group_meta,
                        dedup=dedup,
                        removed_index=removed_index,
                        confidence=confidence,
                        include_unassigned_psms=include_unassigned_psms,
                    )
                    feat_recs.extend(cf_feats)
                    psm_recs.extend(cf_psms)
                    _collect_psm_scan_formats(cf.getPeptideIdentifications(), scan_formats, enabled=want_psm)
                seen.unassigned = True
                if want_psm and include_unassigned_psms:
                    for pid in cm.getUnassignedPeptideIdentifications():
                        psm_recs.extend(psm_records_for_pid(pid, resolve_run, seen, enzyme=enzyme, confidence=confidence))
                        _collect_psm_scan_formats([pid], scan_formats)
            if pg is not None:
                pg.add_source(cm, map_info)
                accumulate_consensus_map(cm, map_info, resolve_run, pg.maps)
        dedup.log()
        feat_recs = _finalize_in_memory_records(feat_recs, psm_recs, seen, dedup)
        if feat_recs:
            written["feature"] = _write_view(
                FeatureWriter,
                out / f"{output_prefix}.feature.parquet",
                feat_recs,
                creator=creator,
                compression=compression,
                identity_composite=_FEATURE_IDENTITY_COMPOSITE,
            )
        if want_psm and psm_recs:
            written["psm"] = _write_view(
                PsmWriter,
                out / f"{output_prefix}.psm.parquet",
                psm_recs,
                creator=creator,
                compression=compression,
                identity_composite=_PSM_IDENTITY_COMPOSITE,
                scan_format=_uniform_scan_format(scan_formats),
            )
        if pg is not None:
            recs = pg.build(sdrf_path, pg_top)
            if recs:
                written["pg"] = _write_view(
                    PgWriter, out / f"{output_prefix}.pg.parquet", recs, creator=creator, compression=compression
                )
        return written
