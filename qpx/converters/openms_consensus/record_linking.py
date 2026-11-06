"""Consensus-feature records and deterministic PSM/feature links."""

from __future__ import annotations

import hashlib
import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Optional

from qpx.converters.openms_consensus.feature_dedup import feature_ids as _feature_ids

_NO_LINK = (0,)


def _key_digest(key: tuple) -> bytes:
    """128-bit digest of a PSM identity key: a tenth of the memory of the tuple itself."""
    return hashlib.blake2b(repr(key).encode(), digest_size=16).digest()


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
        digest = _key_digest(key)
        if digest not in self._keys:
            return False
        if not self.unassigned and digest in self._claimable:
            self.dropped_rows.add(self._claimable.pop(digest))
            return False
        return True

    def add(self, key) -> None:
        """Register an emitted key at the next output row."""
        digest = _key_digest(key)
        if self.unassigned and self.before_features:
            self._claimable[digest] = self.rows
        elif not self.unassigned:
            self._links[digest] = (None, self.rows)
        self._keys.add(digest)
        self.rows += 1

    def offer(self, key, rank: tuple, feature_id: Optional[int]) -> None:
        """Record a candidate link for an emitted key; keep it when it outranks the current one."""
        digest = _key_digest(key)
        entry = self._links.get(digest)
        if entry is None:
            return
        current, row = entry
        if current is None:
            self._links[digest] = (rank, row)
        elif rank > current:
            self._links[digest] = (rank, row)
            self.feature_id_patch[row] = feature_id

    def start_file(self) -> None:
        """Reset the per-file state before the next input."""
        # Runs never span inputs, so no later file can offer a link for these keys.
        self._claimable.clear()
        self._links.clear()
        self.unassigned = False
        self.before_features = True


@dataclass
class _RecordContext:
    """Run, protein and identity state shared by one input map's records."""

    map_info: dict
    protein_groups: tuple
    resolve_run: object
    seen: object
    enzyme: object
    dedup: object
    ms1: object = None

    @classmethod
    def from_map(cls, cm, structures, seen, dedup, sdrf_path):
        """Resolve one input's map columns, groups, runs and search enzyme."""
        from qpx.converters.openms_consensus.channels import Ms1Channels
        from qpx.converters.openms_consensus.feature_adapter import column_runs, feature_map_info, resolve_enzyme
        from qpx.converters.openms_consensus.pg_adapter import protein_group_maps
        from qpx.converters.openms_consensus.psm_adapter import _run_resolver

        map_info = feature_map_info(cm)
        protein_groups = protein_group_maps(cm) if "feature" in structures or "pg" in structures else (None, None)
        ms1 = Ms1Channels(cm, column_runs(cm))
        return cls(map_info, protein_groups, _run_resolver(cm), seen, resolve_enzyme(cm, sdrf_path), dedup, ms1)


def _cf_psm_records(cf, context, removed, include_unassigned_psms, duplicates, confidence):
    """Collect kept and recovered PSMs belonging to one consensus feature."""
    from qpx.converters.openms_consensus.psm_adapter import _cf_element_runs, psm_records_for_pid

    # Multi-run isobaric PIDs carry a local id_merge_index; the feature's
    # element runs disambiguate which run they belong to.
    cf_runs = _cf_element_runs(cf, context.map_info)
    pids = list(cf.getPeptideIdentifications())
    if include_unassigned_psms:
        pids.extend(removed)
    records = []
    for pid in pids:
        records.extend(
            psm_records_for_pid(
                pid,
                context.resolve_run,
                context.seen,
                enzyme=context.enzyme,
                cf_runs=cf_runs,
                duplicates=duplicates,
                confidence=confidence,
                ms1=context.ms1,
            )
        )
    return records


def _finalize_cf_records(cf, cf_feats, cf_psms, context, structures, duplicates):
    """Assign feature identities, link PSMs and retain preferred feature rows."""
    want_feature, want_psm = "feature" in structures, "psm" in structures
    feature_ids = _feature_ids(cf_feats) if cf_feats and (want_psm or context.dedup is not None) else None
    # Stamp psm.feature_id only when both views are emitted (the FK references a
    # feature row written in this dataset). feature.psm_ids is the computed inverse.
    if want_feature and want_psm:
        _link_psm_feature(
            cf_feats,
            cf_psms,
            feature_ids,
            links=context.seen if duplicates is not None else None,
            duplicates=duplicates or (),
            quality=float(cf.getQuality()),
        )
    if context.dedup is not None and cf_feats:
        cf_feats = context.dedup.keep(cf_feats, feature_ids, float(cf.getQuality()))
    return cf_feats if want_feature else [], cf_psms


def _cf_feature_psm_records(cf, context, structures, *, removed_index=None, confidence=None, include_unassigned_psms=True):
    """Feature + PSM records for one consensus feature, linked and deduplicated."""
    from qpx.converters.openms_consensus.feature_adapter import feature_records_for_cf, removed_same_peptide_ids

    want_feature, want_psm, want_pg = "feature" in structures, "psm" in structures, "pg" in structures
    removed = removed_same_peptide_ids(cf, removed_index)
    cf_feats = (
        feature_records_for_cf(
            cf,
            context.map_info,
            context.protein_groups[0],
            enzyme=context.enzyme,
            group_meta=context.protein_groups[1],
            resolve_run=context.resolve_run,
            removed_pids=removed,
            confidence=confidence,
            ms1=context.ms1,
        )
        if want_feature or want_pg
        else []
    )
    cf_psms: list[dict] = []
    link = want_feature and want_psm
    duplicates: Optional[list] = [] if link and isinstance(context.seen, _PsmKeys) else None
    if want_psm:
        cf_psms = _cf_psm_records(cf, context, removed, include_unassigned_psms, duplicates, confidence)
    return _finalize_cf_records(cf, cf_feats, cf_psms, context, structures, duplicates)
