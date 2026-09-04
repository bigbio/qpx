"""consensusXML -> QPX psm records.

Each ``PeptideIdentification`` (assigned to a consensus feature or unassigned) is
one spectrum match. We emit one psm record per hit with PK
``[peptidoform, charge, run_file_name, scan]``. The run is resolved from the
identification's global ``map_index`` (→ the map's file), the parent consensus
feature's element runs, or — for an unassigned ID in a merged multi-run
consensusXML — its ``id_merge_index`` (→ the i-th merged MS run, see
:func:`_merge_index_runs`); falling back to the sole run of a single-run file.
"""

from __future__ import annotations

import hashlib
import logging
import math
import re

from qpx.converters.openms_consensus.feature_adapter import (
    _run_stem,
    feature_map_info,
    load_consensus_map,
    localization_scores,
    to_modifications,
    to_proforma,
)
from qpx.converters.openms_consensus.feature_adapter import (
    mass_error_ppm as _mass_error_ppm,
)
from qpx.converters.openms_consensus.feature_adapter import (
    pep_of as _pep_of,
)
from qpx.core.cleavage import count_missed_cleavages

_log = logging.getLogger(__name__)

# Thermo/Bruker/single-peak-list nativeIDs expose the spectrum ordinal directly.
_SCAN_RE = re.compile(r"(?:scan|index|spectrum)=(\d+)", re.IGNORECASE)
# Sciex WIFF nativeIDs (``sample=.. period=.. cycle=.. experiment=..``) carry no
# scan/index/spectrum token; the cycle is the acquisition ordinal (scan-equivalent).
_CYCLE_RE = re.compile(r"cycle=(\d+)", re.IGNORECASE)
# scan is a list<int32>; keep any surrogate within the signed 32-bit range.
_INT32_MASK = 0x7FFFFFFF


def _surrogate_scan(spectrum_ref: str) -> int:
    """Deterministic int32 surrogate ordinal for a nativeID with no scan token.

    Derived only from the spectrum reference string, so the same spectrum maps to
    the same value in both the pyopenms and streaming paths. Used as a last resort
    (instead of dropping the PSM) for nativeID schemes that expose no recognizable
    ordinal at all.
    """
    digest = hashlib.blake2s(str(spectrum_ref).encode("utf-8"), digest_size=4).digest()
    return int.from_bytes(digest, "big") & _INT32_MASK


def _scan_of(spectrum_ref: str) -> list[int]:
    """Parse the scan number(s) from a spectrum reference into a list<int>.

    Recognizes the common ``scan=``/``index=``/``spectrum=`` tokens (Thermo,
    Bruker, single peak lists, Waters), falls back to the Sciex ``cycle=``
    ordinal, and finally to a deterministic surrogate for nativeID schemes with
    no recognizable ordinal. A completely empty reference returns ``[]`` (the
    caller skips those PSMs). Shared with the feature adapter so psm.scan and
    feature.scan stay consistent.
    """
    ref = str(spectrum_ref or "")
    if not ref:
        return []
    scans = [int(m) for m in _SCAN_RE.findall(ref)]
    if scans:
        return scans
    cycles = [int(m) for m in _CYCLE_RE.findall(ref)]
    if cycles:
        return cycles
    return [_surrogate_scan(ref)]


def _protein_accessions(hits) -> list[str] | None:
    """Return the distinct protein evidence carried by all collapsed hits."""
    accessions: list[str] = []
    seen: set[str] = set()
    for hit in hits:
        for evidence in hit.getPeptideEvidences():
            raw = evidence.getProteinAccession()
            accession = raw.decode() if isinstance(raw, (bytes, bytearray)) else str(raw or "")
            if accession and accession not in seen:
                seen.add(accession)
                accessions.append(accession)
    return accessions or None


def _append_unique_score(additional_scores: list[dict], name: str, value: float, higher_better: bool) -> None:
    """Append a score, disambiguating the name if it already occurs.

    Colliding hits from different search engines share the identification score
    type, so a plain name would repeat; suffix ``_2``, ``_3``, ... keeps every
    engine's score without overwriting.
    """
    existing = {s["score_name"] for s in additional_scores}
    unique = name
    n = 2
    while unique in existing:
        unique = f"{name}_{n}"
        n += 1
    additional_scores.append({"score_name": unique, "score_value": value, "higher_better": higher_better})


def _cf_element_runs(cf, map_info: dict[int, tuple[str, str]]) -> set[str]:
    """Runs the consensus feature's elements map to (positive intensity).

    Mirrors the feature adapter's run attribution (``feature_records_for_cf``
    skips non-positive elements), so a PSM record lands on the same
    ``run_file_name`` as the feature record it links to.
    """
    runs: set[str] = set()
    for sub in cf.getFeatureList():
        if float(sub.getIntensity()) <= 0:
            continue
        run = map_info.get(sub.getMapIndex(), (None, None))[0]
        if run is not None:
            runs.add(run)
    return runs


def _merge_index_runs(cm) -> list[str]:
    """Run stems indexed by ``id_merge_index`` (the merged-run ordering).

    When OpenMS merges several identification runs into one consensusXML, each
    moved PeptideIdentification is tagged with ``id_merge_index`` = the position
    of its ORIGINAL MS run in the merge order — NOT a map-column index. OpenMS
    records that order as the merged ProteinIdentification's primary MS run paths
    (the ``spectra_data`` StringList), so concatenating those paths across the
    ProteinIdentifications, in document order, yields ``id_merge_index -> run``.

    Returns an empty list when no primary MS run path is recorded (e.g. a bare
    single-run file), in which case the resolver falls back to the sole run.
    """
    runs: list[str] = []
    for prot in cm.getProteinIdentifications():
        get_paths = getattr(prot, "getPrimaryMSRunPath", None)
        if get_paths is None:
            continue
        paths: list = []
        get_paths(paths)
        for p in paths:
            runs.append(_run_stem(p.decode() if isinstance(p, (bytes, bytearray)) else str(p)))
    return runs


def _run_resolver(cm):
    """Build a callable mapping a PeptideIdentification to its run_file_name.

    ``map_index`` is a global map-column index (label-free: one map per run, so it
    is authoritative). ``id_merge_index``, however, is a per-run MERGE index in a
    multi-run consensusXML (isobaric TMT/iTRAQ, or any merged file) — it selects
    the i-th original MS run, not the i-th map column — so it is resolved through
    the merged-run ordering (:func:`_merge_index_runs`), never the map-column dict.
    Callers with a consensus feature pass ``cf_runs`` (the runs its
    positive-intensity elements map to); with exactly one such run it is
    authoritative for that feature's PIDs.
    """
    headers = cm.getColumnHeaders()
    map_run = {idx: _run_stem(headers[idx].filename) for idx in headers}
    distinct = sorted(set(map_run.values()))
    sole_run = distinct[0] if len(distinct) == 1 else None
    merge_runs = _merge_index_runs(cm)

    def resolve(pid, cf_runs=None) -> str | None:
        # Global map index, when present, is always authoritative.
        if pid.metaValueExists("map_index"):
            run = map_run.get(int(pid.getMetaValue("map_index")))
            if run:
                return run
        # Assigned PID: fall back to the consensus feature's single element run.
        if cf_runs and len(cf_runs) == 1:
            return next(iter(cf_runs))
        # Unassigned PID in a merged multi-run consensusXML: id_merge_index selects
        # the i-th original MS run (merge order), resolved via the merged
        # ProteinIdentification's primary MS run paths — NOT the map-column dict.
        if pid.metaValueExists("id_merge_index"):
            idx = int(pid.getMetaValue("id_merge_index"))
            if 0 <= idx < len(merge_runs):
                return merge_runs[idx]
        return sole_run

    return resolve


def consensus_psms_to_records(consensusxml_path: str | None = None, cm=None) -> list[dict]:
    """Return QPX psm record dicts extracted from a consensusXML.

    Pass either ``consensusxml_path`` (loaded here) or an already-loaded ``cm``.
    """
    cm = cm if cm is not None else load_consensus_map(consensusxml_path)
    resolve_run = _run_resolver(cm)
    records: list[dict] = []
    seen: set[tuple] = set()
    for pid in cm.getUnassignedPeptideIdentifications():
        records.extend(psm_records_for_pid(pid, resolve_run, seen))
    map_info = feature_map_info(cm)
    for cf in cm:
        cf_runs = _cf_element_runs(cf, map_info)
        for pid in cf.getPeptideIdentifications():
            records.extend(psm_records_for_pid(pid, resolve_run, seen, cf_runs=cf_runs))
    return records


def _psm_additional_scores(primary, hits, score, score_type, score_is_qvalue, higher_better):
    """Assemble one PSM's ``additional_scores`` and its localisation site scores.

    Extracted from :func:`psm_records_for_pid` to keep that function within the
    project's complexity limits; the behaviour is unchanged.

    Returns ``(additional_scores, site_scores)`` — the latter feeds
    :func:`to_modifications` so per-site localisation scores land on the
    modification they belong to.
    """
    additional_scores = []
    if score is not None:
        # Route the identification score into additional_scores: the psm schema
        # has no dedicated q-value column.
        name = "q-value" if score_is_qvalue else (score_type or "search_score")
        additional_scores.append({"score_name": name, "score_value": score, "higher_better": higher_better})
    # Preserve the other colliding hits' (i.e. the other engines') search scores
    # so a comet+msgf merged spectrum does not lose either engine's score.
    for other in hits:
        if other is primary:
            continue
        other_score = float(other.getScore()) if other.getScore() is not None else None
        if other_score is not None:
            base = "q-value" if score_is_qvalue else (score_type or "search_score")
            _append_unique_score(additional_scores, base, other_score, higher_better)
    if primary.metaValueExists("consensus_support"):
        additional_scores.append(
            {
                "score_name": "consensus_support",
                "score_value": float(primary.getMetaValue("consensus_support")),
                "higher_better": True,
            }
        )
    loc_scores, site_scores = localization_scores(primary)
    if loc_scores:
        additional_scores.extend(loc_scores)
    return additional_scores, site_scores


def psm_records_for_pid(pid, resolve_run, seen: set[tuple], cf_runs=None, enzyme=None) -> list[dict]:
    """PSM records for one PeptideIdentification (deduped via the shared ``seen`` set).

    ``cf_runs`` is the consensus feature's element-run set (passed for assigned
    PIDs); it lets :func:`_run_resolver` attribute a PID whose ``id_merge_index``
    is only a local per-run channel index to the correct run.
    """
    run = resolve_run(pid, cf_runs=cf_runs)
    if run is None:
        return []
    spectrum_ref = pid.getSpectrumReference() if hasattr(pid, "getSpectrumReference") else ""
    if not spectrum_ref and pid.metaValueExists("spectrum_reference"):
        spectrum_ref = pid.getMetaValue("spectrum_reference")
    scan = _scan_of(spectrum_ref)
    if not scan:
        # No spectrum reference at all (or nothing keyable): nothing to key on, skip.
        _log.debug("Skipping consensusXML PSM with empty spectrum_reference: %r", spectrum_ref)
        return []
    obs_mz = float(pid.getMZ()) if pid.getMZ() else 0.0
    # When the identification score IS the q-value, the hit score is the peptide
    # q-value (OpenMS FDR output); otherwise it is a search score.
    score_type = str(pid.getScoreType() or "")
    score_is_qvalue = score_type.lower() in ("q-value", "qvalue", "fdr")
    # Group hits by identity key. A PeptideIdentification usually holds one hit
    # (quantms ships Percolator-merged single-hit PIDs), but comet+msgf merged
    # spectra can carry several hits for the SAME peptidoform/charge/scan. Those
    # collisions must resolve to one row (lowest PEP), keeping the other engines'
    # search scores; genuinely different peptidoforms map to different keys and are
    # all emitted as distinct PSMs.
    groups: dict[tuple, list] = {}
    order: list[tuple] = []
    for hit in pid.getHits():
        peptidoform = to_proforma(hit.getSequence())
        charge = int(hit.getCharge() or 0)
        key = (peptidoform, charge, run, tuple(scan))
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(hit)

    records: list[dict] = []
    higher_better = bool(pid.isHigherScoreBetter())
    for key in order:
        if key in seen:
            continue
        seen.add(key)
        hits = groups[key]
        # Keep the lowest-PEP hit; hits without a PEP sort last (worst). ``min`` is
        # stable, so ties (and the common single-hit case) keep the first hit,
        # preserving existing output exactly.
        primary = min(hits, key=lambda h: _pep_of(h) if _pep_of(h) is not None else math.inf)
        seq_obj = primary.getSequence()
        peptidoform = to_proforma(seq_obj)
        charge = int(primary.getCharge() or 0)
        calc_mz = float(seq_obj.getMZ(charge)) if charge > 0 else None
        is_decoy = primary.metaValueExists("target_decoy") and "decoy" in str(primary.getMetaValue("target_decoy")).lower()
        pep = _pep_of(primary)
        score = float(primary.getScore()) if primary.getScore() is not None else None
        additional_scores, site_scores = _psm_additional_scores(primary, hits, score, score_type, score_is_qvalue, higher_better)
        modifications = to_modifications(seq_obj, site_scores)
        records.append(
            {
                "sequence": seq_obj.toUnmodifiedString(),
                "peptidoform": peptidoform,
                "modifications": modifications,
                "charge": charge,
                "run_file_name": run,
                "scan": scan,
                "rt": float(pid.getRT()) if pid.getRT() else None,
                "calculated_mz": calc_mz,
                "observed_mz": obs_mz,
                "mass_error_ppm": _mass_error_ppm(calc_mz, obs_mz),
                "missed_cleavages": (count_missed_cleavages(seq_obj.toUnmodifiedString(), enzyme) if enzyme else None),
                "posterior_error_probability": pep,
                "additional_scores": additional_scores or None,
                "is_decoy": is_decoy,
                "protein_accessions": _protein_accessions(hits),
            }
        )
    return records
