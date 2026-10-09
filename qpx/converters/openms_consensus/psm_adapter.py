"""consensusXML -> QPX psm records.

Each ``PeptideIdentification`` (assigned to a consensus feature or unassigned) is
one spectrum match. We emit one psm record per hit with PK
``[peptidoform, charge, run_file_name, scan]``. The run is resolved from the
identification's global ``map_index`` (→ the map's file) unless it is a copy of
merged identifications whose ``id_merge_index`` places the spectrum in another
map run, the parent consensus feature's element runs, or — for an unassigned ID
in a merged multi-run consensusXML — its ``id_merge_index`` (→ the i-th merged
MS run, see :func:`_merge_index_runs`); falling back to the sole run of a
single-run file.
"""

from __future__ import annotations

import hashlib
import logging
import math
from collections import Counter

from qpx.converters.openms_consensus.feature_adapter import (
    _map_label,
    column_runs,
    feature_map_info,
    load_consensus_map,
    localization_scores,
    primary_run_stems,
    to_modifications,
    to_proforma,
)
from qpx.converters.openms_consensus.feature_adapter import (
    mass_error_ppm as _mass_error_ppm,
)
from qpx.converters.openms_consensus.feature_adapter import (
    pep_of as _pep_of,
)
from qpx.converters.openms_consensus.feature_adapter import (
    qvalue_meta_of as _qvalue_meta_of,
)
from qpx.converters.openms_consensus.protein_groups import identification_identifier
from qpx.converters.utils import safe_float
from qpx.core.cleavage import count_missed_cleavages
from qpx.core.scan import scan_format_from_native_id, scan_from_native_id

_log = logging.getLogger(__name__)

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
    """Parse numeric native ID components with the shared QPX scan convention."""
    ref = str(spectrum_ref or "")
    if not ref:
        return []
    scans = scan_from_native_id(ref)
    if scans:
        return scans
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


def _identifier_of(identification) -> str:
    """Identification run identifier, or ``""`` for objects that carry none."""
    try:
        return identification_identifier(identification) or ""
    except AttributeError:
        return ""


def _merge_index_runs(cm) -> dict[str, list[str]]:
    """Run stems indexed by ``id_merge_index``, per ProteinIdentification identifier.

    When OpenMS merges several identification runs, each moved PeptideIdentification
    is tagged with ``id_merge_index`` = the position of its ORIGINAL MS run in the
    primary MS run paths (the ``spectra_data`` StringList) of the ProteinIdentification
    it references — NOT a map-column index. Resolving through the PID's own
    ProteinIdentification (matched by identifier) mirrors OpenMS'
    ``IdentifierMSRunMapper``. Identifiers without recorded run paths are omitted.
    """
    runs: dict[str, list[str]] = {}
    for prot in cm.getProteinIdentifications():
        stems = primary_run_stems(prot)
        if stems:
            runs.setdefault(_identifier_of(prot), []).extend(stems)
    return runs


def _run_resolver(cm):
    """Build a callable mapping a PeptideIdentification to its run_file_name.

    ``map_index`` is a global map-column index (label-free: one map per run, so it
    is authoritative) — except when merged identifications were copied into every
    run's map (FeatureFinderIdentification fed a group-merged idXML, then linked).
    Then several ProteinIdentifications record the same multi-run ``spectra_data``
    (OpenMS' ``IdentifierMSRunMapper`` rejects that layout, so quantms output never
    has it), ``map_index`` names the map the copy landed in, and ``id_merge_index``
    names the spectrum's run. That run wins when it is one of the label-free map
    runs; when the two agree nothing changes.

    ``id_merge_index`` is otherwise only used without a ``map_index``: through the
    PID's own ProteinIdentification (:func:`_merge_index_runs`), falling back to
    all run paths in document order when the own list cannot hold the index
    (e.g. one ProteinIdentification per run). Callers with a consensus feature
    pass ``cf_runs`` (the runs its positive-intensity elements map to); with
    exactly one such run it is authoritative for that feature's PIDs.
    """
    headers = cm.getColumnHeaders()
    map_run = column_runs(cm)
    map_runs = set(map_run.values())
    label_free_maps = {idx for idx in headers if _map_label(getattr(headers[idx], "label", "")) == "LFQ"}
    sole_run = next(iter(map_runs)) if len(map_runs) == 1 else None
    runs_by_identifier = _merge_index_runs(cm)
    run_lists = [primary_run_stems(prot) for prot in cm.getProteinIdentifications()]
    all_runs = [run for runs in run_lists for run in runs]
    list_counts = Counter(tuple(runs) for runs in run_lists if len(runs) > 1)
    copied_lists = {runs for runs, count in list_counts.items() if count > 1}

    def merged_run(pid) -> tuple[str | None, bool]:
        """The id_merge_index run, and whether it comes from a copied multi-run list."""
        if not pid.metaValueExists("id_merge_index"):
            return None, False
        idx = int(pid.getMetaValue("id_merge_index"))
        own = runs_by_identifier.get(_identifier_of(pid), [])
        if 0 <= idx < len(own):
            return own[idx], tuple(own) in copied_lists
        return (all_runs[idx] if 0 <= idx < len(all_runs) else None), False

    def resolve(pid, cf_runs=None) -> str | None:
        origin, copied_ids = merged_run(pid)
        if pid.metaValueExists("map_index"):
            idx = int(pid.getMetaValue("map_index"))
            run = map_run.get(idx)
            if run:
                moved = copied_ids and idx in label_free_maps and origin in map_runs
                return origin if moved else run
        # Assigned PID: fall back to the consensus feature's single element run.
        if cf_runs and len(cf_runs) == 1:
            return next(iter(cf_runs))
        # Unassigned PID in a merged multi-run consensusXML: id_merge_index selects
        # the original MS run, NOT a map column.
        if origin:
            return origin
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


def _psm_additional_scores(primary, hits, score, score_type, score_is_qvalue, higher_better, peptide_level=False):
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
    if not score_is_qvalue and not peptide_level:
        # A search-score primary (e.g. COMET:xcorr) leaves the Percolator q-value in the hit meta values.
        qvalue = _qvalue_meta_of(primary)
        if qvalue is not None:
            additional_scores.append({"score_name": "q-value", "score_value": qvalue, "higher_better": False})
    # Preserve the other colliding hits' (i.e. the other engines') search scores
    # so a comet+msgf merged spectrum does not lose either engine's score.
    for other in hits:
        if other is primary:
            continue
        other_score = safe_float(other.getScore())
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


def psm_records_for_pid(
    pid, resolve_run, seen: set[tuple], cf_runs=None, enzyme=None, duplicates=None, confidence=None
) -> list[dict]:
    """PSM records for one PeptideIdentification (deduped via the shared ``seen`` set).

    ``duplicates``, when a list, receives ``(key, rt, observed_mz)`` for every hit
    whose key was already emitted, so the caller can weigh this copy's feature link.

    ``confidence`` (a :class:`PeptideLevelConfidence`) marks peptide-level Percolator
    runs: their PEP/q-value belong to the peptide, so the PSM PEP stays null and the
    peptide q-value goes to ``additional_scores`` as ``peptide_qvalue``.

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
            if duplicates is not None:
                duplicates.append((key, float(pid.getRT()) if pid.getRT() else None, obs_mz))
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
        peptide_level = confidence is not None and confidence.applies(pid)
        pep = None if peptide_level else _pep_of(primary)
        score = safe_float(primary.getScore())
        additional_scores, site_scores = _psm_additional_scores(
            primary, hits, score, score_type, score_is_qvalue, higher_better, peptide_level=peptide_level
        )
        if peptide_level:
            peptide_qvalue = confidence.of(pid, to_proforma(seq_obj))[1]
            if peptide_qvalue is not None:
                additional_scores.append({"score_name": "peptide_qvalue", "score_value": peptide_qvalue, "higher_better": False})
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


def collect_psm_scan_formats(identifications, formats: set[str | None], *, enabled=True) -> None:
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


def uniform_scan_format(formats: set[str | None]) -> str | None:
    """Return the declaration only after the full input confirms one known format."""
    return next(iter(formats)) if len(formats) == 1 else None


def declare_psm_scan_format(writer, formats: set[str | None]) -> None:
    """Finalize a streamed PSM writer's confirmed declaration before it closes."""
    scan_format = uniform_scan_format(formats)
    if writer is not None and scan_format is not None:
        writer.set_scan_format(scan_format)
