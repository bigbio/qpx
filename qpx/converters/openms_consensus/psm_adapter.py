"""consensusXML -> QPX psm records.

Each ``PeptideIdentification`` (assigned to a consensus feature or unassigned) is
one spectrum match. We emit one psm record per hit with PK
``[run_file_name, scan, peptidoform, charge]``. The run is resolved from the
identification's global ``map_index`` (→ the map's file) unless it is a copy of
merged identifications whose ``id_merge_index`` places the spectrum in another
map run, the parent consensus feature's element runs, or — for an unassigned ID
in a merged multi-run consensusXML — its ``id_merge_index`` (→ the i-th merged
MS run, see :func:`runs._merge_index_runs`); falling back to the sole run of a
single-run file.
"""

from __future__ import annotations

import logging
import math
from typing import NamedTuple

from qpx.converters.openms_consensus.channels import label_state_params, matched_sequence
from qpx.converters.openms_consensus.feature_adapter import (
    is_decoy_hit,
    localization_scores,
    to_modifications,
    to_proforma,
)
from qpx.converters.openms_consensus.feature_adapter import (
    pep_of as _pep_of,
)
from qpx.converters.openms_consensus.feature_adapter import (
    qvalue_meta_of as _qvalue_meta_of,
)
from qpx.converters.openms_consensus.runs import _scan_of
from qpx.converters.utils import canonical_accession, safe_float
from qpx.converters.utils import mass_error_ppm as _mass_error_ppm
from qpx.core.cleavage import count_missed_cleavages
from qpx.core.scan import scan_format_from_native_id
from qpx.core.scores import normalize_score_name

_log = logging.getLogger(__name__)


def _protein_accessions(hits) -> list[str] | None:
    """Return the distinct protein evidence carried by all collapsed hits."""
    accessions: list[str] = []
    seen: set[str] = set()
    for hit in hits:
        for evidence in hit.getPeptideEvidences():
            raw = evidence.getProteinAccession()
            accession = canonical_accession(raw.decode() if isinstance(raw, (bytes, bytearray)) else str(raw or ""))
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
        # Route the identification score into additional_scores; a q-value named
        # qvalue fills psm.qvalue when written (scores.md).
        name = "qvalue" if score_is_qvalue else (normalize_score_name(score_type) or "search_score")
        additional_scores.append({"score_name": name, "score_value": score, "higher_better": higher_better})
    if not score_is_qvalue and not peptide_level:
        # A search-score primary (e.g. COMET:xcorr) leaves the Percolator q-value in the hit meta values.
        qvalue = _qvalue_meta_of(primary)
        if qvalue is not None:
            additional_scores.append({"score_name": "qvalue", "score_value": qvalue, "higher_better": False})
    # Preserve the other colliding hits' (i.e. the other engines') search scores
    # so a comet+msgf merged spectrum does not lose either engine's score.
    for other in hits:
        if other is primary:
            continue
        other_score = safe_float(other.getScore())
        if other_score is not None:
            base = "qvalue" if score_is_qvalue else (normalize_score_name(score_type) or "search_score")
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


class _PsmContext(NamedTuple):
    """Spectrum location and score interpretation shared by one PID's hits."""

    run: str
    scan: list[int]
    obs_mz: float
    score_type: str
    score_is_qvalue: bool
    higher_better: bool
    ms1: object = None


def _psm_hit_groups(pid, run, scan) -> tuple[dict, list]:
    """Group colliding hits while retaining their first-seen order."""
    groups: dict[tuple, list] = {}
    order: list[tuple] = []
    for hit in pid.getHits():
        peptidoform = to_proforma(matched_sequence(hit))
        charge = int(hit.getCharge() or 0)
        key = (peptidoform, charge, run, tuple(scan))
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(hit)

    return groups, order


def _psm_confidence_scores(pid, primary, hits, context, confidence, seq_obj):
    """Build source-scoped confidence and additional score fields."""
    peptide_level = confidence is not None and confidence.applies(pid)
    pep = None if peptide_level else _pep_of(primary)
    score = safe_float(primary.getScore())
    if peptide_level and context.score_is_qvalue:
        # Under peptide-level FDR a q-value score is the peptide's, 1.0 on all but its best PSM:
        # the PSM gets it as peptide_qvalue below, not as a score of its own or of its colliding hits.
        score, hits = None, [primary]
    additional_scores, site_scores = _psm_additional_scores(
        primary, hits, score, context.score_type, context.score_is_qvalue, context.higher_better, peptide_level=peptide_level
    )
    if peptide_level:
        peptide_qvalue = confidence.of(pid, to_proforma(seq_obj))[1]
        if peptide_qvalue is not None:
            additional_scores.append({"score_name": "peptide_qvalue", "score_value": peptide_qvalue, "higher_better": False})
    return pep, additional_scores, site_scores


def _psm_record(pid, hits, context, enzyme, confidence) -> dict:
    """Build one collapsed PSM from its peptide, spectrum and score fields."""
    # Keep the lowest-PEP hit; hits without a PEP sort last (worst). ``min`` is
    # stable, so ties (and the common single-hit case) keep the first hit,
    # preserving existing output exactly.
    primary = min(hits, key=lambda h: _pep_of(h) if _pep_of(h) is not None else math.inf)
    # An MS1-labelled PSM is the partner the spectrum was matched with; confidence
    # belongs to the peptide identity, which is the hit's own sequence.
    seq_obj = matched_sequence(primary)
    peptidoform = to_proforma(seq_obj)
    charge = int(primary.getCharge() or 0)
    calc_mz = float(seq_obj.getMZ(charge)) if charge > 0 else None
    is_decoy = is_decoy_hit(primary)
    pep, additional_scores, site_scores = _psm_confidence_scores(pid, primary, hits, context, confidence, primary.getSequence())
    modifications = to_modifications(seq_obj, site_scores)
    record = {
        "sequence": seq_obj.toUnmodifiedString(),
        "peptidoform": peptidoform,
        "modifications": modifications,
        "charge": charge,
        "run_file_name": context.run,
        "scan": context.scan,
        "rt": float(pid.getRT()) if pid.getRT() else None,
        "calculated_mz": calc_mz,
        "observed_mz": context.obs_mz,
        "mass_error_ppm": _mass_error_ppm(calc_mz, context.obs_mz),
        "missed_cleavages": (count_missed_cleavages(seq_obj.toUnmodifiedString(), enzyme) if enzyme else None),
        "posterior_error_probability": pep,
        "additional_scores": additional_scores or None,
        "is_decoy": is_decoy,
        "protein_accessions": _protein_accessions(hits),
    }
    if context.ms1:
        record["label"] = context.ms1.identified_label(pid, primary, context.run)
        record["cv_params"] = label_state_params(primary)
    return record


def psm_records_for_pid(
    pid, resolve_run, seen: set[tuple], cf_runs=None, enzyme=None, duplicates=None, confidence=None, ms1=None
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

    ``ms1`` (:class:`~qpx.converters.openms_consensus.channels.Ms1Channels`) names the
    channel of an MS1-labelled PSM.
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
    # When the identification score is a q-value, the hit score is that q-value
    # (OpenMS FDR output), the peptide's under peptide-level FDR; otherwise it is a search score.
    score_type = str(pid.getScoreType() or "")
    score_is_qvalue = score_type.lower() in ("q-value", "qvalue", "fdr")
    # Group hits by identity key. A PeptideIdentification usually holds one hit
    # (quantms ships Percolator-merged single-hit PIDs), but comet+msgf merged
    # spectra can carry several hits for the SAME peptidoform/charge/scan. Those
    # collisions must resolve to one row (lowest PEP), keeping the other engines'
    # search scores; genuinely different peptidoforms map to different keys and are
    # all emitted as distinct PSMs.
    groups, order = _psm_hit_groups(pid, run, scan)

    records: list[dict] = []
    higher_better = bool(pid.isHigherScoreBetter())
    context = _PsmContext(run, scan, obs_mz, score_type, score_is_qvalue, higher_better, ms1)
    for key in order:
        if key in seen:
            if duplicates is not None:
                duplicates.append((key, float(pid.getRT()) if pid.getRT() else None, obs_mz))
            continue
        seen.add(key)
        records.append(_psm_record(pid, groups[key], context, enzyme, confidence))
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
