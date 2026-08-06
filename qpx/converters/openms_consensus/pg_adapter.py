"""consensusXML -> QPX pg (protein group) records.

Protein groups come from the OpenMS protein-inference graph: indistinguishable
protein groups when present, else each protein hit as a singleton group. Peptide
counts come from the peptide→protein evidence links.

**Interim protein intensity**: the consensusXML has no protein-level abundance
(that lived in the mzTab), so until OpenMS ``-out_qpx`` provides the authoritative
number we roll one up ourselves — the **unnormalized sum of the group's unique
peptides** (quantms ``unique_peptides`` policy, no normalization), per
``(protein group, grouped_runs unit, label)``. It is *not* the authoritative
quant, so every quantified row carries a ``quantification_method`` cv_param; set
``top=3`` to mirror the ProteomicsLFQ/IsobaricWorkflow default instead of summing
all peptides. Rows stay ``intensity``-null where a group has no unique-peptide
signal. One row per ``(protein group, grouped_runs, label)`` — the quantification
unit(s) from the SDRF when given, otherwise the whole run set as a single unit.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Optional

from qpx.converters.channel_labels import fraction_groups_from_sdrf
from qpx.converters.openms_consensus.feature_adapter import (
    _map_label,
    _run_stem,
    load_consensus_map,
    to_proforma,
)
from qpx.converters.openms_consensus.psm_adapter import _run_resolver

_GENE_RE = re.compile(r"GN=([^\s]+)")


class _ProteinMaps:
    """Peptide/feature/run evidence indexed by accession, with the inverse maps.

    ``acc_to_*`` power the per-group unions; ``pep_to_accs`` / ``feat_to_accs`` are
    the inverse (peptide/feature -> the accessions it maps to) used to decide which
    peptides/features are *unique* to a group (map only to proteins in it) — the
    same ``unique_peptides`` policy the intensity rollup uses.
    """

    def __init__(self):
        self.acc_to_pep: dict[str, set[str]] = defaultdict(set)
        self.acc_to_runs: dict[str, set[str]] = defaultdict(set)
        self.acc_to_feat: dict[str, set[tuple]] = defaultdict(set)
        self.pep_to_accs: dict[str, set[str]] = defaultdict(set)
        self.feat_to_accs: dict[tuple, set[str]] = defaultdict(set)


def _evidence_accession(ev) -> str | None:
    acc = ev.getProteinAccession()
    if not acc:
        return None
    return acc.decode() if isinstance(acc, bytes) else str(acc)


def _collect_pid(pid, runs: set[str], m: _ProteinMaps) -> None:
    """Fold one PeptideIdentification's evidence into the accession maps."""
    for hit in pid.getHits():
        seq_obj = hit.getSequence()
        seq = seq_obj.toUnmodifiedString()
        feat = (to_proforma(seq_obj), int(hit.getCharge() or 0))
        for ev in hit.getPeptideEvidences():
            acc = _evidence_accession(ev)
            if acc:
                m.acc_to_pep[acc].add(seq)
                m.acc_to_feat[acc].add(feat)
                m.acc_to_runs[acc] |= runs
                m.pep_to_accs[seq].add(acc)
                m.feat_to_accs[feat].add(acc)


def _cf_runs(cf, map_run: dict[int, str]) -> set[str]:
    runs = {map_run.get(sub.getMapIndex()) for sub in cf.getFeatureList() if sub.getIntensity() > 0}
    runs.discard(None)
    return runs


def accumulate_cf_maps(cf, map_run: dict[int, str], m: _ProteinMaps) -> None:
    """Fold one consensus feature's assigned IDs into the accession maps (per-cf)."""
    runs = _cf_runs(cf, map_run)
    for pid in cf.getPeptideIdentifications():
        _collect_pid(pid, runs, m)


def accumulate_unassigned_maps(pid, resolve_run, m: _ProteinMaps) -> None:
    """Fold one unassigned PeptideIdentification into the accession maps."""
    run = resolve_run(pid)
    _collect_pid(pid, {run} if run else set(), m)


def _protein_maps(cm) -> _ProteinMaps:
    """Index peptide/feature/run evidence by accession (see :class:`_ProteinMaps`)."""
    headers = cm.getColumnHeaders()
    map_run = {i: _run_stem(headers[i].filename) for i in headers}
    resolve_run = _run_resolver(cm)
    m = _ProteinMaps()
    for cf in cm:  # assigned IDs: runs are the consensus feature's member maps
        accumulate_cf_maps(cf, map_run, m)
    for pid in cm.getUnassignedPeptideIdentifications():  # run from map_index or id_merge_index (merge order)
        accumulate_unassigned_maps(pid, resolve_run, m)
    return m


def _is_decoy_accession(acc: str) -> bool:
    return str(acc).upper().startswith(("DECOY", "REV_", "RANDOM_"))


def _is_contaminant(acc: str) -> bool:
    return "CONTAM" in str(acc).upper()


def _acc_str(acc) -> str:
    return acc.decode() if isinstance(acc, bytes) else str(acc)


def _protein_hit_meta(prot) -> tuple[dict[str, bool], dict[str, float], dict[str, str]]:
    """Per-accession decoy flag, q-value (when the protein score IS a q-value), gene."""
    score_is_qvalue = str(prot.getScoreType() or "").lower() in ("q-value", "qvalue", "fdr")
    acc_decoy: dict[str, bool] = {}
    acc_qvalue: dict[str, float] = {}
    acc_gene: dict[str, str] = {}
    for hit in prot.getHits():
        acc = _acc_str(hit.getAccession())
        if hit.metaValueExists("target_decoy"):
            acc_decoy[acc] = "decoy" in str(hit.getMetaValue("target_decoy")).lower()
        if score_is_qvalue and hit.getScore() is not None:
            acc_qvalue[acc] = float(hit.getScore())
        gene = _GENE_RE.search(str(hit.getDescription() or "")) if hasattr(hit, "getDescription") else None
        if gene:
            acc_gene[acc] = gene.group(1)
    return acc_decoy, acc_qvalue, acc_gene


def _build_groups(prot) -> list[list[str]]:
    """Protein groups: indistinguishable groups plus a singleton for every hit
    not covered by one.

    ``getIndistinguishableProteins()`` only returns the *grouped* proteins, so a
    hit that OpenMS left ungrouped would be lost if we returned those groups
    alone. Append singletons for the uncovered hits (this also reproduces the
    all-singletons behaviour when there are no indistinguishable groups).
    """
    groups: list[list[str]] = []
    covered: set[str] = set()
    for grp in prot.getIndistinguishableProteins():
        accs = [_acc_str(a) for a in grp.accessions]
        if accs:
            groups.append(accs)
            covered.update(accs)
    for hit in prot.getHits():
        acc = _acc_str(hit.getAccession())
        if acc not in covered:
            groups.append([acc])
            covered.add(acc)
    return groups


def _merge_protein_ids(cm) -> tuple[dict[str, bool], dict[str, float], dict[str, str], list[list[str]]]:
    """Merge every ProteinIdentification's metadata + groups from the inference graph.

    A merged multi-run consensusXML carries one ProteinIdentification per run, so
    process them all: merge per-accession decoy/gene, keep the best (min) q-value,
    and concatenate groups while deduplicating equivalent accession sets.
    """
    acc_decoy: dict[str, bool] = {}
    acc_qvalue: dict[str, float] = {}
    acc_gene: dict[str, str] = {}
    groups: list[list[str]] = []
    seen_groups: set[frozenset] = set()
    for prot in cm.getProteinIdentifications():
        decoy, qvalue, gene = _protein_hit_meta(prot)
        acc_decoy.update(decoy)
        acc_gene.update(gene)
        for acc, qv in qvalue.items():
            acc_qvalue[acc] = min(acc_qvalue[acc], qv) if acc in acc_qvalue else qv
        for grp in _build_groups(prot):
            key = frozenset(grp)
            if key not in seen_groups:
                seen_groups.add(key)
                groups.append(grp)
    return acc_decoy, acc_qvalue, acc_gene, groups


def accession_to_anchor(cm) -> dict[str, str]:
    """Map each accession to its protein-group leader (the pg ``anchor_protein``).

    OpenMS defines the group leader as the group's first accession; the pg view
    uses it as ``anchor_protein``. Sharing this map with the feature adapter lets
    a feature stamp the SAME anchor as its protein group, so the documented
    feature->pg join on ``anchor_protein`` holds for multi-accession groups
    (a peptide's evidence order alone does not identify the leader).
    """
    _, _, _, groups = _merge_protein_ids(cm)
    return {acc: grp[0] for grp in groups for acc in grp}


def _map_info(cm) -> dict[int, tuple[str, str]]:
    """Map index -> (run_file_name, channel label), matching the feature adapter.

    Isobaric channels are detected from the map label (``tmt6plex_126`` ->
    ``TMT126``); everything else is ``LFQ``. ``experiment_type`` is not used — it
    is ``"label-free"`` even for real quantms TMT output.
    """
    headers = cm.getColumnHeaders()
    return {i: (_run_stem(headers[i].filename), _map_label(headers[i].label)) for i in headers}


def _peptide_intensities(cm, map_info: dict[int, tuple[str, str]]) -> dict[tuple[str, str, str], float]:
    """Return ``(peptide, run, label) -> summed feature intensity``.

    Peptides are keyed by unmodified sequence (matching :func:`_protein_maps`, so
    the same peptide->accession map drives the unique-to-group test). Feature
    intensities are summed per (peptide, run, label) — charge states / peptidoforms
    of one sequence roll up together.
    """
    pep_intensity: dict[tuple[str, str, str], float] = defaultdict(float)
    for cf in cm:
        accumulate_cf_intensity(cf, map_info, pep_intensity)
    return pep_intensity


def accumulate_cf_intensity(cf, map_info: dict[int, tuple[str, str]], pep_intensity: dict) -> None:
    """Sum one consensus feature's per-map intensities into ``(seq, run, label)`` (per-cf)."""
    pids = cf.getPeptideIdentifications()
    if not pids or not pids[0].getHits():
        return
    seq = pids[0].getHits()[0].getSequence().toUnmodifiedString()
    for sub in cf.getFeatureList():
        inten = float(sub.getIntensity())
        if inten <= 0:
            continue
        run, label = map_info.get(sub.getMapIndex(), (None, None))
        if run is not None:
            pep_intensity[(seq, run, label)] += inten


def _protein_intensity(group_peps, group_accs, unit, label, pep_intensity, pep_accs, top) -> Optional[float]:
    """Interim protein-group intensity for one ``(unit, label)``.

    Sum of the group's **unique** peptides — those mapping only to proteins in the
    group (the quantms ``unique_peptides`` policy) — with **no normalization**.
    ``top > 0`` keeps only the N most-abundant peptides (quantms/ProteomicsLFQ
    default is 3); ``top = 0`` sums all. Returns ``None`` when the group has no
    unique-peptide signal in this unit+label (stays identification-only).
    """
    abundances = []
    for pep in group_peps:
        if not pep_accs.get(pep, set()).issubset(group_accs):
            continue  # shared outside the group -> excluded by unique_peptides
        ab = sum(pep_intensity.get((pep, run, label), 0.0) for run in unit)
        if ab > 0:
            abundances.append(ab)
    if not abundances:
        return None
    abundances.sort(reverse=True)
    if top and top > 0:
        abundances = abundances[:top]
    return sum(abundances)


def consensus_protein_groups_to_records(
    consensusxml_path: str | None = None,
    sdrf_path: Optional[str] = None,
    cm=None,
    top: int = 0,
) -> list[dict]:
    """Return QPX pg record dicts: one per (protein group, unit, label).

    ``intensity`` is an **interim, unnormalized total** — the sum of the group's
    unique-peptide feature intensities for that unit+label (quantms
    ``unique_peptides`` policy, no normalization), until OpenMS ``-out_qpx``
    provides the authoritative protein quant. ``top`` bounds the peptides used
    (``0`` = all; set to 3 to mirror the ProteomicsLFQ/IsobaricWorkflow default).
    Each quantified row is stamped with a ``quantification_method`` cv_param so it
    is never mistaken for the authoritative number. Pass either
    ``consensusxml_path`` (loaded here) or an already-loaded ``cm``.
    """
    cm = cm if cm is not None else load_consensus_map(consensusxml_path)
    map_info = _map_info(cm)
    m = _protein_maps(cm)
    pep_intensity = _peptide_intensities(cm, map_info)
    return build_pg_records(cm, map_info, m, pep_intensity, sdrf_path, top)


def build_pg_records(cm, map_info, m: _ProteinMaps, pep_intensity: dict, sdrf_path, top: int) -> list[dict]:
    """Build pg records from already-accumulated maps + peptide intensities.

    Separated from the accumulation so both the multi-pass adapter and the
    single-pass streaming driver share the exact record-building logic.
    """
    all_runs = sorted({run for run, _ in map_info.values()})
    # grouped_runs unit per run, from the SDRF (fractions grouped); else one unit.
    run_to_grouped = fraction_groups_from_sdrf(sdrf_path)
    units = {tuple(v) for v in run_to_grouped.values()} if run_to_grouped else {tuple(all_runs)}
    # One pg row per (protein group, unit, label): the isobaric channels, or "LFQ".
    labels = sorted({label for _, label in map_info.values()})

    acc_decoy, acc_qvalue, acc_gene, groups = _merge_protein_ids(cm)
    quant_method = "unnormalized_unique_peptide_sum" if not top else f"unnormalized_unique_peptide_top{top}_sum"

    records: list[dict] = []
    for accs in groups:
        anchor = accs[0]
        group_accs = set(accs)
        peptide_seqs: set[str] = set()
        feats: set[tuple] = set()
        for acc in accs:
            peptide_seqs |= m.acc_to_pep.get(acc, set())
            feats |= m.acc_to_feat.get(acc, set())
        # total = every peptide/feature of the group; unique = those mapping only
        # to proteins in this group (same policy as the intensity rollup).
        n_pep_total = len(peptide_seqs)
        n_feat_total = len(feats)
        n_pep_unique = sum(1 for p in peptide_seqs if m.pep_to_accs.get(p, set()).issubset(group_accs))
        n_feat_unique = sum(1 for ft in feats if m.feat_to_accs.get(ft, set()).issubset(group_accs))
        # Prefer the target_decoy meta; fall back to the accession prefix.
        is_decoy = all(acc_decoy.get(a, _is_decoy_accession(a)) for a in accs)
        qvals = [acc_qvalue[a] for a in accs if a in acc_qvalue]
        global_qvalue = min(qvals) if qvals else None
        genes = [acc_gene[a] for a in accs if a in acc_gene] or None
        # Only the quantification units where this group was actually identified
        # (its peptides appear in a run of that unit) — not every unit.
        group_runs: set[str] = set()
        for acc in accs:
            group_runs |= m.acc_to_runs.get(acc, set())
        group_units = [unit for unit in units if group_runs.intersection(unit)] or list(units)
        for unit in group_units:
            for label in labels:
                intensity = _protein_intensity(peptide_seqs, group_accs, unit, label, pep_intensity, m.pep_to_accs, top)
                cv_params = [{"cv_name": "quantification_method", "cv_value": quant_method}] if intensity is not None else None
                records.append(
                    {
                        "pg_accessions": list(accs),
                        "anchor_protein": anchor,
                        "grouped_runs": list(unit),
                        "label": label,
                        # interim unnormalized total; null when the group has no
                        # unique-peptide signal in this unit+label (see cv_params).
                        "intensity": intensity,
                        "global_qvalue": global_qvalue,
                        "is_decoy": is_decoy,
                        "contaminant": any(_is_contaminant(a) for a in accs),
                        "gg_accessions": genes,
                        "gg_names": genes,
                        "peptide_counts": {"unique_sequences": n_pep_unique, "total_sequences": n_pep_total},
                        "feature_counts": {"unique_features": n_feat_unique, "total_features": n_feat_total},
                        "peptides": [{"protein_name": a, "peptide_count": len(m.acc_to_pep.get(a, set()))} for a in accs],
                        "cv_params": cv_params,
                    }
                )
    return records
