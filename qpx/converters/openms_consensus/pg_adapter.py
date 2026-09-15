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
from qpx.converters.openms_consensus.protein_groups import ProteinGroupIndex, identification_identifier
from qpx.converters.openms_consensus.psm_adapter import _run_resolver
from qpx.converters.utils import is_contaminant_accession, safe_float, uniprot_entry_name
from qpx.core.protein_sequence import average_molecular_weight_kda

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
    return is_contaminant_accession(acc)


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
        score = safe_float(hit.getScore())
        if score_is_qvalue and score is not None:
            acc_qvalue[acc] = score
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


def _protein_molecular_weight(sequences: set[str]) -> float | None:
    """Average mass in kDa for one agreed, unmodified protein sequence."""
    if len(sequences) != 1:
        return None
    # Shared with the FASTA-based protein-properties transform, so both agree.
    return average_molecular_weight_kda(next(iter(sequences)))


def _anchor_properties(coverages: set[float], probabilities: set[float], sequences: set[str]) -> dict:
    """Keep an anchor's known properties only when its source records agree."""
    coverage = next(iter(coverages)) if len(coverages) == 1 else None
    scores = None
    if len(probabilities) == 1:
        scores = [{"score_name": "posterior_probability", "score_value": next(iter(probabilities)), "higher_better": True}]
    return {
        "sequence_coverage": coverage,
        "additional_scores": scores,
        "molecular_weight": _protein_molecular_weight(sequences),
    }


def _protein_properties(cm) -> dict[str, dict]:
    """Index coverage, posterior probability and theoretical mass by accession.

    ProteinHit properties describe the representative protein, not a group-wide
    aggregate. Mass is computed from the complete unmodified ProteinHit sequence.
    Missing, unknown or conflicting properties do not supply a value; other group
    members are never used as a substitute.
    """
    coverages: dict[str, set[float]] = defaultdict(set)
    probabilities: dict[str, set[float]] = defaultdict(set)
    sequences: dict[str, set[str]] = defaultdict(set)
    for prot in cm.getProteinIdentifications():
        for hit in prot.getHits():
            acc = _acc_str(hit.getAccession())
            sequence = hit.getSequence()
            if sequence:
                sequences[acc].add(sequence)
            coverage = hit.getCoverage()
            if 0 <= coverage <= 100:
                coverages[acc].add(float(coverage))
            if hit.metaValueExists("Posterior Probability_score"):
                probability = safe_float(hit.getMetaValue("Posterior Probability_score"))
                if probability is not None and 0 <= probability <= 1:
                    probabilities[acc].add(probability)
    return {
        acc: _anchor_properties(coverages[acc], probabilities[acc], sequences[acc])
        for acc in coverages.keys() | probabilities.keys() | sequences.keys()
    }


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


def _group_qvalue_and_genes(
    accs, acc_qvalue: dict[str, float], acc_gene: dict[str, str]
) -> tuple[float | None, list[str] | None]:
    """A protein group's global q-value (best member) and gene names.

    The single definition behind both ``pg.global_qvalue``/``pg.gg_names`` and the
    feature view's ``pg_global_qvalue``/``gg_names``, so the two views cannot
    disagree about the same group.
    """
    qvals = [acc_qvalue[a] for a in accs if a in acc_qvalue]
    genes = [acc_gene[a] for a in accs if a in acc_gene] or None
    return (min(qvals) if qvals else None), genes


def protein_group_maps(cm) -> tuple[ProteinGroupIndex, dict[tuple[str, ...], tuple[float | None, list[str] | None]]]:
    """Group lookup by membership/source, plus group -> (global q-value, genes).

    Both come from one ``_merge_protein_ids`` pass. The second map lets the
    feature view carry the group's confidence and gene names: they were computed
    here for pg and then discarded, leaving ``feature.pg_global_qvalue`` and
    ``feature.gg_names`` null on every OpenMS dataset while pg had them for every
    group (the DIA-NN converter fills both).
    """
    _, acc_qvalue, acc_gene, groups = _merge_protein_ids(cm)
    group_map = ProteinGroupIndex.from_groups(groups)
    groups_by_identification: dict[str, list[list[str]]] = defaultdict(list)
    for prot in cm.getProteinIdentifications():
        identifier = identification_identifier(prot)
        if identifier:
            groups_by_identification[identifier].extend(_build_groups(prot))
    group_map.by_identification = {
        identifier: ProteinGroupIndex.from_groups(source_groups) for identifier, source_groups in groups_by_identification.items()
    }
    group_meta = {tuple(grp): _group_qvalue_and_genes(grp, acc_qvalue, acc_gene) for grp in groups}
    return group_map, group_meta


def accession_to_group(cm) -> dict[str, list[str]]:
    """Map unambiguous accessions to their full protein-group membership.

    group[0] remains the producer's leader. For shared accessions and groups from
    multiple identification runs, use the full index in ``protein_group_maps``.
    """
    group_map, _ = protein_group_maps(cm)
    return group_map.unambiguous_accessions()


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


def pg_units_and_labels(map_info, sdrf_path) -> tuple[set[tuple[str, ...]], list[str]]:
    """Return the ``(grouped_runs units, labels)`` a pg build spans for this map.

    Used by :func:`build_pg_records` (which emits one row per (group, unit,
    label)) to decide how runs are grouped into quantification units and which
    labels exist. A unit is the set of raw files aggregated together (SDRF
    fractions grouped; else all runs as one unit); labels are the isobaric
    channels or ``LFQ``.
    """
    all_runs = sorted({run for run, _ in map_info.values()})
    # grouped_runs unit per run, from the SDRF (fractions grouped); else one unit.
    run_to_grouped = fraction_groups_from_sdrf(sdrf_path)
    units = {tuple(v) for v in run_to_grouped.values()} if run_to_grouped else {tuple(all_runs)}
    labels = sorted({label for _, label in map_info.values()})
    return units, labels


def build_pg_records(cm, map_info, m: _ProteinMaps, pep_intensity: dict, sdrf_path, top: int) -> list[dict]:
    """Build pg records from already-accumulated maps + peptide intensities.

    Separated from the accumulation so both the multi-pass adapter and the
    single-pass streaming driver share the exact record-building logic.
    """
    # One pg row per (protein group, unit, label): the isobaric channels, or "LFQ".
    units, labels = pg_units_and_labels(map_info, sdrf_path)

    acc_decoy, acc_qvalue, acc_gene, groups = _merge_protein_ids(cm)
    properties = _protein_properties(cm)
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
        global_qvalue, genes = _group_qvalue_and_genes(accs, acc_qvalue, acc_gene)
        # Entry names from the ``db|ACC|NAME`` accessions, aligned with pg_accessions;
        # null unless every member has one, so a partial list never misaligns.
        member_names = [uniprot_entry_name(a) for a in accs]
        pg_names = member_names if all(member_names) else None
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
                        "pg_names": pg_names,
                        "anchor_protein": anchor,
                        **properties.get(anchor, {}),
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
