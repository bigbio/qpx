"""consensusXML -> QPX pg (protein group) records.

Protein groups come from the OpenMS protein-inference graph: indistinguishable
protein groups when present, else each protein hit as a singleton group. Peptide
counts come from the peptide→protein evidence links.

**Protein intensity**: the abundance OpenMS recorded for the group
(``fraction_group_level_abundance``, see :mod:`.protein_quant`) when the input has
any. Files without one (quantms feeds ProteinQuantifier separately) get a roll-up
instead: the **unnormalized sum of the group's unique peptides** (quantms
``unique_peptides`` policy, no normalization), per ``(protein group, grouped_runs
unit, label)``, optionally over the ``top`` most abundant peptides. A
``quantification_method`` cv_param names the source of every quantified row, and
the roll-up stays available as an additional intensity beside OpenMS' value. Rows
stay ``intensity``-null where a group has no signal. One row per ``(protein group,
grouped_runs, label)``; a quantification unit comes from the SDRF, else from the
``fraction_group`` of the map list, else it is one run.
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from typing import Optional

from qpx.converters.channel_labels import fraction_groups_from_sdrf
from qpx.converters.openms_consensus.feature_adapter import is_decoy_hit, to_proforma
from qpx.converters.openms_consensus.protein_groups import ProteinGroupIndex, identification_identifier
from qpx.converters.openms_consensus.protein_quant import fraction_group_units, quant_by_unit
from qpx.converters.utils import is_contaminant_accession, safe_float, uniprot_entry_name
from qpx.core.protein_sequence import average_molecular_weight_kda

_log = logging.getLogger(__name__)
_GENE_RE = re.compile(r"GN=([^\s]+)")
_OPENMS_QUANT_METHOD = "openms_protein_abundance"


class _ProteinMaps:
    """Peptide/feature/run evidence indexed by accession, with the inverse maps.

    ``acc_to_*`` power the per-group unions; ``pep_to_accs`` / ``feat_to_accs`` are
    the inverse (peptide/feature -> the accessions it maps to) used to decide which
    peptides/features are *unique* to a group (map only to proteins in it) — the
    same ``unique_peptides`` policy the intensity rollup uses. ``acc_to_runs``
    holds every run that can carry a row (quantified or identified);
    ``acc_to_psm_runs`` only the runs of the accession's own PSMs, its
    identification evidence.
    """

    def __init__(self):
        self.acc_to_pep: dict[str, set[str]] = defaultdict(set)
        self.acc_to_runs: dict[str, set[str]] = defaultdict(set)
        self.acc_to_psm_runs: dict[str, set[str]] = defaultdict(set)
        self.acc_to_feat: dict[str, set[tuple]] = defaultdict(set)
        self.pep_to_accs: dict[str, set[str]] = defaultdict(set)
        self.feat_to_accs: dict[tuple, set[str]] = defaultdict(set)


class ProteinGroupAccumulator:
    """Protein-group evidence and peptide intensities accumulated over every input map.

    Also serves as the protein-identification source for ``build_pg_records``, so
    groups, decoy flags, q-values and properties span all inputs.
    """

    def __init__(self):
        self.maps = _ProteinMaps()
        self.pep_intensity: dict = defaultdict(float)
        self.map_info: dict = {}
        self.fraction_units: set[tuple[str, ...]] = set()
        self.protein_quant: dict = {}
        self._protein_ids: list = []
        self._inputs = 0

    def add_source(self, cm, map_info, path: str) -> None:
        """Register one input's protein identifications, map columns, units and abundances."""
        self.map_info.update({(self._inputs, idx): info for idx, info in map_info.items()})
        self._protein_ids.extend(cm.getProteinIdentifications())
        self._inputs += 1
        runs = {idx: run for idx, (run, _label) in map_info.items()}
        labels = {idx: label for idx, (_run, label) in map_info.items()}
        self.fraction_units.update(fraction_group_units(cm, runs).values())
        self.protein_quant.update(quant_by_unit(cm, runs, labels, path))

    def getProteinIdentifications(self) -> list:  # pylint: disable=invalid-name
        """ProteinIdentifications of every input (pyopenms accessor name)."""
        return self._protein_ids

    def build(self, sdrf_path, top) -> list[dict]:
        """pg records over all inputs."""
        return build_pg_records(
            self, self.map_info, self.maps, self.pep_intensity, sdrf_path, top, self.fraction_units, self.protein_quant
        )


def _evidence_accession(ev) -> str | None:
    acc = ev.getProteinAccession()
    if not acc:
        return None
    return acc.decode() if isinstance(acc, bytes) else str(acc)


def _collect_pid(pid, runs: set[str], psm_run: str | None, m: _ProteinMaps) -> None:
    """Fold one PeptideIdentification's evidence into the accession maps.

    ``runs`` are the runs its consensus feature quantified; ``psm_run`` is the
    PSM's own run (``None`` when it cannot be resolved).
    """
    psm_runs = {psm_run} if psm_run else set()
    for hit in pid.getHits():
        seq_obj = hit.getSequence()
        seq = seq_obj.toUnmodifiedString()
        feat = (to_proforma(seq_obj), int(hit.getCharge() or 0))
        for ev in hit.getPeptideEvidences():
            acc = _evidence_accession(ev)
            if acc:
                m.acc_to_pep[acc].add(seq)
                m.acc_to_feat[acc].add(feat)
                m.acc_to_runs[acc] |= runs | psm_runs
                m.acc_to_psm_runs[acc] |= psm_runs
                m.pep_to_accs[seq].add(acc)
                m.feat_to_accs[feat].add(acc)


def _cf_runs(cf, map_run: dict[int, str]) -> set[str]:
    runs = {map_run.get(sub.getMapIndex()) for sub in cf.getFeatureList() if sub.getIntensity() > 0}
    runs.discard(None)
    return runs


def accumulate_cf_maps(cf, map_run: dict[int, str], resolve_run, m: _ProteinMaps) -> None:
    """Fold one consensus feature's assigned IDs into the accession maps (per-cf)."""
    runs = _cf_runs(cf, map_run)
    for pid in cf.getPeptideIdentifications():
        # The PSM's own run, resolved as the PSM view does.
        _collect_pid(pid, runs, resolve_run(pid, cf_runs=runs), m)


def accumulate_unassigned_maps(pid, resolve_run, m: _ProteinMaps) -> None:
    """Fold one unassigned PeptideIdentification into the accession maps."""
    _collect_pid(pid, set(), resolve_run(pid), m)


def accumulate_consensus_map(cm, map_info, resolve_run, m: _ProteinMaps) -> None:
    """Fold all identification evidence into the shared maps, including duplicate peaks."""
    map_run = {idx: run for idx, (run, _label) in map_info.items()}
    for cf in cm:
        accumulate_cf_maps(cf, map_run, resolve_run, m)
    for pid in cm.getUnassignedPeptideIdentifications():
        accumulate_unassigned_maps(pid, resolve_run, m)


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
            acc_decoy[acc] = is_decoy_hit(hit)
        score = safe_float(hit.getScore())
        if score_is_qvalue and score is not None:
            acc_qvalue[acc] = score
        gene = _GENE_RE.search(str(hit.getDescription() or ""))
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


def _protein_intensity(group_peps, group_accs, unit, label, pep_intensity, pep_accs, top) -> Optional[float]:
    """Unnormalized protein-group roll-up for one ``(unit, label)``.

    Sum of the group's **unique** peptides — those mapping only to proteins in the
    group (the quantms ``unique_peptides`` policy) — with **no normalization**.
    ``top > 0`` keeps only the N most-abundant peptides; ``top = 0`` sums all. Returns ``None`` when the group has no
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


def pg_units_and_labels(map_info, sdrf_path, fraction_units=()) -> tuple[set[tuple[str, ...]], list[str]]:
    """Return the ``(grouped_runs units, labels)`` a pg build spans for this map.

    Used by :func:`build_pg_records` (which emits one row per (group, unit,
    label)) to decide how runs are grouped into quantification units and which
    labels exist. A unit is the set of raw files aggregated together: from the
    SDRF, else the map list's ``fraction_group``, else each run on its own.
    """
    labels = sorted({label for _, label in map_info.values()})
    run_to_grouped = fraction_groups_from_sdrf(sdrf_path)
    if run_to_grouped:
        return {tuple(v) for v in run_to_grouped.values()}, labels
    all_runs = sorted({run for run, _ in map_info.values()})
    covered = {run for unit in fraction_units for run in unit}
    if len(all_runs) > 1 and covered != set(all_runs):
        _log.warning(
            "No SDRF and no complete fraction_group in the map list: each of the %d runs becomes its own "
            "pg quantification unit. Pass an SDRF if runs are fractions of one sample.",
            len(all_runs),
        )
    return set(fraction_units) | {(run,) for run in all_runs if run not in covered}, labels


def _group_sequence_counts(accs, m: _ProteinMaps) -> tuple[set, set, dict]:
    """Count all and group-unique peptide sequences and features."""
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
    return (
        peptide_seqs,
        group_accs,
        {
            "peptide_counts": {"unique_sequences": n_pep_unique, "total_sequences": n_pep_total},
            "feature_counts": {"unique_features": n_feat_unique, "total_features": n_feat_total},
        },
    )


def _pg_group_fields(accs, m, acc_decoy, acc_qvalue, acc_gene, properties) -> dict:
    """Build protein identity, properties and confidence shared by all units."""
    anchor = accs[0]
    # Prefer the target_decoy meta; fall back to the accession prefix.
    is_decoy = all(acc_decoy.get(a, _is_decoy_accession(a)) for a in accs)
    global_qvalue, genes = _group_qvalue_and_genes(accs, acc_qvalue, acc_gene)
    # Entry names from the ``db|ACC|NAME`` accessions, aligned with pg_accessions;
    # null unless every member has one, so a partial list never misaligns.
    member_names = [uniprot_entry_name(a) for a in accs]
    pg_names = member_names if all(member_names) else None
    return {
        "pg_accessions": list(accs),
        "pg_names": pg_names,
        "anchor_protein": anchor,
        **properties.get(anchor, {}),
        "global_qvalue": global_qvalue,
        "is_decoy": is_decoy,
        "contaminant": any(_is_contaminant(a) for a in accs),
        "gg_accessions": genes,
        "gg_names": genes,
        "peptides": [{"protein_name": a, "peptide_count": len(m.acc_to_pep.get(a, set()))} for a in accs],
    }


def _quant_fields(openms, rollup, method: str, label: str, reported: bool) -> dict:
    """Intensity of one pg row: OpenMS' abundance in an input that reports them, else the roll-up."""
    if not reported:
        return {
            "intensity": rollup,
            "cv_params": [{"cv_name": "quantification_method", "cv_value": method}] if rollup is not None else None,
        }
    openms = openms or {}
    extra = [{"intensity_name": method, "intensity_value": rollup}] if rollup is not None else []
    extra += [{"intensity_name": name, "intensity_value": value} for name, value in openms.get("ratios", {}).items()]
    abundance = openms.get("abundance")
    return {
        "intensity": abundance,
        "cv_params": [{"cv_name": "quantification_method", "cv_value": _OPENMS_QUANT_METHOD}] if abundance is not None else None,
        "additional_intensities": [{"label": label, "intensities": extra}] if extra else None,
    }


def build_pg_records(
    cm, map_info, m: _ProteinMaps, pep_intensity: dict, sdrf_path, top: int, fraction_units=(), protein_quant=None
) -> list[dict]:
    """Build pg records from already-accumulated maps + peptide intensities.

    Separated from the accumulation so both the multi-pass adapter and the
    streaming driver share the exact record-building logic.
    """
    # One pg row per (protein group, unit, label): the isobaric or MS1 channels, or "LFQ".
    units, labels = pg_units_and_labels(map_info, sdrf_path, fraction_units)
    protein_quant = protein_quant or {}
    reported_units = {key[1] for key in protein_quant}
    acc_decoy, acc_qvalue, acc_gene, groups = _merge_protein_ids(cm)
    properties = _protein_properties(cm)
    quant_method = "unnormalized_unique_peptide_sum" if not top else f"unnormalized_unique_peptide_top{top}_sum"
    records: list[dict] = []
    for accs in groups:
        peptide_seqs, group_accs, counts = _group_sequence_counts(accs, m)
        group_fields = _pg_group_fields(accs, m, acc_decoy, acc_qvalue, acc_gene, properties)
        # Only the quantification units where this group was identified or
        # quantified — not every unit. A group with no attributable run gets no row.
        group_runs: set[str] = set()
        psm_runs: set[str] = set()
        for acc in accs:
            group_runs |= m.acc_to_runs.get(acc, set())
            psm_runs |= m.acc_to_psm_runs.get(acc, set())
        group_units = [unit for unit in units if group_runs.intersection(unit)]
        for unit in group_units:
            for label in labels:
                rollup = _protein_intensity(peptide_seqs, group_accs, unit, label, pep_intensity, m.pep_to_accs, top)
                unit_key = tuple(sorted(unit))
                openms = protein_quant.get((frozenset(accs), unit_key, label))
                # Without a PSM of the group in this unit, only a measured quantity
                # (e.g. a transferred feature) is evidence: no null row.
                if rollup is None and not (openms and openms["abundance"]) and not psm_runs.intersection(unit):
                    continue
                records.append(
                    {
                        **group_fields,
                        "grouped_runs": list(unit),
                        "label": label,
                        **counts,
                        **_quant_fields(openms, rollup, quant_method, label, unit_key in reported_units),
                    }
                )
    return records
