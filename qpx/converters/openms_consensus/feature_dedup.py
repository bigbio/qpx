"""Select measured peaks once for Feature rows and protein-group quantification."""

from __future__ import annotations

import logging
from dataclasses import dataclass

from qpx.writers.feature import derive_feature_ids

_log = logging.getLogger(__name__)

# OpenMS producer identity (bigbio/qpx#229), also passed to FeatureWriter.
FEATURE_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name", "rt", "scan", "observed_mz", "consensus_rt")
_PROTEIN_FIELDS = ("anchor_protein", "pg_accessions", "unique", "gg_accessions", "gg_names", "pg_global_qvalue")
_EVIDENCE_FIELDS = ("scan", "id_run_file_name", "pg_positions", "posterior_error_probability", "peptide_qvalue", *_PROTEIN_FIELDS)


def feature_ids(records: list[dict]) -> list[int]:
    """Derive IDs from exactly the values persisted by FeatureWriter."""
    return derive_feature_ids(records, FEATURE_IDENTITY_COMPOSITE)


def _merge_evidence(target: dict, source: dict) -> None:
    """Recover compatible own-run evidence without changing the selected quantity.

    Confidence stays with the selected direct identification, as a PEP/q-value
    pair. A transferred row takes the donor's direct evidence instead. Distinct
    inferred groups are not combined into a guessed protein assignment.
    """
    run = target["run_file_name"]
    if source.get("id_run_file_name") != run:
        return
    if target.get("id_run_file_name") != run:
        target.update({name: source.get(name) for name in _EVIDENCE_FIELDS})
        return
    if any(target.get(name) != source.get(name) for name in _PROTEIN_FIELDS):
        return
    target["scan"] = sorted(set(target.get("scan") or ()) | set(source.get("scan") or ()))
    positions = {
        (pos["protein_accession"], pos["start"], pos["end"]) for rec in (target, source) for pos in rec.get("pg_positions") or ()
    }
    target["pg_positions"] = [
        {"protein_accession": acc, "start": start, "end": end} for acc, start, end in sorted(positions)
    ] or None


@dataclass
class _PeakChoice:
    """Compact winner state; never retain the source XML or pyopenms objects."""

    rank: tuple
    row: int
    fid: int
    record: dict


class FeatureDeduplicator:
    """Keep the feature with more runs, then higher quality (ties keep the first).

    Duplicate rows share an ID or a measured peak. Identification evidence is
    consolidated separately from this quantitative choice. ``row_patches`` and
    ``superseded`` address original output positions, including already flushed
    batches. PSM IDs resolve directly to the final winner, without redirect cycles.
    An optional peptide accumulator receives only the winning peak contributions.
    """

    def __init__(self, pep_intensity=None):
        self._best: dict[tuple, _PeakChoice] = {}
        self.superseded: set[int] = set()
        self.redirect: dict[int, _PeakChoice] = {}
        self.row_patches: dict[int, dict] = {}
        self.dropped = 0
        self.rows = 0
        self.pep_intensity = pep_intensity

    @staticmethod
    def _peak(rec: dict) -> tuple:
        intensities = tuple(sorted((i["label"], i["intensity"]) for i in rec.get("intensities") or ()))
        return (rec["run_file_name"], rec["peptidoform"], rec["charge"], rec["rt"], intensities)

    def _compact(self, rec: dict) -> dict:
        fields = (*FEATURE_IDENTITY_COMPOSITE, *_EVIDENCE_FIELDS)
        if self.pep_intensity is not None:
            fields += ("sequence", "intensities")
        return {name: rec.get(name) for name in fields}

    def _accumulate(self, record: dict, sign: int) -> None:
        if self.pep_intensity is not None:
            for value in record["intensities"]:
                key = (record["sequence"], record["run_file_name"], value["label"])
                self.pep_intensity[key] += sign * value["intensity"]

    def _combine(self, choice: _PeakChoice, rec: dict, fid: int, rank: tuple) -> bool:
        self.dropped += 1
        self.redirect[choice.fid] = self.redirect[fid] = choice
        replace = rank > choice.rank
        if replace:
            self.superseded.add(choice.row)
            self.row_patches.pop(choice.row, None)
            self._accumulate(choice.record, -1)
            merged = self._compact(rec)
            _merge_evidence(merged, choice.record)
            choice.record = merged
            choice.rank, choice.row = rank, self.rows
            self._accumulate(merged, 1)
        else:
            _merge_evidence(choice.record, rec)
        choice.fid = feature_ids([choice.record])[0]
        self._best[("id", choice.fid)] = choice
        self.row_patches[choice.row] = {**choice.record, "feature_id": choice.fid}
        return replace

    def keep(self, records: list[dict], ids: list[int], quality: float) -> list[dict]:
        """Register candidates and return rows to emit; later evidence is patched on close."""
        rank = (len(records), quality)
        kept = []
        for rec, fid in zip(records, ids):
            keys = (("id", fid), ("peak", self._peak(rec)))
            choice = next((self._best[key] for key in keys if key in self._best), None)
            if choice is None:
                choice = _PeakChoice(rank, self.rows, fid, self._compact(rec))
                self._accumulate(choice.record, 1)
                emit = True
            else:
                emit = self._combine(choice, rec, fid, rank)
            for key in keys:
                self._best[key] = choice
            if emit:
                self.rows += 1
                kept.append(rec)
        return kept

    def resolve(self, fid):
        """Final persisted feature ID for an original candidate ID."""
        choice = self.redirect.get(fid)
        return choice.fid if choice is not None else fid

    def log(self) -> None:
        """Warn once about duplicate measured peaks."""
        if self.dropped:
            _log.warning(
                "Dropped %d duplicate feature row(s): one peak in one run reported by two consensus features "
                "(typically isobaric targets linked to one peak); kept the consensus feature with more runs, "
                "then the higher quality, preserving compatible direct identification evidence.",
                self.dropped,
            )
