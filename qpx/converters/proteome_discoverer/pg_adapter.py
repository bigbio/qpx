"""Proteome Discoverer ``Proteins.txt`` to the pg view: one row per file group and channel."""

from __future__ import annotations

import re

import pandas as pd

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.channel_labels import normalize_label
from qpx.converters.mappings import get_field_mappings
from qpx.converters.proteome_discoverer.reports import split_accessions
from qpx.converters.utils import safe_float
from qpx.writers.pg import PgWriter

_PG_MAP = get_field_mappings("proteome_discoverer", "pg")
# "Abundance: F16: 127N, Sample, CORT", "Abundance F16 127N Sample CORT", "Abundances (Grouped): F1, 126" and,
# without a channel, "Abundance: F1: Sample". Grouped by a study factor rather than a file group, a column does not match.
_ABUNDANCE = re.compile(
    r"^(Abundance|Abundances \(?(?:Normalized|Scaled|Grouped)\)?):? (F\d+)(?:[:,]? (\d{3}[NCD]?)\b)?(?:[:,]? \w.*)?$"
)
_GROUPED = "Abundances Grouped"
_ADDITIONAL_NAMES = {"Abundances Normalized": "normalized_abundance", "Abundances Scaled": "scaled_abundance"}


def read_proteins(path: str) -> pd.DataFrame:
    """The ``Proteins.txt`` export, with all values as text."""
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)


def master_groups(proteins: pd.DataFrame) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Members of each master protein's group (itself first) and its gene, as far as the export says."""
    resolved = resolve_columns(_PG_MAP, set(proteins.columns))
    accession, group_ids, gene = (resolved.get(field) for field in ("anchor_protein", "protein_group_ids", "gg_names"))
    masters = proteins[proteins["Master"] == "IsMasterProtein"]
    members_by_group: dict[str, list[str]] = {}
    if group_ids:
        for name, ids in zip(proteins[accession], proteins[group_ids].fillna("")):
            for group in split_accessions(ids):
                members_by_group.setdefault(group, []).append(name)
    groups, genes = {}, {}
    for _, row in masters.iterrows():
        master = row[accession]
        others = [name for group in split_accessions(row[group_ids]) for name in members_by_group[group]] if group_ids else []
        groups[master] = list(dict.fromkeys([master, *others]))
        if gene and isinstance(row[gene], str):
            genes[master] = row[gene]
    return groups, genes


def _abundance_columns(columns) -> dict[str, tuple[str, str, str | None]]:
    """Column -> (kind, file group, canonical label or None for a label-free column)."""
    found = {}
    for column in columns:
        match = _ABUNDANCE.match(column)
        if match:
            kind = match.group(1).replace("(", "").replace(")", "")
            found[column] = (kind, match.group(2), normalize_label(match.group(3)) if match.group(3) else None)
    return found


class ProteomeDiscovererPgAdapter(BaseConverter):
    """Write the master proteins, one row per quantified file group and channel."""

    _resolved: dict[str, str] = {}
    _groups: dict[str, list[str]] = {}
    _file_groups: dict[str, list[str]] = {}

    def convert(self, proteins: pd.DataFrame, output_path: str, runs: dict[str, str]) -> None:
        self._resolved = resolve_columns(_PG_MAP, set(proteins.columns))
        self._groups, _ = master_groups(proteins)
        columns = _abundance_columns(proteins.columns)
        # PD writes the per-file abundance only on request; the grouped one is then all a file group has.
        primary = {column: spec for column, spec in columns.items() if spec[0] == "Abundance"} or {
            column: spec for column, spec in columns.items() if spec[0] == _GROUPED
        }
        additional = {
            (spec[1], spec[2]): (column, _ADDITIONAL_NAMES[spec[0]])
            for column, spec in columns.items()
            if spec[0] in _ADDITIONAL_NAMES
        }
        self._file_groups = _runs_by_file_group(runs)
        unknown = sorted({spec[1] for spec in primary.values()} - set(self._file_groups))
        if unknown:
            raise ValueError(
                f"No raw runs for Proteome Discoverer file group(s) {', '.join(unknown)}; pass the InputFiles export"
            )
        masters = proteins[proteins["Master"] == "IsMasterProtein"].to_dict("records")
        with PgWriter(output_path, creator="proteome_discoverer", compression=self._compression) as writer:
            records = self._transform_rows_with_skip(
                masters, lambda row: self._records(row, primary, additional), "Proteome Discoverer PG"
            )
            if records:
                self._track_scores(records)
                writer.write_batch(records)

    def _records(self, row: dict, primary: dict, additional: dict) -> list[dict]:
        common = self._common_fields(row)
        if not primary:
            # An identification-only export: the group is known, its quantity is not.
            runs = sorted(run for group in self._file_groups.values() for run in group)
            return [
                {
                    **common,
                    "grouped_runs": runs,
                    "intensities": [{"label": "LFQ", "intensity": None}],
                    "additional_intensities": None,
                }
            ]
        records = []
        for column, (_, file_group, channel) in primary.items():
            value = safe_float(row.get(column))
            if value is None or value <= 0:
                continue
            label = channel or "LFQ"
            records.append(
                {
                    **common,
                    "grouped_runs": self._file_groups[file_group],
                    "intensities": [{"label": label, "intensity": value}],
                    "additional_intensities": _additional_abundances(row, additional, file_group, channel, label),
                }
            )
        return records

    def _common_fields(self, row: dict) -> dict:
        r = self._resolved
        master = row[r["anchor_protein"]]
        members = self._groups.get(master, [master])
        peptides = _count(row, r, "peptide_count_total")
        unique = _count(row, r, "peptide_count_unique")
        gene = row.get(r.get("gg_names"))
        description = row.get(r.get("pg_names"))
        length = _count(row, r, "sequence_length")
        scores = [
            {"score_name": "sum_pep_score", "score_value": value, "higher_better": True}
            for value in [safe_float(row.get(r.get("sum_pep_score")))]
            if value is not None
        ]
        return {
            "pg_accessions": members,
            "pg_names": [description] if isinstance(description, str) else None,
            "gg_accessions": [gene] if isinstance(gene, str) else None,
            "gg_names": [gene] if isinstance(gene, str) else None,
            "gg_qvalue": None,
            "anchor_protein": master,
            "global_qvalue": None,
            "pg_qvalue": safe_float(row.get(r.get("pg_qvalue"))),
            "is_decoy": False,
            "contaminant": row.get(r.get("contaminant")) == "True",
            "peptides": [{"protein_name": master, "peptide_count": peptides}] if peptides is not None else None,
            "peptide_counts": (
                None if peptides is None or unique is None else {"unique_sequences": unique, "total_sequences": peptides}
            ),
            "feature_counts": None,
            "sequence_coverage": safe_float(row.get(r.get("sequence_coverage"))),
            "molecular_weight": safe_float(row.get(r.get("molecular_weight"))),
            "cv_params": [{"cv_name": "sequence_length", "cv_value": str(length)}] if length is not None else None,
            "additional_scores": scores or None,
        }


def _additional_abundances(row: dict, additional: dict, file_group: str, channel: str | None, label: str) -> list[dict] | None:
    column_and_name = additional.get((file_group, channel))
    if column_and_name is None:
        return None
    column, name = column_and_name
    value = safe_float(row.get(column))
    if value is None:
        return None
    return [{"label": label, "intensities": [{"intensity_name": name, "intensity_value": value}]}]


def _count(row: dict, resolved: dict, field: str) -> int | None:
    value = safe_float(row.get(resolved.get(field)))
    return None if value is None else int(value)


def _runs_by_file_group(runs: dict[str, str]) -> dict[str, list[str]]:
    """``F16`` -> the runs of its fractions ``F16.1``, ``F16.2``, in fraction order."""
    groups: dict[str, list[str]] = {}
    for file_id, run in sorted(runs.items(), key=lambda item: int(item[0].partition(".")[2] or 0)):
        groups.setdefault(file_id.split(".")[0], []).append(run)
    return groups
