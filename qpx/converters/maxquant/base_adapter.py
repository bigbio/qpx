"""Shared base adapter for MaxQuant converters.

Provides common data-loading helpers used by
:class:`MaxQuantFeatureAdapter` and :class:`MaxQuantPgAdapter`.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.channel_labels import label_unimod, normalize_label, plex_for_columns, read_sdrf_labels
from qpx.converters.maxquant.constants import FixedMods, parameters_fixed_mods, sdrf_fixed_mods, summary_fixed_mods
from qpx.converters.maxquant.protein_fields import parse_genes
from qpx.converters.utils import canonical_accession, safe_float
from qpx.core.files import run_file_stem

_REPORTER_COLUMN = re.compile(r"Reporter intensity (\d+)(?: .+)?")


class MaxQuantBaseAdapter(BaseConverter):
    """Base adapter with MaxQuant-specific loading utilities.

    Subclasses inherit ``_load_sdrf()`` so the SDRF loading logic is
    defined once rather than duplicated across the feature and PG adapters.
    """

    # Fixed modifications written into peptidoforms; set by ``convert``.
    _fixed: FixedMods = ()
    # SILAC multiplet letter -> label, lightest first; set by ``convert``.
    _multiplet: dict[str, str] = {}

    @staticmethod
    def _norm_uniprot_id(pid: str) -> str:
        """Normalise a MaxQuant UniProt protein ID to its bare accession.

        Thin wrapper over :func:`qpx.converters.utils.canonical_accession`::

            sp|P55011|S12A2_HUMAN  →  P55011
            tr|A0A3B3IS91|..._HUMAN  →  A0A3B3IS91
            P55011  →  P55011  (returned as-is)
        """
        return canonical_accession(pid)

    @classmethod
    def _protein_list(cls, cell) -> list[str]:
        """Normalised accessions of a ``;``-separated MaxQuant protein cell.

        This is the ``pg_accessions`` list the pg view writes for a
        ``proteinGroups.txt`` row, so the feature view reuses it to carry the
        same membership.
        """
        if cell is None or pd.isna(cell):
            return []
        return [cls._norm_uniprot_id(a.strip()) for a in str(cell).split(";") if a.strip()]

    @classmethod
    def _group_anchor(cls, majority_cell, accessions: list[str]) -> str:
        """The group's representative: first ``Majority protein IDs`` entry, else first member."""
        if majority_cell is not None and pd.notna(majority_cell) and majority_cell:
            return cls._norm_uniprot_id(str(majority_cell).split(";", maxsplit=1)[0])
        return accessions[0] if accessions else ""

    def _build_protein_group_maps(self, protein_groups_path: Optional[str]) -> dict:
        """Build protein-group lookup maps from proteinGroups.txt.

        Every map value is ``{"accessions", "anchor", "qvalue", "genes"}``,
        where ``accessions`` and ``anchor`` are exactly what the pg view writes
        for the group. Returns a dict with three sub-maps, and ``columns``,
        the proteinGroups.txt columns ``accessions`` and ``anchor`` come from:
            ``by_id``: proteinGroups ``id`` -> value. evidence.txt names each
                row's groups in ``Protein group IDs``, so this is the primary,
                exact key.
            ``by_group``: full protein-group membership (a ``frozenset`` of
                normalised accessions) -> value. Used, like ``by_accession``,
                only for rows without usable ``Protein group IDs``.
            ``by_accession``: single accession -> value, but *only* for
                accessions that belong to exactly one protein group.
                Ambiguous razor accessions (shared across groups) are excluded
                so a feature cannot inherit another group's q-value/genes via
                its leader.
        """
        empty = {"by_id": {}, "by_group": {}, "by_accession": {}, "columns": {}}
        if not protein_groups_path:
            return empty
        df = self._conn.execute(
            "SELECT * FROM read_csv_auto($1, delim='\t', header=true, auto_detect=true, null_padding=true)",
            [protein_groups_path],
        ).df()
        acc_col, qval_col, gene_col = self._detect_pg_columns(df)
        if not acc_col:
            raise ValueError(f"{protein_groups_path} has no 'Protein IDs' or 'Majority protein IDs' column")
        maps = self._build_pg_lookup(df, acc_col, qval_col, gene_col)
        anchor_col = "Majority protein IDs" if "Majority protein IDs" in df.columns else acc_col
        maps["columns"] = {"pg_accessions": acc_col, "anchor_protein": anchor_col}
        self.logger.info(
            "Built protein group maps: %d group(s), %d group id(s), %d unambiguous accession(s)",
            len(maps["by_group"]),
            len(maps["by_id"]),
            len(maps["by_accession"]),
        )
        return maps

    @staticmethod
    def _detect_pg_columns(df: pd.DataFrame) -> tuple:
        """Detect accession, q-value, and gene name columns in proteinGroups."""
        acc_col = next(
            (c for c in ["Protein IDs", "Majority protein IDs"] if c in df.columns),
            None,
        )
        qval_col = next((c for c in ["Q-value", "q-value", "Q.value"] if c in df.columns), None)
        gene_col = next((c for c in ["Gene names", "Gene Names"] if c in df.columns), None)
        return acc_col, qval_col, gene_col

    @staticmethod
    def _canonical_group_key(accessions) -> frozenset:
        """Order-independent key for a protein-group membership."""
        return frozenset(a.strip() for a in accessions if a and a.strip())

    def _build_pg_lookup(self, df: pd.DataFrame, acc_col: str, qval_col: Optional[str], gene_col: Optional[str]) -> dict:
        """Build the ``by_id``/``by_group``/``by_accession`` lookup maps (see caller)."""
        by_id: dict[int, dict] = {}
        by_group: dict[frozenset, dict] = {}
        acc_to_groups: dict[str, set] = {}
        acc_value: dict[str, dict] = {}

        def values(column: Optional[str]) -> list:
            return df[column].tolist() if column and column in df.columns else [None] * len(df)

        qval_vals = pd.to_numeric(df[qval_col], errors="coerce").tolist() if qval_col else [None] * len(df)
        columns = zip(
            values(acc_col), qval_vals, values(gene_col), values("Fasta headers"), values("Majority protein IDs"), values("id")
        )
        for acc_raw, qval, gene_raw, fasta_raw, majority_raw, group_id in columns:
            accs = self._protein_list(acc_raw)
            if not accs:
                continue
            value = {
                "accessions": accs,
                "anchor": self._group_anchor(majority_raw, accs),
                "qvalue": safe_float(qval),
                "genes": parse_genes(gene_raw, fasta_raw),
            }
            if pd.notna(group_id):
                by_id[int(group_id)] = value
            key = self._canonical_group_key(accs)
            by_group.setdefault(key, value)
            for acc in value["accessions"]:
                acc_to_groups.setdefault(acc, set()).add(key)
                acc_value.setdefault(acc, value)

        # Keep only accessions that map to a single protein group for the
        # leader fallback; ambiguous accessions would otherwise collide.
        by_accession = {acc: acc_value[acc] for acc, keys in acc_to_groups.items() if len(keys) == 1}
        return {"by_id": by_id, "by_group": by_group, "by_accession": by_accession}

    def _lookup_pg_metadata(self, pg_maps: Optional[dict], pg_acc_list: list[str], anchor_protein: str) -> tuple:
        """Resolve ``(pg_global_qvalue, genes)`` for a feature's protein group.

        Keyed on the full protein-group membership first; only falls back to the
        leader accession when that leader is unambiguous (belongs to one group).
        """
        if not pg_maps:
            return None, None
        entry = pg_maps.get("by_group", {}).get(self._canonical_group_key(pg_acc_list))
        if entry is None and anchor_protein:
            entry = pg_maps.get("by_accession", {}).get(anchor_protein)
        if entry is None:
            return None, None
        return entry.get("qvalue"), entry.get("genes")

    def _fixed_mods(self, sdrf_path: Optional[str], source_path: str) -> FixedMods:
        """Fixed modifications as ``(accession, residues, terminus)``, as :func:`to_proforma` takes them.

        MaxQuant leaves fixed modifications and isobaric labels out of
        ``Modified sequence`` while ``Mass`` and ``m/z`` include the fixed ones.
        They are read from the SDRF ``MT=fixed`` entries, else from the
        ``parameters.txt`` or ``summary.txt`` MaxQuant writes next to
        *source_path*; an isobaric label none lists is added from the SDRF
        channel labels.
        """
        fixed: FixedMods = ()
        if sdrf_path:
            from qpx.core.sdrf import SDRFHandler

            fixed = sdrf_fixed_mods(SDRFHandler(sdrf_path).get_modification_parameters())
        folder = Path(source_path)
        fixed = (
            fixed
            or parameters_fixed_mods(folder.with_name("parameters.txt"))
            or summary_fixed_mods(folder.with_name("summary.txt"))
        )
        if not fixed:
            self.logger.warning(
                "Neither the SDRF nor MaxQuant's parameters.txt or summary.txt declares a fixed modification that can "
                "be placed, so peptidoforms carry none; MaxQuant leaves fixed modifications out of 'Modified sequence' "
                "although calculated_mz includes them"
            )
        label = label_unimod(read_sdrf_labels(sdrf_path))
        if label is not None and all(accession != f"UNIMOD:{label}" for accession, _, _ in fixed):
            fixed += ((f"UNIMOD:{label}", "", "N"), (f"UNIMOD:{label}", "K", ""))
        return fixed

    @staticmethod
    def _load_sdrf(
        sdrf_path: Optional[str], columns: Iterable[str] = ()
    ) -> tuple[dict[str, set[str]], str, list[tuple[str, int]]]:
        """Load SDRF and return ``(run_labels, experiment_type, channels)``.

        When *sdrf_path* is ``None`` the method returns LFQ defaults.

        Args:
            sdrf_path: Filesystem path to the SDRF file, or ``None``.
            columns: Columns of the MaxQuant table; its reporter columns size the plex.

        Returns:
            A 3-tuple of:
              - *run_labels*: ``{run file stem -> canonical labels the SDRF assigns a sample}``
              - *experiment_type*: e.g. ``"LFQ"``, ``"TMT11"``
              - *channels*: ``(canonical label, 0-based reporter index)`` per declared
                isobaric channel, in plex order (empty for LFQ)
        """
        if not sdrf_path:
            return {}, "LFQ", []

        from qpx.core.sdrf import SDRFHandler

        handler = SDRFHandler(sdrf_path)
        experiment_type = handler.get_experiment_type_from_sdrf()
        table = handler.sdrf_table
        run_labels: dict[str, set[str]] = {}
        for data_file, label in zip(table["comment[data file]"], table["comment[label]"]):
            if pd.notna(data_file) and pd.notna(label):
                run_labels.setdefault(run_file_stem(str(data_file)), set()).add(normalize_label(str(label)))

        channels: list[tuple[str, int]] = []
        family = {"TMT": "TMT", "ITRAQ": "iTRAQ"}.get(re.sub(r"\d+", "", experiment_type))
        if family:
            labels = set().union(*run_labels.values())
            # MaxQuant numbers reporter columns by plex position, N-isomer before C-isomer.
            n_columns = len({match.group(1) for column in columns if (match := _REPORTER_COLUMN.fullmatch(column))})
            index = {label: position - 1 for position, label in plex_for_columns(family, labels, n_columns).items()}
            channels = sorted(((label, index[label]) for label in labels), key=lambda channel: channel[1])
        return run_labels, experiment_type, channels
