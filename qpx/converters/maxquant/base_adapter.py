"""Shared base adapter for MaxQuant converters.

Provides common data-loading helpers used by
:class:`MaxQuantFeatureAdapter` and :class:`MaxQuantPgAdapter`.
"""

from __future__ import annotations

from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.utils import strip_uniprot_prefix


class MaxQuantBaseAdapter(BaseConverter):
    """Base adapter with MaxQuant-specific loading utilities.

    Subclasses inherit ``_load_sdrf()`` so the SDRF loading logic is
    defined once rather than duplicated across the feature and PG adapters.
    """

    @staticmethod
    def _norm_uniprot_id(pid: str) -> str:
        """Normalise a MaxQuant UniProt protein ID to its bare accession.

        Thin wrapper over :func:`qpx.converters.utils.strip_uniprot_prefix`::

            sp|P55011|S12A2_HUMAN  →  P55011
            tr|A0A3B3IS91|..._HUMAN  →  A0A3B3IS91
            P55011  →  P55011  (returned as-is)
        """
        return strip_uniprot_prefix(pid)

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

    def _load_sdrf(self, sdrf_path: Optional[str]) -> tuple[dict, str, list]:
        """Load SDRF and return ``(sample_map, experiment_type, tmt_channels)``.

        When *sdrf_path* is ``None`` the method returns sensible LFQ defaults.

        Args:
            sdrf_path: Filesystem path to the SDRF file, or ``None``.

        Returns:
            A 3-tuple of:
              - *sample_map*: ``{run_file -> sample_accession}``
              - *experiment_type*: e.g. ``"LFQ"``, ``"TMT11"``
              - *tmt_channels*: sorted list of TMT channel labels (empty for LFQ)
        """
        if not sdrf_path:
            return {}, "LFQ", []

        from qpx.core.sdrf import SDRFHandler

        handler = SDRFHandler(sdrf_path)
        sample_map = handler.get_sample_map_run()
        experiment_type = handler.get_experiment_type_from_sdrf()

        tmt_channels: list[str] = []
        if experiment_type and ("TMT" in experiment_type.upper() or "ITRAQ" in experiment_type.upper()):
            labels = handler.sdrf_table.get("comment[label]")
            if labels is not None:
                from qpx.converters.maxquant.constants import TMT_LABEL_TO_MQ_COL

                # Normalize to canonical uppercase form so TMT_LABEL_TO_MQ_COL lookups
                # are robust to inconsistent SDRF casing / leading-trailing whitespace.
                raw_labels = [str(lbl).strip().upper() for lbl in labels.unique() if lbl and str(lbl).strip()]
                _MAX_IDX = max(TMT_LABEL_TO_MQ_COL.values()) + 1
                tmt_channels = sorted(raw_labels, key=lambda lbl: TMT_LABEL_TO_MQ_COL.get(lbl, _MAX_IDX))

        return sample_map, experiment_type or "LFQ", tmt_channels
