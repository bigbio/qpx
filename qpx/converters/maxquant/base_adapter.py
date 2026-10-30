"""Shared base adapter for MaxQuant converters.

Provides common data-loading helpers used by
:class:`MaxQuantFeatureAdapter` and :class:`MaxQuantPgAdapter`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from qpx.converters.base import BaseConverter
from qpx.converters.channel_labels import label_unimod, read_sdrf_labels
from qpx.converters.maxquant.constants import FixedMods, parameters_fixed_mods, sdrf_fixed_mods
from qpx.converters.ptm import NTERM
from qpx.converters.utils import strip_uniprot_prefix


class MaxQuantBaseAdapter(BaseConverter):
    """Base adapter with MaxQuant-specific loading utilities.

    Subclasses inherit ``_load_sdrf()`` so the SDRF loading logic is
    defined once rather than duplicated across the feature and PG adapters.
    """

    # Fixed modifications written into peptidoforms; set by ``convert``.
    _fixed: FixedMods = ()

    @staticmethod
    def _norm_uniprot_id(pid: str) -> str:
        """Normalise a MaxQuant UniProt protein ID to its bare accession.

        Thin wrapper over :func:`qpx.converters.utils.strip_uniprot_prefix`::

            sp|P55011|S12A2_HUMAN  →  P55011
            tr|A0A3B3IS91|..._HUMAN  →  A0A3B3IS91
            P55011  →  P55011  (returned as-is)
        """
        return strip_uniprot_prefix(pid)

    @staticmethod
    def _fixed_mods(sdrf_path: Optional[str], source_path: str) -> FixedMods:
        """Fixed modifications as ``(site, "UNIMOD:n")`` pairs; *site* is a residue or ``"N-term"``.

        MaxQuant leaves fixed modifications and isobaric labels out of
        ``Modified sequence``. They are read from the SDRF ``MT=fixed`` entries,
        else from the ``parameters.txt`` MaxQuant writes next to *source_path*;
        an isobaric label neither lists is added from the SDRF channel labels.
        """
        fixed: FixedMods = ()
        if sdrf_path:
            from qpx.converters.sdrf import SdrfConverter
            from qpx.core.sdrf import SDRFHandler

            fixed = sdrf_fixed_mods(SdrfConverter.parse_modification_params(SDRFHandler(sdrf_path).sdrf_table))
        fixed = fixed or parameters_fixed_mods(Path(source_path).with_name("parameters.txt"))
        label = label_unimod(read_sdrf_labels(sdrf_path))
        if label is not None and all(tag != f"UNIMOD:{label}" for _, tag in fixed):
            fixed += ((NTERM, f"UNIMOD:{label}"), ("K", f"UNIMOD:{label}"))
        return fixed

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
