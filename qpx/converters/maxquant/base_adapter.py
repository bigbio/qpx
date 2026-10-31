"""Shared base adapter for MaxQuant converters.

Provides common data-loading helpers used by
:class:`MaxQuantFeatureAdapter` and :class:`MaxQuantPgAdapter`.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.channel_labels import label_unimod, normalize_label, read_sdrf_labels, resolve_channel_labels
from qpx.converters.maxquant.constants import FixedMods, parameters_fixed_mods, sdrf_fixed_mods, summary_fixed_mods
from qpx.converters.ptm import NTERM
from qpx.converters.utils import canonical_accession
from qpx.core.files import run_file_stem


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

    @staticmethod
    def _fixed_mods(sdrf_path: Optional[str], source_path: str) -> FixedMods:
        """Fixed modifications as ``(site, "UNIMOD:n")`` pairs; *site* is a residue or ``"N-term"``.

        MaxQuant leaves fixed modifications and isobaric labels out of
        ``Modified sequence``. They are read from the SDRF ``MT=fixed`` entries,
        else from the ``parameters.txt`` or ``summary.txt`` MaxQuant writes next
        to *source_path*; an isobaric label none lists is added from the SDRF
        channel labels.
        """
        fixed: FixedMods = ()
        if sdrf_path:
            from qpx.converters.sdrf import SdrfConverter
            from qpx.core.sdrf import SDRFHandler

            fixed = sdrf_fixed_mods(SdrfConverter.parse_modification_params(SDRFHandler(sdrf_path).sdrf_table))
        folder = Path(source_path)
        fixed = (
            fixed
            or parameters_fixed_mods(folder.with_name("parameters.txt"))
            or summary_fixed_mods(folder.with_name("summary.txt"))
        )
        label = label_unimod(read_sdrf_labels(sdrf_path))
        if label is not None and all(tag != f"UNIMOD:{label}" for _, tag in fixed):
            fixed += ((NTERM, f"UNIMOD:{label}"), ("K", f"UNIMOD:{label}"))
        return fixed

    @staticmethod
    def _load_sdrf(sdrf_path: Optional[str]) -> tuple[dict[str, set[str]], str, list[tuple[str, int]]]:
        """Load SDRF and return ``(run_labels, experiment_type, channels)``.

        When *sdrf_path* is ``None`` the method returns LFQ defaults.

        Args:
            sdrf_path: Filesystem path to the SDRF file, or ``None``.

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
            index = {label: position - 1 for position, label in resolve_channel_labels(family, labels).items()}
            channels = sorted(((label, index[label]) for label in labels), key=lambda channel: channel[1])
        return run_labels, experiment_type, channels
