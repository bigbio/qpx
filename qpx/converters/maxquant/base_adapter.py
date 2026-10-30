"""Shared base adapter for MaxQuant converters.

Provides common data-loading helpers used by
:class:`MaxQuantFeatureAdapter` and :class:`MaxQuantPgAdapter`.
"""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.maxquant.channels import channel_label, reporter_channels
from qpx.converters.maxquant.constants import FixedModification
from qpx.converters.utils import strip_uniprot_prefix

# SDRF ``PP`` values a fixed modification can be placed for from the peptide
# alone; protein termini need the peptide's protein position, which
# evidence.txt and msms.txt do not give.
_TERMINUS = {"anywhere": "", "any n-term": "N", "any c-term": "C"}

# Unimod names of quantification labels, which MaxQuant configures as labels,
# not as fixed modifications: isobaric tags (TMT, iTRAQ) and the isotopic
# ``Label:`` modifications of SILAC.
_LABEL_NAME = re.compile(r"(?i)tmt|itraq|^label:")

_NO_FIXED_MODS = (
    "Peptidoforms carry no fixed modification because %s; MaxQuant leaves fixed modifications "
    "out of 'Modified sequence' although calculated_mz includes them"
)


class MaxQuantBaseAdapter(BaseConverter):
    """Base adapter with MaxQuant-specific loading utilities.

    Subclasses inherit ``_load_sdrf()`` so the SDRF loading logic is
    defined once rather than duplicated across the feature and PG adapters.
    """

    # Fixed modifications added to peptidoforms (see _load_fixed_modifications).
    _fixed_mods: tuple[FixedModification, ...] = ()

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

    def _load_sdrf(self, sdrf_path: Optional[str], reporter_columns: int = 0) -> tuple[dict, str, dict[str, int]]:
        """Load SDRF and return ``(sample_map, experiment_type, channels)``.

        When *sdrf_path* is ``None`` the method returns LFQ defaults, and warns
        when the MaxQuant table has reporter columns that nothing maps.

        Args:
            sdrf_path: Filesystem path to the SDRF file, or ``None``.
            reporter_columns: Number of reporter channels in the MaxQuant table
                (:func:`~qpx.converters.maxquant.channels.reporter_column_count`).

        Returns:
            A 3-tuple of:
              - *sample_map*: ``{"<run stem>-<label>" -> sample_accession}``, with
                the label spelled as the *channels* keys (``{run stem -> ...}`` for LFQ)
              - *experiment_type*: e.g. ``"LFQ"``, ``"TMT11"``
              - *channels*: ``{label: 0-based reporter column}`` in column
                order (empty for LFQ); see
                :func:`~qpx.converters.maxquant.channels.reporter_channels`
        """
        if not sdrf_path:
            if reporter_columns:
                self.logger.warning(
                    "The MaxQuant table has %d reporter intensity columns, but without an SDRF no channel maps "
                    "to them; they are not converted and the MS1 intensity is written as LFQ",
                    reporter_columns,
                )
            return {}, "LFQ", {}

        from qpx.core.sdrf import SDRFHandler

        handler = SDRFHandler(sdrf_path)
        experiment_type, channels = reporter_channels(handler.sdrf_table.get("comment[label]"), reporter_columns)
        return handler.get_sample_map_run(label_key=channel_label), experiment_type, channels

    def _load_fixed_modifications(self, sdrf_path: Optional[str]) -> tuple[FixedModification, ...]:
        """The SDRF's fixed modifications (``MT=fixed``), as :func:`to_proforma` takes them.

        MaxQuant leaves fixed modifications out of ``Modified sequence`` while
        ``Mass`` and ``m/z`` include them, so the peptidoform needs them to match
        ``calculated_mz``. Each one is taken by its SDRF accession (``AC=UNIMOD:n``;
        without one, the Unimod accession pyopenms has for its ``NT`` name), on
        its target residues (``TA``, comma-separated) or peptide terminus (``PP``).
        Quantification labels (TMT, iTRAQ, SILAC) are left out: MaxQuant
        quantifies them as labels, reports ``Mass`` and ``m/z`` without the
        isobaric ones, and carries a SILAC label on the labelled rows only.
        """
        if not sdrf_path:
            self.logger.warning(_NO_FIXED_MODS, "no SDRF was given")
            return ()

        from qpx.core.sdrf import SDRFHandler

        params = SDRFHandler(sdrf_path).get_modification_parameters() or []
        fixed_mods = tuple(mod for mod in map(self._fixed_modification, filter(lambda p: p["fixed"], params)) if mod)
        if not fixed_mods:
            self.logger.warning(_NO_FIXED_MODS, "the SDRF declares none that can be placed")
        return fixed_mods

    def _fixed_modification(self, param: dict) -> Optional[FixedModification]:
        """One SDRF fixed modification as ``(accession, residues, terminus)``, or ``None`` to skip it.

        Quantification labels are skipped, whether the SDRF names them (``NT``)
        or only gives their accession: MaxQuant reports ``Mass`` and ``m/z``
        without isobaric tags, and a SILAC label is carried by the labelled
        rows only.
        """
        unimod = _sdrf_unimod(param)
        if _is_quantification_label(param.get("name"), unimod):
            self.logger.info(
                "SDRF fixed modification %s is a quantification label (TMT, iTRAQ or SILAC), which MaxQuant does "
                "not report as a modification; peptidoforms omit it",
                param.get("name") or unimod,
            )
            return None
        target = _fixed_modification_target(param)
        if unimod is None or target is None:
            self.logger.warning("Cannot place SDRF fixed modification %r on MaxQuant peptides; peptidoforms omit it", param)
            return None
        return (unimod, *target)


@lru_cache(maxsize=1)
def _unimod_names() -> dict[str, frozenset[str]]:
    """``{"UNIMOD:n": names}`` of the Unimod modifications in the pyopenms modification database."""
    import pyopenms as oms

    database = oms.ModificationsDB()
    names: dict[str, set[str]] = {}
    for index in range(database.getNumberOfModifications()):
        modification = database.getModification(index)
        names.setdefault(str(modification.getUniModAccession()).upper(), set()).add(str(modification.getId()))
    return {accession: frozenset(found) for accession, found in names.items() if accession}


@lru_cache(maxsize=1)
def _unimod_accessions() -> dict[str, str]:
    """``{name (case-folded): "UNIMOD:n"}`` for the names the pyopenms database gives one Unimod accession."""
    found: dict[str, set[str]] = {}
    for accession, names in _unimod_names().items():
        for name in names:
            found.setdefault(name.casefold(), set()).add(accession)
    return {name: min(accessions) for name, accessions in found.items() if len(accessions) == 1}


def _sdrf_unimod(param: dict) -> Optional[str]:
    """The ``UNIMOD:n`` accession of an SDRF modification: its ``AC``, else the Unimod accession of its ``NT``."""
    accession = re.fullmatch(r"(?i)\s*unimod:(\d+)\s*", param.get("accession") or "")
    if accession:
        return f"UNIMOD:{accession.group(1)}"
    return _unimod_accessions().get((param.get("name") or "").strip().casefold())


def _is_quantification_label(name: Optional[str], unimod: Optional[str]) -> bool:
    """Whether an SDRF modification is a TMT, iTRAQ or SILAC label, by its ``NT`` or its Unimod name."""
    if name and _LABEL_NAME.search(name):
        return True
    return unimod is not None and any(_LABEL_NAME.search(known) for known in _unimod_names().get(unimod, ()))


_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTUVWY")


def _fixed_modification_target(param: dict) -> Optional[tuple[str, str]]:
    """``(residues, terminus)`` of an SDRF fixed modification, or ``None`` when it cannot be placed.

    ``TA`` lists the target residues (comma-separated); a terminal ``PP``
    without a residue (``TA=N-term``, ``TA=X``) applies to any terminal residue.
    """
    terminus = _TERMINUS.get(str(param.get("position") or "Anywhere").strip().lower())
    tokens = str(param.get("target_amino_acid") or "").upper().replace(" ", "").split(",")
    residues = "".join(token for token in tokens if token in _RESIDUES)
    if terminus is None or (terminus == "" and not residues):
        return None
    return residues, terminus
