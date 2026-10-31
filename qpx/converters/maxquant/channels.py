"""SDRF channel labels -> MaxQuant ``Reporter intensity N`` columns.

MaxQuant numbers the reporter columns of an isobaric run by channel position
in the plex (mass order, N before C isomer): column ``k`` (0-based) is the
``k + 1``-th channel. Which label sits at a position depends on the plex
(``TMT131`` is channel 6 of TMT6 but channel 10 of TMT10; ``ITRAQ114`` is
channel 1 of iTRAQ4 but channel 2 of iTRAQ8), so the plex is resolved from the
SDRF labels with the shared sdrf-pipelines channel map, and the number of
reporter columns in the MaxQuant table chooses between plexes that hold the
same labels.
"""

from __future__ import annotations

import logging
import re
from typing import Iterable, Optional

import pandas as pd
from sdrf_pipelines.converters.channel_map import CHANNEL_LABELS

from qpx.converters.channel_labels import experiment_type_from_labels, normalize_label, resolve_channel_labels

logger = logging.getLogger(__name__)

# MaxQuant and SDRF spellings of a reporter channel: TMT127N, TMT6plex-127,
# TMTpro-134C, iTRAQ114, iTRAQ4plex-114, "iTRAQ reagent 114".
_REPORTER_LABEL = re.compile(r"(TMT|ITRAQ)(?:\d+PLEX|PRO)?[\s_-]*(?:REAGENT)?[\s_-]*(\d{3}[NC]?)")

_REPORTER_COLUMN = re.compile(r"Reporter intensity (\d+)(?: .+)?")


def reporter_column_count(columns: Iterable[str]) -> int:
    """How many reporter channels a MaxQuant table has (``Reporter intensity N [experiment]``)."""
    return len({match.group(1) for column in columns if (match := _REPORTER_COLUMN.fullmatch(str(column)))})


def _canonical_reporter_label(label: str) -> str:
    """The sdrf-pipelines spelling of a reporter channel label (``iTRAQ4plex-114`` -> ``ITRAQ114``)."""
    match = _REPORTER_LABEL.fullmatch(label.strip().upper())
    return normalize_label(f"{match.group(1)}{match.group(2)}") if match else normalize_label(label)


def channel_label(value) -> str:
    """An SDRF ``comment[label]`` value spelled as the run view writes it, which keys the channels."""
    return normalize_label(str(value))


def _sdrf_labels(values: Optional[Iterable]) -> set[str]:
    """The distinct ``comment[label]`` values, spelled as the run view writes them."""
    texts = (str(value).strip() for value in (values if values is not None else ()) if pd.notna(value))
    return {channel_label(text) for text in texts if text}


def _resolve_plex(family: str, labels: set[str], n_columns: int) -> dict[int, str]:
    """``{1-based position: label}`` of the plex holding ``labels``, sized ``n_columns`` when one is."""
    plex = resolve_channel_labels(family, labels)
    if not plex or not n_columns or len(plex) == n_columns:
        return plex
    prefix = family.lower()
    for name, channels in CHANNEL_LABELS.items():
        if name.startswith(prefix) and len(channels) == n_columns and labels <= set(channels.values()):
            return channels
    logger.warning(
        "The MaxQuant table has %d reporter columns, but the SDRF labels %s are channels of a %d-channel %s plex; "
        "each is read from the column of its position in that plex, so check the labels if the run used another plex",
        n_columns,
        sorted(labels),
        len(plex),
        family,
    )
    return plex


def _isobaric_family(labels: set[str]) -> Optional[str]:
    """``"TMT"`` or ``"iTRAQ"`` for isobaric labels, ``None`` for label-free ones; SILAC raises."""
    family = experiment_type_from_labels(labels)
    if family == "SILAC":
        raise ValueError(
            f"SDRF labels {sorted(labels)} describe a SILAC experiment, which convert maxquant does not support "
            "yet (bigbio/qpx#308)"
        )
    return family if family in ("TMT", "iTRAQ") else None


def reporter_channels(sdrf_labels: Optional[Iterable], n_columns: int = 0) -> tuple[str, dict[str, int]]:
    """Resolve the SDRF labels to ``(experiment_type, {label: 0-based reporter column})``.

    ``experiment_type`` is ``"TMT<n>"`` or ``"ITRAQ<n>"`` for an ``n``-channel
    plex, with the channels ordered by column, or ``"LFQ"`` with no channels.
    Partial plexes are resolved like full ones. ``n_columns`` is the number of
    reporter columns in the MaxQuant table (0 when unknown).

    Raises:
        ValueError: for SILAC labels, which this converter does not support
            yet, and for isobaric labels that do not all belong to one known
            plex (rather than reading them from the wrong columns).
    """
    labels = _sdrf_labels(sdrf_labels)
    canonical = {label: _canonical_reporter_label(label) for label in labels}
    family = _isobaric_family(set(canonical.values()))
    if family is None:
        return "LFQ", {}
    position = {label: index for index, label in _resolve_plex(family, set(canonical.values()), n_columns).items()}
    unknown = sorted(label for label, channel in canonical.items() if channel not in position)
    if unknown:
        raise ValueError(f"SDRF labels {unknown} are not channels of a {family} plex that holds all of {sorted(labels)}")
    ordered = sorted(labels, key=lambda label: position[canonical[label]])
    return f"{family.upper()}{len(position)}", {label: position[canonical[label]] - 1 for label in ordered}
