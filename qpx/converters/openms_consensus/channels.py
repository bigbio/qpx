"""Quantification channels of a consensusXML: isobaric, MS1-labelled or label-free.

The MS1-label roles follow OpenMS' ``ArrowIOHelpers::qpxIntensityLabels``, so this
converter and the native ``-out_qpx`` writer name the same channel the same way.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Optional

import pyopenms as oms

from qpx.converters.channel_labels import normalize_label

_CHANNEL_RE = re.compile(r"(tmt|itraq)\d*plex?_(\d+[NC]?)", re.IGNORECASE)
_MASS_TAG_RE = re.compile(r"^(dimethyl|mtraq)(\d+)$", re.IGNORECASE)
_LIGHT, _MEDIUM, _HEAVY = "SILAC light", "SILAC medium", "SILAC heavy"
_MEDIUM_MARKERS = ("arg6", "lys4", "lys6")
_HEAVY_MARKERS = ("arg10", "lys8")
_MS1_LABEL_PREFIXES = ("SILAC", "DIMETHYL", "MTRAQ")

# The label state MS1LabeledWorkflow records on every identification.
LABELED_SEQUENCE = "MS1Label:labeled_sequence"
REMOVED_LABELS = "MS1Label:removed_labels"
CHANNEL = "MS1Label:channel"


def _canonical_channel(label: Optional[str]) -> str:
    """Canonicalize an OpenMS map label to the QPX channel label.

    ``tmt6plex_126`` -> ``TMT126``; ``itraq4plex_114`` -> ``ITRAQ114``. Falls back
    to :func:`normalize_label`, then the raw string, so unknown labels are kept
    rather than dropped.
    """
    if not label:
        return "LFQ"
    m = _CHANNEL_RE.search(str(label))
    if m:
        family = "TMT" if m.group(1).lower() == "tmt" else "ITRAQ"
        return normalize_label(f"{family}{m.group(2).upper()}")
    return normalize_label(str(label))


def _map_label(raw_label) -> str:
    """Resolve an isobaric consensusXML map label to its channel, anything else to ``LFQ``.

    Detection is by the label itself, NOT ``experiment_type``: quantms/OpenMS
    write ``experiment_type="label-free"`` even for isobaric runs (the maps still
    carry ``tmt6plex_126`` / ``itraq4plex_114`` labels; label-free maps carry the
    literal ``"label-free"``). MS1 labels depend on the other columns of the run:
    see :func:`column_labels`.
    """
    return _canonical_channel(raw_label) if _CHANNEL_RE.search(str(raw_label or "")) else "LFQ"


def _silac_role(label: str) -> Optional[str]:
    """Role a map label names; ``None`` when it names none or contradicts itself."""
    lower = label.lower()
    if lower in ("", "no_label", "silac light", "light"):
        return _LIGHT
    medium = lower in ("silac medium", "medium") or any(marker in lower for marker in _MEDIUM_MARKERS)
    heavy = lower in ("silac heavy", "heavy") or any(marker in lower for marker in _HEAVY_MARKERS)
    if medium == heavy:
        return None
    return _MEDIUM if medium else _HEAVY


def _silac_labels(source: str, raw: list[str]) -> list[str]:
    """Roles of one run's columns; a two-plex calls its non-light channel heavy, as QPX does."""
    roles = [_silac_role(label) for label in raw]
    valid = None not in roles and len(raw) in (2, 3) and roles.count(_LIGHT) == 1
    if valid and len(raw) == 2:
        roles = [_LIGHT if role == _LIGHT else _HEAVY for role in roles]
    if not valid or len(set(roles)) != len(raw):
        raise ValueError(f"Run {source!r} has SILAC labels {raw} that are not an unambiguous two- or three-plex")
    return roles


def _source_labels(source: str, raw: list[str]) -> list[str]:
    """QPX labels of the columns of one source run."""
    if any(_CHANNEL_RE.search(label) for label in raw):
        return [_map_label(label) for label in raw]
    tags = [_MASS_TAG_RE.match(label) for label in raw]
    if any(tags):
        if not all(tags) or len({tag.group(0).lower() for tag in tags}) != len(raw):
            raise ValueError(f"Run {source!r} mixes mass-tag labels with others: {raw}")
        return [f"{tag.group(1).upper()}{tag.group(2)}" for tag in tags]
    if any("silac" in label.lower() or _silac_role(label) in (_MEDIUM, _HEAVY) for label in raw):
        return _silac_labels(source, raw)
    if len(set(raw)) > 1:
        # One intensity per channel is the contract: these cannot share one label.
        raise ValueError(f"Run {source!r} has channel labels QPX cannot represent: {raw}")
    return ["LFQ"] * len(raw)


def _columns_by_source(headers) -> dict[str, list[int]]:
    sources: dict[str, list[int]] = defaultdict(list)
    for idx in sorted(headers):
        sources[headers[idx].filename].append(idx)
    return sources


def column_labels(cm) -> dict[int, str]:
    """Map index -> QPX channel label. Raises ``ValueError`` for labels QPX cannot name."""
    headers = cm.getColumnHeaders()
    labels: dict[int, str] = {}
    for source, columns in _columns_by_source(headers).items():
        labels.update(zip(columns, _source_labels(source, [str(headers[idx].label or "") for idx in columns])))
    return labels


def column_channel(header, position: int) -> int:
    """1-based design channel of a column: its ``channel_id``, else its place in the run."""
    if header.metaValueExists("channel_id"):
        return int(header.getMetaValue("channel_id")) + 1
    return position + 1


def channels_by_column(cm) -> dict[int, int]:
    """Map index -> 1-based design channel within its source run."""
    headers = cm.getColumnHeaders()
    return {
        idx: column_channel(headers[idx], position)
        for columns in _columns_by_source(headers).values()
        for position, idx in enumerate(columns)
    }


def is_ms1_label(label: str) -> bool:
    """Whether a QPX channel label is an MS1 label (SILAC, dimethyl, mTRAQ)."""
    return str(label).startswith(_MS1_LABEL_PREFIXES)


def ratio_suffix(label: str) -> str:
    """Channel suffix used in ``ratio_<num>_<den>``: ``h`` for SILAC heavy, ``4`` for DIMETHYL4."""
    return label.rsplit(" ", 1)[-1][0] if label.startswith("SILAC") else re.sub(r"\D", "", label)


def ratio_name(numerator: str, reference: str, normalized: bool = False) -> str:
    """``ratio_<num>_<den>`` for two labels, optionally ``_normalized``."""
    name = f"ratio_{ratio_suffix(numerator)}_{ratio_suffix(reference)}"
    return f"{name}_normalized" if normalized else name


def list_value(value) -> list[str]:
    """A list-typed meta value as strings: pyopenms gives a list, the streaming reader ``"[a, b]"``."""
    if isinstance(value, str):
        text = value.strip().removeprefix("[").removesuffix("]")
        return [item.strip() for item in text.split(",")] if text else []
    return [item.decode() if isinstance(item, (bytes, bytearray)) else str(item) for item in value]


def matched_sequence(hit):
    """The peptidoform a spectrum was matched with (MS1 labels included), else the hit's sequence."""
    if hit.metaValueExists(LABELED_SEQUENCE):
        labeled = str(hit.getMetaValue(LABELED_SEQUENCE))
        if labeled:
            return oms.AASequence.fromString(labeled)
    return hit.getSequence()


def label_state_params(hit) -> list[dict] | None:
    """The label state of a hit as ``cv_params``; ``None`` when it records none."""
    params = [
        {"cv_name": key, "cv_value": str(hit.getMetaValue(key))}
        for key in (LABELED_SEQUENCE, REMOVED_LABELS, CHANNEL)
        if hit.metaValueExists(key)
    ]
    return params or None


class Ms1Channels:
    """Channel lookups of an MS1-labelled consensusXML; falsy for any other input."""

    def __init__(self, cm, runs: dict[int, str]):
        labels = column_labels(cm)
        channels = channels_by_column(cm)
        self.by_map = {idx: label for idx, label in labels.items() if is_ms1_label(label)}
        self.by_run_channel = {(runs[idx], channels[idx]): label for idx, label in self.by_map.items()}

    def __bool__(self) -> bool:
        return bool(self.by_map)

    def identified_label(self, pid, hit, run: str) -> Optional[str]:
        """Channel of the identified partner: ``MS1Label:channel``, else the identification's map."""
        channel = int(hit.getMetaValue(CHANNEL)) if hit.metaValueExists(CHANNEL) else 0
        if channel > 0:
            return self.by_run_channel.get((run, channel))
        if pid.metaValueExists("map_index"):
            return self.by_map.get(int(pid.getMetaValue("map_index")))
        return None

    def run_labels(self, run: str) -> list[str]:
        """The labels of one run, lightest first."""
        return [label for (other, _), label in sorted(self.by_run_channel.items()) if other == run]
