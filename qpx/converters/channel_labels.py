"""
Shared quantification channel-label resolution.

The converters put the *same* canonical reporter labels into
``intensities[].label`` so TMT/iTRAQ/LFQ results are consistent and join to
``sample``/``run`` the same way whichever tool wrote them.

The canonical vocabulary is the shared sdrf-pipelines ``channel_map`` (single
source of truth for label spellings, incl. the TMT10 ch10 = ``TMT131`` vs
TMT11 ch10 = ``TMT131N`` distinction). The plex is resolved, most-solid first,
from the SDRF ``comment[label]`` set (ground truth) and only then from the
channel indices present in the data.
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from typing import Iterable, Optional

import pandas as pd
from sdrf_pipelines.converters.channel_map import CHANNEL_LABELS
from sdrf_pipelines.converters.channel_map import normalize_label as _normalize_label
from sdrf_pipelines.converters.openms.utils import infer_itraqplex, infer_tmtplex

from qpx.core.files import run_file_stem

logger = logging.getLogger(__name__)

_MS1_LABEL = re.compile(r"(?i)^(?:silac|dimethyl|mtraq)")

__all__ = [
    "read_sdrf_labels",
    "experiment_runs_from_sdrf",
    "fraction_groups_from_sdrf",
    "experiment_type_from_labels",
    "resolve_channel_labels",
    "normalize_label",
]

_FAMILY_PREFIX = {"TMT": "tmt", "ITRAQ": "itraq", "iTRAQ": "itraq"}
_CANONICAL_LABELS = {label.casefold(): label for channel_labels in CHANNEL_LABELS.values() for label in channel_labels.values()}
_SILAC_LABELS = set(CHANNEL_LABELS["silac3plex"].values())
# TMTpro-134C, iTRAQ4plex-114, "iTRAQ reagent 114"
_REPORTER_SPELLING = re.compile(r"(?i)(TMT|iTRAQ)(?:(\d+)plex|pro)?[\s_-]*(?:reagent)?[\s_-]*(\d{3}[NC]?)")
_BARE_REPORTER = re.compile(r"\d{3}[NC]?D?")
_FIRST_TMT_CHANNEL = 126
_ONTOLOGY_NAME = re.compile(r"(?:^|;)\s*NT\s*=\s*([^;]+)", re.IGNORECASE)


def _reporter_channel(text: str) -> str:
    """``iTRAQ4plex-114`` -> ``ITRAQ114``; a channel its named plex does not have keeps its spelling.

    A bare ``126`` names one family: iTRAQ reporters stop at 121, TMT starts at 126.
    """
    if _BARE_REPORTER.fullmatch(text):
        return f"{'TMT' if int(text[:3]) >= _FIRST_TMT_CHANNEL else 'ITRAQ'}{text}"
    spelling = _REPORTER_SPELLING.fullmatch(text)
    if spelling is None:
        return text
    family, plex_size, channel = spelling.groups()
    label = f"{family}{channel}".upper()
    plex = CHANNEL_LABELS.get(f"{family}{plex_size}plex".lower())
    return text if plex is not None and label not in plex.values() else label


def normalize_label(label: str) -> str:
    """Return the canonical spelling of a known SDRF channel label."""
    text = label.strip()
    match = _ONTOLOGY_NAME.search(text)
    if match:
        text = match.group(1).strip()
    normalized = _normalize_label(_reporter_channel(text))
    return _CANONICAL_LABELS.get(normalized.casefold(), normalized)


def read_sdrf_labels(sdrf_path: Optional[str]) -> Optional[set[str]]:
    """
    Return the distinct ``comment[label]`` values declared in the SDRF.

    This is the experiment's ground-truth channel set. Returns ``None`` when no
    SDRF is available or it has no label column.
    """
    if not sdrf_path:
        return None
    try:
        df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError):
        # Unreadable, malformed, or otherwise invalid SDRF -> no ground truth;
        # fall back to count-based inference rather than aborting the conversion.
        return None
    col = next((c for c in df.columns if c.strip().lower() == "comment[label]"), None)
    if col is None:
        return None
    labels = {str(v).strip() for v in df[col].dropna() if str(v).strip()}
    return labels or None


def experiment_runs_from_sdrf(sdrf_path: Optional[str]) -> Optional[dict[str, list[str]]]:
    """Map each SDRF ``source name`` to its raw run files (``comment[data file]`` stems).

    The SDRF is the tool-agnostic ground truth for which raw files form one
    quantification unit: rows sharing a ``source name`` are the fractions of one
    sample. This is the fallback source for pg ``grouped_runs`` when a converter's
    own experiment->runs result mapping (FragPipe ``experiment_annotation.tsv`` /
    MaxQuant ``evidence.txt``) is absent. Files are returned distinct and in SDRF
    row order (fraction order preserved, not sorted). Returns ``None`` when no
    SDRF is available or it lacks ``source name`` / ``comment[data file]``.
    """
    if not sdrf_path:
        return None
    try:
        df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError):
        return None
    cols = {c.strip().lower(): c for c in df.columns}
    source_col = cols.get("source name")
    file_col = cols.get("comment[data file]")
    if source_col is None or file_col is None:
        return None
    mapping: dict[str, list[str]] = {}
    for source, file_value in zip(df[source_col], df[file_col]):
        if pd.isna(source) or pd.isna(file_value):
            continue
        experiment = str(source).strip()
        run = run_file_stem(str(file_value))
        if not experiment or not run:
            continue
        runs = mapping.setdefault(experiment, [])
        if run not in runs:
            runs.append(run)
    return mapping or None


def _sdrf_unit_columns(df) -> Optional[dict[str, Optional[str]]]:
    """Resolve the SDRF columns that define a quantification unit, or None."""
    cols = {c.strip().lower(): c for c in df.columns}
    file_col = cols.get("comment[data file]")
    source_col = cols.get("source name")
    if file_col is None or source_col is None:
        return None
    return {
        "file": file_col,
        "source": source_col,
        "label": cols.get("comment[label]"),
        "fraction": cols.get("comment[fraction identifier]"),
        "brep": next((cols[k] for k in cols if "biological replicate" in k), None),
        "trep": next((cols[k] for k in cols if "technical replicate" in k), None),
    }


def _safe_fraction(value: object) -> int:
    try:
        return int(str(value).strip()) if value is not None else 0
    except (ValueError, TypeError):
        return 0


def _collect_run_signatures(df, c: dict[str, Optional[str]]):
    """Per run stem: its ``(source,label,brep,trep)`` signature set + fraction number."""
    signature: dict[str, set] = defaultdict(set)
    fraction: dict[str, int] = {}
    first_seen: list[str] = []
    for row in df.to_dict("records"):
        raw = row.get(c["file"])
        if raw is None or not str(raw).strip():
            continue
        run = run_file_stem(str(raw))
        if run not in signature:
            first_seen.append(run)
        signature[run].add(tuple(str(row.get(c[k], "")).strip() if c[k] else "" for k in ("source", "label", "brep", "trep")))
        if run not in fraction:
            fraction[run] = _safe_fraction(row.get(c["fraction"])) if c["fraction"] else 0
    return signature, fraction, first_seen


def fraction_groups_from_sdrf(sdrf_path: Optional[str]) -> Optional[dict[str, list[str]]]:
    """Map each SDRF raw file to the ``grouped_runs`` of its quantification unit.

    A quantification unit is the set of raw files that differ **only in fraction**:
    they carry the same set of ``(source name, label, biological replicate,
    technical replicate)`` assignments. For LFQ, a sample's fractions group
    together; for TMT/iTRAQ, all fraction files of one labelled mixture group
    together (they share the same channel→sample assignments). Technical
    replicates (which differ in ``technical replicate``) stay in separate units.

    Returns a dict mapping each ``run_file_name`` stem to its ``grouped_runs``
    list (the unit's files, distinct, ordered by fraction identifier), or ``None``
    when no SDRF or the required columns are absent.
    """
    if not sdrf_path:
        return None
    try:
        df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError):
        return None
    cols = _sdrf_unit_columns(df)
    if cols is None:
        return None
    signature, fraction, first_seen = _collect_run_signatures(df, cols)
    if not signature:
        return None

    groups: dict[frozenset, list[str]] = defaultdict(list)
    for run in first_seen:
        groups[frozenset(signature[run])].append(run)

    run_to_grouped: dict[str, list[str]] = {}
    for members in groups.values():
        ordered: list[str] = []
        for run in sorted(members, key=lambda r: (fraction.get(r, 0), first_seen.index(r))):
            if run not in ordered:
                ordered.append(run)
        for run in ordered:
            run_to_grouped[run] = ordered
    return run_to_grouped or None


def experiment_type_from_labels(sdrf_labels: Optional[set[str]]) -> str:
    """Classify TMT, iTRAQ, SILAC or LFQ from the SDRF channel labels."""
    if not sdrf_labels:
        return "LFQ"
    upper = {label.upper() for label in sdrf_labels}
    if any("TMT" in label for label in upper):
        return "TMT"
    if any("ITRAQ" in label for label in upper):
        return "iTRAQ"
    if any(normalize_label(label) in _SILAC_LABELS for label in sdrf_labels):
        return "SILAC"
    return "LFQ"


# Unimod accession of the reagent each isobaric plex is labelled with.
_PLEX_UNIMOD = {"tmt6plex": 737, "tmt10plex": 737, "tmt11plex": 737, "tmt16plex": 2016, "tmt18plex": 2016}
_PLEX_UNIMOD |= {"itraq4plex": 214, "itraq8plex": 730}


def label_unimod(sdrf_labels: Optional[set[str]]) -> Optional[int]:
    """UNIMOD accession of the isobaric reagent the SDRF channel labels imply, or ``None``."""
    family = _FAMILY_PREFIX.get(experiment_type_from_labels(sdrf_labels))
    if family is None:
        return None
    canonical = {normalize_label(label) for label in sdrf_labels if label}
    infer = infer_tmtplex if family == "tmt" else infer_itraqplex
    try:
        return _PLEX_UNIMOD.get(infer(canonical))
    except ValueError:
        return None


def plex_for_columns(family: str, labels: set[str], n_columns: int) -> dict[int, str]:
    """``{1-based position: label}`` of the plex holding *labels*, sized *n_columns* when one is.

    A label's position depends on the plex: ``ITRAQ114`` is channel 1 of iTRAQ4
    and channel 2 of iTRAQ8, so four labels in an 8-column table read 8-plex columns.
    A table with one column per label and no plex of that size was set up with
    these channels only, which take its columns in plex order: TMT128N to TMT131C
    of TMT11 in eight columns read columns 1 to 8.
    """
    plex = resolve_channel_labels(family, labels)
    if not plex or not n_columns or len(plex) == n_columns:
        return plex
    for name, channels in CHANNEL_LABELS.items():
        if name.startswith(_FAMILY_PREFIX[family]) and len(channels) == n_columns and labels <= set(channels.values()):
            return channels
    if len(labels) == n_columns:
        return dict(enumerate((label for _, label in sorted(plex.items()) if label in labels), start=1))
    logger.warning(
        "The table has %d reporter columns, but the SDRF labels %s are channels of a %d-channel %s plex; "
        "each is read from the column of its position in that plex, so check the labels if the run used another plex",
        n_columns,
        sorted(labels),
        len(plex),
        family,
    )
    return plex


def _channel_indices(channel_indices: Optional[Iterable]) -> set[int]:
    """Integer channel ordinals in ``channel_indices``; values that are not numbers are skipped."""
    indices: set[int] = set()
    for value in () if channel_indices is None else channel_indices:
        try:
            indices.add(int(float(value)))
        except (TypeError, ValueError):
            continue
    return indices


def resolve_channel_labels(
    experiment_type: str,
    sdrf_labels: Optional[set[str]] = None,
    channel_indices: Optional[Iterable] = None,
) -> dict[int, str]:
    """
    Resolve the ``{1-based channel index -> canonical label}`` map for a run.

    Resolution order (most-solid first):

    1. From the SDRF ``comment[label]`` set via sdrf-pipelines
       ``infer_tmtplex`` / ``infer_itraqplex`` — robust even when a channel is
       empty in the data.
    2. Fallback: the smallest plex of the reagent family whose channel count
       covers the highest index present in ``channel_indices``.

    Returns ``{}`` for non-isobaric or unresolvable input (caller falls back to
    the raw value / ``"LFQ"``).
    """
    prefix = _FAMILY_PREFIX.get(experiment_type)
    if prefix is None:
        return {}

    # 1. SDRF-declared plex (ground truth).
    if sdrf_labels:
        canonical = {normalize_label(label) for label in sdrf_labels if label}
        infer = infer_tmtplex if prefix == "tmt" else infer_itraqplex
        try:
            return CHANNEL_LABELS[infer(canonical)]
        except (ValueError, KeyError):
            pass

    # 2. Fallback: infer plex from the channel indices present in the data.
    # NB: channel_indices may be a numpy array (df[col].values), so test for
    # None explicitly — `channel_indices or ()` raises on array truthiness.
    indices = _channel_indices(channel_indices)
    if not indices:
        return {}
    max_index = max(indices)

    family = {plex: labels for plex, labels in CHANNEL_LABELS.items() if plex.startswith(prefix)}
    for plex in sorted(family, key=lambda p: len(family[p])):
        if len(family[plex]) >= max_index:
            return family[plex]
    return {}


def is_ms1_label(text: str) -> bool:
    """A SILAC, dimethyl or mTRAQ label names its quantity already, whatever the experiment type says."""
    return _MS1_LABEL.match(text) is not None
