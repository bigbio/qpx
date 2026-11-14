"""Spectronaut report layouts: pivot columns, channel quantities and the views a report can produce."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Optional

import pyarrow.parquet as pq
from sdrf_pipelines.converters.channel_map import CHANNEL_LABELS

from qpx.converters.base import resolve_columns
from qpx.converters.channel_labels import normalize_label
from qpx.converters.mappings import get_field_mappings
from qpx.core.constants import FEATURE, PG
from qpx.core.query import _escape_sql_string

# A pivot report has one column per run and quantity: "[1] run.raw.EG.Channel1Quantity".
_PIVOT_COLUMN = re.compile(r"\[\d+\] (?P<run>.+?)\.(?P<metric>(?:PG|PEP|EG|FG)\.[^.]+)")
_CHANNEL_COLUMN = re.compile(r"EG\.Channel(\d+)Quantity")
_RUN_COLUMN = "R.FileName"
_CHARGE_COLUMN = "FG.Charge"
_PRECURSOR_ID_COLUMN = "EG.PrecursorId"
# Spectronaut writes this in place of a quantity its filters removed.
_FILTERED = "Filtered"

_SILAC_ORDER = list(CHANNEL_LABELS["silac3plex"].values())
_FEATURE_FIELDS = ("sequence", "modified_sequence", "charge", "run_file_name", "pg_accessions")


def report_header(path: str) -> list[str]:
    """The column names of a TSV or Parquet report, without reading its rows."""
    if str(path).endswith(".parquet"):
        return pq.read_schema(path).names
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        return handle.readline().rstrip("\r\n").split("\t")


def pivot_columns(columns: Iterable[str]) -> dict[str, dict[str, str]]:
    """Run -> metric -> column of a pivot report; empty for the long layout."""
    runs: dict[str, dict[str, str]] = {}
    for column in columns:
        if match := _PIVOT_COLUMN.fullmatch(column):
            runs.setdefault(match["run"], {})[match["metric"]] = column
    return runs


def long_columns(columns: Iterable[str]) -> set[str]:
    """The columns a report has in the long layout, which is what a pivot report melts to."""
    columns = list(columns)
    runs = pivot_columns(columns)
    if not runs:
        return set(columns)
    pivoted = {column for metrics in runs.values() for column in metrics.values()}
    melted = set(columns) - pivoted | {_RUN_COLUMN} | {metric for metrics in runs.values() for metric in metrics}
    return melted | {_CHARGE_COLUMN} if _PRECURSOR_ID_COLUMN in melted else melted


def channel_columns(columns: Iterable[str]) -> dict[int, str]:
    """Channel number -> ``EG.Channel<N>Quantity`` column, in channel order."""
    found = {int(match[1]): column for column in columns if (match := _CHANNEL_COLUMN.fullmatch(column))}
    return dict(sorted(found.items()))


def channel_labels(channels: Iterable[int], sdrf_labels: Optional[set[str]]) -> dict[int, str]:
    """Channel number -> label, lightest first; the report itself does not name its channels.

    The channels take the SDRF's SILAC labels in plex order, or ``SILAC
    light``/``medium``/``heavy`` without an SDRF. A single channel without
    SILAC labels is label-free. An SDRF that disagrees with the report raises.
    """
    channels = list(channels)
    declared = {normalize_label(label) for label in sdrf_labels or ()}
    if len(channels) == 1 and declared <= {"LFQ"}:
        return {channels[0]: "LFQ"}
    if not declared:
        if len(channels) > len(_SILAC_ORDER):
            raise ValueError(f"The report has {len(channels)} channels; pass an SDRF that declares their labels")
        return dict(zip(channels, ["SILAC light", "SILAC heavy"] if len(channels) == 2 else _SILAC_ORDER))
    labels = [label for label in _SILAC_ORDER if label in declared]
    if len(labels) != len(channels):
        raise ValueError(f"The report has {len(channels)} quantity channels, but the SDRF declares labels {sorted(declared)}")
    return dict(zip(channels, labels))


def unavailable_views(columns: Iterable[str]) -> dict[str, str]:
    """View -> why the report (in the long layout) cannot produce it."""
    columns = set(columns)
    reasons: dict[str, str] = {}
    feature = resolve_columns(get_field_mappings("spectronaut", FEATURE), columns)
    missing = [field for field in _FEATURE_FIELDS if field not in feature]
    if "intensity" not in feature and not channel_columns(columns):
        missing.append("intensity (FG.Quantity or EG.Channel<N>Quantity)")
    if missing:
        reasons[FEATURE] = f"the report has no column for {', '.join(missing)}"
    if channel_columns(columns):
        reasons[PG] = "the report quantifies channels, and PG.Quantity is not per channel"
    else:
        pg = resolve_columns(get_field_mappings("spectronaut", PG), columns)
        missing = [field for field in ("pg_accessions", "intensity", "run_file_name") if field not in pg]
        if missing:
            reasons[PG] = f"the report has no column for {', '.join(missing)}"
    return reasons


def _identifier(column: str) -> str:
    """A quoted SQL identifier; report columns are not valid unquoted names."""
    return '"' + column.replace('"', '""') + '"'


def _melted_metric(column: Optional[str], metric: str) -> str:
    """One run's value of a metric, with Spectronaut's ``Filtered`` as null."""
    if column is None:
        return f"NULL AS {_identifier(metric)}"
    text = f"CAST({_identifier(column)} AS VARCHAR)"
    if metric.endswith("Quantity"):
        return f"NULLIF(TRY_CAST(replace({text}, ',', '.') AS DOUBLE), 'NaN'::DOUBLE) AS {_identifier(metric)}"
    return f"NULLIF({text}, '{_FILTERED}') AS {_identifier(metric)}"


def melted_report_sql(columns: Iterable[str], source: str) -> str:
    """SQL melting the pivot table *source* to one row per precursor and run; unquantified rows are dropped."""
    columns = list(columns)
    runs = pivot_columns(columns)
    pivoted = {column for metrics in runs.values() for column in metrics.values()}
    fixed = [column for column in columns if column not in pivoted and column != _CHARGE_COLUMN]
    fixed_sql = [_identifier(column) for column in fixed]
    if _CHARGE_COLUMN in columns:
        fixed_sql.append(f"TRY_CAST({_identifier(_CHARGE_COLUMN)} AS BIGINT) AS {_identifier(_CHARGE_COLUMN)}")
    elif _PRECURSOR_ID_COLUMN in columns:
        # "_PEPTIDE_.2": the precursor id ends in the charge.
        charge = f"regexp_extract({_identifier(_PRECURSOR_ID_COLUMN)}, '\\.(\\d+)$', 1)"
        fixed_sql.append(f"TRY_CAST({charge} AS BIGINT) AS {_identifier(_CHARGE_COLUMN)}")

    metrics = sorted({metric for run_metrics in runs.values() for metric in run_metrics})
    quantified = " OR ".join(f"{_identifier(metric)} IS NOT NULL" for metric in metrics if metric.endswith("Quantity"))
    selects = []
    for run, run_metrics in runs.items():
        parts = [*fixed_sql, f"'{_escape_sql_string(run)}' AS {_identifier(_RUN_COLUMN)}"]
        parts.extend(_melted_metric(run_metrics.get(metric), metric) for metric in metrics)
        selects.append(f"SELECT {', '.join(parts)} FROM {source}")
    union = " UNION ALL ".join(selects)
    return f"SELECT * FROM ({union}) WHERE {quantified}" if quantified else union
