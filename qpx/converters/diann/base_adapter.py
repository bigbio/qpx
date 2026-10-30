"""Shared base adapter for DIA-NN converters.

Provides common data-loading helpers used by both
:class:`DiannFeatureAdapter` and :class:`DiannPgAdapter`.
"""

from __future__ import annotations

import os
import re

import duckdb
import pandas as pd
from sdrf_pipelines.converters.diann.constants import PLEXDIA_REGISTRY

from qpx.converters.base import BaseConverter
from qpx.converters.channel_labels import normalize_label, read_sdrf_labels
from qpx.core.sql import escape_path, sql_build, validate_identifier

# Multiplier from on-disk file size to an estimated in-memory (materialized,
# columnar) size. Parquet is compressed on disk so it inflates more when
# materialized; a TSV is ~1x on disk but string/typing overhead still inflates
# it in memory.
_MATERIALIZE_SIZE_FACTOR = {"parquet": 5.0, "tsv": 2.5}
# Materialize only when the estimate fits within this fraction of memory_limit,
# leaving headroom for query working memory (memory_limit is not a strict bound
# on process RSS in DuckDB).
_MATERIALIZE_LIMIT_FRACTION = 0.5
# When the memory_limit is unknown/unlimited, only materialize clearly small
# reports where memory pressure is a non-issue.
_MATERIALIZE_UNKNOWN_LIMIT_CAP = 512 * 1024 * 1024  # 512 MiB estimated

_SIZE_UNITS = {
    "B": 1,
    "KB": 1000,
    "MB": 1000**2,
    "GB": 1000**3,
    "TB": 1000**4,
    "KIB": 1024,
    "MIB": 1024**2,
    "GIB": 1024**3,
    "TIB": 1024**4,
}


def _parse_duckdb_size(value: str) -> int | None:
    """Parse a DuckDB size setting (e.g. ``'14.3 GiB'``, ``'8192MB'``) into bytes.

    Returns ``None`` for an unset/unlimited value (``'-1'`` / empty) or an
    unrecognized format, so callers treat it as "budget unknown".
    """
    text = value.strip()
    if not text or text == "-1":
        return None
    match = re.fullmatch(r"([\d.]+)\s*([KMGT]i?B|B)?", text, re.IGNORECASE)
    if not match:
        return None
    return int(float(match.group(1)) * _SIZE_UNITS.get((match.group(2) or "B").upper(), 1))


# DIA-NN channel name and mass offset of each canonical plexDIA label, from the
# sdrf-pipelines registry that builds DIA-NN's --channels (MTRAQ4 -> "4", SILAC
# heavy -> "H"), and the mass offset of each channel name.
_PLEXDIA_CHANNELS = {
    label: (channel["channel_name"], channel["masses"][0])
    for registry in PLEXDIA_REGISTRY.values()
    for label, channel in registry["channels"].items()
}
_CHANNEL_MASSES = dict(_PLEXDIA_CHANNELS.values())

# Raw-file extensions stripped from a report ``Run`` to give ``run_file_name``.
RUN_EXTENSIONS = "mzML|raw|d|wiff|htrms"


def run_name_sql(run_expr: str) -> str:
    """SQL for a report ``Run`` value without its raw-file extension."""
    return f"regexp_replace(CAST({run_expr} AS VARCHAR), '(?i)\\.({RUN_EXTENSIONS})$', '')"


def decoy_protein_group_sql(pg_expr: str, decoy_expr: str) -> str:
    """SQL for a report row's ``;``-separated protein group, ``DECOY_``-prefixed on decoy rows.

    DIA-NN ``--report-decoys`` reports a decoy precursor under its target's
    protein group and flags it only in ``Decoy``. Prefixing every accession gives
    the decoy rows their own pg group and links their features to it, not to the
    target group.
    """
    prefixed = (
        f"array_to_string(list_transform(list_filter(string_split({pg_expr}, ';'), lambda x: trim(x) <> ''), "
        "lambda x: 'DECOY_' || trim(x)), ';')"
    )
    return f"CASE WHEN CAST({decoy_expr} AS BOOLEAN) THEN {prefixed} ELSE {pg_expr} END"


def _populated_channel_values(raw_channels: list[object], channel_col: str) -> list[object]:
    """Return non-empty channel values, rejecting partially missing labels."""
    populated = [channel for channel in raw_channels if channel is not None and str(channel).strip()]
    if populated and len(populated) != len(raw_channels):
        raise ValueError(f"DIA-NN {channel_col} contains an empty channel identifier")
    return populated


def _declared_channel_label(raw_text: str, normalized: str, declared_by_key: dict[str, str]) -> str:
    """The SDRF label of a report channel: the label itself, or the plexDIA label DIA-NN names by that channel."""
    key = normalized.casefold()
    if key in declared_by_key:
        return declared_by_key[key]
    aliases = [
        label
        for label in declared_by_key.values()
        if label in _PLEXDIA_CHANNELS and _PLEXDIA_CHANNELS[label][0].casefold() == key
    ]
    if len(aliases) == 1:
        return aliases[0]
    raise ValueError(
        f"DIA-NN channel {raw_text!r} is not declared in SDRF comment[label]"
        + (f" (it names {', '.join(sorted(aliases))})" if aliases else "")
        + "; channel identifiers must match so intensities can resolve to run samples"
    )


def _channel_label_rows(raw_channels: list[object], declared_labels: set[str] | None) -> list[tuple[str, str, float | None]]:
    """Resolve report channels to distinct canonical SDRF labels and their mass offsets.

    A channel matches an SDRF label of the same name, or the declared plexDIA
    label DIA-NN calls by that channel name (``MTRAQ4`` for ``4``, ``SILAC heavy``
    for ``H``). The mass offset orders the partners of a multiplet; it is null
    for a channel outside the plexDIA labels.
    """
    declared_by_key = {normalize_label(label).casefold(): normalize_label(label) for label in declared_labels or ()}
    rows: list[tuple[str, str, float | None]] = []
    canonical_to_raw: dict[str, str] = {}
    for raw_channel in raw_channels:
        raw_text = str(raw_channel).strip()
        normalized = normalize_label(raw_text)
        if declared_by_key:
            normalized = _declared_channel_label(raw_text, normalized, declared_by_key)

        canonical_key = normalized.casefold()
        previous = canonical_to_raw.get(canonical_key)
        if previous is not None and previous != raw_text:
            raise ValueError(f"DIA-NN channels {previous!r} and {raw_text!r} both canonicalize to {normalized!r}")
        canonical_to_raw[canonical_key] = raw_text
        mass = _PLEXDIA_CHANNELS[normalized][1] if normalized in _PLEXDIA_CHANNELS else _CHANNEL_MASSES.get(raw_text)
        rows.append((raw_text, normalized, mass))
    return rows


class DiaNNBaseAdapter(BaseConverter):
    """Base adapter with DIA-NN-specific loading utilities.

    Subclasses inherit ``_load_diann_report()`` so the logic is defined
    once rather than duplicated across the feature and PG adapters.
    """

    def _load_diann_report(self, path: str) -> None:
        """Load a DIA-NN report into DuckDB, choosing TABLE vs VIEW by memory.

        Materializing the report as a ``TABLE`` parses the source once so every
        downstream scan is a cheap columnar read; a ``VIEW`` instead re-parses
        the file on each scan (~hundreds of x slower per scan on a mid-size
        report). But a ``TABLE`` holds the report in memory, which is unsafe for
        the multi-GB reports DIA-NN can produce (median public report ~570 MB,
        some > 10 GB). So the report is materialized only when its estimated
        in-memory size fits safely within DuckDB's configured ``memory_limit``,
        falls back to a ``VIEW`` otherwise, and a materialization that hits the
        limit is rolled back to a ``VIEW``. ``memory_limit`` is not a strict RSS
        bound, so unknown sizes take the conservative (VIEW) path.

        Supports both TSV (``report.tsv``) and Parquet (``report.parquet``).
        If the ``report`` relation already exists it is skipped.
        """
        if self._table_exists("report"):
            self.logger.debug("report already loaded -- skipping reload")
            return
        safe_path = escape_path(path)
        if path.endswith(".parquet"):
            reader = sql_build("read_parquet('$path')", path=safe_path)
        else:
            reader = sql_build(
                "read_csv_auto('$path', delim='\t', header=true, auto_detect=true, null_padding=true)",
                path=safe_path,
            )

        materialized = False
        if self._should_materialize_report(path):
            try:
                self._conn.execute(sql_build("CREATE TABLE report AS SELECT * FROM $reader", reader=reader))
                materialized = True
            except duckdb.OutOfMemoryException:
                self._conn.execute("DROP TABLE IF EXISTS report")
                self.logger.warning("DIA-NN report exceeded memory_limit while materializing; falling back to a lazy VIEW")
        if not materialized:
            self._conn.execute(sql_build("CREATE VIEW report AS SELECT * FROM $reader", reader=reader))

        if materialized:
            count = self._conn.execute("SELECT COUNT(*) FROM report").fetchall()[0][0]
            self.logger.info("DIA-NN report table created (%s rows)", format(count, ","))
        else:
            self.logger.info("DIA-NN report view created")

    def _should_materialize_report(self, path: str) -> bool:
        """Whether the report likely fits safely within DuckDB's ``memory_limit``.

        Conservative by design: estimates the in-memory size from the on-disk
        file size, requires it to fit within a fraction of the limit, and
        declines for compressed or unknown-size inputs, or an unknown/unlimited
        limit above a small absolute cap. A gzip file's compressed size cannot
        safely bound the memory needed for its expanded rows.
        """
        if path.casefold().endswith(".gz"):
            return False
        try:
            file_bytes = os.path.getsize(path)
        except OSError:
            return False  # unknown size -> conservative VIEW
        factor = _MATERIALIZE_SIZE_FACTOR["parquet" if path.endswith(".parquet") else "tsv"]
        estimated = file_bytes * factor
        limit = self._duckdb_memory_limit_bytes()
        if limit is None:
            return estimated <= _MATERIALIZE_UNKNOWN_LIMIT_CAP
        return estimated <= _MATERIALIZE_LIMIT_FRACTION * limit

    def _duckdb_memory_limit_bytes(self) -> int | None:
        """DuckDB's configured ``memory_limit`` in bytes, or ``None`` if unknown."""
        try:
            raw = self._conn.execute("SELECT current_setting('memory_limit')").fetchall()[0][0]
        except duckdb.Error:
            return None
        return _parse_duckdb_size(str(raw))

    def _register_channel_labels(
        self,
        report_columns: set[str],
        sdrf_path: str | None,
    ) -> str | None:
        """Register canonical DIA-NN channel labels and return the source column."""
        channel_col = next((col for col in ("Channel", "Label") if col in report_columns), None)
        if channel_col is None:
            return None

        quoted_channel = validate_identifier(channel_col)
        raw_channels = [
            row[0]
            for row in self._conn.execute(
                sql_build(
                    """
                    SELECT DISTINCT CAST($channel AS VARCHAR)
                    FROM report
                    ORDER BY CAST($channel AS VARCHAR)
                    """,
                    channel=quoted_channel,
                )
            ).fetchall()
        ]
        raw_channels = _populated_channel_values(raw_channels, channel_col)
        if not raw_channels:
            self.logger.info("DIA-NN %s is empty; treating the report as label-free", channel_col)
            return None

        declared_labels = read_sdrf_labels(sdrf_path)
        rows = _channel_label_rows(raw_channels, declared_labels)

        self._conn.execute("DROP TABLE IF EXISTS diann_channel_labels")
        self._conn.from_df(
            pd.DataFrame(
                rows,
                columns=["channel_value", "canonical_label", "channel_mass"],
            ).astype({"channel_mass": "float64"})
        ).create("diann_channel_labels")
        return channel_col
