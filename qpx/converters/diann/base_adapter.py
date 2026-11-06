"""Shared base adapter for DIA-NN converters.

Provides common data-loading helpers used by both
:class:`DiannFeatureAdapter` and :class:`DiannPgAdapter`.
"""

from __future__ import annotations

import os
import tempfile
from typing import Iterator

import pandas as pd
from sdrf_pipelines.converters.diann.constants import PLEXDIA_REGISTRY

from qpx.converters.base import BaseConverter
from qpx.converters.channel_labels import normalize_label, read_sdrf_labels
from qpx.core.files import RUN_FILE_SUFFIX_PATTERN
from qpx.core.sql import escape_path, sql_build, validate_identifier

# DIA-NN channel name and mass offset of each canonical plexDIA label, from the
# sdrf-pipelines registry that builds DIA-NN's --channels (MTRAQ4 -> "4", SILAC
# heavy -> "H"), and the mass offset of each channel name.
_PLEXDIA_CHANNELS = {
    label: (channel["channel_name"], channel["masses"][0])
    for registry in PLEXDIA_REGISTRY.values()
    for label, channel in registry["channels"].items()
}
_CHANNEL_MASSES = dict(_PLEXDIA_CHANNELS.values())

# Report rows per conversion batch of runs. The feature step of a 5.4 M-row
# DIA-NN 2.5 report peaks at 4.3 GB with it, against 6.9 GB in one batch.
BATCH_ROWS = 1_000_000

# Raw-file extensions stripped from a report ``Run`` to give ``run_file_name``.


def run_name_sql(run_expr: str) -> str:
    """SQL for a report ``Run`` value without its raw-file extension."""
    return f"regexp_replace(CAST({run_expr} AS VARCHAR), '{RUN_FILE_SUFFIX_PATTERN}', '')"


def decoy_protein_group_sql(pg_expr: str, decoy_expr: str) -> str:
    """SQL for a report row's ``;``-separated protein group, ``DECOY_``-prefixed on decoy rows.

    DIA-NN ``--report-decoys`` reports a decoy precursor under its target's
    protein group and flags it only in ``Decoy``. Prefixing every accession gives
    the decoy rows their own pg group and links their features to it, not to the
    target group.
    """
    prefixed = (
        f"array_to_string(list_transform(list_filter(string_split({pg_expr}, ';'), lambda x: trim(x) <> ''), "
        "lambda x: 'DECOY_' || regexp_replace(trim(x), '^(?:sp|tr)\\|([^|]+)(?:\\|.*)?$', '\\1')), ';')"
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

    # Temporary Parquet copy of a TSV report, removed on close.
    _report_copy: str | None = None

    def _load_diann_report(self, path: str) -> None:
        """Expose a DIA-NN report to DuckDB as the ``report`` view.

        A Parquet report is read in place: DuckDB scans only the columns and row
        groups a query needs, as fast as a materialized table and in less memory.
        A TSV (or gzipped TSV) report is first converted once to a temporary
        Parquet file in DuckDB's ``temp_directory``, since a view over the text
        would re-parse it for every batch of runs. The copy is removed when the
        adapter closes. If the ``report`` relation already exists it is skipped.
        """
        if self._table_exists("report"):
            self.logger.debug("report already loaded -- skipping reload")
            return
        if not path.endswith(".parquet"):
            temp_dir = self._conn.execute("SELECT current_setting('temp_directory')").fetchone()[0]
            os.makedirs(temp_dir, exist_ok=True)
            handle, self._report_copy = tempfile.mkstemp(prefix="qpx_diann_report_", suffix=".parquet", dir=temp_dir)
            os.close(handle)
            self._conn.execute(
                sql_build(
                    "COPY (SELECT * FROM read_csv_auto('$path', delim='\t', header=true, auto_detect=true, "
                    "null_padding=true)) TO '$copy' (FORMAT parquet)",
                    path=escape_path(path),
                    copy=escape_path(self._report_copy),
                )
            )
            self.logger.info("DIA-NN report copied to %s", self._report_copy)
            path = self._report_copy
        self._conn.execute(sql_build("CREATE VIEW report AS SELECT * FROM read_parquet('$path')", path=escape_path(path)))
        self.logger.info("DIA-NN report view created")

    def close(self):
        """Close the connection and remove the temporary copy of a TSV report."""
        super().close()
        if self._report_copy:
            os.remove(self._report_copy)
            self._report_copy = None

    def _run_batches(self, run_col: str, runs: list[str], batch_rows: int) -> Iterator[list[str]]:
        """Split *runs*, values of the report column *run_col*, into consecutive batches of about *batch_rows* rows.

        A batch's memory follows its row count, which a fixed number of runs does
        not bound: DIA-NN runs range from thousands to hundreds of thousands of
        precursors.
        """
        counts = dict(self._conn.execute(sql_build("SELECT $run, count(*) FROM report GROUP BY ALL", run=run_col)).fetchall())
        batch: list[str] = []
        rows = 0
        for run in runs:
            batch.append(run)
            rows += counts.get(run, 0)
            if rows >= batch_rows:
                yield batch
                batch, rows = [], 0
        if batch:
            yield batch

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
