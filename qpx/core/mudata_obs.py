"""Run and sample metadata for MuData observations."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from qpx.core.files import run_file_stem
from qpx.core.sample_sql import UNIT_SEPARATOR
from qpx.core.sql import sql_build

if TYPE_CHECKING:
    from qpx.core.engine import DuckDBEngine


_RUN_OBS_SQL = """
                SELECT r.run_file_name,
                       r.run_accession,
                       rs.sample_accession,
                       rs.label,
                       rs.biological_replicate,
                       rs.technical_replicate,
                       $optional
                FROM run r, UNNEST(r.samples) AS _t(rs)
                """

# Columns of run.parquet that describe one raw file rather than its samples.
RUN_LEVEL_COLUMNS = ("run_accession", "fraction", "instrument")

# Sample-channel columns on which the member runs of a quantification unit can
# disagree, e.g. a producer that groups re-injections of a sample into one unit.
# Such a unit leaves them empty instead of taking one member's value.
_REPLICATE_COLUMNS = ("biological_replicate", "technical_replicate")


def _has_view(engine: DuckDBEngine, name: str) -> bool:
    """Return True if the dataset registered the view *name*."""
    rows = engine.execute("SELECT 1 FROM information_schema.tables WHERE table_name = ? LIMIT 1", [name]).fetchall()
    return len(rows) > 0


def _run_sample_lookups(df: pd.DataFrame) -> tuple:
    """Index run/sample records by run name, stem and label."""
    exact: dict[tuple[str, str], dict] = {}
    stem_exact: dict[tuple[str, str], dict] = {}
    first_by_run: dict[str, dict] = {}
    first_by_stem: dict[str, dict] = {}
    labels_by_run: dict[str, set[str]] = {}
    labels_by_stem: dict[str, set[str]] = {}
    for record in df.to_dict("records"):
        run_name = str(record["run_file_name"])
        stem = run_file_stem(run_name)
        label = str(record["label"])
        exact.setdefault((run_name, label), record)
        stem_exact.setdefault((stem, label), record)
        first_by_run.setdefault(run_name, record)
        first_by_stem.setdefault(stem, record)
        labels_by_run.setdefault(run_name, set()).add(label)
        labels_by_stem.setdefault(stem, set()).add(label)
    return exact, stem_exact, first_by_run, first_by_stem, labels_by_run, labels_by_stem


class _RunSampleLookup:
    """Resolve a ``(run, intensity label)`` pair to its ``run.samples`` record."""

    def __init__(self, df: pd.DataFrame) -> None:
        (
            self._exact,
            self._stem_exact,
            self._first_by_run,
            self._first_by_stem,
            self._labels_by_run,
            self._labels_by_stem,
        ) = _run_sample_lookups(df)

    def resolve(self, run_name: str, label: str) -> dict | None:
        """Return the sample record of *label* in *run_name*, or None."""
        stem = run_file_stem(run_name)
        record = self._exact.get((run_name, label)) or self._stem_exact.get((stem, label))
        # Fall back to the run's sole sample when the intensity label does not
        # match any SDRF sample label. This is the label-free case (intensity
        # label "LFQ" vs sample label "label free sample"), detected by the run
        # carrying at most one distinct SDRF label. For multiplexed runs (more
        # than one distinct label / channel) we never guess, so each channel
        # keeps its own real sample_accession.
        if record is None and len(self._labels_by_run.get(run_name, ())) <= 1:
            record = self._first_by_run.get(run_name)
        if record is None and len(self._labels_by_stem.get(stem, ())) <= 1:
            record = self._first_by_stem.get(stem)
        return record

    def run_record(self, run_name: str) -> dict:
        """Return any record of *run_name* (for its run-level columns), or ``{}``."""
        return self._first_by_run.get(run_name) or self._first_by_stem.get(run_file_stem(run_name)) or {}


def _run_sample_records(engine: DuckDBEngine) -> pd.DataFrame:
    """Return one row per run sample channel, or an empty frame without a run view."""
    if not _has_view(engine, "run"):
        return pd.DataFrame()
    # fraction and instrument are optional run columns: an absent one reads as null.
    optional = engine.optional_columns("run", "r", ("fraction", "instrument"))
    return engine.execute(sql_build(_RUN_OBS_SQL, optional=optional)).fetchdf()


def _without_run_lookup(observation_keys: pd.DataFrame, all_labels: bool) -> pd.DataFrame:
    """Observations of a dataset without run metadata: just their keys."""
    if all_labels:
        return observation_keys.copy()
    return pd.DataFrame(index=observation_keys.index)


def _finish_obs(rows: list[dict], index: pd.Index) -> pd.DataFrame:
    """Build obs, filling NaN in string columns to avoid h5py serialization errors."""
    obs = pd.DataFrame(rows, index=index)
    str_cols = obs.select_dtypes(include=["object"]).columns
    obs[str_cols] = obs[str_cols].fillna("")
    return obs


def _load_run_obs(
    engine: DuckDBEngine,
    observation_keys: pd.DataFrame,
    all_labels: bool,
) -> pd.DataFrame:
    """Load run/sample metadata for matrix observations."""
    df = _run_sample_records(engine)
    if df.empty:
        return _without_run_lookup(observation_keys, all_labels)

    lookup = _RunSampleLookup(df)
    rows: list[dict] = []
    for run_name, label in zip(observation_keys["run_file_name"], observation_keys["intensity_label"]):
        run_name, label = str(run_name), str(label)
        values = dict(lookup.resolve(run_name, label) or {})
        values.pop("run_file_name", None)
        if all_labels:
            values = {"run_file_name": run_name, "intensity_label": label, **values}
        rows.append(values)
    return _finish_obs(rows, observation_keys.index)


def _joined_run_columns(lookup: _RunSampleLookup, runs: list[str]) -> dict[str, str]:
    """Join each run-level column over the member runs of a quantification unit."""
    records = [lookup.run_record(run) for run in runs]
    joined = {}
    for column in RUN_LEVEL_COLUMNS:
        values = (record.get(column) for record in records)
        joined[column] = UNIT_SEPARATOR.join(dict.fromkeys(str(v) for v in values if v is not None and not pd.isna(v)))
    return joined


def _unit_sample_record(lookup: _RunSampleLookup, runs: list[str], label: str) -> dict:
    """Return the sample record of *label* in a unit's runs, without disputed replicate numbers."""
    records = [record for record in (lookup.resolve(run, label) for run in runs) if record]
    if not records:
        return {}
    values = dict(records[0])
    for column in _REPLICATE_COLUMNS:
        if len({record.get(column) for record in records}) > 1:
            values[column] = float("nan")
    return values


def _load_unit_obs(
    engine: DuckDBEngine,
    unit_keys: pd.DataFrame,
    all_labels: bool,
) -> pd.DataFrame:
    """Load sample metadata for observations that are quantification units.

    ``unit_keys`` has ``grouped_runs`` (the unit's member runs joined by
    :data:`UNIT_SEPARATOR`) and ``intensity_label``. Member runs are fractions of
    one sample context, so the sample is resolved through the first member with
    a record for the label (docs/spec/data-model.md); ``biological_replicate``
    and ``technical_replicate`` are left empty when the members disagree.
    Run-level columns (``run_accession``, ``fraction``, ``instrument``) describe
    single raw files and are joined over the members.
    """
    df = _run_sample_records(engine)
    if df.empty:
        return _without_run_lookup(unit_keys, all_labels)

    lookup = _RunSampleLookup(df)
    rows: list[dict] = []
    for unit, label in zip(unit_keys["grouped_runs"], unit_keys["intensity_label"]):
        unit, label = str(unit), str(label)
        runs = unit.split(UNIT_SEPARATOR)
        values = _unit_sample_record(lookup, runs, label)
        values.pop("run_file_name", None)
        values.update(_joined_run_columns(lookup, runs))
        if all_labels:
            values = {"grouped_runs": unit, "intensity_label": label, **values}
        rows.append(values)
    return _finish_obs(rows, unit_keys.index)
