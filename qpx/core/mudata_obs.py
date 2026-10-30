"""Run and sample metadata for MuData observations."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from qpx.core.files import run_file_stem
from qpx.core.sql import sql_build, validate_identifier

if TYPE_CHECKING:
    from qpx.core.engine import DuckDBEngine


_RUN_OBS_SQL = """
                SELECT r.run_file_name,
                       r.run_accession,
                       rs.sample_accession,
                       rs.label,
                       rs.biological_replicate,
                       rs.technical_replicate,
                       $optional$characteristics
                FROM run r
                CROSS JOIN UNNEST(r.samples) AS _t(rs)
                $sample_join
                """


def _sample_characteristics(engine: DuckDBEngine) -> tuple[str, str]:
    """Select list and join adding the sample.parquet characteristics, when the dataset has one."""
    if not _has_view(engine, "sample"):
        return "", ""
    columns = [
        name
        for name, column_type, *_ in engine.execute("DESCRIBE sample").fetchall()
        if column_type == "VARCHAR" and name != "sample_accession"
    ]
    select = "".join(f", s.{validate_identifier(name)}" for name in columns)
    return select, "LEFT JOIN sample s ON s.sample_accession = rs.sample_accession"


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


def _load_run_obs(
    engine: DuckDBEngine,
    observation_keys: pd.DataFrame,
    all_labels: bool,
) -> pd.DataFrame:
    """Load run/sample metadata for matrix observations."""
    if _has_view(engine, "run"):
        # fraction and instrument are optional run columns: an absent one reads as null.
        optional = engine.optional_columns("run", "r", ("fraction", "instrument"))
        characteristics, sample_join = _sample_characteristics(engine)
        df = engine.execute(
            sql_build(
                _RUN_OBS_SQL,
                optional=optional,
                characteristics=characteristics,
                sample_join=sample_join,
            )
        ).fetchdf()
    else:
        df = pd.DataFrame()

    if df.empty:
        if all_labels:
            return observation_keys.copy()
        return pd.DataFrame(index=observation_keys.index)

    exact, stem_exact, first_by_run, first_by_stem, labels_by_run, labels_by_stem = _run_sample_lookups(df)

    rows: list[dict] = []
    for _, key in observation_keys.iterrows():
        run_name = str(key["run_file_name"])
        stem = run_file_stem(run_name)
        label = str(key["intensity_label"])
        record = exact.get((run_name, label))
        if record is None:
            record = stem_exact.get((stem, label))
        # Fall back to the run's sole sample when the intensity label does not
        # match any SDRF sample label. This is the label-free case (intensity
        # label "LFQ" vs sample label "label free sample"), detected by the run
        # carrying at most one distinct SDRF label. For multiplexed runs (more
        # than one distinct label / channel) we never guess, so each channel
        # keeps its own real sample_accession.
        if record is None and len(labels_by_run.get(run_name, ())) <= 1:
            record = first_by_run.get(run_name)
        if record is None and len(labels_by_stem.get(stem, ())) <= 1:
            record = first_by_stem.get(stem)

        values = dict(record) if record is not None else {}
        values.pop("run_file_name", None)
        if all_labels:
            values = {
                "run_file_name": run_name,
                "intensity_label": label,
                **values,
            }
        rows.append(values)

    obs = pd.DataFrame(rows, index=observation_keys.index)

    # Fill remaining NaN in string columns to avoid h5py serialization errors
    str_cols = obs.select_dtypes(include=["object"]).columns
    obs[str_cols] = obs[str_cols].fillna("")
    return obs
