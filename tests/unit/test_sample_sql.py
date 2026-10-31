"""Sample resolution SQL of the Dataset helpers and API views.

The helpers resolve feature and pg quantities to samples through
``run.samples[]``. Unnesting ``run.samples`` and ``feature.intensities`` inside
one ``FROM`` list with ``feature`` and ``run`` planned a feature x run cross
product, and the views joined every channel of a TMT feature to every channel
of its run before filtering. These tests pin the results against an
independent reference and bound the plan without timing anything.
"""

from __future__ import annotations

import json

import pandas as pd
import pytest

from qpx.dataset import Dataset
from qpx.writers import FeatureWriter, PgWriter, RunWriter
from tests.conftest import make_feature_record, make_pg_record, make_run_record

TMT16 = [
    "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N",
    "TMT130C", "TMT131N", "TMT131C", "TMT132N", "TMT132C", "TMT133N", "TMT133C", "TMT134N",
]  # fmt: skip
PLEXES = ("plexA", "plexB")
FRACTIONS = (1, 2)
FEATURES_PER_RUN = 6
LFQ_RUNS = 120


def _tmt_runs():
    """Two TMT16 plexes, each measured in two fractions (four runs)."""
    runs = []
    for plex in PLEXES:
        for fraction in FRACTIONS:
            record = make_run_record(run_accession=f"{plex}_F{fraction}", run_file_name=f"{plex}_F{fraction}")
            record["fraction"] = str(fraction)
            record["samples"] = [
                {
                    "sample_accession": f"{plex}_{label}",
                    "label": label,
                    "biological_replicate": 1,
                    "technical_replicate": 1,
                }
                for label in TMT16
            ]
            runs.append(record)
    return runs


def _tmt_features():
    features = []
    for run_index, run in enumerate(_tmt_runs()):
        for k in range(FEATURES_PER_RUN):
            intensities = [
                {"label": label, "intensity": float(1000 * run_index + 10 * k + channel)} for channel, label in enumerate(TMT16)
            ]
            features.append(
                make_feature_record(
                    sequence=f"PEPTIDE{k}K",
                    peptidoform=f"PEPTIDE{k}K",
                    charge=2 + k % 2,
                    is_decoy=k == FEATURES_PER_RUN - 1,
                    run_file_name=run["run_file_name"],
                    anchor_protein=f"P{k % 3}",
                    scan=[100 + k],
                    intensities=intensities,
                )
            )
    return features


def _tmt_pg():
    rows = []
    for plex in PLEXES:
        for k in range(3):
            record = make_pg_record(
                anchor_protein=f"P{k}",
                is_decoy=k == 2,
                intensities=[{"label": label, "intensity": float(100 * k + c)} for c, label in enumerate(TMT16)],
            )
            record["grouped_runs"] = [f"{plex}_F{fraction}" for fraction in FRACTIONS]
            rows.append(record)
    return rows


@pytest.fixture(name="tmt_dataset")
def _tmt_dataset(tmp_path):
    """Fractionated TMT16 dataset written through the qpx writers."""
    with RunWriter(tmp_path / "t.run.parquet") as writer:
        writer.write_batch(_tmt_runs())
    with FeatureWriter(tmp_path / "t.feature.parquet") as writer:
        writer.write_batch(_tmt_features())
    with PgWriter(tmp_path / "t.pg.parquet") as writer:
        writer.write_batch(_tmt_pg())
    return tmp_path


@pytest.fixture(name="lfq_dataset")
def _lfq_dataset(tmp_path):
    """120 label-free runs whose features carry DIA-NN's ``raw`` label, not ``LFQ``."""
    runs, features = [], []
    for i in range(LFQ_RUNS):
        run = make_run_record(run_accession=f"assay_{i:03d}", run_file_name=f"run_{i:03d}")
        run["samples"][0].update(sample_accession=f"S{i:03d}", label="LFQ")
        runs.append(run)
        for k in range(2):
            features.append(
                make_feature_record(
                    sequence=f"PEP{k}K",
                    peptidoform=f"PEP{k}K",
                    run_file_name=run["run_file_name"],
                    scan=[k],
                    intensities=[{"label": "raw", "intensity": float(i * 10 + k)}],
                )
            )
    with RunWriter(tmp_path / "l.run.parquet") as writer:
        writer.write_batch(runs)
    with FeatureWriter(tmp_path / "l.feature.parquet") as writer:
        writer.write_batch(features)
    return tmp_path


def _reference_feature_samples(features, runs) -> pd.DataFrame:
    """One row per (feature, intensity, sample): the documented resolution rule."""
    by_run = {run["run_file_name"]: run["samples"] for run in runs}
    rows = []
    for feature in features:
        samples = by_run[feature["run_file_name"]]
        for entry in feature["intensities"]:
            for sample in samples:
                if entry["label"] == sample["label"] or len(samples) == 1:
                    rows.append(
                        {
                            "sequence": feature["sequence"],
                            "peptidoform": feature["peptidoform"],
                            "is_decoy": feature["is_decoy"],
                            "sample_accession": sample["sample_accession"],
                            "label": sample["label"],
                            "intensity": entry["intensity"],
                        }
                    )
    return pd.DataFrame(rows)


def _sorted(df: pd.DataFrame, columns) -> pd.DataFrame:
    return df[list(columns)].sort_values(list(columns)).reset_index(drop=True)


def test_feature_helpers_match_the_resolution_rule(tmt_dataset):
    reference = _reference_feature_samples(_tmt_features(), _tmt_runs())
    with Dataset(tmt_dataset) as ds:
        view = ds.peptide_view.intensity().to_df()
        quantified = ds.feature.for_quantification(ds).to_df()
        abundance = ds.intensity("peptide").to_df()

    columns = ("sequence", "sample_accession", "label", "intensity")
    pd.testing.assert_frame_equal(_sorted(view, columns), _sorted(reference, columns), check_dtype=False)
    pd.testing.assert_frame_equal(_sorted(quantified, columns), _sorted(reference, columns), check_dtype=False)
    assert set(quantified["fraction"]) == {"1", "2"}

    targets = reference[~reference["is_decoy"]]
    expected = targets.groupby(["sample_accession", "peptidoform"], as_index=False)["intensity"].sum()
    expected = expected.rename(columns={"peptidoform": "feature_id"})
    columns = ("sample_accession", "feature_id", "intensity")
    pd.testing.assert_frame_equal(_sorted(abundance, columns), _sorted(expected, columns), check_dtype=False)


def test_pg_helpers_resolve_each_unit_once_per_sample(tmt_dataset):
    """A pg row of a two-fraction unit is returned once for its channel's sample."""
    with Dataset(tmt_dataset) as ds:
        abundance = ds.intensity("protein").to_df()
        view = ds.protein_view.intensity().to_df()

    # 2 plexes x 3 groups x 16 channels; the decoy group drops out of intensity().
    assert len(view) == 2 * 3 * 16
    assert len(abundance) == 2 * 2 * 16
    row = view[(view["protein_accession"] == "P1") & (view["sample_accession"] == "plexB_TMT127N")]
    assert row["intensity"].tolist() == [101.0]
    assert row["label"].tolist() == ["TMT127N"]


def test_single_sample_runs_take_any_feature_label(lfq_dataset):
    with Dataset(lfq_dataset) as ds:
        abundance = ds.intensity("peptide").to_df()
        view = ds.peptide_view.intensity().to_df()

    assert len(abundance) == len(view) == LFQ_RUNS * 2
    assert set(view["label"]) == {"LFQ"}
    assert abundance.set_index(["sample_accession", "feature_id"]).loc[("S007", "PEP1K"), "intensity"] == 71.0


def _helper_sql(ds, name):
    return {
        "intensity_peptide": lambda: ds._abundance_sql("peptide"),
        "intensity_protein": lambda: ds._abundance_sql("protein"),
        "peptide_view": lambda: ds.peptide_view.intensity()._relation.sql_query(),
        "protein_view": lambda: ds.protein_view.intensity()._relation.sql_query(),
        "for_quantification": lambda: ds.feature.for_quantification(ds)._relation.sql_query(),
    }[name]()


def _operator_rows(node, out):
    """Collect every operator's output row count from a JSON profile (DuckDB 1.5 and 2.0 keys)."""
    if isinstance(node, dict):
        for key in ("operator_cardinality", "intermediate_rows"):
            if key in node:
                out.append(int(node[key]))
        for child in node.get("children", []) + node.get("operator", []):
            _operator_rows(child, out)
    elif isinstance(node, list):
        for child in node:
            _operator_rows(child, out)
    return out


def _operator_names(node, out):
    """Collect the operator names of a JSON plan."""
    if isinstance(node, dict):
        if "name" in node:
            out.append(node["name"])
        for child in node.get("children", []):
            _operator_names(child, out)
    elif isinstance(node, list):
        for child in node:
            _operator_names(child, out)
    return out


def _max_operator_rows(ds, sql) -> int:
    profile = ds._engine.execute("EXPLAIN (ANALYZE, FORMAT JSON) " + sql).fetchall()[0][1]
    rows = _operator_rows(json.loads(profile), [])
    assert rows, "the profile reported no operator row counts"
    return max(rows)


FEATURE_HELPERS = ("intensity_peptide", "peptide_view", "for_quantification")
PG_HELPERS = ("intensity_protein", "protein_view")


@pytest.mark.parametrize("name", FEATURE_HELPERS + PG_HELPERS)
def test_plans_have_no_cross_product(tmt_dataset, lfq_dataset, name):
    for path in (tmt_dataset, lfq_dataset):
        with Dataset(path) as ds:
            if name in PG_HELPERS and ds.pg is None:
                continue
            plan = json.loads(ds._engine.execute("EXPLAIN (FORMAT JSON) " + _helper_sql(ds, name)).fetchall()[0][1])
        assert "CROSS_PRODUCT" not in _operator_names(plan, [])


@pytest.mark.parametrize("name", FEATURE_HELPERS)
def test_feature_plans_never_multiply_channels(tmt_dataset, name):
    """No operator emits more rows than the feature intensities themselves.

    Joining a TMT16 feature to its run before matching the label produced
    16 x 16 rows per feature; the label is a join key now.
    """
    with Dataset(tmt_dataset) as ds:
        assert _max_operator_rows(ds, _helper_sql(ds, name)) <= len(_tmt_features()) * len(TMT16)


@pytest.mark.parametrize("name", PG_HELPERS)
def test_pg_plans_never_multiply_members(tmt_dataset, name):
    """No operator emits more rows than the pg view: units are resolved once, not per row."""
    n_pg_rows = len(_tmt_pg()) * len(TMT16)
    with Dataset(tmt_dataset) as ds:
        assert _max_operator_rows(ds, _helper_sql(ds, name)) <= n_pg_rows
