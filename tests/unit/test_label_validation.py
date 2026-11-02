"""Cross-view checks beside run and label resolution: PSM labels, grouped-run fraction order and feature-to-pg coverage."""

import pytest

from qpx.core.data.identity import derive_id
from qpx.core.label_validation import pg_link_coverage_issues
from qpx.dataset import Dataset
from qpx.writers import DatasetWriter, FeatureWriter, PgWriter, PsmWriter, RunWriter
from tests.conftest import make_dataset_record, make_feature_record, make_pg_record, make_psm_record, make_run_record


def _write_dataset(ds_dir, *, features, run_labels, psms=()):
    """Write feature (+ psm) records with one run_01 whose samples carry *run_labels*."""
    ds_dir.mkdir(parents=True, exist_ok=True)
    run = make_run_record()
    run["samples"] = [
        {"sample_accession": accession, "label": label, "biological_replicate": 1, "technical_replicate": 1}
        for accession, label in run_labels
    ]
    with RunWriter(ds_dir / "exp.run.parquet") as w:
        w.write_batch([run])
    with DatasetWriter(ds_dir / "exp.dataset.parquet") as w:
        w.write_batch([make_dataset_record()])
    with FeatureWriter(ds_dir / "exp.feature.parquet") as w:
        w.write_batch(list(features))
    if psms:
        with PsmWriter(ds_dir / "exp.psm.parquet") as w:
            w.write_batch(list(psms))
    return ds_dir


def _issues(ds_dir, structure, check, *, strict=False):
    with Dataset(ds_dir) as ds:
        results = ds.validate(strict=strict)
    return [issue for issue in results[structure].issues if issue.check == check]


TMT_RUN = [("S1", "TMT126"), ("S2", "TMT127N")]


def test_resolving_labels_raise_no_issue(tmp_path):
    """Labels that each name one sample of a known run pass every run and label check."""
    feature = make_feature_record(intensities=[{"label": "TMT126", "intensity": 1.0}, {"label": "TMT127N", "intensity": 2.0}])
    ds_dir = _write_dataset(tmp_path, features=[feature], run_labels=TMT_RUN)
    with Dataset(ds_dir) as ds:
        results = ds.validate(strict=True)
    assert results["feature"].issues == []


def test_psm_label_outside_identified_labels(tmp_path):
    """A PSM's label must be one of the labels its feature was identified in."""
    feature = make_feature_record(
        intensities=[{"label": "SILAC light", "intensity": 1.0}, {"label": "SILAC heavy", "intensity": 2.0}]
    )
    feature["identified_labels"] = ["SILAC light"]
    psm = make_psm_record()
    psm["label"] = "SILAC heavy"
    psm["feature_id"] = derive_id([feature["peptidoform"], feature["charge"], feature["run_file_name"], feature["rt"]])
    ds_dir = _write_dataset(
        tmp_path,
        features=[feature],
        run_labels=[("S1", "SILAC light"), ("S2", "SILAC heavy")],
        psms=[psm],
    )
    [issue] = _issues(ds_dir, "psm", "label_not_identified")
    assert issue.severity == "warning"
    [strict_issue] = _issues(ds_dir, "psm", "label_not_identified", strict=True)
    assert strict_issue.severity == "error"


@pytest.mark.parametrize(("grouped_runs", "flagged"), [(["run_02", "run_10"], False), (["run_10", "run_02"], True)])
def test_grouped_runs_follow_the_fraction_order(tmp_path, grouped_runs, flagged):
    """grouped_runs is stored in fraction order, compared as numbers (fraction 10 after fraction 2)."""
    runs = []
    for name, fraction in (("run_02", "2"), ("run_10", "10")):
        run = make_run_record(run_accession=name, run_file_name=name)
        run["fraction"] = fraction
        runs.append(run)
    pg = make_pg_record()
    pg["grouped_runs"] = grouped_runs
    with RunWriter(tmp_path / "exp.run.parquet") as w:
        w.write_batch(runs)
    with PgWriter(tmp_path / "exp.pg.parquet") as w:
        w.write_batch([pg])
    issues = _issues(tmp_path, "pg", "grouped_runs_fraction_order")
    assert [issue.severity for issue in issues] == (["warning"] if flagged else [])


def test_low_feature_to_pg_coverage_warns():
    """Two producers naming groups differently leave most features without a pg row."""
    issues = pg_link_coverage_issues({"n_features_linked": 90, "n_features_without_pg": 10})
    assert [(issue.check, issue.severity) for issue in issues] == [("pg_link_coverage", "warning")]
    assert "90.0%" in issues[0].message
    assert pg_link_coverage_issues({"n_features_linked": 95, "n_features_without_pg": 5}) == []
    assert pg_link_coverage_issues({"n_features_linked": 0, "n_features_without_pg": 0}) == []


def test_semantic_validation_reports_features_without_a_pg_row(tmp_path):
    ds_dir = _write_dataset(tmp_path, features=[make_feature_record()], run_labels=TMT_RUN)
    with PgWriter(ds_dir / "exp.pg.parquet") as w:
        w.write_batch([make_pg_record(anchor_protein="Q99999", pg_accessions=["Q99999"])])

    with Dataset(ds_dir) as ds:
        issues = ds.validate(semantic=True)["feature"].issues

    assert "pg_link_coverage" in {issue.check for issue in issues}
