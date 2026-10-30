"""Dataset validation of run references and label-to-sample resolution across views."""

import pytest

from qpx.core.data.identity import derive_id
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


def test_feature_naming_an_unknown_run_is_an_error(tmp_path):
    """Orphaned feature runs resolve to no sample, as in FragPipe output with renamed runs."""
    features = [make_feature_record(run_file_name="NOT_A_RUN", charge=charge) for charge in (1, 2, 3)]
    ds_dir = _write_dataset(tmp_path, features=features, run_labels=TMT_RUN)
    [issue] = _issues(ds_dir, "feature", "dangling_run")
    assert issue.severity == "error"
    assert "3 feature row(s)" in issue.message
    assert "'NOT_A_RUN'" in issue.message


def test_psm_naming_an_unknown_run_is_an_error(tmp_path):
    """The run reference check covers PSMs as well."""
    ds_dir = _write_dataset(
        tmp_path,
        features=[make_feature_record()],
        run_labels=TMT_RUN,
        psms=[make_psm_record(run_file_name="NOT_A_RUN")],
    )
    [issue] = _issues(ds_dir, "psm", "dangling_run")
    assert issue.severity == "error"


def test_repeated_label_in_one_row_is_an_error(tmp_path):
    """A label listed twice in one row would be counted twice by every consumer."""
    feature = make_feature_record(intensities=[{"label": "TMT126", "intensity": 1.0}] * 3)
    ds_dir = _write_dataset(tmp_path, features=[feature], run_labels=TMT_RUN)
    [issue] = _issues(ds_dir, "feature", "duplicate_label")
    assert issue.severity == "error"


@pytest.mark.parametrize(
    ("run_labels", "label", "check", "severity"),
    [
        (TMT_RUN, "TMT128N", "unresolved_label", "error"),
        ([("S1", "TMT126"), ("S2", "TMT126")], "TMT126", "ambiguous_label_sample", "error"),
        ([("S1", "LFQ")], "raw", "label_single_sample_fallback", "warning"),
    ],
)
def test_label_must_resolve_to_one_sample(tmp_path, run_labels, label, check, severity):
    """A (run, label) pair resolves to exactly one sample; a single-sample run only falls back."""
    feature = make_feature_record(intensities=[{"label": label, "intensity": 1.0}])
    ds_dir = _write_dataset(tmp_path, features=[feature], run_labels=run_labels)
    [issue] = _issues(ds_dir, "feature", check)
    assert issue.severity == severity
    assert f"('run_01', '{label}')" in issue.message


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
