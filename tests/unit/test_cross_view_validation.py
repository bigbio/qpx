"""Cross-view checks of Dataset.validate() that resolve rows through the run view.

The pg checks look up each distinct grouped run, and resolve each distinct
(grouped runs, label) unit through ``run.samples``, once rather than once per pg
row. feature, psm and mz rows must name a run of the run view, and each feature
label must name exactly one sample of its run. The bounded-memory tests run the
checks on fractionated TMT data with 32 MB of DuckDB memory and no temp space,
where joining every row to its runs' channels runs out of memory.
"""

import contextlib

import duckdb
import pyarrow.parquet as pq
import pytest

from qpx.core.data.schema import ValidationResult
from qpx.core.dataset_validation import _GROUPED_RUN_SAMPLE_MAPPING_SQL, _skipped_check_reason
from qpx.dataset import Dataset
from qpx.writers import DatasetWriter, FeatureWriter, MzWriter, PgWriter, PsmWriter, RunWriter
from tests.conftest import make_dataset_record, make_feature_record, make_pg_record, make_psm_record, make_run_record

TMT = ("TMT126", "TMT127N", "TMT127C", "TMT128N")


def _run_name(plex: int, fraction: int) -> str:
    return f"plex{plex}_F{fraction:02d}"


def _tmt_runs(n_plex: int, n_frac: int) -> list[dict]:
    """Fractionated TMT design: every fraction of a plex carries the plex's channels."""
    runs = []
    for plex in range(n_plex):
        for fraction in range(n_frac):
            run = make_run_record(run_accession=_run_name(plex, fraction), run_file_name=_run_name(plex, fraction))
            run["fraction"] = str(fraction + 1)
            run["samples"] = [
                {"sample_accession": f"S{plex}_{label}", "label": label, "biological_replicate": 1, "technical_replicate": 1}
                for label in TMT
            ]
            runs.append(run)
    return runs


def _tmt_pg(n_plex: int, n_frac: int, n_groups: int) -> list[dict]:
    """One pg record per (protein group, plex); the writer flattens it to one row per label."""
    records = []
    for plex in range(n_plex):
        for group in range(n_groups):
            record = make_pg_record(anchor_protein=f"P{group:05d}")
            record["grouped_runs"] = [_run_name(plex, fraction) for fraction in range(n_frac)]
            record["intensities"] = [{"label": label, "intensity": 1000.0 + i} for i, label in enumerate(TMT)]
            records.append(record)
    return records


def _write_dataset(ds_dir, *, runs, pg):
    ds_dir.mkdir(parents=True, exist_ok=True)
    with RunWriter(ds_dir / "exp.run.parquet") as writer:
        writer.write_batch(runs)
    with PgWriter(ds_dir / "exp.pg.parquet") as writer:
        writer.write_batch(pg)
    with DatasetWriter(ds_dir / "exp.dataset.parquet") as writer:
        writer.write_batch([make_dataset_record()])
    return ds_dir


def _bad_pg(anchor: str, grouped_runs: list[str], label: str) -> dict:
    record = make_pg_record(anchor_protein=anchor)
    record["grouped_runs"] = grouped_runs
    record["intensities"] = [{"label": label, "intensity": 1.0}]
    return record


def _issues(result, check):
    return [issue for issue in result.issues if issue.check == check]


@contextlib.contextmanager
def _low_memory(ds_dir):
    """Open *ds_dir* with 32 MB of DuckDB memory, one thread and no temp space to spill to."""
    with Dataset(ds_dir, duckdb_memory="32MB", duckdb_threads=1) as ds:
        ds._engine.execute("SET max_temp_directory_size='0KiB'")
        yield ds


@pytest.fixture(scope="module")
def large_tmt_dir(tmp_path_factory):
    """2 plexes x 12 fractions x 4,000 protein groups: 32,000 pg rows in 8 (grouped runs, label) units."""
    return _write_dataset(tmp_path_factory.mktemp("large") / "tmt", runs=_tmt_runs(2, 12), pg=_tmt_pg(2, 12, 4000))


@pytest.mark.parametrize("strict", [False, True])
def test_fractionated_tmt_pg_resolves_every_label_to_one_sample(tmp_path, strict):
    """2 plexes x 12 fractions: each (fraction set, channel) unit names one sample of its plex."""
    ds_dir = _write_dataset(tmp_path / "tmt", runs=_tmt_runs(2, 12), pg=_tmt_pg(2, 12, 40))
    with Dataset(ds_dir) as ds:
        assert ds.pg.count() == 2 * 40 * len(TMT)
        pg = ds.validate(structures=["pg"], strict=strict)["pg"]

    assert not _issues(pg, "ambiguous_grouped_run_sample")
    assert not [issue for issue in pg.issues if issue.check.endswith("check_skipped")]
    assert pg.is_valid


@pytest.mark.parametrize("strict", [False, True])
def test_ambiguous_pg_sample_mapping_is_reported_once_per_unit(tmp_path, strict):
    """pg rows sharing an ambiguous unit give one issue that counts them.

    The unit is the set of grouped runs, so a row listing the same runs in
    another order, or one of them twice, belongs to it.
    """
    spanning = [_run_name(0, 0), _run_name(1, 0)]
    bad = [
        _bad_pg("BAD1", spanning, "TMT126"),
        _bad_pg("BAD2", spanning[::-1], "TMT126"),
        _bad_pg("BAD3", [*spanning, spanning[0]], "TMT126"),
    ]
    ds_dir = _write_dataset(tmp_path / "ambiguous", runs=_tmt_runs(2, 3), pg=_tmt_pg(2, 3, 5) + bad)
    with Dataset(ds_dir) as ds:
        pg = ds.validate(structures=["pg"], strict=strict)["pg"]

    ambiguous = _issues(pg, "ambiguous_grouped_run_sample")
    assert len(ambiguous) == 1
    assert ambiguous[0].severity == "error"
    assert ambiguous[0].message.startswith("3 pg row(s) with label 'TMT126'")
    assert "2 samples (['S0_TMT126', 'S1_TMT126'])" in ambiguous[0].message
    assert "e.g. PG 'BAD1'" in ambiguous[0].message
    assert not pg.is_valid


@pytest.mark.parametrize("strict", [False, True])
def test_pg_label_without_a_sample_is_an_error(tmp_path, strict):
    """A label no run of the unit carries resolves to no sample."""
    bad = [_bad_pg("BAD1", [_run_name(0, 0), _run_name(0, 1)], "TMT131")]
    ds_dir = _write_dataset(tmp_path / "unknown-label", runs=_tmt_runs(1, 2), pg=_tmt_pg(1, 2, 3) + bad)
    with Dataset(ds_dir) as ds:
        pg = ds.validate(structures=["pg"], strict=strict)["pg"]

    ambiguous = _issues(pg, "ambiguous_grouped_run_sample")
    assert len(ambiguous) == 1
    assert "label 'TMT131'" in ambiguous[0].message
    assert "resolve to 0 samples ([])" in ambiguous[0].message
    assert not pg.is_valid


@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("check", ["_check_grouped_runs_referential", "_check_grouped_runs_sample_mapping"])
def test_pg_run_checks_run_in_bounded_memory(large_tmt_dir, check, strict):
    """The checks read pg once and resolve its 24 runs and 8 units, so 32 MB without spill is enough.

    Unnesting every pg row into its 12 runs (384,000 rows, each looked up in
    run) or joining them to their 4 channels (1.5 million rows) ran out of
    memory here under DuckDB 1.5 and 2.0; resolving the distinct runs and
    units needs less than half of the limit.
    """
    with _low_memory(large_tmt_dir) as ds:
        assert ds.pg.count() == 32_000
        result = ValidationResult(structure="pg")
        getattr(ds, check)(result, strict=strict)

    assert result.issues == []


@pytest.mark.parametrize(("strict", "severity"), [(False, "warning"), (True, "error")])
def test_sample_mapping_check_that_runs_out_of_memory_says_so(tmp_path, monkeypatch, strict, severity):
    """A resource failure is reported as such, not as schema drift."""
    ds_dir = _write_dataset(tmp_path / "oom", runs=_tmt_runs(1, 2), pg=_tmt_pg(1, 2, 3))
    with Dataset(ds_dir) as ds:
        execute = ds._engine.execute

        def out_of_memory(sql, params=None):
            if sql == _GROUPED_RUN_SAMPLE_MAPPING_SQL:
                raise duckdb.OutOfMemoryException("Out of Memory Error: failed to pin block")
            return execute(sql, params)

        monkeypatch.setattr(ds._engine, "execute", out_of_memory)
        pg = ds.validate(structures=["pg"], strict=strict)["pg"]

    skipped = _issues(pg, "sample_mapping_check_skipped")
    assert len(skipped) == 1
    assert skipped[0].severity == severity
    assert "ran out of memory or temp space" in skipped[0].message
    assert "schema drift" not in skipped[0].message
    assert pg.is_valid is not strict


@pytest.mark.parametrize(
    ("exc", "reason"),
    [
        (duckdb.OutOfMemoryException("x"), "ran out of memory or temp space"),
        (duckdb.IOException("x"), "I/O error"),
        (duckdb.BinderException("x"), "possible schema drift"),
        (duckdb.CatalogException("x"), "possible schema drift"),
        (duckdb.ConversionException("x"), "query error"),
    ],
)
def test_skipped_check_reason_names_the_failure(exc, reason):
    """Only binder and catalog errors point at the data's schema."""
    assert reason in _skipped_check_reason(exc)


# ---------------------------------------------------------------------------
# feature / psm / mz run_file_name -> run.run_file_name
# ---------------------------------------------------------------------------


def _mz_record(scan: int, run_file_name: str) -> dict:
    return {
        "id": f"spectrum_{scan}",
        "run_file_name": run_file_name,
        "scan": [scan],
        "ms_level": 2,
        "centroid": True,
        "scan_start_time": 1.0,
        "inverse_ion_mobility": None,
        "ion_injection_time": None,
        "total_ion_current": 10.0,
        "precursors": None,
        "mz": [100.0],
        "intensity": [5.0],
        "cv_params": None,
    }


def _write_run_keyed_dataset(ds_dir, *, feature_runs, psm_runs, mz_runs, with_run=True):
    """A one-plex TMT dataset whose feature, psm and mz rows name the given runs."""
    ds_dir.mkdir(parents=True, exist_ok=True)
    if with_run:
        with RunWriter(ds_dir / "exp.run.parquet") as writer:
            writer.write_batch(_tmt_runs(1, 2))
    features = []
    for i, run in enumerate(feature_runs):
        feature = make_feature_record(sequence=f"PEPTIDE{i}K", peptidoform=f"PEPTIDE{i}K", run_file_name=run)
        feature["intensities"] = [{"label": label, "intensity": 100.0} for label in TMT]
        features.append(feature)
    with FeatureWriter(ds_dir / "exp.feature.parquet") as writer:
        writer.write_batch(features)
    with PsmWriter(ds_dir / "exp.psm.parquet") as writer:
        writer.write_batch([make_psm_record(run_file_name=run, scan=[i]) for i, run in enumerate(psm_runs)])
    with MzWriter(ds_dir / "exp.mz.parquet") as writer:
        writer.write_batch([_mz_record(i, run) for i, run in enumerate(mz_runs)])
    with DatasetWriter(ds_dir / "exp.dataset.parquet") as writer:
        writer.write_batch([make_dataset_record()])
    return ds_dir


@pytest.mark.parametrize("strict", [False, True])
def test_feature_psm_and_mz_rows_naming_an_unknown_run_are_errors(tmp_path, strict):
    """A row whose run_file_name is not a run resolves to no run metadata or sample."""
    known = _run_name(0, 0)
    ds_dir = _write_run_keyed_dataset(
        tmp_path / "dangling",
        feature_runs=[known, "NOT_A_RUN", "NOT_A_RUN"],
        psm_runs=[known, "GHOST_RUN"],
        mz_runs=["GHOST_RUN", known, "OTHER_RUN"],
    )
    with Dataset(ds_dir) as ds:
        results = ds.validate(structures=["feature", "psm", "mz"], strict=strict)

    expected = {
        "feature": "2 feature row(s) name 1 run_file_name value(s) that are not in run.parquet ('NOT_A_RUN')",
        "psm": "1 psm row(s) name 1 run_file_name value(s) that are not in run.parquet ('GHOST_RUN')",
        "mz": "2 mz row(s) name 2 run_file_name value(s) that are not in run.parquet ('GHOST_RUN', 'OTHER_RUN')",
    }
    for view, message in expected.items():
        dangling = _issues(results[view], "dangling_run_file_name")
        assert len(dangling) == 1, view
        assert dangling[0].severity == "error"
        assert dangling[0].column == "run_file_name"
        assert dangling[0].message.startswith(message)
        assert not results[view].is_valid


@pytest.mark.parametrize("strict", [False, True])
def test_rows_naming_known_runs_pass_the_run_check(tmp_path, strict):
    runs = [_run_name(0, 0), _run_name(0, 1)]
    ds_dir = _write_run_keyed_dataset(tmp_path / "known", feature_runs=runs, psm_runs=runs, mz_runs=runs)
    with Dataset(ds_dir) as ds:
        results = ds.validate(structures=["feature", "psm", "mz"], strict=strict)

    for view in ("feature", "psm", "mz"):
        assert not _issues(results[view], "dangling_run_file_name")
        assert not [issue for issue in results[view].issues if issue.check.endswith("check_skipped")]


def test_dangling_run_message_counts_every_run_but_lists_five(tmp_path):
    known = _run_name(0, 0)
    ds_dir = _write_run_keyed_dataset(
        tmp_path / "many", feature_runs=[f"exp{i}" for i in range(7)], psm_runs=[known], mz_runs=[known]
    )
    with Dataset(ds_dir) as ds:
        feature = ds.validate(structures=["feature"])["feature"]

    (dangling,) = _issues(feature, "dangling_run_file_name")
    assert "7 feature row(s) name 7 run_file_name value(s)" in dangling.message
    assert "('exp0', 'exp1', 'exp2', 'exp3', 'exp4', ...)" in dangling.message


def test_mz_without_run_file_name_is_not_run_checked(tmp_path):
    """mz.run_file_name is optional: leaving it out is neither dangling nor a skipped check."""
    known = _run_name(0, 0)
    ds_dir = _write_run_keyed_dataset(tmp_path / "mz-no-run", feature_runs=[known], psm_runs=[known], mz_runs=["GHOST_RUN"])
    mz_path = ds_dir / "exp.mz.parquet"
    pq.write_table(pq.read_table(mz_path).drop_columns(["run_file_name"]), mz_path)
    with Dataset(ds_dir) as ds:
        mz = ds.validate(structures=["mz"], strict=True)["mz"]

    assert not _issues(mz, "dangling_run_file_name")
    assert not [issue for issue in mz.issues if issue.check.endswith("check_skipped")]


def test_run_check_needs_a_run_view(tmp_path):
    """Without run.parquet there is nothing to resolve run names against."""
    ds_dir = _write_run_keyed_dataset(
        tmp_path / "no-run", feature_runs=["NOT_A_RUN"], psm_runs=["NOT_A_RUN"], mz_runs=["NOT_A_RUN"], with_run=False
    )
    with Dataset(ds_dir) as ds:
        results = ds.validate(strict=True)

    assert not [issue for result in results.values() for issue in result.issues if issue.check == "dangling_run_file_name"]


@pytest.mark.parametrize(("strict", "severity"), [(False, "warning"), (True, "error")])
def test_run_check_that_cannot_run_is_reported(tmp_path, strict, severity):
    """A run view without run_file_name makes the check fail explicitly instead of passing."""
    known = _run_name(0, 0)
    ds_dir = _write_run_keyed_dataset(tmp_path / "drift", feature_runs=[known], psm_runs=[known], mz_runs=[known])
    with Dataset(ds_dir) as ds:
        ds._engine.execute("CREATE OR REPLACE VIEW run AS SELECT 1 AS wrong_column")
        psm = ds.validate(structures=["psm"], strict=strict)["psm"]

    skipped = _issues(psm, "referential_check_skipped")
    assert len(skipped) == 1
    assert skipped[0].severity == severity
    assert "psm run_file_name referential check could not run (possible schema drift)" in skipped[0].message


# ---------------------------------------------------------------------------
# feature (run_file_name, label) -> exactly one run.samples entry
# ---------------------------------------------------------------------------


def _write_feature_dataset(ds_dir, runs, features):
    """*features* are (run_file_name, [labels]) pairs; each becomes one feature row."""
    ds_dir.mkdir(parents=True, exist_ok=True)
    with RunWriter(ds_dir / "exp.run.parquet") as writer:
        writer.write_batch(runs)
    records = []
    for i, (run, labels) in enumerate(features):
        record = make_feature_record(sequence=f"PEPTIDE{i}K", peptidoform=f"PEPTIDE{i}K", run_file_name=run)
        record["intensities"] = [{"label": label, "intensity": 100.0} for label in labels]
        records.append(record)
    with FeatureWriter(ds_dir / "exp.feature.parquet") as writer:
        writer.write_batch(records)
    return ds_dir


def _lfq_runs(labels_per_run: dict[str, list[str]]) -> list[dict]:
    """Runs whose samples carry the given channel labels, one sample per label."""
    runs = []
    for run_name, labels in labels_per_run.items():
        run = make_run_record(run_accession=run_name, run_file_name=run_name)
        run["samples"] = [
            {"sample_accession": f"{run_name}_{label}", "label": label, "biological_replicate": 1, "technical_replicate": 1}
            for label in labels
        ]
        runs.append(run)
    return runs


@pytest.mark.parametrize("strict", [False, True])
def test_fractionated_tmt_features_resolve_every_label(tmp_path, strict):
    features = [(_run_name(plex, fraction), list(TMT)) for plex in range(2) for fraction in range(3) for _ in range(4)]
    ds_dir = _write_feature_dataset(tmp_path / "tmt", _tmt_runs(2, 3), features)
    with Dataset(ds_dir) as ds:
        feature = ds.validate(structures=["feature"], strict=strict)["feature"]

    assert not _issues(feature, "ambiguous_label_sample")
    assert not _issues(feature, "duplicate_intensity_label")
    assert not [issue for issue in feature.issues if issue.check.endswith("check_skipped")]


@pytest.mark.parametrize("strict", [False, True])
def test_feature_run_and_label_checks_run_in_bounded_memory(tmp_path, strict):
    """The checks group 24,000 fractionated TMT features to 24 runs and 96 (run, label) pairs.

    Joining every intensity to all channels of its run instead (384,000 rows)
    ran out of memory at this limit under DuckDB 2.0.
    """
    features = [(_run_name(plex, fraction), list(TMT)) for plex in range(2) for fraction in range(12) for _ in range(1000)]
    ds_dir = _write_feature_dataset(tmp_path / "tmt", _tmt_runs(2, 12), features)
    with _low_memory(ds_dir) as ds:
        result = ValidationResult(structure="feature")
        ds._check_run_file_names("feature", result, strict=strict)
        ds._check_feature_label_samples(result, strict=strict)

    assert result.issues == []


@pytest.mark.parametrize("strict", [False, True])
def test_feature_label_without_a_sample_is_an_error(tmp_path, strict):
    """A channel the run does not carry resolves to no sample: the quantity is orphaned."""
    run = _run_name(0, 0)
    ds_dir = _write_feature_dataset(
        tmp_path / "unknown-label", _tmt_runs(1, 1), [(run, ["TMT126", "TMT131"]), (run, ["TMT131"]), (run, ["TMT126"])]
    )
    with Dataset(ds_dir) as ds:
        feature = ds.validate(structures=["feature"], strict=strict)["feature"]

    (issue,) = _issues(feature, "ambiguous_label_sample")
    assert issue.severity == "error"
    assert issue.column == "intensities"
    assert issue.message.startswith(
        f"2 feature intensities in 1 (run_file_name, label) pair(s) match no sample_channel label of their run (('{run}', 'TMT131'))"
    )
    assert not feature.is_valid


@pytest.mark.parametrize("strict", [False, True])
def test_lfq_features_of_silac_runs_are_reported(tmp_path, strict):
    """Features labelled LFQ in runs whose channels are SILAC light and heavy name no sample."""
    runs = _lfq_runs({"A_F1": ["SILAC light", "SILAC heavy"], "A_F2": ["SILAC light", "SILAC heavy"]})
    ds_dir = _write_feature_dataset(tmp_path / "silac", runs, [("A_F1", ["LFQ"]), ("A_F2", ["LFQ"]), ("A_F2", ["LFQ"])])
    with Dataset(ds_dir) as ds:
        feature = ds.validate(strict=strict)["feature"]

    (issue,) = _issues(feature, "ambiguous_label_sample")
    assert "3 feature intensities in 2 (run_file_name, label) pair(s) match no sample_channel label" in issue.message
    assert not feature.is_valid


def test_feature_label_must_match_a_single_sample_run_exactly(tmp_path):
    """run.samples[].label is the canonical channel: a single-sample run does not excuse a mismatch.

    The feature->pg softlink joins on the exact label, so a feature labelled
    'raw' in an 'LFQ' run links to no pg row.
    """
    ds_dir = _write_feature_dataset(tmp_path / "lfq", _lfq_runs({"R1": ["LFQ"]}), [("R1", ["raw"]), ("R1", ["LFQ"])])
    with Dataset(ds_dir) as ds:
        feature = ds.validate(structures=["feature"])["feature"]

    (issue,) = _issues(feature, "ambiguous_label_sample")
    assert "1 feature intensities in 1 (run_file_name, label) pair(s) match no sample_channel label" in issue.message
    assert "('R1', 'raw')" in issue.message


@pytest.mark.parametrize("strict", [False, True])
def test_feature_label_matching_several_samples_is_an_error(tmp_path, strict):
    runs = _lfq_runs({"R1": ["TMT126", "TMT127N"]})
    runs[0]["samples"][1]["label"] = "TMT126"  # two samples on one channel
    ds_dir = _write_feature_dataset(tmp_path / "ambiguous", runs, [("R1", ["TMT126"])])
    with Dataset(ds_dir) as ds:
        feature = ds.validate(structures=["feature"], strict=strict)["feature"]

    (issue,) = _issues(feature, "ambiguous_label_sample")
    assert issue.severity == "error"
    assert "match the labels of several samples of their run (('R1', 'TMT126') -> ['R1_TMT126', 'R1_TMT127N'])" in issue.message
    assert not feature.is_valid


@pytest.mark.parametrize("strict", [False, True])
def test_feature_row_listing_a_label_twice_is_an_error(tmp_path, strict):
    """Each label is one sample's quantity; listing it three times triples that sample's value."""
    run = _run_name(0, 0)
    ds_dir = _write_feature_dataset(
        tmp_path / "labels-x3", _tmt_runs(1, 1), [(run, ["TMT126"] * 3), (run, ["TMT126", "TMT127N"]), (run, list(TMT) * 2)]
    )
    with Dataset(ds_dir) as ds:
        feature = ds.validate(structures=["feature"], strict=strict)["feature"]

    (issue,) = _issues(feature, "duplicate_intensity_label")
    assert issue.severity == "error"
    assert issue.message.startswith("2 feature row(s) list a label more than once in intensities")
    assert not _issues(feature, "ambiguous_label_sample")
    assert not feature.is_valid


def test_labels_of_an_unknown_run_are_left_to_the_run_check(tmp_path):
    """A row of an unknown run is reported once, as dangling, not also as an unresolved label."""
    ds_dir = _write_feature_dataset(tmp_path / "unknown-run", _tmt_runs(1, 1), [("NOT_A_RUN", list(TMT))])
    with Dataset(ds_dir) as ds:
        feature = ds.validate(structures=["feature"])["feature"]

    assert len(_issues(feature, "dangling_run_file_name")) == 1
    assert not _issues(feature, "ambiguous_label_sample")


@pytest.mark.parametrize(("strict", "severity"), [(False, "warning"), (True, "error")])
def test_feature_label_check_that_cannot_run_is_reported(tmp_path, strict, severity):
    """A run view without samples makes the label check fail explicitly instead of passing."""
    run = _run_name(0, 0)
    ds_dir = _write_feature_dataset(tmp_path / "drift", _tmt_runs(1, 1), [(run, list(TMT))])
    with Dataset(ds_dir) as ds:
        ds._engine.execute("CREATE OR REPLACE VIEW run AS SELECT 'plex0_F00' AS run_file_name")
        feature = ds.validate(structures=["feature"], strict=strict)["feature"]

    skipped = _issues(feature, "sample_mapping_check_skipped")
    assert len(skipped) == 1
    assert skipped[0].severity == severity
    assert "feature label sample-mapping check could not run (possible schema drift)" in skipped[0].message
