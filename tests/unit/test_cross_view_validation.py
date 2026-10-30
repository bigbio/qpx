"""Cross-view checks of Dataset.validate() that resolve rows through the run view.

The pg label -> sample check resolves each distinct (grouped runs, label) unit
through ``run.samples``, once per unit rather than once per pg row.
"""

import contextlib

import duckdb
import pytest

from qpx.core.data.schema import ValidationResult
from qpx.core.dataset_validation import _GROUPED_RUN_SAMPLE_MAPPING_SQL, _skipped_check_reason
from qpx.dataset import Dataset
from qpx.writers import DatasetWriter, PgWriter, RunWriter
from tests.conftest import make_dataset_record, make_pg_record, make_run_record

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
