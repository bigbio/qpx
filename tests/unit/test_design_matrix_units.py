"""Dataset.intensity() and design_matrix() keep quantification units apart.

A quantity belongs to one quantification unit: a pg ``grouped_runs`` set and
label, or for peptides the runs that are fractions of one sample context. The
matrix used to sum every value of a sample, so the three technical replicates
of each PXD001819 sample came out three times too large and a TMT reference
channel present in every plex was summed across plexes.
"""

from __future__ import annotations

import math

import pandas as pd
import pyarrow.parquet as pq
import pytest

from qpx.dataset import Dataset
from qpx.writers import FeatureWriter, PgWriter, RunWriter
from tests.conftest import make_feature_record, make_pg_record, make_run_record


def _run(name, sample, label="LFQ", technical_replicate=1, fraction=None):
    record = make_run_record(run_accession=name, run_file_name=name)
    record["fraction"] = fraction
    record["samples"] = [
        {
            "sample_accession": sample,
            "label": label,
            "biological_replicate": 1,
            "technical_replicate": technical_replicate,
        }
    ]
    return record


def _write(path, runs, pgs, features):
    with RunWriter(path / "d.run.parquet") as writer:
        writer.write_batch(runs)
    with PgWriter(path / "d.pg.parquet") as writer:
        writer.write_batch(pgs)
    with FeatureWriter(path / "d.feature.parquet") as writer:
        writer.write_batch(features)


@pytest.fixture(name="replicates")
def _replicates(tmp_path):
    """S1 measured in three technical replicates (R1-R3), S2 once; S3 in two fractions."""
    runs = [_run(f"R{i}", "S1", technical_replicate=i) for i in (1, 2, 3)]
    runs += [_run("R4", "S2"), _run("F1", "S3", fraction="1"), _run("F2", "S3", fraction="2")]
    pgs, features = [], []
    for run, value in (("R1", 10.0), ("R2", 20.0), ("R3", 60.0), ("R4", 5.0)):
        pgs.append(make_pg_record(anchor_protein="P1", run_file_name=run, intensities=[{"label": "LFQ", "intensity": value}]))
        features.append(make_feature_record(run_file_name=run, intensities=[{"label": "LFQ", "intensity": value}]))
    pgs.append(make_pg_record(anchor_protein="P2", run_file_name="R4", intensities=[{"label": "LFQ", "intensity": 9.0}]))
    unit = make_pg_record(anchor_protein="P1", intensities=[{"label": "LFQ", "intensity": 7.0}])
    unit["grouped_runs"] = ["F2", "F1"]
    pgs.append(unit)
    for run, value in (("F1", 3.0), ("F2", 4.0)):
        features.append(make_feature_record(run_file_name=run, intensities=[{"label": "LFQ", "intensity": value}]))
    _write(tmp_path, runs, pgs, features)
    return tmp_path


@pytest.fixture(name="bridge")
def _bridge(tmp_path):
    """Two TMT plexes sharing a POOL reference channel (TMT126)."""
    runs, pgs, features = [], [], []
    for plex, pool in (("plex1", 1000.0), ("plex2", 1100.0)):
        run = _run(plex, "POOL", label="TMT126")
        run["samples"].append(
            {"sample_accession": f"{plex}_A", "label": "TMT127N", "biological_replicate": 1, "technical_replicate": 1}
        )
        runs.append(run)
        channels = [{"label": "TMT126", "intensity": pool}, {"label": "TMT127N", "intensity": 500.0}]
        pgs.append(make_pg_record(anchor_protein="P1", run_file_name=plex, intensities=channels))
        features.append(make_feature_record(run_file_name=plex, intensities=channels))
    _write(tmp_path, runs, pgs, features)
    return tmp_path


@pytest.mark.parametrize("level", ["protein", "peptide"])
def test_repeated_units_are_not_summed_by_default(replicates, level):
    with Dataset(replicates) as ds:
        with pytest.raises(ValueError, match=r"1 samples have more than one quantification unit \(e.g. 'S1'\)"):
            ds.design_matrix(level)


def test_reference_channel_in_every_plex_is_not_summed(bridge):
    with Dataset(bridge) as ds:
        with pytest.raises(ValueError, match="reference channel repeated across plexes"):
            ds.design_matrix()
        median = ds.design_matrix(agg="median")
        units = ds.design_matrix(by="unit")

    assert median.loc["POOL", "P1"] == 1050.0
    assert units.loc[("POOL", "plex1|TMT126"), "P1"] == 1000.0
    assert units.loc[("POOL", "plex2|TMT126"), "P1"] == 1100.0
    assert units.loc[("plex2_A", "plex2|TMT127N"), "P1"] == 500.0


def test_per_fraction_pg_rows_are_named_in_the_error(tmp_path):
    """CDAP and Spectronaut write a pg row per raw file, so each fraction is a protein unit of its own."""
    runs = [_run("F1", "S1", fraction="1"), _run("F2", "S1", fraction="2")]
    pgs = [
        make_pg_record(anchor_protein="P1", run_file_name=run, intensities=[{"label": "LFQ", "intensity": value}])
        for run, value in (("F1", 3.0), ("F2", 4.0))
    ]
    features = [make_feature_record(run_file_name=run, intensities=[{"label": "LFQ", "intensity": 1.0}]) for run in ("F1", "F2")]
    _write(tmp_path, runs, pgs, features)
    with Dataset(tmp_path) as ds:
        with pytest.raises(ValueError, match=r"fractions whose pg rows name only their own run .* agg='sum'"):
            ds.design_matrix()
        total = ds.design_matrix(agg="sum")
        peptides = ds.design_matrix("peptide")

    assert total.loc["S1", "P1"] == 7.0
    # The peptide units already join the two fractions.
    assert peptides.loc["S1", "PEPTIDEK"] == 2.0


@pytest.mark.parametrize(("level", "fraction_value"), [("protein", 7.0), ("peptide", 7.0)])
def test_aggregations_combine_a_samples_units(replicates, level, fraction_value):
    with Dataset(replicates) as ds:
        mean = ds.design_matrix(level, agg="mean")
        total = ds.design_matrix(level, agg="sum")
        largest = ds.design_matrix(level, agg="max")

    column = mean.columns[0]
    assert mean.loc["S1", column] == 30.0
    assert total.loc["S1", column] == 90.0
    assert largest.loc["S1", column] == 60.0
    # Fractions of one sample are one unit: the pg value once, the peptide summed over fractions.
    assert mean.loc["S3", column] == fraction_value


def test_by_unit_keeps_one_row_per_unit(replicates):
    with Dataset(replicates) as ds:
        matrix = ds.design_matrix("peptide", by="unit")

    assert list(matrix.index) == [
        ("S1", "R1|LFQ"),
        ("S1", "R2|LFQ"),
        ("S1", "R3|LFQ"),
        ("S2", "R4|LFQ"),
        ("S3", "F1;F2|LFQ"),
    ]
    assert matrix["PEPTIDEK"].tolist() == [10.0, 20.0, 60.0, 5.0, 7.0]


def test_intensity_names_the_quantification_unit(replicates):
    with Dataset(replicates) as ds:
        protein = ds.intensity("protein").to_df()
        peptide = ds.intensity("peptide").to_df()

    assert list(protein.columns) == ["sample_accession", "feature_id", "intensity", "quantification_unit"]
    assert set(protein["quantification_unit"]) == {"F1;F2|LFQ", "R1|LFQ", "R2|LFQ", "R3|LFQ", "R4|LFQ"}
    assert set(peptide["quantification_unit"]) == set(protein["quantification_unit"])


def test_missing_cells_stay_nan_by_default(replicates):
    """fillna defaults to None: a protein not measured in a sample is NaN, not 0."""
    with Dataset(replicates) as ds:
        matrix = ds.design_matrix(agg="mean")
        filled = ds.design_matrix(agg="mean", fillna=0.0)

    assert matrix.loc["S2", "P2"] == 9.0
    assert math.isnan(matrix.loc["S1", "P2"])
    assert filled.loc["S1", "P2"] == 0.0


def test_parquet_pivot_matches_the_in_memory_matrix(bridge, tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("matrices")
    out = out_dir / "m.parquet"
    with Dataset(bridge) as ds:
        expected = ds.design_matrix(agg="mean")
        assert ds.design_matrix(agg="mean", output_path=out) == out
        by_unit = ds.design_matrix(by="unit", output_path=out_dir / "u.parquet")

    on_disk = pd.read_parquet(out).set_index("sample_accession").sort_index()
    pd.testing.assert_frame_equal(on_disk, expected.sort_index(), check_names=False, check_dtype=False)
    assert pq.read_table(by_unit).num_rows == 4
    assert {"sample_accession", "quantification_unit", "P1"} <= set(pq.read_table(by_unit).column_names)


def test_parquet_pivot_refuses_repeated_units_without_agg(bridge, tmp_path_factory):
    out = tmp_path_factory.mktemp("matrices") / "m.parquet"
    with Dataset(bridge) as ds:
        with pytest.raises(ValueError, match="more than one quantification unit"):
            ds.design_matrix(output_path=out)
    assert not out.exists()


@pytest.mark.parametrize(("kwargs", "message"), [({"agg": "first"}, "agg must be"), ({"by": "run"}, "by must be")])
def test_invalid_options_are_rejected(bridge, kwargs, message):
    with Dataset(bridge) as ds:
        with pytest.raises(ValueError, match=message):
            ds.design_matrix(**kwargs)
