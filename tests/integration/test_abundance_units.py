"""Quantification units, labels and analyte keys of Dataset.intensity() and design_matrix()."""

import logging

import pyarrow.parquet as pq
import pytest

from qpx.dataset import Dataset
from qpx.writers import FeatureWriter, PgWriter, RunWriter
from tests.conftest import make_feature_record, make_pg_record, make_run_record


def _run(name, samples):
    """A run whose samples are (sample_accession, label, technical_replicate) triples."""
    run = make_run_record(run_accession=name, run_file_name=name)
    run["samples"] = [
        {"sample_accession": sample, "label": label, "biological_replicate": 1, "technical_replicate": replicate}
        for sample, label, replicate in samples
    ]
    return run


def _pg(accessions, grouped_runs, label, intensity):
    pg = make_pg_record(
        anchor_protein=accessions[0], pg_accessions=accessions, intensities=[{"label": label, "intensity": intensity}]
    )
    pg["grouped_runs"] = grouped_runs
    return pg


def _feature(**kwargs):
    return make_feature_record(run_file_name="r1", **kwargs)


def _write(path, runs, pgs=(), features=()):
    with RunWriter(path / "x.run.parquet") as w:
        w.write_batch(runs)
    if pgs:
        with PgWriter(path / "x.pg.parquet") as w:
            w.write_batch(list(pgs))
    if features:
        with FeatureWriter(path / "x.feature.parquet") as w:
            w.write_batch(list(features))
    return path


@pytest.fixture(name="replicates")
def _replicates(tmp_path):
    """S1 injected twice (technical replicates 1 and 2), S2 once."""
    runs = [_run("r1", [("S1", "LFQ", 1)]), _run("r2", [("S1", "LFQ", 2)]), _run("r3", [("S2", "LFQ", 1)])]
    pgs = [_pg(["P1"], ["r1"], "LFQ", 100.0), _pg(["P1"], ["r2"], "LFQ", 300.0), _pg(["P1"], ["r3"], "LFQ", 50.0)]
    return _write(tmp_path, runs, pgs)


def test_technical_replicates_are_combined_by_agg(replicates):
    """A sample's replicate units are averaged by default, not summed."""
    with Dataset(replicates) as ds:
        assert ds.design_matrix("protein").loc["S1", "P1"] == 200.0
        assert ds.design_matrix("protein", agg="sum").loc["S1", "P1"] == 400.0
        units = ds.design_matrix("protein", by="unit")
    assert units.loc[("S1", "LFQ", "r1"), "P1"] == 100.0
    assert units.loc[("S1", "LFQ", "r2"), "P1"] == 300.0


def test_output_path_applies_agg_and_rows(replicates, tmp_path):
    """The out-of-core pivot groups and aggregates like the in-memory one."""
    with Dataset(replicates) as ds:
        samples = pq.read_table(ds.design_matrix("protein", output_path=tmp_path / "s.parquet")).to_pylist()
        units = pq.read_table(ds.design_matrix("protein", by="unit", output_path=tmp_path / "u.parquet"))
    assert {row["sample_accession"]: row["P1"] for row in samples} == {"S1": 200.0, "S2": 50.0}
    assert units.num_rows == 3
    assert {"sample_accession", "label", "unit"} <= set(units.column_names)


def test_bridge_channel_is_one_unit_per_plex(tmp_path):
    """A pool in two plexes is two units; the fractions of one plex form one unit."""
    plex_a = [("POOL", "TMT126", 1), ("A1", "TMT127N", 1)]
    plex_b = [("POOL", "TMT126", 1), ("B1", "TMT127N", 1)]
    runs = [_run("a_f1", plex_a), _run("a_f2", plex_a), _run("b_f1", plex_b)]
    pgs = [_pg(["P1"], ["a_f1", "a_f2"], "TMT126", 5000.0), _pg(["P1"], ["b_f1"], "TMT126", 5000.0)]
    features = [
        make_feature_record(run_file_name=run, intensities=[{"label": "TMT126", "intensity": value}])
        for run, value in (("a_f1", 10.0), ("a_f2", 30.0), ("b_f1", 20.0))
    ]
    with Dataset(_write(tmp_path, runs, pgs, features)) as ds:
        proteins = ds.intensity("protein").to_df()
        assert sorted(proteins["unit"]) == ["a_f1", "b_f1"]
        assert ds.design_matrix("protein").loc["POOL", "P1"] == 5000.0
        peptides = ds.design_matrix("peptide", by="unit")
    # Fractions a_f1 and a_f2 add up within their unit.
    assert peptides.loc[("POOL", "TMT126", "a_f1"), "PEPTIDEK"] == 40.0
    assert peptides.loc[("POOL", "TMT126", "b_f1"), "PEPTIDEK"] == 20.0


def test_groups_sharing_a_leader_stay_apart(tmp_path):
    """Two groups led by the same protein are two analytes."""
    runs = [_run("r1", [("S1", "LFQ", 1)])]
    pgs = [_pg(["P1"], ["r1"], "LFQ", 1358.62), _pg(["P2", "P1"], ["r1"], "LFQ", 4187.25)]
    pgs[1]["anchor_protein"] = "P1"
    features = []
    for peptide, members in (("ONLYK", ["P1"]), ("SHAREDK", ["P2", "P1"])):
        feature = _feature(peptidoform=peptide, sequence=peptide, anchor_protein="P1")
        feature["pg_accessions"] = [{"accession": accession, "start": None, "end": None} for accession in members]
        features.append(feature)
    with Dataset(_write(tmp_path, runs, pgs, features)) as ds:
        matrix = ds.design_matrix("protein")
        view = ds.protein_view.intensity().to_df()
        assert ds.pg.by_protein("P1").count() == 2
        assert ds.pg.by_group(["P1", "P2"]).to_df()["intensity"].tolist() == [pytest.approx(4187.25)]
        assert ds.feature.by_protein("P1").count() == 2
        assert ds.feature.by_group(["P1"]).to_df()["peptidoform"].tolist() == ["ONLYK"]
    assert matrix.loc["S1", "P1"] == pytest.approx(1358.62)
    assert matrix.loc["S1", "P1;P2"] == pytest.approx(4187.25)
    assert sorted(view["protein_group"]) == ["P1", "P1;P2"]
    assert set(view["protein_accession"]) == {"P1"}


def test_peptide_key(tmp_path):
    """Peptides are keyed by peptidoform unless a sequence or precursor key is asked for."""
    runs = [_run("r1", [("S1", "LFQ", 1)])]
    features = [
        _feature(peptidoform="PEPTIDEK", intensities=[{"label": "LFQ", "intensity": 1.0}]),
        _feature(peptidoform="PEPTIDEK", charge=3, intensities=[{"label": "LFQ", "intensity": 2.0}]),
        _feature(peptidoform="PEPT[UNIMOD:21]IDEK", intensities=[{"label": "LFQ", "intensity": 4.0}]),
    ]
    with Dataset(_write(tmp_path, runs, features=features)) as ds:
        assert ds.design_matrix("peptide").loc["S1"].to_dict() == {"PEPT[UNIMOD:21]IDEK": 4.0, "PEPTIDEK": 3.0}
        assert ds.design_matrix("peptide", key="sequence").loc["S1", "PEPTIDEK"] == 7.0
        assert ds.design_matrix("peptide", key="precursor").loc["S1", "PEPTIDEK/3"] == 2.0
        with pytest.raises(ValueError, match="key"):
            ds.intensity("protein", key="sequence")


def test_two_labels_of_one_sample(tmp_path, caplog):
    """Pulsed SILAC: light and heavy of one sample are kept apart by unit rows and flagged by sample rows."""
    runs = [_run("r1", [("S1", "SILAC light", 1), ("S1", "SILAC heavy", 1)])]
    feature = _feature(intensities=[{"label": "SILAC light", "intensity": 1000.0}, {"label": "SILAC heavy", "intensity": 10.0}])
    with Dataset(_write(tmp_path, runs, features=[feature])) as ds:
        units = ds.design_matrix("peptide", by="unit")
        with caplog.at_level(logging.WARNING, logger="qpx.core.dataset_queries"):
            samples = ds.design_matrix("peptide")
    assert units.loc[("S1", "SILAC light", "r1"), "PEPTIDEK"] == 1000.0
    assert units.loc[("S1", "SILAC heavy", "r1"), "PEPTIDEK"] == 10.0
    assert samples.loc["S1", "PEPTIDEK"] == 505.0
    assert "several labels within one unit" in caplog.text


def test_unmatched_labels_of_a_single_sample_run_do_not_add_up(tmp_path):
    """Only a lone intensity falls back to the run's single sample; a light/heavy pair matching no label is not summed."""
    runs = [_run("r1", [("S1", "LFQ", 1)])]
    pair = _feature(intensities=[{"label": "L", "intensity": 1e6}, {"label": "H", "intensity": 4e6}])
    lone = _feature(peptidoform="ANOTHERK", sequence="ANOTHERK", intensities=[{"label": "raw", "intensity": 5.0}])
    with Dataset(_write(tmp_path, runs, features=[pair, lone])) as ds:
        rows = ds.intensity("peptide").to_df()
    assert rows[["feature_id", "intensity"]].values.tolist() == [["ANOTHERK", 5.0]]
