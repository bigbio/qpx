"""MuData observations, decoys and duplicate cells.

A protein group is quantified once per quantification unit (a ``grouped_runs``
set of fractions) and label. The proteins modality used to explode
``grouped_runs`` into one observation per fraction, each holding a copy of the
unit's value, so a three-fraction TMT design summed to three times the pg
total. It also summed cells that received two intensities without saying so.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

pytest.importorskip("mudata")

from qpx.dataset import Dataset  # noqa: E402
from qpx.mudata import build_mudata  # noqa: E402
from qpx.writers import FeatureWriter, PgWriter, RunWriter  # noqa: E402
from tests.conftest import make_feature_record, make_pg_record, make_run_record  # noqa: E402

LABELS = ("TMT126", "TMT127N")


def _fractionated_tmt(path):
    """Two TMT plexes, each measured in two fractions; one pg row per plex unit and label."""
    runs, features, pgs = [], [], []
    for plex in ("plexA", "plexB"):
        for fraction in (1, 2):
            run = make_run_record(run_accession=f"{plex}_F{fraction}", run_file_name=f"{plex}_F{fraction}")
            run["fraction"] = str(fraction)
            run["samples"] = [
                {"sample_accession": f"{plex}_{label}", "label": label, "biological_replicate": 1, "technical_replicate": 1}
                for label in LABELS
            ]
            runs.append(run)
            features.append(
                make_feature_record(
                    run_file_name=run["run_file_name"],
                    intensities=[{"label": label, "intensity": 10.0} for label in LABELS],
                )
            )
        for anchor, decoy in (("P1", False), ("DECOY_P1", True)):
            record = make_pg_record(
                anchor_protein=anchor,
                is_decoy=decoy,
                intensities=[{"label": label, "intensity": 100.0 + i} for i, label in enumerate(LABELS)],
            )
            record["grouped_runs"] = [f"{plex}_F2", f"{plex}_F1"]
            pgs.append(record)
    with RunWriter(path / "t.run.parquet") as writer:
        writer.write_batch(runs)
    with FeatureWriter(path / "t.feature.parquet") as writer:
        writer.write_batch(features)
    with PgWriter(path / "t.pg.parquet") as writer:
        writer.write_batch(pgs)


def _build(path, **kwargs):
    with Dataset(path, file_prefix="t") as ds:
        return build_mudata(ds, **kwargs)


def test_protein_observations_are_quantification_units(tmp_path):
    _fractionated_tmt(tmp_path)
    proteins = _build(tmp_path, modalities=["proteins"]).mod["proteins"]

    assert list(proteins.obs_names) == [
        "plexA_F1;plexA_F2|TMT126",
        "plexA_F1;plexA_F2|TMT127N",
        "plexB_F1;plexB_F2|TMT126",
        "plexB_F1;plexB_F2|TMT127N",
    ]
    obs = proteins.obs
    assert list(obs["grouped_runs"]) == ["plexA_F1;plexA_F2"] * 2 + ["plexB_F1;plexB_F2"] * 2
    assert list(obs["sample_accession"]) == ["plexA_TMT126", "plexA_TMT127N", "plexB_TMT126", "plexB_TMT127N"]
    assert set(obs["fraction"]) == {"1;2"}
    assert "run_file_name" not in obs.columns
    # Each pg value appears once: 2 plexes x 2 groups x (100 + 101).
    np.testing.assert_allclose(proteins.X.toarray().sum(), 804.0)


def test_unit_observations_leave_disputed_replicate_numbers_empty(tmp_path):
    """A unit grouping two injections of a sample has no single technical replicate."""
    runs = []
    for injection in (1, 2):
        run = make_run_record(run_accession=f"inj{injection}", run_file_name=f"inj{injection}")
        run["samples"][0]["technical_replicate"] = injection
        runs.append(run)
    record = make_pg_record()
    record["grouped_runs"] = ["inj1", "inj2"]
    with RunWriter(tmp_path / "t.run.parquet") as writer:
        writer.write_batch(runs)
    with PgWriter(tmp_path / "t.pg.parquet") as writer:
        writer.write_batch([record])

    obs = _build(tmp_path, modalities=["proteins"], intensity_label="TMT126").mod["proteins"].obs

    assert list(obs.index) == ["inj1;inj2"]
    assert obs["sample_accession"].tolist() == ["SAMPLE_01"]
    assert obs["biological_replicate"].tolist() == [1]
    assert obs["technical_replicate"].isna().all()


def test_single_label_protein_observations_keep_run_names(tmp_path):
    """An unfractionated unit is named by its run, as before."""
    run = make_run_record()
    with RunWriter(tmp_path / "t.run.parquet") as writer:
        writer.write_batch([run])
    with PgWriter(tmp_path / "t.pg.parquet") as writer:
        writer.write_batch([make_pg_record()])

    proteins = _build(tmp_path, modalities=["proteins"], intensity_label="TMT126").mod["proteins"]

    assert list(proteins.obs_names) == ["run_01"]
    assert list(proteins.obs["sample_accession"]) == ["SAMPLE_01"]
    assert list(proteins.obs["run_accession"]) == ["assay_01"]


def _duplicate_precursor(path):
    """One precursor quantified by two features (e.g. MS/MS and MBR) in one run."""
    run = make_run_record()
    features = [
        make_feature_record(scan=[1], intensities=[{"label": "TMT126", "intensity": 100.0}]),
        make_feature_record(scan=[2], intensities=[{"label": "TMT126", "intensity": 300.0}]),
    ]
    with RunWriter(path / "t.run.parquet") as writer:
        writer.write_batch([run])
    with FeatureWriter(path / "t.feature.parquet") as writer:
        writer.write_batch(features)


def test_duplicate_cells_are_reported(tmp_path, caplog):
    _duplicate_precursor(tmp_path)
    with caplog.at_level(logging.WARNING, logger="qpx.mudata"):
        precursors = _build(tmp_path, modalities=["precursors"]).mod["precursors"]

    assert list(precursors.var_names) == ["PEPTIDEK|2"]
    assert precursors.X.toarray().tolist() == [[400.0]]
    assert precursors.uns["qpx_duplicates"] == {"aggregation": "sum", "cells": 1}
    assert "1 cells received more than one intensity" in caplog.text


def test_duplicate_cells_can_be_averaged(tmp_path):
    _duplicate_precursor(tmp_path)
    precursors = _build(tmp_path, modalities=["precursors"], duplicates="mean").mod["precursors"]

    np.testing.assert_allclose(precursors.X.toarray(), [[200.0]])
    assert precursors.uns["qpx_duplicates"] == {"aggregation": "mean", "cells": 1}


def test_unknown_duplicate_aggregation_is_rejected(tmp_path):
    _duplicate_precursor(tmp_path)
    with pytest.raises(ValueError, match="duplicates must be one of"):
        _build(tmp_path, modalities=["precursors"], duplicates="max")


def test_unit_observations_survive_the_h5mu_round_trip(tmp_path):
    import mudata as mu

    _fractionated_tmt(tmp_path)
    mdata = _build(tmp_path)
    mdata.write(str(tmp_path / "t.h5mu"))
    proteins = mu.read_h5mu(tmp_path / "t.h5mu").mod["proteins"]

    assert proteins.n_obs == 4
    assert proteins.uns["qpx_duplicates"]["cells"] == 0


def test_pg_rows_without_grouped_runs_are_not_an_observation(tmp_path):
    """A pg row with an empty grouped_runs belongs to no unit, so it adds no '' observation."""
    orphan = make_pg_record(anchor_protein="P2")
    orphan["grouped_runs"] = []
    with RunWriter(tmp_path / "t.run.parquet") as writer:
        writer.write_batch([make_run_record()])
    with PgWriter(tmp_path / "t.pg.parquet") as writer:
        writer.write_batch([make_pg_record(), orphan])

    single = _build(tmp_path, modalities=["proteins"], intensity_label="TMT126").mod["proteins"]
    expanded = _build(tmp_path, modalities=["proteins"], all_intensity_labels=True).mod["proteins"]

    assert list(single.obs_names) == ["run_01"]
    assert list(single.var_names) == ["P12345"]
    assert list(expanded.obs_names) == ["run_01|TMT126"]
    assert list(expanded.obs["grouped_runs"]) == ["run_01"]
