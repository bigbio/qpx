"""Peptide abundance is keyed by peptidoform and keeps channels apart.

The peptide long form and matrix were keyed by the unmodified sequence and
summed every label a single-sample run took: a heavy-labelled standard was added
to its light peptide (1e6 + 4e6 + 2e6 = 7e6), and a DIA-NN plexDIA feature whose
L and H channels both fell back to the run's one sample came out as L + H.
Whether labels are separate units is decided per quantification unit, over all
its runs, so a peptide seen in only one channel, or a fraction that saw only one
label, keeps that channel.
"""

from __future__ import annotations

import pytest

from qpx.dataset import Dataset
from qpx.writers import FeatureWriter, RunWriter
from tests.conftest import make_feature_record, make_run_record


def _run(run_file_name="run_01", sample="SAMPLE_01", sample_label="LFQ", fraction=None):
    run = make_run_record(run_accession=f"assay_{run_file_name}", run_file_name=run_file_name)
    run["samples"][0].update(sample_accession=sample, label=sample_label)
    run["fraction"] = fraction
    return run


def _write(path, features, sample_label="LFQ", runs=None):
    with RunWriter(path / "p.run.parquet") as writer:
        writer.write_batch(runs or [_run(sample_label=sample_label)])
    with FeatureWriter(path / "p.feature.parquet") as writer:
        writer.write_batch(features)


def test_heavy_standard_is_its_own_peptidoform(tmp_path):
    lfq = [{"label": "LFQ", "intensity": value} for value in (1e6, 4e6, 2e6)]
    _write(
        tmp_path,
        [
            make_feature_record(charge=2, scan=[1], intensities=[lfq[0]]),
            make_feature_record(charge=3, scan=[2], intensities=[lfq[1]]),
            make_feature_record(peptidoform="PEPTIDEK[UNIMOD:259]", scan=[3], intensities=[lfq[2]]),
        ],
    )
    with Dataset(tmp_path) as ds:
        matrix = ds.design_matrix("peptide")
        long = ds.intensity("peptide").to_df()

    # Charge states of one peptidoform are summed; the heavy form is not.
    assert matrix.loc["SAMPLE_01", "PEPTIDEK"] == 5e6
    assert matrix.loc["SAMPLE_01", "PEPTIDEK[UNIMOD:259]"] == 2e6
    assert set(long["quantification_unit"]) == {"run_01|LFQ"}


def test_labels_of_a_single_sample_run_are_not_summed(tmp_path):
    channels = [{"label": "L", "intensity": 2.46e6}, {"label": "H", "intensity": 5.49e6}]
    _write(tmp_path, [make_feature_record(intensities=channels)])
    with Dataset(tmp_path) as ds:
        long = ds.intensity("peptide").to_df().sort_values("quantification_unit")
        with pytest.raises(ValueError, match="several labels of one sample"):
            ds.design_matrix("peptide")
        units = ds.design_matrix("peptide", by="unit")

    assert long["quantification_unit"].tolist() == ["run_01|H", "run_01|L"]
    assert long["intensity"].tolist() == [5.49e6, 2.46e6]
    assert units.loc[("SAMPLE_01", "run_01|L"), "PEPTIDEK"] == 2.46e6


def test_single_label_features_still_fall_back_to_the_sole_sample(tmp_path):
    """DIA-NN's generic feature label still reaches the run's one sample, under its label."""
    _write(tmp_path, [make_feature_record(intensities=[{"label": "raw", "intensity": 10.0}])])
    with Dataset(tmp_path) as ds:
        long = ds.intensity("peptide").to_df()

    assert long[["sample_accession", "feature_id", "intensity", "quantification_unit"]].values.tolist() == [
        ["SAMPLE_01", "PEPTIDEK", 10.0, "run_01|LFQ"]
    ]


def test_single_channel_peptides_keep_their_channel_in_a_multi_label_run(tmp_path):
    """L/H is decided per unit: CCK (H only) and DDK (L only) are not filed under the sample label."""
    light, heavy = ({"label": label, "intensity": value} for label, value in (("L", 1.0), ("H", 10.0)))
    _write(
        tmp_path,
        [
            make_feature_record(sequence="AAK", peptidoform="AAK", scan=[1], intensities=[light, heavy]),
            make_feature_record(sequence="CCK", peptidoform="CCK", scan=[2], intensities=[{"label": "H", "intensity": 20.0}]),
            make_feature_record(sequence="DDK", peptidoform="DDK", scan=[3], intensities=[{"label": "L", "intensity": 3.0}]),
            # A label-free run in the same dataset keeps its sample label.
            make_feature_record(
                sequence="AAK", peptidoform="AAK", run_file_name="run_02", intensities=[{"label": "raw", "intensity": 5.0}]
            ),
        ],
        runs=[_run(), _run("run_02", "SAMPLE_02")],
    )
    with Dataset(tmp_path) as ds:
        long = ds.intensity("peptide").to_df().sort_values(["quantification_unit", "feature_id"])
        units = ds.design_matrix("peptide", by="unit")

    assert long[["feature_id", "intensity", "quantification_unit"]].values.tolist() == [
        ["AAK", 10.0, "run_01|H"],
        ["CCK", 20.0, "run_01|H"],
        ["AAK", 1.0, "run_01|L"],
        ["DDK", 3.0, "run_01|L"],
        ["AAK", 5.0, "run_02|LFQ"],
    ]
    assert units.index.tolist() == [("SAMPLE_01", "run_01|H"), ("SAMPLE_01", "run_01|L"), ("SAMPLE_02", "run_02|LFQ")]


def test_a_fraction_with_one_label_keeps_the_channels_of_its_unit(tmp_path):
    """L/H is decided over the unit's fractions: F2 saw only L, and its L joins F1's L."""
    light, heavy = ({"label": label, "intensity": value} for label, value in (("L", 1.0), ("H", 10.0)))
    _write(
        tmp_path,
        [
            make_feature_record(sequence="AAK", peptidoform="AAK", run_file_name="F1", scan=[1], intensities=[light, heavy]),
            make_feature_record(
                sequence="CCK", peptidoform="CCK", run_file_name="F2", scan=[2], intensities=[{"label": "L", "intensity": 3.0}]
            ),
            make_feature_record(
                sequence="AAK", peptidoform="AAK", run_file_name="F2", scan=[3], intensities=[{"label": "L", "intensity": 5.0}]
            ),
        ],
        runs=[_run("F1", "S", "label free sample", "1"), _run("F2", "S", "label free sample", "2")],
    )
    with Dataset(tmp_path) as ds:
        long = ds.intensity("peptide").to_df().sort_values(["quantification_unit", "feature_id"])

    assert long[["feature_id", "intensity", "quantification_unit"]].values.tolist() == [
        ["AAK", 10.0, "F1;F2|H"],
        ["AAK", 6.0, "F1;F2|L"],
        ["CCK", 3.0, "F1;F2|L"],
    ]
