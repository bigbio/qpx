"""Decoys, q-value filters, labels and missing structures in the built-in summary views."""

import pytest

from qpx.dataset import Dataset
from tests.conftest import make_feature_record
from tests.integration.test_abundance_units import _pg, _run, _write


def _feature(peptide, members, intensities, is_decoy=False):
    feature = make_feature_record(
        run_file_name="plex", sequence=peptide, peptidoform=peptide, anchor_protein=members[0], intensities=intensities
    )
    feature["pg_accessions"] = [{"accession": accession, "start": None, "end": None} for accession in members]
    feature["is_decoy"] = is_decoy
    return feature


@pytest.fixture(name="plex")
def _plex(tmp_path):
    """One TMT run of samples A (TMT126) and B (TMT127N); P1 leads two groups; one decoy group and feature."""
    runs = [_run("plex", [("A", "TMT126", 1), ("B", "TMT127N", 1)])]
    pgs = [_pg(["P1"], ["plex"], "TMT126", 10.0), _pg(["P1", "P2"], ["plex"], "TMT126", 20.0)]
    pgs[1]["anchor_protein"] = "P1"
    pgs[1]["global_qvalue"] = 0.05
    decoy = _pg(["REV__P9"], ["plex"], "TMT127N", 5.0)
    decoy["is_decoy"] = True
    features = [
        _feature("ONLYK", ["P1"], [{"label": "TMT126", "intensity": 1.0}]),
        _feature("SHAREDK", ["P1", "P2"], [{"label": "TMT126", "intensity": 2.0}]),
        _feature("KEDITPEP", ["REV__P9"], [{"label": "TMT127N", "intensity": 3.0}], is_decoy=True),
    ]
    return _write(tmp_path, runs, [*pgs, decoy], features)


def test_summaries_count_target_groups(plex):
    """Run summaries leave decoys out and count the two groups led by P1 as two."""
    with Dataset(plex) as ds:
        for summary in (ds.identification_summary.summary(), ds.run_summary.summary()):
            row = summary.to_df().iloc[0]
            assert (row["n_proteins"], row["n_peptides"], row["n_features"]) == (2, 2, 2)
        metrics = ds.qc_view.metrics().to_df().iloc[0]
        assert (metrics["unique_proteins"], metrics["unique_peptides"]) == (2, 2)


def test_protein_and_peptide_views_leave_decoys_out(plex):
    """protein_view filters q-values only when asked; neither view returns decoys."""
    with Dataset(plex) as ds:
        proteins = ds.protein_view.intensity().to_df()
        confident = ds.protein_view.intensity(q_value_threshold=0.01).to_df()
        one_group = ds.protein_view.intensity().filter("protein_group = 'P1;P2'").to_df()
        peptides = ds.peptide_view.intensity().to_df()
    assert sorted(proteins["protein_group"]) == ["P1", "P1;P2"]
    assert confident["protein_group"].tolist() == ["P1"]
    assert one_group["intensity"].tolist() == [20.0]
    assert sorted(peptides["sequence"]) == ["ONLYK", "SHAREDK"]


def test_per_sample_counts_follow_the_label(plex):
    """A channel with no quantity gets no peptides or proteins of the other channel."""
    with Dataset(plex) as ds:
        peptides = ds.sample_summary.peptides_per_sample().to_df()
        proteins = ds.sample_summary.proteins_per_sample().to_df()
    assert peptides.values.tolist() == [["A", 2]]
    assert proteins.values.tolist() == [["A", 2]]


def test_missing_structures_are_named(tmp_path):
    """A view over an absent structure says which one is missing."""
    with Dataset(_write(tmp_path, [_run("r1", [("S1", "LFQ", 1)])])) as ds:
        with pytest.raises(ValueError, match="ProteinView requires the pg structure"):
            ds.protein_view.intensity()
        with pytest.raises(ValueError, match="requires the feature structure"):
            ds.run_summary.summary()
        with pytest.raises(ValueError, match="requires the psm structure"):
            ds.modification_view.frequency()
