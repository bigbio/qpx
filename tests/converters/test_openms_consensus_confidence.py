"""Confidence and feature-metadata regressions for the consensusXML converter."""

from unittest.mock import Mock

import pyarrow.parquet as pq
import pytest

from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter
from qpx.converters.openms_consensus.feature_adapter import (
    feature_map_info,
    feature_records_for_cf,
    load_consensus_map,
    mass_error_ppm,
    pep_of,
    qvalue_of,
)
from qpx.converters.openms_consensus.psm_adapter import consensus_psms_to_records
from tests.converters.test_openms_consensus import (
    _TMT_CONSENSUSXML,
    _write_multirun_confidence_consensusxml,
)


def _scored_hit(score=None, meta=None):
    """Return a PeptideHit-like mock carrying a score and optional metadata."""
    values = meta or {}
    hit = Mock()
    hit.getScore.return_value = score
    hit.metaValueExists.side_effect = values.__contains__
    hit.getMetaValue.side_effect = values.__getitem__
    return hit


@pytest.mark.parametrize("streaming", [False, True])
def test_openms_consensus_confidence_is_attributed_per_run(tmp_path, streaming):
    """Each run receives confidence from its own PeptideIdentification."""
    cx = tmp_path / "multirun-confidence.consensusXML"
    _write_multirun_confidence_consensusxml(cx)
    written = OpenMSConsensusConverter().convert(
        str(cx),
        str(tmp_path / ("stream" if streaming else "mem")),
        output_prefix="t",
        structures=("feature",),
        streaming=streaming,
    )

    records = {record["run_file_name"]: record for record in pq.read_table(written["feature"]).to_pylist()}
    assert records["run_01"]["posterior_error_probability"] == pytest.approx(0.001)
    assert records["run_01"]["peptide_qvalue"] == pytest.approx(0.001)
    assert records["run_02"]["posterior_error_probability"] == pytest.approx(0.02)
    assert records["run_02"]["peptide_qvalue"] == pytest.approx(0.02)


def test_pep_of_reads_openms_meta_keys():
    """PEP is read from whichever meta key OpenMS used."""
    assert pep_of(_scored_hit(meta={"Posterior Error Probability_score": 0.004})) == 0.004
    assert pep_of(_scored_hit(meta={"PEP": 0.02})) == 0.02
    assert pep_of(_scored_hit(meta={"pep": 0.03})) == 0.03
    assert pep_of(_scored_hit()) is None


def test_qvalue_of_only_accepts_qvalue_score_types():
    """A raw search score must never be written to peptide_qvalue.

    Once it is in the column a consumer cannot tell an FDR from a Percolator SVM
    score, so anything that is not declared a q-value has to stay null
    (bigbio/qpx#284).
    """
    hit = _scored_hit(score=0.008)
    for score_type in ("q-value", "Q-Value", "qvalue", "FDR"):
        assert qvalue_of(hit, score_type) == 0.008, score_type
    for score_type in ("Posterior Error Probability", "expect", "svm_score", "", None):
        assert qvalue_of(hit, score_type) is None, score_type


def test_qvalue_of_returns_none_without_a_score():
    """A declared q-value score type remains null when the hit has no score."""
    assert qvalue_of(_scored_hit(score=None), "q-value") is None


def test_mass_error_ppm_is_computed_from_the_two_mz_values():
    """mass_error_ppm is derived when both m/z inputs are present."""
    assert mass_error_ppm(456.5589294433594, 456.5606994628906) == pytest.approx(3.876, abs=1e-2)
    assert mass_error_ppm(500.0, 499.995) == pytest.approx(-10.0, abs=1e-6)


def test_mass_error_ppm_is_none_only_when_not_measurable():
    """None means absent or unmeasurable, never a real zero."""
    assert mass_error_ppm(456.55, 456.55) == 0.0
    assert mass_error_ppm(456.55, None) is None
    assert mass_error_ppm(None, 456.55) is None
    assert mass_error_ppm(0.0, 456.55) is None


def test_mass_error_ppm_rejects_non_positive_mz():
    """A non-positive m/z is missing data, not a measurement."""
    assert mass_error_ppm(456.55, 0.0) is None
    assert mass_error_ppm(0.0, 456.55) is None
    assert mass_error_ppm(456.55, -1.0) is None
    assert mass_error_ppm(-456.55, 456.55) is None


def test_zero_charge_feature_has_no_mass_error(tmp_path):
    """A missing charge cannot produce a theoretical mass error."""
    path = tmp_path / "zero-charge.consensusXML"
    path.write_text(_TMT_CONSENSUSXML.replace('charge="2"', 'charge="0"'))
    consensus_map = load_consensus_map(str(path))
    map_info = feature_map_info(consensus_map)

    records = feature_records_for_cf(
        list(consensus_map)[0],
        map_info,
        group_map={"P12345": ["P12345"]},
    )

    assert records
    assert all(record["charge"] == 0 for record in records)
    assert all(record["mass_error_ppm"] is None for record in records)


def test_zero_charge_psm_has_no_mass_error(tmp_path):
    """A missing charge cannot produce a theoretical PSM mass error."""
    path = tmp_path / "zero-charge.consensusXML"
    path.write_text(_TMT_CONSENSUSXML.replace('charge="2"', 'charge="0"'))

    records = consensus_psms_to_records(str(path))

    assert records
    assert all(record["charge"] == 0 for record in records)
    assert all(record["calculated_mz"] is None for record in records)
    assert all(record["mass_error_ppm"] is None for record in records)


def test_unique_is_unknown_without_a_resolved_group(tmp_path):
    """No resolved group means unknown, not unique."""
    path = tmp_path / "in.consensusXML"
    path.write_text(_TMT_CONSENSUSXML)
    consensus_map = load_consensus_map(str(path))
    map_info = feature_map_info(consensus_map)
    consensus_feature = list(consensus_map)[0]

    unresolved = feature_records_for_cf(consensus_feature, map_info, group_map={})
    assert unresolved, "expected at least one feature record"
    for record in unresolved:
        assert record["anchor_protein"] == "P12345", "anchor still resolves from peptide evidence"
        assert record["unique"] is None, "unique must be unknown when no group resolved"

    resolved = feature_records_for_cf(consensus_feature, map_info, group_map={"P12345": ["P12345"]})
    assert all(record["unique"] is True for record in resolved), "a resolved group of one is unique"
