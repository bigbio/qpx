"""Preserve protein properties and derive mass from the group's own anchor."""

from copy import deepcopy

import pyarrow.parquet as pq
import pytest
from defusedxml.ElementTree import fromstring, tostring

from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter
from tests.converters.test_openms_consensus import _TMT_CONSENSUSXML


def _set_properties(hit, coverage, probability):
    hit.set("coverage", str(coverage))
    score = hit.find("UserParam[@name='Posterior Probability_score']")
    if score is None:
        score = hit.makeelement("UserParam", {"type": "float", "name": "Posterior Probability_score"})
        hit.append(score)
    score.set("value", str(probability))


def _protein_root(grouped=False):
    root = fromstring(_TMT_CONSENSUSXML)
    protein_id = root.find(".//ProteinIdentification")
    protein_id.set("score_type", "q-value")
    protein_id.set("higher_score_better", "false")
    hit = protein_id.find("ProteinHit")
    hit.set("score", "0.01")
    _set_properties(hit, 80, 0.95)
    if grouped:
        member = fromstring('<ProteinHit id="PH_1" accession="P67890" score="0.02" sequence=""/>')
        _set_properties(member, 40, 0.5)
        protein_id.append(member)
        protein_id.append(fromstring('<UserParam type="string" name="indistinguishable_proteins_0" value="0,PH_0,PH_1"/>'))
    return root


def _pg_rows(root, tmp_path, streaming):
    source = tmp_path / "properties.consensusXML"
    source.write_bytes(tostring(root, encoding="utf-8", xml_declaration=True))
    written = OpenMSConsensusConverter().convert(str(source), str(tmp_path / "out"), structures=("pg",), streaming=streaming)
    return pq.read_table(written["pg"]).to_pylist()


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("grouped", [False, True])
@pytest.mark.parametrize(("coverage", "probability"), [(80, 0.95), (0, 0), (100, 1)])
def test_pg_preserves_anchor_properties(tmp_path, streaming, grouped, coverage, probability):
    """Recorded anchor values, including zero, survive both conversion paths."""
    root = _protein_root(grouped)
    _set_properties(root.find(".//ProteinHit"), coverage, probability)
    rows = _pg_rows(root, tmp_path, streaming)

    assert len(rows) == 2
    for row in rows:
        assert row["anchor_protein"] == "P12345"
        assert row["sequence_coverage"] == pytest.approx(coverage)
        assert row["additional_scores"] == [
            {"score_name": "posterior_probability", "score_value": pytest.approx(probability), "higher_better": True}
        ]
        assert row["global_qvalue"] == pytest.approx(0.01)
        assert row["molecular_weight"] is None


@pytest.mark.parametrize("streaming", [False, True])
def test_pg_does_not_borrow_missing_properties_from_group_member(tmp_path, streaming):
    """A well-annotated member cannot supply the anchor's unknown properties."""
    root = _protein_root(grouped=True)
    anchor = root.find(".//ProteinHit")
    del anchor.attrib["coverage"]
    anchor.remove(anchor.find("UserParam[@name='Posterior Probability_score']"))
    rows = _pg_rows(root, tmp_path, streaming)

    assert rows
    assert all(row["sequence_coverage"] is None and row["additional_scores"] is None for row in rows)


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("conflict", ["coverage", "probability", None])
def test_pg_repeated_anchor_records_must_agree(tmp_path, streaming, conflict):
    """Duplicate inference records retain only the properties that agree."""
    root = _protein_root()
    duplicate = deepcopy(root.find("IdentificationRun"))
    duplicate.set("id", "PI_1")
    hit = duplicate.find(".//ProteinHit")
    hit.set("id", "PH_1")
    _set_properties(hit, 40 if conflict == "coverage" else 80, 0.5 if conflict == "probability" else 0.95)
    root.insert(1, duplicate)
    rows = _pg_rows(root, tmp_path, streaming)

    assert rows
    for row in rows:
        expected_coverage = None if conflict == "coverage" else pytest.approx(80)
        assert row["sequence_coverage"] == expected_coverage
        if conflict == "probability":
            assert row["additional_scores"] is None
        else:
            assert row["additional_scores"][0]["score_value"] == pytest.approx(0.95)


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("grouped", [False, True])
@pytest.mark.parametrize(("sequence", "expected_kda"), [("ACD", 0.3073251), ("MUOJ", 0.6497068)])
def test_pg_molecular_weight_uses_anchor_sequence(tmp_path, streaming, grouped, sequence, expected_kda):
    """Use the anchor's average neutral mass in kDa, independent of its peptides."""
    root = _protein_root(grouped)
    root.find(".//ProteinHit").set("sequence", sequence)
    if grouped:
        root.find(".//ProteinHit[@id='PH_1']").set("sequence", "WWWW")
    rows = _pg_rows(root, tmp_path, streaming)

    assert len(rows) == 2
    for row in rows:
        assert row["anchor_protein"] == "P12345"
        assert row["molecular_weight"] == pytest.approx(expected_kda)


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("sequence", [None, "", "ACX", "ABZ", "C(Carbamidomethyl)PEPTIDEK"])
def test_pg_molecular_weight_unavailable_for_unknown_anchor_sequence(tmp_path, streaming, sequence):
    """Missing, ambiguous or modified sequences cannot borrow a member's mass."""
    root = _protein_root(grouped=True)
    anchor = root.find(".//ProteinHit")
    if sequence is None:
        del anchor.attrib["sequence"]
    else:
        anchor.set("sequence", sequence)
    root.find(".//ProteinHit[@id='PH_1']").set("sequence", "ACD")
    rows = _pg_rows(root, tmp_path, streaming)

    assert rows
    assert all(row["molecular_weight"] is None for row in rows)


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize(("duplicate_sequence", "expected_kda"), [("", 0.3073251), ("ACD", 0.3073251), ("ADC", None)])
def test_pg_molecular_weight_requires_consistent_anchor_sequence(tmp_path, streaming, duplicate_sequence, expected_kda):
    """Conflicting sequences remain unknown even when their masses are equal."""
    root = _protein_root()
    root.find(".//ProteinHit").set("sequence", "ACD")
    duplicate = deepcopy(root.find("IdentificationRun"))
    duplicate.set("id", "PI_1")
    hit = duplicate.find("ProteinIdentification/ProteinHit")
    hit.set("id", "PH_1")
    hit.set("sequence", duplicate_sequence)
    root.insert(1, duplicate)
    rows = _pg_rows(root, tmp_path, streaming)

    assert rows
    for row in rows:
        if expected_kda is None:
            assert row["molecular_weight"] is None
        else:
            assert row["molecular_weight"] == pytest.approx(expected_kda)
        assert row["sequence_coverage"] == pytest.approx(80)
