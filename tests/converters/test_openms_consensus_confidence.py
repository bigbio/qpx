"""Confidence and feature-metadata regressions for the consensusXML converter."""

from copy import deepcopy
from unittest.mock import Mock
from xml.etree import ElementTree as ET

import pyarrow.parquet as pq
import pytest
from defusedxml.ElementTree import fromstring

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


@pytest.mark.parametrize(("score", "expected"), [(None, None), ("", None), ("0", 0.0), ("0.125", 0.125)])
def test_streaming_missing_scores_are_not_zero(tmp_path, score, expected):
    """A native hit's default zero must not turn missing confidence into zero FDR."""
    root = fromstring(_TMT_CONSENSUSXML)
    for identification in root.findall(".//ProteinIdentification") + root.findall(".//PeptideIdentification"):
        identification.set("score_type", "q-value")
    for hit in root.findall(".//ProteinHit") + root.findall(".//PeptideHit"):
        if score is None:
            del hit.attrib["score"]
        else:
            hit.set("score", score)
    source = tmp_path / "scores.consensusXML"
    ET.ElementTree(root).write(source, encoding="utf-8", xml_declaration=True)
    written = OpenMSConsensusConverter().convert(
        str(source), str(tmp_path / "out"), structures=("feature", "psm", "pg"), streaming=True
    )

    for view, field in (("feature", "peptide_qvalue"), ("pg", "global_qvalue")):
        rows = pq.read_table(written[view]).to_pylist()
        assert rows
        assert all(row[field] == expected for row in rows)
    psms = pq.read_table(written["psm"]).to_pylist()
    assert psms
    scores = [
        score["score_value"] for row in psms for score in row["additional_scores"] or [] if score["score_name"] == "q-value"
    ]
    assert scores == ([expected] * len(psms) if expected is not None else [])


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


# Shared-leader groups [A,B] and [A,C] with protein q-values and GN= genes.
# Distinct q-values per member make "best member" observable, and the shared
# leader means keying the annotation on anchor_protein would give both groups
# the same values — the bug class of bigbio/qpx#266.
def _annotated_shared_leader_consensusxml():
    from tests.converters.test_openms_consensus import _SHARED_LEADER_CONSENSUSXML

    xml = _SHARED_LEADER_CONSENSUSXML.replace(
        '<ProteinIdentification score_type="" higher_score_better="true"',
        '<ProteinIdentification score_type="q-value" higher_score_better="false"',
    )
    # Real quantms consensusXML carries the FASTA header as a "Description"
    # UserParam, not a description attribute — pyopenms ignores the attribute,
    # so a fixture using it would test a file OpenMS never writes.
    for index, (acc, score, gene) in enumerate((("A", "0.004", "GENEA"), ("B", "0.001", "GENEB"), ("C", "0.009", "GENEC"))):
        xml = xml.replace(
            f'<ProteinHit id="PH_{index}" accession="{acc}" score="0" sequence=""></ProteinHit>',
            f'<ProteinHit id="PH_{index}" accession="{acc}" score="{score}" sequence="">'
            f'<UserParam type="string" name="Description" value="Protein {acc} OS=Homo sapiens GN={gene} PE=1"/>'
            "</ProteinHit>",
        )
    assert xml.count('name="Description"') == 3
    return xml


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("reverse_evidence", [False, True])
@pytest.mark.parametrize("reverse_groups", [False, True])
def test_feature_carries_its_protein_groups_qvalue_and_genes(tmp_path, streaming, reverse_evidence, reverse_groups):
    """feature.pg_global_qvalue / gg_names were null on every OpenMS dataset.

    The converter computed both for the pg view and discarded them, so on
    PXD000612 pg carried them for all 10,105 groups while the feature view had
    0% — although every feature's group exists in pg. Each feature must now carry
    exactly its own group's values, including when groups share a leader.
    """
    import duckdb

    cx = tmp_path / "annotated_shared_leader.consensusXML"
    root = fromstring(_annotated_shared_leader_consensusxml())
    if reverse_evidence:
        for hit in root.findall(".//PeptideHit"):
            hit.set("protein_refs", " ".join(reversed(hit.get("protein_refs").split())))
    if reverse_groups:
        groups = root.findall(".//ProteinIdentification/UserParam")
        first, second = groups
        first_value, second_value = first.get("value"), second.get("value")
        first.set("value", second_value)
        second.set("value", first_value)
    cx.write_text(ET.tostring(root, encoding="unicode"))
    written = OpenMSConsensusConverter().convert(
        str(cx),
        str(tmp_path / ("stream" if streaming else "mem")),
        output_prefix="t",
        structures=("feature", "pg"),
        streaming=streaming,
    )
    con = duckdb.connect()
    features = {
        seq: (qv, genes)
        for seq, qv, genes in con.execute(
            "SELECT sequence, pg_global_qvalue, gg_names FROM read_parquet($1)", [str(written["feature"])]
        ).fetchall()
    }
    assert features["PEPTIDEK"][0] == pytest.approx(0.001)
    assert sorted(features["PEPTIDEK"][1]) == ["GENEA", "GENEB"]
    assert features["ELVISLIVK"][0] == pytest.approx(0.004)
    assert sorted(features["ELVISLIVK"][1]) == ["GENEA", "GENEC"]

    # And the two views agree group-for-group.
    disagreements = con.execute(
        """
        SELECT count(*) FROM read_parquet($1) f
        JOIN read_parquet($2) p
          ON list_sort(list_transform(f.pg_accessions, x -> x.accession)) = list_sort(p.pg_accessions)
        WHERE f.pg_global_qvalue IS DISTINCT FROM p.global_qvalue
           OR list_sort(f.gg_names) IS DISTINCT FROM list_sort(p.gg_names)
        """,
        [str(written["feature"]), str(written["pg"])],
    ).fetchone()[0]
    assert disagreements == 0
    con.close()


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("references", ["PH_0", "PH_1 PH_2", "PH_2 PH_1", ""])
def test_ambiguous_feature_group_keeps_quantification(tmp_path, streaming, references):
    """Ambiguous protein evidence leaves annotations null without losing intensities."""
    root = fromstring(_annotated_shared_leader_consensusxml())
    for hit in root.findall(".//PeptideHit"):
        hit.set("protein_refs", references)
    cx = tmp_path / "ambiguous.consensusXML"
    cx.write_text(ET.tostring(root, encoding="unicode"))
    written = OpenMSConsensusConverter().convert(
        str(cx), str(tmp_path / "out"), output_prefix="t", structures=("feature", "pg"), streaming=streaming
    )
    features = pq.read_table(written["feature"]).to_pylist()
    assert len(features) == 2
    assert sorted(row["intensities"][0]["intensity"] for row in features) == [1000.0, 3000.0]
    for row in features:
        for field in ("pg_accessions", "pg_global_qvalue", "gg_accessions", "gg_names", "unique"):
            assert row[field] is None, field
        assert row["anchor_protein"] == ("A" if references == "PH_0" else None)


def _separate_identification_groups_xml(reverse_identifications=False):
    """One consensus feature, two runs, shared A assigned to different source groups."""
    root = fromstring(_annotated_shared_leader_consensusxml())
    first = root.find("IdentificationRun")
    second = deepcopy(first)
    first.set("date", "2026-09-13T00:00:00")
    second.set("date", "2026-09-14T00:00:00")
    second.set("id", "PI_1")
    root.insert(1, second)
    for run, excluded, group_number in [(first, "C", "1"), (second, "B", "0")]:
        proteins = run.find("ProteinIdentification")
        proteins.remove(proteins.find(f"ProteinHit[@accession='{excluded}']"))
        proteins.remove(proteins.find(f"UserParam[@name='indistinguishable_proteins_{group_number}']"))
    for hit in second.findall(".//ProteinHit"):
        hit.set("id", hit.get("id").replace("PH_", "SECOND_"))
    group = second.find(".//ProteinIdentification/UserParam")
    group.set("name", "indistinguishable_proteins_0")
    group.set("value", group.get("value").replace("PH_", "SECOND_"))
    maps = root.find("mapList")
    maps.set("count", "2")
    ET.SubElement(maps, "map", id="1", name="run_02.mzML", unique_id="2", label="label-free", size="1")
    elements = root.find("consensusElementList")
    cf, other = list(elements)
    for number, element in enumerate((cf, other)):
        pid = element.find("PeptideIdentification")
        pid.set("identification_run_ref", f"PI_{number}")
        pid.find("PeptideHit").set("protein_refs", "PH_0" if number == 0 else "SECOND_0")
        pid.find("PeptideHit").set("sequence", "PEPTIDEK")
        ET.SubElement(pid, "UserParam", type="int", name="map_index", value=str(number))
    sub = other.find("groupedElementList/element")
    sub.set("map", "1")
    cf.find("groupedElementList").append(sub)
    cf.append(other.find("PeptideIdentification"))
    elements.remove(other)
    if reverse_identifications:
        root.remove(second)
        root.insert(0, second)
    return ET.tostring(root, encoding="unicode")


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("reverse_identifications", [False, True])
@pytest.mark.parametrize("merge_index", [False, True])
def test_feature_group_uses_its_runs_identification(tmp_path, streaming, reverse_identifications, merge_index):
    """Each run resolves its group from its own source identification."""
    root = fromstring(_separate_identification_groups_xml(reverse_identifications))
    if merge_index:
        indices = {}
        for index, identification in enumerate(root.findall("IdentificationRun")):
            identifier = identification.get("id")
            indices[identifier] = str(index)
            run = "run_01" if identifier == "PI_0" else "run_02"
            ET.SubElement(
                identification.find("ProteinIdentification"),
                "UserParam",
                type="stringList",
                name="spectra_data",
                value=f"[{run}.mzML]",
            )
        for pid in root.findall(".//PeptideIdentification"):
            mapping = pid.find("UserParam[@name='map_index']")
            mapping.set("name", "id_merge_index")
            mapping.set("value", indices[pid.get("identification_run_ref")])
    cx = tmp_path / "separate_identifications.consensusXML"
    cx.write_text(ET.tostring(root, encoding="unicode"))
    written = OpenMSConsensusConverter().convert(
        str(cx), str(tmp_path / "out"), output_prefix="t", structures=("feature", "pg"), streaming=streaming
    )
    features = {row["run_file_name"]: row for row in pq.read_table(written["feature"]).to_pylist()}
    assert set(features) == {"run_01", "run_02"}
    for run, member, qvalue, intensity in [("run_01", "B", 0.001, 1000.0), ("run_02", "C", 0.004, 3000.0)]:
        row = features[run]
        assert [acc["accession"] for acc in row["pg_accessions"]] == ["A", member]
        assert row["pg_global_qvalue"] == pytest.approx(qvalue)
        assert row["gg_names"] == row["gg_accessions"] == ["GENEA", f"GENE{member}"]
        assert row["intensities"][0]["intensity"] == intensity


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("add_matching_pid", [False, True])
def test_sequence_conflict_clears_run_protein_fields(tmp_path, streaming, add_matching_pid):
    """A conflicting identification must not inherit another run's protein group."""
    root = fromstring(_separate_identification_groups_xml())
    cf = root.find("consensusElementList/consensusElement")
    pid = cf.findall("PeptideIdentification")[1]
    if add_matching_pid:
        cf.append(deepcopy(pid))
    pid.find("PeptideHit").set("sequence", "ANOTHERK")
    cx = tmp_path / "sequence_conflict.consensusXML"
    cx.write_text(ET.tostring(root, encoding="unicode"))

    written = OpenMSConsensusConverter().convert(
        str(cx), str(tmp_path / "out"), output_prefix="t", structures=("feature", "pg"), streaming=streaming
    )
    features = {row["run_file_name"]: row for row in pq.read_table(written["feature"]).to_pylist()}
    assert set(features) == {"run_01", "run_02"}
    matched = features["run_01"]
    assert [entry["accession"] for entry in matched["pg_accessions"]] == ["A", "B"]
    assert matched["pg_global_qvalue"] == pytest.approx(0.001)
    assert matched["gg_names"] == matched["gg_accessions"] == ["GENEA", "GENEB"]
    for field in ("anchor_protein", "pg_accessions", "unique", "pg_global_qvalue", "gg_accessions", "gg_names"):
        assert features["run_02"][field] is None, field
    for run, intensity in (("run_01", 1000.0), ("run_02", 3000.0)):
        assert features[run]["intensities"][0]["intensity"] == intensity


@pytest.mark.parametrize("streaming", [False, True])
def test_pg_names_come_from_uniprot_entry_names(tmp_path, streaming):
    """pg.pg_names was null on every OpenMS dataset although the accessions carry
    the entry name (``sp|ACC|NAME``); DIA-NN fills it (bigbio/qpx#300)."""
    import duckdb

    from tests.converters.test_openms_consensus import _SHARED_LEADER_CONSENSUSXML

    xml = _SHARED_LEADER_CONSENSUSXML
    for acc, full in (("A", "sp|P11111|AAA_HUMAN"), ("B", "sp|P22222|BBB_HUMAN"), ("C", "sp|P33333|CCC_HUMAN")):
        xml = xml.replace(f'accession="{acc}"', f'accession="{full}"')
    cx = tmp_path / "named.consensusXML"
    cx.write_text(xml)
    written = OpenMSConsensusConverter().convert(
        str(cx), str(tmp_path / ("stream" if streaming else "mem")), output_prefix="t", structures=("pg",), streaming=streaming
    )

    names = {
        tuple(sorted(accs)): name_list
        for accs, name_list in duckdb.connect()
        .execute("SELECT DISTINCT pg_accessions, pg_names FROM read_parquet($1)", [str(written["pg"])])
        .fetchall()
    }
    assert names[("sp|P11111|AAA_HUMAN", "sp|P22222|BBB_HUMAN")] in (["AAA_HUMAN", "BBB_HUMAN"], ["BBB_HUMAN", "AAA_HUMAN"])
    assert names[("sp|P11111|AAA_HUMAN", "sp|P33333|CCC_HUMAN")] in (["AAA_HUMAN", "CCC_HUMAN"], ["CCC_HUMAN", "AAA_HUMAN"])


@pytest.mark.parametrize("streaming", [False, True])
def test_pg_names_stay_null_for_bare_accessions(tmp_path, streaming):
    """A bare accession has no entry name; pg_names must not echo pg_accessions."""
    import duckdb

    from tests.converters.test_openms_consensus import _SHARED_LEADER_CONSENSUSXML

    cx = tmp_path / "bare.consensusXML"
    cx.write_text(_SHARED_LEADER_CONSENSUSXML)
    written = OpenMSConsensusConverter().convert(
        str(cx), str(tmp_path / ("stream" if streaming else "mem")), output_prefix="t", structures=("pg",), streaming=streaming
    )

    assert duckdb.connect().execute("SELECT count(pg_names) FROM read_parquet($1)", [str(written["pg"])]).fetchone()[0] == 0
