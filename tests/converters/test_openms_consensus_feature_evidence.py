"""Preserve peptide evidence and direct identification origins in Feature rows."""

from copy import deepcopy

import pyarrow.parquet as pq
import pytest
from defusedxml.ElementTree import fromstring, tostring

from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter
from qpx.converters.openms_consensus.feature_adapter import consensus_features_to_records, load_consensus_map
from qpx.converters.openms_consensus.pg_adapter import protein_group_maps
from qpx.converters.openms_consensus.streaming import StreamingConsensusMap
from tests.converters.test_openms_consensus import _TMT_CONSENSUSXML, _write_multirun_confidence_consensusxml


def _feature_rows(root, tmp_path, streaming, structures=("feature",)):
    source = tmp_path / "evidence.consensusXML"
    source.write_bytes(tostring(root, encoding="utf-8", xml_declaration=True))
    written = OpenMSConsensusConverter().convert(str(source), str(tmp_path / "out"), structures=structures, streaming=streaming)
    return pq.read_table(written["feature"]).to_pylist()


def _position_root():
    root = fromstring(_TMT_CONSENSUSXML)
    root.find(".//PeptideHit").attrib.update(start="1", end="8")
    return root


def _multirun_root(tmp_path):
    source = tmp_path / "multirun.consensusXML"
    _write_multirun_confidence_consensusxml(source)
    root = fromstring(source.read_text())
    hits = root.findall(".//PeptideHit")
    hits[0].attrib.update(start="1", end="8")
    hits[1].attrib.update(start="21", end="28")
    return root


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("structures", [("feature",), ("feature", "psm")])
def test_feature_preserves_positions_and_direct_run(tmp_path, streaming, structures):
    """Both readers retain evidence regardless of whether PSMs are exported."""
    rows = _feature_rows(_position_root(), tmp_path, streaming, structures)

    assert len(rows) == 1
    assert rows[0]["pg_positions"] == [{"protein_accession": "P12345", "start": 2, "end": 9}]
    assert rows[0]["id_run_file_name"] == "run_01"


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
def test_feature_keeps_multiple_positions_aligned_with_proteins(tmp_path, streaming, reverse):
    """Repeated occurrences survive without inventing positions for other members."""
    root = _position_root()
    protein_id = root.find(".//ProteinIdentification")
    for index, accession in [(1, "P67890"), (2, "P99999")]:
        protein_id.append(fromstring(f'<ProteinHit id="PH_{index}" accession="{accession}" score="0" sequence=""/>'))
    protein_id.append(fromstring('<UserParam type="string" name="indistinguishable_proteins_0" value="0,PH_0,PH_1,PH_2"/>'))
    mappings = [("PH_1", 10, 17), ("PH_0", 0, 7), ("PH_1", 30, 37), ("PH_1", 10, 17)]
    if reverse:
        mappings.reverse()
    hit = root.find(".//PeptideHit")
    hit.set("protein_refs", " ".join(ref for ref, _, _ in mappings))
    hit.set("start", " ".join(str(start) for _, start, _ in mappings))
    hit.set("end", " ".join(str(end) for _, _, end in mappings))
    rows = _feature_rows(root, tmp_path, streaming)

    assert rows[0]["pg_positions"] == [
        {"protein_accession": "P12345", "start": 1, "end": 8},
        {"protein_accession": "P67890", "start": 11, "end": 18},
        {"protein_accession": "P67890", "start": 31, "end": 38},
    ]
    assert {entry["accession"] for entry in rows[0]["pg_accessions"]} == {"P12345", "P67890", "P99999"}


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("coordinates", [{}, {"start": "-1", "end": "-1"}, {"start": "-1", "end": "8"}])
def test_feature_unknown_positions_remain_null(tmp_path, streaming, coordinates):
    """Absent and partially unknown coordinates are not converted into real positions."""
    root = fromstring(_TMT_CONSENSUSXML)
    root.find(".//PeptideHit").attrib.update(coordinates)
    rows = _feature_rows(root, tmp_path, streaming)

    assert rows[0]["pg_positions"] is None
    assert rows[0]["id_run_file_name"] == "run_01"


@pytest.mark.parametrize("streaming", [False, True])
def test_feature_positions_and_identification_are_attributed_per_run(tmp_path, streaming):
    """Each run receives its own evidence without clearing a consistent protein group."""
    rows = _feature_rows(_multirun_root(tmp_path), tmp_path, streaming)
    by_run = {row["run_file_name"]: row for row in rows}

    for run, start in [("run_01", 2), ("run_02", 22)]:
        assert by_run[run]["id_run_file_name"] == run
        assert by_run[run]["pg_positions"] == [{"protein_accession": "P12345", "start": start, "end": start + 7}]
        assert by_run[run]["anchor_protein"] == "P12345"


@pytest.mark.parametrize("streaming", [False, True])
def test_transferred_feature_has_positions_but_no_invented_identification_run(tmp_path, streaming):
    """Sequence positions can transfer while the actual origin of the best PSM is unknown."""
    root = _multirun_root(tmp_path)
    consensus = root.find(".//consensusElement")
    consensus.remove(consensus.findall("PeptideIdentification")[1])
    rows = _feature_rows(root, tmp_path, streaming)
    by_run = {row["run_file_name"]: row for row in rows}

    assert by_run["run_01"]["id_run_file_name"] == "run_01"
    assert by_run["run_02"]["id_run_file_name"] is None
    assert by_run["run_02"]["scan"] == []
    assert by_run["run_02"]["pg_positions"] == by_run["run_01"]["pg_positions"]


@pytest.mark.parametrize("streaming", [False, True])
def test_ambiguous_source_run_remains_null(tmp_path, streaming):
    """Multiple possible source runs cannot supply a direct identification origin."""
    root = _multirun_root(tmp_path)
    for pid in root.findall(".//PeptideIdentification"):
        pid.remove(pid.find("UserParam[@name='map_index']"))
    rows = _feature_rows(root, tmp_path, streaming)

    assert len(rows) == 2
    assert all(row["id_run_file_name"] is None for row in rows)


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("with_spectrum", [False, True])
def test_identification_origin_requires_a_spectrum_reference(tmp_path, streaming, with_spectrum):
    """Peptide annotations alone cannot identify a run containing a supporting PSM."""
    root = _position_root()
    consensus = root.find(".//consensusElement")
    pid = consensus.find("PeptideIdentification")
    del pid.attrib["spectrum_reference"]
    if with_spectrum:
        supporting = deepcopy(pid)
        supporting.set("spectrum_reference", "scan=43")
        consensus.append(supporting)
    rows = _feature_rows(root, tmp_path, streaming)

    assert rows[0]["id_run_file_name"] == ("run_01" if with_spectrum else None)
    assert rows[0]["pg_positions"] == [{"protein_accession": "P12345", "start": 2, "end": 9}]


@pytest.mark.parametrize("streaming", [False, True])
def test_conflicting_peptide_assignment_has_no_identification_origin(tmp_path, streaming):
    """An identification of another peptide cannot identify the exported feature."""
    root = _position_root()
    consensus = root.find(".//consensusElement")
    conflicting = deepcopy(consensus.find("PeptideIdentification"))
    conflicting.set("spectrum_reference", "scan=43")
    conflicting.find("PeptideHit").set("sequence", "ELVISLIVK")
    consensus.append(conflicting)
    rows = _feature_rows(root, tmp_path, streaming)

    assert rows[0]["id_run_file_name"] is None
    assert rows[0]["pg_positions"] is None


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("entry_point", ["converter", "adapter"])
@pytest.mark.parametrize("conflicting_peptide", [False, True])
def test_identification_origin_uses_the_recorded_merge_order(tmp_path, streaming, entry_point, conflicting_peptide):
    """Merge order resolves source runs while conflicting direct evidence stays null."""
    root = _multirun_root(tmp_path)
    root.find(".//ProteinIdentification").append(
        fromstring('<UserParam type="stringList" name="spectra_data" value="[run_02.mzML,run_01.mzML]"/>')
    )
    for index, pid in enumerate(root.findall(".//PeptideIdentification")):
        mapping = pid.find("UserParam[@name='map_index']")
        mapping.set("name", "id_merge_index")
        mapping.set("value", str(1 - index))
    if conflicting_peptide:
        consensus = root.find(".//consensusElement")
        conflicting = deepcopy(consensus.find("PeptideIdentification"))
        conflicting.set("spectrum_reference", "scan=44")
        conflicting.find("PeptideHit").set("sequence", "ELVISLIVK")
        consensus.append(conflicting)
    if entry_point == "converter":
        rows = _feature_rows(root, tmp_path, streaming, structures=("feature", "psm"))
    else:
        source = tmp_path / "adapter.consensusXML"
        source.write_bytes(tostring(root, encoding="utf-8", xml_declaration=True))
        consensus = StreamingConsensusMap(str(source)) if streaming else load_consensus_map(str(source))
        group_map, group_meta = protein_group_maps(consensus)
        rows = consensus_features_to_records(cm=consensus, group_map=group_map, group_meta=group_meta)

    by_run = {row["run_file_name"]: row for row in rows}
    assert by_run["run_01"]["id_run_file_name"] == (None if conflicting_peptide else "run_01")
    assert by_run["run_02"]["id_run_file_name"] == "run_02"
    assert by_run["run_01"]["pg_positions"] == (
        None if conflicting_peptide else [{"protein_accession": "P12345", "start": 2, "end": 9}]
    )
    assert by_run["run_02"]["pg_positions"][0]["start"] == 22
    assert by_run["run_01"]["scan"] == ([42, 44] if conflicting_peptide else [42])
    assert by_run["run_02"]["scan"] == [43]
    for run, confidence in [("run_01", 0.001), ("run_02", 0.02)]:
        assert by_run[run]["peptide_qvalue"] == pytest.approx(confidence)
        assert by_run[run]["posterior_error_probability"] == pytest.approx(confidence)
    if entry_point == "converter":
        psms = pq.read_table(tmp_path / "out" / "openms.psm.parquet").to_pylist()
        for psm in psms:
            assert psm["feature_id"] == by_run[psm["run_file_name"]]["feature_id"]


@pytest.mark.parametrize("streaming", [False, True])
def test_shared_peptide_across_groups_keeps_positions_on_each_protein(tmp_path, streaming):
    """A peptide whose evidence spans two inferred groups gets no group (that
    would be a guess), but its coordinates on each protein are recorded facts
    and each pg_positions entry names its own protein, so none are dropped."""
    root = _position_root()
    protein_id = root.find(".//ProteinIdentification")
    # A second protein that OpenMS left as its own group: P12345 and P67890 are
    # two singleton groups, and the peptide maps to both.
    protein_id.append(fromstring('<ProteinHit id="PH_1" accession="P67890" score="0" sequence=""/>'))
    hit = root.find(".//PeptideHit")
    hit.set("protein_refs", "PH_0 PH_1")
    hit.set("start", "0 10")
    hit.set("end", "7 17")

    rows = _feature_rows(root, tmp_path, streaming)

    assert rows[0]["pg_accessions"] is None
    assert rows[0]["anchor_protein"] is None
    assert rows[0]["pg_positions"] == [
        {"protein_accession": "P12345", "start": 1, "end": 8},
        {"protein_accession": "P67890", "start": 11, "end": 18},
    ]
