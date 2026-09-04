"""consensusXML -> QPX converter (interim path).

Reads a tiny literal consensusXML fixture (no pyopenms construction) so the
extraction runs in CI without a large fixture file and without depending on the
pyopenms setter APIs, which differ across versions.
"""

import duckdb
import pytest

pytest.importorskip("pyopenms")
import pyopenms as oms  # noqa: E402

import qpx  # noqa: E402
from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter  # noqa: E402
from qpx.converters.openms_consensus.feature_adapter import to_proforma  # noqa: E402


@pytest.mark.parametrize(
    ("openms_seq", "expected"),
    [
        ("PEPTIDEK", "PEPTIDEK"),
        (".(TMT6plex)THSQEEM(Oxidation)QHMQR", "[UNIMOD:737]-THSQEEM[UNIMOD:35]QHMQR"),
        ("C(Carbamidomethyl)PEPTIDEK", "C[UNIMOD:4]PEPTIDEK"),
        ("PEPTIDER.(Amidated)", "PEPTIDER-[UNIMOD:2]"),
    ],
)
def test_to_proforma(openms_seq, expected):
    assert to_proforma(oms.AASequence.fromString(openms_seq)) == expected


# A 2-channel isobaric consensusXML written as a literal fixture: one peptide,
# 2 TMT channels (126/127), 1 protein. Kept as text — not built through pyopenms
# constructors — so the test exercises only the *read* path and is immune to
# pyopenms-version drift in the setter APIs (e.g. list vs PeptideIdentificationList).
# NOTE: experiment_type is "label-free" ON PURPOSE — real quantms IsobaricWorkflow
# output stamps TMT/iTRAQ runs "label-free" while the maps carry tmt6plex_* labels,
# so the channels must be detected from the map labels, not experiment_type.
_TMT_CONSENSUSXML = """<?xml version="1.0" encoding="ISO-8859-1"?>
<consensusXML version="1.7" experiment_type="label-free"
  xsi:noNamespaceSchemaLocation="https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/ConsensusXML_1_7.xsd"
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <IdentificationRun id="PI_0" date="0000-00-00T00:00:00" search_engine="" search_engine_version="">
    <SearchParameters db="" db_version="" taxonomy="" mass_type="monoisotopic" charges=""
      enzyme="unknown_enzyme" missed_cleavages="0" precursor_peak_tolerance="0"
      precursor_peak_tolerance_ppm="false" peak_mass_tolerance="0"
      peak_mass_tolerance_ppm="false">
    </SearchParameters>
    <ProteinIdentification score_type="" higher_score_better="true" significance_threshold="0">
      <ProteinHit id="PH_0" accession="P12345" score="0" sequence="">
      </ProteinHit>
    </ProteinIdentification>
  </IdentificationRun>
  <mapList count="2">
    <map id="0" name="run_01.mzML" unique_id="1" label="tmt6plex_126" size="1">
    </map>
    <map id="1" name="run_01.mzML" unique_id="2" label="tmt6plex_127" size="1">
    </map>
  </mapList>
  <consensusElementList>
    <consensusElement id="e_0" quality="0.0" charge="2">
      <centroid rt="100.123456" mz="450.251234" it="0.0"/>
      <groupedElementList>
        <element map="0" id="0" rt="100.123456" mz="450.251234" it="1000.0"/>
        <element map="1" id="1" rt="100.123456" mz="450.251234" it="2000.0"/>
      </groupedElementList>
      <PeptideIdentification identification_run_ref="PI_0" score_type=""
        higher_score_better="true" significance_threshold="0" MZ="450.26" RT="100"
        spectrum_reference="controllerType=0 controllerNumber=1 scan=42">
        <PeptideHit score="0" sequence="PEPTIDEK" charge="2" protein_refs="PH_0">
          <UserParam type="string" name="target_decoy" value="target"/>
          <UserParam type="float" name="Posterior Error Probability_score" value="1.0e-03"/>
        </PeptideHit>
      </PeptideIdentification>
    </consensusElement>
  </consensusElementList>
</consensusXML>
"""


def _write_tmt_consensusxml(path):
    """Write the literal 2-channel TMT consensusXML fixture to ``path``."""
    path.write_text(_TMT_CONSENSUSXML)


def test_converter_rejects_invalid_structure_before_creating_output(tmp_path):
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="Unknown.*bogus"):
        OpenMSConsensusConverter().convert(
            str(tmp_path / "missing.consensusXML"),
            str(out),
            structures=("bogus",),
        )
    assert not out.exists()


def test_converter_requires_sdrf_for_requested_metadata(tmp_path):
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="SDRF.*run"):
        OpenMSConsensusConverter().convert(
            str(tmp_path / "missing.consensusXML"),
            str(out),
            structures=("run",),
        )
    assert not out.exists()


def test_converter_writes_only_requested_sdrf_structure(tmp_path):
    sdrf = tmp_path / "test.sdrf.tsv"
    sdrf.write_text(
        "source name\tcharacteristics[organism]\tcharacteristics[organism part]\t"
        "comment[data file]\tcomment[label]\n"
        "S1\tHomo sapiens\tliver\trun_01.raw\tlabel free sample\n"
    )
    out = tmp_path / "out"
    written = OpenMSConsensusConverter().convert(
        str(tmp_path / "unused.consensusXML"),
        str(out),
        sdrf_path=str(sdrf),
        structures=("run",),
    )

    assert set(written) == {"run"}
    assert written["run"].exists()
    assert not (out / "openms.sample.parquet").exists()


def test_pg_peptide_counts_are_per_protein(monkeypatch):
    from qpx.converters.openms_consensus import pg_adapter

    class Header:
        filename = "run_01.mzML"
        label = ""

    class ConsensusMap:
        @staticmethod
        def getColumnHeaders():
            return {0: Header()}

        @staticmethod
        def getExperimentType():
            return "label-free"

        @staticmethod
        def getProteinIdentifications():
            return [object()]

    protein_maps = pg_adapter._ProteinMaps()
    protein_maps.acc_to_pep.update({"P1": {"PEPA", "PEPB"}, "P2": {"PEPB"}})
    protein_maps.acc_to_runs.update({"P1": {"run_01"}, "P2": {"run_01"}})
    protein_maps.acc_to_feat.update({"P1": {("PEPA", 2), ("PEPB", 2)}, "P2": {("PEPB", 2)}})
    protein_maps.pep_to_accs.update({"PEPA": {"P1"}, "PEPB": {"P1", "P2"}})
    protein_maps.feat_to_accs.update({("PEPA", 2): {"P1"}, ("PEPB", 2): {"P1", "P2"}})
    monkeypatch.setattr(pg_adapter, "_protein_maps", lambda _cm: protein_maps)
    monkeypatch.setattr(
        pg_adapter,
        "_merge_protein_ids",
        lambda _cm: ({}, {}, {}, [["P1", "P2"]]),
    )
    monkeypatch.setattr(pg_adapter, "_peptide_intensities", lambda _cm, _map_info: {})

    records = pg_adapter.consensus_protein_groups_to_records(cm=ConsensusMap())

    assert records[0]["peptide_counts"]["unique_sequences"] == 2
    assert records[0]["peptides"] == [
        {"protein_name": "P1", "peptide_count": 2},
        {"protein_name": "P2", "peptide_count": 1},
    ]


def test_pg_uses_id_merge_index_after_unmapped_map_index():
    """An unassigned PID with an unmapped map_index resolves its run via
    id_merge_index -> the i-th merged MS run (spectra_data order), not the
    map-column dict."""
    from qpx.converters.openms_consensus.pg_adapter import _protein_maps

    class Header:
        def __init__(self, filename):
            self.filename = filename

    class ProteinIdentification:
        @staticmethod
        def getPrimaryMSRunPath(output):
            # Merge order: id_merge_index 0 -> run_01, 1 -> run_02.
            output.extend([b"run_01.mzML", b"run_02.mzML"])

    class ConsensusMap:
        @staticmethod
        def getColumnHeaders():
            return {0: Header("run_01.mzML"), 1: Header("run_02.mzML")}

        @staticmethod
        def getProteinIdentifications():
            return [ProteinIdentification()]

        @staticmethod
        def getUnassignedPeptideIdentifications():
            return [PeptideIdentification()]

        def __iter__(self):
            return iter(())

    class Sequence:
        @staticmethod
        def toUnmodifiedString():
            return "PEPTIDE"

        @staticmethod
        def toUniModString():
            return "PEPTIDE"

    class Evidence:
        @staticmethod
        def getProteinAccession():
            return "P1"

    class Hit:
        @staticmethod
        def getSequence():
            return Sequence()

        @staticmethod
        def getCharge():
            return 2

        @staticmethod
        def getPeptideEvidences():
            return [Evidence()]

    class PeptideIdentification:
        values = {"map_index": 99, "id_merge_index": 1}

        def metaValueExists(self, key):
            return key in self.values

        def getMetaValue(self, key):
            return self.values[key]

        @staticmethod
        def getHits():
            return [Hit()]

    protein_maps = _protein_maps(ConsensusMap())
    assert protein_maps.acc_to_pep == {"P1": {"PEPTIDE"}}
    assert protein_maps.acc_to_runs == {"P1": {"run_02"}}
    assert protein_maps.acc_to_feat == {"P1": {("PEPTIDE", 2)}}


def test_consensus_psms_use_parent_feature_run_context(tmp_path):
    """Assigned PSMs inherit the run represented by their parent feature."""
    from qpx.converters.openms_consensus.psm_adapter import consensus_psms_to_records

    consensusxml = _TMT_CONSENSUSXML.replace(
        '<map id="0" name="run_01.mzML"',
        '<map id="0" name="run_A.mzML"',
    ).replace(
        '<map id="1" name="run_01.mzML"',
        '<map id="1" name="run_B.mzML"',
    )
    consensusxml = consensusxml.replace(
        'map="0" id="0" rt="100.123456" mz="450.251234" it="1000.0"', 'map="0" id="0" rt="100.123456" mz="450.251234" it="0.0"'
    )
    consensusxml = consensusxml.replace(
        "        </PeptideHit>\n      </PeptideIdentification>",
        '        </PeptideHit>\n        <UserParam type="int" name="id_merge_index" value="0"/>\n      </PeptideIdentification>',
    )
    path = tmp_path / "multi_run.consensusXML"
    path.write_text(consensusxml)

    records = consensus_psms_to_records(str(path))

    assert len(records) == 1
    assert records[0]["run_file_name"] == "run_B"


# A merged 2-plex isobaric consensusXML: two input MS runs (plexA, plexB), each
# contributing two TMT channels (maps 0-1 = plexA, 2-3 = plexB). The merged
# ProteinIdentification records the merge order as ``spectra_data`` = [plexA, plexB],
# so an unassigned PID's ``id_merge_index`` selects the i-th run. The unassigned
# PID below carries ``id_merge_index=1`` -> it belongs to plexB (the SECOND run),
# NOT plexA — the regression for issue #243.
_TWO_PLEX_ISOBARIC_CONSENSUSXML = """<?xml version="1.0" encoding="ISO-8859-1"?>
<consensusXML version="1.7" experiment_type="label-free"
  xsi:noNamespaceSchemaLocation="https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/ConsensusXML_1_7.xsd"
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <IdentificationRun id="PI_0" date="0000-00-00T00:00:00" search_engine="" search_engine_version="">
    <SearchParameters db="" db_version="" taxonomy="" mass_type="monoisotopic" charges=""
      enzyme="unknown_enzyme" missed_cleavages="0" precursor_peak_tolerance="0"
      precursor_peak_tolerance_ppm="false" peak_mass_tolerance="0"
      peak_mass_tolerance_ppm="false">
    </SearchParameters>
    <ProteinIdentification score_type="" higher_score_better="true" significance_threshold="0">
      <ProteinHit id="PH_0" accession="P12345" score="0" sequence=""></ProteinHit>
      <UserParam type="stringList" name="spectra_data" value="[plexA.mzML,plexB.mzML]"/>
    </ProteinIdentification>
  </IdentificationRun>
  <mapList count="4">
    <map id="0" name="plexA.mzML" unique_id="1" label="tmt6plex_126" size="1"></map>
    <map id="1" name="plexA.mzML" unique_id="2" label="tmt6plex_127" size="1"></map>
    <map id="2" name="plexB.mzML" unique_id="3" label="tmt6plex_126" size="1"></map>
    <map id="3" name="plexB.mzML" unique_id="4" label="tmt6plex_127" size="1"></map>
  </mapList>
  <consensusElementList>
    <consensusElement id="e_0" quality="0.0" charge="2">
      <centroid rt="100.0" mz="450.25" it="0.0"/>
      <groupedElementList>
        <element map="0" id="0" rt="100.0" mz="450.25" it="1000.0"/>
        <element map="1" id="1" rt="100.0" mz="450.25" it="2000.0"/>
      </groupedElementList>
      <PeptideIdentification identification_run_ref="PI_0" score_type="" higher_score_better="true"
        significance_threshold="0" MZ="450.26" RT="100"
        spectrum_reference="controllerType=0 controllerNumber=1 scan=42">
        <PeptideHit score="0" sequence="PEPTIDEK" charge="2" protein_refs="PH_0">
          <UserParam type="string" name="target_decoy" value="target"/>
        </PeptideHit>
      </PeptideIdentification>
    </consensusElement>
  </consensusElementList>
  <UnassignedPeptideIdentification identification_run_ref="PI_0" score_type="" higher_score_better="true"
    significance_threshold="0" MZ="500.30" RT="200"
    spectrum_reference="controllerType=0 controllerNumber=1 scan=77">
    <PeptideHit score="0" sequence="ELVISLIVEK" charge="2" protein_refs="PH_0">
      <UserParam type="string" name="target_decoy" value="target"/>
    </PeptideHit>
    <UserParam type="int" name="id_merge_index" value="1"/>
  </UnassignedPeptideIdentification>
</consensusXML>
"""


def test_unassigned_isobaric_psm_run_from_id_merge_index(tmp_path):
    """An unassigned PSM in a merged multi-run isobaric consensusXML is attributed
    to the run its id_merge_index selects (the SECOND run), not the first map's
    run — the regression for issue #243."""
    from qpx.converters.openms_consensus.psm_adapter import consensus_psms_to_records

    path = tmp_path / "two_plex.consensusXML"
    path.write_text(_TWO_PLEX_ISOBARIC_CONSENSUSXML)

    records = consensus_psms_to_records(str(path))

    unassigned = [r for r in records if r["peptidoform"] == "ELVISLIVEK"]
    assert len(unassigned) == 1
    # id_merge_index=1 -> spectra_data[1] = plexB. The bug attributed it to plexA.
    assert unassigned[0]["run_file_name"] == "plexB"
    assert list(unassigned[0]["scan"]) == [77]


def test_unassigned_isobaric_pg_run_from_id_merge_index(tmp_path):
    """The corrected run flows into the pg adapter's acc_to_runs for an unassigned
    PID (via accumulate_unassigned_maps)."""
    from qpx.converters.openms_consensus.feature_adapter import load_consensus_map
    from qpx.converters.openms_consensus.pg_adapter import _protein_maps

    path = tmp_path / "two_plex_pg.consensusXML"
    path.write_text(_TWO_PLEX_ISOBARIC_CONSENSUSXML)

    m = _protein_maps(load_consensus_map(str(path)))
    # ELVISLIVEK (unassigned, id_merge_index=1) contributes plexB, not plexA.
    assert "plexB" in m.acc_to_runs["P12345"]
    assert "plexA" in m.acc_to_runs["P12345"]  # from the assigned PEPTIDEK feature


def test_unassigned_isobaric_run_streaming_matches_pyopenms(tmp_path):
    """The streaming reader resolves the unassigned run identically to pyopenms."""
    import json

    import pyarrow.parquet as pq

    cx = tmp_path / "two_plex_stream.consensusXML"
    cx.write_text(_TWO_PLEX_ISOBARIC_CONSENSUSXML)

    def convert(streaming):
        out = tmp_path / ("stream" if streaming else "pyopenms")
        return OpenMSConsensusConverter().convert(
            str(cx), str(out), output_prefix="d", structures=("feature", "psm", "pg"), streaming=streaming
        )

    wp, ws = convert(False), convert(True)

    def canon(path):
        return sorted(json.dumps(r, sort_keys=True, default=str) for r in pq.read_table(str(path)).to_pylist())

    for view in ("feature", "psm", "pg"):
        assert canon(wp[view]) == canon(ws[view]), f"{view} differs between pyopenms and streaming"
    psm_accessions = pq.read_table(str(wp["psm"]), columns=["protein_accessions"]).column(0).to_pylist()
    assert psm_accessions and all(accessions == ["P12345"] for accessions in psm_accessions)


def test_consensus_psm_multihit_keeps_lowest_pep_and_merges_scores(tmp_path):
    """A PID with two hits colliding on the identity key resolves to one PSM: the
    lowest-PEP hit is kept and the other (engine's) search score is preserved."""
    from qpx.converters.openms_consensus.psm_adapter import consensus_psms_to_records

    # First hit: score 1.5, PEP 5.0e-02. Second hit (inserted): score 2.7, PEP
    # 1.0e-03 -> the second must win (lower PEP) even though it is not first.
    xml = _TMT_CONSENSUSXML.replace('score="0" sequence="PEPTIDEK"', 'score="1.5" sequence="PEPTIDEK"')
    xml = xml.replace('value="1.0e-03"', 'value="5.0e-02"')
    xml = xml.replace(
        "    </ProteinIdentification>",
        '      <ProteinHit id="PH_1" accession="P67890" score="0" sequence=""></ProteinHit>\n    </ProteinIdentification>',
    )
    second_hit = (
        '        <PeptideHit score="2.7" sequence="PEPTIDEK" charge="2" protein_refs="PH_1">\n'
        '          <UserParam type="string" name="target_decoy" value="target"/>\n'
        '          <UserParam type="float" name="Posterior Error Probability_score" value="1.0e-03"/>\n'
        "        </PeptideHit>"
    )
    xml = xml.replace(
        "        </PeptideHit>\n      </PeptideIdentification>",
        f"        </PeptideHit>\n{second_hit}\n      </PeptideIdentification>",
    )
    assert second_hit in xml  # guard: the fixture edit actually applied
    path = tmp_path / "multihit.consensusXML"
    path.write_text(xml)

    records = consensus_psms_to_records(str(path))

    assert len(records) == 1  # collapsed to a single PSM, not two
    rec = records[0]
    assert rec["peptidoform"] == "PEPTIDEK"
    assert rec["posterior_error_probability"] == pytest.approx(1.0e-03)  # lower PEP kept
    score_values = [s["score_value"] for s in (rec["additional_scores"] or [])]
    assert any(v == pytest.approx(2.7) for v in score_values)  # kept hit's own search score
    assert any(v == pytest.approx(1.5) for v in score_values)  # the other engine's search score preserved
    assert rec["protein_accessions"] == ["P12345", "P67890"]


def test_consensus_psm_distinct_peptidoforms_both_emitted(tmp_path):
    """Two *different* peptidoforms in one PID are distinct PSMs and both survive."""
    from qpx.converters.openms_consensus.psm_adapter import consensus_psms_to_records

    second_hit = (
        '        <PeptideHit score="0" sequence="PEPTIDER" charge="2" protein_refs="PH_0">\n'
        '          <UserParam type="string" name="target_decoy" value="target"/>\n'
        '          <UserParam type="float" name="Posterior Error Probability_score" value="2.0e-03"/>\n'
        "        </PeptideHit>"
    )
    xml = _TMT_CONSENSUSXML.replace(
        "        </PeptideHit>\n      </PeptideIdentification>",
        f"        </PeptideHit>\n{second_hit}\n      </PeptideIdentification>",
    )
    path = tmp_path / "twopep.consensusXML"
    path.write_text(xml)

    records = consensus_psms_to_records(str(path))

    assert {r["peptidoform"] for r in records} == {"PEPTIDEK", "PEPTIDER"}


def test_consensus_psm_sciex_nativeid_not_dropped(tmp_path):
    """A Sciex WIFF nativeID (no scan token, cycle-based) is kept, keyed by cycle."""
    from qpx.converters.openms_consensus.psm_adapter import consensus_psms_to_records

    xml = _TMT_CONSENSUSXML.replace(
        'spectrum_reference="controllerType=0 controllerNumber=1 scan=42"',
        'spectrum_reference="sample=1 period=1 cycle=123 experiment=2"',
    )
    path = tmp_path / "sciex.consensusXML"
    path.write_text(xml)

    records = consensus_psms_to_records(str(path))

    assert len(records) == 1  # not silently dropped
    assert list(records[0]["scan"]) == [123]  # cycle is the scan-equivalent ordinal


def test_consensus_psm_unknown_nativeid_uses_surrogate_scan(tmp_path):
    """A nativeID with no recognizable ordinal falls back to a deterministic
    int32 surrogate rather than being dropped."""
    from qpx.converters.openms_consensus.psm_adapter import consensus_psms_to_records

    xml = _TMT_CONSENSUSXML.replace(
        'spectrum_reference="controllerType=0 controllerNumber=1 scan=42"',
        'spectrum_reference="sample=1 period=1 experiment=2"',
    )
    path = tmp_path / "exotic.consensusXML"
    path.write_text(xml)

    records = consensus_psms_to_records(str(path))
    repeated = consensus_psms_to_records(str(path))

    assert len(records) == 1  # not dropped
    scan = list(records[0]["scan"])
    assert len(scan) == 1 and 0 <= scan[0] <= 0x7FFFFFFF  # deterministic, int32-safe
    assert repeated[0]["scan"] == records[0]["scan"]


def test_streaming_retains_all_protein_identification_run_paths(tmp_path):
    """Multiple ProteinIdentification blocks contribute to id_merge_index in
    document order; the streaming reader must not retain only the last block."""
    from qpx.converters.openms_consensus.psm_adapter import _merge_index_runs
    from qpx.converters.openms_consensus.streaming import StreamingConsensusMap

    xml = """<?xml version="1.0" encoding="UTF-8"?>
<consensusXML>
  <IdentificationRun id="PI_0">
    <ProteinIdentification score_type="">
      <ProteinHit id="PH_0" accession="P1" score="0"/>
      <UserParam type="stringList" name="spectra_data" value="[run_A.mzML]"/>
    </ProteinIdentification>
  </IdentificationRun>
  <IdentificationRun id="PI_1">
    <ProteinIdentification score_type="">
      <ProteinHit id="PH_1" accession="P2" score="0"/>
      <UserParam type="stringList" name="spectra_data" value="[run_B.mzML]"/>
    </ProteinIdentification>
  </IdentificationRun>
  <mapList count="0"/>
  <consensusElementList/>
</consensusXML>
"""
    path = tmp_path / "multiple-identification-runs.consensusXML"
    path.write_text(xml)

    consensus_map = StreamingConsensusMap(str(path))
    assert len(consensus_map.getProteinIdentifications()) == 2
    assert _merge_index_runs(consensus_map) == ["run_A", "run_B"]


def _write_multi_reference_consensusxml(path):
    """Write one ConsensusFeature supported by two spectrum references."""
    second_pid = """
      <PeptideIdentification identification_run_ref="PI_0" score_type=""
        higher_score_better="true" significance_threshold="0" MZ="450.26" RT="100"
        spectrum_reference="controllerType=0 controllerNumber=1 scan=43">
        <PeptideHit score="0" sequence="PEPTIDEK" charge="2" protein_refs="PH_0">
          <UserParam type="string" name="target_decoy" value="target"/>
        </PeptideHit>
      </PeptideIdentification>"""
    path.write_text(_TMT_CONSENSUSXML.replace("\n    </consensusElement>", f"{second_pid}\n    </consensusElement>"))


def _write_multirun_confidence_consensusxml(path):
    """Write one ConsensusFeature with a different confidence in each run."""
    xml = _TMT_CONSENSUSXML.replace(
        'name="run_01.mzML" unique_id="2" label="tmt6plex_127"',
        'name="run_02.mzML" unique_id="2" label="label-free"',
    ).replace('label="tmt6plex_126"', 'label="label-free"')
    xml = xml.replace(
        'PeptideIdentification identification_run_ref="PI_0" score_type=""',
        'PeptideIdentification identification_run_ref="PI_0" score_type="q-value"',
        1,
    ).replace(
        '<PeptideHit score="0" sequence="PEPTIDEK"',
        '<PeptideHit score="1.0e-03" sequence="PEPTIDEK"',
        1,
    )
    first_pid_end = """        </PeptideHit>
      </PeptideIdentification>"""
    pids = """        </PeptideHit>
        <UserParam type="int" name="map_index" value="0"/>
      </PeptideIdentification>
      <PeptideIdentification identification_run_ref="PI_0" score_type="q-value"
        higher_score_better="false" significance_threshold="0" MZ="450.26" RT="100"
        spectrum_reference="controllerType=0 controllerNumber=1 scan=43">
        <PeptideHit score="2.0e-02" sequence="PEPTIDEK" charge="2" protein_refs="PH_0">
          <UserParam type="string" name="target_decoy" value="target"/>
          <UserParam type="float" name="Posterior Error Probability_score" value="2.0e-02"/>
        </PeptideHit>
        <UserParam type="int" name="map_index" value="1"/>
      </PeptideIdentification>"""
    path.write_text(xml.replace(first_pid_end, pids, 1))


def test_streaming_matches_pyopenms(tmp_path):
    """The low-memory streaming reader produces the same parquet as pyopenms."""
    import json

    import pyarrow.parquet as pq

    cx = tmp_path / "test.consensusXML"
    _write_tmt_consensusxml(cx)

    def convert(streaming):
        out = tmp_path / ("stream" if streaming else "pyopenms")
        return OpenMSConsensusConverter().convert(
            str(cx), str(out), output_prefix="d", structures=("feature", "psm", "pg"), streaming=streaming
        )

    wp, ws = convert(False), convert(True)

    def canon(path):
        return sorted(json.dumps(r, sort_keys=True, default=str) for r in pq.read_table(str(path)).to_pylist())

    for view in ("feature", "psm", "pg"):
        assert canon(wp[view]) == canon(ws[view]), f"{view} differs between pyopenms and streaming"


def test_channel_sdrf_consistency_check(tmp_path):
    """Channels read from the consensusXML maps are checked against SDRF comment[label]."""
    from qpx.converters.openms_consensus.feature_adapter import check_channels_vs_sdrf, load_consensus_map

    cx = tmp_path / "test.consensusXML"
    _write_tmt_consensusxml(cx)  # channels: TMT126, TMT127
    cm = load_consensus_map(str(cx))

    # Matching SDRF -> no warnings.
    ok = tmp_path / "ok.sdrf.tsv"
    ok.write_text("comment[label]\nTMT126\nTMT127\n")
    assert check_channels_vs_sdrf(cm, str(ok)) == []

    # Mismatched SDRF -> flags both directions (TMT127 only in consensusXML,
    # TMT131 only in the SDRF).
    bad = tmp_path / "bad.sdrf.tsv"
    bad.write_text("comment[label]\nTMT126\nTMT131\n")
    msgs = check_channels_vs_sdrf(cm, str(bad))
    assert any("TMT127" in m and "consensusXML but not" in m for m in msgs)
    assert any("TMT131" in m and "SDRF comment[label] but not" in m for m in msgs)


def test_consensusxml_to_qpx_feature_has_channels_and_interim_pg_intensity(tmp_path):
    cx = tmp_path / "test.consensusXML"
    _write_tmt_consensusxml(cx)
    out = tmp_path / "out"
    written = OpenMSConsensusConverter().convert(str(cx), str(out), output_prefix="t", structures=("feature", "psm", "pg"))
    con = duckdb.connect()

    feat = con.execute(
        "SELECT peptidoform, charge, run_file_name, intensities FROM read_parquet($1)",
        [str(written["feature"])],
    ).fetchall()
    assert len(feat) == 1
    pep, charge, run, intensities = feat[0]
    assert pep == "PEPTIDEK" and charge == 2 and run == "run_01"
    assert len(intensities) == 2  # exactly the 2 channels, no duplicate rows
    labels = {e["label"]: e["intensity"] for e in intensities}
    assert labels == {"TMT126": 1000.0, "TMT127": 2000.0}  # both channels, canonicalized, quant kept

    psm = con.execute(
        "SELECT peptidoform, scan FROM read_parquet($1)",
        [str(written["psm"])],
    ).fetchall()
    assert len(psm) == 1  # one spectrum match, not collapsed/duplicated
    assert psm[0][0] == "PEPTIDEK" and list(psm[0][1]) == [42]

    pg = con.execute(
        "SELECT anchor_protein, label, intensity, cv_params, peptide_counts, feature_counts FROM read_parquet($1)",
        [str(written["pg"])],
    ).fetchall()
    assert len(pg) == 2  # one row per channel, no duplicate protein-group rows
    assert all(anchor == "P12345" for anchor, *_ in pg)
    # interim: one row per channel; intensity is the unnormalized sum of the
    # group's unique peptides for that channel (PEPTIDEK is unique to P12345, so
    # the protein total == its per-channel feature intensity).
    by_label = {row[1]: row[2] for row in pg}
    assert by_label == {"TMT126": 1000.0, "TMT127": 2000.0}
    for _, _, intensity, cv, pep_counts, feat_counts in pg:
        names = {p["cv_name"]: p["cv_value"] for p in (cv or [])}
        assert intensity is not None and names.get("quantification_method") == "unnormalized_unique_peptide_sum"
        # PEPTIDEK is unique to the single-protein group -> unique == total == 1.
        assert pep_counts == {"unique_sequences": 1, "total_sequences": 1}
        assert feat_counts == {"unique_features": 1, "total_features": 1}


def test_feature_pg_accessions_populated_full_membership(tmp_path):
    """feature.pg_accessions carries the full protein-group membership (not just the
    leader), and matches the pg view's pg_accessions for that group (bigbio/qpx#266)."""
    cx = tmp_path / "test.consensusXML"
    _write_tmt_consensusxml(cx)
    out = tmp_path / "out"
    written = OpenMSConsensusConverter().convert(str(cx), str(out), output_prefix="t", structures=("feature", "pg"))
    con = duckdb.connect()
    feat = con.execute(
        "SELECT anchor_protein, pg_accessions FROM read_parquet($1)",
        [str(written["feature"])],
    ).fetchall()
    pg = con.execute(
        "SELECT pg_accessions FROM read_parquet($1)",
        [str(written["pg"])],
    ).fetchall()
    con.close()

    assert len(feat) == 1
    anchor, pg_accessions = feat[0]
    assert pg_accessions is not None  # populated, not omitted
    accs = [p["accession"] for p in pg_accessions]
    assert accs == ["P12345"] and anchor == "P12345"  # the full membership, incl. the leader
    # the feature's membership equals the pg view's membership for that group
    assert accs == list(pg[0][0])


# Two DISTINCT indistinguishable groups that SHARE a leading protein A: [A, B] and
# [A, C] (the bigbio/qpx#240 class). Each peptide/feature maps to one group (its
# first protein evidence is the distinguishing member B / C). Anchor alone (A) does
# not identify the group; the full pg_accessions membership does.
_SHARED_LEADER_CONSENSUSXML = """<?xml version="1.0" encoding="ISO-8859-1"?>
<consensusXML version="1.7" experiment_type="label-free"
  xsi:noNamespaceSchemaLocation="https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/ConsensusXML_1_7.xsd"
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <IdentificationRun id="PI_0" date="0000-00-00T00:00:00" search_engine="" search_engine_version="">
    <SearchParameters db="" db_version="" taxonomy="" mass_type="monoisotopic" charges=""
      enzyme="unknown_enzyme" missed_cleavages="0" precursor_peak_tolerance="0"
      precursor_peak_tolerance_ppm="false" peak_mass_tolerance="0"
      peak_mass_tolerance_ppm="false">
    </SearchParameters>
    <ProteinIdentification score_type="" higher_score_better="true" significance_threshold="0">
      <ProteinHit id="PH_0" accession="A" score="0" sequence=""></ProteinHit>
      <ProteinHit id="PH_1" accession="B" score="0" sequence=""></ProteinHit>
      <ProteinHit id="PH_2" accession="C" score="0" sequence=""></ProteinHit>
      <UserParam type="string" name="indistinguishable_proteins_0" value="0,PH_0,PH_1"/>
      <UserParam type="string" name="indistinguishable_proteins_1" value="0,PH_0,PH_2"/>
    </ProteinIdentification>
  </IdentificationRun>
  <mapList count="1">
    <map id="0" name="run_01.mzML" unique_id="1" label="label-free" size="2">
    </map>
  </mapList>
  <consensusElementList>
    <consensusElement id="e_0" quality="0.0" charge="2">
      <centroid rt="100.0" mz="450.25" it="0.0"/>
      <groupedElementList>
        <element map="0" id="0" rt="100.0" mz="450.25" it="1000.0"/>
      </groupedElementList>
      <PeptideIdentification identification_run_ref="PI_0" score_type=""
        higher_score_better="true" significance_threshold="0" MZ="450.26" RT="100"
        spectrum_reference="scan=42">
        <PeptideHit score="0" sequence="PEPTIDEK" charge="2" protein_refs="PH_1 PH_0">
          <UserParam type="string" name="target_decoy" value="target"/>
        </PeptideHit>
      </PeptideIdentification>
    </consensusElement>
    <consensusElement id="e_1" quality="0.0" charge="2">
      <centroid rt="200.0" mz="500.25" it="0.0"/>
      <groupedElementList>
        <element map="0" id="1" rt="200.0" mz="500.25" it="3000.0"/>
      </groupedElementList>
      <PeptideIdentification identification_run_ref="PI_0" score_type=""
        higher_score_better="true" significance_threshold="0" MZ="500.26" RT="200"
        spectrum_reference="scan=43">
        <PeptideHit score="0" sequence="ELVISLIVK" charge="2" protein_refs="PH_2 PH_0">
          <UserParam type="string" name="target_decoy" value="target"/>
        </PeptideHit>
      </PeptideIdentification>
    </consensusElement>
  </consensusElementList>
</consensusXML>
"""


@pytest.mark.parametrize("streaming", [False, True])
def test_feature_pg_accessions_disambiguate_shared_leader(tmp_path, streaming):
    """Two distinct groups sharing a leader ([A,B] and [A,C]) give the two features
    DISTINCT feature.pg_accessions even though both share anchor_protein A, and each
    feature's pg_accessions matches its pg row's pg_accessions (bigbio/qpx#266)."""
    cx = tmp_path / "shared_leader.consensusXML"
    cx.write_text(_SHARED_LEADER_CONSENSUSXML)
    out = tmp_path / ("stream" if streaming else "mem")
    written = OpenMSConsensusConverter().convert(
        str(cx), str(out), output_prefix="t", structures=("feature", "pg"), streaming=streaming
    )
    con = duckdb.connect()
    feat = con.execute(
        "SELECT sequence, anchor_protein, pg_accessions FROM read_parquet($1)",
        [str(written["feature"])],
    ).fetchall()
    pg = con.execute(
        "SELECT anchor_protein, pg_accessions FROM read_parquet($1)",
        [str(written["pg"])],
    ).fetchall()
    con.close()

    feat_membership = {seq: [p["accession"] for p in (pg_accessions or [])] for seq, _, pg_accessions in feat}
    feat_anchor = {seq: anchor for seq, anchor, _ in feat}
    # Both features share the leader as anchor_protein ...
    assert feat_anchor == {"PEPTIDEK": "A", "ELVISLIVK": "A"}
    # ... but pg_accessions distinguishes their groups (the #240 ambiguity resolved).
    assert feat_membership["PEPTIDEK"] == ["A", "B"]
    assert feat_membership["ELVISLIVK"] == ["A", "C"]
    assert feat_membership["PEPTIDEK"] != feat_membership["ELVISLIVK"]

    # Each feature's membership matches exactly one pg row's pg_accessions.
    pg_memberships = [list(members) for _, members in pg]
    assert sorted(pg_memberships) == [["A", "B"], ["A", "C"]]
    assert feat_membership["PEPTIDEK"] in pg_memberships
    assert feat_membership["ELVISLIVK"] in pg_memberships


# A single protein group P shared across TWO label-free runs (run_01, run_02) so the
# pg row's grouped_runs is multi-run ([run_01, run_02]); plus a peptide with NO protein
# evidence (protein_refs="") whose feature has no protein group and therefore no pg row.
_MULTIRUN_CONSENSUSXML = """<?xml version="1.0" encoding="ISO-8859-1"?>
<consensusXML version="1.7" experiment_type="label-free"
  xsi:noNamespaceSchemaLocation="https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/ConsensusXML_1_7.xsd"
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <IdentificationRun id="PI_0" date="0000-00-00T00:00:00" search_engine="" search_engine_version="">
    <SearchParameters db="" db_version="" taxonomy="" mass_type="monoisotopic" charges=""
      enzyme="unknown_enzyme" missed_cleavages="0" precursor_peak_tolerance="0"
      precursor_peak_tolerance_ppm="false" peak_mass_tolerance="0" peak_mass_tolerance_ppm="false">
    </SearchParameters>
    <ProteinIdentification score_type="" higher_score_better="true" significance_threshold="0">
      <ProteinHit id="PH_0" accession="P99999" score="0" sequence=""></ProteinHit>
    </ProteinIdentification>
  </IdentificationRun>
  <mapList count="2">
    <map id="0" name="run_01.mzML" unique_id="1" label="label-free" size="1"></map>
    <map id="1" name="run_02.mzML" unique_id="2" label="label-free" size="1"></map>
  </mapList>
  <consensusElementList>
    <consensusElement id="e_0" quality="0.0" charge="2">
      <centroid rt="100.0" mz="450.25" it="0.0"/>
      <groupedElementList><element map="0" id="0" rt="100.0" mz="450.25" it="1000.0"/></groupedElementList>
      <PeptideIdentification identification_run_ref="PI_0" score_type="" higher_score_better="true"
        significance_threshold="0" MZ="450.26" RT="100" spectrum_reference="scan=1">
        <PeptideHit score="0" sequence="PEPTIDEK" charge="2" protein_refs="PH_0">
          <UserParam type="string" name="target_decoy" value="target"/>
        </PeptideHit>
      </PeptideIdentification>
    </consensusElement>
    <consensusElement id="e_1" quality="0.0" charge="2">
      <centroid rt="110.0" mz="450.25" it="0.0"/>
      <groupedElementList><element map="1" id="1" rt="110.0" mz="450.25" it="2000.0"/></groupedElementList>
      <PeptideIdentification identification_run_ref="PI_0" score_type="" higher_score_better="true"
        significance_threshold="0" MZ="450.26" RT="110" spectrum_reference="scan=2">
        <PeptideHit score="0" sequence="PEPTIDEK" charge="2" protein_refs="PH_0">
          <UserParam type="string" name="target_decoy" value="target"/>
        </PeptideHit>
      </PeptideIdentification>
    </consensusElement>
    <consensusElement id="e_2" quality="0.0" charge="2">
      <centroid rt="200.0" mz="600.25" it="0.0"/>
      <groupedElementList><element map="0" id="2" rt="200.0" mz="600.25" it="500.0"/></groupedElementList>
      <PeptideIdentification identification_run_ref="PI_0" score_type="" higher_score_better="true"
        significance_threshold="0" MZ="600.26" RT="200" spectrum_reference="scan=3">
        <PeptideHit score="0" sequence="NOPROTEINR" charge="2" protein_refs="">
          <UserParam type="string" name="target_decoy" value="target"/>
        </PeptideHit>
      </PeptideIdentification>
    </consensusElement>
  </consensusElementList>
</consensusXML>
"""


from tests.converters.pg_ids_roundtrip import assert_softlink_valid, open_converted, pg_ids_by_sequence  # noqa: E402


@pytest.mark.parametrize("streaming", [False, True])
def test_feature_pg_softlink_tmt_multilabel(tmp_path, streaming):
    """The computed softlink links each feature to real pg rows with matching
    pg_accessions + run + a carried label; a TMT feature (two channels) links to
    one pg row per channel/label, never over-linking (bigbio/qpx#269)."""
    cx = tmp_path / "tmt.consensusXML"
    _write_tmt_consensusxml(cx)
    out = tmp_path / ("stream" if streaming else "mem")
    OpenMSConsensusConverter().convert(str(cx), str(out), output_prefix="t", structures=("feature", "pg"), streaming=streaming)
    with open_converted(out, prefix="t") as ds:
        feat, pg, link = assert_softlink_valid(ds)

    assert len(feat) == 1  # one run, channels folded into intensities
    fid = feat[0][0]
    linked = {(pid, label) for f, pid, label in link if f == fid}
    # TMT126 + TMT127 -> one pg row per channel, so the feature links to BOTH.
    assert len(linked) == 2
    assert {label for _pid, label in linked} == {"TMT126", "TMT127"}
    assert {pid for pid, _label in linked} == {pg_id for pg_id, *_ in pg}


@pytest.mark.parametrize("streaming", [False, True])
def test_feature_pg_softlink_shared_leader_distinct(tmp_path, streaming):
    """Shared-leader ([A,B] vs [A,C]): the two features link to the CORRECT distinct
    pg rows (not each other's) — membership match is set-wise (bigbio/qpx#269)."""
    cx = tmp_path / "shared_leader.consensusXML"
    cx.write_text(_SHARED_LEADER_CONSENSUSXML)
    out = tmp_path / ("stream" if streaming else "mem")
    OpenMSConsensusConverter().convert(str(cx), str(out), output_prefix="t", structures=("feature", "pg"), streaming=streaming)
    with open_converted(out, prefix="t") as ds:
        feat, pg, link = assert_softlink_valid(ds)

    pg_id_by_membership = {tuple(sorted(memb)): pg_id for pg_id, memb, *_ in pg}
    ids_by_seq = pg_ids_by_sequence(feat, link)
    assert ids_by_seq["PEPTIDEK"] == [pg_id_by_membership[("A", "B")]]
    assert ids_by_seq["ELVISLIVK"] == [pg_id_by_membership[("A", "C")]]
    # The two features do NOT link to each other's pg row.
    assert set(ids_by_seq["PEPTIDEK"]).isdisjoint(ids_by_seq["ELVISLIVK"])


@pytest.mark.parametrize("streaming", [False, True])
def test_feature_pg_softlink_multirun_and_no_pg_row(tmp_path, streaming):
    """Multi-run grouped_runs: a feature's run is a member of a multi-run grouped_runs
    set, so both runs' features link to the same pg_id; and a feature whose group has
    NO pg row (no protein evidence) produces no softlink edge (bigbio/qpx#269)."""
    cx = tmp_path / "multirun.consensusXML"
    cx.write_text(_MULTIRUN_CONSENSUSXML)
    out = tmp_path / ("stream" if streaming else "mem")
    OpenMSConsensusConverter().convert(str(cx), str(out), output_prefix="t", structures=("feature", "pg"), streaming=streaming)
    with open_converted(out, prefix="t") as ds:
        feat, pg, link = assert_softlink_valid(ds)

    # One pg row: P99999 quantified over the two-run unit.
    assert len(pg) == 1
    pg_id, _memb, grouped_runs, _label = pg[0]
    assert sorted(grouped_runs) == ["run_01", "run_02"]  # multi-run grouped_runs

    # feature_id -> (sequence, run) so we can key the multi-run assertions.
    fmeta = {fid: (seq, run) for fid, seq, run, _memb, _labels in feat}
    ids_by_key: dict[tuple, set] = {}
    for fid, pid, _label in link:
        ids_by_key.setdefault(fmeta[fid], set()).add(pid)
    # Both runs' PEPTIDEK features link to the SAME single multi-run pg row.
    assert ids_by_key[("PEPTIDEK", "run_01")] == {pg_id}
    assert ids_by_key[("PEPTIDEK", "run_02")] == {pg_id}
    # The peptide with no protein evidence has no pg row -> no softlink edge.
    assert ("NOPROTEINR", "run_01") not in ids_by_key


def test_label_free_consensusxml_uses_lfq_labels(tmp_path):
    cx = tmp_path / "label_free.consensusXML"
    xml = _TMT_CONSENSUSXML.replace('label="tmt6plex_126"', 'label="label-free"')
    xml = xml.replace('label="tmt6plex_127"', 'label="label-free"')
    cx.write_text(xml)
    out = tmp_path / "out"
    written = OpenMSConsensusConverter().convert(
        str(cx),
        str(out),
        output_prefix="lfq",
        structures=("feature", "pg"),
    )
    con = duckdb.connect()
    feature_labels = con.execute(
        "SELECT UNNEST(intensities).label FROM read_parquet($1)",
        [str(written["feature"])],
    ).fetchall()
    pg_labels = con.execute(
        "SELECT label FROM read_parquet($1)",
        [str(written["pg"])],
    ).fetchall()
    con.close()

    assert {label for (label,) in feature_labels} == {"LFQ"}
    assert {label for (label,) in pg_labels} == {"LFQ"}


@pytest.mark.parametrize("streaming", [False, True])
def test_openms_consensus_cross_refs_resolve(tmp_path, streaming):
    """psm.feature_id is the authoritative FK: it is materialized and every value
    resolves to a real feature.feature_id (both the in-memory and streaming paths).
    feature.psm_ids is NOT materialized — the feature->psms mapping is recovered on
    read via Dataset.link_feature_psm() (bigbio/qpx#267)."""
    cx = tmp_path / "test.consensusXML"
    _write_tmt_consensusxml(cx)
    out = tmp_path / ("stream" if streaming else "mem")
    written = OpenMSConsensusConverter().convert(
        str(cx), str(out), output_prefix="t", structures=("feature", "psm"), streaming=streaming
    )
    con = duckdb.connect()

    feature_ids = {r[0] for r in con.execute(f"SELECT feature_id FROM read_parquet('{written['feature']}')").fetchall()}

    # psm.feature_id is the authoritative producer assignment: materialized, and
    # every non-null value resolves to a real feature (guards a silently null FK).
    linked = con.execute(
        f"SELECT psm_id, feature_id FROM read_parquet('{written['psm']}') WHERE feature_id IS NOT NULL"
    ).fetchall()
    assert linked, "expected at least one psm.feature_id to be populated"
    for _pid, fid in linked:
        assert fid in feature_ids, f"psm.feature_id {fid} does not resolve to a feature.feature_id"

    # feature.psm_ids is the computed inverse — qpx does NOT materialize it.
    psm_ids_col = con.execute(f"SELECT psm_ids FROM read_parquet('{written['feature']}')").fetchall()
    assert all(ids is None for (ids,) in psm_ids_col), "feature.psm_ids must not be materialized (computed on read)"
    con.close()

    # Dataset.link_feature_psm() is the inverse of psm.feature_id; grouping it by
    # feature_id recovers the feature->psms mapping we no longer persist.
    expected: dict[int, set[int]] = {}
    for pid, fid in linked:
        expected.setdefault(fid, set()).add(pid)
    with qpx.Dataset(str(out), file_prefix="t", structures=["feature", "psm"]) as ds:
        recovered: dict[int, set[int]] = {}
        for fid, pid in ds.link_feature_psm().fetchall():
            recovered.setdefault(fid, set()).add(pid)
    assert recovered == expected


@pytest.mark.parametrize("streaming", [False, True])
def test_openms_consensus_identity_retains_all_spectrum_references(tmp_path, streaming):
    """The Feature identity keeps every supporting scan and the parent consensus RT."""
    import pyarrow.parquet as pq

    cx = tmp_path / "multi.consensusXML"
    _write_multi_reference_consensusxml(cx)
    written = OpenMSConsensusConverter().convert(
        str(cx),
        str(tmp_path / ("stream" if streaming else "mem")),
        output_prefix="t",
        structures=("feature",),
        streaming=streaming,
    )
    table = pq.read_table(written["feature"])
    assert table.column("scan").to_pylist() == [[42, 43]]
    assert table.column("consensus_rt").to_pylist() == pytest.approx([100.123456])
    assert table.schema.metadata[b"identity_composite"] == (b"peptidoform,charge,run_file_name,rt,scan,observed_mz,consensus_rt")
