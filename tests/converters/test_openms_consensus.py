"""consensusXML -> QPX converter (interim path).

Reads a tiny literal consensusXML fixture (no pyopenms construction) so the
extraction runs in CI without a large fixture file and without depending on the
pyopenms setter APIs, which differ across versions.
"""

import duckdb
import pytest

pytest.importorskip("pyopenms")
import pyopenms as oms  # noqa: E402

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
    from qpx.converters.openms_consensus.pg_adapter import _protein_maps

    class Header:
        def __init__(self, filename):
            self.filename = filename

    class ConsensusMap:
        @staticmethod
        def getColumnHeaders():
            return {0: Header("run_01.mzML"), 1: Header("run_02.mzML")}

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


def test_consensus_psm_multihit_keeps_lowest_pep_and_merges_scores(tmp_path):
    """A PID with two hits colliding on the identity key resolves to one PSM: the
    lowest-PEP hit is kept and the other (engine's) search score is preserved."""
    from qpx.converters.openms_consensus.psm_adapter import consensus_psms_to_records

    # First hit: score 1.5, PEP 5.0e-02. Second hit (inserted): score 2.7, PEP
    # 1.0e-03 -> the second must win (lower PEP) even though it is not first.
    xml = _TMT_CONSENSUSXML.replace('score="0" sequence="PEPTIDEK"', 'score="1.5" sequence="PEPTIDEK"')
    xml = xml.replace('value="1.0e-03"', 'value="5.0e-02"')
    second_hit = (
        '        <PeptideHit score="2.7" sequence="PEPTIDEK" charge="2" protein_refs="PH_0">\n'
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

    assert len(records) == 1  # not dropped
    scan = list(records[0]["scan"])
    assert len(scan) == 1 and 0 <= scan[0] <= 0x7FFFFFFF  # deterministic, int32-safe


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
    """Every psm.feature_id and every feature.psm_ids element resolves to a real
    sibling id — the converter-populated FK matches the writer-derived PK (both
    the in-memory and streaming paths)."""
    cx = tmp_path / "test.consensusXML"
    _write_tmt_consensusxml(cx)
    out = tmp_path / ("stream" if streaming else "mem")
    written = OpenMSConsensusConverter().convert(
        str(cx), str(out), output_prefix="t", structures=("feature", "psm"), streaming=streaming
    )
    con = duckdb.connect()

    feature_ids = {r[0] for r in con.execute(f"SELECT feature_id FROM read_parquet('{written['feature']}')").fetchall()}
    psm_ids = {r[0] for r in con.execute(f"SELECT psm_id FROM read_parquet('{written['psm']}')").fetchall()}

    # At least one link was actually populated (guards against a silently null FK).
    linked_feature_ids = [
        r[0] for r in con.execute(f"SELECT feature_id FROM read_parquet('{written['psm']}')").fetchall() if r[0] is not None
    ]
    assert linked_feature_ids, "expected at least one psm.feature_id to be populated"
    for fid in linked_feature_ids:
        assert fid in feature_ids, f"psm.feature_id {fid} does not resolve to a feature.feature_id"

    referenced_psm_ids = [
        pid
        for (ids,) in con.execute(f"SELECT psm_ids FROM read_parquet('{written['feature']}')").fetchall()
        for pid in (ids or [])
    ]
    assert referenced_psm_ids, "expected at least one feature.psm_ids element to be populated"
    for pid in referenced_psm_ids:
        assert pid in psm_ids, f"feature.psm_ids element {pid} does not resolve to a psm.psm_id"


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
