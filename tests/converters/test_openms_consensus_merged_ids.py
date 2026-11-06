"""consensusXML layouts produced by FeatureFinderIdentification on group-merged IDs.

Covers run attribution of identification copies, unnamed single-run maps,
duplicate consensus features and multi-file input, on both readers.
"""

import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

pytest.importorskip("pyopenms")

import qpx  # noqa: E402
from qpx.cli.main import qpx_main  # noqa: E402
from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter  # noqa: E402
from qpx.converters.openms_consensus.feature_adapter import column_runs, load_consensus_map  # noqa: E402
from qpx.converters.openms_consensus.streaming import StreamingConsensusMap  # noqa: E402
from tests.converters.test_openms_consensus import (  # noqa: E402
    _MULTIRUN_CONSENSUSXML,
    _TWO_PLEX_ISOBARIC_CONSENSUSXML,
    _write_multirun_confidence_consensusxml,
)

_HEADER = """<?xml version="1.0" encoding="ISO-8859-1"?>
<consensusXML version="1.7" experiment_type="label-free">
"""


def _identification_run(run_id, spectra_data, proteins=("P1",)):
    hits = "".join(
        f'      <ProteinHit id="{run_id}_PH_{i}" accession="{acc}" score="0" sequence=""></ProteinHit>\n'
        for i, acc in enumerate(proteins)
    )
    spectra = f'      <UserParam type="stringList" name="spectra_data" value="[{spectra_data}]"/>\n' if spectra_data else ""
    return (
        f'  <IdentificationRun id="{run_id}" date="0000-00-00T00:00:00" search_engine="Percolator" search_engine_version="3">\n'
        '    <ProteinIdentification score_type="" higher_score_better="true" significance_threshold="0">\n'
        f"{hits}{spectra}"
        "    </ProteinIdentification>\n"
        "  </IdentificationRun>\n"
    )


def _pid(tag, run_ref, scan, hits, map_index=None, merge_index=None, rt=100.0):
    meta = ""
    if map_index is not None:
        meta += f'    <UserParam type="int" name="map_index" value="{map_index}"/>\n'
    if merge_index is not None:
        meta += f'    <UserParam type="int" name="id_merge_index" value="{merge_index}"/>\n'
    peptide_hits = "".join(
        f'    <PeptideHit score="0.9" sequence="{sequence}" charge="2" protein_refs="{run_ref}_PH_0">\n'
        '      <UserParam type="string" name="target_decoy" value="target"/>\n'
        "    </PeptideHit>\n"
        for sequence in hits
    )
    return (
        f'  <{tag} identification_run_ref="{run_ref}" score_type="q-value" higher_score_better="false"'
        f' significance_threshold="0" MZ="450.26" RT="{rt}" spectrum_reference="controllerType=0 controllerNumber=1 scan={scan}">\n'
        f"{peptide_hits}{meta}  </{tag}>\n"
    )


def _maps(names):
    rows = "".join(
        f'    <map id="{i}" name="{name}" unique_id="{i + 1}" label="" size="1"></map>\n' for i, name in enumerate(names)
    )
    return f'  <mapList count="{len(names)}">\n{rows}  </mapList>\n'


def _element(index, elements, pids, quality=1.0, mz=450.25, rt=100.0):
    grouped = "".join(f'      <element map="{m}" id="{index}{m}" rt="{e_rt}" mz="{mz}" it="{it}"/>\n' for m, e_rt, it in elements)
    return (
        f'  <consensusElement id="e_{index}" quality="{quality}" charge="2">\n'
        f'    <centroid rt="{rt}" mz="{mz}" it="0.0"/>\n'
        f"    <groupedElementList>\n{grouped}    </groupedElementList>\n"
        f"{''.join(pids)}"
        "  </consensusElement>\n"
    )


def _merged_copies_xml():
    """Two runs whose maps each carry a copy of every group-merged identification.

    Spectrum scan=10 comes from run_A, scan=20 and scan=30 from run_B (``id_merge_index``
    into the identical ``spectra_data`` of both ProteinIdentifications). Unassigned
    identifications precede the consensus features, as OpenMS writes them.
    """
    assigned = "PEPTIDEK"
    unassigned = [
        _pid("UnassignedPeptideIdentification", "PI_1", 10, [assigned], map_index=1, merge_index=0),
        _pid("UnassignedPeptideIdentification", "PI_0", 30, ["ELVISLIVEK"], map_index=0, merge_index=1, rt=300),
        _pid("UnassignedPeptideIdentification", "PI_1", 30, ["ELVISLIVEK"], map_index=1, merge_index=1, rt=300),
    ]
    element = _element(
        0,
        [(0, 99.0, 1000.0), (1, 101.0, 2000.0)],
        [
            _pid("PeptideIdentification", "PI_0", 10, [assigned], map_index=0, merge_index=0).replace("\n  ", "\n    "),
            _pid("PeptideIdentification", "PI_0", 20, [assigned], map_index=0, merge_index=1).replace("\n  ", "\n    "),
            _pid("PeptideIdentification", "PI_1", 20, [assigned], map_index=1, merge_index=1).replace("\n  ", "\n    "),
        ],
    )
    return (
        _HEADER
        + _identification_run("PI_0", "run_A.mzML,run_B.mzML")
        + _identification_run("PI_1", "run_A.mzML,run_B.mzML")
        + "".join(unassigned)
        + _maps(["run_A.mzML", "run_B.mzML"])
        + f"  <consensusElementList>\n{element}  </consensusElementList>\n"
        + "</consensusXML>\n"
    )


def _convert(tmp_path, source, streaming, structures=("feature", "psm"), name="out", **kwargs):
    out = tmp_path / f"{name}_{'stream' if streaming else 'memory'}"
    written = OpenMSConsensusConverter().convert(source, str(out), structures=structures, streaming=streaming, **kwargs)
    return out, written


def _rows(written, view):
    return pq.read_table(written[view]).to_pylist()


def _assert_valid(out, structures):
    with qpx.open_dataset(out) as dataset:
        results = dataset.validate(structures=list(structures), strict=True)
    invalid = {name: [issue.message for issue in result.issues] for name, result in results.items() if not result.is_valid}
    assert not invalid


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_merged_id_copies_collapse_to_their_spectrum_run(tmp_path, streaming):
    """Each spectrum yields one PSM in the run its id_merge_index names, linked to that run's feature."""
    source = tmp_path / "merged.consensusXML"
    source.write_text(_merged_copies_xml())

    out, written = _convert(tmp_path, str(source), streaming)

    psms = {(row["scan"][0], row["run_file_name"]): row for row in _rows(written, "psm")}
    features = {row["run_file_name"]: row for row in _rows(written, "feature")}
    assert set(psms) == {(10, "run_A"), (20, "run_B"), (30, "run_B")}
    assert features["run_A"]["scan"] == [10]
    assert features["run_B"]["scan"] == [20]
    assert psms[(10, "run_A")]["feature_id"] == features["run_A"]["feature_id"]
    assert psms[(20, "run_B")]["feature_id"] == features["run_B"]["feature_id"]
    assert psms[(30, "run_B")]["feature_id"] is None
    _assert_valid(out, ("feature", "psm"))


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_single_merged_identification_keeps_map_index(tmp_path, streaming):
    """One ProteinIdentification (quantms ProteomicsLFQ layout): map_index stays authoritative
    even where id_merge_index points elsewhere, as in OpenMS 3.1 ProteomicsLFQ output."""
    xml = _merged_copies_xml()
    start, end = xml.index('  <IdentificationRun id="PI_1"'), xml.index("  <UnassignedPeptideIdentification")
    source = tmp_path / "single_identification.consensusXML"
    source.write_text(xml[:start] + xml[end:].replace("PI_1", "PI_0"))

    _, written = _convert(tmp_path, str(source), streaming, structures=("psm",))

    assert sorted((row["scan"][0], row["run_file_name"]) for row in _rows(written, "psm")) == [
        (10, "run_A"),
        (10, "run_B"),
        (20, "run_A"),
        (20, "run_B"),
        (30, "run_A"),
        (30, "run_B"),
    ]


def _lfq_fixtures(tmp_path):
    """Existing label-free fixtures, plus quantms-style merged-run annotations of them."""
    confidence = tmp_path / "confidence.consensusXML"
    _write_multirun_confidence_consensusxml(confidence)
    base = confidence.read_text()
    fixtures = {"confidence": base, "multirun": _MULTIRUN_CONSENSUSXML, "isobaric": _TWO_PLEX_ISOBARIC_CONSENSUSXML}
    for label, order in [("agreeing", "run_01.mzML,run_02.mzML"), ("reversed", "run_02.mzML,run_01.mzML")]:
        xml = base.replace(
            "    </ProteinIdentification>",
            f'      <UserParam type="stringList" name="spectra_data" value="[{order}]"/>\n    </ProteinIdentification>',
        )
        # Each PID gets id_merge_index 0 and 1 in turn next to its map_index.
        for index in (0, 1):
            xml = xml.replace(
                f'<UserParam type="int" name="map_index" value="{index}"/>',
                f'<UserParam type="int" name="map_index" value="{index}"/>'
                f'<UserParam type="int" name="id_merge_index" value="{index}"/>',
                1,
            )
        fixtures[label] = xml
    return fixtures


# Run assignment of the fixtures above before identification copies were recognised.
_EXPECTED_LFQ_RUNS = {
    "confidence": [("PEPTIDEK", 42, "run_01"), ("PEPTIDEK", 43, "run_02")],
    "agreeing": [("PEPTIDEK", 42, "run_01"), ("PEPTIDEK", 43, "run_02")],
    "reversed": [("PEPTIDEK", 42, "run_01"), ("PEPTIDEK", 43, "run_02")],
    "multirun": [("NOPROTEINR", 3, "run_01"), ("PEPTIDEK", 1, "run_01"), ("PEPTIDEK", 2, "run_02")],
    "isobaric": [("ELVISLIVEK", 77, "plexB"), ("PEPTIDEK", 42, "plexA")],
}


@pytest.mark.parametrize("fixture", sorted(_EXPECTED_LFQ_RUNS))
@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_quantms_style_psm_runs_are_unchanged(tmp_path, fixture, streaming):
    """ProteomicsLFQ/IsobaricWorkflow-style maps keep their established PSM run assignment."""
    source = tmp_path / f"{fixture}.consensusXML"
    source.write_text(_lfq_fixtures(tmp_path)[fixture])

    _, written = _convert(tmp_path, str(source), streaming, structures=("feature", "psm"))
    runs = sorted((row["peptidoform"], row["scan"][0], row["run_file_name"]) for row in _rows(written, "psm"))

    assert runs == _EXPECTED_LFQ_RUNS[fixture]


def _single_unnamed_map_xml(spectra_data="IP_lung_1T.mzML"):
    element = _element(
        0, [(0, 100.0, 1000.0)], [_pid("PeptideIdentification", "PI_0", 7, ["PEPTIDEK"], map_index=0).replace("\n  ", "\n    ")]
    )
    return (
        _HEADER
        + _identification_run("PI_0", spectra_data)
        + _pid("UnassignedPeptideIdentification", "PI_0", 8, ["ELVISLIVEK"], rt=200)
        + _maps([""])
        + f"  <consensusElementList>\n{element}  </consensusElementList>\n"
        + "</consensusXML>\n"
    )


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_unnamed_single_map_takes_the_identification_run(tmp_path, streaming):
    """A FileConverter-promoted single-run map has no column filename; its run comes from spectra_data."""
    source = tmp_path / "single.consensusXML"
    source.write_text(_single_unnamed_map_xml())

    _, written = _convert(tmp_path, str(source), streaming, structures=("feature", "psm", "pg"))

    assert {row["run_file_name"] for row in _rows(written, "psm")} == {"IP_lung_1T"}
    assert {row["run_file_name"] for row in _rows(written, "feature")} == {"IP_lung_1T"}
    assert {tuple(row["grouped_runs"]) for row in _rows(written, "pg")} == {("IP_lung_1T",)}


@pytest.mark.parametrize("reader", [load_consensus_map, StreamingConsensusMap], ids=["memory", "streaming"])
def test_unnamed_map_stays_unnamed_without_a_single_identification_run(tmp_path, reader):
    source = tmp_path / "ambiguous.consensusXML"
    source.write_text(_single_unnamed_map_xml("run_A.mzML,run_B.mzML"))

    assert column_runs(reader(str(source))) == {0: ""}


def _duplicate_features_xml(better_first):
    """Two consensus features identical in every feature identity column; the lower-quality
    one also carries a second, distinct peptidoform for the same spectrum."""
    winner = _element(
        0,
        [(0, 100.0, 5000.0)],
        [_pid("PeptideIdentification", "PI_0", 7, ["VLFHNLPSL"], map_index=0).replace("\n  ", "\n    ")],
        quality=2.0,
    )
    loser = _element(
        1,
        [(0, 100.0, 1000.0)],
        [_pid("PeptideIdentification", "PI_0", 7, ["VLFHNLPSL", "AVSSFIFLR"], map_index=0).replace("\n  ", "\n    ")],
        quality=1.0,
    )
    elements = winner + loser if better_first else loser + winner
    return (
        _HEADER
        + _identification_run("PI_0", "run_01.mzML")
        + _maps(["run_01.mzML"])
        + f"  <consensusElementList>\n{elements}  </consensusElementList>\n"
        + "</consensusXML>\n"
    )


@pytest.mark.parametrize("better_first", [True, False], ids=["better-first", "better-last"])
@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_duplicate_consensus_features_keep_the_higher_quality_one(tmp_path, streaming, better_first, caplog):
    source = tmp_path / "duplicates.consensusXML"
    source.write_text(_duplicate_features_xml(better_first))

    out, written = _convert(tmp_path, str(source), streaming, structures=("feature", "psm", "pg"))

    features = _rows(written, "feature")
    assert len(features) == 1
    assert features[0]["intensities"][0]["intensity"] == 5000.0
    psms = _rows(written, "psm")
    assert sorted(row["peptidoform"] for row in psms) == ["AVSSFIFLR", "VLFHNLPSL"]
    assert {row["feature_id"] for row in psms} == {features[0]["feature_id"]}
    assert "Dropped 1 duplicate feature row(s)" in caplog.text
    _assert_valid(out, ("feature", "psm"))


def _renamed(xml, mapping):
    for old, new in mapping.items():
        xml = xml.replace(old, new)
    return xml


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_several_consensusxml_inputs_form_one_dataset(tmp_path, streaming):
    first, second = tmp_path / "group1.consensusXML", tmp_path / "group2.consensusXML"
    first.write_text(_merged_copies_xml())
    second.write_text(_renamed(_merged_copies_xml(), {"run_A": "run_C", "run_B": "run_D", '"P1"': '"P2"'}))

    out, written = _convert(tmp_path, [str(first), first.parent / second.name], streaming, structures=("feature", "psm", "pg"))

    runs = {row["run_file_name"] for row in _rows(written, "psm")}
    assert runs == {row["run_file_name"] for row in _rows(written, "feature")} == {"run_A", "run_B", "run_C", "run_D"}
    assert len(_rows(written, "psm")) == 6
    pg = _rows(written, "pg")
    assert {row["anchor_protein"] for row in pg} == {"P1", "P2"}
    assert {tuple(row["grouped_runs"]) for row in pg} == {("run_A",), ("run_B",), ("run_C",), ("run_D",)}
    provenance = pq.read_table(out / "openms.provenance.parquet").to_pylist()
    qpx_step = next(step for step in provenance if step["tool_name"] == "qpx")
    assert qpx_step["parameters"] == [
        {"key": "consensusxml", "value": "group1.consensusXML"},
        {"key": "consensusxml", "value": "group2.consensusXML"},
    ]
    _assert_valid(out, ("feature", "psm"))


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_single_input_matches_a_one_element_list(tmp_path, streaming):
    source = tmp_path / "merged.consensusXML"
    source.write_text(_merged_copies_xml())

    _, single = _convert(tmp_path, str(source), streaming, name="single")
    _, listed = _convert(tmp_path, [str(source)], streaming, name="listed")

    for view in ("feature", "psm"):
        assert _rows(single, view) == _rows(listed, view)


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_inputs_sharing_a_run_are_rejected(tmp_path, streaming):
    first, second = tmp_path / "group1.consensusXML", tmp_path / "group2.consensusXML"
    first.write_text(_merged_copies_xml())
    second.write_text(_renamed(_merged_copies_xml(), {"run_A": "run_C"}))
    out = tmp_path / "out"

    with pytest.raises(ValueError, match="'run_B' appears in both"):
        OpenMSConsensusConverter().convert(
            [str(first), str(second)], str(out), structures=("feature", "psm"), streaming=streaming
        )
    assert not list(out.glob("*.parquet"))


def test_cli_accepts_consensusxml_list(tmp_path):
    first, second = tmp_path / "group1.consensusXML", tmp_path / "group2.consensusXML"
    first.write_text(_merged_copies_xml())
    second.write_text(_renamed(_merged_copies_xml(), {"run_A": "run_C", "run_B": "run_D"}))
    out = tmp_path / "out"
    args = ["convert", "openms-consensus", "--structures", "psm", "--no-mudata", "--output-folder", str(out)]

    result = CliRunner().invoke(qpx_main, [*args, "--consensusxml", f"{first},{second}"])
    assert result.exit_code == 0, result.output
    assert len(pq.read_table(out / "openms.psm.parquet")) == 6

    rejected = CliRunner().invoke(qpx_main, [*args, "--consensusxml", f"{first},{first}"])
    assert rejected.exit_code != 0
    assert "appears in both" in rejected.output

    missing = CliRunner().invoke(qpx_main, [*args, "--consensusxml", f"{first},{tmp_path / 'absent.consensusXML'}"])
    assert missing.exit_code != 0
    assert "does not exist" in missing.output


def _removed_from(pid_xml, feature_uid):
    """Mark an unassigned identification as removed from a consensus feature, as IDConflictResolver does."""
    tag = "</UnassignedPeptideIdentification>"
    return pid_xml.replace(tag, f'    <UserParam type="string" name="feature_id" value="{feature_uid}"/>\n  {tag}')


def _conflict_resolved_xml():
    """IDConflictResolver kept run_A's identification on e_7; run_B's own one is unassigned.

    e_8 is a genuine transfer: run_B's feature carries no identification from run_B.
    """
    peptide = "PEPTIDEK"
    unassigned = [
        _removed_from(_pid("UnassignedPeptideIdentification", "PI_1", 20, [peptide], map_index=1, merge_index=1), 7),
        _removed_from(_pid("UnassignedPeptideIdentification", "PI_1", 40, ["ELVISLIVEK"], map_index=1, merge_index=1), 7),
    ]
    kept = _pid("PeptideIdentification", "PI_0", 10, [peptide], map_index=0, merge_index=0).replace("\n  ", "\n    ")
    transfer = _pid("PeptideIdentification", "PI_0", 50, ["SAMPLER"], map_index=0, merge_index=0, rt=500)
    elements = _element(7, [(0, 99.0, 1000.0), (1, 101.0, 2000.0)], [kept]) + _element(
        8, [(0, 499.0, 1000.0), (1, 501.0, 2000.0)], [transfer.replace("\n  ", "\n    ")], rt=500
    )
    return (
        _HEADER
        + _identification_run("PI_0", "run_A.mzML,run_B.mzML")
        + _identification_run("PI_1", "run_A.mzML,run_B.mzML")
        + "".join(unassigned)
        + _maps(["run_A.mzML", "run_B.mzML"])
        + f"  <consensusElementList>\n{elements}  </consensusElementList>\n"
        + "</consensusXML>\n"
    )


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_conflict_resolved_own_run_identification_is_recovered(tmp_path, streaming):
    """A run's own identification removed by IDConflictResolver still identifies that run's feature."""
    source = tmp_path / "resolved.consensusXML"
    source.write_text(_conflict_resolved_xml())

    out, written = _convert(tmp_path, str(source), streaming)

    features = {(row["peptidoform"], row["run_file_name"]): row for row in _rows(written, "feature")}
    recovered = features[("PEPTIDEK", "run_B")]
    assert recovered["id_run_file_name"] == "run_B"
    assert recovered["scan"] == [20]
    assert features[("PEPTIDEK", "run_A")]["scan"] == [10]
    transferred = features[("SAMPLER", "run_B")]
    assert transferred["id_run_file_name"] is None
    assert transferred["scan"] == []

    psms = {row["scan"][0]: row for row in _rows(written, "psm")}
    assert psms[20]["feature_id"] == recovered["feature_id"]
    assert psms[40]["feature_id"] is None
    _assert_valid(out, ("feature", "psm"))


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_no_include_unassigned_psms_excludes_recovered_source_ids(tmp_path, streaming):
    """The opt-out excludes recovered source-unassigned PSMs but keeps Feature metadata recovery."""
    source = tmp_path / "resolved-no-unassigned.consensusXML"
    source.write_text(_conflict_resolved_xml())

    _, written = _convert(tmp_path, str(source), streaming, include_unassigned_psms=False)

    psms = _rows(written, "psm")
    assert sorted(row["scan"][0] for row in psms) == [10, 50]
    features = {(row["peptidoform"], row["run_file_name"]): row for row in _rows(written, "feature")}
    assert features[("PEPTIDEK", "run_B")]["scan"] == [20]


def _percolator_run(run_id, spectra_data, peptide_level):
    """An identification run as PercolatorAdapter writes it, with its FDR level in the search parameters."""
    flag = "1" if peptide_level else "0"
    return _identification_run(run_id, spectra_data).replace(
        "    <ProteinIdentification",
        '    <SearchParameters db="db.fasta" db_version="" taxonomy="" mass_type="monoisotopic" charges="2:3"'
        ' enzyme="unspecific cleavage" missed_cleavages="0" precursor_peak_tolerance="0"'
        ' precursor_peak_tolerance_ppm="true" peak_mass_tolerance="0.01" peak_mass_tolerance_ppm="false" >\n'
        f'      <UserParam type="int" name="Percolator:peptide_level_fdrs" value="{flag}"/>\n'
        "    </SearchParameters>\n"
        "    <ProteinIdentification",
        1,
    )


def _xcorr_pid(tag, scan, sequence, qvalue, pep, rt=100.0, indent="  "):
    """A PSM whose primary score is COMET:xcorr, with Percolator q-value and PEP as hit meta values."""
    return (
        f'{indent}<{tag} identification_run_ref="PI_0" score_type="COMET:xcorr" higher_score_better="true"'
        f' significance_threshold="0" MZ="450.26" RT="{rt}" spectrum_reference="controllerType=0 controllerNumber=1 scan={scan}">\n'
        f'{indent}  <PeptideHit score="2.5" sequence="{sequence}" charge="2" protein_refs="PI_0_PH_0">\n'
        f'{indent}    <UserParam type="string" name="target_decoy" value="target"/>\n'
        f'{indent}    <UserParam type="float" name="q-value" value="{qvalue}"/>\n'
        f'{indent}    <UserParam type="float" name="MS:1001491" value="{qvalue}"/>\n'
        f'{indent}    <UserParam type="float" name="MS:1001493" value="{pep}"/>\n'
        f"{indent}  </PeptideHit>\n"
        f'{indent}  <UserParam type="int" name="map_index" value="0"/>\n'
        f"{indent}</{tag}>\n"
    )


def _percolator_xml(peptide_level):
    """PEPTIDEK: best PSM (scan 1) unassigned, a second PSM (scan 2) on the feature. ELVISLIVEK: only a
    non-representative PSM on its feature, so peptide-level mode has nothing but the 1.0 placeholder."""
    second = 1.0 if peptide_level else 0.004
    elements = _element(
        0, [(0, 100.0, 1000.0)], [_xcorr_pid("PeptideIdentification", 2, "PEPTIDEK", second, second, indent="    ")]
    ) + _element(
        1,
        [(0, 300.0, 2000.0)],
        [_xcorr_pid("PeptideIdentification", 3, "ELVISLIVEK", second, second, rt=300.0, indent="    ")],
        mz=460.0,
        rt=300.0,
    )
    return (
        _HEADER
        + _percolator_run("PI_0", "run_01.mzML", peptide_level)
        + _xcorr_pid("UnassignedPeptideIdentification", 1, "PEPTIDEK", 0.001, 0.01)
        + _maps(["run_01.mzML"])
        + f"  <consensusElementList>\n{elements}  </consensusElementList>\n"
        + "</consensusXML>\n"
    )


def _score(row, name):
    return next((s["score_value"] for s in row["additional_scores"] or () if s["score_name"] == name), None)


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_peptide_level_percolator_confidence(tmp_path, streaming):
    """Peptide-level q-value/PEP come from the peptide's best PSM; the 1.0 placeholders are never reported."""
    source = tmp_path / "percolator.consensusXML"
    source.write_text(_percolator_xml(peptide_level=True))

    out, written = _convert(tmp_path, str(source), streaming)

    features = {row["peptidoform"]: row for row in _rows(written, "feature")}
    assert features["PEPTIDEK"]["peptide_qvalue"] == pytest.approx(0.001)
    assert features["PEPTIDEK"]["posterior_error_probability"] == pytest.approx(0.01)
    assert features["ELVISLIVEK"]["peptide_qvalue"] is None
    assert features["ELVISLIVEK"]["posterior_error_probability"] is None
    psms = {(row["peptidoform"], row["scan"][0]): row for row in _rows(written, "psm")}
    assert all(row["posterior_error_probability"] is None for row in psms.values())
    assert _score(psms[("PEPTIDEK", 1)], "peptide_qvalue") == pytest.approx(0.001)
    assert _score(psms[("PEPTIDEK", 2)], "peptide_qvalue") == pytest.approx(0.001)
    assert _score(psms[("ELVISLIVEK", 3)], "peptide_qvalue") is None
    assert all(_score(row, "q-value") is None for row in psms.values())
    _assert_valid(out, ("feature", "psm"))


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_psm_level_percolator_confidence(tmp_path, streaming):
    """PSM-level q-value/PEP stay per PSM and are not reported as a peptide-level q-value."""
    source = tmp_path / "percolator.consensusXML"
    source.write_text(_percolator_xml(peptide_level=False))

    _, written = _convert(tmp_path, str(source), streaming)

    features = {row["peptidoform"]: row for row in _rows(written, "feature")}
    assert features["PEPTIDEK"]["peptide_qvalue"] is None
    assert features["PEPTIDEK"]["posterior_error_probability"] == pytest.approx(0.004)
    psms = {(row["peptidoform"], row["scan"][0]): row for row in _rows(written, "psm")}
    assert psms[("PEPTIDEK", 1)]["posterior_error_probability"] == pytest.approx(0.01)
    assert _score(psms[("PEPTIDEK", 1)], "q-value") == pytest.approx(0.001)
    assert _score(psms[("PEPTIDEK", 2)], "q-value") == pytest.approx(0.004)


def _shared_peak_xml():
    """VLFHNLPSL's run_01 peak sits in two consensus features: one spanning both runs, one
    holding only that peak (different consensus RT). The PSM is attached to the singleton."""
    spanning = _element(0, [(0, 100.0, 5000.0), (1, 102.0, 4000.0)], [], rt=101.0)
    spanning = spanning.replace(
        "  </consensusElement>",
        _pid("PeptideIdentification", "PI_0", 9, ["VLFHNLPSL"], map_index=1).replace("\n  ", "\n    ") + "  </consensusElement>",
    )
    singleton = _element(
        1,
        [(0, 100.0, 5000.0)],
        [_pid("PeptideIdentification", "PI_0", 7, ["VLFHNLPSL"], map_index=0).replace("\n  ", "\n    ")],
        quality=3.0,
        rt=100.0,
    )
    return (
        _HEADER
        + _identification_run("PI_0", "")
        + _maps(["run_01.mzML", "run_02.mzML"])
        + f"  <consensusElementList>\n{spanning}{singleton}  </consensusElementList>\n"
        + "</consensusXML>\n"
    )


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_shared_peak_is_counted_once(tmp_path, streaming, caplog):
    """One peak in one run reported by two consensus features keeps the row of the one with more runs."""
    source = tmp_path / "shared_peak.consensusXML"
    source.write_text(_shared_peak_xml())

    out, written = _convert(tmp_path, str(source), streaming)

    features = _rows(written, "feature")
    assert sorted(row["run_file_name"] for row in features) == ["run_01", "run_02"]
    assert {row["consensus_rt"] for row in features} == {101.0}
    run_01 = next(row["feature_id"] for row in features if row["run_file_name"] == "run_01")
    psm = next(row for row in _rows(written, "psm") if row["scan"] == [7])
    assert psm["feature_id"] == run_01
    assert "Dropped 1 duplicate feature row(s)" in caplog.text
    _assert_valid(out, ("feature", "psm"))


def _two_candidate_xml(near_first):
    """Two consensus features of one peptide in run_01, at RT 150 and 101, each holding a copy of
    the scan=10 PSM (precursor RT 100). The link must go to the nearer one whatever the order."""
    copy = _pid("PeptideIdentification", "PI_0", 10, ["PEPTIDEK"], map_index=0).replace("\n  ", "\n    ")
    far = _element(0, [(0, 150.0, 1000.0)], [copy], quality=5.0, rt=150.0)
    near = _element(1, [(0, 101.0, 900.0)], [copy], quality=1.0, rt=101.0)
    elements = near + far if near_first else far + near
    return (
        _HEADER
        + _identification_run("PI_0", "")
        + _maps(["run_01.mzML"])
        + f"  <consensusElementList>\n{elements}  </consensusElementList>\n"
        + "</consensusXML>\n"
    )


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_psm_feature_link_does_not_depend_on_feature_order(tmp_path, streaming):
    links = {}
    for near_first in (True, False):
        source = tmp_path / f"order_{near_first}.consensusXML"
        source.write_text(_two_candidate_xml(near_first))
        _, written = _convert(tmp_path, str(source), streaming, name=f"order_{near_first}")
        features = {row["feature_id"]: row["rt"] for row in _rows(written, "feature")}
        (psm,) = _rows(written, "psm")
        links[near_first] = features[psm["feature_id"]]
    assert links == {True: 101.0, False: 101.0}
