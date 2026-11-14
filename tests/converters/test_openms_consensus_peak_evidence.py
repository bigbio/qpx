"""Regression coverage for confidence scopes and duplicate measured-peak evidence."""

from copy import deepcopy

import pytest
from defusedxml import ElementTree

from qpx.converters.openms_consensus import conversion
from qpx.converters.openms_consensus.conversion import _stream_feature_psm
from qpx.core import files
from qpx.core.data.identity import derive_id
from tests.converters.test_openms_consensus_merged_ids import (
    _HEADER,
    _assert_valid,
    _convert,
    _element,
    _maps,
    _percolator_run,
    _percolator_xml,
    _rows,
    _score,
    _shared_peak_xml,
    _xcorr_pid,
)


def _scoped_confidence_xml(second_value, same_source):
    """The same peptide in two runs, analysed either together or independently."""
    source = "PI_0" if same_source else "PI_1"
    second = _xcorr_pid("PeptideIdentification", 2, "PEPTIDEK", second_value, second_value)
    second = second.replace("PI_0", source).replace('name="map_index" value="0"', 'name="map_index" value="1"')
    header = _HEADER + _percolator_run("PI_0", "run_A.mzML,run_B.mzML" if same_source else "run_A.mzML", True)
    if not same_source:
        header += _percolator_run("PI_1", "run_B.mzML", True)
    elements = _element(0, [(0, 100, 1000)], [_xcorr_pid("PeptideIdentification", 1, "PEPTIDEK", 0.001, 0.01)])
    elements += _element(1, [(1, 100, 2000)], [second])
    return (
        header + _maps(["run_A.mzML", "run_B.mzML"]) + f"<consensusElementList>{elements}</consensusElementList></consensusXML>"
    )


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("same_source", [False, True])
@pytest.mark.parametrize("second_value", [0.2, 1.0])
def test_confidence_stays_with_its_fdr_source(tmp_path, streaming, same_source, second_value):
    """Only a shared analysis may propagate confidence across runs; placeholders stay null otherwise."""
    source = tmp_path / "sources.consensusXML"
    source.write_text(_scoped_confidence_xml(second_value, same_source))
    out, written = _convert(tmp_path, str(source), streaming)
    features = {rec["run_file_name"]: rec for rec in _rows(written, "feature")}
    expected_q = 0.001 if same_source else (second_value if second_value < 1 else None)
    expected_pep = 0.01 if same_source else expected_q
    assert features["run_A"]["peptide_qvalue"] == pytest.approx(0.001)
    assert features["run_B"]["peptide_qvalue"] == pytest.approx(expected_q)
    assert features["run_B"]["posterior_error_probability"] == pytest.approx(expected_pep)
    psm = next(rec for rec in _rows(written, "psm") if rec["run_file_name"] == "run_B")
    assert _score(psm, "peptide_qvalue") == pytest.approx(expected_q)
    _assert_valid(out, ("feature", "psm"))


def _shared_peak_source(path, reverse=False, both_identified=False):
    """Vary duplicate order and whether the quantitative winner already has an own-run PSM."""
    root = ElementTree.fromstring(_shared_peak_xml())
    elements = root.find("consensusElementList")
    if both_identified:
        pid = ElementTree.fromstring(_xcorr_pid("PeptideIdentification", 8, "VLFHNLPSL", 0.02, 0.03))
        elements[0].append(pid)
    if reverse:
        elements[:] = list(reversed(elements))
    path.write_text(ElementTree.tostring(root, encoding="unicode"))


def _flush_every_record(monkeypatch):
    """Exercise updates after the earlier candidate has already reached Parquet."""

    def write_small_batches(*args, **kwargs):
        kwargs["batch"] = 1
        return _stream_feature_psm(*args, **kwargs)

    monkeypatch.setattr(conversion, "_stream_feature_psm", write_small_batches)


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("structures", [("feature", "psm", "pg"), ("feature",), ("pg",), ("psm", "pg")])
def test_shared_peak_preserves_evidence_and_pg_quantity(tmp_path, monkeypatch, streaming, reverse, structures):
    """The 5000 peak contributes once, retains scan 7, and links to its recomputed persisted ID."""
    source = tmp_path / "shared.consensusXML"
    _shared_peak_source(source, reverse)
    _flush_every_record(monkeypatch)
    out, written = _convert(tmp_path, str(source), streaming, structures=structures, include_unassigned_psms=False)
    if "feature" in structures:
        features = _rows(written, "feature")
        assert len(features) == 2
        own = next(rec for rec in features if rec["run_file_name"] == "run_01")
        assert own["consensus_rt"] == 101
        assert own["scan"] == [7]
        assert own["id_run_file_name"] == "run_01"
        assert own["peptide_qvalue"] == pytest.approx(0.9)
        expected = derive_id(
            [own[key] for key in ("peptidoform", "charge", "run_file_name", "rt", "scan", "observed_mz", "consensus_rt")]
        )
        assert own["feature_id"] == expected
        if "psm" in structures:
            psm = next(rec for rec in _rows(written, "psm") if rec["scan"] == [7])
            assert psm["feature_id"] == expected
    if "pg" in structures:
        assert sorted(rec["intensity"] for rec in _rows(written, "pg")) == [4000, 5000]
    if "feature" not in structures and "psm" in structures:
        assert all(rec["feature_id"] is None for rec in _rows(written, "psm"))
    _assert_valid(out, structures)


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
def test_duplicate_direct_ids_keep_scans_and_confidence_together(tmp_path, monkeypatch, streaming, reverse):
    """Combine matching own-run scans while retaining the quantitative winner's confidence pair."""
    source = tmp_path / "direct.consensusXML"
    _shared_peak_source(source, reverse, both_identified=True)
    _flush_every_record(monkeypatch)
    out, written = _convert(tmp_path, str(source), streaming)
    own = next(rec for rec in _rows(written, "feature") if rec["run_file_name"] == "run_01")
    assert own["scan"] == [7, 8]
    assert own["posterior_error_probability"] == pytest.approx(0.03)
    assert own["peptide_qvalue"] is None  # The selected PSM's primary score is xcorr.
    psms = [rec for rec in _rows(written, "psm") if rec["run_file_name"] == "run_01"]
    assert {tuple(rec["scan"]) for rec in psms} == {(7,), (8,)}
    assert {rec["feature_id"] for rec in psms} == {own["feature_id"]}
    _assert_valid(out, ("feature", "psm"))


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("structures", [("feature", "pg"), ("pg",)])
def test_pg_replaces_a_superseded_quantity(tmp_path, streaming, reverse, structures):
    """Identical Feature IDs with different quantities use the higher-quality candidate exactly once."""
    root = ElementTree.fromstring(_shared_peak_xml())
    elements = root.find("consensusElementList")
    earlier = elements[1]
    better = deepcopy(earlier)
    better.set("id", "e_2")
    better.set("quality", "5")
    better.find("groupedElementList/element").set("it", "7000")
    elements[:] = [better, earlier] if reverse else [earlier, better]
    source = tmp_path / "quantities.consensusXML"
    source.write_text(ElementTree.tostring(root, encoding="unicode"))
    out, written = _convert(tmp_path, str(source), streaming, structures=structures)
    assert [rec["intensity"] for rec in _rows(written, "pg")] == [7000]
    if "feature" in structures:
        assert [rec["intensities"][0]["intensity"] for rec in _rows(written, "feature")] == [7000]
    _assert_valid(out, structures)


@pytest.mark.parametrize("peptide_level", [False, True])
@pytest.mark.parametrize("unassigned_last", [False, True])
def test_streaming_uses_one_context_pass_and_one_output_pass(tmp_path, monkeypatch, peptide_level, unassigned_last):
    """Recovered IDs and peptide confidence share a pre-pass, even with trailing unassigned IDs."""
    source = tmp_path / "passes.consensusXML"
    root = ElementTree.fromstring(_percolator_xml(peptide_level))
    if unassigned_last:
        pid = root.find("UnassignedPeptideIdentification")
        root.remove(pid)
        root.append(pid)
    source.write_text(ElementTree.tostring(root, encoding="unicode"))
    parse = files.iterparse
    completed = []

    def count_complete_parses(*args, **kwargs):
        yield from parse(*args, **kwargs)
        completed.append(True)

    monkeypatch.setattr(files, "iterparse", count_complete_parses)
    _, written = _convert(tmp_path, str(source), True)
    assert len(completed) == 2
    feature = next(rec for rec in _rows(written, "feature") if rec["peptidoform"] == "PEPTIDEK")
    assert feature["posterior_error_probability"] == pytest.approx(0.01 if peptide_level else 0.004)
