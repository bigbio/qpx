"""Ontology versions describe the actual mapped resource, not unrelated terms."""

import pytest

from qpx.converters.sdrf import SdrfConverter
from qpx.core import scores
from qpx.core.ontology import PublicOntology


@pytest.fixture(name="loaded_ms_version", params=["4.1.235", None])
def _loaded_ms_version(request, tmp_path, monkeypatch):
    """Load a local ontology with an explicit version or no version metadata."""
    header = f"data-version: {request.param}\n" if request.param else ""
    path = tmp_path / "psi-ms.obo"
    path.write_text(
        header
        + '\n[Term]\nid: MS:1001492\nname: percolator:score\ndef: "Percolator score." []\n'
        + '\n[Term]\nid: MS:1000016\nname: scan start time\ndef: "Scan start time." []\n'
        + '\n[Term]\nid: MS:1001251\nname: Trypsin\ndef: "Trypsin cleavage." []\n'
    )
    with PublicOntology.from_obo(path) as ontology:
        monkeypatch.setattr(scores, "_ontology", ontology)
        monkeypatch.setattr(scores, "_lookup_cache", {})
        monkeypatch.setattr("qpx.core.ontology.PublicOntology", lambda *_args, **_kwargs: ontology)
        yield request.param


def test_field_versions_belong_only_to_mapped_ontology(loaded_ms_version):
    """Unmapped fields cannot inherit the PSI-MS registry version."""
    entries = scores.field_ontology_entries(
        view="feature",
        resolved_mappings={"rt": "RT", "lfq": "Precursor.Normalised", "diann_cscore": "CScore"},
        tool_name="DIA-NN",
    )
    by_field = {entry["field_name"]: entry for entry in entries}

    assert by_field["rt"]["ontology_version"] == loaded_ms_version
    for field in ("lfq", "diann_cscore"):
        assert by_field[field]["ontology_accession"] is None
        assert by_field[field]["ontology_source"] is None
        assert by_field[field]["ontology_version"] is None
    assert all(entry["ontology_version"] == loaded_ms_version for entry in scores.field_ontology_entries())


def test_score_versions_belong_only_to_mapped_ontology(loaded_ms_version):
    """Tool-specific scores keep a null version while PSI-MS uses the loaded file."""
    entries = scores.score_ontology_entries({"percolator_score", "diann_cscore", "unmapped_custom_score"})
    by_field = {entry["field_name"]: entry for entry in entries}

    assert by_field["percolator_score"]["ontology_source"] == "MS"
    assert by_field["percolator_score"]["ontology_version"] == loaded_ms_version
    assert by_field["diann_cscore"]["ontology_source"] is None
    assert by_field["diann_cscore"]["ontology_version"] is None
    assert "unmapped_custom_score" not in by_field


def test_run_versions_use_the_loaded_resource(loaded_ms_version, tmp_path):
    """SDRF run ontology mappings use the same known-or-null version convention."""
    sdrf = tmp_path / "input.sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\tcomment[cleavage agent details]\n"
        "sample1\trun1.raw\tlabel free sample\tNT=Trypsin;AC=MS:1001251\n"
    )
    with SdrfConverter(duckdb_threads=6) as converter:
        converter.convert(str(sdrf), run_output=str(tmp_path / "run.parquet"))
        entries = converter.run_ontology_entries()

    assert len(entries) == 1
    assert entries[0]["ontology_accession"] == "MS:1001251"
    assert entries[0]["ontology_version"] == loaded_ms_version


def test_unavailable_ontology_does_not_invent_a_version(monkeypatch):
    """Known mappings remain usable when their resource cannot be loaded."""

    def unavailable():
        raise FileNotFoundError("PSI-MS resource unavailable")

    monkeypatch.setattr(scores, "_get_ontology", unavailable)
    monkeypatch.setattr(scores, "_lookup_cache", {})
    entries = scores.field_ontology_entries(view="feature", resolved_mappings={"rt": "RT", "lfq": "Precursor.Normalised"})
    entries.extend(scores.score_ontology_entries({"qvalue", "diann_cscore"}))

    assert any(entry["ontology_accession"] == "MS:1002354" for entry in entries)
    assert all(entry["ontology_version"] is None for entry in entries)
