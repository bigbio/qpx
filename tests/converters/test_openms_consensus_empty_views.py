"""Empty quantification views must not register missing or stale output files."""

import pyarrow.parquet as pq
import pytest
from defusedxml.ElementTree import fromstring, tostring

from qpx.converters.openms_consensus import converter
from qpx.dataset import Dataset
from tests.converters.test_openms_consensus import _TMT_CONSENSUSXML


@pytest.fixture(name="streaming", params=[False, True])
def streaming_mode(request):
    """Exercise both consensusXML readers."""
    return request.param


@pytest.fixture(name="empty_view", params=["feature", "pg"])
def empty_view_name(request):
    """The requested quantification view with no exportable records."""
    return request.param


@pytest.fixture(name="empty_view_xml")
def empty_view_xml_file(tmp_path, empty_view):
    """Retain valid PSMs and the other quantification view."""
    root = fromstring(_TMT_CONSENSUSXML)
    if empty_view == "feature":
        consensus = root.find(".//consensusElement")
        pid = consensus.find("PeptideIdentification")
        consensus.remove(pid)
        pid.tag = "UnassignedPeptideIdentification"
        root.append(pid)
    else:
        proteins = root.find(".//ProteinIdentification")
        proteins.remove(proteins.find("ProteinHit"))
        root.find(".//PeptideHit").attrib.pop("protein_refs")
    source = tmp_path / "empty_view.consensusXML"
    source.write_bytes(tostring(root, encoding="utf-8"))
    return source


@pytest.fixture(name="existing_view")
def existing_view_file(tmp_path, empty_view, request):
    """A valid quantification file from an earlier export."""
    original = request.getfixturevalue(f"{empty_view}_parquet")
    out = tmp_path / "out"
    out.mkdir()
    path = out / f"rerun.{empty_view}.parquet"
    path.write_bytes(original.read_bytes())
    return path


@pytest.mark.parametrize("include_others", [False, True])
def test_empty_quantification_view_is_skipped(tmp_path, empty_view_xml, empty_view, streaming, include_others, caplog):
    """Return and describe only emitted views, including a completely empty export."""
    structures = ("feature", "psm", "pg") if include_others else (empty_view,)
    out = tmp_path / "out"
    written = converter.OpenMSConsensusConverter().convert(
        str(empty_view_xml), str(out), structures=structures, streaming=streaming
    )

    assert set(written) == set(structures) - {empty_view}
    assert not (out / f"openms.{empty_view}.parquet").exists()
    assert sum(f"No exportable {empty_view.upper()} records" in message for message in caplog.messages) == 1
    assert all(pq.read_metadata(path).num_rows > 0 for path in written.values())
    if written:
        provenance = pq.read_table(out / "openms.provenance.parquet").to_pylist()
        assert set(provenance[0]["output_views"]) == set(written)
        assert all(empty_view not in step["output_views"] for step in provenance)
        assert (out / "openms.dataset.parquet").exists()
    else:
        assert not list(out.iterdir())


def test_empty_quantification_rerun_removes_stale_view(existing_view, empty_view_xml, empty_view, streaming):
    """Dataset must not rediscover the old view or lose another output prefix."""
    out = existing_view.parent
    previous = existing_view.read_bytes()
    other_view = out / f"other.{empty_view}.parquet"
    other_view.write_bytes(previous)
    with Dataset(out, file_prefix="rerun", duckdb_threads=24) as dataset:
        assert getattr(dataset, empty_view) is not None

    written = converter.OpenMSConsensusConverter().convert(
        str(empty_view_xml), str(out), output_prefix="rerun", structures=("feature", "psm", "pg"), streaming=streaming
    )

    assert set(written) == {"feature", "psm", "pg"} - {empty_view}
    assert not existing_view.exists()
    assert other_view.read_bytes() == previous
    with Dataset(out, file_prefix="rerun", duckdb_threads=24) as dataset:
        assert getattr(dataset, empty_view) is None
        assert dataset.psm is not None
    provenance = pq.read_table(out / "rerun.provenance.parquet").to_pylist()
    assert all(empty_view not in step["output_views"] for step in provenance)


def test_unrequested_quantification_view_is_preserved(existing_view, empty_view_xml, empty_view, streaming, caplog):
    """A PSM-only export preserves quantification files outside its requested views."""
    previous = existing_view.read_bytes()
    converter.OpenMSConsensusConverter().convert(
        str(empty_view_xml), str(existing_view.parent), output_prefix="rerun", structures=("psm",), streaming=streaming
    )

    assert existing_view.read_bytes() == previous
    assert f"No exportable {empty_view.upper()} records" not in caplog.text


@pytest.mark.parametrize("keep_feature", [False, True])
def test_empty_rerun_removes_only_orphaned_metadata(tmp_path, streaming, ontology_parquet, keep_feature):
    """Clear stale metadata only when no same-prefix data remains after a rerun."""
    source = tmp_path / "input.consensusXML"
    source.write_text(_TMT_CONSENSUSXML, encoding="utf-8")
    out = tmp_path / "out"
    conv = converter.OpenMSConsensusConverter()
    conv.convert(str(source), str(out), output_prefix="rerun", structures=("feature", "psm", "pg"), streaming=streaming)
    (out / "rerun.ontology.parquet").write_bytes(ontology_parquet.read_bytes())
    before = {path.name: path.read_bytes() for path in out.iterdir()}
    other = {name.replace("rerun.", "other.", 1): data for name, data in before.items()}
    for name, data in other.items():
        (out / name).write_bytes(data)
    root = fromstring(_TMT_CONSENSUSXML)
    for parent in root.iter():
        for child in list(parent):
            if child.tag in ("PeptideIdentification", "ProteinHit"):
                parent.remove(child)
    source.write_bytes(tostring(root, encoding="utf-8"))
    structures = ("psm", "pg") if keep_feature else ("feature", "psm", "pg")

    written = conv.convert(str(source), str(out), output_prefix="rerun", structures=structures, streaming=streaming)

    assert not written
    expected = {name: data for name, data in before.items() if name not in ("rerun.psm.parquet", "rerun.pg.parquet")}
    assert {path.name: path.read_bytes() for path in out.glob("rerun.*")} == (expected if keep_feature else {})
    assert {path.name: path.read_bytes() for path in out.glob("other.*")} == other
    with Dataset(out, file_prefix="rerun", duckdb_threads=24) as dataset:
        assert bool(dataset.available_structures) == keep_feature


def test_failed_conversion_preserves_quantification_view(existing_view, empty_view_xml, empty_view, streaming, monkeypatch):
    """Clear an old empty view only after core conversion has completed successfully."""
    previous = existing_view.read_bytes()

    def fail(*_args, **_kwargs):
        raise RuntimeError("conversion failed")

    if streaming:
        monkeypatch.setattr(converter, "_convert_streaming", fail)
    else:
        monkeypatch.setattr(converter.OpenMSConsensusConverter, "_convert_in_memory", fail)
    with pytest.raises(RuntimeError, match="conversion failed"):
        converter.OpenMSConsensusConverter().convert(
            str(empty_view_xml), str(existing_view.parent), output_prefix="rerun", structures=(empty_view,), streaming=streaming
        )
    assert existing_view.read_bytes() == previous
