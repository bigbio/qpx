"""Empty PSM exports warn and leave outputs consistent in both readers."""

import logging

import pyarrow.parquet as pq
import pytest
from defusedxml.ElementTree import fromstring, tostring

from qpx.converters.openms_consensus import converter
from qpx.converters.openms_consensus.converter import _stream_feature_psm
from qpx.dataset import Dataset
from tests.converters.test_openms_consensus import _TMT_CONSENSUSXML

_LOGGER = "qpx.converters.openms_consensus.converter"
_WARNING = "No exportable PSM records"


@pytest.fixture(name="streaming", params=[False, True], ids=["memory", "streaming"])
def streaming_mode(request):
    """Exercise each reader explicitly."""
    return request.param


@pytest.fixture(name="no_psm_xml")
def no_psm_xml_file(tmp_path):
    """Keep feature and protein evidence but remove the spectrum reference."""
    root = fromstring(_TMT_CONSENSUSXML)
    for pid in root.iter("PeptideIdentification"):
        pid.attrib.pop("spectrum_reference")
    path = tmp_path / "no_psm.consensusXML"
    path.write_bytes(tostring(root, encoding="utf-8"))
    return path


@pytest.fixture(name="existing_psm")
def existing_psm_file(tmp_path, psm_parquet):
    """A valid PSM output from an earlier conversion using a custom prefix."""
    out = tmp_path / "out"
    out.mkdir()
    path = out / "rerun.psm.parquet"
    path.write_bytes(psm_parquet.read_bytes())
    return path


@pytest.mark.parametrize("structures", [("psm",), ("feature", "psm", "pg")])
def test_empty_psm_is_skipped(tmp_path, no_psm_xml, streaming, structures, caplog):
    """Missing spectra must not break metadata or create an empty PSM file."""
    out = tmp_path / "out"
    with caplog.at_level(logging.WARNING, logger=_LOGGER):
        written = converter.OpenMSConsensusConverter().convert(
            str(no_psm_xml), str(out), structures=structures, streaming=streaming
        )

    assert set(written) == set(structures) - {"psm"}
    assert not (out / "openms.psm.parquet").exists()
    assert sum(_WARNING in message for message in caplog.messages) == 1
    assert all(pq.read_metadata(path).num_rows > 0 for path in written.values())
    if written:
        provenance = pq.read_table(out / "openms.provenance.parquet").to_pylist()
        assert set(provenance[0]["output_views"]) == set(written)
        assert all("psm" not in step["output_views"] for step in provenance)
        assert (out / "openms.dataset.parquet").exists()
    else:
        assert not list(out.iterdir())


def test_no_psm_request_does_not_warn(existing_psm, no_psm_xml, streaming, caplog):
    """Feature-only exports neither warn about nor delete an unrequested PSM view."""
    previous = existing_psm.read_bytes()
    converter.OpenMSConsensusConverter().convert(
        str(no_psm_xml), str(existing_psm.parent), output_prefix="rerun", structures=("feature",), streaming=streaming
    )
    assert _WARNING not in caplog.text
    assert existing_psm.read_bytes() == previous


@pytest.mark.parametrize("structures", [("psm",), ("feature", "psm")])
def test_empty_psm_removes_existing_file(existing_psm, no_psm_xml, streaming, structures):
    """Dataset discovery must not expose stale matches after an empty PSM rerun."""
    out = existing_psm.parent
    previous = existing_psm.read_bytes()
    other_psm = out / "other.psm.parquet"
    other_psm.write_bytes(previous)
    with Dataset(out, file_prefix="rerun", duckdb_threads=24) as dataset:
        assert dataset.psm is not None

    written = converter.OpenMSConsensusConverter().convert(
        str(no_psm_xml), str(out), output_prefix="rerun", structures=structures, streaming=streaming
    )

    assert set(written) == set(structures) - {"psm"}
    assert not existing_psm.exists()
    assert other_psm.read_bytes() == previous
    with Dataset(out, file_prefix="rerun", duckdb_threads=24) as dataset:
        assert dataset.psm is None
        assert (dataset.feature is not None) == ("feature" in structures)
    if written:
        provenance = pq.read_table(out / "rerun.provenance.parquet").to_pylist()
        assert all("psm" not in step["output_views"] for step in provenance)


def test_failed_conversion_preserves_existing_psm(existing_psm, no_psm_xml, streaming, monkeypatch):
    """A failed conversion must retain the previous PSM output."""
    previous = existing_psm.read_bytes()

    def fail(*_args, **_kwargs):
        raise RuntimeError("conversion failed")

    if streaming:
        monkeypatch.setattr(converter, "_convert_streaming", fail)
    else:
        monkeypatch.setattr(converter.OpenMSConsensusConverter, "_convert_in_memory", fail)
    with pytest.raises(RuntimeError, match="conversion failed"):
        converter.OpenMSConsensusConverter().convert(
            str(no_psm_xml), str(existing_psm.parent), output_prefix="rerun", structures=("psm",), streaming=streaming
        )
    assert existing_psm.read_bytes() == previous


@pytest.mark.parametrize("include_unassigned", [True, False])
def test_unassigned_only_psms(tmp_path, streaming, include_unassigned, caplog):
    """Unlinked PSMs remain valid; explicitly excluding all of them warns."""
    root = fromstring(_TMT_CONSENSUSXML)
    element = root.find(".//consensusElement")
    pid = element.find("PeptideIdentification")
    element.remove(pid)
    pid.tag = "UnassignedPeptideIdentification"
    root.append(pid)
    source = tmp_path / "unassigned.consensusXML"
    source.write_bytes(tostring(root, encoding="utf-8"))
    out = tmp_path / "out"

    written = converter.OpenMSConsensusConverter().convert(
        str(source),
        str(out),
        structures=("psm",),
        streaming=streaming,
        include_unassigned_psms=include_unassigned,
    )

    if include_unassigned:
        rows = pq.read_table(written["psm"]).to_pylist()
        assert len(rows) == 1
        assert rows[0]["feature_id"] is None
        assert _WARNING not in caplog.text
    else:
        assert not written
        assert not list(out.iterdir())
        assert _WARNING in caplog.text


def test_nonempty_psm_is_registered_after_flush(tmp_path, streaming, monkeypatch, caplog):
    """A fully flushed streaming buffer still represents a nonempty PSM view."""

    def flush_every_record(*args, **kwargs):
        kwargs["batch"] = 1
        return _stream_feature_psm(*args, **kwargs)

    monkeypatch.setattr(converter, "_stream_feature_psm", flush_every_record)
    source = tmp_path / "input.consensusXML"
    source.write_text(_TMT_CONSENSUSXML, encoding="utf-8")
    out = tmp_path / "out"

    written = converter.OpenMSConsensusConverter().convert(
        str(source), str(out), structures=("feature", "psm", "pg"), streaming=streaming
    )

    assert set(written) == {"feature", "psm", "pg"}
    assert pq.read_metadata(written["psm"]).num_rows == 1
    provenance = pq.read_table(out / "openms.provenance.parquet").to_pylist()
    assert "psm" in provenance[0]["output_views"]
    assert _WARNING not in caplog.text
