"""ConsensusXML PSM format metadata reflects native IDs across the full input."""

from copy import deepcopy

import pyarrow.parquet as pq
import pytest
from defusedxml.ElementTree import fromstring, tostring

from qpx.converters.openms_consensus import converter
from qpx.writers.psm import PsmWriter
from tests.converters.test_openms_consensus import _TMT_CONSENSUSXML
from tests.converters.test_openms_consensus_native_id import _convert_references


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize(
    ("references", "expected"),
    [
        (["scan=42", "controllerType=0 controllerNumber=1 scan=43"], b"scan"),
        (["index=0", "index=1"], b"index"),
        (["frame=120 scan=475", "frame=121 scan=475"], b"nativeId"),
        (["scan=42", "index=0"], None),
        (["scan=42", "uuid=opaque"], None),
    ],
)
def test_public_conversion_declares_only_known_uniform_psm_format(tmp_path, streaming, references, expected):
    """Both readers declare one confirmed format, without labelling Feature arrays."""
    written = _convert_references(tmp_path, references, streaming)
    metadata = pq.read_schema(written["psm"]).metadata
    assert metadata.get(b"scan_format") == expected
    assert b"scan_format" not in pq.read_schema(written["feature"]).metadata


def test_streaming_declares_format_after_flushing_batches(tmp_path, monkeypatch):
    """A single streaming pass may confirm the format after row groups were written."""
    from qpx.converters.openms_consensus.streaming import StreamingConsensusMap

    class SmallBatchPsmWriter(PsmWriter):
        """Force multiple real Parquet batches from a small fixture."""

        def __init__(self, *args, **kwargs):
            super().__init__(*args, batch_size=1, **kwargs)

    passes = []
    original = StreamingConsensusMap.iter_all

    def tracked_pass(stream):
        passes.append(1)
        yield from original(stream)

    monkeypatch.setattr(converter, "PsmWriter", SmallBatchPsmWriter)
    monkeypatch.setattr(StreamingConsensusMap, "iter_all", tracked_pass)
    written = _convert_references(tmp_path, ["index=0", "index=1"], streaming=True)
    parquet = pq.ParquetFile(written["psm"])
    assert len(passes) == 1
    assert parquet.metadata.num_row_groups == 2
    expected = b"index" if hasattr(pq.ParquetWriter, "add_key_value_metadata") else None
    assert parquet.metadata.metadata.get(b"scan_format") == expected
    assert parquet.schema_arrow.metadata.get(b"scan_format") == expected


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("include_unassigned", [False, True])
def test_unassigned_identifications_participate_when_included(tmp_path, streaming, include_unassigned):
    """An included unassigned index prevents a scan-only file declaration."""
    source = tmp_path / "unassigned.consensusXML"
    root = fromstring(_TMT_CONSENSUSXML)
    unassigned = deepcopy(root.find(".//PeptideIdentification"))
    unassigned.tag = "UnassignedPeptideIdentification"
    unassigned.set("spectrum_reference", "index=99")
    root.append(unassigned)
    source.write_bytes(tostring(root, encoding="utf-8", xml_declaration=True))
    written = converter.OpenMSConsensusConverter().convert(
        str(source),
        str(tmp_path / "out"),
        structures=("psm",),
        streaming=streaming,
        include_unassigned_psms=include_unassigned,
    )
    assert pq.read_schema(written["psm"]).metadata.get(b"scan_format") == (None if include_unassigned else b"scan")


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize(
    ("second", "expected"),
    [(["scan=44"], b"scan"), (["index=0"], None), (["uuid=opaque"], None)],
)
def test_multiple_inputs_declare_one_format_for_all_files(tmp_path, streaming, second, expected):
    """A later input file neither replaces nor hides an earlier file's format."""
    paths = []
    for run, references in (("run_01", ["scan=42", "scan=43"]), ("run_02", second)):
        root = fromstring(_TMT_CONSENSUSXML.replace("run_01", run))
        feature = root.find(".//consensusElement")
        template = feature.find("PeptideIdentification")
        feature.remove(template)
        for reference in references:
            identification = deepcopy(template)
            identification.set("spectrum_reference", reference)
            feature.append(identification)
        path = tmp_path / f"{run}.consensusXML"
        path.write_bytes(tostring(root, encoding="utf-8", xml_declaration=True))
        paths.append(str(path))
    written = converter.OpenMSConsensusConverter().convert(paths, str(tmp_path / "out"), structures=("psm",), streaming=streaming)
    assert pq.read_schema(written["psm"]).metadata.get(b"scan_format") == expected
