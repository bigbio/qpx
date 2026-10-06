"""Native spectrum components survive conversion and distinguish PSM identities."""

from copy import deepcopy

import pyarrow.parquet as pq
import pytest
from defusedxml.ElementTree import fromstring, tostring

from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter
from qpx.core.scan import scan_from_native_id
from tests.converters.test_openms_consensus import _TMT_CONSENSUSXML


@pytest.mark.parametrize(
    ("native_id", "expected"),
    [
        ("scan=42", [42]),
        ("index=0", [0]),
        ("spectrum=42", [42]),
        ("controllerType=0 controllerNumber=1 scan=42", [42]),
        ("controllerType=5 controllerNumber=1 scan=7", [5, 1, 7]),
        ("frame=120 scan=475", [120, 475]),
        ("frame=120 scan=475 precursor=3", [120, 475, 3]),
        ("frame=120 windowGroup=2 scan=475", [120, 2, 475]),
        ("merged=0 frame=120 scanStart=4 scanEnd=8", [0, 120, 4, 8]),
        ("function=10 process=1 scan=345", [10, 1, 345]),
        ("sample=1 period=1 cycle=2740 experiment=10", [1, 1, 2740, 10]),
        ("SCAN=42", [42]),
        ("uuid=opaque", []),
        ("", []),
    ],
)
def test_scan_encoding_preserves_native_components(native_id, expected):
    """The documented vendor encodings retain order and repeated components."""
    assert scan_from_native_id(native_id) == expected


def _convert_references(tmp_path, references, streaming):
    """Convert one feature with the requested supporting spectrum references."""
    root = fromstring(_TMT_CONSENSUSXML)
    feature = root.find(".//consensusElement")
    template = feature.find("PeptideIdentification")
    feature.remove(template)
    for reference in references:
        identification = deepcopy(template)
        identification.set("spectrum_reference", reference)
        feature.append(identification)
    source = tmp_path / "native.consensusXML"
    source.write_bytes(tostring(root, encoding="utf-8", xml_declaration=True))
    return OpenMSConsensusConverter().convert(
        str(source), str(tmp_path / "out"), structures=("feature", "psm"), streaming=streaming
    )


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize(
    "references",
    [
        ("frame=120 scan=475", "frame=121 scan=475"),
        ("function=1 process=1 scan=345", "function=2 process=1 scan=345"),
        ("sample=1 period=1 cycle=123 experiment=1", "sample=1 period=1 cycle=123 experiment=2"),
        ("controllerType=5 controllerNumber=1 scan=7", "controllerType=5 controllerNumber=2 scan=7"),
    ],
)
def test_multicomponent_psms_are_not_collapsed(tmp_path, streaming, references):
    """Different native spectra with the same scan ordinal retain separate PSMs."""
    written = _convert_references(tmp_path, (*references, references[0]), streaming)
    psms = pq.read_table(written["psm"]).to_pylist()
    features = pq.read_table(written["feature"]).to_pylist()
    expected = [scan_from_native_id(reference) for reference in references]

    assert len(psms) == 2
    assert len({row["psm_id"] for row in psms}) == 2
    assert [row["scan"] for row in psms] == expected
    assert len(features) == 1
    assert features[0]["scan"] == expected[0] + expected[1]
    assert {row["feature_id"] for row in psms} == {features[0]["feature_id"]}


@pytest.mark.parametrize("streaming", [False, True])
def test_feature_and_psm_keep_repeated_native_components(tmp_path, streaming):
    """A single Sciex identification has identical full Feature and PSM scans."""
    written = _convert_references(tmp_path, ["sample=1 period=1 cycle=123 experiment=2"], streaming)
    for view in ("feature", "psm"):
        assert pq.read_table(written[view]).column("scan").to_pylist() == [[1, 1, 123, 2]]
