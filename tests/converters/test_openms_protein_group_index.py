"""ProteinGroupIndex.resolve: per-run precedence and the merged-index fallback."""

import pytest

from qpx.converters.openms_consensus.protein_groups import ProteinGroupIndex


@pytest.fixture
def index():
    merged = ProteinGroupIndex.from_groups([["A", "B"], ["C"]])
    merged.by_identification = {
        "PI_0": ProteinGroupIndex.from_groups([["A", "B"]]),
        "PI_EMPTY": ProteinGroupIndex(),
    }
    return merged


def test_known_run_resolves_from_its_own_groups(index):
    assert index.resolve({"A"}, "PI_0") == ("A", "B")


def test_unknown_run_falls_back_to_the_merged_index(index):
    """An identification naming a run the index does not know used to null the
    protein attribution of every feature from that run."""
    assert index.resolve({"C"}, "PI_UNKNOWN") == ("C",)


def test_run_without_protein_inference_falls_back_to_the_merged_index(index):
    assert index.resolve({"C"}, "PI_EMPTY") == ("C",)


def test_known_run_that_excluded_the_protein_stays_authoritative(index):
    """PI_0 inferred groups and did not include C; the merged index must not override it."""
    assert index.resolve({"C"}, "PI_0") is None


def test_no_identifier_uses_the_merged_index(index):
    assert index.resolve({"C"}, "") == ("C",)


def test_empty_evidence_resolves_to_nothing(index):
    assert index.resolve(set(), "PI_0") is None
    assert index.resolve(set(), "PI_UNKNOWN") is None
