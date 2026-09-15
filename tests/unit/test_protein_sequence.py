"""Shared protein-sequence calculations."""

import pytest

from qpx.core.protein_sequence import average_molecular_weight_kda, peptide_occurrences, sequence_coverage_percent


def test_occurrences_are_one_based_inclusive_and_include_overlaps():
    assert peptide_occurrences("PEP", "MPEPTIDEPEP") == [(2, 4), (9, 11)]
    assert peptide_occurrences("AA", "AAA") == [(1, 2), (2, 3)]


def test_occurrences_are_exact_so_i_and_l_differ():
    assert peptide_occurrences("PEPLIDE", "MPEPIIDE") == []


@pytest.mark.parametrize("peptide, protein", [(None, "MPEP"), ("", "MPEP"), ("PEP", None), ("PEP", "")])
def test_occurrences_empty_for_missing_inputs(peptide, protein):
    assert peptide_occurrences(peptide, protein) == []


def test_coverage_is_the_union_of_spans():
    # residues 2-4 and 3-6 overlap -> 5 of 10 residues covered
    assert sequence_coverage_percent("ABCDEFGHIJ", [(2, 4), (3, 6)]) == pytest.approx(50.0)


def test_coverage_null_when_nothing_maps_or_protein_unknown():
    assert sequence_coverage_percent("ABCDEF", []) is None
    assert sequence_coverage_percent(None, [(1, 2)]) is None


def test_coverage_clips_spans_to_the_protein():
    assert sequence_coverage_percent("ABCD", [(3, 9)]) == pytest.approx(50.0)


def test_molecular_weight_matches_the_openms_consensus_rule():
    # Glycine x 10 average mass = 10 * 57.0513 + 18.0153 = 588.53 Da
    assert average_molecular_weight_kda("GGGGGGGGGG") == pytest.approx(0.58853, abs=1e-4)


@pytest.mark.parametrize("sequence", [None, "", "PEPBTIDE", "PEP(Oxidation)K"])
def test_molecular_weight_null_without_a_defined_mass(sequence):
    assert average_molecular_weight_kda(sequence) is None
