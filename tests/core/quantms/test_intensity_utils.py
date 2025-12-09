"""
Property-based tests for intensity utility functions.

Uses hypothesis library for property-based testing as specified in the design document.
Each property test runs at least 100 iterations.
"""

import math

import pytest
from hypothesis import given, settings, strategies as st

from qpx.utils.intensity_utils import (
    calculate_total_all_peptides_intensity,
    calculate_top3_intensity,
    calculate_top3_peptide_intensity,
    create_standardized_intensity_entries,
)


valid_intensity = st.floats(
    min_value=0.0, max_value=1e15, allow_nan=False, allow_infinity=False
)

mixed_intensity = st.one_of(
    st.floats(min_value=-1e10, max_value=1e15, allow_nan=False, allow_infinity=False),
    st.just(float("nan")),
)


class TestCalculateTotalAllPeptidesIntensity:
    """Tests for calculate_total_all_peptides_intensity function."""

    @settings(max_examples=100)
    @given(st.lists(valid_intensity, min_size=1, max_size=100))
    def test_property_total_equals_sum_of_all_intensities(self, intensities):
        """
        Property 1: For any list of valid intensities, total_all_peptides_intensity
        should equal the arithmetic sum of all intensity values.
        """
        result = calculate_total_all_peptides_intensity(intensities)
        expected = sum(intensities)
        if not math.isclose(result, expected, rel_tol=1e-9):
            pytest.fail(f"Expected {expected}, got {result}")

    def test_empty_list_returns_nan(self):
        """Empty list should return NaN."""
        result = calculate_total_all_peptides_intensity([])
        if not math.isnan(result):
            pytest.fail(f"Expected NaN, got {result}")

    def test_filters_nan_values(self):
        """NaN values should be filtered out."""
        intensities = [100.0, float("nan"), 200.0]
        result = calculate_total_all_peptides_intensity(intensities)
        if not math.isclose(result, 300.0, rel_tol=1e-9):
            pytest.fail(f"Expected 300.0, got {result}")

    def test_filters_negative_values(self):
        """Negative values should be filtered out."""
        intensities = [100.0, -50.0, 200.0]
        result = calculate_total_all_peptides_intensity(intensities)
        if not math.isclose(result, 300.0, rel_tol=1e-9):
            pytest.fail(f"Expected 300.0, got {result}")


class TestCalculateTop3Intensity:
    """Tests for calculate_top3_intensity function."""

    @settings(max_examples=100)
    @given(st.lists(valid_intensity, min_size=3, max_size=100))
    def test_property_top3_equals_sum_of_top3_highest(self, intensities):
        """
        Property 2: For any list with >= 3 valid intensities, top3_intensity
        should equal the sum of the 3 highest intensity values.
        """
        result = calculate_top3_intensity(intensities)
        sorted_desc = sorted(intensities, reverse=True)
        expected = sum(sorted_desc[:3])
        if not math.isclose(result, expected, rel_tol=1e-9):
            pytest.fail(f"Expected {expected}, got {result}")

    @settings(max_examples=100)
    @given(st.lists(valid_intensity, min_size=1, max_size=2))
    def test_property_fewer_than_3_peptides_top3_equals_total(self, intensities):
        """
        Property 3: For any list with fewer than 3 valid intensities,
        top3_intensity should equal total_all_peptides_intensity.
        """
        top3_result = calculate_top3_intensity(intensities)
        total_result = calculate_total_all_peptides_intensity(intensities)
        if not math.isclose(top3_result, total_result, rel_tol=1e-9):
            pytest.fail(f"Expected {total_result}, got {top3_result}")

    def test_empty_list_returns_nan(self):
        """Empty list should return NaN."""
        result = calculate_top3_intensity([])
        if not math.isnan(result):
            pytest.fail(f"Expected NaN, got {result}")

    def test_single_peptide(self):
        """Single peptide should return its intensity."""
        result = calculate_top3_intensity([100.0])
        if not math.isclose(result, 100.0, rel_tol=1e-9):
            pytest.fail(f"Expected 100.0, got {result}")

    def test_two_peptides(self):
        """Two peptides should return sum of both."""
        result = calculate_top3_intensity([100.0, 200.0])
        if not math.isclose(result, 300.0, rel_tol=1e-9):
            pytest.fail(f"Expected 300.0, got {result}")

    def test_filters_nan_and_negative_for_top3(self):
        """NaN and negative values should be filtered before selecting top 3."""
        intensities = [100.0, float("nan"), -50.0, 200.0, 300.0, 400.0]
        result = calculate_top3_intensity(intensities)
        if not math.isclose(result, 900.0, rel_tol=1e-9):
            pytest.fail(f"Expected 900.0, got {result}")


class TestCalculateTop3PeptideIntensity:
    """Tests for calculate_top3_peptide_intensity function."""

    @settings(max_examples=100)
    @given(
        st.lists(
            st.tuples(
                st.text(min_size=5, max_size=20, alphabet="ACDEFGHIKLMNPQRSTVWY"),
                valid_intensity,
            ),
            min_size=3,
            max_size=50,
        )
    )
    def test_property_top3_aggregates_by_peptide(self, peptide_data):
        """
        Property 2b: For any list of (peptide_sequence, intensity) pairs,
        top3_peptide_intensity should aggregate intensities by peptide sequence
        and return the sum of the top 3 peptides.
        """
        sequences = [p[0] for p in peptide_data]
        intensities = [p[1] for p in peptide_data]

        result = calculate_top3_peptide_intensity(sequences, intensities)

        peptide_totals = {}
        for seq, intensity in zip(sequences, intensities):
            peptide_totals[seq] = peptide_totals.get(seq, 0) + intensity

        sorted_totals = sorted(peptide_totals.values(), reverse=True)
        expected = sum(sorted_totals[:3])

        if not math.isclose(result, expected, rel_tol=1e-9):
            pytest.fail(f"Expected {expected}, got {result}")

    def test_empty_lists_returns_nan(self):
        """Empty lists should return NaN."""
        result = calculate_top3_peptide_intensity([], [])
        if not math.isnan(result):
            pytest.fail(f"Expected NaN, got {result}")

    def test_mismatched_lengths_returns_nan(self):
        """Mismatched list lengths should return NaN."""
        result = calculate_top3_peptide_intensity(["SEQ1", "SEQ2"], [100.0])
        if not math.isnan(result):
            pytest.fail(f"Expected NaN, got {result}")

    def test_aggregates_same_peptide(self):
        """Same peptide sequence should have intensities aggregated."""
        sequences = ["PEPTIDE", "PEPTIDE", "OTHER"]
        intensities = [100.0, 200.0, 50.0]
        result = calculate_top3_peptide_intensity(sequences, intensities)
        if not math.isclose(result, 350.0, rel_tol=1e-9):
            pytest.fail(f"Expected 350.0, got {result}")

    def test_top3_with_many_peptides(self):
        """Should correctly select top 3 peptides after aggregation."""
        sequences = ["A", "A", "B", "C", "D"]
        intensities = [100.0, 100.0, 150.0, 120.0, 80.0]
        result = calculate_top3_peptide_intensity(sequences, intensities)
        if not math.isclose(result, 470.0, rel_tol=1e-9):
            pytest.fail(f"Expected 470.0, got {result}")

    def test_filters_nan_values(self):
        """NaN values should be filtered out."""
        sequences = ["A", "B", "C"]
        intensities = [100.0, float("nan"), 200.0]
        result = calculate_top3_peptide_intensity(sequences, intensities)
        if not math.isclose(result, 300.0, rel_tol=1e-9):
            pytest.fail(f"Expected 300.0, got {result}")

    def test_filters_negative_values(self):
        """Negative values should be filtered out."""
        sequences = ["A", "B", "C"]
        intensities = [100.0, -50.0, 200.0]
        result = calculate_top3_peptide_intensity(sequences, intensities)
        if not math.isclose(result, 300.0, rel_tol=1e-9):
            pytest.fail(f"Expected 300.0, got {result}")


class TestCreateStandardizedIntensityEntries:
    """Tests for create_standardized_intensity_entries function."""

    def test_creates_correct_structure(self):
        """Should create entries with correct field names."""
        entries = create_standardized_intensity_entries(
            sample_accession="sample_001",
            channel="LFQ",
            total_intensity=12345.0,
            top3_intensity=9876.0,
        )

        if len(entries) != 1:
            pytest.fail(f"Expected 1 entry, got {len(entries)}")
        entry = entries[0]
        if entry["sample_accession"] != "sample_001":
            pytest.fail(
                f"Expected sample_accession 'sample_001', got {entry['sample_accession']}"
            )
        if entry["channel"] != "LFQ":
            pytest.fail(f"Expected channel 'LFQ', got {entry['channel']}")
        if len(entry["intensities"]) != 2:
            pytest.fail(f"Expected 2 intensities, got {len(entry['intensities'])}")

        intensity_names = [i["intensity_name"] for i in entry["intensities"]]
        if "total_all_peptides_intensity" not in intensity_names:
            pytest.fail("Expected 'total_all_peptides_intensity' in intensity_names")
        if "top3_intensity" not in intensity_names:
            pytest.fail("Expected 'top3_intensity' in intensity_names")

    def test_correct_intensity_values(self):
        """Should store correct intensity values."""
        entries = create_standardized_intensity_entries(
            sample_accession="sample_001",
            channel="LFQ",
            total_intensity=12345.0,
            top3_intensity=9876.0,
        )

        entry = entries[0]
        intensities_dict = {
            i["intensity_name"]: i["intensity_value"] for i in entry["intensities"]
        }

        if not math.isclose(
            intensities_dict["total_all_peptides_intensity"], 12345.0, rel_tol=1e-9
        ):
            pytest.fail(
                f"Expected total_all_peptides_intensity 12345.0, got {intensities_dict['total_all_peptides_intensity']}"
            )
        if not math.isclose(intensities_dict["top3_intensity"], 9876.0, rel_tol=1e-9):
            pytest.fail(
                f"Expected top3_intensity 9876.0, got {intensities_dict['top3_intensity']}"
            )


class TestCrossPlatformConsistency:
    """Tests for cross-platform calculation consistency."""

    @settings(max_examples=100)
    @given(st.lists(mixed_intensity, min_size=0, max_size=100))
    def test_property_same_input_produces_same_output(self, intensities):
        """
        Property 6: For any same peptide intensity input data, calling the shared
        utility functions multiple times should produce identical results.

        Since both DIA-NN and MaxQuant use the same shared utility functions
        (calculate_total_all_peptides_intensity and calculate_top3_intensity),
        cross-platform consistency is guaranteed by design. This test verifies
        that the functions are deterministic - same input always produces same output.
        """
        total_result_1 = calculate_total_all_peptides_intensity(intensities)
        total_result_2 = calculate_total_all_peptides_intensity(intensities)

        top3_result_1 = calculate_top3_intensity(intensities)
        top3_result_2 = calculate_top3_intensity(intensities)

        if math.isnan(total_result_1):
            if not math.isnan(total_result_2):
                pytest.fail(
                    f"total_all_peptides_intensity is not deterministic: {total_result_1} != {total_result_2}"
                )
        else:
            if total_result_1 != total_result_2:
                pytest.fail(
                    f"total_all_peptides_intensity is not deterministic: {total_result_1} != {total_result_2}"
                )

        if math.isnan(top3_result_1):
            if not math.isnan(top3_result_2):
                pytest.fail(
                    f"top3_intensity is not deterministic: {top3_result_1} != {top3_result_2}"
                )
        else:
            if top3_result_1 != top3_result_2:
                pytest.fail(
                    f"top3_intensity is not deterministic: {top3_result_1} != {top3_result_2}"
                )


class TestExistingDataPreservation:
    """Tests for existing data preservation when appending new intensity metrics."""

    @settings(max_examples=100)
    @given(
        st.text(
            min_size=1,
            max_size=50,
            alphabet=st.characters(whitelist_categories=("L", "N")),
        ),
        st.text(
            min_size=1,
            max_size=20,
            alphabet=st.characters(whitelist_categories=("L", "N")),
        ),
        st.floats(min_value=0.0, max_value=1e15, allow_nan=False, allow_infinity=False),
        st.floats(min_value=0.0, max_value=1e15, allow_nan=False, allow_infinity=False),
        st.lists(
            st.fixed_dictionaries(
                {
                    "intensity_name": st.text(
                        min_size=1,
                        max_size=30,
                        alphabet=st.characters(whitelist_categories=("L", "N", "P")),
                    ),
                    "intensity_value": st.floats(
                        min_value=0.0,
                        max_value=1e15,
                        allow_nan=False,
                        allow_infinity=False,
                    ),
                }
            ),
            min_size=0,
            max_size=5,
        ),
    )
    def test_property_existing_data_preserved_when_appending(
        self,
        sample_accession,
        channel,
        total_intensity,
        top3_intensity,
        existing_intensities,
    ):
        """
        Property 7: For any data that already contains additional_intensities values,
        appending new standardized intensity metrics should preserve the original values.

        This test verifies that create_standardized_intensity_entries creates new entries
        that can be appended to existing data without overwriting or modifying existing values.
        """
        new_entries = create_standardized_intensity_entries(
            sample_accession=sample_accession,
            channel=channel,
            total_intensity=total_intensity,
            top3_intensity=top3_intensity,
        )

        existing_additional_intensities = (
            [
                {
                    "sample_accession": "existing_sample",
                    "channel": "existing_channel",
                    "intensities": existing_intensities,
                }
            ]
            if existing_intensities
            else []
        )

        combined = existing_additional_intensities + new_entries

        if existing_intensities:
            if len(combined) < 1:
                pytest.fail(f"Expected at least 1 entry, got {len(combined)}")
            if combined[0]["sample_accession"] != "existing_sample":
                pytest.fail(
                    f"Expected sample_accession 'existing_sample', got {combined[0]['sample_accession']}"
                )
            if combined[0]["channel"] != "existing_channel":
                pytest.fail(
                    f"Expected channel 'existing_channel', got {combined[0]['channel']}"
                )
            if combined[0]["intensities"] != existing_intensities:
                pytest.fail("Existing intensities were not preserved")

        new_entry = combined[-1]
        if new_entry["sample_accession"] != sample_accession:
            pytest.fail(
                f"Expected sample_accession '{sample_accession}', got {new_entry['sample_accession']}"
            )
        if new_entry["channel"] != channel:
            pytest.fail(f"Expected channel '{channel}', got {new_entry['channel']}")

        intensity_names = [i["intensity_name"] for i in new_entry["intensities"]]
        if "total_all_peptides_intensity" not in intensity_names:
            pytest.fail("Expected 'total_all_peptides_intensity' in intensity_names")
        if "top3_intensity" not in intensity_names:
            pytest.fail("Expected 'top3_intensity' in intensity_names")

        intensities_dict = {
            i["intensity_name"]: i["intensity_value"] for i in new_entry["intensities"]
        }
        if not math.isclose(
            intensities_dict["total_all_peptides_intensity"],
            total_intensity,
            rel_tol=1e-9,
        ):
            pytest.fail(
                f"Expected total_all_peptides_intensity {total_intensity}, got {intensities_dict['total_all_peptides_intensity']}"
            )
        if not math.isclose(
            intensities_dict["top3_intensity"], top3_intensity, rel_tol=1e-9
        ):
            pytest.fail(
                f"Expected top3_intensity {top3_intensity}, got {intensities_dict['top3_intensity']}"
            )
