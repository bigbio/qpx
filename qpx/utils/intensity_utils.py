"""
Utility functions for standardized protein intensity calculations.

This module provides shared functions for calculating total_all_peptides_intensity
and top3_intensity metrics across different proteomics platforms (DIA-NN, MaxQuant).
"""

import math
from typing import List


def calculate_total_all_peptides_intensity(peptide_intensities: List[float]) -> float:
    """
    Calculate the sum of all valid peptide intensities.

    Args:
        peptide_intensities: List of peptide intensity values

    Returns:
        Sum of all valid (non-NaN, non-negative) intensities.
        Returns NaN if the list is empty or contains no valid values.
    """
    if not peptide_intensities:
        return float("nan")

    valid_intensities = [
        intensity
        for intensity in peptide_intensities
        if not math.isnan(intensity) and intensity >= 0
    ]

    return sum(valid_intensities) if valid_intensities else float("nan")


def calculate_top3_intensity(peptide_intensities: List[float]) -> float:
    """
    Calculate the sum of the top 3 highest peptide intensities.

    Note: This function is deprecated for cases where peptidoform-to-peptide
    aggregation is needed. Use calculate_top3_peptide_intensity instead.

    Args:
        peptide_intensities: List of peptide intensity values (already aggregated by peptide)

    Returns:
        Sum of the top 3 highest valid (non-NaN, non-negative) intensities.
        If fewer than 3 valid peptides exist, returns the sum of all valid intensities.
        Returns NaN if the list is empty or contains no valid values.
    """
    if not peptide_intensities:
        return float("nan")

    valid_intensities = [
        intensity
        for intensity in peptide_intensities
        if not math.isnan(intensity) and intensity >= 0
    ]

    if not valid_intensities:
        return float("nan")

    # Sort in descending order and take top 3
    sorted_intensities = sorted(valid_intensities, reverse=True)
    top3 = sorted_intensities[:3]

    return sum(top3)


def calculate_top3_peptide_intensity(
    peptide_sequences: List[str], intensities: List[float]
) -> float:
    """
    Calculate the sum of the top 3 most intense peptides, considering all peptidoforms.

    This function first aggregates intensities by peptide sequence (summing all
    peptidoforms for each peptide), then selects the top 3 peptides by total intensity.

    Definition: "Intensity from the top 3 most intense peptides per protein group
    and RAW file (considering all peptidoforms)."

    Args:
        peptide_sequences: List of peptide sequences (e.g., Stripped.Sequence)
        intensities: List of corresponding intensity values (e.g., Precursor.Quantity)

    Returns:
        Sum of the top 3 highest peptide intensities (after aggregating peptidoforms).
        If fewer than 3 peptides exist, returns the sum of all peptide intensities.
        Returns NaN if the lists are empty or contain no valid values.
    """
    if not peptide_sequences or not intensities:
        return float("nan")

    if len(peptide_sequences) != len(intensities):
        return float("nan")

    # Aggregate intensities by peptide sequence
    peptide_totals: dict = {}
    for seq, intensity in zip(peptide_sequences, intensities):
        if seq is None or (isinstance(intensity, float) and math.isnan(intensity)):
            continue
        if intensity < 0:
            continue

        if seq in peptide_totals:
            peptide_totals[seq] += intensity
        else:
            peptide_totals[seq] = intensity

    if not peptide_totals:
        return float("nan")

    # Sort peptides by total intensity and take top 3
    sorted_peptides = sorted(peptide_totals.values(), reverse=True)
    top3 = sorted_peptides[:3]

    return sum(top3)


def create_standardized_intensity_entries(
    sample_accession: str,
    channel: str,
    total_intensity: float,
    top3_intensity: float,
) -> List[dict]:
    """
    Create standardized intensity entries for additional_intensities array.

    Args:
        sample_accession: Sample identifier
        channel: Channel identifier (e.g., "LFQ")
        total_intensity: Calculated total_all_peptides_intensity value
        top3_intensity: Calculated top3_intensity value

    Returns:
        List of dictionaries with standardized intensity entries,
        each containing sample_accession, channel, and intensities array.
    """
    return [
        {
            "sample_accession": sample_accession,
            "channel": channel,
            "intensities": [
                {
                    "intensity_name": "total_all_peptides_intensity",
                    "intensity_value": total_intensity,
                },
                {
                    "intensity_name": "top3_intensity",
                    "intensity_value": top3_intensity,
                },
            ],
        }
    ]
