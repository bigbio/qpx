"""Protein-sequence calculations shared by converters and transforms.

One definition of each property, so a value computed from a producer's recorded
sequence (OpenMS consensusXML) and one computed from a FASTA agree exactly.
"""

from __future__ import annotations

# U/O have defined masses; J denotes the isobaric I/L pair. B/Z/X do not.
_MASS_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWYJUO")


def average_molecular_weight_kda(sequence: str | None) -> float | None:
    """Theoretical average mass of an unmodified protein sequence, in kDa.

    Includes the terminal water. Null for an empty sequence or one containing a
    residue without a defined mass (B/Z/X) or modification notation.
    """
    if not sequence or not set(sequence).issubset(_MASS_RESIDUES):
        return None
    from pyopenms import AASequence  # pylint: disable=no-name-in-module

    return AASequence.fromString(sequence).getAverageWeight() / 1000.0


def peptide_occurrences(peptide: str | None, protein: str | None) -> list[tuple[int, int]]:
    """Every occurrence of ``peptide`` in ``protein`` as one-based inclusive spans.

    Overlapping occurrences are all reported (``AA`` in ``AAA`` occurs twice). The
    match is exact: I and L are distinct residues here, as in the search database.
    """
    if not peptide or not protein:
        return []
    spans: list[tuple[int, int]] = []
    start = protein.find(peptide)
    while start != -1:
        spans.append((start + 1, start + len(peptide)))
        start = protein.find(peptide, start + 1)
    return spans


def sequence_coverage_percent(protein: str | None, spans) -> float | None:
    """Percent of ``protein`` residues covered by the union of one-based ``spans``.

    Null when the protein is unknown or nothing maps to it: zero coverage from
    an absent peptide list is not a measurement.
    """
    if not protein:
        return None
    covered = bytearray(len(protein))
    for start, end in spans:
        first, last = max(start, 1) - 1, min(end, len(protein))
        if first < last:
            covered[first:last] = b"\x01" * (last - first)
    n_covered = sum(covered)
    if not n_covered:
        return None
    return 100.0 * n_covered / len(protein)
