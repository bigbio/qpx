"""Unit tests for converter utility methods added in recent commits.

Tests cover:
- BaseConverter._escape_path (SQL path injection prevention)
- MaxQuantPgAdapter._extract_protein_names / _extract_gene_names
- MaxQuantFeatureAdapter._detect_pg_columns / _build_pg_lookup / _lookup_pg_metadata
"""

import pandas as pd
import pytest

from qpx.converters.base import BaseConverter
from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter
from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter
from qpx.converters.utils import parse_uniprot_id, strip_uniprot_prefix
from qpx.core.files import run_file_stem


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("run.part.1.raw", "run.part.1"),
        ("file:///data/run.part.1.mzML.gz", "run.part.1"),
        ("file:///data/run.part.1.mzXML", "run.part.1"),
        (r"C:\data\run.part.1.wiff", "run.part.1"),
        ("run.part.1.dia", "run.part.1"),
        ("run.part.1", "run.part.1"),
    ],
)
def test_run_file_stem_preserves_meaningful_dots(value, expected):
    assert run_file_stem(value) == expected


# ---------------------------------------------------------------------------
# BaseConverter._escape_path
# ---------------------------------------------------------------------------


class TestEscapePath:
    """Validate SQL path escaping."""

    def test_normal_path(self):
        result = BaseConverter._escape_path("/data/results.txt")
        if result != "/data/results.txt":
            raise AssertionError(f"Expected '/data/results.txt', got {result!r}")

    def test_path_with_single_quote(self):
        result = BaseConverter._escape_path("O'Brien/data.txt")
        if result != "O''Brien/data.txt":
            raise AssertionError(f"Expected \"O''Brien/data.txt\", got {result!r}")

    def test_path_with_multiple_quotes(self):
        result = BaseConverter._escape_path("a'b'c.txt")
        if result != "a''b''c.txt":
            raise AssertionError(f"Expected \"a''b''c.txt\", got {result!r}")

    def test_windows_path(self):
        result = BaseConverter._escape_path("C:\\Users\\test\\file.txt")
        if result != "C:\\Users\\test\\file.txt":
            raise AssertionError(f"Expected windows path unchanged, got {result!r}")

    def test_empty_path(self):
        result = BaseConverter._escape_path("")
        if result != "":
            raise AssertionError(f"Expected empty string, got {result!r}")


# ---------------------------------------------------------------------------
# MaxQuantPgAdapter._extract_gene_names
# ---------------------------------------------------------------------------


class TestExtractGeneNames:
    """Validate GN= tag extraction from FASTA headers."""

    def test_single_gene(self):
        result = MaxQuantPgAdapter._extract_gene_names("sp|P12345|PROT_HUMAN Protein X OS=Homo sapiens GN=TP53 PE=1 SV=1")
        if result != ["TP53"]:
            raise AssertionError(f"Unexpected result: {result!r}")

    def test_multiple_headers(self):
        result = MaxQuantPgAdapter._extract_gene_names("sp|P12345|X GN=TP53 PE=1;sp|Q99999|Y GN=BRCA1 PE=1")
        if result != ["TP53", "BRCA1"]:
            raise AssertionError(f"Unexpected result: {result!r}")

    def test_no_gene_name(self):
        result = MaxQuantPgAdapter._extract_gene_names("sp|P12345|PROT_HUMAN Protein X OS=Homo sapiens PE=1 SV=1")
        if result is not None:
            raise AssertionError(f"Expected None, got {result!r}")

    def test_duplicate_gene_names(self):
        result = MaxQuantPgAdapter._extract_gene_names("sp|P12345|X GN=TP53 PE=1;sp|Q99999|Y GN=TP53 PE=1")
        if result != ["TP53"]:
            raise AssertionError(f"Unexpected result: {result!r}")

    def test_empty_string(self):
        result = MaxQuantPgAdapter._extract_gene_names("")
        if result is not None:
            raise AssertionError(f"Expected None, got {result!r}")


# ---------------------------------------------------------------------------
# MaxQuantFeatureAdapter._detect_pg_columns
# ---------------------------------------------------------------------------


class TestDetectPgColumns:
    """Validate column detection in proteinGroups DataFrame."""

    def test_standard_columns(self):
        df = pd.DataFrame({"Protein IDs": ["P1"], "Q-value": [0.01], "Gene names": ["TP53"]})
        acc, qval, gene = MaxQuantFeatureAdapter._detect_pg_columns(df)
        if acc != "Protein IDs":
            raise AssertionError(f"Expected 'Protein IDs', got {acc!r}")
        if qval != "Q-value":
            raise AssertionError(f"Expected 'Q-value', got {qval!r}")
        if gene != "Gene names":
            raise AssertionError(f"Expected 'Gene names', got {gene!r}")

    def test_alternative_columns(self):
        df = pd.DataFrame({"Majority protein IDs": ["P1"], "q-value": [0.01], "Gene Names": ["TP53"]})
        acc, qval, gene = MaxQuantFeatureAdapter._detect_pg_columns(df)
        if acc != "Majority protein IDs":
            raise AssertionError(f"Expected 'Majority protein IDs', got {acc!r}")
        if qval != "q-value":
            raise AssertionError(f"Expected 'q-value', got {qval!r}")
        if gene != "Gene Names":
            raise AssertionError(f"Expected 'Gene Names', got {gene!r}")

    def test_missing_columns(self):
        df = pd.DataFrame({"some_col": [1]})
        acc, qval, gene = MaxQuantFeatureAdapter._detect_pg_columns(df)
        if acc is not None:
            raise AssertionError(f"Expected None for acc, got {acc!r}")
        if qval is not None:
            raise AssertionError(f"Expected None for qval, got {qval!r}")
        if gene is not None:
            raise AssertionError(f"Expected None for gene, got {gene!r}")


# ---------------------------------------------------------------------------
# MaxQuantFeatureAdapter._build_pg_lookup / _lookup_pg_metadata
# ---------------------------------------------------------------------------


class TestBuildPgLookup:
    """Validate membership-keyed protein-group q-value/gene maps."""

    def _build(self, df, qval_col="Q-value", gene_col="Gene names"):
        with MaxQuantFeatureAdapter() as adapter:
            return adapter._build_pg_lookup(df, "Protein IDs", qval_col, gene_col)

    def test_group_keyed_qvalue(self):
        df = pd.DataFrame({"Protein IDs": ["P1;P2", "P3"], "Q-value": [0.01, 0.05], "Gene names": [None, None]})
        maps = self._build(df)
        by_group = maps["by_group"]
        if by_group[frozenset({"P1", "P2"})]["qvalue"] != 0.01:
            raise AssertionError(f"Unexpected group qvalue: {by_group[frozenset({'P1', 'P2'})]!r}")
        if by_group[frozenset({"P3"})]["qvalue"] != 0.05:
            raise AssertionError(f"Unexpected group qvalue: {by_group[frozenset({'P3'})]!r}")

    def test_ambiguous_accession_excluded_from_leader_fallback(self):
        # P1 appears in two distinct groups → must NOT be in by_accession.
        df = pd.DataFrame(
            {
                "Protein IDs": ["P1;P2", "P1;P3"],
                "Q-value": [0.01, 0.20],
                "Gene names": ["GENEA", "GENEB"],
            }
        )
        maps = self._build(df)
        if "P1" in maps["by_accession"]:
            raise AssertionError("Ambiguous accession P1 must be excluded from the leader fallback")
        # P2 and P3 are unique to one group each.
        if maps["by_accession"]["P2"]["qvalue"] != 0.01:
            raise AssertionError(f"Unexpected P2 qvalue: {maps['by_accession']['P2']!r}")
        if maps["by_accession"]["P3"]["qvalue"] != 0.20:
            raise AssertionError(f"Unexpected P3 qvalue: {maps['by_accession']['P3']!r}")

    def test_no_qval_col(self):
        df = pd.DataFrame({"Protein IDs": ["P1"], "Gene names": [None]})
        maps = self._build(df, qval_col=None)
        if maps["by_group"][frozenset({"P1"})]["qvalue"] is not None:
            raise AssertionError("Expected None qvalue when no q-value column present")

    def test_genes_from_gene_column(self):
        df = pd.DataFrame({"Protein IDs": ["P1;P2"], "Q-value": [0.01], "Gene names": ["TP53;BRCA1"]})
        maps = self._build(df)
        genes = maps["by_group"][frozenset({"P1", "P2"})]["genes"]
        if genes != ["TP53", "BRCA1"]:
            raise AssertionError(f"Unexpected genes: {genes!r}")

    def test_genes_fasta_fallback(self):
        df = pd.DataFrame(
            {
                "Protein IDs": ["P1"],
                "Q-value": [0.01],
                "Gene names": [None],
                "Fasta headers": ["sp|P1|X GN=TP53 PE=1"],
            }
        )
        maps = self._build(df)
        genes = maps["by_group"][frozenset({"P1"})]["genes"]
        if genes != ["TP53"]:
            raise AssertionError(f"Unexpected genes: {genes!r}")


class TestLookupPgMetadata:
    """Validate feature-side lookup keyed on full membership, leader as fallback."""

    def test_full_membership_wins_over_leader(self):
        # Razor accession P1 is shared by two groups with different q-values.
        # A feature whose membership is {P1,P3} must get the {P1,P3} group's
        # q-value/genes, not the first group's, even though P1 is the leader.
        df = pd.DataFrame(
            {
                "Protein IDs": ["P1;P2", "P1;P3"],
                "Q-value": [0.01, 0.20],
                "Gene names": ["GENEA", "GENEB"],
            }
        )
        with MaxQuantFeatureAdapter() as adapter:
            maps = adapter._build_pg_lookup(df, "Protein IDs", "Q-value", "Gene names")
            qvalue, genes = adapter._lookup_pg_metadata(maps, ["P1", "P3"], "P1")
        if qvalue != 0.20:
            raise AssertionError(f"Expected 0.20 for the {{P1,P3}} membership, got {qvalue!r}")
        if genes != ["GENEB"]:
            raise AssertionError(f"Expected ['GENEB'] for the {{P1,P3}} membership, got {genes!r}")

    def test_leader_fallback_only_when_unambiguous(self):
        df = pd.DataFrame(
            {
                "Protein IDs": ["P1;P2", "P3;P4"],
                "Q-value": [0.01, 0.05],
                "Gene names": ["GENEA", "GENEB"],
            }
        )
        with MaxQuantFeatureAdapter() as adapter:
            maps = adapter._build_pg_lookup(df, "Protein IDs", "Q-value", "Gene names")
            # Membership doesn't match any group, but leader P3 is unambiguous.
            qvalue, genes = adapter._lookup_pg_metadata(maps, ["P3"], "P3")
        if qvalue != 0.05 or genes != ["GENEB"]:
            raise AssertionError(f"Expected unambiguous leader fallback, got {(qvalue, genes)!r}")


# ---------------------------------------------------------------------------
# UniProt accession helpers (shared parser)
# ---------------------------------------------------------------------------


class TestParseUniprotId:
    """``db|ACCESSION|NAME`` -> ``(accession, name)`` with graceful fallbacks."""

    @pytest.mark.parametrize(
        ("entry", "expected"),
        [
            ("sp|P12345|PROT_HUMAN", ("P12345", "PROT_HUMAN")),
            ("tr|A0A3B3IS91|X_HUMAN", ("A0A3B3IS91", "X_HUMAN")),
            ("A|B", ("B", "B")),  # 2 fields: accession doubles as name
            ("P12345", ("P12345", "P12345")),  # no pipes: entry is both
            ("A|", ("", "")),  # trailing empty field preserved
            ("db|X|Y|Z", ("X", "Y")),  # >3 fields: first two used
        ],
    )
    def test_parse(self, entry, expected):
        assert parse_uniprot_id(entry) == expected


class TestStripUniprotPrefix:
    """Only ``sp|``/``tr|``-prefixed ids are stripped; everything else is verbatim."""

    @pytest.mark.parametrize(
        ("entry", "expected"),
        [
            ("sp|P55011|S12A2_HUMAN", "P55011"),
            ("tr|A0A3B3IS91|X_HUMAN", "A0A3B3IS91"),
            ("P55011", "P55011"),  # bare accession untouched
            ("CON__P1", "CON__P1"),  # contaminant prefix, no sp|/tr| -> verbatim
            ("REV__sp|P1|N", "REV__sp|P1|N"),  # decoy: not a leading sp|/tr| -> verbatim
            ("xyz|abc", "xyz|abc"),  # pipe but not a UniProt db prefix -> verbatim
            ("", ""),
        ],
    )
    def test_strip(self, entry, expected):
        assert strip_uniprot_prefix(entry) == expected
