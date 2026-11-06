"""Tests for PublicOntology registry and OBO parser."""

import tempfile

import pytest

from qpx.core.obo import (
    _normalize_name,
    extract_obo_version,
    parse_obo,
    read_parquet_ontology_metadata,
    terms_to_arrow,
    write_terms_parquet,
)
from qpx.core.ontology import PublicOntology

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Minimal OBO snippet for offline testing
_MINI_OBO = """\
format-version: 1.2
data-version: 4.1.235
ontology: ms

[Term]
id: MS:1001143
name: PSM-level search engine specific statistic
def: "A score or other statistic." [PSI:MS]

[Term]
id: MS:1001492
name: percolator:score
def: "The Percolator score." [PSI:MS]
is_a: MS:1001143 ! PSM-level search engine specific statistic
relationship: has_order MS:1002108 ! higher score is better
synonym: "percolator score" EXACT []

[Term]
id: MS:1001493
name: percolator:PEP
def: "The Percolator posterior error probability." [PSI:MS]
is_a: MS:1001143 ! PSM-level search engine specific statistic
relationship: has_order MS:1002109 ! lower score is better

[Term]
id: MS:1002523
name: Q Exactive HF
def: "Thermo Scientific Q Exactive HF." [PSI:MS]

[Term]
id: MS:9999999
name: obsolete term
def: "This term is obsolete." [PSI:MS]
is_obsolete: true
"""


@pytest.fixture
def mini_obo_terms():
    """Parse the mini OBO snippet."""
    return parse_obo(_MINI_OBO, source="MS")


@pytest.fixture
def mini_parquet(mini_obo_terms, tmp_path):
    """Write mini OBO terms to a temp Parquet file."""
    path = tmp_path / "test_ms.parquet"
    write_terms_parquet(mini_obo_terms, path, source="test_ms", version="4.1.235")
    return path


@pytest.fixture
def mini_ontology(mini_parquet):
    """Create a PublicOntology from the mini Parquet file."""
    return PublicOntology.from_parquet(mini_parquet)


# ---------------------------------------------------------------------------
# Tests (6 total)
# ---------------------------------------------------------------------------


def test_normalize_name():
    """Name normalization handles colons, dots, mixed case, special chars."""
    cases = [
        ("percolator:score", "percolator_score"),
        ("Comet.xcorr", "comet_xcorr"),
        ("MSFragger:Hyperscore", "msfragger_hyperscore"),
        ("andromeda_score", "andromeda_score"),
        ("x!tandem:expect", "x_tandem_expect"),
    ]
    for input_name, expected in cases:
        assert _normalize_name(input_name) == expected


def test_parse_obo(mini_obo_terms):
    """Comprehensive OBO parsing: count, fields, scores, parents, synonyms, obsolete."""
    assert len(mini_obo_terms) == 5
    terms_by_acc = {t.accession: t for t in mini_obo_terms}

    # Basic fields
    assert "MS:1001492" in terms_by_acc
    perc = terms_by_acc["MS:1001492"]
    assert perc.name == "percolator:score"
    assert perc.normalized_name == "percolator_score"
    assert "Percolator score" in perc.definition

    # Score attributes
    assert perc.is_score is True
    assert perc.higher_better is True
    assert terms_by_acc["MS:1001493"].higher_better is False
    assert terms_by_acc["MS:1002523"].is_score is False

    # Hierarchy and synonyms
    assert "MS:1001143" in perc.is_a
    assert "percolator score" in perc.synonyms

    # Obsolete
    assert terms_by_acc["MS:9999999"].is_obsolete is True
    assert perc.is_obsolete is False

    # Version extraction
    assert extract_obo_version(_MINI_OBO) == "4.1.235"
    assert extract_obo_version("remark: version: 1.2.3\n[Term]") == "1.2.3"
    assert extract_obo_version("[Term]\nid: X:1") == "unknown"


def test_parquet_round_trip(mini_parquet, mini_obo_terms):
    """Write/read metadata and Arrow table schema."""
    meta = read_parquet_ontology_metadata(mini_parquet)
    assert meta["qpx_ontology_source"] == "test_ms"
    assert meta["qpx_ontology_version"] == "4.1.235"

    table = terms_to_arrow(mini_obo_terms)
    assert "accession" in table.schema.names
    assert "normalized_name" in table.schema.names
    assert "is_score" in table.schema.names
    assert "is_obsolete" in table.schema.names
    assert len(table) == 5


def test_public_ontology_api(mini_ontology):
    """Core PublicOntology API: version, len, search, lookup, validate, standardize, scores, df, repr."""
    ont = mini_ontology

    # version and length
    assert ont.version == "4.1.235"
    assert len(ont) == 4  # 5 terms, 1 obsolete

    # search
    df = ont.search("percolator")
    assert len(df) >= 2
    assert "percolator:score" in df["name"].values

    # lookup by accession
    term = ont.lookup("MS:1001492")
    assert term is not None
    assert term["name"] == "percolator:score"
    assert ont.lookup("MS:0000000") is None

    # lookup by name (case insensitive)
    term = ont.lookup_name("percolator:score")
    assert term["accession"] == "MS:1001492"
    assert ont.lookup_name("Percolator:Score") is not None

    # lookup by normalized name
    term = ont.lookup_normalized("percolator_score")
    assert term["accession"] == "MS:1001492"

    # validate
    result = ont.validate(["percolator_score", "unknown_score", "percolator_pep"])
    assert result["percolator_score"]
    assert not result["unknown_score"]
    assert result["percolator_pep"]
    assert len(ont.validate([])) == 0

    # standardize
    assert ont.standardize("percolator:score") == "MS:1001492"
    assert ont.standardize("nonexistent_term") is None

    # scores
    scores_df = ont.scores()
    assert len(scores_df) >= 2
    accessions = scores_df["accession"].tolist()
    assert "MS:1001492" in accessions
    assert "MS:1001493" in accessions

    # df
    assert len(ont.df()) == 4
    assert len(ont.df(include_obsolete=True)) == 5

    # repr
    r = repr(ont)
    assert "PublicOntology" in r
    assert "4.1.235" in r


def test_from_obo_and_context_manager(mini_parquet):
    """from_obo creates ontology; context manager works."""
    with tempfile.NamedTemporaryFile(suffix=".obo", mode="w", delete=False) as f:
        f.write(_MINI_OBO)
        f.flush()
        onto = PublicOntology.from_obo(f.name, source="MS")
        assert onto.version == "4.1.235"
        assert len(onto) == 4
        onto.close()

    with PublicOntology.from_parquet(mini_parquet) as onto:
        assert len(onto) == 4


def test_scores_integration():
    """Integration with scores.py: lookup, higher_better, builtin, ontology entries."""
    from qpx.core.scores import is_higher_better, lookup_score, score_ontology_entries

    info = lookup_score("percolator_score")
    assert info is not None
    assert info["ontology_accession"] == "MS:1001492"
    assert is_higher_better("percolator_score") is True

    info = lookup_score("diann_qvalue")
    assert info is not None
    assert info["higher_better"] is False

    entries = score_ontology_entries({"percolator_score", "diann_qvalue"}, view="psm")
    assert len(entries) == 2
    names = {e["field_name"] for e in entries}
    assert "percolator_score" in names
    assert "diann_qvalue" in names
    for e in entries:
        assert "ontology_version" in e


@pytest.mark.parametrize(
    ("score_name", "expected"),
    [
        ("mascot_expectation_value", False),
        ("psm_level_q_value", False),
        ("percolator_q_value", False),
        ("ms_gf_specevalue", False),
        ("x_tandem_expect", False),
        ("comet_expectation_value", False),
        ("posterior_error_probability", False),
        ("mascot_score", True),
        ("ms_gf_rawscore", True),
        ("some_unknown_metric", None),
    ],
)
def test_score_direction_is_never_guessed(score_name, expected):
    """E-values and q-values are lower-is-better; an unknown direction stays None (#154)."""
    from qpx.core.scores import is_higher_better

    assert is_higher_better(score_name) is expected
