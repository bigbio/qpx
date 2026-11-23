"""Tests for field_ontology_entries() with source provenance."""


def test_field_ontology_entries():
    """field_ontology_entries: with source provenance, missing CV, backward compat."""
    from qpx.core.scores import field_ontology_entries

    # With source provenance
    resolved = {"intensity": "Precursor.Quantity", "rt": "RT"}
    entries = field_ontology_entries(
        view="feature",
        resolved_mappings=resolved,
        tool_name="DIA-NN",
    )
    rt_entries = [e for e in entries if e["field_name"] == "rt"]
    assert len(rt_entries) == 1
    assert rt_entries[0]["source_column_name"] == "RT"
    assert rt_entries[0]["source_tool"] == "DIA-NN"
    assert rt_entries[0]["ontology_accession"] == "MS:1000894"
    assert rt_entries[0]["view"] == "feature"

    # The rt of a PSM is the start time of its MS2 scan.
    entries = field_ontology_entries(view="psm", resolved_mappings={"rt": "Retention time"}, tool_name="MaxQuant")
    assert entries[0]["ontology_accession"] == "MS:1000016"

    # A shared field name is resolved according to its view: PG global q-value
    # must not inherit precursor-level Global.Q.Value semantics.
    entries = field_ontology_entries(
        view="pg",
        resolved_mappings={"global_qvalue": "Global.PG.Q.Value"},
        tool_name="DIA-NN",
    )
    assert entries[0]["ontology_accession"] == "MS:1002373"
    assert entries[0]["source_column_name"] == "Global.PG.Q.Value"

    # Missing CV still written
    resolved = {"lfq": "Precursor.Normalised"}
    entries = field_ontology_entries(
        view="feature",
        resolved_mappings=resolved,
        tool_name="DIA-NN",
    )
    lfq_entries = [e for e in entries if e["field_name"] == "lfq"]
    assert len(lfq_entries) == 1
    assert lfq_entries[0]["source_column_name"] == "Precursor.Normalised"
    assert lfq_entries[0]["ontology_accession"] is None

    # Backward compat (no resolved_mappings)
    entries = field_ontology_entries(view="psm")
    field_names = {e["field_name"] for e in entries}
    assert "posterior_error_probability" in field_names
    assert "rt" in field_names
    for e in entries:
        assert "source_column_name" in e
        assert "source_tool" in e


def test_field_terms_follow_the_view():
    """A column takes the term of its view, never a PSI-MS name match: is_decoy
    is a decoy peptide only for peptides, and the mz id is no PEFF "ID"."""
    from qpx.core.scores import field_ontology_entries

    entries = field_ontology_entries(
        view="pg", resolved_mappings={"is_decoy": "Reverse", "intensity": "Intensity"}, tool_name="MaxQuant"
    )
    assert [entry["ontology_accession"] for entry in entries] == [None, None]
    assert field_ontology_entries(view="mz", resolved_mappings={"id": "id"})[0]["ontology_accession"] is None
    terms = {
        ("xic", "rt"): "MS:1000595",
        ("xic", "intensity"): "MS:1000515",
        ("fragment", "observed_mz"): "MS:1001225",
        ("psm", "observed_mz"): "MS:1000744",
        ("protein", "sequence"): None,
        ("mz", "ms_level"): "MS:1000511",
        ("mz", "scan_start_time"): "MS:1000016",
        ("psm", "sequence"): "MS:1000888",
    }
    for (view, field), accession in terms.items():
        entry = field_ontology_entries(view=view, resolved_mappings={field: field})[0]
        assert entry["ontology_accession"] == accession, (view, field)


def test_dedupe_ontology_entries_merges_complementary_values():
    """Duplicate score and field entries retain CV metadata and provenance."""
    from qpx.converters.ontology_util import dedupe_ontology_entries

    entries = [
        {
            "field_name": "qvalue",
            "view": "feature",
            "ontology_name": "q-value",
            "ontology_accession": "MS:1002354",
            "source_column_name": None,
            "source_tool": None,
        },
        {
            "field_name": "qvalue",
            "view": "feature",
            "ontology_name": "DIA-NN:Q.Value",
            "ontology_accession": "MS:1002354",
            "source_column_name": "Q.Value",
            "source_tool": "DIA-NN",
        },
        {"field_name": "qvalue", "view": "pg", "source_column_name": "PG.Q.Value"},
        {"field_name": "rt", "view": "feature", "source_column_name": "RT"},
    ]
    deduped = dedupe_ontology_entries(entries)

    keys = [(e["field_name"], e["view"]) for e in deduped]
    assert keys == [("qvalue", "feature"), ("qvalue", "pg"), ("rt", "feature")]
    qv_feature = next(e for e in deduped if e["field_name"] == "qvalue" and e["view"] == "feature")
    assert qv_feature["ontology_accession"] == "MS:1002354"
    assert qv_feature["ontology_name"] == "q-value"
    assert qv_feature["source_column_name"] == "Q.Value"
    assert qv_feature["source_tool"] == "DIA-NN"


def test_dedupe_ontology_entries_rejects_conflicts():
    """A primary key cannot silently select between conflicting source columns."""
    import pytest

    from qpx.converters.ontology_util import dedupe_ontology_entries

    entries = [
        {"field_name": "qvalue", "view": "feature", "source_column_name": "Q.Value"},
        {"field_name": "qvalue", "view": "feature", "source_column_name": "Lib.Q.Value"},
    ]
    with pytest.raises(ValueError, match="Conflicting ontology entries"):
        dedupe_ontology_entries(entries)

    entries = [
        {"field_name": "tool_score", "view": "feature", "ontology_name": "Tool score A"},
        {"field_name": "tool_score", "view": "feature", "ontology_name": "Tool score B"},
    ]
    with pytest.raises(ValueError, match="Conflicting ontology entries"):
        dedupe_ontology_entries(entries)


def test_dedupe_ontology_entries_requires_primary_key_fields():
    """Malformed entries cannot be silently collapsed under a null key."""
    import pytest

    from qpx.converters.ontology_util import dedupe_ontology_entries

    with pytest.raises(ValueError, match="non-null field_name and view"):
        dedupe_ontology_entries([{"field_name": "qvalue"}])
