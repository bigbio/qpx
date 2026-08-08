"""Shared helpers for converter ontology writing."""

from __future__ import annotations


def dedupe_ontology_entries(entries: list[dict]) -> list[dict]:
    """Merge complementary entries sharing an ontology primary key.

    Score discovery supplies CV metadata while converter mappings supply source
    provenance. Conflicting identifiers or provenance are rejected.
    """
    conflict_fields = {
        "ontology_accession",
        "ontology_source",
        "ontology_version",
        "source_column_name",
        "source_tool",
    }
    by_key: dict[tuple[str, str], dict] = {}
    deduped: list[dict] = []
    for entry in entries:
        field_name = entry.get("field_name")
        view = entry.get("view")
        if field_name is None or view is None:
            raise ValueError("Ontology entries must contain non-null field_name and view")
        key = (field_name, view)
        current = by_key.get(key)
        if current is not None:
            for field, value in entry.items():
                if value is None:
                    continue
                existing = current.get(field)
                if existing is None:
                    current[field] = value
                elif (
                    existing != value
                    and field == "ontology_name"
                    and current.get("ontology_accession") is None
                    and entry.get("ontology_accession") is None
                ):
                    raise ValueError(f"Conflicting ontology entries for {key}: {field}={existing!r} versus {value!r}")
                elif existing != value and field in conflict_fields:
                    raise ValueError(f"Conflicting ontology entries for {key}: {field}={existing!r} versus {value!r}")
            continue
        merged = dict(entry)
        by_key[key] = merged
        deduped.append(merged)
    return deduped
