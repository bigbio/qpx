"""Shared helpers for converter ontology writing."""

from __future__ import annotations

from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from qpx.core.ontology_fields import emitted_field_names
from qpx.core.scores import (
    _FIELD_CV_MAP,
    field_ontology_entries,
    score_ontology_entries,
    unresolved_field_entry,
)

_NON_DATA_VIEWS = frozenset({"ontology", "provenance", "dataset"})


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


def _entry_for_emitted_name(name: str, view: str, columns: set[str]) -> dict | None:
    """Return the row of an emitted name without one; schema columns need one only when a CV term names them."""
    if name in columns:
        if name not in _FIELD_CV_MAP:
            return None
        return field_ontology_entries(view, {name: None}, None)[0]
    entries = score_ontology_entries({name}, view)
    return entries[0] if entries else unresolved_field_entry(name, view)


def reconcile_ontology_entries(entries: list[dict], output_folder: Path, prefix: str) -> list[dict]:
    """Make the ontology rows describe the names the views actually emit.

    A row whose view or name was never written is dropped, and an emitted score,
    CV parameter, intensity name or CV-mapped column without a row gets one, so
    the view stays true however the converter assembled its mappings.
    """
    emitted: dict[str, set[str]] = {}
    with duckdb.connect() as con:
        for path in sorted(Path(output_folder).glob(f"{prefix}.*.parquet")):
            view = path.name[len(prefix) + 1 : -len(".parquet")]
            if view not in _NON_DATA_VIEWS:
                relation = "read_parquet('" + str(path).replace("'", "''") + "')"
                emitted[view] = emitted_field_names(con, relation)

    kept = [entry for entry in entries if entry["field_name"] in emitted.get(entry["view"], ())]
    described = {(entry["field_name"], entry["view"]) for entry in kept}
    for view, names in sorted(emitted.items()):
        columns = set(pq.read_schema(Path(output_folder) / f"{prefix}.{view}.parquet").names)
        for name in sorted(names):
            if (name, view) not in described and (entry := _entry_for_emitted_name(name, view, columns)):
                kept.append(entry)
    return kept
