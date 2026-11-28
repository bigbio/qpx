"""Shared helpers for converter ontology writing."""

from __future__ import annotations

from pathlib import Path

import pyarrow.parquet as pq

from qpx.core.builtin_scores import BUILTIN_SCORES
from qpx.core.engine import DuckDBEngine
from qpx.core.ontology_fields import emitted_field_names
from qpx.core.scores import (
    field_cv_term,
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
        if field_cv_term(view, name) is None:
            return None
        return field_ontology_entries(view, {name: None}, None)[0]
    entries = score_ontology_entries({name}, view)
    return entries[0] if entries else unresolved_field_entry(name, view)


def _current_entry(entry: dict, columns: set[str]) -> dict | None:
    """The row a new conversion writes for the name and source of *entry*, None for no row.

    *entry* itself when QPX defines no term for the name: a score or CV parameter
    whose term came from the source file or the PSI-MS OBO.
    """
    name, view, source = entry["field_name"], entry["view"], entry["source_column_name"]
    if name in columns:
        if source is None:
            return _entry_for_emitted_name(name, view, columns)
        return field_ontology_entries(view, {name: source}, entry["source_tool"])[0]
    if name not in BUILTIN_SCORES:
        return entry
    return {**score_ontology_entries({name}, view)[0], "source_column_name": source, "source_tool": entry["source_tool"]}


def refresh_ontology_entries(entries: list[dict], output_folder: Path, prefix: str) -> list[dict]:
    """Give the rows of an older file the terms QPX defines now, keeping their source column and tool."""
    columns: dict[str, set[str]] = {}
    refreshed = []
    for entry in entries:
        view = entry["view"]
        if view not in columns:
            path = Path(output_folder) / f"{prefix}.{view}.parquet"
            columns[view] = set(pq.read_schema(path).names) if path.is_file() else set()
        current = _current_entry(entry, columns[view])
        if current is not None:
            refreshed.append(current)
    return refreshed


def write_view_ontology(dataset: Path, prefix: str, views: dict[str, Path], target: Path) -> None:
    """Write to ``target`` the ontology view of ``dataset`` with rows describing ``views``.

    For views written after the ontology view, such as an mz view or the protein view
    of a FASTA; ``views`` maps each view name to its Parquet file.
    """
    from qpx.writers.ontology import OntologyWriter

    entries = describe_views(pq.read_table(Path(dataset) / f"{prefix}.ontology.parquet").to_pylist(), views)
    with OntologyWriter(target, creator="qpx") as writer:
        writer.write_batch(entries)


def describe_views(entries: list[dict], views: dict[str, Path]) -> list[dict]:
    """Make the rows of ``views`` (view name to Parquet file) describe the names those files emit.

    A row of one of these views whose name its file does not emit is dropped, and an
    emitted score, CV parameter, intensity name or CV-mapped column without a row gets
    one. The rows of other views are kept as they are.
    """
    emitted: dict[str, set[str]] = {}
    with DuckDBEngine() as engine:
        for view, path in views.items():
            engine.register_parquet(view, path)
            emitted[view] = emitted_field_names(engine, view)

    kept = [entry for entry in entries if entry["view"] not in emitted or entry["field_name"] in emitted[entry["view"]]]
    described = {(entry["field_name"], entry["view"]) for entry in kept}
    for view, names in sorted(emitted.items()):
        columns = set(pq.read_schema(views[view]).names)
        for name in sorted(names):
            if (name, view) not in described and (entry := _entry_for_emitted_name(name, view, columns)):
                kept.append(entry)
    return kept


def reconcile_ontology_entries(entries: list[dict], output_folder: Path, prefix: str) -> list[dict]:
    """Make the ontology rows describe the names the views actually emit.

    A row whose view or name was never written is dropped, and an emitted score,
    CV parameter, intensity name or CV-mapped column without a row gets one, so
    the view stays true however the converter assembled its mappings.
    """
    views = {}
    for path in sorted(Path(output_folder).glob(f"{prefix}.*.parquet")):
        view = path.name[len(prefix) + 1 : -len(".parquet")]
        if view not in _NON_DATA_VIEWS:
            views[view] = path
    return describe_views([entry for entry in entries if entry["view"] in views], views)
