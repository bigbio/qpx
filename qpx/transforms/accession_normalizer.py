"""Normalize protein accessions between UniProt full format and bare accession.

Forward:  sp|P04114|APOB_HUMAN  →  P04114
Reverse:  P04114  →  sp|P04114|APOB_HUMAN  (requires FASTA database)

Handles prefixes: sp|, tr|, CONTAM_sp|, ENTRAP_sp|, DECOY_sp|, REV__sp|, rev_sp|
Decoy, contaminant and entrapment markers are preserved. The forward direction
moves entry names into a null ``pg_names`` so they are not lost.
"""

from __future__ import annotations

import logging
import os
import re
import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from qpx.converters.utils import canonical_accession, uniprot_entry_name

logger = logging.getLogger(__name__)

# Matches:  sp|ACCESSION|NAME  or  tr|ACCESSION|NAME
# Also:     CONTAM_sp|CONTAM_ACCESSION|...  or  ENTRAP_sp|ENTRAP_ACCESSION|...
_FULL_RE = re.compile(
    r"^(?P<prefix>(?:CONTAM_|ENTRAP_|DECOY_|REV__|rev_)?(?:sp|tr))\|"
    r"(?P<accession>[^|]+)\|"
    r"(?P<name>.+)$"
)

_DECOY_MARKER = re.compile(r"^(?:DECOY_|REV__|rev_)")

_SCALAR_ACCESSION_COLUMNS = ("anchor_protein", "protein_accession", "accession", "family")
# List column -> the struct field holding the protein (a list<string> column uses its values).
# sequence_edits[].accession is a variant id, not a protein, so it is left alone.
_LIST_ACCESSION_COLUMNS = {"pg_accessions": "accession", "protein_accessions": "accession", "sequence_edits": "protein_accession"}
ACCESSION_VIEWS = ("feature", "pg", "psm", "pepmap", "site", "proteoform")


def to_bare(accession: str) -> str:
    """Convert a full UniProt identifier to its bare accession.

    Only ``sp|``/``tr|`` headers whose accession is a UniProtKB accession are stripped; any other
    identifier, including other ``|``-separated database formats, is returned unchanged.

    sp|P04114|APOB_HUMAN           → P04114
    CONTAM_sp|CONTAM_P02768|...    → CONTAM_P02768
    REV__sp|P04114|APOB_HUMAN      → REV__P04114
    tr|A0A123|A0A123_HUMAN         → A0A123
    DECOY_sp|P04114|APOB_HUMAN     → DECOY_P04114
    P04114                         → P04114  (no-op)
    """
    return canonical_accession(accession)


def to_full(accession: str, lookup: dict[str, str]) -> str:
    """Convert bare accession to full UniProt format using a FASTA-derived lookup.

    P04114 → sp|P04114|APOB_HUMAN  (if in lookup)
    REV__P04114 → REV__sp|P04114|APOB_HUMAN
    P04114 → P04114                (if not in lookup, returned as-is)
    """
    if accession is None:
        return accession
    # Already in full format?
    if _FULL_RE.match(accession):
        return accession
    if accession in lookup:
        return lookup[accession]
    marker = _DECOY_MARKER.match(accession)
    if marker and accession[marker.end() :] in lookup:
        return marker.group() + lookup[accession[marker.end() :]]
    return accession


def parse_fasta_headers(fasta_path: str | Path) -> dict[str, str]:
    """Parse FASTA file and return {bare_accession: full_header} mapping.

    Reads headers like:
        >sp|P04114|APOB_HUMAN Apolipoprotein B-100 ...
        >tr|A0A123|A0A123_HUMAN ...
        >CONTAM_sp|CONTAM_P02768|CONTAM_ALBU_HUMAN ...

    Returns:
        {"P04114": "sp|P04114|APOB_HUMAN",
         "A0A123": "tr|A0A123|A0A123_HUMAN",
         "CONTAM_P02768": "CONTAM_sp|CONTAM_P02768|CONTAM_ALBU_HUMAN"}
    """
    lookup: dict[str, str] = {}
    fasta_path = Path(fasta_path)
    with open(fasta_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            # Extract the identifier (first whitespace-delimited token after >)
            header = line[1:].split()[0]
            m = _FULL_RE.match(header)
            if m:
                lookup[m.group("accession")] = header
    logger.info("Parsed %d protein headers from %s", len(lookup), fasta_path.name)
    return lookup


def normalize_parquet(
    parquet_path: Path,
    output_path: Path,
    direction: str = "forward",
    fasta_lookup: dict[str, str] | None = None,
) -> dict:
    """Normalize accessions in a QPX parquet file.

    The accession columns are scanned first and each distinct accession is
    transformed once. A file with nothing to change is copied (or, in place,
    left untouched); otherwise it is rewritten in bounded batches with the QPX
    Parquet encoding.

    Args:
        parquet_path: Input feature, pg, psm, pepmap, site or proteoform Parquet file.
        output_path: Where to write the normalized file.
        direction: "forward" (full → bare) or "reverse" (bare → full).
        fasta_lookup: Required for reverse. {bare_accession: full_header}.

    Returns:
        Dict with stats: {"rows", "accessions_changed", "distinct_changed", "file"}.
    """
    from qpx.transforms.protein_properties import _rewrite_view

    if direction == "reverse" and not fasta_lookup:
        raise ValueError("FASTA lookup required for reverse normalization")

    transform_fn = to_bare if direction == "forward" else lambda a: to_full(a, fasta_lookup)

    parquet_path, output_path = Path(parquet_path), Path(output_path)
    parquet = pq.ParquetFile(str(parquet_path))
    total_rows = parquet.metadata.num_rows
    accession_counts, accessions = _scan_accessions(parquet)
    changes = {accession: new for accession in accessions if accession and (new := transform_fn(accession)) != accession}
    n_changed = sum(count for accession, count in accession_counts.items() if accession in changes)

    if not changes:
        if output_path.resolve() != parquet_path.resolve():
            shutil.copyfile(parquet_path, output_path)
    else:
        staging = output_path.with_name(f".{output_path.name}.normalizing")
        try:
            keep_names = direction == "forward"
            _rewrite_view(
                parquet_path, staging, lambda table: _normalize_table(table, changes, keep_names), _identity_footer(parquet)
            )
            os.replace(staging, output_path)
        finally:
            staging.unlink(missing_ok=True)

    logger.info(
        "Normalized %s → %s (%d rows, %d accession changes)",
        parquet_path.name,
        output_path.name,
        total_rows,
        n_changed,
    )
    return {"rows": total_rows, "accessions_changed": n_changed, "distinct_changed": len(changes), "file": str(output_path)}


def accession_view_files(dataset_dir: Path, prefix: str) -> list[Path]:
    """The accession-bearing view files: each flat ``<prefix>.<view>.parquet``, else every part of a
    partitioned ``<view>/`` directory."""
    dataset_dir = Path(dataset_dir)
    files: list[Path] = []
    for view in ACCESSION_VIEWS:
        flat = dataset_dir / f"{prefix}.{view}.parquet"
        if flat.is_file():
            files.append(flat)
        elif (dataset_dir / view).is_dir():
            files.extend(sorted((dataset_dir / view).rglob("*.parquet")))
    return files


def normalize_dataset_accessions(dataset_dir: Path, prefix: str) -> int:
    """Strip every view of a dataset to bare accessions in place; returns the files rewritten."""
    return sum(bool(normalize_parquet(f, f)["distinct_changed"]) for f in accession_view_files(dataset_dir, prefix))


def _identity_footer(parquet: pq.ParquetFile) -> dict[bytes, bytes] | None:
    """Footer keys for a rewrite that changes accessions: ids hashed from them no longer re-derive."""
    declared = (parquet.schema_arrow.metadata or {}).get(b"identity_composite", b"").decode()
    if set(declared.split(",")) & {*_SCALAR_ACCESSION_COLUMNS, *_LIST_ACCESSION_COLUMNS}:
        return {b"identity_source": b"provided"}
    return None


def _scan_accessions(parquet: pq.ParquetFile) -> tuple[dict[str, int], set[str]]:
    """Count the first scalar accession column (else the first list column's members) and collect
    every distinct accession."""
    names = parquet.schema_arrow.names
    scalar_columns = [name for name in _SCALAR_ACCESSION_COLUMNS if name in names]
    list_columns = [name for name in _LIST_ACCESSION_COLUMNS if name in names]
    columns = [*scalar_columns, *list_columns]
    accession_counts: dict[str, int] = {}
    accessions: set[str] = set()
    if not columns:
        return accession_counts, accessions
    for batch in parquet.iter_batches(columns=columns):
        for name in scalar_columns:
            _add_value_counts(batch.column(name), accessions, accession_counts if name == scalar_columns[0] else None)
        for name in list_columns:
            members = _member_accessions(batch.column(name), _LIST_ACCESSION_COLUMNS[name])
            if members is None:
                continue
            if not scalar_columns and name == list_columns[0]:
                _add_value_counts(members, accessions, accession_counts)
            else:
                accessions.update(pc.call_function("unique", [members]).to_pylist())
    accessions.discard(None)
    return accession_counts, accessions


def _add_value_counts(values: pa.Array, accessions: set[str], counts: dict[str, int] | None) -> None:
    """Collect the distinct non-null *values*, and add their occurrences to *counts* when given."""
    for entry in pc.call_function("value_counts", [values]).to_pylist():
        if entry["values"] is not None:
            accessions.add(entry["values"])
            if counts is not None:
                counts[entry["values"]] = counts.get(entry["values"], 0) + entry["counts"]


def _member_accessions(column: pa.Array, field: str) -> pa.Array | None:
    """The accession strings inside an accession list column, or None for another layout."""
    if _is_list_of_strings(column.type):
        return column.values
    if _is_list_of_structs_with(column.type, field):
        return column.values.field(field)
    return None


def _names_from_members(table: pa.Table) -> pa.Table:
    """Fill a null ``pg_names`` with the entry names of ``db|ACC|NAME`` members, before those are stripped."""
    index = table.schema.get_field_index("pg_names")
    names = [
        name if name is not None else _entry_names(members)
        for name, members in zip(table.column(index).to_pylist(), table.column("pg_accessions").to_pylist())
    ]
    return table.set_column(index, table.schema.field(index), pa.array(names, type=table.schema.field(index).type))


def _entry_names(members: list[str] | None) -> list[str] | None:
    """Entry names aligned with *members*, or None unless every member has one."""
    names = [uniprot_entry_name(member) for member in members or []]
    return names if names and all(names) else None


def _normalize_table(table: pa.Table, changes: dict[str, str], keep_names: bool = False) -> pa.Table:
    """Apply *changes* to scalar accession columns and accession list members of one batch."""
    if keep_names and {"pg_names", "pg_accessions"} <= set(table.schema.names):
        if _is_list_of_strings(table.schema.field("pg_accessions").type):
            table = _names_from_members(table)
    for name in _SCALAR_ACCESSION_COLUMNS:
        if name in table.schema.names:
            index = table.schema.get_field_index(name)
            column = table.column(index)
            table = table.set_column(
                index, table.schema.field(index), pa.chunked_array([_mapped(c, changes) for c in column.chunks], column.type)
            )
    for name, field in _LIST_ACCESSION_COLUMNS.items():
        if name not in table.schema.names:
            continue
        index = table.schema.get_field_index(name)
        column = table.column(index)
        if _is_list_of_strings(column.type) or _is_list_of_structs_with(column.type, field):
            chunks = [_normalize_members(chunk, changes, field) for chunk in column.chunks]
            table = table.set_column(index, table.schema.field(index), pa.chunked_array(chunks, column.type))
    return table


def _normalize_members(chunk: pa.Array, changes: dict[str, str], field: str) -> pa.Array:
    """Rebuild one accession list chunk with its member accessions mapped."""
    values = chunk.values
    if _is_list_of_strings(chunk.type):
        mapped = _mapped(values, changes)
    else:
        children = [
            _mapped(values.field(i), changes) if child.name == field else values.field(i) for i, child in enumerate(values.type)
        ]
        mapped = pa.StructArray.from_arrays(children, fields=list(values.type), mask=values.is_null())
    return type(chunk).from_arrays(chunk.offsets, mapped, type=chunk.type, mask=chunk.is_null())


def _mapped(values: pa.Array, changes: dict[str, str]) -> pa.Array:
    """Map each distinct string through *changes*, keeping unmapped values and nulls."""
    encoded = pc.call_function("dictionary_encode", [values])
    dictionary = [changes.get(value, value) for value in encoded.dictionary.to_pylist()]
    return pa.array(dictionary, type=values.type).take(encoded.indices)


def _is_list_of_strings(arrow_type: pa.DataType) -> bool:
    """Check if type is list<string> or large_list<string>."""
    if not (pa.types.is_list(arrow_type) or pa.types.is_large_list(arrow_type)):
        return False
    return pa.types.is_string(arrow_type.value_type) or pa.types.is_large_string(arrow_type.value_type)


def _is_list_of_structs_with(arrow_type: pa.DataType, field: str) -> bool:
    """Check if type is list<struct{<field>: string, ...}>."""
    if not (pa.types.is_list(arrow_type) or pa.types.is_large_list(arrow_type)):
        return False
    vt = arrow_type.value_type
    if not pa.types.is_struct(vt):
        return False
    return any(f.name == field for f in vt)
