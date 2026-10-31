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


def to_bare(accession: str) -> str:
    """Convert full UniProt accession to bare form.

    sp|P04114|APOB_HUMAN           → P04114
    CONTAM_sp|CONTAM_P02768|...    → CONTAM_P02768
    REV__sp|P04114|APOB_HUMAN      → REV__P04114
    tr|A0A123|A0A123_HUMAN         → A0A123
    P04114                         → P04114  (no-op)
    """
    if accession is None:
        return accession
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
        parquet_path: Input feature, pg, site or proteoform Parquet file.
        output_path: Where to write the normalized file.
        direction: "forward" (full → bare) or "reverse" (bare → full).
        fasta_lookup: Required for reverse. {bare_accession: full_header}.

    Returns:
        Dict with stats: {"rows", "accessions_changed", "file"}.
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
            _rewrite_view(parquet_path, staging, lambda table: _normalize_table(table, changes, keep_names))
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
    return {"rows": total_rows, "accessions_changed": n_changed, "file": str(output_path)}


def _scan_accessions(parquet: pq.ParquetFile) -> tuple[dict[str, int], set[str]]:
    """Count the first scalar accession column and collect every distinct accession."""
    names = parquet.schema_arrow.names
    scalar_columns = [name for name in _SCALAR_ACCESSION_COLUMNS if name in names]
    columns = [name for name in (*_SCALAR_ACCESSION_COLUMNS, "pg_accessions") if name in names]
    accession_counts: dict[str, int] = {}
    accessions: set[str] = set()
    if not columns:
        return accession_counts, accessions
    for batch in parquet.iter_batches(columns=columns):
        for name in scalar_columns:
            for entry in pc.call_function("value_counts", [batch.column(name)]).to_pylist():
                if entry["values"] is not None:
                    accessions.add(entry["values"])
                    if name == scalar_columns[0]:
                        accession_counts[entry["values"]] = accession_counts.get(entry["values"], 0) + entry["counts"]
        if "pg_accessions" in columns:
            members = _member_accessions(batch.column("pg_accessions"))
            if members is not None:
                accessions.update(pc.call_function("unique", [members]).to_pylist())
    accessions.discard(None)
    return accession_counts, accessions


def _member_accessions(column: pa.Array) -> pa.Array | None:
    """The accession strings inside a ``pg_accessions`` list, or None for another layout."""
    if _is_list_of_strings(column.type):
        return column.values
    if _is_list_of_structs_with_accession(column.type):
        return column.values.field("accession")
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
    """Apply *changes* to scalar accession columns and the ``pg_accessions`` members of one batch."""
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
    if "pg_accessions" in table.schema.names:
        index = table.schema.get_field_index("pg_accessions")
        column = table.column(index)
        if _is_list_of_strings(column.type) or _is_list_of_structs_with_accession(column.type):
            chunks = [_normalize_members(chunk, changes) for chunk in column.chunks]
            table = table.set_column(index, table.schema.field(index), pa.chunked_array(chunks, column.type))
    return table


def _normalize_members(chunk: pa.Array, changes: dict[str, str]) -> pa.Array:
    """Rebuild one ``pg_accessions`` list chunk with its member accessions mapped."""
    values = chunk.values
    if _is_list_of_strings(chunk.type):
        mapped = _mapped(values, changes)
    else:
        children = [
            _mapped(values.field(i), changes) if field.name == "accession" else values.field(i)
            for i, field in enumerate(values.type)
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


def _is_list_of_structs_with_accession(arrow_type: pa.DataType) -> bool:
    """Check if type is list<struct{accession: string, ...}>."""
    if not (pa.types.is_list(arrow_type) or pa.types.is_large_list(arrow_type)):
        return False
    vt = arrow_type.value_type
    if not pa.types.is_struct(vt):
        return False
    return any(f.name == "accession" for f in vt)
