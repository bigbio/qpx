"""The protein view: what a sequence database says of the proteins a dataset names."""

from __future__ import annotations

import re
from pathlib import Path

import pyarrow.parquet as pq

_UNIPROT_FIELD = re.compile(r"\b(OS|OX|GN|PE|SV)=")


def fasta_header_fields(header: str) -> dict:
    """Entry name, description, gene, organism and taxonomy id of a FASTA header.

    ``sp|P02768|ALBU_HUMAN Serum albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2``.
    A header of another form yields its description, and nothing it does not state.
    """
    identifier, _, rest = header.strip().partition(" ")
    parts = identifier.split("|")
    fields = {"entry_name": parts[2] if len(parts) >= 3 and parts[2] else None}
    found = _UNIPROT_FIELD.search(rest)
    fields["description"] = (rest[: found.start()] if found else rest).strip() or None
    pieces = _UNIPROT_FIELD.split(rest)
    tagged = dict(zip(pieces[1::2], (piece.strip() for piece in pieces[2::2])))
    taxonomy = tagged.get("OX", "")
    return {
        **fields,
        "organism": tagged.get("OS"),
        "gene_name": tagged.get("GN"),
        "taxonomy_id": int(taxonomy) if taxonomy.isdigit() else None,
    }


_ACCESSION_QUERIES = {
    "pg": "SELECT unnest(pg_accessions) AS accession FROM read_parquet('$source')",
    "feature": "SELECT unnest(pg_accessions).accession AS accession FROM read_parquet('$source')",
}


def _dataset_accessions(dataset_dir: Path, prefix: str) -> list[str]:
    """Every protein accession of the dataset's pg and feature views."""
    from qpx.core.engine import create_connection
    from qpx.core.sql import escape_path, sql_build

    queries = []
    for view, query in _ACCESSION_QUERIES.items():
        path = dataset_dir / f"{prefix}.{view}.parquet"
        if path.is_file():
            queries.append(sql_build(query, source=escape_path(str(path))))
    if not queries:
        return []
    con = create_connection()
    try:
        union = " UNION ".join(queries)
        rows = con.sql(f"SELECT accession FROM ({union}) WHERE accession IS NOT NULL ORDER BY accession").fetchall()
    finally:
        con.close()
    return [accession for (accession,) in rows]


def write_protein_view(dataset_dir: Path, prefix: str, fasta, staging: Path) -> bool:
    """Write the protein view of the dataset's accessions the FASTA has; values a producer recorded stay."""
    from qpx.writers.protein import ProteinWriter

    existing_path = dataset_dir / f"{prefix}.protein.parquet"
    existing = {row["accession"]: row for row in pq.read_table(existing_path).to_pylist()} if existing_path.is_file() else {}
    rows = dict(existing)
    for accession in _dataset_accessions(dataset_dir, prefix):
        entry = fasta.entry(accession)
        if entry is None:
            continue
        header, sequence = entry
        filled = {**fasta_header_fields(header), "accession": accession, "length": len(sequence), "sequence": sequence}
        known = rows.get(accession, {})
        rows[accession] = {name: known[name] if known.get(name) is not None else value for name, value in filled.items()}
    if not rows:
        return False
    with ProteinWriter(staging / f"{prefix}.protein.parquet") as writer:
        writer.write_batch([rows[accession] for accession in sorted(rows)])
    return True
