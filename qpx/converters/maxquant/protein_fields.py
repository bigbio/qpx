"""Parsing of MaxQuant protein cells: gene names and protein group ids."""

from __future__ import annotations

import re
from typing import Optional

import pandas as pd

from qpx.core.protein_annotation import fasta_header_fields
from qpx.writers.protein import ProteinWriter


def _genes_from_headers(fasta_raw: str) -> Optional[list[str]]:
    """Extract distinct gene symbols from FASTA headers."""
    genes = []
    for hdr in str(fasta_raw).split(";"):
        m = re.search(r"GN=([^\s;]+)", hdr)
        if m and m.group(1) not in genes:
            genes.append(m.group(1))
    return genes or None


def parse_genes(gene_raw, fasta_raw) -> Optional[list[str]]:
    """Parse gene names from a proteinGroups row (Gene names, else Fasta headers)."""
    genes: list[str] | None = None
    if gene_raw is not None and pd.notna(gene_raw) and str(gene_raw).strip():
        genes = [g.strip() for g in str(gene_raw).split(";") if g.strip()] or None
    if not genes and fasta_raw is not None and pd.notna(fasta_raw) and str(fasta_raw).strip():
        genes = _genes_from_headers(fasta_raw)
    return genes


def group_ids(value) -> list[int]:
    """Parse an evidence ``Protein group IDs`` cell (``"12"`` or ``"12;40"``)."""
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return []
    text = value if isinstance(value, str) else str(int(value))
    return [int(part) for part in text.split(";") if part.strip()]


def _aligned_cell(cell, count: int) -> list:
    """The tokens of a ``;``-joined cell when there is one per protein, else nothing to align."""
    tokens = str(cell).split(";") if pd.notna(cell) else []
    return tokens if len(tokens) == count else [None] * count


def _protein_rows(table: pd.DataFrame) -> dict[str, dict]:
    rows: dict[str, dict] = {}
    for cells in table.to_dict("records"):
        accessions = str(cells["Protein IDs"]).split(";")
        headers = _aligned_cell(cells.get("Fasta headers"), len(accessions))
        lengths = _aligned_cell(cells.get("Sequence lengths"), len(accessions))
        for accession, header, length in zip(accessions, headers, lengths):
            if accession.startswith("REV__") or accession in rows:
                continue
            fields = fasta_header_fields(header.lstrip(">")) if header else {}
            rows[accession] = {**fields, "accession": accession, "length": int(length) if length else None}
    return rows


def write_protein_view(protein_groups_path: str, output_path: str) -> None:
    """The protein view of what ``proteinGroups.txt`` states: its FASTA headers and sequence lengths."""
    columns = {"Protein IDs", "Fasta headers", "Sequence lengths"}
    table = pd.read_csv(protein_groups_path, sep="\t", usecols=lambda name: name in columns, dtype=str)
    rows = _protein_rows(table)
    if rows:
        with ProteinWriter(output_path) as writer:
            writer.write_batch([rows[accession] for accession in sorted(rows)])
