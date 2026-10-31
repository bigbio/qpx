"""Parsing of MaxQuant protein cells: gene names and protein group ids."""

from __future__ import annotations

import re
from typing import Optional

import pandas as pd


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
