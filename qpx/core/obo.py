"""OBO ontology parser and Parquet I/O for CV term files.

Downloads and parses OBO (Open Biomedical Ontology) files, converts
them to Arrow tables, and writes versioned Parquet files.

Usage::

    from qpx.core.obo import parse_obo, download_obo, write_terms_parquet

    text = download_obo(PSI_MS_OBO_URL)
    version = extract_obo_version(text)
    terms = parse_obo(text, source="MS")
    write_terms_parquet(terms, "psi_ms.parquet", source="psi_ms", version=version)
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pyarrow as pa
import pyarrow.parquet as pq

from qpx.core.cv_terms import (
    CV_HIGHER_BETTER,
    CV_LOWER_BETTER,
    SCORE_PARENT_IDS,
)
from qpx.core.http import safe_urlopen

logger = logging.getLogger(__name__)

# Default OBO source URL
PSI_MS_OBO_URL = "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"

# Arrow schema for ontology Parquet files
TERMS_SCHEMA = pa.schema(
    [
        pa.field("accession", pa.string()),
        pa.field("name", pa.string()),
        pa.field("normalized_name", pa.string()),
        pa.field("definition", pa.string()),
        pa.field("source", pa.string()),
        pa.field("is_score", pa.bool_()),
        pa.field("higher_better", pa.bool_()),
        pa.field("is_a", pa.list_(pa.string())),
        pa.field("synonyms", pa.list_(pa.string())),
        pa.field("is_obsolete", pa.bool_()),
    ]
)


@dataclass(frozen=True)
class CVTerm:
    """A controlled vocabulary term parsed from an OBO file."""

    accession: str
    name: str
    normalized_name: str
    definition: str
    source: str  # "MS" or "PRIDE"
    is_score: bool
    higher_better: Optional[bool]  # None if direction unknown
    is_a: tuple[str, ...] = ()
    synonyms: tuple[str, ...] = ()
    is_obsolete: bool = False


# ---------------------------------------------------------------------------
# OBO parsing
# ---------------------------------------------------------------------------


def _normalize_name(raw_name: str) -> str:
    """Normalize a CV term name to snake_case."""
    name = raw_name.lower()
    name = re.sub(r"[^a-z0-9]+", "_", name)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")


def extract_obo_version(text: str) -> str:
    """Extract the data-version from an OBO file header.

    Returns ``"unknown"`` if no version line is found.
    """
    for line in text.split("\n", 50):  # version is near the top
        if line.startswith("data-version:"):
            return line.split(":", 1)[1].strip()
        if line.startswith("remark: version:"):
            return line.split("version:", 1)[1].strip()
    return "unknown"


def parse_obo(text: str, source: str = "MS") -> list[CVTerm]:
    """Parse an OBO-format string into a list of ``CVTerm`` objects.

    All ``[Term]`` blocks are processed, including obsolete terms
    (flagged via ``is_obsolete=True``). ``[Typedef]`` blocks are skipped.
    """
    terms: list[CVTerm] = []
    blocks = text.split("\n[Term]\n")

    for block in blocks[1:]:  # skip preamble before first [Term]
        accession = None
        name = None
        definition = ""
        is_a_parents: list[str] = []
        higher_better: Optional[bool] = None
        synonyms: list[str] = []
        is_obsolete = False

        for line in block.split("\n"):
            line = line.strip()
            if line.startswith("["):
                break  # hit next section type (e.g. [Typedef])
            if line.startswith("id: "):
                accession = line[4:].strip()
            elif line.startswith("name: "):
                name = line[6:].strip()
            elif line.startswith("def: "):
                m = re.match(r'^def:\s*"(.+?)"\s*\[', line)
                if m:
                    definition = m.group(1)
            elif line.startswith("is_a: "):
                parent_id = line[6:].split("!")[0].strip()
                is_a_parents.append(parent_id)
            elif line.startswith("relationship: has_order "):
                order_acc = line.split("has_order ")[1].split("!")[0].strip()
                if order_acc == CV_HIGHER_BETTER:
                    higher_better = True
                elif order_acc == CV_LOWER_BETTER:
                    higher_better = False
            elif line.startswith("synonym: "):
                m = re.match(r'^synonym:\s*"(.+?)"', line)
                if m:
                    synonyms.append(m.group(1))
            elif line.startswith("is_obsolete: true"):
                is_obsolete = True

        if accession and name:
            is_score = bool(set(is_a_parents) & SCORE_PARENT_IDS)
            terms.append(
                CVTerm(
                    accession=accession,
                    name=name,
                    normalized_name=_normalize_name(name),
                    definition=definition,
                    source=source,
                    is_score=is_score,
                    higher_better=higher_better,
                    is_a=tuple(is_a_parents),
                    synonyms=tuple(synonyms),
                    is_obsolete=is_obsolete,
                )
            )

    return terms


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------


def download_obo(url: str, timeout: int = 30) -> str:
    """Download an OBO file from a URL."""
    logger.info("Downloading OBO from %s ...", url)
    with safe_urlopen(url, timeout=timeout) as resp:
        return resp.read().decode("utf-8")


# ---------------------------------------------------------------------------
# Arrow / Parquet I/O
# ---------------------------------------------------------------------------


def terms_to_arrow(terms: list[CVTerm]) -> pa.Table:
    """Convert a list of CVTerm objects to a PyArrow Table."""
    arrays = {
        "accession": pa.array([t.accession for t in terms], type=pa.string()),
        "name": pa.array([t.name for t in terms], type=pa.string()),
        "normalized_name": pa.array([t.normalized_name for t in terms], type=pa.string()),
        "definition": pa.array([t.definition for t in terms], type=pa.string()),
        "source": pa.array([t.source for t in terms], type=pa.string()),
        "is_score": pa.array([t.is_score for t in terms], type=pa.bool_()),
        "higher_better": pa.array([t.higher_better for t in terms], type=pa.bool_()),
        "is_a": pa.array([list(t.is_a) for t in terms], type=pa.list_(pa.string())),
        "synonyms": pa.array([list(t.synonyms) for t in terms], type=pa.list_(pa.string())),
        "is_obsolete": pa.array([t.is_obsolete for t in terms], type=pa.bool_()),
    }
    return pa.table(arrays, schema=TERMS_SCHEMA)


def write_terms_parquet(
    terms: list[CVTerm],
    path: str | Path,
    source: str,
    version: str,
) -> Path:
    """Write CV terms to a Parquet file with ontology metadata in the footer.

    Args:
        terms: Parsed CV terms.
        path: Output Parquet file path.
        source: Ontology source identifier (e.g., ``"psi_ms"``).
        version: Ontology version string.

    Returns:
        Path to the written file.
    """
    from datetime import datetime, timezone

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    table = terms_to_arrow(terms)

    metadata = {
        b"qpx_ontology_source": source.encode(),
        b"qpx_ontology_version": version.encode(),
        b"qpx_built_at": datetime.now(timezone.utc).isoformat().encode(),
        b"qpx_num_terms": str(len(terms)).encode(),
    }
    existing = table.schema.metadata or {}
    merged = {**existing, **metadata}
    table = table.replace_schema_metadata(merged)

    pq.write_table(table, path, compression="zstd")
    logger.info(
        "Wrote %d terms to %s (source=%s, version=%s)",
        len(terms),
        path,
        source,
        version,
    )
    return path


def read_parquet_ontology_metadata(path: str | Path) -> dict[str, str]:
    """Read qpx_ontology_* metadata from a Parquet file footer."""
    pf = pq.read_metadata(path)
    meta = pf.schema.to_arrow_schema().metadata or {}
    result = {}
    for key, value in meta.items():
        k = key.decode() if isinstance(key, bytes) else key
        if k.startswith("qpx_ontology_"):
            result[k] = value.decode() if isinstance(value, bytes) else value
    return result
