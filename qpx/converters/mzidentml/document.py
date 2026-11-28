"""One-pass reading of an mzIdentML file.

``SpectrumIdentificationResult`` elements, by far the largest part of a file, are
yielded one at a time and freed. The sections they need come before them in the
schema's element order, and the sections after them stay in the tree, so one parse
serves the PSM, protein-group and provenance views.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any, Iterator

try:
    from lxml import etree
except ImportError:
    etree: Any = None

# Decoy prefixes for files that mark no PeptideEvidence as a decoy (Kojak writes DECOY1_).
DECOY_PREFIXES = ("DECOY", "decoy_", "rev_")


def is_true(value: str) -> bool:
    """Read an xsd:boolean attribute, which a writer may spell true or 1 (OpenPepXL)."""
    return value.lower() in ("true", "1")


class MzidDocument:
    """The tree of a streamed mzIdentML file, without its results and, once released, its sequences."""

    def __init__(self) -> None:
        self.root = None
        self.ns = ""
        # DBSequence id -> accession, which outlives the SequenceCollection.
        self.db_accessions: dict[str, str] = {}
        # SpectrumIdentificationList id -> SpectraData references of its results, in file order.
        self.spectra_refs: dict[str, list[str]] = {}
        # DBSequence ids of decoys: those a decoy PeptideEvidence refers to, or by prefix in a file that marks none.
        self.decoy_sequences: set[str] = set()

    def attach(self, root) -> None:
        """Take the tree as parsed so far: everything up to the first result."""
        self.root = root
        self.ns = root.nsmap.get(None, "")
        accessions = {d.get("id"): d.get("accession", d.get("id")) for d in root.iter(f"{{{self.ns}}}DBSequence")}
        self.decoy_sequences = {
            evidence.get("dBSequence_ref")
            for evidence in root.iter(f"{{{self.ns}}}PeptideEvidence")
            if is_true(evidence.get("isDecoy", "false"))
        } or {ref for ref, accession in accessions.items() if accession.startswith(DECOY_PREFIXES)}
        # A decoy never shares its target's accession: one that does (xiFDR) gets DECOY_, as DIA-NN's decoys do.
        self.db_accessions = {
            ref: f"DECOY_{accession}" if ref in self.decoy_sequences and not accession.startswith(DECOY_PREFIXES) else accession
            for ref, accession in accessions.items()
        }

    def release_sequences(self) -> None:
        """Free the SequenceCollection; a caller keeps what it needs from it first."""
        for collection in self.root.iter(f"{{{self.ns}}}SequenceCollection"):
            collection.clear()


def stream_results(path: str | Path, document: MzidDocument) -> Iterator:
    """Yield each SpectrumIdentificationResult of *path* and free it once the caller moves on."""
    if etree is None:
        raise ImportError("lxml is required for mzIdentML conversion. Install it with: pip install qpx[mzidentml]")
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as source:
        context = etree.iterparse(source, events=("end",), tag="{*}SpectrumIdentificationResult")
        for _, result in context:
            if document.root is None:
                document.attach(result.getroottree().getroot())
            parent = result.getparent()
            refs = document.spectra_refs.setdefault(parent.get("id", ""), [])
            if result.get("spectraData_ref") not in refs:
                refs.append(result.get("spectraData_ref"))
            yield result
            result.clear()
            # The parser may still point at the current element, so only earlier siblings go.
            while result.getprevious() is not None:
                del parent[0]
        if document.root is None:
            document.attach(context.root)
