"""Fill protein properties from a FASTA when the producer did not record them.

Some producers never report a protein's sequence (DIA-NN), and some consensusXML
files carry protein hits without it (MSV000085836). Given the FASTA used for the
search, this module fills, for target rows only and only where the value is
null:

- ``pg.molecular_weight``   — average mass of the anchor protein, in kDa
- ``pg.sequence_coverage``  — percent of the anchor covered by the dataset's
  target peptides whose protein evidence includes it
- ``feature.pg_positions``  — every one-based occurrence of the peptide in each
  member of its protein group

A value the producer recorded is never overwritten. The FASTA is optional: a
protein absent from it (a DIA-NN internal decoy, a contaminant from another
database, an isoform not in this FASTA) keeps a null value, and the result reports
how many rows could not be matched so a wrong FASTA is visible rather than
silently partial.
"""

from __future__ import annotations

import gzip
import hashlib
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from qpx.core.protein_sequence import (
    average_molecular_weight_kda,
    peptide_occurrences,
    sequence_coverage_percent,
)

logger = logging.getLogger(__name__)

# Union of the decoy conventions used across qpx converters and quantms databases.
_DECOY_PREFIXES = ("DECOY", "REV_", "RANDOM_", "XXX_")


def _is_decoy_identifier(identifier: str) -> bool:
    """True for a FASTA entry whose identifier, or its accession field, is a decoy."""
    upper = identifier.upper()
    if upper.startswith(_DECOY_PREFIXES):
        return True
    parts = identifier.split("|")
    return len(parts) >= 2 and parts[1].upper().startswith(_DECOY_PREFIXES)


class FastaSequences:
    """Protein sequences from a FASTA, looked up by full identifier or accession.

    ``>sp|P12345|NAME ...`` is reachable as ``sp|P12345|NAME`` (OpenMS) and as
    ``P12345`` (DIA-NN). Decoy entries are skipped: a quantms target-decoy database
    carries ``DECOY_sp|P12345|NAME``, whose accession field would otherwise collide
    with the target and erase it as a conflict. A key that maps to two different
    sequences is dropped rather than guessed.
    """

    def __init__(self) -> None:
        self._sequences: dict[str, str] = {}
        self._conflicting: set[str] = set()
        self.entries = 0
        self.decoy_entries = 0

    @classmethod
    def from_path(cls, path: str | Path) -> FastaSequences:
        """Index a FASTA (plain or ``.gz``); multi-line sequences are joined."""
        index = cls()
        opener = gzip.open if str(path).endswith(".gz") else open
        identifier: str | None = None
        chunks: list[str] = []
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    index._add(identifier, chunks)
                    header = line[1:].strip()
                    identifier = header.split(None, 1)[0] if header else ""
                    chunks = []
                elif identifier is not None:
                    chunks.append(line.strip())
        index._add(identifier, chunks)
        if index._conflicting:
            logger.warning(
                "%d FASTA identifier(s) map to different sequences and are left unresolved, e.g. %s",
                len(index._conflicting),
                sorted(index._conflicting)[:3],
            )
        return index

    def _add(self, identifier: str | None, chunks: list[str]) -> None:
        if not identifier:
            return
        self.entries += 1
        if _is_decoy_identifier(identifier):
            self.decoy_entries += 1
            return
        sequence = "".join(chunks).upper().rstrip("*")
        if not sequence:
            return
        keys = {identifier}
        parts = identifier.split("|")
        if len(parts) >= 2 and parts[1]:
            keys.add(parts[1])
        for key in keys:
            if key in self._conflicting:
                continue
            known = self._sequences.get(key)
            if known is None:
                self._sequences[key] = sequence
            elif known != sequence:
                del self._sequences[key]
                self._conflicting.add(key)

    def __len__(self) -> int:
        return len(self._sequences)

    def get(self, accession: str | None) -> str | None:
        """The sequence for ``accession``, or None when absent or ambiguous."""
        if not accession or accession in self._conflicting:
            return None
        sequence = self._sequences.get(accession)
        if sequence is None and "|" in accession:
            parts = accession.split("|")
            if len(parts) >= 2 and parts[1] not in self._conflicting:
                sequence = self._sequences.get(parts[1])
        return sequence

    def canonical_accession(self, accession: str) -> str:
        """Share an evidence key only for an unambiguous target FASTA accession."""
        if _is_decoy_identifier(accession):
            return accession
        parts = accession.split("|")
        if len(parts) >= 2 and parts[1] in self._sequences:
            return parts[1]
        return accession


@dataclass
class ProteinPropertiesReport:  # pylint: disable=too-many-instance-attributes
    """What was filled, and what could not be matched to the FASTA (plain counters)."""

    fasta_entries: int = 0
    fasta_decoy_entries: int = 0
    pg_rows: int = 0
    pg_rows_eligible: int = 0
    pg_coverage_filled: int = 0
    pg_molecular_weight_filled: int = 0
    pg_anchors_not_in_fasta: int = 0
    feature_rows: int = 0
    feature_rows_eligible: int = 0
    feature_positions_filled: int = 0
    feature_rows_without_match: int = 0
    unmatched_examples: list[str] = field(default_factory=list)
    coverage_evidence: list[str] = field(default_factory=list)

    @property
    def pg_match_rate(self) -> float | None:
        """Share of eligible target pg rows whose anchor was found in the FASTA."""
        if not self.pg_rows_eligible:
            return None
        return 1.0 - self.pg_anchors_not_in_fasta / self.pg_rows_eligible


def _note_unmatched(report: ProteinPropertiesReport, accession: str) -> None:
    if len(report.unmatched_examples) < 5 and accession not in report.unmatched_examples:
        report.unmatched_examples.append(accession)


def _accession_lists(column: pa.ChunkedArray) -> list[tuple[str, ...] | None]:
    """Per-row accession tuples from a ``list<struct<accession,...>>`` column.

    Vectorised: flattening the struct field avoids materialising every struct
    as a Python dict, which dominates the cost on a 20M-row feature view.
    """
    import pyarrow.compute as pc

    out: list[tuple[str, ...] | None] = []
    for chunk in column.chunks:
        if len(chunk) == 0:
            continue
        offsets = chunk.offsets.to_numpy()
        accessions = pc.struct_field(chunk.flatten(), "accession").to_pylist()  # pylint: disable=no-member
        validity = chunk.is_valid().to_pylist()
        base = offsets[0]
        for row, valid in enumerate(validity):
            if not valid:
                out.append(None)
                continue
            lo, hi = offsets[row] - base, offsets[row + 1] - base
            out.append(tuple(a for a in accessions[lo:hi] if a))
    return out


# Evidence column -> SQL expression yielding that row's protein accessions. Fixed
# templates with the file bound as a parameter: no query text is built from input.
_EVIDENCE_QUERIES = {
    "protein_accessions": (
        "SELECT sequence, list_distinct(protein_accessions) FROM read_parquet($1) "
        "WHERE {decoy} sequence IS NOT NULL AND protein_accessions IS NOT NULL GROUP BY ALL"
    ),
    "pg_accessions": (
        "SELECT sequence, list_distinct(list_transform(pg_accessions, x -> x.accession)) FROM read_parquet($1) "
        "WHERE {decoy} sequence IS NOT NULL AND pg_accessions IS NOT NULL GROUP BY ALL"
    ),
    "pg_positions": (
        "SELECT sequence, list_distinct(list_transform(pg_positions, x -> x.protein_accession)) FROM read_parquet($1) "
        "WHERE {decoy} sequence IS NOT NULL AND pg_positions IS NOT NULL GROUP BY ALL"
    ),
}
_TARGET_ONLY = "NOT coalesce(is_decoy, false) AND"


def _evidence_query(column: str, has_decoy_flag: bool) -> str:
    """The fixed evidence query for ``column``, target-restricted when the view has ``is_decoy``."""
    template = _EVIDENCE_QUERIES[column]
    return template.replace("{decoy}", _TARGET_ONLY if has_decoy_flag else "")


def _view_source(dataset_dir: Path, prefix: str, view: str) -> str | None:
    """A flat ``<prefix>.<view>.parquet`` file, else a partitioned ``<view>/`` directory glob.

    ``Dataset`` reads either layout, so evidence must too: a partitioned PSM view
    that is skipped leaves coverage built from group membership alone, the
    undercount PSM evidence exists to prevent.
    """
    flat = dataset_dir / f"{prefix}.{view}.parquet"
    if flat.is_file():
        return str(flat)
    directory = dataset_dir / view
    if directory.is_dir() and any(directory.rglob("*.parquet")):
        return str(directory / "**" / "*.parquet")
    return None


def _source_schema_names(source: str) -> set[str]:
    """Column names of a flat file or of the first file behind a partition glob."""
    if "*" not in source:
        return set(pq.read_schema(source).names)
    first = next(Path(source.split("**", 1)[0]).rglob("*.parquet"))
    return set(pq.read_schema(first).names)


def _peptide_protein_candidates(
    feature_path: str | Path | None, psm_path: str | Path | None = None, used: list[str] | None = None
) -> dict[str, set[str]]:
    """Map each target peptide sequence to every protein its evidence names.

    Sources, unioned: the PSM view's ``protein_accessions`` (every protein an
    identification maps to, shared peptides included), the feature view's group
    memberships, and positions a producer already recorded.

    Coverage needs the PSM evidence. Group membership alone drops every peptide
    shared across groups — on PXD000612 that put ACTB at 4.5% coverage against the
    84.3% OpenMS recorded, since actin peptides are mostly shared with other
    actins. Positions do not use these candidates; they stay tied to the group.
    """
    import duckdb

    candidates: dict[str, set[str]] = defaultdict(set)
    con = duckdb.connect()
    try:
        for view, source in (("feature", feature_path), ("psm", psm_path)):
            if source is None or ("*" not in str(source) and not Path(source).is_file()):
                continue
            if _collect_evidence(con, str(source), candidates) and used is not None:
                used.append(view)
    finally:
        con.close()
    return candidates


def _collect_evidence(con, source: str, candidates: dict[str, set[str]]) -> bool:
    """Add one view's peptide -> protein evidence to ``candidates``; True if it had any."""
    columns = _source_schema_names(source)
    if "sequence" not in columns:
        return False
    found = False
    for column in (name for name in _EVIDENCE_QUERIES if name in columns):
        rows = con.execute(_evidence_query(column, "is_decoy" in columns), [source]).fetchall()
        for sequence, accessions in rows:
            candidates[sequence].update(accession for accession in accessions or () if accession)
        found = True
    return found


def _protein_peptides(candidates: dict[str, set[str]], fasta: FastaSequences) -> dict[str, set[str]]:
    proteins: dict[str, set[str]] = defaultdict(set)
    for sequence, accessions in candidates.items():
        for accession in accessions:
            proteins[fasta.canonical_accession(accession)].add(sequence)
    return proteins


def _stamped_schema(schema: pa.Schema) -> tuple[pa.Schema, str]:
    """The source schema with its footer re-stamped (new uuid/date, same identity)."""
    from qpx.writers.base import _stamp_footer_metadata

    metadata = schema.metadata or {}
    compression = metadata.get(b"compression_format", b"zstd").decode() or "zstd"
    stamped = _stamp_footer_metadata(schema.empty_table(), compression).schema
    return stamped, compression


def _rewrite_view(source: Path, destination: Path, fill_batch) -> None:
    """Rewrite a view row group by row group, keeping its schema and footer identity."""
    from qpx.writers.base import parquet_write_options

    parquet = pq.ParquetFile(source)
    schema, compression = _stamped_schema(parquet.schema_arrow)
    # The same encoding every QPX writer uses (byte-stream-split on rt/mz leaves,
    # per-column dictionaries, zstd level, format 2.6). A bare compression= kept
    # the footer identity but silently re-encoded the file with pyarrow defaults.
    with pq.ParquetWriter(str(destination), schema, **parquet_write_options(schema, compression)) as writer:
        for group in range(parquet.num_row_groups):
            table = parquet.read_row_group(group)
            writer.write_table(fill_batch(table).cast(schema))


def _fill_pg_rows(anchors, decoys, columns, fasta: FastaSequences, protein_peptides, report: ProteinPropertiesReport) -> None:
    """Fill missing target-row properties using the same accession keys as the evidence."""
    cache: dict[str, dict[str, float | None]] = {}
    for row, anchor in enumerate(anchors):
        missing = [name for name, values in columns.items() if values[row] is None]
        if decoys[row] or not anchor or not missing:
            continue
        report.pg_rows_eligible += 1
        sequence = fasta.get(anchor)
        if sequence is None:
            report.pg_anchors_not_in_fasta += 1
            _note_unmatched(report, anchor)
            continue
        key = fasta.canonical_accession(anchor)
        if key not in cache:
            spans = [span for peptide in protein_peptides.get(key, ()) for span in peptide_occurrences(peptide, sequence)]
            cache[key] = {
                "sequence_coverage": sequence_coverage_percent(sequence, spans),
                "molecular_weight": average_molecular_weight_kda(sequence),
            }
        for name in missing:
            columns[name][row] = cache[key][name]


def _fill_pg_table(table: pa.Table, fasta: FastaSequences, protein_peptides, report: ProteinPropertiesReport) -> pa.Table:
    names = table.schema.names
    report.pg_rows += table.num_rows
    columns = {name: table.column(name).to_pylist() for name in ("sequence_coverage", "molecular_weight") if name in names}
    if not columns:
        return table
    anchors = table.column("anchor_protein").to_pylist()
    decoys = table.column("is_decoy").to_pylist() if "is_decoy" in names else [False] * table.num_rows
    _fill_pg_rows(anchors, decoys, columns, fasta, protein_peptides, report)
    for name, values in columns.items():
        index = names.index(name)
        column_field = table.schema.field(index)
        filled = pa.array(values, type=column_field.type)
        count = table.column(name).null_count - filled.null_count
        if name == "sequence_coverage":
            report.pg_coverage_filled += count
        else:
            report.pg_molecular_weight_filled += count
        table = table.set_column(index, column_field, filled)
    return table


def _fill_feature_table(table: pa.Table, fasta: FastaSequences, report: ProteinPropertiesReport) -> pa.Table:
    names = table.schema.names
    report.feature_rows += table.num_rows
    if "pg_positions" not in names or "pg_accessions" not in names:
        return table
    positions_column = table.column("pg_positions")
    if positions_column.null_count == 0:
        return table
    positions = positions_column.to_pylist()
    sequences = table.column("sequence").to_pylist()
    decoys = table.column("is_decoy").to_pylist() if "is_decoy" in names else [False] * table.num_rows
    groups = _accession_lists(table.column("pg_accessions"))
    cache: dict[tuple, list[dict] | None] = {}
    for row, existing in enumerate(positions):
        group = groups[row]
        if existing is not None or decoys[row] or not sequences[row] or not group:
            continue
        report.feature_rows_eligible += 1
        key = (sequences[row], group)
        if key not in cache:
            found: list[dict] = []
            for accession in group:
                protein = fasta.get(accession)
                for start, end in peptide_occurrences(sequences[row], protein):
                    found.append({"protein_accession": accession, "start": start, "end": end})
            cache[key] = found or None
        if cache[key] is None:
            report.feature_rows_without_match += 1
            continue
        positions[row] = cache[key]
        report.feature_positions_filled += 1
    index = names.index("pg_positions")
    return table.set_column(index, table.schema.field(index), pa.array(positions, type=table.schema.field(index).type))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _provenance_step(report: ProteinPropertiesReport, fasta_path: Path, views: list[str], step_order: int) -> dict:
    from qpx._version import __version__

    parameters = {
        "fasta": fasta_path.name,
        "fasta_sha256": _sha256(fasta_path),
        "fill_policy": "null values on target rows only; producer values are never overwritten",
        "pg_sequence_coverage_filled": report.pg_coverage_filled,
        "pg_molecular_weight_filled": report.pg_molecular_weight_filled,
        "pg_anchors_not_in_fasta": report.pg_anchors_not_in_fasta,
        "feature_pg_positions_filled": report.feature_positions_filled,
        "feature_rows_without_match": report.feature_rows_without_match,
    }
    return {
        "step_order": step_order,
        "step_category": "annotation",
        "step_name": "protein_properties_from_fasta",
        "tool_name": "qpx",
        "tool_version": __version__,
        "parameters": [{"key": key, "value": str(value)} for key, value in parameters.items()],
        "output_views": views,
    }


def _write_provenance(source: Path | None, destination: Path, step: dict) -> None:
    from qpx.writers.provenance import ProvenanceWriter

    steps = pq.read_table(source).to_pylist() if source is not None and source.is_file() else []
    step["step_order"] = max((s.get("step_order") or 0 for s in steps), default=0) + 1
    with ProvenanceWriter(destination, creator="qpx") as writer:
        writer.write_batch([*steps, step])


def annotate_protein_properties(
    dataset_dir: str | Path,
    prefix: str,
    fasta_path: str | Path,
    staging: str | Path,
) -> tuple[ProteinPropertiesReport, list[str]]:
    """Write FASTA-annotated pg/feature (and provenance) views for a dataset into ``staging``.

    Returns the report and the file names written, which the caller moves into
    place. Nothing in ``dataset_dir`` is modified here.
    """
    dataset_dir, staging, fasta_path = Path(dataset_dir), Path(staging), Path(fasta_path)
    fasta = FastaSequences.from_path(fasta_path)
    report = ProteinPropertiesReport(fasta_entries=fasta.entries, fasta_decoy_entries=fasta.decoy_entries)
    written: list[str] = []
    views: list[str] = []

    pg_path = dataset_dir / f"{prefix}.pg.parquet"
    feature_path = dataset_dir / f"{prefix}.feature.parquet"
    candidates = _peptide_protein_candidates(
        _view_source(dataset_dir, prefix, "feature"),
        _view_source(dataset_dir, prefix, "psm"),
        used=report.coverage_evidence,
    )

    if pg_path.is_file():
        protein_peptides = _protein_peptides(candidates, fasta)
        name = pg_path.name
        _rewrite_view(pg_path, staging / name, lambda t: _fill_pg_table(t, fasta, protein_peptides, report))
        written.append(name)
        views.append("pg")
    if feature_path.is_file():
        name = feature_path.name
        _rewrite_view(feature_path, staging / name, lambda t: _fill_feature_table(t, fasta, report))
        written.append(name)
        views.append("feature")

    if written:
        provenance = dataset_dir / f"{prefix}.provenance.parquet"
        name = provenance.name
        _write_provenance(provenance, staging / name, _provenance_step(report, fasta_path, views, 0))
        written.append(name)

    match_rate = report.pg_match_rate
    if match_rate is not None and match_rate < 0.5:
        logger.warning(
            "Only %.1f%% of target protein groups were found in %s; is this the FASTA used for the search? "
            "Unmatched examples: %s",
            100 * match_rate,
            fasta_path.name,
            report.unmatched_examples,
        )
    return report, written
