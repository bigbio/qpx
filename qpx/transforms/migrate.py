"""Bring a QPX dataset from an older specification to the current one, without its source files."""

from __future__ import annotations

import logging
import os
import re
import shutil
import uuid
from functools import lru_cache
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from qpx.converters.channel_labels import normalize_label
from qpx.converters.ontology_util import dedupe_ontology_entries, reconcile_ontology_entries, refresh_ontology_entries
from qpx.converters.ptm import _add_nterm_dash, _unimod_title_of, name_to_unimod, normalize_peptidoform
from qpx.converters.utils import canonical_accession
from qpx.core.engine import create_connection
from qpx.core.scores import field_ontology_entries
from qpx.core.semantic_validation import _CANONICAL, _TAGS
from qpx.core.sql import sql_build, validate_identifier
from qpx.core.struct_lists import file_layouts, restore_table
from qpx.dataset import Dataset
from qpx.transforms.accession_normalizer import _LIST_ACCESSION_COLUMNS, _SCALAR_ACCESSION_COLUMNS, ACCESSION_VIEWS, _entry_names
from qpx.transforms.mantissa import MANTISSA_BITS_FOOTER_KEY
from qpx.transforms.producer_rows import _diann_feature_row, _irt_row, _msnet_row, _pg_row_fix, _raw_msnet_row
from qpx.transforms.protein_properties import _rewrite_view
from qpx.transforms.relayout import ID_REFERENCES, ORDINAL_IDS, _quoted
from qpx.version import QPX_SPEC_VERSION, parse_spec_version
from qpx.writers.xic import XIC_FOOTER_KEYS

logger = logging.getLogger(__name__)

_BATCH_ROWS = 50_000
# Derived again from the rewritten rows, and computed links that are not stored any more.
_DROPPED_ID_COLUMNS = ("feature_id", "psm_id", "pg_id", "psm_ids", "pg_ids")
# The columns that hold another view's ids and follow them; the feature's own links are dropped.
_ID_REFERENCES = {view: columns for view, columns in ID_REFERENCES.items() if view != "feature"}
# The feature composite each producer declares (docs/spec/identity.md), by footer `creator`.
_PRODUCER_FEATURE_COMPOSITES = {
    "diann": ("peptidoform", "charge", "run_file_name", "is_decoy", "is_library_decoy"),
    "spectronaut": ("peptidoform", "charge", "run_file_name"),
    "cdap": ("peptidoform", "charge", "run_file_name"),
    "maxquant": ("peptidoform", "charge", "run_file_name", "rt_start", "rt_stop", "scan"),
    "quantms-msstats": ("run_file_name", "peptidoform", "charge", "rt", "scan"),
    "fragpipe": ("quantification_unit_id", "peptidoform", "charge", "compensation_voltage"),
}
# These producers supply PSM ids that cannot be re-derived from the persisted columns.
_PRODUCER_PSM_IDS = {"mzidentml"}
# Non-conformant label-free labels of 1.0.2: an erratum, not an older vocabulary.
_LABEL_ERRATA = {"raw": "LFQ", "precursor_area": "LFQ"}
_FILE_LEVEL_DATASET_COLUMNS = ("file_checksums", "file_row_counts", "file_sizes_bytes", "total_structures", "packaged_at")
_NOT_MIGRATED_SUFFIXES = (".h5mu", ".h5ad")
# The accession lists of normalize-accessions, and feature.pg_positions, which names proteins too.
_ACCESSION_LISTS = {**_LIST_ACCESSION_COLUMNS, "pg_positions": "protein_accession"}
# What the semantic validator removes from a peptidoform to compare it with its sequence.
_MODIFICATION_MARKS = re.compile(rf"{_TAGS}|[-?]")
# Run names old converters spelled apart from run.parquet: with the suffix of a quantms idXML, without the
# dissociation suffix of a dataset split by it, or without a "%".
_IDXML_RUN_SUFFIX = "_consensus_fdr_filter.idXML"
_SPLIT_RUN_SUFFIXES = ("_CID", "_HCD")
_RUN_NAMES_SQL = "SELECT DISTINCT $column FROM read_parquet($source)"
# openms-consensus spelled the PSM q-value score `q-value`; score names are snake_case (ontology.md).
_PSM_SCORE_NAMES = {"q-value": "qvalue"}
# The DIA-NN scores of QPX 1.0.2 that a migration moves to the column they repeat (producer_rows.py).
_MOVED_SCORES = {("feature", "qvalue"): "precursor_qvalue", ("pg", "qvalue"): "pg_qvalue"}
# The first specification whose DIA-NN observed_mz is the measured MS1 apex rather than Precursor.Mz.
_MEASURED_DIANN_MZ = (1, 3)
# The rows that repeat an id, except the first, or for a feature the first copy holding its group's
# q-value: QPX 1.0.2 repeated a quantms feature once per q-value of its first protein.
_REPEATED_ROWS_SQL = (
    "SELECT position FROM (SELECT v.file_row_number AS position, row_number() OVER "
    "(PARTITION BY v.$id ORDER BY $preferred DESC, v.file_row_number) AS copy "
    "FROM read_parquet($target, file_row_number=true) AS v $join "
    "WHERE v.$id IN (SELECT $id FROM read_parquet($target) GROUP BY $id HAVING count(*) > 1)) "
    "WHERE copy > 1 ORDER BY position"
)
_GROUP_QVALUE_JOIN = (
    "LEFT JOIN (SELECT list_sort(pg_accessions) AS members, min(global_qvalue) AS qvalue "
    "FROM read_parquet($pg) GROUP BY members) AS g "
    "ON list_sort([member.accession FOR member IN v.pg_accessions]) = g.members"
)


def _label(value):
    return value if value is None else normalize_label(_LABEL_ERRATA.get(value, value))


def _entries(values):
    return values and [{**entry, "label": _label(entry["label"])} for entry in values]


def _footer(path: Path) -> dict[bytes, bytes]:
    return pq.read_schema(path).metadata or {}


def _batches(path: Path):
    """The records of *path*, with its struct lists read back as lists."""
    source = pq.ParquetFile(path)
    layouts = file_layouts(source.metadata.metadata)
    try:
        with create_connection() as con:
            for batch in source.iter_batches(_BATCH_ROWS):
                yield restore_table(con, pa.Table.from_batches([batch]), layouts).to_pylist()
    finally:
        source.close()


def _keep_schema_fields(writer, records: list[dict]) -> list[dict]:
    """Drop columns the current schema no longer has, naming them once."""
    names = set(writer.arrow_schema.names)
    dropped = {key for record in records for key in record} - names
    if dropped:
        logger.warning("Dropping columns outside the schema: %s", sorted(dropped))
    return [{key: value for key, value in record.items() if key in names} for record in records]


@lru_cache(maxsize=65536)
def _accession(value: str) -> str:
    return canonical_accession(value)


def _member(value, field: str):
    """A list member with a bare accession: a string, or a struct whose *field* holds it."""
    if isinstance(value, str):
        return _accession(value)
    if value and value.get(field):
        return {**value, field: _accession(value[field])}
    return value


def _canonical_accessions(row: dict) -> dict:
    """Bare UniProt accessions, as the converters write them; a null pg_names takes the members' entry names."""
    row = dict(row)
    members = row.get("pg_accessions")
    if "pg_names" in row and row["pg_names"] is None and members and isinstance(members[0], str):
        row["pg_names"] = _entry_names(members)
    for column in _SCALAR_ACCESSION_COLUMNS:
        if row.get(column):
            row[column] = _accession(row[column])
    for column, field in _ACCESSION_LISTS.items():
        if row.get(column):
            row[column] = [_member(value, field) for value in row[column]]
    return row


def _listed_scores(position: dict) -> dict:
    """The scores of a position as a list: writer 1.2.1.dev5 wrote one struct, with an empty name for none."""
    scores = position.get("scores")
    if isinstance(scores, dict):
        return {**position, "scores": [scores] if scores["score_name"] else None}
    return position


def _modification(entry: dict) -> dict:
    """The Unimod title and accession the converters write, where the old entry named only one of them.

    Writer 1.2.1.dev5 wrote a missing accession as an empty string.
    """
    name, accession = entry["name"], entry["accession"] or None
    positions = entry.get("positions") and [_listed_scores(position) for position in entry["positions"]]
    if accession is None and name:
        number = name_to_unimod(name, (positions or [{}])[0].get("amino_acid"))
        accession = None if number is None else f"UNIMOD:{number}"
    if name and name.startswith("UNIMOD:"):
        name = _unimod_title_of(name)
    return {**entry, "name": name, "accession": accession, "positions": positions}


def _peptidoform(peptidoform: str) -> str:
    """The canonical peptidoform; MSNet wrote ProForma tags after the mzTab dot of the N-terminus, .[Acetyl]PEPTIDE."""
    if peptidoform.startswith("."):
        peptidoform = _add_nterm_dash(peptidoform.removeprefix("."))
    return normalize_peptidoform(peptidoform)


@lru_cache(maxsize=65536)
def _residues(peptidoform: str) -> str | None:
    """The sequence of a canonical peptidoform; quantms 1.0.2 TMT features spelled modification names into it."""
    return _MODIFICATION_MARKS.sub("", peptidoform) if _CANONICAL.fullmatch(peptidoform) else None


def _repaired_sequence(peptidoform: str, sequence: str | None) -> dict:
    """The sequence read from the peptidoform where they differ, and no missed cleavages counted on the old one."""
    residues = _residues(peptidoform)
    return {"sequence": residues, "missed_cleavages": None} if residues and residues != sequence else {}


def _identity_row(record: dict) -> dict:
    """The rewritten peptidoform, modifications and labels of a psm or feature row, without its old id."""
    row = {key: value for key, value in record.items() if key not in _DROPPED_ID_COLUMNS}
    if record.get("feature_id") is not None and "psm_id" in record:
        row["feature_id"] = record["feature_id"]
    if row.get("peptidoform"):
        row["peptidoform"] = _peptidoform(row["peptidoform"])
        row.update(_repaired_sequence(row["peptidoform"], row.get("sequence")))
    if row.get("modifications"):
        row["modifications"] = [_modification(entry) for entry in row["modifications"]]
    if "consensus_support" in row:
        row = _msnet_row(row)
    for column in ("calculated_mz", "observed_mz"):
        # 0.0 was the old sentinel for an m/z the producer did not report
        if row.get(column) == 0.0:
            row[column] = None
    for column in ("intensities", "additional_intensities"):
        if column in row:
            row[column] = _entries(row[column])
    if row.get("label") is not None:
        row["label"] = _label(row["label"])
    return row


def _psm_row(record: dict) -> dict:
    """A psm row rewritten like a feature row, with its scores under their snake_case names."""
    row = _identity_row(_raw_msnet_row(record) if "reference_file_name" in record else record)
    if row.get("additional_scores"):
        row["additional_scores"] = [
            {**score, "score_name": _PSM_SCORE_NAMES.get(score["score_name"], score["score_name"])}
            for score in row["additional_scores"]
        ]
    return row


def _flat_pg_rows(record: dict) -> list[dict]:
    """One pg row per label of a pre-1.1 row, keyed by a one-run grouped_runs."""
    row = {key: value for key, value in record.items() if key not in _DROPPED_ID_COLUMNS}
    if "run_file_name" in row:
        row["grouped_runs"] = [row.pop("run_file_name")]
    if "intensities" not in row:
        if row.get("label") is not None:
            row["label"] = _label(row["label"])
        return [row]
    entries = row.pop("intensities")
    additional = row.pop("additional_intensities", None) or []
    if not entries:
        return [{**row, "label": None, "intensity": None, "additional_intensities": None}]
    return [
        {
            **row,
            "label": _label(entry["label"]),
            "intensity": entry["intensity"],
            "additional_intensities": [
                {**item, "label": _label(item["label"])} for item in additional if item["label"] == entry["label"]
            ]
            or None,
        }
        for entry in entries
    ]


def _run_row(record: dict) -> dict:
    return {**record, "samples": _entries(record.get("samples"))}


def _dataset_row(record: dict) -> dict:
    """File-level columns describe the old files, and this is a new analysis id."""
    row = {key: value for key, value in record.items() if key not in _FILE_LEVEL_DATASET_COLUMNS}
    return {**row, "analysis_id": str(uuid.uuid4())}


def _mz_row(record: dict) -> dict:
    scan = record.get("scan")
    return {**record, "scan": [scan] if isinstance(scan, int) else scan}


_ROW_MIGRATIONS = {
    "feature": lambda record: [_identity_row(record)],
    "psm": lambda record: [_psm_row(record)],
    "pg": _flat_pg_rows,
    "run": lambda record: [_run_row(record)],
    "dataset": lambda record: [_dataset_row(record)],
    "mz": lambda record: [_mz_row(record)],
}


def _holds_each_label_once(row: dict) -> bool:
    labels = [entry["label"] for entry in row.get("intensities") or []]
    return len(labels) == len(set(labels))


def _migrated_pairs(view: str, footer: dict[bytes, bytes], records: list[dict]) -> list[tuple[dict, dict]]:
    """Keep each source row beside its migrated rows, including producer-specific identities.

    A feature holding a label twice is two features merged into one row, which cannot be told apart: it is left out.
    """
    if view == "feature" and footer.get(b"creator") == b"diann":
        measured_mz = (parse_spec_version(footer.get(b"qpx_version", b"").decode()) or (0, 0)) >= _MEASURED_DIANN_MZ
        pairs = [(record, _identity_row(_diann_feature_row(record, measured_mz))) for record in records]
    else:
        migrate_row = _ROW_MIGRATIONS.get(view, lambda record: [record])
        pairs = [(record, row) for record in records for row in migrate_row(record)]
    pg_fix = _pg_row_fix(footer) if view == "pg" else None
    if pg_fix is not None:
        pairs = [(record, pg_fix(row)) for record, row in pairs]
    if view in ACCESSION_VIEWS:
        pairs = [(record, _canonical_accessions(row)) for record, row in pairs]
    if view == "feature":
        pairs = [(record, _irt_row(row)) for record, row in pairs if _holds_each_label_once(row)]
    return pairs


def _writer_options(view: str, footer: dict[bytes, bytes]) -> dict:
    creator = footer.get(b"creator", b"qpx").decode()
    options = {"creator": creator}
    if b"scan_format" in footer:
        options["scan_format"] = footer[b"scan_format"].decode()
    if view == "feature":
        declared = footer.get(b"identity_composite")
        composite = tuple(declared.decode().split(",")) if declared else _PRODUCER_FEATURE_COMPOSITES.get(creator)
        # Files written before the decoy flags joined the DIA-NN composite declare its first three fields.
        if creator == "diann":
            composite = _PRODUCER_FEATURE_COMPOSITES[creator]
        if composite is not None:
            options["identity_composite"] = composite
    return options


def _extra_column_options(writer_class, source: Path) -> dict:
    """The writer option keeping the source columns outside the schema, such as the SDRF characteristics of a sample."""
    schema = writer_class.view_schema()
    known = set(schema.get_arrow_schema().names)
    extra = [name for name in pq.read_schema(source).names if name not in known]
    refused = schema.undeclared_columns(extra)
    kept = [name for name in extra if name not in refused]
    if not kept:
        return {}
    logger.info("Keeping the %s columns outside the schema as strings: %s", schema.view_name, kept)
    return {"extra_columns": kept}


def _migrate_ontology(source: Path, target: Path) -> int:
    """One row per field and view, naming only fields the migrated views emit, with the current terms."""
    writer_class, suffix = Dataset.writer_for("ontology")
    prefix = target.name.removesuffix(suffix)
    records = refresh_ontology_entries(dedupe_ontology_entries(pq.read_table(source).to_pylist()), target.parent, prefix)
    entries = _follow_moved_scores(records, reconcile_ontology_entries(records, target.parent, prefix))
    with writer_class(target, **_writer_options("ontology", _footer(source))) as writer:
        writer.write_batch(entries)
    return len(entries)


def _follow_moved_scores(records: list[dict], entries: list[dict]) -> list[dict]:
    """Give the row of a column that took the values of a moved score the source of the score's dropped row."""
    kept = {(entry["view"], entry["field_name"]) for entry in entries}
    moved = {
        (record["view"], _MOVED_SCORES[key]): record
        for record in records
        if (key := (record["view"], record["field_name"])) in _MOVED_SCORES and key not in kept and record["source_column_name"]
    }
    followed = []
    for entry in entries:
        score = moved.get((entry["view"], entry["field_name"]))
        if score is not None and entry["source_column_name"] is None:
            entry = field_ontology_entries(
                entry["view"], {entry["field_name"]: score["source_column_name"]}, score["source_tool"]
            )[0]
        followed.append(entry)
    return followed


def _followed(value, new_ids: dict[int, int], field: str):
    """*value* with the new ids: one id, a list of them, or a list of structs whose *field* holds one."""
    if field:
        # The field cannot be null, so an id the map lacks stays.
        return [{**item, field: new_ids.get(item[field], item[field])} for item in value]
    if isinstance(value, list):
        return [new_ids.get(item) for item in value]
    return new_ids.get(value)


def _follow_ids(pairs: list[tuple[dict, dict]], references: dict[str, dict[int, int]], keep_psm_ids: bool) -> None:
    """Give each migrated row the new ids its references hold; mzIdentML PSM ids stay."""
    for record, row in pairs:
        if keep_psm_ids:
            row["psm_id"] = record["psm_id"]
        for reference, new_ids in references.items():
            column, _, field = reference.partition(".")
            if row.get(column) is not None:
                row[column] = _followed(row[column], new_ids, field)


def _kept_footer(footer: dict[bytes, bytes]) -> dict[bytes, bytes]:
    """The footer keys that still hold: the values stay rounded, and the traces keep their selection."""
    return {key: footer[key] for key in (MANTISSA_BITS_FOOTER_KEY, *XIC_FOOTER_KEYS) if key in footer}


def _run_names(folder: Path) -> set[str]:
    paths = folder.glob(f"*{Dataset.writer_for('run')[1]}")
    return {name for path in paths for name in pq.read_table(path, columns=["run_file_name"]).column(0).to_pylist()}


def _old_spellings(run: str) -> set[str]:
    """The spellings of *run* that old converters wrote in psm files."""
    return {run + _IDXML_RUN_SUFFIX, run.replace("%", ""), *(run.removesuffix(s) for s in _SPLIT_RUN_SUFFIXES)} - {run}


def _run_spellings(view: str, source: Path) -> dict[str, str]:
    """The spelling in run.parquet of the run names of a psm *source* it lacks, when each matches one run.

    A matched run the file names already would merge two runs, so then no name is respelled.
    """
    runs = _run_names(source.parent) if view == "psm" else set()
    if not runs:
        return {}
    column = "reference_file_name" if "reference_file_name" in pq.read_schema(source).names else "run_file_name"
    sql = sql_build(_RUN_NAMES_SQL, column=validate_identifier(column), source=_quoted(source))
    with create_connection() as con:
        names = {name for (name,) in con.sql(sql).fetchall()}
    missing = names - runs - {None}
    if not missing:
        return {}
    matches: dict[str, set[str]] = {}
    for run in runs:
        for spelling in _old_spellings(run):
            matches.setdefault(spelling, set()).add(run)
    spellings = {name: next(iter(matches[name])) for name in missing if len(matches.get(name, ())) == 1}
    if len(spellings) < len(missing) or len(set(spellings.values())) < len(spellings) or set(spellings.values()) & names:
        logger.warning("%d run_file_name values of %s are not in run.parquet", len(missing), source.name)
        return {}
    logger.warning("Respelling %d run_file_name values of %s as run.parquet does", len(spellings), source.name)
    return spellings


def _respell_runs(pairs: list[tuple[dict, dict]], spellings: dict[str, str]) -> None:
    if not spellings:
        return
    for _, row in pairs:
        row["run_file_name"] = spellings.get(row["run_file_name"], row["run_file_name"])


def _migrate_view(view: str, source: Path, target: Path, ids: dict[str, dict[int, int]], *, track: bool) -> int:
    """Migrate the rows of *source* to *target*, following the views already in *ids*.

    With *track*, the view's old ids are added to *ids* with their new ones, matched
    by row position: the writer keeps the order of its rows.
    """
    if view == "ontology":
        return _migrate_ontology(source, target)
    writer_class, _ = Dataset.writer_for(view)
    footer = _footer(source)
    keep_ids = view == "psm" and footer.get(b"creator", b"").decode() in _PRODUCER_PSM_IDS
    references = {column: ids.get(referenced, {}) for column, referenced in _ID_REFERENCES.get(view, {}).items()}
    spellings = _run_spellings(view, source)
    old_ids = []
    rows = source_rows = 0
    options = {**_writer_options(view, footer), **_extra_column_options(writer_class, source)}
    with writer_class(target, **options) as writer:
        # Added before the rows, the keys are in the Arrow schema too.
        if kept := _kept_footer(footer):
            writer.add_footer_metadata(kept)
        for records in _batches(source):
            pairs = _migrated_pairs(view, footer, records)
            _respell_runs(pairs, spellings)
            _follow_ids(pairs, references, keep_ids)
            if track:
                old_ids.extend(record.get(ORDINAL_IDS[view]) for record, _ in pairs)
            writer.write_batch(_keep_schema_fields(writer, [row for _, row in pairs]))
            rows += len(pairs)
            source_rows += len(records)
    if view == "feature" and source_rows > rows:
        logger.warning("Dropped %d features of %s that hold a label twice", source_rows - rows, source.name)
    if track:
        new_ids = pq.read_table(target, columns=[ORDINAL_IDS[view]]).column(0).to_pylist()
        ids.setdefault(view, {}).update((old, new) for old, new in zip(old_ids, new_ids) if old is not None)
    return rows


def _repeated_rows(target: Path, id_field: str) -> list[int]:
    """The positions of the rows of *target* that repeat an id, without the copy that stays."""
    feature_suffix, pg_suffix = Dataset.writer_for("feature")[1], Dataset.writer_for("pg")[1]
    pg = target.with_name(target.name.removesuffix(feature_suffix) + pg_suffix)
    grouped = target.name.endswith(feature_suffix) and pg.is_file()
    sql = sql_build(
        _REPEATED_ROWS_SQL,
        id=validate_identifier(id_field),
        target=_quoted(target),
        preferred="g.members IS NOT NULL AND v.pg_global_qvalue IS NOT DISTINCT FROM g.qvalue" if grouped else "true",
        join=sql_build(_GROUP_QVALUE_JOIN, pg=_quoted(pg)) if grouped else "",
    )
    with create_connection() as con:
        return [position for (position,) in con.sql(sql).fetchall()]


def _drop_repeated_ids(target: Path) -> int:
    """Keep one row per id of *target*, returning the rows dropped."""
    id_field = _footer(target)[b"primary_key"].decode()
    positions = _repeated_rows(target, id_field)
    if not positions:
        return 0
    logger.warning("Dropping %d rows of %s that repeat a %s", len(positions), target.name, id_field)
    dropped = np.asarray(positions, dtype=np.int64)
    offset = 0

    def without_repeats(table: pa.Table) -> pa.Table:
        nonlocal offset
        index = np.arange(offset, offset + table.num_rows, dtype=np.int64)
        offset += table.num_rows
        return table.filter(pa.array(np.isin(index, dropped, invert=True)))

    staging = target.with_name(f".{target.name}.deduplicated")
    try:
        _rewrite_view(target, staging, without_repeats)
        os.replace(staging, target)
    finally:
        staging.unlink(missing_ok=True)
    return len(positions)


def _without_repeated_ids(output: Path, files: dict[str, list[Path]], written: dict[str, int]) -> dict[str, int]:
    """The rows left per view once the repeated feature and psm ids are dropped.

    Runs after every view: a feature copy is chosen by its pg group, and a dropped
    copy shares its id with the kept one, which the references follow.
    """
    dropped = {
        view: sum(_drop_repeated_ids(output / path.name) for path in files[view])
        for view in ("feature", "psm")
        if view in written
    }
    return {view: rows - dropped.get(view, 0) for view, rows in written.items()}


def _append_migration_step(target: Path, source_version: str) -> None:
    """Record the migration as the last provenance step."""
    from qpx._version import __version__
    from qpx.writers.provenance import ProvenanceWriter

    rows = pq.read_table(target).to_pylist()
    step = {
        "step_order": len(rows) + 1,
        "step_category": "format_conversion",
        "step_name": "qpxc migrate",
        "tool_name": "qpx",
        "tool_version": __version__,
        "tool_uri": None,
        "parameters": [{"key": "source_qpx_version", "value": source_version}],
        "config": None,
        "output_views": None,
    }
    with ProvenanceWriter(target) as writer:
        writer.write_batch([*rows, step])


def migrate_dataset(source: Path, output: Path) -> dict[str, int]:
    """Write *source* under the current specification to the new folder *output*.

    Returns the rows written per view. Views keep their file names. Ids are derived
    again, except for mzIdentML psm_id values that cannot be reproduced from the
    stored columns. The columns that hold ids, such as `psm.feature_id` and
    `site.pg_id`, follow them, and labels, peptidoforms, sequences, modification
    names and UniProt accessions take their canonical spelling. DIA-NN feature ids use the
    current composite and library-decoy flags, with the old caret prefix removed
    from their peptidoforms; files before specification 1.3 lose their theoretical
    `observed_mz`. Rows that repeat an id keep one copy, and a feature holding a
    label twice is dropped. Sample and run keep their columns outside the schema.
    MSNet psm files take the columns scripts/convert_msnet.py writes, also those it
    never converted, and a null is_decoy there becomes False.
    psm run names an old converter spelled apart from run.parquet take its spelling
    when each matches one run.
    A fixed modification the old file left out cannot be recovered: re-convert
    from the source files to get one.
    """
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    names = Dataset.structure_names()
    first, last = ("feature", "psm"), ("ontology",)
    order = [*first, *[name for name in names if name not in (*first, *last)], *last]
    files = {view: sorted(source.glob(f"*{Dataset.writer_for(view)[1]}")) for view in order}
    # Only the ids some view of this dataset refers to are kept, old to new.
    tracked = {referenced for view in order if files[view] for referenced in _ID_REFERENCES.get(view, {}).values()}
    written: dict[str, int] = {}
    ids: dict[str, dict[int, int]] = {}
    migrated_files = set()
    source_version = "unknown"
    for view in order:
        for path in files[view]:
            source_version = _footer(path).get(b"qpx_version", source_version.encode()).decode()
            target = output / path.name
            written[view] = written.get(view, 0) + _migrate_view(view, path, target, ids, track=view in tracked)
            migrated_files.add(path)
    written = _without_repeated_ids(output, files, written)
    for path in sorted(source.iterdir()):
        if path in migrated_files or not path.is_file():
            continue
        if path.suffix in _NOT_MIGRATED_SUFFIXES:
            logger.warning("%s is not migrated: build it again from the migrated views", path.name)
        else:
            shutil.copy2(path, output / path.name)
    for path in sorted(output.glob(f"*{Dataset.writer_for('provenance')[1]}")):
        _append_migration_step(path, source_version)
    logger.info("Migrated %s from QPX %s to %s", source, source_version, QPX_SPEC_VERSION)
    return written
