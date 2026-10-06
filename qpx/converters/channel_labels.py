"""
Shared quantification channel-label resolution for QuantMS/OpenMS QPX.

Both QPX paths for OpenMS output — the consensusXML converter
(:mod:`qpx.converters.openms_consensus`) and the enrichment of OpenMS ``-out_qpx``
(:mod:`qpx.converters.openms`) — must put the *same* canonical reporter labels
into ``intensities[].label`` so TMT/iTRAQ/LFQ results are consistent and join
to ``sample``/``run`` the same way DIA-NN (quantmsdiann) output does.

The canonical vocabulary is the shared sdrf-pipelines ``channel_map`` (single
source of truth for label spellings, incl. the TMT10 ch10 = ``TMT131`` vs
TMT11 ch10 = ``TMT131N`` distinction). The plex is resolved, most-solid first,
from the SDRF ``comment[label]`` set (ground truth) and only then from the
channel indices present in the data.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Iterable, Optional

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from defusedxml.common import DefusedXmlException
from defusedxml.ElementTree import ParseError, iterparse
from sdrf_pipelines.converters.channel_map import CHANNEL_LABELS
from sdrf_pipelines.converters.channel_map import normalize_label as _normalize_label
from sdrf_pipelines.converters.openms.utils import infer_itraqplex, infer_tmtplex

from qpx.core.files import run_file_stem
from qpx.writers.base import parquet_write_options

__all__ = [
    "read_sdrf_labels",
    "experiment_runs_from_sdrf",
    "fraction_groups_from_sdrf",
    "experiment_type_from_labels",
    "resolve_channel_labels",
    "parse_consensusxml_maplist",
    "channel_labels_from_consensusxml",
    "fraction_groups_from_consensusxml",
    "relabel_intensities_table",
    "relabel_intensities_parquet",
    "normalize_label",
]

_FAMILY_PREFIX = {"TMT": "tmt", "ITRAQ": "itraq", "iTRAQ": "itraq"}
_CANONICAL_LABELS = {label.casefold(): label for channel_labels in CHANNEL_LABELS.values() for label in channel_labels.values()}
_SILAC_LABELS = set(CHANNEL_LABELS["silac3plex"].values())
_ONTOLOGY_NAME = re.compile(r"(?:^|;)\s*NT\s*=\s*([^;]+)", re.IGNORECASE)


def normalize_label(label: str) -> str:
    """Return the canonical spelling of a known SDRF channel label."""
    text = label.strip()
    match = _ONTOLOGY_NAME.search(text)
    if match:
        text = match.group(1).strip()
    normalized = _normalize_label(text)
    return _CANONICAL_LABELS.get(normalized.casefold(), normalized)


def read_sdrf_labels(sdrf_path: Optional[str]) -> Optional[set[str]]:
    """
    Return the distinct ``comment[label]`` values declared in the SDRF.

    This is the experiment's ground-truth channel set. Returns ``None`` when no
    SDRF is available or it has no label column.
    """
    if not sdrf_path:
        return None
    try:
        df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError):
        # Unreadable, malformed, or otherwise invalid SDRF -> no ground truth;
        # fall back to count-based inference rather than aborting the conversion.
        return None
    col = next((c for c in df.columns if c.strip().lower() == "comment[label]"), None)
    if col is None:
        return None
    labels = {str(v).strip() for v in df[col].dropna() if str(v).strip()}
    return labels or None


def experiment_runs_from_sdrf(sdrf_path: Optional[str]) -> Optional[dict[str, list[str]]]:
    """Map each SDRF ``source name`` to its raw run files (``comment[data file]`` stems).

    The SDRF is the tool-agnostic ground truth for which raw files form one
    quantification unit: rows sharing a ``source name`` are the fractions of one
    sample. This is the fallback source for pg ``grouped_runs`` when a converter's
    own experiment->runs result mapping (FragPipe ``experiment_annotation.tsv`` /
    MaxQuant ``evidence.txt``) is absent. Files are returned distinct and in SDRF
    row order (fraction order preserved, not sorted). Returns ``None`` when no
    SDRF is available or it lacks ``source name`` / ``comment[data file]``.
    """
    if not sdrf_path:
        return None
    try:
        df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError):
        return None
    cols = {c.strip().lower(): c for c in df.columns}
    source_col = cols.get("source name")
    file_col = cols.get("comment[data file]")
    if source_col is None or file_col is None:
        return None
    mapping: dict[str, list[str]] = {}
    for source, file_value in zip(df[source_col], df[file_col]):
        if pd.isna(source) or pd.isna(file_value):
            continue
        experiment = str(source).strip()
        run = run_file_stem(str(file_value))
        if not experiment or not run:
            continue
        runs = mapping.setdefault(experiment, [])
        if run not in runs:
            runs.append(run)
    return mapping or None


def _sdrf_unit_columns(df) -> Optional[dict[str, Optional[str]]]:
    """Resolve the SDRF columns that define a quantification unit, or None."""
    cols = {c.strip().lower(): c for c in df.columns}
    file_col = cols.get("comment[data file]")
    source_col = cols.get("source name")
    if file_col is None or source_col is None:
        return None
    return {
        "file": file_col,
        "source": source_col,
        "label": cols.get("comment[label]"),
        "fraction": cols.get("comment[fraction identifier]"),
        "brep": next((cols[k] for k in cols if "biological replicate" in k), None),
        "trep": next((cols[k] for k in cols if "technical replicate" in k), None),
    }


def _safe_fraction(value: object) -> int:
    try:
        return int(str(value).strip()) if value is not None else 0
    except (ValueError, TypeError):
        return 0


def _collect_run_signatures(df, c: dict[str, Optional[str]]):
    """Per run stem: its ``(source,label,brep,trep)`` signature set + fraction number."""
    signature: dict[str, set] = defaultdict(set)
    fraction: dict[str, int] = {}
    first_seen: list[str] = []
    for row in df.to_dict("records"):
        raw = row.get(c["file"])
        if raw is None or not str(raw).strip():
            continue
        run = run_file_stem(str(raw))
        if run not in signature:
            first_seen.append(run)
        signature[run].add(tuple(str(row.get(c[k], "")).strip() if c[k] else "" for k in ("source", "label", "brep", "trep")))
        if run not in fraction:
            fraction[run] = _safe_fraction(row.get(c["fraction"])) if c["fraction"] else 0
    return signature, fraction, first_seen


def fraction_groups_from_sdrf(sdrf_path: Optional[str]) -> Optional[dict[str, list[str]]]:
    """Map each SDRF raw file to the ``grouped_runs`` of its quantification unit.

    A quantification unit is the set of raw files that differ **only in fraction**:
    they carry the same set of ``(source name, label, biological replicate,
    technical replicate)`` assignments. For LFQ, a sample's fractions group
    together; for TMT/iTRAQ, all fraction files of one labelled mixture group
    together (they share the same channel→sample assignments). Technical
    replicates (which differ in ``technical replicate``) stay in separate units.

    Returns a dict mapping each ``run_file_name`` stem to its ``grouped_runs``
    list (the unit's files, distinct, ordered by fraction identifier), or ``None``
    when no SDRF or the required columns are absent.
    """
    if not sdrf_path:
        return None
    try:
        df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError):
        return None
    cols = _sdrf_unit_columns(df)
    if cols is None:
        return None
    signature, fraction, first_seen = _collect_run_signatures(df, cols)
    if not signature:
        return None

    groups: dict[frozenset, list[str]] = defaultdict(list)
    for run in first_seen:
        groups[frozenset(signature[run])].append(run)

    run_to_grouped: dict[str, list[str]] = {}
    for members in groups.values():
        ordered: list[str] = []
        for run in sorted(members, key=lambda r: (fraction.get(r, 0), first_seen.index(r))):
            if run not in ordered:
                ordered.append(run)
        for run in ordered:
            run_to_grouped[run] = ordered
    return run_to_grouped or None


def experiment_type_from_labels(sdrf_labels: Optional[set[str]]) -> str:
    """Classify TMT, iTRAQ, SILAC or LFQ from the SDRF channel labels."""
    if not sdrf_labels:
        return "LFQ"
    upper = {label.upper() for label in sdrf_labels}
    if any("TMT" in label for label in upper):
        return "TMT"
    if any("ITRAQ" in label for label in upper):
        return "iTRAQ"
    if any(normalize_label(label) in _SILAC_LABELS for label in sdrf_labels):
        return "SILAC"
    return "LFQ"


def resolve_channel_labels(
    experiment_type: str,
    sdrf_labels: Optional[set[str]] = None,
    channel_indices: Optional[Iterable] = None,
) -> dict[int, str]:
    """
    Resolve the ``{1-based channel index -> canonical label}`` map for a run.

    Resolution order (most-solid first):

    1. From the SDRF ``comment[label]`` set via sdrf-pipelines
       ``infer_tmtplex`` / ``infer_itraqplex`` — robust even when a channel is
       empty in the data.
    2. Fallback: the smallest plex of the reagent family whose channel count
       covers the highest index present in ``channel_indices``.

    Returns ``{}`` for non-isobaric or unresolvable input (caller falls back to
    the raw value / ``"LFQ"``).
    """
    prefix = _FAMILY_PREFIX.get(experiment_type)
    if prefix is None:
        return {}

    # 1. SDRF-declared plex (ground truth).
    if sdrf_labels:
        canonical = {normalize_label(label) for label in sdrf_labels if label}
        infer = infer_tmtplex if prefix == "tmt" else infer_itraqplex
        try:
            return CHANNEL_LABELS[infer(canonical)]
        except (ValueError, KeyError):
            pass

    # 2. Fallback: infer plex from the channel indices present in the data.
    # NB: channel_indices may be a numpy array (df[col].values), so test for
    # None explicitly — `channel_indices or ()` raises on array truthiness.
    indices: set[int] = set()
    for value in () if channel_indices is None else channel_indices:
        try:
            indices.add(int(float(value)))
        except (TypeError, ValueError):
            continue
    if not indices:
        return {}
    max_index = max(indices)

    family = {plex: labels for plex, labels in CHANNEL_LABELS.items() if plex.startswith(prefix)}
    for plex in sorted(family, key=lambda p: len(family[p])):
        if len(family[plex]) >= max_index:
            return family[plex]
    return {}


# OpenMS ConsensusMap ``<map>`` UserParams that describe the experimental-design
# grouping of each run. ``fraction_group`` is OpenMS's replicate/fraction
# grouping key: runs that share a ``fraction_group`` are fractions of the same
# quantification unit (see the OpenMS experimental design). ``fraction`` and
# ``sample_name`` are captured alongside for provenance.
_MAP_DESIGN_PARAMS = ("fraction_group", "fraction", "sample_name")


def parse_consensusxml_maplist(consensusxml_path: str) -> dict[int, dict[str, str]]:
    """Parse a consensusXML's leading ``<mapList>`` exactly once.

    Both channel-label resolution and fraction-group capture need the same
    leading ``<mapList>`` block, so this single raw pass is the shared source
    for :func:`channel_labels_from_consensusxml` and
    :func:`fraction_groups_from_consensusxml` — avoiding two iterparse passes
    over the (possibly tens-of-GB) file.

    Returns ``{map_index (0-based ``id``) -> {"label", "name", and any of
    ``fraction_group``/``fraction``/``sample_name`` that were declared}}`` in
    document order. Only the leading ``mapList`` is read (defused ``iterparse``,
    bounded memory), and ``{}`` is returned on a missing or malformed file.
    """
    maps: dict[int, dict[str, str]] = {}
    try:
        in_map_list = False
        current_id: Optional[int] = None
        current: Optional[dict[str, str]] = None
        for event, element in iterparse(consensusxml_path, events=("start", "end")):
            tag = element.tag.rsplit("}", 1)[-1]
            if event == "start" and tag == "mapList":
                in_map_list = True
            elif in_map_list and event == "start" and tag == "map":
                current_id = int(element.attrib["id"])
                current = {
                    "label": element.attrib.get("label", "").strip(),
                    "name": element.attrib.get("name", "").strip(),
                }
            elif in_map_list and event == "start" and tag == "UserParam" and current is not None:
                param_name = element.attrib.get("name", "")
                if param_name in _MAP_DESIGN_PARAMS:
                    current.update({param_name: element.attrib.get("value", "")})
            elif in_map_list and event == "end" and tag == "map":
                if current_id is not None and current is not None:
                    maps[current_id] = current
                current_id = None
                current = None
            elif event == "end" and tag == "mapList":
                break
            if event == "end":
                element.clear()
    except (OSError, ParseError, DefusedXmlException, KeyError, ValueError):
        return {}
    return maps


def channel_labels_from_consensusxml(
    consensusxml_path: str,
    experiment_type: str,
    sdrf_labels: Optional[set[str]] = None,
    maplist: Optional[dict[int, dict[str, str]]] = None,
) -> dict[int, str]:
    """
    Resolve ``{1-based channel index -> canonical label}`` from a consensusXML.

    The OpenMS ConsensusMap ``ColumnHeaders`` give the **authoritative** channel
    count and order (map-index), independent of which channels happen to carry
    data — so this is more robust than inferring the plex from the data alone.
    The canonical spelling still comes from the SDRF-declared plex + channel_map;
    a ColumnHeader ``label`` is used only as a fallback for a position the plex
    map does not cover. Only the leading ``mapList`` is parsed, so even very
    large consensusXML files are handled with bounded memory. Returns ``{}``
    when the file is unavailable or malformed.

    Pass a pre-parsed *maplist* (from :func:`parse_consensusxml_maplist`) to
    reuse a single parse when the caller also needs the fraction groups.
    """
    maps = maplist if maplist is not None else parse_consensusxml_maplist(consensusxml_path)
    headers = {map_index: entry.get("label", "") for map_index, entry in maps.items()}
    if not headers:
        return {}

    n_channels = len(headers)
    # Plex resolved with the consensusXML channel count as the index range.
    plex = resolve_channel_labels(experiment_type, sdrf_labels, range(1, n_channels + 1))
    labels: dict[int, str] = {}
    for map_index in sorted(headers):
        position = map_index + 1  # ColumnHeaders are 0-based; qpx channels are 1-based
        raw = headers[map_index]
        labels[position] = plex.get(position) or (normalize_label(raw) if raw else str(position))
    return labels


def fraction_groups_from_consensusxml(
    consensusxml_path: str,
    maplist: Optional[dict[int, dict[str, str]]] = None,
) -> dict[str, dict[str, str]]:
    """
    Extract the per-map experimental-design grouping from a consensusXML.

    OpenMS writes ``<UserParam name="fraction_group">``,
    ``<UserParam name="fraction">`` and ``<UserParam name="sample_name">`` into
    every ``<mapList>`` ``<map>`` entry for LFQ (label-free) runs. This is the
    grouping key that makes ``grouped_runs`` correct for fractionated LFQ, yet it
    is otherwise dropped on the OpenMS ``-out_qpx`` path.

    Returns ``{map_name -> {"fraction_group", "fraction", "sample_name"}}`` keyed
    by the map ``name`` attribute (the run file name; falls back to the 0-based
    ``id`` when a map has no name). Values are the raw string values from the
    XML. Only maps that actually declare a ``fraction_group`` are returned, so
    isobaric / pre-design consensusXML files (no such UserParams) yield ``{}``.
    Like :func:`channel_labels_from_consensusxml`, only the leading ``mapList``
    is parsed with a bounded, defused-XML ``iterparse`` — even very large files
    are handled with constant memory, and ``{}`` is returned on a missing or
    malformed file.

    Pass a pre-parsed *maplist* (from :func:`parse_consensusxml_maplist`) to
    reuse a single parse when the caller also needs the channel labels.
    """
    maps = maplist if maplist is not None else parse_consensusxml_maplist(consensusxml_path)
    groups: dict[str, dict[str, str]] = {}
    for map_index in sorted(maps):
        entry = maps[map_index]
        design = {name: entry[name] for name in _MAP_DESIGN_PARAMS if name in entry}
        if not design.get("fraction_group"):
            continue
        key = entry.get("name") or str(map_index)
        groups[key] = design
    return groups


def _relabel_entries(rows, channel_labels: dict[int, str], is_lfq: bool):
    """
    Relabel the ``label`` field of a list<struct> intensities column.

    OpenMS ``-out_qpx`` writes the *filename* into feature ``intensities[].label``
    and a bare 1-based index (``"1"``, ``"2"``, ...) into pg ``intensities[].label``.
    Recover the channel index from a numeric label, else from the array position
    (channels are emitted in index order), then map to the canonical label.
    LFQ collapses to the single ``"LFQ"`` label.
    """
    out = []
    for row in rows:
        if row is None:
            out.append(None)
            continue
        new_row = []
        for position, entry in enumerate(row):
            label = entry.get("label")
            if is_lfq:
                new_label = "LFQ"
            else:
                text = str(label) if label is not None else ""
                index = int(text) if text.isdigit() else position + 1
                new_label = channel_labels.get(index, normalize_label(text) if text else text)
            updated = dict(entry)
            updated["label"] = new_label
            new_row.append(updated)
        out.append(new_row)
    return out


def _relabel_scalar_labels(values, channel_labels: dict[int, str], is_lfq: bool):
    """Relabel a flat scalar ``label`` column (pg is flattened since QPX 1.1).

    OpenMS writes a bare 1-based index (``"1"``, ``"2"``, ...) into the pg
    ``label`` column; map the numeric index to its canonical channel label. LFQ
    collapses to ``"LFQ"``. Null labels (identification-only rows) pass through.
    Non-numeric labels are normalized but otherwise kept.
    """
    out = []
    for value in values:
        if value is None:
            out.append(None)
            continue
        if is_lfq:
            out.append("LFQ")
            continue
        text = str(value)
        if text.isdigit():
            out.append(channel_labels.get(int(text), normalize_label(text)))
        else:
            out.append(normalize_label(text))
    return out


def _resolve_parquet_compression(
    parquet: pq.ParquetFile,
    source_metadata: dict[bytes, bytes],
    compression: str | None,
) -> str:
    """Resolve the output codec from an override, QPX metadata, or the file."""
    if compression is not None:
        return compression

    declared = source_metadata.get(b"compression_format")
    if declared:
        return declared.decode()

    if parquet.num_row_groups:
        row_group = parquet.metadata.row_group(0)
        codecs = {row_group.column(index).compression.lower() for index in range(row_group.num_columns)}
        if len(codecs) == 1:
            codec = codecs.pop()
            return "none" if codec == "uncompressed" else codec
    return "snappy"


def _append_cv_param_column(
    table,
    run_column: str,
    cv_param_name: str,
    resolver,
) -> tuple[object, int]:
    """Append a ``{cv_name, cv_value}`` param to each row whose run resolves.

    ``resolver(run_value)`` returns the cv_value string for a row's
    ``run_column`` value (or ``None`` to leave the row untouched); existing
    ``cv_params`` are preserved and the param is not duplicated. Returns the
    (possibly updated) table and the number of rows annotated.
    """
    if run_column not in table.column_names or "cv_params" not in table.column_names:
        return table, 0
    run_values = table.column(run_column).to_pylist()
    cv_values = table.column("cv_params").to_pylist()
    new_cv = []
    annotated = 0
    for run_value, cv_params in zip(run_values, cv_values):
        cv_value = resolver(run_value)
        params = list(cv_params) if cv_params else []
        if cv_value is not None and not any((p or {}).get("cv_name") == cv_param_name for p in params):
            params.append({"cv_name": cv_param_name, "cv_value": str(cv_value)})
            annotated += 1
        new_cv.append(params or None)
    field_index = table.schema.get_field_index("cv_params")
    cv_array = pa.array(new_cv, type=table.column("cv_params").type)
    return table.set_column(field_index, "cv_params", cv_array), annotated


def relabel_intensities_table(
    table: pa.Table,
    channel_labels: dict[int, str],
    is_lfq: bool,
    columns: tuple[str, ...] = ("intensities", "label", "additional_intensities"),
    *,
    relabel: bool = True,
    cv_param_name: str | None = None,
    run_column: str | None = None,
    cv_param_resolver=None,
) -> tuple[pa.Table, int]:
    """Relabel one Arrow table and return it with the annotated-row count."""
    if relabel:
        for column in columns:
            if column not in table.column_names:
                continue
            field_index = table.schema.get_field_index(column)
            original = table.column(column)
            if pa.types.is_list(original.type) or pa.types.is_large_list(original.type):
                values = _relabel_entries(original.to_pylist(), channel_labels, is_lfq)
            elif column == "label" and pa.types.is_string(original.type):
                values = _relabel_scalar_labels(original.to_pylist(), channel_labels, is_lfq)
            else:
                continue
            table = table.set_column(field_index, column, pa.array(values, type=original.type))

    do_cv_param = cv_param_resolver is not None and run_column is not None and cv_param_name is not None
    if not do_cv_param:
        return table, 0
    return _append_cv_param_column(table, run_column, cv_param_name, cv_param_resolver)


def relabel_intensities_parquet(
    src_path: str,
    dst_path: str,
    channel_labels: dict[int, str],
    is_lfq: bool,
    columns: tuple[str, ...] = ("intensities", "label", "additional_intensities"),
    compression: str | None = None,
    *,
    relabel: bool = True,
    cv_param_name: str | None = None,
    run_column: str | None = None,
    cv_param_resolver=None,
) -> int:
    """
    Rewrite ``src_path`` to ``dst_path`` in a single streaming pass.

    Two transforms are applied together so each file is rewritten only once:

    * **relabel** (default) — canonicalize ``intensities[].label`` /
      ``additional_intensities[].label`` (see :func:`_relabel_entries`). Set
      ``relabel=False`` to leave intensity labels untouched (e.g. a cv_param-only
      pass when no channel evidence is available).
    * **cv_param annotation** — when ``cv_param_resolver``, ``run_column`` and
      ``cv_param_name`` are all supplied, append a ``{cv_name, cv_value}`` param
      to each row whose ``run_column`` value resolves via the resolver (used to
      stamp OpenMS's ``fraction_group`` onto pg/feature rows).

    Streams by row group so memory stays bounded on large tables. Non-target
    columns pass through untouched. QPX compression and selective
    BYTE_STREAM_SPLIT encoding are retained; *compression* overrides the source
    codec when supplied. Returns the number of rows annotated with a cv_param.
    """
    do_cv_param = cv_param_resolver is not None and run_column is not None and cv_param_name is not None

    parquet = pq.ParquetFile(src_path)
    source_metadata = dict(parquet.schema_arrow.metadata or {})
    output_compression = _resolve_parquet_compression(parquet, source_metadata, compression)

    if compression is not None or b"compression_format" in source_metadata:
        source_metadata[b"compression_format"] = output_compression.encode()
    output_schema = parquet.schema_arrow.with_metadata(source_metadata or None)
    writer = pq.ParquetWriter(
        dst_path,
        output_schema,
        **parquet_write_options(output_schema, output_compression),
    )
    annotated = 0
    try:
        for group in range(parquet.num_row_groups):
            table = parquet.read_row_group(group)
            table, group_annotated = relabel_intensities_table(
                table,
                channel_labels,
                is_lfq,
                columns,
                relabel=relabel,
                cv_param_name=cv_param_name if do_cv_param else None,
                run_column=run_column,
                cv_param_resolver=cv_param_resolver,
            )
            annotated += group_annotated
            table = table.replace_schema_metadata(output_schema.metadata)
            writer.write_table(table)
    finally:
        writer.close()
    return annotated
