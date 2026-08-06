"""OpenMS QPX converter — enriches ProteomicsLFQ ``-out_qpx`` output.

OpenMS ``ProteomicsLFQ -out_qpx`` produces the core ``psm.parquet``,
``feature.parquet``, and ``pg.parquet`` tables. This converter upgrades them to
the current QPX schema, enriches them, and generates the missing metadata tables
(``run``, ``sample``, ``ontology``, ``provenance``, ``dataset``) from an SDRF.
"""

from __future__ import annotations

import logging
import tempfile
from functools import partial
from pathlib import Path
from typing import Optional

import pyarrow as pa
import pyarrow.parquet as pq

from qpx._version import __version__
from qpx.converters.channel_labels import (
    channel_labels_from_consensusxml,
    experiment_type_from_labels,
    fraction_groups_from_consensusxml,
    parse_consensusxml_maplist,
    read_sdrf_labels,
    relabel_intensities_table,
    resolve_channel_labels,
)
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.converters.sdrf import SdrfConverter
from qpx.core.constants import FEATURE, ONTOLOGY, PG, PSM, RUN, SAMPLE
from qpx.core.data.loader import load_schema
from qpx.core.scores import score_ontology_entries
from qpx.writers import FeatureWriter, PgWriter, PsmWriter

logger = logging.getLogger(__name__)

# Mapping from QPX view name to the schema loader name
_VIEW_SCHEMAS = {
    PSM: "psm",
    FEATURE: "feature",
    PG: "pg",
}

_OPENMS_IDENTITY_COMPOSITES = {
    FEATURE: ("run_file_name", "peptidoform", "charge", "rt", "scan", "observed_mz"),
    PSM: ("run_file_name", "scan", "peptidoform", "charge"),
}

_VIEW_WRITERS = {
    FEATURE: FeatureWriter,
    PSM: PsmWriter,
    PG: PgWriter,
}

# There is no PSI-MS CV term for OpenMS's experimental-design ``fraction_group``
# grouping key, so it is captured under a local quantms CV accession placeholder
# (``quantms:fraction_group``) until a PSI-MS term is minted (bigbio/qpx#221,
# OpenMS#9817). The qpx cv_param struct carries only (cv_name, cv_value) and has no
# accession slot, so the value is recorded via ``cv_name`` and round-trips
# unambiguously under that name.
FRACTION_GROUP_CV_NAME = "fraction_group"

# Run-file extensions stripped when matching consensusXML ``<map name>`` values
# against pg ``grouped_runs`` / feature ``run_file_name`` (which -out_qpx may
# write with or without an extension).
_RUN_EXTENSIONS = (".mzml", ".raw", ".mzxml", ".d", ".wiff", ".dia")


def _run_key_candidates(name: str | None) -> list[str]:
    """Return match candidates for a run/map filename (raw, basename, stem)."""
    if not name:
        return []
    text = str(name).strip()
    if not text:
        return []
    candidates = [text]
    base = text.rsplit("/", 1)[-1].rsplit("\\", 1)[-1]
    if base and base not in candidates:
        candidates.append(base)
    lower = base.lower()
    for ext in _RUN_EXTENSIONS:
        if lower.endswith(ext):
            stem = base[: -len(ext)]
            if stem and stem not in candidates:
                candidates.append(stem)
            break
    return candidates


def _fraction_group_lookup(fraction_groups: dict[str, dict[str, str]]) -> dict[str, str]:
    """Flatten ``{map_name -> design}`` into ``{run_key -> fraction_group}``.

    Each map name contributes several keys (raw, basename, extension-stripped)
    so a ``grouped_runs`` / ``run_file_name`` value matches regardless of whether
    -out_qpx kept the run-file extension.
    """
    lookup: dict[str, str] = {}
    for name, design in fraction_groups.items():
        fraction_group = design.get("fraction_group")
        if not fraction_group:
            continue
        for key in _run_key_candidates(name):
            lookup.setdefault(key, fraction_group)
    return lookup


def _fraction_group_for(run_names, lookup: dict[str, str]) -> str | None:
    """First fraction_group matching any of ``run_names`` (str or iterable)."""
    if run_names is None:
        return None
    names = [run_names] if isinstance(run_names, str) else list(run_names)
    for name in names:
        for key in _run_key_candidates(name):
            if key in lookup:
                return lookup[key]
    return None


def _discover_parquet(qpx_dir: Path) -> dict[str, Path]:
    """Discover QPX parquet files in a directory.

    Looks for files matching ``*.psm.parquet``, ``*.feature.parquet``,
    and ``*.pg.parquet``.  Returns a dict mapping view name to path.
    """
    found: dict[str, Path] = {}
    for view in (PSM, FEATURE, PG):
        matches = sorted(qpx_dir.glob(f"*.{view}.parquet"))
        if matches:
            found[view] = matches[0]
            if len(matches) > 1:
                logger.warning(
                    "Multiple %s parquet files found in %s; using %s",
                    view,
                    qpx_dir,
                    matches[0].name,
                )
    return found


def _collect_score_names(table_path: Path) -> set[str]:
    """Read ``additional_scores`` from a parquet file and collect score names."""
    table = pq.read_table(str(table_path), columns=["additional_scores"])
    names: set[str] = set()
    col = table.column("additional_scores")
    for row in col.to_pylist():
        if row:
            for score in row:
                if score and "name" in score and score["name"]:
                    names.add(score["name"])
    return names


def _validate_core(discovered: dict[str, Path]) -> None:
    """Validate each discovered parquet file against its QPX schema.

    Uses ``strict=False`` so the convert path persists source data as-produced:
    a duplicate primary key (and a null in a non-nullable column) is a warning,
    not a blocking error, while the format stabilises. Missing columns and type
    mismatches remain errors regardless of ``strict``. Null primary keys are
    still fatal — the writer rejects them at close time before this runs. The
    ``qpxc validate --strict`` audit/CI path stays strict and is unaffected.
    """
    for view, path in discovered.items():
        schema = load_schema(_VIEW_SCHEMAS[view])
        table = pq.read_table(str(path))
        result = schema.validate_full(table, strict=False)
        if not result.is_valid:
            errors = "; ".join(i.message for i in result.errors)
            raise ValueError(f"Validation failed for {path.name}: {errors}")
        logger.info(
            "%s: %s (%d rows, %d errors, %d warnings)",
            view,
            path.name,
            table.num_rows,
            len(result.errors),
            len(result.warnings),
        )


# Column on each view that identifies the row's run(s) for fraction_group
# lookup: pg groups several runs (a list), feature names a single run (scalar).
_FRACTION_GROUP_RUN_COLUMN = {PG: "grouped_runs", FEATURE: "run_file_name"}


def _source_identity_composite(parquet: pq.ParquetFile, view: str) -> tuple[str, ...] | None:
    """Return a declared source composite, or the OpenMS composite for legacy data."""
    metadata = parquet.schema_arrow.metadata or {}
    declared = metadata.get(b"identity_composite")
    if declared:
        return tuple(declared.decode().split(","))
    return _OPENMS_IDENTITY_COMPOSITES.get(view)


def _add_legacy_grouped_runs(table: pa.Table) -> pa.Table:
    """Map the pre-1.1 scalar pg run column to grouped_runs."""
    if "grouped_runs" in table.column_names or "run_file_name" not in table.column_names:
        return table
    values = [[value] if value is not None else [] for value in table.column("run_file_name").to_pylist()]
    return table.append_column("grouped_runs", pa.array(values, type=pa.list_(pa.string())))


def _rewrite_core_file(
    view: str,
    src_path: Path,
    dst: Path,
    channel_labels: dict[int, str],
    is_lfq: bool | None,
    compression: str,
    fraction_group_lookup: dict[str, str],
) -> tuple[Path, int, int]:
    """Upgrade one OpenMS core file while streaming by Parquet row group."""
    parquet = pq.ParquetFile(src_path)
    metadata = parquet.schema_arrow.metadata or {}
    writer_class = _VIEW_WRITERS[view]
    writer_kwargs = {
        "creator": metadata.get(b"creator", b"OpenMS").decode(),
        "software_provider": metadata.get(b"software_provider", b"OpenMS").decode(),
        "compression": compression,
    }
    identity_composite = _source_identity_composite(parquet, view)
    if identity_composite:
        writer_kwargs["identity_composite"] = identity_composite
    if view == PSM and b"scan_format" in metadata:
        writer_kwargs["scan_format"] = metadata[b"scan_format"].decode()

    with tempfile.NamedTemporaryFile(prefix=f".{dst.name}.", suffix=".tmp", dir=dst.parent, delete=False) as temp_file:
        temp_path = Path(temp_file.name)
    annotated = 0
    rows = 0
    fraction_group_resolver = partial(_fraction_group_for, lookup=fraction_group_lookup)
    try:
        with writer_class(temp_path, **writer_kwargs) as writer:
            for group in range(parquet.num_row_groups):
                table = parquet.read_row_group(group)
                if view == PG:
                    table = _add_legacy_grouped_runs(table)
                run_column = _FRACTION_GROUP_RUN_COLUMN.get(view)
                table, group_annotated = relabel_intensities_table(
                    table,
                    channel_labels,
                    bool(is_lfq),
                    relabel=view in (FEATURE, PG) and is_lfq is not None,
                    cv_param_name=FRACTION_GROUP_CV_NAME if fraction_group_lookup else None,
                    run_column=run_column,
                    cv_param_resolver=fraction_group_resolver,
                )
                annotated += group_annotated
                rows += table.num_rows
                if view != PG or "intensities" not in table.column_names:
                    table = writer.align_table_to_schema(table)
                writer.write_table(table)
        parquet.close()
    except Exception:
        parquet.close()
        temp_path.unlink(missing_ok=True)
        raise
    return temp_path, rows, annotated


def _copy_core(
    discovered: dict[str, Path],
    output_folder: Path,
    output_prefix: str,
    channel_labels: Optional[dict[int, str]] = None,
    is_lfq: bool | None = None,
    compression: str = "zstd",
    fraction_group_lookup: Optional[dict[str, str]] = None,
) -> dict[str, Path]:
    """
    Upgrade and copy core parquet files to the output directory.

    ``feature`` and ``pg`` carry ``intensities[].label`` — OpenMS ``-out_qpx``
    writes the run filename (feature) or a bare channel index (pg) there, so
    those two are relabeled with canonical channel labels when the experiment
    type is known (``is_lfq=None`` preserves their source labels: no SDRF
    evidence). When a consensusXML ``fraction_group`` design is available, that
    cv_param is stamped onto pg/feature rows in the same streaming pass. Legacy
    pre-1.1 tables are projected onto the current schemas and receive mandatory
    IDs before they are atomically installed at the destination.
    """
    fraction_group_lookup = fraction_group_lookup or {}
    output_paths: dict[str, Path] = {}
    staged: dict[str, tuple[Path, Path, Path, int, int]] = {}
    try:
        for view, src_path in discovered.items():
            dst = output_folder / f"{output_prefix}.{view}.parquet"
            temp_path, rows, annotated = _rewrite_core_file(
                view,
                src_path,
                dst,
                channel_labels or {},
                is_lfq,
                compression,
                fraction_group_lookup,
            )
            staged[view] = (temp_path, src_path, dst, rows, annotated)

        _validate_core({view: values[0] for view, values in staged.items()})
        for view, (temp_path, src_path, dst, rows, annotated) in staged.items():
            temp_path.replace(dst)
            output_paths[view] = dst
            logger.info(
                "Rewrote %s -> %s (%d source rows; stamped fraction_group on %d row(s))",
                src_path.name,
                dst.name,
                rows,
                annotated,
            )
    except Exception:
        for temp_path, *_ in staged.values():
            temp_path.unlink(missing_ok=True)
        raise
    return output_paths


class OpenMSConverter(BaseOrchestrator):
    """Orchestrate OpenMS QPX enrichment to a full QPX dataset.

    Usage::

        converter = OpenMSConverter(
            qpx_dir="path/to/openms_qpx_output",
            sdrf_path="metadata.sdrf.tsv",
        )
        converter.convert(
            output_folder="./qpx_full",
            output_prefix="openms",
        )
    """

    def __init__(
        self,
        qpx_dir: str | Path,
        sdrf_path: str | Path | None = None,
        consensusxml_path: str | Path | None = None,
        compression: str = "zstd",
    ):
        """Initialize the OpenMS QPX converter.

        Parameters
        ----------
        qpx_dir : str or Path
            Directory containing OpenMS ``-out_qpx`` parquet files.
        sdrf_path : str, Path or None
            Optional SDRF metadata file for sample/run generation.
        consensusxml_path : str, Path or None
            Optional OpenMS ``.consensusXML`` (the ``-out_cxml`` companion of
            ``-out_qpx``). When given, its ColumnHeaders provide the
            authoritative channel count/order for relabeling isobaric channels;
            otherwise the plex is resolved from the SDRF + data indices.
        compression : str
            Parquet compression codec (default ``zstd``).

        """
        self.qpx_dir = Path(qpx_dir)
        self.sdrf_path = str(sdrf_path) if sdrf_path else None
        self.consensusxml_path = str(consensusxml_path) if consensusxml_path else None
        self._compression = compression

    def convert(
        self,
        output_folder: str | Path,
        output_prefix: str = "openms",
        project_accession: str | None = None,
    ) -> None:
        """Run the full enrichment pipeline.

        1. Discover core QPX parquet files
        2. Upgrade, enrich, and validate core files in the output folder
        3. Convert SDRF to run.parquet + sample.parquet
        4. Collect score names for ontology
        5. Write ontology.parquet, provenance.parquet, dataset.parquet
        """
        logger.warning(
            "`qpxc convert openms` (over the OpenMS -out_qpx parquet folder) is DEPRECATED. "
            "OpenMS -out_qpx mis-assigns every PSM's run_file_name to the first run (OpenMS#9872) "
            "and emits duplicate PSMs (OpenMS#9871). Use `qpxc convert openms-consensus` "
            "(reads the consensusXML directly, resolving the run per PSM) instead. This path will "
            "be reconsidered once OpenMS ships an -out_qpx that carries the correct per-PSM run."
        )

        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)

        discovered = self.discover_and_validate()

        # Resolve canonical channel labels so the OpenMS -out_qpx run-filename /
        # bare-index labels become TMT126.. / LFQ, consistent with the
        # consensusXML and DIA-NN (quantmsdiann) QPX paths. Prefer the ColumnHeaders
        # (authoritative channel count/order) when available; otherwise resolve
        # from the SDRF-declared plex.
        sdrf_labels = read_sdrf_labels(self.sdrf_path)
        experiment_type = experiment_type_from_labels(sdrf_labels) if sdrf_labels else None

        # Parse the consensusXML leading <mapList> ONCE and derive both the
        # channel labels and the experimental-design fraction_group grouping from
        # that single pass (see parse_consensusxml_maplist) — the file may be tens
        # of GB, so it must not be read twice.
        maplist = parse_consensusxml_maplist(self.consensusxml_path) if self.consensusxml_path else {}

        channel_labels = {}
        if maplist and experiment_type:
            channel_labels = channel_labels_from_consensusxml(
                self.consensusxml_path, experiment_type, sdrf_labels, maplist=maplist
            )
            if channel_labels:
                logger.info("Resolved %d channels from consensusXML", len(channel_labels))
        if not channel_labels and experiment_type:
            channel_labels = resolve_channel_labels(experiment_type, sdrf_labels)
        is_lfq = experiment_type == "LFQ" if experiment_type else None
        if experiment_type is None:
            logger.info("No SDRF channel labels available; preserving OpenMS intensity labels")

        # OpenMS's experimental-design ``fraction_group`` (the replicate/fraction
        # grouping key) is stamped as a cv_param on pg + feature rows during the
        # copy/relabel rewrite below (same streaming pass). Best-effort: empty for
        # isobaric / pre-design consensusXML files (no fraction_group UserParams)
        # and when the file is absent. See bigbio/qpx#221, OpenMS#9817.
        fraction_groups = fraction_groups_from_consensusxml(self.consensusxml_path, maplist=maplist) if maplist else {}
        fraction_group_lookup = _fraction_group_lookup(fraction_groups) if fraction_groups else {}
        if fraction_group_lookup:
            logger.info("Resolved fraction_group for %d run(s) from consensusXML", len(fraction_groups))

        output_paths = _copy_core(
            discovered,
            output_folder,
            output_prefix,
            channel_labels,
            is_lfq,
            self._compression,
            fraction_group_lookup,
        )

        ontology_entries = self._convert_sdrf(output_folder, output_prefix)
        ontology_entries.extend(self._collect_ontology(output_paths))

        self._write_ontology(output_folder, output_prefix, ontology_entries)
        structures = list(discovered.keys())
        provenance_records = self._build_provenance(structures, experiment_type)
        self._write_provenance(output_folder, output_prefix, provenance_records)
        self._write_dataset(
            output_folder,
            output_prefix,
            project_accession,
            software_name=provenance_records[0]["tool_name"],
            software_version=None,
            provenance_records=provenance_records,
        )
        self._write_mudata(output_folder, output_prefix)
        logger.info("OpenMS QPX enrichment complete -> %s", output_folder)

    def discover_and_validate(self) -> dict[str, Path]:
        """Discover core files; validation follows their schema upgrade."""
        discovered = _discover_parquet(self.qpx_dir)
        if not discovered:
            raise FileNotFoundError(
                f"No QPX parquet files (*.psm.parquet, *.feature.parquet, *.pg.parquet) found in {self.qpx_dir}"
            )
        logger.info(
            "Discovered %d core QPX file(s): %s",
            len(discovered),
            ", ".join(f"{v}={p.name}" for v, p in discovered.items()),
        )
        return discovered

    def _convert_sdrf(
        self,
        output_folder: Path,
        output_prefix: str,
    ) -> list[dict]:
        """Convert SDRF to run + sample, returning ontology entries."""
        if not self.sdrf_path:
            logger.warning("No SDRF path provided — skipping sample/run conversion")
            return []
        with SdrfConverter(compression=self._compression) as sdrf_conv:
            sdrf_conv.convert(
                sdrf_path=self.sdrf_path,
                sample_output=str(output_folder / f"{output_prefix}.sample.parquet"),
                run_output=str(output_folder / f"{output_prefix}.run.parquet"),
            )
            entries = list(sdrf_conv.run_ontology_entries())
        logger.info("SDRF conversion complete (sample + run)")
        return entries

    @staticmethod
    def _collect_ontology(output_paths: dict[str, Path]) -> list[dict]:
        """Collect score-based ontology entries from core tables."""
        entries: list[dict] = []
        all_scores: set[str] = set()
        for view, path in output_paths.items():
            try:
                scores = _collect_score_names(path)
            except (KeyError, pa.ArrowInvalid):
                logger.debug("No additional_scores in %s", path.name)
                continue
            if scores:
                entries.extend(score_ontology_entries(scores, view=view))
                all_scores |= scores
        if all_scores:
            logger.info("Discovered scores: %s", sorted(all_scores))
        return entries

    @staticmethod
    def _build_provenance(structures: list[str], experiment_type: str | None = None) -> list[dict]:
        """Build provenance records for OpenMS + QPX conversion."""
        is_isobaric = experiment_type in {"TMT", "iTRAQ"}
        step_name = "isobaric_quantification" if is_isobaric else "label_free_quantification"
        tool_name = "OpenMS/IsobaricWorkflow" if is_isobaric else "OpenMS/ProteomicsLFQ"
        return [
            {
                "step_order": 1,
                "step_category": "quantification",
                "step_name": step_name,
                "tool_name": tool_name,
                "tool_version": None,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": structures,
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "openms_qpx_enrichment",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": [SAMPLE, RUN, ONTOLOGY],
            },
        ]
