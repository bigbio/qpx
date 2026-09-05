"""Orchestrate consensusXML (+ SDRF) -> a QPX dataset (feature/psm/pg [+ run/sample]).

Interim quantms path while OpenMS ``-out_qpx`` is pre-1.1. pg carries an interim
unnormalized unique-peptide-sum intensity (until ``-out_qpx`` provides the real
quant). run/sample come from the SDRF (reusing :class:`SdrfConverter`) when an
SDRF is provided; when it is, the consensusXML channels are also checked against
the SDRF ``comment[label]`` and mismatches are logged as warnings.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pyarrow as pa
import pyarrow.parquet as pq

from qpx._version import __version__
from qpx.converters.openms_consensus.feature_adapter import (
    check_channels_vs_sdrf,
    load_consensus_map,
)
from qpx.converters.openms_consensus.pg_adapter import (
    accession_to_group,
    consensus_protein_groups_to_records,
)
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.core.constants import FEATURE, ONTOLOGY, PG, PSM, RUN, SAMPLE
from qpx.core.data import FeatureSchema
from qpx.core.data.identity import derive_id
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter

_log = logging.getLogger(__name__)

_STRUCTURE_ALL = ("feature", "psm", "pg", "run", "sample")

# Producer-specific feature identity composite for the openms-consensus path
# (measured keys agreed in bigbio/qpx#229). Passed to the FeatureWriter so the
# derived feature_id hashes exactly these columns (all present in feature.yaml)
# instead of the schema default.
_FEATURE_IDENTITY_COMPOSITE = (
    "peptidoform",
    "charge",
    "run_file_name",
    "rt",
    "scan",
    "observed_mz",
    "consensus_rt",
)
# The psm view's schema identity_composite (psm.yaml) — passed to the PsmWriter
# so the derived psm_id hashes exactly these columns.
_PSM_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name", "scan")


def _derive_persisted_record_id(record: dict, composite: tuple[str, ...], schema: pa.Schema) -> int:
    """Derive an id from the values exactly as they will be stored by Arrow."""
    values = [pa.scalar(record.get(name), type=schema.field(name).type).as_py() for name in composite]
    return derive_id(values)


def _link_psm_feature(feature_records: list[dict], psm_records: list[dict]) -> None:
    """Stamp ``psm.feature_id`` in place for one consensus feature.

    ``psm.feature_id`` is the authoritative producer assignment — which consensus
    feature a PSM belongs to. The id is derived with the SAME composite the
    FeatureWriter uses, so this converter-populated foreign key matches the
    writer-derived ``feature_id`` byte-for-byte (``derive_id`` is deterministic).
    A consensus feature yields one feature record per run; each PSM links to the
    feature record of its own run. PSMs whose run has no feature record keep
    ``feature_id`` null (resolves #182).

    The inverse ``feature.psm_ids`` is NOT materialized here: it is the pure
    inverse of ``psm.feature_id`` (group PSMs by their ``feature_id``) and is
    computed on read via :meth:`qpx.dataset.Dataset.link_feature_psm`.
    """
    feature_schema = FeatureSchema.get_arrow_schema()
    feat_id_by_run: dict[str, int] = {
        rec["run_file_name"]: _derive_persisted_record_id(rec, _FEATURE_IDENTITY_COMPOSITE, feature_schema)
        for rec in feature_records
    }
    for prec in psm_records:
        feat_id = feat_id_by_run.get(prec.get("run_file_name"))
        if feat_id is not None:
            prec["feature_id"] = feat_id


# Above this consensusXML size, auto-select the streaming reader: the pyopenms
# in-memory load needs ~0.8x the file in RAM, so ~4 GB (≈3 GB RAM) is a safe point
# to switch to the low-memory path on a normal node.
_STREAM_THRESHOLD_BYTES = 4 * 1024**3


def _should_stream(path: str) -> bool:
    try:
        return Path(path).stat().st_size > _STREAM_THRESHOLD_BYTES
    except OSError:
        return False


def _cf_feature_psm_records(cf, map_info, group_map, resolve_run, seen, *, want_feature, want_psm, enzyme=None):
    """Feature + PSM records for one consensus feature, cross-linked when both views
    are emitted. Shared by the streaming and in-memory paths so their output matches.
    """
    from qpx.converters.openms_consensus.feature_adapter import feature_records_for_cf
    from qpx.converters.openms_consensus.psm_adapter import _cf_element_runs, psm_records_for_pid

    cf_feats = feature_records_for_cf(cf, map_info, group_map, enzyme=enzyme) if want_feature else []
    cf_psms: list[dict] = []
    if want_psm:
        # Multi-run isobaric PIDs carry a local id_merge_index; the feature's
        # element runs disambiguate which run they belong to.
        cf_runs = _cf_element_runs(cf, map_info)
        for pid in cf.getPeptideIdentifications():
            cf_psms.extend(psm_records_for_pid(pid, resolve_run, seen, enzyme=enzyme, cf_runs=cf_runs))
    # Stamp psm.feature_id only when both views are emitted (the FK references a
    # feature row written in this dataset). feature.psm_ids is the computed inverse.
    if want_feature and want_psm:
        _link_psm_feature(cf_feats, cf_psms)
    return cf_feats, cf_psms


def _write_view(writer_cls, path, records, *, creator, compression, identity_composite=None):
    """Write ``records`` to a view parquet via its writer (empty -> header only)."""
    kwargs = {"creator": creator, "compression": compression}
    if identity_composite is not None:
        kwargs["identity_composite"] = identity_composite
    with writer_cls(str(path), **kwargs) as w:
        if records:
            w.write_batch(records)
    return path


def _stream_feature_psm(
    cm,
    fw,
    pw,
    *,
    map_info,
    group_map,
    resolve_run,
    maps,
    pep_intensity,
    map_run,
    seen,
    batch,
    include_unassigned_psms=True,
    enzyme=None,
):
    """One ordered element/unassigned pass: write feature/psm in batches and
    accumulate the pg maps in place (the streaming path's inner loop)."""
    from qpx.converters.openms_consensus.pg_adapter import (
        accumulate_cf_intensity,
        accumulate_cf_maps,
        accumulate_unassigned_maps,
    )
    from qpx.converters.openms_consensus.psm_adapter import psm_records_for_pid

    feat_buf: list[dict] = []
    psm_buf: list[dict] = []
    for kind, obj in cm.iter_all():
        if kind == "element":
            cf_feats, cf_psms = _cf_feature_psm_records(
                obj,
                map_info,
                group_map,
                resolve_run,
                seen,
                want_feature=fw is not None,
                want_psm=pw is not None,
                enzyme=enzyme,
            )
            feat_buf.extend(cf_feats)
            psm_buf.extend(cf_psms)
            if maps is not None:
                accumulate_cf_maps(obj, map_run, maps)
                accumulate_cf_intensity(obj, map_info, pep_intensity)
        else:  # unassigned peptide identification
            if pw is not None and include_unassigned_psms:
                # Unassigned PSMs map to no feature -> feature_id stays null.
                psm_buf.extend(psm_records_for_pid(obj, resolve_run, seen, enzyme=enzyme))
            if maps is not None:
                # Protein inference always sees every identification, whether or
                # not the PSM rows are emitted: dropping evidence would change the
                # protein groups, which is not what this option is for.
                accumulate_unassigned_maps(obj, resolve_run, maps)
        if fw is not None and len(feat_buf) >= batch:
            fw.write_batch(feat_buf)
            feat_buf = []
        if pw is not None and len(psm_buf) >= batch:
            pw.write_batch(psm_buf)
            psm_buf = []
    if fw is not None and feat_buf:
        fw.write_batch(feat_buf)
    if pw is not None and psm_buf:
        pw.write_batch(psm_buf)


def _convert_streaming(
    consensusxml_path,
    out,
    output_prefix,
    structures,
    sdrf_path,
    creator,
    pg_top,
    compression="zstd",
    include_unassigned_psms=True,
) -> dict:
    """Single-pass, low-memory feature/psm/pg from a streamed consensusXML.

    One ordered ``iter_all()`` pass over the elements + unassigned IDs feeds the
    same per-element builders the in-memory adapters use (so output is identical),
    writing feature/psm in batches and accumulating the pg maps; pg records are
    built at the end. The whole map is never held, and the file is parsed once.
    """
    from collections import defaultdict
    from contextlib import ExitStack

    from qpx.converters.openms_consensus.feature_adapter import _run_stem, feature_map_info, resolve_enzyme
    from qpx.converters.openms_consensus.pg_adapter import (
        _ProteinMaps,
        accession_to_group,
        build_pg_records,
    )
    from qpx.converters.openms_consensus.psm_adapter import _run_resolver
    from qpx.converters.openms_consensus.streaming import StreamingConsensusMap

    cm = StreamingConsensusMap(consensusxml_path)
    if sdrf_path:
        for msg in check_channels_vs_sdrf(cm, sdrf_path):
            _log.warning("consensusXML/SDRF channel mismatch: %s", msg)

    map_info = feature_map_info(cm)
    headers = cm.getColumnHeaders()
    map_run = {i: _run_stem(headers[i].filename) for i in headers}
    want_feature, want_psm, want_pg = ("feature" in structures, "psm" in structures, "pg" in structures)
    group_map = accession_to_group(cm) if want_feature else None
    resolve_run = _run_resolver(cm) if want_psm or want_pg else None
    maps = _ProteinMaps() if want_pg else None
    pep_intensity: dict = defaultdict(float) if want_pg else {}
    seen: set = set()
    written: dict[str, Path] = {}

    with ExitStack() as stack:
        fw = pw = None
        if want_feature:
            fw = stack.enter_context(
                FeatureWriter(
                    str(out / f"{output_prefix}.feature.parquet"),
                    creator=creator,
                    identity_composite=_FEATURE_IDENTITY_COMPOSITE,
                    compression=compression,
                )
            )
        if want_psm:
            pw = stack.enter_context(
                PsmWriter(
                    str(out / f"{output_prefix}.psm.parquet"),
                    creator=creator,
                    identity_composite=_PSM_IDENTITY_COMPOSITE,
                    compression=compression,
                )
            )
        _stream_feature_psm(
            cm,
            fw,
            pw,
            map_info=map_info,
            group_map=group_map,
            resolve_run=resolve_run,
            maps=maps,
            pep_intensity=pep_intensity,
            map_run=map_run,
            seen=seen,
            batch=100_000,
            include_unassigned_psms=include_unassigned_psms,
            enzyme=resolve_enzyme(cm, sdrf_path),
        )
        if fw is not None:
            written["feature"] = out / f"{output_prefix}.feature.parquet"
        if pw is not None:
            written["psm"] = out / f"{output_prefix}.psm.parquet"

    if want_pg:
        records = build_pg_records(cm, map_info, maps, pep_intensity, sdrf_path, pg_top)
        written["pg"] = _write_view(
            PgWriter, out / f"{output_prefix}.pg.parquet", records, creator=creator, compression=compression
        )
    return written


def _validate_structures(structures: tuple[str, ...], sdrf_path: Optional[str]) -> None:
    """Reject unknown structure names and run/sample requested without an SDRF."""
    unknown = [s for s in structures if s not in _STRUCTURE_ALL]
    if unknown:
        raise ValueError(f"Unknown structure(s) {unknown}; valid values are {list(_STRUCTURE_ALL)}")
    metadata = set(structures).intersection({"run", "sample"})
    if metadata and not sdrf_path:
        raise ValueError(f"An SDRF is required to write {sorted(metadata)}")


def _collect_score_names(table_path: Path) -> set[str]:
    """Read ``additional_scores`` score names from a written parquet file."""
    table = pq.read_table(str(table_path), columns=["additional_scores"])
    names: set[str] = set()
    for row in table.column("additional_scores").to_pylist():
        if row:
            for score in row:
                if score and score.get("score_name"):
                    names.add(score["score_name"])
    return names


def _consensus_is_isobaric(consensusxml_path: str) -> bool:
    """Return whether a consensusXML's map labels indicate an isobaric (TMT/iTRAQ) run.

    Reads only the small column-map header via the streaming reader, so it never
    loads the whole ConsensusMap into memory — the full pyopenms load would defeat
    the streaming path and can exhaust RAM on multi-GB files.
    """
    from qpx.converters.openms_consensus.feature_adapter import consensus_channels
    from qpx.converters.openms_consensus.streaming import StreamingConsensusMap

    try:
        return bool(consensus_channels(StreamingConsensusMap(consensusxml_path)))
    except (OSError, ValueError, KeyError, AttributeError, RuntimeError) as exc:
        # Best-effort provenance hint: any read/parse failure -> assume label-free.
        _log.warning("could not determine isobaric labeling from %s: %s", consensusxml_path, exc)
        return False


def _write_sdrf_metadata(
    output_folder: Path,
    output_prefix: str,
    sdrf_path: str,
    requested: set[str],
    compression: str = "zstd",
) -> tuple[dict[str, Path], list[dict]]:
    """Write the requested SDRF-backed run/sample structures.

    Returns ``(paths, run_ontology_entries)`` — the latter collected from the
    same ``SdrfConverter`` instance (its parse state populates the run ontology),
    so callers need not re-parse the SDRF.
    """
    metadata = requested.intersection({"run", "sample"})
    if not metadata:
        return {}, []

    from qpx.converters.sdrf import SdrfConverter

    paths = {
        "sample": output_folder / f"{output_prefix}.sample.parquet",
        "run": output_folder / f"{output_prefix}.run.parquet",
    }
    with SdrfConverter(compression=compression) as sdrf_converter:
        sdrf_converter.convert(
            sdrf_path=sdrf_path,
            sample_output=str(paths["sample"]) if "sample" in metadata else None,
            run_output=str(paths["run"]) if "run" in metadata else None,
        )
        run_ontology = list(sdrf_converter.run_ontology_entries())
    return {name: paths[name] for name in metadata}, run_ontology


class OpenMSConsensusConverter(BaseOrchestrator):  # pylint: disable=too-few-public-methods
    """consensusXML + SDRF -> QPX views.

    A single-entry orchestrator (``convert``) — the interim counterpart to the
    other converter classes, kept as a class for call-site symmetry with them.
    """

    def convert(
        self,
        consensusxml_path: str,
        output_folder: str,
        output_prefix: str = "openms",
        sdrf_path: Optional[str] = None,
        structures: tuple[str, ...] = _STRUCTURE_ALL,
        creator: str = "openms-consensus",
        pg_top: int = 0,
        streaming: Optional[bool] = None,
        project_accession: Optional[str] = None,
        compression: str = "zstd",
        include_unassigned_psms: bool = True,
    ) -> dict[str, Path]:
        """Write the requested QPX views and return ``{structure: parquet path}``.

        ``include_unassigned_psms`` (default ``True``) controls whether
        unassigned PeptideIdentifications from the source reach ``psm.parquet``.
        They are retained to preserve identification evidence; pass ``False``
        to exclude them.

        ``feature_id`` records a link in the exported dataset, not quantification
        status. It is only populated when both feature and PSM views are emitted;
        PSM-only output leaves it null even for assigned identifications.
        ``PSM.with_feature()`` and ``PSM.without_feature()`` filter these links.

        The feature view always omits consensus features without peptide hits.
        Protein inference always uses every identification, regardless of this
        option.

        ``structures`` selects which of feature/psm/pg/run/sample to emit. pg
        carries an interim unnormalized unique-peptide-sum intensity; ``pg_top``
        bounds the peptides used (0 = all; 3 mirrors ProteomicsLFQ/IsobaricWorkflow).

        ``streaming`` picks the consensusXML reader: ``None`` (default) auto-selects
        — the low-memory streaming reader for files above ``_STREAM_THRESHOLD_BYTES``
        (which pyopenms would otherwise load whole into ~0.8x-file RAM), else the
        faster in-memory pyopenms load. ``True``/``False`` forces the choice.

        ``project_accession`` (e.g. ``PXD001819``) is stamped into dataset.parquet.

        ``compression`` is the Parquet codec for every written view (default zstd).
        """
        self._compression = compression
        _validate_structures(structures, sdrf_path)
        requested = set(structures)

        out = Path(output_folder)
        out.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}

        if {"feature", "psm", "pg"}.intersection(structures):
            use_stream = streaming if streaming is not None else _should_stream(consensusxml_path)
            if use_stream:
                # Low-memory path: a single ordered pass over the consensusXML builds
                # feature/psm/pg together and writes in batches, so the whole map is
                # never in memory and the file is parsed once (not once per view).
                written.update(
                    _convert_streaming(
                        consensusxml_path,
                        out,
                        output_prefix,
                        structures,
                        sdrf_path,
                        creator,
                        pg_top,
                        compression,
                        include_unassigned_psms,
                    )
                )
            else:
                # In-memory path: pyopenms loads the map once (fast for smaller files);
                # the adapters iterate it cheaply. Output is identical either way.
                written.update(
                    self._convert_in_memory(
                        consensusxml_path,
                        out,
                        output_prefix,
                        structures,
                        sdrf_path,
                        creator,
                        pg_top,
                        compression,
                        include_unassigned_psms,
                    )
                )

        sdrf_paths, run_ontology = _write_sdrf_metadata(out, output_prefix, sdrf_path, requested, compression)
        written.update(sdrf_paths)

        # Metadata tables (ontology / provenance / dataset) — written when the
        # core + run/sample structures they describe are present. Best-effort:
        # ontology entries only materialise when a PSI-MS term resolves.
        if written:
            self._write_metadata(out, output_prefix, written, consensusxml_path, run_ontology, requested, project_accession)

        return written

    def _write_metadata(
        self,
        out: Path,
        output_prefix: str,
        written: dict[str, Path],
        consensusxml_path: str,
        run_ontology: list[dict],
        requested: set[str],
        project_accession: Optional[str] = None,
    ) -> None:
        """Write ontology/provenance/dataset metadata tables (best-effort).

        Mirrors the ``openms`` enrichment path: ontology collects the run-level
        terms from the SDRF plus the score names discovered in the written core
        tables; provenance records the consensusXML -> QPX conversion steps; the
        dataset table carries the project-level metadata.
        """
        ontology_entries = list(run_ontology)
        for view in (PSM, FEATURE, PG):
            path = written.get(view)
            if not path:
                continue
            try:
                names = _collect_score_names(path)
            except (KeyError, pa.ArrowInvalid):
                continue
            if names:
                from qpx.core.scores import score_ontology_entries

                ontology_entries.extend(score_ontology_entries(names, view=view))  # noqa: PERF401
        self._write_ontology(out, output_prefix, ontology_entries)

        structures = sorted(requested & {PSM, FEATURE, PG, RUN, SAMPLE})
        provenance_records = self._build_provenance(structures, consensusxml_path)
        self._write_provenance(out, output_prefix, provenance_records)
        if provenance_records:
            self._write_dataset(
                out,
                output_prefix,
                project_accession,
                software_name=provenance_records[0]["tool_name"],
                software_version=None,
                provenance_records=provenance_records,
            )

    @staticmethod
    def _build_provenance(structures: list[str], consensusxml_path: str) -> list[dict]:
        """Provenance records: consensusXML parsing + QPX conversion."""
        is_isobaric = _consensus_is_isobaric(consensusxml_path)
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
                "output_views": [v for v in (FEATURE, PSM, PG) if v in structures],
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "openms_consensus_qpx_conversion",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": [v for v in (SAMPLE, RUN) if v in structures] + [ONTOLOGY],
            },
        ]

    @staticmethod
    def _convert_in_memory(
        consensusxml_path,
        out,
        output_prefix,
        structures,
        sdrf_path,
        creator,
        pg_top,
        compression="zstd",
        include_unassigned_psms=True,
    ) -> dict:
        """feature/psm/pg via the in-memory pyopenms map (loaded once, iterated cheaply)."""
        from qpx.converters.openms_consensus.feature_adapter import feature_map_info, resolve_enzyme
        from qpx.converters.openms_consensus.psm_adapter import _run_resolver, psm_records_for_pid

        cm = load_consensus_map(consensusxml_path)
        if sdrf_path:
            for msg in check_channels_vs_sdrf(cm, sdrf_path):
                _log.warning("consensusXML/SDRF channel mismatch: %s", msg)
        written: dict[str, Path] = {}
        want_feature, want_psm = "feature" in structures, "psm" in structures
        if want_feature or want_psm:
            # Build feature and psm per consensus feature (mirroring the streaming
            # path's element loop) so their cross-refs can be linked identically:
            # assigned PSMs first (cf order), then the unassigned PSMs.
            map_info = feature_map_info(cm)
            # Share the full protein-group membership so feature.anchor_protein and
            # feature.pg_accessions match pg (unambiguous even for shared leaders).
            group_map = accession_to_group(cm) if want_feature else None
            resolve_run = _run_resolver(cm) if want_psm else None
            seen: set = set()
            enzyme = resolve_enzyme(cm, sdrf_path)
            feat_recs: list[dict] = []
            psm_recs: list[dict] = []
            for cf in cm:
                cf_feats, cf_psms = _cf_feature_psm_records(
                    cf,
                    map_info,
                    group_map,
                    resolve_run,
                    seen,
                    want_feature=want_feature,
                    want_psm=want_psm,
                    enzyme=enzyme,
                )
                feat_recs.extend(cf_feats)
                psm_recs.extend(cf_psms)
            if want_psm and include_unassigned_psms:
                for pid in cm.getUnassignedPeptideIdentifications():
                    psm_recs.extend(psm_records_for_pid(pid, resolve_run, seen, enzyme=enzyme))
            if want_feature:
                written["feature"] = _write_view(
                    FeatureWriter,
                    out / f"{output_prefix}.feature.parquet",
                    feat_recs,
                    creator=creator,
                    compression=compression,
                    identity_composite=_FEATURE_IDENTITY_COMPOSITE,
                )
            if want_psm:
                written["psm"] = _write_view(
                    PsmWriter,
                    out / f"{output_prefix}.psm.parquet",
                    psm_recs,
                    creator=creator,
                    compression=compression,
                    identity_composite=_PSM_IDENTITY_COMPOSITE,
                )
        if "pg" in structures:
            recs = consensus_protein_groups_to_records(sdrf_path=sdrf_path, cm=cm, top=pg_top)
            written["pg"] = _write_view(
                PgWriter, out / f"{output_prefix}.pg.parquet", recs, creator=creator, compression=compression
            )
        return written
