"""Orchestrate consensusXML (+ SDRF) -> a QPX dataset (feature/psm/pg [+ run/sample]).

Interim quantms path while OpenMS ``-out_qpx`` is pre-1.1. pg carries the protein
abundance OpenMS recorded in the file, else an unnormalized unique-peptide sum.
run/sample come from the SDRF (reusing :class:`SdrfConverter`) when an
SDRF is provided; when it is, the consensusXML channels are also checked against
the SDRF ``comment[label]`` and mismatches are logged as warnings.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from pathlib import Path
from typing import Optional, Union

import pyarrow as pa
import pyarrow.parquet as pq

from qpx.converters.openms_consensus.conversion import (
    _ConversionConfig,
    _convert_in_memory,
    _convert_streaming,
    _OutputConfig,
)
from qpx.converters.openms_consensus.provenance import _build_provenance
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.core.constants import DATASET, FEATURE, ONTOLOGY, PG, PROVENANCE, PSM, RUN, SAMPLE
from qpx.core.scores import score_ontology_entries

_log = logging.getLogger(__name__)
_STRUCTURE_ALL = ("feature", "psm", "pg", "run", "sample")


def _as_paths(consensusxml_path) -> list[str]:
    """Normalise one path or a sequence of paths to a non-empty list of strings."""
    if isinstance(consensusxml_path, (str, os.PathLike)):
        return [str(consensusxml_path)]
    paths = [str(p) for p in consensusxml_path]
    if not paths:
        raise ValueError("At least one consensusXML input is required")
    return paths


# Above this decompressed consensusXML size, auto-select the streaming reader. A whole
# in-memory conversion peaks at ~7x the file in RAM (~7 GB at this size); streaming
# peaks at ~2.3x and takes ~1.7x as long.
_STREAM_THRESHOLD_BYTES = 1024**3


def _should_stream(path: str) -> bool:
    from qpx.converters.openms_consensus.streaming import exceeds_xml_size

    try:
        return exceeds_xml_size(path, _STREAM_THRESHOLD_BYTES)
    except OSError:
        return False


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


def _remove_orphaned_metadata(output_folder: Path, output_prefix: str) -> None:
    """Remove metadata after an empty rerun only if no same-prefix data remains."""
    if any((output_folder / f"{output_prefix}.{view}.parquet").is_file() for view in _STRUCTURE_ALL):
        return
    for view in (DATASET, ONTOLOGY, PROVENANCE):
        (output_folder / f"{output_prefix}.{view}.parquet").unlink(missing_ok=True)


class OpenMSConsensusConverter(BaseOrchestrator):  # pylint: disable=too-few-public-methods
    """consensusXML + SDRF -> QPX views.

    A single-entry orchestrator (``convert``) — the interim counterpart to the
    other converter classes, kept as a class for call-site symmetry with them.
    """

    def convert(
        self,
        consensusxml_path: Union[str, os.PathLike, Sequence[Union[str, os.PathLike]]],
        output_folder: str,
        output_prefix: str = "openms",
        sdrf_path: Optional[str] = None,
        structures: tuple[str, ...] = _STRUCTURE_ALL,
        creator: str = "openms-consensus",
        streaming: Optional[bool] = None,
        project_accession: Optional[str] = None,
        compression: str = "zstd",
        include_unassigned_psms: bool = True,
    ) -> dict[str, Path]:
        """Write the requested QPX views and return ``{structure: parquet path}``.

        ``consensusxml_path`` is one consensusXML or a sequence of them. Several
        inputs (e.g. one per sample group) become ONE dataset: each file keeps its
        own column -> run mapping and identification metadata, rows go to the same
        views, pg spans all inputs and provenance lists every input. Two inputs
        must not contain the same run (``ValueError``).

        ``include_unassigned_psms`` (default ``True``) controls whether
        unassigned PeptideIdentifications from the source reach ``psm.parquet``.
        They are retained to preserve identification evidence; pass ``False``
        to exclude them.

        Requested feature/PSM/PG views with no exportable records are skipped
        with a warning and excluded from output metadata. After successful core
        conversion, remove any existing file for each skipped view and this
        output prefix so Dataset cannot discover stale records. Other requested
        views are still exported; an entirely empty export returns an empty dict.
        When no same-prefix core or run/sample views remain, their orphaned
        ontology/provenance/dataset metadata is also removed.

        ``feature_id`` records a link in the exported dataset, not quantification
        status. It is only populated when both feature and PSM views are emitted;
        PSM-only output leaves it null even for assigned identifications.
        ``PSM.with_feature()`` and ``PSM.without_feature()`` filter these links.

        The feature view always omits consensus features without peptide hits.
        Protein inference always uses every identification, regardless of this
        option.

        ``structures`` selects which of feature/psm/pg/run/sample to emit. pg
        carries the protein abundance OpenMS recorded in the file, else an
        unnormalized unique-peptide sum. The sum stays an additional intensity
        beside OpenMS' value.

        ``streaming`` picks the consensusXML reader: ``None`` (default) auto-selects
        — the low-memory streaming reader when any input is above ``_STREAM_THRESHOLD_BYTES``
        (an in-memory conversion peaks at ~7x the file in RAM, streaming at ~2.3x), else
        the faster in-memory pyopenms load. ``True``/``False`` forces the choice.

        ``project_accession`` (e.g. ``PXD001819``) is stamped into dataset.parquet.

        ``compression`` is the Parquet codec for every written view (default zstd).
        """
        self._compression = compression
        _validate_structures(structures, sdrf_path)
        consensusxml_paths = _as_paths(consensusxml_path)
        requested = set(structures)

        out = Path(output_folder)
        out.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}

        config = _ConversionConfig(
            _OutputConfig(out, output_prefix, creator, compression), structures, sdrf_path, include_unassigned_psms
        )

        if {"feature", "psm", "pg"}.intersection(structures):
            use_stream = streaming if streaming is not None else any(_should_stream(p) for p in consensusxml_paths)
            if use_stream:
                # A shared identification pre-pass precedes batched output of all
                # views; the complete XML map is never retained in memory.
                written.update(_convert_streaming(consensusxml_paths, config))
            else:
                # In-memory path: pyopenms loads the map once (fast for smaller files);
                # the adapters iterate it cheaply. Both paths write the same rows; a few
                # float values (PEP, mass error) can differ in the last bit.
                written.update(_convert_in_memory(consensusxml_paths, config))

        for view in (FEATURE, PSM, PG):
            if view in requested and view not in written:
                (out / f"{output_prefix}.{view}.parquet").unlink(missing_ok=True)
                _log.warning("No exportable %s records; skipping %s Parquet output.", view.upper(), view.upper())

        sdrf_paths, run_ontology = _write_sdrf_metadata(out, output_prefix, sdrf_path, requested, compression)
        written.update(sdrf_paths)

        # Metadata tables (ontology / provenance / dataset) — written when the
        # core + run/sample structures they describe are present. Best-effort:
        # ontology entries only materialise when a PSI-MS term resolves.
        if written:
            self._write_metadata(config, written, consensusxml_paths, run_ontology, project_accession)
        else:
            _remove_orphaned_metadata(out, output_prefix)

        return written

    def _write_metadata(
        self,
        config: _ConversionConfig,
        written: dict[str, Path],
        consensusxml_paths: list[str],
        run_ontology: list[dict],
        project_accession: Optional[str] = None,
    ) -> None:
        """Write ontology/provenance/dataset metadata tables (best-effort).

        Mirrors the ``openms`` enrichment path: ontology collects the run-level
        terms from the SDRF plus the score names discovered in the written core
        tables; provenance records the consensusXML -> QPX conversion steps; the
        dataset table carries the project-level metadata.
        """
        out, output_prefix = config.output.folder, config.output.prefix
        requested = set(config.structures)
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
                ontology_entries.extend(score_ontology_entries(names, view=view))  # noqa: PERF401
        self._write_ontology(out, output_prefix, ontology_entries)

        structures = sorted(requested.intersection(written, {PSM, FEATURE, PG, RUN, SAMPLE}))
        provenance_records, software_name, software_version = _build_provenance(structures, consensusxml_paths)
        self._write_provenance(out, output_prefix, provenance_records)
        # No provenance_records here: a None version would take the first step's version (the search engine).
        self._write_dataset(
            out,
            output_prefix,
            project_accession,
            software_name=software_name,
            software_version=software_version,
        )
