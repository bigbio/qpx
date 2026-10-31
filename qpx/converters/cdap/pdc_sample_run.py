"""Build CPTAC/PDC ``sample.parquet`` and ``run.parquet`` from PDC metadata.

CDAP ``.psm`` files carry the run (``FileName``) but no biological sample
metadata, and for isobaric (TMT/iTRAQ) studies no channel -> sample mapping.
For 79% of CPTAC studies the reporter intensities in ``feature.parquet`` are
therefore uninterpretable on their own. The missing link lives only in the PDC
GraphQL API and is recovered here:

    run_file_name --fileMetadata--> plex (study_run_metadata)
    plex          --studyExperimentalDesign--> {channel: aliquot}
    aliquot_id    --biospecimenPerStudy--> {organism, tissue, disease, case}

The reporter channels use the **canonical SDRF labels** (``TMT126``,
``ITRAQ114``, ``LFQ`` for label-free), as ``feature.parquet`` does, so that each
``run.samples[].label`` joins directly to an intensity label.

This builder lives at the ``pdc2qpx`` layer (it needs network), not inside the
offline :class:`~qpx.converters.cdap.converter.CdapConverter`.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from qpx.converters.cdap.constants import CANONICAL_LABELS, CHANNEL_DEFS
from qpx.pdc.metadata import (
    DEFAULT_METADATA_THREADS,
    ITRAQ4_CHANNEL_FIELDS,
    TMT_CHANNEL_FIELDS,
    fetch_biospecimen,
    fetch_experimental_design,
    map_runs_to_plex,
)
from qpx.writers.run import RunWriter
from qpx.writers.sample import SampleWriter

logger = logging.getLogger(__name__)

_LFQ_LABEL = "LFQ"
_SAMPLE_EXTRA_COLUMNS = ["sample_type"]
_DEFAULT_ORGANISM = "Homo sapiens"
_NOT_REPORTED = "Not Reported"


def _is_reference(aliquot_submitter_id: Optional[str]) -> bool:
    """Return True for pooled internal-reference channels (no real case)."""
    text = (aliquot_submitter_id or "").lower()
    return "reference" in text or "pooled" in text


def _plex_channel_map(plex: dict) -> list[tuple[str, str, str]]:
    """Map one plex to ``[(label, aliquot_id, aliquot_submitter_id)]``.

    Aligns the plex's non-empty PDC channel fields (canonical mass order) with
    the CDAP channel labels for its ``experiment_type``. Label-free plexes yield
    a single ``LFQ`` channel. Returns ``[]`` (with a warning) for experiment
    types CDAP does not model (e.g. TMT16) or on a count mismatch, so a wrong
    mapping is never emitted.
    """
    if plex.get("label_free"):
        entry = plex["label_free"][0]
        return [(_LFQ_LABEL, entry["aliquot_id"], entry.get("aliquot_submitter_id") or entry["aliquot_id"])]

    experiment_type = plex.get("experiment_type")
    cdap_labels = CHANNEL_DEFS.get(experiment_type)
    if not cdap_labels:
        logger.warning(
            "Plex %s has experiment_type %r not modelled by CDAP; skipping its channels",
            plex.get("study_run_metadata_submitter_id"),
            experiment_type,
        )
        return []

    canonical = ITRAQ4_CHANNEL_FIELDS if experiment_type.startswith("iTRAQ") else TMT_CHANNEL_FIELDS
    present = [field for field in canonical if plex.get(field)]

    expected_labels = list(cdap_labels)
    # TMT11 pilots may also carry the legacy TMTXX 132N/C channels; CDAP appends
    # them for TMT11 (see CdapBaseAdapter._detect_channels), so mirror that here
    # instead of skipping the whole plex on the resulting count mismatch.
    if experiment_type == "TMT11" and (plex.get("tmt_132n") or plex.get("tmt_132c")):
        expected_labels.extend(CHANNEL_DEFS["TMTXX"])

    if len(present) != len(expected_labels):
        logger.warning(
            "Plex %s: %d PDC channels but %d CDAP labels for %s; skipping its channels",
            plex.get("study_run_metadata_submitter_id"),
            len(present),
            len(expected_labels),
            experiment_type,
        )
        return []

    mapping: list[tuple[str, str, str]] = []
    for field, cdap_label in zip(present, expected_labels):
        entry = plex[field][0]
        mapping.append(
            (CANONICAL_LABELS[cdap_label], entry["aliquot_id"], entry.get("aliquot_submitter_id") or entry["aliquot_id"])
        )
    return mapping


def _sample_record(aliquot_id: str, aliquot_submitter_id: str, biospecimen: dict[str, dict]) -> dict:
    """Build one ``sample.parquet`` record from a biospecimen aliquot."""
    if _is_reference(aliquot_submitter_id):
        return {
            "sample_accession": aliquot_submitter_id,
            "organism": _DEFAULT_ORGANISM,
            "organism_part": _NOT_REPORTED,
            "disease": None,
            "individual": None,
            "sample_type": "Internal Reference",
        }
    meta = biospecimen.get(aliquot_id) or {}
    return {
        "sample_accession": aliquot_submitter_id,
        "organism": meta.get("taxon") or _DEFAULT_ORGANISM,
        "organism_part": meta.get("primary_site") or _NOT_REPORTED,
        "disease": meta.get("disease_type"),
        "individual": meta.get("case_submitter_id"),
        "sample_type": meta.get("sample_type"),
    }


def _write_sample_run(
    sample_path: Path,
    run_path: Path,
    sample_records: list[dict],
    run_records: list[dict],
    compression: str,
) -> None:
    """Write the sample/run parquet via temp files + atomic rename.

    A writer emits no file for an empty batch (a run with no mappable channels),
    so only rename temps that were actually written. Temp files are removed on
    any failure, so a crash never clobbers a previously-good parquet from an
    earlier run.
    """
    sample_tmp = sample_path.with_name(sample_path.name + ".tmp")
    run_tmp = run_path.with_name(run_path.name + ".tmp")
    try:
        with SampleWriter(
            sample_tmp,
            creator="pdc",
            extra_columns=_SAMPLE_EXTRA_COLUMNS,
            compression=compression,
        ) as writer:
            writer.write_batch(sample_records)
        with RunWriter(run_tmp, creator="pdc", compression=compression) as writer:
            writer.write_batch(run_records)
        if sample_tmp.exists():
            sample_tmp.replace(sample_path)
        if run_tmp.exists():
            run_tmp.replace(run_path)
    finally:
        for tmp in (sample_tmp, run_tmp):
            if tmp.exists():
                tmp.unlink()


def _pdc_sample_records(plex_channels: dict, biospecimen: dict) -> dict[str, dict]:
    """Build one sample record per referenced aliquot."""
    # Sample rows: one per referenced aliquot, keyed by aliquot_submitter_id.
    sample_records: dict[str, dict] = {}
    for channels in plex_channels.values():
        for _label, aliquot_id, aliquot_submitter_id in channels:
            if aliquot_submitter_id not in sample_records:
                sample_records[aliquot_submitter_id] = _sample_record(aliquot_id, aliquot_submitter_id, biospecimen)

    return sample_records


def _pdc_run_records(run_to_plex: dict, plex_channels: dict) -> tuple[list[dict], int]:
    """Build run records and count runs without channel mappings."""
    # Run rows: one per run file, samples = its plex's channels.
    run_records: list[dict] = []
    runs_without_channels = 0
    for run_file_name, info in sorted(run_to_plex.items()):
        channels = plex_channels.get(info["plex"], [])
        if not channels:
            runs_without_channels += 1
        samples = [
            {
                "sample_accession": aliquot_submitter_id,
                "label": label,
                "biological_replicate": None,
                "technical_replicate": None,
            }
            for label, _aliquot_id, aliquot_submitter_id in channels
        ]
        run_records.append(
            {
                "run_accession": run_file_name,
                "run_file_name": run_file_name,
                "file_name": info.get("file_name"),
                "samples": samples,
                "fraction": str(info["fraction"]) if info.get("fraction") is not None else None,
                "instrument": info.get("instrument"),
                "enzymes": None,
                "dissociation_method": None,
                "modification_parameters": None,
            }
        )

    return run_records, runs_without_channels


def run_units_from_plexes(run_to_plex: dict[str, dict]) -> dict[str, list[str]]:
    """Map each run to the runs of its plex: PDC's own quantification unit."""
    plexes: dict[str, list[str]] = {}
    for run, info in sorted(run_to_plex.items()):
        plexes.setdefault(info["plex"], []).append(run)
    return {run: plexes[info["plex"]] for run, info in run_to_plex.items()}


def build_sample_run_from_pdc(
    pdc_study_id: str,
    output_folder: str | Path,
    prefix: str,
    *,
    compression: str = "zstd",
    threads: int = DEFAULT_METADATA_THREADS,
    run_to_plex: Optional[dict[str, dict]] = None,
) -> Optional[dict[str, Path]]:
    """Build ``<prefix>.sample.parquet`` and ``<prefix>.run.parquet`` from PDC.

    Returns ``{"sample": path, "run": path}`` on success or ``None`` when the
    run -> plex map is empty (study has no harmonized metadata / no raw files).
    Network errors propagate; the caller (pdc2qpx) wraps this in a warn-skip
    guard so metadata failure never breaks the main conversion. ``run_to_plex``
    reuses a :func:`~qpx.pdc.metadata.map_runs_to_plex` result already fetched.
    """
    output_folder = Path(output_folder)

    design = fetch_experimental_design(pdc_study_id)
    biospecimen = fetch_biospecimen(pdc_study_id)
    if run_to_plex is None:
        run_to_plex = map_runs_to_plex(pdc_study_id, threads=threads)

    if not run_to_plex:
        logger.warning("No run->plex mapping for %s; skipping sample/run views", pdc_study_id)
        return None

    # plex submitter id -> channel mapping
    plex_channels = {
        plex["study_run_metadata_submitter_id"]: _plex_channel_map(plex)
        for plex in design
        if plex.get("study_run_metadata_submitter_id")
    }

    sample_records = _pdc_sample_records(plex_channels, biospecimen)

    run_records, runs_without_channels = _pdc_run_records(run_to_plex, plex_channels)

    if runs_without_channels:
        logger.warning(
            "%d/%d runs in %s have no channel mapping (unsupported experiment type or missing plex)",
            runs_without_channels,
            len(run_records),
            pdc_study_id,
        )

    sample_path = output_folder / f"{prefix}.sample.parquet"
    run_path = output_folder / f"{prefix}.run.parquet"
    _write_sample_run(sample_path, run_path, list(sample_records.values()), run_records, compression)

    logger.info(
        "Wrote %d samples and %d runs from PDC metadata for %s",
        len(sample_records),
        len(run_records),
        pdc_study_id,
    )
    return {"sample": sample_path, "run": run_path}
