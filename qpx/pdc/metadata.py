"""Minimal PDC GraphQL client for sample/run metadata.

CDAP ``.psm`` files identify the run (``FileName``) but carry no biological
sample metadata and, for isobaric studies, no channel -> sample mapping. That
mapping lives only in the PDC GraphQL API. This module exposes the three
queries needed to rebuild it, plus a threaded helper to map each run file to
its plex (``study_run_metadata``):

    run_file_name --fileMetadata--> study_run_metadata_submitter_id (plex)
    plex          --studyExperimentalDesign--> {channel: aliquot}
    aliquot_id    --biospecimenPerStudy--> {organism, tissue, disease, ...}

All queries pass ``acceptDUA: true`` (the harmonized metadata is open access).
Read-only; no pridepy dependency.
"""

from __future__ import annotations

import json
import logging
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from typing import Optional

from qpx.core.files import run_file_stem
from qpx.core.http import safe_urlopen

logger = logging.getLogger(__name__)

PDC_API = "https://pdc.cancer.gov/graphql"
_USER_AGENT = "qpx-pdc2qpx"
_TIMEOUT_SECONDS = 120
# Global rule: PDC file-metadata lookups fan out at 24 threads by default.
DEFAULT_METADATA_THREADS = 24

# Reporter-ion channel fields in canonical mass-ascending order. The same order
# CDAP's CHANNEL_DEFS use, so a plex's non-empty channels align positionally
# with the CDAP channel labels (see pdc_sample_run.py).
TMT_CHANNEL_FIELDS: tuple[str, ...] = (
    "tmt_126",
    "tmt_127n",
    "tmt_127c",
    "tmt_128n",
    "tmt_128c",
    "tmt_129n",
    "tmt_129c",
    "tmt_130n",
    "tmt_130c",
    "tmt_131",
    "tmt_131c",
    "tmt_132n",
    "tmt_132c",
    "tmt_133n",
    "tmt_133c",
    "tmt_134n",
    "tmt_134c",
    "tmt_135n",
)
ITRAQ4_CHANNEL_FIELDS: tuple[str, ...] = ("itraq_114", "itraq_115", "itraq_116", "itraq_117")
ALL_CHANNEL_FIELDS: tuple[str, ...] = ("label_free",) + TMT_CHANNEL_FIELDS + ITRAQ4_CHANNEL_FIELDS


def post_graphql(query: str, variables: Optional[dict] = None, *, api_url: str = PDC_API) -> dict:
    """POST a GraphQL query to PDC and return the ``data`` object.

    Raises ``RuntimeError`` if the response carries ``errors``.
    """
    payload = {"query": query, "variables": variables or {}}
    request = urllib.request.Request(
        api_url,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "User-Agent": _USER_AGENT},
    )
    with safe_urlopen(request, timeout=_TIMEOUT_SECONDS) as response:
        body = json.load(response)
    if body.get("errors"):
        raise RuntimeError(f"PDC GraphQL returned errors: {body['errors']}")
    return body.get("data") or {}


_FILES_PER_STUDY_QUERY = """
query StudyFiles($studyId: String!) {
  filesPerStudy(pdc_study_id: $studyId, acceptDUA: true) {
    file_id
    file_name
    data_category
  }
}
"""

_FILE_METADATA_QUERY = """
query FileMeta($fileId: String!) {
  fileMetadata(file_id: $fileId, acceptDUA: true) {
    file_name
    study_run_metadata_submitter_id
    fraction_number
    instrument
  }
}
"""

_CHANNEL_SELECTION = "\n".join(f"    {field} {{ aliquot_id aliquot_submitter_id }}" for field in ALL_CHANNEL_FIELDS)
_EXPERIMENTAL_DESIGN_QUERY = f"""
query Design($studyId: String!) {{
  studyExperimentalDesign(pdc_study_id: $studyId, acceptDUA: true) {{
    study_run_metadata_submitter_id
    experiment_type
{_CHANNEL_SELECTION}
  }}
}}
"""

_BIOSPECIMEN_QUERY = """
query Biospecimen($studyId: String!) {
  biospecimenPerStudy(pdc_study_id: $studyId, acceptDUA: true) {
    aliquot_id
    aliquot_submitter_id
    sample_submitter_id
    sample_type
    disease_type
    primary_site
    case_submitter_id
    taxon
  }
}
"""


def fetch_study_files(pdc_study_id: str) -> list[dict]:
    """Return all files for a study (``file_id``, ``file_name``, ``data_category``)."""
    data = post_graphql(_FILES_PER_STUDY_QUERY, {"studyId": pdc_study_id})
    return data.get("filesPerStudy") or []


def fetch_experimental_design(pdc_study_id: str) -> list[dict]:
    """Return the per-plex experimental design (channel -> aliquot mapping)."""
    data = post_graphql(_EXPERIMENTAL_DESIGN_QUERY, {"studyId": pdc_study_id})
    return data.get("studyExperimentalDesign") or []


def fetch_biospecimen(pdc_study_id: str) -> dict[str, dict]:
    """Return ``{aliquot_id: biospecimen_record}`` for a study."""
    data = post_graphql(_BIOSPECIMEN_QUERY, {"studyId": pdc_study_id})
    rows = data.get("biospecimenPerStudy") or []
    return {row["aliquot_id"]: row for row in rows if row.get("aliquot_id")}


def fetch_file_metadata(file_id: str) -> Optional[dict]:
    """Return the first ``fileMetadata`` record for a file id, or ``None``."""
    data = post_graphql(_FILE_METADATA_QUERY, {"fileId": file_id})
    records = data.get("fileMetadata") or []
    return records[0] if records else None


def _strip_raw_extension(file_name: str) -> str:
    """Drop a trailing acquisition-file suffix to match CDAP ``run_file_name``."""
    return run_file_stem(file_name)


def map_runs_to_plex(pdc_study_id: str, *, threads: int = DEFAULT_METADATA_THREADS) -> dict[str, dict]:
    """Map each raw run file to its plex via ``fileMetadata``.

    PDC has no per-study bulk ``fileMetadata`` query, so each raw file is looked
    up individually. Lookups fan out across *threads* workers (metadata-only,
    fast). Returns ``{run_file_name: {plex, file_name, fraction, instrument}}``
    keyed by the extension-stripped file name (matching CDAP ``run_file_name``).
    """
    raw_files = [f for f in fetch_study_files(pdc_study_id) if f.get("data_category") == "Raw Mass Spectra"]
    if not raw_files:
        logger.warning("No raw spectra files found for %s; run->plex map will be empty", pdc_study_id)
        return {}

    def _lookup(entry: dict) -> Optional[tuple[str, dict]]:
        meta = fetch_file_metadata(entry["file_id"])
        if not meta or not meta.get("study_run_metadata_submitter_id"):
            return None
        original_name = meta.get("file_name") or entry.get("file_name") or ""
        run_file_name = _strip_raw_extension(original_name)
        return run_file_name, {
            "plex": meta["study_run_metadata_submitter_id"],
            "file_name": original_name or None,
            "fraction": meta.get("fraction_number"),
            "instrument": meta.get("instrument"),
        }

    run_to_plex: dict[str, dict] = {}
    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        for result in pool.map(_lookup, raw_files):
            if result is not None:
                run_to_plex[result[0]] = result[1]
    logger.info("Mapped %d/%d raw run files to plexes for %s", len(run_to_plex), len(raw_files), pdc_study_id)
    return run_to_plex
