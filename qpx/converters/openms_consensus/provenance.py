"""Search, processing and conversion provenance from consensusXML headers."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from qpx._version import __version__
from qpx.converters.orchestrator import build_provenance_step
from qpx.core.constants import FEATURE, ONTOLOGY, PG, PSM, RUN, SAMPLE
from qpx.core.scores import normalize_score_name

_log = logging.getLogger("qpx.converters.openms_consensus.converter")

# quantms-rescoring writes one of these labels instead of the engine names when it merges
# several search engines ("quantms-rescoring" before June 2026). The file does not name the engines.
_MERGED_ENGINE_LABELS = frozenset({"quantms-consensus-rescoring", "quantms-rescoring"})


def _read_consensus_header(consensusxml_path: str):
    """Parse only the consensusXML header with the streaming reader; None when unreadable.

    The full pyopenms load would defeat the streaming path and can exhaust RAM on multi-GB files.
    """
    from qpx.converters.openms_consensus.streaming import StreamingConsensusMap

    try:
        return StreamingConsensusMap(consensusxml_path)
    except (OSError, ValueError, KeyError, AttributeError, RuntimeError) as exc:
        _log.warning("could not read the consensusXML header of %s: %s", consensusxml_path, exc)
        return None


def _read_provenance_header(consensusxml_path: str):
    """Read the search engines and data processing of a consensusXML; None when unreadable."""
    from qpx.converters.openms_consensus.streaming import read_provenance_header

    try:
        return read_provenance_header(consensusxml_path)
    except (OSError, ValueError, KeyError, AttributeError, RuntimeError) as exc:
        _log.warning("could not read the consensusXML header of %s: %s", consensusxml_path, exc)
        return None


def _consensus_is_isobaric(header) -> bool:
    """Return whether a parsed header's map labels indicate an isobaric (TMT/iTRAQ) run."""
    from qpx.converters.openms_consensus.channels import _map_label

    try:
        return any(_map_label(column.label) != "LFQ" for column in header.getColumnHeaders().values())
    except (ValueError, KeyError, AttributeError) as exc:
        _log.warning("could not determine isobaric labeling: %s", exc)
        return False


def _inputs_are_isobaric(consensusxml_paths: list[str]) -> bool:
    """Whether any input's map labels are isobaric.

    The map labels follow the unassigned identifications, so this parses the whole header.
    """
    headers = (_read_consensus_header(path) for path in consensusxml_paths)
    return any(_consensus_is_isobaric(header) for header in headers if header is not None)


def _provenance_sources(headers) -> tuple[list, list, list]:
    """Deduplicate original engines, rescoring engines and processing tools."""
    original = list(dict.fromkeys(engine for header in headers for engine in header.original_search_engines))
    original_names = {name for name, _ in original}
    all_engines = list(dict.fromkeys(engine for header in headers for engine in header.search_engines))
    engines = [engine for engine in all_engines if engine[0] not in original_names]
    # One step per (tool, version, action): ProteomicsLFQ records alignment and feature grouping separately.
    software = list(dict.fromkeys(tool for header in headers for tool in header.data_processing))

    return original, engines, software


def _search_steps(original, engines, psm_views) -> list[dict]:
    """Record original searches before subsequent scoring engines."""
    steps = [
        build_provenance_step(
            step_category="rescoring" if name in _MERGED_ENGINE_LABELS else "database_search",
            step_name="engine_consensus" if name in _MERGED_ENGINE_LABELS else "sequence_search",
            tool_name=name,
            tool_version=version,
            output_views=psm_views,
        )
        for name, version in original
    ]
    category, step_name = ("rescoring", "psm_rescoring") if original else ("database_search", "sequence_search")
    steps += [
        build_provenance_step(
            step_category=category,
            step_name=step_name,
            tool_name=name,
            tool_version=version,
            output_views=psm_views,
        )
        for name, version in engines
    ]
    return steps


def _quantification_steps(software, core_views, consensusxml_paths) -> tuple[list[dict], str, Optional[str]]:
    """Quantification steps from the recorded processing tools, or one OpenMS workflow step when none is recorded."""
    steps = [
        build_provenance_step(
            step_category="quantification",
            step_name=normalize_score_name(action) or "quantification",
            tool_name=name,
            tool_version=version,
            output_views=core_views,
        )
        for name, version, action in software
    ]
    if software:
        versions = {version for _, version, _ in software}
        software_name = "OpenMS"
        software_version = next(iter(versions)) if len(versions) == 1 and "" not in versions else None
    else:
        is_isobaric = _inputs_are_isobaric(consensusxml_paths)
        software_name = "OpenMS/IsobaricWorkflow" if is_isobaric else "OpenMS/ProteomicsLFQ"
        software_version = None
        fallback_step = "isobaric_quantification" if is_isobaric else "label_free_quantification"
        steps.append(
            build_provenance_step(
                step_category="quantification",
                step_name=fallback_step,
                tool_name=software_name,
                output_views=core_views,
            )
        )
    return steps, software_name, software_version


def _build_provenance(structures: list[str], consensusxml_paths: list[str]) -> tuple[list[dict], str, Optional[str]]:
    """Provenance records and dataset software (name, version) for the consensusXML -> QPX conversion.

    Search engines and processing tools come from the consensusXML headers of every input.
    The qpx conversion step is always last.
    """
    headers = [header for header in (_read_provenance_header(path) for path in consensusxml_paths) if header is not None]
    psm_views = [PSM] if PSM in structures else []
    core_views = [v for v in (FEATURE, PSM, PG) if v in structures]
    original, engines, software = _provenance_sources(headers)

    steps = _search_steps(original, engines, psm_views)
    quantification, software_name, software_version = _quantification_steps(software, core_views, consensusxml_paths)
    steps.extend(quantification)
    steps.append(
        build_provenance_step(
            step_category="format_conversion",
            step_name="openms_consensus_qpx_conversion",
            tool_name="qpx",
            tool_version=__version__,
            parameters=[{"key": "consensusxml", "value": Path(path).name} for path in consensusxml_paths],
            output_views=[v for v in (SAMPLE, RUN) if v in structures] + [ONTOLOGY],
        )
    )
    for order, step in enumerate(steps, start=1):
        step["step_order"] = order
    return steps, software_name, software_version
