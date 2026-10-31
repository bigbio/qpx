"""Score normalization and CV term lookup.

All QPX score names use snake_case — no colons, dots, spaces, or mixed case.
This module provides:
    - ``normalize_score_name()``: convert raw tool score names to snake_case
    - ``is_higher_better()``: determine score direction
    - ``score_ontology_entries()``: generate ontology.parquet rows for discovered scores

Score CV term information comes from two sources:
    1. A small built-in cache for scores not in any OBO file (tool-specific
       scores like DIA-NN, MSFragger, QPX-internal metrics).
    2. The PSI-MS controlled vocabulary via ``PublicOntology("psi_ms")``,
       backed by a Parquet file with three-layer resolution
       (repo checkout > auto-updated cache > bundled fallback).
"""

from __future__ import annotations

import logging
import re
from typing import Optional

from qpx.core.builtin_scores import BUILTIN_SCORES

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Score name normalization
# ---------------------------------------------------------------------------

# Characters that get replaced by underscore
_NON_IDENT = re.compile(r"[^a-z0-9]+")
# Collapse multiple underscores and strip leading/trailing
_MULTI_UNDER = re.compile(r"_+")


def normalize_score_name(raw_name: str) -> str:
    """Normalize a raw score name to snake_case.

    Lowercases, replaces colons / dots / spaces / special chars with ``_``,
    collapses runs of underscores, and strips leading/trailing underscores.

    Examples::

        >>> normalize_score_name("percolator:score")
        'percolator_score'
        >>> normalize_score_name("MSFragger:Hyperscore")
        'msfragger_hyperscore'
        >>> normalize_score_name("Comet:xcorr")
        'comet_xcorr'
        >>> normalize_score_name("x!tandem:expect")
        'x_tandem_expect'
        >>> normalize_score_name("MS-GF:RawScore")
        'ms_gf_rawscore'
        >>> normalize_score_name("andromeda_score")
        'andromeda_score'
    """
    name = raw_name.lower()
    name = _NON_IDENT.sub("_", name)
    name = _MULTI_UNDER.sub("_", name)
    return name.strip("_")


# ---------------------------------------------------------------------------
# PublicOntology-backed CV lookup (lazy singleton)
# ---------------------------------------------------------------------------

_ontology = None  # type: Optional["PublicOntology"]


def _get_ontology():
    """Return the lazily-loaded ``PublicOntology`` singleton for PSI-MS.

    Uses ``auto_update=False`` to avoid network calls during conversions.
    The Parquet file is resolved via three-layer lookup: repo checkout,
    cached download, or bundled fallback.
    """
    global _ontology
    if _ontology is None:
        from qpx.core.ontology import PublicOntology

        _ontology = PublicOntology("psi_ms", auto_update=False)
    return _ontology


def _lookup_from_obo(normalized_name: str) -> Optional[dict]:
    """Try to resolve a score from the PSI-MS PublicOntology."""
    term = _get_ontology().lookup_normalized(normalized_name)
    if term is None:
        return None

    higher_better = term["higher_better"]
    if higher_better is None:
        higher_better = _direction_from_name(normalize_score_name(term["name"]))
    return {
        "ontology_name": term["name"],
        "ontology_accession": term["accession"],
        "ontology_source": term["source"],
        "description": term["definition"],
        "higher_better": higher_better,
    }


_lookup_cache: dict[str, Optional[dict]] = {}


def lookup_score(normalized_name: str) -> Optional[dict]:
    """Look up CV term info for a normalized score name.

    Checks the built-in cache first, then falls back to the PSI-MS OBO file.
    Results are cached in ``_lookup_cache`` so that repeated calls for the
    same score name (common during row-by-row conversion) avoid redundant
    DuckDB queries against the ontology Parquet file.

    Returns ``None`` if the score is not in either source.
    """
    # 0. Check runtime cache (populated on first miss)
    if normalized_name in _lookup_cache:
        return _lookup_cache[normalized_name]

    # 1. Check built-in cache
    entry = BUILTIN_SCORES.get(normalized_name)
    if entry is not None:
        _lookup_cache[normalized_name] = entry
        return entry

    # 2. Check OBO registry
    result = _lookup_from_obo(normalized_name)
    _lookup_cache[normalized_name] = result
    return result


# ---------------------------------------------------------------------------
# Score direction
# ---------------------------------------------------------------------------


# Most PSI-MS score terms carry no has_order relation. Error rates and
# probabilities of being wrong are lower-is-better whatever the engine.
_LOWER_IS_BETTER = re.compile(r"e_?value|expect|q_?value|p_?value|fdr|posterior_error|(?:^|_)pep(?:_|$)")
# Engine match scores and probabilities of being right are higher-is-better.
_HIGHER_IS_BETTER = re.compile(r"score$|xcorr|probability$")


def _direction_from_name(normalized_name: str) -> Optional[bool]:
    """Score direction implied by a score's name, ``None`` when the name does not tell."""
    if _LOWER_IS_BETTER.search(normalized_name):
        return False
    if _HIGHER_IS_BETTER.search(normalized_name):
        return True
    return None


def is_higher_better(normalized_name: str) -> Optional[bool]:
    """Return whether the score with *normalized_name* is higher-is-better.

    Checks built-in cache, then OBO registry. A score whose direction neither
    declares is lower-is-better when it is an e-value, q-value, p-value, PEP
    or FDR, and ``None`` (unknown) otherwise; a consumer must not guess.
    """
    info = lookup_score(normalized_name)
    if info is not None and info["higher_better"] is not None:
        return info["higher_better"]
    return _direction_from_name(normalized_name)


# ---------------------------------------------------------------------------
# Field-level CV term mappings
# ---------------------------------------------------------------------------

# Standard QPX fields with known PSI-MS CV accessions
_FIELD_CV_MAP: dict[str, dict] = {
    "posterior_error_probability": {
        "ontology_name": "posterior error probability",
        "ontology_accession": "MS:1001493",
        "ontology_source": "MS",
        "description": "Posterior error probability (PEP) for PSM correctness",
    },
    "observed_mz": {
        "ontology_name": "selected ion m/z",
        "ontology_accession": "MS:1001475",
        "ontology_source": "MS",
        "description": "Experimentally measured precursor m/z",
    },
    "calculated_mz": {
        "ontology_name": "calculated mass to charge ratio",
        "ontology_accession": "MS:1001234",
        "ontology_source": "MS",
        "description": "Theoretically calculated m/z of the peptide",
    },
    "rt": {
        "ontology_name": "scan start time",
        "ontology_accession": "MS:1000016",
        "ontology_source": "MS",
        "description": "Retention time at peak apex",
    },
    "charge": {
        "ontology_name": "charge state",
        "ontology_accession": "MS:1000041",
        "ontology_source": "MS",
        "description": "Precursor ion charge state",
    },
    "intensity": {
        "ontology_name": "MS1 feature area",
        "ontology_accession": "MS:1002498",
        "ontology_source": "MS",
        "description": "Primary precursor intensity",
    },
    "missed_cleavages": {
        "ontology_name": "number of missed cleavages",
        "ontology_accession": "MS:1003044",
        "ontology_source": "MS",
        "description": "Number of amino acid residue bonds that should have been cleaved by the cleavage agent used, but were not",
    },
    "sequence": {
        "ontology_name": "stripped peptide sequence",
        "ontology_accession": "MS:1000888",
        "ontology_source": "MS",
        "description": "The amino acid sequence of a peptide without modifications",
    },
    "peptidoform": {
        "ontology_name": "peptidoform sequence",
        "ontology_accession": "MS:1000889",
        "ontology_source": "MS",
        "description": "A string representation of an amino acid sequence with modifications",
    },
    "is_decoy": {
        "ontology_name": "decoy peptide",
        "ontology_accession": "MS:1002217",
        "ontology_source": "MS",
        "description": "A peptide from a decoy or randomized database used for FDR estimation",
    },
    "predicted_rt": {
        "ontology_name": "predicted retention time",
        "ontology_accession": "MS:1000897",
        "ontology_source": "MS",
        "description": "Predicted retention time of the peptide (seconds)",
    },
}

# A shared QPX field name can represent a different level in each view. Resolve
# those cases through the corresponding score definition rather than attaching
# precursor-level semantics to a protein-group field.
_VIEW_FIELD_SCORE_ALIASES: dict[tuple[str, str], str] = {
    ("pg", "global_qvalue"): "pg_global_qvalue",
}


# ---------------------------------------------------------------------------
# Ontology entry generation
# ---------------------------------------------------------------------------


def field_ontology_entries(
    view: str = "psm",
    resolved_mappings: dict[str, str] | None = None,
    tool_name: str | None = None,
) -> list[dict]:
    """Generate ontology.parquet rows for field-to-CV mappings.

    When *resolved_mappings* is provided, produces entries with
    ``source_column_name`` and ``source_tool`` populated from the
    actual converter resolution. Fields in the resolved mapping that
    have no CV term in ``_FIELD_CV_MAP`` still produce an entry
    (with ``ontology_accession=None``).

    For backward compatibility, calling without *resolved_mappings*
    emits entries for all ``_FIELD_CV_MAP`` fields with null source info.

    Args:
        view: The QPX view name (e.g. ``"psm"``, ``"feature"``).
        resolved_mappings: QPX field -> resolved tool column name.
        tool_name: Tool name string (e.g. ``"DIA-NN"``).

    Returns:
        List of dicts matching the ``OntologySchema``.
    """
    ontology_version = _get_ontology().version
    if ontology_version == "unknown":
        ontology_version = None

    entries: list[dict] = []

    if resolved_mappings is not None:
        # New path: produce entries for each resolved field
        for field_name, source_column in sorted(resolved_mappings.items()):
            cv_info = _FIELD_CV_MAP.get(field_name)
            if cv_info is None:
                score_name = _VIEW_FIELD_SCORE_ALIASES.get((view, field_name))
                cv_info = lookup_score(score_name) if score_name else _lookup_from_obo(field_name)
            entries.append(
                {
                    "field_name": field_name,
                    "ontology_name": cv_info["ontology_name"] if cv_info else None,
                    "ontology_accession": (cv_info.get("ontology_accession") if cv_info else None),
                    "ontology_source": (cv_info.get("ontology_source") if cv_info else None),
                    "ontology_version": (ontology_version if cv_info and cv_info.get("ontology_source") == "MS" else None),
                    "view": view,
                    "description": (
                        cv_info["description"] if cv_info else f"{tool_name or 'unknown'} {source_column} (no CV term)"
                    ),
                    "source_column_name": source_column,
                    "source_tool": tool_name,
                }
            )
    else:
        # Backward-compatible path: emit _FIELD_CV_MAP entries
        for field_name, info in sorted(_FIELD_CV_MAP.items()):
            entries.append(
                {
                    "field_name": field_name,
                    "ontology_name": info["ontology_name"],
                    "ontology_accession": info["ontology_accession"],
                    "ontology_source": info["ontology_source"],
                    "ontology_version": ontology_version,
                    "view": view,
                    "description": info["description"],
                    "source_column_name": None,
                    "source_tool": None,
                }
            )

    return entries


def score_ontology_entries(
    score_names: set[str],
    view: str = "psm",
) -> list[dict]:
    """Generate ontology.parquet rows for a set of normalized score names.

    Looks up each score in the built-in cache and the PSI-MS OBO registry.
    Scores not found in either source are skipped.

    Args:
        score_names: Set of normalized (snake_case) score names found in data.
        view: The QPX view name (e.g. ``"psm"``, ``"feature"``, ``"pg"``).

    Returns:
        List of dicts matching the ``OntologySchema``.
    """
    # Fetch ontology version for provenance tracking
    ontology_version = _get_ontology().version
    if ontology_version == "unknown":
        ontology_version = None

    entries: list[dict] = []
    for name in sorted(score_names):
        info = lookup_score(name)
        if info is None:
            continue
        entries.append(
            {
                "field_name": name,
                "ontology_name": info["ontology_name"],
                "ontology_accession": info["ontology_accession"],
                "ontology_source": info["ontology_source"],
                "ontology_version": ontology_version if info["ontology_source"] == "MS" else None,
                "view": view,
                "description": info["description"],
                "source_column_name": None,
                "source_tool": None,
            }
        )
    return entries
