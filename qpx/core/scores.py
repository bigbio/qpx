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
# Built-in cache for scores NOT in PSI-MS OBO
# ---------------------------------------------------------------------------
# Only tool-specific scores that have no CV accession, or QPX-internal
# metrics.  Everything with an MS: accession is resolved from the OBO file.

_BUILTIN_SCORES: dict[str, dict] = {
    # --- DIA-NN (no MS accessions yet) ---
    "diann_qvalue": {
        "ontology_name": "DIA-NN:Q.Value",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "DIA-NN run-level q-value (lower is better)",
        "higher_better": False,
    },
    "precursor_qvalue": {
        "ontology_name": "DIA-NN:Q.Value",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "DIA-NN run-level precursor q-value (lower is better)",
        "higher_better": False,
    },
    "diann_global_qvalue": {
        "ontology_name": "DIA-NN:Global.Q.Value",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "DIA-NN global q-value (lower is better)",
        "higher_better": False,
    },
    "diann_cscore": {
        "ontology_name": "DIA-NN:CScore",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "DIA-NN confidence score (higher is better)",
        "higher_better": True,
    },
    # --- MSFragger / FragPipe (no MS accessions) ---
    "msfragger_hyperscore": {
        "ontology_name": "MSFragger:Hyperscore",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "MSFragger hyperscore (higher is better)",
        "higher_better": True,
    },
    "msfragger_expectation": {
        "ontology_name": "MSFragger:Expectation",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "MSFragger expectation value (lower is better)",
        "higher_better": False,
    },
    # --- Sage (no MS accessions yet) ---
    "sage_hyperscore": {
        "ontology_name": "Sage:hyperscore",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "Sage hyperscore (higher is better)",
        "higher_better": True,
    },
    # --- QPX-internal metrics ---
    "parent_ion_fraction": {
        "ontology_name": "parent ion fraction",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "Fraction of target peak intensity in isolation window (higher is better)",
        "higher_better": True,
    },
    "consensus_support": {
        "ontology_name": "consensus support",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "Number of search engines supporting the identification (higher is better)",
        "higher_better": True,
    },
    "precursor_quantification_score": {
        "ontology_name": "precursor quantification score",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "QuantUMS-derived quantification confidence (higher is better)",
        "higher_better": True,
    },
    # --- DIA-NN additional scores ---
    "diann_evidence": {
        "ontology_name": "DIA-NN:Evidence",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "DIA-NN co-elution evidence score (higher is better)",
        "higher_better": True,
    },
    "diann_spectrum_similarity": {
        "ontology_name": "spectral dot product",
        "ontology_accession": "MS:1003304",
        "ontology_source": "MS",
        "description": "Observed vs library MS2 spectrum similarity (higher is better)",
        "higher_better": True,
    },
    "diann_ms1_profile_corr": {
        "ontology_name": "DIA-NN:Ms1.Profile.Corr",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "MS1/MS2 chromatographic profile correlation (higher is better)",
        "higher_better": True,
    },
    "diann_mass_evidence": {
        "ontology_name": "DIA-NN:Mass.Evidence",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "Mass accuracy scoring contribution (higher is better)",
        "higher_better": True,
    },
    "diann_averagine": {
        "ontology_name": "DIA-NN:Averagine",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "Isotope pattern fit against Averagine model (higher is better)",
        "higher_better": True,
    },
    "diann_quantity_quality": {
        "ontology_name": "DIA-NN:Quantity.Quality",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "Quantification quality 1/(1+SD) (higher is better)",
        "higher_better": True,
    },
    "diann_lib_qvalue": {
        "ontology_name": "DIA-NN:Lib.Q.Value",
        "ontology_accession": "MS:1002354",
        "ontology_source": "MS",
        "description": "Library-pass precursor q-value (lower is better)",
        "higher_better": False,
    },
    "diann_lib_pg_qvalue": {
        "ontology_name": "DIA-NN:Lib.PG.Q.Value",
        "ontology_accession": "MS:1002373",
        "ontology_source": "MS",
        "description": "Library-pass protein group q-value (lower is better)",
        "higher_better": False,
    },
    "diann_protein_qvalue": {
        "ontology_name": "protein-level q-value",
        "ontology_accession": "MS:1001869",
        "ontology_source": "MS",
        "description": "Protein-level FDR using proteotypic peptides (lower is better)",
        "higher_better": False,
    },
    "diann_decoy_cscore": {
        "ontology_name": "DIA-NN:Decoy.CScore",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "CScore for paired decoy match (higher is better)",
        "higher_better": True,
    },
    "diann_decoy_evidence": {
        "ontology_name": "DIA-NN:Decoy.Evidence",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "Evidence score for paired decoy match (higher is better)",
        "higher_better": True,
    },
    "diann_translated_quality": {
        "ontology_name": "DIA-NN:Translated.Quality",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "plexDIA translation quality (higher is better)",
        "higher_better": True,
    },
    "diann_translated_qvalue": {
        "ontology_name": "DIA-NN:Translated.Q.Value",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "plexDIA translated q-value (lower is better)",
        "higher_better": False,
    },
    # --- q-value variants used directly by converters ---
    # Names preserve DIA-NN / tool origin for backward compatibility.
    "qvalue": {
        "ontology_name": "q-value (DIA-NN:Q.Value)",
        "ontology_accession": "MS:1002354",
        "ontology_source": "MS",
        "description": "Run-level precursor q-value (DIA-NN Q.Value, lower is better)",
        "higher_better": False,
    },
    "pg_qvalue": {
        "ontology_name": "protein group q-value (DIA-NN:PG.Q.Value)",
        "ontology_accession": "MS:1002373",
        "ontology_source": "MS",
        "description": "Protein group run-level q-value (DIA-NN PG.Q.Value, lower is better)",
        "higher_better": False,
    },
    "protein_global_qvalue": {
        "ontology_name": "protein-level global FDR",
        "ontology_accession": "MS:1001214",
        "ontology_source": "MS",
        "description": "Protein global q-value (lower is better)",
        "higher_better": False,
    },
    "global_qvalue": {
        "ontology_name": "global q-value (DIA-NN:Global.Q.Value)",
        "ontology_accession": "MS:1002350",
        "ontology_source": "MS",
        "description": "Global q-value at experiment level (DIA-NN Global.Q.Value, lower is better)",
        "higher_better": False,
    },
    "pg_global_qvalue": {
        "ontology_name": "protein group global FDR (DIA-NN:Global.PG.Q.Value)",
        "ontology_accession": "MS:1001214",
        "ontology_source": "MS",
        "description": "Protein group global q-value (DIA-NN Global.PG.Q.Value, lower is better)",
        "higher_better": False,
    },
    # --- Phospho site localization scores ---
    "phosphors_site_probability": {
        "ontology_name": "PhosphoRS site probability",
        "ontology_accession": "MS:1001969",
        "ontology_source": "MS",
        "description": "PhosphoRS phosphorylation site probability (higher is better)",
        "higher_better": True,
    },
    "ptmrs_site_probability": {
        "ontology_name": "ptmRS site probability",
        "ontology_accession": "MS:1002536",
        "ontology_source": "MS",
        "description": "ptmRS PTM site probability (higher is better)",
        "higher_better": True,
    },
    "luciphor_site_probability": {
        "ontology_name": "Luciphor site probability",
        "ontology_accession": "MS:1003291",
        "ontology_source": "MS",
        "description": "Luciphor PTM site localization probability (higher is better)",
        "higher_better": True,
    },
    "d_score_site_probability": {
        "ontology_name": "D-Score site probability",
        "ontology_accession": "MS:1002536",
        "ontology_source": "MS",
        "description": "D-Score phospho site localization probability (higher is better)",
        "higher_better": True,
    },
    "phospho_sty_probability": {
        "ontology_name": "MaxQuant Phospho (STY) probability",
        "ontology_accession": "MS:1001982",
        "ontology_source": "MS",
        "description": "MaxQuant per-site phosphorylation probability (higher is better)",
        "higher_better": True,
    },
    # --- MaxQuant Andromeda scores ---
    "andromeda_score": {
        "ontology_name": "Andromeda score",
        "ontology_accession": "MS:1002338",
        "ontology_source": "MS",
        "description": "Andromeda search engine score (higher is better)",
        "higher_better": True,
    },
    "andromeda_delta_score": {
        "ontology_name": "Andromeda delta score",
        "ontology_accession": "MS:1003433",
        "ontology_source": "MS",
        "description": "Andromeda delta score between best and second-best match (higher is better)",
        "higher_better": True,
    },
    # --- Intensity terms ---
    "lfq": {
        "ontology_name": "MaxLFQ intensity",
        "ontology_accession": None,
        "ontology_source": None,
        "description": "Label-free quantification intensity (MaxLFQ algorithm)",
        "higher_better": True,
    },
}


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
    try:
        ontology = _get_ontology()
    except Exception:
        return None

    term = ontology.lookup_normalized(normalized_name)
    if term is None:
        return None

    return {
        "ontology_name": term["name"],
        "ontology_accession": term["accession"],
        "ontology_source": term["source"],
        "description": term["definition"],
        "higher_better": (term["higher_better"] if term["higher_better"] is not None else True),
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
    entry = _BUILTIN_SCORES.get(normalized_name)
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


def is_higher_better(normalized_name: str) -> bool:
    """Return whether the score with *normalized_name* is higher-is-better.

    Checks built-in cache, then OBO registry. Falls back to ``True``
    (higher is better) for completely unknown scores.
    """
    info = lookup_score(normalized_name)
    if info is not None:
        return info["higher_better"]
    return True


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
    try:
        ontology_version = _get_ontology().version
    except Exception:
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
                    "ontology_version": ontology_version,
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


def modification_ontology_entries(
    modifications_meta: dict,
    view: str = "psm",
) -> list[dict]:
    """Generate ontology.parquet entries for discovered modifications.

    Args:
        modifications_meta: Dict keyed by accession, values
            ``(name, list[site], list[position])``.
        view: The QPX view name.

    Returns:
        List of dicts matching the ``OntologySchema``.
    """
    entries: list[dict] = []
    for accession, (name, _sites, _positions) in modifications_meta.items():
        if not accession:
            continue
        source = "UNIMOD" if accession.startswith("UNIMOD:") else "MOD"
        entries.append(
            {
                "field_name": name,
                "ontology_name": name,
                "ontology_accession": accession,
                "ontology_source": source,
                "ontology_version": None,
                "view": view,
                "description": f"Modification: {name} ({accession})",
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
    try:
        ontology_version = _get_ontology().version
    except Exception:
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
                "ontology_version": ontology_version,
                "view": view,
                "description": info["description"],
                "source_column_name": None,
                "source_tool": None,
            }
        )
    return entries
