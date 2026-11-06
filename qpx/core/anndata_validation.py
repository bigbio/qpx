"""Checks of the AnnData views (``.ae.h5ad``, ``.de.h5ad``), which have no Parquet schema to validate against."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from qpx.core.data.schema import ValidationIssue, ValidationResult

_FILE_TYPES = {"ae": "absolute_expression", "de": "differential_expression"}
_REQUIRED_DE_ARRAYS = ("names", "pvals", "pvals_adj")
# Scalars of a contrast, not arrays parallel to ``names``.
_DE_SCALARS = ("condition_test", "condition_reference", "contrast_definition", "statistic_type")
# These tests compare several coefficients at once, so they have no single fold change.
_OMNIBUS_TESTS = ("F", "LRT")
_ENTITY_TYPES = ("protein_group", "gene_group", "peptidoform", "site", "proteoform")
_AE_RECORDED_KEYS = ("log_base", "var_key", "aggregation")
_VAR_KEYS = ("pg_accessions", "anchor_protein")
_STANDARDS_COLUMNS = ("protein", "amount", "unit")
_LOG_TOLERANCE = 1e-5


def _duplicate_index_problems(adata) -> list[str]:
    return [f"{slot} index has duplicates" for slot in ("obs", "var") if not getattr(adata, slot).index.is_unique]


def _dense(matrix) -> np.ndarray:
    return np.asarray(matrix.todense() if hasattr(matrix, "todense") else matrix, dtype="float64")


def _log_base_problems(adata) -> list[str]:
    """``X`` must be the declared logarithm of ``ibaq_raw``."""
    base = adata.uns.get("log_base")
    if base is None or "ibaq_raw" not in adata.layers:
        return []
    raw = _dense(adata.layers["ibaq_raw"])
    values = _dense(adata.X)
    measured = np.isfinite(values) & (raw > 0)
    if measured.any() and np.abs(values[measured] - np.log(raw[measured]) / np.log(float(base))).max() > _LOG_TOLERANCE:
        return [f"X is not the base {base} logarithm of layer ibaq_raw"]
    return []


def _ae_problems(adata) -> list[str]:
    problems = _duplicate_index_problems(adata)
    if adata.n_obs == 0 or adata.n_vars == 0:
        problems.append(f"matrix is empty ({adata.n_obs} samples x {adata.n_vars} proteins)")
    problems.extend(
        f"layer {name} has negative values" for name, matrix in adata.layers.items() if np.nanmin(_dense(matrix), initial=0.0) < 0
    )
    if np.isinf(_dense(adata.X)).any():
        problems.append("X has infinite values; a protein not detected in a sample is NaN")
    problems.extend(_log_base_problems(adata))
    if adata.uns.get("var_key", _VAR_KEYS[0]) not in _VAR_KEYS:
        problems.append(f"uns['var_key'] is not one of {', '.join(_VAR_KEYS)}")
    standards = adata.uns.get("standards")
    if standards is not None and not set(_STANDARDS_COLUMNS) <= set(getattr(standards, "columns", ())):
        problems.append(f"uns['standards'] is not a table with the columns {', '.join(_STANDARDS_COLUMNS)}")
    return problems


def _ae_warnings(adata) -> list[str]:
    return [f"uns has no {key}" for key in _AE_RECORDED_KEYS if key not in adata.uns]


def _contrast_problems(contrast: str, fields: dict) -> list[str]:
    required = list(_REQUIRED_DE_ARRAYS)
    if fields.get("statistic_type") not in _OMNIBUS_TESTS:
        required.append("logfoldchanges")
    if "contrast_definition" not in fields:
        required.extend(("condition_test", "condition_reference"))
    problems = [f"contrast {contrast!r} lacks {name}" for name in required if name not in fields]
    arrays = {name: np.asarray(fields[name]) for name in fields if name not in _DE_SCALARS}
    if len({values.shape[0] for values in arrays.values()}) > 1:
        problems.append(f"contrast {contrast!r} has arrays of different lengths")
    for name in ("pvals", "pvals_adj"):
        if name in arrays:
            values = arrays[name].astype("float64")
            if ((values < 0) | (values > 1)).any():
                problems.append(f"contrast {contrast!r} has {name} outside [0, 1]")
    return problems


def _de_problems(adata) -> list[str]:
    problems = _duplicate_index_problems(adata)
    results = adata.uns.get("de_results")
    if results is None:
        return [*problems, "uns has no de_results"]
    if not isinstance(results, dict):
        return [*problems, f"uns['de_results'] is a {type(results).__name__}; the specification keeps a dict keyed by contrast"]
    for contrast, fields in results.items():
        problems.extend(_contrast_problems(contrast, fields))
    if adata.uns.get("entity_type", _ENTITY_TYPES[0]) not in _ENTITY_TYPES:
        problems.append(f"uns['entity_type'] is not one of {', '.join(_ENTITY_TYPES)}")
    return problems


def validate_anndata_view(path: Path, view: str) -> ValidationResult:
    """Check an AnnData view against its specification page; a file that cannot be read is an error."""
    import anndata

    issues = []

    def issue(check: str, severity: str, message: str) -> None:
        issues.append(ValidationIssue(structure=view, check=check, severity=severity, column=None, message=message))

    try:
        adata = anndata.read_h5ad(path)
    except OSError as exc:
        issue("unreadable", "error", f"{path.name} could not be read: {exc}")
        return ValidationResult(structure=view, issues=issues)
    for problem in _ae_problems(adata) if view == "ae" else _de_problems(adata):
        issue("anndata_layout", "error", f"{path.name}: {problem}")
    for problem in _ae_warnings(adata) if view == "ae" else ():
        issue("anndata_metadata", "warning", f"{path.name}: {problem}")
    if adata.uns.get("file_type") != _FILE_TYPES[view]:
        issue("anndata_file_type", "warning", f"{path.name}: uns['file_type'] is not {_FILE_TYPES[view]!r}")
    return ValidationResult(structure=view, issues=issues)
