"""Checks of the AnnData views (``.ae.h5ad``, ``.de.h5ad``), which have no Parquet schema to validate against."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from qpx.core.data.schema import ValidationIssue, ValidationResult

_FILE_TYPES = {"ae": "absolute_expression", "de": "differential_expression"}
_REQUIRED_DE_ARRAYS = ("names", "logfoldchanges", "pvals", "pvals_adj")
_REQUIRED_DE_CONDITIONS = ("condition_test", "condition_reference")


def _duplicate_index_problems(adata) -> list[str]:
    return [f"{slot} index has duplicates" for slot in ("obs", "var") if not getattr(adata, slot).index.is_unique]


def _ae_problems(adata) -> list[str]:
    problems = _duplicate_index_problems(adata)
    if adata.n_obs == 0 or adata.n_vars == 0:
        problems.append(f"matrix is empty ({adata.n_obs} samples x {adata.n_vars} proteins)")
    if "ibaq_raw" in adata.layers and np.nanmin(np.asarray(adata.layers["ibaq_raw"], dtype="float64"), initial=0.0) < 0:
        problems.append("layer ibaq_raw has negative values")
    return problems


def _contrast_problems(contrast: str, fields: dict) -> list[str]:
    problems = [
        f"contrast {contrast!r} lacks {name}" for name in (*_REQUIRED_DE_ARRAYS, *_REQUIRED_DE_CONDITIONS) if name not in fields
    ]
    arrays = {name: np.asarray(fields[name]) for name in fields if name not in _REQUIRED_DE_CONDITIONS}
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
    if adata.uns.get("file_type") != _FILE_TYPES[view]:
        issue("anndata_file_type", "warning", f"{path.name}: uns['file_type'] is not {_FILE_TYPES[view]!r}")
    return ValidationResult(structure=view, issues=issues)
