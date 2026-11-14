"""Expression and differential modalities of MuData, read from the dataset's ae and de AnnData files."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def _build_expression_adata(ds_path: Path, file_prefix: str | None = None):
    """Auto-detect and load the absolute-expression AnnData: ``.ae``, or the older ``.pe``, as h5ad or zarr.

    Returns None if not found.
    """
    import anndata as ad

    stem = file_prefix or "*"
    for view in ("ae", "pe"):
        candidates = list(ds_path.glob(f"{stem}.{view}.h5ad"))
        if candidates:
            return ad.read_h5ad(candidates[0])
        candidates = list(ds_path.glob(f"{stem}.{view}.zarr"))
        if candidates:
            return ad.read_zarr(candidates[0])

    return None


def _build_differential_adata(ds_path: Path, file_prefix: str | None = None):
    """Auto-detect and load ``.de.h5ad`` or ``.de.zarr``.

    Restructures ``uns["de_results"]`` into matrix form:
    obs = contrasts, var = proteins, X = log2FC,
    layers = {pvals_adj, scores, pvals, se}.

    Returns None if not found.
    """
    import anndata as ad

    # Try h5ad first
    h5ad_pattern = f"{file_prefix}.de.h5ad" if file_prefix else "*.de.h5ad"
    candidates = list(ds_path.glob(h5ad_pattern))
    if not candidates:
        zarr_pattern = f"{file_prefix}.de.zarr" if file_prefix else "*.de.zarr"
        candidates = list(ds_path.glob(zarr_pattern))
        if not candidates:
            return None
        adata = ad.read_zarr(candidates[0])
    else:
        adata = ad.read_h5ad(candidates[0])

    # If it has de_results in uns, restructure into matrix form
    if "de_results" in adata.uns:
        de = adata.uns["de_results"]
        if isinstance(de, dict):
            de = _long_de_results(de)
        if not de.empty:
            return _reshape_de_results(de)

    # Already in matrix form or no de_results — return as-is
    return adata


def _long_de_results(contrasts: dict) -> pd.DataFrame:
    """Flatten the specification's ``{contrast: {names, logfoldchanges, ...}}`` layout into long form."""
    frames = []
    for contrast, fields in contrasts.items():
        frame = pd.DataFrame(
            {"contrast": contrast, "protein": list(fields["names"]), "log2FC": fields.get("logfoldchanges", np.nan)}
        )
        for column in ("pvals_adj", "scores", "pvals", "se"):
            if column in fields:
                frame[column] = fields[column]
        frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _reshape_de_results(de: pd.DataFrame):
    """Reshape long-form DE results into an AnnData matrix.

    Expects columns: contrast, protein, log2FC and optionally pvals_adj, scores,
    pvals, se. A protein a contrast did not test is NaN, because 0 would read as
    no change and a p-value of 0.
    """
    import anndata as ad

    # Identify contrast and protein columns
    contrast_col = "contrast" if "contrast" in de.columns else de.columns[0]
    protein_col = "protein" if "protein" in de.columns else de.columns[1]

    if "log2FC" not in de.columns:
        raise ValueError(f"DE results missing required 'log2FC' column. Found columns: {list(de.columns)}")

    contrasts = pd.Index(sorted(de[contrast_col].unique()), name="contrast")
    proteins = pd.Index(sorted(de[protein_col].unique()), name="protein")

    def matrix(column: str) -> np.ndarray:
        return (
            de.pivot_table(index=contrast_col, columns=protein_col, values=column, aggfunc="first", dropna=False)
            .reindex(index=contrasts, columns=proteins)
            .to_numpy(dtype="float64")
        )

    return ad.AnnData(
        X=matrix("log2FC"),
        obs=pd.DataFrame(index=contrasts),
        var=pd.DataFrame(index=proteins),
        layers={column: matrix(column) for column in ("pvals_adj", "scores", "pvals", "se") if column in de.columns},
    )
