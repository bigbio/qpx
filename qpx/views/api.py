"""Computed views — cross-structure projections using shared DuckDB connection."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from qpx.core.convert import QueryResult
from qpx.core.sample_resolution import RUN_SAMPLES_SQL, feature_samples_sql
from qpx.core.sql import sql_build
from qpx.views.base import BaseView

if TYPE_CHECKING:
    from qpx.dataset import Dataset


# The summary views count targets only, and count protein groups by membership:
# an anchor protein can lead several groups.
_FEATURE_GROUP_SQL = "list_sort(list_transform(f.pg_accessions, lambda x: x.accession))"

_RUN_COUNTS_SQL = f"""
        SELECT f.run_file_name,
               COUNT(DISTINCT {_FEATURE_GROUP_SQL}) AS n_proteins,
               COUNT(DISTINCT f.sequence) AS n_peptides,
               COUNT(*) AS n_features
        FROM feature f
        WHERE NOT f.is_decoy
        GROUP BY f.run_file_name
        """


class ProteinView(BaseView):
    """Protein intensity per sample — joins PG + Run."""

    def intensity(self, q_value_threshold: float | None = None) -> QueryResult:
        """Target group intensity per sample and label.

        ``protein_group`` (the sorted membership joined by ``;``) identifies the
        analyte; ``protein_accession`` is the group's anchor, a display label that
        can lead several groups. *q_value_threshold* keeps the groups whose
        ``global_qvalue`` is at most the threshold or null; by default no
        q-value filter is applied.
        """
        self._require("pg", "run")
        threshold = (
            ""
            if q_value_threshold is None
            else f"AND (p.global_qvalue IS NULL OR p.global_qvalue <= {float(q_value_threshold)!r})"
        )
        # Any member run of grouped_runs resolves the label to the unit's sample.
        sql = """
        WITH run_samples AS ($run_samples)
        SELECT p.anchor_protein AS protein_accession,
               array_to_string(list_sort(p.pg_accessions), ';') AS protein_group,
               rs.sample_accession,
               p.label,
               p.intensity,
               p.global_qvalue,
               p.gg_names AS gene_names
        FROM pg p
        JOIN run_samples rs ON rs.run_file_name = p.grouped_runs[1] AND rs.label = p.label
        WHERE NOT p.is_decoy $threshold
        """
        sql = sql_build(sql, run_samples=RUN_SAMPLES_SQL, threshold=threshold)
        return QueryResult(self._engine.connection, sql)


class PeptideView(BaseView):
    """Peptide intensity per sample — joins Feature + Run."""

    def intensity(self) -> QueryResult:
        """Target peptide intensity per sample and label."""
        self._require("feature", "run")
        sql = sql_build(
            """
        SELECT sequence, peptidoform, charge, anchor_protein, sample_accession, label, intensity
        FROM ($feature_samples)
        WHERE NOT is_decoy
        """,
            feature_samples=feature_samples_sql("feature", ("sequence", "peptidoform", "charge", "anchor_protein", "is_decoy")),
        )
        return QueryResult(self._engine.connection, sql)


class IdentificationSummaryView(BaseView):
    """Per-run summary: protein groups, peptides, features."""

    def summary(self) -> QueryResult:
        """Target protein groups, peptide sequences and features per run."""
        self._require("feature")
        return QueryResult(self._engine.connection, _RUN_COUNTS_SQL)

    def plot(self):
        """Bar chart of identifications per run."""
        from qpx.views.plotting import grouped_bar_chart

        df = self.summary().to_df()
        return grouped_bar_chart(
            df,
            x="run_file_name",
            y_cols=["n_proteins", "n_peptides"],
            title="Identifications per Run",
            xlabel="Run",
            ylabel="Count",
        )


class RunSummaryView(BaseView):
    """Per-run statistics with intensity distribution."""

    def summary(self):
        """Target peptide sequences, protein groups and features per run (cached)."""
        self._require("feature")
        return self._execute_cached("run_summary", _RUN_COUNTS_SQL)

    def plot(self):
        """Bar chart of run-level statistics."""
        from qpx.views.plotting import grouped_bar_chart

        df = self.summary().to_df()
        return grouped_bar_chart(
            df,
            x="run_file_name",
            y_cols=["n_peptides", "n_proteins", "n_features"],
            title="Run Summary",
            xlabel="Run",
            ylabel="Count",
        )


class ModificationView(BaseView):
    """Modification frequency across PSMs."""

    def frequency(self):
        """Target PSMs and peptide sequences per modification."""
        self._require("psm")
        sql = """
        SELECT m.name AS modification_name,
               m.accession AS modification_accession,
               COUNT(*) AS n_psms,
               COUNT(DISTINCT p.sequence) AS n_peptides
        FROM psm p, UNNEST(p.modifications) AS _t(m)
        WHERE p.modifications IS NOT NULL AND NOT p.is_decoy
        GROUP BY m.name, m.accession
        ORDER BY n_psms DESC
        """
        return self._execute_cached("mod_frequency", sql)

    def plot(self, top_n=20):
        """Bar chart of top modifications by PSM count."""
        from qpx.views.plotting import bar_chart

        df = self.frequency().to_df().head(top_n)
        return bar_chart(
            df,
            x="modification_name",
            y="n_psms",
            title=f"Top {top_n} Modifications by PSM Count",
            xlabel="Modification",
            ylabel="PSM Count",
        )


class QualityControlView(BaseView):
    """Dataset-wide QC metrics."""

    def metrics(self):
        """Target peptide sequences, protein groups, runs and features of the dataset (cached)."""
        self._require("feature")
        sql = f"""
        SELECT COUNT(DISTINCT f.sequence) AS unique_peptides,
               COUNT(DISTINCT {_FEATURE_GROUP_SQL}) AS unique_proteins,
               COUNT(DISTINCT f.run_file_name) AS n_runs,
               COUNT(*) AS total_features
        FROM feature f
        WHERE f.is_decoy = false
        """
        return self._execute_cached("qc_metrics", sql)

    def plot(self):
        """Bar chart of dataset QC summary metrics."""
        import pandas as pd

        from qpx.views.plotting import bar_chart

        df = self.metrics().to_df()
        summary = pd.DataFrame(
            {
                "metric": [
                    "Unique Peptides",
                    "Unique Proteins",
                    "Runs",
                    "Total Features",
                ],
                "value": [
                    df["unique_peptides"].iloc[0],
                    df["unique_proteins"].iloc[0],
                    df["n_runs"].iloc[0],
                    df["total_features"].iloc[0],
                ],
            }
        )
        return bar_chart(
            summary,
            x="metric",
            y="value",
            title="Dataset QC Summary",
            rotation=0,
        )


class SampleSummaryView(BaseView):
    """Per-sample and per-run peptide/protein counts."""

    def peptides_per_run(self):
        """Count distinct peptide sequences per run."""
        self._require("feature")
        sql = """
        SELECT f.run_file_name,
               COUNT(DISTINCT f.sequence) AS n_peptides
        FROM feature f
        WHERE f.is_decoy = false
        GROUP BY f.run_file_name
        ORDER BY f.run_file_name
        """
        return self._execute_cached("peptides_per_run", sql)

    def peptides_per_sample(self):
        """Count distinct peptide sequences per sample, resolved through the label of each intensity."""
        self._require("feature", "run")
        sql = sql_build(
            """
        SELECT sample_accession, COUNT(DISTINCT sequence) AS n_peptides
        FROM ($feature_samples)
        WHERE is_decoy = false
        GROUP BY sample_accession
        ORDER BY sample_accession
        """,
            feature_samples=feature_samples_sql("feature", ("sequence", "is_decoy")),
        )
        return self._execute_cached("peptides_per_sample", sql)

    def proteins_per_sample(self):
        """Count distinct protein groups per sample.

        A row resolves to its sample through its label; an identification-only
        row (no label) resolves in a single-sample run.
        """
        self._require("pg", "run")
        sql = sql_build(
            """
        WITH run_samples AS ($run_samples)
        SELECT rs.sample_accession, COUNT(DISTINCT list_sort(p.pg_accessions)) AS n_proteins
        FROM pg p
        JOIN run_samples rs
          ON rs.run_file_name = p.grouped_runs[1]
         AND (rs.label = p.label OR (p.label IS NULL AND rs.single_sample))
        WHERE p.is_decoy = false
        GROUP BY rs.sample_accession
        ORDER BY rs.sample_accession
        """,
            run_samples=RUN_SAMPLES_SQL,
        )
        return self._execute_cached("proteins_per_sample", sql)

    def plot(self, metric: str = "peptides_per_run"):
        """Bar chart for the given metric.

        Parameters
        ----------
        metric : str
            One of ``"peptides_per_run"``, ``"peptides_per_sample"``,
            ``"proteins_per_sample"``.
        """
        from qpx.views.plotting import bar_chart

        dispatch = {
            "peptides_per_run": (
                self.peptides_per_run,
                "run_file_name",
                "n_peptides",
                "Peptides per Run",
                "Run",
                "Distinct Peptides",
            ),
            "peptides_per_sample": (
                self.peptides_per_sample,
                "sample_accession",
                "n_peptides",
                "Peptides per Sample",
                "Sample",
                "Distinct Peptides",
            ),
            "proteins_per_sample": (
                self.proteins_per_sample,
                "sample_accession",
                "n_proteins",
                "Proteins per Sample",
                "Sample",
                "Distinct Proteins",
            ),
        }
        if metric not in dispatch:
            raise ValueError(f"Unknown metric '{metric}'. Choose from: {list(dispatch.keys())}")
        fn, x, y, title, xlabel, ylabel = dispatch[metric]
        df = fn().to_df()
        return bar_chart(df, x=x, y=y, title=title, xlabel=xlabel, ylabel=ylabel)


class AbsoluteExpressionView:
    """AnnData-based view for absolute expression analysis (PCA, etc.).

    Lazily loads the ``.ae.h5ad`` file from the dataset directory.
    """

    def __init__(self, dataset: Dataset):
        self._dataset = dataset
        self._adata = None

    def _load_anndata(self):
        """Find and load the .ae.h5ad file from the dataset directory."""
        if self._adata is not None:
            return self._adata

        try:
            import anndata as ad
        except ImportError as exc:
            raise ImportError(
                "anndata is required for AbsoluteExpressionView; reinstall qpx (anndata is a core dependency)"
            ) from exc

        ds_path = Path(self._dataset.path)
        # Look for <prefix>.ae.h5ad or any .ae.h5ad
        candidates = list(ds_path.glob("*.ae.h5ad"))
        if not candidates:
            raise FileNotFoundError(
                f"No .ae.h5ad file found in {ds_path}. Generate one with save_anndata(view='ae') or a downstream tool."
            )
        self._adata = ad.read_h5ad(candidates[0])
        return self._adata

    def pca(
        self,
        color_by: str | None = None,
        layer: str | None = None,
        n_components: int = 2,
        fillna: str = "zero",
    ):
        """PCA scatter plot of samples.

        Parameters
        ----------
        color_by : str, optional
            obs column to color samples by (e.g. "organism", "disease").
        layer : str, optional
            AnnData layer to use. If None, uses ``.X``.
        n_components : int
            Number of PCA components (default 2).
        fillna : str
            How to handle NaN: ``"zero"`` or ``"median"``.
        """
        import numpy as np
        import pandas as pd
        from sklearn.decomposition import PCA

        from qpx.views.plotting import scatter_plot

        adata = self._load_anndata()
        X = adata.layers[layer] if layer else adata.X
        X = np.array(X, dtype=np.float64)

        # Handle NaN
        if fillna == "median":
            col_medians = np.nanmedian(X, axis=0)
            inds = np.where(np.isnan(X))
            X[inds] = np.take(col_medians, inds[1])
        # Zero-fill any remaining NaN (all-NaN columns after median, or default)
        X = np.nan_to_num(X, nan=0.0)

        pca_model = PCA(n_components=n_components)
        coords = pca_model.fit_transform(X)

        df = pd.DataFrame(
            coords,
            columns=[f"PC{i + 1}" for i in range(n_components)],
            index=adata.obs.index,
        )
        if color_by and color_by in adata.obs.columns:
            df[color_by] = adata.obs[color_by].values

        var_explained = pca_model.explained_variance_ratio_ * 100
        return scatter_plot(
            df,
            x="PC1",
            y="PC2",
            color_by=color_by if color_by and color_by in adata.obs.columns else None,
            title="PCA — Samples",
            xlabel=f"PC1 ({var_explained[0]:.1f}%)",
            ylabel=f"PC2 ({var_explained[1]:.1f}%)",
        )
