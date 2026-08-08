"""Computed views — cross-structure projections using shared DuckDB connection."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from qpx.core.convert import QueryResult
from qpx.views.base import BaseView

if TYPE_CHECKING:
    from qpx.dataset import Dataset


class ProteinView(BaseView):
    """Protein intensity per sample — joins PG + Run."""

    def intensity(self, q_value_threshold: float = 0.01) -> QueryResult:
        q_value_threshold = float(q_value_threshold)
        sql = """
        WITH numbered_pg AS MATERIALIZED (
            SELECT ROW_NUMBER() OVER () AS pg_row_id, *
            FROM pg
        ),
        pg_samples AS (
            SELECT DISTINCT p.pg_row_id, rs.sample_accession, rs.label
            FROM numbered_pg p
            CROSS JOIN UNNEST(list_distinct(p.grouped_runs)) AS _g(run_file_name)
            JOIN run r USING (run_file_name)
            CROSS JOIN UNNEST(r.samples) AS _s(rs)
        )
        SELECT p.anchor_protein AS protein_accession,
               ps.sample_accession,
               p.label,
               p.intensity,
               p.global_qvalue,
               p.gg_names AS gene_names
        FROM numbered_pg p
        JOIN pg_samples ps USING (pg_row_id)
        WHERE p.label = ps.label
          AND (p.global_qvalue IS NULL OR p.global_qvalue <= $1)
        """
        return QueryResult(self._engine.execute(sql, [q_value_threshold]))


class PeptideView(BaseView):
    """Peptide intensity per sample — joins Feature + Run."""

    def intensity(self) -> QueryResult:
        # The label carried by a feature intensity is not always the SDRF channel
        # label carried by run.samples. For label-free (LFQ) runs the converters
        # emit a generic feature label (e.g. DIA-NN's ``'raw'``) while run.samples
        # carries the canonical ``'LFQ'`` channel, so an exact ``i.label = rs.label``
        # join returns 0 rows. A label-free run has exactly one sample, so when the
        # run is single-sample we attribute the intensity to that lone sample
        # regardless of label; multiplexed (TMT/iTRAQ) runs still require an exact
        # label match to disambiguate channels.
        sql = """
        SELECT f.sequence, f.peptidoform, f.charge, f.anchor_protein,
               rs.sample_accession, rs.label AS label, i.intensity
        FROM feature f
        JOIN run r ON f.run_file_name = r.run_file_name,
             UNNEST(r.samples) AS _t1(rs),
             UNNEST(f.intensities) AS _t2(i)
        WHERE (i.label = rs.label OR len(r.samples) = 1)
        """
        return QueryResult(self._engine.execute(sql))


class IdentificationSummaryView(BaseView):
    """Per-run summary: proteins, peptides, features."""

    def summary(self) -> QueryResult:
        sql = """
        SELECT f.run_file_name,
               COUNT(DISTINCT f.anchor_protein) AS n_proteins,
               COUNT(DISTINCT f.sequence) AS n_peptides,
               COUNT(*) AS n_features
        FROM feature f
        GROUP BY f.run_file_name
        """
        return QueryResult(self._engine.execute(sql))

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
        sql = """
        SELECT f.run_file_name,
               COUNT(DISTINCT f.sequence) AS n_peptides,
               COUNT(DISTINCT f.anchor_protein) AS n_proteins,
               COUNT(*) AS n_features
        FROM feature f
        GROUP BY f.run_file_name
        """
        return self._execute_cached("run_summary", sql)

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
        sql = """
        SELECT m.name AS modification_name,
               m.accession AS modification_accession,
               COUNT(*) AS n_psms,
               COUNT(DISTINCT p.sequence) AS n_peptides
        FROM psm p, UNNEST(p.modifications) AS _t(m)
        WHERE p.modifications IS NOT NULL
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
        sql = """
        SELECT COUNT(DISTINCT f.sequence) AS unique_peptides,
               COUNT(DISTINCT f.anchor_protein) AS unique_proteins,
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
        """Count distinct peptide sequences per sample (joins Feature + Run)."""
        sql = """
        SELECT rs.sample_accession,
               COUNT(DISTINCT f.sequence) AS n_peptides
        FROM feature f
        JOIN run r ON f.run_file_name = r.run_file_name,
             UNNEST(r.samples) AS _t(rs)
        WHERE f.is_decoy = false
        GROUP BY rs.sample_accession
        ORDER BY rs.sample_accession
        """
        return self._execute_cached("peptides_per_sample", sql)

    def proteins_per_sample(self):
        """Count distinct proteins per sample (joins PG + Run)."""
        # Pre-explode grouped_runs so the run join is a hash-join on
        # run_file_name instead of a nested-loop list_contains predicate.
        sql = """
        WITH pg_runs AS (
            SELECT p.anchor_protein, p.is_decoy, gr AS run_file_name
            FROM pg p, UNNEST(p.grouped_runs) AS _g(gr)
        )
        SELECT rs.sample_accession,
               COUNT(DISTINCT pg.anchor_protein) AS n_proteins
        FROM pg_runs pg
        JOIN run r ON pg.run_file_name = r.run_file_name,
             UNNEST(r.samples) AS _t(rs)
        WHERE pg.is_decoy = false
        GROUP BY rs.sample_accession
        ORDER BY rs.sample_accession
        """
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
