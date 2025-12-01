"""
Report generation module for QPX projects.

This module provides functionality to generate comprehensive statistical reports for QPX
projects, including analysis of ibaqpy pipeline results with detailed statistics about samples,
proteins, peptides, and intensity distributions.
"""

import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple
from multiprocessing import Pool

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import umap
import pyarrow.parquet as pq
import psutil


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def _process_chunk_for_proteins(args):
    """Extract protein IDs from a file partition"""
    file_path, start_row, num_rows = args

    try:
        parquet_file = pq.ParquetFile(file_path)

        proteins = set()
        batch_size = 50000

        current_row = 0
        for batch in parquet_file.iter_batches(
            batch_size=batch_size, columns=["pg_accessions"]
        ):
            batch_end = current_row + len(batch)

            if batch_end > start_row and current_row < start_row + num_rows:
                batch_start_idx = max(0, start_row - current_row)
                batch_end_idx = min(len(batch), start_row + num_rows - current_row)

                df = batch.to_pandas().iloc[batch_start_idx:batch_end_idx]

                for pg_acc in df["pg_accessions"]:
                    if isinstance(pg_acc, (list, tuple)) and len(pg_acc) > 0:
                        proteins.add(pg_acc[0])
                    else:
                        proteins.add(str(pg_acc))

            current_row = batch_end

            if current_row >= start_row + num_rows:
                break

        return proteins
    except Exception as e:
        logger.warning(f"Error in chunk {start_row}-{start_row+num_rows}: {e}")
        return set()


def _process_chunk_for_filter(args):
    """Calculate unique peptides per protein from a file partition"""
    file_path, start_row, num_rows = args

    try:
        parquet_file = pq.ParquetFile(file_path)

        unique_counts = {}
        batch_size = 50000

        current_row = 0
        for batch in parquet_file.iter_batches(
            batch_size=batch_size, columns=["pg_accessions", "unique", "sequence"]
        ):
            batch_end = current_row + len(batch)

            if batch_end > start_row and current_row < start_row + num_rows:
                batch_start_idx = max(0, start_row - current_row)
                batch_end_idx = min(len(batch), start_row + num_rows - current_row)

                df = batch.to_pandas().iloc[batch_start_idx:batch_end_idx]

                df["_protein"] = df["pg_accessions"].apply(
                    lambda x: x[0] if isinstance(x, list) and len(x) > 0 else str(x)
                )

                unique_df = df[df["unique"] == 1][["_protein", "sequence"]]

                for protein, group in unique_df.groupby("_protein"):
                    sequences = set(group["sequence"].unique())
                    if protein in unique_counts:
                        unique_counts[protein].update(sequences)
                    else:
                        unique_counts[protein] = sequences

            current_row = batch_end

            if current_row >= start_row + num_rows:
                break

        return unique_counts
    except Exception as e:
        logger.warning(f"Error in chunk {start_row}-{start_row+num_rows}: {e}")
        return {}


def _normalize_protein_column(df_chunk):
    """Normalize protein column from pg_accessions."""
    if df_chunk["pg_accessions"].apply(lambda x: isinstance(x, list)).any():
        df_chunk["protein"] = df_chunk["pg_accessions"].apply(
            lambda x: x[0] if isinstance(x, list) and len(x) > 0 else str(x)
        )
    else:
        df_chunk["protein"] = df_chunk["pg_accessions"].astype(str)
    return df_chunk


def _aggregate_sample_intensities(df_chunk, sample_data):
    """Aggregate protein intensities for each sample."""
    for sample in df_chunk["sample_accession"].unique():
        if sample not in sample_data:
            sample_data[sample] = {}

        sample_rows = df_chunk[df_chunk["sample_accession"] == sample]
        for _, row in sample_rows.iterrows():
            protein = row["protein"]
            intensity = row["intensity"]
            if intensity > 0:
                if protein in sample_data[sample]:
                    sample_data[sample][protein] = max(
                        sample_data[sample][protein], intensity
                    )
                else:
                    sample_data[sample][protein] = intensity


def _filter_chunk_by_range(df_chunk, current_row, start_row, batch_end, end_row):
    """Filter chunk to include only rows within the specified range."""
    if current_row < start_row:
        skip_rows = start_row - current_row
        df_chunk = df_chunk.iloc[skip_rows:]

    if batch_end > end_row:
        keep_rows = end_row - max(current_row, start_row)
        df_chunk = df_chunk.iloc[:keep_rows]

    return df_chunk


def _process_parquet_partition(target_samples, file_path, start_row, end_row):
    """Process a parquet file partition."""
    parquet_file = pq.ParquetFile(file_path)
    total_rows = parquet_file.metadata.num_rows

    if end_row is None:
        end_row = total_rows

    sample_data = {}
    current_row = 0

    for batch in parquet_file.iter_batches(
        batch_size=500000,
        columns=["sample_accession", "pg_accessions", "intensity"],
    ):
        batch_size = len(batch)
        batch_end = current_row + batch_size

        if batch_end < start_row:
            current_row = batch_end
            continue

        if current_row >= end_row:
            break

        df_chunk = batch.to_pandas()
        df_chunk = _filter_chunk_by_range(
            df_chunk, current_row, start_row, batch_end, end_row
        )
        df_chunk = df_chunk[df_chunk["sample_accession"].isin(target_samples)]

        if len(df_chunk) == 0:
            current_row = batch_end
            continue

        df_chunk = _normalize_protein_column(df_chunk)
        _aggregate_sample_intensities(df_chunk, sample_data)
        current_row = batch_end

    return sample_data


def _process_csv_partition(target_samples, file_path, start_row, end_row):
    """Process a CSV file partition."""
    chunk_iter = pd.read_csv(
        file_path,
        sep="\t",
        usecols=["SampleID", "ProteinName", "Ibaq"],
        chunksize=100000,
        skiprows=range(1, start_row) if start_row > 0 else None,
        nrows=(end_row - start_row) if end_row else None,
    )

    sample_data = {}
    for df_chunk in chunk_iter:
        df_chunk = df_chunk.rename(
            columns={
                "SampleID": "sample_accession",
                "ProteinName": "protein",
                "Ibaq": "intensity",
            }
        )
        df_chunk = df_chunk[df_chunk["sample_accession"].isin(target_samples)]

        if len(df_chunk) == 0:
            continue

        _aggregate_sample_intensities(df_chunk, sample_data)

    return sample_data


def _process_file_partition(args):
    """Process a partition of the file with specific row ranges."""
    target_samples, file_path, is_parquet, start_row, end_row = args
    try:
        if is_parquet:
            return _process_parquet_partition(
                target_samples, file_path, start_row, end_row
            )
        else:
            return _process_csv_partition(target_samples, file_path, start_row, end_row)
    except Exception as e:
        logger.warning(f"Error processing partition [{start_row}:{end_row}]: {e}")
        return {}


@dataclass
class ProjectStatistics:
    """Container for project statistics."""

    total_samples: int = 0
    samples_per_condition: Dict[str, int] = None
    total_proteins_quantified: int = 0
    proteins_per_sample: Dict[str, int] = None
    proteins_per_condition: Dict[str, int] = None
    proteins_across_samples: Dict[int, int] = None
    intensity_stats: Dict[str, Dict[str, float]] = None
    conditions: List[str] = None
    project_accession: Optional[str] = None
    feature_stats: Dict[str, any] = None
    psm_stats: Dict[str, any] = None

    def __post_init__(self):
        """Initialize mutable default values"""
        if self.samples_per_condition is None:
            self.samples_per_condition = {}
        if self.proteins_per_sample is None:
            self.proteins_per_sample = {}
        if self.proteins_per_condition is None:
            self.proteins_per_condition = {}
        if self.proteins_across_samples is None:
            self.proteins_across_samples = {}
        if self.intensity_stats is None:
            self.intensity_stats = {}
        if self.conditions is None:
            self.conditions = []
        if self.feature_stats is None:
            self.feature_stats = {}
        if self.psm_stats is None:
            self.psm_stats = {}


class ProjectReportGenerator:
    """
    Generate comprehensive statistical reports for QPX projects.

    This class analyzes ibaqpy TSV output files and generates detailed reports
    in multiple formats (text, JSON, HTML) with statistics and visualizations.
    """

    def __init__(
        self,
        ibaq_file: Union[str, Path],
        feature_file: Optional[Union[str, Path]] = None,
        psm_file: Optional[Union[str, Path]] = None,
        min_unique_peptides: int = 2,
        memory_limit_gb: Optional[float] = None,
        n_workers: int = 8,
    ):
        """
        Initialize the report generator.

        Args:
            ibaq_file: Path to ibaq TSV or parquet file (required)
            feature_file: Path to feature parquet file (optional)
            psm_file: Path to PSM parquet file (optional)
            min_unique_peptides: Minimum unique peptides per protein (default: 2, 0 to disable)
            memory_limit_gb: Memory usage limit in GB for batch processing (default: None, use available memory)
            n_workers: Number of parallel workers (default: 8)
        """
        self.ibaq_file = Path(ibaq_file)
        self.feature_file = Path(feature_file) if feature_file else None
        self.psm_file = Path(psm_file) if psm_file else None
        self.min_unique_peptides = min_unique_peptides

        if memory_limit_gb is None:
            available_memory_gb = psutil.virtual_memory().available / (1024**3)
            self.memory_limit_gb = available_memory_gb
        else:
            self.memory_limit_gb = memory_limit_gb

        self.n_workers = max(1, n_workers)

        self.stats = ProjectStatistics()
        self.ibaq_data = None
        self.feature_data = None
        self.psm_data = None

    def compute_statistics(self) -> ProjectStatistics:
        """Compute all statistics from ibaq TSV or parquet file."""

        if not self.ibaq_file.exists():
            raise FileNotFoundError(f"IBAQ file not found: {self.ibaq_file}")

        self._compute_ibaq_statistics()

        if self.feature_file and self.feature_file.exists():
            self._compute_feature_statistics()

        if self.psm_file and self.psm_file.exists():
            self._compute_psm_statistics()

        logger.info("Statistics computation completed")
        return self.stats

    def _compute_ibaq_statistics(self):
        """Compute statistics from IBAQ TSV or parquet file."""
        logger.info(f"Reading ibaq file: {self.ibaq_file}")

        try:
            if self.ibaq_file.suffix == ".parquet":
                parquet_file = pq.ParquetFile(str(self.ibaq_file))
                total_rows = parquet_file.metadata.num_rows

                n_workers = self.n_workers
                chunk_size = (total_rows + n_workers - 1) // n_workers

                tasks = []
                for i in range(n_workers):
                    start_row = i * chunk_size
                    num_rows = min(chunk_size, total_rows - start_row)
                    if num_rows > 0:
                        tasks.append((str(self.ibaq_file), start_row, num_rows))

                logger.info(f"Processing {total_rows:,} rows")

                with Pool(processes=n_workers) as pool:
                    results = pool.map(_process_chunk_for_proteins, tasks)

                all_proteins = set()
                for proteins in results:
                    all_proteins.update(proteins)

                initial_proteins = len(all_proteins)

                if self.min_unique_peptides > 0:
                    with Pool(processes=n_workers) as pool:
                        results = pool.map(_process_chunk_for_filter, tasks)

                    unique_peptides_count = {}
                    for counts in results:
                        for protein, sequences in counts.items():
                            if protein in unique_peptides_count:
                                unique_peptides_count[protein].update(sequences)
                            else:
                                unique_peptides_count[protein] = sequences

                    proteins_to_keep = set(
                        p
                        for p, seqs in unique_peptides_count.items()
                        if len(seqs) >= self.min_unique_peptides
                    )
                    final_proteins = len(proteins_to_keep)

                    logger.info(
                        f"Unique peptides filter (>={self.min_unique_peptides}): {initial_proteins:,} → {final_proteins:,} proteins"
                    )

                    df = pd.read_parquet(self.ibaq_file)
                    df["_protein"] = df["pg_accessions"].apply(
                        lambda x: x[0] if isinstance(x, list) and len(x) > 0 else str(x)
                    )
                    df = df[df["_protein"].isin(proteins_to_keep)].copy()
                    df = df.drop(columns=["_protein"])
                else:
                    df = pd.read_parquet(self.ibaq_file)

                column_mapping = {
                    "pg_accessions": "ProteinName",
                    "sample_accession": "SampleID",
                    "condition": "Condition",
                    "intensity": "Ibaq",
                }
                rename_dict = {
                    k: v for k, v in column_mapping.items() if k in df.columns
                }
                df = df.rename(columns=rename_dict)

                if "ProteinName" in df.columns and df["ProteinName"].dtype == "object":
                    df["ProteinName"] = df["ProteinName"].apply(
                        lambda x: "|".join(x) if isinstance(x, list) else str(x)
                    )
            else:
                df = pd.read_csv(self.ibaq_file, sep="\t")

            self.stats.total_proteins_quantified = df["ProteinName"].nunique()
            self.stats.total_samples = df["SampleID"].nunique()
            self.stats.samples_per_condition = (
                df.groupby("Condition")["SampleID"].nunique().to_dict()
            )
            self.stats.proteins_per_sample = (
                df.groupby("SampleID")["ProteinName"].nunique().to_dict()
            )
            self.stats.proteins_per_condition = (
                df.groupby("Condition")["ProteinName"].nunique().to_dict()
            )
            self.stats.conditions = sorted(df["Condition"].unique().tolist())

            protein_sample_counts = df.groupby("ProteinName")["SampleID"].nunique()
            self.stats.proteins_across_samples = (
                protein_sample_counts.value_counts().sort_index().to_dict()
            )

            for condition in df["Condition"].unique():
                condition_df = df[df["Condition"] == condition]
                if "Ibaq" in df.columns:
                    ibaq_data = condition_df[condition_df["Ibaq"] > 0]["Ibaq"]

                    self.stats.intensity_stats[condition] = {
                        "ibaq_mean": float(ibaq_data.mean()),
                        "ibaq_median": float(ibaq_data.median()),
                        "ibaq_std": float(ibaq_data.std()),
                        "ibaq_min": float(ibaq_data.min()),
                        "ibaq_max": float(ibaq_data.max()),
                    }

        except Exception as e:
            logger.error(f"Error computing IBAQ statistics: {e}")
            raise

    def _compute_feature_statistics(self):
        """Compute statistics from feature parquet file using chunked reading."""
        logger.info(f"Loading feature data from {self.feature_file}")

        parquet_file = pq.ParquetFile(str(self.feature_file))

        total_features = parquet_file.metadata.num_rows
        unique_sequences = set()
        unique_peptidoforms = set()
        unique_proteins = set()
        charge_counts = {}
        pep_values = []

        batch_size = min(500000, max(100000, int(self.memory_limit_gb * 50000)))

        for batch in parquet_file.iter_batches(batch_size=batch_size):
            chunk = batch.to_pandas()

            if "sequence" in chunk.columns:
                unique_sequences.update(chunk["sequence"].dropna().unique())

            if "peptidoform" in chunk.columns:
                unique_peptidoforms.update(chunk["peptidoform"].dropna().unique())

            if "pg_accessions" in chunk.columns:
                for pg in chunk["pg_accessions"].dropna():
                    if isinstance(pg, list):
                        unique_proteins.update(pg)
                    else:
                        unique_proteins.add(str(pg))

            if "precursor_charge" in chunk.columns:
                charge_vc = chunk["precursor_charge"].value_counts()
                for charge, count in charge_vc.items():
                    charge_counts[charge] = charge_counts.get(charge, 0) + count

            if "posterior_error_probability" in chunk.columns:
                pep_values.extend(
                    chunk["posterior_error_probability"].dropna().tolist()
                )

        self.stats.feature_stats = {
            "total_features": total_features,
            "unique_peptides": len(unique_sequences),
            "unique_peptidoforms": len(unique_peptidoforms),
            "unique_proteins": len(unique_proteins),
            "charge_distribution": charge_counts,
            "mean_pep": float(np.mean(pep_values)) if pep_values else 0.0,
            "median_pep": float(np.median(pep_values)) if pep_values else 0.0,
        }

        logger.info(
            f"Feature statistics computed: {self.stats.feature_stats['total_features']} features"
        )

    def _compute_psm_statistics(self):
        """Compute statistics from PSM parquet file using chunked reading."""
        logger.info(f"Loading PSM data from {self.psm_file}")

        parquet_file = pq.ParquetFile(str(self.psm_file))

        total_psms = parquet_file.metadata.num_rows
        unique_sequences = set()
        unique_peptidoforms = set()
        charge_counts = {}
        pep_values = []
        decoy_count = 0

        batch_size = min(500000, max(100000, int(self.memory_limit_gb * 50000)))

        for batch in parquet_file.iter_batches(batch_size=batch_size):
            chunk = batch.to_pandas()

            if "sequence" in chunk.columns:
                unique_sequences.update(chunk["sequence"].dropna().unique())

            if "peptidoform" in chunk.columns:
                unique_peptidoforms.update(chunk["peptidoform"].dropna().unique())

            if "precursor_charge" in chunk.columns:
                charge_vc = chunk["precursor_charge"].value_counts()
                for charge, count in charge_vc.items():
                    charge_counts[charge] = charge_counts.get(charge, 0) + count

            if "posterior_error_probability" in chunk.columns:
                pep_values.extend(
                    chunk["posterior_error_probability"].dropna().tolist()
                )

            if "is_decoy" in chunk.columns:
                decoy_count += int(chunk["is_decoy"].sum())

        self.stats.psm_stats = {
            "total_psms": total_psms,
            "unique_peptides": len(unique_sequences),
            "unique_peptidoforms": len(unique_peptidoforms),
            "charge_distribution": charge_counts,
            "mean_pep": float(np.mean(pep_values)) if pep_values else 0.0,
            "median_pep": float(np.median(pep_values)) if pep_values else 0.0,
            "decoy_count": decoy_count,
        }

        logger.info(
            f"PSM statistics computed: {self.stats.psm_stats['total_psms']} PSMs"
        )

    def generate_text_report(
        self, output_file: Optional[Union[str, Path]] = None
    ) -> str:
        """
        Generate a text-based statistical report.

        Args:
            output_file: Optional path to save the report

        Returns:
            Report text as string
        """
        if self.stats.total_samples == 0:
            self.compute_statistics()

        lines = []
        lines.append("=" * 80)
        lines.append("QPX PROJECT REPORT")
        lines.append("=" * 80)
        lines.append("")

        lines.append("SAMPLE STATISTICS")
        lines.append("-" * 80)
        lines.append(f"Total Samples: {self.stats.total_samples}")
        lines.append(f"Total Conditions: {len(self.stats.conditions)}")

        if self.stats.samples_per_condition:
            lines.append("\nSamples per Condition:")
            for condition, count in sorted(self.stats.samples_per_condition.items()):
                lines.append(f"  {condition}: {count}")

        lines.append("")

        lines.append("PROTEIN STATISTICS")
        lines.append("-" * 80)
        lines.append(
            f"Total Proteins Quantified: {self.stats.total_proteins_quantified}"
        )

        if self.stats.proteins_per_sample:
            lines.append("\nProteins per Sample Statistics:")
            sample_counts = list(self.stats.proteins_per_sample.values())
            lines.append(f"  Mean: {np.mean(sample_counts):.1f}")
            lines.append(f"  Median: {np.median(sample_counts):.1f}")
            lines.append(f"  Min: {np.min(sample_counts)}")
            lines.append(f"  Max: {np.max(sample_counts)}")
            lines.append(f"  Std: {np.std(sample_counts):.1f}")

        lines.append("")

        if self.stats.proteins_per_condition:
            lines.append("PROTEINS PER CONDITION")
            lines.append("-" * 80)
            for condition, count in sorted(self.stats.proteins_per_condition.items()):
                lines.append(f"  {condition}: {count}")

        lines.append("")

        if self.stats.intensity_stats:
            lines.append("INTENSITY DISTRIBUTION STATISTICS")
            lines.append("-" * 80)

            for condition, stats_dict in sorted(self.stats.intensity_stats.items()):
                lines.append(f"\nCondition: {condition}")

                if "mean" in stats_dict:
                    lines.append("  Raw Intensity Statistics:")
                    lines.append(f"    Mean: {stats_dict['mean']:.2e}")
                    lines.append(f"    Median: {stats_dict['median']:.2e}")
                    lines.append(f"    Std: {stats_dict['std']:.2e}")
                    lines.append(f"    Min: {stats_dict['min']:.2e}")
                    lines.append(f"    Max: {stats_dict['max']:.2e}")
                    lines.append(f"    Q25: {stats_dict['q25']:.2e}")
                    lines.append(f"    Q75: {stats_dict['q75']:.2e}")

                if "ibaq_mean" in stats_dict:
                    lines.append("  IBAQ Statistics:")
                    lines.append(f"    Mean: {stats_dict['ibaq_mean']:.2e}")
                    lines.append(f"    Median: {stats_dict['ibaq_median']:.2e}")
                    lines.append(f"    Std: {stats_dict['ibaq_std']:.2e}")
                    lines.append(f"    Min: {stats_dict['ibaq_min']:.2e}")
                    lines.append(f"    Max: {stats_dict['ibaq_max']:.2e}")

        lines.append("")
        lines.append("=" * 80)
        lines.append("END OF REPORT")
        lines.append("=" * 80)

        report_text = "\n".join(lines)

        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(report_text)
            logger.info(f"Report saved to {output_path}")

        return report_text

    def generate_json_report(self, output_file: Union[str, Path]) -> Dict:
        """
        Generate a JSON-formatted report with all statistics.

        Args:
            output_file: Path to save the JSON report

        Returns:
            Dictionary with all statistics
        """
        if self.stats.total_samples == 0:
            self.compute_statistics()

        report_dict = {
            "project_accession": self.stats.project_accession,
            "sample_statistics": {
                "total_samples": self.stats.total_samples,
                "total_conditions": len(self.stats.conditions),
                "conditions": self.stats.conditions,
                "samples_per_condition": self.stats.samples_per_condition,
            },
            "protein_statistics": {
                "total_proteins_quantified": self.stats.total_proteins_quantified,
                "proteins_per_condition": self.stats.proteins_per_condition,
                "proteins_per_sample_summary": {
                    "mean": (
                        float(np.mean(list(self.stats.proteins_per_sample.values())))
                        if self.stats.proteins_per_sample
                        else 0
                    ),
                    "median": (
                        float(np.median(list(self.stats.proteins_per_sample.values())))
                        if self.stats.proteins_per_sample
                        else 0
                    ),
                    "min": (
                        int(np.min(list(self.stats.proteins_per_sample.values())))
                        if self.stats.proteins_per_sample
                        else 0
                    ),
                    "max": (
                        int(np.max(list(self.stats.proteins_per_sample.values())))
                        if self.stats.proteins_per_sample
                        else 0
                    ),
                    "std": (
                        float(np.std(list(self.stats.proteins_per_sample.values())))
                        if self.stats.proteins_per_sample
                        else 0
                    ),
                },
            },
            "intensity_statistics": self.stats.intensity_stats,
        }

        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(report_dict, f, indent=2)

        logger.info(f"JSON report saved to {output_path}")
        return report_dict

    def _filter_chunk_for_dimensionality_reduction(
        self,
        chunk: pd.DataFrame,
        remove_contaminants: bool,
        remove_conditions: Optional[List[str]],
    ) -> pd.DataFrame:
        """Filter a data chunk for dimensionality reduction."""
        chunk = chunk[chunk["intensity"] > 0]

        if remove_contaminants and len(chunk) > 0:
            mask = chunk["pg_accessions"].apply(
                lambda x: (
                    any("CONT" in str(p) for p in x)
                    if isinstance(x, list)
                    else "CONT" in str(x)
                )
            )
            chunk = chunk[~mask]

        if "condition" in chunk.columns and len(chunk) > 0:
            chunk["condition"] = chunk["condition"].astype(str).str.split("|").str[-1]

            if remove_conditions:
                chunk = chunk[~chunk["condition"].isin(remove_conditions)]

        return chunk

    def _read_parquet_batches(
        self, remove_contaminants: bool, remove_conditions: Optional[List[str]]
    ) -> pd.DataFrame:
        """Read and filter parquet data in batches."""
        parquet_file = pq.ParquetFile(str(self.ibaq_file))
        df_chunks = []
        batch_size = 500000

        for batch in parquet_file.iter_batches(
            batch_size=batch_size,
            columns=["sample_accession", "condition", "pg_accessions", "intensity"],
        ):
            chunk = batch.to_pandas()
            chunk = self._filter_chunk_for_dimensionality_reduction(
                chunk, remove_contaminants, remove_conditions
            )
            if len(chunk) > 0:
                df_chunks.append(chunk)

        df = pd.concat(df_chunks, ignore_index=True)
        del df_chunks
        return df

    def _get_valid_samples(
        self, df: pd.DataFrame, min_samples_per_condition: int
    ) -> Tuple[Optional[List[str]], Optional[pd.DataFrame]]:
        """Get valid samples based on minimum samples per condition."""
        if "sample_accession" not in df.columns or "condition" not in df.columns:
            logger.error(
                "Required columns 'sample_accession' and 'condition' not found"
            )
            return None, None

        sample_conditions = df[["sample_accession", "condition"]].drop_duplicates()

        condition_counts = sample_conditions["condition"].value_counts()
        valid_conditions = condition_counts[
            condition_counts >= min_samples_per_condition
        ].index

        if len(valid_conditions) == 0:
            logger.warning(f"No conditions with >= {min_samples_per_condition} samples")
            return None, None

        valid_samples = sample_conditions[
            sample_conditions["condition"].isin(valid_conditions)
        ]["sample_accession"].tolist()

        return valid_samples, sample_conditions

    def _aggregate_partition_results(self, results: List[Dict]) -> Dict:
        """Aggregate results from parallel partitions."""
        sample_data = {}
        for result in results:
            for sample, protein_dict in result.items():
                if sample not in sample_data:
                    sample_data[sample] = {}
                for protein, intensity in protein_dict.items():
                    if protein in sample_data[sample]:
                        sample_data[sample][protein] = max(
                            sample_data[sample][protein], intensity
                        )
                    else:
                        sample_data[sample][protein] = intensity
        return sample_data

    def _build_intensity_matrix(
        self, sample_data: Dict, sample_conditions: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Build intensity matrix from sample data."""
        all_proteins = set()
        for sample_dict in sample_data.values():
            all_proteins.update(sample_dict.keys())

        matrix_data = []
        matrix_index = []
        for sample, protein_dict in sample_data.items():
            matrix_index.append(sample)
            row = [protein_dict.get(protein, 0) for protein in all_proteins]
            matrix_data.append(row)

        matrix = pd.DataFrame(
            matrix_data, index=matrix_index, columns=list(all_proteins)
        )

        labels = sample_conditions[
            sample_conditions["sample_accession"].isin(matrix.index)
        ].set_index("sample_accession")
        labels = labels.loc[matrix.index]

        return matrix, labels

    def _prepare_matrix_for_dimensionality_reduction(
        self,
        remove_contaminants: bool = True,
        remove_conditions: Optional[List[str]] = None,
        min_samples_per_condition: int = 2,
    ) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
        """
        Prepare intensity matrix for dimensionality reduction from parquet data.

        Args:
            remove_contaminants: Remove proteins with "CONT" in name
            remove_conditions: List of conditions to remove (e.g., ["Reference", "not available"])
            min_samples_per_condition: Minimum samples required per condition

        Returns:
            Tuple of (matrix, labels) or (None, None) if insufficient data or not parquet
        """
        if self.ibaq_file.suffix != ".parquet":
            logger.info("Dimensionality reduction only available for parquet files")
            return None, None

        # Read and filter parquet data
        df = self._read_parquet_batches(remove_contaminants, remove_conditions)

        # Get valid samples
        valid_samples, sample_conditions = self._get_valid_samples(
            df, min_samples_per_condition
        )
        if valid_samples is None:
            return None, None

        df = df[df["sample_accession"].isin(valid_samples)].copy()

        if "intensity" not in df.columns:
            logger.error("'intensity' column not found")
            return None, None

        # Normalize protein column
        df = _normalize_protein_column(df)

        # Process partitions in parallel
        samples = df["sample_accession"].unique()
        n_workers = self.n_workers
        target_samples = samples.tolist()

        parquet_file = pq.ParquetFile(str(self.ibaq_file))
        total_rows = parquet_file.metadata.num_rows

        partition_size = total_rows // n_workers
        partitions = []
        for i in range(n_workers):
            start_row = i * partition_size
            end_row = (i + 1) * partition_size if i < n_workers - 1 else None
            partitions.append(
                (target_samples, str(self.ibaq_file), True, start_row, end_row)
            )

        with Pool(processes=n_workers) as pool:
            results = pool.map(_process_file_partition, partitions)

        # Aggregate results and build matrix
        sample_data = self._aggregate_partition_results(results)
        matrix, labels = self._build_intensity_matrix(sample_data, sample_conditions)

        return matrix, labels

    def _compute_pca(
        self, matrix: pd.DataFrame, labels: pd.DataFrame
    ) -> Tuple[pd.DataFrame, np.ndarray]:
        """Compute PCA dimensionality reduction."""
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(matrix)

        pca = PCA(n_components=2, random_state=42)
        X_pca = pca.fit_transform(X_scaled)

        pca_df = pd.DataFrame(X_pca, columns=["PC1", "PC2"], index=matrix.index)
        pca_df["Sample"] = matrix.index
        pca_df["Condition"] = labels["condition"].values

        return pca_df, pca.explained_variance_ratio_

    def _compute_tsne(self, matrix: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
        """Compute t-SNE dimensionality reduction."""
        perplexity = min(30, matrix.shape[0] - 1)

        tsne = TSNE(
            n_components=2, random_state=42, metric="cosine", perplexity=perplexity
        )
        tsne_result = tsne.fit_transform(matrix)

        tsne_df = pd.DataFrame(
            tsne_result, columns=["TSNE1", "TSNE2"], index=matrix.index
        )
        tsne_df["Sample"] = matrix.index
        tsne_df["Condition"] = labels["condition"].values

        return tsne_df

    def _compute_umap(self, matrix: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
        """Compute UMAP dimensionality reduction."""
        # Use 'random' init for small samples to avoid spectral embedding issues
        init_method = "random" if matrix.shape[0] < 10 else "spectral"
        n_neighbors = min(15, matrix.shape[0] - 1)

        reducer = umap.UMAP(
            n_components=2,
            random_state=42,
            n_jobs=1,
            metric="cosine",
            n_neighbors=n_neighbors,
            min_dist=0.1,
            init=init_method,
        )
        umap_result = reducer.fit_transform(matrix)

        umap_df = pd.DataFrame(
            umap_result, columns=["UMAP1", "UMAP2"], index=matrix.index
        )
        umap_df["Sample"] = matrix.index
        umap_df["Condition"] = labels["condition"].values

        return umap_df

    def _create_ibaq_boxplot(self, plt, sns) -> Optional[str]:
        """Create FOT normalised iBAQ box plot."""
        import base64
        import io

        if self.ibaq_file.suffix != ".parquet":
            logger.info("Box plot only available for parquet files")
            return None

        if self.ibaq_data is None:
            self.ibaq_data = pd.read_parquet(self.ibaq_file)

        df = self.ibaq_data.copy()

        if "sample_accession" not in df.columns or "intensity" not in df.columns:
            logger.warning("Required columns not found for box plot")
            return None

        samples = sorted(
            df["sample_accession"].unique(),
            key=lambda x: [
                int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", x)
            ],
        )

        df_fot = df[df["intensity"] > 0].copy()

        # Calculate FOT normalization more efficiently using vectorized operations
        sample_totals = df_fot.groupby("sample_accession")["intensity"].sum()
        df_fot["sample_total"] = df_fot["sample_accession"].map(sample_totals)
        df_fot["fot_intensity"] = (df_fot["intensity"] / df_fot["sample_total"]) * 1e6
        df_fot = df_fot.drop(columns=["sample_total"])

        # Process samples sequentially to avoid memory issues
        data_to_plot = []

        for sample in samples:
            sample_data = df_fot[df_fot["sample_accession"] == sample]["fot_intensity"]
            if len(sample_data) > 0:
                filtered_data = sample_data[sample_data <= 1e9]
                data_to_plot.append(filtered_data)

        num_samples = len(samples)
        fig_height = max(8, num_samples * 0.3)

        fig, ax = plt.subplots(figsize=(10, fig_height))
        fig.patch.set_facecolor("white")

        ax.boxplot(
            data_to_plot,
            vert=False,
            patch_artist=True,
            widths=0.65,
            flierprops=dict(
                marker="o",
                markerfacecolor="none",
                markeredgecolor="gray",
                markersize=2,
                markeredgewidth=0.8,
                alpha=0.6,
                linestyle="none",
            ),
            boxprops=dict(facecolor="#4a90d9", edgecolor="#1f77b4", linewidth=1.5),
            medianprops=dict(color="orange", linewidth=2),
            whiskerprops=dict(color="#1f77b4", linewidth=1.5),
            capprops=dict(color="#1f77b4", linewidth=1.5),
        )

        ax.set_yticks(range(1, num_samples + 1))
        ax.set_yticklabels(samples, fontsize=10)
        ax.set_xlabel("FOT Intensity (ppm)", fontsize=12, fontweight="bold")
        ax.set_ylabel("Sample", fontsize=12, fontweight="bold")
        ax.set_title("FOT normalised iBAQ", fontsize=14, fontweight="bold", pad=20)
        ax.set_xscale("log")

        whisker_min = float("inf")
        whisker_max = 0
        for sample_data in data_to_plot:
            if len(sample_data) > 0:
                q1 = np.percentile(sample_data, 25)
                q3 = np.percentile(sample_data, 75)
                iqr = q3 - q1
                lower_whisker = max(sample_data.min(), q1 - 1.5 * iqr)
                upper_whisker = min(sample_data.max(), q3 + 1.5 * iqr)
                whisker_min = min(whisker_min, lower_whisker)
                whisker_max = max(whisker_max, upper_whisker)

        ax.set_xlim(whisker_min * 0.8, whisker_max * 1.2)
        ax.grid(True, alpha=0.4, axis="x", which="major", linestyle="-", linewidth=0.5)
        ax.set_facecolor("white")

        plt.tight_layout()

        img_buffer = io.BytesIO()
        plt.savefig(
            img_buffer,
            format="png",
            dpi=150,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )
        img_buffer.seek(0)
        img_base64 = base64.b64encode(img_buffer.read()).decode()
        plt.close()

        return f'<img src="data:image/png;base64,{img_base64}" style="max-width: 100%; height: auto;" alt="FOT normalised iBAQ Box Plot">'

    def _setup_plotting_libraries(self):
        """Setup matplotlib and seaborn for plotting."""
        try:
            import matplotlib

            matplotlib.use("Agg")  # Use non-interactive backend
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError:
            logger.error(
                "matplotlib and seaborn are required for HTML report generation"
            )
            raise

        sns.set_style("whitegrid")
        return plt, sns

    def _create_proteins_per_condition_plot(self) -> Optional[str]:
        """Create proteins per condition bar chart."""
        if not self.stats.proteins_per_condition:
            return None

        conditions = list(self.stats.proteins_per_condition.keys())
        counts = list(self.stats.proteins_per_condition.values())

        num_conditions = len(conditions)
        if num_conditions <= 10:
            tick_fontsize = 28
        elif num_conditions <= 20:
            tick_fontsize = 23
        elif num_conditions <= 30:
            tick_fontsize = 20
        else:
            tick_fontsize = 17

        fig_plotly = go.Figure()
        fig_plotly.add_trace(
            go.Bar(
                x=conditions,
                y=counts,
                marker=dict(
                    color=counts,
                    colorscale="Viridis",
                    line=dict(color="black", width=0.5),
                ),
                text=counts,
                textposition="none",
                hovertemplate="<b>%{x}</b><br>Proteins: %{y}<extra></extra>",
                showlegend=False,
            )
        )

        x_min = -0.5
        x_max = len(conditions) - 0.5
        y_min = 0
        y_max = max(counts) * 1.1
        initial_range = min(19.5, len(conditions) - 0.5)

        fig_plotly.update_layout(
            xaxis=dict(
                title="Condition",
                tickangle=-45,
                tickfont=dict(size=tick_fontsize),
                fixedrange=False,
                range=[-0.5, initial_range],
                rangeslider=dict(visible=False),
            ),
            yaxis=dict(
                title="Number of Proteins",
                gridcolor="lightgray",
                fixedrange=True,
                range=[y_min, y_max],
            ),
            hovermode="closest",
            dragmode="pan",
            height=600,
            plot_bgcolor="white",
            paper_bgcolor="white",
            margin=dict(l=80, r=40, t=40, b=150),
        )

        fig_plotly.update_xaxes(range=[-0.5, initial_range], autorange=False)
        fig_plotly.update_yaxes(range=[y_min, y_max], autorange=False)

        config = {
            "scrollZoom": False,
            "displayModeBar": True,
            "modeBarButtonsToRemove": [
                "select2d",
                "lasso2d",
                "autoScale2d",
            ],
            "displaylogo": False,
        }

        plotly_html = fig_plotly.to_html(
            include_plotlyjs="cdn", div_id="proteins_per_condition", config=config
        )

        sorted_cond_items = sorted(
            zip(conditions, counts), key=lambda x: x[1], reverse=True
        )
        sorted_conditions = [x[0] for x in sorted_cond_items]
        sorted_cond_counts = [x[1] for x in sorted_cond_items]

        range_limit_js = f"""
        <script>
        document.addEventListener('DOMContentLoaded', function() {{
            var plot = document.getElementById('proteins_per_condition');
            if (plot) {{
                var originalConditions = {conditions};
                var originalCounts = {counts};
                var sortedConditions = {sorted_conditions};
                var sortedCounts = {sorted_cond_counts};
                var isSorted = false;

                plot.on('plotly_relayout', function(eventdata) {{
                    if (eventdata['xaxis.range[0]'] !== undefined) {{
                        var x0 = eventdata['xaxis.range[0]'];
                        var x1 = eventdata['xaxis.range[1]'];
                        var xMin = {x_min};
                        var xMax = {x_max};

                        if (x0 < xMin || x1 > xMax) {{
                            var update = {{}};
                            update['xaxis.range[0]'] = Math.max(x0, xMin);
                            update['xaxis.range[1]'] = Math.min(x1, xMax);
                            Plotly.relayout(plot, update);
                        }}
                    }}
                }});
                var sortBtn = document.createElement('button');
                sortBtn.textContent = '↓ Sort by Value';
                sortBtn.style.cssText = (
                    'margin: 10px; padding: 8px 16px; background: #667eea; ' +
                    'color: white; border: none; border-radius: 4px; ' +
                    'cursor: pointer; font-size: 14px;'
                );
                sortBtn.onclick = function() {{
                    if (isSorted) {{
                        Plotly.restyle(plot, {{
                            x: [originalConditions],
                            y: [originalCounts],
                            'marker.color': [originalCounts]
                        }});
                        sortBtn.textContent = '↓ Sort by Value';
                        isSorted = false;
                    }} else {{
                        Plotly.restyle(plot, {{
                            x: [sortedConditions],
                            y: [sortedCounts],
                            'marker.color': [sortedCounts]
                        }});
                        sortBtn.textContent = '↑ Sort by Name';
                        isSorted = true;
                    }}
                }};
                plot.parentElement.insertBefore(sortBtn, plot);
            }}
        }});
        </script>
        """
        plotly_html = plotly_html + range_limit_js

        return f"""
    <div class="plot-container">
        <h3>Proteins per Condition ({len(conditions)} Conditions)</h3>
        {plotly_html}
    </div>
    """

    def _create_intensity_distribution_plot(self, plt, sns) -> Optional[str]:
        """Create intensity distribution violin plot."""
        if not self.stats.intensity_stats:
            return None

        fig, ax = plt.subplots(figsize=(12, 6))

        conditions = []
        log_intensities = []

        for condition, stats_dict in self.stats.intensity_stats.items():
            if "mean" in stats_dict and "std" in stats_dict:
                mean = stats_dict["mean"]
                std = stats_dict["std"]
                log_mean = np.log10(mean + 1)
                log_std = np.log10(std + 1) / 3  # Scale down std
                samples = np.random.normal(log_mean, log_std, 100)

                conditions.extend([condition] * len(samples))
                log_intensities.extend(samples)

        if not log_intensities:
            plt.close()
            return None

        plot_df = pd.DataFrame(
            {"Condition": conditions, "Log10(Intensity)": log_intensities}
        )

        sns.violinplot(
            data=plot_df,
            x="Condition",
            y="Log10(Intensity)",
            ax=ax,
            palette="muted",
        )
        ax.set_xlabel("Condition", fontsize=14)
        ax.set_ylabel("Log10(Intensity)", fontsize=14)
        plt.xticks(rotation=45, ha="right", fontsize=16)
        plt.tight_layout()

        html = self._fig_to_html(fig, "Intensity Distribution")
        plt.close()
        return html

    def _create_peptides_distribution_plot(self, plt) -> Optional[str]:
        """Create peptides distribution histogram."""
        if (
            not hasattr(self.stats, "peptides_per_sample")
            or not self.stats.peptides_per_sample
        ):
            return None

        fig, ax = plt.subplots(figsize=(10, 6))
        counts = list(self.stats.peptides_per_sample.values())

        ax.hist(counts, bins=30, color="coral", alpha=0.7, edgecolor="black")
        ax.axvline(
            np.mean(counts),
            color="red",
            linestyle="--",
            linewidth=2,
            label=f"Mean: {np.mean(counts):.0f}",
        )
        ax.axvline(
            np.median(counts),
            color="blue",
            linestyle="--",
            linewidth=2,
            label=f"Median: {np.median(counts):.0f}",
        )
        ax.set_xlabel("Number of Peptides", fontsize=14)
        ax.set_ylabel("Frequency", fontsize=14)
        ax.legend()
        ax.grid(axis="y", alpha=0.3)

        plt.tight_layout()
        html = self._fig_to_html(fig, "Peptides Distribution")
        plt.close()
        return html

    def _create_sample_distribution_pie_chart(self) -> Optional[str]:
        """Create sample distribution pie chart."""
        if not self.stats.samples_per_condition:
            return None

        conditions = list(self.stats.samples_per_condition.keys())
        counts = list(self.stats.samples_per_condition.values())

        num_conditions = len(conditions)
        if num_conditions <= 10:
            pie_fontsize = 30
        elif num_conditions <= 20:
            pie_fontsize = 25
        elif num_conditions <= 30:
            pie_fontsize = 22
        else:
            pie_fontsize = 18

        fig_plotly = go.Figure()
        fig_plotly.add_trace(
            go.Pie(
                labels=conditions,
                values=counts,
                textfont=dict(size=pie_fontsize),
                textposition="auto",
                hovertemplate="<b>%{label}</b><br>Samples: %{value}<br>Percentage: %{percent}<extra></extra>",
            )
        )

        fig_plotly.update_layout(
            height=1000,
            margin=dict(l=60, r=60, t=60, b=60),
            legend=dict(
                font=dict(size=16), x=1.02, y=0.5, xanchor="left", yanchor="middle"
            ),
        )

        config = {
            "scrollZoom": False,
            "displayModeBar": True,
            "modeBarButtonsToRemove": [
                "zoom2d",
                "select2d",
                "lasso2d",
                "zoomIn2d",
                "zoomOut2d",
                "autoScale2d",
            ],
            "displaylogo": False,
        }

        plotly_html = fig_plotly.to_html(
            include_plotlyjs="cdn", div_id="sample_distribution", config=config
        )

        return f"""
    <div class="plot-container">
        <h3>Sample Distribution</h3>
        {plotly_html}
    </div>
    """

    def _add_parquet_specific_plots(self, plots_html: List[str], plt, sns) -> None:
        """Add plots specific to parquet files."""
        if self.ibaq_file.suffix != ".parquet":
            return

        project_id = (
            self.ibaq_file.stem.split("-")[0]
            if "-" in self.ibaq_file.stem
            else "Project"
        )

        boxplot_html = self._create_ibaq_boxplot(plt, sns)
        if boxplot_html:
            plots_html.append(
                f"""
                <div class="plot-container">
                    <h3>FOT normalised iBAQ Distribution</h3>
                    <p>FOT normalized intensity distribution across all samples (ppm).</p>
                    {boxplot_html}
                </div>
            """
            )

        self._add_dimensionality_reduction_plots(plots_html, project_id)

    def _create_proteins_across_samples_plot(self) -> Optional[str]:
        """Create proteins across samples bar chart."""
        if not self.stats.proteins_across_samples:
            return None

        sample_nums = list(self.stats.proteins_across_samples.keys())
        protein_counts = list(self.stats.proteins_across_samples.values())
        sample_labels = [str(n) for n in sample_nums]

        num_points = len(sample_nums)
        if num_points <= 10:
            tick_fontsize = 28
        elif num_points <= 20:
            tick_fontsize = 23
        elif num_points <= 30:
            tick_fontsize = 20
        else:
            tick_fontsize = 17

        fig_plotly = go.Figure()
        fig_plotly.add_trace(
            go.Bar(
                x=sample_labels,
                y=protein_counts,
                marker=dict(
                    color=protein_counts,
                    colorscale="Viridis",
                    line=dict(color="black", width=0.5),
                ),
                text=protein_counts,
                textposition="none",
                hovertemplate="<b>%{x} samples</b><br>Proteins: %{y}<extra></extra>",
                showlegend=False,
            )
        )

        x_min = -0.5
        x_max = len(sample_labels) - 0.5
        y_min = 0
        y_max = max(protein_counts) * 1.1
        initial_range = min(19.5, len(sample_labels) - 0.5)

        fig_plotly.update_layout(
            xaxis=dict(
                title="Number of samples",
                tickangle=-45,
                tickfont=dict(size=tick_fontsize),
                fixedrange=False,
                range=[-0.5, initial_range],
                rangeslider=dict(visible=False),
            ),
            yaxis=dict(
                title="Number of Proteins",
                gridcolor="lightgray",
                fixedrange=True,
                range=[y_min, y_max],
            ),
            hovermode="closest",
            dragmode="pan",
            height=600,
            plot_bgcolor="white",
            paper_bgcolor="white",
            margin=dict(l=80, r=40, t=40, b=150),
        )

        fig_plotly.update_xaxes(range=[-0.5, initial_range], autorange=False)
        fig_plotly.update_yaxes(range=[y_min, y_max], autorange=False)

        config = {
            "scrollZoom": False,
            "displayModeBar": True,
            "modeBarButtonsToRemove": [
                "select2d",
                "lasso2d",
                "autoScale2d",
            ],
            "displaylogo": False,
        }

        plotly_html = fig_plotly.to_html(
            include_plotlyjs="cdn", div_id="proteins_across_samples", config=config
        )

        sorted_items = sorted(
            zip(sample_labels, protein_counts), key=lambda x: x[1], reverse=True
        )
        sorted_labels = [x[0] for x in sorted_items]
        sorted_counts = [x[1] for x in sorted_items]

        range_limit_js = f"""
        <script>
        document.addEventListener('DOMContentLoaded', function() {{
            var plot = document.getElementById('proteins_across_samples');
            if (plot) {{
                var originalLabels = {sample_labels};
                var originalCounts = {protein_counts};
                var sortedLabels = {sorted_labels};
                var sortedCounts = {sorted_counts};
                var isSorted = false;

                plot.on('plotly_relayout', function(eventdata) {{
                    if (eventdata['xaxis.range[0]'] !== undefined) {{
                        var x0 = eventdata['xaxis.range[0]'];
                        var x1 = eventdata['xaxis.range[1]'];
                        var xMin = {x_min};
                        var xMax = {x_max};

                        if (x0 < xMin || x1 > xMax) {{
                            var update = {{}};
                            update['xaxis.range[0]'] = Math.max(x0, xMin);
                            update['xaxis.range[1]'] = Math.min(x1, xMax);
                            Plotly.relayout(plot, update);
                        }}
                    }}
                }});

                var sortBtn = document.createElement('button');
                sortBtn.textContent = '↓ Sort by Value';
                sortBtn.style.cssText = (
                    'margin: 10px; padding: 8px 16px; background: #667eea; ' +
                    'color: white; border: none; border-radius: 4px; ' +
                    'cursor: pointer; font-size: 14px;'
                );
                sortBtn.onclick = function() {{
                    if (isSorted) {{
                        Plotly.restyle(plot, {{
                            x: [originalLabels],
                            y: [originalCounts],
                            'marker.color': [originalCounts]
                        }}).then(function() {{
                            Plotly.relayout(plot, {{
                                'xaxis.range': [{x_min}, Math.min({initial_range}, {x_max})]
                            }});
                        }});
                        sortBtn.textContent = '↓ Sort by Value';
                        isSorted = false;
                    }} else {{
                        Plotly.restyle(plot, {{
                            x: [sortedLabels],
                            y: [sortedCounts],
                            'marker.color': [sortedCounts]
                        }}).then(function() {{
                            Plotly.relayout(plot, {{
                                'xaxis.range': [{x_min}, Math.min({initial_range}, {x_max})]
                            }});
                        }});
                        sortBtn.textContent = '↑ Sort by Sample Number';
                        isSorted = true;
                    }}
                }};
                plot.parentElement.insertBefore(sortBtn, plot);
            }}
        }});
        </script>
        """
        plotly_html = plotly_html + range_limit_js

        return f"""
    <div class="plot-container">
        <h3>Proteins Shared Across Samples ({len(sample_labels)} Groups)</h3>
        {plotly_html}
    </div>
    """

    def _create_proteins_per_sample_plot(self) -> Optional[str]:
        """Create proteins per sample bar chart."""
        if not self.stats.proteins_per_sample:
            return None

        sample_items = sorted(
            self.stats.proteins_per_sample.items(),
            key=lambda x: [
                int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", x[0])
            ],
        )
        sample_labels = [s[0] for s in sample_items]
        sample_counts = [s[1] for s in sample_items]
        mean_val = np.mean(sample_counts)

        fig_plotly = go.Figure()
        fig_plotly.add_trace(
            go.Bar(
                x=sample_labels,
                y=sample_counts,
                marker=dict(
                    color=sample_counts,
                    colorscale="Viridis",
                    line=dict(color="black", width=0.5),
                ),
                text=sample_counts,
                textposition="none",
                hovertemplate="<b>%{x}</b><br>Proteins: %{y}<extra></extra>",
                showlegend=False,
            )
        )

        fig_plotly.add_hline(
            y=mean_val,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Mean: {mean_val:.0f}",
            annotation_position="top right",
        )

        x_min = -0.5
        x_max = len(sample_labels) - 0.5
        y_min = 0
        y_max = max(sample_counts) * 1.1

        initial_display = min(30, len(sample_labels))

        fig_plotly.update_layout(
            xaxis=dict(
                title="Sample ID",
                tickangle=-45,
                tickfont=dict(size=16),
                fixedrange=False,
                range=[-0.5, initial_display - 0.5],
                rangeslider=dict(visible=False),
            ),
            yaxis=dict(
                title="Number of Proteins",
                gridcolor="lightgray",
                fixedrange=True,
                range=[y_min, y_max],
            ),
            hovermode="closest",
            dragmode="pan",
            height=600,
            plot_bgcolor="white",
            paper_bgcolor="white",
            margin=dict(l=80, r=40, t=80, b=200),
        )

        fig_plotly.update_xaxes(range=[-0.5, initial_display - 0.5], autorange=False)

        fig_plotly.update_yaxes(range=[y_min, y_max], autorange=False)

        config = {
            "scrollZoom": False,
            "displayModeBar": True,
            "modeBarButtonsToRemove": [
                "select2d",
                "lasso2d",
                "autoScale2d",
            ],
            "displaylogo": False,
        }

        plotly_html = fig_plotly.to_html(
            include_plotlyjs="cdn", div_id="proteins_per_sample", config=config
        )

        sorted_items = sorted(
            zip(sample_labels, sample_counts), key=lambda x: x[1], reverse=True
        )
        sorted_labels = [x[0] for x in sorted_items]
        sorted_counts = [x[1] for x in sorted_items]

        range_limit_js = f"""
        <script>
        document.addEventListener('DOMContentLoaded', function() {{
            var plot = document.getElementById('proteins_per_sample');
            if (plot) {{
                var originalLabels = {sample_labels};
                var originalCounts = {sample_counts};
                var sortedLabels = {sorted_labels};
                var sortedCounts = {sorted_counts};
                var isSorted = false;

                plot.on('plotly_relayout', function(eventdata) {{
                    if (eventdata['xaxis.range[0]'] !== undefined) {{
                        var x0 = eventdata['xaxis.range[0]'];
                        var x1 = eventdata['xaxis.range[1]'];
                        var xMin = {x_min};
                        var xMax = {x_max};

                        if (x0 < xMin || x1 > xMax) {{
                            var update = {{}};
                            update['xaxis.range[0]'] = Math.max(x0, xMin);
                            update['xaxis.range[1]'] = Math.min(x1, xMax);
                            Plotly.relayout(plot, update);
                        }}
                    }}
                }});

                var sortBtn = document.createElement('button');
                sortBtn.textContent = '↓ Sort by Value';
                sortBtn.style.cssText = (
                    'margin: 10px; padding: 8px 16px; background: #667eea; ' +
                    'color: white; border: none; border-radius: 4px; ' +
                    'cursor: pointer; font-size: 14px;'
                );
                sortBtn.onclick = function() {{
                    if (isSorted) {{
                        Plotly.restyle(plot, {{
                            x: [originalLabels],
                            y: [originalCounts],
                            'marker.color': [originalCounts]
                        }});
                        sortBtn.textContent = '↓ Sort by Value';
                        isSorted = false;
                    }} else {{
                        Plotly.restyle(plot, {{
                            x: [sortedLabels],
                            y: [sortedCounts],
                            'marker.color': [sortedCounts]
                        }});
                        sortBtn.textContent = '↑ Sort by Name';
                        isSorted = true;
                    }}
                }};
                plot.parentElement.insertBefore(sortBtn, plot);
            }}
        }});
        </script>
        """
        plotly_html = plotly_html + range_limit_js

        return f"""
    <div class="plot-container">
        <h3>Proteins per Sample ({len(sample_counts)} Samples)</h3>
        {plotly_html}
    </div>
    """

    def _create_dimensionality_reduction_plot(
        self, result_df: pd.DataFrame, algorithm: str, project_id: str = "Project"
    ) -> str:
        """Create interactive dimensionality reduction plot using Plotly."""
        import plotly.express as px

        if algorithm == "PCA":
            axis_text = "PC"
            title_text = "PCA"
        elif algorithm == "TSNE":
            axis_text = "TSNE"
            title_text = "t-SNE"
        elif algorithm == "UMAP":
            axis_text = "UMAP"
            title_text = "UMAP"
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}")

        sorted_conditions = sorted(result_df["Condition"].unique())
        colors = (
            px.colors.qualitative.Plotly
            + px.colors.qualitative.D3
            + px.colors.qualitative.Set1
        )
        color_map = {
            condition: colors[i % len(colors)]
            for i, condition in enumerate(sorted_conditions)
        }

        fig = go.Figure()

        for condition in sorted_conditions:
            condition_data = result_df[result_df["Condition"] == condition]

            hover_text = []
            for _, row in condition_data.iterrows():
                text = f"<b>Sample:</b> {row['Sample']}<br>"
                text += f"<b>Condition:</b> {row['Condition']}<br>"
                text += f"<b>{axis_text}1:</b> {row[f'{axis_text}1']:.3f}<br>"
                text += f"<b>{axis_text}2:</b> {row[f'{axis_text}2']:.3f}"
                hover_text.append(text)

            fig.add_trace(
                go.Scatter(
                    x=condition_data[f"{axis_text}1"],
                    y=condition_data[f"{axis_text}2"],
                    mode="markers",
                    name=condition,
                    marker=dict(
                        size=10,
                        color=color_map[condition],
                        opacity=0.85,
                        line=dict(width=1, color="white"),
                        symbol="circle",
                    ),
                    hovertext=hover_text,
                    hoverinfo="text",
                )
            )

        fig.update_layout(
            title={
                "text": f"{title_text}({project_id})",
                "x": 0.5,
                "xanchor": "center",
                "font": {"size": 18},
            },
            xaxis_title=f"{axis_text}1",
            yaxis_title=f"{axis_text}2",
            hovermode="closest",
            plot_bgcolor="white",
            paper_bgcolor="white",
            font=dict(size=12),
            legend=dict(
                yanchor="top", y=0.99, xanchor="left", x=1.01, font=dict(size=16)
            ),
            margin=dict(l=60, r=60, t=60, b=60),
            height=1000,
        )

        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor="lightgray")
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor="lightgray")

        config = {
            "displayModeBar": True,
            "modeBarButtonsToRemove": [
                "zoom2d",
                "pan2d",
                "select2d",
                "lasso2d",
                "zoomIn2d",
                "zoomOut2d",
                "autoScale2d",
                "hoverClosestCartesian",
                "hoverCompareCartesian",
                "toggleSpikelines",
            ],
            "modeBarButtonsToAdd": [],
            "displaylogo": False,
            "toImageButtonOptions": {
                "format": "png",
                "filename": f"{project_id}_{algorithm}_plot",
                "height": 1000,
                "width": 1200,
                "scale": 2,
            },
        }

        return fig.to_html(
            include_plotlyjs="cdn", config=config, div_id=f"{algorithm.lower()}_plot"
        )

    def generate_html_report(self, output_file: Union[str, Path]) -> None:
        """Generate an HTML report with interactive visualizations."""
        if self.stats.total_samples == 0:
            self.compute_statistics()

        plt, sns = self._setup_plotting_libraries()
        plots_html = []

        # Add proteins per sample plot
        proteins_per_sample_html = self._create_proteins_per_sample_plot()
        if proteins_per_sample_html:
            plots_html.append(proteins_per_sample_html)

        # Add proteins per condition plot
        proteins_per_condition_html = self._create_proteins_per_condition_plot()
        if proteins_per_condition_html:
            plots_html.append(proteins_per_condition_html)

        # Add proteins across samples plot
        proteins_across_samples_html = self._create_proteins_across_samples_plot()
        if proteins_across_samples_html:
            plots_html.append(proteins_across_samples_html)

        # Add intensity distribution plot
        intensity_html = self._create_intensity_distribution_plot(plt, sns)
        if intensity_html:
            plots_html.append(intensity_html)

        # Add peptides distribution plot
        peptides_html = self._create_peptides_distribution_plot(plt)
        if peptides_html:
            plots_html.append(peptides_html)

        # Add sample distribution pie chart
        sample_dist_html = self._create_sample_distribution_pie_chart()
        if sample_dist_html:
            plots_html.append(sample_dist_html)

        # Add parquet-specific plots (ibaq boxplot and dimensionality reduction)
        self._add_parquet_specific_plots(plots_html, plt, sns)

        self._save_html_report(output_file, plots_html)

    def _add_dimensionality_reduction_plots(
        self, plots_html: List[str], project_id: str
    ) -> None:
        """Add dimensionality reduction plots (PCA, t-SNE, UMAP) to plots list."""
        matrix, labels = self._prepare_matrix_for_dimensionality_reduction(
            remove_contaminants=True,
            remove_conditions=["Reference", "not available"],
            min_samples_per_condition=2,
        )

        if matrix is None or labels is None:
            return

        logger.info("Computing PCA...")
        pca_df, _ = self._compute_pca(matrix, labels)
        pca_html = self._create_dimensionality_reduction_plot(pca_df, "PCA", project_id)
        plots_html.append(
            f"""
            <div class="plot-container">
                <h3>PCA</h3>
                {pca_html}
            </div>
        """
        )

        logger.info("Computing t-SNE...")
        tsne_df = self._compute_tsne(matrix, labels)
        tsne_html = self._create_dimensionality_reduction_plot(
            tsne_df, "TSNE", project_id
        )
        plots_html.append(
            f"""
            <div class="plot-container">
                <h3>t-SNE</h3>
                {tsne_html}
            </div>
        """
        )

        logger.info("Computing UMAP...")
        umap_df = self._compute_umap(matrix, labels)
        umap_html = self._create_dimensionality_reduction_plot(
            umap_df, "UMAP", project_id
        )
        plots_html.append(
            f"""
            <div class="plot-container">
                <h3>UMAP</h3>
                {umap_html}
            </div>
        """
        )

    def _save_html_report(
        self, output_file: Union[str, Path], plots_html: List[str]
    ) -> None:
        """Save HTML report to file."""
        html_content = self._generate_html_content(plots_html)

        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(html_content)

        logger.info(f"HTML report saved to {output_path}")

    def _fig_to_html(self, fig, title: str) -> str:
        """Convert matplotlib figure to HTML img tag with base64 encoding."""
        from io import BytesIO
        import base64

        buf = BytesIO()
        fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")
        buf.close()

        return f"""
        <div class="plot-container">
            <h3>{title}</h3>
            <img src="data:image/png;base64,{img_base64}" alt="{title}">
        </div>
        """

    def _generate_html_content(self, plots_html: List[str]) -> str:
        """Generate complete HTML document with embedded plots and statistics."""
        stats_html = self._generate_stats_tables()

        html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>QPX Project Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
            font-size: 16px;
        }}

        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}

        .header h1 {{
            margin: 0;
            font-size: 3em;
        }}

        .header p {{
            margin: 10px 0 0 0;
            font-size: 1.3em;
            opacity: 0.9;
        }}

        .section {{
            background: white;
            padding: 25px;
            margin-bottom: 25px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}

        .section h2 {{
            color: #667eea;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
            margin-top: 0;
            font-size: 2em;
        }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}

        .stat-card {{
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}

        .stat-card.blue {{
            background: linear-gradient(135deg, #0891b2 0%, #06b6d4 100%);
        }}

        .stat-card.green {{
            background: linear-gradient(135deg, #16a34a 0%, #4ade80 100%);
        }}

        .stat-card.purple {{
            background: linear-gradient(135deg, #2563eb 0%, #60a5fa 100%);
        }}

        .stat-card h3 {{
            margin: 0 0 10px 0;
            font-size: 1.1em;
            opacity: 0.9;
        }}

        .stat-card .value {{
            font-size: 3em;
            font-weight: bold;
            margin: 0;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
            font-size: 15px;
        }}

        th, td {{
            padding: 14px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}

        th {{
            background-color: #667eea;
            color: white;
            font-weight: bold;
            position: sticky;
            top: 0;
            z-index: 10;
            font-size: 16px;
        }}

        tr:hover {{
            background-color: #f5f5f5;
        }}

        .table-scroll-container {{
            max-height: 480px;
            overflow-y: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
        }}

        .table-scroll-container table {{
            margin: 0;
        }}

        .table-scroll-container-horizontal {{
            max-height: 480px;
            overflow-y: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
        }}

        .table-scroll-container-horizontal table {{
            margin: 0;
        }}

        .plot-container {{
            margin: 30px 0;
            text-align: center;
        }}

        .plot-container h3 {{
            color: #667eea;
            margin-bottom: 15px;
            font-size: 2.0em;
        }}

        .plot-container p {{
            font-size: 20px;
            line-height: 1.8;
        }}

        .plot-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}

        .footer {{
            text-align: center;
            padding: 20px;
            color: #666;
            font-size: 1em;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>QPX Project Report</h1>
        <p>Comprehensive Statistical Analysis of Proteomics Data</p>
        {f'<p>Project: {self.stats.project_accession}</p>' if self.stats.project_accession else ''}
    </div>

    <div class="section">
        <h2>Key Statistics</h2>
        <div class="stats-grid">
            <div class="stat-card purple">
                <h3>Samples</h3>
                <p class="value">{self.stats.total_samples}</p>
            </div>
            <div class="stat-card blue">
                <h3>Conditions</h3>
                <p class="value">{len(self.stats.conditions)}</p>
            </div>
            <div class="stat-card green">
                <h3>Proteins{
                    f' <span style="font-size: 0.8em; font-weight: normal;">'
                    f'(≥{self.min_unique_peptides} unique peptides)</span>'
                    if self.min_unique_peptides > 0 else ''
                }</h3>
                <p class="value">{self.stats.total_proteins_quantified}</p>
            </div>
        </div>
    </div>
    """

        feature_stats_html = ""
        if self.stats.feature_stats:
            feature_stats_html = f"""
    <div class="section">
        <h2>Feature Statistics</h2>
        <div class="stats-grid">
            <div class="stat-card blue">
                <h3>Total Features</h3>
                <p class="value">{self.stats.feature_stats['total_features']:,}</p>
            </div>
            <div class="stat-card green">
                <h3>Unique Peptides</h3>
                <p class="value">{self.stats.feature_stats['unique_peptides']:,}</p>
            </div>
            <div class="stat-card purple">
                <h3>Unique Peptidoforms</h3>
                <p class="value">{self.stats.feature_stats['unique_peptidoforms']:,}</p>
            </div>
            <div class="stat-card blue">
                <h3>Mean PEP</h3>
                <p class="value">{self.stats.feature_stats['mean_pep']:.4f}</p>
            </div>
        </div>
     </div>
     """

        psm_stats_html = ""
        if self.stats.psm_stats:
            psm_stats_html = f"""
    <div class="section">
        <h2>PSM Statistics</h2>
        <div class="stats-grid">
            <div class="stat-card blue">
                <h3>Total PSMs</h3>
                <p class="value">{self.stats.psm_stats['total_psms']:,}</p>
            </div>
            <div class="stat-card green">
                <h3>Unique Peptides</h3>
                <p class="value">{self.stats.psm_stats['unique_peptides']:,}</p>
            </div>
            <div class="stat-card purple">
                <h3>Unique Peptidoforms</h3>
                <p class="value">{self.stats.psm_stats['unique_peptidoforms']:,}</p>
            </div>
            <div class="stat-card blue">
                <h3>Mean PEP</h3>
                <p class="value">{self.stats.psm_stats['mean_pep']:.4f}</p>
            </div>
        </div>
    </div>
    """

        html += (
            feature_stats_html
            + psm_stats_html
            + f"""

    {stats_html}

    <div class="section">
        <h2>Visualizations</h2>
        {''.join(plots_html)}
    </div>

    <div class="footer">
        <p>Generated by QPX Report Generator</p>
        <p>For more information, visit <a href="https://github.com/bigbio/QPX">github.com/bigbio/QPX</a></p>
    </div>
</body>
</html>
        """
        )

        return html

    def _generate_stats_tables(self) -> str:
        """Generate HTML tables for detailed statistics."""
        if (
            self.stats.samples_per_condition
            or self.stats.proteins_per_condition
            or self.stats.intensity_stats
        ):
            all_conditions = set()
            if self.stats.samples_per_condition:
                all_conditions.update(self.stats.samples_per_condition.keys())
            if self.stats.proteins_per_condition:
                all_conditions.update(self.stats.proteins_per_condition.keys())
            if self.stats.intensity_stats:
                all_conditions.update(self.stats.intensity_stats.keys())

            table = '<div class="section"><h2>Statistics per Condition</h2>'
            table += '<div class="table-scroll-container-horizontal">'
            table += "<table>"
            table += (
                "<tr><th>Condition</th><th>Samples</th><th>Proteins</th>"
                "<th>Mean Intensity</th><th>Median Intensity</th><th>Std Dev</th>"
                "<th>Min Intensity</th><th>Max Intensity</th></tr>"
            )

            for condition in sorted(all_conditions):
                samples = (
                    self.stats.samples_per_condition.get(condition, 0)
                    if self.stats.samples_per_condition
                    else 0
                )
                proteins = (
                    self.stats.proteins_per_condition.get(condition, 0)
                    if self.stats.proteins_per_condition
                    else 0
                )

                intensity_data = (
                    self.stats.intensity_stats.get(condition, {})
                    if self.stats.intensity_stats
                    else {}
                )
                mean_val = intensity_data.get("ibaq_mean", 0)
                median_val = intensity_data.get("ibaq_median", 0)
                std_val = intensity_data.get("ibaq_std", 0)
                min_val = intensity_data.get("ibaq_min", 0)
                max_val = intensity_data.get("ibaq_max", 0)

                mean_str = f"{mean_val:.2e}" if mean_val > 0 else "N/A"
                median_str = f"{median_val:.2e}" if median_val > 0 else "N/A"
                std_str = f"{std_val:.2e}" if std_val > 0 else "N/A"
                min_str = f"{min_val:.2e}" if min_val > 0 else "N/A"
                max_str = f"{max_val:.2e}" if max_val > 0 else "N/A"

                table += f"""<tr>
                    <td>{condition}</td>
                    <td>{samples}</td>
                    <td>{proteins}</td>
                    <td>{mean_str}</td>
                    <td>{median_str}</td>
                    <td>{std_str}</td>
                    <td>{min_str}</td>
                    <td>{max_str}</td>
                </tr>"""

            table += "</table></div></div>"
            return table

        return ""
