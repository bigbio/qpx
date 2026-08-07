"""Spectronaut PG adapter -- report to pg.parquet.

Unlike DIA-NN, Spectronaut does not produce a separate PG matrix file.
Protein-group quantities (``PG.Quantity``) are present directly in the
report alongside precursor/fragment rows.  We aggregate per
(PG.ProteinGroups, R.FileName) to produce one PG record per protein
group per run.
"""

from __future__ import annotations

import logging
import re
from typing import Optional

import pandas as pd

from qpx.converters.base import resolve_columns
from qpx.converters.mappings import get_field_mappings
from qpx.converters.spectronaut.base_adapter import SpectronautBaseAdapter
from qpx.converters.utils import safe_float
from qpx.core.sql import sql_build, validate_identifier
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)

_EXT_RE = re.compile(r"\.(mzML|raw|d|wiff|htrms)$", re.IGNORECASE)


class SpectronautPgAdapter(SpectronautBaseAdapter):
    """Convert Spectronaut report to ``pg.parquet``.

    Usage::

        with SpectronautPgAdapter() as adapter:
            adapter.convert(
                spectronaut_report="report.tsv",
                output_path="pg.parquet",
            )
    """

    def __init__(self, **kwargs):
        """Initialize PG adapter with empty resolved mappings."""
        super().__init__(**kwargs)
        self._resolved_pg: dict | None = None

    def convert(
        self,
        spectronaut_report: str,
        output_path: str,
        file_num: int = 20,
        creator: str = "spectronaut",
    ) -> None:
        """Run the Spectronaut report -> pg.parquet conversion."""
        self._load_spectronaut_report(spectronaut_report)

        report_cols = self._get_report_columns()
        self._resolved_pg = resolve_columns(get_field_mappings("spectronaut", "pg"), report_cols)

        runs = self._get_unique_runs()

        with PgWriter(output_path, creator=creator, compression=self._compression) as writer:
            for i in range(0, len(runs), file_num):
                batch_runs = runs[i : i + file_num]
                self.logger.info("Processing PG runs %d-%d of %d", i + 1, min(i + file_num, len(runs)), len(runs))
                records = self._process_batch(batch_runs, report_cols)
                if records:
                    writer.write_batch(records)

        self.logger.info("Spectronaut PG conversion complete -> %s", output_path)

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _get_unique_runs(self) -> list[str]:
        """Get sorted list of unique run values from the report."""
        return self._query_distinct_runs(self._resolved_pg["run_file_name"])

    # ------------------------------------------------------------------
    # Batch processing
    # ------------------------------------------------------------------

    @staticmethod
    def _build_pg_select_parts(
        r: dict,
        report_cols: set[str],
    ) -> list[str]:
        """Build SELECT clause parts for PG aggregation query."""
        pg_col = r["pg_accessions"]
        run_col = r["run_file_name"]
        int_col = r["intensity"]

        parts = [
            f'FIRST(r."{pg_col}") AS pg_accessions',
            f"regexp_replace(FIRST(r.\"{run_col}\"), '\\.(mzML|raw|d|wiff|htrms)$', '', 'i') AS run_file_name",
            f'FIRST(CAST(r."{int_col}" AS DOUBLE)) AS pg_quantity',
        ]

        # Optional columns
        qv_col = r.get("qvalue")
        if qv_col and qv_col in report_cols:
            parts.append(f'FIRST(CAST(r."{qv_col}" AS DOUBLE)) AS qvalue')
        pqv_col = r.get("pg_qvalue_runwise")
        if pqv_col and pqv_col in report_cols:
            parts.append(f'FIRST(CAST(r."{pqv_col}" AS DOUBLE)) AS qvalue_runwise')
        pa_col = r.get("pg_protein_accessions")
        if pa_col and pa_col in report_cols:
            parts.append(f'FIRST(r."{pa_col}") AS protein_accessions')
        pn_col = r.get("pg_names")
        if pn_col and pn_col in report_cols:
            parts.append(f'FIRST(r."{pn_col}") AS pg_names')
        gg_col = r.get("gg_accessions")
        if gg_col and gg_col in report_cols:
            parts.append(f'FIRST(r."{gg_col}") AS gg_accessions')

        # Peptide counting columns
        seq_col = "PEP.StrippedSequence"
        prot_col = "PEP.IsProteotypic"
        if seq_col in report_cols:
            parts.append(f'COUNT(DISTINCT r."{seq_col}") AS total_sequences')
        if prot_col in report_cols and seq_col in report_cols:
            parts.append(f'COUNT(DISTINCT CASE WHEN CAST(r."{prot_col}" AS BOOLEAN) THEN r."{seq_col}" END) AS unique_sequences')
        return parts

    def _process_batch(
        self,
        runs: list[str],
        report_cols: set[str],
    ) -> list[dict]:
        """Process a batch of runs for PG quantification."""
        r = self._resolved_pg
        select_parts = self._build_pg_select_parts(r, report_cols)

        select_clause = ",\n                ".join(select_parts)
        placeholders = ", ".join(["?" for _ in runs])
        pg_id = validate_identifier(r["pg_accessions"])
        run_id = validate_identifier(r["run_file_name"])

        stmt = sql_build(
            "SELECT $select FROM report r WHERE r.$run IN ($ph) AND r.$pg IS NOT NULL GROUP BY r.$pg, r.$run",
            select=select_clause,
            run=run_id,
            ph=placeholders,
            pg=pg_id,
        )
        report_df = self._conn.execute(stmt, runs).df()

        if report_df.empty:
            return []

        records: list[dict] = []
        for _, row in report_df.iterrows():
            rec = self._build_pg_record(row)
            if rec:
                records.append(rec)

        return records

    @staticmethod
    def _extract_pg_names_genes(
        row: pd.Series,
    ) -> tuple[list[str] | None, list[str] | None]:
        """Extract protein names and gene accessions from a PG row."""
        pg_names = None
        if "pg_names" in row.index and pd.notna(row["pg_names"]):
            pg_names = str(row["pg_names"]).split(";")
        gg_accessions = None
        if "gg_accessions" in row.index and pd.notna(row["gg_accessions"]):
            gg_accessions = str(row["gg_accessions"]).split(";")
        return pg_names, gg_accessions

    def _build_pg_scores(
        self,
        row: pd.Series,
    ) -> tuple[float | None, float | None, list[dict] | None]:
        """Build additional_scores and extract q-values from a PG row."""
        global_qvalue = None
        if "qvalue" in row.index:
            global_qvalue = safe_float(row["qvalue"])
        pg_qvalue = None
        if "qvalue_runwise" in row.index:
            pg_qvalue = safe_float(row["qvalue_runwise"])

        additional_scores: list[dict] = []
        if global_qvalue is not None:
            additional_scores.append({"score_name": "pg_qvalue", "score_value": global_qvalue, "higher_better": False})
            self._discovered_scores.add("pg_qvalue")
        if pg_qvalue is not None:
            additional_scores.append({"score_name": "pg_qvalue_runwise", "score_value": pg_qvalue, "higher_better": False})
            self._discovered_scores.add("pg_qvalue_runwise")

        return global_qvalue, pg_qvalue, additional_scores or None

    def _build_pg_record(
        self,
        row: pd.Series,
    ) -> Optional[dict]:
        """Build a single PG record from an aggregated row."""
        pg_acc_str = str(row["pg_accessions"])
        pg_accessions = pg_acc_str.split(";")
        anchor_protein = pg_accessions[0] if pg_accessions else ""
        run_file_name = _EXT_RE.sub("", str(row["run_file_name"]))

        pg_names, gg_accessions = self._extract_pg_names_genes(row)

        # Intensities
        pg_quantity = safe_float(row["pg_quantity"]) or 0.0
        intensities = [{"label": "LFQ", "intensity": float(pg_quantity)}]

        # Q-values and scores
        global_qvalue, pg_qvalue, additional_scores = self._build_pg_scores(row)

        # Peptide counts
        total_sequences = int(row["total_sequences"]) if "total_sequences" in row.index else 0
        unique_sequences = int(row["unique_sequences"]) if "unique_sequences" in row.index else 0
        peptide_count = max(total_sequences, 1)
        peptides = [{"protein_name": acc, "peptide_count": peptide_count} for acc in pg_accessions]

        # is_decoy: Spectronaut normally pre-filters decoys; infer from prefix
        is_decoy = any(acc.startswith(("DECOY_", "decoy_", "rev_", "REV_")) for acc in pg_accessions)

        return {
            "pg_accessions": pg_accessions,
            "pg_names": pg_names,
            "gg_accessions": gg_accessions,
            "gg_names": gg_accessions,
            "gg_qvalue": None,
            "anchor_protein": anchor_protein,
            "grouped_runs": [run_file_name],
            "global_qvalue": global_qvalue,
            "pg_qvalue": pg_qvalue,
            "intensities": intensities,
            "additional_intensities": None,
            "is_decoy": is_decoy,
            "contaminant": None,
            "peptides": peptides,
            "peptide_counts": {
                "unique_sequences": unique_sequences,
                "total_sequences": total_sequences,
            },
            "feature_counts": None,
            "sequence_coverage": None,
            "molecular_weight": None,
            "additional_scores": additional_scores,
            "cv_params": None,
        }
