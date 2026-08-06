"""CDAP PG adapter -- aggregate features into per-run protein groups.

Reads the ``feature.parquet`` written by :class:`CdapFeatureAdapter` and rolls
features up to ``(anchor_protein, run_file_name)`` groups.  Reporter-ion or
LFQ intensities are summed across the contributing features per channel; the
peptide / feature counts and score summaries follow the QPX :class:`PgSchema`.
"""

from __future__ import annotations

from typing import Optional

import pandas as pd

from qpx.converters.cdap.base_adapter import CdapBaseAdapter
from qpx.writers.pg import PgWriter


class CdapPgAdapter(CdapBaseAdapter):
    """Convert CDAP-derived feature.parquet into pg.parquet."""

    def convert(
        self,
        feature_path: str,
        output_path: str,
        chunksize: int = 100_000,
        creator: str = "cdap",
    ) -> None:
        """Aggregate the given feature parquet into per-run protein groups."""
        self._load_parquet("features", feature_path)
        actual_cols = self.get_table_columns("features")
        required = {"anchor_protein", "run_file_name", "pg_accessions"}
        missing = required - set(actual_cols)
        if missing:
            raise ValueError(f"feature.parquet missing columns: {sorted(missing)}")

        sql = self._build_aggregation_sql()
        self.logger.info("Transforming CDAP protein groups ...")
        with PgWriter(output_path, creator=creator, compression=self._compression) as writer:
            self._stream_transform_write(sql, writer, chunksize)

        self.logger.info("CDAP PG conversion complete -> %s", output_path)

    # ------------------------------------------------------------------
    # SQL aggregation
    # ------------------------------------------------------------------

    @staticmethod
    def _build_aggregation_sql() -> str:
        """Build a query that aggregates features per (protein group, run).

        The group key is the FULL protein-group membership (the sorted set of
        ``pg_accessions`` accessions), not the leading ``anchor_protein``.
        Keying on the leader alone collapsed two genuinely distinct groups that
        share a leading protein (e.g. ``P1;P2`` and ``P1;P3``) into one, and
        filed shared peptides under whichever accession happened to be listed
        first. This mirrors the pg-identity fix (bigbio/qpx#240) that keys
        ``pg_id`` on ``pg_accessions``; ``anchor_protein`` stays a descriptive
        leader carried through with ``any_value``.
        """
        return (
            "WITH keyed AS (\n"
            "    SELECT\n"
            "        *,\n"
            "        COALESCE(\n"
            "            array_to_string(list_sort(list_transform(pg_accessions, x -> x.accession)), ';'),\n"
            "            anchor_protein\n"
            "        ) AS pg_membership_key\n"
            "    FROM features\n"
            ")\n"
            "SELECT\n"
            "    any_value(anchor_protein) AS anchor_protein,\n"
            "    run_file_name,\n"
            "    any_value(is_decoy) AS any_decoy,\n"
            "    bool_and(is_decoy) AS all_decoy,\n"
            "    count(*) AS feature_count,\n"
            "    count(DISTINCT sequence) AS unique_sequence_count,\n"
            "    count(DISTINCT peptidoform) AS unique_peptidoform_count,\n"
            "    list(sequence) AS sequences,\n"
            "    list(intensities) AS intensities_per_feature,\n"
            "    list(pg_accessions) AS pg_accessions_per_feature,\n"
            "    list(additional_scores) AS additional_scores_per_feature\n"
            "FROM keyed\n"
            "GROUP BY pg_membership_key, run_file_name"
        )

    # ------------------------------------------------------------------
    # Transformation
    # ------------------------------------------------------------------

    def _transform_batch(self, df: pd.DataFrame) -> list[dict]:
        """Build PG records from one aggregation batch."""
        return self._transform_batch_with_skip(df, self._transform_row, "PG")

    def _transform_row(self, row: dict) -> Optional[dict]:
        """Convert one aggregated feature group into a PG record."""
        anchor_protein = row.get("anchor_protein")
        if not isinstance(anchor_protein, str) or not anchor_protein:
            return None
        run_file_name = row.get("run_file_name") or ""

        # Reporter / LFQ intensities: union channels across the features then sum.
        intensities = self._merge_intensities(row.get("intensities_per_feature"))

        # Protein-group accessions: union of every feature's anchored group, with
        # anchor_protein guaranteed to be the first entry.
        pg_accessions = self._collect_pg_accessions(anchor_protein, row.get("pg_accessions_per_feature"))

        # Peptide bookkeeping
        feature_count = int(row.get("feature_count") or 0)
        unique_sequence = int(row.get("unique_sequence_count") or 0)
        unique_peptidoform = int(row.get("unique_peptidoform_count") or 0)
        sequences = row.get("sequences")
        if sequences is None:
            sequences = []
        peptides = [
            {"protein_name": anchor_protein, "peptide_count": unique_sequence},
        ]

        # Best-feature scores rolled up by ``score_name`` (use min for FDR-ish
        # metrics, max for raw scores -- decided by ``higher_better``).
        additional_scores = self._roll_up_scores(row.get("additional_scores_per_feature"))

        is_decoy = bool(row.get("all_decoy"))

        return {
            "pg_accessions": pg_accessions,
            "pg_names": pg_accessions,
            "gg_accessions": None,
            "gg_names": None,
            "gg_qvalue": None,
            "anchor_protein": anchor_protein,
            "grouped_runs": [run_file_name],
            "global_qvalue": None,
            "pg_qvalue": None,
            "intensities": intensities or None,
            "additional_intensities": None,
            "is_decoy": is_decoy,
            "contaminant": False,
            "peptides": peptides,
            "peptide_counts": {
                "unique_sequences": unique_sequence,
                "total_sequences": len(sequences) if sequences is not None else 0,
            },
            "feature_counts": {
                "unique_features": unique_peptidoform,
                "total_features": feature_count,
            },
            "sequence_coverage": None,
            "molecular_weight": None,
            "additional_scores": additional_scores or None,
            "cv_params": None,
        }

    # ------------------------------------------------------------------
    # Helpers: merge nested lists from DuckDB
    # ------------------------------------------------------------------

    @staticmethod
    def _merge_intensities(intensities_per_feature) -> list[dict]:
        """Sum intensities by channel across all features in the group."""
        if intensities_per_feature is None:
            return []
        totals: dict[str, float] = {}
        for feature_intensities in intensities_per_feature:
            if feature_intensities is None or (hasattr(feature_intensities, "__len__") and len(feature_intensities) == 0):
                continue
            for entry in feature_intensities:
                if entry is None:
                    continue
                label = entry.get("label") if isinstance(entry, dict) else None
                value = entry.get("intensity") if isinstance(entry, dict) else None
                if label is None or value is None:
                    continue
                try:
                    totals[label] = totals.get(label, 0.0) + float(value)
                except (TypeError, ValueError):
                    continue
        return [{"label": label, "intensity": value} for label, value in totals.items()]

    @staticmethod
    def _collect_pg_accessions(anchor_protein: str, pg_accessions_per_feature) -> list[str]:
        """Union pg_accessions across features, keeping anchor first."""
        seen: list[str] = [anchor_protein]
        seen_set = {anchor_protein}
        if pg_accessions_per_feature is None:
            return seen
        for entries in pg_accessions_per_feature:
            if entries is None or (hasattr(entries, "__len__") and len(entries) == 0):
                continue
            for entry in entries:
                accession = entry.get("accession") if isinstance(entry, dict) else None
                if not accession or accession in seen_set:
                    continue
                seen_set.add(accession)
                seen.append(accession)
        return seen

    @staticmethod
    def _roll_up_scores(scores_per_feature) -> list[dict]:
        """Reduce per-feature score lists to one entry per ``score_name``.

        Uses the best value -- max for ``higher_better=True``, min otherwise.
        """
        if scores_per_feature is None:
            return []
        bucket: dict[str, dict] = {}
        for scores in scores_per_feature:
            if scores is None or (hasattr(scores, "__len__") and len(scores) == 0):
                continue
            for entry in scores:
                if not isinstance(entry, dict):
                    continue
                name = entry.get("score_name")
                value = entry.get("score_value")
                higher_better = entry.get("higher_better")
                if name is None or value is None:
                    continue
                try:
                    numeric = float(value)
                except (TypeError, ValueError):
                    continue
                existing = bucket.get(name)
                if existing is None:
                    bucket[name] = {
                        "score_name": name,
                        "score_value": numeric,
                        "higher_better": higher_better,
                    }
                    continue
                keep_new = numeric > existing["score_value"] if higher_better else numeric < existing["score_value"]
                if keep_new:
                    existing["score_value"] = numeric
        return list(bucket.values())
