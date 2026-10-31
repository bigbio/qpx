"""Feature data structure — quantified peptide features per MS run."""

from __future__ import annotations

from typing import TYPE_CHECKING

from qpx.core.convert import QueryResult
from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.query import _escape_sql_string
from qpx.core.sample_resolution import feature_samples_sql
from qpx.core.sql import sql_build

FeatureSchema = load_schema("feature")

if TYPE_CHECKING:
    from qpx.dataset import Dataset


class Feature(BaseStructure):
    """Quantified peptide features per MS run."""

    _schema_class = FeatureSchema

    @classmethod
    def from_file(cls, path, **engine_kwargs) -> Feature:
        """Open a standalone feature.parquet file, guarding the pre-1.0 ``channel`` layout."""
        from qpx.version import check_feature_file_compatible

        check_feature_file_compatible(path)
        return super().from_file(path, **engine_kwargs)

    def unique_proteins(self) -> list[str]:
        """Get all unique anchor proteins."""
        stmt = sql_build("SELECT DISTINCT anchor_protein FROM $src", src=self._query.source)
        result = self._engine.execute(stmt)
        return [row[0] for row in result.fetchall()]

    def by_protein(self, protein: str) -> Feature:
        """Filter features by anchor protein: those of every group the protein leads."""
        return self.filter(f"anchor_protein = '{_escape_sql_string(protein)}'")

    def by_group(self, accessions: list[str]) -> Feature:
        """Filter to the features of the protein group whose membership is exactly ``accessions``, in any order."""
        members = ", ".join(f"'{_escape_sql_string(accession)}'" for accession in accessions)
        return self.filter(f"list_sort(list_transform(pg_accessions, lambda x: x.accession)) = list_sort([{members}])")

    def by_run(self, run_file_name: str) -> Feature:
        """Filter features by run file."""
        return self.filter(f"run_file_name = '{_escape_sql_string(run_file_name)}'")

    def peptide_intensities(self) -> QueryResult:
        """
        Flatten intensities: one row per (peptide, label, intensity).
        Useful for downstream analysis and plotting.
        """
        stmt = sql_build(
            """SELECT sequence, peptidoform, charge, anchor_protein,
               run_file_name, i.label, i.intensity
        FROM $src,
             UNNEST(intensities) AS _t(i)""",
            src=self._query.source,
        )
        return QueryResult(self._engine.connection, stmt)

    def intensities_wide(self, labels: list[str] | None = None) -> QueryResult:
        """One row per feature with one intensity column per label; *labels* picks the columns, default all."""
        if labels is None:
            stmt = sql_build(
                "SELECT DISTINCT i.label FROM $src, UNNEST(intensities) AS t(i) WHERE i.label IS NOT NULL ORDER BY 1",
                src=self._query.source,
            )
            labels = [row[0] for row in self._engine.execute(stmt).fetchall()]
        columns = "".join(
            f""", MAX(i.intensity) FILTER (WHERE i.label = '{_escape_sql_string(label)}') AS "{label.replace('"', '""')}\""""
            for label in labels
        )
        stmt = sql_build(
            "SELECT feature_id, peptidoform, charge, anchor_protein, run_file_name$columns "
            "FROM $src, UNNEST(intensities) AS t(i) "
            "GROUP BY feature_id, peptidoform, charge, anchor_protein, run_file_name",
            src=self._query.source,
            columns=columns,
        )
        return QueryResult(self._engine.connection, stmt)

    def for_quantification(self, dataset: Dataset) -> QueryResult:
        """
        Return feature intensities in long format with sample metadata,
        ready for downstream quantification tools.

        Returns one row per (peptidoform, charge, protein, sample, label)
        with columns: sequence, peptidoform, charge, anchor_protein,
        sample_accession, label, intensity, run_file_name,
        biological_replicate, fraction.
        """
        # fraction is an optional run column: an absent one reads as null.
        fraction = self._engine.optional_columns("run", "r", ("fraction",))
        stmt = sql_build(
            """SELECT q.sequence, q.peptidoform, q.charge, q.anchor_protein,
               q.sample_accession, q.label, q.intensity,
               q.run_file_name, q.biological_replicate, $fraction
        FROM ($feature_samples) q
        JOIN run r ON r.run_file_name = q.run_file_name""",
            feature_samples=feature_samples_sql(self._query.source, ("sequence", "peptidoform", "charge", "anchor_protein")),
            fraction=fraction,
        )
        return QueryResult(self._engine.connection, stmt)
