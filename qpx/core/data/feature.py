"""Feature data structure — quantified peptide features per MS run."""

from __future__ import annotations

from typing import TYPE_CHECKING

from qpx.core.convert import QueryResult
from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.query import _escape_sql_string
from qpx.core.sample_sql import feature_sample_ctes
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
        """Filter features by anchor protein."""
        return self.filter(f"anchor_protein = '{_escape_sql_string(protein)}'")

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
        ctes = feature_sample_ctes(
            ("sequence", "peptidoform", "charge", "anchor_protein"),
            source=self._query.source,
            run_extra=(", " + fraction, ", fraction"),
        )
        stmt = (
            ctes
            + """
        SELECT sequence, peptidoform, charge, anchor_protein,
               sample_accession, sample_label AS label, intensity,
               run_file_name, biological_replicate, fraction
        FROM feature_samples"""
        )
        return QueryResult(self._engine.connection, stmt)
