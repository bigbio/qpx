"""PG data structure — Protein Groups with quantification."""

import logging

from qpx.core.convert import QueryResult
from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.data.run import Run
from qpx.core.query import _escape_sql_string
from qpx.core.sql import sql_build

logger = logging.getLogger(__name__)

PgSchema = load_schema("pg")


class PG(BaseStructure):
    """Protein groups with per-run quantification."""

    _schema_class = PgSchema

    @classmethod
    def from_file(cls, path, **engine_kwargs) -> "PG":
        """Open a standalone pg.parquet file, guarding the pre-1.1 layout."""
        from qpx.version import check_pg_file_compatible

        check_pg_file_compatible(path)
        return super().from_file(path, **engine_kwargs)

    def by_protein(self, protein: str) -> "PG":
        """Filter by anchor protein: every group the protein leads, which can be several."""
        return self.filter(f"anchor_protein = '{_escape_sql_string(protein)}'")

    def by_group(self, accessions: list[str]) -> "PG":
        """Filter to the protein group whose membership is exactly ``accessions``, in any order."""
        members = ", ".join(f"'{_escape_sql_string(accession)}'" for accession in accessions)
        return self.filter(f"list_sort(pg_accessions) = list_sort([{members}])")

    def by_run(self, run_file_name: str) -> "PG":
        """Filter to protein groups whose grouped_runs contains the given run file."""
        return self.filter(f"list_contains(grouped_runs, '{_escape_sql_string(run_file_name)}')")

    def join(self, other: BaseStructure, on: str | None = None) -> "PG":
        """Join PG rows to Run by expanding grouped_runs membership."""
        if on is None and isinstance(other, Run):
            return self._join_list_membership(
                other,
                list_column="grouped_runs",
                value_column="run_file_name",
            )
        return super().join(other, on=on)

    def targets_only(self) -> "PG":
        """Filter to target proteins only (exclude decoys)."""
        return self.filter("is_decoy = false")

    def protein_intensities(self) -> QueryResult:
        """One row per (protein, label, intensity).

        Since QPX 1.1 the pg view is flattened, so ``label`` and ``intensity`` are
        scalar columns and no ``UNNEST`` is needed. ID-only rows (null label) are
        excluded.
        """
        stmt = sql_build(
            """SELECT anchor_protein, pg_accessions, gg_names, grouped_runs,
               global_qvalue, label, intensity
        FROM $src
        WHERE label IS NOT NULL""",
            src=self._query.source,
        )
        return QueryResult(self._engine.connection, stmt)
