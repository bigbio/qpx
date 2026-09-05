"""PSM data structure — Peptide Spectrum Matches."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.query import _escape_sql_string

PsmSchema = load_schema("psm")


class PSM(BaseStructure):
    """Peptide Spectrum Matches from identification searches."""

    _schema_class = PsmSchema

    def by_protein(self, protein: str) -> "PSM":
        """Filter PSMs by protein accession (searches within array)."""
        return self.filter(f"list_contains(protein_accessions, '{_escape_sql_string(protein)}')")

    def by_run(self, run_file_name: str) -> "PSM":
        """Filter PSMs by run file."""
        return self.filter(f"run_file_name = '{_escape_sql_string(run_file_name)}'")

    def targets_only(self) -> "PSM":
        """Filter to target PSMs only (exclude decoys)."""
        return self.filter("is_decoy = false")

    def with_feature(self) -> "PSM":
        """Filter to PSMs with a populated ``feature_id``."""
        return self.filter("feature_id IS NOT NULL")

    def without_feature(self) -> "PSM":
        """Filter to PSMs without a populated ``feature_id``.

        A missing link does not establish whether the precursor was quantified.
        OpenMS PSM-only exports also leave assigned PSMs' feature links null.
        """
        return self.filter("feature_id IS NULL")
