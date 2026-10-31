"""FragPipe PSM adapter -- psm.tsv to psm.parquet.

Loads a FragPipe ``psm.tsv`` into DuckDB, transforms via SQL into
``PsmSchema``, and streams results through ``PsmWriter``.

Key schema changes:
    - ``precursor_charge`` (int32) -> ``charge`` (int16)
    - ``is_decoy`` (int32, 0/1) -> ``is_decoy`` (bool)
    - ``scan`` (string) -> ``scan`` (list<int32>)
    - ``source_file`` -> ``run_file_name``
"""

from __future__ import annotations

import logging
from typing import Optional, Tuple

import pandas as pd

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.fragpipe.constants import (
    is_decoy_accession,
    name_glycan,
    place_mass_shift,
    protein_accessions,
    to_proforma,
    unassigned_mass_shift,
)
from qpx.converters.mappings import get_field_mappings
from qpx.converters.ptm import compute_precursor_mz, from_proforma
from qpx.converters.utils import safe_float
from qpx.core.scores import normalize_score_name
from qpx.core.sql import escape_path, sql_build
from qpx.writers.psm import PsmWriter

logger = logging.getLogger(__name__)

# Derive field map from central YAML mappings
_PSM_MAP = get_field_mappings("fragpipe", "psm")


def _parse_spectrum_id(identifier: str) -> Tuple[str, int]:
    """Parse FragPipe Spectrum identifier into (source_file, scan_number).

    Format: ``<source_file>.<scan>.<scan>.<charge>``
    """
    tokens = str(identifier).split(".")
    source_file = tokens[0] if tokens else ""
    scan_number = 0
    if len(tokens) >= 2:
        try:
            scan_number = int(tokens[1].lstrip("0") or "0")
        except ValueError:
            pass
    return source_file, scan_number


class FragPipePsmAdapter(BaseConverter):
    """Convert FragPipe ``psm.tsv`` to ``psm.parquet``.

    Usage::

        with FragPipePsmAdapter() as adapter:
            adapter.convert(
                psm_path="psm.tsv",
                output_path="psm.parquet",
            )
    """

    def convert(
        self,
        psm_path: str,
        output_path: str,
        chunksize: int = 500_000,
        creator: str = "fragpipe",
    ) -> None:
        """Run the psm.tsv -> psm.parquet conversion.

        Args:
            psm_path: Path to FragPipe ``psm.tsv``.
            output_path: Destination Parquet path.
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load psm.tsv into DuckDB
        self._load_psm_tsv(psm_path)

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='fragpipe_psms'"
            ).fetchall()
        }
        self._resolved = resolve_columns(_PSM_MAP, actual_cols)

        # Step 3: Stream and transform
        self.logger.info("Transforming FragPipe PSMs ...")

        with PsmWriter(output_path, creator=creator, compression=self._compression) as writer:
            for batch in self._query_batched("SELECT * FROM fragpipe_psms", chunksize):
                df = batch.to_pandas()
                records = self._transform_batch(df)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"FragPipe PSM conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_psm_tsv(self, path: str) -> None:
        """Load FragPipe psm.tsv into DuckDB."""
        self._conn.execute(
            sql_build(
                """CREATE TABLE fragpipe_psms AS
            SELECT * FROM read_csv_auto('$path',
                delim='\t', header=true, auto_detect=true)""",
                path=escape_path(path),
            )
        )
        count = self._conn.execute("SELECT COUNT(*) FROM fragpipe_psms").fetchall()[0][0]
        self.logger.info(f"Loaded {count:,} FragPipe PSM rows")

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(self, df: pd.DataFrame) -> list[dict]:
        # Pre-extract column arrays for faster per-row access than to_dict("records")
        col_arrays = {col: df[col].values for col in df.columns}
        rows = ({col: vals[i] for col, vals in col_arrays.items()} for i in range(len(df)))
        return self._transform_rows_with_skip(rows, self._transform_row, "FragPipe PSM")

    def _transform_row(self, row) -> Optional[dict]:
        """Transform a single FragPipe psm.tsv row."""
        r = self._resolved  # shorthand for resolved column mappings

        sequence = str(row.get(r.get("sequence", "Peptide"), ""))
        charge = int(row.get(r.get("charge", "Charge"), 0))

        spectrum_fields = self._psm_spectrum_fields(row)
        observed_mz = spectrum_fields["observed_mz"]

        # Protein accessions: the razor protein and every other protein the peptide maps to
        protein_raw = row.get(r.get("pg_accessions", "Protein"))
        accessions = protein_accessions(protein_raw, row.get("Mapped Proteins"))

        # Is decoy (bool) -- decoy prefix on any accession in the RAW protein field
        # (rev_/REV_/DECOY_/decoy_), kept consistent with the feature/pg adapters.
        is_decoy = is_decoy_accession(protein_raw)

        peptidoform, modifications, candidates = self._psm_modifications(row, sequence)

        # Calculated M/Z lacks the mass shift, and a glycan PTM-Shepherd only added to Assigned Modifications
        calculated_mz = compute_precursor_mz(peptidoform, charge) or spectrum_fields["calculated_mz"]

        # Mass error (ppm) — compute from observed and calculated m/z
        if observed_mz and calculated_mz:
            mass_error_ppm = 1e6 * (observed_mz - calculated_mz) / calculated_mz
        else:
            mass_error_ppm = None

        # Ion mobility
        ion_mobility = safe_float(row.get("Ion Mobility"))

        additional_scores, missed_cleavages, cv_params = self._psm_annotations(row)

        return {
            **spectrum_fields,
            "calculated_mz": calculated_mz,
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "modification_candidates": candidates,
            "charge": charge,
            "is_decoy": is_decoy,
            "mass_error_ppm": mass_error_ppm,
            "additional_scores": additional_scores or None,
            "predicted_rt": None,
            "cv_params": cv_params or None,
            "ion_mobility": ion_mobility,
            "missed_cleavages": missed_cleavages,
            "protein_accessions": accessions or None,
            "mz_array": None,
            "intensity_array": None,
            "charge_array": None,
            "ion_type_array": None,
            "ion_mobility_array": None,
        }

    @staticmethod
    def _psm_modifications(row: dict, sequence: str) -> tuple:
        """Peptidoform, modifications and candidates from Assigned Modifications and any unassigned mass shift."""
        assigned_mods_raw = row.get("Assigned Modifications")
        assigned_mods = str(assigned_mods_raw) if pd.notna(assigned_mods_raw) and assigned_mods_raw else ""
        peptidoform = to_proforma(assigned_mods, sequence)
        candidates = None
        shift = unassigned_mass_shift(
            assigned_mods,
            sequence,
            safe_float(row.get("Calculated Peptide Mass")),
            safe_float(row.get("Delta Mass")),
        )
        if shift is not None:
            peptidoform, candidates = place_mass_shift(assigned_mods, sequence, shift, row.get("MSFragger Localization"))
        _, modifications = from_proforma(peptidoform, sequence, meta=None)
        return peptidoform, name_glycan(modifications, row.get("Total Glycan Composition")), candidates

    def _psm_spectrum_fields(self, row: dict) -> dict:
        """Extract source, spectrum, and precursor confidence fields."""
        r = self._resolved
        # Parse Spectrum column for source file and scan
        spectrum_raw = str(row.get("Spectrum", ""))
        source_file, scan_number = _parse_spectrum_id(spectrum_raw)

        # run_file_name is the source file stem
        run_file_name = source_file

        # scan (list<int32>)
        scan = [scan_number] if scan_number > 0 else []

        # m/z
        observed_mz = safe_float(row.get("Calibrated Observed M/Z", row.get(r.get("observed_mz", "Observed M/Z"))))
        calculated_mz = safe_float(row.get(r.get("calculated_mz", "Calculated M/Z")))

        # RT (FragPipe Retention is already in seconds)
        rt = safe_float(row.get(r.get("rt", "Retention")))

        # PEP (PeptideProphet Probability, "Probability" in FragPipe 22; convert to 1-prob)
        pp_prob = safe_float(row.get("PeptideProphet Probability", row.get("Probability")))
        pep = None
        if pp_prob is not None:
            pep = 1.0 - pp_prob if pp_prob <= 1.0 else pp_prob
        return {
            "run_file_name": run_file_name,
            "scan": scan,
            "observed_mz": observed_mz,
            "calculated_mz": calculated_mz,
            "rt": rt,
            "posterior_error_probability": pep,
        }

    @staticmethod
    def _psm_annotations(row: dict) -> tuple:
        """Collect search scores, cleavage count, and enzymatic termini."""
        # Additional scores
        additional_scores = []
        hyperscore = safe_float(row.get("Hyperscore"))
        if hyperscore is not None:
            additional_scores.append(
                {
                    "score_name": normalize_score_name("MSFragger:Hyperscore"),
                    "score_value": hyperscore,
                    "higher_better": True,
                }
            )
        expectation = safe_float(row.get("Expectation"))
        if expectation is not None:
            additional_scores.append(
                {
                    "score_name": normalize_score_name("MSFragger:Expectation"),
                    "score_value": expectation,
                    "higher_better": False,
                }
            )

        # Purity: the precursor's share of the ion abundance in the isolation window
        purity = safe_float(row.get("Purity"))
        if purity is not None:
            additional_scores.append(
                {
                    "score_name": "parent_ion_fraction",
                    "score_value": purity,
                    "higher_better": True,
                }
            )

        # Glycan assignment of PTM-Shepherd
        for column, score_name, higher_better in (
            ("Glycan Score", "glycan_score", True),
            ("Glycan q-value", "glycan_qvalue", False),
        ):
            value = safe_float(row.get(column))
            if value is not None:
                additional_scores.append({"score_name": score_name, "score_value": value, "higher_better": higher_better})

        # Missed cleavages (dedicated field)
        missed = row.get("Number of Missed Cleavages")
        missed_cleavages = int(missed) if pd.notna(missed) else None

        # CV params (enzymatic termini)
        cv_params = []
        termini = row.get("Number of Enzymatic Termini")
        if pd.notna(termini):
            cv_params.append(
                {
                    "cv_name": "number_of_enzymatic_termini",
                    "cv_value": str(int(termini)),
                }
            )
        return additional_scores, missed_cleavages, cv_params
