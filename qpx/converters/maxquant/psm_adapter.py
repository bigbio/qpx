"""MaxQuant PSM adapter -- msms.txt to psm.parquet.

Loads ``msms.txt`` into DuckDB, transforms via SQL into ``PsmSchema``,
and streams results through ``PsmWriter``.

Key schema changes:
    - ``reference_file_name`` -> ``run_file_name``  (Raw file column)
    - ``precursor_charge`` (int32) -> ``charge`` (int16)
    - ``is_decoy`` (int32, 0/1) -> ``is_decoy`` (bool)
    - ``scan`` (string) -> ``scan`` (list<int32>)
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

from qpx.converters.base import resolve_columns
from qpx.converters.mappings import get_field_mappings
from qpx.converters.maxquant.base_adapter import MaxQuantBaseAdapter
from qpx.converters.maxquant.constants import PROTON_MASS, parse_phospho_probabilities, to_proforma
from qpx.converters.ptm import from_proforma
from qpx.converters.utils import mq_flag_to_bool, safe_float
from qpx.writers.psm import PsmWriter

logger = logging.getLogger(__name__)

# Derive expected columns from central mappings
_PSM_MAP = get_field_mappings("maxquant", "psm")

# MaxQuant columns to read from msms.txt (derived from constants)
_MQ_PSM_USECOLS = list({candidates[0] for candidates in _PSM_MAP.values()}) + [
    "Proteins",  # protein accessions (not in FIELD_MAPPINGS)
    "Mass",  # neutral mass for calculated_mz computation
    "PIF",  # parent ion fraction score
    "Masses",  # spectral masses
    "Intensities",  # spectral intensities
    "Number of matches",
    "1/K0",  # ion mobility
    "Phospho (STY) Probabilities",  # per-site phospho localization
    "Phospho (STY) Score Diffs",  # per-site phospho score differences
]

# Optional spectral data columns
_MQ_PSM_SPECTRAL_COLS = [
    "Matches",
    "Fragmentation",
    "Mass analyzer",
    "Type",
]


class MaxQuantPsmAdapter(MaxQuantBaseAdapter):
    """Convert MaxQuant ``msms.txt`` to ``psm.parquet``.

    Usage::

        with MaxQuantPsmAdapter() as adapter:
            adapter.convert(
                msms_path="msms.txt",
                output_path="psm.parquet",
            )
    """

    def convert(
        self,
        msms_path: str,
        output_path: str,
        chunksize: int = 500_000,
        spectral_data: bool = False,
        creator: str = "maxquant",
    ) -> None:
        """Run the msms.txt -> psm.parquet conversion.

        Args:
            msms_path: Path to MaxQuant ``msms.txt``.
            output_path: Destination Parquet path.
            chunksize: Rows per batch.
            spectral_data: Whether to include spectral arrays.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load msms.txt into DuckDB
        self._load_msms(msms_path)

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0]
            for c in self._conn.execute("SELECT column_name FROM information_schema.columns WHERE table_name='msms'").fetchall()
        }
        self._resolved = resolve_columns(_PSM_MAP, actual_cols)

        # Step 3: Stream and transform
        self.logger.info("Transforming MaxQuant PSMs ...")

        with PsmWriter(output_path, creator=creator, compression=self._compression) as writer:
            for batch in self._query_batched("SELECT * FROM msms", chunksize):
                df = batch.to_pandas()
                records = self._transform_batch(df, spectral_data)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"MaxQuant PSM conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_msms(self, path: str) -> None:
        """Load msms.txt into DuckDB."""
        self._conn.execute(
            "CREATE TABLE msms AS SELECT * FROM read_csv_auto($1, header=true, auto_detect=true, delim='\t', null_padding=true)",
            [path],
        )
        count = self._conn.execute("SELECT COUNT(*) FROM msms").fetchone()[0]
        self.logger.info(f"Loaded {count:,} MaxQuant PSM rows from msms.txt")

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(self, df: pd.DataFrame, spectral_data: bool) -> list[dict]:
        """Transform a batch of msms.txt rows into QPX PSM records."""
        records: list[dict] = []
        # Pre-extract column arrays for faster per-row access than to_dict("records")
        col_arrays = {col: df[col].values for col in df.columns}
        n_rows = len(df)
        for i in range(n_rows):
            row = {col: vals[i] for col, vals in col_arrays.items()}
            rec = self._transform_row(row, spectral_data)
            if rec:
                records.append(rec)
        return records

    def _transform_row(self, row, spectral_data: bool) -> Optional[dict]:
        """Transform a single msms.txt row into a QPX PSM record."""
        r = self._resolved  # shorthand for resolved column mappings

        sequence = str(row.get(r.get("sequence", "Sequence"), ""))
        peptidoform = to_proforma(
            str(row.get(r.get("modified_sequence", "Modified sequence"), "")),
        )

        # Parse phospho site localization probabilities (if available)
        phospho_raw = row.get("Phospho (STY) Probabilities")
        site_scores = None
        if pd.notna(phospho_raw) and phospho_raw:
            site_scores = parse_phospho_probabilities(str(phospho_raw))

        # Mirror the feature adapter: an empty/missing "Modified sequence"
        # yields no peptidoform, and the tuple-unpack must degrade to
        # (None, None) rather than a bare None (which raises TypeError and
        # aborts the whole PSM view).
        peptidoform, modifications = (
            from_proforma(
                peptidoform,
                sequence,
                site_scores=site_scores,
            )
            if peptidoform
            else (None, None)
        )
        if peptidoform is None:
            # PsmSchema requires a peptidoform; a row with an empty/missing
            # "Modified sequence" is skipped rather than aborting the whole view.
            self.logger.debug("Skipping PSM row with empty modified sequence (sequence=%r)", sequence)
            return None
        charge_raw = row.get(r.get("charge", "Charge"), 0)
        charge = int(charge_raw) if pd.notna(charge_raw) else 0

        # is_decoy (bool) -- MaxQuant uses '+' for Reverse
        is_decoy = mq_flag_to_bool(row.get(r.get("is_decoy", "Reverse"), ""))

        # Scan (list<int32>)
        scan_raw = row.get(r.get("scan", "Scan number"))
        scan = []
        if pd.notna(scan_raw):
            try:
                scan = [int(scan_raw)]
            except (ValueError, TypeError):
                pass

        # Run file name (strip extension)
        run_file_name = str(row.get(r.get("run_file_name", "Raw file"), ""))

        # m/z values
        observed_mz = safe_float(row.get(r.get("observed_mz", "m/z"))) or 0.0

        # Calculated m/z from neutral mass and charge
        neutral_mass = safe_float(row.get("Mass"))
        if neutral_mass and charge:
            calculated_mz = (neutral_mass + charge * PROTON_MASS) / charge
        else:
            calculated_mz = 0.0

        # RT (minutes → seconds)
        rt = safe_float(row.get(r.get("rt", "Retention time")))
        if rt is not None:
            rt *= 60.0

        # PEP
        pep = safe_float(row.get(r.get("posterior_error_probability", "PEP")))

        # Protein accessions — normalise sp|/tr| prefixes for consistency
        proteins_cell = row.get("Proteins")
        protein_accessions = (
            [self._norm_uniprot_id(a.strip()) for a in str(proteins_cell).split(";") if a.strip()]
            if pd.notna(proteins_cell) and proteins_cell is not None
            else []
        )

        # Mass error (ppm) — direct column from msms.txt
        mass_error_ppm = safe_float(row.get(r.get("mass_error_ppm", "Mass error [ppm]")))

        # Ion mobility
        ion_mobility = safe_float(row.get("1/K0"))

        # Missed cleavages (dedicated field from msms.txt)
        mc_raw = row.get(r.get("missed_cleavages", "Missed cleavages"))
        missed_cleavages = int(mc_raw) if pd.notna(mc_raw) else None

        # Additional scores
        additional_scores = []
        andromeda_score = safe_float(row.get(r.get("andromeda_score", "Score")))
        if andromeda_score is not None:
            additional_scores.append(
                {
                    "score_name": "andromeda_score",
                    "score_value": andromeda_score,
                    "higher_better": True,
                }
            )
        delta_score = safe_float(row.get(r.get("andromeda_delta_score", "Delta score")))
        if delta_score is not None:
            additional_scores.append(
                {
                    "score_name": "andromeda_delta_score",
                    "score_value": delta_score,
                    "higher_better": True,
                }
            )
        pif = safe_float(row.get("PIF"))
        if pif is not None:
            additional_scores.append(
                {
                    "score_name": "parent_ion_fraction",
                    "score_value": pif,
                    "higher_better": True,
                }
            )

        # CV params (from Fragmentation, Mass analyzer, Type)
        cv_params = []
        for cv_col in ["Fragmentation", "Mass analyzer", "Type"]:
            val = row.get(cv_col)
            if pd.notna(val) and val:
                cv_params.append({"cv_name": cv_col.lower(), "cv_value": str(val)})

        # Spectral arrays
        mz_array = None
        intensity_array = None
        charge_array = None
        ion_type_array = None
        if spectral_data:
            masses_raw = row.get("Masses")
            if pd.notna(masses_raw) and masses_raw:
                try:
                    mz_array = [float(x) for x in str(masses_raw).split(";") if x.strip()]
                except ValueError:
                    pass

            intensities_raw = row.get("Intensities")
            if pd.notna(intensities_raw) and intensities_raw:
                try:
                    intensity_array = [float(x) for x in str(intensities_raw).split(";") if x.strip()]
                except ValueError:
                    pass

        return {
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "charge": charge,
            "posterior_error_probability": pep,
            "is_decoy": is_decoy,
            "calculated_mz": calculated_mz,
            "observed_mz": observed_mz,
            "mass_error_ppm": mass_error_ppm,
            "additional_scores": additional_scores or None,
            "predicted_rt": None,
            "run_file_name": run_file_name,
            "cv_params": cv_params or None,
            "scan": scan,
            "rt": rt,
            "ion_mobility": ion_mobility,
            "missed_cleavages": missed_cleavages,
            "protein_accessions": protein_accessions or None,
            "mz_array": mz_array,
            "intensity_array": intensity_array,
            "charge_array": charge_array,
            "ion_type_array": ion_type_array,
            "ion_mobility_array": None,
        }
