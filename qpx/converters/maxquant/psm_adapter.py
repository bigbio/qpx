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
from pathlib import Path
from typing import Optional

import pandas as pd

from qpx.converters.base import resolve_columns
from qpx.converters.localization import localize_modifications
from qpx.converters.mappings import get_field_mappings
from qpx.converters.maxquant.base_adapter import MaxQuantBaseAdapter
from qpx.converters.maxquant.constants import (
    MQ_LABEL_LETTERS,
    FixedMods,
    localization_columns,
    localization_sites,
    mq_calculated_mz,
    mq_observed_mz,
    summary_label_mods,
    to_proforma,
)
from qpx.converters.maxquant.evidence_fields import peak_annotations
from qpx.converters.maxquant.silac import multiplet_labels
from qpx.converters.ptm import from_proforma
from qpx.converters.utils import mass_error_ppm, mq_flag_to_bool, safe_float
from qpx.writers.psm import PsmWriter

logger = logging.getLogger(__name__)

# Derive expected columns from central mappings
_PSM_MAP = get_field_mappings("maxquant", "psm")

# Optional spectral data columns
_MQ_PSM_SPECTRAL_COLS = [
    "Matches",
    "Fragmentation",
    "Mass analyzer",
    "Type",
]


class MaxQuantPsmAdapter(MaxQuantBaseAdapter):
    """Convert MaxQuant ``msms.txt`` to ``psm.parquet``.

    ``psm.feature_id`` is filled from ``Evidence ID`` when the feature view's
    ids are passed to :meth:`convert`. A SILAC PSM carries the label of its
    ``Labeling state`` and that label's modifications.

    Usage::

        with MaxQuantPsmAdapter() as adapter:
            adapter.convert(
                msms_path="msms.txt",
                output_path="psm.parquet",
            )
    """

    # Site-probability columns of the current batch.
    _localization_columns: list[tuple[str, str]] = []

    _evidence_feature_ids: dict[int, int] = {}
    # Raw file -> (label, label modifications) of each SILAC labeling state.
    _states: dict[str, tuple[tuple[str, FixedMods], ...]] = {}

    def convert(
        self,
        msms_path: str,
        output_path: str,
        chunksize: int = 500_000,
        spectral_data: bool = False,
        creator: str = "maxquant",
        sdrf_path: Optional[str] = None,
        evidence_feature_ids: Optional[dict[int, int]] = None,
    ) -> None:
        """Run the msms.txt -> psm.parquet conversion.

        Args:
            msms_path: Path to MaxQuant ``msms.txt``.
            output_path: Destination Parquet path.
            chunksize: Rows per batch.
            spectral_data: Whether to include spectral arrays.
            creator: Creator tag in Parquet metadata.
            sdrf_path: Optional SDRF file declaring fixed modifications and labels.
            evidence_feature_ids: evidence.txt ``id`` -> feature_id of the
                feature view written alongside, linking each PSM through its
                ``Evidence ID``.
        """
        self._fixed = self._fixed_mods(sdrf_path, msms_path)
        self._evidence_feature_ids = evidence_feature_ids or {}
        # Step 1: Load msms.txt into DuckDB
        self._load_msms(msms_path)

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0]
            for c in self._conn.execute("SELECT column_name FROM information_schema.columns WHERE table_name='msms'").fetchall()
        }
        self._resolved = resolve_columns(_PSM_MAP, actual_cols)
        self._states = self._labeling_states(msms_path, actual_cols, sdrf_path)

        # Step 3: Stream and transform
        self.logger.info("Transforming MaxQuant PSMs ...")

        with PsmWriter(output_path, creator=creator, compression=self._compression) as writer:
            for batch in self._query_batched("SELECT * FROM msms", chunksize):
                df = batch.to_pandas()
                records = self._transform_batch(df, spectral_data)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)
            if spectral_data:
                writer.add_footer_metadata({b"spectrum_content": b"matched_peaks"})

        self.logger.info(f"MaxQuant PSM conversion complete -> {output_path}")

    def _labeling_states(self, msms_path: str, columns: set[str], sdrf_path: Optional[str]) -> dict:
        """``Raw file`` -> ``(label, label modifications)`` of each SILAC labeling state.

        SILAC ``msms.txt`` has a ``Labeling state`` column; the labels of each
        state are in the ``summary.txt`` MaxQuant writes beside it.
        """
        run_labels, experiment_type, _ = self._load_sdrf(sdrf_path)
        if "Labeling state" not in columns:
            if experiment_type == "SILAC":
                raise ValueError("The SDRF declares SILAC labels, but msms.txt has no Labeling state column")
            return {}
        summary_path = Path(msms_path).with_name("summary.txt")
        label_mods = summary_label_mods(summary_path)
        if not label_mods:
            raise ValueError(f"SILAC msms.txt needs MaxQuant's summary.txt beside it for the label definitions: {summary_path}")
        # A SILAC pair is "L"/"H", a triplet "L"/"M"/"H".
        letters = {run: "LH" if len(states) == 2 else MQ_LABEL_LETTERS[: len(states)] for run, states in label_mods.items()}
        labels = multiplet_labels("".join(letters.values()), run_labels, experiment_type)
        return {
            run: tuple((labels[letter], mods) for letter, mods in zip(letters[run], states)) for run, states in label_mods.items()
        }

    def _psm_label(self, row: dict) -> tuple[Optional[str], FixedMods]:
        """The SILAC label and label modifications of a PSM; none when its labeling state is unknown (-1)."""
        state = row.get("Labeling state")
        if not self._states or pd.isna(state) or state < 0:
            return None, ()
        return self._states[str(row.get(self._resolved.get("run_file_name", "Raw file")))][int(state)]

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_msms(self, path: str) -> None:
        """Load msms.txt into DuckDB."""
        self._conn.execute(
            "CREATE TABLE msms AS SELECT * FROM read_csv_auto($1, header=true, auto_detect=true, delim='\t', null_padding=true)",
            [path],
        )
        count = self._conn.execute("SELECT COUNT(*) FROM msms").fetchall()[0][0]
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
        self._localization_columns = localization_columns(df.columns)
        for i in range(n_rows):
            row = {col: vals[i] for col, vals in col_arrays.items()}
            rec = self._transform_row(row, spectral_data)
            if rec:
                records.append(rec)
        return records

    def _transform_row(self, row, spectral_data: bool) -> Optional[dict]:
        """Transform a single msms.txt row into a QPX PSM record."""
        r = self._resolved  # shorthand for resolved column mappings

        identity = self._psm_identity(row)
        if identity is None:
            return None
        sequence, peptidoform, modifications, candidates, label = identity
        charge_raw = row.get(r.get("charge", "Charge"), 0)
        charge = int(charge_raw) if pd.notna(charge_raw) else 0

        # is_decoy (bool) -- MaxQuant uses '+' for Reverse
        is_decoy = mq_flag_to_bool(row.get(r.get("is_decoy", "Reverse"), ""))

        spectrum_fields = self._psm_spectrum_fields(row, charge, peptidoform)

        # Protein accessions — normalise sp|/tr| prefixes for consistency
        proteins_cell = row.get("Proteins")
        protein_accessions = (
            [self._norm_uniprot_id(a.strip()) for a in str(proteins_cell).split(";") if a.strip()]
            if pd.notna(proteins_cell) and proteins_cell is not None
            else []
        )

        # Ion mobility
        ion_mobility = safe_float(row.get("1/K0"))

        # Missed cleavages (dedicated field from msms.txt)
        mc_raw = row.get(r.get("missed_cleavages", "Missed cleavages"))
        missed_cleavages = int(mc_raw) if pd.notna(mc_raw) else None

        additional_scores, cv_params = self._psm_annotations(row)

        mz_array, intensity_array = self._spectral_arrays(row, spectral_data)
        ion_types, charges = (
            peak_annotations(row.get("Matches") if pd.notna(row.get("Matches")) else None, len(mz_array))
            if mz_array
            else (None, None)
        )

        evidence_id = row.get("Evidence ID")
        return {
            **spectrum_fields,
            "feature_id": self._evidence_feature_ids.get(int(evidence_id)) if pd.notna(evidence_id) else None,
            "label": label,
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "modification_candidates": candidates,
            "charge": charge,
            "is_decoy": is_decoy,
            "additional_scores": additional_scores or None,
            "predicted_rt": None,
            "cv_params": cv_params or None,
            "ion_mobility": ion_mobility,
            "missed_cleavages": missed_cleavages,
            "protein_accessions": protein_accessions or None,
            "mz_array": mz_array,
            "intensity_array": intensity_array,
            "charge_array": charges,
            "ion_type_array": ion_types,
            "ion_mobility_array": None,
        }

    def _psm_identity(self, row: dict) -> tuple | None:
        """Parse the peptide identity and its site localization."""
        r = self._resolved
        sequence = str(row.get(r.get("sequence", "Sequence"), ""))
        label, label_mods = self._psm_label(row)
        peptidoform = to_proforma(str(row.get(r.get("modified_sequence", "Modified sequence"), "")), self._fixed, label_mods)

        # Mirror the feature adapter: an empty/missing "Modified sequence"
        # yields no peptidoform, and the tuple-unpack must degrade to
        # (None, None) rather than a bare None (which raises TypeError and
        # aborts the whole PSM view).
        peptidoform, modifications = from_proforma(peptidoform, sequence) if peptidoform else (None, None)
        if peptidoform is None:
            # PsmSchema requires a peptidoform; a row with an empty/missing
            # "Modified sequence" is skipped rather than aborting the whole view.
            self.logger.debug("Skipping PSM row with empty modified sequence (sequence=%r)", sequence)
            return None
        sites = localization_sites(row, self._localization_columns, sequence)
        modifications, candidates = localize_modifications(modifications, sequence, sites)
        return sequence, peptidoform, modifications, candidates, label

    def _psm_spectrum_fields(self, row: dict, charge: int, peptidoform: str) -> dict:
        """Extract PSM spectrum, m/z, and retention-time fields."""
        r = self._resolved
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

        calculated_mz = mq_calculated_mz(peptidoform, charge, safe_float(row.get("Mass")))
        observed_mz = mq_observed_mz(
            safe_float(row.get(r.get("observed_mz", "m/z"))), safe_float(row.get(r.get("mass_error_ppm", "Mass error [ppm]")))
        )

        # RT (minutes → seconds)
        rt = safe_float(row.get(r.get("rt", "Retention time")))
        if rt is not None:
            rt *= 60.0

        # PEP
        pep = safe_float(row.get(r.get("posterior_error_probability", "PEP")))
        return {
            "scan": scan,
            "run_file_name": run_file_name,
            "observed_mz": observed_mz,
            "calculated_mz": calculated_mz,
            "mass_error_ppm": mass_error_ppm(calculated_mz, observed_mz),
            "rt": rt,
            "posterior_error_probability": pep,
        }

    def _psm_annotations(self, row: dict) -> tuple[list[dict], list[dict]]:
        """Collect PSM scores and spectrum annotations."""
        r = self._resolved
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
        # MaxQuant writes 0 when it did not compute the fraction.
        pif = safe_float(row.get("PIF"))
        if pif:
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
        return additional_scores, cv_params

    @staticmethod
    def _spectral_arrays(row: dict, spectral_data: bool) -> tuple:
        """m/z and intensity arrays from the ``Masses``/``Intensities`` columns when ``spectral_data`` is set.

        An array that does not parse as numbers stays None.
        """
        mz_array = None
        intensity_array = None
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
        return mz_array, intensity_array
