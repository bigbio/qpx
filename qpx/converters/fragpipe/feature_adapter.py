"""FragPipe Feature adapter -- combined_ion.tsv to feature.parquet.

Loads a FragPipe ``combined_ion.tsv`` into DuckDB,
transforms feature rows into ``FeatureSchema``, and writes through ``FeatureWriter``.
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.fragpipe.constants import is_decoy_accession, to_modifications, to_proforma
from qpx.converters.mappings import get_field_mappings
from qpx.converters.utils import parse_uniprot_id, safe_float
from qpx.core.sql import escape_path, sql_build
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

# Derive field map from central YAML mappings
_FEATURE_MAP = get_field_mappings("fragpipe", "feature")
_FEATURE_IDENTITY_COMPOSITE = (
    "quantification_unit_id",
    "peptidoform",
    "charge",
    "compensation_voltage",
)


def _extract_pg_proteins(
    protein_str: str,
    start: int | None = None,
    end: int | None = None,
) -> list[dict]:
    """Extract pg_protein structs from a FragPipe Protein field.

    Handles formats like ``sp|P12345|PROT_HUMAN`` or plain accession ``P12345``.
    Returns list of {accession, start, end, pre, post}. When the tool does not
    provide positional info, pass None and they will be stored as null.
    """
    if not protein_str:
        return []
    result = []
    for part in protein_str.split(","):
        part = part.strip()
        if not part:
            continue
        acc = parse_uniprot_id(part)[0]
        if acc:
            result.append(
                {
                    "accession": acc,
                    "start": start,
                    "end": end,
                    "pre": None,
                    "post": None,
                }
            )
    return result


class FragPipeFeatureAdapter(BaseConverter):
    """Convert FragPipe ``combined_ion.tsv`` to ``feature.parquet``.

    Usage::

        with FragPipeFeatureAdapter() as adapter:
            adapter.convert(
                feature_path="combined_ion.tsv",
                output_path="feature.parquet",
            )
    """

    def convert(
        self,
        feature_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        psm_path: Optional[str] = None,
        chunksize: int = 500_000,
        creator: str = "fragpipe",
    ) -> None:
        """Run the combined_ion.tsv -> feature.parquet conversion.

        Args:
            feature_path: Path to FragPipe ``combined_ion.tsv``.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF file (not used in current impl).
            psm_path: Optional path to FragPipe ``psm.tsv`` for PSM-level lookups
                (mass error, PEP, scan, is_decoy).
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load feature file into DuckDB
        self._load_feature_file(feature_path)

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='fragpipe_features'"
            ).fetchall()
        }
        self._resolved = resolve_columns(_FEATURE_MAP, actual_cols)

        # Step 3: Detect format and experiment columns
        self._require_ion_format()
        experiments = self._detect_experiment_columns()
        self.logger.info("Detected combined_ion format, experiments: %s", experiments)

        # Step 4: Build PSM lookup if psm.tsv provided
        psm_lookup = self._build_psm_lookup(psm_path) if psm_path else {}

        # Step 5: Stream and transform
        with FeatureWriter(
            output_path,
            creator=creator,
            compression=self._compression,
            identity_composite=_FEATURE_IDENTITY_COMPOSITE,
        ) as writer:
            for batch in self._query_batched("SELECT * FROM fragpipe_features", chunksize):
                df = batch.to_pandas()
                records = self._transform_batch(df, experiments, psm_lookup)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"FragPipe feature conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_feature_file(self, path: str) -> None:
        """Load combined_ion.tsv into DuckDB."""
        self._conn.execute(
            sql_build(
                """CREATE TABLE fragpipe_features AS
            SELECT * FROM read_csv_auto('$path',
                delim='\t', header=true, auto_detect=true)""",
                path=escape_path(path),
            )
        )
        count = self._conn.execute("SELECT COUNT(*) FROM fragpipe_features").fetchone()[0]
        self.logger.info(f"Loaded {count:,} FragPipe feature rows")

    def _require_ion_format(self) -> None:
        """Reject peptide-level input whose quantities cannot be split by charge."""
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name='fragpipe_features'"
        ).fetchall()
        col_names = {c[0] for c in cols}
        if "Charge" in col_names:
            return
        if "Charges" in col_names:
            raise ValueError(
                "FragPipe combined_peptide.tsv is peptide-level: its intensity cannot be assigned to individual "
                "charges without duplicating quantity. Use combined_ion.tsv for Feature conversion."
            )
        raise ValueError("FragPipe Feature input must be combined_ion.tsv with a Charge column")

    def _detect_experiment_columns(self) -> list[str]:
        """Detect per-experiment intensity columns.

        FragPipe feature files have columns like ``<experiment> Intensity``.
        ``<experiment> MaxLFQ Intensity`` is a derived value for the same
        experiment, not a separate experiment.
        """
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name='fragpipe_features'"
        ).fetchall()
        col_names = [c[0] for c in cols]

        experiments = set()
        for col in col_names:
            if col.endswith(" Intensity") and not col.endswith(" MaxLFQ Intensity"):
                experiments.add(col[: -len(" Intensity")])

        return sorted(experiments)

    # ------------------------------------------------------------------
    # PSM lookup
    # ------------------------------------------------------------------

    def _build_psm_lookup(self, psm_path: str) -> dict[tuple, dict]:
        """Build a lookup for feature-level PSM enrichment.

        Loads FragPipe ``psm.tsv`` into a temporary DuckDB table and extracts the
        first-occurrence PSM providing mass error, PEP, scan, and decoy flag.

        Each PSM is registered under two keys so that the feature side (which is
        keyed by *experiment*, not raw file) can always resolve a match:

        * ``(source_file, peptidoform, charge)`` — exact, used when the FragPipe
          experiment name equals the raw-file stem (non-fractionated case);
        * ``(peptidoform, charge)`` — fallback, used when experiment != raw file
          (fractionated runs), where the 3-tuple would otherwise always miss.

        The 2- and 3-element tuples never collide, so a single dict holds both.
        """
        lookup: dict[tuple, dict] = {}

        try:
            self._conn.execute(
                sql_build(
                    """CREATE TEMPORARY TABLE _fp_psm_lookup AS
                SELECT * FROM read_csv_auto('$path',
                    delim='\t', header=true, auto_detect=true)""",
                    path=escape_path(psm_path),
                )
            )
            count = self._conn.execute("SELECT COUNT(*) FROM _fp_psm_lookup").fetchone()[0]
            self.logger.info(f"PSM lookup: loaded {count:,} PSM rows from {psm_path}")

            # Detect columns
            actual_cols = {
                c[0]
                for c in self._conn.execute(
                    "SELECT column_name FROM information_schema.columns WHERE table_name='_fp_psm_lookup'"
                ).fetchall()
            }

            # Observed m/z: prefer calibrated
            obs_col = "Calibrated Observed M/Z" if "Calibrated Observed M/Z" in actual_cols else "Observed M/Z"
            calc_col = "Calculated M/Z"
            has_obs = obs_col in actual_cols
            has_calc = calc_col in actual_cols

            obs_expr = f'TRY_CAST("{obs_col}" AS DOUBLE)' if has_obs else "NULL"
            calc_expr = f'TRY_CAST("{calc_col}" AS DOUBLE)' if has_calc else "NULL"

            # PeptideProphet Probability → PEP = 1 - prob
            has_pp = "PeptideProphet Probability" in actual_cols
            pp_expr = 'TRY_CAST("PeptideProphet Probability" AS DOUBLE)' if has_pp else "NULL"

            # Ion mobility
            has_im = "Ion Mobility" in actual_cols
            im_expr = 'TRY_CAST("Ion Mobility" AS DOUBLE)' if has_im else "NULL"

            # Missed cleavages
            has_mc = "Number of Missed Cleavages" in actual_cols
            mc_expr = 'TRY_CAST("Number of Missed Cleavages" AS INTEGER)' if has_mc else "NULL"

            sql = sql_build(
                """
                SELECT
                    CAST("Spectrum" AS VARCHAR) AS spectrum,
                    CAST("Peptide" AS VARCHAR) AS sequence,
                    COALESCE(CAST("Assigned Modifications" AS VARCHAR), '') AS assigned_mods,
                    CAST("Charge" AS INTEGER) AS charge,
                    $obs_expr AS obs_mz,
                    $calc_expr AS calc_mz,
                    $pp_expr AS pp_prob,
                    CAST("Protein" AS VARCHAR) AS protein,
                    $im_expr AS ion_mobility,
                    $mc_expr AS missed_cleavages
                FROM _fp_psm_lookup""",
                obs_expr=obs_expr,
                calc_expr=calc_expr,
                pp_expr=pp_expr,
                im_expr=im_expr,
                mc_expr=mc_expr,
            )
            rows = self._conn.execute(sql).fetchall()

            for (
                spectrum,
                sequence,
                assigned_mods,
                charge,
                obs_mz,
                calc_mz,
                pp_prob,
                protein,
                ion_mobility,
                missed_cleavages_val,
            ) in rows:
                # Parse source file and scan from Spectrum column
                tokens = str(spectrum).split(".") if spectrum else []
                source_file = tokens[0] if tokens else ""
                scan_number = 0
                if len(tokens) >= 2:
                    try:
                        scan_number = int(tokens[1].lstrip("0") or "0")
                    except ValueError:
                        pass

                sequence = str(sequence) if sequence else ""
                assigned_mods = str(assigned_mods) if assigned_mods else ""
                peptidoform = to_proforma(assigned_mods, sequence)
                charge_str = str(charge) if charge is not None else "0"

                key = (source_file, peptidoform, charge_str)
                fallback_key = (peptidoform, charge_str)
                if key not in lookup:
                    # PEP = 1 - PeptideProphet Probability
                    pep = None
                    if pp_prob is not None:
                        pep = 1.0 - pp_prob if pp_prob <= 1.0 else pp_prob

                    is_decoy = is_decoy_accession(protein)

                    info = {
                        "observed_mz": float(obs_mz) if obs_mz is not None else None,
                        "calculated_mz": (float(calc_mz) if calc_mz is not None else None),
                        "pep": pep,
                        "is_decoy": is_decoy,
                        "scan": [scan_number] if scan_number > 0 else [],
                        "ion_mobility": (float(ion_mobility) if ion_mobility is not None else None),
                        "missed_cleavages": (int(missed_cleavages_val) if missed_cleavages_val is not None else None),
                    }
                    lookup[key] = info
                    # Register the run-agnostic fallback (first occurrence wins) so
                    # the experiment-keyed feature query still resolves when the
                    # experiment name differs from the raw-file stem.
                    lookup.setdefault(fallback_key, info)

            # Clean up temporary table
            self._conn.execute("DROP TABLE IF EXISTS _fp_psm_lookup")
            self.logger.info(f"PSM lookup: {len(lookup):,} unique keys built")

        except Exception as e:
            self.logger.warning(f"Could not build PSM lookup from {psm_path}: {e}")
            try:
                self._conn.execute("DROP TABLE IF EXISTS _fp_psm_lookup")
            except Exception:
                self.logger.debug("Failed to drop _fp_psm_lookup table", exc_info=True)

        return lookup

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(
        self,
        df: pd.DataFrame,
        experiments: list[str],
        psm_lookup: dict[tuple, dict] | None = None,
    ) -> list[dict]:
        records: list[dict] = []
        skipped = 0
        # Pre-extract column arrays for faster per-row access than to_dict("records")
        col_arrays = {col: df[col].values for col in df.columns}
        n_rows = len(df)
        for i in range(n_rows):
            try:
                row = {col: vals[i] for col, vals in col_arrays.items()}
                recs = self._transform_row(row, experiments, psm_lookup)
                records.extend(recs)
            except Exception as e:
                skipped += 1
                self.logger.debug(f"Skipping FragPipe feature row: {e}")
        if skipped:
            total = skipped + len(records)
            self.logger.warning(
                "Skipped %d / %d rows (%.1f%%) in batch",
                skipped,
                total,
                100 * skipped / total if total else 0,
            )
        return records

    def _transform_row(
        self,
        row,
        experiments: list[str],
        psm_lookup: dict[tuple, dict] | None = None,
    ) -> list[dict]:
        """Transform a single row into one or more feature records.

        Expands into one record per experiment with non-zero intensity.
        """
        records: list[dict] = []
        r = self._resolved  # shorthand for resolved column mappings

        # Sequence
        sequence = str(row.get(r.get("sequence", "Peptide Sequence"), ""))

        # Peptidoform -- build ProForma from sequence + Assigned Modifications
        mods_col = r.get("modifications", "Assigned Modifications")
        assigned_mods_raw = row.get(mods_col, "")
        assigned_mods_str = str(assigned_mods_raw) if pd.notna(assigned_mods_raw) and assigned_mods_raw else ""
        peptidoform = to_proforma(assigned_mods_str, sequence)

        # Protein mapping (schema: list<pg_protein>; FragPipe does not provide start/end)
        protein_raw = str(row.get(r.get("pg_accessions", "Protein"), ""))
        pg_accessions = _extract_pg_proteins(protein_raw, start=None, end=None) or None
        anchor_protein = pg_accessions[0]["accession"] if pg_accessions else ""

        # Gene
        gene_raw = row.get(r.get("gg_names", "Gene"), "")
        gg_names = [str(gene_raw)] if pd.notna(gene_raw) and gene_raw else None

        # Modifications (reuse assigned_mods_str already extracted for peptidoform)
        modifications = None
        if assigned_mods_str:
            _, modifications = to_modifications(assigned_mods_str, sequence)

        # M/Z (from feature file — used as fallback)
        mz = safe_float(row.get(r.get("observed_mz", "M/Z"))) or 0.0
        compensation_voltage = safe_float(row.get(r.get("compensation_voltage", "Compensation Voltage"))) or 0.0

        charge = int(row.get(r.get("charge", "Charge"), 0))

        # Expand per experiment
        for experiment in experiments:
            int_col = f"{experiment} Intensity"
            intensity_val = safe_float(row.get(int_col)) or 0.0
            if intensity_val <= 0:
                continue

            intensities = [{"label": "LFQ", "intensity": float(intensity_val)}]

            # PSM lookup: enrich with mass error, PEP, scan, decoy flag. Try the
            # exact (experiment==raw-file) key first, then the run-agnostic
            # (peptidoform, charge) fallback so fractionated experiments still hit.
            psm_info = {}
            if psm_lookup:
                psm_info = psm_lookup.get((experiment, peptidoform, str(charge)))
                if not psm_info:
                    psm_info = psm_lookup.get((peptidoform, str(charge)), {})

            _calc = psm_info.get("calculated_mz")
            _obs = psm_info.get("observed_mz")
            mass_error_ppm = 1e6 * (_obs - _calc) / _calc if _calc and _obs else None

            # Is decoy: prefer PSM lookup; fall back to the decoy prefix on the RAW
            # protein field (before parse_uniprot_id strips e.g. the "rev_" prefix).
            _is_decoy_fallback = is_decoy_accession(protein_raw)

            rec = {
                "sequence": sequence,
                "peptidoform": peptidoform,
                "modifications": modifications,
                "charge": charge,
                "posterior_error_probability": psm_info.get("pep"),
                "is_decoy": psm_info.get("is_decoy", _is_decoy_fallback),
                "calculated_mz": _calc or float(mz),
                "observed_mz": _obs or float(mz),
                "mass_error_ppm": mass_error_ppm,
                "additional_scores": None,
                "predicted_rt": None,
                "run_file_name": experiment,
                "quantification_unit_id": experiment,
                "cv_params": None,
                "scan": psm_info.get("scan", []),
                "rt": None,
                "ion_mobility": psm_info.get("ion_mobility"),
                "compensation_voltage": compensation_voltage,
                "missed_cleavages": psm_info.get("missed_cleavages"),
                "intensities": intensities,
                "additional_intensities": None,
                "pg_accessions": pg_accessions,
                "anchor_protein": anchor_protein,
                "unique": len(pg_accessions) <= 1 if pg_accessions else True,
                "pg_global_qvalue": None,
                "ion_mobility_start": None,
                "ion_mobility_stop": None,
                "gg_accessions": None,
                "gg_names": gg_names,
                "id_run_file_name": experiment,
                "rt_start": None,
                "rt_stop": None,
            }
            records.append(rec)

        return records
