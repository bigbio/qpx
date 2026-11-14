"""FragPipe Feature adapter -- combined_ion.tsv to feature.parquet.

Loads a FragPipe ``combined_ion.tsv`` into DuckDB,
transforms feature rows into ``FeatureSchema``, and writes through ``FeatureWriter``.
"""

from __future__ import annotations

import logging
from typing import NamedTuple, Optional

import pandas as pd

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.fragpipe.constants import (
    is_decoy_accession,
    protein_accessions,
    to_modifications,
    to_proforma,
    unassigned_mass_shift,
)
from qpx.converters.fragpipe.isobaric import channel_intensities, ms1_additional_intensity, require_channel_columns
from qpx.converters.fragpipe.pg_adapter import normalize_run_name, require_label_free, resolve_experiment_runs
from qpx.converters.mappings import get_field_mappings
from qpx.converters.ptm import compute_precursor_mz
from qpx.converters.utils import safe_float
from qpx.core.sql import escape_path, sql_build
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

_PSM_LOOKUP_SQL = """
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
                $mc_expr AS missed_cleavages,
                $mass_expr AS calculated_mass,
                $delta_expr AS delta_mass
            FROM _fp_psm_lookup"""

# Derive field map from central YAML mappings
_FEATURE_MAP = get_field_mappings("fragpipe", "feature")
# IonQuant's per-experiment ``Match Type`` column.
_IDENTIFICATION_TYPES = {"MS/MS": "direct", "MBR": "transferred"}
_FEATURE_IDENTITY_COMPOSITE = (
    "quantification_unit_id",
    "peptidoform",
    "charge",
    "compensation_voltage",
)


class _PsmInfo(NamedTuple):
    """Fields a feature takes from its best PSM; tuples keep a million-PSM lookup small."""

    observed_mz: float | None
    calculated_mz: float | None
    pep: float | None
    is_decoy: bool
    scan: int
    ion_mobility: float | None
    missed_cleavages: int | None


def _psm_rank(pep: float | None, scan: int) -> tuple[float, int]:
    """Order PSMs best first: lowest PEP (a missing one last), then lowest scan."""
    return (pep if pep is not None else float("inf"), scan)


def _pg_proteins(accessions: list[str]) -> list[dict] | None:
    """pg_protein structs of *accessions*; FragPipe gives no positions, so those are null."""
    return [{"accession": acc, "start": None, "end": None, "pre": None, "post": None} for acc in accessions] or None


class FragPipeFeatureAdapter(BaseConverter):
    """Convert FragPipe ``combined_ion.tsv`` to ``feature.parquet``.

    Usage::

        with FragPipeFeatureAdapter() as adapter:
            adapter.convert(
                feature_path="combined_ion.tsv",
                output_path="feature.parquet",
            )
    """

    # Members of each protein group by lead accession, from combined_protein.tsv.
    _protein_groups: dict[str, list[str]] = {}
    # Sequence and every mapped protein of each converted peptidoform.
    _pepmap: dict[str, tuple[str, list[str]]] = {}
    # Channel column -> label of a one-plex ion.tsv; empty for label-free input.
    _channels: dict[str, str] = {}

    def load_channels(self, channels: dict[str, str]) -> None:
        """Read reporter-ion quantities from these columns of a one-plex ion.tsv."""
        self._channels = channels

    def load_protein_groups(self, protein_path: str) -> None:
        """Read the members of each protein group from FragPipe ``combined_protein.tsv``.

        An ion names its group by the lead ``Protein``; the group also holds the
        lead's ``Indistinguishable Proteins``. Features then carry the full
        membership that the pg view writes, so the two join; without this call a
        feature lists only its lead protein.
        """
        columns = ["Protein", "Indistinguishable Proteins"]
        proteins = pd.read_csv(protein_path, sep="\t", dtype=str, usecols=lambda col: col in columns)
        groups = (protein_accessions(lead, others) for lead, others in proteins.reindex(columns=columns).itertuples(index=False))
        self._protein_groups = {members[0]: members for members in groups if members}

    def convert(
        self,
        feature_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        psm_path: Optional[str] = None,
        experiment_to_runs: dict[str, list[str]] | None = None,
        experiment_annotation_path: Optional[str] = None,
        chunksize: int = 500_000,
        creator: str = "fragpipe",
    ) -> None:
        """Run the combined_ion.tsv -> feature.parquet conversion.

        Args:
            feature_path: Path to FragPipe ``combined_ion.tsv``.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF fallback for mapping experiments to raw runs;
                an experiment that no mapping resolves raises ``ValueError``.
            psm_path: Optional path to FragPipe ``psm.tsv`` for PSM-level lookups
                (mass error, PEP, scan, is_decoy).
            experiment_to_runs: Optional explicit mapping from each experiment
                intensity prefix to its member raw runs.
            experiment_annotation_path: Optional FragPipe
                ``experiment_annotation.tsv`` providing that mapping.
            chunksize: Records per batch; a row gives one record per
                experiment, so ``chunksize // experiments`` rows are read at a time.
            creator: Creator tag in Parquet metadata.
        """
        if not self._channels:
            require_label_free(sdrf_path)

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

        # Step 3: Detect format and experiment columns, and resolve each
        # experiment to its raw runs
        self._require_ion_format()
        if self._channels:
            require_channel_columns(self._channels, actual_cols, feature_path)
        experiments = list(experiment_to_runs) if self._channels else self._detect_experiment_columns()
        self.logger.info("Detected ion format, experiments: %s", experiments)
        experiment_runs = resolve_experiment_runs(experiments, experiment_to_runs, experiment_annotation_path, sdrf_path)

        # Step 4: Build PSM lookup if psm.tsv provided
        psm_lookup = self._build_psm_lookup(psm_path) if psm_path else {}
        multi_run = [experiment for experiment, runs in experiment_runs.items() if len(runs) > 1]
        if psm_lookup and multi_run:
            self.logger.warning("Skipping PSM enrichment for multi-run FragPipe experiment(s): %s", ", ".join(multi_run))
        self._pepmap = {}

        # Step 5: Stream and transform
        with FeatureWriter(
            output_path,
            creator=creator,
            identity_composite=_FEATURE_IDENTITY_COMPOSITE,
        ) as writer:
            for batch in self._query_batched("SELECT * FROM fragpipe_features", max(1, chunksize // len(experiments))):
                df = batch.to_pandas()
                records = self._transform_batch(df, experiment_runs, psm_lookup)
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
        count = self._conn.execute("SELECT COUNT(*) FROM fragpipe_features").fetchall()[0][0]
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
        if not experiments:
            raise ValueError(
                "FragPipe combined_ion.tsv has no '<experiment> Intensity' columns; a one-plex ion.tsv holds "
                "its reporter quantities in channel columns, which --annotation-file names"
            )
        return sorted(experiments)

    # ------------------------------------------------------------------
    # PSM lookup
    # ------------------------------------------------------------------

    def _build_psm_lookup(self, psm_path: str) -> dict[tuple, _PsmInfo]:
        """Build a lookup for feature-level PSM enrichment.

        Streams the needed columns of FragPipe ``psm.tsv`` and keeps the best PSM
        (lowest PEP, then lowest scan) providing mass error, PEP, scan, and
        decoy flag, whatever order DuckDB returns the rows in.

        Each PSM is registered under its exact
        ``(source_file, peptidoform, charge)`` key. The experiment-to-source
        mapping is resolved separately from FragPipe metadata; no run-agnostic
        precursor fallback is used because it can attach one run's PSM to every
        experiment in a combined results table.
        """
        lookup: dict[tuple, _PsmInfo] = {}

        self._conn.execute(
            sql_build(
                """CREATE TEMPORARY VIEW _fp_psm_lookup AS
            SELECT * FROM read_csv_auto('$path',
                delim='\t', header=true, auto_detect=true)""",
                path=escape_path(psm_path),
            )
        )

        # Detect columns
        actual_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='_fp_psm_lookup'"
            ).fetchall()
        }

        result = self._conn.execute(self._psm_lookup_query(actual_cols))
        count = 0
        while rows := result.fetchmany(100_000):
            count += len(rows)
            for values in rows:
                self._add_psm_lookup_row(lookup, values)

        self._conn.execute("DROP VIEW IF EXISTS _fp_psm_lookup")
        self.logger.info(
            "PSM lookup: %s unique keys from %s PSM rows of %s", format(len(lookup), ","), format(count, ","), psm_path
        )

        return lookup

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(
        self,
        df: pd.DataFrame,
        experiment_runs: dict[str, list[str]],
        psm_lookup: dict[tuple, _PsmInfo],
    ) -> list[dict]:
        # Pre-extract column arrays for faster per-row access than to_dict("records")
        col_arrays = {col: df[col].values for col in df.columns}
        rows = ({col: vals[i] for col, vals in col_arrays.items()} for i in range(len(df)))
        return self._transform_rows_with_skip(
            rows,
            lambda row: self._transform_row(row, experiment_runs, psm_lookup),
            "FragPipe feature",
        )

    def _transform_row(
        self,
        row,
        experiment_runs: dict[str, list[str]],
        psm_lookup: dict[tuple, _PsmInfo],
    ) -> list[dict]:
        """Transform a single row into one or more feature records.

        Expands into one record per experiment with non-zero intensity. The record
        is keyed by the experiment's first raw run, a member of its pg
        ``grouped_runs``; the experiment stays the ``quantification_unit_id``.
        """
        records: list[dict] = []
        common_fields, mz, protein_raw, mapped = self._feature_identity_fields(row)

        # Expand per experiment
        for experiment, runs in experiment_runs.items():
            intensities, additional_intensities = self._experiment_quantities(row, experiment)
            if intensities is None:
                continue

            # The identifying run, and so its PSM, is known only for a one-run
            # experiment whose feature was identified there by MS/MS rather than
            # transferred by match-between-runs. One combined Feature of a
            # multi-run experiment cannot represent one fraction's PSM.
            match_type = row.get(f"{experiment} Match Type")
            id_run = runs[0] if len(runs) == 1 and match_type != "MBR" else None
            psm = psm_lookup.get((id_run, common_fields["peptidoform"], str(common_fields["charge"])))

            psm_fields = self._feature_psm_fields(psm, mz, protein_raw)

            rec = {
                **common_fields,
                **psm_fields,
                "additional_scores": None,
                "predicted_rt": None,
                "run_file_name": runs[0],
                "quantification_unit_id": experiment,
                "cv_params": None,
                "rt": None,
                "intensities": intensities,
                "additional_intensities": additional_intensities,
                "pg_global_qvalue": None,
                "ion_mobility_start": None,
                "ion_mobility_stop": None,
                "gg_accessions": None,
                "id_run_file_name": id_run,
                "identification_type": _IDENTIFICATION_TYPES.get(match_type),
                "rt_start": None,
                "rt_stop": None,
            }
            records.append(rec)

        if records and mapped is not None:
            self._pepmap.setdefault(common_fields["peptidoform"], (common_fields["sequence"], mapped))
        return records

    def _experiment_quantities(self, row, experiment: str) -> tuple[list[dict] | None, list[dict] | None]:
        """The experiment's intensities, or its channels' with the MS1 intensity beside them."""
        if self._channels:
            intensities = channel_intensities(row, self._channels)
            return intensities, ms1_additional_intensity(row, intensities) if intensities else None
        intensity = safe_float(row.get(f"{experiment} Intensity")) or 0.0
        return ([{"label": "LFQ", "intensity": float(intensity)}] if intensity > 0 else None), None

    def write_pepmap(self, output_path: str, creator: str = "fragpipe") -> bool:
        """Write every protein each converted peptidoform maps to as the pepmap view.

        ``pg_accessions`` of the feature view holds the peptide's protein group;
        ``Protein`` and ``Mapped Proteins`` of combined_ion list every protein the
        peptide maps to. Returns whether a file was written: not for an ion table
        without ``Mapped Proteins`` or without features.
        """
        return self._write_pepmap(output_path, self._pepmap, creator)

    @staticmethod
    def _psm_lookup_query(actual_cols: set[str]) -> str:
        """Build the typed PSM enrichment query."""
        # Observed m/z: prefer calibrated
        obs_col = "Calibrated Observed M/Z" if "Calibrated Observed M/Z" in actual_cols else "Observed M/Z"
        calc_col = "Calculated M/Z"
        has_obs = obs_col in actual_cols
        has_calc = calc_col in actual_cols

        obs_expr = f'TRY_CAST("{obs_col}" AS DOUBLE)' if has_obs else "NULL"
        calc_expr = f'TRY_CAST("{calc_col}" AS DOUBLE)' if has_calc else "NULL"

        # PeptideProphet probability ("Probability" in FragPipe 22) → PEP = 1 - prob
        pp_col = next((col for col in ("PeptideProphet Probability", "Probability") if col in actual_cols), None)
        pp_expr = f'TRY_CAST("{pp_col}" AS DOUBLE)' if pp_col else "NULL"

        # Ion mobility
        has_im = "Ion Mobility" in actual_cols
        im_expr = 'TRY_CAST("Ion Mobility" AS DOUBLE)' if has_im else "NULL"

        # Missed cleavages
        has_mc = "Number of Missed Cleavages" in actual_cols
        mc_expr = 'TRY_CAST("Number of Missed Cleavages" AS INTEGER)' if has_mc else "NULL"

        # Mass shift of an open or mass-offset search
        has_shift = {"Calculated Peptide Mass", "Delta Mass"} <= actual_cols
        mass_expr = 'TRY_CAST("Calculated Peptide Mass" AS DOUBLE)' if has_shift else "NULL"
        delta_expr = 'TRY_CAST("Delta Mass" AS DOUBLE)' if has_shift else "NULL"

        sql = sql_build(
            _PSM_LOOKUP_SQL,
            obs_expr=obs_expr,
            calc_expr=calc_expr,
            pp_expr=pp_expr,
            im_expr=im_expr,
            mc_expr=mc_expr,
            mass_expr=mass_expr,
            delta_expr=delta_expr,
        )
        return sql

    @staticmethod
    def _psm_source_scan(spectrum) -> tuple[str, int]:
        """Run name and scan number from a FragPipe ``Spectrum`` value; the scan is 0 when it does not parse."""
        # Parse source file and scan from Spectrum column
        tokens = str(spectrum).split(".") if spectrum else []
        source_file = normalize_run_name(tokens[0]) if tokens else ""
        scan_number = 0
        if len(tokens) >= 2:
            try:
                scan_number = int(tokens[1].lstrip("0") or "0")
            except ValueError:
                pass
        return source_file, scan_number

    def _add_psm_lookup_row(self, lookup: dict, values: tuple) -> None:
        """Register the PSM for an exact source and precursor key if it beats the one there."""
        (
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
            calculated_mass,
            delta_mass,
        ) = values
        source_file, scan_number = self._psm_source_scan(spectrum)

        sequence = str(sequence) if sequence else ""
        assigned_mods = str(assigned_mods) if assigned_mods else ""
        # A PSM with an unassigned mass shift is another molecule than the ion.
        if unassigned_mass_shift(assigned_mods, sequence, calculated_mass, delta_mass) is not None:
            return
        peptidoform = to_proforma(assigned_mods, sequence)
        charge_str = str(charge) if charge is not None else "0"

        key = (source_file, peptidoform, charge_str)
        # PEP = 1 - PeptideProphet Probability
        pep = None
        if pp_prob is not None:
            pep = 1.0 - pp_prob if pp_prob <= 1.0 else pp_prob
        current = lookup.get(key)
        if current is None or _psm_rank(pep, scan_number) < _psm_rank(current.pep, current.scan):
            lookup[key] = _PsmInfo(
                observed_mz=float(obs_mz) if obs_mz is not None else None,
                calculated_mz=compute_precursor_mz(peptidoform, int(charge_str))
                or (float(calc_mz) if calc_mz is not None else None),
                pep=pep,
                is_decoy=is_decoy_accession(protein),
                scan=scan_number,
                ion_mobility=float(ion_mobility) if ion_mobility is not None else None,
                missed_cleavages=int(missed_cleavages_val) if missed_cleavages_val is not None else None,
            )

    def _feature_identity_fields(self, row: dict) -> tuple:
        """Parse the feature peptide and protein identity fields."""
        r = self._resolved  # shorthand for resolved column mappings

        # Sequence
        sequence = str(row.get(r.get("sequence", "Peptide Sequence"), ""))

        # Peptidoform -- build ProForma from sequence + Assigned Modifications
        mods_col = r.get("modifications", "Assigned Modifications")
        assigned_mods_raw = row.get(mods_col, "")
        assigned_mods_str = str(assigned_mods_raw) if pd.notna(assigned_mods_raw) and assigned_mods_raw else ""
        peptidoform = to_proforma(assigned_mods_str, sequence)

        # Protein group: the lead protein's group members; every protein the
        # peptide maps to (null without Mapped Proteins) decides uniqueness.
        protein_raw = row.get(r.get("pg_accessions", "Protein"))
        leads = protein_accessions(protein_raw)
        anchor_protein = leads[0] if leads else ""
        mapped = protein_accessions(protein_raw, row.get("Mapped Proteins")) if "Mapped Proteins" in row else None

        # Gene
        gene_raw = row.get(r.get("gg_names", "Gene"), "")
        gg_names = [str(gene_raw)] if pd.notna(gene_raw) and gene_raw else None

        # Modifications (reuse assigned_mods_str already extracted for peptidoform)
        modifications = None
        if assigned_mods_str:
            _, modifications = to_modifications(assigned_mods_str, sequence)

        # M/Z of the feature file is the ion's calculated m/z, used without a PSM
        mz = safe_float(row.get(r.get("calculated_mz", "M/Z")))
        compensation_voltage = safe_float(row.get(r.get("compensation_voltage", "Compensation Voltage"))) or 0.0

        charge = int(row.get(r.get("charge", "Charge"), 0))
        return (
            {
                "sequence": sequence,
                "peptidoform": peptidoform,
                "modifications": modifications,
                "charge": charge,
                "compensation_voltage": compensation_voltage,
                "pg_accessions": _pg_proteins(self._protein_groups.get(anchor_protein) or leads),
                "anchor_protein": anchor_protein,
                "unique": len(mapped) == 1 if mapped is not None else None,
                "gg_names": gg_names,
            },
            mz,
            protein_raw,
            mapped,
        )

    @staticmethod
    def _feature_psm_fields(psm: _PsmInfo | None, mz: float | None, protein_raw: str | None) -> dict:
        """Resolve spectrum fields from exact PSM enrichment, else the feature's calculated m/z."""
        if psm is None:
            # Without a PSM, is_decoy falls back to the decoy prefix on the RAW
            # protein field (before parsing strips e.g. the "rev_" prefix).
            psm = _PsmInfo(None, None, None, is_decoy_accession(protein_raw), 0, None, None)
        _calc = psm.calculated_mz
        _obs = psm.observed_mz
        return {
            "posterior_error_probability": psm.pep,
            "is_decoy": psm.is_decoy,
            "calculated_mz": _calc or mz,
            "observed_mz": _obs,
            "mass_error_ppm": 1e6 * (_obs - _calc) / _calc if _calc and _obs else None,
            "scan": [psm.scan] if psm.scan > 0 else [],
            "ion_mobility": psm.ion_mobility,
            "missed_cleavages": psm.missed_cleavages,
        }
