"""MaxQuant Feature adapter -- evidence.txt to feature.parquet.

Loads ``evidence.txt`` into DuckDB, optionally merges SDRF for sample
mapping and TMT channels, transforms into ``FeatureSchema``, and writes
through ``FeatureWriter``.

Key schema changes:
    - ``reference_file_name`` -> ``run_file_name``
    - ``precursor_charge`` (int32) -> ``charge`` (int16)
    - ``is_decoy`` (int32, 0/1) -> ``is_decoy`` (bool)
    - ``scan`` (string) -> ``scan`` (list<int32>)
    - ``start_ion_mobility``/``stop_ion_mobility`` -> ``ion_mobility_start``/``ion_mobility_stop``
    - Intensities struct: ``{sample_accession, channel, intensity}`` ->
      ``{label, intensity}``
"""

from __future__ import annotations

import logging
import re
from typing import Optional

import duckdb
import pandas as pd

from qpx.converters.base import resolve_columns
from qpx.converters.mappings import get_field_mappings
from qpx.converters.maxquant.base_adapter import MaxQuantBaseAdapter
from qpx.converters.maxquant.constants import (
    PROTON_MASS,
    TMT_LABEL_TO_MQ_COL,
    to_proforma,
)
from qpx.converters.ptm import from_proforma
from qpx.converters.utils import mq_flag_to_bool, safe_float, strip_uniprot_prefix
from qpx.core.sql import sql_build, validate_identifier
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

# Derive field map from central mappings
_FEATURE_MAP = get_field_mappings("maxquant", "feature")

# Producer-specific feature identity composite (bigbio/qpx#229). MaxQuant reports
# an RT window per feature, so its measured key uses rt_start/rt_stop (both in
# feature.yaml) rather than a single rt. Passed to FeatureWriter so feature_id
# hashes exactly these columns instead of the schema default.
_FEATURE_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name", "rt_start", "rt_stop", "scan")


class MaxQuantFeatureAdapter(MaxQuantBaseAdapter):
    """Convert MaxQuant ``evidence.txt`` to ``feature.parquet``.

    Usage::

        with MaxQuantFeatureAdapter() as adapter:
            adapter.convert(
                evidence_path="evidence.txt",
                output_path="feature.parquet",
                sdrf_path="sdrf.tsv",
            )
    """

    def convert(
        self,
        evidence_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        protein_groups_path: Optional[str] = None,
        chunksize: int = 500_000,
        creator: str = "maxquant",
        fixed_mod_only: bool = False,
    ) -> None:
        """Run the evidence.txt -> feature.parquet conversion.

        Args:
            evidence_path: Path to MaxQuant ``evidence.txt``.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF file for sample/channel mapping.
            protein_groups_path: Optional ``proteinGroups.txt`` for
                populating ``pg_global_qvalue`` and gene names on features.
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Track which dropped-TMT-channel warnings have already been emitted so
        # the message is logged once, not once per batch.
        self._warned_tmt_channels: set = set()

        # Step 1: Load evidence into DuckDB
        self._load_evidence(evidence_path)

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='evidence'"
            ).fetchall()
        }
        self._resolved = resolve_columns(_FEATURE_MAP, actual_cols)

        # Step 3: Load SDRF for sample mapping
        sample_map, experiment_type, tmt_channels = self._load_sdrf(sdrf_path)

        # Step 3b: Build lookup maps from proteinGroups.txt
        pg_maps = self._build_protein_group_maps(protein_groups_path)

        # Step 4: Stream and transform
        self.logger.info("Transforming MaxQuant features ...")
        evidence_query = self._deduplicated_evidence_query(actual_cols)

        with FeatureWriter(
            output_path,
            creator=creator,
            compression=self._compression,
            identity_composite=_FEATURE_IDENTITY_COMPOSITE,
        ) as writer:
            for batch in self._query_batched(evidence_query, chunksize):
                df = batch.to_pandas()
                records = self._transform_batch(
                    df, sample_map, experiment_type, tmt_channels, pg_maps, fixed_mod_only=fixed_mod_only
                )
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"MaxQuant feature conversion complete -> {output_path}")

    def _deduplicated_evidence_query(self, actual_cols: set[str]) -> str:
        """Collapse rows that share the feature identity, preferring identified evidence.

        The partition key mirrors ``_FEATURE_IDENTITY_COMPOSITE`` -- the columns
        ``feature_id`` is hashed from -- so rows resolving to the same
        ``feature_id`` collapse to one. The previous key partitioned on
        ``Calibrated retention time`` (which is *not* part of the identity), so
        identity-equal rows differing only in calibrated RT survived as
        duplicate primary keys. Only identity columns actually present are used;
        a missing column narrows the partition rather than degrading to a
        no-dedup ``SELECT * FROM evidence`` (which would let MBR duplicates
        through).
        """
        # Map each identity-composite field to its resolved MaxQuant column.
        # "Modified sequence" stands in for ``peptidoform`` (a deterministic
        # function of it); the remaining fields match the composite directly.
        identity_field_columns = [
            self._resolved.get("modified_sequence", "Modified sequence"),
            self._resolved.get("charge", "Charge"),
            self._resolved.get("run_file_name", "Raw file"),
            self._resolved.get("rt_start", "Calibrated retention time start"),
            self._resolved.get("rt_stop", "Calibrated retention time finish"),
            self._resolved.get("scan", "MS/MS scan number"),
        ]
        identity_columns = [column for column in identity_field_columns if column in actual_cols]
        if not identity_columns:
            return "SELECT * FROM evidence"

        partition_by = ", ".join(validate_identifier(column) for column in identity_columns)
        order_by = []
        if "Type" in actual_cols:
            type_column = validate_identifier("Type")
            order_by.append(
                sql_build(
                    "CASE WHEN UPPER(TRIM(CAST($type_column AS VARCHAR))) = 'MULTI-MATCH' THEN 1 ELSE 0 END",
                    type_column=type_column,
                )
            )
        pep_column = self._resolved.get("posterior_error_probability", "PEP")
        if pep_column in actual_cols:
            order_by.append(
                sql_build(
                    "TRY_CAST($pep_column AS DOUBLE) ASC NULLS LAST",
                    pep_column=validate_identifier(pep_column),
                )
            )
        if "id" in actual_cols:
            order_by.append(
                sql_build(
                    "TRY_CAST($id_column AS BIGINT) ASC NULLS LAST",
                    id_column=validate_identifier("id"),
                )
            )
        if not order_by:
            # No preference columns present; keep an arbitrary-but-stable row.
            order_by.append("1")

        return sql_build(
            """
            SELECT *
            FROM evidence
            QUALIFY ROW_NUMBER() OVER (
                PARTITION BY $partition_by
                ORDER BY $order_by
            ) = 1
            """,
            partition_by=partition_by,
            order_by=", ".join(order_by),
        )

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_evidence(self, path: str) -> None:
        """Load evidence.txt into DuckDB."""
        self._conn.execute(
            "CREATE TABLE evidence AS SELECT * FROM read_csv_auto($1, delim='\t', header=true, auto_detect=true, null_padding=true)",
            [path],
        )
        count = self._conn.execute("SELECT COUNT(*) FROM evidence").fetchone()[0]
        self.logger.info(f"Loaded {count:,} MaxQuant evidence rows")

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _build_protein_group_maps(self, protein_groups_path: Optional[str]) -> dict:
        """Build protein-group lookup maps from proteinGroups.txt.

        Returns a dict with two sub-maps:
            ``by_group``: full protein-group membership (a ``frozenset`` of
                normalised accessions) -> ``{"qvalue", "genes"}``. This is the
                primary key, consistent with the pg identity fix which keys on
                the whole ``pg_accessions`` membership rather than the leader.
            ``by_accession``: single accession -> ``{"qvalue", "genes"}``, but
                *only* for accessions that belong to exactly one protein group.
                Ambiguous razor accessions (shared across groups) are excluded
                so a feature cannot inherit another group's q-value/genes via
                its leader.
        """
        empty = {"by_group": {}, "by_accession": {}}
        if not protein_groups_path:
            return empty
        try:
            df = self._conn.execute(
                "SELECT * FROM read_csv_auto($1, delim='\t', header=true, auto_detect=true, null_padding=true)",
                [protein_groups_path],
            ).df()
            acc_col, qval_col, gene_col = self._detect_pg_columns(df)
            if not acc_col:
                return empty
            maps = self._build_pg_lookup(df, acc_col, qval_col, gene_col)
            self.logger.info(
                "Built protein group maps: %d group(s), %d unambiguous accession(s)",
                len(maps["by_group"]),
                len(maps["by_accession"]),
            )
            return maps
        except (FileNotFoundError, pd.errors.ParserError, KeyError, ValueError, duckdb.Error) as e:
            self.logger.warning("Could not build protein group maps: %s", e)
        return empty

    @staticmethod
    def _detect_pg_columns(df: pd.DataFrame) -> tuple:
        """Detect accession, q-value, and gene name columns in proteinGroups."""
        acc_col = next(
            (c for c in ["Protein IDs", "Majority protein IDs"] if c in df.columns),
            None,
        )
        qval_col = next((c for c in ["Q-value", "q-value", "Q.value"] if c in df.columns), None)
        gene_col = next((c for c in ["Gene names", "Gene Names"] if c in df.columns), None)
        return acc_col, qval_col, gene_col

    @staticmethod
    def _canonical_group_key(accessions) -> frozenset:
        """Order-independent key for a protein-group membership."""
        return frozenset(a.strip() for a in accessions if a and a.strip())

    @staticmethod
    def _parse_genes(gene_raw, fasta_raw) -> Optional[list[str]]:
        """Parse gene names from a proteinGroups row (Gene names, else Fasta headers)."""
        genes: list[str] | None = None
        if gene_raw is not None and pd.notna(gene_raw) and str(gene_raw).strip():
            genes = [g.strip() for g in str(gene_raw).split(";") if g.strip()] or None
        if not genes and fasta_raw is not None and pd.notna(fasta_raw) and str(fasta_raw).strip():
            genes = []
            for hdr in str(fasta_raw).split(";"):
                m = re.search(r"GN=([^\s;]+)", hdr)
                if m and m.group(1) not in genes:
                    genes.append(m.group(1))
            genes = genes or None
        return genes

    def _build_pg_lookup(self, df: pd.DataFrame, acc_col: str, qval_col: Optional[str], gene_col: Optional[str]) -> dict:
        """Build the ``by_group``/``by_accession`` lookup maps (see caller)."""
        by_group: dict[frozenset, dict] = {}
        acc_to_groups: dict[str, set] = {}
        acc_value: dict[str, dict] = {}

        acc_vals = df[acc_col].astype(str).tolist()
        qval_vals = pd.to_numeric(df[qval_col], errors="coerce").tolist() if qval_col else [None] * len(df)
        gene_vals = df[gene_col].tolist() if (gene_col and gene_col in df.columns) else [None] * len(df)
        fasta_vals = df["Fasta headers"].tolist() if "Fasta headers" in df.columns else [None] * len(df)

        for acc_raw, qval, gene_raw, fasta_raw in zip(acc_vals, qval_vals, gene_vals, fasta_vals):
            accs = [strip_uniprot_prefix(a.strip()) for a in acc_raw.split(";") if a.strip()]
            if not accs:
                continue
            qvalue = float(qval) if qval is not None and pd.notna(qval) else None
            genes = self._parse_genes(gene_raw, fasta_raw)
            value = {"qvalue": qvalue, "genes": genes}
            key = self._canonical_group_key(accs)
            by_group.setdefault(key, value)
            for acc in accs:
                acc_to_groups.setdefault(acc, set()).add(key)
                acc_value.setdefault(acc, value)

        # Keep only accessions that map to a single protein group for the
        # leader fallback; ambiguous accessions would otherwise collide.
        by_accession = {acc: acc_value[acc] for acc, keys in acc_to_groups.items() if len(keys) == 1}
        return {"by_group": by_group, "by_accession": by_accession}

    def _lookup_pg_metadata(self, pg_maps: Optional[dict], pg_acc_list: list[str], anchor_protein: str) -> tuple:
        """Resolve ``(pg_global_qvalue, genes)`` for a feature's protein group.

        Keyed on the full protein-group membership first; only falls back to the
        leader accession when that leader is unambiguous (belongs to one group).
        """
        if not pg_maps:
            return None, None
        entry = pg_maps.get("by_group", {}).get(self._canonical_group_key(pg_acc_list))
        if entry is None and anchor_protein:
            entry = pg_maps.get("by_accession", {}).get(anchor_protein)
        if entry is None:
            return None, None
        return entry.get("qvalue"), entry.get("genes")

    def _detect_ri_offset(self, columns) -> int:
        """Detect whether reporter intensity columns are 0-indexed or 1-indexed.

        MaxQuant versions and experiment configurations differ:
        - 0-indexed: ``Reporter intensity 0`` … ``Reporter intensity N-1``
        - 1-indexed: ``Reporter intensity 1`` … ``Reporter intensity N``

        Returns 0 or 1 accordingly.
        """
        col_set = set(columns)
        if "Reporter intensity 0" in col_set:
            return 0
        return 1

    def _warn_dropped_tmt_channels(self, columns, experiment_type: str, tmt_channels: list[str], ri_offset: int) -> None:
        """Warn when the SDRF plex declares channels absent from evidence.txt.

        A ``Reporter intensity N`` column missing for a declared channel means
        that channel is silently dropped from every feature; surface it once.
        """
        if getattr(self, "_warned_tmt_channels", None) is None:  # pragma: no cover - defensive
            self._warned_tmt_channels = set()
        exp_type = re.sub(r"\d+", "", experiment_type or "").upper()
        if exp_type not in ("TMT", "ITRAQ") or not tmt_channels:
            return
        col_set = set(columns)
        dropped: list[tuple[str, str]] = []
        for seq_idx, channel_name in enumerate(tmt_channels):
            col_idx = TMT_LABEL_TO_MQ_COL.get(str(channel_name).strip().upper())
            if col_idx is None:
                col_idx = seq_idx
            col_name = f"Reporter intensity {col_idx + ri_offset}"
            if col_name not in col_set:
                dropped.append((str(channel_name), col_name))
        if not dropped:
            return
        key = tuple(sorted(ch for ch, _ in dropped))
        if key in self._warned_tmt_channels:
            return
        self._warned_tmt_channels.add(key)
        logger.warning(
            "SDRF declares %d TMT/iTRAQ channel(s) with no matching reporter intensity "
            "column in evidence.txt; these channels are dropped from all features: %s",
            len(dropped),
            ", ".join(f"{ch} (expected {col})" for ch, col in dropped),
        )

    def _transform_batch(
        self,
        df: pd.DataFrame,
        sample_map: dict,
        experiment_type: str,
        tmt_channels: list[str],
        pg_maps: dict | None = None,
        fixed_mod_only: bool = False,
    ) -> list[dict]:
        """Transform a batch of evidence.txt rows into QPX feature records."""
        records: list[dict] = []
        # Detect reporter intensity column indexing once per batch
        ri_offset = self._detect_ri_offset(df.columns)
        # Warn (once) about SDRF channels with no reporter column in evidence.txt
        self._warn_dropped_tmt_channels(df.columns, experiment_type, tmt_channels, ri_offset)
        # Pre-extract column arrays for faster per-row access than to_dict("records")
        col_arrays = {col: df[col].values for col in df.columns}
        n_rows = len(df)
        for i in range(n_rows):
            row = {col: vals[i] for col, vals in col_arrays.items()}
            rec = self._transform_row(
                row,
                sample_map,
                experiment_type,
                tmt_channels,
                pg_maps,
                ri_offset,
                fixed_mod_only=fixed_mod_only,
            )
            if rec:
                records.append(rec)
        return records

    _FIXED_MODS_ALLOWED = frozenset({"unmodified", "carbamidomethyl (c)"})

    def _transform_row(
        self,
        row,
        sample_map: dict,
        experiment_type: str,
        tmt_channels: list[str],
        pg_maps: dict | None = None,
        ri_offset: int = 0,
        fixed_mod_only: bool = False,
    ) -> Optional[dict]:
        """Transform a single evidence.txt row."""
        r = self._resolved  # shorthand for resolved column mappings

        if fixed_mod_only:
            mods_raw = str(row.get(r.get("modifications", "Modifications"), "Unmodified")).strip()
            if mods_raw.lower() not in self._FIXED_MODS_ALLOWED:
                return None

        sequence = str(row.get(r.get("sequence", "Sequence"), ""))
        peptidoform = to_proforma(
            str(row.get(r.get("modified_sequence", "Modified sequence"), "")),
        )
        _, modifications = from_proforma(peptidoform, sequence) if peptidoform else (None, None)
        charge = int(row.get(r.get("charge", "Charge"), 0))
        run_file_name = str(row.get(r.get("run_file_name", "Raw file"), ""))

        # is_decoy (bool)
        is_decoy = mq_flag_to_bool(row.get(r.get("is_decoy", "Reverse"), ""))

        # Scan (list<int32>)
        scan_raw = row.get(r.get("scan", "MS/MS scan number"))
        scan = []
        if pd.notna(scan_raw):
            try:
                scan = [int(scan_raw)]
            except (ValueError, TypeError):
                pass

        # Calculated m/z from neutral mass and charge
        neutral_mass = safe_float(row.get(r.get("mass", "Mass")))
        if neutral_mass and charge:
            calculated_mz = (neutral_mass + charge * PROTON_MASS) / charge
        else:
            calculated_mz = 0.0

        # m/z
        observed_mz = safe_float(row.get(r.get("observed_mz", "m/z"))) or 0.0

        # RT (minutes → seconds)
        rt = safe_float(row.get(r.get("rt", "Calibrated retention time")))
        if rt is not None:
            rt *= 60.0
        rt_start = safe_float(row.get(r.get("rt_start", "Calibrated retention time start")))
        if rt_start is not None:
            rt_start *= 60.0
        rt_stop = safe_float(row.get(r.get("rt_stop", "Calibrated retention time finish")))
        if rt_stop is not None:
            rt_stop *= 60.0

        # PEP
        pep = safe_float(row.get(r.get("posterior_error_probability", "PEP")))

        # Ion mobility
        ion_mobility = safe_float(row.get(r.get("ion_mobility", "1/K0")))

        # Protein accessions (list of pg_protein structs)
        pg_acc_raw = str(row.get(r.get("pg_accessions", "Leading proteins"), ""))
        pg_acc_list = [self._norm_uniprot_id(a) for a in pg_acc_raw.split(";")] if pg_acc_raw else []
        anchor_protein = self._norm_uniprot_id(str(row.get(r.get("anchor_protein", "Leading razor protein"), "")))
        if not anchor_protein and pg_acc_list:
            anchor_protein = pg_acc_list[0]
        pg_accessions = [{"accession": acc, "start": None, "end": None, "pre": None, "post": None} for acc in pg_acc_list] or None

        # Unique peptide indicator
        unique = len(pg_acc_list) <= 1

        # Protein-group q-value and genes, keyed on the full membership.
        pg_qvalue, pg_genes = self._lookup_pg_metadata(pg_maps, pg_acc_list, anchor_protein)

        # Gene names
        gg_raw = row.get(r.get("gg_names", "Gene names"))
        gg_names = str(gg_raw).split(";") if pd.notna(gg_raw) and gg_raw else None
        if not gg_names:
            gg_names = pg_genes

        # Mass error (ppm) — direct column from evidence.txt
        mass_error_ppm = safe_float(row.get(r.get("mass_error_ppm", "Mass error [ppm]")))

        # Missed cleavages (dedicated field from evidence.txt)
        mc_raw = row.get(r.get("missed_cleavages", "Missed cleavages"))
        missed_cleavages = int(mc_raw) if pd.notna(mc_raw) else None

        # MBR detection: Type == "MULTI-MATCH" means transferred (MBR) feature
        evidence_type = str(row.get("Type", "")).strip().upper()
        id_run_file_name = None if evidence_type == "MULTI-MATCH" else run_file_name

        # Build intensities (new schema: {label, intensity})
        intensities, additional_intensities = self._build_intensities(
            row, run_file_name, sample_map, experiment_type, tmt_channels, ri_offset
        )

        # Additional scores
        additional_scores = []
        andromeda = safe_float(row.get(r.get("andromeda_score", "Score")))
        if andromeda is not None:
            additional_scores.append(
                {
                    "score_name": "andromeda_score",
                    "score_value": andromeda,
                    "higher_better": True,
                }
            )
        delta = safe_float(row.get(r.get("andromeda_delta_score", "Delta score")))
        if delta is not None:
            additional_scores.append(
                {
                    "score_name": "andromeda_delta_score",
                    "score_value": delta,
                    "higher_better": True,
                }
            )

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
            "cv_params": None,
            "scan": scan,
            "rt": rt,
            "ion_mobility": ion_mobility,
            "missed_cleavages": missed_cleavages,
            "intensities": intensities or None,
            "additional_intensities": additional_intensities,
            "pg_accessions": pg_accessions,
            "anchor_protein": anchor_protein,
            "unique": unique,
            "pg_global_qvalue": pg_qvalue,
            "ion_mobility_start": None,
            "ion_mobility_stop": None,
            "gg_accessions": gg_names,
            "gg_names": gg_names,
            "id_run_file_name": id_run_file_name,
            "rt_start": rt_start,
            "rt_stop": rt_stop,
        }

    def _build_intensities(
        self,
        row,
        run_file_name: str,
        sample_map: dict,
        experiment_type: str,
        tmt_channels: list[str],
        ri_offset: int = 0,
    ) -> tuple[list[dict], Optional[list[dict]]]:
        """Build intensities and additional_intensities from evidence row."""
        intensities: list[dict] = []
        additional_intensities: list[dict] = []

        exp_type = re.sub(r"\d+", "", experiment_type).upper()

        if exp_type in ("TMT", "ITRAQ") and tmt_channels:
            # TMT/iTRAQ: one intensity per channel.
            # Use the explicit label→column-index map so that N/C isomer pairs
            # (e.g. TMT127N at col 1, TMT127C at col 2) are always read from
            # their correct MaxQuant 0-indexed reporter intensity column,
            # regardless of how tmt_channels happens to be ordered.
            for seq_idx, channel_name in enumerate(tmt_channels):
                col_idx = TMT_LABEL_TO_MQ_COL.get(str(channel_name).strip().upper())
                if col_idx is None:
                    # Unknown label (e.g. bare TMT6-plex "TMT128" without N/C suffix).
                    # Fall back to sequential position within the channel list, which is
                    # correct for single-mass plexes where no N/C isomers exist.
                    logger.debug(
                        "Unknown TMT/iTRAQ channel label %r — using sequential index %d",
                        channel_name,
                        seq_idx,
                    )
                    col_idx = seq_idx
                col_name = f"Reporter intensity {col_idx + ri_offset}"
                val = row.get(col_name)
                if pd.notna(val) and float(val) > 0:
                    intensities.append({"label": channel_name, "intensity": float(val)})

                    # Corrected intensity
                    corr_col = col_name.replace("Reporter intensity", "Reporter intensity corrected")
                    corr_val = row.get(corr_col)
                    if pd.notna(corr_val):
                        additional_intensities.append(
                            {
                                "label": channel_name,
                                "intensities": [
                                    {
                                        "intensity_name": "corrected_reporter_intensity",
                                        "intensity_value": float(corr_val),
                                    }
                                ],
                            }
                        )
        else:
            # LFQ: prefer MaxQuant-normalised "LFQ intensity" when available,
            # fall back to raw "Intensity" (e.g. for datasets without LFQ algo).
            lfq_val = safe_float(row.get("LFQ intensity"))
            raw_val = safe_float(row.get(_FEATURE_MAP["intensity"][0]))
            intensity_val = lfq_val if (lfq_val is not None and lfq_val > 0) else (raw_val or 0.0)
            label = "LFQ"
            intensities.append({"label": label, "intensity": float(intensity_val)})

        return intensities, additional_intensities or None
