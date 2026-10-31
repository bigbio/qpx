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
from typing import Optional

import duckdb
import pandas as pd

from qpx.converters.base import resolve_columns
from qpx.converters.localization import localize_modifications
from qpx.converters.mappings import get_field_mappings
from qpx.converters.maxquant.base_adapter import MaxQuantBaseAdapter
from qpx.converters.maxquant.constants import (
    localization_columns,
    localization_sites,
    mq_calculated_mz,
    mq_observed_mz,
    to_proforma,
)
from qpx.converters.maxquant.protein_fields import group_ids, parse_genes
from qpx.converters.maxquant.silac import (
    identified_labels,
    identified_states,
    multiplet_intensities,
    multiplet_labels,
    ratio_letters,
)
from qpx.converters.ptm import from_proforma
from qpx.converters.utils import canonical_accession, mass_error_ppm, mq_flag_to_bool, safe_float
from qpx.core.files import run_file_stem
from qpx.core.sql import sql_build, validate_identifier
from qpx.writers.feature import FeatureWriter, derive_feature_ids

logger = logging.getLogger(__name__)

# Derive field map from central mappings
_FEATURE_MAP = get_field_mappings("maxquant", "feature")

_EVIDENCE_DEDUP_SQL = """
            SELECT *$evidence_ids
            FROM evidence
            QUALIFY ROW_NUMBER() OVER (
                PARTITION BY $partition_by
                ORDER BY $order_by
            ) = 1
            ORDER BY $partition_by
            """

# Evidence ids collapsed into each kept evidence row.
_EVIDENCE_IDS = "__evidence_ids"

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

    def __init__(
        self,
        duckdb_memory: str | None = None,
        duckdb_threads: int | None = None,
        conn: duckdb.DuckDBPyConnection | None = None,
        compression: str = "zstd",
    ):
        super().__init__(
            duckdb_memory=duckdb_memory,
            duckdb_threads=duckdb_threads,
            conn=conn,
            compression=compression,
        )
        self._warned_tmt_channels: set[str] = set()
        # Canonical labels the SDRF assigns a sample, per run file stem.
        self._run_labels: dict[str, set[str]] = {}
        # evidence.txt ``id`` -> feature_id, filled when ``convert`` links evidence.
        self.evidence_feature_ids: dict[int, int] = {}
        # SILAC: evidence ``id`` -> labeling states identified by MS/MS.
        self._states_by_evidence: dict[int, list[int]] = {}
        # The batch's ``<mod> Probabilities`` columns and their modification names.
        self._localization_columns: list[tuple[str, str]] = []

    def convert(
        self,
        evidence_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        protein_groups_path: Optional[str] = None,
        chunksize: int = 500_000,
        creator: str = "maxquant",
        link_evidence: bool = False,
    ) -> None:
        """Run the evidence.txt -> feature.parquet conversion.

        For SILAC, the ``msms.txt`` beside *evidence_path* gives the labeling
        state of every MS/MS for ``identified_labels``.

        Args:
            evidence_path: Path to MaxQuant ``evidence.txt``.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF file for sample/channel mapping.
            protein_groups_path: Optional ``proteinGroups.txt`` for
                populating ``pg_global_qvalue`` and gene names on features.
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
            link_evidence: Record each evidence ``id``'s feature_id in
                :attr:`evidence_feature_ids`, for linking msms.txt PSMs.
        """
        self.evidence_feature_ids.clear()
        # Track which dropped-TMT-channel warnings have already been emitted so
        # the message is logged once, not once per batch.
        self._warned_tmt_channels.clear()

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

        # Step 3: Load the experiment type, isobaric channels, SILAC labels and fixed modifications
        channels = self._load_labels(sdrf_path, evidence_path, actual_cols)

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
                records = self._transform_batch(df, channels, pg_maps)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)
                    if link_evidence:
                        self._link_evidence(df, derive_feature_ids(records, _FEATURE_IDENTITY_COMPOSITE))

        self.logger.info(f"MaxQuant feature conversion complete -> {output_path}")

    def _load_labels(self, sdrf_path: Optional[str], evidence_path: str, columns: set[str]) -> list[tuple[str, int]]:
        """Load the SDRF run labels, SILAC labels and fixed modifications; return the isobaric channels."""
        self._run_labels, experiment_type, channels = self._load_sdrf(sdrf_path, columns)
        self._multiplet = multiplet_labels(ratio_letters(columns), self._run_labels, experiment_type)
        self._states_by_evidence = identified_states(self._conn, columns, evidence_path) if self._multiplet else {}
        self._fixed = self._fixed_mods(sdrf_path, evidence_path)
        if not channels and any(column.startswith("Reporter intensity") for column in columns):
            logger.warning(
                "evidence.txt has reporter intensity columns but no SDRF declares the channels; "
                "features carry the MS1 intensity only"
            )
        return channels

    def _link_evidence(self, df: pd.DataFrame, feature_ids: list[int]) -> None:
        """Map every evidence ``id`` behind the batch's features, collapsed duplicates included."""
        if _EVIDENCE_IDS in df.columns:
            id_lists = df[_EVIDENCE_IDS]
        elif "id" in df.columns:
            id_lists = [[value] for value in df["id"]]
        else:
            return
        for evidence_ids, feature_id in zip(id_lists, feature_ids):
            for evidence_id in evidence_ids:
                if evidence_id is not None and pd.notna(evidence_id):
                    self.evidence_feature_ids[int(evidence_id)] = feature_id

    def _identified_labels(self, row: dict) -> Optional[list[str]]:
        """Labels of the SILAC partners identified by MS/MS; ``None`` outside SILAC."""
        if not self._multiplet:
            return None
        evidence_ids = row[_EVIDENCE_IDS] if _EVIDENCE_IDS in row else [row.get("id")]
        return identified_labels(evidence_ids, row.get("Labeling State"), self._multiplet, self._states_by_evidence)

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
                sql_build("TRY_CAST($pep_column AS DOUBLE) ASC NULLS LAST", pep_column=validate_identifier(pep_column))
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

        # Every evidence id a kept row stands for, so msms ``Evidence ID`` can resolve to the feature.
        evidence_ids = (
            f", LIST(TRY_CAST({validate_identifier('id')} AS BIGINT)) OVER (PARTITION BY {partition_by}) AS {_EVIDENCE_IDS}"
            if "id" in actual_cols
            else ""
        )
        return sql_build(
            _EVIDENCE_DEDUP_SQL,
            evidence_ids=evidence_ids,
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
        count = self._conn.execute("SELECT COUNT(*) FROM evidence").fetchall()[0][0]
        self.logger.info(f"Loaded {count:,} MaxQuant evidence rows")

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _build_protein_group_maps(self, protein_groups_path: Optional[str]) -> dict:
        """Build protein-group lookup maps from proteinGroups.txt.

        Returns a dict with three sub-maps:
            ``by_id``: proteinGroups.txt ``id`` -> ``{"accessions", "qvalue",
                "genes"}``, resolving the evidence ``Protein group IDs``.
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
        empty = {"by_id": {}, "by_group": {}, "by_accession": {}}
        if not protein_groups_path:
            return empty
        df = self._conn.execute(
            "SELECT * FROM read_csv_auto($1, delim='\t', header=true, auto_detect=true, null_padding=true)",
            [protein_groups_path],
        ).df()
        acc_col, qval_col, gene_col = self._detect_pg_columns(df)
        if not acc_col:
            raise ValueError(f"{protein_groups_path} has no 'Protein IDs' or 'Majority protein IDs' column")
        maps = self._build_pg_lookup(df, acc_col, qval_col, gene_col)
        self.logger.info(
            "Built protein group maps: %d group id(s), %d group(s), %d unambiguous accession(s)",
            len(maps["by_id"]),
            len(maps["by_group"]),
            len(maps["by_accession"]),
        )
        return maps

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

    def _build_pg_lookup(self, df: pd.DataFrame, acc_col: str, qval_col: Optional[str], gene_col: Optional[str]) -> dict:
        """Build the ``by_id``/``by_group``/``by_accession`` lookup maps (see caller)."""
        by_id: dict[int, dict] = {}
        by_group: dict[frozenset, dict] = {}
        acc_to_groups: dict[str, set] = {}
        acc_value: dict[str, dict] = {}

        def values(column: Optional[str]) -> list:
            return df[column].tolist() if column and column in df.columns else [None] * len(df)

        acc_vals = df[acc_col].astype(str).tolist()
        qval_vals = pd.to_numeric(df[qval_col], errors="coerce").tolist() if qval_col else [None] * len(df)

        for acc_raw, qval, gene_raw, fasta_raw, group_id in zip(
            acc_vals, qval_vals, values(gene_col), values("Fasta headers"), values("id")
        ):
            accs = [canonical_accession(a.strip()) for a in acc_raw.split(";") if a.strip()]
            if not accs:
                continue
            value = {"qvalue": safe_float(qval), "genes": parse_genes(gene_raw, fasta_raw)}
            if pd.notna(group_id):
                by_id[int(group_id)] = {**value, "accessions": accs}
            key = self._canonical_group_key(accs)
            by_group.setdefault(key, value)
            for acc in accs:
                acc_to_groups.setdefault(acc, set()).add(key)
                acc_value.setdefault(acc, value)

        # Keep only accessions that map to a single protein group for the
        # leader fallback; ambiguous accessions would otherwise collide.
        by_accession = {acc: acc_value[acc] for acc, keys in acc_to_groups.items() if len(keys) == 1}
        return {"by_id": by_id, "by_group": by_group, "by_accession": by_accession}

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

    def _warn_dropped_tmt_channels(self, columns, channels: list[tuple[str, int]], ri_offset: int) -> None:
        """Warn when the SDRF plex declares channels absent from evidence.txt.

        A ``Reporter intensity N`` column missing for a declared channel means
        that channel is silently dropped from every feature; surface it once.
        """
        col_set = set(columns)
        dropped = [
            (label, f"Reporter intensity {index + ri_offset}")
            for label, index in channels
            if f"Reporter intensity {index + ri_offset}" not in col_set
        ]
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
        channels: list[tuple[str, int]],
        pg_maps: dict | None = None,
    ) -> list[dict]:
        """Transform a batch of evidence.txt rows into QPX feature records."""
        records: list[dict] = []
        # Detect reporter intensity column indexing once per batch
        ri_offset = self._detect_ri_offset(df.columns)
        # Warn (once) about SDRF channels with no reporter column in evidence.txt
        self._warn_dropped_tmt_channels(df.columns, channels, ri_offset)
        self._localization_columns = localization_columns(df.columns)
        # Pre-extract column arrays for faster per-row access than to_dict("records")
        col_arrays = {col: df[col].values for col in df.columns}
        n_rows = len(df)
        for i in range(n_rows):
            row = {col: vals[i] for col, vals in col_arrays.items()}
            rec = self._transform_row(row, channels, pg_maps, ri_offset)
            if rec:
                records.append(rec)
        return records

    def _transform_row(
        self,
        row,
        channels: list[tuple[str, int]],
        pg_maps: dict | None = None,
        ri_offset: int = 0,
    ) -> Optional[dict]:
        """Transform a single evidence.txt row."""
        r = self._resolved  # shorthand for resolved column mappings

        sequence = str(row.get(r.get("sequence", "Sequence"), ""))
        peptidoform = to_proforma(str(row.get(r.get("modified_sequence", "Modified sequence"), "")), self._fixed)
        _, modifications = from_proforma(peptidoform, sequence) if peptidoform else (None, None)
        sites = localization_sites(row, self._localization_columns, sequence)
        modifications, candidates = localize_modifications(modifications, sequence, sites)
        charge = int(row.get(r.get("charge", "Charge"), 0))
        run_file_name = str(row.get(r.get("run_file_name", "Raw file"), ""))

        # is_decoy (bool)
        is_decoy = mq_flag_to_bool(row.get(r.get("is_decoy", "Reverse"), ""))

        spectrum_fields = self._evidence_spectrum_fields(row, charge, peptidoform)

        pg_accessions, anchor_protein, unique, pg_qvalue, gg_names = self._evidence_protein_fields(row, pg_maps)

        match_fields = self._evidence_match_fields(row, run_file_name)

        # Build intensities (new schema: {label, intensity})
        intensities, additional_intensities = self._build_intensities(row, run_file_name, channels, ri_offset)

        additional_scores = self._evidence_scores(row)

        return {
            **match_fields,
            **spectrum_fields,
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "modification_candidates": candidates,
            "charge": charge,
            "is_decoy": is_decoy,
            "additional_scores": additional_scores or None,
            "predicted_rt": None,
            "run_file_name": run_file_name,
            "cv_params": None,
            "intensities": intensities or None,
            "additional_intensities": additional_intensities,
            "identified_labels": self._identified_labels(row),
            "pg_accessions": pg_accessions,
            "anchor_protein": anchor_protein,
            "unique": unique,
            "pg_global_qvalue": pg_qvalue,
            "ion_mobility_start": None,
            "ion_mobility_stop": None,
            "gg_accessions": gg_names,
            "gg_names": gg_names,
        }

    def write_pepmap(self, output_path: str) -> bool:
        """Write every protein ``Proteins`` lists for each target peptidoform of the converted evidence."""
        columns = self.get_table_columns("evidence")
        r = self._resolved
        sequence, modified = r.get("sequence", "Sequence"), r.get("modified_sequence", "Modified sequence")
        if "Proteins" not in columns:
            return False
        target = "" if "Reverse" not in columns else """WHERE COALESCE("Reverse", '') <> '+'"""
        rows = self._conn.execute(f"""SELECT DISTINCT "{sequence}", "{modified}", "Proteins" FROM evidence {target}""").fetchall()
        pepmap: dict[str, tuple[str, list[str]]] = {}
        for seq, mod, proteins in rows:
            peptidoform = to_proforma(str(mod), self._fixed)
            if proteins and peptidoform:
                pepmap.setdefault(peptidoform, (seq, [self._norm_uniprot_id(p) for p in proteins.split(";")]))
        return self._write_pepmap(output_path, pepmap, "maxquant")

    def _build_intensities(
        self,
        row,
        run_file_name: str,
        channels: list[tuple[str, int]],
        ri_offset: int = 0,
    ) -> tuple[list[dict], Optional[list[dict]]]:
        """Build intensities and additional_intensities from evidence row.

        Every evidence row is an identification, so a label without a quantity
        is kept with a null intensity; MaxQuant writes a missing quantity as 0.
        """
        if channels:
            intensities, additional_intensities = self._reporter_intensities(row, run_file_name, channels, ri_offset)
            return intensities, additional_intensities or None
        if self._multiplet:
            return multiplet_intensities(row, self._multiplet, self._run_labels.get(run_file_stem(run_file_name)))
        # LFQ: prefer MaxQuant-normalised "LFQ intensity" when available,
        # fall back to raw "Intensity" (e.g. for datasets without LFQ algo).
        lfq_val = safe_float(row.get("LFQ intensity"))
        raw_val = safe_float(row.get(_FEATURE_MAP["intensity"][0]))
        intensity_val = lfq_val if (lfq_val is not None and lfq_val > 0) else raw_val
        quantified = intensity_val is not None and intensity_val > 0
        count = safe_float(row.get("MS/MS count"))
        additional = None
        if count is not None:
            additional = [{"label": "LFQ", "intensities": [{"intensity_name": "spectral_count", "intensity_value": count}]}]
        return [{"label": "LFQ", "intensity": intensity_val if quantified else None}], additional

    def _reporter_intensities(self, row, run_file_name: str, channels: list[tuple[str, int]], ri_offset: int) -> tuple:
        """One intensity per channel the SDRF assigns in this run, read from its plex position."""
        intensities: list[dict] = []
        additional_intensities: list[dict] = []
        run_labels = self._run_labels.get(run_file_stem(run_file_name))
        for label, index in channels:
            col_name = f"Reporter intensity {index + ri_offset}"
            if col_name not in row or (run_labels is not None and label not in run_labels):
                continue
            val = safe_float(row.get(col_name))
            quantified = val is not None and val > 0
            intensities.append({"label": label, "intensity": val if quantified else None})
            corr_val = safe_float(row.get(col_name.replace("Reporter intensity", "Reporter intensity corrected")))
            # Beside a missing reporter intensity, a zero is MaxQuant's sentinel.
            if corr_val is not None and (quantified or corr_val != 0):
                additional_intensities.append(
                    {
                        "label": label,
                        "intensities": [{"intensity_name": "reporter_intensity_corrected", "intensity_value": corr_val}],
                    }
                )
        return intensities, additional_intensities

    def _evidence_spectrum_fields(self, row: dict, charge: int, peptidoform: str) -> dict:
        """Extract the spectrum, m/z, and retention-time fields."""
        r = self._resolved
        # Scan (list<int32>)
        scan_raw = row.get(r.get("scan", "MS/MS scan number"))
        scan = []
        if pd.notna(scan_raw):
            try:
                scan = [int(scan_raw)]
            except (ValueError, TypeError):
                pass

        calculated_mz = mq_calculated_mz(peptidoform, charge, safe_float(row.get(r.get("mass", "Mass"))))
        observed_mz = mq_observed_mz(
            safe_float(row.get(r.get("observed_mz", "m/z"))), safe_float(row.get(r.get("mass_error_ppm", "Mass error [ppm]")))
        )

        # RT (minutes → seconds)
        rt = safe_float(row.get(r.get("rt", "Retention time")))
        if rt is not None:
            rt *= 60.0
        rt_aligned = safe_float(row.get(r.get("rt_aligned", "Calibrated retention time")))
        if rt_aligned is not None:
            rt_aligned *= 60.0
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
        return {
            "scan": scan,
            "calculated_mz": calculated_mz,
            "observed_mz": observed_mz,
            "mass_error_ppm": mass_error_ppm(calculated_mz, observed_mz),
            "rt": rt,
            "rt_aligned": rt_aligned,
            "rt_start": rt_start,
            "rt_stop": rt_stop,
            "posterior_error_probability": pep,
            "ion_mobility": ion_mobility,
        }

    def _assigned_group(self, row: dict, pg_maps: dict | None, razor: str) -> tuple[Optional[dict], bool]:
        """The proteinGroups.txt entry MaxQuant assigns an evidence row, and whether its peptide is shared.

        A peptide in one group belongs to it; a peptide shared by several groups
        is assigned to the group of its leading razor protein. ``None`` when the
        groups are not in proteinGroups.txt or no group holds the razor protein.
        """
        by_id = (pg_maps or {}).get("by_id")
        ids = group_ids(row.get("Protein group IDs")) if by_id else []
        groups = [by_id[group_id] for group_id in ids if group_id in by_id]
        if len(ids) == 1:
            return (groups[0] if groups else None), False
        return next((group for group in groups if razor in group["accessions"]), None), len(ids) > 1

    def _evidence_protein_fields(self, row: dict, pg_maps: dict | None) -> tuple:
        """Resolve protein membership, group confidence, and genes."""
        r = self._resolved
        anchor_protein = self._norm_uniprot_id(str(row.get(r.get("anchor_protein", "Leading razor protein"), "")))
        if (pg_maps or {}).get("by_id") and "Protein group IDs" in row:
            pg_acc_list, anchor_protein, unique, pg_qvalue, pg_genes = self._assigned_group_fields(row, pg_maps, anchor_protein)
        else:
            pg_acc_list, anchor_protein, unique, pg_qvalue, pg_genes = self._leading_protein_fields(row, pg_maps, anchor_protein)
        pg_accessions = [{"accession": acc, "start": None, "end": None, "pre": None, "post": None} for acc in pg_acc_list] or None

        # Gene names
        gg_raw = row.get(r.get("gg_names", "Gene names"))
        gg_names = str(gg_raw).split(";") if pd.notna(gg_raw) and gg_raw else None
        if not gg_names:
            gg_names = pg_genes
        return pg_accessions, anchor_protein, unique, pg_qvalue, gg_names

    def _assigned_group_fields(self, row: dict, pg_maps: dict, razor: str) -> tuple:
        """``(members, anchor, unique, qvalue, genes)`` of the group MaxQuant assigns the evidence row."""
        group, shared = self._assigned_group(row, pg_maps, razor)
        if group is None:
            return [], None, None, None, None
        members = group["accessions"]
        anchor = razor if razor in members else members[0]
        return members, anchor, len(members) == 1 and not shared, group["qvalue"], group["genes"]

    def _leading_protein_fields(self, row: dict, pg_maps: dict | None, anchor: str) -> tuple:
        """``(members, anchor, unique, qvalue, genes)`` without proteinGroups.txt group ids.

        The group membership is unknown, so the evidence's leading proteins stand in for it.
        """
        pg_acc_raw = str(row.get(self._resolved.get("pg_accessions", "Leading proteins"), ""))
        members = [self._norm_uniprot_id(a) for a in pg_acc_raw.split(";")] if pg_acc_raw else []
        if not anchor and members:
            anchor = members[0]
        proteins = row.get("Proteins")
        unique = len(proteins.split(";")) == 1 if isinstance(proteins, str) else len(members) <= 1
        # Protein-group q-value and genes, keyed on the full membership.
        qvalue, genes = self._lookup_pg_metadata(pg_maps, members, anchor)
        return members, anchor, unique, qvalue, genes

    def _evidence_scores(self, row: dict) -> list[dict]:
        """Collect Andromeda score fields from an evidence row."""
        r = self._resolved
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
        return additional_scores

    def _evidence_match_fields(self, row: dict, run_file_name: str) -> dict:
        """Extract cleavage count and identification source."""
        r = self._resolved
        # Missed cleavages (dedicated field from evidence.txt)
        mc_raw = row.get(r.get("missed_cleavages", "Missed cleavages"))
        missed_cleavages = int(mc_raw) if pd.notna(mc_raw) else None

        # MBR detection: Type == "MULTI-MATCH" means transferred (MBR) feature
        evidence_type = str(row.get("Type", "")).strip().upper()
        id_run_file_name = None if evidence_type == "MULTI-MATCH" else run_file_name
        return {
            "missed_cleavages": missed_cleavages,
            "id_run_file_name": id_run_file_name,
        }
