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
from qpx.converters.maxquant.channels import reporter_column_count
from qpx.converters.maxquant.constants import (
    PROTON_MASS,
    to_proforma,
)
from qpx.converters.ptm import from_proforma
from qpx.converters.utils import mq_flag_to_bool, safe_float
from qpx.core.sql import sql_build, validate_identifier
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

# Derive field map from central mappings
_FEATURE_MAP = get_field_mappings("maxquant", "feature")

_EVIDENCE_DEDUP_SQL = """
            SELECT *
            FROM evidence
            QUALIFY ROW_NUMBER() OVER (
                PARTITION BY $partition_by
                ORDER BY $order_by
            ) = 1
            """

# Producer-specific feature identity composite (bigbio/qpx#229). MaxQuant reports
# an RT window per feature, so its measured key uses rt_start/rt_stop (both in
# feature.yaml) rather than a single rt. Passed to FeatureWriter so feature_id
# hashes exactly these columns instead of the schema default.
_FEATURE_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name", "rt_start", "rt_stop", "scan")


def _reporter_columns(tmt_channels: dict[str, int], ri_offset: int) -> list[tuple[str, str]]:
    """``(label, evidence.txt column)`` per channel; ``ri_offset`` is the number of the first column."""
    return [(name, f"Reporter intensity {col_idx + ri_offset}") for name, col_idx in tmt_channels.items()]


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
        # proteinGroups.txt columns pg_accessions/anchor_protein are read from
        # when groups resolve through ``Protein group IDs`` (see get_resolved_columns).
        self._group_columns: dict[str, str] = {}

    def get_resolved_columns(self) -> dict[str, str]:
        """QPX field -> source column, as recorded in ontology.parquet.

        When groups resolve through ``Protein group IDs``, ``pg_accessions`` and
        ``anchor_protein`` are read from the proteinGroups.txt columns the pg
        view also maps (``Protein IDs``, ``Majority protein IDs``), not from
        ``Leading proteins`` / ``Leading razor protein``.
        """
        return {**super().get_resolved_columns(), **self._group_columns}

    def convert(
        self,
        evidence_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        protein_groups_path: Optional[str] = None,
        chunksize: int = 500_000,
        creator: str = "maxquant",
    ) -> None:
        """Run the evidence.txt -> feature.parquet conversion.

        Args:
            evidence_path: Path to MaxQuant ``evidence.txt``.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF file for sample/channel mapping.
            protein_groups_path: Optional ``proteinGroups.txt``. With it, each
                feature takes the group its ``Protein group IDs`` name (the
                razor group for a shared peptide): the group's full membership
                in ``pg_accessions``, its representative in ``anchor_protein``,
                ``unique``, ``pg_global_qvalue`` and gene names. Without it,
                ``pg_accessions`` falls back to ``Leading proteins``.
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
        """
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

        # Step 3: Load the experiment type and TMT channels from the SDRF
        _, experiment_type, tmt_channels = self._load_sdrf(sdrf_path, reporter_column_count(actual_cols))
        self._fixed_mods = self._load_fixed_modifications(sdrf_path)

        # Step 3b: Build lookup maps from proteinGroups.txt
        pg_maps = self._build_protein_group_maps(protein_groups_path)
        self._group_columns = self._group_source_columns(pg_maps, actual_cols)

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
                records = self._transform_batch(df, experiment_type, tmt_channels, pg_maps)
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
            _EVIDENCE_DEDUP_SQL,
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

        Every map value is ``{"accessions", "anchor", "qvalue", "genes"}``,
        where ``accessions`` and ``anchor`` are exactly what the pg view writes
        for the group. Returns a dict with three sub-maps, and ``columns``,
        the proteinGroups.txt columns ``accessions`` and ``anchor`` come from:
            ``by_id``: proteinGroups ``id`` -> value. evidence.txt names each
                row's groups in ``Protein group IDs``, so this is the primary,
                exact key.
            ``by_group``: full protein-group membership (a ``frozenset`` of
                normalised accessions) -> value. Used, like ``by_accession``,
                only for rows without usable ``Protein group IDs``.
            ``by_accession``: single accession -> value, but *only* for
                accessions that belong to exactly one protein group.
                Ambiguous razor accessions (shared across groups) are excluded
                so a feature cannot inherit another group's q-value/genes via
                its leader.
        """
        empty = {"by_id": {}, "by_group": {}, "by_accession": {}, "columns": {}}
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
        anchor_col = "Majority protein IDs" if "Majority protein IDs" in df.columns else acc_col
        maps["columns"] = {"pg_accessions": acc_col, "anchor_protein": anchor_col}
        self.logger.info(
            "Built protein group maps: %d group(s), %d group id(s), %d unambiguous accession(s)",
            len(maps["by_group"]),
            len(maps["by_id"]),
            len(maps["by_accession"]),
        )
        return maps

    @staticmethod
    def _group_source_columns(pg_maps: dict, evidence_columns: set[str]) -> dict[str, str]:
        """The proteinGroups.txt columns of the group fields, when ``Protein group IDs`` can resolve groups."""
        return pg_maps.get("columns", {}) if pg_maps.get("by_id") and "Protein group IDs" in evidence_columns else {}

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
            genes = MaxQuantFeatureAdapter._genes_from_headers(fasta_raw)
        return genes

    @staticmethod
    def _column_values(df: pd.DataFrame, column: Optional[str]) -> list:
        """The values of ``column``, or ``None`` for every row when it is absent."""
        return df[column].tolist() if column and column in df.columns else [None] * len(df)

    def _build_pg_lookup(self, df: pd.DataFrame, acc_col: str, qval_col: Optional[str], gene_col: Optional[str]) -> dict:
        """Build the ``by_id``/``by_group``/``by_accession`` lookup maps (see caller)."""
        by_id: dict[str, dict] = {}
        by_group: dict[frozenset, dict] = {}
        acc_to_groups: dict[str, set] = {}
        acc_value: dict[str, dict] = {}

        qval_vals = pd.to_numeric(df[qval_col], errors="coerce").tolist() if qval_col else [None] * len(df)
        columns = zip(
            self._column_values(df, acc_col),
            qval_vals,
            self._column_values(df, gene_col),
            self._column_values(df, "Fasta headers"),
            self._column_values(df, "Majority protein IDs"),
            self._column_values(df, "id"),
        )
        for acc_raw, qval, gene_raw, fasta_raw, majority_raw, group_id in columns:
            value = self._group_value(acc_raw, qval, self._parse_genes(gene_raw, fasta_raw), majority_raw)
            if value is None:
                continue
            group_key = self._group_ids(group_id)
            if len(group_key) == 1:
                by_id[group_key[0]] = value
            key = self._canonical_group_key(value["accessions"])
            by_group.setdefault(key, value)
            for acc in value["accessions"]:
                acc_to_groups.setdefault(acc, set()).add(key)
                acc_value.setdefault(acc, value)

        # Keep only accessions that map to a single protein group for the
        # leader fallback; ambiguous accessions would otherwise collide.
        by_accession = {acc: acc_value[acc] for acc, keys in acc_to_groups.items() if len(keys) == 1}
        return {"by_id": by_id, "by_group": by_group, "by_accession": by_accession}

    def _group_value(self, acc_raw, qval, genes: Optional[list[str]], majority_raw) -> Optional[dict]:
        """One proteinGroups row as a lookup value, or ``None`` when it has no accession."""
        accs = self._protein_list(acc_raw)
        if not accs:
            return None
        qvalue = float(qval) if qval is not None and pd.notna(qval) else None
        return {"accessions": accs, "anchor": self._group_anchor(majority_raw, accs), "qvalue": qvalue, "genes": genes}

    @staticmethod
    def _group_ids(cell) -> list[str]:
        """The proteinGroups ids in a MaxQuant id cell (``"12"``, ``12`` or ``"12;40"``)."""
        if cell is None or pd.isna(cell):
            return []
        if isinstance(cell, float):
            cell = int(cell)
        return [part.strip() for part in str(cell).split(";") if part.strip()]

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

    def _group_from_ids(self, row: dict, pg_maps: Optional[dict]) -> Optional[tuple]:
        """Resolve the row's protein group from ``Protein group IDs``, MaxQuant's own key.

        Returns ``(group, unique)``, where ``group`` is a ``by_id`` value, or
        ``None`` when the ids do not resolve (no proteinGroups.txt, no column, or
        ids it does not contain). A peptide in one group takes that group, and is
        unique when the group has a single protein. A peptide shared by several
        groups takes the group MaxQuant assigned it to, the one whose leading
        protein is the ``Leading razor protein``, and is not unique; when no
        single group matches, the group is unresolved (``(None, None)``).
        """
        by_id = pg_maps.get("by_id") if pg_maps else None
        if not by_id:
            return None
        ids = self._group_ids(row.get("Protein group IDs"))
        groups = [by_id[group_id] for group_id in ids if group_id in by_id]
        if not groups:
            return None
        if len(ids) == 1:
            return groups[0], len(groups[0]["accessions"]) == 1
        razor_group = self._razor_group(row, groups)
        return (razor_group, False) if razor_group else (None, None)

    def _razor_group(self, row: dict, groups: list[dict]) -> Optional[dict]:
        """The one group led by the row's ``Leading razor protein``, else the one containing it."""
        razor = next(iter(self._protein_list(row.get(self._resolved.get("anchor_protein", "Leading razor protein")))), "")
        matches = [group for group in groups if group["anchor"] == razor]
        if len(matches) != 1:
            matches = [group for group in groups if razor in group["accessions"]]
        return matches[0] if len(matches) == 1 else None

    def _group_from_leading_proteins(self, row: dict, pg_maps: Optional[dict]) -> tuple:
        """Fallback group for rows ``Protein group IDs`` cannot resolve: ``Leading proteins``.

        ``Leading proteins`` holds one leader per group, so the group's size is
        unknown here, and ``unique`` means the peptide is in a single group of
        any size.
        """
        r = self._resolved
        pg_acc_raw = str(row.get(r.get("pg_accessions", "Leading proteins"), ""))
        pg_acc_list = [self._norm_uniprot_id(a) for a in pg_acc_raw.split(";")] if pg_acc_raw else []
        anchor_protein = self._norm_uniprot_id(str(row.get(r.get("anchor_protein", "Leading razor protein"), "")))
        if not anchor_protein and pg_acc_list:
            anchor_protein = pg_acc_list[0]
        pg_qvalue, pg_genes = self._lookup_pg_metadata(pg_maps, pg_acc_list, anchor_protein)
        group = {"accessions": pg_acc_list, "anchor": anchor_protein, "qvalue": pg_qvalue, "genes": pg_genes}
        return group, len(pg_acc_list) <= 1

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

    def _warn_dropped_tmt_channels(self, columns, experiment_type: str, tmt_channels: dict[str, int], ri_offset: int) -> None:
        """Warn when the SDRF plex declares channels absent from evidence.txt.

        ``tmt_channels`` maps each label to its 0-based reporter column. A
        ``Reporter intensity N`` column missing for a declared channel means
        that channel is silently dropped from every feature; surface it once.
        """
        if re.sub(r"\d+", "", experiment_type).upper() not in ("TMT", "ITRAQ"):
            return
        col_set = set(columns)
        dropped = [(name, column) for name, column in _reporter_columns(tmt_channels, ri_offset) if column not in col_set]
        key = tuple(sorted(channel_name for channel_name, _ in dropped))
        if not dropped or key in self._warned_tmt_channels:
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
        experiment_type: str,
        tmt_channels: dict[str, int],
        pg_maps: dict | None = None,
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
                experiment_type,
                tmt_channels,
                pg_maps,
                ri_offset,
            )
            if rec:
                records.append(rec)
        return records

    def _transform_row(
        self,
        row,
        experiment_type: str,
        tmt_channels: dict[str, int],
        pg_maps: dict | None = None,
        ri_offset: int = 0,
    ) -> Optional[dict]:
        """Transform a single evidence.txt row."""
        r = self._resolved  # shorthand for resolved column mappings

        sequence = str(row.get(r.get("sequence", "Sequence"), ""))
        peptidoform = to_proforma(
            str(row.get(r.get("modified_sequence", "Modified sequence"), "")),
            self._fixed_mods,
        )
        _, modifications = from_proforma(peptidoform, sequence) if peptidoform else (None, None)
        charge = int(row.get(r.get("charge", "Charge"), 0))
        run_file_name = str(row.get(r.get("run_file_name", "Raw file"), ""))

        # is_decoy (bool)
        is_decoy = mq_flag_to_bool(row.get(r.get("is_decoy", "Reverse"), ""))

        spectrum_fields = self._evidence_spectrum_fields(row, charge)

        pg_accessions, anchor_protein, unique, pg_qvalue, gg_names = self._evidence_protein_fields(row, pg_maps)

        match_fields = self._evidence_match_fields(row, run_file_name)

        # Build intensities (new schema: {label, intensity})
        intensities, additional_intensities = self._build_intensities(row, experiment_type, tmt_channels, ri_offset)

        additional_scores = self._evidence_scores(row)

        return {
            **match_fields,
            **spectrum_fields,
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "charge": charge,
            "is_decoy": is_decoy,
            "additional_scores": additional_scores or None,
            "predicted_rt": None,
            "run_file_name": run_file_name,
            "cv_params": None,
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
        }

    def _build_intensities(
        self,
        row,
        experiment_type: str,
        tmt_channels: dict[str, int],
        ri_offset: int = 0,
    ) -> tuple[list[dict], Optional[list[dict]]]:
        """Build intensities and additional_intensities from evidence row."""
        if re.sub(r"\d+", "", experiment_type).upper() in ("TMT", "ITRAQ") and tmt_channels:
            return self._reporter_intensities(row, tmt_channels, ri_offset)
        return self._lfq_intensities(row), None

    @staticmethod
    def _reporter_intensities(row, tmt_channels: dict[str, int], ri_offset: int) -> tuple[list[dict], Optional[list[dict]]]:
        """One intensity per channel with a positive reporter intensity, read from the channel's column.

        ``tmt_channels`` maps each label to its 0-based MaxQuant reporter column
        (its position in the plex), so N/C isomers and the plexes that share a
        label at different positions read their own column.
        """
        intensities: list[dict] = []
        additional_intensities: list[dict] = []
        for channel_name, col_name in _reporter_columns(tmt_channels, ri_offset):
            val = safe_float(row.get(col_name))
            if not val or val <= 0:
                continue
            intensities.append({"label": channel_name, "intensity": val})
            corr_val = safe_float(row.get(col_name.replace("Reporter intensity", "Reporter intensity corrected")))
            if corr_val is not None:
                additional_intensities.append(
                    {
                        "label": channel_name,
                        "intensities": [{"intensity_name": "corrected_reporter_intensity", "intensity_value": corr_val}],
                    }
                )
        return intensities, additional_intensities or None

    @staticmethod
    def _lfq_intensities(row) -> list[dict]:
        """The label-free intensity of an evidence row, or none without a quantity.

        Prefers MaxQuant-normalised "LFQ intensity" when available and falls back
        to raw "Intensity" (e.g. for datasets without LFQ algo). MaxQuant writes 0
        (or nothing) for evidence without an MS1 quantity, such as MS/MS-only
        rows. The intensity struct cannot hold a null, so such a feature carries
        no intensity entry, like an empty TMT channel.
        """
        lfq_val = safe_float(row.get("LFQ intensity"))
        raw_val = safe_float(row.get(_FEATURE_MAP["intensity"][0]))
        intensity_val = lfq_val if (lfq_val is not None and lfq_val > 0) else raw_val
        return [{"label": "LFQ", "intensity": intensity_val}] if intensity_val else []

    @staticmethod
    def _genes_from_headers(fasta_raw: str) -> Optional[list[str]]:
        """Extract distinct gene symbols from FASTA headers."""
        genes = []
        for hdr in str(fasta_raw).split(";"):
            m = re.search(r"GN=([^\s;]+)", hdr)
            if m and m.group(1) not in genes:
                genes.append(m.group(1))
        genes = genes or None
        return genes

    def _evidence_spectrum_fields(self, row: dict, charge: int) -> dict:
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

        # Calculated m/z from neutral mass and charge
        neutral_mass = safe_float(row.get(r.get("mass", "Mass")))
        if neutral_mass and charge:
            calculated_mz = (neutral_mass + charge * PROTON_MASS) / charge
        else:
            calculated_mz = None

        # m/z
        observed_mz = safe_float(row.get(r.get("observed_mz", "m/z")))

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
        return {
            "scan": scan,
            "calculated_mz": calculated_mz,
            "observed_mz": observed_mz,
            "rt": rt,
            "rt_start": rt_start,
            "rt_stop": rt_stop,
            "posterior_error_probability": pep,
            "ion_mobility": ion_mobility,
        }

    def _evidence_protein_fields(self, row: dict, pg_maps: dict | None) -> tuple:
        """Resolve protein membership, group confidence, and genes."""
        resolved = self._group_from_ids(row, pg_maps)
        group, unique = resolved if resolved is not None else self._group_from_leading_proteins(row, pg_maps)
        group = group or {"accessions": [], "anchor": None, "qvalue": None, "genes": None}
        # Protein accessions (list of pg_protein structs)
        pg_accessions = [
            {"accession": acc, "start": None, "end": None, "pre": None, "post": None} for acc in group["accessions"]
        ] or None

        # Gene names
        gg_raw = row.get(self._resolved.get("gg_names", "Gene names"))
        gg_names = str(gg_raw).split(";") if pd.notna(gg_raw) and gg_raw else None
        if not gg_names:
            gg_names = group["genes"]
        return pg_accessions, group["anchor"], unique, group["qvalue"], gg_names

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
        """Extract mass error, cleavage count, and identification source."""
        r = self._resolved
        # Mass error (ppm) — direct column from evidence.txt
        mass_error_ppm = safe_float(row.get(r.get("mass_error_ppm", "Mass error [ppm]")))

        # Missed cleavages (dedicated field from evidence.txt)
        mc_raw = row.get(r.get("missed_cleavages", "Missed cleavages"))
        missed_cleavages = int(mc_raw) if pd.notna(mc_raw) else None

        # MBR detection: Type == "MULTI-MATCH" means transferred (MBR) feature
        evidence_type = str(row.get("Type", "")).strip().upper()
        id_run_file_name = None if evidence_type == "MULTI-MATCH" else run_file_name
        return {
            "mass_error_ppm": mass_error_ppm,
            "missed_cleavages": missed_cleavages,
            "id_run_file_name": id_run_file_name,
        }
