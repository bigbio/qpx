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
from qpx.converters.maxquant.evidence_fields import identification_fields
from qpx.converters.maxquant.protein_fields import group_ids
from qpx.converters.maxquant.silac import (
    identified_labels,
    identified_states,
    multiplet_intensities,
    multiplet_labels,
    ratio_letters,
)
from qpx.converters.ptm import from_proforma
from qpx.converters.utils import mass_error_ppm, mq_flag_to_bool, safe_float
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
    ):
        super().__init__(
            duckdb_memory=duckdb_memory,
            duckdb_threads=duckdb_threads,
            conn=conn,
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
        link_evidence: bool = False,
    ) -> None:
        """Run the evidence.txt -> feature.parquet conversion.

        For SILAC, the ``msms.txt`` beside *evidence_path* gives the labeling
        state of every MS/MS for ``identified_labels``.

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
        self._group_columns = self._group_source_columns(pg_maps, actual_cols)

        # Step 4: Stream and transform
        self.logger.info("Transforming MaxQuant features ...")
        evidence_query = self._deduplicated_evidence_query(actual_cols)

        with FeatureWriter(
            output_path,
            creator=creator,
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

    @staticmethod
    def _group_source_columns(pg_maps: dict, evidence_columns: set[str]) -> dict[str, str]:
        """The proteinGroups.txt columns of the group fields, when ``Protein group IDs`` can resolve groups."""
        return pg_maps.get("columns", {}) if pg_maps.get("by_id") and "Protein group IDs" in evidence_columns else {}

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
        ids = group_ids(row.get("Protein group IDs"))
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
        unknown here; ``unique`` is whether the peptide maps to a single protein
        (``Proteins``), else whether it has a single leader.
        """
        r = self._resolved
        pg_acc_raw = str(row.get(r.get("pg_accessions", "Leading proteins"), ""))
        pg_acc_list = [self._norm_uniprot_id(a) for a in pg_acc_raw.split(";")] if pg_acc_raw else []
        anchor_protein = self._norm_uniprot_id(str(row.get(r.get("anchor_protein", "Leading razor protein"), "")))
        if not anchor_protein and pg_acc_list:
            anchor_protein = pg_acc_list[0]
        pg_qvalue, pg_genes = self._lookup_pg_metadata(pg_maps, pg_acc_list, anchor_protein)
        group = {"accessions": pg_acc_list, "anchor": anchor_protein, "qvalue": pg_qvalue, "genes": pg_genes}
        proteins = row.get("Proteins")
        return group, len(proteins.split(";")) == 1 if isinstance(proteins, str) else len(pg_acc_list) <= 1

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
        rows = self._conn.execute(
            sql_build(
                'SELECT DISTINCT $sequence, $modified, "Proteins" FROM evidence $target',
                sequence=validate_identifier(sequence),
                modified=validate_identifier(modified),
                target=target,
            )
        ).fetchall()
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
            "ccs": safe_float(row.get("CCS")),
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
            **identification_fields(row, evidence_type),
        }
