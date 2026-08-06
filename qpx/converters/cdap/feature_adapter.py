"""CDAP Feature adapter -- aggregate PSMs into peptidoform+charge+run features.

CDAP only emits PSM-level files; quantitative features are derived by rolling
up PSMs that share ``(peptidoform, charge, run_file_name)``.  Reporter-ion
intensities are summed per channel across the contributing PSMs; representative
fields (m/z, RT, scan, scores) come from the best PSM (lowest spectral
``Evalue``).  Ambiguous top-ranked identifications remain available in
``psm.parquet`` but are excluded from feature quantification because one MS2
reporter signal cannot be assigned independently to multiple peptide matches.

For LFQ studies (no reporter channels) the ``PrecursorArea`` column from the
best PSM is emitted as a single intensity entry with label ``"LFQ"`` so
label-free consumers can infer the quantification category from QPX labels.
"""

from __future__ import annotations

import re
from typing import Optional

import pandas as pd

from qpx.converters.base import resolve_columns
from qpx.converters.cdap.base_adapter import CdapBaseAdapter
from qpx.converters.mappings import get_field_mappings
from qpx.converters.ptm import compute_precursor_mz
from qpx.writers.feature import FeatureWriter

_FEATURE_MAP = get_field_mappings("cdap", "feature")

# Producer-specific feature identity composite (bigbio/qpx#229). CDAP reports one
# feature per precursor per run, so the measured key is just peptidoform + charge +
# run (all in feature.yaml). Passed to FeatureWriter so feature_id hashes exactly
# these columns instead of the schema default.
_FEATURE_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name")

# Reporter-ion meta column suffixes (per channel prefix) skipped during the
# channel auto-detection.  Channels of interest end with the reporter-ion m/z
# token (e.g. ``-127N``, ``-131C``, or ``114``..``117`` for iTRAQ).
_META_SUFFIX_RE = re.compile(r"(Flags|FractionOfTotalAb|TotalAb)$")


class CdapFeatureAdapter(CdapBaseAdapter):
    """Convert CDAP ``*.psm`` files into ``feature.parquet`` via aggregation.

    Usage::

        with CdapFeatureAdapter() as adapter:
            adapter.convert(
                psm_dir="/data/CPTAC/PDC000440",
                output_path="cdap.feature.parquet",
            )
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._resolved: dict = {}
        self._experiment_type: str = "LFQ"
        self._channel_cols: list[str] = []

    def convert(
        self,
        psm_dir: str,
        output_path: str,
        chunksize: int = 200_000,
        creator: str = "cdap",
    ) -> None:
        """Run the feature aggregation for a single CPTAC study."""
        # Step 1: load all PSM rows for the study.
        self._load_psm_directory(psm_dir)

        # Step 2: detect reporter channels from the header.
        actual_cols = self.get_table_columns("psm")
        self._resolved = resolve_columns(_FEATURE_MAP, actual_cols)
        channel_cols, experiment_type = self._detect_channels(actual_cols)
        self.logger.info("Detected %s experiment with %d reporter channels", experiment_type, len(channel_cols))
        self._experiment_type = experiment_type
        self._channel_cols = channel_cols

        # Step 3: stream aggregated groups via DuckDB SQL.
        agg_sql = self._build_aggregation_sql(actual_cols, channel_cols)

        self.logger.info("Transforming CDAP features ...")
        with FeatureWriter(
            output_path,
            creator=creator,
            compression=self._compression,
            identity_composite=_FEATURE_IDENTITY_COMPOSITE,
        ) as writer:
            self._stream_transform_write(agg_sql, writer, chunksize)

        self.logger.info("CDAP feature conversion complete -> %s", output_path)

    # ------------------------------------------------------------------
    # SQL aggregation: rollup PSMs to (peptidoform, charge, run) groups
    # ------------------------------------------------------------------

    def _build_aggregation_sql(
        self,
        actual_cols: set[str],
        channel_cols: list[str],
    ) -> str:
        """Build a DuckDB SQL query that aggregates PSMs into features.

        The query relies on Python-side per-row helpers for protein parsing,
        modification re-parsing, and intensity ``intensity/ppm`` extraction,
        so it returns row arrays for the aggregated PSMs rather than a single
        canonical value where multiple PSMs disagree.
        """
        has_rt = "RTAtPrecursorHalfElution" in actual_cols
        has_precursor_area = "PrecursorArea" in actual_cols

        typed_cols: list[str] = [
            '"PeptideSequence" AS peptide_raw',
            'CAST("QueryCharge" AS INTEGER) AS charge',
            '"FileName" AS file_name',
            '"Protein" AS protein_cell',
            '"ScanNum" AS scan_num',
            'TRY_CAST("QueryPrecursorMz" AS DOUBLE) AS observed_mz',
            'TRY_CAST("OriginalPrecursorMz" AS DOUBLE) AS original_mz',
            'TRY_CAST("PrecursorError(ppm)" AS DOUBLE) AS ppm',
            'TRY_CAST("Evalue" AS DOUBLE) AS evalue',
            'TRY_CAST("MSGFScore" AS DOUBLE) AS msgf',
            'TRY_CAST("DeNovoScore" AS DOUBLE) AS denovo',
            'TRY_CAST("Qvalue" AS DOUBLE) AS qvalue',
            'TRY_CAST("PepQvalue" AS DOUBLE) AS pep_qvalue',
        ]
        typed_cols.append('TRY_CAST("RTAtPrecursorHalfElution" AS DOUBLE) AS rt' if has_rt else "CAST(NULL AS DOUBLE) AS rt")
        if has_precursor_area:
            typed_cols.append('TRY_CAST("PrecursorArea" AS DOUBLE) AS precursor_area')
        # Pre-extract reporter-ion intensities from "value/ppm" format in SQL
        # so Python never calls parse_intensity_cell.
        for col in channel_cols:
            qcol = self._quote_ident(col)
            typed_cols.append(
                f"TRY_CAST(split_part(CAST({qcol} AS VARCHAR), '/', 1) AS DOUBLE) AS {self._safe_channel_alias(col)}"
            )

        # Representative fields must all come from the SAME best PSM row --
        # independent ``arg_min`` per column produced a "Frankenstein" record
        # (scan/observed_mz/rt from different PSMs) on tied or NULL Evalue.
        # Applying one deterministic ORDER BY to every ``first()`` guarantees a
        # single source row. NULLS LAST so a real Evalue is preferred over a
        # missing one; scan/mz/rt break ties reproducibly.
        best_order = "evalue ASC NULLS LAST, scan_num ASC NULLS LAST, observed_mz ASC NULLS LAST, rt ASC NULLS LAST"
        agg_cols: list[str] = [
            "peptide_raw",
            "charge",
            "file_name",
            f"first(scan_num ORDER BY {best_order}) AS best_scan",
            f"first(observed_mz ORDER BY {best_order}) AS best_observed_mz",
            f"first(original_mz ORDER BY {best_order}) AS best_original_mz",
            f"first(ppm ORDER BY {best_order}) AS best_ppm",
            f"first(rt ORDER BY {best_order}) AS best_rt",
            f"first(evalue ORDER BY {best_order}) AS best_evalue",
            f"first(msgf ORDER BY {best_order}) AS best_msgf",
            f"first(denovo ORDER BY {best_order}) AS best_denovo",
            f"first(qvalue ORDER BY {best_order}) AS best_qvalue",
            f"first(pep_qvalue ORDER BY {best_order}) AS best_pep_qvalue",
            "list(protein_cell ORDER BY evalue) AS protein_cells",
            "count(*) AS psm_count",
        ]
        if has_precursor_area:
            agg_cols.append(f"first(precursor_area ORDER BY {best_order}) AS best_precursor_area")
        # SUM pre-extracted intensities directly in SQL.
        for col in channel_cols:
            alias = self._safe_channel_alias(col)
            agg_cols.append(f"SUM({alias}) AS {alias}")

        typed_block = ",\n                ".join(typed_cols)
        agg_block = ",\n            ".join(agg_cols)
        ambiguous_filter = ""
        if "AmbiguousMatch" in actual_cols:
            ambiguous_filter = (
                "AND COALESCE(lower(trim(CAST(\"AmbiguousMatch\" AS VARCHAR))), '') NOT IN ('1', 'true', 't', 'y', 'yes')\n"
            )
        return (
            f"WITH typed AS (\n"
            f"    SELECT\n                {typed_block}\n"
            f"    FROM psm\n"
            f'    WHERE "PeptideSequence" IS NOT NULL '
            f'AND "QueryCharge" IS NOT NULL\n'
            f"    {ambiguous_filter}"
            f")\n"
            f"SELECT\n            {agg_block}\n"
            f"FROM typed\n"
            f"GROUP BY peptide_raw, charge, file_name"
        )

    @staticmethod
    def _quote_ident(name: str) -> str:
        """Quote a DuckDB identifier, escaping embedded double quotes."""
        escaped = name.replace('"', '""')
        return f'"{escaped}"'

    @staticmethod
    def _safe_channel_alias(col: str) -> str:
        """Map a reporter-ion column name to a SQL-safe alias."""
        # ``TMT11-127N`` -> ``ch_TMT11_127N``; ``iTRAQ114`` -> ``ch_iTRAQ114``.
        return "ch_" + re.sub(r"[^A-Za-z0-9]+", "_", col)

    # ------------------------------------------------------------------
    # Transformation
    # ------------------------------------------------------------------

    def _transform_batch(self, df: pd.DataFrame) -> list[dict]:
        """Build feature records from one aggregation batch."""
        channel_cols = self._channel_cols
        channel_aliases = [self._safe_channel_alias(c) for c in channel_cols]
        experiment_type = self._experiment_type

        def _row_fn(row: dict):
            return self._transform_row(row, channel_cols, channel_aliases, experiment_type)

        return self._transform_batch_with_skip(df, _row_fn, "feature")

    def _transform_row(
        self,
        row: dict,
        channel_cols: list[str],
        channel_aliases: list[str],
        experiment_type: str,
    ) -> Optional[dict]:
        """Materialise one aggregated PSM group into a QPX Feature record.

        Numeric fields pre-cast to DOUBLE in SQL may be NaN for unparseable
        values; ``_finite_float`` / ``_optional_float`` normalise them safely.
        """
        identity = self._parse_peptide_identity(row)
        if identity is None:
            return None
        charge = self._int_or_zero(row.get("charge"))
        peptidoform = identity[1]
        run_file_name = self.strip_run_extension(row.get("file_name"))

        scan_val = self._int_or_zero(row.get("best_scan"))
        scan: list[int] = [scan_val] if scan_val else []

        observed_mz = self._finite_float(row.get("best_observed_mz"))
        # ``calculated_mz`` is a THEORETICAL mass. When the peptidoform carries
        # an unparseable token (e.g. a CHEMMOD unknown-mass tag), keep it NULL
        # rather than falling back to a measured precursor m/z, which would
        # silently corrupt the theoretical-mass semantics (bigbio/qpx#250).
        calculated_mz = compute_precursor_mz(peptidoform, charge) if peptidoform and charge else None

        mass_error_ppm = self._optional_float(row.get("best_ppm"))
        rt = self._optional_float(row.get("best_rt"))

        # Protein + decoy resolution (delegated to sub-method).
        pg_accessions, anchor_protein, unique, is_decoy = self._resolve_proteins(row)

        # Intensities (delegated to sub-method).
        intensities = self._collect_intensities(row, channel_cols, channel_aliases)

        # Scores (delegated to sub-method).
        additional_scores = self._collect_feature_scores(row)

        psm_count = self._int_or_zero(row.get("psm_count"))
        cv_params = [
            {"cv_name": "experiment_type", "cv_value": experiment_type},
            {"cv_name": "psm_count", "cv_value": str(psm_count)},
        ]

        mz_values = (calculated_mz, observed_mz, mass_error_ppm, additional_scores)
        rec = self._common_record_fields(identity, charge, is_decoy, mz_values, run_file_name)
        rec.update(
            {
                "cv_params": cv_params,
                "scan": scan,
                "rt": rt,
                "ion_mobility": None,
                "missed_cleavages": None,
                "intensities": intensities or None,
                "additional_intensities": None,
                "pg_accessions": pg_accessions,
                "anchor_protein": anchor_protein,
                "unique": unique,
                "pg_global_qvalue": None,
                "ion_mobility_start": None,
                "ion_mobility_stop": None,
                "gg_accessions": None,
                "gg_names": None,
                "id_run_file_name": run_file_name,
                "rt_start": None,
                "rt_stop": None,
            }
        )
        return rec

    # ------------------------------------------------------------------
    # Sub-methods extracted from _transform_row to reduce CCN
    # ------------------------------------------------------------------

    def _resolve_proteins(self, row: dict) -> tuple[list[dict] | None, str, bool, bool]:
        """Parse protein cells, deduplicate, and resolve decoy status.

        Returns (pg_accessions, anchor_protein, unique, is_decoy).
        """
        protein_records: list[dict] = []
        seen_accessions: set[str] = set()
        decoy_flags: list[bool] = []
        protein_cells = row.get("protein_cells")
        if protein_cells is not None:
            for cell in protein_cells:
                cell_records, cell_is_decoy = self.parse_protein_cell(cell)
                decoy_flags.append(cell_is_decoy)
                for rec in cell_records:
                    if rec["accession"] not in seen_accessions:
                        seen_accessions.add(rec["accession"])
                        protein_records.append(rec)
        is_decoy = bool(decoy_flags) and all(decoy_flags)

        pg_accessions = [
            {
                "accession": rec["accession"],
                "start": None,
                "end": None,
                "pre": rec.get("pre"),
                "post": rec.get("post"),
            }
            for rec in protein_records
        ] or None
        anchor_protein = protein_records[0]["accession"] if protein_records else ""
        unique = len(protein_records) == 1
        return pg_accessions, anchor_protein, unique, is_decoy

    def _collect_intensities(
        self,
        row: dict,
        channel_cols: list[str],
        channel_aliases: list[str],
    ) -> list[dict]:
        """Extract reporter-ion or LFQ intensities from the row."""
        intensities: list[dict] = []
        if channel_cols:
            for col, alias in zip(channel_cols, channel_aliases):
                val = self._optional_float(row.get(alias))
                if val is not None and val > 0:
                    intensities.append({"label": col, "intensity": val})
        elif "best_precursor_area" in row:
            val = self._optional_float(row.get("best_precursor_area"))
            if val is not None and val > 0:
                intensities.append({"label": "LFQ", "intensity": val})
        return intensities

    _FEATURE_SCORE_DEFS: tuple[tuple[str, str, bool], ...] = (
        ("best_evalue", "msgf_specevalue", False),
        ("best_msgf", "msgf_rawscore", True),
        ("best_denovo", "msgf_denovo_score", True),
        ("best_qvalue", "psm_level_qvalue", False),
        ("best_pep_qvalue", "peptide_level_qvalue", False),
    )

    def _collect_feature_scores(self, row: dict) -> list[dict]:
        """Extract best-PSM scores for the feature record."""
        return self._extract_scores_from_columns(row, self._FEATURE_SCORE_DEFS)
