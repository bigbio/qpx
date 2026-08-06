"""CDAP PSM adapter -- *.psm tab-separated files to psm.parquet.

Loads every ``*.psm`` file in a CPTAC study directory into DuckDB and emits
one ``PsmSchema`` record per row.  Reporter-ion intensities are not stored on
PSM rows (PSM schema has no ``intensities`` field); they roll up into the
``feature.parquet`` and ``pg.parquet`` outputs produced by sibling adapters.
"""

from __future__ import annotations

from typing import Optional

import pandas as pd

from qpx.converters.base import resolve_columns
from qpx.converters.cdap.base_adapter import CdapBaseAdapter
from qpx.converters.mappings import get_field_mappings
from qpx.converters.ptm import compute_precursor_mz
from qpx.writers.psm import PsmWriter

_PSM_MAP = get_field_mappings("cdap", "psm")

# CV-param columns surfaced as ``cv_params`` on PSM records.  Keys map to the
# QPX-friendly ``cv_name`` used in the output struct.
_CV_PARAM_COLUMNS: dict[str, str] = {
    "AmbiguousMatch": "ambiguous_match",
    "HCDEnergy": "hcd_collision_energy",
    "PrecursorPurity": "precursor_purity",
    "FractionDecomposition": "fraction_decomposition",
    "PhosphoRSPeptide": "phosphors_peptide",
    "FullyLocalized": "fully_localized",
}

# Score columns surfaced as ``additional_scores``.  ``higher_better`` follows
# the MS-GF+ documentation -- raw scores ascend, statistical metrics descend.
_SCORE_COLUMNS: tuple[tuple[str, str, bool], ...] = (
    ("DeNovoScore", "msgf_denovo_score", True),
    ("MSGFScore", "msgf_rawscore", True),
    ("Evalue", "msgf_specevalue", False),
    ("Qvalue", "psm_level_qvalue", False),
    ("PepQvalue", "peptide_level_qvalue", False),
)


class CdapPsmAdapter(CdapBaseAdapter):
    """Convert CDAP ``*.psm`` files into ``psm.parquet``.

    Usage::

        with CdapPsmAdapter() as adapter:
            adapter.convert(
                psm_dir="/data/CPTAC/PDC000440",
                output_path="cdap.psm.parquet",
            )
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._resolved: dict = {}

    def convert(
        self,
        psm_dir: str,
        output_path: str,
        chunksize: int = 500_000,
        creator: str = "cdap",
    ) -> None:
        """Run the PSM conversion for a single CPTAC study directory.

        Args:
            psm_dir: Directory containing ``*.psm`` files for one study.
            output_path: Destination ``psm.parquet`` path.
            chunksize: Rows per DuckDB→pandas batch.
            creator: Creator tag stored in the Parquet metadata.
        """
        # Step 1: load every .psm into one DuckDB table.
        self._load_psm_directory(psm_dir)

        # Step 2: resolve column mappings against the actual header.
        actual_cols = self.get_table_columns("psm")
        self._resolved = resolve_columns(_PSM_MAP, actual_cols)

        # Step 3: build a SQL projection that pre-casts numerics and computes
        # simple derived columns so the Python row loop does minimal work.
        sql = self._build_typed_sql(actual_cols)

        # Step 4: stream batches and transform.
        self.logger.info("Transforming CDAP PSMs ...")
        with PsmWriter(output_path, creator=creator, compression=self._compression) as writer:
            self._stream_transform_write(sql, writer, chunksize)

        self.logger.info("CDAP PSM conversion complete -> %s", output_path)

    def _build_typed_sql(self, actual_cols: set[str]) -> str:
        """Build a SQL SELECT that pre-casts all numeric/string columns."""
        r = self._resolved
        has_rt = "RTAtPrecursorHalfElution" in actual_cols

        cols: list[str] = [
            f'"{r.get("peptide_raw", "PeptideSequence")}" AS peptide_raw',
            f'CAST("{r.get("charge", "QueryCharge")}" AS INTEGER) AS charge',
            # Pre-clean run file name: strip .raw/.mzML/.mgf extensions
            f"regexp_replace(\"{r.get('run_file_name', 'FileName')}\", '\\.(raw|RAW|mzML|mzml|mgf)$', '') AS run_file_name",
            # Scan as integer
            f'TRY_CAST("{r.get("scan", "ScanNum")}" AS INTEGER) AS scan_num',
            # All numeric fields pre-cast to DOUBLE
            f'TRY_CAST("{r.get("observed_mz", "QueryPrecursorMz")}" AS DOUBLE) AS observed_mz',
            'TRY_CAST("OriginalPrecursorMz" AS DOUBLE) AS original_mz',
            f'TRY_CAST("{r.get("mass_error_ppm", "PrecursorError(ppm)")}" AS DOUBLE) AS mass_error_ppm',
            'TRY_CAST("RTAtPrecursorHalfElution" AS DOUBLE) AS rt' if has_rt else "CAST(NULL AS DOUBLE) AS rt",
            # Protein cell kept as string for Python-side regex parsing
            f'"{r.get("protein", "Protein")}" AS protein_cell',
        ]
        # Score columns
        for col, _, _ in _SCORE_COLUMNS:
            if col in actual_cols:
                cols.append(f'TRY_CAST("{col}" AS DOUBLE) AS "{col}"')
            else:
                cols.append(f'CAST(NULL AS DOUBLE) AS "{col}"')
        # CV param columns (keep as strings)
        for col in _CV_PARAM_COLUMNS:
            if col in actual_cols:
                cols.append(f'"{col}"')

        return f"SELECT\n    {','.join(chr(10) + '    ' + c for c in cols)}\nFROM psm"

    # ------------------------------------------------------------------
    # Transformation
    # ------------------------------------------------------------------

    def _transform_batch(self, df: pd.DataFrame) -> list[dict]:
        """Vectorised per-row transformation of a DuckDB batch."""
        return self._transform_batch_with_skip(df, self._transform_row, "PSM")

    def _transform_row(self, row: dict) -> Optional[dict]:
        """Convert a single ``.psm`` row into a QPX PSM record.

        Column names match the aliases produced by :meth:`_build_typed_sql`.
        Numeric fields are pre-cast to DOUBLE/INTEGER in SQL but may still be
        NaN for unparseable values; use ``_finite_float`` / ``_optional_float``
        / ``_int_or_zero`` helpers to normalise them safely.
        """
        identity = self._parse_peptide_identity(row)
        if identity is None:
            return None
        peptidoform = identity[1]
        charge = self._int_or_zero(row.get("charge"))
        run_file_name = row.get("run_file_name") or ""

        scan_val = self._int_or_zero(row.get("scan_num"))
        scan: list[int] = [scan_val] if scan_val else []

        observed_mz = self._finite_float(row.get("observed_mz"))
        # ``calculated_mz`` is a THEORETICAL mass. When the peptidoform carries
        # an unparseable token (e.g. a CHEMMOD unknown-mass tag), keep it NULL
        # rather than falling back to a measured precursor m/z, which would
        # silently corrupt the theoretical-mass semantics (bigbio/qpx#250).
        calculated_mz = compute_precursor_mz(peptidoform, charge) if peptidoform and charge else None

        mass_error_ppm = self._optional_float(row.get("mass_error_ppm"))
        rt = self._optional_float(row.get("rt"))

        protein_records, is_decoy = self.parse_protein_cell(row.get("protein_cell"))
        protein_accessions = [rec["accession"] for rec in protein_records] or None

        additional_scores = self._extract_scores(row)
        cv_params = self._extract_cv_params(row)

        rec = self._common_record_fields(
            identity,
            charge,
            is_decoy,
            (calculated_mz, observed_mz, mass_error_ppm, additional_scores),
            run_file_name,
        )
        rec.update(
            {
                "cv_params": cv_params or None,
                "scan": scan,
                "rt": rt,
                "ion_mobility": None,
                "missed_cleavages": None,
                "protein_accessions": protein_accessions,
            }
        )
        return rec

    def _extract_scores(self, row: dict) -> list[dict]:
        """Extract additional_scores from pre-cast DOUBLE columns."""
        return self._extract_scores_from_columns(row, _SCORE_COLUMNS)

    def _extract_cv_params(self, row: dict) -> list[dict]:
        """Extract CV params from string columns."""
        params: list[dict] = []
        for col, cv_name in _CV_PARAM_COLUMNS.items():
            value = row.get(col)
            if value is not None and value != "" and str(value) != "nan":
                params.append({"cv_name": cv_name, "cv_value": str(value)})
        return params
