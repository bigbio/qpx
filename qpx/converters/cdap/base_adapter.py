"""Shared base adapter for CDAP converters.

Centralises helpers reused across :class:`CdapPsmAdapter`,
:class:`CdapFeatureAdapter`, and :class:`CdapPgAdapter`:

    * Loading every ``*.psm`` file inside a study directory into one DuckDB
      table while tolerating per-file optional-column differences.
    * Reporter-ion channel auto-detection from a study's column header.
    * Parsing CDAP's ``ACC(pre=X,post=Y);...`` protein notation.
"""

from __future__ import annotations

import math
import re
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Optional

from qpx.converters.base import BaseConverter
from qpx.converters.cdap.constants import (
    CHANNEL_DEFS,
    DECOY_PREFIX,
)
from qpx.converters.cdap.peptidoform import cdap_to_proforma, strip_label_tags
from qpx.converters.ptm import from_proforma
from qpx.core.sql import sql_build, validate_table

if TYPE_CHECKING:
    import pandas as pd

# CDAP protein cell format: ``ACC(pre=X,post=Y)`` -- pre/post are flanking
# residues used for digestion/peptide-position context.  The accession can
# contain pipes, dots, and underscores so we capture greedy up to the opening
# parenthesis (Ensembl entries use ``ENSP|ENST|ENSG|SYMBOL``).
_PROTEIN_TOKEN_RE = re.compile(r"^(.+?)\(pre=([^,]+),post=([^)]+)\)$")

# Channel header detection: token must start with one of the known prefixes.
_CHANNEL_PREFIX_RE = re.compile(r"^(?:TMT\d+|TMTXX|iTRAQ)")


class CdapBaseAdapter(BaseConverter):
    """Base adapter with CDAP-specific loading and parsing utilities."""

    # ------------------------------------------------------------------
    # Loading: many *.psm files -> single DuckDB table
    # ------------------------------------------------------------------

    def _load_psm_directory(self, psm_dir: str, table_name: str = "psm") -> int:
        """Load every ``*.psm`` file in *psm_dir* into a DuckDB table.

        Uses ``union_by_name=true`` so studies with optional columns
        (``PhosphoRSPeptide``, ``nPhospho``, ``FullyLocalized``) merge cleanly
        with files that omit them.  Returns the loaded row count.

        When a shared DuckDB connection is used, the table may already exist
        from a prior adapter; in that case we skip the reload and just log.
        """
        tbl = validate_table(table_name)
        count_sql = sql_build("SELECT COUNT(*) FROM $tbl", tbl=tbl)

        if self._table_exists(table_name):
            count = self._conn.execute(count_sql).fetchone()[0]
            self.logger.info("Reusing %d CDAP PSM rows already in '%s'", count, table_name)
            return count

        psm_dir = str(Path(psm_dir).resolve())
        pattern = f"{psm_dir}/*.psm"
        self._conn.execute(
            sql_build(
                "CREATE OR REPLACE TABLE $tbl AS "
                "SELECT * FROM read_csv($1, delim='\\t', header=true, "
                "auto_detect=true, null_padding=true, union_by_name=true, "
                "ignore_errors=false, sample_size=-1)",
                tbl=tbl,
            ),
            [pattern],
        )
        count = self._conn.execute(count_sql).fetchone()[0]
        self.logger.info("Loaded %d CDAP PSM rows from %s", count, psm_dir)
        return count

    # ------------------------------------------------------------------
    # Channel auto-detection from a study's column header
    # ------------------------------------------------------------------

    def _detect_channels(self, columns: Iterable[str]) -> tuple[list[str], str]:
        """Identify reporter-ion data channels present in *columns*.

        Returns a ``(channel_columns, experiment_type)`` tuple where
        ``experiment_type`` is one of ``"TMT10"``, ``"TMT11"``, ``"TMT18"``,
        ``"iTRAQ4"``, or ``"LFQ"`` (no reporter channels detected).
        """
        col_set = set(columns)
        # Walk the known channel definitions in mass-ascending order so that
        # the returned list preserves the canonical layout that downstream
        # consumers (SDRF resolution, intensity aggregation) rely on.
        for prefix, channel_list in CHANNEL_DEFS.items():
            present = [c for c in channel_list if c in col_set]
            if present and prefix != "TMTXX":
                # If a TMT11 study has TMTXX-132N/C as well, append them.
                if prefix == "TMT11":
                    present.extend(c for c in CHANNEL_DEFS["TMTXX"] if c in col_set)
                return present, prefix
        return [], "LFQ"

    # ------------------------------------------------------------------
    # Value parsers (protein)
    # ------------------------------------------------------------------

    @staticmethod
    @lru_cache(maxsize=100_000)
    def _parse_protein_text(text: str) -> tuple[list[dict], bool]:
        """Parse a non-empty CDAP protein cell into ``(records, is_decoy)``.

        Bounded LRU cache: the same protein cell string recurs across many PSMs
        in a study, but the cache must stay bounded so long-lived / multi-study
        runs do not accumulate memory without limit.
        """
        records: list[dict] = []
        decoy_flags: list[bool] = []
        for token in text.split(";"):
            token = token.strip()
            if not token:
                continue
            match = _PROTEIN_TOKEN_RE.match(token)
            if match:
                accession, pre, post = match.group(1), match.group(2), match.group(3)
            else:
                # Fall back to the raw token (no flanking residue annotation).
                accession, pre, post = token, None, None
            is_decoy_token = accession.startswith(DECOY_PREFIX)
            if is_decoy_token:
                accession = accession[len(DECOY_PREFIX) :]
            decoy_flags.append(is_decoy_token)
            records.append({"accession": accession, "pre": pre, "post": post})

        is_decoy = bool(records) and all(decoy_flags)
        return records, is_decoy

    @classmethod
    def parse_protein_cell(cls, cell: object) -> tuple[list[dict], bool]:
        """Parse CDAP ``Protein`` cell into protein records + decoy flag.

        Returns ``(records, is_decoy)`` where ``records`` is a list of
        ``{"accession": str, "pre": str|None, "post": str|None}`` and
        ``is_decoy`` is True only when *every* listed accession carries the
        :data:`DECOY_PREFIX`.  Parsing of non-empty cells is memoised via a
        bounded LRU cache (see :meth:`_parse_protein_text`).
        """
        if cell is None:
            return [], False
        text = str(cell).strip()
        if not text:
            return [], False
        return cls._parse_protein_text(text)

    @staticmethod
    def strip_run_extension(filename: object) -> str:
        """Return the run stem with the trailing ``.raw``/``.mzML`` removed."""
        if filename is None:
            return ""
        stem = str(filename).strip()
        for ext in (".raw", ".RAW", ".mzML", ".mzml", ".mgf"):
            if stem.endswith(ext):
                return stem[: -len(ext)]
        return stem

    # ------------------------------------------------------------------
    # Shared record field builder
    # ------------------------------------------------------------------

    @staticmethod
    def _common_record_fields(identity, charge, is_decoy, mz_values, run_file_name) -> dict:
        """Build the common field dict shared by PSM and Feature records.

        Args:
            identity: ``(sequence, peptidoform, modifications)`` tuple from
                :meth:`_parse_peptide_identity`.
            charge: Precursor charge state.
            is_decoy: Whether the record is a decoy hit.
            mz_values: ``(calculated_mz, observed_mz, mass_error_ppm,
                additional_scores)`` tuple.
            run_file_name: Source run file name.
        """
        sequence, peptidoform, modifications = identity
        calculated_mz, observed_mz, mass_error_ppm, additional_scores = mz_values
        return {
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "charge": charge,
            "posterior_error_probability": None,
            "is_decoy": is_decoy,
            "calculated_mz": float(calculated_mz) if calculated_mz is not None else None,
            "observed_mz": observed_mz,
            "mass_error_ppm": mass_error_ppm,
            "additional_scores": additional_scores or None,
            "predicted_rt": None,
            "run_file_name": run_file_name,
        }

    # ------------------------------------------------------------------
    # Shared score extraction
    # ------------------------------------------------------------------

    def _extract_scores_from_columns(self, row: dict, score_defs: Iterable[tuple[str, str, bool]]) -> list[dict]:
        """Build additional_scores list from *score_defs*.

        Each element of *score_defs* is ``(column_key, score_name, higher_better)``.
        """
        scores: list[dict] = []
        for col_key, score_name, higher_better in score_defs:
            value = self._optional_float(row.get(col_key))
            if value is not None:
                scores.append({"score_name": score_name, "score_value": value, "higher_better": higher_better})
        return scores

    # ------------------------------------------------------------------
    # Shared peptide identity parsing
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_peptide_identity(row: dict) -> tuple[str, str, object] | None:
        """Parse peptide_raw from *row*, returning (sequence, peptidoform, modifications).

        Returns ``None`` if the row should be skipped (empty/invalid peptide).
        """
        peptide_raw = row.get("peptide_raw")
        if not isinstance(peptide_raw, str) or not peptide_raw:
            return None
        sequence = strip_label_tags(peptide_raw)
        if not sequence:
            return None
        peptidoform = cdap_to_proforma(peptide_raw)
        _, modifications = from_proforma(peptidoform, sequence) if peptidoform else (peptidoform, None)
        return sequence, peptidoform, modifications

    # ------------------------------------------------------------------
    # Stream-transform-write loop (shared across adapters)
    # ------------------------------------------------------------------

    def _stream_transform_write(self, sql: str, writer, chunksize: int) -> None:
        """Execute *sql* in batches, transform, track scores, and write.

        Subclasses override ``_transform_batch`` to provide the actual
        per-batch transformation logic.
        """
        for batch in self._query_batched(sql, chunksize):
            df = batch.to_pandas()
            records = self._transform_batch(df)
            if records:
                self._track_scores(records)
                writer.write_batch(records)

    def _transform_batch(self, df: "pd.DataFrame") -> list[dict]:
        """Override in subclass."""
        raise NotImplementedError

    # ------------------------------------------------------------------
    # Batch transformation with skip-counting (shared across adapters)
    # ------------------------------------------------------------------

    def _transform_batch_with_skip(
        self,
        df: "pd.DataFrame",
        transform_fn,
        row_label: str = "row",
    ) -> list[dict]:
        """Apply *transform_fn* to each row of *df*, tolerating per-row errors.

        Returns the list of successfully transformed records.  Rows that raise
        or return ``None`` are counted and logged as warnings.

        Args:
            df: pandas DataFrame (one DuckDB batch).
            transform_fn: Callable accepting a ``dict`` row and returning a
                record dict or ``None``.
            row_label: Label used in log messages (e.g. "PSM", "feature").
        """
        records: list[dict] = []
        skipped = 0
        col_arrays = {col: df[col].values for col in df.columns}
        n_rows = len(df)
        for i in range(n_rows):
            try:
                row = {col: vals[i] for col, vals in col_arrays.items()}
                rec = transform_fn(row)
                if rec is not None:
                    records.append(rec)
            except (ValueError, TypeError, KeyError, IndexError, ArithmeticError) as exc:
                skipped += 1
                self.logger.debug("Skipping CDAP %s row: %s", row_label, exc)
        if skipped:
            total = skipped + len(records)
            self.logger.warning(
                "Skipped %d / %d %s rows (%.1f%%) in batch",
                skipped,
                total,
                row_label,
                100 * skipped / total if total else 0.0,
            )
        return records

    # ------------------------------------------------------------------
    # NaN-safe value extraction helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _finite_float(value, default: float = 0.0) -> float:
        """Return *value* as float if finite, otherwise *default*.

        DuckDB TRY_CAST may produce ``None`` (SQL NULL) or ``float('nan')``
        for unparseable/missing values; this helper normalises both cases.
        """
        if value is None:
            return default
        try:
            f = float(value)
        except (TypeError, ValueError):
            return default
        if math.isnan(f):
            return default
        return f

    @staticmethod
    def _optional_float(value) -> Optional[float]:
        """Return *value* as float if finite, otherwise ``None``."""
        if value is None:
            return None
        try:
            f = float(value)
        except (TypeError, ValueError):
            return None
        if math.isnan(f):
            return None
        return f

    @staticmethod
    def _int_or_zero(value) -> int:
        """Return *value* as int if valid, otherwise ``0``."""
        if value is None:
            return 0
        try:
            f = float(value)
        except (TypeError, ValueError):
            return 0
        if math.isnan(f):
            return 0
        return int(f)
