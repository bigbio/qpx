"""DIA-NN Feature adapter -- report to feature.parquet.

Performance strategy:
    - DuckDB VIEW over parquet for memory efficiency (no in-memory table)
    - Precursor lookup (peptidoform, mz, missed_cleavages) via DuckDB temp table JOIN
    - All scalar and nested fields (scores, intensities, cv_params) constructed in SQL
    - Only modifications built in Python (requires PTM parsing); pg_accessions
      split with Arrow kernels
    - Arrow tables written directly via write_table() — no Python dict intermediaries
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from functools import partial
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from qpx.converters.base import resolve_columns
from qpx.converters.diann.base_adapter import BATCH_ROWS, DiaNNBaseAdapter, run_name_sql
from qpx.converters.diann.constants import occupancy_sites, to_modifications, to_proforma
from qpx.converters.diann.feature_sql import _SCORE_MAPPINGS, _build_batch_sql, build_pepmap_sql
from qpx.converters.diann.fragment_sql import build_fragment_sql
from qpx.converters.fragments import write_fragment_view
from qpx.converters.localization import LOCALIZATION_SCORE, localize_modifications
from qpx.converters.mappings import get_field_mappings
from qpx.converters.utils import pg_accessions_array
from qpx.core.query import _escape_sql_string
from qpx.core.scan import scan_from_native_id
from qpx.core.sql import sql_build, validate_identifier
from qpx.writers.feature import FeatureWriter
from qpx.writers.pepmap import PepMapWriter

logger = logging.getLogger(__name__)

# Producer-specific feature identity composite (bigbio/qpx#229). DIA-NN reports one
# feature per precursor per run, so the measured key is just peptidoform + charge +
# run (all in feature.yaml). Passed to FeatureWriter so feature_id hashes exactly
# these columns instead of the schema default.
_FEATURE_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name")


def _ms_info_scans(values: pd.Series) -> pd.Series:
    """Scan numbers of an MS-info table: numbers, or one-component native IDs such as ``index=5``."""
    scans = pd.to_numeric(values, errors="coerce")
    native = scans.isna() & values.notna()
    if native.any():
        components = values[native].map(lambda value: scan_from_native_id(str(value)))
        scans[native] = components.map(lambda parts: parts[0] if len(parts) == 1 else None)
    return scans


def _nearest_ms2(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    """Nearest MS2 entry by RT for each feature row.

    DIA spectra acquired together (the windows of one diaPASEF frame) share one
    time; among them the window centre nearest the precursor m/z
    (``precursor_mz``) is the feature's fragment spectrum.
    """
    times = right[["_ms_rt"]].drop_duplicates()
    nearest = pd.merge_asof(left, times, left_on="rt", right_on="_ms_rt", direction="nearest")
    nearest = nearest.merge(right, on="_ms_rt", how="left")
    gap = (nearest["_matched_mz"] - nearest["precursor_mz"].where(nearest["precursor_mz"] > 0)).abs()
    nearest = nearest.assign(_gap=gap).sort_values(["_merge_row", "_gap"], kind="stable", na_position="last")
    return nearest.drop_duplicates("_merge_row").drop(columns="_gap")


class DiannFeatureAdapter(DiaNNBaseAdapter):
    """Convert DIA-NN report to ``feature.parquet``.

    Usage::

        with DiannFeatureAdapter() as adapter:
            adapter.convert(
                diann_report="report.tsv",
                output_path="feature.parquet",
                mzml_info_folder="ms_info/",
                sdrf_path="sdrf.tsv",
            )
    """

    # The report column naming channels of a plexDIA report, read by write_fragments().
    _channel_col: Optional[str] = None

    # Pepmap query of the converted runs, run by write_pepmap().
    _pepmap_query: str | None = None

    def convert(
        self,
        diann_report: str,
        output_path: str,
        mzml_info_folder: Optional[str] = None,
        sdrf_path: Optional[str] = None,
        qvalue_threshold: Optional[float] = None,
        batch_rows: int = BATCH_ROWS,
        creator: str = "diann",
        mod_masses: tuple[tuple[str, float], ...] = (),
    ) -> None:
        """Run the DIA-NN report -> feature.parquet conversion.

        Filtering is opt-in (bigbio/qpx#241). ``qvalue_threshold`` defaults to
        ``None``, which emits every feature the report contains — DIA-NN already
        FDR-filters its main report and all per-row q-value columns (``precursor_qvalue``,
        ``global_qvalue``, ``diann_lib_qvalue``, ``diann_translated_qvalue`` …)
        are carried through unconditionally for downstream filtering. When a
        threshold is provided, the precursor ``Q.Value`` filter is applied.
        ``mod_masses`` holds the ``(name, mass delta)`` modifications declared in
        the DIA-NN log, which resolve modifications the report names.
        Runs are converted in batches of about ``batch_rows`` report rows.
        """
        # 1. Create VIEW over report (lazy — no data loaded into memory)
        self._load_diann_report(diann_report)

        # 2. SDRF for sample mapping and enzyme info
        enzyme_name: str | None = None
        if sdrf_path:
            enzyme_name = self._load_sdrf_enzyme(sdrf_path)

        # 3. Detect report columns (once, reused for all batches)
        report_cols = self._get_report_columns()
        self._resolved = resolve_columns(get_field_mappings("diann", "feature"), report_cols)
        decoy_col = self._resolved.get("decoy")
        has_decoy_col = decoy_col is not None
        channel_col = self._register_channel_labels(report_cols, sdrf_path)
        self._channel_col = channel_col

        # 4. Pre-compute precursors: DuckDB temp table + Python modifications dict
        mods_lookup = self._register_precursor_lookup(
            self._resolved,
            enzyme_name,
            partial(to_proforma, mod_masses=mod_masses),
            partial(to_modifications, mod_masses=mod_masses),
        )
        self._mod_masses = mod_masses

        # 5. Get target Arrow schema
        target_schema = FeatureWriter._schema_class.get_arrow_schema()

        # 6. Track discovered scores upfront (all score columns present in report)
        for diann_col, score_name, _ in _SCORE_MAPPINGS:
            if diann_col in report_cols:
                self._discovered_scores.add(score_name)

        # 7. Build SQL template (once, reused for all batches)
        sql_template = _build_batch_sql(
            self._resolved,
            report_cols,
            has_decoy_col,
            qvalue_threshold,
            channel_col,
            target_schema,
        )

        # 8. Discover runs
        run_names = self._discover_runs(mzml_info_folder)
        self._pepmap_query = self._fill_runs(build_pepmap_sql(self._resolved, report_cols, qvalue_threshold), run_names)

        # 9. Process in batches → Arrow tables → write directly
        with FeatureWriter(
            output_path,
            creator=creator,
            compression=self._compression,
            identity_composite=_FEATURE_IDENTITY_COMPOSITE,
        ) as writer:
            for table in self._feature_batches(run_names, batch_rows, sql_template, mods_lookup, target_schema, mzml_info_folder):
                writer.write_table(table)

        self.logger.info(f"DIA-NN feature conversion complete -> {output_path}")

    def write_fragments(self, output_path: str, qvalue_threshold: Optional[float] = None) -> bool:
        """Write the ``Fragment.Quant.*`` lists of the report as the fragment view; False when it has none."""
        sql = build_fragment_sql(self._resolved, self._get_report_columns(), self._channel_col, qvalue_threshold)
        if sql is None:
            self.logger.warning("The report has no Fragment.Quant.Raw column: no fragment view")
            return False
        rows = write_fragment_view(
            self._conn, sql, output_path, _FEATURE_IDENTITY_COMPOSITE, creator="diann", compression=self._compression
        )
        self.logger.info("Wrote %d fragment rows to %s", rows, output_path)
        return True

    def write_pepmap(self, output_path: str, creator: str = "diann") -> bool:
        """Write the proteins each converted peptidoform maps to (``Protein.Ids``) as the pepmap view.

        ``pg_accessions`` of the feature view holds the peptide's protein group;
        ``Protein.Ids`` lists every protein the peptide maps to. Returns whether a
        file was written: not for a report without ``Protein.Ids`` or features.
        """
        rows = self._conn.execute(self._pepmap_query).fetchall() if self._pepmap_query else []
        if not rows:
            return False
        records = [
            {
                "sequence": sequence,
                "peptidoform": peptidoform,
                "pg_accessions": [{"accession": acc, "start": None, "end": None, "pre": None, "post": None} for acc in accessions]
                or None,
                "is_unique": len(accessions) == 1 if accessions else None,
            }
            for peptidoform, sequence, accessions in rows
        ]
        with PepMapWriter(output_path, creator=creator, compression=self._compression) as writer:
            writer.write_batch(records)
        self.logger.info("Wrote %d pepmap entries to %s", len(records), output_path)
        return True

    def _fill_runs(self, sql: str | None, run_names: list[str]) -> str | None:
        """*sql* restricted to the report runs of *run_names*; None without SQL or runs."""
        runs = self._report_runs(run_names) if sql else []
        return sql.replace("{run_placeholders}", ", ".join(f"'{_escape_sql_string(run)}'" for run in runs)) if runs else None

    def _feature_batches(
        self,
        run_names: list[str],
        batch_rows: int,
        sql_template: str,
        mods_lookup: dict,
        target_schema: pa.Schema,
        mzml_info_folder: Optional[str],
    ) -> Iterator[pa.Table]:
        """Yield assembled feature tables for consecutive batches of runs."""
        run_col = validate_identifier(self._resolved["run_file_name"])
        done = 0
        for batch_runs in self._run_batches(run_col, self._report_runs(run_names), batch_rows):
            self.logger.info(f"Processing runs {done + 1}-{done + len(batch_runs)} of {len(run_names)}")
            done += len(batch_runs)
            table = self._process_batch_arrow(
                batch_runs,
                sql_template,
                mods_lookup,
                target_schema,
                mzml_info_folder,
            )
            if table is not None:
                # DuckDB's parallel GROUP BY returns rows in a different order on each run.
                yield table.sort_by([(column, "ascending") for column in _FEATURE_IDENTITY_COMPOSITE])

    # ------------------------------------------------------------------
    # Data loading helpers
    # ------------------------------------------------------------------

    def _get_report_columns(self) -> set[str]:
        """Return the set of column names in the DuckDB report table/view."""
        return {
            c[0]
            for c in self._conn.execute("SELECT column_name FROM information_schema.columns WHERE table_name='report'").fetchall()
        }

    def _discover_runs(self, mzml_info_folder: Optional[str]) -> list[str]:
        """Discover run names from ms_info files or from the report."""
        if mzml_info_folder:
            info_files = list(Path(mzml_info_folder).glob("*_ms_info.parquet"))
        else:
            info_files = []

        if info_files:
            run_names = [f.stem.replace("_ms_info", "") for f in info_files]
            self.logger.info(f"Found {len(run_names)} MS info files")
            # The feature SQL filters `run IN (run_names)`, so any report Run
            # without a matching *_ms_info.parquet is silently dropped from
            # feature.parquet. Surface the discrepancy (bigbio/qpx#244) so the
            # loss is visible; the drop behaviour itself is unchanged.
            self._warn_on_dropped_runs(set(run_names))
        else:
            run_names = sorted(self._report_run_names())
            self.logger.info(f"Discovered {len(run_names)} runs from report")

        return run_names

    def _report_runs(self, run_names: list[str]) -> list[str]:
        """The report ``Run`` values of *run_names*, which lack the raw-file extension, in that order.

        Batches filter the ``Run`` column itself, which DuckDB pushes into the scan.
        """
        run_col = validate_identifier(self._resolved["run_file_name"])
        by_name: dict[str, list[str]] = {}
        rows = self._conn.execute(
            sql_build("SELECT DISTINCT $run, $name FROM report ORDER BY 1", run=run_col, name=run_name_sql(run_col))
        ).fetchall()
        for raw, name in rows:
            by_name.setdefault(name, []).append(raw)
        return [raw for name in run_names for raw in by_name.get(name, [])]

    def _report_run_names(self) -> set[str]:
        """Return the distinct extension-stripped Run values in the report."""
        run_col = validate_identifier(self._resolved["run_file_name"])
        rows = self._conn.execute(
            sql_build("SELECT DISTINCT $run AS run_file_name FROM report", run=run_name_sql(run_col))
        ).fetchall()
        return {row[0] for row in rows}

    def _warn_on_dropped_runs(self, kept_runs: set[str]) -> None:
        """Warn about report runs with no matching ms_info file.

        When runs are discovered from ``*_ms_info.parquet`` stems, the feature
        SQL keeps only those runs. Any report ``Run`` absent from *kept_runs* is
        dropped from ``feature.parquet`` without reconciliation; log the dropped
        run names so the loss is visible (bigbio/qpx#244).
        """
        dropped = sorted(self._report_run_names() - kept_runs)
        if dropped:
            self.logger.warning(
                "%d report run(s) have no matching *_ms_info.parquet and will be dropped from feature.parquet: %s",
                len(dropped),
                ", ".join(dropped),
            )

    # ------------------------------------------------------------------
    # Batch processing (SQL → Arrow → write)
    # ------------------------------------------------------------------

    def _process_batch_arrow(
        self,
        runs: list[str],
        sql_template: str,
        mods_lookup: dict,
        target_schema: pa.Schema,
        mzml_info_folder: Optional[str],
    ) -> Optional[pa.Table]:
        """Execute SQL batch, add Python-only columns, cast to target schema."""
        placeholders = ", ".join(f"'{_escape_sql_string(r)}'" for r in runs)
        sql = sql_template.replace("{run_placeholders}", placeholders)

        # Execute SQL — DuckDB pushes predicates to parquet, only reads needed data
        table = self._conn.sql(sql).to_arrow_table()

        if len(table) == 0:
            return None

        # --- Merge scan info from mzml files if available ---
        if mzml_info_folder:
            table = self._merge_scan_info_arrow(table, mzml_info_folder)

        # --- Build modifications column from Python lookup ---
        mod_seqs = table.column("_modified_sequence").to_pylist()
        seqs = table.column("sequence").to_pylist()
        charges = table.column("charge").to_pylist()

        mods_list = [mods_lookup.get((str(ms), str(s), int(c))) for ms, s, c in zip(mod_seqs, seqs, charges)]
        mods_list, candidates = self._localized(mods_list, seqs, table.column("_site_occupancy").to_pylist())
        mods_array = pa.array(mods_list, type=target_schema.field("modifications").type)
        candidates_type = target_schema.field("modification_candidates").type
        table = table.append_column("modification_candidates", pa.array(candidates, type=candidates_type))

        # --- Build pg_accessions column ---
        pg_array = pg_accessions_array(table.column("_pg_group"), target_schema.field("pg_accessions").type)

        # feature.pg_ids is left null: the feature->pg association is computed on
        # read via the Dataset softlink (bigbio/qpx#269), not materialized here.

        # --- Drop helper columns ---
        table = table.drop("_modified_sequence")
        table = table.drop("_pg_group")
        table = table.drop("_site_occupancy")

        return self._assemble_feature_table(table, target_schema, mods_array, pg_array)

    def _localized(self, modifications: list, sequences: list, occupancies: list) -> tuple[list, list]:
        """Modifications with their site scores, and the candidates, from each row's occupancy probabilities."""
        localized = [
            localize_modifications(mods, sequence, occupancy_sites(occupancy, sequence, self._mod_masses))
            if occupancy and "{" in occupancy
            else (mods, None)
            for mods, sequence, occupancy in zip(modifications, sequences, occupancies)
        ]
        candidates = [row_candidates for _, row_candidates in localized]
        if any(candidates):
            self._discovered_scores.add(LOCALIZATION_SCORE)
        return [mods for mods, _ in localized], candidates

    def _merge_scan_info_arrow(self, table: pa.Table, mzml_info_folder: str) -> pa.Table:
        """Merge MS scan info from mzml files into the Arrow table.

        Only the columns the merge reads go through pandas; the merged ``scan``
        column replaces the table's, with rows grouped by run.
        """
        df = table.select(["run_file_name", "rt", "calculated_mz", "observed_mz", "scan", "_ms2_index"]).to_pandas()
        rows: list[np.ndarray] = []
        scans: list = []

        for run_name, group_df in df.groupby("run_file_name", sort=False):
            ms_info_path = next(
                Path(mzml_info_folder).glob(f"*{run_name}_ms_info.parquet"),
                None,
            )
            if ms_info_path:
                group_df = self._merge_scan_info(group_df, ms_info_path)
            rows.append(group_df.index.to_numpy())
            scans.extend(group_df["scan"])

        merged = table.take(pa.array(np.concatenate(rows), type=pa.int64()))
        scan_field = merged.schema.field("scan")
        return merged.set_column(merged.schema.get_field_index("scan"), scan_field, pa.array(scans, type=scan_field.type))

    def _merge_scan_info(self, run_data: pd.DataFrame, ms_info_path: Path) -> pd.DataFrame:
        """Fill feature scans from one run's MS-info table.

        A DIA-NN ``MS2.Scan`` (``Ms2.Scan`` in 2.x) position selects that run's MS2 spectrum in
        acquisition order. Rows without one take the nearest MS2 entry by RT (by
        window among entries sharing that time), never past the measured
        acquisition range. The window is chosen by the theoretical m/z, or the
        measured one when the peptidoform has no computable mass. The table's
        precursor m/z is an isolation-window centre for DIA spectra, so it never
        stands in for a feature's observed m/z.
        """
        columns = [name for name in ("rt", "scan", "precursor_mz", "ms_level") if name in pq.read_schema(ms_info_path).names]
        target = pd.read_parquet(ms_info_path, columns=columns)
        run_data = run_data.copy()
        run_data["_merge_row"] = range(len(run_data))
        run_data["rt"] = run_data["rt"].astype("float64")
        run_data["precursor_mz"] = (
            pd.to_numeric(run_data["calculated_mz"], errors="coerce")
            .fillna(pd.to_numeric(run_data["observed_mz"], errors="coerce"))
            .astype("float64")
        )
        target["scan"] = _ms_info_scans(target["scan"])
        target["precursor_mz"] = pd.to_numeric(target["precursor_mz"], errors="coerce")
        target["rt"] = target["rt"].astype("float64")

        ms2 = target[target["ms_level"] == 2] if "ms_level" in target else target[target["precursor_mz"].notna()]
        positions = pd.Series(ms2["scan"].to_numpy())
        indexed_scan = pd.to_numeric(run_data["_ms2_index"], errors="coerce").astype("Int64").map(positions)
        has_index = indexed_scan.notna()
        run_data.loc[has_index, "scan"] = indexed_scan[has_index].map(lambda value: [int(value)])

        missing_scan = ~has_index & run_data["rt"].notna()
        rt_target = target.dropna(subset=["rt", "scan", "precursor_mz"]).sort_values("rt")
        if missing_scan.any() and not rt_target.empty:
            left = run_data.loc[missing_scan, ["_merge_row", "rt", "precursor_mz"]].sort_values("rt")
            right = rt_target[["rt", "scan", "precursor_mz"]].rename(
                columns={"rt": "_ms_rt", "scan": "_matched_scan", "precursor_mz": "_matched_mz"}
            )
            nearest = _nearest_ms2(left, right)
            outside_acquisition = (nearest["rt"] < right["_ms_rt"].min()) | (nearest["rt"] > right["_ms_rt"].max())
            nearest.loc[outside_acquisition, "_matched_scan"] = None

            matched_scan = run_data["_merge_row"].map(nearest.set_index("_merge_row")["_matched_scan"])
            has_match = matched_scan.notna()
            run_data.loc[has_match, "scan"] = matched_scan[has_match].map(lambda value: [int(value)])

        return run_data.drop(columns=["_merge_row", "precursor_mz"])
