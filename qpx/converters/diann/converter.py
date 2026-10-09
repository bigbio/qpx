"""DIA-NN orchestrator — composes Feature and PG adapters."""

from __future__ import annotations

import logging
import re
from pathlib import Path

from qpx._version import __version__
from qpx.converters.base import resolve_columns
from qpx.converters.diann.feature_adapter import DiannFeatureAdapter
from qpx.converters.diann.pg_adapter import DiannPgAdapter
from qpx.converters.mappings import get_field_mappings, get_tool_meta
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.core.constants import FEATURE, ONTOLOGY, PG, RUN, SAMPLE
from qpx.core.scores import field_ontology_entries, score_ontology_entries

logger = logging.getLogger(__name__)


class DiaNNConverter(BaseOrchestrator):
    """Orchestrate full DIA-NN conversion to QPX format."""

    _DIANN_VERSION_RE = re.compile(r"DIA-NN\s+([\d.]+)")

    def __init__(
        self,
        report_path,
        sdrf_path=None,
        max_memory=None,
        max_cpus=None,
        compression: str = "zstd",
        diann_log: str | None = None,
    ):
        self.report_path = str(report_path)
        self.sdrf_path = str(sdrf_path) if sdrf_path else None
        self._memory = max_memory or "16GB"
        self._threads = max_cpus or 4
        self._compression = compression
        self._diann_version = self._parse_diann_version(diann_log)
        self._ontology_entries: list[dict] = []
        self._resolved_mappings_by_view: dict[str, dict] = {}

    def convert_features(
        self,
        mzml_info_folder=None,
        qvalue_threshold=None,
        output_folder=".",
        output_prefix=None,
        protein_file=None,
        partitions=None,
        batch_size=100,
    ):
        output_folder = Path(output_folder)
        prefix = output_prefix or "diann"
        with DiannFeatureAdapter(
            duckdb_memory=self._memory,
            duckdb_threads=self._threads,
            compression=self._compression,
        ) as adapter:
            adapter.convert(
                diann_report=self.report_path,
                output_path=str(output_folder / f"{prefix}.feature.parquet"),
                mzml_info_folder=str(mzml_info_folder) if mzml_info_folder else None,
                sdrf_path=self.sdrf_path,
                qvalue_threshold=qvalue_threshold,
            )
            self._ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=FEATURE))
            cols = adapter.get_table_columns("report")
            self._resolved_mappings_by_view[FEATURE] = resolve_columns(get_field_mappings("diann", "feature"), cols)
        logger.info("DIA-NN feature conversion complete")

    def convert_pg(
        self,
        pg_matrix_path,
        output_folder=".",
        output_prefix=None,
        batch_size=100,
        standardized_intensities=False,
        qvalue_threshold=None,
    ):
        output_folder = Path(output_folder)
        prefix = output_prefix or "diann"
        with DiannPgAdapter(
            duckdb_memory=self._memory,
            duckdb_threads=self._threads,
            compression=self._compression,
        ) as adapter:
            adapter.convert(
                diann_report=self.report_path,
                pg_matrix_path=str(pg_matrix_path),
                sdrf_path=self.sdrf_path,
                output_path=str(output_folder / f"{prefix}.pg.parquet"),
                qvalue_threshold=qvalue_threshold,
            )
            self._ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=PG))
            cols = adapter.get_table_columns("report")
            self._resolved_mappings_by_view[PG] = resolve_columns(get_field_mappings("diann", "pg"), cols)
        logger.info("DIA-NN PG conversion complete")

    def convert_mz(
        self,
        ams_info_folder,
        output_folder=".",
        output_prefix=None,
        matched_scans=None,
    ):
        """Generate the ``mz`` view from mzML-derived AMS_info spectra.

        Writes only the feature-matched MS2 scans. When ``matched_scans`` (a set of
        ``(run_file_name, scan)`` tuples) is omitted it is derived from the DIA-NN
        report's ``Run`` / ``MS2.Scan`` columns. Peak arrays are losslessly encoded
        into ``mz_blob`` / ``intensity_blob``.
        """
        from qpx.converters.diann.mz_adapter import DiannMzAdapter

        output_folder = Path(output_folder)
        prefix = output_prefix or "diann"
        if matched_scans is None:
            matched_scans = self._extract_matched_scans()
        adapter = DiannMzAdapter(compression=self._compression)
        adapter.convert(
            ams_info_folder=str(ams_info_folder),
            output_path=str(output_folder / f"{prefix}.mz.parquet"),
            matched_scans=matched_scans,
        )
        logger.info("DIA-NN mz conversion complete")

    def _extract_matched_scans(self) -> set[tuple[str, int]] | None:
        """Return the ``(run_file_name, scan)`` pairs with a DIA feature match.

        Resolves the report's run and MS2 scan columns via the DIA-NN column
        mapping, then returns the distinct pairs. Returns ``None`` when the columns
        cannot be resolved (callers then keep all scans).
        """
        import duckdb

        from qpx.core.sql import escape_path

        con = duckdb.connect()
        try:
            if self.report_path.endswith(".parquet"):
                reader = f"read_parquet('{escape_path(self.report_path)}')"
            else:
                reader = f"read_csv_auto('{escape_path(self.report_path)}', delim='\t', header=true, auto_detect=true)"
            cols = {row[0] for row in con.execute(f"DESCRIBE SELECT * FROM {reader}").fetchall()}
            resolved = resolve_columns(get_field_mappings("diann", "feature"), cols)
            run_col = resolved.get("run_file_name")
            scan_col = resolved.get("ms2_scan")
            if run_col is None or scan_col is None:
                logger.warning("could not resolve Run / MS2.Scan columns; keeping all scans")
                return None
            rows = con.execute(f'SELECT DISTINCT "{run_col}", "{scan_col}" FROM {reader}').fetchall()
            matched = {(str(run), int(scan)) for run, scan in rows if run is not None and scan is not None}
            logger.info("matched %s (run_file_name, scan) from DIA-NN report", len(matched))
            return matched
        finally:
            con.close()

    def convert_sdrf(
        self,
        output_folder: str | Path,
        prefix: str = "diann",
    ) -> None:
        """Convert SDRF to sample.parquet and run.parquet."""
        output_folder = Path(output_folder)
        if not self.sdrf_path:
            logger.warning("No SDRF path provided — skipping sample/run conversion")
            return
        try:
            from qpx.converters.sdrf import SdrfConverter

            with SdrfConverter() as sdrf_conv:
                sdrf_conv.convert(
                    sdrf_path=self.sdrf_path,
                    sample_output=str(output_folder / f"{prefix}.sample.parquet"),
                    run_output=str(output_folder / f"{prefix}.run.parquet"),
                )
                self._ontology_entries.extend(sdrf_conv.run_ontology_entries())
            logger.info("SDRF conversion complete (sample + run)")
        except Exception as exc:
            logger.warning("SDRF conversion skipped (incomplete SDRF?): %s", exc)
            for suffix in (".sample.parquet", ".run.parquet"):
                corrupt = output_folder / f"{prefix}{suffix}"
                if corrupt.exists():
                    corrupt.unlink()
                    logger.debug("Removed corrupt %s", corrupt)

    def write_ontology(self, output_folder: str | Path, prefix: str = "diann") -> None:
        """Write combined ontology.parquet with all accumulated entries."""
        entries = list(self._ontology_entries)
        for view_name, mappings in self._resolved_mappings_by_view.items():
            entries.extend(
                field_ontology_entries(
                    view=view_name,
                    resolved_mappings=mappings,
                    tool_name=get_tool_meta("diann")["tool_name"],
                )
            )
        self._write_ontology(Path(output_folder), prefix, entries)

    def write_provenance(self, output_folder: str | Path, prefix: str = "diann") -> None:
        """Write provenance.parquet with DIA-NN + QPX conversion steps."""
        records = [
            {
                "step_order": 1,
                "step_category": "quantification",
                "step_name": "precursor_quantification",
                "tool_name": "DIA-NN",
                "tool_version": self._diann_version,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": [FEATURE, PG],
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "diann_to_qpx",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": [
                    {"key": "report_path", "value": Path(self.report_path).name},
                ],
                "config": None,
                "output_views": [FEATURE, PG, SAMPLE, RUN, ONTOLOGY],
            },
        ]
        self._write_provenance(Path(output_folder), prefix, records)

    def write_dataset(
        self,
        output_folder: str | Path,
        prefix: str = "diann",
        project_accession: str | None = None,
    ) -> None:
        """Write dataset.parquet with project-level metadata."""
        self._write_dataset(
            Path(output_folder),
            prefix,
            project_accession,
            software_name="DIA-NN",
            software_version=self._diann_version,
        )

    @classmethod
    def _parse_diann_version(cls, log_path: str | None) -> str | None:
        """Extract the first DIA-NN version found in a summary log."""
        if not log_path:
            return None
        try:
            with open(log_path, encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    match = cls._DIANN_VERSION_RE.search(line)
                    if match:
                        return match.group(1)
        except OSError:
            logger.debug("Could not read DIA-NN log: %s", log_path)
        return None
