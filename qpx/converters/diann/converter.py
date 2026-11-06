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
from qpx.converters.sdrf import convert_sdrf_views
from qpx.core.constants import FEATURE, MZ, ONTOLOGY, PG, RUN, SAMPLE, XIC
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
        self._memory = max_memory
        self._threads = max_cpus
        self._compression = compression
        self._diann_version = self._parse_diann_version(diann_log)
        self._ontology_entries: list[dict] = []
        self._resolved_mappings_by_view: dict[str, dict] = {}
        # Views written by this run, in write order; provenance lists only these.
        self._written_views: list[str] = []

    def convert_features(
        self,
        mzml_info_folder=None,
        qvalue_threshold=None,
        output_folder=".",
        output_prefix=None,
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
        self._written_views.append(FEATURE)
        logger.info("DIA-NN feature conversion complete")

    def convert_pg(
        self,
        pg_matrix_path,
        output_folder=".",
        output_prefix=None,
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
        self._written_views.append(PG)
        logger.info("DIA-NN PG conversion complete")

    def convert_mz(
        self,
        mzml_statistics_folder,
        output_folder=".",
        output_prefix=None,
        matched_scans=None,
        feature_parquet=None,
    ):
        """Generate the ``mz`` view from mzML-statistics MS2 spectra.

        Writes only the feature-matched MS2 scans (the PSM fragment spectra).
        ``matched_scans`` (a set of ``(run_file_name, scan)`` tuples) may be
        supplied directly; otherwise it is read from ``feature_parquet``, the
        generated feature view (default ``{prefix}.feature.parquet`` in
        ``output_folder``), whose ``scan`` column links each feature to its MS2
        scan. The view is skipped when that feature view cannot be read.
        Peak arrays are written to the ``mz`` / ``intensity`` float32 list columns.
        """
        from qpx.converters.diann.mz_adapter import DiannMzAdapter

        output_folder = Path(output_folder)
        prefix = output_prefix or "diann"
        if not mzml_statistics_folder or not list(Path(mzml_statistics_folder).glob("*_ms2_info.parquet")):
            logger.warning("no *_ms2_info.parquet in %s; skipping m/z view", mzml_statistics_folder)
            return None
        if matched_scans is None:
            matched_scans = self._extract_matched_scans_from_feature(
                feature_parquet or output_folder / f"{prefix}.feature.parquet"
            )
            if matched_scans is None:
                return None
        adapter = DiannMzAdapter(compression=self._compression)
        written = adapter.convert(
            mzml_statistics_folder=str(mzml_statistics_folder),
            output_path=str(output_folder / f"{prefix}.mz.parquet"),
            matched_scans=matched_scans,
        )
        # The writer creates no file when no spectrum matched.
        if written:
            self._written_views.append(MZ)
        logger.info("DIA-NN mz conversion complete")
        return written

    def convert_xic(self, xic_paths, output_folder=".", output_prefix=None, **kwargs):
        """Convert DIA-NN XIC exports and record the written view for provenance."""
        from qpx.converters.diann.xic_adapter import convert_diann_xic

        prefix = output_prefix or "diann"
        traces = convert_diann_xic(
            xic_paths,
            Path(output_folder) / f"{prefix}.xic.parquet",
            duckdb_memory=self._memory,
            duckdb_threads=self._threads,
            **kwargs,
        )
        self._written_views.append(XIC)
        return traces

    def _extract_matched_scans_from_feature(self, feature_parquet) -> set[tuple[str, int]] | None:
        """Return the ``(run_file_name, scan)`` pairs from a generated feature view.

        The feature view's ``scan`` column (filled by the DIA-NN feature adapter,
        e.g. via RT matching in ``_merge_scan_info``) is the authoritative
        feature-to-MS2-scan link. Returns ``None`` if the columns are absent.
        """
        import pyarrow as pa
        import pyarrow.compute as pc
        import pyarrow.parquet as pq

        try:
            table = pq.read_table(str(feature_parquet), columns=["run_file_name", "scan"])
        except (FileNotFoundError, KeyError, pa.ArrowInvalid, pa.ArrowIOError):
            logger.warning("could not read run_file_name/scan from %s; skipping m/z view", feature_parquet)
            return None
        scans = table.column("scan")
        pairs = (
            pa.table(
                {
                    "run": table.column("run_file_name").take(pc.call_function("list_parent_indices", [scans])),
                    "scan": pc.call_function("list_flatten", [scans]),
                }
            )
            .filter(pc.field("scan").is_valid())
            .group_by(["run", "scan"])
            .aggregate([])
        )
        matched = set(zip(pairs.column("run").to_pylist(), pairs.column("scan").to_pylist()))
        logger.info("matched %s (run_file_name, scan) from feature view", len(matched))
        return matched

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
        if convert_sdrf_views(self.sdrf_path, output_folder, prefix, self._ontology_entries):
            self._written_views.extend([SAMPLE, RUN])
            logger.info("SDRF conversion complete (sample + run)")

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
        if self._write_ontology(Path(output_folder), prefix, entries) is not None:
            self._written_views.append(ONTOLOGY)

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
                "output_views": [view for view in self._written_views if view in (FEATURE, PG)],
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
                "output_views": list(self._written_views),
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
