"""Spectronaut orchestrator — composes Feature and PG adapters."""

from __future__ import annotations

import logging
from pathlib import Path

from qpx._version import __version__
from qpx.converters.base import resolve_columns
from qpx.converters.mappings import get_field_mappings, get_tool_meta
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.converters.sdrf import convert_sdrf_views
from qpx.converters.spectronaut.feature_adapter import SpectronautFeatureAdapter
from qpx.converters.spectronaut.pg_adapter import SpectronautPgAdapter
from qpx.converters.spectronaut.report_layout import long_columns, report_header, unavailable_views
from qpx.core.constants import FEATURE, ONTOLOGY, PG, RUN, SAMPLE
from qpx.core.scores import field_ontology_entries, score_ontology_entries

logger = logging.getLogger(__name__)


class SpectronautConverter(BaseOrchestrator):
    """Orchestrate full Spectronaut conversion to QPX format."""

    def __init__(
        self,
        report_path,
        sdrf_path=None,
        max_memory=None,
        max_cpus=None,
        compression: str = "zstd",
    ):
        """Initialize converter with report path and DuckDB settings."""
        self.report_path = str(report_path)
        self.sdrf_path = str(sdrf_path) if sdrf_path else None
        self._memory = max_memory
        self._threads = max_cpus
        self._compression = compression
        self._ontology_entries: list[dict] = []
        self._resolved_mappings_by_view: dict[str, dict] = {}

    def available_views(self) -> list[str]:
        """The views this report can produce, warning about the others; raises when it produces none."""
        unavailable = unavailable_views(long_columns(report_header(self.report_path)))
        for view, reason in unavailable.items():
            logger.warning("Skipping the %s view: %s", view, reason)
        views = [view for view in (FEATURE, PG) if view not in unavailable]
        if not views:
            raise ValueError(f"{self.report_path} is not a convertible Spectronaut report: {'; '.join(unavailable.values())}")
        return views

    def convert_features(
        self,
        qvalue_threshold=None,
        output_folder=".",
        output_prefix=None,
    ):
        """Convert Spectronaut features to feature.parquet."""
        output_folder = Path(output_folder)
        prefix = output_prefix or "spectronaut"
        with SpectronautFeatureAdapter(
            duckdb_memory=self._memory,
            duckdb_threads=self._threads,
            compression=self._compression,
        ) as adapter:
            adapter.convert(
                spectronaut_report=self.report_path,
                output_path=str(output_folder / f"{prefix}.feature.parquet"),
                sdrf_path=self.sdrf_path,
                qvalue_threshold=qvalue_threshold,
            )
            self._ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=FEATURE))
            cols = adapter.get_table_columns("report")
            self._resolved_mappings_by_view[FEATURE] = resolve_columns(get_field_mappings("spectronaut", "feature"), cols)
        logger.info("Spectronaut feature conversion complete")

    def convert_pg(
        self,
        output_folder=".",
        output_prefix=None,
    ):
        """Convert Spectronaut protein groups to pg.parquet."""
        output_folder = Path(output_folder)
        prefix = output_prefix or "spectronaut"
        with SpectronautPgAdapter(
            duckdb_memory=self._memory,
            duckdb_threads=self._threads,
            compression=self._compression,
        ) as adapter:
            adapter.convert(
                spectronaut_report=self.report_path,
                output_path=str(output_folder / f"{prefix}.pg.parquet"),
            )
            self._ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=PG))
            cols = adapter.get_table_columns("report")
            self._resolved_mappings_by_view[PG] = resolve_columns(get_field_mappings("spectronaut", "pg"), cols)
        logger.info("Spectronaut PG conversion complete")

    def convert_sdrf(
        self,
        output_folder: str | Path,
        prefix: str = "spectronaut",
    ) -> None:
        """Convert SDRF to sample.parquet and run.parquet."""
        output_folder = Path(output_folder)
        if not self.sdrf_path:
            logger.warning("No SDRF path provided — skipping sample/run conversion")
            return
        if convert_sdrf_views(self.sdrf_path, output_folder, prefix, self._ontology_entries):
            logger.info("SDRF conversion complete (sample + run)")

    def write_ontology(
        self,
        output_folder: str | Path,
        prefix: str = "spectronaut",
    ) -> None:
        """Write combined ontology.parquet."""
        entries = list(self._ontology_entries)
        tool_name = get_tool_meta("spectronaut")["tool_name"]
        for view_name, mappings in self._resolved_mappings_by_view.items():
            entries.extend(
                field_ontology_entries(
                    view=view_name,
                    resolved_mappings=mappings,
                    tool_name=tool_name,
                )
            )
        self._write_ontology(Path(output_folder), prefix, entries)

    def write_provenance(
        self,
        output_folder: str | Path,
        prefix: str = "spectronaut",
    ) -> None:
        """Write provenance.parquet."""
        records = [
            {
                "step_order": 1,
                "step_category": "quantification",
                "step_name": "precursor_quantification",
                "tool_name": "Spectronaut",
                "tool_version": None,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": list(self._resolved_mappings_by_view),
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "spectronaut_to_qpx",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": [
                    {
                        "key": "report_path",
                        "value": Path(self.report_path).name,
                    },
                ],
                "config": None,
                "output_views": [*self._resolved_mappings_by_view, SAMPLE, RUN, ONTOLOGY],
            },
        ]
        self._write_provenance(Path(output_folder), prefix, records)

    def write_dataset(
        self,
        output_folder: str | Path,
        prefix: str = "spectronaut",
        project_accession: str | None = None,
    ) -> None:
        """Write dataset.parquet."""
        self._write_dataset(
            Path(output_folder),
            prefix,
            project_accession,
            software_name="Spectronaut",
        )
