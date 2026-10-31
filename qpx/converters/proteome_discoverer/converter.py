"""Proteome Discoverer orchestrator: PSM, feature and pg views from the text exports."""

from __future__ import annotations

import logging
from pathlib import Path

from qpx._version import __version__
from qpx.converters.mappings import get_tool_meta
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.converters.proteome_discoverer.feature_adapter import ProteomeDiscovererFeatureAdapter
from qpx.converters.proteome_discoverer.pg_adapter import ProteomeDiscovererPgAdapter, master_groups, read_proteins
from qpx.converters.proteome_discoverer.psm_adapter import ProteomeDiscovererPsmAdapter
from qpx.converters.proteome_discoverer.reports import runs_by_file_id
from qpx.core.constants import FEATURE, ONTOLOGY, PEPMAP, PG, PSM, RUN, SAMPLE
from qpx.core.scores import field_ontology_entries, score_ontology_entries

logger = logging.getLogger(__name__)
_TOOL = "proteome_discoverer"


class ProteomeDiscovererConverter(BaseOrchestrator):
    """Orchestrate the conversion of Proteome Discoverer ``PSMs.txt`` and ``Proteins.txt``."""

    def __init__(self, output_directory, max_memory=None, max_cpus=None, compression: str = "zstd"):
        self._output_dir = Path(output_directory)
        self._adapter_options = {"duckdb_memory": max_memory, "duckdb_threads": max_cpus, "compression": compression}
        self._compression = compression
        self._resolved_mappings_by_view: dict[str, dict] = {}

    def convert(self, psm_file, proteins_file=None, *, input_files=None, sdrf_file=None, output_prefix=None, **options):
        """Write the views the exports support; *input_files* names the raw run of each ``File ID``."""
        out = self._output_dir
        out.mkdir(parents=True, exist_ok=True)
        prefix = output_prefix or _TOOL
        runs = runs_by_file_id(input_files) if input_files else None
        ontology_entries: list[dict] = []
        produced = [PSM]
        groups, genes = {}, {}
        proteins = read_proteins(str(proteins_file)) if proteins_file else None
        if proteins is not None:
            groups, genes = master_groups(proteins)

        with ProteomeDiscovererPsmAdapter(**self._adapter_options) as psm_adapter:
            psm_adapter.load_psms(str(psm_file), runs)
            psm_adapter.convert(str(out / f"{prefix}.psm.parquet"), rank1_only=options.get("rank1_only", False))
            ontology_entries.extend(score_ontology_entries(psm_adapter.get_discovered_scores(), view=PSM))
            self._resolved_mappings_by_view[PSM] = psm_adapter.get_resolved_columns()
            psm_runs = runs or psm_adapter.runs_by_file_id()

        with ProteomeDiscovererFeatureAdapter(**self._adapter_options) as feature_adapter:
            feature_adapter.load_psms(str(psm_file), runs)
            feature_adapter.load_protein_groups(groups, genes)
            if feature_adapter.convert(str(out / f"{prefix}.feature.parquet")):
                produced.append(FEATURE)
                ontology_entries.extend(score_ontology_entries(feature_adapter.get_discovered_scores(), view=FEATURE))
                if feature_adapter.write_pepmap(str(out / f"{prefix}.pepmap.parquet")):
                    produced.append(PEPMAP)

        if proteins is not None:
            with ProteomeDiscovererPgAdapter(**self._adapter_options) as pg_adapter:
                pg_adapter.convert(proteins, str(out / f"{prefix}.pg.parquet"), psm_runs)
                ontology_entries.extend(score_ontology_entries(pg_adapter.get_discovered_scores(), view=PG))
                self._resolved_mappings_by_view[PG] = pg_adapter.get_resolved_columns()
            produced.append(PG)

        if sdrf_file:
            ontology_entries.extend(self.convert_sdrf(sdrf_file, prefix))
        self._write_metadata(out, prefix, produced, ontology_entries, options.get("project_accession"))

    def convert_sdrf(self, sdrf_file, prefix: str) -> list[dict]:
        """Write the sample and run views of an SDRF; returns their ontology terms."""
        from qpx.converters.sdrf import SdrfConverter

        with SdrfConverter() as conv:
            conv.convert(
                sdrf_path=str(sdrf_file),
                sample_output=str(self._output_dir / f"{prefix}.sample.parquet"),
                run_output=str(self._output_dir / f"{prefix}.run.parquet"),
            )
            return conv.run_ontology_entries()

    def _write_metadata(self, out, prefix, produced, ontology_entries, project_accession):
        tool_name = get_tool_meta(_TOOL)["tool_name"]
        for view_name, mappings in self._resolved_mappings_by_view.items():
            ontology_entries.extend(field_ontology_entries(view=view_name, resolved_mappings=mappings, tool_name=tool_name))
        self._write_ontology(out, prefix, ontology_entries)
        self._write_provenance(out, prefix, self._provenance(produced))
        self._write_dataset(out, prefix, project_accession, software_name=tool_name, software_version=None)

    @staticmethod
    def _provenance(structures: list[str]) -> list[dict]:
        common = {"tool_version": None, "tool_uri": None, "parameters": None, "config": None}
        return [
            {
                **common,
                "step_order": 1,
                "step_category": "database_search",
                "step_name": "spectrum_identification",
                "tool_name": "Proteome Discoverer",
                "output_views": structures,
            },
            {
                **common,
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "proteome_discoverer_to_qpx",
                "tool_name": "qpx",
                "tool_version": __version__,
                "output_views": structures + [SAMPLE, RUN, ONTOLOGY],
            },
        ]
