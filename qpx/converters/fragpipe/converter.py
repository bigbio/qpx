"""FragPipe orchestrator — composes PSM, Feature, and PG adapters."""

import logging
from pathlib import Path

from qpx._version import __version__
from qpx.converters.fragpipe.pg_adapter import FragPipePgAdapter
from qpx.converters.fragpipe.psm_adapter import FragPipePsmAdapter
from qpx.converters.mappings import get_tool_meta
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.core.constants import FEATURE, ONTOLOGY, PEPMAP, PG, PSM, RUN, SAMPLE
from qpx.core.scores import field_ontology_entries, score_ontology_entries

logger = logging.getLogger(__name__)

# Rows per batch read from each FragPipe table.
_BATCH_SIZE = 1_000_000


class FragPipeConverter(BaseOrchestrator):
    """Orchestrate full FragPipe conversion to QPX format."""

    def __init__(self, output_directory=None, max_memory=None, max_cpus=None, compression: str = "zstd"):
        self._output_dir = Path(output_directory) if output_directory else None
        self._memory = max_memory
        self._cpus = max_cpus
        self._compression = compression
        self._resolved_mappings_by_view: dict[str, dict] = {}

    def convert(
        self,
        psm_file=None,
        ion_file=None,
        pg_file=None,
        sdrf_file=None,
        output_prefix=None,
        output_folder=None,
        project_accession=None,
        experiment_to_runs=None,
        experiment_annotation_file=None,
    ):
        out = Path(output_folder) if output_folder else self._output_dir
        if out is None:
            raise ValueError("output_folder or output_directory required")
        out.mkdir(parents=True, exist_ok=True)
        prefix = output_prefix or "fragpipe"

        ontology_entries: list[dict] = []
        produced_structures: list[str] = []

        if psm_file:
            self._convert_psms(psm_file, out / f"{prefix}.psm.parquet", ontology_entries)
            produced_structures.append(PSM)
            logger.info("FragPipe PSM conversion complete")

        if ion_file:
            produced_structures += self._convert_features(
                ion_file,
                (psm_file, pg_file, sdrf_file),
                (out, prefix),
                (experiment_to_runs, experiment_annotation_file),
                ontology_entries,
            )
            logger.info("FragPipe feature conversion complete")

        if pg_file:
            self._convert_pg(
                pg_file,
                sdrf_file,
                out / f"{prefix}.pg.parquet",
                (experiment_to_runs, experiment_annotation_file),
                ontology_entries,
            )
            produced_structures.append(PG)
            logger.info("FragPipe PG conversion complete")

        if sdrf_file:
            self._convert_sdrf(sdrf_file, out / f"{prefix}.sample.parquet", out / f"{prefix}.run.parquet", ontology_entries)

        self._write_conversion_metadata(out, prefix, produced_structures, project_accession, ontology_entries)

    def _build_provenance_records(self, structures: list[str]) -> list[dict]:
        """Build provenance records for FragPipe + QPX conversion steps."""
        return [
            {
                "step_order": 1,
                "step_category": "database_search",
                "step_name": "spectrum_identification",
                "tool_name": "FragPipe/MSFragger",
                "tool_version": None,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": structures,
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "fragpipe_to_qpx",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": structures + [SAMPLE, RUN, ONTOLOGY],
            },
        ]

    def _convert_psms(self, psm_file, output_path, ontology_entries):
        """Write PSMs and collect their ontology and source mappings."""
        with FragPipePsmAdapter(duckdb_memory=self._memory, duckdb_threads=self._cpus, compression=self._compression) as adapter:
            adapter.convert(
                psm_path=str(psm_file),
                output_path=str(output_path),
                chunksize=_BATCH_SIZE,
            )
            ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=PSM))
            self._resolved_mappings_by_view[PSM] = adapter.get_resolved_columns()

    def _convert_features(self, ion_file, sources, output, mapping, ontology_entries) -> list[str]:
        """Write features, and the pepmap of their mapped proteins, using the supplied experiment mapping.

        The protein groups of ``pg_file`` give each feature its full group
        membership. Returns the views written.
        """
        psm_file, pg_file, sdrf_file = sources
        out, prefix = output
        experiment_to_runs, experiment_annotation_file = mapping
        from qpx.converters.fragpipe.feature_adapter import FragPipeFeatureAdapter

        with FragPipeFeatureAdapter(
            duckdb_memory=self._memory, duckdb_threads=self._cpus, compression=self._compression
        ) as adapter:
            if pg_file:
                adapter.load_protein_groups(str(pg_file))
            adapter.convert(
                feature_path=str(ion_file),
                output_path=str(out / f"{prefix}.feature.parquet"),
                sdrf_path=str(sdrf_file) if sdrf_file else None,
                psm_path=str(psm_file) if psm_file else None,
                experiment_to_runs=experiment_to_runs,
                experiment_annotation_path=(str(experiment_annotation_file) if experiment_annotation_file else None),
                chunksize=_BATCH_SIZE,
            )
            ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=FEATURE))
            self._resolved_mappings_by_view[FEATURE] = adapter.get_resolved_columns()
            if adapter.write_pepmap(str(out / f"{prefix}.pepmap.parquet")):
                return [FEATURE, PEPMAP]
        return [FEATURE]

    def _convert_pg(self, pg_file, sdrf_file, output_path, mapping, ontology_entries):
        """Write protein groups using the supplied experiment mapping."""
        experiment_to_runs, experiment_annotation_file = mapping
        with FragPipePgAdapter(duckdb_memory=self._memory, duckdb_threads=self._cpus, compression=self._compression) as adapter:
            adapter.convert(
                protein_path=str(pg_file),
                output_path=str(output_path),
                chunksize=_BATCH_SIZE,
                experiment_to_runs=experiment_to_runs,
                experiment_annotation_path=(str(experiment_annotation_file) if experiment_annotation_file else None),
                sdrf_path=str(sdrf_file) if sdrf_file else None,
            )
            ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=PG))
            self._resolved_mappings_by_view[PG] = adapter.get_resolved_columns()

    def _convert_sdrf(self, sdrf_file, sample_path, run_path, ontology_entries):
        """Write SDRF views and collect run ontology terms."""
        from qpx.converters.sdrf import SdrfConverter

        with SdrfConverter() as conv:
            conv.convert(
                sdrf_path=str(sdrf_file),
                sample_output=str(sample_path),
                run_output=str(run_path),
            )
            ontology_entries.extend(conv.run_ontology_entries())
        logger.info("SDRF conversion complete")

    def _write_conversion_metadata(self, out, prefix, produced_structures, project_accession, ontology_entries):
        """Write field ontology, provenance, and dataset metadata."""
        for view_name, mappings in self._resolved_mappings_by_view.items():
            ontology_entries.extend(
                field_ontology_entries(
                    view=view_name,
                    resolved_mappings=mappings,
                    tool_name=get_tool_meta("fragpipe")["tool_name"],
                )
            )

        self._write_ontology(out, prefix, ontology_entries)
        self._write_provenance(out, prefix, self._build_provenance_records(produced_structures))
        self._write_dataset(
            out,
            prefix,
            project_accession,
            software_name="FragPipe/MSFragger",
            software_version=None,
        )
