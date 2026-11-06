"""MaxQuant orchestrator — composes PSM, Feature, and PG adapters."""

import logging
from pathlib import Path

from qpx._version import __version__
from qpx.converters.channel_labels import read_sdrf_labels
from qpx.converters.mappings import get_tool_meta
from qpx.converters.maxquant.channels import reporter_channels
from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter
from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter
from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.converters.sdrf import convert_sdrf_views
from qpx.core.constants import FEATURE, ONTOLOGY, PG, PSM, RUN, SAMPLE
from qpx.core.scores import field_ontology_entries, score_ontology_entries

logger = logging.getLogger(__name__)

# Rows per batch read from each MaxQuant table.
_BATCH_SIZE = 100_000


def _views_with_input(msms_file, evidence_file, protein_groups_file) -> list[str]:
    """The views whose MaxQuant input file was given."""
    inputs = ((PSM, msms_file), (FEATURE, evidence_file), (PG, protein_groups_file))
    return [view for view, path in inputs if path]


def _reads_reporter_columns(structures, evidence_file, protein_groups_file) -> bool:
    """Whether a view that maps SDRF labels to reporter columns (feature, pg) will be written."""
    return bool({FEATURE, PG} & set(structures) & set(_views_with_input(None, evidence_file, protein_groups_file)))


class MaxQuantConverter(BaseOrchestrator):
    """Orchestrate full MaxQuant conversion to QPX format."""

    def __init__(self, max_memory=None, max_cpus=None, compression: str = "zstd"):
        self._memory = max_memory
        self._cpus = max_cpus
        self._compression = compression
        self._resolved_mappings_by_view: dict[str, dict] = {}

    def convert(
        self,
        output_folder,
        msms_file=None,
        evidence_file=None,
        protein_groups_file=None,
        sdrf_file=None,
        output_prefix=None,
        structures=None,
        spectral_data=False,
        project_accession=None,
    ):
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        prefix = output_prefix or "maxquant"

        if structures is None:
            structures = _views_with_input(msms_file, evidence_file, protein_groups_file)

        ontology_entries: list[dict] = []

        check_labels = _reads_reporter_columns(structures, evidence_file, protein_groups_file)
        self._write_sdrf_views(sdrf_file, output_folder, prefix, ontology_entries, check_labels)

        if PSM in structures and msms_file:
            self._convert_psms(msms_file, sdrf_file, output_folder / f"{prefix}.psm.parquet", spectral_data, ontology_entries)
            logger.info("MaxQuant PSM conversion complete")

        if FEATURE in structures and evidence_file:
            self._convert_features(
                evidence_file, protein_groups_file, sdrf_file, output_folder / f"{prefix}.feature.parquet", ontology_entries
            )
            logger.info("MaxQuant feature conversion complete")

        if PG in structures and protein_groups_file:
            self._convert_pg(
                protein_groups_file, evidence_file, sdrf_file, output_folder / f"{prefix}.pg.parquet", ontology_entries
            )
            logger.info("MaxQuant PG conversion complete")

        self._write_conversion_metadata(output_folder, prefix, structures, project_accession, ontology_entries)

    @staticmethod
    def _write_sdrf_views(sdrf_file, output_folder: Path, prefix: str, ontology_entries: list[dict], check_labels: bool) -> None:
        """Write sample/run metadata from the SDRF (if it has the required columns).

        With ``check_labels`` (a feature or pg view will be written), labels the
        adapters cannot map to reporter columns (SILAC, labels of no known
        plex) fail here, before any view is written. The PSM view reads no
        reporter columns, so a PSM-only conversion does not check them.
        """
        if not sdrf_file:
            return
        if check_labels:
            reporter_channels(read_sdrf_labels(str(sdrf_file)))
        if convert_sdrf_views(str(sdrf_file), output_folder, prefix, ontology_entries):
            logger.info("SDRF conversion complete")

    def _write_provenance(self, output_folder: Path, prefix: str, structures: list[str]) -> None:
        """Write provenance.parquet with MaxQuant + QPX conversion steps."""
        records = [
            {
                "step_order": 1,
                "step_category": "database_search",
                "step_name": "spectrum_identification",
                "tool_name": "MaxQuant",
                "tool_version": None,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": structures,
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "maxquant_to_qpx",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": structures + [SAMPLE, RUN, ONTOLOGY],
            },
        ]
        super()._write_provenance(output_folder, prefix, records)

    def _convert_psms(self, msms_file, sdrf_file, output_path, spectral_data, ontology_entries):
        """Write PSMs and collect their ontology and source mappings."""
        with MaxQuantPsmAdapter(
            duckdb_memory=self._memory,
            duckdb_threads=self._cpus,
            compression=self._compression,
        ) as adapter:
            adapter.convert(
                msms_path=str(msms_file),
                output_path=str(output_path),
                chunksize=_BATCH_SIZE,
                spectral_data=spectral_data,
                sdrf_path=str(sdrf_file) if sdrf_file else None,
            )
            ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=PSM))
            self._resolved_mappings_by_view[PSM] = adapter.get_resolved_columns()

    def _convert_features(self, evidence_file, protein_groups_file, sdrf_file, output_path, ontology_entries):
        """Write features and collect their ontology and source mappings."""
        with MaxQuantFeatureAdapter(
            duckdb_memory=self._memory,
            duckdb_threads=self._cpus,
            compression=self._compression,
        ) as adapter:
            adapter.convert(
                evidence_path=str(evidence_file),
                output_path=str(output_path),
                sdrf_path=str(sdrf_file) if sdrf_file else None,
                protein_groups_path=(str(protein_groups_file) if protein_groups_file else None),
                chunksize=_BATCH_SIZE,
            )
            ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=FEATURE))
            self._resolved_mappings_by_view[FEATURE] = adapter.get_resolved_columns()

    def _convert_pg(self, protein_groups_file, evidence_file, sdrf_file, output_path, ontology_entries):
        """Write protein groups and collect their source mappings."""
        with MaxQuantPgAdapter(
            duckdb_memory=self._memory,
            duckdb_threads=self._cpus,
            compression=self._compression,
        ) as adapter:
            adapter.convert(
                protein_groups_path=str(protein_groups_file),
                output_path=str(output_path),
                sdrf_path=str(sdrf_file) if sdrf_file else None,
                evidence_path=str(evidence_file) if evidence_file else None,
                chunksize=_BATCH_SIZE,
            )
            ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=PG))
            self._resolved_mappings_by_view[PG] = adapter.get_resolved_columns()

    def _write_conversion_metadata(self, output_folder, prefix, structures, project_accession, ontology_entries):
        """Write field ontology, provenance, and dataset metadata."""
        for view_name, mappings in self._resolved_mappings_by_view.items():
            ontology_entries.extend(
                field_ontology_entries(
                    view=view_name,
                    resolved_mappings=mappings,
                    tool_name=get_tool_meta("maxquant")["tool_name"],
                )
            )

        self._write_ontology(output_folder, prefix, ontology_entries)
        self._write_provenance(output_folder, prefix, structures)
        self._write_dataset(
            output_folder,
            prefix,
            project_accession,
            software_name="MaxQuant",
            software_version=None,
        )
