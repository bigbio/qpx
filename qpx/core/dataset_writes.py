"""Dataset write-back: PRIDE metadata, QPX views and AnnData files."""

from __future__ import annotations

from pathlib import Path

import pyarrow.parquet as pq

from qpx._version import __version__
from qpx.version import QPX_SPEC_VERSION


class _DatasetWriteMixin:
    """Write files into a local dataset directory."""

    # --- PRIDE enrichment ---
    def enrich_from_pride(self, project_accession: str | None = None) -> dict:
        """Fetch PRIDE metadata and update dataset.parquet.

        Retrieves project title, description, and PubMed ID from the
        PRIDE REST API and writes them into the dataset metadata.

        Parameters
        ----------
        project_accession : str or None
            PRIDE/ProteomeXchange accession. If *None*, reads from
            the existing dataset.parquet.

        Returns
        -------
        dict
            The fetched PRIDE metadata dict.

        Raises
        ------
        ValueError
            If no accession is provided and none is found in dataset.parquet,
            or if the accession is not found in PRIDE.
        """
        self._require_local("enrich_from_pride")

        # Resolve accession
        if project_accession is None:
            if self.dataset_meta is None:
                raise ValueError("No project_accession provided and no dataset.parquet found.")
            meta_df = self.dataset_meta.to_df()
            project_accession = meta_df["project_accession"].iloc[0]
            if not project_accession or project_accession == "unknown":
                raise ValueError("No valid project_accession in dataset.parquet. Provide one explicitly.")

        from qpx.core.pride import fetch_pride_metadata

        metadata = fetch_pride_metadata(project_accession)

        # Read existing dataset record, update fields, and rewrite
        if self.dataset_meta is not None:
            meta_df = self.dataset_meta.to_df()
            record = meta_df.iloc[0].to_dict()
        else:
            record = {
                "project_accession": project_accession,
                "creation_date": None,
            }

        record["project_title"] = metadata["project_title"]
        record["project_description"] = metadata["project_description"]
        record["pubmed_id"] = metadata["pubmed_id"]

        self.save_structure([record], "dataset")
        self.refresh()

        return metadata

    # --- Write-back ---
    # Writer registry: name → (WriterClassName, file_suffix)
    _WRITER_REGISTRY = {
        "psm": ("PsmWriter", ".psm.parquet"),
        "feature": ("FeatureWriter", ".feature.parquet"),
        "pg": ("PgWriter", ".pg.parquet"),
        "mz": ("MzWriter", ".mz.parquet"),
        "sample": ("SampleWriter", ".sample.parquet"),
        "run": ("RunWriter", ".run.parquet"),
        "dataset": ("DatasetWriter", ".dataset.parquet"),
        "ontology": ("OntologyWriter", ".ontology.parquet"),
        "provenance": ("ProvenanceWriter", ".provenance.parquet"),
        "pepmap": ("PepMapWriter", ".pepmap.parquet"),
        "site": ("SiteWriter", ".site.parquet"),
        "proteoform": ("ProteoformWriter", ".proteoform.parquet"),
        "protein": ("ProteinWriter", ".protein.parquet"),
    }

    def _existing_files(self, structure: str) -> list[tuple[str, Path]]:
        """(structure, single file) of the local files of *structure*, then of the other structures."""
        names = [structure, *(name for name in self._structures if name != structure)]
        return [
            (name, Path(path))
            for name in names
            if name in self._structures
            for path in self._structures[name].file_paths[:1]
            if Path(path).is_file()
        ]

    def _dataset_prefix(self, structure: str) -> str:
        """The prefix of the dataset's files, so a write-back replaces a view rather than adding a shard beside it."""
        if self._file_prefix:
            return self._file_prefix
        for name, path in self._existing_files(structure):
            suffix = self._STRUCTURE_REGISTRY[name][1]
            if path.name.endswith(suffix):
                return path.name[: -len(suffix)]
        return Path(self.path).name

    def _source_writer_options(self, structure: str) -> dict:
        """Writer options that keep what the structure's footer declares, and its ids.

        Hashing other columns than the source's composite would change every id
        and leave ``psm.feature_id`` dangling.
        """
        options: dict = {"override_provided_ids": False}
        files = [path for name, path in self._existing_files(structure) if name == structure]
        metadata = pq.read_schema(files[0]).metadata or {} if files else {}
        for key, option in (
            (b"creator", "creator"),
            (b"software_provider", "software_provider"),
            (b"scan_format", "scan_format"),
        ):
            if key in metadata:
                options[option] = metadata[key].decode()
        if b"identity_composite" in metadata:
            options["identity_composite"] = metadata[b"identity_composite"].decode().split(",")
        return options

    @classmethod
    def writer_for(cls, structure: str) -> tuple[type, str]:
        """The writer class and file suffix of *structure*."""
        import qpx.writers as writers_mod

        writer_name, suffix = cls._WRITER_REGISTRY[structure]
        return getattr(writers_mod, writer_name), suffix

    def save_structure(
        self,
        data,
        structure: str,
        prefix: str | None = None,
    ) -> Path:
        """Write validated Parquet back into the dataset directory.

        Rows keep their ids; a row without one gets it from the identity composite
        the existing file declares, which the new footer declares too.

        Only supported for local datasets. For S3-backed datasets, use the
        appropriate writer (e.g. ``FeatureWriter(path)``) with an explicit
        local output path.

        Parameters
        ----------
        data : list[dict] | pd.DataFrame | pa.Table
            Records to write. Schema-validated by the QPX writer.
        structure : str
            Structure name (e.g., "feature", "pg", "sample").
        prefix : str | None
            File name prefix. Defaults to the prefix of the dataset's files, which
            replaces an existing view. Another prefix adds a file that is read
            together with it.

        Returns
        -------
        Path
            Path to the written Parquet file.
        """
        import pyarrow as pa

        self._require_local("save_structure")
        if structure not in self._WRITER_REGISTRY:
            raise ValueError(f"Unknown structure '{structure}'. Valid: {list(self._WRITER_REGISTRY.keys())}")

        writer_cls, suffix = self.writer_for(structure)

        path = Path(self.path)
        file_prefix = prefix or self._dataset_prefix(structure)
        output_path = path / f"{file_prefix}{suffix}"

        with writer_cls(output_path, **self._source_writer_options(structure)) as writer:
            if isinstance(data, list):
                writer.write_batch(data)
            elif isinstance(data, pa.Table):
                # Project onto the writer's schema (adds any absent columns such
                # as a derived id or optional cross-refs, then casts to normalize
                # nullability/types). write_table stamps the derived id.
                writer.write_table(writer.align_table_to_schema(data))
            else:
                writer.write_dataframe(data)

        return output_path

    # AnnData view type -> canonical file suffix (before .h5ad)
    _ANNDATA_VIEWS = {"ae", "de"}

    def save_anndata(
        self,
        adata,
        name: str | None = None,
        *,
        view: str | None = None,
    ) -> Path:
        """Write an AnnData object to the dataset directory as ``.h5ad``.

        Only supported for local datasets. For S3-backed datasets, write to
        a local path with ``adata.write_h5ad(path)``.

        Parameters
        ----------
        adata : anndata.AnnData
            The AnnData object to save.
        name : str, optional
            Base name for the file (without extension), e.g., "de_results".
            When omitted the dataset prefix is used.
        view : str, optional
            AnnData view type: ``"ae"`` (absolute expression) or ``"de"``
            (differential expression).  When provided the output file
            follows the QPX naming convention
            ``<prefix>.<view>.h5ad`` (e.g., ``PXD000000.ae.h5ad``).
            Ignored when *name* is explicitly given.

        Returns
        -------
        Path
            Path to the written .h5ad file.
        """
        self._require_local("save_anndata")
        path = Path(self.path)
        if name is not None:
            output_path = path / f"{name}.h5ad"
        elif view is not None:
            if view not in self._ANNDATA_VIEWS:
                raise ValueError(f"Unknown AnnData view '{view}'. Choose from: {sorted(self._ANNDATA_VIEWS)}")
            output_path = path / f"{self._dataset_prefix('dataset')}.{view}.h5ad"
        else:
            raise ValueError("Either 'name' or 'view' must be provided")
        adata.uns["qpx_version"] = QPX_SPEC_VERSION
        adata.uns["writer_version"] = __version__
        adata.write_h5ad(output_path)
        return output_path
