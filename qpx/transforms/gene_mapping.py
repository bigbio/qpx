"""Gene mapping transform — FASTA protein-to-gene annotation.

This transform maps protein accessions to gene names using the ``GN=`` field of
UniProt FASTA headers. The FASTA supplied by the caller is the only source: no
network service is queried, and no optional dependency is required.

The FASTA header format expected (UniProt):
    >sp|P12345|PROT_HUMAN Some description GN=BRCA1 PE=1 SV=2

The 'GN=' field is extracted as the gene name (symbol).

Usage:
    mapping = GeneMappingTransform(fasta_path="proteins.fasta")

    # Annotate a Feature or PG structure
    annotated_df = mapping.annotate_dataframe(feature_df)

    # Annotate via Dataset (returns new DataFrame, does not modify originals)
    annotated_df = mapping.annotate_dataset_features(dataset)

    # Write annotated Parquet
    mapping.write_annotated_features(dataset, "output.feature.parquet")
"""

from __future__ import annotations

import gzip
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)

_GENE_NAME_RE = re.compile(r"\bGN=(\S+)")

# "both" keys the map by accession AND entry name; it was previously reachable
# only by passing an unrecognised value, which made a typo silently change mode.
_MAP_BY_CHOICES = frozenset({"accession", "name", "both"})


def _parse_fasta_header(header: str) -> tuple[str, str, Optional[str]]:
    """Split one FASTA header line into (accession, entry name, gene name).

    Handles the UniProt ``>db|ACCESSION|NAME description`` form, the two-field
    ``>db|ACCESSION`` form, and a bare ``>IDENTIFIER``. The gene name is the
    ``GN=`` field when present, else None.
    """
    header = header[1:] if header.startswith(">") else header
    header = header.strip()
    identifier = header.split(None, 1)[0] if header else ""
    parts = identifier.split("|")

    if len(parts) >= 3:
        accession, name = parts[1], parts[2]
    elif len(parts) == 2:
        accession = name = parts[1]
    else:
        accession = name = identifier

    match = _GENE_NAME_RE.search(header)
    return accession, name, match.group(1) if match else None


def _parse_gene_names_from_fasta(
    fasta_path: str,
    map_by: str = "accession",
) -> dict[str, set[str]]:
    """
    Parse FASTA file and build a protein -> gene name mapping.

    Extracts gene names from the 'GN=' field in UniProt FASTA headers. Only
    header lines are read (sequences are skipped), and ``.gz`` files are opened
    transparently.

    Args:
        fasta_path: Path to the FASTA file (optionally gzip-compressed).
        map_by: How to key the mapping:
            - "accession": use UniProt accession (e.g., P12345)
            - "name": use UniProt entry name (e.g., PROT_HUMAN)

    Returns:
        Dict mapping protein identifier to set of gene names.
    """
    gene_map: dict[str, set[str]] = defaultdict(set)

    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    with opener(fasta_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            accession, name, gene_name = _parse_fasta_header(line)

            if map_by == "accession":
                keys = (accession,)
            elif map_by == "name":
                keys = (name,)
            else:
                # Map by both accession and name for maximum matching
                keys = (accession, name)
            for key in keys:
                gene_map[key].add(gene_name)

    # Count identifiers that actually carry a gene, not every header parsed. The
    # old count included GN-less entries, so a FASTA yielding no genes at all
    # still reported "Parsed gene names for N" — the opposite of the truth, on
    # the one line a caller would check before overwriting annotations.
    with_gene = sum(1 for names in gene_map.values() if any(n is not None for n in names))
    logger.info(
        "Parsed gene names for %d/%d protein identifiers from %s",
        with_gene,
        len(gene_map),
        fasta_path,
    )
    return gene_map


def _normalize_protein_list(accessions) -> Optional[list]:
    """Normalize a protein accessions value into a list suitable for _resolve_gene_names.

    - None / float NaN  → None
    - str               → [str]
    - list / array      → list(item)
    """
    if accessions is None or (isinstance(accessions, float) and accessions != accessions):
        return None
    if isinstance(accessions, str):
        return [accessions]
    return list(accessions)


def _extract_accession_from_raw(raw) -> Optional[str]:
    """Extract a bare accession string from a raw protein entry (str or struct dict)."""
    if raw is None or isinstance(raw, str):
        return raw
    return raw.get("accession")


def _resolve_gene_names(
    protein_accessions: list[str],
    gene_map: dict[str, set[str]],
) -> Optional[list[str]]:
    """
    Resolve gene names for a list of protein accessions.

    Args:
        protein_accessions: List of protein accessions from pg_accessions.
        gene_map: Protein-to-gene mapping from _parse_gene_names_from_fasta.

    Returns:
        List of gene names, or None if no mappings found.
    """
    if protein_accessions is None:
        return None

    gene_names = []
    for raw in protein_accessions:
        accession = _extract_accession_from_raw(raw)
        if accession is None:
            continue
        # Try full accession first, then extract short UniProt ID (sp|P12345|NAME → P12345)
        if accession not in gene_map:
            parts = accession.split("|")
            if len(parts) >= 2:
                accession = parts[1]
        if accession in gene_map:
            names = gene_map[accession]
            for name in names:
                if name is not None and name not in gene_names:
                    gene_names.append(name)

    return gene_names if gene_names else None


def _view_batches(file_paths):
    """Tables of bounded size read from a view's Parquet files, in order, with struct lists read back as lists."""
    from qpx.core.engine import create_connection
    from qpx.core.struct_lists import file_layouts, restore_table
    from qpx.transforms.protein_properties import _rewrite_batch_rows

    with create_connection() as con:
        for path in file_paths:
            parquet = pq.ParquetFile(path)
            layouts = file_layouts(parquet.metadata.metadata)
            for group in range(parquet.num_row_groups):
                for batch in parquet.iter_batches(batch_size=_rewrite_batch_rows(parquet, group), row_groups=[group]):
                    yield restore_table(con, pa.Table.from_batches([batch]), layouts)


def _annotate_table(table: pa.Table, gene_map: dict[str, set[str]], names_type: pa.DataType) -> tuple[pa.Table, int]:
    """*table* with gene names resolved from ``pg_accessions``, and the number of rows resolved.

    Existing ``gg_names`` stay where the FASTA resolves none.
    """
    # Rows share few protein groups: key each group by the dictionary codes of its
    # accessions, which no separator can occur in, resolve one row per distinct key
    # and take the names back to every row.
    groups = table.column("pg_accessions").combine_chunks()
    # feature lists accession structs, pg bare accessions.
    accessions = (
        pc.call_function("struct_field", [groups.values], pc.StructFieldOptions("accession"))
        if pa.types.is_struct(groups.type.value_type)
        else groups.values
    )
    codes = pc.fill_null(accessions.dictionary_encode().indices, -1)
    code_lists = type(groups).from_arrays(groups.offsets, codes, mask=groups.is_null())
    keys = pc.call_function("binary_join", [pc.cast(code_lists, pa.list_(pa.string())), pa.scalar(",")])
    encoded = pc.fill_null(keys, "").dictionary_encode()
    representatives = np.unique(encoded.indices.to_numpy(), return_index=True)[1]
    mapped = pa.array(
        [_resolve_gene_names(_normalize_protein_list(group), gene_map) for group in groups.take(representatives).to_pylist()],
        type=names_type,
    ).take(encoded.indices)
    names = mapped
    if "gg_names" in table.column_names:
        names = pc.call_function("if_else", [mapped.is_valid(), mapped, table.column("gg_names").cast(names_type)])
        table = table.drop_columns("gg_names")
    table = table.append_column("gg_names", names)
    return table, table.num_rows - mapped.null_count


class GeneMappingTransform:
    """
    Map protein accessions to gene names from a FASTA database.

    This transform reads UniProt FASTA headers to extract gene symbols (GN= field)
    and annotates QPX Feature or PG data structures with a gg_names column. The
    FASTA is the only source of truth, so annotation is offline and reproducible.

    Usage:
        mapping = GeneMappingTransform(fasta_path="proteins.fasta")

        # Get the raw gene map
        gene_map = mapping.gene_map

        # Annotate a DataFrame (adds a gg_names column)
        annotated_df = mapping.annotate_dataframe(df, protein_col="pg_accessions")

        # Annotate a Dataset's features
        feature_df = mapping.annotate_dataset_features(dataset)

        # Write annotated features to a new Parquet file
        mapping.write_annotated_features(dataset, "output.feature.parquet")
    """

    def __init__(
        self,
        fasta_path: Union[str, Path],
        map_by: str = "accession",
    ):
        """
        Initialize the gene mapping transform.

        Args:
            fasta_path: Path to the UniProt FASTA file (optionally gzip-compressed).
            map_by: Mapping strategy ("accession" or "name").
        """
        self._fasta_path = Path(fasta_path)
        if not self._fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        if map_by not in _MAP_BY_CHOICES:
            raise ValueError(f"map_by must be one of {sorted(_MAP_BY_CHOICES)}, got {map_by!r}")
        self._map_by = map_by

        # Lazily computed
        self._gene_map: Optional[dict[str, set[str]]] = None
        # Share of rows the last annotate_dataframe call resolved from the FASTA,
        # so callers can surface a zero-yield run instead of silently writing one.
        self._last_mapped_share: Optional[float] = None

    @property
    def last_mapped_share(self) -> Optional[float]:
        """Percentage of rows the most recent annotation resolved from the FASTA."""
        return self._last_mapped_share

    @property
    def gene_map(self) -> dict[str, set[str]]:
        """Protein-to-gene name mapping (lazy-loaded from FASTA)."""
        if self._gene_map is None:
            self._gene_map = _parse_gene_names_from_fasta(
                str(self._fasta_path),
                map_by=self._map_by,
            )
        return self._gene_map

    def annotate_dataframe(
        self,
        df: pd.DataFrame,
        protein_col: str = "pg_accessions",
    ) -> pd.DataFrame:
        """
        Add a gg_names column to a DataFrame.

        The protein_col should contain either:
        - A list of protein accessions (e.g., from QPX pg_accessions column)
        - A single protein accession string

        A FASTA carries no genomic accessions, so ``gg_accessions`` is left
        exactly as the converter wrote it (and created as NULL when absent).

        Args:
            df: Input DataFrame with protein identifiers.
            protein_col: Column name containing protein accessions.

        Returns:
            DataFrame with an added gg_names column.
        """
        if protein_col not in df.columns:
            raise ValueError(f"Column '{protein_col}' not found in DataFrame.")

        result = df.copy()
        gene_map = self.gene_map

        # Map protein accessions to gene names
        mapped = result[protein_col].apply(
            lambda accessions: _resolve_gene_names(
                _normalize_protein_list(accessions),
                gene_map,
            )
        )

        # A FASTA the accession is absent from says nothing about that protein's
        # gene — it is not evidence the gene is unknown. Overwriting regardless
        # wiped every gg_names a converter had already written (DIA-NN reports
        # carry them) whenever the FASTA lacked GN= fields, and --in-place made
        # that unrecoverable. Only rows that actually resolved are replaced.
        if "gg_names" in result.columns:
            existing = result["gg_names"]
            result["gg_names"] = [
                new_names if new_names is not None else old_names for new_names, old_names in zip(mapped, existing, strict=True)
            ]
        else:
            result["gg_names"] = mapped

        if "gg_accessions" not in result.columns:
            result["gg_accessions"] = None

        self._record_mapped_share(int(mapped.notna().sum()), len(result))
        return result

    def _record_mapped_share(self, n_mapped: int, n_rows: int) -> None:
        """Log and keep the share of rows resolved from the FASTA."""
        mapped_share = n_mapped / n_rows * 100 if n_rows else 0.0
        logger.info(
            "Mapped gene names from FASTA for %d/%d rows (%.1f%%)",
            n_mapped,
            n_rows,
            mapped_share,
        )
        self._last_mapped_share = mapped_share

    def _write_annotated_view(self, structure, writer_class, output_path: Path) -> None:
        """Write *structure* with gene names to *output_path*, keeping its identity recipe.

        The view's Parquet files are read in bounded batches, so only one batch is
        in memory. As in :meth:`annotate_dataframe`, existing names stay where the
        FASTA resolves none.
        """
        gene_map = self.gene_map
        source_composite = structure.file_metadata.get("identity_composite")
        identity_composite = tuple(source_composite.split(",")) if source_composite else None
        n_mapped = n_rows = 0
        with writer_class(output_path, override_provided_ids=False, identity_composite=identity_composite) as writer:
            names_type = writer.arrow_schema.field("gg_names").type
            for table in _view_batches(structure.file_paths):
                table, mapped = _annotate_table(table, gene_map, names_type)
                writer.write_table(writer.align_table_to_schema(table))
                n_mapped += mapped
                n_rows += table.num_rows
            if not n_rows:
                # An empty view still gets its (empty) file.
                writer.write_table(writer.align_table_to_schema(pq.read_schema(structure.file_paths[0]).empty_table()))
        self._record_mapped_share(n_mapped, n_rows)

    def annotate_dataset_features(
        self,
        dataset,
    ) -> pd.DataFrame:
        """
        Annotate a Dataset's Feature data with gene names.

        Materializes the Feature data as a DataFrame, adds gene annotations,
        and returns the annotated DataFrame.

        Args:
            dataset: A qpx.Dataset with feature data.

        Returns:
            Annotated Feature DataFrame with gg_names.
        """
        if dataset.feature is None:
            raise ValueError("Dataset does not contain feature data.")

        feature_df = dataset.feature.to_df()
        return self.annotate_dataframe(
            feature_df,
            protein_col="pg_accessions",
        )

    def write_annotated_features(
        self,
        dataset,
        output_path: Union[str, Path],
    ) -> Path:
        """
        Write gene-annotated Feature data to a new Parquet file.

        Uses the FeatureWriter to produce a schema-validated output file.

        Args:
            dataset: A qpx.Dataset with feature data.
            output_path: Path for the output .feature.parquet file.

        Returns:
            Path to the written file.
        """
        from qpx.writers.feature import FeatureWriter

        if dataset.feature is None:
            raise ValueError("Dataset does not contain feature data.")
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        self._write_annotated_view(dataset.feature, FeatureWriter, output_path)

        logger.info("Wrote gene-annotated features to %s", output_path)
        return output_path

    def annotate_dataset_pg(
        self,
        dataset,
    ) -> pd.DataFrame:
        """
        Annotate a Dataset's protein-group data with gene names.

        Args:
            dataset: A qpx.Dataset with pg data.

        Returns:
            Annotated PG DataFrame with gg_names.
        """
        if dataset.pg is None:
            raise ValueError("Dataset does not contain protein group data.")

        pg_df = dataset.pg.to_df()
        return self.annotate_dataframe(
            pg_df,
            protein_col="pg_accessions",
        )

    def write_annotated_pg(
        self,
        dataset,
        output_path: Union[str, Path],
    ) -> Path:
        """
        Write gene-annotated protein-group data to a new Parquet file.

        Uses the PgWriter to produce a schema-validated output file, preserving
        the source file's identity recipe so pg_id values do not change.

        Args:
            dataset: A qpx.Dataset with pg data.
            output_path: Path for the output .pg.parquet file.

        Returns:
            Path to the written file.
        """
        from qpx.writers.pg import PgWriter

        if dataset.pg is None:
            raise ValueError("Dataset does not contain protein group data.")
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        self._write_annotated_view(dataset.pg, PgWriter, output_path)

        logger.info("Wrote gene-annotated protein groups to %s", output_path)
        return output_path
