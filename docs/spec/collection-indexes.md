# Collection Indexes

!!! warning "Design proposal — not implemented"

    This entire page describes proposed index formats and APIs. QPX does not
    currently implement `qpx.open_collection()`, `coll.index()`,
    `coll.build_index()`, or `coll.build_all_indexes()`. The schemas, SQL, and
    code examples below are design sketches, not supported package behavior.
    For current cross-dataset queries, use the
    [DatasetCollection API](collection.md#current-api).

An **index** is a materialized, pre-computed data structure that enables fast search across all datasets in a [collection](collection.md). Unlike the collection itself (which is programmatic), indexes are **persisted as Parquet files** because the alternative -- scanning the full collection for every query -- is impractical at scale.

## Why Indexes?

Consider a collection with dozens of datasets. Without an index, searching for a peptide means scanning every PSM file; finding a differentially expressed protein means opening every DE view. With an index, the query reads only a small, targeted partition:

```python
# Without index: scans ALL parquet files (~minutes)
coll.sql("SELECT * FROM psm WHERE sequence = 'PEPTIDEK'")

# With index: reads only the relevant partition (~milliseconds)
coll.index("peptide").search("PEPTIDEK")
```

An index trades **one-time build cost** (minutes to hours) for **instant queries** (milliseconds). Think of it as a library catalog: the catalog is small and fast to search, and it tells you which book (dataset) to open.

## Index Types

QPX defines a standard framework for collection indexes. Each index type materializes data from a specific QPX view into a Hive-partitioned Parquet structure optimized for its search pattern.

### Peptide Index

Maps peptide sequences to the datasets that contain them. Built from `psm.parquet` across all datasets.

| | |
|---|---|
| **Source view** | PSM |
| **Stored at** | `_index/peptide/` |
| **Partitioned by** | First two amino acids of the sequence (~400 partitions) |

```text
_index/peptide/
├── sequence_prefix=AA/data.parquet
├── sequence_prefix=AC/data.parquet
├── ...
├── sequence_prefix=PE/data.parquet    # "PEPTIDEK" lives here
└── sequence_prefix=YY/data.parquet
```

#### Schema

Each row represents one (peptide, peptidoform, dataset) combination:

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | string | Unmodified peptide amino acid sequence |
| `peptidoform` | string | Peptide with modifications in ProForma notation |
| `project_accession` | string | Dataset identifier (e.g., `PXD000561`) |
| `charge_states` | array[int16] | Distinct precursor charge states observed |
| `spectra_count` | int32 | Number of PSMs for this peptide in this dataset |
| `best_pep` | float64 | Best (lowest) posterior error probability across all PSMs |
| `protein_accessions` | array[string] | All proteins this peptide maps to in this dataset |

#### Querying

```python
coll = qpx.open_collection("/data/my_collection/")

# Search by exact sequence
results = coll.index("peptide").search("PEPTIDEK")

# Search by prefix (wildcard)
results = coll.index("peptide").search_prefix("PEPTID")

# Search by peptidoform (with modifications)
results = coll.index("peptide").search_peptidoform("PEP(Phospho)TIDEK")
```

#### Build Algorithm

```sql
SELECT
    sequence, peptidoform, project_accession,
    LIST(DISTINCT charge)                     AS charge_states,
    COUNT(*)                                  AS spectra_count,
    MIN(posterior_error_probability)           AS best_pep,
    LIST(DISTINCT UNNEST(protein_accessions))  AS protein_accessions,
    LEFT(sequence, 2)                         AS sequence_prefix
FROM psm
WHERE is_decoy = false
GROUP BY sequence, peptidoform, project_accession
```

### Protein Index

Maps protein accessions to the datasets that identified or quantified them. Built from `pg.parquet` across all datasets.

| | |
|---|---|
| **Source view** | Protein Group (PG) |
| **Stored at** | `_index/protein/` |
| **Partitioned by** | First two characters of the anchor protein accession (~300-400 partitions) |

```text
_index/protein/
├── accession_prefix=A0/data.parquet
├── accession_prefix=P0/data.parquet   # "P04637" (TP53) lives here
├── accession_prefix=Q9/data.parquet
└── ...
```

#### Schema

Each row represents one (protein, dataset) combination:

| Field | Type | Description |
|-------|------|-------------|
| `anchor_protein` | string | Representative protein accession for the group |
| `protein_accessions` | array[string] | All accessions in the protein group |
| `project_accession` | string | Dataset identifier |
| `gg_names` | array[string] | Gene names associated with the protein group |
| `global_qvalue` | float64 | Protein-level global q-value |
| `num_peptides` | int32 | Number of distinct peptides supporting the protein group |
| `num_runs` | int32 | Number of MS runs where the protein was quantified |

#### Querying

```python
# Find all datasets that identified TP53
results = coll.index("protein").search("P04637")

# Search by gene name (requires scanning gene name column within partitions)
results = coll.index("protein").search_gene("TP53")

# Find all datasets with a protein family (prefix search)
results = coll.index("protein").search_prefix("P046")
```

#### Build Algorithm

```sql
SELECT
    anchor_protein, protein_accessions, project_accession,
    gg_names, global_qvalue,
    num_peptides,
    COUNT(DISTINCT run_file_name) AS num_runs,
    LEFT(anchor_protein, 2)       AS accession_prefix
FROM pg
GROUP BY anchor_protein, protein_accessions, project_accession,
         gg_names, global_qvalue, num_peptides
```

### Differential Expression Index

Maps proteins to their differential expression results across studies. Built from differential expression (DE) views. Enables meta-analysis queries like "which proteins are consistently up-regulated in cancer studies?"

| | |
|---|---|
| **Source view** | Differential Expression (DE) |
| **Stored at** | `_index/de/` |
| **Partitioned by** | First two characters of the protein accession |

```text
_index/de/
├── accession_prefix=A0/data.parquet
├── accession_prefix=P0/data.parquet
└── ...
```

#### Schema

Each row represents one (protein, contrast, dataset) combination:

| Field | Type | Description |
|-------|------|-------------|
| `protein_accession` | string | Protein accession |
| `gene_name` | string, null | Gene name |
| `project_accession` | string | Dataset identifier |
| `contrast` | string | Comparison label (e.g., "disease_vs_control") |
| `log2_fold_change` | float64 | Log2 fold change |
| `adj_pvalue` | float64 | Adjusted p-value (BH or similar) |
| `regulation` | string | "up", "down", or "ns" (not significant) |
| `organisms` | array[string] | Organisms in this study (from sample metadata) |

#### Querying

```python
# Which studies show TP53 as differentially expressed?
results = coll.index("de").search("P04637")
#   project_accession | contrast              | log2fc | adj_pvalue | regulation
#   PXD000561         | tumor_vs_normal       | 2.3    | 1.2e-08    | up
#   PXD002137         | treatment_vs_control  | -0.8   | 0.03       | down

# Meta-analysis: proteins consistently up-regulated across 3+ studies
results = coll.index("de").search_consistent(
    min_studies=3,
    direction="up",
    adj_pvalue_cutoff=0.05,
)

# Find all DE results for a gene
results = coll.index("de").search_gene("BRCA1")
```

#### Build Algorithm

```sql
SELECT
    protein_accession, gene_name, project_accession,
    contrast, log2_fold_change, adj_pvalue,
    CASE
        WHEN adj_pvalue < 0.05 AND log2_fold_change > 0 THEN 'up'
        WHEN adj_pvalue < 0.05 AND log2_fold_change < 0 THEN 'down'
        ELSE 'ns'
    END AS regulation,
    s.organisms,
    LEFT(protein_accession, 2) AS accession_prefix
FROM de
JOIN (
    SELECT project_accession, LIST(DISTINCT organism) AS organisms
    FROM sample
    GROUP BY project_accession
) s USING (project_accession)
```

### Future Index Types

The framework supports additional index types following the same pattern:

| Index | Source view | Partition key | Use case |
|-------|-----------|--------------|----------|
| **spectrum** | MZ | Binned precursor m/z | "Find spectra matching this m/z and charge" |
| **modification** | PSM/Feature | Modification name | "Which datasets have phosphorylation data?" |
| **sample** | Sample | Organism/tissue | "Find all liver tissue datasets" |

## Building Indexes

```python
coll = qpx.open_collection("/data/my_collection/")

# Build a specific index
coll.build_index("peptide")
coll.build_index("protein")
coll.build_index("de")

# Build all available indexes
coll.build_all_indexes()
```

The build process for each index:

1. Reads the source view from every dataset in the collection.
2. Applies index-specific filters (e.g., `is_decoy = false` for peptide index).
3. Aggregates per the index schema's grouping key.
4. Partitions by the index's partition column.
5. Writes Hive-partitioned Parquet files to `_index/{type}/`.

### Rebuilding

Indexes can go stale when datasets are added or removed. The rebuild process is idempotent -- it drops the existing index and recreates it from scratch:

```python
# Rebuild after adding new datasets
coll.build_index("peptide", force=True)

# Check if index is stale (compares dataset list vs index metadata)
coll.index("peptide").is_stale()
```

### Index Metadata

Each index stores a metadata file at `_index/{type}/_metadata.json`:

```json
{
    "index_type": "peptide",
    "source_view": "psm",
    "qpx_version": "1.1",
    "built_at": "2026-04-08T14:30:00Z",
    "datasets_included": ["PXD000561", "PXD002137"],
    "total_entries": 14400000,
    "partitions": 400,
    "compression": "zstd"
}
```

The `datasets_included` list enables staleness detection: if the current collection has datasets not in this list, the index is stale.

## Storage

Indexes are lightweight compared to the source data:

| Index type | Typical size | Notes |
|-----------|-------------|-------|
| Peptide | ~0.05% of PSM data | Aggregates millions of PSMs into unique peptide-dataset pairs |
| Protein | ~0.01% of PG data | One row per protein group per dataset |
| DE | ~0.02% of DE data | One row per protein per contrast per dataset |

## Conventions

- Indexes live under `_index/` at the collection root.
- The `_` prefix signals "infrastructure, not a dataset" -- collection discovery skips it.
- Each index type gets its own subdirectory: `_index/peptide/`, `_index/protein/`, `_index/de/`.
- Partition column names use the `{field}_prefix` naming pattern.
- Index Parquet files use ZSTD compression.
- Index metadata is stored as `_metadata.json` (not Parquet) for easy inspection.
- Each index declares its `source_view` in metadata, linking it to the QPX view it materializes.
