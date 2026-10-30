# API Views

API views are **not** stored as standalone Parquet files. They are derived, on-demand views that the API computes at query time by joining and aggregating the primary data views (`psm.parquet`, `feature.parquet`, `pg.parquet`) with the metadata views (`run.parquet`, `sample.parquet`).

!!! info "Programmable, not prescriptive"
    Unlike the primary data views which have fixed YAML schemas, API views are **programmable**. The API can expose different fields, aggregations, and summaries depending on the use case. There is no rigid response schema -- what you retrieve is determined by your query.

## Concept

The QPX primary views (PSM, Feature, Protein Group) store the full, detailed data. API views sit on top of these and provide **convenience summaries** tailored to specific needs:

- A **peptide-level summary** might aggregate features across runs into one abundance value per peptide per sample.
- A **protein-level summary** might roll up protein group data into protein counts, abundances, and identification metrics per sample.
- A **custom view** might combine fields from multiple primary views to answer a specific question.

The key insight is that the primary data model is rich enough to derive many different summaries. Rather than defining a fixed schema for each possible summary, the API allows users to request exactly the data they need.

## Examples of what you can retrieve

Because views are programmable, different users can request different slices of the data:

| Use case | What you might request |
|----------|----------------------|
| Quick protein report | Protein accessions, gene names, number of peptides per protein |
| Abundance comparison | Protein accessions, abundance per sample, global q-values |
| Peptide coverage | Peptide sequences, modifications, which proteins they map to |
| Identification summary | Number of proteins, number of peptides, number of PSMs per sample |
| Sample overview | Per-sample counts and total intensities |

## How views are derived

API views are computed by joining and aggregating the primary QPX data. For example:

**Protein-level data** can be derived from the PG view by:

1. Joining `pg.parquet` with `run.parquet` to resolve sample-channel mappings.
2. Selecting the `anchor_protein` as the representative for each protein group.
3. Reading the per-label intensity directly from the flattened scalar `pg.intensity` (since QPX 1.1 the PG view has one row per label — there is no `intensities` list).
4. Optionally aggregating peptide/PSM counts, scores, or gene annotations.

```sql
-- Example: deriving protein abundance per sample from pg + run
SELECT DISTINCT
       pg.anchor_protein AS protein_accession,
       rs.sample_accession,
       pg.label,
       pg.intensity AS abundance,
       pg.global_qvalue,
       pg.gg_names
FROM 'PXD014414.pg.parquet' pg
CROSS JOIN UNNEST(pg.grouped_runs) AS grouped(run_file_name)
JOIN 'PXD014414.run.parquet' r USING (run_file_name)
CROSS JOIN UNNEST(r.samples) AS samples(rs)
WHERE pg.label = rs.label;   -- flattened: scalar pg.label/pg.intensity, no intensities list
```

`DISTINCT` is required because several member runs can be fractions of the
same sample. The protein-group quantity is already aggregated across those
fractions and must be returned only once per sample and label.

**Peptide-level data** can be derived from the Feature view by:

1. Joining `feature.parquet` with `run.parquet` to resolve sample-channel mappings.
2. Aggregating features across files and channels into a single abundance per sample.
3. Optionally joining with `pg.parquet` for gene annotations.

```sql
-- Example: deriving peptide abundance per sample from feature + run
SELECT f.sequence,
       f.peptidoform,
       rs.sample_accession,
       i.intensity AS abundance,
       pg.gg_names
FROM 'PXD014414.feature.parquet' f,
     'PXD014414.run.parquet' r,
     UNNEST(r.samples) AS rs,
     UNNEST(f.intensities) AS i,
     'PXD014414.pg.parquet' pg
WHERE f.run_file_name = r.run_file_name
  AND i.label = rs.label
  AND f.anchor_protein = pg.anchor_protein
  AND list_contains(pg.grouped_runs, f.run_file_name);
```

## Relationship to primary views

```mermaid
graph LR
    FEAT[feature.parquet] --> API["API Views<br/>(on demand)"]
    PG[pg.parquet] --> API
    RUN[run.parquet] --> API
    SAMPLE[sample.parquet] --> API
    API --> PEP[Peptide summaries]
    API --> PROT[Protein summaries]
    API --> CUSTOM[Custom queries]

    style API fill:#e8f5e9,stroke-dasharray: 5 5
    style PEP fill:#f5f5f5,stroke-dasharray: 5 5
    style PROT fill:#f5f5f5,stroke-dasharray: 5 5
    style CUSTOM fill:#f5f5f5,stroke-dasharray: 5 5
```

!!! tip "Use primary views for full detail"
    API views are for convenience and fast lookups. For full detail -- per-file intensities, scan references, protein group inference context, spectral data -- always use the primary views directly: [PSM](psm.md), [Feature](feature.md), [Protein Group](pg.md), [Mass Spectra](mz.md).

## Built-in API Views

The Python API provides several pre-built views accessible as properties on the `Dataset` object:

| View | Property | Data Method | Description |
|------|----------|-------------|-------------|
| Identification Summary | `ds.identification_summary` | `.summary()` | Per-run protein, peptide, and feature counts |
| Run Summary | `ds.run_summary` | `.summary()` | Per-run statistics (cached) |
| Modification View | `ds.modification_view` | `.frequency()` | Modification frequency across PSMs |
| QC View | `ds.qc_view` | `.metrics()` | Dataset-wide quality control metrics |
| Protein View | `ds.protein_view` | `.intensity()` | Protein abundance per sample (joins PG + Run) |
| Peptide View | `ds.peptide_view` | `.intensity()` | Peptide abundance per sample (joins Feature + Run) |

### Plotting

Each summary view has a `.plot()` method that returns a `matplotlib.Figure` for quick visualization. Matplotlib is imported lazily -- it is only required when `.plot()` is called.

```python
import qpx

ds = qpx.open_dataset("PXD014414/")

# Identification summary bar chart
fig = ds.identification_summary.plot()
fig.savefig("identifications.png")

# Run summary grouped bar chart
fig = ds.run_summary.plot(figsize=(14, 7))
fig.savefig("run_summary.png")

# QC metrics overview
fig = ds.qc_view.plot()
fig.savefig("qc_overview.png")

# Modification frequency (top 20 by PSM count)
fig = ds.modification_view.plot(top_n=20)
fig.savefig("modifications.png")
```

!!! tip "Customization"
    Each `.plot()` method accepts a `figsize` tuple. The `ModificationView.plot()` also accepts `top_n` to control how many modifications are shown. For more advanced customization, retrieve the data with `.summary().to_df()` and build your own plots.

### Abundance matrices

`ds.intensity(level)` returns the long form of the target (non-decoy) abundances, one row per sample, feature and [quantification unit](data-model.md#quantification-units-fractions-one-grouped_runs-one-quantity-per-label), with the columns `sample_accession`, `feature_id`, `intensity` and `quantification_unit`. At protein level a unit is a pg `grouped_runs` set and label, and `feature_id` is the `anchor_protein`. At peptide level the features of a peptide are summed within a unit: the charge states, and the runs that are fractions of one sample (runs with the same `run.samples` entries). The unit is named by its runs, sorted and joined by `;`, then `|` and the channel label: `F1;F2;F3|LFQ`, `plex1|TMT126`.

`ds.design_matrix(level)` pivots the long form into samples × features. A sample measured in more than one unit, such as a technical replicate or a reference channel present in every plex, has several values per feature, and so, at protein level, does a unit in which distinct groups share an `anchor_protein`. The matrix does not add them up: the default raises `ValueError`, `agg="mean"` (or `"median"`, `"sum"`, `"min"`, `"max"`) combines them, and `by="unit"` keeps one row per `(sample_accession, quantification_unit)`. The CDAP and Spectronaut converters write one pg row per raw file (`grouped_runs` names only that run), so at protein level each fraction of a sample is a unit of its own that holds part of the sample's quantity: use `agg="sum"` there (peptide units already join the fractions). Missing cells stay `NaN` unless `fillna` is given. `output_path` writes the same matrix to Parquet with DuckDB, out of core.

```python
ds = qpx.open_dataset("PXD001819/")   # 9 samples x 3 technical replicates
ds.design_matrix(agg="mean")          # 9 rows: replicate means
ds.design_matrix(by="unit")           # 27 rows: one per replicate
```

## Dataset Collections

For multi-dataset analysis, `DatasetCollection` registers multiple datasets in a single DuckDB engine. Each structure is indexed by dataset number (e.g., `feature_0`, `feature_1`).

### Virtual Mode (Cross-Dataset SQL)

```python
import qpx

ds1 = qpx.open_dataset("PXD014414/")
ds2 = qpx.open_dataset("PXD016999/")

coll = qpx.DatasetCollection([ds1, ds2])

# Query across datasets using indexed table names
result = coll.sql("""
    SELECT 'PXD014414' AS dataset, COUNT(*) AS n_features FROM feature_0
    UNION ALL
    SELECT 'PXD016999' AS dataset, COUNT(*) AS n_features FROM feature_1
""")
print(result.to_df())

# See what structures are available per dataset
print(coll.structure_names)
# {0: ['feature', 'pg', 'sample', 'run', 'dataset'], 1: ['feature', 'pg', ...]}
```

### Physical Merge

Merge matching structures from multiple datasets into a new directory. A `source_dataset` column is added to distinguish origins.

```python
coll.merge("merged_output/", structures=["feature", "pg"])

# Open the merged dataset
merged = qpx.open_dataset("merged_output/")
df = merged.feature.to_df()
print(df["source_dataset"].value_counts())
```

!!! note "Common structures only"
    When no `structures` argument is passed to `merge()`, only structures present in **all** datasets are merged.
