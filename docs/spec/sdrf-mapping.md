# SDRF Conversion

When SDRF is ingested into QPX, bracket-style column names (e.g. `characteristics[organism]`) are mapped to clean `snake_case` field names (e.g. `organism`) and the data is split into separate Parquet files -- one for biological samples and one for data acquisition runs. This page documents the column mapping registry, the multi-valued metadata convention, and the end-to-end conversion process.

---

## SDRF Column Name Mapping

When SDRF is converted to QPX, bracket-style column names are mapped to `snake_case` field names according to the following registry. This mapping is deterministic and applied automatically during ingestion.

### Sample Columns (→ sample.parquet)

Sample columns describe biological properties of the specimen. These are deduplicated and stored in [sample.parquet](sample.md).

| SDRF Column | QPX Field Name | Target File |
|---|---|---|
| `source name` | `sample_accession` | sample.parquet |
| `characteristics[organism]` | `organism` | sample.parquet |
| `characteristics[organism part]` | `organism_part` | sample.parquet |
| `characteristics[disease]` | `disease` | sample.parquet |
| `characteristics[cell line]` | `cell_line` | sample.parquet |
| `characteristics[cell type]` | `cell_type` | sample.parquet |
| `characteristics[sex]` | `sex` | sample.parquet |
| `characteristics[age]` | `age` | sample.parquet |
| `characteristics[developmental stage]` | `developmental_stage` | sample.parquet |
| `characteristics[ancestry category]` | `ancestry` | sample.parquet |
| `characteristics[individual]` | `individual` | sample.parquet |

### Run Columns (→ run.parquet)

Run columns describe technical and instrument-related properties of the data acquisition. These are stored in [run.parquet](run.md).

| SDRF Column | QPX Field Name | Target File |
|---|---|---|
| `comment[data file]` (stem) | `run_accession`, `run_file_name` | run.parquet |
| `comment[data file]` | `run_file_name` | run.parquet |
| `comment[label]` | label (inside `samples` list) | run.parquet |
| `comment[fraction identifier]` | `fraction` | run.parquet |
| `characteristics[biological replicate]` | `biological_replicate` (inside `samples` list) | run.parquet |
| `comment[technical replicate]` | `technical_replicate` (inside `samples` list) | run.parquet |
| `comment[instrument]` | `instrument` (plain string; the term is in `additional_terms`) | run.parquet |
| `comment[cleavage agent details]` | `enzymes` (list of plain strings; the terms are in `additional_terms`) | run.parquet |
| `comment[proteomics data acquisition method]` | `acquisition_method` (plain string; also resolved into `additional_terms` as ONTOLOGY_TERM) | run.parquet |
| `comment[dissociation method]` | `dissociation_method` (plain string; the term is in `additional_terms`) | run.parquet |
| `comment[modification parameters]` | `modification_parameters` (as list of MODIFICATION) | run.parquet |

### Factor Columns

Factor columns capture experimental design variables from `factor value[X]` columns in the SDRF. Each one becomes a `factor_X` column of `sample.parquet` (lower-case, with runs of other characters replaced by an underscore), so a factor stays distinguishable from a characteristic of the same name. `qpxc transform quantify --factor X` chooses the factor used as the condition.

| SDRF Column | QPX Field Name | View |
|---|---|---|
| `factor value[disease]` | `factor_disease` | sample.parquet |
| `factor value[organism part]` | `factor_organism_part` | sample.parquet |

!!! note
    Factor columns are open-ended: any `factor value[X]` column maps to `factor_X`. A planned `factor.parquet` view may later add the factor's type and unit. Replicate identifiers (`biological_replicate`, `technical_replicate`) are stored in the `samples` struct of [run.parquet](run.md), so replicate structure is already self-contained there.

---

## Multi-valued Metadata

In rare cases a single SDRF row may describe a sample derived from multiple organisms (e.g. a spike-in experiment) or multiple cell lines (pooled samples). QPX handles this natively using **list columns**.

### Convention

Biological property columns in `sample.parquet` are strings. When several SDRF columns give different values for one sample (for example `characteristics[organism]` and a second organism column), the values are joined with `"; "`:

```python
# Single-valued (most samples)
organism = "Homo sapiens"

# Multi-valued (spike-in / pooled samples)
organism = "Homo sapiens; Saccharomyces cerevisiae"
```

### SQL Queries

```sql
-- All unique organisms in the project
SELECT DISTINCT UNNEST(string_split(organism, '; ')) AS organism
FROM 'PXD014414.sample.parquet';

-- Samples with a specific disease
SELECT sample_accession, disease
FROM 'PXD014414.sample.parquet'
WHERE list_contains(string_split(disease, '; '), 'breast cancer');
```

---

## SDRF Conversion Process

At ingestion time, the SDRF TSV is split into `sample.parquet` and `run.parquet`:

1. **Load** the SDRF TSV file
2. **Rename** columns using the field registry above
3. **Extract sample-level columns** -- deduplicate rows by `sample_accession` (= `source name`) to produce `sample.parquet`
4. **Extract run-level columns** -- one row per `run_file_name`, with the `samples` list linking to sample accessions, labels, and replicate numbers, to produce `run.parquet`
5. **Parse ontology terms** -- instrument, enzyme, dissociation and acquisition cells become plain strings on the run, and their resolved terms go to `additional_terms`
6. **Preserve provenance** -- the SDRF the conversion was given is copied to `{prefix}.sdrf.tsv`, and its SHA-256 is recorded in `provenance.parquet`

A blank required `organism` or `organism_part` is written as `not available` with a warning, so one incomplete row does not lose the views. Run-level `comment[...]` columns that have no QPX field stay available in the SDRF copy.

### Project Layout

```text
PXD014414/
  PXD014414.sdrf.tsv             # Original SDRF (provenance)
  PXD014414.sample.parquet       # Biological samples (deduplicated)
  PXD014414.run.parquet          # Data acquisition runs with ontology terms
  PXD014414.dataset.parquet      # Project-level metadata
  PXD014414.feature.parquet
  PXD014414.psm.parquet
  ...
```

---

## Related Pages

- [Sample Metadata](sample.md) -- `sample.parquet`
- [Run Metadata](run.md) -- `run.parquet`
- [Dataset Metadata](dataset.md) -- `dataset.parquet`
- [QPX Format Overview](index.md)
