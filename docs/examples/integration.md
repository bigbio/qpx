# Integration Examples

Examples of integrating qpx with other tools and frameworks.

## Nextflow Pipeline Integration

Real production Nextflow pipeline from the QPX project for processing MaxQuant PSM data.

### Pipeline Features

- Converts Thermo RAW files to mzML format using ThermoRawFileParser
- Processes MaxQuant PSM data to QPX parquet format
- Extracts spectral information from mzML files
- Supports containerization (Docker/Singularity)
- Configurable resource allocation for HPC environments
- Automatic error handling and completion reporting

### Run the Pipeline

```bash
# Using Docker
nextflow run nf-mq-psm.nf -profile docker

# Using Singularity
nextflow run nf-mq-psm.nf -profile singularity

# With custom parameters
nextflow run nf-mq-psm.nf \
    --raw_dir ./raw_files \
    --mzml_dir ./mzml_output \
    --output_dir ./results \
    --msms_file ./maxquant/msms.txt \
    --chunksize 2000000 \
    -profile docker
```

**Full pipeline source**: [nextflow/nf-mq-psm/](https://github.com/bigbio/qpx/tree/main/nextflow/nf-mq-psm)

## Python API Usage

Use qpx programmatically in Python scripts. The central entry point is `qpx.open_dataset()`, which discovers all Parquet structures in a dataset directory and registers them in a DuckDB engine for fast analytical queries.

### Opening and Exploring a Dataset

```python
import qpx

# Open a dataset directory — auto-discovers all Parquet structures
ds = qpx.open_dataset("PXD014414/")

# See which structures are available
print(ds.available_structures)
# ['feature', 'pg', 'sample', 'run', 'dataset']

# Query individual structures
print(f"Features: {ds.feature.count()}")
print(f"Protein groups: {ds.pg.count()}")
print(f"Samples: {ds.sample.count()}")

# Convert to pandas DataFrame
feature_df = ds.feature.to_df()
print(feature_df.head())

# Filter with SQL-like predicates
high_charge = ds.feature.filter("charge > 2")
print(f"Features with charge > 2: {high_charge.count()}")

# Run arbitrary SQL
result = ds.sql("SELECT sequence, COUNT(*) AS n FROM feature GROUP BY sequence ORDER BY n DESC LIMIT 10")
print(result.to_df())
```

### Built-in Views and Plotting

QPX provides pre-built views for common summaries, each with a `.plot()` method for quick visualization.

```python
import qpx

ds = qpx.open_dataset("PXD014414/")

# Identification summary (per-run protein/peptide counts)
summary = ds.identification_summary.summary()
print(summary.to_df())

fig = ds.identification_summary.plot()
fig.write_html("identifications.html")

# Run summary (per-run statistics)
fig = ds.run_summary.plot()
fig.write_html("run_summary.html")

# Modification frequency across PSMs
fig = ds.modification_view.plot(top_n=20)
fig.write_html("modifications.html")

# Dataset QC metrics
fig = ds.qc_view.plot()
fig.write_html("qc_overview.html")
```

### Multi-Dataset Analysis

Use `DatasetCollection` to query across multiple datasets or physically merge them.

```python
import qpx

ds1 = qpx.open_dataset("PXD014414/")
ds2 = qpx.open_dataset("PXD016999/")

# Virtual mode — cross-dataset SQL
coll = qpx.DatasetCollection([ds1, ds2])
result = coll.sql("SELECT source_dataset, COUNT(*) AS n_features FROM feature GROUP BY source_dataset")
print(result.to_df())

# Physical merge — combine structures into a new directory
coll.merge("merged_output/", structures=["feature", "pg"])

merged = qpx.open_dataset("merged_output/")
df = merged.feature.to_df()
print(df["source_dataset"].value_counts())
```

### Dataset Integrity

Compute and verify checksums, row counts, and file sizes for a dataset.

```python
import qpx

ds = qpx.open_dataset("PXD014414/")

# Compute integrity metadata
integrity = ds.compute_integrity()
print(f"Files tracked: {integrity['total_structures']}")
print(f"Packaged at: {integrity['packaged_at']}")

# Verify an existing dataset
result = ds.verify_integrity()
if result["errors"]:
    print("Integrity errors:", result["errors"])
else:
    print("Dataset integrity verified.")
```

---

[← Back to Examples Overview](index.md) | [View QC Examples →](qc.md)
