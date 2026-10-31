# Serialization & Parquet Details

QPX uses three serialization formats, each chosen for its suitability to the data it stores: Apache Parquet for columnar proteomics data and metadata, AnnData (`.h5ad`) for the standalone expression views, and MuData (`.h5mu`) for the combined multi-modal expression container.

## Supported Serialization Formats

| Format | Extension | Used By | Description |
|--------|-----------|---------|-------------|
| Apache Parquet | `.parquet` | PSM, Feature, PG, Peptide, Protein, MZ, Dataset, Sample, Run | Columnar storage with compression |
| AnnData (HDF5) | `.h5ad` | AE, DE | Matrix-form expression data (scverse ecosystem) |
| MuData (HDF5) | `.h5mu` | Combined expression container | Multi-modal bundle of the `precursors`, `proteins`, `expression`, and `differential` modalities plus a `varp["feature_mapping"]` precursor↔protein adjacency (see [MuData](anndata.md#mudata-h5mu)) |

!!! note
    The choice of format for each view is driven by its access patterns. Data views and metadata use Parquet for column-level filtering and large-scale analytical queries. Expression views use AnnData (`.h5ad`) for native interoperability with the scverse ecosystem (scanpy, scvi-tools, muon). The original SDRF (`.sdrf.tsv`) is preserved alongside the Parquet files for provenance but is not a QPX output format.

## Parquet Format

Apache Parquet is a columnar storage format designed for efficient analytical processing. Unlike row-oriented formats (CSV, TSV, mzTab), Parquet stores data column by column, which provides significant advantages for proteomics workflows:

- **File metadata and column metadata**: Each Parquet file includes a footer that describes the schema, row groups, column statistics (min/max values), and custom key-value metadata. This allows query engines to skip irrelevant data without reading the full file.
- **Column-oriented design**: Analytical queries that access a subset of columns (e.g., retrieving only `sequence` and `charge` from a PSM file) read only the relevant column chunks, reducing I/O dramatically compared to row-based formats.
- **Broad ecosystem support**: Parquet is natively supported by PyArrow, DuckDB, Polars, Apache Spark, pandas, R (via the `arrow` package), and many other tools. This makes QPX files immediately usable across languages and frameworks without custom parsers.

## Compression

Parquet's columnar layout achieves substantial compression ratios because values within a column tend to have similar data types and distributions. The following benchmarks illustrate the storage savings and write performance for real-world proteomics datasets:

| Project | Type | Original Size (GB) | Parquet Size (MB) | PSM Write Time (s) | Feature Write Time (s) |
|---------|------|--------------------:|-------------------:|--------------------:|-----------------------:|
| PXD046440 | MaxQuant | 48 | 337/343 | 985.27 | 678.47 |
| PXD016999 | mzTab | 160 | 155/228 | 539.00 | 3554.53 |
| PXD019909 | DIA-NN | 1.9 | 195 | -- | 229.48 |

!!! tip "Compression algorithms"
    QPX supports the following compression algorithms for Parquet files:

    - **Zstd** -- Best balance of compression ratio and speed; broad cross-language support (default)
    - **Snappy** -- Fast compression/decompression with moderate ratio
    - **Gzip** -- Higher compression ratio, slower speed
    - **LZO** -- Optimized for decompression speed
    - **None** -- No compression

    In addition to block-level compression, Parquet applies encoding schemes at the column level:

    - **RLE (Run-Length Encoding)** -- Efficient for columns with many repeated values (e.g., `charge`). Boolean columns like `is_decoy` use native bit-packing (1 bit per value) with optional RLE for maximum compression
    - **Dictionary encoding** -- Replaces repeated string values with integer codes (e.g., `run_file_name`, `sequence`)

## Parquet Features

QPX leverages the following Parquet capabilities:

- **Columnar Storage**: Each column is stored independently, enabling query engines to read only the columns needed for a given analysis. This reduces I/O for analytical queries that access a subset of fields.
- **Efficient Compression**: Column-level compression with Zstd, Snappy, Gzip, or LZO, combined with RLE and dictionary encoding, achieves compression ratios of 10x--100x on typical proteomics datasets.
- **Schema Evolution**: Columns can be added, removed, or modified across QPX versions without breaking compatibility with existing files. Readers that encounter unknown columns can safely ignore them.
- **Complex Data Types**: Parquet natively supports nested structs, arrays, and maps. QPX uses these for fields like `modifications` (array of structs), `additional_scores` (array of structs), and spectral arrays (`mz_array`, `intensity_array`).

## Hive Partitioning

QPX supports **Hive-style partitioning** for large datasets. This splits a single Parquet file into a directory tree where each subdirectory encodes a partition column value (e.g., `run_file_name=run_01/`). This is particularly useful for datasets with 100+ runs where a single Parquet file would exceed several GB.

Partitioning applies to the **per-run** views — `psm` and `feature` — which carry a scalar `run_file_name`. The **`pg` (protein-group) view is not Hive-partitionable**: since QPX 1.1 it is keyed on `grouped_runs` (`list<string>`), and Hive partitioning encodes each partition value as a directory name, which cannot represent a list. Write `pg` as a single Parquet file. (`write_partitioned` raises an actionable error if asked to partition a `pg` table.)

### Directory Structure

The default partition column is `run_file_name` (for `psm` / `feature`):

```text
PXD004683/
├── exp.pg.parquet                   # Single-file structure (small)
├── exp.sample.parquet
├── exp.run.parquet
├── exp.dataset.parquet
└── feature/                         # Partitioned structure (large)
    ├── run_file_name=run_01.raw/
    │   └── part-1f3a9c2e-0.parquet
    ├── run_file_name=run_02.raw/
    │   └── part-1f3a9c2e-0.parquet
    └── run_file_name=run_03.raw/
        └── part-1f3a9c2e-0.parquet
```

### Auto-Discovery

When opening a dataset, QPX checks for single Parquet files first. If no single file is found for a structure, it falls back to checking for a Hive-partitioned directory with the same name (e.g., `feature/`).

```python
import qpx

# Automatically discovers both single files and partitioned directories
ds = qpx.open_dataset("PXD004683/")
print(ds.feature.count())  # Works regardless of storage layout
```

### Writing Partitioned Data

```python
from qpx.writers import FeatureWriter
import pyarrow.parquet as pq

# Read existing data
table = pq.read_table("exp.feature.parquet")

# Write as Hive-partitioned; ids, the identity composite and the footer are kept
FeatureWriter.write_partitioned(table, "feature/", partition_cols=["run_file_name"])
```

A directory that already holds Parquet files raises by default. `mode="append"`
adds to it, each call under its own file names, so repeated calls stream a large
table into one directory. `mode="replace_partitions"` replaces the partitions the
table writes and keeps the others. `qpxc transform quantify --feature-path` also
accepts a directory partitioned by `run_file_name`.

No converter writes partitioned output: partition a view after the conversion.
A run name read back from a partition path stays a string, even if it looks like
a number or a date.

### DuckDB Integration

DuckDB reads Hive-partitioned directories natively:

```sql
-- The partition column (run_file_name) is automatically available
SELECT run_file_name, COUNT(*) AS n_features
FROM read_parquet('feature/**/*.parquet', hive_partitioning=true)
GROUP BY run_file_name;
```

!!! warning "Practical limits of Parquet metadata"
    Parquet stores all file-level and column-level metadata in the file footer. When reading a partitioned dataset, the footer of every partition file must be loaded into memory to reconstruct the full schema and row group statistics. For datasets with thousands of partitions, this metadata overhead can become significant. Consider using fewer, larger partitions rather than many small ones.

## Writing from Python

Each view has a writer in `qpx.writers`. Rows are written to a staging file that
replaces the destination when the writer closes without an error.

```python
from qpx.writers import FeatureWriter

with FeatureWriter("exp.feature.parquet", identity_composite=("peptidoform", "charge", "run_file_name")) as writer:
    writer.write_batch(records)
```

- A record key outside the schema raises when the record is written. Keys that
  start with an underscore are ignored, for hints between the steps of a converter.
- A value the schema cannot hold raises when its batch is written, naming the
  column, the row of the batch and the value. That can be at `close()`, which
  writes the last batch.
- Non-fatal findings, such as scans that violate the declared `scan_format`, are
  logged once per file.
- `FeatureWriter` re-derives identified feature ids by default, and keeps a provided
  one as a `provided_feature_id` cv_param. Other writers keep provided ids. Pass
  `override_provided_ids=False` to keep the feature ids, for instance when
  `psm.feature_id` already refers to them.
- `FeatureWriter.derive_ids(table, identity_composite=...)` returns the table with
  the ids the writer would derive, so that `psm.feature_id` can be filled first.

`Dataset.save_structure()` replaces a view of a dataset: it writes under the
prefix of the dataset's files, keeps the rows' ids and takes the identity composite,
`creator` and `scan_format` from the existing file. Another `prefix` adds a file
that is read together with the view, with a warning.

## Cloud / S3 Access

QPX datasets can be opened directly from S3-compatible cloud storage without downloading files locally. This uses DuckDB's `httpfs` extension under the hood.

### Opening a Dataset from S3

```python
import qpx

ds = qpx.open_dataset(
    "s3://my-bucket/datasets/PXD014414/",
    s3_config={
        "region": "us-east-1",
        "access_key_id": "AKIA...",
        "secret_access_key": "...",
    },
)

# Use exactly like a local dataset
print(ds.feature.count())
df = ds.feature.filter("charge > 2").to_df()
```

### Configuration Options

| Parameter | Description |
|-----------|-------------|
| `region` | AWS region (e.g., `"us-east-1"`) |
| `access_key_id` | AWS access key ID |
| `secret_access_key` | AWS secret access key |
| `endpoint` | Custom S3 endpoint for MinIO or other S3-compatible services |
| `anonymous` | Set `True` for public buckets (no credentials needed) |

### S3-Compatible Services

For MinIO or other S3-compatible services, set the `endpoint` parameter:

```python
ds = qpx.open_dataset(
    "s3://proteomics/PXD014414/",
    s3_config={
        "endpoint": "minio.example.com:9000",
        "access_key_id": "minioadmin",
        "secret_access_key": "minioadmin",
    },
)
```

!!! tip "Performance"
    Cloud access adds network latency to every query. For interactive analysis, consider downloading the dataset locally first. Cloud access is most useful for automated pipelines and one-off queries against large remote datasets.

## File-Level Metadata

Every QPX Parquet file includes metadata as key-value pairs stored in the Parquet file footer. These fields identify the file's origin, format version, and provenance:

| Key | Description |
|-----|-------------|
| `qpx_version` | Version of the QPX specification used to generate the file |
| `software_provider` | Software that wrote the file (`qpx` for the qpx writers) |
| `writer_version` | Version of the qpx package that wrote the file |
| `project_accession` | ProteomeXchange or project accession (e.g., `PXD012345`) |
| `project_title` | Human-readable title of the project |
| `scan_format` | Format of scan identifiers: `scan`, `index`, or `nativeId` |
| `primary_key` | Primary-key column for identity-bearing views (`feature_id`, `psm_id`, or `pg_id`) |
| `identity_composite` | Comma-separated producer properties used when QPX derives the row ID |
| `creator` | Tool whose output the data was converted from (e.g. `maxquant`, `diann`) |
| `file_type` | QPX view type (e.g., `psm_file`, `feature_file`, `pg_file`) |
| `creation_date` | ISO 8601 date when the file was created |
| `compression_format` | Compression algorithm used: `zstd`, `snappy`, `gzip`, `lzo`, or `none` |

For QPX-derived identity columns, the ordered `identity_composite` values are
encoded as compact canonical JSON with escaped strings and JSON `null` values.
QPX hashes those bytes with BLAKE2b using an 8-byte digest and stores the result
as a signed, big-endian `int64`. `grouped_runs` is sorted for identity hashing;
ordered lists such as `scan` retain their order.

One identity composite value names one row identity: rows with equal composite
values carry the same id, also when the producer supplies its own ids.

`qpxc validate` checks the footer of every local file: an unparseable
`qpx_version` or another major version than the reader's, a `file_type` of another
view, an `identity_composite` naming fields outside the view and a `primary_key`
other than the view's are errors; a missing `qpx_version`, `file_type`,
`identity_composite` or `primary_key` fails strict validation.

!!! tip "Writing Parquet with metadata in Python"
    ```python
    import pyarrow as pa
    import pyarrow.parquet as pq

    # Build your table
    table = pa.table({"sequence": ["PEPTIDER", "ANOTHERPEPTIDE"]})

    # Define file-level metadata
    file_metadata = {
        "qpx_version": "1.3",
        "software_provider": "qpx",
        "creator": "quantms",
        "project_accession": "PXD012345",
        "file_type": "psm_file",
        "creation_date": "2021-01-01",
    }

    # Merge with existing schema metadata and write
    existing_metadata = table.schema.metadata or {}
    merged_metadata = {**existing_metadata, **{k.encode(): v.encode() for k, v in file_metadata.items()}}
    table = table.replace_schema_metadata(merged_metadata)
    pq.write_table(table, "psm_data.parquet")
    ```

!!! tip "Reading file metadata in Python"
    ```python
    import pyarrow.parquet as pq

    parquet_file = pq.ParquetFile("psm_data.parquet")
    metadata = parquet_file.schema_arrow.metadata
    for key, value in metadata.items():
        print(f"{key.decode()}: {value.decode()}")
    ```

## See Also

- [File Naming](file-naming.md) -- file extension and naming conventions for QPX files
- [Versioning](versioning.md) -- version format and backward compatibility rules
- [QPX Format Overview](index.md) -- overview of all QPX views and their relationships
