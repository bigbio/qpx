# Dataset Metadata

`dataset.parquet` is a **single-row** Parquet file capturing project-level metadata for a QPX dataset -- project identity and software provenance. It does **not** store per-run or per-sample metadata -- those live in [run.parquet](run.md) and [sample.parquet](sample.md), respectively.

See the full YAML schema in [`dataset.yaml`](schemas/dataset.yaml).

!!! info "Aggregated biology is derived, not stored"
    Fields like organisms, organism parts, diseases, and cell lines are **not** stored in `dataset.parquet`. These are derived at query time from [sample.parquet](sample.md):

    ```sql
    SELECT DISTINCT UNNEST(organism) AS organism FROM 'PXD014414.sample.parquet';
    ```

    This eliminates the risk of project metadata drifting out of sync with the actual sample annotations. Similarly, instrument, enzyme, and modification information lives in [run.parquet](run.md) since it can vary per run.

---

## PyArrow Schema

```python
import pyarrow as pa

dataset_schema = pa.schema(
    [
        # Project identity
        pa.field("project_accession", pa.string()),
        pa.field("project_title", pa.string(), nullable=True),
        pa.field("project_description", pa.string(), nullable=True),
        pa.field("pubmed_id", pa.string(), nullable=True),
        # Software
        pa.field("software_name", pa.string(), nullable=True),
        pa.field("software_version", pa.string(), nullable=True),
        # Provenance
        pa.field("creation_date", pa.string()),
        # Integrity / Packaging
        pa.field("file_checksums", pa.map_(pa.string(), pa.string()), nullable=True),
        pa.field("file_row_counts", pa.map_(pa.string(), pa.int64()), nullable=True),
        pa.field("file_sizes_bytes", pa.map_(pa.string(), pa.int64()), nullable=True),
        pa.field("total_structures", pa.int32(), nullable=True),
        pa.field("packaged_at", pa.string(), nullable=True),
    ]
)
```

### Field Reference

#### Core Identity

| Field | Description | Type | Required |
|---|---|---|---|
| `project_accession` | Project accession identifier (e.g. `PXD014414`) | string | yes |
| `project_title` | Title of the project | string | no |
| `project_description` | Description of the project | string | no |
| `pubmed_id` | PubMed ID associated with the project | string | no |
| `analysis_id` | Names this analysis, new for each conversion and merge ([Row Identity](identity.md#scope)) | string | no |

The specification version is not a column: it is the `qpx_version` footer key of every file ([Versioning](versioning.md#version-in-files)).

#### Software Provenance

| Field | Description | Type | Required |
|---|---|---|---|
| `software_name` | Name of the software that generated the data (e.g. `quantms`) | string | no |
| `software_version` | Version of the software | string | no |

#### File Provenance

| Field | Description | Type | Required |
|---|---|---|---|
| `creation_date` | ISO 8601 date when the dataset was created | string | yes |

#### Integrity / Packaging

These fields enable dataset validation after packaging or transfer. They are populated by calling `ds.compute_integrity()` and verified with `ds.verify_integrity()`.

| Field | Description | Type | Required |
|---|---|---|---|
| `file_checksums` | SHA-256 hex digests keyed by file name (relative to dataset root) | map&lt;string, string&gt; | no |
| `file_row_counts` | Row counts keyed by file name | map&lt;string, int64&gt; | no |
| `file_sizes_bytes` | File sizes in bytes keyed by file name | map&lt;string, int64&gt; | no |
| `total_structures` | Number of Parquet structures in this dataset | int32 | no |
| `packaged_at` | ISO 8601 timestamp when integrity was computed | string | no |

!!! example "Computing and verifying integrity"
    ```python
    import qpx

    ds = qpx.open_dataset("PXD014414/")

    # Compute integrity (checksums, row counts, file sizes)
    integrity = ds.compute_integrity()

    # Save integrity data into dataset.parquet
    meta_df = ds.dataset_meta.to_df()
    meta_dict = meta_df.iloc[0].to_dict()
    meta_dict.update(integrity)
    ds.save_structure([meta_dict], "dataset", prefix="PXD014414")
    ds.refresh()

    # Later: verify the dataset is intact
    result = ds.verify_integrity()
    if result["errors"]:
        print("Integrity check failed:", result["errors"])
    else:
        print("Dataset integrity verified")
    ```

!!! tip "Reading dataset.parquet in Python"
    ```python
    import pyarrow.parquet as pq

    dataset = pq.read_table("PXD014414.dataset.parquet")
    row = dataset.to_pydict()

    print(row["project_accession"])  # ['PXD014414']
    print(row["software_name"])  # ['quantms']
    ```

!!! tip "Reading dataset.parquet with DuckDB"
    ```sql
    SELECT project_accession, software_name, software_version
    FROM 'PXD014414.dataset.parquet';
    ```

---

## Related Pages

- [Sample Metadata](sample.md) -- `sample.parquet`
- [Run Metadata](run.md) -- `run.parquet`
- [SDRF Conversion](sdrf-mapping.md) -- column mapping and ingestion process
- [QPX Format Overview](index.md)
