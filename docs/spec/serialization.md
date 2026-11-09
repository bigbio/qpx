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

## Physical Layout

The physical layout is not normative: any row order and row-group size is a valid QPX file. It does decide what a selective read costs.

QPX writers use row groups of at most 131,072 rows (fewer when rows are large, at most about 256 MB uncompressed per group), a page index, and for each column the encoding its first row group measures smallest: a dictionary, `PLAIN`, `BYTE_STREAM_SPLIT` (floats) or `DELTA_BINARY_PACKED` (integers); a string keeps its dictionary unless `PLAIN` is more than 20% smaller. No fixed choice per column suits every file: TMT reporter intensities are smallest byte-stream-split, and PEP is smallest as a dictionary in DIA-NN files but byte-stream-split in OpenMS ones. Published 1.1 files rewritten this way took 2–14% less space (MSV000085836 feature 674.6 to 598.9 MB, pg 65.1 to 55.8 MB; PXD017199 feature 1,596.8 to 1,562.4 MB). On a 783k-row dataset, sorting with these row groups made protein lookups 3–4× faster and read 5–7× fewer remote bytes than unsorted 1 Mi-row groups, for about 3% more bytes on disk.

For a published dataset, `qpxc transform relayout` sorts the views and records the order in `sorting_columns` and in the footer key `qpx.sort_order`:

| View | Sorted by | Good for |
|---|---|---|
| `feature` | `anchor_protein`, `peptidoform`, `charge`, `run_file_name` | Looking up a protein or a peptide |
| `pg` | `anchor_protein`, `grouped_runs`, `label` | Looking up a protein |
| `psm` | `peptidoform`, `charge`, `run_file_name`, `scan` | Looking up a peptide |
| `mz` | `run_file_name`, `scan` | Reading a run, looking up a scan |

`sorting_columns` cannot name a list column, so for `pg` it records `anchor_protein` only. relayout also numbers the ids of feature, pg and psm 1, 2, 3, … in their new row order, ties broken by the old id, and rewrites every column that holds them: `pg_ids` and `psm_ids` of feature, `feature_id` of psm and fragment, `pg_id` and `best_psm_id` of site and proteoform, and `feature_ids` and `feature_scores[].feature_id` of xic. A reference to a row the view lacks keeps its value. The footer then states `identity_source` `provided`, since the ids no longer derive from the composite (see [Row Identity](identity.md)); so does the footer of fragment, whose ids keep their values but were hashed from the old `feature_id`. Like any id they mean something only within their file ([Row Identity](identity.md#scope)), and every relayouted dataset numbers from 1, so the ids of two datasets collide: joined on `psm_id`, the relayouted psm views of PD PXD046576 and PXD081613 match 16,359 rows, their converter outputs none. In a `DatasetCollection`, join on `source_dataset` and the id, and merge two relayouted datasets with `on_collision='namespace'`. Other values do not change, apart from the [rounded measurements](#rounded-measurements).

The same input gets the same numbers, but a change to the input moves nearly all of them: PXD081613 converted without one of its runs kept 60 of 104,322 psm ids and 25 of 41,444 feature ids after relayout, where the hashed ids all stayed. Keep links to rows on the composite or a natural key. Hashing the composite of a relayouted feature, pg or psm row gives back the id it had before, if that id was derived, as composite columns are not rounded. A view relayout does not sort keeps its ids; they derive from its composite as before unless the composite holds a rewritten reference, as fragment's does.

Rows of one peptide sitting together compress better, and a hashed id is eight random bytes per row that no encoding shrinks, while ordinals delta-encode to almost nothing. On MSV000085836 (TMT, 6.2 million features) sorting took feature from 674.6 to 430.7 MB and psm from 390.7 to 245.4 MB, and numbering the ids took them on to 381.1 and 150.7 MB and pg from 55.8 to 17.0 MB, in 38 s. With hashed ids, pg sorted by its default keys stays at the 55.8 MB of its file order, which keeps a group's runs together; sorted by `anchor_protein` and `label` alone it measured 62.2 MB.

Choose by the reads you expect, because one order serves one kind of query well. On that sorted feature view, one protein's rows took 0.06 s to read instead of 2.12 s, and one run's 2.23 s instead of 0.22 s (DuckDB, 24 threads, warm cache). Sorting by protein costs per-run reads their pruning; Hive partitioning by run is fastest per run, and measured 57% more bytes and 3.5× slower for proteins. Partition by a batch of runs rather than one file per run, and only for a dataset mostly read run by run.

The sort also changes what a query that relies on row order returns. Without `ORDER BY`, `LIMIT` and `first()` can return other rows after relayout: none of the first 10 rows of PD PXD081613 psm or DIA-NN feature stayed the same. `ORDER BY` on an id now follows the sort keys, where hashed ids gave a pseudo-random order. A query that needs an order should name it.

### Struct lists

Unlike the rest of the layout, this changes the type a reader sees, so it is part of the specification from 1.3. relayout writes a name-value list, `additional_scores`, `cv_params` or `additional_intensities`, whose rows all list the same names in the same order as a struct with one field per name: `additional_scores.diann_cscore`, `additional_intensities.LFQ.maxlfq`. The values of one name then share a column and its encoding, where the scores of a DIA-NN feature interleaved in one column compress poorly. A null list is a null struct. A value is stored as `float` when every one of its values reads back unchanged, else as `double`; a `cv_value` string becomes a number only when it reads back as the same string, so a `proteotypic` of `"1"` stays a string. The footer key `qpx.struct_lists` gives each such column's names in order as JSON, with the `higher_better` of each score:

```json
{"additional_scores": [{"name": "diann_cscore", "higher_better": true}],
 "cv_params": [{"name": "irt"}],
 "additional_intensities": [{"name": "LFQ", "names": ["maxlfq", "ms1_area"]}]}
```

QPX reads such a column back as its list, each file by its own footer (`Dataset`, `DatasetCollection`, `qpxc query` and `qpxc validate`, MuData, `merge()`, gene-map, migrate), so the API and SQL through QPX do not change, also for a view whose files differ in layout or store a field as different types; such a view reads its files in groups of one layout and storage and unions them, so its rows follow the groups rather than the file order. A file that holds such a struct without the footer key, as a DuckDB `COPY` or a pandas round trip of a relayouted file leaves it, reads as the struct, and QPX refuses to read it in one view with files that hold the column as a list or with a footer layout. A column keeps its list when its layout varies, such as a feature `additional_intensities` without `ms1_area` on some rows, when a name holds a character such as `-` (OpenMS psm `q-value`), or when two names differ only in case, which DuckDB reads as two other names. On PXD071205 (DIA-NN, 6.8 million features) struct lists took the sorted feature view from 406.1 to 311.2 MB, and relayout from 32 to 17 s and from 9.0 to 5.0 GB of memory, as DuckDB sorts structs more cheaply than lists. Run relayout last: gene-map and migrate rewrite a view with its lists and without the recorded sort order.

DuckDB, Polars or PyArrow opening the file see the struct; a reader that expects the list must turn the struct back first. In DuckDB a value is a field, `additional_scores.precursor_qvalue` (DIA-NN feature) or `additional_intensities.LFQ.spectral_count` (PD feature), and a name that is not a plain identifier is quoted, as MaxQuant psm `cv_params."mass analyzer"`. List functions such as `list_filter` and `len` fail on the struct, and so does `unnest` in a `FROM` clause or under a field access; at the root of a `SELECT`, `unnest` raises no error but returns the fields as columns, one row per row of the file instead of one per entry. Field access fails on the list. Which form a column has depends on the rows of the file, not on its converter: the PD converter gave a struct `additional_scores` for PXD046576 psm and a list for PXD081613 psm. Check `DESCRIBE` or the footer before writing SQL for one form. In the struct, the `higher_better` of a score is only in the footer, whose values DuckDB returns as BLOBs:

```sql
SELECT decode(value) FROM parquet_kv_metadata('d.feature.parquet') WHERE key = 'qpx.struct_lists';
```

The flat q-value and PEP columns ([Typed q-value columns](scores.md#typed-q-value-columns)) read the same in both forms. Field access is also the fast form: on DIA-NN feature (195,691 rows) `additional_scores.precursor_qvalue < 0.001` took 3.24 ms, where `list_filter` on the list took 114.5 ms (DuckDB 2.0.0, 4 cores, warm cache; 3.95 and 93.6 ms with 1.5.0). Through QPX, which reads the list, it took 85 ms.

Read several relayouted files in one `read_parquet` with `union_by_name = true`. Without it DuckDB 1.4.4 to 2.0.0 take the struct type of the first file: fields only in later files are dropped, with no error unless a query names one, and a field the first file stores as `float` reads as `float`. After the PD PXD046576 feature view, the 3,131,056 score values of the DIA-NN one read as null. With it the fields merge, and `float` with `double` reads as `double`. DuckDB 1.1.3 stops with a type mismatch in both modes. A list file and a struct file of one column cannot be read in one `read_parquet` at all, and the error comes only once the column is read, after `DESCRIBE` and `count(*)` have worked. Relayouting every file does not prevent this, as each file's rows decide its form: the relayouted PD PXD046576 and PXD081613 psm files hold `additional_scores` as a struct and as a list. `DESCRIBE` of the combined read shows the list type in either file order, so check each file's `DESCRIBE` or footer before reading files together, or read them through QPX (`Dataset`, `DatasetCollection`), which restores each file by its own footer.

### Rounded measurements

relayout rounds the measured values of the views it sorts to the precision of the measurement, as the bits below it are noise that no encoding compresses. Intensities and mass errors keep 12 mantissa bits, a relative error of at most 2^-13 (1.2 × 10^-4); retention times, ion mobility and CCS keep 16, at most 2^-17, and so do the iRT and ion mobility terms of `cv_params` once it is a struct. Values are rounded to nearest with ties to even, as numcodecs' `BitRound` does, so without bias, and stay ordinary floats. On MSV000085836 the log2 ratio between repeated PSMs of one precursor varies by 0.22 (MAD), 3,000 times the error 12 bits add to it (RMS 7.1 × 10^-5). Scores, q-values, PEP, m/z and every other value keep all their bits, so the rows that pass a threshold do not change, and so does a column of the view's identity composite, whose values name its rows (`rt` and `consensus_rt` of MSV000085836 feature). The footer key `qpx.mantissa_bits` gives the bits of each rounded column as JSON, a path naming the column and its struct fields:

```json
{"intensities.intensity": 12, "mass_error_ppm": 12, "rt": 16}
```

Rounding took the sorted MSV000085836 feature view from 381.1 to 294.6 MB, psm from 150.7 to 142.8 MB and pg to 11.7 MB, and PXD071205 feature from 311.2 to 263.1 MB, without slowing relayout. The dataset's h5mu, which holds feature and pg values, is built again. gene-map and migrate keep the rounded values and the footer key.

An exact comparison with a rounded value stops matching. Within PD PXD081613 a join of psm to feature on run, peptidoform, charge and `rt` matched 43,121 rows before relayout and 336 after, because psm `rt` is rounded and feature `rt`, a column of its composite, is not. An exact filter on the iRT 55.2095 of DIA-NN feature matched 22 rows before and none after; a range of ±0.0001 around it matched the 22 in both. Compare a rounded column with a relative tolerance, 2^-17 where it keeps 16 bits and 2^-13 where it keeps 12 (add 2^-24 to compare a `cv_params` term stored as `float` with its original string), or join on ids and natural keys. Through QPX a struct `cv_params` reads back as the strings of its rounded values: DIA-NN iRT `52.8228` becomes `52.822754`.

### Size of the long pg format

`pg` is long: one row per label, each with its own `pg_id`, so a TMT file is larger than the nested form it replaced (a factor of 2.34 on a TMT dataset, nearly all of it the hashed ids: without them the long form is 5% smaller). A value costs 14–21 bytes in a realistic single-cell matrix. Sorting by `anchor_protein` took 34% off the same TMT file, and `qpxc transform relayout`, which also numbers the ids, took MSV000085836 pg from 55.8 to 17.0 MB. Only the primary quantity is a row; `additional_intensities` cannot be pivoted the same way.

### Reader compatibility

Writers choose each column's encoding from its data (see [Physical Layout](#physical-layout)), so any float column may be `BYTE_STREAM_SPLIT` and any integer column `DELTA_BINARY_PACKED`; the others are `RLE_DICTIONARY` or `PLAIN`. Quantities and scores sit in nested `list<struct>` columns, or in structs after relayout ([Struct lists](#struct-lists)).

| Reader | QPX files | Checked with |
|---|---|---|
| PyArrow, pandas with `engine="pyarrow"` | Reads everything | 24.0, read in the test suite |
| DuckDB | Reads everything | 1.1.3, 1.4.4, 1.5.0, 1.5.6, 2.0.0 |
| Polars | Reads everything, nested columns included | 1.41 |
| R `arrow` | Reads everything, nested columns included | 25.0 |
| Apache Spark | Reads everything, nested columns included, with `spark.sql.parquet.enableVectorizedReader=false`. The default vectorized reader **fails** on `BYTE_STREAM_SPLIT` columns, which published 1.1 files already have for `rt` and m/z: `Unsupported encoding: BYTE_STREAM_SPLIT` | 4.2.1 |
| fastparquet | **Fails** on `BYTE_STREAM_SPLIT` columns: `NotImplementedError: Encoding 9` (it reads `DELTA_BINARY_PACKED`). Even without that encoding it returns `None` for every `list<struct>` value, silently; strings and scalars read | 2026.3 |

Each reader read every column of a relayouted DIA-NN feature view (struct lists) and a TMT feature and psm view (`list<struct>`), and gave DuckDB's row counts and sums. To read with Spark, turn its vectorized reader off: `spark.conf.set("spark.sql.parquet.enableVectorizedReader", "false")`. To read with fastparquet, rewrite through PyArrow without the encoding, and expect nested columns to be lost.

DuckDB read 101 files without an error in every version listed: 98 written by OpenMS's QPX export, by the DIA-NN, PD, MaxQuant, FragPipe, OpenMS and Spectronaut converters, and by relayout from them, and 3 synthetic files that put NaN, infinities, -0.0, denormals and the int64 extremes in `BYTE_STREAM_SPLIT` and `DELTA_BINARY_PACKED` columns, flat and nested. All 4,018 leaf columns were bit-identical to PyArrow 26. Every version skips row groups by their min/max statistics on `BYTE_STREAM_SPLIT` floats, `DELTA_BINARY_PACKED` ids, dictionary strings and struct fields. None uses the page index, which takes 0.06–0.49% of each real file over 1 MB; other readers may use it. Several relayouted files of one view need `union_by_name = true`, and DuckDB 1.4.4 or later when their structs differ (1.1.3 fails; see [Struct lists](#struct-lists)).

### Comparing two conversions

Converters write rows in a fixed order, but two conversions of the same input do not produce identical bytes: the footer holds a uuid and a date, and `dataset.parquet` a packaging time. Compare rows and ids, not checksums.

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
