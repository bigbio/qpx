# Versioning

QPX uses a simple versioning scheme to track specification changes and ensure that tools can correctly interpret files produced by different versions of the format.

## Version Format

QPX versions follow a **`{major}.{minor}`** format -- for example, **1.0**.

This is a two-component scheme without patch numbers. The version applies to the entire QPX specification, not to individual views.

## Rules

- **Major releases** (`X.0`): Reserved for large-scale, backward-incompatible restructuring of the data model once the format is declared stable.
- **Minor updates** (`1.X`): The active line of development while the format stabilises. Minor updates normally add backward-compatible fields, views, or metadata keys, but **may also carry backward-incompatible changes** (removing or renaming a field, changing a field type) during the pre-2.0 stabilisation period. Every such change is called out in the [changelog](#changelog) below, and tools should check the `qpx_version` in the file footer and consult the changelog rather than assuming minor updates are always compatible.

!!! warning "Pre-2.0 stabilisation"
    While QPX is on the `1.x` line the specification is still stabilising, so a
    minor release may change the data model in a backward-incompatible way. Pin
    the `qpx_version` you target and read the changelog before upgrading. Once the
    format is declared stable the strict "minor = additive only" rule takes over
    and breaking changes move to a major (`2.0`) release.

!!! note
    QPX does not use patch versions. Bug fixes to the specification text that do not change the data model are tracked in the documentation changelog but do not increment the version number.

## Version in Files

All serialized QPX views (PSM, Feature, PG, MZ, and the other Parquet and AnnData files) include a `qpx_version` field in their file-level metadata identifying which version of the QPX specification generated the file. In Parquet files it is a key-value pair in the file footer; in AnnData (`.h5ad`) and MuData (`.h5mu`) files it is stored under `uns`. Peptide- and protein-level summaries are derived [API views](views.md) computed on demand from these files, not standalone serialized files.

```python
import pyarrow.parquet as pq

# Writing a file with version metadata
metadata = {"qpx_version": "1.1"}
```

!!! tip "Reading the version from a file"
    ```python
    import pyarrow.parquet as pq

    parquet_file = pq.ParquetFile("experiment.psm.parquet")
    schema_metadata = parquet_file.schema_arrow.metadata
    qpx_version = schema_metadata.get(b"qpx_version", b"unknown").decode()
    print(f"QPX version: {qpx_version}")
    ```

## Current Version

The current QPX specification version is **1.1**.

This version defines all core serialized views (PSM, Feature, PG, MZ), the derived [API views](views.md) (Peptide, Protein) computed on demand from them, expression views (Absolute, Differential), and metadata views (SDRF, Project).

### Changelog

- **1.1.2** (on-disk spec version stays `1.1`): two related fixes to identity
  handling.
    - **PG identity now keys on full group membership** (backward-incompatible for
      `pg_id` values): the PG schema-default `identity_composite` changed from
      `[anchor_protein, grouped_runs, label]` to
      `[pg_accessions, grouped_runs, label]`. Keying on only the group **leader**
      (`anchor_protein`) gave two genuinely distinct protein groups that share a
      leading protein the **same** `pg_id`, and split a single group with
      inconsistent per-precursor annotations into duplicate-identity rows. The
      identity now hashes the full `pg_accessions` membership (order-independently,
      like `grouped_runs`), so distinct groups get distinct ids. `anchor_protein`
      remains a descriptive field.
    - **Duplicate primary keys on the write/convert path are now a warning, not a
      hard error.** While the format stabilises, writers and converters persist
      source data as-produced and log a warning on a duplicate identity id rather
      than crashing (a duplicate usually signals a converter identity bug worth
      investigating). A **NULL primary key remains a fatal error**, and
      `qpxc validate --strict` (the audit/CI path) still reports a duplicate
      primary key as an **error**.
    - **No in-place migration.** These changes affect derived `pg_id` values and
      identity handling; regenerate qpx files by **re-converting from the original
      source files** — there is no on-disk migration.
- **1.1** (backward-incompatible, pre-2.0 stabilisation): the PG view now keys on
  `grouped_runs` (`list<string>`) instead of the scalar `run_file_name`. A
  protein-group quantity applies to the set of raw files (fractions) aggregated
  into one quantification unit, not a single file; the sample is resolved via
  `(any file in grouped_runs, label) -> run.samples[]`. This removes and retypes
  a field, so files written by 1.1 are **not** readable by strict 1.0 tooling —
  shipped as a minor under the pre-2.0 stabilisation rule above. The PG
  `intensities` list (`list<struct<label,intensity>>`) is also **flattened** into
  scalar `label` + `intensity` columns — one row per label. `label` participates
  in the schema-default identity composite `[anchor_protein, grouped_runs,
  label]`; `label`/`intensity` are null for identification-only groups. This
  removes and retypes columns, so pg files must be regenerated under >= 1.1.
- **1.1** (backward-incompatible): Feature, PSM, and PG now use mandatory opaque
  `feature_id`, `psm_id`, and `pg_id` columns as their primary keys. A producer ID
  is preserved when supplied; otherwise QPX derives the ID from the
  footer-declared `identity_composite`. The schema-default Feature composite
  changes from `[sequence, charge, run_file_name, anchor_protein]` to
  `[peptidoform, charge, run_file_name, rt]`, and the PSM composite changes from
  `[sequence, charge, run_file_name, scan]` to `[peptidoform, charge,
  run_file_name, scan]`. Converters may declare a different producer-specific
  Feature composite. Regenerate Feature, PSM, and PG files under >= 1.1.
- **1.1** (OpenMS bridge): OpenMS and QPX share the *fraction_group* concept — a
  group of fraction raw files that together quantify one protein, which is exactly
  what `grouped_runs` encodes. Rather than change OpenMS's Arrow schema, the QPX
  OpenMS converter stamps each pg (and feature) row with a `fraction_group`
  cv_param carrying OpenMS's fraction_group number (read from the consensusXML
  experimental design), so the grouping is preserved across both sides via a
  shared cv term. No OpenMS schema change is required.
- **1.0**: initial specification.

## Software Provider

In addition to the specification version, every QPX file records the software that generated the data. This is stored in the `software_provider` metadata field in the Parquet file footer. The value identifies the tool name and version that produced the file:

```json
{
  "software_provider": {
    "name": "quantms",
    "version": "1.3.0"
  }
}
```

!!! info "Why track the software provider?"
    Tracking the generating software enables reproducibility and debugging. If a downstream analysis produces unexpected results, the software provider metadata allows users to determine which tool version created the input data and whether known issues apply.

The `software_provider` field is a free-text string in the Parquet metadata. The JSON structure shown above is the recommended format, but tools may also use a flat string such as `"QuantMS 1.3.0"`.

## Version Compatibility Matrix

| Reader Version | File Version | Compatible? | Notes |
|---------------|-------------|-------------|-------|
| 1.0 | 1.0 | Yes | Exact match |
| 1.1 | 1.1 | Yes | Exact match |
| 1.1 | 1.0 | No\* | 1.0 core files lack the mandatory identity-ID columns, and 1.0 PG uses scalar `run_file_name`; re-convert the dataset. |
| 1.0 | 1.1 | No | 1.1 adds mandatory identity-ID columns and replaces PG `run_file_name` with `grouped_runs`; a strict 1.0 reader cannot read it. |
| 2.0 | 1.0/1.1 | Maybe | Major version change; check migration guide |
| 1.x | 2.0 | No | Reader cannot handle a newer major version |

\* The 1.x line carries backward-incompatible changes during pre-2.0 stabilisation (see the Rules above), so even within `1.x` a reader must check `qpx_version` and the changelog — it is not safe to assume minor updates are compatible.

## See Also

- [Serialization](serialization.md) -- serialization formats and file-level metadata details
- [QPX Format Overview](index.md) -- overview of all QPX views and the current specification status
