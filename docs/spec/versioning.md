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

All serialized QPX views (PSM, Feature, PG, MZ, and the other Parquet and AnnData files) include a `qpx_version` field in their file-level metadata identifying which version of the QPX specification generated the file. In Parquet files it is a key-value pair in the file footer; in AnnData (`.h5ad`) and MuData (`.h5mu`) files it is stored under `uns`. Files written by the qpx 1.0.x releases carry the package version instead, such as `1.0.2`; readers compare the first two components only, so such a file is specification 1.0. Peptide- and protein-level summaries are derived [API views](views.md) computed on demand from these files, not standalone serialized files.

```python
import pyarrow.parquet as pq

# Writing a file with version metadata
metadata = {"qpx_version": "1.3"}
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

The current QPX specification version is **1.3**.

This version defines all core serialized views (PSM, Feature, PG, MZ), the derived [API views](views.md) (Peptide, Protein) computed on demand from them, expression views (Absolute, Differential), and metadata views (SDRF, Project).

### Changelog

- **1.3** (backward-incompatible for identities):
    - **Canonical peptidoform.** Every producer writes the
      [canonical form](peptidoform.md#canonical-form): fixed modifications and
      labels inline, `UNIMOD:n` accessions in upper case, and signed
      four-decimal mass shifts instead of `CHEMMOD:`. `calculated_mz` is the
      theoretical m/z of that peptidoform. `feature_id` and `psm_id` hash the
      peptidoform, so they change wherever a 1.2 producer left fixed
      modifications out (MaxQuant) or spelled a modification differently.
      A mass shift of unknown position is written `[+129.0804]?PEPTIDE`.
      Re-convert from the source files.
    - **Missing values.** An `intensities` entry may have a null `intensity`:
      the feature was identified in that label but not quantified
      ([Intensities](intensities.md#missing-values)). `psm.observed_mz` is
      nullable. `feature.unique` has one definition: a peptide of a one-protein
      group that maps to no other group.
    - **MS1 labels.** SILAC, dimethyl and mTRAQ multiplets have a normative
      encoding ([MS1 labels](intensities.md#ms1-labels-silac-dimethyl-mtraq)):
      one feature row per multiplet, with the lightest partner's peptidoform
      and one `intensities` entry per channel label, the optional
      `feature.identified_labels` and `psm.label`, labelled PSM peptidoforms,
      and `ratio_<num>_<den>` ratios on the numerator label.
    - **Value ranges.** The YAML schemas declare `minimum` and `maximum` for
      q-values and posterior error probabilities (0 to 1), sequence coverage
      (0 to 100) and primary intensities (at least 0); NaN and infinite values
      are outside every range.
    - **Identity composites.** One identity composite value names one row
      identity, also for producer-supplied ids.
- **1.2** (backward-compatible): adds the [site](site.md) and
  [proteoform](proteoform.md) views, the optional `modification_candidates` column
  of the PSM and Feature views ([Modifications](modifications.md)) and the optional
  `sequence_edits` column of the pepmap view ([Sequence Edits](sequence-edits.md)).
  Run, MZ and PSM files may leave out the columns the specification does not
  require, and `mz.ion_injection_time` may be null. Nothing is removed or retyped,
  so 1.1 files are valid 1.2 files.
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
      `qpxc validate` (the audit/CI path, which always validates strictly) and
      `Dataset.validate(strict=True)` still report a duplicate primary key as an
      **error**; `Dataset.validate()` reports it as a warning.
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

In addition to the specification version, every QPX file records the software behind it in two footer fields: `software_provider` names the software that wrote the file (`qpx` for the qpx writers, which also stamp their package version as `writer_version`), and `creator` names the tool whose output the data was converted from (`maxquant`, `diann`, `quantms`, ...):

```json
{
  "software_provider": "qpx",
  "writer_version": "1.2.0",
  "creator": "diann"
}
```

!!! info "Why track the software provider?"
    Tracking the generating software enables reproducibility and debugging. If a downstream analysis produces unexpected results, the software provider metadata allows users to determine which tool version created the input data and whether known issues apply.

Both fields are free-text strings in the Parquet metadata.

## Version Compatibility Matrix

| Reader Version | File Version | Compatible? | Notes |
|---------------|-------------|-------------|-------|
| 1.0 | 1.0 | Yes | Exact match |
| 1.1 | 1.1 | Yes | Exact match |
| 1.1 | 1.0 | No\* | 1.0 core files lack the mandatory identity-ID columns, and 1.0 PG uses scalar `run_file_name`; re-convert the dataset. |
| 1.0 | 1.1 | No | 1.1 adds mandatory identity-ID columns and replaces PG `run_file_name` with `grouped_runs`; a strict 1.0 reader cannot read it. |
| 1.3 | 1.2 | Yes\* | Same columns; identities of rewritten peptidoforms differ, so do not merge 1.2 and 1.3 files of one dataset |
| 1.2 | 1.3 | Yes\* | Same columns; see the 1.3 changelog for value changes |
| 1.2 | 1.1 | Yes | 1.2 only adds views and optional columns |
| 1.1 | 1.2 | Partly | A 1.1 reader ignores the new views and optional columns, but 1.1 validation reports run, MZ and PSM files that leave out columns 1.1 requires. |
| 2.0 | 1.0–1.3 | Maybe | Major version change; check migration guide |
| 1.x | 2.0 | No | Reader cannot handle a newer major version |

\* The 1.x line carries backward-incompatible changes during pre-2.0 stabilisation (see the Rules above), so even within `1.x` a reader must check `qpx_version` and the changelog — it is not safe to assume minor updates are compatible.

## See Also

- [Serialization](serialization.md) -- serialization formats and file-level metadata details
- [QPX Format Overview](index.md) -- overview of all QPX views and the current specification status
