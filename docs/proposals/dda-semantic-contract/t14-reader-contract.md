# T14. Reader contract: how consumers must interpret DDA files

[← T13](t13-metadata-conformance.md) · [Overview](index.md) · [Strict profile](strict-profile.md) · [Annex overview](annex-overview.md)

!!! note "Annex: producer-neutral analysis"
    Part of the annex of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text, and they cover every producer. For OpenMS-written files the [strict profile](strict-profile.md) is the contract and wins where the two differ. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** A normative `consumer-guide.md` merges all topic reader rules into one default profile, implemented by qpx's read API: decoys and unidentified rows out; accepted rank-1 rows identify; transfers are never evidence and never combined with a measured own-run entry of the same precursor, run and label; signals and fractions are summed; the group is the protein; samples resolve exactly via `run.samples`; absent keys, null declaring columns and missing provenance parameters mean unknown; filters never silently return empty or unfiltered results; outputs record the policy applied.

## Questions a reader has
- Which rows drop by default: decoys, null `is_decoy`, rank > 1, failing, unidentified, transferred, contaminants?
- What does q ≤ t return without a q-value at that level?
- If a direct and a transferred row share precursor, run and label, which counts?
- Rows, fractions, charges on one key: sum, max or error?
- Is 0 missing? Is null 0? Is an imputed value measured?
- Group or anchor? What does `unique` mean?
- `run.samples` or SDRF? Normalised labels?
- How are spectra, PSMs, peptides, groups, spectral counts counted?
- What does an absent footer key mean?

## What happens today
- Readers MUST filter decoys (docs/spec/pg.md:99-103). The group, never its anchor, is the analyte (pg.md:66-67, data-model.md:159).
- Labels match exactly (data-model.md:239-246), yet views.md:96-97 accepts any label of a single-sample run.
- The protein is the group in views.md:153 and the anchor in views.md:212. The code uses the group (qpx/core/dataset_queries.py:90) [[CRITIC-66](findings-pg.md#critic-66)].
- anndata.md:153 documents summing a MaxQuant MS/MS feature with its MBR copy [[VALIDATOR-19](findings-feature.md#validator-19)].
- scores.md:38 calls `confident()` "the portable way to filter".
- intensities.md:121-133 defines null, 0 and imputed values, but sets no reader rule.
- The footer recipes (psm.md:275, feature.md:291, serialization.md:343, versioning.md:42), `file_metadata` (qpx/core/data/base.py:220-222 → qpx/core/parquet_io.py:11-15, first shard only) and qpx/cli/info.py:352 read `schema_arrow.metadata`. They miss keys added after the first row group (qpx/writers/base.py:359-365), e.g. `identity_source` (:462-463) [[INFRA-21](findings-cross.md#infra-21), [CRITIC-61](findings-cross.md#critic-61)].
- No page says what an absent key means.

| Producer | Behaviour | Findings |
|---|---|---|
| MaxQuant, FragPipe | Filtered upstream, but with no typed q-value, a null `pass_threshold` and no filter key. `feature.confident()` keeps 0 of 243,821 MaxQuant rows | [[MAXQUANT-05](findings-psm.md#maxquant-05), -26, [FRAGPIPE-02](findings-psm.md#fragpipe-02), [VALIDATOR-07](findings-feature.md#validator-07)] |
| mzIdentML (Mascot) | 11,517 rows have `pass_threshold = true` and a null `qvalue`. `psm.confident()` keeps 0 of 13,386 | [[MZID_CDAP-M1](findings-psm.md#mzid_cdap-m1), [VALIDATOR-06](findings-psm.md#validator-06)] |
| MaxQuant | 9,272 target (peptidoform, charge, run) keys (9,337 with decoys) have both a measured direct row and a measured transferred row. In 799 of them the transferred value is larger | [[MOKUME-09](findings-feature.md#mokume-09)] |

The qpx and mokume reader defects are listed under Reader rules [[INFRA-M1](findings-cross.md#infra-m1), [VALIDATOR-17](findings-feature.md#validator-17), -18, -19, -25, -M3, [MOKUME-01](findings-feature.md#mokume-01), -04, -06, -07, -09, -12, -15, -25, -29].

## Why it matters
Producer declarations ([T01](t01-psm-rows.md)–[T13](t13-metadata-conformance.md)) pay off only if readers act on them identically. Today the same call returns an empty set, a PSM-level cut or a precursor-level cut, depending on the producer. Readers also disagree on SUM vs MAX, on transfers, and on null flags. The failures are silent: empty results, MBR-inflated cells, and collapsed channels.

## Options considered
| Option | For | Against |
|---|---|---|
| A. Rules stay spread over topics | No new page | No checklist; defaults unreconciled |
| B. Strict defaults (q ≤ 0.01, no transfers, no contaminants) | Conservative | Drops MBR quantities the producer's protein values include; unusable without q-values; hides rows [#241](https://github.com/bigbio/qpx/issues/241) keeps |
| C. Pass-through (decoys only) | Simple | Keeps double and rank-inflated counts |
| **D. Guide + default profile + tri-state filters + policy record** | One meaning per call; testable | Changes documented API results |
| E. Refuse files with unknown keys | Forces declaration | 1.3 files unreadable |

D, with warnings and upper bounds from E and opt-in filters from B.

## Recommendation
**RDR-1 (MUST)** `docs/spec/consumer-guide.md` is the normative reader contract, added to [T13](t13-metadata-conformance.md) META-12 (producer-guide.md stays informative). It lists every topic reader rule (the `<prefix>-R<n>` rules of [T01](t01-psm-rows.md)–[T13](t13-metadata-conformance.md)) and RDR-2–14; a topic rule wins over the guide. qpx/core/data, qpx/views, qpx/core/dataset_queries.py, qpx/mudata.py and `qpxc transform` conform.

**RDR-2 (MUST)** Footer keys are read per META-R1–R3: key-value metadata, per file, shards combined by META-4; absent or off-vocabulary means `unknown` (off-vocabulary warns); never inferred from data (FLT-R1). `file_metadata` returns every registry key, `unknown` filled in. Declarations that are not footer keys come from their columns (`rank`, `precursor_index`, `identification_type`, `identified_in_unit`, `imputed_labels`, `is_imputed`), from `identity_composite`, or from provenance step parameters of the row's source (META-4); the table at the end of Reader rules says what a missing one means.

**RDR-3 (MUST)** A threshold filter (q-value at an explicit level, `pass_threshold`, `match_qvalue`) classes rows as pass, fail or unassessable: null, and not settled by `fdr_filter` or a provenance-stated `pass_threshold` basis (FLT-R3, CNF-R1). It returns the pass rows, reports the unassessable count and the keys used, and raises an error if every row is unassessable. When a q-value column is null on every row and `fdr_filter` does not settle it, FDR was not estimated at that level ([T04](t04-confidence-fields.md) CNF-2) and filtering at that level is not available: the error says so and names the levels that have q-values, and the reader filters at a level that exists. It never reads a null as 1.0 and never drops every row silently. On a transferred row, a column [T04](t04-confidence-fields.md) CNF-8 leaves null is not unassessable: the row passes through to the transfer policy and is reported apart (CNF-R5, MBR-R5). A null flag follows its owner rule: `rank` → PSM-R3 (upper bounds), `is_decoy` → PRT-R2 (invalid file), `contaminant` → PRT-R3 (unknown), `identification_type` → MBR-R2.

**RDR-4 (MUST)** The default profile is below. Caller overrides are recorded (RDR-13).

| Aspect | Identification counts, evidence | Quantities | Source |
|---|---|---|---|
| `is_decoy` true | excluded | excluded | FLT-R2 |
| `is_decoy` null | error (lenient: excluded and reported) | same | PRT-R2 |
| Rank | `rank = 1` (1.3: the score); none → upper bound | – | PSM-R1, -R3 |
| Confidence | FLT-R4; unknown state → non-decoy rank 1, as an upper bound | only on request, at an explicit level | CNF-R4 |
| `transferred` | excluded | kept, de-duplicated (RDR-5), masked | MBR-R3–R5 |
| Null type | own-run only when provenance records `match_between_runs = false` or `transferred_rows_removed = true`; otherwise unknown, count reported | kept | MBR-R2 |
| Unidentified: `identification_type = 'unidentified'`, tested as `peptidoform IS NULL` (below 1.4 also a `''` or `'None'` peptidoform, or charge ≤ 0) | excluded | excluded | UNI-R1, -R2 |
| Null entry | identified, not quantified | no cell, never 0 | UNI-R5 |
| Contaminant | included; true and null counts reported | kept, with a `contaminant` column | PRT-R3 |
| Imputed | – | excluded or masked | INT-R3 |
| 0 | – | measured (MAY be read as censored, recorded) | INT-R2 |
| Protein | `pg_group_id` or sorted membership | same | PRT-R1 |

Whether to use transferred quantities is the reader's choice. The default keeps them, because they improve quantification; filtering them out is optional (MBR-R5).

In 1.3 psm and feature files, contamination comes from the linked pg row; otherwise it is unknown.

**RDR-5 (MUST)** Combining entries of one peptidoform, charge, run (or unit) and label applies [T08](t08-match-between-runs.md) MBR-R4 first: a transferred entry is dropped beside a measured `direct`/`secondary` entry and kept beside an unquantified one; it is never summed, averaged or MAXed with a measured own-run entry. This replaces anndata.md:153.

**RDR-6 (MUST)** Aggregation order:
1. Rows of one (peptidoform, charge, run, label) are summed unless `identity_composite` gives the `precursor_run` or `precursor_unit` grain ([T05](t05-feature-rows.md) FRW-3, FRW-R1/R2); under those grains a repeated key is an error. This happens after RDR-5, and for reporter rows after CHI-R6's `quant_scan` de-duplication wherever `quant_scan` is set.
2. Fractions of a unit are summed. A multi-run `grouped_runs` row belongs to its unit (FRW-R3/R4). pg follows PGV-R4.
3. Roll-ups over charges, peptidoforms and proteins use the consumer's method, which is recorded (qpx sums). Reporter roll-ups across peptidoforms follow CHI-R6.
4. Units or replicates of one sample combine only by a named aggregation (`design_matrix` raises, qpx/core/dataset_queries.py:538-558).

Picking one row (MAX, best row) in steps 1–2 is a selection, not the quantity. It is allowed only as a named, recorded option.

**RDR-7 (MUST)** The protein is the group: `pg_group_id` or sorted distinct `pg_accessions`. `anchor_protein` is for display; an anchor-keyed output declares `var_key = anchor_protein` (PRT-R1). Group-specific = `group_unique`, single-protein group = `unique`, proteotypic = `pepmap.is_unique`; in 1.3 files group specificity is unknown (PRT-R4). Null `pg_accessions` is unattributed; a pepmap attribution is a declared inference (PRT-R5). Protein FDR comes from `pg.global_qvalue`, `pg.pg_qvalue` or pg `pass_threshold` with its basis (PGV-R3). `feature.pg_global_qvalue` equals `pg.global_qvalue` (CNF-1).

**RDR-8 (MUST)** A (`run_file_name`, label) pair resolves through `run.samples` by exact equality (data-model.md:239-246; SPC-R4).
- No lone-sample fallback, case folding, positional channel aliases or name parsing. Use `samples[].mixture`, `channel_role`, the replicate fields and `run.fraction`; take conditions from sample.parquet `factor_*` (sdrf-mapping.md:49).
- A pair matching no sample, or several, is excluded and reported (SILAC and LFQ alike).
- An SDRF matches runs via `run_file_stem` of `comment[data file]` and MAY add columns. Contradicting `run.samples` is an error unless the user overrides it, and the override is recorded.
- Without a run view, the shipped `{prefix}.sdrf.tsv` (file-naming.md:71) stands in. Without either, LFQ MAY use sample = run (recorded), and multiplexed data raises.
- `qpxc migrate` fixes legacy `raw` and `precursor_area` labels (qpx/transforms/migrate.py:41).

**RDR-9 (MUST)** Counts run over accepted rows (FLT-R4):
- **Spectra**: distinct (`run_file_name`, `scan`), or `spectrum_reference` where `scan` is empty (CHI-R1).
- **PSMs**: per PSM-R2 and CHI-R2.
- **Peptidoforms and peptides**: distinct, over accepted PSMs or `direct`/`secondary` features; never `''` or null.
- **Quantified features**: target, identified rows with a non-null entry (FRW-R6); transfers reported apart.
- **Groups per unit**: target, accepted pg rows with `identified_in_unit = true`; a null value gives an upper bound (PGV-R1).
- **Group spectral count**: distinct spectra among accepted PSMs whose linked feature carries the group, complete in 1.4 files except in runs listed in provenance `skipped_runs_psm` (FID-9, FID-R2). Counts via `protein_accessions` or pepmap count shared spectra in every group and say so. Tool counts (`spectral_count`, `unique_spectral_count`, `total_spectral_count`; pg.md:233-235) stay the tool's. Summing feature `spectral_count` over co-identified rows overcounts (CHI-R5).

Under an `unknown` filter state, counts are upper bounds (FLT-R5).

**RDR-10 (MUST)** Null is not 0, and an unstored MuData cell is not 0 (anndata.md:101). A stored 0 is measured. Imputed values stay flagged. Values the reader fills in appear only beside a mask.

**RDR-11 (MUST)** Readers check `qpx_version` and `file_type` (META-R4). Below 1.4, every 1.4 key is `unknown` (META-2), and the legacy clauses apply:
- PSM-R1: the score `rank`;
- [T02](t02-chimeric-spectra.md): `direct` may include secondary-only rows;
- CNF-R6, FRW-R7, FID-R7, MBR-R7 and SPC-R5;
- UNI-R2: a feature row with a null, `''` or `'None'` peptidoform or a charge ≤ 0 is unidentified;
- PRT-R4: below 1.1.3, a null `unique` with one member means unique.

**RDR-12 (MUST)** K (PSM-R4), the feature key under the `precursor_run`/`precursor_unit` grain of `identity_composite` (FRW-R1) and `pg_id` are unique only after strict validation. Otherwise a reader drops exact duplicates, reports keys whose rows differ, and never silently joins on or sums them.

**RDR-13 (SHOULD)** Derived outputs record a JSON `qpx_read_policy`: in `uns` of `.h5mu`, `.ae.h5ad` and `.de.h5ad`, and in provenance or a sidecar for tables such as mokume's TSVs. Fields:
- `decoys`, `rank`, `contaminants`, `imputed`, `zero`;
- `confidence{basis, level, threshold}`;
- `transferred` (`keep_dedup`, `drop` or `match_qvalue<=t`);
- `aggregation{within_key, fractions, units}`;
- `protein_key`, `sample_source`;
- `unassessable`, `contract_version`.

Source keys follow META-R7. `qpx_include_decoys` and `qpx_duplicates` remain as aliases.

**RDR-14 (SHOULD)** `tests/examples/reader_contract/` holds tiny datasets with expected counts and matrices. qpx CI runs them, and mokume can too. The cases are:
- direct and transferred entries on one key, with the direct entry quantified and unquantified;
- a null q-value under `fdr_filter` `psm:0.01`;
- rank-2 rows and tied rank-1 rows;
- labels that differ only in case;
- a 0 and an imputed value;
- groups sharing an anchor;
- a multi-run unit;
- a late footer key.

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| No new Parquet footer keys; `consumer-guide.md`; recipes read the KV footer; views.md:212 → group | docs | – | – | No |
| `qpx_read_policy` | `uns`; table sidecar | JSON | SHOULD | No |
| `transferred_mask` layer (pg: `identified_in_unit = false`) | MuData | sparse 0/1 | SHOULD with transfers | No |
| `contaminant` column | MuData var, `intensity()` | bool, null | SHOULD | No |
| `file_metadata` reads every shard's KV footer | qpx API | – | – | No |
| `confident()` per CNF-R4/RDR-3; `Feature.targets_only()`; `transferred=`, `contaminants=` | qpx API | – | – | **Yes**: signature, results |
| Drop the views.md:96-97 fallback; RDR-5 replaces the anndata.md:153 and `ds.intensity` sums | API, docs | – | – | **Yes**: legacy-label, MBR results |
| RDR-9 counts; MuData obs keyed on the unit (FRW-R4); profile in `for_quantification` | qpx API | – | – | **Yes**: counts, obs names, rows |
| A null `is_decoy` raises | MuData, quantify, site_quant | – | – | **Yes**, for invalid files only |

## Producer changes
| Producer | Change needed |
|---|---|
| qpx writers | Keys visible in `ARROW:schema` (META-3, FLT-13): `_add_footer_metadata` (qpx/writers/base.py:359-365) and the close-time `identity_source` stamp (:462-463) re-serialize it, like `_set_scan_format` (:367-383) |
| All converters | Write the [T03](t03-filtering-state.md), [T07](t07-unidentified-features.md), [T09](t09-intensities-missing-values.md), [T12](t12-protein-groups.md) and [T13](t13-metadata-conformance.md) footer keys, and the columns and provenance parameters [T01](t01-psm-rows.md)–[T12](t12-protein-groups.md) use instead of keys; until then default counts are upper bounds |
| MaxQuant, FragPipe, mzIdentML | Declare acceptance (`pass_threshold`, `fdr_filter`; FLT-9); `applied` (e.g. `psm:applied`) when the threshold is unknown |
| `qpxc migrate` | Write only the keys that META-4 lists |

## Reader rules
RDR-1–14 above are the normative reader rules. Current violations:

| Rule | Where | Violation | Fix |
|---|---|---|---|
| RDR-2 | qpx/core/parquet_io.py:11-15 (`file_metadata`, qpx/core/data/base.py:220-222); qpx/cli/info.py:352; recipes | Arrow metadata, first shard only | KV, every shard |
| RDR-3 | qpx/core/data/psm.py:23-25 | `qvalue <= t` only (0/13,386 Mascot); no rank | CNF-R4; raise |
| RDR-3 | qpx/core/data/feature.py:52-54 | COALESCEs levels; drops transfers | explicit level; CNF-R5 |
| RDR-3 | qpx/views/api.py:45-48 | Null `global_qvalue` passes | unassessable |
| RDR-4 | qpx/core/data/feature.py:91-115 | `for_quantification` keeps decoys, transfers and null entries; no `targets_only()` | profile |
| RDR-4 | qpx/mudata.py:279-280; qpx/cli/transform/quantify.py:305; qpx/transforms/site_quant.py:25 | Null `is_decoy` read as target | raise |
| RDR-4 | qpx views, MuData | `contaminant` never read | column |
| RDR-5 | qpx/core/dataset_queries.py:124-137; qpx/mudata.py:691; quantify.py:261-272 | Direct + transferred summed (12,672 cells, MaxQuant full) | MBR-R4; mask |
| RDR-6 | qpx/mudata.py:178-189 | obs keyed on `run_file_name` | unit |
| RDR-7 | qpx/core/data/feature.py:33-41; qpx/views/api.py:80-85 | Anchor-keyed helpers and export | `pg_group_id` |
| RDR-8 | qpx/core/sample_sql.py:79-84; qpx/core/mudata_obs.py:92-106 | Lone-sample and stem fallback [[VALIDATOR-17](findings-feature.md#validator-17)] | exact; report |
| RDR-8 | quantify.py:315-321 | External SDRF, else `run::label` | `run.samples` |
| RDR-9 | qpx/views/api.py:139-151 | `ModificationView` counts every rank and threshold: Deamidated 4,323 vs 926 accepted [[VALIDATOR-M3](findings-psm.md#validator-m3)] | accepted rows |
| RDR-9 | qpx/views/api.py:21-29, 172-183, 247-268; qpx/converters/summary.py:109-131 | Count transfers, ID-only rows, rank > 1, every pg row (via `grouped_runs[1]`), decoy accessions | RDR-9, or label as row counts |
| RDR-10 | quantify.py:269; qpx/core/data/schema.py:306; qpx/core/dataset_queries.py:124-137 | Zeros dropped or read as "no quantity"; imputed values unflagged | keep zeros; mask |
| RDR-2/3 | mokume:python/mokume/io/feature.py:185-197 | No filter key read, though `parquet_kv_metadata` is used (mokume:python/mokume/io/_qpx.py:35-38) | read `fdr_filter` and the [T03](t03-filtering-state.md)/T08 provenance parameters |
| RDR-3 | mokume:python/mokume/io/_qpx.py:13-46; mokume:python/mokume/preprocessing/filters/peptide.py:403-409; mokume:rust/crates/mokume-pipeline/src/lib.rs:6966-6972 | Typed column used only when it has no nulls; nulls dropped silently; score filters drop transfers [[MOKUME-21](findings-feature.md#mokume-21), -22, -M3] | tri-state |
| RDR-4 | mokume:python/mokume/core/contaminants.py:8; mokume:python/mokume/pipeline/stages.py:436-438; mokume:rust/crates/mokume-pipeline/src/lib.rs:8461, 8902-8909 | Substring matching; null decoy read as target | typed flags |
| RDR-5 | mokume:python/mokume/io/feature.py:560-643; mokume:rust/crates/mokume-io/src/qpx.rs:493-524; mokume:rust/crates/mokume-core/src/config.rs:172-180 | Type ignored (Python) or inferred from `id_run_file_name` (Rust); `keep` unrecorded | read; de-duplicate; record |
| RDR-5/6 | mokume:python/mokume/preprocessing/aggregation.py:263-322; mokume:rust/crates/mokume-pipeline/src/lib.rs:2761-2772 | MAX over rows, transfers, fractions | RDR-5, then sum |
| RDR-6 | mokume:python/mokume/io/feature.py:679-686 | Unit never read [[MOKUME-11](findings-feature.md#mokume-11)] | unit |
| RDR-7 | mokume:python/mokume/io/feature.py:225-230, 295, 1018-1030; mokume:python/mokume/preprocessing/aggregation.py:167; mokume:python/mokume/pipeline/stages.py:1771-1776 | `unique` as group-specific; null membership dropped; order and anchor as keys | `group_unique`; `pg_group_id` |
| RDR-8 | mokume:python/mokume/io/feature.py:44-79, 605-643, 743-746; mokume:rust/crates/mokume-io/src/sdrf.rs:108-146 | `run.samples` ignored; lowercasing; label as run; mixture from names; TMT collapses without SDRF [[MOKUME-12](findings-cross.md#mokume-12), -13, -14, -28] | RDR-8 |
| RDR-9 | mokume:rust/crates/mokume-io/src/psm.rs:14-25; mokume:rust/crates/mokume-pipeline/src/spectral_count.rs:158-208 | No rank or threshold; counts rows; no label | RDR-9 |
| RDR-4/7/8 | mokume:python/mokume/tissuemap/loader.py:51, 56, 203-262, 360 | Anchor key; null decoy kept; positional aliases; TMT126/131 fallback; first sample | `channel_role` |
| INT-R4 | mokume:rust/crates/mokume-pipeline/src/tmt.rs:88-114 | `--tmt-impurity` ignores the input's correction state [[MOKUME-26](findings-feature.md#mokume-26)] | check `intensity_quantity` (`reporter_corrected`) |

**When a declaration is absent, null or `unknown`, a reader MUST:**

| Declaration (owner) | Behaviour |
|---|---|
| `fdr_filter` ([T03](t03-filtering-state.md)) | Treat a null typed value as unassessable (FLT-R3); report counts as upper bounds (FLT-R5) |
| `unquantified_policy` ([T07](t07-unidentified-features.md), [T12](t12-protein-groups.md)) | Never read a missing row as "unquantified" or "not identified" (FRW-R6, PGV-R2) |
| `intensity_quantity`, `intensity_processing` ([T09](t09-intensities-missing-values.md)) | Warn before pooling files (INT-R1); never re-correct silently, and warn on `reporter` (INT-R4) |
| `scan_format` ([T11](t11-spectrum-precursor-fields.md)) | Up to 1.3, absent means `scan` (SPC-7). From 1.4, `unknown` means join on `spectrum_reference`, or on `scan` only if `rt` agrees (SPC-R2) |
| `acquisition_method` ([T13](t13-metadata-conformance.md)) | Use `run.acquisition_method`. Otherwise apply DDA rules only if the caller asserts DDA, and record it (META-R5) |
| `rank` column ([T01](t01-psm-rows.md) PSM-4) | Count `rank = 1`; with a null rank, counts are upper bounds (PSM-R3) |
| `precursor_index`, `quant_scan` ([T02](t02-chimeric-spectra.md) CHI-10, CHI-11) | Read a null index as unknown, never 1 (CHI-R7); count spectra by scan (CHI-R1); de-duplicate reporter roll-ups on `quant_scan` wherever it is set (CHI-R6) |
| `identity_composite` ([T05](t05-feature-rows.md) FRW-3) | Derive the grain; a missing composite is the `peak` default. Apply RDR-5, then sum rows sharing a natural key, except under `precursor_run` or `precursor_unit` (FRW-R1/R2) |
| PSM links, provenance `skipped_runs_psm` ([T06](t06-identifications-per-feature.md) FID-9) | Treat links in 1.4 files as complete except in listed runs; there, never read a missing link as "not supporting" (FID-R2) |
| Unidentified rows in the feature file, provenance `skipped_unidentified_features` ([T07](t07-unidentified-features.md) UNI-5) | Compute an identification rate only when the file holds unidentified rows and none were skipped (UNI-R3) |
| `identification_type`, `identified_in_unit`, provenance `match_between_runs` step ([T08](t08-match-between-runs.md) MBR-8, MBR-9) | Report null types as unknown unless provenance settles them (MBR-R2); a null `identified_in_unit` gives upper bounds (PGV-R1); when `mbr_fdr` is `none`, `unknown` or absent, warn and offer to drop transfers (MBR-R6) |
| `imputed_labels`, `is_imputed` ([T09](t09-intensities-missing-values.md) INT-9) | Read a null flag as unknown, not measured: exclude or mask (INT-R3) |
| Provenance `pass_threshold_level`, `pass_threshold_fdr` ([T03](t03-filtering-state.md) FLT-8) | Treat `pass_threshold = true` as an FDR cut at t only if they state an FDR level and a value ≤ t (FLT-R4) |
| Provenance `decoy_policy`, `contaminant_policy` ([T03](t03-filtering-state.md) FLT-10, FLT-11) | Exclude `is_decoy = true` anyway, and filter contaminants by the typed flag (FLT-R2, FLT-R6); estimate FDR only under `decoy_policy = retained` |
| Provenance `include_unassigned_psms` ([T03](t03-filtering-state.md) FLT-12) | Under `false`, never read a PSM missing from the psm view as "not identified" |
| Provenance `psm_qvalue_scope`, `peptide_qvalue_scope` ([T04](t04-confidence-fields.md) CNF-3) | Report no experiment-wide FDR (CNF-R3) |
| Provenance `observed_mz_calibration` ([T11](t11-spectrum-precursor-fields.md) SPC-2) | Treat `unknown` as a calibration state of its own: never pool its mass errors with `raw` or `recalibrated` files (SPC-R3) |
| Provenance `pg_grouping` ([T12](t12-protein-groups.md) PGV-6) | Pivot `pg_group_id` across units only under `dataset_inference` or `mapping_set`; otherwise warn (PGV-R5) |

## Open questions for maintainers
1. *Resolved by the maintainer:* whether to use transferred quantities is the reader's choice. The default keeps them, de-duplicated and masked, because they improve quantification and MaxQuant and IonQuant protein values include MBR. Filtering them out is optional (RDR-4, [T08](t08-match-between-runs.md) MBR-R5; mokume#117).
2. Should matrices drop contaminants by default, as the usual Perseus workflow does?
3. Should the lone-sample fallback survive one release with a deprecation warning?
4. Should FRW-R2/R3 allow a recorded best-row or MAX selection (mokume TMT `best-row`), as RDR-6 does?
