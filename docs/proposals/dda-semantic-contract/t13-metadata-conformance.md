# T13. Metadata, validation and conformance: making the contract checkable

[← T12](t12-protein-groups.md) · [Overview](index.md) · [T14 →](t14-reader-contract.md)

!!! note "Proposal"
    Part of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** One registry (`footer.yaml`) defines the Parquet key-value footer as each file's declaration, adding exactly five keys (`acquisition_method`, `fdr_filter`, `unquantified_policy`, `intensity_quantity`, `intensity_processing`); other declarations are columns, fixed rules or provenance parameters; rewrites carry every key by its rule, merge keeps every input's provenance steps, and none relabels old rows with a new `qpx_version`; writers and `qpxc validate` share one check set with error, warning and info severities and `core`/`dda` profiles chosen by `--profile` (files claim none), so no allowed state fails `--conformance`; any change of meaning moves `qpx_version`.

## Questions a reader has
- Which footer keys exist, where are they read, and what does an absent or unrecognised value mean?
- Do merge, migrate, relayout and write-back keep the declarations and the version?
- Is this file DDA or DIA, and which rules apply?
- What is "valid" when writer, `Dataset.validate()` and `--conformance` disagree? Can a file claim the DDA contract?
- Can meanings change under one `qpx_version`? Which page wins on a contradiction?
- Do tests and fixtures check the contract?

## What happens today
**Spec**
- `serialization.md:283-296` lists 12 keys, two never written (`project_accession`, `project_title`); `uuid`, `identity_hash_version`, `identity_source`, `spectrum_content`, `acquisition_start`, `xic_*` and `qpx.*` live in prose (`identity.md:20-21`, `serialization.md:63,72,80,92`, `psm.md:290`, `mz.md:154,162`, `xic.md:21`); `psm.md:256-268` and `feature.md:272-284` repeat 7 keys [[INFRA-22](findings-cross.md#infra-22)][[SPEC-36](findings-cross.md#spec-36)][[CRITIC-61](findings-cross.md#critic-61)].
- The four reading recipes use `schema_arrow.metadata` (`psm.md:275`, `feature.md:291`, `serialization.md:343`, `versioning.md:42`).
- Acquisition lives only in the optional free-text `run.acquisition_method` from the SDRF (`run.md:125`; fixture value `NT=Data-Independent Acquisition;AC=NCIT:C161786`); the psm view is "primarily DDA" (`psm.md:282-283`) [[EXTERNAL-33](findings-cross.md#external-33)][[CRITIC-59](findings-cross.md#critic-59)]; `provenance.md:105` shows `acquisition_method = lfq`.
- Writers "validate on close()" (`producer-guide.md:60`); a duplicate PK is a write-time warning but a `qpxc validate` error (`versioning.md:120-127`); `--conformance` adds `--semantic` and `--check-identity` and fails on any warning (`qpx/cli/validate.py:112-113,125-126`).
- Only id derivation must move the version (`identity.md:69`); `pg_id` changed under one version (`versioning.md:131-133`) [[CRITIC-63](findings-cross.md#critic-63)]; minors may break before 2.0 (`versioning.md:14`).
- "Required: no" means "must be present" (`data-model.md:302-307`), yet view pages mark YAML-optional columns "no" [[SPEC-31](findings-cross.md#spec-31)].

| Component | Behaviour | Findings |
|---|---|---|
| Writers | `_add_footer_metadata` (`qpx/writers/base.py:359-365`) does not re-serialize `ARROW:schema`, which only `_set_scan_format` (`:367-383`) does, so `identity_source`, added at close (`:462-463`), is missing from `schema_arrow` in nearly every output. The table path (`:416`) persists range, vocabulary and sign violations as warnings; the record path checks only nulls, PK, pg references and `scan_format`. No converter calls `Dataset.validate`; only the golden dataset runs `--conformance` (`tests/conformance/test_conformance.py:20`) | [[INFRA-21](findings-cross.md#infra-21)][[VALIDATOR-M4](findings-cross.md#validator-m4)][[INFRA-12](findings-cross.md#infra-12)][[VALIDATOR-01](findings-cross.md#validator-01)] |
| Footer content | No key for filtering, rank, decoys, MBR, intensity kind or acquisition. MaxQuant, FragPipe, PD and CDAP write no `scan_format`; OpenMS and mzIdentML only for uniform references (`qpx/converters/mzidentml/psm_records.py:38-43`). Provenance lacks tool versions and MaxQuant sidecars | [[VALIDATOR-31](findings-cross.md#validator-31)][[INFRA-20](findings-cross.md#infra-20)][[MZID_CDAP-20](findings-cross.md#mzid_cdap-20)][[MAXQUANT-06](findings-cross.md#maxquant-06)][[FRAGPIPE-29](findings-cross.md#fragpipe-29)][[PD-31](findings-cross.md#pd-31)][[MAXQUANT-M3](findings-cross.md#maxquant-m3)] |
| Validator | `footer_issues` checks only version, file type and identity keys (`qpx/core/data/footer_validation.py:77-92`). `identity_not_checked` (`:112-114`) alone fails PXD054720 mzIdentML under `--conformance`, though its composites are unique (`:153` separates cross-links). `pg_link_coverage` warns on allowed null-membership features. Without a run view, run and label checks are skipped silently (`qpx/core/dataset_validation.py:193-194`). OpenMS fixtures fail on `q-value` (`qpx/converters/openms_consensus/psm_adapter.py:224`), MaxQuant simple on a dangling run | [[VALIDATOR-02](findings-psm.md#validator-02)][[VALIDATOR-03](findings-psm.md#validator-03)][[VALIDATOR-24](findings-feature.md#validator-24)][[VALIDATOR-29](findings-cross.md#validator-29)] |
| Rewrites | **merge** reads `pq.read_schema` (`qpx/core/dataset_merge.py:240`), drops `scan_format` and `spectrum_content`, writes `creator` `qpx` (`:252-267`) and checks no version. **migrate** reads `read_schema` (`qpx/transforms/migrate.py:54-55`), keeps `creator`, `scan_format`, the feature composite (`:150-160`), `qpx.mantissa_bits` and xic keys (`:195-197`), and would re-derive mzIdentML ids because `identity_source` is invisible (`:211`, `:81`; code reading). **write-back** keeps four keys (`qpx/core/dataset_writes.py:123-133`); the OpenMS `-out_qpx` upgrade and run annotation also read `schema_arrow` (META-R1). **relayout** reads the key-value footer (`qpx/transforms/relayout.py:136-143`). All stamp the running `qpx_version` (`qpx/writers/partitioned.py:60-67`) | [[VALIDATOR-30](findings-cross.md#validator-30)] |
| Docs | `converter-coverage.md:181,285,303`, `tool-mappings.md:17-30,105-109` and `ontology.md:244,290-291` contradict the code. `producer-guide.md:48` names a non-existent `imputation` flag. `views.md` contradicts itself, and `scores.md` contradicts `types.yaml` | [[CRITIC-M6](findings-cross.md#critic-m6)][[MOKUME-M4](findings-cross.md#mokume-m4)][[CRITIC-33](findings-feature.md#critic-33)][[INFRA-04](findings-feature.md#infra-04)][[FRAGPIPE-24](findings-cross.md#fragpipe-24)][[SPEC-33](findings-cross.md#spec-33)][[CRITIC-65](findings-cross.md#critic-65)][[SPEC-M6](findings-cross.md#spec-m6)][[SPEC-34](findings-cross.md#spec-34)][[VALIDATOR-26](findings-cross.md#validator-26)][[SPEC-30](findings-cross.md#spec-30)][[FRAGPIPE-M3](findings-cross.md#fragpipe-m3)] |
| Tests, fixtures | PD tests pin `unique=True` on a 2-member group, `normalized_abundance` and an identification-only `LFQ` pg row (`tests/converters/test_proteome_discoverer_converter.py:151-169,197-227`). `psm_empty_scan` fails only through a stale id. No case covers an allowed null-membership feature (`feature_no_group` keeps its anchor, so it rightly fails). MaxQuant fixtures hold MBR (1,469 / 104,491 `MULTI-MATCH` rows) and second peptides (30 msms / 5,110 evidence `MULTI-SECPEP` rows); PXD054720 has ranks 1-12, one engine, no decoys. No PD, FragPipe, two-engine, chimeric-search or IonQuant MBR export | [[PD-33](findings-cross.md#pd-33)][[VALIDATOR-13](findings-psm.md#validator-13)] |
| Readers | qpx reads `schema_arrow` of the first shard and stamps its own version on exports; mokume never checks `qpx_version`, `file_type` or acquisition (sites under Reader rules) | [[MOKUME-23](findings-cross.md#mokume-23)] |

## Why it matters
[T01](t01-psm-rows.md)-[T12](t12-protein-groups.md) move most semantics into footer keys, typed columns and provenance parameters (P-B), which work only if readers find them, rewrites keep them and a validator checks them. Today a reader who follows the spec misses `identity_source`; a merge reads zero-based indices as scan numbers; relayout or merge by a newer qpx labels 1.3 rows as 1.4; `--conformance` fails spec-sanctioned mzIdentML output; and nothing says whether DDA rules (one row per hypothesis, `scan` = the triggering MS2) apply.

## Options considered
| Option | For | Against |
|---|---|---|
| Provenance parameters | Exist (`provenance.md:65`) | Free text, unvalidated, per dataset not per file [[CRITIC-60](findings-cross.md#critic-60)] |
| Registry of footer keys | Travel with the file. Precedents: xic keys, mzTab MTD, mzIdentML `Threshold` | About 30 keys, each needing a carry rule |
| Readers use `schema_arrow` | Current recipes | Misses late keys; DuckDB and Rust read key-value pairs |
| `acquisition_method` only on `run` | No new key | Optional free text; needs an SDRF; absent from psm-only exports |
| Profiles `core`/`dda` plus an optional claim | Testable; DIA and `unknown` still pass `core` (P-G); mirrors mzTab `mzTab-mode` | A claim can be false |
| Meaning changes via `writer_version` | No spec churn | Dev builds (`0.1.dev877`) map to no meaning |

## Recommendation
**META-1 (MUST) One registry.**
- `qpx/core/data/schemas/footer.yaml` (mirrored in `docs/spec/schemas/`) defines each key's views, value type, first `qpx_version`, owner, carry rule and `dda` requirement. `serialization.md` renders it as the only footer table; delete `psm.md:256-278` and `feature.md:272-294`; drop `project_*`.
- Every 1.4 enum also allows `unknown`. Private keys use `x_<namespace>_<name>`. "The T13 profile" in [T01](t01-psm-rows.md)-[T12](t12-protein-groups.md) means `dda`.
- 1.4 adds exactly five keys: `acquisition_method`, `fdr_filter`, `unquantified_policy`, `intensity_quantity` and `intensity_processing` (the existing `scan_format` keeps its [T11](t11-spectrum-precursor-fields.md) SPC-7 extension). Every other [T01](t01-psm-rows.md)-[T12](t12-protein-groups.md) declaration is a column, a fixed rule or a provenance step parameter, defined by its owner; the registry lists none of them.

Blank Carry cell: kept when all sources agree, else `unknown`.

| Key | Views | Values (owner) | Carry | `dda` |
|---|---|---|---|---|
| `qpx_version`, `file_type` | all | `<major>.<minor>`; `<view>_file` | META-4 | required |
| `creator`, `software_provider`, `writer_version`, `creation_date`, `uuid`, `compression_format` | all | free text (`creator` = converter id); date; UUID; codec | `creator`: sorted list; rest new | `creator` |
| `primary_key`, `identity_composite`, `identity_hash_version`, `identity_source` | id views | identity.md | recomputed | required |
| `scan_format` | psm, feature, mz | [T11](t11-spectrum-precursor-fields.md) SPC-7 | agreed value; else `mixed` (psm, every row with `spectrum_reference`); else `unknown` | psm; feature with a `scan` |
| `spectrum_content`; `acquisition_start`; `qpx.*`; `xic_*`, `transition_format` | psm, mz; mz; any; xic | psm.md:290, mz.md:154; JSON; JSON; xic.md:21 | default; union; recomputed; dropped | psm with arrays |
| `acquisition_method` | psm, feature, pg | `dda`, `dia`, `prm`, `mixed` (T13) | differing known values → `mixed` | `dda`, `mixed` |
| `fdr_filter` | psm, feature, pg | `none`, or `<level>:<value>` entries joined by `;` ([T03](t03-filtering-state.md) FLT-1) | weakest value per level (FLT-13) | present |
| `unquantified_policy` | feature, pg | `retained`, `removed` ([T07](t07-unidentified-features.md) UNI-6, [T12](t12-protein-groups.md) PGV-4) | | present |
| `intensity_quantity`, `intensity_processing` | feature, pg, site, unidentified_feature | [T09](t09-intensities-missing-values.md) INT-2 | | present (feature, pg) |

**META-2 (MUST) Reading values.**
- An absent key, `unknown` and an unrecognised value (also a V1 error) all mean unknown.
- A key counts only in files whose `qpx_version` is at least its first version: every 1.4 key is unknown in a 1.3 file. An absent `scan_format` follows [T11](t11-spectrum-precursor-fields.md) SPC-7 (`scan` up to 1.3).

**META-3 (MUST) Where keys live.**
- The footer is the Parquet key-value metadata. qpx writers keep `ARROW:schema` metadata equal to it: they stamp keys before the first row group or re-serialize the schema on every later add (`_add_footer_metadata`, `qpx/writers/base.py:359-365`, adopts `_set_scan_format`, `:367-383`).
- In a 1.4 file an `ARROW:schema` value that differs is an error; a key missing there is not checked.

**META-4 (MUST) Rewrites.**
- merge, relayout, write-back, the OpenMS `-out_qpx` upgrade and the annotation transforms read the key-value footer and apply the carry rules through one shared function, warning on every `unknown` outcome. A transform that removes rows updates the keys it affects ([T03](t03-filtering-state.md) FLT-5, FLT-13).
- merge keeps the provenance steps of every input, each tagged with its `source_dataset` (as `qpx/core/dataset_merge.py:369-391` does today), because several [T01](t01-psm-rows.md)-[T12](t12-protein-groups.md) declarations are provenance step parameters or options ([T01](t01-psm-rows.md) PSM-4; [T02](t02-chimeric-spectra.md) CHI-10; [T03](t03-filtering-state.md) FLT-8, FLT-10–12; [T04](t04-confidence-fields.md) CNF-3; [T06](t06-identifications-per-feature.md) FID-8; [T07](t07-unidentified-features.md) UNI-5; [T08](t08-match-between-runs.md) MBR-8, MBR-9; [T11](t11-spectrum-precursor-fields.md) SPC-2; [T12](t12-protein-groups.md) PGV-6).
- Only converters (the OpenMS upgrade included) and `migrate` write a new `qpx_version`. Any other rewrite refuses a source at another version and names `qpxc migrate`.
- `migrate` writes a 1.4 key only when certain ([T09](t09-intensities-missing-values.md) INT-3) and keeps the source's provenance steps. It writes no other key: `rank1_only` and `include_unassigned_psms` stay provenance options, and rows [T07](t07-unidentified-features.md) UNI-8 moves go to `unidentified_feature`, whose presence means retained (UNI-5).

**META-5 (MUST) `acquisition_method`.**
- Values: `dda` (MS:1003221), `dia` (MS:1003215), `prm` (MS:1003907), `mixed`, `unknown`.
- Source, first that applies: (1) the SDRF `comment[proteomics data acquisition method]` of every run, mapped by accession or name (MS:1003221, NCIT:C161785 → `dda`; MS:1003215, NCIT:C161786, DIA subtypes → `dia`; MS:1003907 → `prm`; unmappable → `unknown` with a warning); (2) the converter option `--acquisition-method`; (3) an input that exists for one method only (Producer changes); (4) `unknown`.
- `mixed` requires a mappable `run.acquisition_method` for every run; a footer contradicting a run value is an error.
- DDA-scoped rules ([T01](t01-psm-rows.md), [T02](t02-chimeric-spectra.md), [T06](t06-identifications-per-feature.md), and any rule that says so) apply to `dda` files and `dda` runs under `mixed`; under `unknown` they run as info.

**META-6 (MUST) Conformance profiles** (not [T14](t14-reader-contract.md)'s reader profile).
- `core`: everything `--conformance` checks today (strict, semantic, identity) [[VALIDATOR-33](findings-cross.md#validator-33)], plus V1, V2, V4-V8.
- `dda`: `core` plus V3: every key marked in META-1 is present (`unknown` counts, except for `acquisition_method`); a run view exists when rows carry intensities; owners' profile-level errors apply.
- The validator picks the profile: `qpxc validate --profile dda`. No file claims one (META-7).

**META-7 (MUST) No claims.** No footer key names a profile, and a file makes no conformance claim. `--conformance` checks `core` unless `--profile` names another; whoever needs `dda` runs `qpxc validate --profile dda`. **(MAY)** A writer given `profile=` raises on every file-local rule of that profile; it stamps nothing.

**META-8 (MUST) Severity.**
- **error**: a violated MUST. As a transition, a 1.4 rule its owner grades "warning outside the profile" is an error only under `dda`.
- **warning**: a violated SHOULD, or such a 1.4 rule outside `dda`.
- **info**: a check that could not run (no run view, unknown key or acquisition, provided ids).
- A check may be an error or warning only if every state that fires it violates a rule; a heuristic that an allowed state can fire is info.
- `qpxc validate` fails on errors, `--conformance` also on warnings, neither on info (`qpx/cli/validate.py:125-126`).
- Hence `identity_not_checked` becomes info ([T01](t01-psm-rows.md) PSM-10 checks the composite for every `identity_source`), `pg_link_coverage` becomes info and skips null-membership and decoy features, and a skipped run check is info outside `dda`.

**META-9 (MUST) Write paths.**
- Records, tables and DataFrames get the same file-local checks, and each raises on range, vocabulary, NaN/inf, sign and footer-value violations, naming check and column.
- Vocabularies come from a YAML `enum:` attribute, not `qpx/core/data/validation.py:70-74`.
- A duplicate PK stays a write-time warning (`versioning.md:120-127`) unless the writer is given a profile (META-7).

| # | Validator check | Severity |
|---|---|---|
| V1 | Off-vocabulary or ill-typed key; differing `ARROW:schema` value (1.4) | error |
| V2 | Missing or bad version, file type or identity keys; other major | error (existing); newer minor: warning |
| V3 | `--profile dda`: presence of the keys META-1 marks (of the new ones, `fdr_filter` on psm, feature and pg; `unquantified_policy`, `intensity_quantity` and `intensity_processing` on feature and pg); `acquisition_method` `dda` or `mixed`; run view | error under `dda` |
| V4 | `acquisition_method` against `run`; `mixed` without per-run values | error |
| V5 | Declarations against rows (footer keys, `rank`, `precursor_index`, `identity_composite`, provenance parameters) and the fixed rules that replace keys: PSM-4, CHI-10, CHI-11, FLT-14, CNF-3, FRW-12, FID-8, FID-9, UNI-9, MBR-12, INT-14, PGV-13 | owner's grade |
| V6 | Row invariants: PSM-1, 7; CHI-1, 6; CNF-6; FID-1, 5, 6, 7, 10; INT-5; [T10](t10-protein-mapping.md) anchor and `unique`; SPC-3, 7, 8 | owner's grade; DDA scope (META-5) |
| V7 | `identity_not_checked`, `pg_link_coverage`, skipped run checks, DDA checks under `unknown` | info |
| V8 | Undeclared non-`x_` column; unregistered footer key (except `ARROW:schema`, `pandas`, `x_*`) or one newer than the file | warning |

**META-10 (MUST) Versioning.**
- `qpx_version` moves, with a changelog entry, when the meaning of a typed column, a footer key or value, a vocabulary, or the rows a view holds under a declared key value changes (extends `identity.md:69`).
- A converter fix that changes values but not meanings keeps the version; the changelog names its first release.
- Removing a vocabulary value is breaking; adding one is not (META-R3).

**META-11 (SHOULD) Source provenance.** The source tool's step records `tool_version` whenever the source states it [[MAXQUANT-06](findings-cross.md#maxquant-06)][[FRAGPIPE-29](findings-cross.md#fragpipe-29)][[PD-31](findings-cross.md#pd-31)], and every input read, sidecars included, its SHA-256 [[MAXQUANT-M3](findings-cross.md#maxquant-m3)].

**META-12 (MUST) Page precedence.**
- On structure (presence, type, nullability, range, vocabulary) the YAML schemas and `footer.yaml` win; on meaning, the view pages, data-model, identity, serialization, scores, intensities, peptidoform, scan, versioning and `consumer-guide.md` ([T14](t14-reader-contract.md) RDR-1).
- converter-coverage, tool-mappings, ontology's mapping tables, producer-guide and views are informative: they link to rules and state none.
- CI generates view field tables from YAML, so optional columns read "optional" [[SPEC-31](findings-cross.md#spec-31)], and `converter-coverage.md` from fixture conversions ("not verified" without a fixture). `tool-mappings.md` drops its tool tables.
- Fix `producer-guide.md:48` (`imputed_labels`, `is_imputed`) and `provenance.md:105,107` (`lfq` is a label type). `source_column_name` names the column read, or `derived`.

**META-13 (MUST) Tests.**
- Every check has a negative fixture, and a positive fixture per allowed state it touches. Every fixture conversion passes `--conformance` at the profile its test names (`--profile`).
- Tests that pin divergent behaviour change with the rule: the PD tests above; `psm_empty_scan` re-derives its id so only [T11](t11-spectrum-precursor-fields.md) SPC-8 fails it; a new positive case has a shared peptide with null anchor, membership and `pg_group_id` [[VALIDATOR-24](findings-feature.md#validator-24)].

**META-14 (SHOULD) Fixtures.**
- The golden dataset becomes DDA, in LFQ and TMT variants with every 1.4 key known: rank-2 rows, a chimeric secondary, a transferred feature, decoys, an identification-only feature and an `unidentified_feature` row ([T07](t07-unidentified-features.md)).
- Real exports: PD multi-engine TMT and Minora LFQ; FragPipe with IonQuant MBR; a Sage or MSFragger DDA+ search; a two-engine mzIdentML.

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| `footer.yaml` with the existing keys plus five new ones (`acquisition_method`, `fdr_filter`, `unquantified_policy`, `intensity_quantity`, `intensity_processing`); generated table; `project_*` dropped | schemas, docs | – | – | No |
| `acquisition_method` | psm, feature, pg | `dda`, `dia`, `prm`, `mixed`, `unknown` | SHOULD; `dda`/`mixed` under `dda` | No |
| merge keeps every input's provenance steps (META-4) | provenance | `source_dataset` per step | MUST | No (today's behaviour) |
| YAML `enum:` attribute; recipes read `metadata.metadata` | schemas, docs | – | – | No |
| Severity `info`, `--profile`, two checks to info | CLI | – | – | No; `--conformance` outcomes change |
| Writers raise on value violations | writers | – | – | **Yes** for producers that persisted invalid values |
| Rewrites refuse other versions; merged `creator` becomes a list | merge, relayout, write-back | – | – | **Yes** for callers rewriting files of another version |
| Absent key read by `qpx_version` | all | – | – | **Yes** for `scan_format` ([T11](t11-spectrum-precursor-fields.md)) |

## Producer changes
| Producer | Change needed |
|---|---|
| All converters | Stamp keys early or by the re-serializing add; add `--acquisition-method` |
| CDAP, OpenMS consensusXML and `-out_qpx` | `dda` unless the SDRF differs (DDA-only workflows); OpenMS writes `psm_qvalue`, not `q-value` ([T04](t04-confidence-fields.md)) |
| MaxQuant, PD | SDRF or option, else `unknown` (MaxDIA, PD DIA exist; Q3); record the tool version |
| FragPipe | SDRF or option; `dda` for isobaric input; else `unknown`; record MSFragger and IonQuant versions |
| mzIdentML, quantms-MSstats | SDRF or option, else `unknown` |
| DIA-NN, Spectronaut | `dia` |
| Rewrites | META-4 at the sites in META-R1, plus `qpx/writers/partitioned.py:60-67` and `migrate.py:195-197,211`; merge and migrate keep every input's provenance steps |
| Validator, CLI | Registry checks, `info`, `--profile` (no claim key); `_passes` skips info; V5 checks that read provenance run at dataset level |
| Third-party writers | Write the key-value footer; an Arrow copy is optional but must match |

## Reader rules
- **META-R1 (MUST)** Read footers from the key-value metadata ([T14](t14-reader-contract.md) RDR-2).
  - Violations: `qpx/core/parquet_io.py:11-15`, `qpx/cli/info.py:352`, `qpx/core/dataset_merge.py:240`, `qpx/core/dataset_writes.py:123`, `qpx/transforms/migrate.py:54-55`, `qpx/converters/openms/converter.py:184,211`, `qpx/converters/channel_labels.py:683`, and the four recipes.
  - mokume complies (mokume:python/mokume/io/_qpx.py:36).
- **META-R2 (MUST)** Read declarations per file and combine shards by the carry rules. Violation: `file_metadata` reads the first shard (`qpx/core/data/base.py:52-53,221-222`).
- **META-R3 (MUST)** Absent, `unknown` and unrecognised values mean unknown; never assume a convenient state such as "decoys removed". Violation: mokume reads no declaration [[MOKUME-23](findings-cross.md#mokume-23)].
- **META-R4 (MUST)** Check `qpx_version` and `file_type` first. Apply a rule only to the versions it covers, since minors may break. Violation: mokume checks neither.
- **META-R5 (MUST)** Apply DDA rules only to `dda` files and to `dda` runs under `mixed`. Violation: mokume takes `data_type` from its caller (mokume:rust/python/mokume/agentic/profiler.py:26-31), not the footer.
- **META-R6 (SHOULD)** Files make no conformance claim (META-7). A reader that needs a guarantee validates with `--profile dda`.
- **META-R7 (MUST)** Exports copy each source's registry keys and `qpx_version` into `uns`. Violations: `qpx/mudata.py:400-401` and `qpx/core/dataset_writes.py:243-244` stamp the reader's version and no keys.
- **META-R8 (MUST)** Ignore unregistered keys. An unregistered key never changes a column's meaning.

## Open questions for maintainers
1. Should `dda` reject `unknown` for `fdr_filter` and `unquantified_policy` once every qpx converter fills them ([T07](t07-unidentified-features.md) Q4), and require provenance parameters such as `decoy_policy` ([T03](t03-filtering-state.md) Q6)?
2. Wide-window DDA (DDA+, CHIMERYS): own value, or `dda` plus a cv_param [[EXTERNAL-33](findings-cross.md#external-33)]?
3. How can the MaxQuant and PD converters detect DIA output without an SDRF?
4. Should `qpxc validate` record its result (profile, outcome, qpx version) as a provenance step, now that files make no claim?
5. Should a `dia` profile reuse `dda`'s key-presence checks, and `dataset` carry `acquisition_method`?
6. Should merge split sources by `scan_format` rather than write `unknown`?
7. Should exports (META-R7) also copy the provenance parameters readers act on, such as `pass_threshold_level` ([T03](t03-filtering-state.md) FLT-R4) and `pg_grouping` ([T12](t12-protein-groups.md) PGV-R5)?
