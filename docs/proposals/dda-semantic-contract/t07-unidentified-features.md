# T07. Features without an identification

[← T06](t06-identifications-per-feature.md) · [Overview](index.md) · [T08 →](t08-match-between-runs.md)

!!! note "Proposal"
    Part of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** The feature view holds identified precursors only: every row has a canonical peptidoform, a sequence and a charge ≥ 1, and `identification_type` gets no `unidentified` value. MS1 features without an identification MAY go into a new optional view, `unidentified_feature`, with ids derived from a declared composite. Two feature-footer keys declare what was kept: `unidentified_feature_policy` and `unquantified_feature_policy`.

## Questions a reader has
- Can a feature row lack an identification? How do I recognise one: `identification_type`, an empty `peptidoform`, charge 0?
- What do its required `peptidoform`, `sequence`, `charge`, `is_decoy` and `scan` hold, and how is its `feature_id` built?
- No unidentified features in a file: did the tool detect none, or did the converter drop them?
- An identified precursor has no feature row in a run. Was it not identified there, or not quantified?
- Can I compute an identification rate, or tell missing-not-at-random (MNAR) from "not identified", from the feature view?

## What happens today
Spec:
- Unidentified Features need a producer id, which QPX namespaces by run (docs/spec/feature.md:7-10, [#229](https://github.com/bigbio/qpx/issues/229); [#232](https://github.com/bigbio/qpx/issues/232) deferred a null-peptidoform form).
- Yet `sequence`, `peptidoform` and `charge` are required (feature.md:51-55), and required means non-null (qpx/core/data/loader.py:150). So are `is_decoy` (:60) and `scan` (:139).
- feature.md:145 lists `direct`, `transferred`, `secondary` or null, but feature.yaml:210-214 has no enum. Only qpx/core/data/validation.py:70-71 checks it: `write_table` warns and `write_batch` does not check (qpx/writers/base.py:414-419,599-617) [[INFRA-11](findings-feature.md#infra-11)][[SPEC-07](findings-feature.md#spec-07)].
- psm.md:175: the OpenMS consensus converter omits features without hits. scope.md:21: "Every quantity joins to a peptide or protein accession."
- Mirror case: a null entry means identified, not quantified; no entry means not identified (producer-guide.md:43-48; intensities.md:123-126). That implies a row ([[FRAGPIPE-16](findings-feature.md#fragpipe-16)] and [[PD-12](findings-feature.md#pd-12)] call the drops divergences), but only pg requires one (pg.md:126-135, [#280](https://github.com/bigbio/qpx/issues/280)).
- [#378](https://github.com/bigbio/qpx/issues/378) (open) has not decided between "same view, `unidentified`" and "separate file".

Writer and validator (re-run):
- `FeatureWriter._should_override_provided_id` treats a falsy peptidoform as unidentified (qpx/writers/feature.py:49-51). Tests write `''` with charge 2 (tests/unit/test_identity.py:140-205); a `''`/`''`/0 probe is written and validates [[INFRA-19](findings-feature.md#infra-19)].
- A provided id is recorded as `provided` (base.py:222-226), then replaced by `derive_id([run_file_name, id])` (feature.py:71-86). One identified plus two unidentified rows give `identity_source=mixed`, so `--check-identity` skips the file (qpx/core/data/footer_validation.py:112-114); unidentified rows alone give `provided`, though their ids are hashes [[INFRA-18](findings-feature.md#infra-18)].
- Without a producer id the writer derives one from the composite; the fallback and duplicates only warn (feature.py:64-66,93-100; test_identity.py:169-189).
- The native OpenMS route sets `override_provided_ids=False` (qpx/converters/openms/converter.py:216-217), so FeatureWriter keeps raw OpenMS ids (feature.py:68-69), as feature.md:14-17 allows: feature.md:8's namespacing is not universal.
- `--semantic` fails `''` as `peptidoform_malformed` (qpx/core/semantic_validation.py:23,85-86) [[VALIDATOR-11](findings-feature.md#validator-11)]. Charge 0 passes, and the m/z check skips charge ≤ 0 (:50) [[VALIDATOR-M6](findings-cross.md#validator-m6)].

| Producer | Unidentified MS1 features | Identified, not quantified | Findings |
|---|---|---|---|
| OpenMS consensusXML | Dropped and not counted when the first PID has no hit (qpx/converters/openms_consensus/feature_adapter.py:714-716) | No row: only runs with an element > 0 get one (:348-369). ProteomicsLFQ_1: 19 precursor-runs have a PSM but no row | [[OPENMS-10](findings-feature.md#openms-10)][[OPENMS-M2](findings-feature.md#openms-m2)][[OPENMS-13](findings-feature.md#openms-13)] |
| OpenMS `-out_qpx`; `qpxc convert openms` (deprecated, qpx/cli/convert/openms.py:16) | Kept in feature. On the OpenMS branch measured in [#378](https://github.com/bigbio/qpx/issues/378), 216 of 304 rows are seeds: null `identification_type`, a grouping key reduced to `(charge, consensus_rt)` that collides, charge 0 linked to charge z. The upgrade keeps every row (converter.py:261-302) | No fixture | [[OPENMS-M5](findings-cross.md#openms-m5)] |
| MaxQuant | Never: allPeptides.txt is not read | Kept: MSMS rows get a null entry (feature_adapter.py:504-531). Full fixture: 12,758 null-only direct rows; 8,731 of them sit beside a quantified row of the same precursor-run. TMT fixture: 1,702 direct and 1,469 transferred. A NULL `Modified sequence` becomes `'None'` (:441) | [[MAXQUANT-30](findings-feature.md#maxquant-30) corrected][[MAXQUANT-09](findings-feature.md#maxquant-09)][[MAXQUANT-12](findings-feature.md#maxquant-12)][[MAXQUANT-M1](findings-cross.md#maxquant-m1)] |
| FragPipe | IonQuant reports none | Dropped: LFQ intensity 0 (feature_adapter.py:332-333,376-377); TMT plex with every channel 0 (isobaric.py:44-50 via :374). Rows are per experiment | [[FRAGPIPE-16](findings-feature.md#fragpipe-16)] (isobaric.py:255-261 is stale) |
| Proteome Discoverer | Not read | Dropped: LFQ abundance null or 0 (feature_adapter.py:80-83); no channel quantified (:89-92,99-100); every PSM "Excluded by Method" (:76). Omitting `Rejected` (:72) is correct | [[PD-12](findings-feature.md#pd-12) corrected] |
| CDAP | n/a | Kept: 6 of 94 rows have only null entries | [[MZID_CDAP-23](findings-feature.md#mzid_cdap-23)] |
| quantms MSstats | n/a | Aborts on a null intensity | [[INFRA-14](findings-feature.md#infra-14)] |
| mzIdentML | No feature view | – | – |

Readers mishandle legacy rows (UNI-R2): on a probe, `QualityControlView.metrics()` reports 2 peptides and 3 features for 1 identified feature.

## Why it matters
- One OpenMS analysis gives feature views with or without unidentified rows, by route; "not detected", "unidentified" and "dropped" look alike.
- The sentinels `''`, `'None'` and charge 0 enter peptide counts, MuData variables and id hashes; one such row disables id verification for the file.
- FragPipe, PD and OpenMS make "identified, not quantified" look like "not identified"; MaxQuant and CDAP keep the difference, so MNAR-aware imputation gets producer-dependent input.
- Unidentified features are real tool output [[EXTERNAL-15](findings-feature.md#external-15)]: MaxQuant allPeptides.txt, OpenMS `targeted_only=false` seeds, Dinosaur, Biosaur2, PSI QC term MS:4000093. mzTab-M keeps them in SMF (empty `SME_ID_REFS`) because identifications live in SME [[EXTERNAL-14](findings-feature.md#external-14)]. A QPX feature row carries its identification like an mzTab 1.0 PEP row (`sequence` mandatory), so SMF's counterpart is an identification-free view.

## Options considered
| Option | For | Against |
|---|---|---|
| A. Forbid unidentified features | Smallest; matches scope.md:21 | Reverses [#229](https://github.com/bigbio/qpx/issues/229); loses OpenMS seeds, allPeptides and QC data |
| B. Same view, `identification_type='unidentified'`, nullable `peptidoform`/`sequence`/`charge`/`is_decoy` | One file; what `-out_qpx` does; a row flag (P-C) | Breaks 4 required columns; every reader must filter (qpx MuData, views and mokume fail today); ~25 columns meaningless; contradicts [T05](t05-feature-rows.md) FRW-1, [T06](t06-identifications-per-feature.md) FID-1 |
| C. Same view with `''`/0 sentinels (today's writer) | No schema change | Violates P-D; fails `--semantic`; ids unverifiable |
| **D. Separate optional view, identified-only feature view, footer key** | Additive; existing readers stay correct; membership in the view marks the row; ids verifiable | One more view; `-out_qpx` must split |
| Mirror M1: identification-only rows MUST exist | Strongest | Infeasible for quantms MSstats, which cannot tell whether NA means identified (P-G) |
| **Mirror M2: declare the policy; SHOULD keep the rows** | Feasible, testable | Readers handle two modes |

## Recommendation
Definitions:
- An *unidentified feature* is an MS1 signal (peak, isotope pattern or multiplet) reported in a run without a link to any peptidoform, by MS/MS or by transfer. A signal that cross-run linking ties to a peptidoform is a feature row (`transferred` per [T08](t08-match-between-runs.md) MBR-1/2).
- A *precursor-run* is a `(peptidoform, charge)` in one run, or in one unit for unit rows ([T05](t05-feature-rows.md) FRW-6).
- An *identification-only row* has no non-null intensity entry ([T05](t05-feature-rows.md) FRW-5, [T09](t09-intensities-missing-values.md) INT-6).

**UNI-1 (MUST)** Every feature row is an identified precursor.
- Basic tier (`qpxc validate`): `peptidoform` and `sequence` are non-empty, and `charge` ≥ 1 ([T11](t11-spectrum-precursor-fields.md) SPC-6).
- `peptidoform` is canonical and equals `sequence` without its tags ([T11](t11-spectrum-precursor-fields.md) SPC-15). semantic_validation.py:85-88 already checks this and catches `'None'`; SPC-15 moves the check from `--semantic` to the basic tier.

Errors for `qpx_version` ≥ 1.4; warnings pointing to `qpxc migrate` before. UNI-1 replaces the unidentified sentence of feature.md:7-10 and settles [#378](https://github.com/bigbio/qpx/issues/378) for feature. It reverses [#229](https://github.com/bigbio/qpx/issues/229)/#232's placement of unidentified features in feature: their required columns cannot be filled truthfully.

**UNI-2 (MUST)** [T08](t08-match-between-runs.md) MBR-1's enum gets no `unidentified` value. A null `identification_type` means an identified row of unstated origin, never "unidentified".

**UNI-3 (MUST)** Unidentified features appear only in the optional view `unidentified_feature` (`{PREFIX}.unidentified_feature.parquet`, `file_type` `unidentified_feature_file`), never in feature or pg.

| Column | Type | Meaning |
|---|---|---|
| `feature_id` | int64, PK | Derived (UNI-4); `psm.feature_id` never points here |
| `source_feature_id` | string, null | The tool's own feature id |
| `run_file_name` | string, required | Must exist in `run` |
| `charge` | int16, null | ≥ 1; null if not assigned, never 0 |
| `observed_mz` | float32, null | Monoisotopic m/z from the feature finder |
| `rt`, `rt_start`, `rt_stop`, `rt_aligned` | float32, null | As [T05](t05-feature-rows.md) FRW-7 to FRW-9 |
| `ion_mobility`, `ion_mobility_start`, `ion_mobility_stop`, `compensation_voltage` | float32, null | As in feature |
| `intensities`, `additional_intensities` | as feature | `LFQ`, or the run's MS1-label channels |
| `consensus_feature_id` | int64, null | [T06](t06-identifications-per-feature.md) FID-11; a value never occurs in both views |
| `additional_scores`, `cv_params` | as feature | |

The view has no peptidoform, protein, decoy, q-value or `match_*` columns. Until open question 3 is decided, it holds label-free and MS1-label runs only. Its footer carries [T09](t09-intensities-missing-values.md)'s `intensity_quantity` and `intensity_normalization` (INT-2) and the identity keys of UNI-4.

**UNI-4 (MUST)** Ids in `unidentified_feature` are derived from the declared `identity_composite`, and `identity_source` is `derived`.
- A tool with its own feature id (OpenMS unique id, Dinosaur index, allPeptides row number) writes it to `source_feature_id` and declares `[run_file_name, source_feature_id]`.
- Otherwise the composite is `[run_file_name, charge, observed_mz, rt]`; if that collides, the producer writes a deterministic ordinal ([#229](https://github.com/bigbio/qpx/issues/229)) to `source_feature_id` and uses the first composite.
- A duplicate id is a write error.

This keeps [#229](https://github.com/bigbio/qpx/issues/229)'s producer id, verifiable by `--check-identity`.

**UNI-5 (MUST)** The feature footer key `unidentified_feature_policy` takes one of these values; absent means `unknown`.
- `retained`: every unidentified feature the producer obtained is in the dataset's `unidentified_feature` file.
- `removed`: the source has unidentified features and none were exported; SHOULD be counted as `skipped_rows_unidentified_feature` ([T03](t03-filtering-state.md) FLT-6).
- `not_applicable`: the source has none (the tool quantifies only identified or transferred precursors, or the input contains none).
- `unknown`.

Values follow [T03](t03-filtering-state.md) FLT-12; carriage follows [T13](t13-metadata-conformance.md) META-4. Writers SHOULD write this key and UNI-6's; the [T13](t13-metadata-conformance.md) profile requires both.

**UNI-6 (MUST)** The feature footer key `unquantified_feature_policy` takes one of these values; absent means `unknown`.
- `retained`: every target precursor-run with an accepted own-run identification ([T08](t08-match-between-runs.md) MBR-1, under [T03](t03-filtering-state.md) filters) has at least one feature row; if none is quantified, its rows are identification-only.
- `removed`: every row has a non-null entry.
- `unknown`.

Transfers are out of scope ([T08](t08-match-between-runs.md) MBR-3 forbids transferred rows without a measured entry); decoys follow [T03](t03-filtering-state.md) `decoy_policy`. The key is the feature mirror of [T03](t03-filtering-state.md) FLT-12; [T12](t12-protein-groups.md) PGV-4 mirrors it for pg.

**UNI-7 (SHOULD)** A producer whose source states per run that a precursor was identified writes `retained`. "Identified" follows the source: a rejected PSM (PD `Rejected`) identifies nothing; a PSM only excluded from quantification (PD `Excluded by Method`) identifies the precursor, with null entries. The tool's reason MAY go in the cv_param `unquantified_reason`.

**UNI-8 (MUST)** Every qpxc command that writes a 1.4 feature file from 1.3 input (`migrate`, and the deprecated `convert openms`) handles rows with a null or `''` peptidoform. It:
- moves them to `unidentified_feature` and writes `retained`;
- makes a charge ≤ 0 null;
- fills `source_feature_id` from the cv_param `provided_feature_id`, or else from the 1.3 `feature_id` when the file says `identity_source=provided` (feature.md:14-17);
- nulls any `psm.feature_id` that points at a moved row, with a warning.

**UNI-9 (MUST)** The validator reports:

| Check | Severity |
|---|---|
| UNI-1 violation | Error (≥ 1.4); warning before |
| Policy value outside its vocabulary | Error |
| `unidentified_feature_policy=retained` without the companion file, or a companion file under `removed` or `not_applicable` | Error, dataset level |
| A companion file with that key absent or `unknown` | Warning |
| `unquantified_feature_policy=removed` and an identification-only row | Error |
| `unquantified_feature_policy=retained`, a psm view, and an accepted target PSM ([T03](t03-filtering-state.md) FLT-R4) whose precursor-run has no feature row (units match through `grouped_runs`) | Warning |
| `unidentified_feature`: a duplicate id, an unknown run, a charge < 1, or a `consensus_feature_id` that also occurs in feature | Error |

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| UNI-1 checks. Text changes: feature.md:7-10, psm.md:175, scope.md:12,21 (scope gains the new view as an exception) | feature | – | – | **Yes for producers**: OpenMS `-out_qpx`, FeatureWriter's `''` path, MaxQuant `'None'`. 1.3 files only warn |
| New view `unidentified_feature` (schemas/unidentified_feature.yaml) | new | UNI-3 | Optional | No |
| `unidentified_feature_policy` | feature footer | `retained`, `removed`, `not_applicable`, `unknown` | SHOULD; MUST in [T13](t13-metadata-conformance.md) | No |
| `unquantified_feature_policy` | feature footer | `retained`, `removed`, `unknown` | SHOULD; MUST in [T13](t13-metadata-conformance.md) | No |
| cv_param `unquantified_reason` | feature | string | MAY | No |
| FeatureWriter: remove the falsy-peptidoform branch (feature.py:44-100); add `UnidentifiedFeatureWriter` | writers | – | – | API: `''` rows raise |
| UNI-8 migration | feature, new | – | – | **Yes** for migrated `-out_qpx` bundles: rows leave feature |

## Producer changes
| Producer | Keys (unidentified / unquantified) | Change needed |
|---|---|---|
| OpenMS consensusXML | `retained` or `removed` (`not_applicable` if the map has none) / `removed` until fixed, then `retained` | Test every PID for a hit, not only the first (feature_adapter.py:715). Write hitless consensus features to the new view (map position → `consensus_feature_id`, handle unique id → `source_feature_id`), or keep dropping them with a count. Add identification-only rows for precursor-runs with an accepted identification but no row (runs without a positive element; unassigned PIDs) and link them ([T06](t06-identifications-per-feature.md) FID-8) |
| OpenMS `-out_qpx` (upstream) | `retained` / as above | Seeds go to the new view with null charge, never 0 |
| MaxQuant | `removed` / `retained` | MAY read allPeptides.txt rows `Identified` leaves unmarked (check the mapping: coxdocs words it per MS/MS). Drop TMT transfers without a reporter value ([T08](t08-match-between-runs.md) MBR-3). Skip and count NULL `Modified sequence` rows; never write `'None'` |
| FragPipe | `not_applicable` / `retained` after the fix | Write null entries for an LFQ ion with `<exp> Intensity` 0 and Match Type MS/MS or spectral count > 0, and for an all-zero TMT plex |
| Proteome Discoverer | `not_applicable` for reporters; `unknown` for LFQ (the PSM export does not say whether Minora detected features independently) / `retained` after the fix | Keep "Excluded by Method" and unquantified precursors as null entries with `unquantified_reason` |
| CDAP | `not_applicable` / `retained` | converter-coverage.md:194 "Null" reads as a null list; CDAP writes null entries |
| quantms MSstats | `not_applicable` for TMT, `unknown` for LFQ (`targeted_only` unknown) / `unknown` | Write a null entry for NA instead of aborting |
| mzIdentML | none (no feature view) | – |
| qpx core | – | UNI-1 checks, `Dataset` discovery of the new view, `UnidentifiedFeatureWriter`, UNI-8 |

## Reader rules
**UNI-R1 (MUST)** A feature row is an identified precursor. Read unidentified features only from `unidentified_feature`. Never put them into roll-ups, MuData variables or peptide counts.

**UNI-R2 (MUST)** In files with `qpx_version` < 1.4, exclude as unidentified any row with a null or `''` `peptidoform`, a `charge` ≤ 0, or `identification_type` `unidentified`. Current violations:
- qpx/mudata.py:183 makes such rows the variable `'|0'` and combines them per `duplicates`.
- qpx/views/api.py:24-25,176-179 count `''` as a peptide and count these rows as features.
- `Feature.for_quantification` passes them through (qpx/core/data/feature.py:91-115).
- mokume Rust accepts `''` and errors on null (mokume:rust/crates/mokume-io/src/qpx.rs:438-439). It reads null `identification_type` with null `id_run_file_name` as transferred (qpx.rs:493-523), and counts `'unidentified'` as identified (:501-504) [[MOKUME-30](findings-feature.md#mokume-30)].
- mokume Python drops them only through the unique filter, and only when `require_unique` is set (mokume:python/mokume/io/feature.py:227-229).

**UNI-R3 (MUST)** Compute an identification rate only under `unidentified_feature_policy=retained`: `direct` and `secondary` rows with a non-null entry, over those plus transferred rows with a non-null entry plus `unidentified_feature` rows. Under other values there is no denominator.

**UNI-R4 (MUST)** Judge a precursor-run as a whole, never a single row, under `unquantified_feature_policy`:
- `retained`: no row means not identified there under the declared filters. A precursor-run whose rows are all identification-only is identified but not quantified, an MNAR candidate. An identification-only row beside a quantified one (8,731 MaxQuant `MSMS` rows) adds nothing.
- `removed` or `unknown`: take identification from accepted PSMs ([T03](t03-filtering-state.md) FLT-R4). Never read a missing row as "not identified".

**UNI-R5 (MUST)** An identification-only row is not a quantified feature: leave it out of feature counts and intensity statistics, and never fill its entries with 0. In legacy files ignore transferred rows with no non-null entry (MaxQuant TMT `MULTI-MATCH`), which [T08](t08-match-between-runs.md) MBR-3 forbids.

Violations: qpx/views/api.py:25,179 count such rows. MuData stores no null cells (qpx/mudata.py:104-106,188; anndata.md:146), so the state exists only in feature. mokume Python's `intensity > 0` (mokume:python/mokume/io/feature.py:212) is correct.

## Open questions for maintainers
1. Should the view be named `unidentified_feature` or `ms1_feature`?
2. Should the view carry the MS/MS scans acquired on each feature (MaxQuant `MSMS Scan Numbers`)?
3. Isobaric runs: should unidentified features be excluded, get a null label, or carry only `ms1_intensity` in `additional_intensities`?
4. Should UNI-6 `retained` become MUST in the [T13](t13-metadata-conformance.md) profile once FragPipe, PD and OpenMS comply?
5. Should MuData gain an `identified_mask` layer, so that "identified, not quantified" survives export?
