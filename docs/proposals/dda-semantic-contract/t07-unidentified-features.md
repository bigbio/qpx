# T07. Features without an identification

[← T06](t06-identifications-per-feature.md) · [Overview](index.md) · [Strict profile](strict-profile.md) · [Annex overview](annex-overview.md) · [T08 →](t08-match-between-runs.md)

!!! note "Annex: producer-neutral analysis"
    Part of the annex of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text, and they cover every producer. For OpenMS-written files the [strict profile](strict-profile.md) is the contract and wins where the two differ. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** MS1 features without an identification stay in the feature view, in the same file as identified ones. A producer MAY write them; each is marked `identification_type = 'unidentified'` and has a null `peptidoform` and `sequence`, and readers ignore these rows by default. Every other row has a canonical peptidoform, a sequence and a charge ≥ 1. There is no separate view, no footer key and no policy parameter for unidentified features: their rows show whether they are present. The provenance parameter `unquantified_rows` states whether identified, unquantified precursors keep a row.

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
| OpenMS `-out_qpx`; `qpxc convert openms` (deprecated, qpx/cli/convert/openms.py:16) | Kept in feature: at HEAD (verified) as rows with `peptidoform = ""`. On the OpenMS branch measured in [#378](https://github.com/bigbio/qpx/issues/378), 216 of 304 rows are seeds: null `identification_type`, a grouping key reduced to `(charge, consensus_rt)` that collides, charge 0 linked to charge z. The upgrade keeps every row (converter.py:261-302) | No fixture. At HEAD (verified) intensity entries are non-nullable, so no null entry can be written | [[OPENMS-M5](findings-cross.md#openms-m5)] |
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
- Unidentified features are real tool output [[EXTERNAL-15](findings-feature.md#external-15)]: MaxQuant allPeptides.txt, OpenMS `targeted_only=false` seeds, Dinosaur, Biosaur2, PSI QC term MS:4000093. mzTab-M keeps them in SMF, beside identified features, with empty `SME_ID_REFS` [[EXTERNAL-14](findings-feature.md#external-14)]. A QPX feature row carries its identification itself, like an mzTab 1.0 PEP row, so an unidentified row needs null identification columns and a marker.

## Options considered
| Option | For | Against |
|---|---|---|
| A. Forbid unidentified features | Smallest; matches scope.md:21 | Reverses [#229](https://github.com/bigbio/qpx/issues/229); loses OpenMS seeds, allPeptides and QC data |
| **B. Same view, optional rows marked `identification_type = 'unidentified'`, null `peptidoform`/`sequence`, nullable `charge`** | One file; keeps [#229](https://github.com/bigbio/qpx/issues/229)/#232; a row flag (P-C); producers that write none change nothing | Three required columns become nullable; every reader must exclude the rows (qpx MuData, views and mokume fail today); about 25 columns are null on them |
| C. Same view with `''`/0 sentinels (today's writer; what `-out_qpx` does at HEAD) | No schema change | Violates P-D; fails `--semantic`; ids unverifiable |
| D. Separate optional view, identified-only feature view | Additive; existing readers stay correct; membership in the view marks the row | One more view and file; `-out_qpx` must split; reverses [#229](https://github.com/bigbio/qpx/issues/229)/#232. Rejected by the maintainer: unidentified features must be storable in the same file |
| Mirror M1: identification-only rows MUST exist | Strongest | Infeasible for quantms MSstats, which cannot tell whether NA means identified (P-G) |
| **Mirror M2: declare the policy; SHOULD keep the rows** | Feasible, testable | Readers handle two modes |

## Recommendation
Definitions:
- An *unidentified feature* is an MS1 signal (peak, isotope pattern or multiplet) reported in a run without a link to any peptidoform, by MS/MS or by transfer. A signal that cross-run linking ties to a peptidoform is a feature row (`transferred` per [T08](t08-match-between-runs.md) MBR-1/2).
- A *precursor-run* is a `(peptidoform, charge)` in one run, or in one unit for unit rows ([T05](t05-feature-rows.md) FRW-6).
- An *identification-only row* has no non-null intensity entry ([T05](t05-feature-rows.md) FRW-5, [T09](t09-intensities-missing-values.md) INT-6).

**UNI-1 (MUST)** A feature row is an identified precursor, or an unidentified feature marked as such. A producer MAY write unidentified features as rows of the feature view; none is required to.
- A row is unidentified exactly when `identification_type = 'unidentified'` (UNI-2). Its columns follow UNI-3.
- Every other row is an identified precursor: `peptidoform` and `sequence` are non-null, and `charge` ≥ 1 ([T11](t11-spectrum-precursor-fields.md) SPC-6). The validator checks this at the basic tier (`qpxc validate`).
- `''`, `'None'` and charge 0 are errors on every row.
- `peptidoform` is canonical and equals `sequence` without its tags ([T11](t11-spectrum-precursor-fields.md) SPC-15). semantic_validation.py:85-88 already checks this and catches `'None'`; SPC-15 moves the check from `--semantic` to the basic tier.
- The psm view is unchanged: `peptidoform` is never null, and `charge` ≥ 1.

Errors for `qpx_version` ≥ 1.4; warnings pointing to `qpxc migrate` before. UNI-1 keeps [#229](https://github.com/bigbio/qpx/issues/229)/#232: unidentified features live in the feature view (feature.md:7-10). It adds the marker and replaces the `''` and 0 sentinels by nulls. It settles [#378](https://github.com/bigbio/qpx/issues/378) as: same view, flagged, optional, ignored by readers.

**UNI-2 (MUST)** [T08](t08-match-between-runs.md) MBR-1's enum is `direct`, `secondary`, `transferred`, `unidentified` and null. `unidentified` marks an unidentified feature. A null `identification_type` means an identified row of unstated origin, never "unidentified".

**UNI-3 (MUST)** Unidentified features appear only as rows of the feature view marked `identification_type = 'unidentified'`. There is no separate view. On such a row:

| Column | Value |
|---|---|
| `peptidoform`, `sequence` | null, never `''` or `'None'` |
| `charge` | The feature's charge when the feature finder assigned one, else null; never 0 |
| `is_decoy` | false |
| `scan` | `[]` |
| `id_run_file_name`, PEP, every q-value, `match_*` | null |
| `contaminant`, `pg_accessions`, `anchor_protein`, `pg_group_id`, `unique`, `group_unique`, `identification_conflict` | null |
| `intensities`, `rt`, `rt_start`, `rt_stop`, `rt_aligned`, `observed_mz`, ion mobility, `consensus_feature_id` | The same rules as on other rows ([T05](t05-feature-rows.md) FRW-7 to FRW-9, [T06](t06-identifications-per-feature.md) FID-11, [T09](t09-intensities-missing-values.md), [T11](t11-spectrum-precursor-fields.md)) |

- `peptidoform`, `sequence` and `charge` become nullable in the feature schema, for these rows only (UNI-1).
- An unidentified row is never evidence: no `psm.feature_id` points to it, and it never enters pg evidence, counts, pepmap or protein quantities. Unidentified features have no pg rows.
- Until open question 3 is decided, unidentified rows are written for label-free and MS1-label runs only.

**UNI-4 (MUST)** Ids of unidentified rows are unchanged from [#229](https://github.com/bigbio/qpx/issues/229): the producer supplies an id (OpenMS unique id, Dinosaur index, allPeptides row number), and QPX namespaces it by run (feature.md:7-10). The original id stays in the cv_param `provided_feature_id`.
- The writer recognises an unidentified row by `identification_type = 'unidentified'` or a null peptidoform, not by a falsy `''` peptidoform (feature.py:49-51).
- These ids count as derived from `[run_file_name, provided id]`. A file with such rows then no longer reports `identity_source` `mixed` or `provided` [[INFRA-18](findings-feature.md#infra-18)], and `--check-identity` can verify it.
- A duplicate id is a write error.

**UNI-5 (MUST)** No footer key and no policy parameter states what happened to unidentified features: the rows themselves show whether they are present. Writing them is optional (UNI-1).
- A file that holds unidentified rows holds every unidentified feature the producer obtained.
- A producer that drops unidentified features the source has MAY record their count as the provenance parameter `skipped_unidentified_features` on its `format_conversion` step ([T03](t03-filtering-state.md) FLT-6).
- A file without unidentified rows and without that count says nothing about them.

**UNI-6 (MUST)** The provenance parameter `unquantified_rows`, on the step that wrote the QPX files (`format_conversion`, or the tool's own export step), takes one of these values; absent means unknown.
- `retained`: every target precursor-run with an accepted own-run identification ([T08](t08-match-between-runs.md) MBR-1, under [T03](t03-filtering-state.md) filters) has at least one feature row; if none is quantified, its rows are identification-only.
- `removed`: every row has a non-null entry.

Transfers are out of scope ([T08](t08-match-between-runs.md) MBR-3 forbids transferred rows without a measured entry); decoys follow provenance `decoy_policy` ([T03](t03-filtering-state.md) FLT-10). The parameter is the feature mirror of [T03](t03-filtering-state.md) FLT-12's provenance option `include_unassigned_psms`; [T12](t12-protein-groups.md) PGV-4 gives the same parameter its pg meaning. One value covers the feature and pg views the step wrote; a producer for which the two differ records none. A merge keeps every input's steps ([T13](t13-metadata-conformance.md) META-4). Writers SHOULD record the parameter; the [T13](t13-metadata-conformance.md) profile requires it.

**UNI-7 (SHOULD)** A producer whose source states per run that a precursor was identified records `unquantified_rows = retained`. "Identified" follows the source: a rejected PSM (PD `Rejected`) identifies nothing; a PSM only excluded from quantification (PD `Excluded by Method`) identifies the precursor, with null entries. The tool's reason MAY go in the cv_param `unquantified_reason`.

**UNI-8 (MUST)** Every qpxc command that writes a 1.4 feature file from 1.3 input (`migrate`, and the deprecated `convert openms`) handles rows with a null, `''` or `'None'` peptidoform. It rewrites them in place; nothing moves to another file. It:
- makes `peptidoform` and `sequence` null;
- makes a charge ≤ 0 null;
- sets `identification_type = 'unidentified'`, and the other columns as UNI-3 requires;
- nulls any `psm.feature_id` that points at such a row, with a warning.

A row that has a real peptidoform but a charge ≤ 0 keeps its peptidoform. The command reports it; it fails validation until the source is converted again.

**UNI-9 (MUST)** The validator reports:

| Check | Severity |
|---|---|
| UNI-1 violation: a row that is not `unidentified` with a null `peptidoform` or `sequence` or a null `charge`; on any row, `''`, `'None'` or a charge < 1 | Error (≥ 1.4); warning before |
| An `unidentified` row that breaks UNI-3: a non-null `peptidoform`, `sequence` or other column that UNI-3 sets to null, `is_decoy = true`, or a non-empty `scan` | Error |
| A `psm.feature_id` that points to an `unidentified` row | Error, dataset level |
| `unquantified_rows` value outside its vocabulary | Error, dataset level |
| `unquantified_rows = removed` and an identification-only row | Error, dataset level |
| `unquantified_rows = retained`, a psm view, and an accepted target PSM ([T03](t03-filtering-state.md) FLT-R4) whose precursor-run has no feature row (units match through `grouped_runs`) | Warning |
| A duplicate id among `unidentified` rows | Error |

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| `peptidoform`, `sequence` and `charge` become nullable; UNI-1 checks at the basic tier. Text changes: feature.md:7-10,51-55, psm.md:175, scope.md:12,21 (unidentified rows are the stated exception) | feature | null only on `unidentified` rows; elsewhere a canonical peptidoform and `charge` ≥ 1 | Non-null on every row that is not `unidentified` | **Yes for readers**: these columns can be null, and predicates must exclude the rows (UNI-R1). **Yes for producers**: OpenMS `-out_qpx`, FeatureWriter's `''` path, MaxQuant `'None'`. 1.3 files only warn |
| `identification_type` gains `unidentified` (UNI-2; [T08](t08-match-between-runs.md) MBR-1) | feature | `direct`, `secondary`, `transferred`, `unidentified`, null | On every unidentified row | No: a new value ([T13](t13-metadata-conformance.md) META-10) |
| Unidentified rows (UNI-3); no new view, no footer key | feature | UNI-3 | Optional (MAY) | No |
| `skipped_unidentified_features` (UNI-5) | provenance (`format_conversion` step) | count | MAY | No |
| `unquantified_rows` (pg: [T12](t12-protein-groups.md) PGV-4) | provenance (the step that wrote the QPX files) | `retained`, `removed` | SHOULD; MUST in [T13](t13-metadata-conformance.md) | No |
| cv_param `unquantified_reason` | feature | string | MAY | No |
| FeatureWriter: recognise unidentified rows by `identification_type = 'unidentified'` or a null peptidoform, not by the falsy-peptidoform branch (feature.py:44-100); record their ids as derived from `[run_file_name, provided id]` (UNI-4) | writers | – | – | API: `''` rows raise |
| UNI-8 migration | feature | – | – | **Yes** for migrated `-out_qpx` bundles: `''` and charge 0 become null and the rows are marked `unidentified`; they stay in feature |

## Producer changes
| Producer | Unidentified features (UNI-5) / provenance `unquantified_rows`, as the feature view allows (UNI-6) | Change needed |
|---|---|---|
| OpenMS consensusXML | Unidentified rows, or a `skipped_unidentified_features` count / `removed` until fixed, then `retained` | Test every PID for a hit, not only the first (feature_adapter.py:715). Write hitless consensus features as unidentified rows (map position → `consensus_feature_id`, handle unique id → the producer id of UNI-4), or keep dropping them with a count. Add identification-only rows for precursor-runs with an accepted identification but no row (runs without a positive element; unassigned PIDs) and link them ([T06](t06-identifications-per-feature.md) FID-8) |
| OpenMS `-out_qpx` (upstream) | Unidentified rows (optional) / as above | Verified at HEAD: unidentified features are feature rows with `peptidoform = ""`, and intensity entries are non-nullable. Seeds, when OpenMS writes them, stay in feature as rows marked `identification_type = 'unidentified'`, with null `peptidoform` and `sequence`, and the feature's charge or null when none was assigned, never 0 (UNI-3); entries become nullable, so identification-only rows can be written. See the [OpenMS worklist](openms-worklist.md). |
| MaxQuant | None (allPeptides.txt not read), or unidentified rows once it is read / `retained` | MAY read allPeptides.txt rows `Identified` leaves unmarked (check the mapping: coxdocs words it per MS/MS). Drop TMT transfers without a reporter value ([T08](t08-match-between-runs.md) MBR-3). Skip and count NULL `Modified sequence` rows; never write `'None'` |
| FragPipe | None (IonQuant reports none) / `retained` after the fix | Write null entries for an LFQ ion with `<exp> Intensity` 0 and Match Type MS/MS or spectral count > 0, and for an all-zero TMT plex |
| Proteome Discoverer | None (not read; the PSM export does not say whether Minora detected features independently) / `retained` after the fix | Keep "Excluded by Method" and unquantified precursors as null entries with `unquantified_reason` |
| CDAP | None / `retained` | converter-coverage.md:194 "Null" reads as a null list; CDAP writes null entries |
| quantms MSstats | None (for LFQ, `targeted_only` is unknown) / not recorded (unknown) | Write a null entry for NA instead of aborting |
| mzIdentML | none (no feature view) | – |
| qpx core | – | Nullable `peptidoform`, `sequence` and `charge` in the feature schema; UNI-1 and UNI-3 checks; FeatureWriter per UNI-4; read helpers exclude unidentified rows by default (UNI-R1); UNI-8 |

## Reader rules
**UNI-R1 (MUST)** Ignore unidentified rows by default. Every reader predicate on the feature view excludes them; the portable test is `peptidoform IS NOT NULL` (equivalently `identification_type IS DISTINCT FROM 'unidentified'`). Never put them into roll-ups, MuData variables, peptide or feature counts, pg evidence or protein quantities. A reader uses them only on request (QC, UNI-R3).

**UNI-R2 (MUST)** In files with `qpx_version` < 1.4, read as unidentified, and so exclude by default, any row with a null, `''` or `'None'` `peptidoform`, a `charge` ≤ 0, or `identification_type` `unidentified`. Current violations:
- qpx/mudata.py:183 makes such rows the variable `'|0'` and combines them per `duplicates`.
- qpx/views/api.py:24-25,176-179 count `''` as a peptide and count these rows as features.
- `Feature.for_quantification` passes them through (qpx/core/data/feature.py:91-115).
- mokume Rust accepts `''` and errors on null (mokume:rust/crates/mokume-io/src/qpx.rs:438-439). It reads null `identification_type` with null `id_run_file_name` as transferred (qpx.rs:493-523), and counts `'unidentified'` as identified (:501-504) [[MOKUME-30](findings-feature.md#mokume-30)].
- mokume Python drops them only through the unique filter, and only when `require_unique` is set (mokume:python/mokume/io/feature.py:227-229).

**UNI-R3 (MUST)** Compute an identification rate only when the feature file holds unidentified rows and provenance counts no `skipped_unidentified_features` (UNI-5): `direct` and `secondary` rows with a non-null entry, over those plus transferred rows with a non-null entry plus unidentified rows. Otherwise there is no denominator.

**UNI-R4 (MUST)** Judge a precursor-run as a whole, never a single row, under the provenance parameter `unquantified_rows`:
- `retained`: no row means not identified there under the declared filters. A precursor-run whose rows are all identification-only is identified but not quantified, an MNAR candidate. An identification-only row beside a quantified one (8,731 MaxQuant `MSMS` rows) adds nothing.
- `removed` or absent (unknown): take identification from accepted PSMs ([T03](t03-filtering-state.md) FLT-R4). Never read a missing row as "not identified".

**UNI-R5 (MUST)** An identification-only row is not a quantified feature: leave it out of feature counts and intensity statistics, and never fill its entries with 0. In legacy files ignore transferred rows with no non-null entry (MaxQuant TMT `MULTI-MATCH`), which [T08](t08-match-between-runs.md) MBR-3 forbids.

Violations: qpx/views/api.py:25,179 count such rows. MuData stores no null cells (qpx/mudata.py:104-106,188; anndata.md:146), so the state exists only in feature. mokume Python's `intensity > 0` (mokume:python/mokume/io/feature.py:212) is correct.

## Open questions for maintainers
1. *Resolved by the maintainer:* unidentified features stay in the feature view, in the same file as identified ones. They are marked `identification_type = 'unidentified'`, optional for producers, and ignored by readers by default (UNI-1 to UNI-3, UNI-R1). This settles [#378](https://github.com/bigbio/qpx/issues/378); there is no separate view to name.
2. Should an unidentified row carry the MS/MS scans acquired on its feature (MaxQuant `MSMS Scan Numbers`)? UNI-3 sets `scan = []`.
3. Isobaric runs: should unidentified features be excluded, get a null label, or carry only `ms1_intensity` in `additional_intensities`?
4. Should UNI-6 `retained` become MUST in the [T13](t13-metadata-conformance.md) profile once FragPipe, PD and OpenMS comply?
5. Should MuData gain an `identified_mask` layer, so that "identified, not quantified" survives export?
6. A merge of datasets with and without unidentified rows yields a file whose unidentified rows cover only some runs, yet their presence means every one was written (UNI-5). Should merge drop these rows, or provenance list the runs they cover?
