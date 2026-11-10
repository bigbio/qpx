# T06. Identifications per feature: the identifying PSM, conflicts, feature.scan and psm.feature_id

[← T05](t05-feature-rows.md) · [Overview](index.md) · [Strict profile](strict-profile.md) · [Annex overview](annex-overview.md) · [T07 →](t07-unidentified-features.md)

!!! note "Annex: producer-neutral analysis"
    Part of the annex of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text, and they cover every producer. For OpenMS-written files the [strict profile](strict-profile.md) is the contract and wins where the two differ. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** An identified DDA feature row has one peptidoform and one charge, 0..n *supporting PSMs* linked by `psm.feature_id`, and at most one *identifying PSM* from its own run: the one the source names, else lowest PEP, then q-value, then scan. Only that PSM feeds `scan` (one spectrum), `id_run_file_name`, PEP and PSM scores. A quantified signal feeds at most one row; competing peptidoforms set `identification_conflict`. Run rows of one cross-run feature share `consensus_feature_id`. A producer writing both views links every supporting PSM; provenance `skipped_runs_psm` names runs it could not link.

## Questions a reader has
- How many identifications stand behind a feature row? Which PSM gave it `scan`, `id_run_file_name`, PEP and scores, and how do I find it?
- Is `feature.scan` one spectrum or a list, and in which run's numbering?
- When two peptidoforms claim one peak, do I get two full-intensity rows, one row, or none?
- Which PSMs carry `psm.feature_id`; does a null link mean "not quantified"? What does `spectral_count` count?
- How are the run rows of one cross-run feature grouped?

## What happens today
**Spec.**
- `id_run_file_name` "may differ from `run_file_name`" (`feature.md:138`), yet is null for a transferred row (`:145`) [[SPEC-04](findings-feature.md#spec-04)] [[CRITIC-27](findings-feature.md#critic-27)] [[VALIDATOR-M5](findings-feature.md#validator-m5)].
- `scan` is "the best PSM's" (`feature.md:139,301`); "best" is defined only per converter (FragPipe lowest PEP, converter-coverage.md:223-224; PD lowest q-value, `:324`) [[SPEC-05](findings-feature.md#spec-05)] [[CRITIC-30](findings-feature.md#critic-30)] [[PD-05](findings-cross.md#pd-05)] [[EXTERNAL-16](findings-feature.md#external-16)]. OpenMS may concatenate every supporting spectrum (`feature.md:26-28`, `scan.md:90-92,123-124`) [[CRITIC-29](findings-feature.md#critic-29)] [[VALIDATOR-12](findings-feature.md#validator-12)].
- `psm.feature_id` is the authoritative, partial link (`psm.md:75,147-161`); only MS1-label partners may differ in peptidoform (`intensities.md:102-105`), and nothing constrains run, charge, rank or decoy status [[VALIDATOR-15](findings-cross.md#validator-15)].
- `psm.md:286` (no PEP on features) contradicts `feature.md:56` [[EXTERNAL-M4](findings-feature.md#external-m4)]. `spectral_count` is "MS/MS scans of the evidence row … or the precursor's PSMs" (`intensities.md:187`). mz joins feature on `run_file_name + scan` (`data-model.md:166`).
- Nothing covers several peptidoforms on one signal [[SPEC-06](findings-feature.md#spec-06)] [[CRITIC-23](findings-feature.md#critic-23)] [[EXTERNAL-14](findings-feature.md#external-14)]. `converter-coverage.md:285` wrongly says MaxQuant fills no links [[SPEC-M5](findings-psm.md#spec-m5)].

| Producer | Behaviour | Findings |
|---|---|---|
| MaxQuant | Scan (`Best MS/MS`) and PEP come from the evidence row; that MS/MS had the lowest PEP in 75,493 of 75,493 rows ([#378](https://github.com/bigbio/qpx/issues/378)). MQ 1.5 `MS/MS Scan Number` is not matched (empty scan). `id_run_file_name` is nulled only for `MULTI-MATCH` (qpx/converters/maxquant/feature_adapter.py:657), so `MULTI-MATCH-MSMS` rows are `transferred` with run, scan and PEP. Links from `Evidence ID`, unchecked (psm_adapter.py:229-232). 312 isotope patterns (same run, charge, m/z, `Intensity`) carry 2–3 sequences, each a `MULTI-MSMS` row with its own MS/MS and the full intensity (re-checked) | [[MAXQUANT-01](findings-feature.md#maxquant-01)] [[MAXQUANT-08](findings-feature.md#maxquant-08)] [[MAXQUANT-25](findings-psm.md#maxquant-25)] [[MAXQUANT-27](findings-cross.md#maxquant-27)] [[SPEC-M1](findings-feature.md#spec-m1)] |
| FragPipe | qpx's pick: lowest PEP, then scan (qpx/converters/fragpipe/feature_adapter.py:72-74,476); key omits the FAIMS CV (:341). Multi-run experiments: no PSM, null run (:340), first run as `run_file_name` (:350). One-run rows keep the run when no PSM is found. No links; `Spectral Count` not read | [[FRAGPIPE-09](findings-psm.md#fragpipe-09)] [[FRAGPIPE-12](findings-feature.md#fragpipe-12)] [[FRAGPIPE-13](findings-feature.md#fragpipe-13)] [[FRAGPIPE-15](findings-feature.md#fragpipe-15)] [[FRAGPIPE-17](findings-feature.md#fragpipe-17)] [[SPEC-M2](findings-feature.md#spec-m2)] |
| Proteome Discoverer | Rank-1 PSMs of a precursor; pick lowest q-value, then scan (qpx/converters/proteome_discoverer/feature_adapter.py:16-20); PD's peptide groups pick by PEP (660 of 3,026 differ, [#378](https://github.com/bigbio/qpx/issues/378)). psm de-duplication breaks ties differently (psm_adapter.py:39-40). `spectral_count` counts rows (two engines count twice), absent on TMT (feature_adapter.py:70,84). No links | [[PD-03](findings-cross.md#pd-03)] [[PD-04](findings-feature.md#pd-04)] [[PD-05](findings-cross.md#pd-05)] [[PD-06](findings-cross.md#pd-06)] [[PD-M4](findings-feature.md#pd-m4)] |
| CDAP | Lowest SpecEValue, full tie order (qpx/converters/cdap/feature_adapter.py:341-372). `AmbiguousMatch` ties are dropped from features (:116-120) and from pg. PSM count in cv_param `psm_count` (:384-388). No links | [[MZID_CDAP-15](findings-feature.md#mzid_cdap-15)] [[MZID_CDAP-16](findings-feature.md#mzid_cdap-16)] [[MZID_CDAP-29](findings-psm.md#mzid_cdap-29)] [[CRITIC-M3](findings-psm.md#critic-m3)] |
| OpenMS (qpxc) | Peptide from `pids[0]` (qpx/converters/openms_consensus/feature_adapter.py:724); in PSM-level mode PEP and q from the run's first PID (:325-345). `scan` concatenates every PID of the run, other peptides included (:303-322), re-flattened by `sorted(set())` when duplicate peaks merge (feature_dedup.py:48). Every PSM of the run is linked, decoys and other peptidoforms included (record_linking.py:66-70). A conflicting run gets null type and run but keeps PEP (feature_adapter.py:533-534) | [[OPENMS-06](findings-feature.md#openms-06)] [[OPENMS-07](findings-feature.md#openms-07)] [[OPENMS-08](findings-feature.md#openms-08)] [[OPENMS-11](findings-feature.md#openms-11)] [[OPENMS-20](findings-psm.md#openms-20)] [[OPENMS-21](findings-cross.md#openms-21)] |
| OpenMS `-out_qpx` | HEAD (verified): peptide, PEP, `scan` and `id_run_file_name` come from the best hit across all runs of the consensus feature, so a row can name a donor run. The unmerged [#378](https://github.com/bigbio/qpx/issues/378) branch adds a per-run identifying PSM and own-run links ([#378](https://github.com/bigbio/qpx/issues/378) rules 1–3 and 5); no `consensus_feature_id` yet. See the [OpenMS worklist](openms-worklist.md). | [#378](https://github.com/bigbio/qpx/issues/378) |
| quantms-MSstats | `id_run_file_name = run_file_name` always, type null (qpx/converters/quantms_msstats/feature_adapter.py:525); LFQ `Reference` names only the file, so `scan` is empty. Wrong for 18 of 85 rows in ProteomicsLFQ_3 ([#378](https://github.com/bigbio/qpx/issues/378)) | [[INFRA-27](findings-feature.md#infra-27)] [[INFRA-M2](findings-feature.md#infra-m2)] |

## Why it matters
- Feature PEP, scan and scores come from different PSMs per producer; joins to psm return another peptide; OpenMS scans cannot be split into spectra.
- A peak credited to two peptidoforms is counted twice in rollups; a dropped tie loses a quantified peptide silently.
- "No link" and "not written" look alike, which breaks spectral counting [[MOKUME-25](findings-psm.md#mokume-25)]; inferring MBR from a null `id_run_file_name` misclassifies rows [[MOKUME-07](findings-feature.md#mokume-07)].

## Options considered
| Option | For | Against |
|---|---|---|
| Hardlink `feature.best_psm_id` | Direct join | Withdrawn in [#378](https://github.com/bigbio/qpx/issues/378): needs the same rule plus another id |
| **Rule over existing `scan`/`id_run_file_name` plus links ([#378](https://github.com/bigbio/qpx/issues/378))** | No new id; works for MS1 labels via links | Needs links or [T01](t01-psm-rows.md)'s unique psm key |
| "Best" = highest engine score | MaxQuant wording | Not comparable across engines |
| **Source's pick, else PEP → q → documented score → scan → peptidoform** ([#378](https://github.com/bigbio/qpx/issues/378) rule 4, made total) | Deterministic; matches MaxQuant and PD | CDAP reaches SpecEValue only after q-value |
| Conflict: full-intensity row per peptidoform | No choice | Double counting; [T05](t05-feature-rows.md) FRW-2 forbids it |
| Conflict: rows share a group id (mzTab-M) | Keeps every candidate | Every reader must de-duplicate |
| Conflict: drop the signal (FlashLFQ, CDAP) | Safe | Loses signal silently |
| **Conflict: one flagged row keeps the signal; losers optional, no intensity** | No double count, nothing hidden | New column |
| `scan` = all supporting spectra (as mzTab 1.0 `spectra_ref`) | Complete | Not one native ID, or a breaking type |
| **`scan` = one spectrum; others via links** | scan_format checkable; joins mz | OpenMS output changes |
| Group by `(peptidoform, charge, consensus_rt)` | No column | Collides and splits ([#378](https://github.com/bigbio/qpx/issues/378) rule 8) |
| **Integer `consensus_feature_id`** | Covers every measured case | New column |

## Recommendation
Scope: identified DDA rows, that is rows of runs whose `run.acquisition_method` is `dda` ([T13](t13-metadata-conformance.md) META-5). DIA keeps scan.md:111. An unidentified row ([T07](t07-unidentified-features.md) UNI-3) has no peptidoform, no identifying or supporting PSM and no link; of the rules below only FID-11 applies to it.

Definitions:
- A *quantified signal* is the MS1 peak or multiplet, or the reporter spectrum, that gives a row its intensities.
- A *supporting PSM* is a rank-1 ([T01](t01-psm-rows.md) PSM-R1) psm row whose spectrum the producer used to identify or quantify the row, with the row's `charge`, `peptidoform` (MS1-label partner: `psm.label ∈ identified_labels`) and decoy status, in `run_file_name` (unit rows: `grouped_runs`, [T05](t05-feature-rows.md) FRW-6).

**FID-1 (MUST)** A row has one peptidoform, one charge, 0..n supporting PSMs and at most one identifying PSM.
- A per-run `direct` or `secondary` row has exactly one identifying PSM in the source; the file may fail to locate it (FID-5, FID-6).
- `transferred` rows ([T08](t08-match-between-runs.md)) have no identifying or supporting PSM. Rows of a multi-run unit ([T05](t05-feature-rows.md) FRW-6) name no identifying PSM: empty `scan`, null `id_run_file_name`, PEP and `psm_qvalue`. A one-run unit row follows the per-run rules.
- A row whose run has PSMs only of competing peptidoforms has no identifying PSM and is never `direct`: `transferred` where the tool says so ([T08](t08-match-between-runs.md) MBR-1, MBR-2; OpenMS: the run is listed in the consensus feature's `mbr_transfer_map_indices`), else null; FID-10 applies.

**FID-2 (MUST)** The identifying PSM is a supporting PSM of the row's own run, never a donor-run PSM. Donors are reached through `consensus_feature_id` (FID-11) and `match_*` ([T08](t08-match-between-runs.md) MBR-5).

**FID-3 (MUST)** The identifying PSM is:
1. the PSM the source names: MaxQuant `Best MS/MS` (`MS/MS scan number`); in the survivor's run, the PID kept by OpenMS IDConflictResolver;
2. otherwise the first supporting PSM of the run (on a `direct` row, [T02](t02-chimeric-spectra.md) `precursor_index = 1` first) by: lowest PEP; lowest `psm.qvalue`; the engine's primary score in its direction, where converter-coverage.md documents one (CDAP SpecEValue); lowest `scan`, component-wise; `peptidoform` in byte order.

Nulls sort last. [T01](t01-psm-rows.md) PSM-2 and [T02](t02-chimeric-spectra.md) index assignment use the same PEP → q-value prefix.

**FID-4 (MUST)** `scan`, `id_run_file_name`, `posterior_error_probability`, [T04](t04-confidence-fields.md)'s `psm_qvalue`, PSM-level `additional_scores` ([T04](t04-confidence-fields.md) CNF-8), `peptidoform` and `charge` equal the identifying PSM's (an MS1 multiplet keeps the lightest partner's peptidoform, [T05](t05-feature-rows.md) FRW-10); without one they are null or empty. `rt`, `rt_start`, `rt_stop`, `observed_mz`, intensities and precursor, peptide or global q-values are feature-level ([T04](t04-confidence-fields.md) CNF-2, [T05](t05-feature-rows.md) FRW-8, [T09](t09-intensities-missing-values.md), [T11](t11-spectrum-precursor-fields.md) SPC-1).

**FID-5 (MUST)** `scan` holds exactly one spectrum, the identifying PSM's, under the file's `scan_format` ([T11](t11-spectrum-precursor-fields.md)); `[]` when none is located. Supporting spectra are reached only through `psm.feature_id`.

**FID-6 (MUST)** `id_run_file_name` is `run_file_name` or null, and non-null if and only if `scan` is non-empty. It is null when there is no identifying PSM or its spectrum cannot be located; `identification_type` still states what the tool reported.

**FID-7 (MUST)** Links.
- A row has 0..n linked PSMs; a PSM links to at most one row, of which it is a supporting PSM: no other peptidoform (MS1-label partners aside), charge or run (unit rows aside), no decoy PSM on a target row, no `rank > 1` ([T01](t01-psm-rows.md) PSM-7), no transferred row, no unidentified row ([T07](t07-unidentified-features.md) UNI-3).
- When several rows qualify, the source's assignment wins (MaxQuant `Evidence ID`). If it names several ([T01](t01-psm-rows.md) PSM-2 merged twins), link the row the PSM identifies; otherwise use the psm.md:158-161 order.

**FID-8 (MUST)** A producer that writes both views writes `psm.feature_id` on every supporting PSM it exports (including when it builds features from PSMs). The validator checks it.
- A producer that cannot cover every run (e.g. a single FragPipe psm.tsv) records the uncovered runs in provenance `skipped_runs_psm`.
- Each non-empty `feature.scan` resolves to exactly one psm row by FID-R1, unless a declared psm filter (a cut stated by provenance `psm_fdr`, [T03](t03-filtering-state.md) FLT-1) removed it.
- `--no-include-unassigned-psms` keeps identifying PSMs [[OPENMS-21](findings-cross.md#openms-21)].

**FID-9 (MUST)** No footer key states link completeness. By FID-8, links in a 1.4 file are complete (every exported supporting PSM is linked to a row it supports), except in runs listed in provenance `skipped_runs_psm`. A psm-only export (mzIdentML) has no links.

**FID-10 (MUST)** Conflicts: the case of [T05](t05-feature-rows.md) FRW-2 (a signal feeds at most one row) in which several peptidoforms claim one MS1 signal; reporter spectra follow [T02](t02-chimeric-spectra.md) CHI-11. A *conflict* exists when the source assigns the row's signal to several peptidoforms, or every supporting PSM has a rank-1 tie with another peptidoform ([T01](t01-psm-rows.md) PSM-6). Then:
- The signal stays on the source's resolution, else on the peptidoform whose best PSM sorts first under FID-3, over the whole `consensus_feature_id` when there is one.
- That row and every loser row written set `identification_conflict = true`. Losers MAY appear with null intensities (intensities.md:123); otherwise their PSMs keep a null `feature_id`.
- A producer that drops the signal counts the rows in provenance ([T03](t03-filtering-state.md) FLT-6).

`false` = checked, none found; null = not checked. Co-identifications ([T02](t02-chimeric-spectra.md)) and MS1-label partners are not conflicts under (b); one MS1 signal credited to two co-identified peptidoforms is a conflict under (a) ([T02](t02-chimeric-spectra.md) CHI-11).

**FID-11 (SHOULD)** `consensus_feature_id` (int64, nullable; owned here, answering [T05](t05-feature-rows.md) question 3) carries one integer on every run row split from one cross-run entity of the source: an OpenMS ConsensusFeature (map position allowed) or an IonQuant `combined_ion` row; null elsewhere.
- Unique within the file; relayout may renumber it consistently, merge namespaces it, migrate keeps it.
- Unidentified rows ([T07](t07-unidentified-features.md) UNI-3) carry it under the same rule. The rows of one id are either all unidentified or all identified, because a signal that cross-run linking ties to a peptidoform is `transferred`, not unidentified ([T07](t07-unidentified-features.md)).
- Within one id, rows with a non-null intensity share one peptidoform, at most one per `run_file_name`.

**FID-12 (MUST)** `spectral_count` is named and placed by [T09](t09-intensities-missing-values.md) INT-11. In T06 terms it counts the distinct `(run_file_name, scan)` of the row's supporting PSMs that pass the file's filters. It is `0` on transferred rows ([T08](t08-match-between-runs.md) MBR-3) and never counts psm rows, engine hits or summed reporter spectra ([T09](t09-intensities-missing-values.md) `quant_spectrum_count`).

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| `identification_conflict` | feature | bool, nullable | Optional; SHOULD in [T13](t13-metadata-conformance.md) `dda` | No |
| `consensus_feature_id` | feature | int64, nullable; add to `_NAMESPACED_IDS` (`qpx/core/dataset_merge.py:41`), relayout, migrate | Optional | No |
| `skipped_runs_psm` (FID-8) | psm provenance | list of `run_file_name` | When some runs are not linked | No |
| `scan` = one spectrum; check feature scan_format (`qpx/core/data/schema.py:134-139`, [T11](t11-spectrum-precursor-fields.md)); edit `feature.md:26-28,139,301`, `scan.md:90-92,120,123-124`, `feature.yaml:98-101` | feature | `list<int32>` unchanged | – | **Yes**: openms-consensus arrays and derived ids change |
| `id_run_file_name` ∈ {`run_file_name`, null}; drop "(may differ)" (`feature.md:138`) | feature | – | – | **Yes** for pre-[#378](https://github.com/bigbio/qpx/issues/378) OpenMS `-out_qpx` files with donor runs (FID-R7) |
| PEP and scores per FID-4 | feature | – | – | Values change (PD, OpenMS); OpenMS peptide-level PEP null per [T04](t04-confidence-fields.md) |
| Delete the `psm.md:286` sentence (with [T04](t04-confidence-fields.md)); `data-model.md:166` → `id_run_file_name + scan` | docs | – | – | No |
| Validator ([T13](t13-metadata-conformance.md) V5, V6; errors under `dda`, else warnings): FID-1, 5, 6, 7, 10, 11 invariants; FID-8 links and resolution outside runs in `skipped_runs_psm` (warning when provenance `psm_fdr` states a cut). Warning: a linked same-run PSM with lower PEP than the identifying one | psm, feature | – | – | New findings |

## Producer changes
| Producer | Change needed |
|---|---|
| MaxQuant | Match `MS/MS Scan Number` (MQ 1.5). `MULTI-MATCH-MSMS`: empty scan, null run and PEP ([T08](t08-match-between-runs.md) MBR-1, MBR-3). Shared isotope patterns (same `Raw file`, `Charge`, `m/z`, `Intensity`, different `Modified sequence`): intensity on the FID-3 winner, null on the others, all flagged; elsewhere `false`. Check run, charge, peptidoform on `Evidence ID` links. Fix converter-coverage.md:285 |
| FragPipe | Key the lookup by CV; null run when no PSM is found; unit rows per [T05](t05-feature-rows.md) FRW-6, linking member-run PSMs. Fill links; runs without a psm.tsv go to provenance `skipped_runs_psm`; read `<exp> Spectral Count`; `consensus_feature_id` from the `combined_ion` row. FID-3's fallback is today's pick (question 1) |
| Proteome Discoverer | FID-3 for the feature pick, [T01](t01-psm-rows.md) PSM-2 for psm de-duplication. Fill links. `spectral_count` = distinct scans, TMT included ([T09](t09-intensities-missing-values.md)) |
| CDAP | Order compatible (SpecEValue documented). Fill links. `AmbiguousMatch` ties: one flagged row, losers optional. `psm_count` → `spectral_count` or `quant_spectrum_count` ([T09](t09-intensities-missing-values.md) INT-12) |
| OpenMS (qpxc) | Peptidoform per ConsensusFeature: resolver survivor, else FID-3 first over its PSMs. Per run: FID-3 pick of that peptidoform; `scan` = its spectrum (drop feature_dedup.py:48); runs with only competing PSMs per FID-1/FID-10, with `transferred` taken from `mbr_transfer_map_indices` ([T08](t08-match-between-runs.md) MBR-2). Link supporting PSMs only; keep identifying PSMs under `--no-include-unassigned-psms`. Write `consensus_feature_id`, feature `scan_format` |
| OpenMS `-out_qpx` | Merge the [#378](https://github.com/bigbio/qpx/issues/378) branch (per-run identifying PSM, own-run links), so that no row names a donor run (FID-2). Take `transferred` from `mbr_transfer_map_indices` ([T08](t08-match-between-runs.md) MBR-2). Add `consensus_feature_id`, `identification_conflict`. See the [OpenMS worklist](openms-worklist.md). |
| quantms-MSstats | `id_run_file_name` only where `Reference` gives a scan (isobaric); LFQ null |
| DIA-NN, Spectronaut | None (DIA, no psm view) |

## Reader rules
**FID-R1 (MUST)** Find the identifying PSM with `psm.feature_id = feature.feature_id AND psm.run_file_name = feature.id_run_file_name AND psm.scan = feature.scan`; without a link, by `(id_run_file_name, scan, peptidoform, charge)` (not for MS1-label rows). Never `(run, scan)` alone: ambiguous in 440 of 5,843 LFQ and 10,661 of 43,121 TMT PD features ([#378](https://github.com/bigbio/qpx/issues/378)).

**FID-R2 (MUST)** Supporting spectra are the linked PSMs; count distinct `(run_file_name, scan)` among accepted PSMs ([T03](t03-filtering-state.md) FLT-R4). A missing link means "not supporting" in 1.4 files (FID-9), except in runs listed in provenance `skipped_runs_psm`; there, and in older files, use `spectral_count` or report unknown.

**FID-R3 (MUST)** A null `id_run_file_name` or empty `scan` means only "no located identifying PSM"; never infer a transfer from it ([T08](t08-match-between-runs.md) MBR-R1).

**FID-R4 (MUST)** Feature PEP and `psm_qvalue` describe one spectrum, never a peptide or precursor FDR ([T04](t04-confidence-fields.md)).

**FID-R5 (MUST)** Never add intensities of rows sharing `(consensus_feature_id, run_file_name)`. Keep `identification_conflict = true` rows distinguishable (a filter or a reported count); never read a null flag as `false`.

**FID-R6 (MUST)** Read `feature.scan` as one spectrum under `scan_format`; join mz on `(id_run_file_name, scan)`.

**FID-R7 (MUST)** In files with `qpx_version` < 1.4, read a feature `scan` without a feature `scan_format` (openms-consensus) as possibly several spectra and use links; read an `id_run_file_name` other than `run_file_name` as unknown, never as a donor.

Current violations:
- No qpx helper finds the identifying PSM; `link_feature_psm` returns links without role (`qpx/core/dataset_queries.py:315-332`); `psms_without_feature` reads a null link as "not tied to a quantified consensus feature" (`:337-339`), against psm.md:147-150.
- The validator checks links only for existence, `psm_ids` sync and `psm.label ∈ identified_labels` (`qpx/core/label_validation.py:259-305`), never feature `scan_format` or `id_run_file_name` (`dataset_validation.py:248-266` covers `run_file_name` only).
- mokume infers a transfer from a null or foreign `id_run_file_name` when the type is null (mokume:rust/crates/mokume-io/src/qpx.rs:493-524), against FID-R3.
- mokume counts linked `psm_id`s without rank, acceptance or scan de-duplication and ignores unlinked PSMs (mokume:rust/crates/mokume-pipeline/src/spectral_count.rs:158-195), against FID-R2.

## Open questions for maintainers
1. FragPipe: keep the FID-3 fallback (lowest PEP), or rank by Expectation like FragPipe (212 ties differ, [#378](https://github.com/bigbio/qpx/issues/378))?
2. Does a `MULTI-MATCH-MSMS` MS/MS ever identify the row's peptide at PSM FDR? Then [T08](t08-match-between-runs.md) makes it `direct` ([T08](t08-match-between-runs.md) question 1).
3. Should loser rows be MUST rather than MAY, to keep their evidence?
4. Should co-eluting localisation isomers count as conflicts (common in phospho), only sequence-level conflicts, or an mzTab-M-style ambiguity code?
5. Final name of `consensus_feature_id` (open in [#378](https://github.com/bigbio/qpx/issues/378)).
