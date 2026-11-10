# T03. Filtering state: FDR level and threshold, pass_threshold, decoys, contaminants and unassigned PSMs

[← T02](t02-chimeric-spectra.md) · [Overview](index.md) · [T04 →](t04-confidence-fields.md)

!!! note "Proposal"
    Part of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** Producers may keep or drop rows, but every psm, feature and pg footer states what survived: one `fdr_filter_<level>` key per level (psm, precursor, peptide, protein), `decoy_policy`, `contaminant_policy`, and on psm `unassigned_psm_policy`. `pass_threshold` stays the source's own acceptance flag, as mzIdentML passThreshold is, with its basis declared in `pass_threshold_level` and `pass_threshold_fdr`. qpxc still removes no rows by default.

## Questions a reader has
- Is this file FDR-filtered? At which level and threshold, on which statistic, by the tool or by qpxc?
- What does `pass_threshold = true` certify, and what do `false` and `null` mean? Is `true` "the spectrum's identification"?
- No decoys: removed, none passed, or no target-decoy search? Are contaminants present?
- Are PSMs that the source assigned to no feature in the psm file?
- Which PSMs count as accepted for pg evidence and counts?

## What happens today
**Spec.**
- `docs/spec/psm.md:66`: `pass_threshold` is the engine's or validator's mark, null "when it does not say"; no threshold is recorded, and `psm.yaml:133-136` lacks `nullable`. `producer-guide.md:11` ("Keep rows above your FDR threshold only if you say so with `pass_threshold`") reads two ways.
- `psm.md:45`, `feature.md:60`: `is_decoy = false` also covers no target-decoy search. `pg.md:99-103`: pg may hold decoys, consumers MUST filter ([#342](https://github.com/bigbio/qpx/issues/342)); psm and feature say nothing on retention. Contaminants are flagged on pg only (`producer-guide.md:34`).
- `pg.md:126-128`: evidence needs a PSM "that passed the FDR threshold", which is undefined.
- No footer table has a filtering key (`psm.md:256-268`, `feature.md:272-284`, `serialization.md:279-296`); precedent `xic.md:21` (`xic_selection`, `xic_qvalue_threshold`).
- `provenance.md:3` promises "what FDR was used?", but only qpxc options are recorded (`:65`); `psm_fdr`/`protein_fdr` are examples (`:159-162,273`).

**Decided issues.**
- [#241](https://github.com/bigbio/qpx/issues/241): no filtering by default, and with a threshold each view filtered at its own level (CHANGELOG.md:358). [#342](https://github.com/bigbio/qpx/issues/342) (f15bd44, merged in PR [#349](https://github.com/bigbio/qpx/issues/349); CHANGELOG.md:265) replaced the second half: DIA-NN `--qvalue-threshold` selects precursors by `Q.Value` for feature and pg, and pg q-values are never cut. [#342](https://github.com/bigbio/qpx/issues/342) also keeps decoy groups. [#229](https://github.com/bigbio/qpx/issues/229) left feature decoys open.
- [#303](https://github.com/bigbio/qpx/issues/303) dropped OpenMS unassigned PSMs by default, and [#307](https://github.com/bigbio/qpx/issues/307) restored them. Its row marker was then dropped as derivable (66281f3).

| Producer | Behaviour | Findings |
|---|---|---|
| MaxQuant | msms and evidence cut at PSM FDR, proteinGroups at protein FDR; Reverse and `CON__` rows kept; `pass_threshold` and `qvalue` null on 1,000/1,000 msms rows; parameters.txt read only for fixed mods (`qpx/converters/maxquant/constants.py:137-144`) | [[MAXQUANT-05](findings-psm.md#maxquant-05)] [[MAXQUANT-06](findings-cross.md#maxquant-06)] [[MAXQUANT-22](findings-cross.md#maxquant-22)] |
| FragPipe | psm.tsv and combined_protein.tsv are cut by Philosopher and decoy-free (`pg.md:365`). `pass_threshold` is null. Provenance has no filtering step (`qpx/converters/fragpipe/converter.py:95-120`) | [[FRAGPIPE-01](findings-psm.md#fragpipe-01)] [[FRAGPIPE-08](findings-psm.md#fragpipe-08)] [[FRAGPIPE-29](findings-cross.md#fragpipe-29)] |
| Proteome Discoverer | `pass_threshold = (Confidence == 'High')`, null without the column (`qpx/converters/proteome_discoverer/psm_adapter.py:68`). Export filters and validator node are unknown. `is_decoy` is always false (`:60`) | [[PD-08](findings-psm.md#pd-08)] [[PD-09](findings-cross.md#pd-09)] [[PD-23](findings-cross.md#pd-23)] |
| mzIdentML | SII `passThreshold` is copied (`qpx/converters/mzidentml/psm_adapter.py:304`; missing → false, only in schema-invalid files). SIP `<Threshold>` is never read; fixtures hold Mascot:SigThreshold 0.01 homology, MS:1003337 crosslinked PSM FDR 0.05, and MS:1001494. PXD054720: 11,517/13,386 true, 1,810 of them rank ≥ 2; no `qvalue` | [[MZID_CDAP-10](findings-psm.md#mzid_cdap-10)] [[INFRA-M3](findings-psm.md#infra-m3)] [[CRITIC-M2](findings-psm.md#critic-m2)] |
| OpenMS consensusXML | No cut: decoys, all hits and unassigned PIDs are written (`qpx/converters/openms_consensus/converter.py:140`). `significance_threshold` is unread. PIDs lacking a run or spectrum reference vanish silently (`psm_adapter.py:348-358`) | [[OPENMS-03](findings-psm.md#openms-03)] [[OPENMS-22](findings-cross.md#openms-22)] [[OPENMS-M4](findings-psm.md#openms-m4)] |
| OpenMS -out_qpx | The upgrade keeps only `creator`, `software_provider`, `scan_format` and `identity_composite` (`qpx/converters/openms/converter.py:209-223`) | [[OPENMS-31](findings-cross.md#openms-31)] |
| CDAP | Looks pre-cut (max Qvalue 0.0089, 10/100 decoys); nothing recorded | [[MZID_CDAP-11](findings-psm.md#mzid_cdap-11)] |
| qpxc filters | DDA converters have none. DIA-NN `--qvalue-threshold` cuts `Q.Value < t` in feature, pg aggregation, fragment and pepmap (`qpx/converters/diann/pg_adapter.py:461-473`), and pepmap drops decoys (`feature_sql.py:250`). Spectronaut cuts feature only, and nothing without `EG.Qvalue` (`qpx/converters/spectronaut/feature_adapter.py:386-390`), yet provenance records the threshold (`qpx/cli/convert/postprocess.py:189-202`) | [[INFRA-32](findings-cross.md#infra-32)] |

**Infrastructure.** Writers stamp no filtering key (`qpx/writers/base.py:114-135`) [[INFRA-20](findings-cross.md#infra-20)]; late keys miss `ARROW:schema` (`:359-365`) [[INFRA-21](findings-cross.md#infra-21)]. The validator ignores filtering (`qpx/core/data/footer_validation.py:77-92`) [[VALIDATOR-31](findings-cross.md#validator-31)]. Of the declaration keys, merge and migrate keep only XIC keys (`qpx/core/dataset_merge.py:252-263`, `qpx/transforms/migrate.py:195-197`); relayout keeps all (`qpx/transforms/relayout.py:136-143`). Failed rows are only logged (`qpx/converters/base.py:211-240`) [[INFRA-30](findings-cross.md#infra-30)].

## Why it matters
- **Same schema, different populations:** MaxQuant psm holds PSM-FDR survivors and passing reverse hits; Mascot mzid 3,679 rank > 1 rows; OpenMS every hit, decoy and unassigned PID.
- **`confident()` returns nothing for pre-filtered files:** 0/1,000 MaxQuant PSMs, 0/13,386 Mascot PSMs and 0/243,821 MaxQuant features [[VALIDATOR-06](findings-psm.md#validator-06)] [[VALIDATOR-07](findings-feature.md#validator-07)].
- **`pass_threshold = true` has three meanings:** Mascot homology p < 0.01 (rank 2 included), PD High confidence, 5% crosslink FDR [[INFRA-25](findings-psm.md#infra-25)].
- **Decoys:** FDR re-estimation and atlas protein FDR ([#342](https://github.com/bigbio/qpx/issues/342)) need them; a decoy-free file looks like one whose decoys failed [[EXTERNAL-31](findings-psm.md#external-31)] [[CRITIC-09](findings-cross.md#critic-09)].
- **No defined "passed".** pg evidence, spectral counts (`mokume:rust/crates/mokume-pipeline/src/spectral_count.rs:158-195`) and `ModificationView` (4,323 vs 926 Deamidated PSMs [[VALIDATOR-M3](findings-psm.md#validator-m3)]) need one.

## Options considered
| Option | For | Against |
|---|---|---|
| Provenance text only | No change | Unvalidated; mzTab and mzIdentML declare it per file [[EXTERNAL-03](findings-cross.md#external-03)] |
| qpxc enforces 1% | Uniform | Reverses [#241](https://github.com/bigbio/qpx/issues/241)/#342; loses rescoring inputs |
| `fdr_filter_level` / `_threshold` / `_metric` (registry draft) | 3 keys | One level per file; MaxQuant and Philosopher cut at several |
| **`fdr_filter_<level>`** | One value per level, like mzTab `false_discovery_rate`; cumulative | 4 keys; a non-FDR cut collapses to `applied` |
| `pass_threshold` = "passes 1% PSM FDR" | One meaning | Only PD could fill it; departs from PSI |
| **`pass_threshold` = source acceptance, basis declared** | PSI-aligned; always fillable or null | Readers check the basis |
| **Decoys declared, never removed by qpx** | Feasible; [#342](https://github.com/bigbio/qpx/issues/342) | Readers filter `is_decoy` everywhere |

## Recommendation
**FLT-1 (MUST)** psm, feature and pg footers use only these values:

| Key | Views | Values |
|---|---|---|
| `fdr_filter_psm`, `fdr_filter_precursor`, `fdr_filter_peptide`, `fdr_filter_protein` | psm, feature, pg | decimal in (0, 1], `applied`, `none`, `unknown` |
| `decoy_policy`, `contaminant_policy` | psm, feature, pg | `retained`, `removed`, `not_searched`, `unknown` |
| `unassigned_psm_policy` | psm | `retained`, `removed`, `not_applicable`, `unknown` |
| `pass_threshold_level` | psm; pg with `pass_threshold` ([T12](t12-protein-groups.md) PGV-9) | `psm`, `precursor`, `peptide`, `protein`, `other`, `unknown` |
| `pass_threshold_fdr` | as above | decimal in (0, 1], `unknown` |

An absent or off-vocabulary value means `unknown`. The [T13](t13-metadata-conformance.md) profile requires every key of the view (`unknown` allowed), but the `pass_threshold_*` keys only when `pass_threshold` is non-null.

**FLT-2 (MUST)** `fdr_filter_<level>` describes this file's rows, cumulatively over every cut by the source and by qpxc.
- **Decimal x:** a row is present only if its entity at that level (PSM; precursor = peptidoform + charge; peptide, [T04](t04-confidence-fields.md); protein group) had a source-estimated FDR or q-value ≤ x. On feature and pg rows, a finer level limits the identifications behind evidence, counts and quantities ([T12](t12-protein-groups.md) PGV-3); a multiplet is cut per channel. Decoys are cut like targets. x bounds the rows; it does not promise every passing entity is present.
- **`applied`:** a cut whose FDR value is unknown, or a cut on another statistic (score, p-value, PEP, probability, confidence class).
- **`none`:** no confidence cut at that level.
- **Several cuts:** write the strictest decimal, else `applied`.

**FLT-3 (MUST)** On `transferred` feature rows, `fdr_filter_*` describes the donor identification; the transfer's cut is [T08](t08-match-between-runs.md) `mbr_fdr`. Rows of [T07](t07-unidentified-features.md)'s `unidentified_feature` view are outside these keys.

**FLT-4** Provenance traces each cut.
- **(MUST)** A qpxc cut is recorded as `qvalue_threshold` and `qvalue_threshold_level` in `format_conversion`.
- **(SHOULD)** An upstream cut is an `fdr_filtering` step with `psm_fdr`, `precursor_fdr`, `peptide_fdr` and `protein_fdr`. An optional `<level>_fdr_term` gives the PSI-MS term: MS:1002350, MS:1001364, MS:1001214 or MS:1002369.
- **(SHOULD)** The scope of each cut (`run` or `experiment`, [T04](t04-confidence-fields.md) CNF-3 vocabulary) goes in the same step as `<level>_fdr_scope`. No footer key carries it in 1.4. A reader that needs the scope and finds none follows [T04](t04-confidence-fields.md) CNF-R3.

**FLT-5 (MUST)** qpxc removes no psm, feature or pg rows by default ([#241](https://github.com/bigbio/qpx/issues/241)), whether by confidence, decoy or contaminant status.
- An option that removes rows updates the keys of every view it changed. DIA-NN `--qvalue-threshold t` sets `fdr_filter_precursor = t` on feature and pg and leaves `fdr_filter_protein` alone ([#342](https://github.com/bigbio/qpx/issues/342)); Spectronaut sets it on feature only. `--rank1-only` sets [T01](t01-psm-rows.md) `psm_rank_policy`; `--no-include-unassigned-psms` sets `unassigned_psm_policy = removed`.
- An option that was not applied changes no key, and provenance records `<option>_applied = false`.

**FLT-6 (MUST)** Rows dropped other than by a declared filter are counted in provenance (`skipped_rows_<view>`) and warned about.

**FLT-7 (MUST)** `pass_threshold` values:
- **`true`:** the source accepted the row under the criterion that `pass_threshold_level` and `pass_threshold_fdr` declare.
- **`false`:** the source evaluated the row and it failed.
- **`null`:** the source states no acceptance: no flag, or the mzIdentML Threshold is MS:1001494 "no threshold" (a deliberate deviation: mzIdentML requires the attribute, but a `true` without a threshold certifies nothing).

The flag ignores rank, decoy status and `precursor_index` ([T02](t02-chimeric-spectra.md)), so `true` never means "the spectrum's identification" (FLT-R4).

**FLT-8 (MUST)** If any `pass_threshold` is non-null, the footer carries `pass_threshold_level`, plus `pass_threshold_fdr` unless the level is `other` or `unknown`. The mzIdentML SIP `<Threshold>` maps MS:1002350, MS:1002354 and MS:1003337 to `psm`, MS:1001364 and MS:1001868 to `peptide`, and other terms (Mascot:SigThreshold) to `other`. Several FDR terms give the finest level; disagreeing protocols give `unknown`. The cvParams go into the search step's provenance.

**FLT-9** Failing and pre-cut rows:
- **(MUST)** A row the source marks as failing (mzIdentML `passThreshold="false"`, PD Medium or Low) keeps `false`, never null. This replaces `producer-guide.md:11`.
- **(SHOULD)** If the source wrote only rows passing a psm, precursor or peptide cut and flags none (MaxQuant msms, Philosopher psm.tsv), write `true` on every row, with the finest cut's level (`other` if not an FDR) and its `fdr_filter_<level>` value (`unknown` for `applied`).

**FLT-10 (MUST)** `decoy_policy` values:
- **`retained`:** every decoy row the source reported at this grain is present, except those cut under FLT-2. Zero decoys then means the source reported none.
- **`removed`:** no decoy is present.
- **`not_searched`:** no target-decoy search; `is_decoy = false` then means "not decoy-labelled".

qpx converters MUST NOT remove decoys from psm, feature or pg; producers SHOULD retain them. Decoy feature rows are allowed, settling [#229](https://github.com/bigbio/qpx/issues/229)'s open point.

**FLT-11 (MUST)** `contaminant_policy` uses the same values, for rows flagged by the typed `contaminant` ([T10](t10-protein-mapping.md) PRT-6). `retained` does not claim that contaminants were searched. qpx converters MUST NOT remove contaminants.

**FLT-12 (MUST)** `unassigned_psm_policy` values:
- **`retained`:** identifications the source assigned to no feature are present.
- **`removed`:** they were dropped.
- **`not_applicable`:** the source has none: no features (mzIdentML), or every identification linked (MaxQuant msms `Evidence ID`).

Link resolution is [T06](t06-identifications-per-feature.md) `psm_feature_link`; no row marker (66281f3).

**FLT-13 (MUST)** Carrying the keys:
- Writers keep `ARROW:schema` equal to the key-value footer ([T13](t13-metadata-conformance.md) META-3).
- Relayout, merge, migrate and the OpenMS `-out_qpx` upgrade carry the keys through [T13](t13-metadata-conformance.md) META-4's shared function; a transform that removes rows updates them (FLT-5).
- Carry rule owned here: a merge writes the weakest `fdr_filter_*` value, `unknown` (or absent) over `none` over `applied` over the largest decimal. The other keys follow [T13](t13-metadata-conformance.md)'s default rule (kept when all sources agree, else `unknown` with a warning).

**FLT-14 (MUST)** The validator reports as errors: (a) values outside FLT-1; (b) `is_decoy = true` rows under `decoy_policy ∈ {removed, not_searched}`; (c) `contaminant = true` rows under `contaminant_policy = removed`. As warnings: (d) a typed q-value of level L ([T04](t04-confidence-fields.md)) above a decimal `fdr_filter_L` (scopes may differ); (e) a non-null `pass_threshold` without `pass_threshold_level` (an error under the `dda` profile, [T13](t13-metadata-conformance.md) META-8).

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| `fdr_filter_psm/_precursor/_peptide/_protein` | psm, feature, pg | FLT-1 | SHOULD; MUST under `dda` | No |
| `decoy_policy`, `contaminant_policy` | psm, feature, pg | FLT-10/11 | SHOULD; MUST under `dda` | No |
| `unassigned_psm_policy` | psm | FLT-12 | SHOULD; MUST under `dda` | No |
| `pass_threshold_level`, `pass_threshold_fdr` | psm (pg via [T12](t12-protein-groups.md)) | FLT-8 | MUST if `pass_threshold` non-null | No |
| `pass_threshold` semantics (FLT-7/9); `psm.yaml:133-136` adds `nullable: true` | psm | bool, null | optional | **Yes for mzIdentML:** rows under SIP "no threshold" become null |
| `producer-guide.md:11` → FLT-9; `pg.md:126-128` → FLT-R4; `psm.md:270-276` → FLT-R1 | – | text | – | No |
| `qvalue_threshold_level`, `<level>_fdr`, `<level>_fdr_term`, `<option>_applied`, `skipped_rows_<view>` | provenance | string | FLT-4/5/6 | No |

## Producer changes
| Producer | Change needed |
|---|---|
| Writers, CLI | `FilterState` like `XicExtraction` (`qpx/writers/xic.py:9-31`), stamped by every writer. Repeatable `--filter-state KEY=VALUE` lets pipelines (quantms) declare unreadable state, marked user-declared in provenance |
| MaxQuant | `fdr_filter_psm` (all views) from parameters.txt `PSM FDR`; pg `fdr_filter_protein` from `Protein FDR`. Without the file, `applied` (row names to confirm; no fixture has it). Other levels `unknown`. Decoys and contaminants `retained`. FLT-9 `true` (`psm`). Unassigned `not_applicable` |
| FragPipe | `fdr_filter_*` from fragpipe.workflow if supplied, else psm `applied`, rest `unknown`. `decoy_policy = removed`, or `unknown` if a decoy row appears. Contaminants `retained` ([T10](t10-protein-mapping.md) keeps `contam_`). FLT-9 `true`. Unassigned `retained` |
| Proteome Discoverer | `pass_threshold_level = psm` and `pass_threshold_fdr` = strict target FDR only when a target-decoy validator is declared, else both `unknown` (the Fixed Value validator uses scores). `fdr_filter_*` `unknown`. Decoys `removed` (to confirm). Unassigned `retained` |
| mzIdentML | FLT-7/8. Decoys `retained` if any `isDecoy="true"` or the SearchDatabase has MS:1001197 (target+decoy), else `unknown`. Unassigned `not_applicable` |
| OpenMS consensusXML | Decoys `retained`. Unassigned follows the option. `fdr_filter_*` from header IDFilter/ProteomicsLFQ parameters, else `unknown`. MAY derive `pass_threshold` from a non-zero `significance_threshold` (`other`). Count the drops at `psm_adapter.py:348-358` |
| OpenMS -out_qpx | OpenMS stamps the keys; `qpxc convert openms` copies them (`converter.py:209-223`) |
| CDAP, quantms-MSstats | Decoys `retained` if decoy rows are written, else `unknown`. The rest `unknown` unless declared. CDAP unassigned `retained` |
| DIA-NN, Spectronaut | FLT-5. DIA-NN main report: `fdr_filter_precursor = applied` (`qpx/cli/convert/diann.py:100-108`) |
| Validator | FLT-14 in `footer_validation.py` |

## Reader rules
**FLT-R1 (MUST)** Read the keys per [T13](t13-metadata-conformance.md) META-R1–R3 (key-value footer, per file; absent means `unknown`). Never infer filtering state from data: no decoys does not mean removed, and max q ≤ 0.01 does not mean filtered.

**FLT-R2 (MUST)** Exclude `is_decoy = true` in every view before counting or quantifying. Estimate FDR from targets and decoys only when `decoy_policy = retained` and that level's cut is `none` or a decimal ≥ the estimated threshold.

**FLT-R3 (MUST)** A confidence filter at level L and threshold t uses the typed q-value of L ([T04](t04-confidence-fields.md)). A null passes when `fdr_filter_L` is a decimal ≤ t; otherwise the row is unassessable, is excluded and counted, and the reader raises if every row is ([T14](t14-reader-contract.md) RDR-3). Never return a silently empty or unfiltered result.

**FLT-R4 (MUST)** An accepted PSM has `is_decoy = false`, rank 1 ([T01](t01-psm-rows.md) PSM-R1; null rank per PSM-R3), and either `pass_threshold = true` or a null `pass_threshold` with `fdr_filter_psm`, `_precursor` or `_peptide` a decimal or `applied`. `true` equals an FDR filter at t only if `pass_threshold_level` is an FDR level and `pass_threshold_fdr` ≤ t. With a null `pass_threshold` and no declared cut, the reader applies FLT-R3 at a stated threshold, or reports non-decoy rank-1 counts as upper bounds. This predicate governs counts, spectral counts, modification statistics and pg evidence.

**FLT-R5 (SHOULD)** Counts and exported matrices carry the file's filter keys. Under `unknown` or `none` they are upper bounds.

**FLT-R6 (MUST)** Remove contaminants by the typed flag ([T10](t10-protein-mapping.md) PRT-R3), never by accession substring.

**Current violations:**
- `qpx/core/data/psm.py:23-25` reads only `qvalue` (R3, R4), and `feature.py:52-54` COALESCEs levels (R3).
- `qpx/views/api.py:19-29,139-151,172-183` count unfiltered rows (R4, R5), and `:45-48` keeps a null `global_qvalue` (R3).
- `qpx/mudata.py:382-401,758` copies no filter key into `uns` (R5). Its default decoy exclusion conforms to R2.
- `psm.md:270-276` recommends `schema_arrow.metadata` (R1). Merge and migrate drop the keys (FLT-13).
- mokume reads no filter key, although it already queries `parquet_kv_metadata` (`mokume:python/mokume/io/_qpx.py:36`).
- mokume spectral counting ignores `pass_threshold` and rank (`mokume:rust/crates/mokume-io/src/psm.rs:14-25`; R4).
- mokume `PeptideFDRFilter` raises only on an all-null column, and otherwise drops null rows silently (`mokume:python/mokume/preprocessing/filters/peptide.py:403-409`; R3).
- mokume removes contaminants and decoys by substring (`mokume:python/mokume/io/feature.py:185-188`; R2, R6).

## Open questions for maintainers
1. *Resolved in the overview:* threshold scope (run vs experiment) stays in provenance as `<level>_fdr_scope` (FLT-4). A footer key is added only if a reader needs it.
2. Spectronaut `--qvalue-threshold` cuts features but pg still comes from every `PG.Quantity`, unlike DIA-NN ([#342](https://github.com/bigbio/qpx/issues/342)). Align?
3. Proteome Discoverer: do real exports lack decoys, and can the validator node and FDR values be read?
4. CDAP: is its PSM cut documented? The fixture suggests q ≤ 0.01.
5. pepmap population [[MZID_CDAP-M2](findings-cross.md#mzid_cdap-m2)] [[CRITIC-57](findings-cross.md#critic-57)]: [T10](t10-protein-mapping.md) PRT-14 sets it (target peptidoforms of rank-1 psm rows and of feature rows). Should pepmap also carry these keys, or inherit those of the psm file?
