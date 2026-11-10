# T03. Filtering state: FDR level and threshold, pass_threshold, decoys, contaminants and unassigned PSMs

[← T02](t02-chimeric-spectra.md) · [Overview](index.md) · [Strict profile](strict-profile.md) · [Annex overview](annex-overview.md) · [T04 →](t04-confidence-fields.md)

!!! note "Annex: producer-neutral analysis"
    Part of the annex of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text, and they cover every producer. For OpenMS-written files the [strict profile](strict-profile.md) is the contract and wins where the two differ. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** Producers may keep or drop rows, but provenance states the confidence cuts behind the psm, feature and pg rows, one parameter per level (`psm_fdr`, `precursor_fdr`, `peptide_fdr`, `protein_fdr`); decoy, contaminant and unassigned-PSM handling are provenance parameters too. `pass_threshold` stays the source's own acceptance flag, as mzIdentML passThreshold is, with its basis in the ontology properties `level` and `fdr` of its `ontology.parquet` row. No footer key is added. qpxc still removes no rows by default.

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
| OpenMS -out_qpx | The native writer (verified at HEAD) stamps no filtering key, and `ProteomicsLFQ -out_qpx` and `IsobaricWorkflow -out_qpx` write no provenance view. The upgrade keeps only `creator`, `software_provider`, `scan_format` and `identity_composite` (`qpx/converters/openms/converter.py:209-223`) | [[OPENMS-31](findings-cross.md#openms-31)] |
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
| One footer key `fdr_filter` with a value per level | Travels with the file | Consumers (qpxui, the quantms.io portal) read tables with SQL, not footers; needs carrying through merge, relayout and migrate |
| **Fixed provenance parameters, one per level (`<level>_fdr`)** | One value per level, like mzTab `false_discovery_rate`; cumulative; read with SQL; no footer plumbing | A file copied out of its dataset loses the statement; validation is cross-file; a non-FDR cut collapses to `applied` |
| `pass_threshold` = "passes 1% PSM FDR" | One meaning | Only PD could fill it |
| **`pass_threshold` = source acceptance, basis declared** | Always fillable or null | Readers check the basis |
| **Decoys declared, never removed by qpx** | Feasible; [#342](https://github.com/bigbio/qpx/issues/342) | Readers filter `is_decoy` everywhere |

## Recommendation
**FLT-1 (MUST)** Provenance states the confidence cuts behind the psm, feature and pg rows, one parameter per level:

| Parameter | Step (`step_category`) | Values |
|---|---|---|
| `psm_fdr`, `precursor_fdr`, `peptide_fdr`, `protein_fdr` | `fdr_filtering`, or the search step that cut | a decimal in (0, 1], `applied` or `none` |

An absent parameter means unknown, and so does an off-vocabulary value. `qpxc validate --profile dda` ([T13](t13-metadata-conformance.md)) requires all four parameters. No footer key is added: these parameters, and decoy, contaminant and unassigned-PSM handling (FLT-10–12), are in provenance; the basis of `pass_threshold` is in ontology properties (FLT-8).

**FLT-2 (MUST)** The `<level>_fdr` parameters describe the dataset's rows, cumulatively over every cut by the source and by qpxc. `psm_fdr`, `precursor_fdr` and `peptide_fdr` bound psm and feature rows and the identifications behind pg rows; `protein_fdr` bounds pg rows only ([T12](t12-protein-groups.md) PGV-3). Per level:
- **Decimal x:** a row is present only if its entity at that level (PSM; precursor = peptidoform + charge; peptide, [T04](t04-confidence-fields.md); protein group) had a source-estimated FDR or q-value ≤ x. On feature and pg rows, a finer level limits the identifications behind evidence, counts and quantities ([T12](t12-protein-groups.md) PGV-3); a multiplet is cut per channel. Decoys are cut like targets. x bounds the rows; it does not promise every passing entity is present.
- **`applied`:** a cut whose FDR value is unknown, or a cut on another statistic (score, p-value, PEP, probability, confidence class).
- **`none`:** no confidence cut at that level.
- **Absent:** whether that level was cut is unknown.
- **Several cuts:** when several steps of one source dataset cut a level (a qpxc cut recorded per FLT-4 included), the strictest decimal holds, else `applied`.

The parameters state the levels at which rows were cut, not the levels at which FDR was estimated. A source may estimate FDR at one level only: OpenMS sometimes estimates only peptide-level FDR (for example in immunopeptidomics), or only PSM-level FDR, and the same holds for every producer. The typed q-value column of a level that was not estimated is null on every row of the file ([T04](t04-confidence-fields.md) CNF-2). This is not a defect and not unknown.

**FLT-3 (MUST)** On `transferred` feature rows, the `<level>_fdr` parameters describe the donor identification; the transfer's cut is the `mbr_fdr` parameter of the provenance `match_between_runs` step ([T08](t08-match-between-runs.md)). Unidentified feature rows ([T07](t07-unidentified-features.md) UNI-3) are outside the meaning of the `<level>_fdr` parameters: no identification stands behind them.

**FLT-4** Provenance traces each cut.
- **(MUST)** A qpxc cut is recorded as `qvalue_threshold` and `qvalue_threshold_level` in `format_conversion`. It counts as a cut at that level (FLT-2).
- **(SHOULD)** An upstream cut is an `fdr_filtering` step that carries the FLT-1 parameter of its level; a cut made inside the search tool is stated on the search step. An optional `<level>_fdr_term` gives the PSI-MS term: MS:1002350, MS:1001364, MS:1001214 or MS:1002369.
- **(SHOULD)** The scope of each cut (`run` or `experiment`) goes in the same step as `<level>_fdr_scope`. A reader that needs the scope and finds none follows [T04](t04-confidence-fields.md) CNF-R3.

**FLT-5 (MUST)** qpxc removes no psm, feature or pg rows by default ([#241](https://github.com/bigbio/qpx/issues/241)), whether by confidence, decoy or contaminant status.
- An option that removes rows is recorded in provenance. DIA-NN `--qvalue-threshold t` is recorded as `qvalue_threshold = t` with `qvalue_threshold_level = precursor`; it cuts feature and pg and states no protein cut ([#342](https://github.com/bigbio/qpx/issues/342)). Spectronaut cuts feature only, which a dataset-wide parameter cannot say (open question 2). `--rank1-only` is recorded as `rank1_only` ([T01](t01-psm-rows.md) PSM-5); `--no-include-unassigned-psms` as `include_unassigned_psms = false` (FLT-12).
- An option that was not applied changes no parameter, and provenance records `<option>_applied = false`.

**FLT-6 (MUST)** Rows dropped other than by a declared filter are counted in provenance (`skipped_rows_<view>`) and warned about.

**FLT-7 (MUST)** `pass_threshold` values:
- **`true`:** the source accepted the row under the criterion that the ontology properties `level` and `fdr` of the `pass_threshold` row declare (FLT-8).
- **`false`:** the source evaluated the row and it failed.
- **`null`:** the source states no acceptance: no flag, or the mzIdentML Threshold is MS:1001494 "no threshold" (a `true` without a threshold certifies nothing).

The flag ignores rank, decoy status and `precursor_index` ([T02](t02-chimeric-spectra.md)), so `true` never means "the spectrum's identification" (FLT-R4).

**FLT-8 (MUST)** If any `pass_threshold` is non-null, the `ontology.parquet` row of `pass_threshold` for that view states its basis in two properties ([T13](t13-metadata-conformance.md) META-1): `level` (`psm`, `precursor`, `peptide`, `protein` or `other`; absent means unknown), plus `fdr` (a decimal in (0, 1] or `unknown`) unless the level is `other` or absent. This holds for pg `pass_threshold` too ([T12](t12-protein-groups.md) PGV-9). The mzIdentML SIP `<Threshold>` maps MS:1002350, MS:1002354 and MS:1003337 to `psm`, MS:1001364 and MS:1001868 to `peptide`, and other terms (Mascot:SigThreshold) to `other`. Several FDR terms give the finest level; disagreeing protocols leave `level` absent. The cvParams go into the source's search or validation step in provenance.

**FLT-9** Failing and pre-cut rows:
- **(MUST)** A row the source marks as failing (mzIdentML `passThreshold="false"`, PD Medium or Low) keeps `false`, never null. This replaces `producer-guide.md:11`.
- **(SHOULD)** If the source wrote only rows passing a psm, precursor or peptide cut and flags none (MaxQuant msms, Philosopher psm.tsv), write `true` on every row, with the property `level` = the finest cut's level (`other` if not an FDR) and `fdr` = that level's `<level>_fdr` value (`unknown` for `applied`).

**FLT-10 (MUST)** `decoy_policy`, a parameter of the provenance `format_conversion` or search step, takes these values (absent means unknown):
- **`retained`:** every decoy row the source reported at this grain is present, except those cut under FLT-2. Zero decoys then means the source reported none.
- **`removed`:** no decoy is present.
- **`not_searched`:** no target-decoy search; `is_decoy = false` then means "not decoy-labelled".

qpx converters MUST NOT remove decoys from psm, feature or pg; producers SHOULD retain them. Decoy feature rows are allowed, settling [#229](https://github.com/bigbio/qpx/issues/229)'s open point.

**FLT-11 (MUST)** `contaminant_policy`, a provenance parameter, uses the same values, for rows flagged by the typed `contaminant` ([T10](t10-protein-mapping.md) PRT-6). `retained` does not claim that contaminants were searched. qpx converters MUST NOT remove contaminants.

**FLT-12 (MUST)** Identifications the source assigned to no feature (unassigned PSMs) are stated by the existing psm provenance option `include_unassigned_psms`:
- **`true`:** they are present.
- **`false`:** they were dropped.
- **Not recorded:** the converter has no such option and keeps them (FLT-5), or the source has none: no features (mzIdentML), or every identification linked (MaxQuant msms `Evidence ID`).

Link resolution follows [T06](t06-identifications-per-feature.md) FID-8 (every supporting PSM linked; uncovered runs recorded as provenance `skipped_runs_psm`); no row marker (66281f3).

**FLT-13 (MUST)** Carrying the statements:
- No footer key carries filtering state, so relayout, write-back and migrate have none to carry. A transform that removes rows records the cut in a provenance step (FLT-5).
- A merge keeps the provenance steps of every input, each tagged with its `source_dataset` ([T13](t13-metadata-conformance.md) META-4), so the parameters of FLT-1 and FLT-10–12 survive.
- Carry rule owned here: for a merged dataset the value of a level is the weakest over the inputs: unknown (an input without the parameter) over `none` over `applied` over the largest decimal. Readers compute it from the tagged steps (FLT-R1); the merge rewrites none of them.
- The ontology properties of FLT-8 follow [T13](t13-metadata-conformance.md) META-4: where the inputs differ, the merged `pass_threshold` row drops `level` and writes `fdr = unknown`, with a warning.

**FLT-14 (MUST)** The validator reports as errors: (a) a `<level>_fdr` value outside FLT-1, or a `level` or `fdr` property outside FLT-8; (b) `is_decoy = true` rows under provenance `decoy_policy ∈ {removed, not_searched}`; (c) `contaminant = true` rows under provenance `contaminant_policy = removed`. As warnings: (d) a typed q-value of level L ([T04](t04-confidence-fields.md)) above a decimal `<L>_fdr` (scopes may differ); (e) a non-null `pass_threshold` whose ontology row has no `level` property (an error under the `dda` profile, [T13](t13-metadata-conformance.md) META-8).

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| Footer | – | none new | – | – |
| `psm_fdr`, `precursor_fdr`, `peptide_fdr`, `protein_fdr` | provenance (`fdr_filtering` step, or the search step that cut) | FLT-1 | SHOULD; MUST under `dda` | No |
| `pass_threshold` semantics (FLT-7/9); `psm.yaml:133-136` adds `nullable: true` | psm | bool, null | optional | **Yes for mzIdentML:** rows under SIP "no threshold" become null |
| `producer-guide.md:11` → FLT-9; `pg.md:126-128` → FLT-R4; `psm.md:270-276` → [T13](t13-metadata-conformance.md) META-R1 | – | text | – | No |
| `qvalue_threshold_level`, `<level>_fdr_scope`, `<level>_fdr_term`, `<option>_applied`, `skipped_rows_<view>` | provenance | string | FLT-4/5/6 | No |
| `decoy_policy`, `contaminant_policy` | provenance (`format_conversion` or search step) | FLT-10/11 | SHOULD | No |
| Properties `level`, `fdr` on the `pass_threshold` row | ontology (psm; pg via [T12](t12-protein-groups.md)) | FLT-8 | MUST if `pass_threshold` non-null | No |

## Producer changes
| Producer | Change needed |
|---|---|
| Writers, CLI | `FilterState` like `XicExtraction` (`qpx/writers/xic.py:9-31`), written by every converter as provenance parameters (FLT-1, FLT-10–12) and as the ontology properties of FLT-8. Repeatable `--filter-state KEY=VALUE` lets pipelines (quantms) declare unreadable state, marked user-declared in provenance |
| MaxQuant | `psm_fdr` from parameters.txt `PSM FDR` and `protein_fdr` from `Protein FDR`. Without the file, `applied` (row names to confirm; no fixture has it). Other levels absent (unknown). Provenance decoys and contaminants `retained`. FLT-9 `true` (`level = psm`) |
| FragPipe | `<level>_fdr` from fragpipe.workflow if supplied, else `psm_fdr = applied` with the other levels absent. Provenance `decoy_policy = removed`, or absent if a decoy row appears; contaminants `retained` ([T10](t10-protein-mapping.md) keeps `contam_`). FLT-9 `true` |
| Proteome Discoverer | Ontology properties `level = psm` and `fdr` = strict target FDR only when a target-decoy validator is declared, else neither (the Fixed Value validator uses scores). No `<level>_fdr` parameter (unknown). Provenance decoys `removed` (to confirm) |
| mzIdentML | FLT-7/8. Provenance decoys `retained` if any `isDecoy="true"` or the SearchDatabase has MS:1001197 (target+decoy), else absent |
| OpenMS consensusXML | Provenance decoys `retained`. `<level>_fdr` from header IDFilter/ProteomicsLFQ parameters, else absent. MAY derive `pass_threshold` from a non-zero `significance_threshold` (`level = other`). Count the drops at `psm_adapter.py:348-358` |
| OpenMS -out_qpx | OpenMS writes a provenance view with the `<level>_fdr` and the other parameters; `qpxc convert openms` keeps it (`converter.py:209-223`). Verified at HEAD: `ProteomicsLFQ -out_qpx` and `IsobaricWorkflow -out_qpx` write only feature, psm and pg (no run, sample or provenance view), so none of these statements exists. See the [OpenMS worklist](openms-worklist.md). |
| CDAP, quantms-MSstats | Provenance decoys `retained` if decoy rows are written, else absent. The rest absent unless declared |
| DIA-NN, Spectronaut | FLT-5. DIA-NN main report: `precursor_fdr = applied` (`qpx/cli/convert/diann.py:100-108`) |
| Validator | FLT-14 as dataset-level checks: (a)–(d) read provenance, (a) and (e) read `ontology.parquet` |

## Reader rules
**FLT-R1 (MUST)** Read the `<level>_fdr` parameters and those of FLT-10–12 from provenance, and the basis of `pass_threshold` from its ontology row (FLT-8); absent means unknown. For one source dataset the strictest stated cut of a level holds (FLT-2); for a merged dataset, the weakest over its sources (FLT-13). Never infer filtering state from data: no decoys does not mean removed, and max q ≤ 0.01 does not mean filtered.

**FLT-R2 (MUST)** Always exclude `is_decoy = true` in every view before counting or quantifying. Estimate FDR from targets and decoys only when provenance states `decoy_policy = retained` and the `<level>_fdr` of that level is `none` or a decimal ≥ the estimated threshold.

**FLT-R3 (MUST)** A confidence filter at level L and threshold t uses the typed q-value of L ([T04](t04-confidence-fields.md)). A null passes when provenance states `<L>_fdr` as a decimal ≤ t; otherwise the row is unassessable, and is excluded and counted. When the q-value column of L is null on every row and `<L>_fdr` does not settle it, FDR was not estimated at L ([T04](t04-confidence-fields.md) CNF-2) and filtering at L is not available: the reader says so (a qpx helper raises an error that names the levels with q-values, [T14](t14-reader-contract.md) RDR-3) and filters at a level that exists. It never reads a null as 1.0 and never drops every row silently. Never return a silently empty or unfiltered result.

**FLT-R4 (MUST)** An accepted PSM has `is_decoy = false`, rank 1 ([T01](t01-psm-rows.md) PSM-R1; null rank per PSM-R3), and either `pass_threshold = true` or a null `pass_threshold` with a `psm_fdr`, `precursor_fdr` or `peptide_fdr` that is a decimal or `applied`. `true` equals an FDR filter at t only if the ontology row of `pass_threshold` states an FDR `level` and an `fdr` ≤ t. With a null `pass_threshold` and no declared cut, the reader applies FLT-R3 at a stated threshold, or reports non-decoy rank-1 counts as upper bounds. This predicate governs counts, spectral counts, modification statistics and pg evidence.

**FLT-R5 (SHOULD)** Counts and exported matrices carry the dataset's `<level>_fdr` values. Where a level is absent or `none`, they are upper bounds.

**FLT-R6 (MUST)** Remove contaminants by the typed flag ([T10](t10-protein-mapping.md) PRT-R3), never by accession substring.

**Current violations:**
- `qpx/core/data/psm.py:23-25` reads only `qvalue` (R3, R4), and `feature.py:52-54` COALESCEs levels (R3).
- `qpx/views/api.py:19-29,139-151,172-183` count unfiltered rows (R4, R5), and `:45-48` keeps a null `global_qvalue` (R3).
- `qpx/mudata.py:382-401,758` copies no filtering state into `uns` (R5). Its default decoy exclusion conforms to R2.
- `psm.md:270-276` recommends `schema_arrow.metadata` ([T13](t13-metadata-conformance.md) META-R1).
- mokume reads no filtering state from provenance; it queries only `parquet_kv_metadata` (`mokume:python/mokume/io/_qpx.py:36`; R1).
- mokume spectral counting ignores `pass_threshold` and rank (`mokume:rust/crates/mokume-io/src/psm.rs:14-25`; R4).
- mokume `PeptideFDRFilter` raises only on an all-null column, and otherwise drops null rows silently (`mokume:python/mokume/preprocessing/filters/peptide.py:403-409`; R3).
- mokume removes contaminants and decoys by substring (`mokume:python/mokume/io/feature.py:185-188`; R2, R6).

## Open questions for maintainers
1. *Resolved in the overview:* threshold scope (run vs experiment) stays in provenance as `<level>_fdr_scope` (FLT-4).
2. Spectronaut `--qvalue-threshold` cuts features but pg still comes from every `PG.Quantity`, unlike DIA-NN ([#342](https://github.com/bigbio/qpx/issues/342)). Align? Until then the dataset-wide `qvalue_threshold` cannot say that only the feature view was cut.
3. Proteome Discoverer: do real exports lack decoys, and can the validator node and FDR values be read?
4. CDAP: is its PSM cut documented? The fixture suggests q ≤ 0.01.
5. pepmap population [[MZID_CDAP-M2](findings-cross.md#mzid_cdap-m2)] [[CRITIC-57](findings-cross.md#critic-57)]: [T10](t10-protein-mapping.md) PRT-14 sets it (target peptidoforms of rank-1 psm rows and of feature rows). Do the provenance `<level>_fdr` parameters describe the pepmap rows as well?
6. Should `qpxc validate --profile dda` require the provenance `decoy_policy` and `contaminant_policy`, as it requires the `<level>_fdr` parameters?
