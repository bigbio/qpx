# T08. Match-between-runs and ID transfer versus identification in the same run

[← T07](t07-unidentified-features.md) · [Overview](index.md) · [Strict profile](strict-profile.md) · [Annex overview](annex-overview.md) · [T09 →](t09-intensities-missing-values.md)

!!! note "Annex: producer-neutral analysis"
    Part of the annex of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text, and they cover every producer. For OpenMS-written files the [strict profile](strict-profile.md) is the contract and wins where the two differ. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** `identification_type` alone says where a row's identity came from: `direct`/`secondary` = an accepted identification from a spectrum of this run (or unit); `transferred` = a signal measured here whose identity came only from other runs, whatever MS/MS exists here. Transferred rows carry no PSM and no run-scope confidence; their own confidence lives only in `match_*` (acceptor − donor, aligned clock). pg gains `identified_in_unit`. Type and `identified_in_unit` are filled whenever the source states them (null = unknown); whether transfer ran, and its FDR (`mbr_fdr`), go in the provenance `match_between_runs` step, not the footer. Readers never count transfers as identifications or protein evidence, nor combine them with an own-run quantity of the same precursor.

## Questions a reader has
- How do I tell "identified by MS/MS in this run" from "quantified through a transfer"? What does a null type mean?
- What are a MaxQuant `MULTI-MATCH-MSMS` row and a transferred donor-run second peptide?
- On a transferred row: is `id_run_file_name` the donor or null? What do `scan`, PEP, q-values, scores and `spectral_count` hold?
- What do the `match_*` fields mean (sign, units, donor)?
- Was the transfer FDR-controlled (IonQuant, OpenMS PIP-ECHO) or an uncontrolled list?
- One precursor is direct and transferred in one run: which quantity counts?
- Is a transferred value imputed (Spectronaut)? Are transfer-only pg rows protein evidence (andes#218)?
- Were transfers kept, dropped or never made? How do I weight them (mokume#117)?

## What happens today
**Spec.**
- `docs/spec/feature.md:145` defines `direct`, `transferred`, `secondary` and null ("tool does not say"), and a null `id_run_file_name` on transferred rows; `feature.md:138` calls it the best PSM's run that "may differ from run_file_name" (`feature.yaml:207-209`), inviting a donor run [[SPEC-04](findings-feature.md#spec-04)] [[CRITIC-27](findings-feature.md#critic-27)] [[VALIDATOR-M5](findings-feature.md#validator-m5)].
- `match_*` are the tool's "score and q-value for the match" and the "difference from its donor", with no sign, clock or donor (`feature.md:146-147`) [[MAXQUANT-11](findings-feature.md#maxquant-11)] [[EXTERNAL-11](findings-feature.md#external-11)]; `match_qvalue` is float32 without range (`feature.yaml:220-224`) [[SPEC-12](findings-feature.md#spec-12)].
- `docs/spec/producer-guide.md:26` asks for an empty `scan`; nothing covers PEP, q-values or scores on transferred rows [[SPEC-09](findings-feature.md#spec-09)] [[MOKUME-M3](findings-feature.md#mokume-m3)], or transfer versus imputation (`docs/spec/intensities.md:131`) [[SPEC-08](findings-feature.md#spec-08)] [[CRITIC-26](findings-feature.md#critic-26)].
- Only the vocabulary is checked (`qpx/core/data/validation.py:71`): a warning on `write_table` (`qpx/writers/base.py:416`), nothing on `write_batch` (`:599-617`), no YAML enum [[VALIDATOR-09](findings-feature.md#validator-09)] [[INFRA-11](findings-feature.md#infra-11)] [[VALIDATOR-M4](findings-cross.md#validator-m4)].
- `docs/spec/pg.md:131-135` allows MBR-only unit rows, marked only by a MaxQuant cv_param (`pg.md:291`) [[SPEC-25](findings-pg.md#spec-25)] [[CRITIC-49](findings-pg.md#critic-49)]. `docs/spec/anndata.md:153` sums an MS/MS and an MBR feature of one precursor by default [[EXTERNAL-24](findings-feature.md#external-24)] [[VALIDATOR-19](findings-feature.md#validator-19)].
- No footer key declares the transfer state [[INFRA-20](findings-cross.md#infra-20)] [[VALIDATOR-31](findings-cross.md#validator-31)] [[CRITIC-60](findings-cross.md#critic-60)]. PSI-MS and mzTab 1.0 have no MBR term [[EXTERNAL-23](findings-feature.md#external-23)].

| Producer | Behaviour | Findings |
|---|---|---|
| MaxQuant | `Type` mapped (qpx/converters/maxquant/evidence_fields.py:10-17); run nulled only for exact `MULTI-MATCH` (qpx/converters/maxquant/feature_adapter.py:655-657), so `MULTI-MATCH-MSMS` rows are `transferred` with run, scan and PEP (injected test; no fixture). Unmapped Types: null, untraced (`cv_params` None, :474). maxquant_full: 104,491/243,821 rows transferred, all null run, scan, PEP; `match_score` = the donor's Andromeda score (coxdocs; ≤ best donor `Score` in 99.4%); `Match q-value` empty; `match_rt_delta` already ×60 (evidence_fields.py:22-28). Dedup demotes only `MULTI-MATCH` (feature_adapter.py:264). pg: 5,687 "By matching" rows marked only by cv_params (pg_adapter.py:403-406) | [[MAXQUANT-01](findings-feature.md#maxquant-01)] [[MAXQUANT-02](findings-feature.md#maxquant-02)] [[MAXQUANT-11](findings-feature.md#maxquant-11)] [[MAXQUANT-20](findings-pg.md#maxquant-20)] [[MAXQUANT-M2](findings-feature.md#maxquant-m2)] [[SPEC-M1](findings-feature.md#spec-m1)] [[EXTERNAL-10](findings-feature.md#external-10)] [[EXTERNAL-11](findings-feature.md#external-11)] |
| FragPipe | `Match Type` `MS/MS` → direct, `MBR` → transferred, else null (qpx/converters/fragpipe/feature_adapter.py:51,361); transferred rows: null run, scan, PEP; no `match_*` or IonQuant FDR. MBR-only pg rows: spectral count 0, no flag | [[FRAGPIPE-14](findings-feature.md#fragpipe-14)] [[FRAGPIPE-15](findings-feature.md#fragpipe-15)] [[FRAGPIPE-23](findings-pg.md#fragpipe-23)] [[FRAGPIPE-29](findings-cross.md#fragpipe-29)] [[EXTERNAL-12](findings-feature.md#external-12)] |
| OpenMS consensus | `transferred` = no own-run PID and every PID's run known (qpx/converters/openms_consensus/feature_adapter.py:542-547); OpenMS's own statement, the ConsensusFeature meta value `mbr_transfer_map_indices` (with `mbr_transfer_qvalues` under PIP-ECHO), is read nowhere. A run whose own competing PID IDConflictResolver removed also becomes `transferred`, with no conflict marker (fx/idconflict). Same-peak copies merge onto the direct evidence (feature_dedup.py:33-50) | [[OPENMS-11](findings-feature.md#openms-11)] [[OPENMS-12](findings-feature.md#openms-12)] |
| Proteome Discoverer, CDAP | Hard-coded `direct` (qpx/converters/proteome_discoverer/feature_adapter.py:139; qpx/converters/cdap/feature_adapter.py:223-224), correct for PSM-built rows. PD pg may include Feature Mapper quantities (unverified); `Found in Sample` unread | [[PD-11](findings-cross.md#pd-11)] [[MZID_CDAP-22](findings-feature.md#mzid_cdap-22)] |
| quantms-MSstats | Null type; `id_run_file_name = run_file_name` (qpx/converters/quantms_msstats/feature_adapter.py:525), though ProteomicsLFQ requantifies runs without own IDs | [[INFRA-27](findings-feature.md#infra-27)] [[INFRA-M2](findings-feature.md#infra-m2)] |
| Spectronaut | `EG.IsImputed` → `transferred`, `imputed_labels`, null run (qpx/converters/spectronaut/feature_adapter.py:337-343) | [[SPEC-08](findings-feature.md#spec-08)] [[MOKUME-20](findings-feature.md#mokume-20)] |

**Re-checked** (audit run on the maxquant_full fixture): 9,337 (peptidoform, charge, run) keys, decoys included, have measured direct and transferred rows (9,272 target keys; 9,639 with `secondary`); transferred is larger in 799 (867) [[MOKUME-09](findings-feature.md#mokume-09)]. `Match time difference` rises with acceptor − median other-run donor calibrated RT (Spearman +0.45, n = 98,546; +0.62 where the acceptor run has no own ID), `Match m/z difference` with `Mass error [Da]` (+0.71): acceptor − donor, aligned clock.

## Why it matters
- Transfers are 43% of MaxQuant feature rows; per-run counts include them (+62% peptides) [[VALIDATOR-20](findings-feature.md#validator-20)].
- One word covers MBR, imputation (Spectronaut), a heuristic (OpenMS) and silence (quantms).
- IonQuant cuts ion FDR; MaxQuant and FlashLFQ v1.0 report uncontrolled lists. In one PIP-ECHO benchmark the FDP at nominal 1% was 3.4% (IonQuant), 2.0% (MaxQuant), 6.1% (FlashLFQ) [[EXTERNAL-13](findings-feature.md#external-13)]. Readers cannot see the regime.
- mokume#117: 1.9–2.4× more MaxLFQ false positives with transfers. Score filters silently drop all transfers [[MOKUME-M3](findings-feature.md#mokume-m3)]; MAX can pick the MBR copy [[MOKUME-09](findings-feature.md#mokume-09)]. Transfer-only pg rows look like protein evidence, against andes#218.

## Options considered
| Option | For | Against |
|---|---|---|
| Infer transfer from null run, empty scan or null PEP | No column | Null for other reasons [[MOKUME-07](findings-feature.md#mokume-07)]; [T06](t06-identifications-per-feature.md) FID-R3 forbids |
| **Typed `identification_type` defined by accepted own-run identification** | One test for all producers | Per-tool mapping |
| `MULTI-MATCH-MSMS` = `direct` or a new value | Has an MS/MS | Typed MATCH, so presumably that MS/MS did not identify it (undocumented, [EXTERNAL-10](findings-feature.md#external-10)); P-I |
| **`MULTI-MATCH-MSMS` = `transferred`, raw Type kept** | One meaning | Scan moves to cv_params |
| `id_run_file_name` = donor | Provenance | Two meanings; [T06](t06-identifications-per-feature.md) FID-6; donor rarely known |
| **`match_score` = transfer score only; donor scores to `additional_scores`** | P-A | MaxQuant values move |
| Producer drops transferred duplicates | Clean | Silent; split peaks lost |
| **Keep both; readers apply in-run precedence** | No column | Every reader implements it |
| pg `identification_type` vs **pg `identified_in_unit`** | Mirrors feature / exactly pg.md:126-135 | Mixed quantities make "direct" overclaim / new column |

## Recommendation
*Accepted own-run identification*: a PSM (peptide-centric tools: a run-level precursor identification) of the row's peptidoform and charge, from a spectrum of `run_file_name` or of a run of the row's unit ([T05](t05-feature-rows.md) FRW-6), accepted by the source ([T03](t03-filtering-state.md)). A *transfer* (MBR, ID transfer, requantification, PIP-ECHO) assigns a peptidoform to this run's signal using only other runs' identifications (*donors*).

**MBR-1 (MUST)** `identification_type` (enum per [T07](t07-unidentified-features.md) UNI-2; other values are errors on every write path):
- `direct`: an accepted own-run identification supports the signal, and not all are `secondary`;
- `secondary`: [T02](t02-chimeric-spectra.md) CHI-9;
- `transferred`: measured in this run, no accepted own-run identification of its peptidoform, whatever unidentified, sub-threshold or differently identified MS/MS the run has. `MULTI-MATCH-MSMS` and a transferred donor-run second peptide are `transferred`;
- `unidentified`: an MS1 feature without any identification, with a null `peptidoform` ([T07](t07-unidentified-features.md) UNI-1, UNI-3);
- null: an identified row of unknown origin.

Any accepted own-run identification wins over transfer.

**MBR-2 (MUST)** The value comes from the tool's per-row statement: MaxQuant `Type`, IonQuant `Match Type`, OpenMS `mbr_transfer_map_indices`, FlashLFQ `Peak Detection Type`.
- `mbr_transfer_map_indices` is a ConsensusFeature meta value that lists the runs the identification was transferred to. It is written after FDR filtering by IDConflictResolver / ProteomicsLFQ, or by PIP-ECHO.
- An unmappable source value gives null and MUST be kept as cv_param `<tool>_<column>` (e.g. `maxquant_type`); mapped values SHOULD be kept.
- Without such a statement, write `transferred` only when every identification of the run the tool could have attached to the signal is visible, assigned or not (OpenMS consensusXML without the meta value: the consensus feature's PIDs plus those IDConflictResolver removed), and none is of the row's peptidoform and charge; if some name another peptidoform, also set `identification_conflict = true` ([T06](t06-identifications-per-feature.md) FID-1, FID-10). Otherwise null.
- Never from imputation, a missing PSM, a null run or an empty scan.
- A producer whose rows all rest on own-run identifications writes `direct`/`secondary`, not null.
- So null means only unknown; no footer key qualifies it (MBR-8).

**MBR-3 (MUST)** A transferred row:
- has no identifying or supporting PSM, so `scan = []`, null `id_run_file_name`, no `psm.feature_id` link ([T06](t06-identifications-per-feature.md) FID-1, FID-5 to FID-7) and `spectral_count`, where written, 0 (FID-12); a `MULTI-MATCH-MSMS` scan MAY go to cv_params;
- has the confidence columns [T04](t04-confidence-fields.md) CNF-8 sets: null PEP, `psm_qvalue`, `precursor_qvalue` and run-scope `peptide_qvalue`, no `psm`-level `additional_scores`; experiment-scope q-values MAY be filled;
- has null `identified_labels`;
- carries the `rt`, `observed_mz` and intensities of the matched signal here ([T05](t05-feature-rows.md) FRW-8, [T11](t11-spectrum-precursor-fields.md) SPC-1);
- has at least one non-null entry whose label is not in `imputed_labels`; other labels per [T09](t09-intensities-missing-values.md) INT-5. A transfer without a measured value (isobaric MBR rows without reporters, andes#148 zeros, [T09](t09-intensities-missing-values.md) INT-7) is not written and is counted in `skipped_rows_feature` ([T03](t03-filtering-state.md) FLT-6).

**MBR-4 (MUST)** `match_*` are null on rows that are not `transferred`.
- `match_qvalue` (float64 in [0, 1], [T04](t04-confidence-fields.md) CNF-9): q-value of transferring this precursor ion into this run (MaxQuant `Match q-value`; OpenMS `mbr_transfer_qvalues`, the per-transfer q-values under PIP-ECHO), never the donor's.
- `match_score`: the tool's score of the transfer itself, higher better. Donor-identification scores go to `additional_scores`: MaxQuant `Match score` → `maxquant_match_score`.
- `match_rt_delta` (s) = acceptor `rt_aligned` − donor-derived aligned RT; `match_mz_delta` (Th) = acceptor measured m/z − donor-expected m/z. Unknown sign → null. MaxQuant complies.
- A transfer PEP or probability goes to `<tool>_match_pep`/`<tool>_match_probability`, level `transfer` ([T04](t04-confidence-fields.md) CNF-5).

**MBR-5 (MUST)** `id_run_file_name` never names a donor ([T06](t06-identifications-per-feature.md) FID-2, FID-6). With `consensus_feature_id` ([T06](t06-identifications-per-feature.md) FID-11), candidate donors are the `direct`/`secondary` rows sharing it.

**MBR-6 (MUST)** `imputed_labels`/`is_imputed` say whether a value was measured, `identification_type` only where the identity came from; neither sets the other. A transferred quantity comes from the quantification tool's own algorithm, so it counts as measured unless the tool itself marks it as imputed ([T09](t09-intensities-missing-values.md) INT-9).

**MBR-7 (MUST)** Rows of different types sharing `(peptidoform, charge, run_file_name, compensation_voltage)` are distinct signals and stay separate; two rows known to be one peak become one row with the own-run type ([T05](t05-feature-rows.md) FRW-2). De-duplication ranks `direct` > `secondary` > `transferred` by `identification_type`, never by the raw tool string. The OpenMS strict profile goes further and requires that OpenMS does not write a transferred row beside a quantified direct row of the same precursor and run ([strict profile](strict-profile.md)).

**MBR-8 (MUST)** No footer key declares the transfer state; rows and provenance carry it:
- `identification_type` is filled whenever the tool states the origin or every row rests on own-run identifications (MBR-2); pg `identified_in_unit` whenever the source states it (MBR-10). Null = unknown.
- Whether the tool made transfers goes in the provenance `match_between_runs` step: parameter `match_between_runs` = `true`/`false` (as `docs/spec/provenance.md:187`), with the settings of MBR-9.
- A qpxc option that drops transfers records provenance option `transferred_rows_removed = true` ([T03](t03-filtering-state.md) FLT-5).

Unless provenance records `match_between_runs = false`, identified pg rows may include transferred signal. Merges keep every input's provenance steps ([T13](t13-metadata-conformance.md) META-4).

**MBR-9 (MUST)** The `match_between_runs` step parameter `mbr_fdr` follows [T03](t03-filtering-state.md) FLT-2 at the transfer level (a precursor ion in a run): decimal x = every transfer had a tool-estimated FDR ≤ x (IonQuant `--ionfdr`), so non-null `match_qvalue` ≤ x; `applied` = a non-FDR cut (RT and m/z windows are not cuts); `none` = no error control; `unknown`. Peptide- and protein-level transfer FDRs and tolerances are parameters of the same step (`mbr_peptide_fdr`, `mbr_protein_fdr`, `mbr_rt_tolerance`; precedent `docs/spec/provenance.md:187-188`). Per-row `match_qvalue` stays (MBR-4).

**MBR-10 (MUST)** New pg column `identified_in_unit` (bool, nullable), filled whenever the source states it: `true` = evidence in the unit as [T12](t12-protein-groups.md) PGV-2 defines it; `false` = the row exists only for a non-null quantity; null = unknown. [T12](t12-protein-groups.md) PGV-1 adopts it.

**MBR-11 (MUST)** Decoy transfers made for error control (FlashLFQ decoy-peptide or random-RT peaks) are never written as target rows; decoy-peptide transfers MAY appear with `is_decoy = true` ([T03](t03-filtering-state.md) FLT-10).

**MBR-12 (MUST)** Off-vocabulary types (any value other than `direct`, `secondary`, `transferred`, `unidentified` or null) are errors on every write path and in the validator. For `qpx_version` ≥ 1.4 the following are errors under the `dda` profile and warnings otherwise ([T13](t13-metadata-conformance.md) META-8): MBR-3 violations (non-empty `scan`, non-null `id_run_file_name`, PEP or run-scope q-values, a linked PSM, no measured entry); `match_*` on non-transferred rows; `match_qvalue` outside [0, 1] or above a decimal provenance `mbr_fdr`; transferred rows when provenance records `match_between_runs = false` or `transferred_rows_removed = true`; `identified_in_unit = false` with null intensity.

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| `identification_type` enum, MBR-1 text (`feature.md:145`, `feature.yaml:210-214`), checked on both write paths | feature | `direct`, `transferred`, `secondary`, `unidentified` ([T07](t07-unidentified-features.md) UNI-2), null | Non-null when the tool states the origin or every row is own-run (MBR-2); `unidentified` on every unidentified row | Yes for writers with other values (none in qpx); `MULTI-MATCH-MSMS` rows lose scan, run, PEP |
| `match_*` definitions (`feature.yaml:215-234`) | feature | `match_qvalue` float64 [0, 1] | Optional | **Yes**: physical type ([T04](t04-confidence-fields.md) CNF-9) |
| `match_score` = transfer score; MaxQuant value → `maxquant_match_score` | feature | – | – | **Yes**, semantic: 1.3 MaxQuant files hold the donor score (MBR-R7) |
| `match_between_runs` step: `match_between_runs`, `mbr_fdr`, `mbr_peptide_fdr`, `mbr_protein_fdr`, `mbr_rt_tolerance`; option `transferred_rows_removed` (MBR-8, MBR-9) | provenance | `match_between_runs`: `true`, `false`; `mbr_fdr`: decimal, `applied`, `none`, `unknown` | Optional | No |
| `identified_in_unit` | pg | bool, nullable | Non-null whenever the source states it (MBR-10) | No |
| Edit `feature.md:138`, `producer-guide.md:26`, `pg.md:131-135,291`, `anndata.md:153`, `converter-coverage.md:55,199-203`, `provenance.md:187-188` | docs | – | – | No |

## Producer changes
| Producer | Change needed |
|---|---|
| MaxQuant | Map Types by token: `*MATCH*` transferred, `*SECPEP` secondary, other `*MSMS*` direct, else null with warning; keep `maxquant_type`. Null run, scan, PEP by type (feature_adapter.py:655-657); dedup by type (:264). Move `Match score`. Drop isobaric MBR rows without reporters, counted. pg: "By MS/MS" → `true`, "By matching" → `false`. Provenance `match_between_runs` from parameters.txt; `mbr_fdr` from parameters.txt if stated, else `unknown` |
| FragPipe | Provenance `mbr_fdr` = `ionquant.ionfdr` from a supplied fragpipe.workflow when MBR is on, else `unknown`. pg per [T12](t12-protein-groups.md) PGV-2 |
| OpenMS (qpxc; `-out_qpx` upstream) | A run listed in the consensus feature's `mbr_transfer_map_indices` → transferred; `match_qvalue` from `mbr_transfer_qvalues` (PIP-ECHO). Without the meta value, MBR-2's fallback using the removed-PID index (feature_adapter.py:573-613): the IDConflictResolver case becomes `transferred` + `identification_conflict`. pg: `identified_in_unit` = a group PSM in the unit (pg_adapter.py:453-465). Native writer (verified): HEAD writes no `identification_type` or `match_*`, and its pg evidence counts every run the feature was quantified in, transfers included; the unmerged [#378](https://github.com/bigbio/qpx/issues/378) branch adds `identification_type` (direct, transferred, null) and `match_qvalue` (float32, from PIP-ECHO). See the [OpenMS worklist](openms-worklist.md). |
| Proteome Discoverer | Features keep `direct`. pg: `identified_in_unit` from `Found in Sample` (`Peak Found` → `false`), else null |
| CDAP; mzIdentML (pg only) | Provenance `match_between_runs = false`; pg `identified_in_unit = true` |
| quantms-MSstats | Null run ([T06](t06-identifications-per-feature.md)). LFQ: null type. Isobaric: `direct`, provenance `match_between_runs = false` |
| Spectronaut, DIA-NN | Stop deriving the type from `EG.IsImputed` (spectronaut/feature_adapter.py:342); keep `imputed_labels`; null type until DIA rules exist |

## Reader rules
**MBR-R1 (MUST)** "Identified in this run" is `identification_type IN ('direct','secondary')` ([T02](t02-chimeric-spectra.md) CHI-R4); "transferred" is `'transferred'`. Never infer either from run, scan, PEP or null q-values ([T06](t06-identifications-per-feature.md) FID-R3).

**MBR-R2 (MUST)** Treat a null type as own-run only when provenance records `match_between_runs = false` or `transferred_rows_removed = true`; otherwise report the null-row count and the assumption made. Reject off-vocabulary values.

**MBR-R3 (MUST)** Identification counts, protein evidence and inference use only `direct`/`secondary` rows and pg rows with `identified_in_unit = true` (andes#218).

**MBR-R4 (MUST)** Before combining rows into one (peptidoform, charge, run or unit, label) value, drop the transferred entry of a label whenever a `direct` or `secondary` row of the same key has a measured entry (non-null, not imputed) for that label; never sum, average or MAX across them. A transferred entry beside an own-run row whose entry for that label is null is kept. [T14](t14-reader-contract.md) RDR-5 applies this rule in the default read profile.

**MBR-R5 (MUST)** Choose and record the transfer policy (keep, drop, `match_qvalue ≤ t`) separately from confidence filters, which never drop transfers merely for null values ([T04](t04-confidence-fields.md) CNF-R5). Whether to use transferred quantities is the reader's choice. The default keeps them, because they improve quantification ([T14](t14-reader-contract.md) RDR-4); filtering them out is optional.

**MBR-R6 (SHOULD)** When provenance `mbr_fdr` is `none`, `unknown` or absent, warn and offer to drop transfers. Transfer weights (mokume#117) use only `identification_type`, `match_qvalue` and provenance `mbr_fdr`, never `match_score` or donor scores.

**MBR-R7 (MUST)** In files with `qpx_version` < 1.4: `match_score` from `creator = maxquant` is the donor's Andromeda score; a `transferred` row is transferred even with scan, run or PEP; Spectronaut `transferred` means imputed (use `imputed_labels`).

Current violations:
- qpx reads the type only in validation (`qpx/core/data/validation.py:71`). Run and QC counts include transfers (`qpx/views/api.py:21-29,172-183`), against R3.
- `ds.intensity` sums every row (`qpx/core/dataset_queries.py:124-137`); `build_mudata` defaults to `duplicates="sum"` (`qpx/mudata.py:685-691`), against R4. `Feature.confident()` drops transfers (`qpx/core/data/feature.py:52-54`), against R5.
- mokume Python never reads the type. mokume Rust falls back to `id_run_file_name` and counts off-vocabulary values as identified (mokume:rust/crates/mokume-io/src/qpx.rs:501-523), against R1, R2.
- mokume takes the MAX across direct and transferred rows (mokume:python/mokume/preprocessing/aggregation.py:316-320; mokume:rust/crates/mokume-pipeline/src/lib.rs:2761-2772), against R4.
- mokume score filters drop transfers (mokume:python/mokume/io/feature.py:266-270; mokume:rust/crates/mokume-pipeline/src/lib.rs:6966-6972), against R5; the Rust default `keep` (mokume:rust/crates/mokume-core/src/config.rs:172-180) conforms once recorded.

## Open questions for maintainers
1. Does a `MULTI-MATCH-MSMS` MS/MS ever identify its own peptidoform at PSM FDR? Then MBR-1 makes it `direct` (also [T06](t06-identifications-per-feature.md) question 2).
2. No producer fills a transfer-quality `match_score`: deprecate it at the next breaking release (P-I), or add typed donor fields ([EXTERNAL-11](findings-feature.md#external-11))?
3. Keep MaxQuant isobaric MBR rows with only `ms1_intensity` ([T09](t09-intensities-missing-values.md) INT-11), relaxing MBR-3?
4. `build_mudata` default: in-run precedence (R4) or an error on duplicates?
5. Which MaxQuant versions control transfer FDR (coxdocs defines `Match q-value`; the fixture's is empty), and which parameter says so? DIA library transfer needs its own definition.
