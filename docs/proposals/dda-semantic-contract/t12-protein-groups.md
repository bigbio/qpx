# T12. Protein-group rows: existence, evidence, counts, q-values, decoys and quantities

[← T11](t11-spectrum-precursor-fields.md) · [Overview](index.md) · [Strict profile](strict-profile.md) · [Annex overview](annex-overview.md) · [T13 →](t13-metadata-conformance.md)

!!! note "Annex: producer-neutral analysis"
    Part of the annex of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text, and they cover every producer. For OpenMS-written files the [strict profile](strict-profile.md) is the contract and wins where the two differ. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** A pg row is one group in one tool quantification unit for one label; it exists when the group has rank-1, non-transferred evidence there or a tool quantity, `identified_in_unit` says which, every count and q-value column has one fixed scope, and provenance records retention (`unquantified_rows`) and grouping (`pg_grouping`).

## Questions a reader has
- Does no row mean "not identified", a null intensity "identified, not quantified"? Does every channel of an identified plex get a row?
- Does a quantified row rest on MS/MS in this unit, or only on match-between-runs (andes#218)?
- Which threshold, ranks, ties, chimeric co-identifications, decoys and transfers make evidence?
- Is `grouped_runs` a raw file, a sample's fractions or the dataset? Is its order meaningful? Are groups the same in every unit ([#51](https://github.com/bigbio/qpx/issues/51))?
- Do counts cover the unit or the analysis, sequences or PSMs, identified or transferred peptides?
- Is `pg_qvalue` unit-level and `global_qvalue` experiment-level everywhere? Where is the PEP? How are `is_decoy` and `contaminant` set for mixed groups; does `false` mean "checked"?
- Which quantity is `intensity`, and can it differ by row?

## What happens today
Spec: `grouped_runs` is one unit's fractions, stored in fraction order (MUST) (docs/spec/pg.md:3,40,414; data-model.md:29-30); views.md:214 wrongly says CDAP writes per raw file (qpx/converters/cdap/pg_adapter.py:161). Existence follows [#280](https://github.com/bigbio/qpx/issues/280)/#363 (pg.md:116-146); evidence is "a PSM that passed the FDR threshold" (:126-128), threshold unrecorded [[SPEC-15](findings-pg.md#spec-15)]. MBR-only rows are allowed (:131-135), marked only by a MaxQuant cv_param (:291). Counts are "in this file" (:81-86), yet intensities.md:206 allows experiment-wide counts. global = experiment, pg = run (:92-93; pg.yaml:52-61); ontology.md:238 maps `global_qvalue` to the list-level MS:1001214. FragPipe guidance contradicts itself (pg.md:405-406 vs converter-coverage.md:299-300). Mixed-group decoys and null `contaminant` are undefined (pg.md:94-95). `pg_names` = entry names (:35), descriptions in the example (:181). [#51](https://github.com/bigbio/qpx/issues/51) is not in the spec; intensities.md:207 only implies it.

| Producer | Behaviour | Findings |
|---|---|---|
| MaxQuant | `pg_qvalue` copies the experiment `Q-value` (qpx/converters/maxquant/pg_adapter.py:484; 27,072/27,072 fixture rows); `PEP` unread. 5,687 `By matching` rows carry intensity, flagged by cv_param only (:403-406), and all have `Peptides <exp>` ≥ 1: per-experiment counts include transfers. Null rows need MS/MS count > 0 (:509,564,589); 489 `By MS/MS` rows have 0. `peptides` falls back to the group total per member (:495-497) | [[MAXQUANT-18](findings-pg.md#maxquant-18)] [[MAXQUANT-19](findings-pg.md#maxquant-19)] [[MAXQUANT-20](findings-pg.md#maxquant-20)] |
| FragPipe | Primary = first positive of `Intensity`, `Total Intensity`, `MaxLFQ Intensity` per cell (qpx/converters/fragpipe/pg_adapter.py:71). Runs sorted `[runB_F1, runB_F10, runB_F2]` (:96; isobaric.py:39). `pg_qvalue` = 1 − probability (:418). One Description as `pg_names` (:437). Leader decides `is_decoy` (:442). TMT drops null channels (:382) | [[FRAGPIPE-20](findings-pg.md#fragpipe-20)] [[FRAGPIPE-21](findings-pg.md#fragpipe-21)] [[FRAGPIPE-22](findings-pg.md#fragpipe-22)] [[FRAGPIPE-23](findings-pg.md#fragpipe-23)] [[FRAGPIPE-M1](findings-pg.md#fragpipe-m1)] [[FRAGPIPE-M2](findings-pg.md#fragpipe-m2)] |
| Proteome Discoverer | Cells ≤ 0 skipped (qpx/converters/proteome_discoverer/pg_adapter.py:113). Identification-only: one `LFQ` row over all runs, sorted (:101). `Exp. q-value: Combined` → `pg_qvalue` (:148). Experiment-wide counts; `peptides` = master only (:151); `contaminant` false without the column (:150) | [[PD-13](findings-pg.md#pd-13)] [[PD-14](findings-pg.md#pd-14)] [[PD-15](findings-pg.md#pd-15)] [[PD-16](findings-pg.md#pd-16)] [[PD-23](findings-cross.md#pd-23)] |
| OpenMS consensusXML | Evidence from hits of any rank, decoys, and the consensus feature's runs (qpx/converters/openms_consensus/pg_adapter.py:111,120). Global counts per unit row (:366-387). `global_qvalue` = min member-protein q over all ProteinIdentifications (:260-297). Native `-out_qpx` files pass through; verified at HEAD: OpenMS's own abundance per (fraction group, label); identification-only rows with null label and intensity; evidence counts every run the feature was quantified in, transfers included; `is_decoy` from the anchor only; `contaminant` and `pg_qvalue` always null; the group probability is reported as a q-value only when a protein FDR ran. See the [OpenMS worklist](openms-worklist.md). | [[OPENMS-23](findings-pg.md#openms-23)] [[OPENMS-25](findings-pg.md#openms-25)] [[OPENMS-26](findings-pg.md#openms-26)] [[OPENMS-31](findings-cross.md#openms-31)] |
| mzIdentML | One row per PAG over every run (qpx/converters/mzidentml/pg_adapter.py:139,240); no run view; passThreshold ignored; `peptide_count` = PeptideHypothesis count (:188), up to 3,496; `is_decoy` = any prefix (:231); PDH scores lack member index: 8/123 PXD054720 rows repeat `mascot_score`; ids and counts stored as scores | [[MZID_CDAP-07](findings-pg.md#mzid_cdap-07)] [[MZID_CDAP-27](findings-pg.md#mzid_cdap-27)] [[MZID_CDAP-28](findings-cross.md#mzid_cdap-28)] [[SPEC-30](findings-cross.md#spec-30)] |
| CDAP | Mapping-set groups; `total_sequences` counts features (qpx/converters/cdap/pg_adapter.py:171); `peptides` = anchor only (:145); `pg_names` = accessions (:156); `contaminant` False on 810/810 rows (:167); PSM/peptide q-values in scores (:150); tied PSMs give no evidence | [[MZID_CDAP-16](findings-feature.md#mzid_cdap-16)] [[MZID_CDAP-25](findings-pg.md#mzid_cdap-25)] [[MZID_CDAP-26](findings-pg.md#mzid_cdap-26)] [[MZID_CDAP-M3](findings-cross.md#mzid_cdap-m3)] [[MZID_CDAP-M5](findings-pg.md#mzid_cdap-m5)] |

DIA-NN (reference) has unit and global group q-values but picks `PG.Quantity` or MaxLFQ per row (qpx/converters/diann/pg_adapter.py:215). Mixed decoy groups: any member (Spectronaut, mzIdentML), all (DIA-NN :238, OpenMS, CDAP), leader (FragPipe) [[SPEC-M3](findings-pg.md#spec-m3)]. `contaminant`: constant false (CDAP), prefix-parsed (FragPipe, DIA-NN, OpenMS) [[SPEC-M4](findings-pg.md#spec-m4)]. No producer writes a filter or grouping key [[INFRA-20](findings-cross.md#infra-20)]. quantms-MSstats writes no pg view (converter-coverage.md:293).

## Why it matters
- `proteins_per_sample` counts 5,687 transfer-only MaxQuant rows as identified [[CRITIC-49](findings-pg.md#critic-49)].
- PD and FragPipe TMT make "identified, not quantified" look like "not identified" [[PD-13](findings-pg.md#pd-13)] [[FRAGPIPE-M1](findings-pg.md#fragpipe-m1)].
- Unit-level protein FDR filters read experiment values (MaxQuant, PD) or probabilities (FragPipe) [[SPEC-24](findings-pg.md#spec-24)].
- "≥ 2 peptides in this sample" reads PSM counts (mzIdentML), feature counts (CDAP), global counts (OpenMS, PD) or transfer-inclusive counts (MaxQuant) [[CRITIC-54](findings-pg.md#critic-54)].
- CDAP contaminants count as sample proteins; `targets_only()` keeps different mixed groups per producer [[CRITIC-53](findings-pg.md#critic-53)].

## Options considered
| Option | For | Against |
|---|---|---|
| MBR-only flag: cv_param / pg `identification_type` / **`identified_in_unit` ([T08](t08-match-between-runs.md) MBR-10)** | – / mirrors feature / exactly the [#280](https://github.com/bigbio/qpx/issues/280) split | untyped / "direct" overclaims mixed quantities / new column |
| Identified-unquantified rows: MUST ([#280](https://github.com/bigbio/qpx/issues/280)) / **declared; `retained` for qpx converters** | decided / P-B, filtered files stay readable | filtered files non-conforming / relaxes [#280](https://github.com/bigbio/qpx/issues/280) elsewhere |
| Count scope: footer key / extra columns / **fixed per column** | flexible / complete / P-A, testable | breaks P-A / bloat / null where a tool lacks the scope |
| Global groups: MUST ([#51](https://github.com/bigbio/qpx/issues/51)) / **SHOULD + `pg_grouping` in provenance** | simple / feasible | merges and consensusXML of OpenMS workflows with per-run inference cannot comply / readers consult provenance |

## Recommendation
**PGV-1 (MUST)** Row existence per group, unit and label:

| Evidence in unit (PGV-2) | Tool quantity | Row | `identified_in_unit` | `intensity` |
|---|---|---|---|---|
| yes | yes | MUST | true | value |
| yes | no | per PGV-4 | true | null |
| no | transferred, mapped (PD `Peak Found`), imputed ([T09](t09-intensities-missing-values.md) INT-9) | MUST, except transfers removed under provenance `transferred_rows_removed = true` ([T08](t08-match-between-runs.md) MBR-8) | false | value |
| unknown | yes | MUST | null | value |
| no or unknown | no | MUST NOT | – | – |

A unit with evidence has one row per label with a sample in its `run.samples` (without a sample map, per reported channel, [T09](t09-intensities-missing-values.md) INT-5): every channel of an identified plex. `label` is null only without label context (no run view: mzIdentML), then with null `intensity` [[VALIDATOR-28](findings-pg.md#validator-28)]. `identified_in_unit` is [T08](t08-match-between-runs.md) MBR-10.

**PGV-2 (MUST)** A group has evidence in a unit when an *evidence identification* from a raw file of the unit maps to a member, by a unique or shared peptide, assigned to a feature or not (pg.md:126-128); under `mapping_set`, only to the group equal to its mapped entries. An evidence identification is rank 1 ([T01](t01-psm-rows.md) PSM-R1; every rank-1 tie and chimeric co-identification, PSM-6, [T02](t02-chimeric-spectra.md) CHI-4), not transferred (andes#218), not marked failing (`pass_threshold` ≠ false, [T03](t03-filtering-state.md) FLT-7), and target for target groups, decoy for decoy groups. This is the producer's rule. It counts a null `pass_threshold` as evidence, so provenance MUST state the cut that stands behind such rows (PGV-3). A reader that re-derives evidence from psm rows uses [T03](t03-filtering-state.md) FLT-R4, which agrees with this rule whenever PGV-3's `<level>_fdr` parameters are stated. A producer MAY use the tool's per-unit statement instead: MaxQuant `Identification type <exp>` (`By MS/MS` true, `By matching` false; without the column and without MBR, `Peptides <exp>` > 0), FragPipe `<exp> Total Spectral Count` > 0, PD `Found in Sample` (`Peak Found` false; to confirm). An unidentified feature row ([T07](t07-unidentified-features.md) UNI-3) is never evidence, and it enters no pg row, count or quantity.

**PGV-3 (MUST)** For pg rows, the provenance parameters `psm_fdr`, `precursor_fdr` and `peptide_fdr` ([T03](t03-filtering-state.md) FLT-1) state the cut on PGV-2 evidence and PGV-7 counts, and `protein_fdr` the cut on groups present (FLT-2). They are pg.md:127's "FDR threshold". qpxc applies no protein cut ([#342](https://github.com/bigbio/qpx/issues/342)).

**PGV-4 (MUST)** The provenance parameter `unquantified_rows` (absent = unknown) of [T07](t07-unidentified-features.md) UNI-6 has this pg meaning: `retained` (every group, unit and label with evidence has a row), `removed` (every row has a non-null intensity). One value covers the feature and pg views the step wrote (UNI-6). qpx converters record `retained` ([#280](https://github.com/bigbio/qpx/issues/280)/#363 kept), or nothing in the PGV-5 unattributable case; `removed` serves files filtered later or by other producers (P-B), a declared relaxation of [#280](https://github.com/bigbio/qpx/issues/280).

**PGV-5 (MUST)** `grouped_runs` holds exactly the raw files the tool combined into the row's quantity (MaxQuant or FragPipe experiment, PD file group, OpenMS fraction group or SDRF sample, CDAP sample); a row without a tool quantity covers one unit per qpx/core/sample_sql.py:30-37, or one raw file without a run view. Members SHOULD share one `run.samples` set; a tool unit pooling samples (MaxQuant experiment of replicates) is written as reported and resolves to no single sample. A proper subset of a unit is a *partial unit*. A producer unable to attribute evidence to units writes null `identified_in_unit` on quantified rows and no unquantified row, counting the groups left out ([T03](t03-filtering-state.md) FLT-6). Order SHOULD follow `run.fraction`, else source order, never a lexicographic re-sort; this relaxes pg.md:414, as FragPipe manifests carry no fraction.

**PGV-6 (MUST)** `pg_grouping`, a parameter of the provenance `protein_inference` step, applies to pg, and to feature when `pg_accessions` is written (absent = unknown): `dataset_inference` (one inference over the file, so groups agree across units, [#51](https://github.com/bigbio/qpx/issues/51): MaxQuant, FragPipe, PD, mzIdentML, and always native OpenMS output, because ProteomicsLFQ merges to one identification run before inference and the native exporter requires that single run), `partial_inference` (merged files; consensusXML of other OpenMS workflows with per-run ProteinIdentifications, read by the Python `openms-consensus` converter), `mapping_set` (a group is the entries one peptide maps to: CDAP, DIA-NN `--no-prot-inf`). Producers SHOULD reach `dataset_inference`; [#51](https://github.com/bigbio/qpx/issues/51) becomes a recorded property. A merge keeps every input's `protein_inference` step ([T13](t13-metadata-conformance.md) META-4); a file holding steps of different analyses counts as `partial_inference` unless all are `mapping_set`.

**PGV-7 (MUST)** Fixed count scopes, distinct counts over PGV-2 evidence:

| Field | Scope | Counts |
|---|---|---|
| `peptide_counts.total_` / `unique_sequences` | row's unit, all labels | sequences mapping to a member / only to this group ([T10](t10-protein-mapping.md) PRT-9) |
| `feature_counts.total_` / `unique_features` | row's unit, all labels | (peptidoform, charge), same split |
| `peptides[]` | analysis | one entry per member, set-equal to `pg_accessions`: sequences mapping to it; null when unknown |
| `sequence_coverage` | analysis | anchor protein |

Under `mapping_set` total equals unique. A tool count of another scope or basis (MaxQuant `Peptides <exp>`) is null; it MAY go to cv_params under an explicit name.

**PGV-8 (MUST)** Q-values ([T04](t04-confidence-fields.md) CNF-1/2/5): `global_qvalue` is the group-level q-value (MS:1002373) over the analysis, equal on all rows of a group from one analysis. `pg_qvalue` is that statistic within the row's unit (DIA-NN `PG.Q.Value`), null when the tool has only an analysis value, never a copy. Each is filled exactly when the source estimated protein-group FDR at that scope ([T04](t04-confidence-fields.md) CNF-2). `global_qvalue` is therefore nullable for every producer, OpenMS included: when no protein-level FDR was estimated (a tool may stop at PSM- or peptide-level FDR, as OpenMS sometimes does), it is null on every row. This is not a defect and not `unknown`. 1 − probability, member-protein minima and lower-level values go to `pg_probability`, `protein_qvalue`/`protein_global_qvalue` by scope, and `pg_pep`. Other pg scores are group statistics, or summaries under a non-reserved name whose ontology.parquet property `level` states their level (`top_peptide_probability`: `peptide`), never reserved PSM-, precursor- or peptide-level names ([T04](t04-confidence-fields.md) CNF-5). A score name appears once per row (PDH scores: anchor only); counts and identifiers are not scores. Delete the per-protein sentence of scores.md:87.

**PGV-9 (SHOULD)** Optional pg `pass_threshold` (bool, null), [T03](t03-filtering-state.md) FLT-7 at group level, where the source states acceptance (mzIdentML PAG MS:1002415; PD `Protein FDR Confidence` High, to confirm), with the ontology properties `level = protein_group` and `fdr` on the pg `pass_threshold` row (FLT-8).

**PGV-10 (MUST)** `is_decoy` per [T10](t10-protein-mapping.md) PRT-5 (all members decoy); `contaminant` per PRT-6 (tool flag, else any member matching a recorded marker; null = unknown; never constant). Markers stay on accessions (PRT-3), so no aggregation merges decoy and target groups [[MZID_CDAP-06](findings-pg.md#mzid_cdap-06)].

**PGV-11 (MUST)** `intensity` follows [T09](t09-intensities-missing-values.md): one source column per file (INT-1); kind in the ontology properties `quantity` and `processing` of the pg `intensity` row (INT-2), superseding cv_param `quantification_method` (INT-3); a QPX roll-up only when the tool provides no protein quantity in the whole file, declared `processing = sum` and never mixed with tool quantities (INT-4); `lfq` → `maxlfq` (INT-12). Such a roll-up stands for the tool quantity in PGV-1, and its rows cover the unit of a row without a tool quantity (PGV-5). An identified row may include transferred signal ([T08](t08-match-between-runs.md) MBR-8).

**PGV-12 (MUST)** `pg_names` is null or one entry name per member, aligned with `pg_accessions`; descriptions MAY go to cv_param `protein_description`. Fix pg.md:181.

**PGV-13 (MUST)** Validator, `qpx_version` ≥ 1.4. Errors ([T13](t13-metadata-conformance.md) profile): null intensity with `identified_in_unit` ≠ true (covers [T08](t08-match-between-runs.md) MBR-12); under `unquantified_rows = retained`, a group and unit with evidence missing a sample label; `peptides` not set-equal to members; `unique_* > total_*`; `total_features < total_sequences`; a repeated or reserved lower-level score name; `pg_names` length ≠ members. Warnings: `global_qvalue` varying within a group of one source; `pg_qvalue = global_qvalue` on every row of a multi-unit file; order against numeric `run.fraction`; members with different `run.samples`. Info: partial units ([T13](t13-metadata-conformance.md) META-8).

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| `identified_in_unit` ([T08](t08-match-between-runs.md)) | pg | bool, null | non-null whenever the source states it ([T08](t08-match-between-runs.md) MBR-10) | No |
| `pass_threshold` | pg | bool, null | optional | No |
| `unquantified_rows` (shared with feature, [T07](t07-unidentified-features.md) UNI-6); no footer key | provenance (the step that wrote the QPX files) | `retained`, `removed` | present in [T13](t13-metadata-conformance.md) profile; `retained` for qpx converters | No |
| `pg_grouping` | provenance `protein_inference` step (pg, feature) | `dataset_inference`, `partial_inference`, `mapping_set` | SHOULD | No |
| `pg_qvalue` = unit level (pg.md:93, pg.yaml:61, ontology.md:239); `global_qvalue` → MS:1002373 (ontology.md:238) | pg | float64 | – | **Yes**, semantic: MaxQuant, PD, FragPipe values move or become null |
| `peptides` nullable, one entry per member | pg | list, null | – | **Yes**: readers must accept null |
| Count scopes; all-member `is_decoy`; null `contaminant`; one row per tool unit; aligned `pg_names`; one score name per row | pg | – | – | **Yes**, semantic for existing MaxQuant, FragPipe, PD, OpenMS, mzIdentML, CDAP, Spectronaut files |
| `grouped_runs` order MUST → SHOULD | pg | – | – | Relaxation; position no longer means fraction |

## Producer changes
| Producer | Change needed |
|---|---|
| MaxQuant | `identified_in_unit` and null rows from `Identification type` (pg_adapter.py:509,564,589); `pg_qvalue` null, `PEP` → `pg_pep`; unit counts from evidence.txt non-`MATCH` rows, else null; coverage `Sequence coverage [%]`; `peptides` null without `Peptide counts (all)`; names from `Fasta headers`; provenance `unquantified_rows = retained` and `pg_grouping = dataset_inference` |
| FragPipe | Source or `run.fraction` order (pg_adapter.py:96; isobaric.py:39); `pg_qvalue` null, probability → `pg_probability`; null rows for annotated channels (:382); `identified_in_unit` from `Total Spectral Count`; names from `Indistinguishable Proteins` headers, else null; all-member `is_decoy` |
| Proteome Discoverer | Evidence per file group from `Found in Sample` or PSMs.txt `File ID`; null rows for unquantified channels (pg_adapter.py:113); per-file-group rows with run-view labels, not one all-run `LFQ` row (:101); q-value → `global_qvalue`; unit counts from PSMs.txt, else null; `peptides` per member from its Proteins.txt row; null `contaminant` without the column; PGV-9 |
| OpenMS consensusXML | Evidence from rank-1 hits in the PSM's own run (pg_adapter.py:111,120); counts per unit (:366-387); member minimum → protein-level name; provenance `pg_grouping` from the ProteinIdentification count: one (always the case for ProteomicsLFQ output) → `dataset_inference`; several (other workflows) → `partial_inference` |
| OpenMS `-out_qpx` | OpenMS writes `identified_in_unit` and the provenance parameters `unquantified_rows` and `pg_grouping`, always `dataset_inference`; `qpxc convert openms` keeps them. Verified at HEAD: evidence counts every run the feature was quantified in, transfers included (PGV-2 needs own-run identifications only); `is_decoy` comes from the anchor only (PGV-10: all members); `contaminant` and `pg_qvalue` are always null; identification-only rows have a null label (PGV-1: one row per label); no provenance view. See the [OpenMS worklist](openms-worklist.md). |
| mzIdentML | One row per PAG and raw file with evidence (SII → SpectraData); `pass_threshold` from MS:1002415; `peptides` = distinct sequences of each PDH's evidence SIIs (MS:1001097 is unthresholded); member decoys from PeptideEvidence `isDecoy` ([T10](t10-protein-mapping.md)); anchor PDH scores only; ids and counts to cv_params; provenance `pg_grouping = dataset_inference` |
| CDAP | Distinct-sequence counts; `peptides` per member; `pg_names` null; `contaminant` null unless a list is given; PSM/peptide q-values dropped or renamed (PGV-8); AmbiguousMatch PSMs are evidence; provenance `pg_grouping = mapping_set` |
| quantms-MSstats | None: no pg view |
| Writers, validator | Record `unquantified_rows`; merge keeps every input's provenance steps ([T13](t13-metadata-conformance.md) META-4; PGV-4, PGV-6); PGV-13 |

## Reader rules
**PGV-R1 (MUST)** Count a group as identified in a unit only on rows with `identified_in_unit = true`, `is_decoy = false` and protein acceptance (PGV-R3), never from row existence; null `identified_in_unit` gives upper bounds. Violation: `proteins_per_sample` counts every row (qpx/views/api.py:247-268).

**PGV-R2 (MUST)** A missing row means "not identified" only when provenance states `unquantified_rows = retained`. Exports dropping null-intensity rows carry identification state separately. Violation: MuData proteins keep non-null intensities only (qpx/mudata.py:233-240).

**PGV-R3 (MUST)** Filter protein FDR on `global_qvalue`, `pg_qvalue` per unit, or `pass_threshold` with its basis in the ontology properties of its row ([T03](t03-filtering-state.md) FLT-8), never on scores. A null q-value is unassessable unless provenance states `protein_fdr` as a decimal ≤ t ([T03](t03-filtering-state.md) FLT-R3, [T14](t14-reader-contract.md) RDR-3). When the q-value column is null on every row, protein-level FDR was not estimated (PGV-8) and filtering at that level is not available: the reader says so and filters at a level that exists. It never reads a null as 1.0 and never drops every row silently. Violations: `ProteinView.intensity` passes nulls (qpx/views/api.py:48); mokume never reads pg.parquet, using `feature.pg_global_qvalue` with a `protein_qvalue` fallback (mokume:python/mokume/io/feature.py:601-603; mokume:rust/crates/mokume-io/src/qpx.rs:30; mokume:python/mokume/preprocessing/filters/protein.py:193-228).

**PGV-R4 (MUST)** A row quantifies its `grouped_runs` set and label. Fractions come from `run.fraction`, never position. Sum partial units only when the pg `intensity` row has the ontology property `processing = sum`. A row resolving to several samples for its label is reported, never copied to each. Violations: `pg_unit_samples` (qpx/core/sample_sql.py:93-103), behind `ProteinView` and `ds.intensity`, copies it; views.md:214 advises `agg="sum"` for CDAP per-fraction rows that no longer exist.

**PGV-R5 (MUST)** Pivot by `pg_group_id` across units only when provenance `pg_grouping` is `dataset_inference` or `mapping_set` (PGV-6); otherwise warn. Violation: qpx/core/dataset_queries.py:88-98 and qpx/mudata.py:230 pivot unconditionally.

**PGV-R6 (MUST)** Minimum-peptide filters use the unit's `peptide_counts` or the analysis `peptides`. Violation: mokume DEqMS counts sequences per `anchor_protein` over all feature rows, decoys and transfers included (mokume:python/mokume/pipeline/stages.py:1772-1776).

**PGV-R7 (MUST)** Decoys and contaminants per [T10](t10-protein-mapping.md) PRT-R2/R3; summaries report null `contaminant` as unknown. Violation: qpx views ignore `contaminant` (qpx/views/api.py:32-67,247-268).

## Open questions for maintainers
1. Does FragPipe `<exp> Total Peptides` count MBR-only peptides, as MaxQuant `Peptides <exp>` does? Then FragPipe unit counts become null too.
2. Confirm PD `Found in Sample` and `Protein FDR Confidence: Combined` on real exports; no PD fixture exists.
3. Tool units pooling samples (MaxQuant experiment of replicates): keep with a warning (PGV-5), or refuse unless the SDRF agrees?
4. Should unattributable protein lists get a dataset-scope view instead of no rows?
5. [#282](https://github.com/bigbio/qpx/issues/282): a feature flag "used for this pg quantity"? Rename `pg_qvalue` to `pg_unit_qvalue` in 2.0?
