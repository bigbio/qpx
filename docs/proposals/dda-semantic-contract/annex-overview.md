# Annex overview: producer-neutral analysis

[Overview](index.md) · [Strict profile](strict-profile.md) · [Consumers](consumers.md) · [OpenMS work list](openms-worklist.md)

!!! warning "Status: proposal, not normative"
    This is the annex of the [DDA semantic contract proposal](index.md). It analyses qpx `dev` at commit `32cbbdc` (specification 1.3) for every producer and proposes rules for specification 1.4. For OpenMS-written files the [strict profile](strict-profile.md) is the contract, and it wins where the two differ. The rules here are the relaxations other producers need.

This page was the first overview of the proposal. It covers every DDA converter (MaxQuant, FragPipe, Proteome Discoverer, OpenMS consensusXML, mzIdentML, CDAP, quantms-MSstats) on equal terms, which is why many rules allow null or `unknown`.

## What the audit found

Some examples. Each number was measured on the repository's own fixtures with qpx `dev`.

- `feature.confident()` keeps **0 of 243,821** MaxQuant features, and `psm.confident()` keeps **0 of 13,386** Mascot PSMs. Both inputs were filtered upstream (MaxQuant at its PSM FDR, Mascot at its p < 0.01 homology threshold), but no typed q-value was written, and 11,517 of the Mascot rows carry `pass_threshold = true`. ([VALIDATOR-06](findings-psm.md#validator-06), [VALIDATOR-07](findings-feature.md#validator-07))
- A Mascot mzIdentML file gives 13,386 PSM rows for 9,707 spectra. 3,679 rows have rank > 1, and 1,810 of those pass the threshold. `ModificationView` reports 4,323 deamidation sites (on 3,279 rows), against 926 sites (on 827 rows) when only accepted rank-1 rows are counted. ([MZID_CDAP-01](findings-psm.md#mzid_cdap-01), [VALIDATOR-M3](findings-psm.md#validator-m3))
- 43% of MaxQuant feature rows are match-between-runs transfers. In 9,272 target (peptidoform, charge, run) keys a measured direct row and a measured transferred row coexist, and in 799 of them the transferred value is larger. The MuData export sums 12,672 duplicate cells by default, mostly direct-plus-transferred pairs, while mokume takes the maximum. ([VALIDATOR-19](findings-feature.md#validator-19), [MOKUME-09](findings-feature.md#mokume-09))
- `feature.precursor_qvalue` holds a run-level precursor FDR for DIA-NN, the best PSM's q-value for Proteome Discoverer and CDAP, and nothing for MaxQuant, FragPipe and OpenMS. ([INFRA-02](findings-feature.md#infra-02))
- The meaning of MaxQuant `unique` depends on whether proteinGroups.txt was supplied. mokume reads `unique` as "group-specific" and drops 43,823 measured MaxQuant rows. ([MOKUME-01](findings-feature.md#mokume-01), [MOKUME-02](findings-feature.md#mokume-02))
- One MaxQuant isotope pattern is credited in full to two or three peptidoforms (312 cases). On 83,413 rows, `rt` lies outside `[rt_start, rt_stop]` because the window is on another clock. ([MAXQUANT-08](findings-feature.md#maxquant-08), [MAXQUANT-10](findings-feature.md#maxquant-10))
- No footer says whether a file was FDR-filtered, at which level or threshold, whether decoys or transfers were kept, or what kind of intensity it holds. Only the xic view declares its selection. ([INFRA-20](findings-cross.md#infra-20))

In total the audit produced **443 verified findings**: 268 confirmed, 119 confirmed with a correction, and 56 added by the verifiers. None was refuted. They are collected in [Findings: PSM](findings-psm.md), [feature](findings-feature.md), [protein group](findings-pg.md) and [cross-view](findings-cross.md). The 174 converter and tooling defects among them are in the [worklist](converter-worklist.md).

## The four questions

### How many identifications per feature are expected?

- **One precursor ion per row.** A feature row is one quantity of one precursor ion (one canonical peptidoform at one charge ≥ 1) in one run, or in one declared multi-run unit ([FRW-1](t05-feature-rows.md)). A row never holds several peptidoforms.
- **Supporting PSMs.** A row has 0..n *supporting PSMs*, linked by `psm.feature_id` ([FID-7](t06-identifications-per-feature.md), [FID-8](t06-identifications-per-feature.md)). These are rank-1 PSMs of the row's peptidoform (MS1-label partners aside), charge and decoy status, from its run (unit rows: any run in `grouped_runs`). Transferred rows have none. A producer that writes both views links every supporting PSM, and records gaps in provenance ([FID-8](t06-identifications-per-feature.md), [FID-9](t06-identifications-per-feature.md)).
- **At most one identifying PSM.** It comes from the row's own run ([FID-1](t06-identifications-per-feature.md), [FID-2](t06-identifications-per-feature.md)). If the tool names one, that PSM is used. Otherwise, among the run's supporting PSMs (on a `direct` row, those with `precursor_index = 1` first), the order is lowest PEP, then lowest q-value, then the engine's documented primary score, then scan, then peptidoform ([FID-3](t06-identifications-per-feature.md)).
- **What it supplies.** Only the identifying PSM supplies `scan` (exactly one spectrum), `id_run_file_name` (equal to `run_file_name`, or null), PEP, the new `psm_qvalue` and the PSM-level scores ([FID-4](t06-identifications-per-feature.md) to [FID-6](t06-identifications-per-feature.md), [CNF-8](t04-confidence-fields.md)).
- **Rows without one.** Transferred rows, rows of a multi-run unit and rows whose run has PSMs only of competing peptidoforms have no identifying PSM.
- **Rows per precursor and run.** How many rows one precursor gets in one run follows from the file's `identity_composite`. With exactly `peptidoform, charge, run_file_name` there is at most one row; with signal locators in the composite, rows sharing a key are distinct signals that readers sum ([FRW-3](t05-feature-rows.md), [FRW-4](t05-feature-rows.md)). Each quantified signal feeds at most one row ([FRW-2](t05-feature-rows.md)).
- **Conflicts.** When several peptidoforms claim one signal, at most one row keeps the intensity: the source's choice, else the FID-3 winner, or none if the producer drops the signal and counts it. That row and any loser rows written (MAY, with null intensities) set the new flag `identification_conflict` ([FID-10](t06-identifications-per-feature.md)).

### What if the spectrum is chimeric?

- **One peptidoform per row.** No ProForma `+` joiner is allowed in `peptidoform` ([CHI-1](t02-chimeric-spectra.md)). Co-identifications are separate psm rows that share `(run_file_name, scan)`.
- **`precursor_index`.** A new typed psm column tells co-identifications from alternatives ([CHI-3](t02-chimeric-spectra.md)):
    - rows with the same index are hypotheses for one precursor (lower ranks, ties, other charges, other engines);
    - rows with different indices are different co-isolated precursors;
    - index 1 is the precursor the tool designates, else the one whose best rank-1 row sorts first under FID-3; it is a representative, not necessarily the isolation target;
    - null means unknown, never 1.
- **Rank is per precursor ion.** Each co-identified precursor's best hit therefore has rank 1 ([PSM-3](t01-psm-rows.md), [CHI-4](t02-chimeric-spectra.md)). `observed_mz` and `charge` describe the co-isolated ion the row claims ([CHI-7](t02-chimeric-spectra.md)).
- **Feature type.** A feature is `secondary` only when every supporting PSM in its run is a secondary PSM ([CHI-9](t02-chimeric-spectra.md)).
- **Isobaric data.** The reporter ions of a chimeric spectrum feed at most one feature ([CHI-11](t02-chimeric-spectra.md)).
- **No footer key.** A null `precursor_index` means unknown. Removing secondary PSMs is recorded in provenance ([CHI-10](t02-chimeric-spectra.md)).
- **Counting.** Readers count spectra as distinct `(run, scan)` and PSMs as distinct `(run, scan, precursor_index)` ([CHI-R1](t02-chimeric-spectra.md), [CHI-R2](t02-chimeric-spectra.md)).

### What if there is no identification for a feature?

There are three cases:

1. **An MS1 signal that nothing identifies is not a feature row.**
    - Every feature row has a canonical peptidoform, a sequence and a charge ≥ 1, checked at the basic validation tier ([UNI-1](t07-unidentified-features.md)). `identification_type` gets no `unidentified` value ([UNI-2](t07-unidentified-features.md)).
    - Such signals MAY go to a new optional view, `unidentified_feature`, with derived ids ([UNI-3](t07-unidentified-features.md), [UNI-4](t07-unidentified-features.md)).
    - The presence of that view says they were kept. Features that were not exported are counted in provenance ([UNI-5](t07-unidentified-features.md)).
    - This reverses the placement agreed in [#229](https://github.com/bigbio/qpx/issues/229)/[#232](https://github.com/bigbio/qpx/issues/232) and settles the open point in [#378](https://github.com/bigbio/qpx/issues/378). Maintainers need to confirm this.
2. **A row with an identity but no identifying PSM in its own run** is a transferred row, a multi-run unit row, or a row whose run has only PSMs of competing peptidoforms (`transferred` or null type, [FID-1](t06-identifications-per-feature.md)). A `direct` row whose identifying spectrum cannot be located looks the same ([FID-6](t06-identifications-per-feature.md)). In all these rows `scan` is empty and `id_run_file_name` is null, and only `identification_type` tells them apart. Readers never infer a transfer from those nulls ([FID-R3](t06-identifications-per-feature.md)).
3. **An identified precursor without a quantity** is an identification-only row whose entries have a null intensity ([INT-6](t09-intensities-missing-values.md)). The footer key `unquantified_policy` says whether such rows were kept, so a missing row means "not identified" only when the producer says so ([UNI-6](t07-unidentified-features.md), [UNI-R4](t07-unidentified-features.md)).

### How are MBR and identification in the same run told apart?

- **`identification_type` alone gives the origin** ([MBR-1](t08-match-between-runs.md)).
    - `direct` or `secondary`: an accepted identification from a spectrum of this run (or of a run of the row's unit).
    - `transferred`: the signal was measured here, but its identity came only from other runs, whatever MS/MS this run has. MaxQuant `MULTI-MATCH-MSMS` rows are `transferred`.
    - The value comes from the tool's per-row statement. Without one, a producer writes `transferred` only when every identification of the run is visible and none is of the row's peptidoform and charge, and otherwise null; a producer whose rows all rest on own-run identifications writes `direct` or `secondary`. The value never comes from a null run, scan or PEP ([MBR-2](t08-match-between-runs.md)).
- **A transferred row has no identifying PSM** ([MBR-3](t08-match-between-runs.md)): `scan` is empty; `id_run_file_name`, PEP, `psm_qvalue` and run-scope q-values are null; no PSM links to it; `spectral_count` is 0.
- **`match_*` is the transfer's own confidence** ([MBR-4](t08-match-between-runs.md)). `match_qvalue` is a float64 q-value in [0, 1]. The deltas are acceptor minus donor, on the aligned clock. MaxQuant's `Match score` is the donor's identification score and moves to `maxquant_match_score`.
- **Direct and transferred rows of one precursor stay separate.** A direct and a transferred row of the same precursor in the same run are distinct signals ([MBR-7](t08-match-between-runs.md)).
- **No combining.** When readers combine rows, a transferred entry is dropped wherever a `direct` or `secondary` row of the same key and label has a measured (non-null, non-imputed) value. It is never summed, averaged or maxed with that value ([MBR-R4](t08-match-between-runs.md), [RDR-5](t14-reader-contract.md)).
- **No evidence from transfers.** Transfers never count as identifications or protein evidence ([MBR-R3](t08-match-between-runs.md)). The pg view gains `identified_in_unit`, so transfer-only protein rows are visible ([MBR-10](t08-match-between-runs.md), [PGV-1](t12-protein-groups.md)).
- **No footer key.** A null `identification_type` means unknown. Whether MBR was performed, and its FDR, are recorded in provenance; each row keeps its `match_qvalue` ([MBR-8](t08-match-between-runs.md), [MBR-9](t08-match-between-runs.md)).

## Decisions by topic

Each topic has its own page with the current state, options, numbered rules, schema changes, producer changes, reader rules and open questions.

| Topic | Decision in one line |
|---|---|
| [T01 PSM rows](t01-psm-rows.md) | One row per hypothesis per spectrum, unique on K = `(run_file_name, scan, peptidoform, charge)`, with multi-engine hits merged. A typed `rank` (per precursor ion; ties share 1; null means unknown). Rank stays out of the identity, and only rank-1 rows are identifications. |
| [T02 Chimeric spectra](t02-chimeric-spectra.md) | One peptidoform per row, no `+`. Typed `psm.precursor_index`. `secondary` defined. A chimeric spectrum's reporters feed at most one feature. |
| [T03 Filtering state](t03-filtering-state.md) | Every psm, feature and pg footer states its FDR cuts in one key, `fdr_filter`. Decoy, contaminant and unassigned-PSM handling go to provenance. `pass_threshold` is the source's acceptance, with its basis in provenance. qpxc still removes no rows by default. |
| [T04 Confidence fields](t04-confidence-fields.md) | The column name fixes the level. New `feature.psm_qvalue`, so a best-PSM q-value is never promoted to precursor or peptide level. Typed columns are filled only from level-explicit names. Feature PEP is the identifying PSM's. `confident()` takes an explicit level and never returns a silently empty set. |
| [T05 Feature rows](t05-feature-rows.md) | A row is one quantity of one precursor ion in one run or declared unit. A signal feeds at most one row. The `identity_composite` says whether keys repeat. New `grouped_runs` for unit rows. All times are in seconds: `rt`, `rt_start` and `rt_stop` on the run's acquisition clock, `rt_aligned` and `consensus_rt` on the tool's aligned clock. |
| [T06 Identifications per feature](t06-identifications-per-feature.md) | 0..n supporting PSMs and at most one identifying PSM from the row's own run, chosen by a total order. `scan` is one spectrum. New `identification_conflict` and `consensus_feature_id`. Links are mandatory when a producer writes both views. |
| [T07 Unidentified features](t07-unidentified-features.md) | Feature rows are always identified. Unidentified MS1 features go to an optional `unidentified_feature` view. The footer key `unquantified_policy` says whether identification-only rows were kept. |
| [T08 Match-between-runs](t08-match-between-runs.md) | `identification_type` alone gives the origin. Transfers carry no PSM and no run-level confidence. `match_*` is defined precisely. pg gains `identified_in_unit`. MBR settings go to provenance. |
| [T09 Intensities and missing values](t09-intensities-missing-values.md) | The primary is the least-processed quantity of one kind per file, declared by `intensity_quantity` and `intensity_processing`. One entry per sample label (null = not quantified, 0 = measured zero). New `quant_spectrum_count`. A QPX roll-up is primary only when the tool gives no protein quantity. |
| [T10 Protein mapping](t10-protein-mapping.md) | `protein_accessions` is the set of mapped entries. `pg_accessions` is one group's full membership or null. `anchor_protein` is a representative, never a key. New `group_unique`, from which `unique` is derived. A peptide is a decoy only if every entry is. Typed `contaminant` on psm and feature. Markers are kept verbatim. |
| [T11 Spectrum and precursor fields](t11-spectrum-precursor-fields.md) | `observed_mz` is the measured m/z of the claimed ion, with its calibration recorded in provenance. `mass_error_ppm` is computed from the row; a typed `isotope_error` is added. `charge ≥ 1`. `scan_format` is required, with no hashes, and `spectrum_reference` sits beside `scan`. Ion mobility is 1/K0. One run-name function. The canonical peptidoform is validated. |
| [T12 Protein groups](t12-protein-groups.md) | Row existence follows one table. Evidence is rank-1, non-transferred and not failing. Counts have fixed scopes. `pg_qvalue` is unit-level and `global_qvalue` analysis-level. A group is a decoy only if every member is. `unquantified_policy`; the grouping method goes to provenance. |
| [T13 Metadata and conformance](t13-metadata-conformance.md) | One footer-key registry with five new keys, read from the key-value metadata. Rewrites carry every key by a rule. `acquisition_method`. Validation profiles `core` and `dda` with error, warning and info severities. Any change of meaning moves `qpx_version`. |
| [T14 Reader contract](t14-reader-contract.md) | A normative `consumer-guide.md` and a default read profile. Decoys are excluded. Accepted rank-1 rows identify. Transfers are never evidence and are never combined with a measured own-run value. Signals and fractions are summed. The group is the protein. Samples resolve exactly through `run.samples`. Absent keys mean unknown. Outputs record the policy applied. |

## Principles

These principles govern every section:

| | Principle |
|---|---|
| P-A | **One meaning per typed column.** If a producer cannot supply exactly that meaning, it writes null and puts its near-equivalent in `additional_scores`, `additional_intensities` or `cv_params` under a name that states its level. |
| P-B | **Declare, don't filter silently.** Producers may keep or drop rows, but the file declares what was kept: in the footer when a reader's default behaviour depends on it, otherwise in provenance. Footer keys use closed vocabularies; an absent key means `unknown`, and readers assume neither way. The xic keys are the precedent. |
| P-C | **Flags on rows over hidden policy.** Rank-n alternatives, co-identifications, transfers, imputed values, decoys and contaminants are marked by typed columns. |
| P-D | **Null means unknown or not applicable.** No sentinels: no 0 m/z, no 1.0 q-value placeholders, no `''` or `'None'` peptidoforms, no hashes read as scans. |
| P-E | **PSI terms as references only.** A PSI-MS term or a ProForma, mzIdentML or mzTab definition is cited where it pins down a meaning. Matching PSI formats is not a goal and never decides a rule. |
| P-F | **Additive and versioned.** Rules target 1.4. New optional columns and footer keys come first, with a stricter conformance profile. Every breaking change is flagged. |
| P-G | **Feasible for producers.** Every MUST can be met from the files converters actually read. Where a tool lacks the information, the rule says what to write. |
| P-H | **Readers are part of the contract.** Every section ends with reader rules and the places where qpx and mokume violate them. |
| P-I | **Minimal.** Use the smallest rule set that removes the ambiguity. |

## Proposed footer keys

The footer holds only what a reader needs to interpret the rows: a key stays in the footer only if a reader's default behaviour depends on its value. That leaves five new keys. Every vocabulary also allows `unknown`, and an absent key means `unknown`. [T13 META-1](t13-metadata-conformance.md) holds the full registry, including the existing keys and the carry rules for merge and migrate.

| Key | Views | Values | Owner |
|---|---|---|---|
| `acquisition_method` | psm, feature, pg | `dda`, `dia`, `prm`, `mixed` | [META-5](t13-metadata-conformance.md) |
| `fdr_filter` | psm, feature, pg | `none`, or a `;`-separated list of `<level>:<value>` such as `psm:0.01;protein:0.01`. The level is `psm`, `precursor`, `peptide` or `protein`; the value is a threshold in (0, 1], `applied` (a cut on another statistic) or `unknown`. A level not listed had no cut. | [FLT-1](t03-filtering-state.md), [FLT-2](t03-filtering-state.md) |
| `unquantified_policy` | feature, pg | `retained`, `removed` | [UNI-6](t07-unidentified-features.md), [PGV-4](t12-protein-groups.md) |
| `intensity_quantity` | feature, pg, site, unidentified_feature | `ms1_area`, `ms1_apex`, `ms1_trigger`, `ms1`, `reporter_uncorrected`, `reporter_corrected`, `reporter` (correction unknown), `reporter_sn`, `fragment`, `none` | [INT-2](t09-intensities-missing-values.md) |
| `intensity_processing` | feature, pg, site, unidentified_feature | `none`, `sum`, `topn`, `normalized`, `maxlfq`, `other`; the most processed step applied | [INT-2](t09-intensities-missing-values.md) |
| `scan_format` (existing, extended) | psm, feature, mz | adds `mixed` (psm only) and `unknown`; required from 1.4 on psm, and on feature when any `scan` is non-empty | [SPC-7](t11-spectrum-precursor-fields.md) |

### Declarations that are not footer keys

An earlier draft of this proposal had about 30 new footer keys. The others are replaced as follows:

| Earlier key | Replaced by |
|---|---|
| `psm_rank_policy` | The typed `rank` column; null means unknown ([PSM-4](t01-psm-rows.md)). |
| `chimeric_policy` | The typed `precursor_index` column; null means unknown. Removing secondary PSMs is recorded in provenance ([CHI-10](t02-chimeric-spectra.md)). |
| `feature_grain` | The `identity_composite`. With exactly `peptidoform, charge, run_file_name` there is at most one row per key; otherwise rows sharing a key are distinct signals ([FRW-3](t05-feature-rows.md)). |
| `intensity_imputation` | Row encoding: `imputed_labels = []` means nothing imputed, null means unknown ([INT-9](t09-intensities-missing-values.md)). |
| `mbr_policy` | `identification_type`, filled whenever the tool states the origin (null means unknown), and pg `identified_in_unit` ([MBR-8](t08-match-between-runs.md)). |
| `psm_feature_link` | Links are mandatory when a producer writes both views; gaps are recorded in provenance ([FID-8](t06-identifications-per-feature.md), [FID-9](t06-identifications-per-feature.md)). |
| `chimeric_quant_policy` | A fixed rule ([CHI-11](t02-chimeric-spectra.md)). |
| `conformance_profile` | No footer key. A dataset claims the [strict profile](strict-profile.md#footer) with a provenance parameter of the same name; `qpxc validate --profile` checks ([META-7](t13-metadata-conformance.md)). |
| `fdr_filter_psm`, `_precursor`, `_peptide`, `_protein` | Merged into `fdr_filter`. |
| `unquantified_feature_policy`, `unquantified_pg_policy` | Merged into `unquantified_policy`. |
| `intensity_normalization`, `intensity_rollup`, `reporter_correction` | Merged into `intensity_processing` and the reporter values of `intensity_quantity`. |
| `unassigned_psm_policy`, `decoy_policy`, `contaminant_policy`, `pass_threshold_level`, `pass_threshold_fdr`, `psm_qvalue_scope`, `peptide_qvalue_scope`, `mbr_fdr`, `observed_mz_calibration`, `pg_grouping`, `unidentified_feature_policy` | Provenance parameters with the same names and values. Merges keep the provenance of every input ([META-4](t13-metadata-conformance.md)). |

The cost is that some states can no longer be read from the file alone. For example, "secondary PSMs were removed" looks like "no chimeric search", and "unassigned PSMs were dropped" looks like "there were none". Readers who need these consult provenance.

## New and changed columns

| Column | View | Type | Owner | Breaking? |
|---|---|---|---|---|
| `rank` | psm | int16, null, ≥ 1 | [PSM-3](t01-psm-rows.md) | No. The score `rank` also stays in `additional_scores` during 1.4 ([CNF-6](t04-confidence-fields.md)) |
| `precursor_index` | psm | int16, null, ≥ 1 | [CHI-3](t02-chimeric-spectra.md) | No |
| `spectrum_reference` | psm | string, null; stands in for an empty `scan` in K | [SPC-9](t11-spectrum-precursor-fields.md) | Ids of empty-scan rows change |
| `isotope_error` | psm | int16, null | [SPC-5](t11-spectrum-precursor-fields.md) | No |
| `compensation_voltage` | psm | float32, null | [SPC-12](t11-spectrum-precursor-fields.md) | No |
| `contaminant` | psm, feature | bool, null | [PRT-6](t10-protein-mapping.md) | No |
| `pass_threshold` | psm (nullable), pg (new) | bool, null | [FLT-7](t03-filtering-state.md), [PGV-9](t12-protein-groups.md) | mzIdentML "no threshold" rows become null |
| `psm_qvalue` | feature | float64 [0, 1], null | [CNF-1](t04-confidence-fields.md) | No; receives values that 1.3 wrote at the wrong level |
| `grouped_runs` | feature | list<string>, null | [FRW-6](t05-feature-rows.md) | No |
| `identification_conflict` | feature | bool, null | [FID-10](t06-identifications-per-feature.md) | No |
| `consensus_feature_id` | feature | int64, null | [FID-11](t06-identifications-per-feature.md) | No |
| `group_unique` | feature | bool, null | [PRT-9](t10-protein-mapping.md) | No; `unique` becomes derived and deprecated |
| `quant_spectrum_count` | feature | int32, null, ≥ 1 | [INT-10](t09-intensities-missing-values.md) | No |
| `identification_type` | feature | YAML enum, checked on every write path | [MBR-1](t08-match-between-runs.md) | Writers with other values |
| `match_qvalue` | feature | float32 → float64 [0, 1] | [CNF-9](t04-confidence-fields.md), [MBR-4](t08-match-between-runs.md) | **Physical type** |
| `identified_in_unit` | pg | bool, null | [MBR-10](t08-match-between-runs.md) | No |
| `level` | ontology | enum | [CNF-5](t04-confidence-fields.md) | No |
| New view `unidentified_feature` | – | see T07 | [UNI-3](t07-unidentified-features.md) | No |

## Breaking changes

Additive changes come first. The following changes alter values or meanings and need `qpx_version` 1.4 with `qpxc migrate` support ([META-10](t13-metadata-conformance.md)):

- **Ids.**
    - mzIdentML `psm_id` is derived from K, not from SII ids ([PSM-10](t01-psm-rows.md)).
    - Feature ids change for PD and quantms label-free (composite), MaxQuant (window values), OpenMS (`observed_mz`, `scan`) and FragPipe (null CV) ([FRW-4](t05-feature-rows.md)).
    - Runs named by `run_file_stem` change PD `.RAW` and FragPipe dotted-run ids ([SPC-13](t11-spectrum-precursor-fields.md)).
    - psm rows with an empty or hashed `scan` are keyed by `spectrum_reference`, so their ids change ([SPC-9](t11-spectrum-precursor-fields.md)).
    - mzIdentML open-search masses stay four-decimal deltas, which changes those peptidoforms and ids ([SPC-14](t11-spectrum-precursor-fields.md)).
- **Meaning of existing columns.**
    - PD and CDAP `precursor_qvalue` and OpenMS PSM-level `peptide_qvalue` move to `feature.psm_qvalue` ([CNF-2](t04-confidence-fields.md)).
    - `pg_qvalue` is unit-level only: MaxQuant, PD and FragPipe values move or become null ([PGV-8](t12-protein-groups.md)).
    - `match_score` is a transfer score, not MaxQuant's donor score ([MBR-4](t08-match-between-runs.md)).
    - The primary intensity is "least processed, one kind per file, declared" instead of "raw" ([INT-1](t09-intensities-missing-values.md)).
    - MaxQuant windows move to the acquisition clock ([FRW-7](t05-feature-rows.md)).
    - `MULTI-MATCH-MSMS` rows lose scan, run and PEP ([MBR-3](t08-match-between-runs.md)).
    - `match_qvalue` becomes float64 in [0, 1], a physical type change ([CNF-9](t04-confidence-fields.md)). OpenMS peptide-level feature PEP becomes null ([CNF-8](t04-confidence-fields.md)).
    - mzIdentML rows under "no threshold" get a null `pass_threshold` ([FLT-7](t03-filtering-state.md)).
    - `consensus_support` is an engine count: the OpenMS fraction moves to `consensusid_support`, and engine prefixes replace `_2` suffixes ([PSM-2](t01-psm-rows.md)).
    - `id_run_file_name` is `run_file_name` or null, never a donor run ([FID-6](t06-identifications-per-feature.md)).
    - Target rows hold no decoy entries in `protein_accessions` ([PRT-2](t10-protein-mapping.md)).
    - `scan_format` gains `mixed` and `unknown`, and from 1.4 an absent key means `unknown` ([SPC-7](t11-spectrum-precursor-fields.md)).
    - pg: a QPX roll-up may fill the primary only when the tool gives none, declared as `intensity_processing = sum` ([INT-4](t09-intensities-missing-values.md)). Counts get fixed scopes, `peptides` becomes nullable, `is_decoy` needs every member, `contaminant` may be null, rows follow tool units, and `pg_names` align with members ([PGV-5](t12-protein-groups.md), [PGV-7](t12-protein-groups.md), [PGV-10](t12-protein-groups.md), [PGV-12](t12-protein-groups.md)).
- **Rows.**
    - Unidentified rows leave the feature view ([UNI-1](t07-unidentified-features.md)).
    - K-duplicate PSMs are merged ([PSM-2](t01-psm-rows.md)).
    - OpenMS `feature.scan` holds one spectrum ([FID-5](t06-identifications-per-feature.md)).
    - Identification-only entries and rows are added where producers dropped them ([INT-5](t09-intensities-missing-values.md), [UNI-6](t07-unidentified-features.md)).
- **Names.**
    - `lfq` becomes `maxlfq` ([INT-12](t09-intensities-missing-values.md)).
    - Bare `qvalue`, `q_value`, `q-value` and `fdr` score names are deprecated: a validator warning, and an error under the `dda` profile ([CNF-5](t04-confidence-fields.md)).
    - Typed-column scores (q-values, `rank`) stay in `additional_scores` as well during 1.4. Their removal is deferred to a later minor version ([CNF-6](t04-confidence-fields.md)).
- **Validation.**
    - `charge = 0`, `''` and `'None'` peptidoforms, and the ProForma `+` joiner become errors ([UNI-1](t07-unidentified-features.md), [SPC-6](t11-spectrum-precursor-fields.md), [CHI-1](t02-chimeric-spectra.md)).
    - Writers raise on value violations instead of persisting warnings ([META-9](t13-metadata-conformance.md)).
- **qpx API.**
    - `Feature.confident()` takes an explicit level, `PSM.confident()` respects rank, and both report unassessable rows and raise when every row is unassessable, instead of silently returning an empty set ([CNF-R4](t04-confidence-fields.md), [RDR-3](t14-reader-contract.md)).
    - MuData no longer sums direct and transferred cells ([RDR-5](t14-reader-contract.md)).
    - Lone-sample label fallbacks are removed ([RDR-8](t14-reader-contract.md)).
    - Counts follow RDR-9, MuData `obs` is keyed on the unit, `for_quantification` applies the default profile, and a null `is_decoy` raises ([RDR-4](t14-reader-contract.md), [RDR-9](t14-reader-contract.md), [FRW-R4](t05-feature-rows.md)).
- **Rewrites.**
    - merge, relayout and write-back refuse a source at another `qpx_version`, and a merged `creator` becomes a list ([META-4](t13-metadata-conformance.md)).

## Cross-section decisions taken here

The harmonization pass left ten conflicts between sections open. This overview decides them as follows.

| Conflict | Decision |
|---|---|
| `rank` both as a typed column and as a score ([#377](https://github.com/bigbio/qpx/issues/377)) | Treat it like the q-values: the writer fills `psm.rank` from the score and keeps the entry during 1.4 ([CNF-6](t04-confidence-fields.md)). In 1.3 files the score is the rank. |
| Scope (run vs experiment) of `fdr_filter` thresholds | No new footer key. The scope goes in the provenance filtering step as `<level>_fdr_scope` ([FLT-4](t03-filtering-state.md)). Readers that need it and find none follow [CNF-R3](t04-confidence-fields.md). |
| `PSM.confident()` versus the acceptance predicate | Two helpers with two meanings. `confident(t)` is a q-value filter at PSM level ([CNF-R4](t04-confidence-fields.md)). `accepted()` is the source's acceptance ([FLT-R4](t03-filtering-state.md)). Counts use `accepted()`. |
| pg evidence: the producer counts `pass_threshold ≠ false`, the reader needs a declared cut | Keep the producer rule ([PGV-2](t12-protein-groups.md)) and require the pg footer to declare the cut behind it ([PGV-3](t12-protein-groups.md)). The two rules then agree. |
| Multi-run unit rows (FragPipe) and the identifying PSM | Unit rows name no identifying PSM. `id_run_file_name` stays "`run_file_name` or null" ([FID-1](t06-identifications-per-feature.md), [FRW-6](t05-feature-rows.md)). Supporting PSMs from member runs may still link. |
| Merging files with different `scan_format` | Write `mixed` if every psm row has a `spectrum_reference`, else `unknown` with a warning ([SPC-7](t11-spectrum-precursor-fields.md), [META-4](t13-metadata-conformance.md)). |
| PSM counts when `precursor_index` is null | Count once per spectrum and report the coverage. The count is a lower bound on identified precursors ([PSM-R2](t01-psm-rows.md)). |
| Should the `dda` profile require rank-1-only PSM files, or reject `unknown` footer values? | No, not in 1.4. Requiring rank-1-only files would reverse [#241](https://github.com/bigbio/qpx/issues/241); rank is a typed column and readers count `rank = 1`. Rejecting `unknown` should wait until every qpx converter fills the keys. |
| MaxQuant `MULTI-MATCH-MSMS` | Keep it as `transferred` ([MBR-1](t08-match-between-runs.md)) until MaxQuant documentation or a fixture shows that its MS/MS identifies the peptide at PSM FDR. |
| Footer keys of `unidentified_feature` | `intensity_quantity`, `intensity_processing` and the identity keys. Isobaric runs are out of scope for this view in 1.4 ([UNI-3](t07-unidentified-features.md)). |

## Producer status today

This is a summary of the per-section producer tables, not a test result. ✓ means it already conforms or nearly does. ~ means it partly conforms. ✗ means it diverges or the information is missing. – means not applicable.

| Producer | Rank | Chimeric | Filter state | Q-value level | Feature rows | Identifying PSM | ID-only / unidentified | MBR | Intensities | Proteins | Spectrum fields | pg |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| MaxQuant | ~ top hits, no rank | ✗ second peptides only in cv_params; duplicate ids | ✗ cut, not declared | ✗ no q-values; pg copy of experiment value | ~ peaks, windows on another clock | ✓ `Best MS/MS` | ✓ keeps ID-only rows | ~ `MULTI-MATCH-MSMS`, unmapped types | ~ per-row column choice; TMT without SDRF as LFQ | ✗ `unique` depends on input; contaminant flag unread | ✗ `observed_mz` derived | ✗ |
| FragPipe | ✗ unknown | ✗ | ✗ cut, not declared | ✗ PEP = 1 − probability | ✗ experiment rows on first run | ~ qpx pick; CV ignored | ✗ unquantified dropped | ✓ `Match Type` | ✗ pg mixes columns per cell | ✗ `contam_` stripped | ✗ dotted runs truncated | ✗ |
| Proteome Discoverer | ~ all ranks as scores | ✗ shared reporters | ~ `Confidence` only | ✗ PSM q as precursor q | ~ composite holds picked `rt` | ~ q-value pick | ✗ dropped | ~ always direct | ✗ MAX over PSMs, kind varies by export; double sums | ✗ `unique`, anchor, decoys | ~ First Scan only | ✗ |
| OpenMS consensusXML | ✗ rank unread | ✗ | ~ unfiltered, undeclared | ✗ PSM q as peptide q; 1.0 placeholders | ~ composite includes `scan` | ✗ first PID; scans concatenated | ✗ dropped uncounted | ~ heuristic | ✗ channels omitted; QPX pg sums | ~ decoy rule differs by view | ✗ hashes, 0.0 m/z, 15-decimal masses | ✗ |
| OpenMS native (`-out_qpx`) | ~ best hit only, no `rank` column | ~ separate rows, unmarked | ✗ not declared | ✗ PSM q-value not written | ~ | ✗ best hit across runs (fixed on the #378 branch) | ✗ empty peptidoforms | ✗ no type (the #378 branch adds it) | ~ measured, undeclared; entries cannot be null | ✗ mapped proteins as `pg_accessions` | ✗ aligned `rt`, theoretical m/z | ~ own abundance; evidence counts transfers |
| mzIdentML | ~ all ranks as scores | ✗ ties vs co-identifications | ~ `passThreshold`, threshold dropped | ✗ p-value as PEP | – | – | – | – | – | ✗ target+decoy as decoy | ✗ vendor components dropped | ✗ |
| CDAP | ~ ties only in cv_params | ✓ one precursor | ✗ cut, not declared | ✗ PSM q as precursor q | ~ PSM roll-up | ✓ SpecEValue order | ✓ keeps ID-only rows | ✓ direct | ✗ QPX pg sums | ✗ `XXX_` stripped; contaminant constant | ✗ 0 sentinels | ✗ |
| quantms-MSstats | – | – | ~ undeclared | – | ✓ | ~ claims own run | ✗ aborts on NA | ~ null type | ~ undeclared | ✗ `unique` scope mismatch | ✗ mass shifts | – |

Some current behaviour should be kept and written into the contract: the canonical peptidoform rule (and FragPipe's and PD's conformance to it), IonQuant `Match Type` mapping, MaxQuant's SILAC lightest-partner rule, the null-versus-zero rule, all mapped proteins on PSMs, and the rules the validator already enforces ([SPEC-37](findings-cross.md#spec-37), [FRAGPIPE-26](findings-psm.md#fragpipe-26), [FRAGPIPE-14](findings-feature.md#fragpipe-14), [MAXQUANT-31](findings-feature.md#maxquant-31), [EXTERNAL-30](findings-feature.md#external-30), [EXTERNAL-32](findings-psm.md#external-32), [VALIDATOR-33](findings-cross.md#validator-33)).

## Suggested order of work

This is the producer-neutral order from the first draft. The current order, which puts OpenMS and the consumers first, is on the [overview](index.md#order-of-work).

1. **Plumbing, no change of meaning.** Add the footer registry, writers that keep the Arrow schema equal to the key-value footer, and readers and rewrites that read the key-value footer and carry keys by rule ([META-1](t13-metadata-conformance.md) to [META-4](t13-metadata-conformance.md)). Everything else depends on declarations surviving a merge.
2. **Reader fixes that are wrong under any reading of 1.3.** These are in qpx's read API, the MuData export and mokume ([T14](t14-reader-contract.md)):
    - `confident()` mixing levels or returning empty sets;
    - direct and transferred cells summed;
    - a null `is_decoy` read as target;
    - anchor-keyed protein outputs;
    - zeros read as missing;
    - lone-sample label fallbacks.
3. **Additive 1.4 columns, keys and the `unidentified_feature` view**, with the writer's allowlist fill ([CNF-6](t04-confidence-fields.md)).
4. **Converter work**, by producer, from the [worklist](converter-worklist.md).
5. **Breaking semantic changes** behind `qpx_version` 1.4, with `qpxc migrate`.
6. **`consumer-guide.md`, the `dda` profile and the reader-contract fixtures** ([RDR-1](t14-reader-contract.md), [RDR-14](t14-reader-contract.md), [META-13](t13-metadata-conformance.md), [META-14](t13-metadata-conformance.md)). mokume can then run the same fixtures.

The fixtures lack real Proteome Discoverer exports, FragPipe with IonQuant match-between-runs, a DDA+ or Sage chimeric search and a two-engine mzIdentML file. Several rules are marked "to confirm" until such exports exist.

## Open questions for maintainers

Each section lists its own open questions. These are the ones that shape the contract:

1. **Unidentified features.** Should they go to a separate view, as proposed, which reverses [#229](https://github.com/bigbio/qpx/issues/229)/[#232](https://github.com/bigbio/qpx/issues/232)? Or should they stay in `feature` with nullable identity columns? This also settles [#378](https://github.com/bigbio/qpx/issues/378). ([T07](t07-unidentified-features.md))
2. **Transfers in quantities.** The default read profile keeps transferred quantities, de-duplicated and masked, because MaxQuant and IonQuant protein values include them. Should they be dropped by default, or weighted (mokume#117)? ([T14](t14-reader-contract.md), [T08](t08-match-between-runs.md))
3. **Conflicts.** Should loser rows of an identification conflict be MUST rather than MAY? Do co-eluting localisation isomers count as conflicts ([T06](t06-identifications-per-feature.md))?
4. **`peptide_qvalue`.** Is it per stripped sequence (MS:1001868) or per peptidoform? Fixing one meaning for the column (probably per peptidoform) avoids a footer key ([CRITIC-34](findings-feature.md#critic-34)).
5. **2.0 clean-up.**
    - Retype `feature.pg_accessions` to `list<string>`.
    - Remove `unique`.
    - Move imputed values out of the primary slot.
    - Rename `feature.global_qvalue` to `precursor_global_qvalue` and `pg_qvalue` to `pg_unit_qvalue`.
6. **Facts to confirm with tool authors or real exports.**
    - MaxQuant `MULTI-MATCH-MSMS`.
    - FragPipe DDA+ and CHIMERYS per-precursor m/z.
    - PD `Found in Sample` and decoy content.

### Loose ends from the OpenMS corrections

Correcting the OpenMS facts in the topic sections left three points open:

1. **Feature id on the consensusXML route.** [FRW-4](t05-feature-rows.md) builds the OpenMS feature id from `rt`, `observed_mz` and `consensus_rt`. For aligned ProteomicsLFQ output `rt` and `observed_mz` are now null, which leaves only `consensus_rt`. Either `rt_aligned` replaces `rt` there, or this route stays on its current ids.
2. **Detecting mixed protein quantities.** [INT-4](t09-intensities-missing-values.md) forbids a file that mixes tool quantities and QPX sums in `pg.intensity`, but nothing on a row says which kind it holds, so the validator cannot check it.
3. **Feature m/z fallback.** [SPC-1](t11-spectrum-precursor-fields.md) lets a feature fall back to its identifying PSM's m/z, while the OpenMS producer rows say null. One of the two should change.

## Residual items

The completeness pass checked every finding against the sections. The items below are not yet written into a section. Each has a suggested decision and the section that would own it.

| Item | Suggested decision | Owner |
|---|---|---|
| A psm view may cover only some runs (one FragPipe psm.tsv per experiment) [FRAGPIPE-25](findings-cross.md#fragpipe-25) | The psm view covers every run that has direct or secondary features. Otherwise list the uncovered runs in provenance (`skipped_runs_psm`); the validator warns. | T03, T06 |
| Engine attribution of merged multi-engine rows [CRITIC-12](findings-psm.md#critic-12) | Optional `psm.search_engines` (list of engine ids equal to the score prefixes). `consensus_support` equals its length. | T01 |
| Unit of `peptide_qvalue` [CRITIC-34](findings-feature.md#critic-34) | Fix one meaning for `feature.peptide_qvalue` (probably per peptidoform); values of the other kind go under a named score. No footer key. | T04 |
| Peptide-level confidence on psm rows [EXTERNAL-25](findings-psm.md#external-25) | No typed column. Reserved names stay in `additional_scores`, and CNF-6's removal applies only where a typed column exists. | T04 |
| FDR competition unit when several rank-1 rows share a spectrum [EXTERNAL-08](findings-psm.md#external-08) | Q-values are FDR-valid only for the best row of each `(run, scan, precursor_index)` group. Provenance records the competition unit. | T04 |
| pepmap filter state [MZID_CDAP-M2](findings-cross.md#mzid_cdap-m2), [CRITIC-57](findings-cross.md#critic-57) | Population is already set by [PRT-14](t10-protein-mapping.md) (target peptidoforms of rank-1 psm rows and of feature rows). Residual: pepmap's footer copies the psm filter keys (T03 open question 5). | T03, T10 |
| What a group member is (MaxQuant subset proteins) [CRITIC-52](findings-pg.md#critic-52) | Define membership in the spec (indistinguishable proteins only, or with subsets); converters that cannot meet it say so. No footer key. | T12 |
| Feature traceability to tool row ids [MAXQUANT-M2](findings-feature.md#maxquant-m2) | Optional `feature.source_feature_id` (string). Merged ids are joined with `;`. | T06 |
| Absent optional column vs all-null column; stale 1.3 changelog [CRITIC-62](findings-cross.md#critic-62), [SPEC-32](findings-cross.md#spec-32) | An absent optional column means unknown on every row. The changelog lists every 1.3 column. | T13 |
| Label-field invariants [VALIDATOR-16](findings-psm.md#validator-16), [VALIDATOR-22](findings-feature.md#validator-22) | `psm.label` resolves to exactly one sample. `identified_labels`, `representative_label` and `channel_details` labels are subsets of the row's intensity labels (INT-8 and INT-9 already require this of `additional_intensities` and `imputed_labels`). | T09 |
| Precursor purity definitions [EXTERNAL-21](findings-feature.md#external-21) | INT-10 already maps `parent_ion_fraction` to MS:1003167. Residual: PD writes 1 − interference; other definitions go to `<tool>_pif`. | T09 |
| Consumers dropping labels they do not know [MOKUME-27](findings-feature.md#mokume-27) | Readers accept the whole canonical label vocabulary from `run.samples`, and raise on unknown labels. | T14 |

Questions that neither the findings nor the sections cover, with suggested decisions:

| Question | Suggested decision | Owner |
|---|---|---|
| Group-specific FDR (open-search mass-shift bins, glyco classes) | Recorded in the provenance validation step, with the strata. No footer key. | T04 |
| Are fractions of one sample "identified in the sample" after a cross-fraction transfer? | Sample-level identification = any direct or secondary row in any run of the sample. | T14 |
| Enzyme, specificity, missed cleavages, I/L handling of the protein mapping | Recorded in the provenance search step. `pepmap.is_unique` is defined over that mapping. | T10 |
| Database of each accession (target, contaminant, variant) | The provenance search step lists every FASTA with its role, version and checksum. | T10 |
| How confident is a site in `peptidoform`? (SPC-15 already fixes the site as the final localisation.) | The peptidoform holds the best placement, not a confidence claim. Probabilities live in `modifications`, and the site view declares its threshold. | T11 |
| Cross-link FDR levels; cross-links in features | CSM FDR maps to psm. Residue-pair and PPI FDR go to provenance. No cross-links in feature in 1.4. | T03, T05 |
| DDA-PASEF: feature ion mobility, what one spectrum is | Feature ion mobility follows FRW-8. For TDF, the PASEF precursor is the spectrum. CHI-5 also compares mobility. | T05, T11 |
| Fragment annotation notation | `ion_type_array` uses an mzPAF subset. | T11 |
| Re-quantified MS1-label partners | A measured entry flagged `channel_details[].evidence = requantified`, not in `identified_labels`. | T09 |
| Combined designs (SILAC × TMT) | Out of scope for 1.4. Converters reject them with a clear error. | T09 |
| Spike-in standards (iRT, UPS) | Ordinary targets, optionally listed in provenance, never flagged as contaminants. | T10 |
| Decoy construction and competition mode | Recorded in provenance. FDR re-estimation is allowed only when these are known. | T03 |
| Same-run charge-state transfers | Identification is per (peptidoform, charge), so these are `transferred`. | T08 |
| Clock of `predicted_rt` | Seconds on the run clock. Library and iRT values go to `<tool>_irt`. | T05 |
| One vendor file with several acquisitions | A defined stem suffix, or explicit rejection in 1.4. | T11 |

## How this proposal was produced

1. **Audit.** Eleven independent audits covered the OpenMS, MaxQuant, FragPipe, Proteome Discoverer, mzIdentML and CDAP converters, the shared writer and identity code, the specification text, the validator and read API, mokume, external tools and PSI standards, and a first-principles list of reader questions. Each audit converted the repository's fixtures where useful and cited code by file and line.
2. **Verification.** A separate skeptic tried to refute every finding against the code and could add findings that the audit missed.
3. **Design.** Fourteen topic sections were drafted against the principles above, using the findings, the specification and the GitHub history. A second agent reviewed each section adversarially, checking its facts, producer feasibility, unflagged breaking changes and testability, and corrected it in place.
4. **Consolidation.** A harmonization pass removed most conflicts between sections; this overview decides the ten it left open. A completeness pass checked every finding and the four questions above. Three final reviewers checked this overview against the sections, the findings and the code.

Known limits:

- The OpenMS `-out_qpx` C++ writer was not part of this audit. It was read later; see the [OpenMS work list](openms-worklist.md).
- Proteome Discoverer and FragPipe behaviour comes from the converter code and synthetic cases, not from real exports. The repository has no fixture for either.
- Statements about external tools rely on their documentation, which is sometimes silent (MaxQuant `MULTI-MATCH-MSMS`).

## Files

- Fourteen topic sections: [T01](t01-psm-rows.md) to [T14](t14-reader-contract.md).
- Verified findings: [PSM](findings-psm.md), [feature](findings-feature.md), [protein group](findings-pg.md), [cross-view, footer and documentation](findings-cross.md).
- [Converter and tooling worklist](converter-worklist.md), 174 items, each mapped to the rule that fixes it.

## References

- ProForma 2.0: [specification](https://github.com/HUPO-PSI/ProForma) and [paper](https://arxiv.org/abs/2109.11352) (chimeric `+` notation, §7.2).
- mzIdentML 1.2: [schema documentation](https://peptideatlas.org/PSI/schemas/mzIdentML/1.2/mzIdentML1.2.0_combined.html) (`rank`, `passThreshold`, final PSM list MS:1002439) and [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5500760/).
- mzTab: [1.0 specification](https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/specification_document-releases/1_0-Proteomics-Release/mzTab_format_specification.pdf) and [2.0-M specification](https://hupo-psi.github.io/mzTab/2_0-metabolomics-release/mzTab_format_specification_2_0-M_release.html).
- Universal Spectrum Identifier: [1.0 specification](https://raw.githubusercontent.com/HUPO-PSI/usi/master/Specification/USI_FullSpecification_1.0.0_final.pdf).
- PSI-MS controlled vocabulary: [OLS](https://www.ebi.ac.uk/ols4/ontologies/ms).
- SDRF-Proteomics: [specification](https://github.com/bigbio/proteomics-sample-metadata).
- MaxQuant: [output tables](https://cox-labs.github.io/coxdocs/output_tables.html).
- FragPipe: [output files](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html). IonQuant: [repository](https://github.com/Nesvilab/IonQuant) and [MBR FDR method](https://pmc.ncbi.nlm.nih.gov/articles/PMC8131922).
- MSFragger-DDA+: [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11978857/). CHIMERYS: [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12074992/). Sage: [results documentation](https://sage-docs.vercel.app/docs/results/search).
- FlashLFQ: [output formats](https://github.com/smith-chem-wisc/FlashLFQ/wiki/Output-Formats). PIP-ECHO, match-between-runs error control: [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12488043/).
- OpenMS: [IDConflictResolver](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_IDConflictResolver.html), [FeatureFinderIdentification](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderIdentification.html), [ProteomicsLFQ](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_ProteomicsLFQ.html).
- MaxLFQ: [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4159666).
