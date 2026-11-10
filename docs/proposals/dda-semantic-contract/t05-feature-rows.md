# T05. Feature rows: grain, cardinality, quantification units, labels and retention-time fields

[← T04](t04-confidence-fields.md) · [Overview](index.md) · [Strict profile](strict-profile.md) · [Annex overview](annex-overview.md) · [T06 →](t06-identifications-per-feature.md)

!!! note "Annex: producer-neutral analysis"
    Part of the annex of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text, and they cover every producer. For OpenMS-written files the [strict profile](strict-profile.md) is the contract and wins where the two differ. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** A feature row is one quantity of one precursor ion in one run or one declared multi-run unit, or one unidentified MS1 signal marked as such ([T07](t07-unidentified-features.md)), and a signal feeds at most one row. No footer key declares the grain: the declared `identity_composite` says whether `(peptidoform, charge, run_file_name)` can repeat (exactly those columns: at most one row per key; a leading `quantification_unit_id`: one row per unit key; any other composite: distinct signals, which readers sum). A new column `grouped_runs` names a unit row's runs. `rt`, `rt_start` and `rt_stop` are seconds on the run's acquisition clock, and `rt` is the MS1 peak's time, else the identifying spectrum's.

## Questions a reader has
- What is one row: a peak, a precursor per run, one reporter spectrum, or an experiment-level quantity?
- Can `(peptidoform, charge, run_file_name[, label])` repeat? If it does, are the rows parts to add, or alternatives?
- Is `run_file_name` always the measured run? What does `quantification_unit_id` refer to?
- Are charge states, FAIMS CVs and SILAC partners separate rows?
- Is `rt` the apex or the MS2 time? Which clock and unit do `rt_start` and `rt_stop` use?
- Which `identity_composite` may a producer declare?

## What happens today
Spec:
- A row is "one quantified chromatographic peak" (docs/spec/data-model.md:64; [#217](https://github.com/bigbio/qpx/issues/217)), "per upstream Feature entity" (data-model.md:74), "One row per precursor, charge and run" (docs/spec/producer-guide.md:12), "aggregates one or more PSMs" (docs/spec/feature.md:301), or "(peptidoform, charge, run, rt)" (docs/guide/convert.md:1094).
- Only the composite implies the grain (docs/spec/identity.md:27-35), and the validator checks only it (qpx/core/data/footer_validation.py:157-172).
- `run_file_name` is "the run file name that contains the feature" (feature.md:89), and a row is one run (docs/spec/intensities.md:136; docs/spec/anndata.md:140). `quantification_unit_id` is a free string (feature.md:90).
- `rt_start` and `rt_stop` have no unit or clock (feature.md:69-70; qpx/core/data/schemas/feature.yaml:255-260). `rt` is given only as seconds (feature.md:65), while `xic.rt` is raw run time (docs/spec/xic.md:18).
- An MS1 multiplet is one row with the lightest partner's peptidoform (intensities.md:84-92; docs/spec/peptidoform.md:20-23). convert.md:1038-1040 contradicts this for dimethyl and mTRAQ. [#308](https://github.com/bigbio/qpx/issues/308) settled the encoding.
- Duplicate cells are summed (anndata.md:153; docs/spec/views.md:212). Nothing says whether such rows are parts or alternatives.

| Producer | Behaviour | Findings |
|---|---|---|
| MaxQuant | One row per evidence row; rows equal on the composite collapse (qpx/converters/maxquant/feature_adapter.py:236-294). maxquant_full:<br>- 14,086 keys have several rows (max 101), and 12,738 have ≥ 2 quantified rows.<br>- 9,272 target keys (9,337 with decoys) have a quantified direct row and a quantified transferred row (9,639 with decoys, counting `secondary` rows).<br>- 2,594 keys have ≥ 2 quantified direct rows (median apex span 21 min).<br>- 312 isotope patterns give two peptidoforms their full intensity.<br>- 12,758 MS/MS-only rows carry a constant 60 s window.<br>`rt` is raw but the window is calibrated (:574-585): `rt` falls outside it in 83,413 of 243,821 rows, `rt_aligned` in none. | [[MAXQUANT-07](findings-feature.md#maxquant-07)] [[MAXQUANT-08](findings-feature.md#maxquant-08)] [[MAXQUANT-09](findings-feature.md#maxquant-09)] [[MAXQUANT-10](findings-feature.md#maxquant-10)] [[SPEC-01](findings-feature.md#spec-01)] [[VALIDATOR-10](findings-feature.md#validator-10)] [[VALIDATOR-M2](findings-feature.md#validator-m2)] |
| FragPipe | One row per experiment, peptidoform, charge and CV (qpx/converters/fragpipe/feature_adapter.py:52-57,330-365). `run_file_name` is the first run after a lexicographic sort (qpx/converters/fragpipe/pg_adapter.py:96), so `grouped_runs` comes out `[F1, F10, F2]`. On multi-run rows `scan` is `[]` and `id_run_file_name` is null (:340). `rt` is never read. CV is 0.0 when absent (:519). | [[FRAGPIPE-10](findings-feature.md#fragpipe-10)] [[FRAGPIPE-11](findings-feature.md#fragpipe-11)] [[FRAGPIPE-17](findings-feature.md#fragpipe-17)] [[FRAGPIPE-20](findings-pg.md#fragpipe-20)] [[FRAGPIPE-27](findings-cross.md#fragpipe-27)] [[SPEC-M2](findings-feature.md#spec-m2)] [[MOKUME-11](findings-feature.md#mokume-11)] |
| Proteome Discoverer | One row per run, raw sequence, raw modification text and charge (qpx/converters/proteome_discoverer/feature_adapter.py:16-20). TMT sums the reporters, and LFQ takes the largest abundance. `rt` is `Apex RT`, else the chosen PSM's MS2 time (:108-111). The default composite contains that `rt` (:43). Two spellings of one modification give two rows. PIF comes from one PSM. | [[PD-29](findings-feature.md#pd-29)] [[PD-M1](findings-feature.md#pd-m1)] [[PD-18](findings-feature.md#pd-18)] [[PD-19](findings-feature.md#pd-19)] [[PD-M4](findings-feature.md#pd-m4)] |
| OpenMS consensusXML | One row per consensus feature × run (qpx/converters/openms_consensus/feature_adapter.py:734-757): a peak for LFQ, a reporter spectrum for isobaric data. The composite adds rt, scan, observed_mz and consensus_rt (feature_dedup.py:15). LFQ `observed_mz` is the consensus centroid, or 0.0 (:645). Dimethyl and mTRAQ peptidoforms lack their label (:638-639). `rt_aligned` is never set, though alignment is declared: ProteomicsLFQ transforms element and identification RTs in place and stores no `original_RT`, so the converter's `rt` is then an aligned time. | [[OPENMS-09](findings-feature.md#openms-09)] [[OPENMS-15](findings-feature.md#openms-15)] [[OPENMS-29](findings-feature.md#openms-29)] [[OPENMS-M3](findings-feature.md#openms-m3)] |
| quantms MSstats | Composite run, peptidoform, charge, rt, scan (qpx/converters/quantms_msstats/feature_adapter.py:28): one row per PSM scan for isobaric inputs. Label-free and DIA-NN inputs lack `RetentionTime`, so one row per precursor and run; conflicting duplicates abort (:425-437). | [[INFRA-26](findings-feature.md#infra-26)] [[INFRA-29](findings-feature.md#infra-29)] [[CRITIC-39](findings-feature.md#critic-39)] |
| CDAP | A PSM rollup per raw sequence, charge and file (qpx/converters/cdap/feature_adapter.py:116-131). `rt` is the best PSM's `RTAtPrecursorHalfElution`, and its 0.00 sentinel is kept (:321). | [[MZID_CDAP-14](findings-feature.md#mzid_cdap-14)] [[MZID_CDAP-17](findings-psm.md#mzid_cdap-17)] |
| DIA-NN, Spectronaut | One row per peptidoform, charge and run. The window is converted from minutes to seconds (converter-coverage.md:54; qpx/converters/spectronaut/feature_adapter.py:54). | [[INFRA-26](findings-feature.md#infra-26)] |
| mzIdentML | Writes no feature view. | – |

## Why it matters
"Sum the rows of a precursor in a run" gives a different answer for each producer:
- MuData sums 12,672 cells on maxquant_full, mostly direct+transferred pairs [[VALIDATOR-19](findings-feature.md#validator-19)].
- A FragPipe experiment's quantity lands on its first fraction.
- For fractions, mokume takes the MAX and qpx takes the SUM [[MOKUME-10](findings-cross.md#mokume-10)].
- MaxQuant windows and `rt` use different clocks, which breaks XIC extraction and peak widths.
- PD `feature_id` values depend on the PSM the converter picks.

## Options considered
| Option | For | Against |
|---|---|---|
| A. Grain implied by the composite (status quo) | No change | Readers cannot tell parts from alternatives |
| B. `precursor_run` for everyone | Simple pivots | QPX invents sum or pick rules; loses per-peak origin; infeasible for FragPipe |
| C. Peak only (literal [#217](https://github.com/bigbio/qpx/issues/217)) | Matches [#217](https://github.com/bigbio/qpx/issues/217) | PD, CDAP, DIA-NN and TMT rows are not peaks |
| **D. Signal definition, no double counting, grain read from the identity composite** | Every native entity fits; testable | Amends [#217](https://github.com/bigbio/qpx/issues/217) |
| U1. Null `run_file_name`, or unit quantities only in pg | Pure per-run rows | Breaking (non-nullable column), or loses FragPipe precursor quantities |
| **U2. First member run + new `grouped_runs`** | Additive; joins pg | `run_file_name` is a placeholder |
| W1. Footer key for the window's clock, or aligned-window columns | No value change | Violates P-A, or two columns for one producer |
| **W2. Shift by the row's offset `rt − rt_aligned`, else null** | No new column | Derived value; small error where the calibration curve bends inside a peak |

## Recommendation
**FRW-1 (MUST)** A feature row is one quantity of one precursor ion measured in one run, or in one declared unit (FRW-6), or one unidentified MS1 signal marked as such ([T07](t07-unidentified-features.md) UNI-1). This amends [#217](https://github.com/bigbio/qpx/issues/217) ("a feature is a chromatographic peak"): peaks remain the rows of `peak` producers only.
- On an identified row, the ion is one canonical peptidoform at one charge ≥ 1 ([T07](t07-unidentified-features.md) UNI-1). An MS1 multiplet counts as one ion (FRW-10).
- An unidentified row has a null peptidoform and `identification_type = 'unidentified'` ([T07](t07-unidentified-features.md) UNI-3). It has no natural key, so it stands outside the grain of FRW-3 and its checks (FRW-12); its id follows [T07](t07-unidentified-features.md) UNI-4.
- Labels are `intensities` entries ([T09](t09-intensities-missing-values.md) INT-5), never part of the row's key.
- The signal is an MS1 peak or multiplet, or one or more reporter spectra, or none (FRW-5).
- Charge states and modified forms are separate rows. Rows the tool separates by FAIMS CV carry that CV. A quantity spanning several CVs has a null `compensation_voltage`.

This replaces data-model.md:64,74, feature.md:301, producer-guide.md:12 and convert.md:1094.

**FRW-2 (MUST)** A quantified signal (an MS1 peak or a reporter spectrum) gives intensities to at most one row.
- Rows that share `(peptidoform, charge, run_file_name, compensation_voltage)` are therefore distinct signals. A direct and a transferred copy of one peak are forbidden ([T08](t08-match-between-runs.md) MBR-7).
- One signal claimed by several peptidoforms follows [T06](t06-identifications-per-feature.md) FID-10. The reporters of a chimeric spectrum follow [T02](t02-chimeric-spectra.md) CHI-11.
- A producer that cannot see which signals feed a tool's quantity (tool-aggregated evidence, `ion.tsv`) keeps the tool's rows. It lists the gap in converter-coverage ([T02](t02-chimeric-spectra.md) CHI-11) and leaves `identification_conflict` null ([T06](t06-identifications-per-feature.md)).

**FRW-3 (MUST)** No footer key declares the grain. It follows from the feature footer's `identity_composite` (FRW-4), and FRW-2 makes this safe. A file without an `identity_composite` has the schema default `peptidoform, charge, run_file_name, rt`, a `peak` composite. A merge keeps the grain, because it writes one file per composite (qpx/core/dataset_merge.py:311-316).

| Grain | Composite | Rows per natural key | Producers |
|---|---|---|---|
| `precursor_run` | Exactly `peptidoform, charge, run_file_name`, plus `compensation_voltage` when the tool separates CVs | At most one per `(peptidoform, charge, run_file_name, compensation_voltage)` | DIA-NN, Spectronaut, CDAP, PD; quantms label-free and DIA-NN inputs |
| `peak` | Any other | Any number. Each is one MS1 peak or multiplet, or one MS/MS spectrum whose precursor peak the tool did not detect (MaxQuant `MSMS`) | MaxQuant; OpenMS LFQ and MS1 labels |
| `spectrum` | Any other | Any number. Each is one reporter quantification spectrum | quantms and OpenMS isobaric |
| `precursor_unit` | Starts with `quantification_unit_id` | At most one per `(quantification_unit_id, peptidoform, charge, compensation_voltage)` | FragPipe |

The composite does not tell `peak` from `spectrum`; readers treat both alike (FRW-R2).

**FRW-4 (MUST)** The `identity_composite` depends on the grain:

| Grain | Composite |
|---|---|
| `precursor_run` | Exactly `peptidoform, charge, run_file_name`, plus `compensation_voltage` when the tool separates CVs |
| `precursor_unit` | `quantification_unit_id, peptidoform, charge, compensation_voltage` |
| `peak`, `spectrum` | Those three columns plus the tool-reported locators of the signal: MaxQuant `rt_start, rt_stop, scan`; OpenMS `rt, observed_mz, consensus_rt`; quantms `rt, scan` |

- The composite never holds a value chosen by the converter's best-PSM rule, such as PD's `rt` today or OpenMS `scan` under [T06](t06-identifications-per-feature.md) FID-5.
- Rows are grouped on the canonical peptidoform.
- `compensation_voltage` is null without FAIMS, never 0 ([T11](t11-spectrum-precursor-fields.md) SPC-12).

**FRW-5 (MUST)** An identification-only row has no non-null intensity entry ([T09](t09-intensities-missing-values.md) INT-6).
- Under a `precursor_*` grain, it is the only row of its key.
- Under `peak` or `spectrum`, a producer MAY write one per identifying MS/MS spectrum (MaxQuant `MSMS`).
- "Identified, not quantified" is decided per key and label, over all rows of the key.

**FRW-6 (MUST)** Under `precursor_unit`, every row:
- carries the new column `grouped_runs` (`list<string>`), with the meaning and order rule of `pg.grouped_runs` ([T12](t12-protein-groups.md) PGV-5); a one-run unit has one element;
- has `run_file_name = grouped_runs[0]`, a representative member rather than the measuring run;
- has a `quantification_unit_id` that maps 1:1 to its `grouped_runs` set within the file.

When a unit has several runs, its rows have null `rt`, `rt_aligned`, `rt_start` and `rt_stop`, and name no identifying PSM: `scan = []`, null `id_run_file_name`, PEP and `psm_qvalue` ([T06](t06-identifications-per-feature.md) FID-1, FID-6). A one-run unit row follows the per-run rules. Other grains leave `grouped_runs` null. Producers SHOULD write per-run rows whenever the tool reports them.

**FRW-7 (MUST)** Units and clocks:
- `rt`, `rt_start`, `rt_stop`, `rt_aligned`, `consensus_rt` and `predicted_rt` are in seconds.
- `rt`, `rt_start` and `rt_stop` use the run's acquisition clock: 60 × `mz.scan_start_time`, which is also the clock of `xic.rt`.
- `rt_aligned` and `consensus_rt` use the tool's aligned clock.
- A time the tool reports only on its aligned clock is moved to the run clock by the row's offset `rt − rt_aligned` when both are known. Otherwise it is null.
- OpenMS: ProteomicsLFQ transforms feature RTs, the RTs of their identifications and of unassigned identifications in place when it aligns, and stores no `original_RT`. A consensusXML it wrote with alignment therefore holds only aligned times: they go to `rt_aligned`, and `rt` is null. Only the native writer can supply acquisition-clock times, by carrying the original RT per feature handle.

**FRW-8 (MUST)** `rt`:
- When the tool reports a time for the precursor's MS1 peak in this run, `rt` is that time (apex, centroid or half-elution). On a transferred row this is the matched peak's time, never the donor's.
- Otherwise it is the time of the identifying spectrum ([T06](t06-identifications-per-feature.md) FID-3), equal to its `psm.rt`. This covers summed and one-spectrum reporter rows.
- A time known only on the aligned clock never fills `rt` (FRW-7): in consensusXML written by ProteomicsLFQ with alignment, `rt` and `psm.rt` ([T11](t11-spectrum-precursor-fields.md) SPC-11) are null.
- A 0.0 sentinel is null.

**FRW-9 (MUST)** `rt_start` and `rt_stop` bound the integrated MS1 signal. They are null when nothing was integrated or no bounds are reported, and they are never placeholders (MaxQuant `MSMS`). When all three are set, `rt_start ≤ rt ≤ rt_stop`. Every time is ≥ 0. OpenMS: ProteomicsLFQ clears the feature hulls, so the consensusXML route cannot recover the bounds and writes null; the native writer currently writes `rt ± width/2`, not the hull bounds.

**FRW-10 (MUST)** MS1-label multiplets:
- There is one row per multiplet per run (intensities.md:84). This is the encoding agreed in [#308](https://github.com/bigbio/qpx/issues/308): the light peptidoform on the feature row, one intensity entry per channel, and PSMs keeping their labelled peptidoform.
- `peptidoform` and `calculated_mz` belong to the partner of the run's lightest channel, with any dimethyl or mTRAQ label inline.
- `observed_mz` and `mass_error_ppm` describe that same partner, and are null when it was not measured.
- `rt` is the lightest measured partner's time (partners co-elute).
- `channel_details` SHOULD list every measured partner.
- No partner also appears as its own row.

**FRW-11 (MUST)** A row that sums several quantification spectra follows [T09](t09-intensities-missing-values.md) INT-10: its spectrum count is in `quant_spectrum_count`, and `parent_ion_fraction` and `quant_scan` are null.

**FRW-12 (MUST)** The validator reports the following:

| Check | Severity |
|---|---|
| A repeated natural key under `precursor_run` or `precursor_unit`, among identified rows | Error under `dda`, else warning ([T13](t13-metadata-conformance.md) META-8) |
| A `precursor_unit` row that breaks FRW-6 (null `grouped_runs`, `run_file_name ≠ grouped_runs[0]`, unit id not 1:1) | Error |
| `rt` outside `[rt_start, rt_stop]`, or a negative time | Error under `dda`, else warning ([T13](t13-metadata-conformance.md) META-8) |
| Two intensity-bearing rows equal in natural key, `rt` and `intensities` | Error |
| A composite that breaks FRW-4 | Warning |

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| Grain derived from `identity_composite` (FRW-3, FRW-4); no new footer key | feature footer | – | – | No |
| `grouped_runs` | feature | `list<string>`, nullable, optional | Under `precursor_unit` | No |
| Text for `rt`, the window and `run_file_name` (feature.md:65-70,89; feature.yaml:255-260) | feature | doc | – | **Meaning change** for 1.3 MaxQuant windows (calibrated). Applies from `qpx_version` 1.4; `qpxc migrate` shifts 1.3 windows by `rt − rt_aligned`, both in the file |
| `compensation_voltage` null without FAIMS (feature.md:76; feature.yaml:128-132) | feature | doc | – | **Yes**: FragPipe values and `feature_id` |
| Composites and locator values | feature footer | – | – | **`feature_id` changes**: PD and quantms label-free (composite), MaxQuant (window values), OpenMS (`observed_mz`; `rt` null when the file was aligned; `scan` dropped) |

## Producer changes
| Producer | Change needed |
|---|---|
| MaxQuant | Grain `peak` (current composite). Set `rt_start` to calibrated start + (`rt` − `rt_aligned`), and `rt_stop` likewise (feature_adapter.py:574-585). Null the window on `MSMS` rows. Shared isotope patterns per [T06](t06-identifications-per-feature.md) FID-10. Fix tool-mappings.md:50-51 |
| FragPipe | Grain `precursor_unit` (current composite). `grouped_runs` on every row, ordered per [T12](t12-protein-groups.md) PGV-5, never re-sorted (pg_adapter.py:96). Null CV without FAIMS (feature_adapter.py:519). Apex RT of one-run units, if `combined_ion` has one (unverified) |
| Proteome Discoverer | Grain `precursor_run`; composite `peptidoform, charge, run_file_name`; group on the canonical peptidoform (feature_adapter.py:16-20). Summed TMT rows: null PIF, `quant_spectrum_count` ([T09](t09-intensities-missing-values.md)) |
| OpenMS consensusXML | Grain `peak`, or `spectrum` for isobaric data. Lightest dimethyl or mTRAQ label inline (:638-639), and fix convert.md:1038-1040. `observed_mz` is the run element's m/z, else null (:645); the ProteomicsLFQ element m/z is theoretical, so it is null there ([T11](t11-spectrum-precursor-fields.md) SPC-1). Composite per FRW-4. When the file declares alignment, element RTs are aligned times: they go to `rt_aligned`, and `rt` is null (FRW-7, FRW-8). `rt_start` and `rt_stop` are null, because the hulls are cleared (FRW-9) |
| OpenMS `-out_qpx` | Upstream declares a `peak` or `spectrum` composite per FRW-4; `qpxc convert openms` already copies `identity_composite` (qpx/converters/openms/converter.py:182-220). Verified at HEAD: no `rt_aligned` or `consensus_rt` column (the unmerged [#378](https://github.com/bigbio/qpx/issues/378) branch adds `consensus_rt`); with alignment `rt` is the aligned time; `rt_start`/`rt_stop` are `rt ± width/2`. Needed: carry the original RT per feature handle, so that `rt` is on the acquisition clock and the aligned time goes to `rt_aligned` (FRW-7), and keep the hull bounds (FRW-9). See the [OpenMS worklist](openms-worklist.md). |
| quantms MSstats | `spectrum` for isobaric inputs. `precursor_run` with composite `peptidoform, charge, run_file_name` for label-free and DIA-NN inputs |
| CDAP | Grain `precursor_run` (current composite). `rt` 0.00 becomes null (feature_adapter.py:321) |
| DIA-NN, Spectronaut | None: the current composite gives `precursor_run`. mzIdentML: none |
| Writers | None: the grain follows from the `identity_composite` writers already stamp (FRW-3) |

## Reader rules
**FRW-R1 (MUST)** Derive the grain from `identity_composite` (FRW-3), and read every other composite as `peak`. Treat `(peptidoform, charge, run_file_name)` as unique only under `precursor_run`.

**FRW-R2 (MUST)** Under `peak` or `spectrum`, the rows of one key are distinct signals. A reader that collapses them per key and label:
1. first applies [T08](t08-match-between-runs.md) MBR-R4 (drop a transferred entry beside a measured `direct` or `secondary` entry of the same label) and, to reporter rows, [T02](t02-chimeric-spectra.md) CHI-R6;
2. sums the non-null entries, as MaxQuant peptides.txt sums "all isotopic clusters" [[EXTERNAL-17](findings-feature.md#external-17)];
3. reports how many cells it collapsed.

Under `precursor_run`, duplicates are an error. Current violations:
- qpx/core/dataset_queries.py:124-137 sums without a filter or a count.
- `build_mudata` cannot exclude transferred rows (qpx/mudata.py:179-190,685-691).
- mokume keeps the MAX, as if the rows were redundant (mokume:python/mokume/preprocessing/aggregation.py:295-322; mokume:rust/crates/mokume-pipeline/src/lib.rs:2761-2772).

**FRW-R3 (MUST)** Rows from the fractions of one unit are disjoint parts of one sample, so readers add them. mokume takes the MAX (mokume:python/mokume/preprocessing/aggregation.py:263-292; mokume:python/mokume/pipeline/stages.py:128-140; mokume:rust/crates/mokume-pipeline/src/lib.rs:6205-6211).

**FRW-R4 (MUST)** A row whose `grouped_runs` holds several runs belongs to its unit, not to `run_file_name`. Current violations:
- MuData precursor `obs` uses `run_file_name` (qpx/mudata.py:179-190).
- The run summary counts rows per `run_file_name` (qpx/views/api.py:21-28).
- mokume never reads `quantification_unit_id` (mokume:python/mokume/io/feature.py:679-686).

**FRW-R5 (MUST)** A label is a channel, never a run name. mokume's "label = run" layout (mokume:python/mokume/io/feature.py:605-643; mokume:rust/crates/mokume-io/src/qpx.rs:528-539) is allowed only behind a legacy flag.

**FRW-R6 (MUST)** Quantified-feature counts include only identified rows with a non-null entry; unidentified rows are excluded ([T07](t07-unidentified-features.md) UNI-R1). Identification counts follow [T08](t08-match-between-runs.md) MBR-R3. A missing feature row means "not identified" only when provenance states `unquantified_rows = retained` ([T07](t07-unidentified-features.md) UNI-R4). qpx/views/api.py:25 counts MaxQuant `MSMS` rows as features.

**FRW-R7 (MUST)** Compare times across runs only through `rt_aligned`, and join spectra on 60 × `mz.scan_start_time`. Link multiplet PSMs through `psm.feature_id`. Read the window clock from FRW-7 only in 1.4 or migrated files.

## Open questions for maintainers
1. *Resolved from the OpenMS source (default branch at 0707508):* ProteomicsLFQ transforms feature RTs, the RTs of their identifications and of unassigned identifications in place when it aligns, and does not store `original_RT`. So in consensusXML written by ProteomicsLFQ with alignment, element and identification RTs are on the aligned clock [[OPENMS-M3](findings-feature.md#openms-m3)]: the Python converter writes them to `rt_aligned` and leaves `rt` and `psm.rt` ([T11](t11-spectrum-precursor-fields.md)) null. Only the native writer can supply acquisition-clock times, by carrying the original RT per feature handle (FRW-7, FRW-8).
2. Should MaxQuant `MSMS` rows collapse to one per key? That would drop 8,731 rows that sit beside a quantified row, and it interacts with [T07](t07-unidentified-features.md)'s provenance parameter `unquantified_rows`.
3. Should `quantification_unit_id` be canonical (runs joined by `;`, as views.md:212 names units) instead of the tool's experiment name? FragPipe ids would change.
4. Does FragPipe `combined_ion` report a per-experiment apex RT [[FRAGPIPE-17](findings-feature.md#fragpipe-17)]?
