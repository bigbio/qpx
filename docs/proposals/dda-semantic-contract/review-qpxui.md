# Review: qpxui as a consumer

[Overview](index.md) · [Strict profile](strict-profile.md) · [Consumers](consumers.md) · [OpenMS work list](openms-worklist.md)

!!! note "Review report"
    One of three reviews of the [DDA semantic contract proposal](index.md), made by reading bigbio/qpxui at `ddeda63`. It is kept as written by the reviewer. It predates the restructuring around the strict profile, so some of its recommendations are already applied. Rule ids refer to the topic sections T01 to T14.

Nothing was edited or pushed. qpxui citations are relative to the qpxui repository root (HEAD `ddeda63`). Rule ids come from the topic sections T01 to T14 and the reduced-footer revision.

## Bottom line

- **qpxui reads columns and almost nothing else.** Its only footer key is `qpx.sort_order`. It never reads `qpx_version`, `scan_format`, `rank`, `pass_threshold`, `identification_type`, `imputed_labels`, `pg_group_id`, `grouped_runs`, provenance, pepmap or SDRF, and it applies no q-value threshold anywhere. Its one row filter is `NOT coalesce(is_decoy, false)`.
- **Its largest problem today is outside the proposal.** qpxui reads the pg view in the QPX 1.0 layout (`run_file_name` plus an `intensities` list). On a pg file of 1.1 or later (`grouped_runs`, `label`, scalar `intensity`) it silently loses the pg quantities, per-run protein rows and per-run group counts.
- **Nothing in the proposal hard-breaks qpxui.** It reads no ids, looks up no score or intensity names, and guards every optional column. Several changes do silently empty a chart or shift a number.
- **The proposal fixes most of qpxui's wrong counts, but only if qpxui adopts about six predicates.** Some proposed defaults need provenance, which a SQL-over-range-requests UI will not read.
- **qpxui was not among the audited readers.** The proposal's reader audit covers the qpx API and mokume only.

## 1. How qpxui consumes QPX

**Views.** Reports use `feature`, `psm`, `pg`, `run`, `sample`, `dataset` (`src/lib/qpx.ts:13-17`). It takes the first parquet file per view name (`src/lib/report.ts:37`). Any other view is only browsable as a table. π-MSNet `*-MSNet.parquet` files are mapped to `psm` (`src-tauri/src/qpx/mod.rs:46-77`).

**Columns.** The full list is in `src/lib/qpx.ts:30-95`. Struct fields read:
- `intensities[].label` / `.intensity`
- `modifications[].name` / `.accession` / `.positions[].position` / `.amino_acid`
- `run.samples[].sample_accession` / `.label`
- `pg.peptide_counts.total_sequences` / `.unique_sequences`, `pg.feature_counts.total_features`
- `pg_accessions[].accession` (struct or string form, `src/lib/report.ts:226,233`)
- `cv_params` in three encodings, with name aliases (`src/lib/qpx.ts:115-142`); names used are `dissociation method`, `normalized collision energy`, `mass analyzer` (`src/components/SpectrumViewer.vue:129,134,226`)

**Score names.** None. `additional_scores` is only rendered generically (`src/lib/format.ts:33-39`, `src/components/StructuredView.vue:80-85`). `additional_intensities` is never read.

| # | Assumption | Where |
|---|---|---|
| A1 | Target row = `NOT coalesce(is_decoy,false)`; a null `is_decoy` is a target | `src/lib/report.ts:213-216` |
| A2 | "PSMs" and "Features" KPIs are `count(*)` of all rows, decoys included; from the parquet footer on large views | `src/lib/queries/header.ts:34-35,44-45`; `src/components/ReportView.vue:146-152` |
| A3 | Every target row is an identification: no rank, threshold or MBR filter. Peptides = distinct `sequence` over the evidence view, which is feature if present, else psm | `src/lib/report.ts:104-106`; `src/lib/queries/header.ts:41-43`; `src/lib/queries/qc.ts:25-36`; `src/lib/queries/overview.ts:56-110`; `src/lib/queries/ptms.ts:20-60` |
| A4 | Protein = `anchor_protein`; else `protein_accessions[1]`; else `pg_accessions[1]`. Groups are keyed, counted and matrixed by anchor | `src/lib/report.ts:222-228`; `src/lib/queries/header.ts:36-39`; `src/lib/queries/proteins.ts:17-31`; `src/lib/queries/quant.ts:23` |
| A5 | pg has scalar `run_file_name` and an `intensities` list | `src/lib/queries/quant.ts:16-20`; `src/lib/queries/proteins.ts:21,114-129`; `src/lib/queries/qc.ts:46-50`; `src/components/report/ReportProteins.vue:125` |
| A6 | Feature rows are disjoint signals: all rows of a (protein or peptidoform, sample) are summed, across charges, fractions, replicates, direct and transferred rows | `src/lib/queries/quant.ts:32-42`; `src/components/report/PeptidoformProfile.vue:121-129` |
| A7 | Intensity ≤ 0 or null is missing; log scale only | `src/lib/queries/quant.ts:41`; `src/lib/queries/profile.ts:47`; `src/lib/queries/qc.ts:76`; `src/lib/report.ts:285-287` |
| A8 | Sample = `run.samples` joined on (run, label). Labels are upper-cased, spaces stripped, `TMT131N`→`TMT131`; a single-sample run takes any label; duplicates resolve by `min(sample_accession)`. The QC tab uses exact equality instead | `src/lib/report.ts:334-348`; `src/lib/queries/qc.ts:78-82` |
| A9 | Contaminant = `pg.contaminant`, with null read as false, or the accession containing `CONTAMINANT` or `CON__` | `src/lib/queries/quant.ts:26-28` |
| A10 | `unique` is shown per peptidoform as `bool_or(unique)` | `src/lib/queries/proteins.ts:80`; `src/components/report/ReportProteins.vue:96` |
| A11 | One PSM = (peptidoform, run, first scan element, charge), `LIMIT 1`; list scans matched with `list_contains`; USI always `…:scan:N:…` | `src/lib/queries/spectrum.ts:38-68`; `src/lib/queries/peptides.ts:79-82`; `src/components/SpectrumViewer.vue:74-78` |
| A12 | `rt`, `rt_start`, `rt_stop` are seconds on one clock | `src/lib/queries/profile.ts:33,38`; `src/lib/queries/overview.ts:103`; `src/lib/peptidoformCharts.ts:150-185` |
| A13 | `peptide_qvalue` is a peptide-level q-value; PEP is comparable across rows (`min` = "Best PEP") | `src/lib/queries/overview.ts:105-110`; `src/lib/queries/summary.ts:52-55`; `src/lib/queries/peptides.ts:61` |
| A14 | `sequence` equals the stripped peptidoform, with a sampling workaround when it does not | `src/lib/report.ts:61-77` |
| A15 | Feature and pg are sorted by `anchor_protein` when relayouted | `src/lib/report.ts:175-180`; `src/components/report/ReportPeptides.vue:30-35` |

Fractions are only displayed (`src/lib/queries/qc.ts:55`) and are summed implicitly per sample.

## 2. Do QPX 1.3 producers satisfy these today?

Producer figures are from the proposal's `findings-*.md`; I did not re-measure them.

| # | Status | Where qpxui is wrong or misleading |
|---|---|---|
| A5 | **Broken for every pg file ≥ 1.1** | pg schema at `qpx:qpx/core/data/schemas/pg.yaml:46-72`. The Quantification tab falls back to feature sums per anchor; MaxLFQ and tool protein quantities are never shown. No per-run pg table, no `n_pg` in QC, no pg intensity in the protein list. It works on the portal's published 1.0-layout files, per `quantms-portal-backend:src/quantms_portal_backend/apps/coexpression/README.md:76`. |
| A3 | Violated | **Rank:** Mascot mzIdentML has 13,386 rows for 9,707 spectra (MZID_CDAP-01); PD and OpenMS write every hit (PD-01, OPENMS-01). The PSM KPI, peptide counts and PTM tab are inflated; the same defect in qpx's `ModificationView` gives 4,323 vs 926 deamidation sites. **MBR:** MaxQuant averages 9,021 features per run against 5,156 non-transferred, and 6,128 vs 3,779 peptides (VALIDATOR-20), so per-run "Peptides" and "Proteins" in QC and the peptidoform "Runs" column are inflated. **Failing rows** (mzIdentML, PD `pass_threshold=false`) and unfiltered OpenMS hits (OPENMS-03) count as identifications. |
| A2 | Misleading | The KPIs include decoys (MaxQuant, OpenMS, mzIdentML retain them), rank > 1 rows, transfers and 12,758 MaxQuant `MSMS` identification-only rows (MAXQUANT-09). |
| A6 | Violated | MaxQuant: 9,272 keys hold both a measured direct and a measured transferred row, and both are summed. 312 peaks are credited in full to two or three peptidoforms (MAXQUANT-08). PD sums a spectrum twice with two engines (PD-04). FragPipe puts a whole experiment on its first run (FRAGPIPE-10), so per-run charts are wrong. Technical replicates of a sample are summed too. |
| A4 | Violated | Distinct groups sharing an anchor are merged (CRITIC-66). Shared peptides have a null anchor for OpenMS (OPENMS-27) and drop out of the Proteins and Quantification tabs. `protein_accessions[1]` is order-dependent and can be a decoy entry on a target row (CRITIC-19). MaxQuant without proteinGroups.txt and PD put several groups' leaders in `pg_accessions` (MAXQUANT-28, PD-20). |
| A8 | Mostly satisfied by 1.3 converters | The fallbacks exist for the older published data, where TMT labels are channel indices (same README, line 78). They hide real errors: MaxQuant TMT without SDRF is written as `LFQ` (MAXQUANT-15), PD writes entries for channels without a sample (PD-M5), and PD `.RAW` runs fail the run join (PD-24). |
| A9 | Violated | Feature and psm have no `contaminant` column. The name test misses `CONTAM_` and `contam_`. FragPipe strips the marker (MOKUME-M2); PD and CDAP flags are constant false (PD-23, MZID_CDAP-M3). |
| A1 | Holds for conformant files | mzIdentML flags target+decoy peptides as decoy (MZID_CDAP-04), so qpxui drops some targets. PD `is_decoy` is always false. |
| A10 | Producer-dependent | MaxQuant's meaning depends on whether proteinGroups.txt was given (MOKUME-02); PD and FragPipe use the mapped-protein count (PD-21, FRAGPIPE-18). |
| A11 | Violated | K duplicates exist (MAXQUANT-03, MZID_CDAP-08), so `LIMIT 1` is arbitrary. `scan[0]` with `list_contains` is wrong for nativeId arrays such as `[frame, scan]`. Index or hashed scans produce wrong USIs (OPENMS-17, MZID_CDAP-19). No producer besides OpenMS and mzIdentML declares `scan_format`. |
| A12 | Violated | MaxQuant: on 83,413 rows `rt` lies outside the window (MAXQUANT-10), so the retention chart draws windows that miss their point. CDAP writes `rt=0` sentinels (MZID_CDAP-17). |
| A13 | Violated | OpenMS `peptide_qvalue` is PSM-level with 1.0 placeholders (INFRA-03). MaxQuant, FragPipe and PD leave it null, so the chart is empty. PEP is 1 − probability for FragPipe and a p-value for Mascot. |
| A14 | Workaround needed | Also `'None'` peptidoforms (MAXQUANT-M1) and `''` or charge-0 rows are counted. |
| A7 | Satisfied in practice | Converters write sentinel 0 as null, and quantms-MSstats zeros are missing values. |

## 3. The proposal from qpxui's point of view

### (a) Rules that directly fix a qpxui problem

- **Counts:** PSM-1/2 (unique K), PSM-3 (`rank`), FLT-7/9 (`pass_threshold`), FLT-1/2 (`fdr_filter`), CHI-R1 (spectra = distinct run and scan), RDR-9.
- **MBR:** MBR-1/2/3 (`identification_type` reliable), MBR-R3/R4 with RDR-5, MBR-10 and PGV-1 (`identified_in_unit`).
- **Rows:** UNI-1 and SPC-6/15 (no `''`, `'None'` or charge-0 rows; `sequence` equals the stripped peptidoform, so the A14 workaround can go), FRW-2 and FID-10 (a signal feeds one row), FRW-6 (`grouped_runs` on unit rows), INT-5/6/7.
- **Proteins:** PRT-2 (no decoy entries on target rows), PRT-6 (`contaminant` on psm and feature; `src/lib/queries/quant.ts:28` picks it up with no code change), PRT-7/8/13 (group membership and `pg_group_id`), PRT-9 (`group_unique`), PGV-7/8.
- **Spectra and time:** FID-5, SPC-7/8/9, SPC-11, SPC-13, FRW-7/9.
- **Labels and units:** INT-2 (`intensity_quantity`, `intensity_processing`). qpxui median-normalises by default (`src/components/report/ReportQuant.vue:70`) without knowing whether the data is already MaxLFQ or normalized.
- **Confidence:** CNF-1/2 (`psm_qvalue`), CNF-7/8.

### (b) Changes that break qpxui or force a change

I found no hard break: no SQL errors, no id reads, no score-name lookups. CNF-6 (scores leaving `additional_scores`), `lfq`→`maxlfq`, the id changes, the `match_qvalue` retype and the new `unidentified_feature` view (parsed as its own view name) are harmless.

Silent effects and required edits:

| Change | Effect on qpxui | Site |
|---|---|---|
| CNF-2: OpenMS PSM-level value leaves `peptide_qvalue` | The "Peptide q-value" chart and the summary go empty for quantms data; needs a `psm_qvalue` or `precursor_qvalue` option | `src/lib/queries/overview.ts:108-110`; `src/lib/queries/summary.ts:52-55` |
| CNF-8, MBR-3: PEP null on transfers and OpenMS peptide-level PEP | "Best PEP" and the PEP histogram get sparser | `src/lib/queries/peptides.ts:61`; `src/lib/queries/overview.ts:105` |
| INT-4: QPX roll-ups leave the pg primary (OpenMS, CDAP → `intensity_quantity=none`) | No primary protein quantity for the commonest producers; needs a by-name lookup in `additional_intensities`, which has three encodings | `src/lib/queries/quant.ts:16-20` |
| PRT-8: unassigned shared peptides get a null anchor (PD no longer `masters[0]`) | More peptides vanish from the Proteins and Quantification tabs | `src/lib/report.ts:222-246` |
| PRT-1: `protein_accessions` sorted by code point | The "protein" of psm-only datasets changes | `src/lib/report.ts:225` |
| PRT-10: `unique` deprecated, removed in 2.0 | The column disappears; switch to `group_unique` | `src/lib/queries/proteins.ts:80` |
| INT-5, UNI-6: identification-only rows added | "Features" KPIs and per-run counts grow; split identified from quantified | `src/lib/queries/header.ts:44`; `src/lib/queries/qc.ts:27` |
| SPC-9: empty `scan` with `spectrum_reference` | The PSM key loses `scan`, so the wrong spectrum is shown; add `spectrum_reference` to `PsmKey` | `src/lib/queries/spectrum.ts:8-13,54-60` |
| SPC-7: `scan_format` required | The USI must use `scan`, `index` or `nativeId`; the full scan array must be matched, not its first element | `src/components/SpectrumViewer.vue:74-78`; `src/lib/queries/spectrum.ts:45,51-52` |
| PGV-8: `pg_qvalue` null for MaxQuant, PD, FragPipe | The per-unit q-value column is blank, because it prefers `pg_qvalue` | `src/lib/queries/proteins.ts:122` |
| RDR-8: no label folding or lone-sample fallback | `sampleMapSql` violates it; gate the fallback on `qpx_version` < 1.4 | `src/lib/report.ts:337-346` |

Other reader rules qpxui violates today:
- RDR-4 / PRT-R2: null `is_decoy` read as target (`src/lib/report.ts:215`).
- RDR-5: direct and transferred summed (`src/lib/queries/quant.ts:38`).
- RDR-7 / PRT-R1: anchor as key.
- FLT-R6 / PRT-R3: contaminants by substring.
- PGV-R1: `n_pg` counts every pg row (`src/lib/queries/qc.ts:47`).
- INT-R3: imputed values not masked.
- UNI-R2: unidentified rows not excluded.
- RDR-11: `qpx_version` never read.

### (c) Irrelevant or over-engineered for this consumer

- `precursor_index` and CHI-3 to CHI-11. qpxui only needs "spectra = distinct (run, scan)".
- FID-3 total order, `consensus_feature_id`, `identification_conflict`, `quant_spectrum_count`, the `match_*` definitions, `isotope_error` (qpxui already clips the ppm histogram, `src/lib/report.ts:289-297`), ion mobility.
- The `unidentified_feature` view.
- Every provenance parameter: `pass_threshold_level`/`_fdr`, `*_qvalue_scope`, `mbr_fdr`, `decoy_policy`, `contaminant_policy`, `observed_mz_calibration`, `skipped_runs_psm`.
- Validation profiles and severities (META-7/8), `qpx_read_policy` (RDR-13; a footnote is enough), and "raise when every row is unassessable" (RDR-3; a UI shows a count).
- META-3/4 matter only indirectly: qpxui consumes relayouted and merged portal files, so keys must survive those rewrites.

### (d) What qpxui needs that the proposal does not give

1. **Defaults decidable from the file alone.** Three reader rules need provenance:
   - MBR-R2: is a null type own-run?
   - PGV-R5: may groups be pivoted across units? This is the Quantification tab's default behaviour, so by the revision's own test `pg_grouping` qualifies as a footer key; otherwise PGV-R5 should become "pivot and warn".
   - FLT-R4: the basis of `pass_threshold`.
2. **A feature-level acceptance flag.** `pass_threshold` exists on psm and pg only. On an unfiltered file (`fdr_filter=none`, OpenMS) a UI must guess which q-value column is filled.
3. **Pre-computed counts.** qpxui defers distinct counts on large streamed views (`src/components/ReportView.vue:146-152`). The portal computes its own with different rules (`quantms-portal-backend:src/quantms_portal_backend/qpx_stats.py:110-120`: group-keyed, contaminant groups excluded, peptides from psm). The same dataset therefore shows different numbers in three places.
4. **A displayable protein quantity.** After INT-4 there is no rule for which `additional_intensities` name to show when the primary is null.
5. **Sample on the row.** Resolving pg `(any of grouped_runs, label)` → sample, with the PGV-R4 several-samples check, is awkward SQL. An optional writer-resolved `sample_accession` on pg rows would make the matrix one `SELECT`.
6. **One encoding for name-value lists.** `qpx.struct_lists` leaves three layouts for `cv_params`, `additional_scores` and `additional_intensities`. CNF-6 helps by moving q-values and rank to typed columns.
7. **A cv_params name registry for spectrum metadata** (dissociation method, NCE, analyzer); qpxui carries aliases (`src/lib/qpx.ts:115-118`).
8. **A pg-layout statement in the reader contract.** The 1.0 versus ≥ 1.1 pg difference (A5) is the largest practical trap, and RDR-11 does not mention it.
9. **Reference SQL.** RDR-14 gives fixtures and expected numbers but no queries.

## 4. Minimal consumer-first subset, in priority order

1. **Decoys:** `is_decoy` non-null with the all-entries rule; no decoy entries on target rows (PRT-2, PRT-4/5, PRT-R2).
2. **PSM rows:** unique K and a typed `rank` (PSM-1/2/3). Count `rank = 1`.
3. **Acceptance as a column:** FLT-7, with FLT-9 promoted from SHOULD to MUST for qpx converters, so accepted = `NOT is_decoy AND rank = 1 AND pass_threshold`. `fdr_filter` (FLT-1/2) then only drives a "filtered at … / unfiltered, upper bound" badge.
4. **`identification_type`:** filled whenever known (MBR-1/2/3). Identified = `IN ('direct','secondary')`. Apply MBR-R4 before summing. Null = "origin unknown", shown as a count with no provenance lookup.
5. **Feature rows are identified precursors:** UNI-1, SPC-6/15. Missing values per INT-5/6/7.
6. **Samples:** exact `(run, label)` → `run.samples`, one run-name function, `grouped_runs` on unit rows (RDR-8, SPC-13, FRW-6, PGV-5).
7. **Proteins:** group on `pg_group_id`, display `anchor_protein`, typed `contaminant` in all three views (PRT-6/7/11/13).
8. **Quantity declaration:** `intensity_quantity` and `intensity_processing` (INT-2), plus a fixed display fallback when the pg primary is null.
9. **Confidence levels by column name:** CNF-1/2 with `psm_qvalue`, CNF-7/8.
10. **Spectra and time:** `scan_format` required, `spectrum_reference`, one spectrum in `feature.scan`, seconds on the acquisition clock (SPC-7/9/11, FID-5, FRW-7/9).
11. **pg evidence:** `identified_in_unit` and fixed count and q-value scopes (MBR-10, PGV-1/7/8).
12. **Version gating:** `qpx_version` and the key-value footer (META-2/3, RDR-11), including the pg layout.

Simplifications I recommend so a UI can do this in SQL:
- Add reference DuckDB queries to `consumer-guide.md` for the RDR-9 counts and the RDR-5/6 matrix.
- Have the writer stamp the RDR-9 counts in the footer or `dataset.parquet`.
- Allow QPX sums in the pg primary under `intensity_processing = sum`, or fix a fallback order.
- Keep every default free of provenance: FLT-9 as MUST, and PGV-R5 as a warning or `pg_grouping` back in the footer.

## Not verified

- I ran no qpxui code and no queries; DuckDB is not installed here. All behaviour above is from reading the source.
- The published portal data layout (1.0 pg, TMT index labels) is taken from the portal backend README, not from files. I do not know which producer writes the portal's quantms data; the OpenMS `-out_qpx` writer is unaudited in the proposal too.
- π-MSNet psm files are a producer the proposal does not cover.
- Whether MaxQuant transferred rows carry a donor PEP, which would affect "Best PEP", was not checked.
- The commit `15afbca1` cited in `src/lib/qpx.ts:3` is not in the local qpx clones.
- I read the rule lists of T01–T13 and T14 in full, but not every producer table.
