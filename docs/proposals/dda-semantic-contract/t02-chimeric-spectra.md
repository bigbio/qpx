# T02. Chimeric and co-isolated spectra

[← T01](t01-psm-rows.md) · [Overview](index.md) · [Strict profile](strict-profile.md) · [Annex overview](annex-overview.md) · [T03 →](t03-filtering-state.md)

!!! note "Annex: producer-neutral analysis"
    Part of the annex of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text, and they cover every producer. For OpenMS-written files the [strict profile](strict-profile.md) is the contract and wins where the two differ. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** One peptidoform per row, no ProForma `+`; a typed `psm.precursor_index` (1 = primary, ≥ 2 = co-isolated) with rank 1 per precursor; feature `secondary` = identified here only by secondary PSMs; a chimeric spectrum's reporters feed at most one feature (a plain MUST); no footer keys: a null index means unknown, and removed secondaries are recorded in provenance.

## Questions a reader has
- Two PSM rows share `run_file_name` and `scan`: co-isolated peptides (both may be right), alternatives or ties (at most one right), or engine duplicates?
- Does `rank = 1` drop co-identifications? Do `observed_mz` and `charge` describe the selected precursor or the co-isolated ion?
- What does feature `identification_type = 'secondary'` mean, and who writes it?
- May `peptidoform` hold `PEPTIDE/2+SEQUENCE/3`?
- In TMT data, may one spectrum's reporters feed several features?
- How do I count spectra, PSMs and spectral counts?

## What happens today
In the spec:
- "each spectrum yields one or more peptide identifications" (docs/spec/psm.md:25); "One row per match" (docs/spec/producer-guide.md:11). The composite allows one row per peptidoform and spectrum (psm.md:12-15; serialization.md:304-305).
- `rank` is a recommended score with no scope (psm.md:291; scores.md:58). Nothing covers chimeric PSMs.
- `secondary` exists only on features (feature.md:145; qpx/core/data/schemas/feature.yaml:210-214). Only its vocabulary is checked (qpx/core/data/validation.py:71), and only MaxQuant writes it (converter-coverage.md:199).
- `PEPTIDE/2+SEQUENCE/3` is "passed through unchanged" (scope.md:36,41).

| Producer | Behaviour | Findings |
|---|---|---|
| MaxQuant | `MULTI-SECPEP` PSMs are marked only by cv_param `type` (qpx/converters/maxquant/psm_adapter.py:348-353). A first and second peptide with equal peptidoform and charge duplicate `psm_id` (reproduced; 9 on PXD003133, [#378](https://github.com/bigbio/qpx/issues/378)); maxquant_full evidence has 35 such pairs. Ids 167037/167038 (scan 10817) share m/z and intensity: one MS1 cluster reported twice. Only exact `MULTI-SECPEP` maps to feature `secondary` (qpx/converters/maxquant/evidence_fields.py:10-17). Dedup ranks only `MULTI-MATCH` last, then by PEP (qpx/converters/maxquant/feature_adapter.py:230-295), so a lower-PEP SECPEP twin would make the row `secondary`. maxquant_full: 5,110 of 243,822 rows are SECPEP; 4,525 of 5,066 SECPEP scans are named by no other evidence row (evidence names only its best MS/MS); 44 scans have two SECPEP rows | [[MAXQUANT-03](findings-psm.md#maxquant-03)] [[MAXQUANT-29](findings-feature.md#maxquant-29)] [[EXTERNAL-05](findings-psm.md#external-05)] [[SPEC-18](findings-psm.md#spec-18)] |
| FragPipe | One PSM per psm.tsv row; no rank or co-isolation read or derived (qpx/converters/fragpipe/psm_adapter.py:144-199) | [[FRAGPIPE-03](findings-psm.md#fragpipe-03)] |
| Proteome Discoverer | All ranks unless `--rank1-only` (qpx/converters/proteome_discoverer/psm_adapter.py:37-40); `PSM Ambiguity` only in cv_params (:76-82). A feature sums the reporters of all rank-1 PSMs of a precursor, excluding only `Rejected` and Quan Info `Excluded by…` (feature_adapter.py:16-20,64-76); always `direct` (:139). Two `Ambiguous` PSMs on scan 100 give two features with TMT126=10 each, anchored on P12345 and P77777 (audit test case `ambiguous_same_scan`). [#378](https://github.com/bigbio/qpx/issues/378): 758 of 15,582 PXD046576 LFQ spectra (Sequest HT) have rank-1 PSMs of several peptides, e.g. scan 58857, LGMELLEK 2+ and EFQILSLHGQTK 3+; 358 spectra in a TMT set | [[PD-03](findings-cross.md#pd-03)] [[PD-04](findings-feature.md#pd-04)] [[INFRA-24](findings-psm.md#infra-24)] |
| OpenMS consensus | Every hit of every PID is a row; rank is not read (qpx/converters/openms_consensus/psm_adapter.py:260-273,370-381). PIDs on one spectrum with different MZ go unmarked: MS1LabeledWorkflow_6 scan 11873 gives EYNVTIDAAESTVNK, TALGYYWALTEASGK and its Lys8 partner. A feature takes `pids[0]` (feature_adapter.py:724) | [[OPENMS-01](findings-psm.md#openms-01)] [[OPENMS-02](findings-psm.md#openms-02)] [[OPENMS-20](findings-psm.md#openms-20)] |
| mzIdentML | Every SII becomes a row; rank is a score (qpx/converters/mzidentml/psm_records.py:46-73,140). Ties and co-identifications look alike, although the per-SII `experimentalMassToCharge` is read (:150-160) | [[MZID_CDAP-01](findings-psm.md#mzid_cdap-01)] [[MZID_CDAP-02](findings-psm.md#mzid_cdap-02)] [[EXTERNAL-01](findings-psm.md#external-01)] |
| CDAP | `AmbiguousMatch` ties are rows marked only by a cv_param (qpx/converters/cdap/psm_adapter.py:26-27) and are left out of features: "one MS2 reporter signal cannot be assigned independently" (qpx/converters/cdap/feature_adapter.py:7-9,116-120) | [[MZID_CDAP-16](findings-feature.md#mzid_cdap-16)] [[CRITIC-M3](findings-psm.md#critic-m3)] |
| Sage, CHIMERYS, DDA+ (no qpx converter) | Sage `chimera` puts co-identifications into `rank`. CHIMERYS and DDA+ report several PSMs per spectrum with no primary. CHIMERYS flags I/L ties as ambiguous. DDA+ re-assigns precursor m/z. PSI-MS has no chimeric term, and USI `+` is DIA-only | [[EXTERNAL-04](findings-psm.md#external-04)] [[EXTERNAL-06](findings-psm.md#external-06)] [[EXTERNAL-07](findings-psm.md#external-07)] [[EXTERNAL-22](findings-feature.md#external-22)] |

## Why it matters
- Chimeric spectra are common: 4.9 % of PD LFQ spectra, 2.1 % of MaxQuant evidence rows.
- `rank = 1` drops Sage co-identifications; keeping every row counts alternatives. Spectral counts double-count and `(run, scan)` joins fan out.
- In TMT data one reporter vector enters two peptides and proteins; CDAP and PD disagree.
- A tie inside a chimeric spectrum (CHIMERYS I/L, PD `Ambiguous`) must count once, a second precursor twice.

## Options considered
| Option | For | Against |
|---|---|---|
| A. Status quo | No change | Co-identifications look like alternatives, ties or duplicates |
| B. One row per spectrum with a `+`-joined `peptidoform` | ProForma 2.0 §7.2 | Breaks the per-peptide typed `charge`, scores, proteins and `feature_id`; no qpx converter writes it |
| C. Boolean `is_chimeric` | Smallest | No representative; ties look like second precursors |
| D. Co-identifications in `rank` (Sage; per-SIR mzIdentML) | No column | Rank filters drop valid identifications [[EXTERNAL-06](findings-psm.md#external-06)] |
| E. Role `none`/`primary`/`secondary` | Readable | Two `secondary` rows may be one precursor's tie or two precursors (CHIMERYS, PD, MaxQuant); readers fall back on m/z |
| **E′. `precursor_index` with per-query rank ([T01](t01-psm-rows.md))** | Groups hypotheses by precursor; `= 1` / `> 1` filters; ties count once | New optional column; the producer picks index 1 when the tool does not |
| **F1. Keep and define `secondary`** | Non-breaking; mokume already reads it as in-run | Mixes origin and role ([EXTERNAL-22](findings-feature.md#external-22)), which cannot combine within one run |
| F2. New `feature.chimeric_role` | Orthogonal | Deprecates a 1.3 value |
| Q1. Shared reporters | No loss | Rollups double-count |
| **Q2. One reporter spectrum, at most one feature** | Rollups safe; CDAP and OpenMS isobaric comply | Secondary features lose reporters |
| Q3. Apportion reporters | Uses all signal | Invents intensities (P-A) |

## Recommendation
**CHI-1 (MUST)** One peptidoform per row; co-identifications are separate rows sharing `run_file_name` and `scan`. `psm.peptidoform` and `feature.peptidoform` hold no chimeric joiner: with `[...]` and `{...}` contents removed, no `+` remains. scope.md:41 moves `PEPTIDE/2+SEQUENCE/3` from "passed through" to "invalid" (cross-links stay); the validator checks files of `qpx_version` ≥ 1.4. Readers MAY build `A/2+B/3` for export.

**CHI-2 (MUST)** Definitions:
- A *precursor* is one precursor ion (m/z peak) of an MS2 spectrum.
- A *co-identification* identifies another precursor of the same spectrum. Co-isolated SILAC partners count.
- An *alternative* is another hypothesis for the same precursor: a lower rank, a rank-1 tie, another charge, or another engine ([T01](t01-psm-rows.md)).

**CHI-3 (MUST)** `psm.precursor_index` (int16, ≥ 1, nullable) names the spectrum's precursor that the row explains. On one `(run_file_name, scan)`, equal indices are hypotheses for one precursor; different indices are different precursors.
- Index 1 is the precursor the tool designates (MaxQuant first peptide, Sage chimera hit 1), else the one whose best rank-1 row sorts first under [T06](t06-identifications-per-feature.md) FID-3. Others take 2, 3, … in the tool's order, else FID-3 order. Index 1 is a representative, not a claim that the instrument targeted it.
- A non-chimeric spectrum has only index 1. Null means unknown.
- A *secondary PSM* has index > 1; other sections use this term.

**CHI-4 (MUST)** Rank follows [T01](t01-psm-rows.md) PSM-3: each co-isolated precursor is its own query, so its best hypothesis has rank 1. A per-spectrum tool rank goes to `<engine>_rank` (`sage_rank`). A per-spectrum rank (Sage, mzIdentML per SIR) would make rank filters drop correct identifications [[EXTERNAL-06](findings-psm.md#external-06)].

**CHI-5 (MUST)** The index comes from the tool's marker where one exists (Producer changes). Otherwise the producer groups a spectrum's rank-1 rows by `observed_mz`: rows within the search's precursor tolerance (10 ppm if unknown) are one precursor. A rank > 1 row takes the index of the group its `observed_mz` matches, else null. Without `observed_mz` the index is null.

**CHI-6 (MUST)** Invariants:
- Rank-1 rows of different indices on one spectrum SHOULD differ in `observed_mz` by more than 10 ppm; the validator warns otherwise.
- Index 1 may be absent (not reported, or removed by a declared [T03](t03-filtering-state.md) filter); readers must not assume it.
- K ([T01](t01-psm-rows.md) PSM-1) stays unique. [T01](t01-psm-rows.md) PSM-2 merging happens before indexing: a secondary hit with the primary's peptidoform and charge merges into it and keeps the lower index (MaxQuant 167037/167038).

**CHI-7 (MUST)** `observed_mz` and `charge` describe the claimed ion ([T11](t11-spectrum-precursor-fields.md) SPC-1, SPC-6): for a secondary PSM the co-isolated ion, not the selected precursor, which MAY go to cv_param `MS:1000744`.

**CHI-8 (SHOULD)** Secondary rows keep the source's `qvalue`, PEP and `pass_threshold` ([T04](t04-confidence-fields.md) CNF-10), with the same meaning as on index-1 rows. A converter that knows its source's FDR left co-identifications out nulls `qvalue` and `pass_threshold` there and keeps the tool value under an engine-prefixed name (P-A). This binds converter authors; validators cannot test it.

**CHI-9 (MUST)** For every producer, feature `identification_type = 'secondary'` means MS/MS of this run (or of its unit's runs, [T05](t05-feature-rows.md)) identified the row and every supporting PSM there ([T06](t06-identifications-per-feature.md)) is secondary. Without psm-level indices the tool's type decides (MaxQuant Type ending in `SECPEP`). Any index-1 supporting PSM makes the row `direct`. Secondary rows keep [T06](t06-identifications-per-feature.md)'s `direct` invariants (FID-1, FID-5, FID-6). De-duplication follows [T08](t08-match-between-runs.md) MBR-7 (`direct` > `secondary` > `transferred`).

**CHI-10 (MUST)** No footer key declares chimeric handling; `precursor_index` alone carries it.
- Every rank-1 psm row has an index whenever the producer knows it (tool marker or CHI-5); a source that reports one precursor per spectrum gives index 1. Otherwise the index is null (unknown). Features follow CHI-9.
- Converters keep every reported co-identification by default ([#241](https://github.com/bigbio/qpx/issues/241)), apart from declared [T03](t03-filtering-state.md) filters.
- A producer that removes rows or features with index > 1 records the provenance option `chimeric_secondaries_removed = true`.

**CHI-11 (MUST)** The reporter ions of a chimeric spectrum (its MS2 reporters, or the MS3 `quant_scan` taken from it) give intensities to at most one feature row: an index-1 hypothesis's (ties follow [T06](t06-identifications-per-feature.md) FID-10), or none. This is the chimeric case of [T05](t05-feature-rows.md) FRW-2. A secondary feature MAY keep its row with null intensities (intensities.md:123; [T07](t07-unidentified-features.md)).

No footer key declares it. A producer that cannot ensure it, because it cannot tell which spectra fed its intensities (tool-aggregated reporter sums), says so in converter-coverage. A pg computed by qpx complies when its features do; a tool's pg is listed in converter-coverage unless the tool documents the rule. Readers de-duplicate on `quant_scan` whenever it is set (CHI-R6). MS1 co-identifications are separate peaks; one peak given to two peptidoforms [[MAXQUANT-08](findings-feature.md#maxquant-08)] is a [T06](t06-identifications-per-feature.md) FID-10 conflict.

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| `precursor_index` | psm | int16, null, ≥ 1; null = unknown | Optional; non-null on rank-1 rows whenever the producer knows the precursor (CHI-10) | No (additive) |
| `secondary` defined (CHI-9) | feature | existing value | – | No. In 1.3 files, `direct` may include secondary-only rows |
| No `+` joiner | psm, feature | validator error, `qpx_version` ≥ 1.4 | – | Yes for third-party writers. No qpx converter emits it (qpx/converters/ptm.py:611-614 passes it through) |
| Provenance option `chimeric_secondaries_removed` | psm, feature | `true` | When rows or features with index > 1 were removed (CHI-10) | No |
| cv_param `MS:1000744` | psm | selected precursor m/z | MAY | No |
| Validator | psm, feature | CHI-1, CHI-6; CHI-10: no index > 1 when provenance records `chimeric_secondaries_removed`; CHI-11: no two intensity-bearing rows share a non-null `(run_file_name, quant_scan)` | – | No |

## Producer changes
| Producer | Change needed |
|---|---|
| MaxQuant | PSM: non-SECPEP index 1; SECPEP 2, 3, … in FID-3 order; twins merge into index 1. Feature: every `*SECPEP` Type without `MATCH` is `secondary` (a transferred second peptide is `transferred`, [T08](t08-match-between-runs.md) MBR-1); dedup ranks type before PEP (feature_adapter.py:264). TMT: null `secondary` reporters; listed in converter-coverage as not ensuring CHI-11 until a fixture shows how evidence reporters are summed |
| FragPipe | CHI-5 on `Observed M/Z` (no tool marker). ion.tsv reporter sums: listed in converter-coverage as not ensuring CHI-11 |
| Proteome Discoverer | CHI-5 on `m/z [Da]`; `Ambiguous` ties share an index. A chimeric scan's reporters go only to an index-1 feature ([PD-03](findings-cross.md#pd-03)), once across engines ([PD-04](findings-feature.md#pd-04)) |
| OpenMS consensus | CHI-5 on PID MZ; a PID's hits share its index. Isobaric: one peptidoform per consensus feature, by FID-3 not `pids[0]` |
| OpenMS `-out_qpx` | Verified at HEAD: several identifications on one spectrum with different peptidoforms stay separate rows, each with its own m/z (CHI-1, CHI-7); no `precursor_index` column. IsobaricWorkflow credits a spectrum's reporters to one peptidoform, which meets CHI-11. Needed: `precursor_index` per CHI-5. See the [OpenMS worklist](openms-worklist.md). |
| mzIdentML | CHI-5 on `experimentalMassToCharge` |
| CDAP | MS-GF+ searches one precursor: index 1. Ties are left out of features, which meets CHI-11 |
| quantms/MSstats | No psm view. TMT: CHI-11 holds if no `(run_file_name, scan)` feeds two features; else keep the FID-3 first, or list the gap in converter-coverage |
| Third-party (Sage, CHIMERYS) | Sage `chimera=true`: hit n gets index n, rank 1, `sage_rank = n`; otherwise index 1 and Sage's rank. CHIMERYS: CHI-5; I/L ties share an index |
| DIA-NN, Spectronaut | None (DIA) |

## Reader rules
**CHI-R1 (MUST)** Identified spectra are distinct `(run_file_name, scan)` values, or `spectrum_reference` where `scan` is empty ([T11](t11-spectrum-precursor-fields.md) SPC-9), over accepted rank-1 target rows ([T01](t01-psm-rows.md) PSM-R2). Never count them as rows.

**CHI-R2 (MUST)** PSMs are identified precursors, counted per [T01](t01-psm-rows.md) PSM-R2 as distinct `(run_file_name, scan, precursor_index)`: a tie counts once, each secondary precursor counts. For one precursor per spectrum, take `precursor_index = 1`; this drops spectra whose index 1 is absent (CHI-6).

**CHI-R3 (MUST)** Never key PSMs on `(run_file_name, scan)`. Join on `feature_id`, else on K.

**CHI-R4 (MUST)** "Identified in this run" is `identification_type IN ('direct','secondary')` ([T08](t08-match-between-runs.md) MBR-R1).

**CHI-R5 (MUST)** A spectral count (peptide, protein, group) counts distinct `(run_file_name, scan)` among accepted PSMs; summing feature `spectral_count` overcounts a scan whose co-identifications share a group.

**CHI-R6 (MUST)** Before rolling up reporter intensities across peptidoforms, de-duplicate on `(run_file_name, quant_scan)` wherever `quant_scan` is set (keep the lowest index of the identifying PSM, then FID-3). Rows without `quant_scan` rely on CHI-11.

**CHI-R7 (MUST)** A null or absent `precursor_index` is unknown, never 1; such rows may be co-identifications. Rows with index > 1 were removed only if provenance records `chimeric_secondaries_removed = true`.

Violations:
- qpx counts no spectra: `n_psms` = `COUNT(*)` (qpx/converters/summary.py:109); `ModificationView.frequency()` counts rows (qpx/views/api.py:139-151); `PSM.confident()` ignores rank and index (qpx/core/data/psm.py:23-25).
- No qpx reader reads `identification_type`, so CHI-R4 is never applied.
- Isobaric features are summed unchecked by `intensity('peptide')` (qpx/core/dataset_queries.py:124-137), `build_mudata(duplicates='sum')` (qpx/mudata.py:691) and `transform quantify` (qpx/cli/transform/quantify.py:261-272). The PD tie credits TMT126=10 to both P12345 and P77777.
- mokume counts `psm_id` rows per group (mokume:rust/crates/mokume-pipeline/src/spectral_count.rs:235-249; CHI-R5), aborts on MaxQuant's duplicate SECPEP ids (:210-231), and reads `scan` but no rank or index (mokume:rust/crates/mokume-io/src/psm.rs:14-25). It keys TMT rows on `(peptidoform, charge, plex)` (mokume:rust/crates/mokume-pipeline/src/tmt.rs:116-124), never on `quant_scan` (CHI-R6). Its reading of `secondary` as in-run (mokume:rust/crates/mokume-io/src/qpx.rs:501-504) conforms to CHI-R4.

## Open questions for maintainers
1. MaxQuant TMT: does an evidence row's reporter sum include MS/MS where it was a second peptide? A TMT fixture decides whether MaxQuant meets CHI-11 or is listed in converter-coverage.
2. Do FragPipe DDA+ `Observed M/Z` and PD-with-CHIMERYS `m/z [Da]` give each co-identification its own precursor ([FRAGPIPE-03](findings-psm.md#fragpipe-03))? If not, CHI-5 merges them into ties and the producer must leave the index null.
3. Should index 1 go to the lowest PEP, which every producer can compute, or to the precursor nearest the isolation target, which needs mzML?
