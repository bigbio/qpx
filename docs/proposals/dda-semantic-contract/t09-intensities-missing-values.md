# T09. Intensity semantics and missing values

[← T08](t08-match-between-runs.md) · [Overview](index.md) · [T10 →](t10-protein-mapping.md)

!!! note "Proposal"
    Part of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** The primary holds the least-processed quantity the input exports, from one column per file, declared by two footer keys (`intensity_quantity`, `intensity_processing`). A feature row has one entry per sample label of its run (null = not quantified, 0 = measured zero, no entry = no sample); imputation is marked per row (`[]`/`false` = none, null = unknown). Isobaric sums state their spectrum count, and QPX roll-ups are never primary.

## Questions a reader has
- Is the primary an MS1 apex, area or trigger intensity, a reporter intensity or S/N, a sum, MaxLFQ? Raw, normalized or corrected? Can two files be pooled?
- What do a null entry, a 0, a missing entry, an all-null row and a missing row mean? Which values were imputed?
- How many spectra does an isobaric sum hold, and whose PIF and `quant_scan` apply?
- Which `intensity_name` means what? What did QPX compute? Quantified vs used for quantification ([#282](https://github.com/bigbio/qpx/issues/282))?

## What happens today
- **The spec contradicts itself.**
  - The primary is "the raw measurement" (docs/spec/intensities.md:19,22,219), and QPX computes no processed values (:3).
  - Yet `pg.intensity` is "the primary quantity of the producer", named only by an optional per-row cv_param `quantification_method` (docs/spec/pg.md:54-60,110), and converter-coverage.md:107 accepts an OpenMS QPX sum [[OPENMS-24](findings-pg.md#openms-24)].
  - There is no vocabulary for apex, area or S/N [[CRITIC-38](findings-feature.md#critic-38)] [[EXTERNAL-17](findings-feature.md#external-17)] [[EXTERNAL-20](findings-feature.md#external-20)].
  - The provenance keys `isotope_corrected`/`normalization` are shared by all views of a step (intensities.md:211-213).
- **Missing values** (intensities.md:121-129): null entry = identified, not quantified; no entry = not identified or no sample; 0 = measured. producer-guide.md:48 names a non-existent `imputation` flag [[INFRA-M4](findings-feature.md#infra-m4)]; a null `imputed_labels` means "none, or not said" (qpx/core/data/schemas/feature.yaml:239); `intensity_value` has no missing-value rule (types.yaml:37-40) [[CRITIC-42](findings-cross.md#critic-42)].
- **Names.** MaxLFQ is `lfq` for MaxQuant, `maxlfq` elsewhere (pg.md:228-229) [[MAXQUANT-17](findings-pg.md#maxquant-17)]; examples use `LFQ`/`iBAQ` (pg.md:206-207) and repeat a label (intensities.md:158-176) [[CRITIC-43](findings-cross.md#critic-43)]; `ms1_intensity` sits under the "first quantified channel" (:189).
- PIF and `quant_scan` describe one spectrum (docs/spec/feature.md:154-156); nothing counts summed spectra.

| Producer | Primary | Missing values | Findings |
|---|---|---|---|
| MaxQuant | LFQ: `LFQ intensity` if > 0, else `Intensity`, per row (qpx/converters/maxquant/feature_adapter.py:521-525; latent). TMT: plain `Reporter intensity N`, summed over up to 9 MS/MS; the best MS/MS's `PIF` on 7,666/9,301 maxquant_simple rows (evidence_fields.py:30); MS1 `Intensity` dropped. Site: corrected values (sites.py:19,111) under `isotope_corrected=false` (converter.py:122). In MQ 1.5 plain = corrected | 0 → null; SDRF channels only (feature_adapter.py:537-540). Without an SDRF, TMT is written as `LFQ`. All-null rows: 1,469 MBR, 161 `MSMS`, 1,541 `MULTI-MSMS` | [[MAXQUANT-09](findings-feature.md#maxquant-09),12..17,21] [[EXTERNAL-18](findings-feature.md#external-18),19] |
| FragPipe | `<exp> Intensity`, which IonQuant normalizes by default. pg: the first positive of `Intensity`, `Total Intensity` and `MaxLFQ Intensity`, per cell (qpx/converters/fragpipe/pg_adapter.py:71) | 0 → no row (feature_adapter.py:331-333,376-377). TMT all-zero → no row (isobaric.py:44-50). TMT pg drops null channels (pg_adapter.py:381-383) | [[FRAGPIPE-16](findings-feature.md#fragpipe-16),17,19,22,M1,M4] |
| Proteome Discoverer | LFQ: the MAX of `Precursor Abundance`, else `Intensity` (qpx/converters/proteome_discoverer/feature_adapter.py:69). TMT: the SUM over PSM rows (:66), so two engines double it. PIF from one PSM (:114-116,140). pg: `Abundance`, else `Abundances (Grouped)` (pg_adapter.py:74-76) | An entry per reporter column, whatever the SDRF (feature_adapter.py:86-89). Unquantified rows dropped (:91-92), and pg cells ≤ 0 too (pg_adapter.py:113) | [[PD-04](findings-feature.md#pd-04),12,13,17,18,19,M4,M5] [[EXTERNAL-M5](findings-pg.md#external-m5)] |
| CDAP | LFQ: the best PSM's `PrecursorArea`. TMT: the SUM over PSMs (qpx/converters/cdap/feature_adapter.py:375-379). The count is only the cv_param `psm_count` (:384-388). pg: a QPX sum (pg_adapter.py:188-207) | 0 → null entry | [[MZID_CDAP-15](findings-feature.md#mzid_cdap-15),23,24,M4] |
| OpenMS consensus | `it`; a channel linked twice keeps the max (qpx/converters/openms_consensus/feature_adapter.py:358-368). pg: the OpenMS abundance, else a QPX sum (pg_adapter.py:323-340,414-429) | ≤ 0 → no entry (:360-361), so the identified MYFINGERK has 4/10 | [[OPENMS-13](findings-feature.md#openms-13),14,24] |
| OpenMS -out_qpx | Passthrough; no fixture | Unverified | [[OPENMS-31](findings-cross.md#openms-31)] |
| quantms-MSstats | `Intensity` | NA aborts; 0 kept (qpx/converters/quantms_msstats/feature_adapter.py:394,405-406) | [[INFRA-13](findings-feature.md#infra-13),14,29] |
| DIA-NN (pg) | `PG.Quantity`, else `PG.MaxLFQ`, per row (qpx/converters/diann/pg_adapter.py:215-220) | | [[SPEC-22](findings-cross.md#spec-22)] |

## Why it matters
One column holds apex heights, areas, reporter S/N, sums over 1-9 spectra, normalized values and MaxLFQ, and `merge()` checks only the pg cv_param (qpx/core/dataset_merge.py:276-291). "Identified, not quantified" has three encodings [[INFRA-13](findings-feature.md#infra-13)] [[SPEC-23](findings-feature.md#spec-23)], so missingness models depend on the producer. mokume may correct reporters twice [[MOKUME-26](findings-feature.md#mokume-26)]. A null imputation flag also means "unknown".

## Options considered
| Option | For | Against |
|---|---|---|
| Raw only, else null | Strict P-A | FragPipe LFQ and DIA-NN 2.x pg lose every quantity |
| **Least-processed, per file, declared** | Decidable from columns (P-G); one kind per file | Readers must read the footer |
| Declare per row or per provenance step | Exists | Lost on merge; one step serves several views |
| **One entry per sample label of the run** | Testable against `run.samples` | OpenMS, PD and FragPipe change |
| Imputed values to additional | Primary = measured | Breaks 1.3 and mokume [[MOKUME-19](findings-feature.md#mokume-19)] |
| **Typed spectrum count** | Testable | One column |
| QPX roll-ups as pg primary | pg never empty | Contradicts intensities.md:3,196 |

## Recommendation
**INT-1 (MUST)** `feature.intensities[].intensity` and `pg.intensity` hold the least-processed quantity of one kind that the converter's input exports for the view. For features this is an MS1 or reporter signal; for pg, it is the tool's own protein quantity, and a tool's own top-N or sum roll-up counts as that base quantity.
- Order, least processed first:
  1. measured signal;
  2. isotope-corrected reporters, only when no uncorrected value exists;
  3. normalized values;
  4. cross-run estimates (MaxLFQ).
- The column is chosen once per file from column presence and settings, never per row or cell. A row without a value there gets null.
- Ratios, logs, iBAQ and values labelled scaled (PD `Abundances (Scaled)`) are never primary. PD `Abundances (Grouped)` does not state its scaling. It is primary only without a per-file `Abundance`, declared `intensity_processing = unknown`.
- A stated deviation from P-A: the column means "the tool's base quantity", and INT-2 declares its kind; raw-only would empty FragPipe LFQ features and DIA-NN 2.x pg.
- Replaces intensities.md:3,17-22,219, pg.md:54-60,110 and anndata.md:147. The site primary keeps site.md:112 and only declares itself (INT-2).

**INT-2 (MUST)** Footer keys. An absent key means `unknown`, and any other value is an error. The keys are SHOULD on 1.4 files and MUST on feature and pg under the [T13](t13-metadata-conformance.md) profile ([T13](t13-metadata-conformance.md) V3). They are carried by [T13](t13-metadata-conformance.md) META-4's default rule: when sources disagree, the result is `unknown`, with a warning.

| Key | Views | Values |
|---|---|---|
| `intensity_quantity` | feature, pg, site, unidentified_feature ([T07](t07-unidentified-features.md)) | `ms1_area` (MS:1001844), `ms1_apex` (MS:1001843), `ms1_trigger` (survey-scan precursor, MS:1003085), `ms1` (kind unstated), `reporter_uncorrected`, `reporter_corrected` (isotope-impurity correction not applied / applied; MS:1001847), `reporter` (MS:1001847, correction state unknown), `reporter_sn` (no PSI term), `fragment` (DIA), `none` (every intensity null), `unknown` |
| `intensity_processing` | feature, pg, site, unidentified_feature | `none` (the measured signal, not normalized), `sum` (pg: a plain sum of measured signals), `topn`, `normalized`, `maxlfq`, `other`, `unknown`. If several apply, write the most processed: `maxlfq` > `normalized` > `topn` > `sum` > `none` |

On pg, `intensity_quantity` names the signal the roll-up combines. Details go in provenance: the normalization method (`normalization`), the peptides summed (razor, unique or all) and N for top-N. Imputation has no key (INT-9). PSI-MS accessions were checked in OLS.

**INT-3 (MUST)** The keys supersede the provenance keys `isotope_corrected`/`normalization` and the cv_param `quantification_method`. Writers keep the provenance keys as tool settings and MAY keep the cv_param.

`qpxc migrate` writes a key only when it is certain:
- from a cv_param that every quantified row shares (`PG.MaxLFQ` → `intensity_processing = maxlfq`);
- from OpenMS IsobaricAnalyzer parameters.

It never uses MaxQuant `isotope_corrected=false`, which is wrong for the site view and MQ 1.5 [[MAXQUANT-14](findings-cross.md#maxquant-14)] [[EXTERNAL-19](findings-feature.md#external-19)]. It applies INT-4 to primaries marked `unnormalized_unique_peptide_sum` and to pg files with `creator = cdap`.

**INT-4 (MUST)** A quantity QPX computes (the OpenMS and CDAP pg sums) never fills a primary column. It goes to `additional_intensities` as `unnormalized_unique_peptide_sum` with a provenance step (intensities.md:196), and a pg file without a tool quantity declares `intensity_quantity = none`. This resolves pg.md:240 and converter-coverage.md:107 [[OPENMS-24](findings-pg.md#openms-24)].

**INT-5 (MUST)** Every feature row has a non-null `intensities` list with exactly one entry per label that has a sample in the row's run, or, for unit rows (`identity_composite` starting with `quantification_unit_id`, [T05](t05-feature-rows.md) FRW-4), in its unit's runs (`run.samples`), and no other entry. LFQ has one `LFQ` entry.
- Without a sample map, write one entry per channel the source reports for the run.
- An isobaric source without a channel map is an error, never `LFQ` [[MAXQUANT-15](findings-feature.md#maxquant-15)].
- Identification state belongs to `identification_type` and [T07](t07-unidentified-features.md), not to the entry set. Transferred rows follow [T08](t08-match-between-runs.md) MBR-3.
- Replaces intensities.md:125-126 and producer-guide.md:45.
- The pg softlink (pg.md:73) then reaches every label row of a group.

**INT-6 (MUST)** States of a (row, label) pair:

| State | Encoding |
|---|---|
| Quantified | A finite value ≥ 0 (0 per INT-7) |
| Not quantified in this row | An entry with null `intensity` |
| No sample in the run | No entry |
| Not reported in the run | No row; [T07](t07-unidentified-features.md) `unquantified_policy` says why |

A null entry beside a non-null one means the label was not detected (a channel below detection, a missing SILAC partner). An all-null row is identification-only ([T05](t05-feature-rows.md) FRW-5), such as a MaxQuant `MSMS` row, whose synthetic 60 s window [T05](t05-feature-rows.md) FRW-9 nulls [[MAXQUANT-09](findings-feature.md#maxquant-09)]. pg keeps pg.md:116-146.

**INT-7 (MUST)** A primary `0.0` is a value measured as zero. A converter writes a source 0 as null unless the source documents 0 as measured, and converter-coverage records which applies; MaxQuant, FragPipe, PD, CDAP, OpenMS, DIA-NN and quantms-MSstats (open question 5) write null. NaN, ±inf, negatives and an empty list are write errors. A transfer whose only value is such a 0 (andes#148) is left out ([T08](t08-match-between-runs.md) MBR-3).

**INT-8 (MUST)** In `additional_intensities`, each label appears at most once and only labels of the row's entries (pg: the row's label). An absent name means not reported; sentinels are dropped, so a 0 is real (a count). Values are finite, and negative only under a registered log-scale name (INT-12).

**INT-9 (MUST)** Imputed values stay entries, flagged by `imputed_labels` or `is_imputed` (intensities.md:131).
- `imputed_labels` ⊆ the labels of non-null entries, and `is_imputed = true` ⇒ a non-null `intensity`.
- `imputed_labels = []` means nothing in the row was imputed; pg `is_imputed = false` means measured. A null flag means unknown. Producers that know write `[]`/`false`.
- Imputation never sets `identification_type` ([T08](t08-match-between-runs.md) MBR-6). Fix producer-guide.md:48.

**INT-10 (MUST)** The new feature column `quant_spectrum_count` (int32, null, ≥ 1) counts the distinct spectra `(run_file_name, scan)` whose reporter ions are summed into `intensities`. It is null on MS1 rows and when unknown.
- Sums run over distinct spectra, never over PSM rows [[PD-04](findings-feature.md#pd-04)], after [T02](t02-chimeric-spectra.md) CHI-11.
- At 1, `quant_scan` is that spectrum and `parent_ion_fraction` its purity (MS:1003167 [[EXTERNAL-21](findings-feature.md#external-21)]). Both MAY be null when the producer cannot tell which spectrum it was.
- Above 1, both are null. An aggregate purity goes to `additional_scores` with a tool prefix (`maxquant_reporter_pif`).
- `quant_ms_level` stays set when every spectrum shares the level.
- Choosing one MS1 peak among several (PD max, CDAP best PSM) is [T05](t05-feature-rows.md) FRW-3's grain.

**INT-11 (MUST)** A row-scope name is a precursor value in a multi-label row (`ms1_intensity`, or `spectral_count` on isobaric rows). It is written once, on `intensities[0].label`; this replaces intensities.md:189. `spectral_count` counts what [T06](t06-identifications-per-feature.md) FID-12 defines: the distinct spectra of the row's supporting PSMs in its run (its unit, for unit rows, [T05](t05-feature-rows.md) FRW-4) that pass the file's filters, whether or not their reporters were summed.

**INT-12 (MUST)** Names are lower snake_case, from a registry in intensities.md that gives each name's level, scope (label or row), scale (linear or log2) and sources.
- Renames: `lfq` → `maxlfq`; `LFQ`/`iBAQ` (pg.md examples) → `maxlfq`/`ibaq`. Readers accept the old names in 1.x.
- Fix the repeated-label example (intensities.md:158-176).

**INT-13 (SHOULD)** [#282](https://github.com/bigbio/qpx/issues/282): "quantified" = a non-null entry. "Used for quantification" is never inferred; a source that states it (PD `Quan Usage`) MAY write the cv_param `used_for_quantification`. No typed column until two producers can fill it.

**INT-14 (MUST)** Validator and writers; severities follow [T13](t13-metadata-conformance.md) META-8:

| Check | Severity |
|---|---|
| Off-vocabulary footer value; imputation flag on a null or absent entry; repeated or foreign additional label; count < 1, or > 1 with `quant_scan`/PIF; PIF ∉ [0,1]; `quant_ms_level` ∉ {2,3} [[VALIDATOR-22](findings-feature.md#validator-22)]; non-null intensity under `none`; null-label pg row with an intensity [[VALIDATOR-28](findings-pg.md#validator-28)]; case variant of a registered name (1.4 files) | Error |
| Non-finite or negative primary; non-finite additional value; negative linear-scale value | Error, also on write (qpx/writers/base.py:416 [[VALIDATOR-M4](findings-cross.md#validator-m4)]) |
| Null or empty list; entry without a sample; sample label without an entry | Error under the profile, else warning |
| Missing `intensity_quantity` or `intensity_processing`; a 0.0 primary from an INT-7 converter | Warning (missing keys are errors under the profile) |
| `no_pg_quantity` (qpx/core/data/schema.py:297-313) under `none`; an unregistered name | Info |

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| `quant_spectrum_count` | feature | int32, null, ≥ 1 | SHOULD on reporter rows | No |
| `intensity_quantity`, `intensity_processing` | feature, pg, site, unidentified_feature footers | INT-2 | SHOULD; MUST on feature and pg under the profile | No |
| Imputation per row: `imputed_labels = []` / `is_imputed = false` = none, null = unknown (INT-9) | feature, pg | Constraint | Producers that know write `[]`/`false` | No |
| "Raw" promise withdrawn; one column per file (INT-1) | feature, pg | Constraint | Yes | Semantic: 1.3 files read as `unknown`; some FragPipe and DIA-NN pg cells → null |
| Complete label set, no null list (INT-5) | feature | Constraint | Yes | Semantic: OpenMS, FragPipe TMT and PD gain null entries |
| No QPX roll-up as primary (INT-4) | pg | Constraint | Yes | Semantic: OpenMS without abundance and CDAP pg → null |
| Renames (INT-12) | feature, pg | Registry | – | Yes for `lfq` readers; `migrate` renames |
| `quantification_method`; provenance intensity keys | pg, provenance | Superseded as declarations; provenance keeps the details (INT-2) | – | No |

## Producer changes
| Producer | Change needed |
|---|---|
| MaxQuant | **Feature:** always `Intensity` (remove feature_adapter.py:521-525); `ms1`, or apex/area from a supplied mqpar.xml.<br>**TMT:** `Reporter intensity not corrected N` when present (MQ 1.5), else the plain column; `reporter_uncorrected` if a `corrected` column exists, else `reporter`. Count = largest `Reporter intensity count N`, null at 0. `quant_scan`/PIF only at `MS/MS count` = 1 (25 maxquant_simple rows have reporter count 1 but 2 MS/MS). `Reporter PIF` → `maxquant_reporter_pif`; MS1 `Intensity` → `ms1_intensity`; `MS/MS count` → `spectral_count`. Fail without a channel map.<br>**Other:** site `reporter_corrected`; pg processing `sum` (razor peptides, in provenance), `lfq` → `maxlfq`; `imputed_labels = []`, `is_imputed = false` |
| FragPipe | IonQuant settings from an optional fragpipe.workflow (as [T08](t08-match-between-runs.md)): `normalization` → `intensity_processing` `normalized`/`none`; `intensitymode` 0/1/2 → `ms1_apex`/`ms1_area`/`ms1`; without it `unknown`/`ms1`. Null entries for all-zero ions ([T07](t07-unidentified-features.md)). TMT: `reporter`, count null; reject ratio TMT-Integrator reports [[FRAGPIPE-M4](findings-pg.md#fragpipe-m4)]. pg: `Intensity`, else `Total Intensity`, else MaxLFQ (`maxlfq`), per file; otherwise processing `unknown`; null rows for annotated TMT channels |
| Proteome Discoverer | SDRF labels only. Sum over distinct `(run, scan)`, with the count; `spectral_count` for TMT too; PIF null on sums. `ms1` (Precursor Abundance) or `ms1_trigger` (Intensity); `reporter`. pg: `Abundance` (not normalized; sum or top-N unstated) or Grouped, processing `unknown`; null rows for identified channels; imputation flags null |
| CDAP | Count from distinct spectra; drop `psm_count`. LFQ `ms1_area`; TMT `reporter`. pg sum → additional; pg declared `none`. Imputation: `imputed_labels = []` |
| OpenMS consensus | Null entries for labels without an element > 0 (feature_adapter.py:358-361). Keys from IsobaricAnalyzer and ProteinQuantifier parameters. LFQ `ms1`; isobaric `reporter_corrected` or `reporter_uncorrected`, count 1. pg declared `none` without an abundance |
| OpenMS -out_qpx | OpenMS writes the keys and the null entries. `qpxc convert openms` derives keys only from a companion consensusXML, and adds null entries only with an SDRF |
| quantms-MSstats | NA → null entry, instead of aborting; 0 → null. `ms1` or `reporter`; isobaric count 1 |
| mzIdentML | pg: `intensity_quantity = none` |
| DIA-NN (pg) | `PG.Quantity` per file, else `PG.MaxLFQ` declared `intensity_processing = maxlfq` |

## Reader rules
- **INT-R1 (MUST)** Pool primaries across files only when `intensity_quantity` and `intensity_processing` are equal and known (neither `unknown` nor `reporter`, whose correction state is unknown); otherwise warn. Violations: `merge()` (qpx/core/dataset_merge.py:276-291); mokume, which reads only the `qpx.struct_lists` footer.
- **INT-R2 (MUST)** Read missing values per INT-6. A null is never 0, and a 0 is never "not quantified", though a reader MAY treat it as censored. Violations, all treating zeros as missing: qpx/cli/transform/quantify.py:269; `no_pg_quantity` (qpx/core/data/schema.py:306); mokume (mokume:python/mokume/io/feature.py:212,749; low impact [[MOKUME-18](findings-feature.md#mokume-18)]).
- **INT-R3 (MUST)** Exclude imputed values from measured statistics, or carry a mask; a null flag is unknown, not measured (INT-9). Violations: `ds.intensity` (qpx/core/dataset_queries.py:124-137), qpx/core/data/feature.py:56-116, quantify.py:261-272; mokume ignores `pg.is_imputed`. Compliant: MuData `imputed_mask`, mokume's `imputed_labels` filter (mokume:python/mokume/io/_qpx.py:6-10).
- **INT-R4 (MUST)** Never re-correct `reporter_corrected` values, and warn on `reporter`. Violation: mokume:rust/crates/mokume-pipeline/src/tmt.rs:88-114.
- **INT-R5 (MUST)** Never attribute a row-scope name to a label's sample. Violation: qpx/mudata.py:192-206.
- **INT-R6 (SHOULD)** Compare reporter rows whose spectrum counts differ only through within-row channel proportions, never as absolute sums.

## Open questions for maintainers
1. Should a 0 reporter channel in a quantified spectrum be written as 0 (censored) instead of null?
2. Is the isotope state (`summed`/`monoisotopic`) needed, as a provenance parameter [[EXTERNAL-17](findings-feature.md#external-17)]?
3. Should the DIA-NN 2.x MaxLFQ-only pg be null instead of declared `maxlfq`?
4. Should imputed values leave the primary slot in 2.0?
5. Do quantms MSstats inputs use 0 for censored values (`censoredInt`)?
6. Should `used_for_quantification` become typed ([#282](https://github.com/bigbio/qpx/issues/282))?
7. Do MaxQuant `<mod>Sites.txt` tables carry uncorrected reporters, so the site primary can follow INT-1?
