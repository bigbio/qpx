# T11. Spectrum and precursor fields: m/z, charge, scan, retention time, ion mobility, run names and the canonical peptidoform

[← T10](t10-protein-mapping.md) · [Overview](index.md) · [T12 →](t12-protein-groups.md)

!!! note "Proposal"
    Part of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** `observed_mz` is the measured m/z of the claimed ion, calibration declared; `calculated_mz` and `mass_error_ppm` are computed from the row, plus a typed `isotope_error`; `charge` ≥ 1; `scan` has a declared `scan_format`, no hashes, and a verbatim `spectrum_reference` beside it; `psm.rt` is on the acquisition clock; `ion_mobility` is 1/K0; one stem function names runs; the canonical peptidoform is validated.

## Questions a reader has
- Is `observed_mz` the selected ion, monoisotopic, recalibrated or derived? Can it be 0? Is `mass_error_ppm` isotope-corrected?
- What do charge 0, null or negative mean?
- Is `scan` a scan number, an index, native-ID components or a hash? How do I reach the spectrum?
- Which clock and unit has `psm.rt`? Which quantity is `ion_mobility`? Does CV 0 mean null?
- Is a run spelled alike everywhere, and a molecule one peptidoform string?

## What happens today
Spec:
- `observed_mz` is "experimental observed" (docs/spec/psm.md:47); ontology.md:167,197 say "selected ion m/z" but cite MS:1002234, "selected precursor m/z". Calibration, isotope state and charge source are unstated [[SPEC-20](findings-psm.md#spec-20)][[CRITIC-13](findings-cross.md#critic-13)][[CRITIC-14](findings-cross.md#critic-14)][[EXTERNAL-M3](findings-psm.md#external-m3)].
- "`charge = 0` validates" (scope.md:19); psm.yaml:33-36 has no range [[CRITIC-16](findings-cross.md#critic-16)][[VALIDATOR-M6](findings-cross.md#validator-m6)].
- `psm.rt` is in seconds (psm.md:50), `mz.scan_start_time` in minutes (mz.md:30) [[CRITIC-17](findings-cross.md#critic-17)].
- `ion_mobility` has no unit (psm.md:67, feature.md:72); ontology.md:173 says 1/K0; feature.md:76 prescribes CV `0` for non-FAIMS units [[CRITIC-18](findings-cross.md#critic-18)].
- `scan_format` "can" be declared (scan.md:51; producer-guide.md:23 says "set"); mixed or unknown formats stay undeclared (scan.md:64-67); absent means scan numbers (:80-81); `nativeId` has 2–4 components (:61); features are exempt (:90-92) [[EXTERNAL-26](findings-psm.md#external-26)][[VALIDATOR-12](findings-feature.md#validator-12)].
- scope.md:27 claims ProForma Base but passes Base constructs through (:37-46) [[EXTERNAL-28](findings-cross.md#external-28)]; masses are written only "when no Unimod entry matches" (:33).

| Producer | Behaviour | Findings |
|---|---|---|
| MaxQuant | msms `m/z` is theoretical and omits isobaric labels (tool-mappings.md:26; equals pyopenms on the fixture). psm `observed_mz` = QPX `calculated_mz` + m/z·ppm (qpx/converters/maxquant/psm_adapter.py:290-294), hiding peptidoform errors; feature uses m/z·(1+ppm) (feature_adapter.py:569-571). Missing fixed mods only warn (base_adapter.py:197-202); protein-terminal ones are skipped (constants.py:129-133); NULL `Modified sequence` becomes `'None'` | [[MAXQUANT-23](findings-psm.md#maxquant-23)][[MAXQUANT-24](findings-cross.md#maxquant-24)][[MAXQUANT-M1](findings-cross.md#maxquant-m1)][[VALIDATOR-M1](findings-cross.md#validator-m1)] |
| FragPipe | `Calibrated Observed M/Z` (psm_adapter.py:239); isotope errors stay in `mass_error_ppm` (constants.py:113-140). Run and scan are the first two dot-tokens (psm_adapter.py:45-58): `run.D.v2` → `run`, scan `[]`. CV: string cv_param on psm (:317-319), `0.0` on feature (feature_adapter.py:519) | [[FRAGPIPE-04](findings-psm.md#fragpipe-04)][[FRAGPIPE-06](findings-psm.md#fragpipe-06)][[FRAGPIPE-26](findings-psm.md#fragpipe-26)] |
| Proteome Discoverer | Mass error copied from ΔM (reports.py:118); Comet's +1 to +3 Da offsets sit in `observed_mz` (converter-coverage.md:319). First Scan only (reports.py:90,116). Run from a case-sensitive regex (:99) or `run_file_stem` (:49) | [[PD-24](findings-cross.md#pd-24)][[PD-25](findings-psm.md#pd-25)][[PD-26](findings-psm.md#pd-26)] |
| OpenMS consensusXML, quantms | psm `observed_mz` 0.0 when absent (psm_adapter.py:359; streaming.py:408). Feature `observed_mz` is the consensus centroid (feature_adapter.py:645); ProteomicsLFQ elements carry FeatureFinderIdentification's *theoretical* m/z, so all 68 fixture rows have `mass_error_ppm` < 1e-9 (PSMs: median −0.3 ppm). Hash scans (psm_adapter.py:49-69); `spectrum=N` undeclared (qpx/core/scan.py:17-31). Mass shifts become 15-decimal residue masses (`K[+54.0106]`→`K[182.105597942799989]`), also in quantms (quantms_msstats/feature_adapter.py:104), which parses only `scan=` (:346) | [[OPENMS-15](findings-feature.md#openms-15)][[OPENMS-16](findings-psm.md#openms-16)][[OPENMS-17](findings-psm.md#openms-17)][[OPENMS-28](findings-psm.md#openms-28)] |
| mzIdentML | One number kept (psm_records.py:253-268): `frame=120 scan=475`→`[475]`; `spectrum=` and bare integers declared `scan`; mixed schemes undeclared (:38-43). `unknown modification` masses mapped to Unimod (qpx/converters/ptm.py:250-256) | [[MZID_CDAP-19](findings-psm.md#mzid_cdap-19)][[MZID_CDAP-20](findings-cross.md#mzid_cdap-20)] |
| CDAP | `rt` = `RTAtPrecursorHalfElution`, 0.00 kept (cdap/psm_adapter.py:118,169; 8/100). Missing values → charge 0, m/z 0.0 (:155,161). Mass error copied (:168). 3-decimal masses padded (ptm.py:217-219) | [[MZID_CDAP-17](findings-psm.md#mzid_cdap-17)][[MZID_CDAP-18](findings-psm.md#mzid_cdap-18)][[MZID_CDAP-31](findings-psm.md#mzid_cdap-31)] |
| Validator, merge | Names without `:` and masses with any decimals pass (qpx/core/semantic_validation.py:25,68-81); `'None'`, `''` fail only under `--semantic` (:85). Mass error: fixed 50 ppm warning (:146-152); charge ≤ 0 skipped (:50). `scan_format` checked for single psm files only (qpx/core/data/schema.py:134-157; base.py:173-178); merge drops it (qpx/core/dataset_merge.py:252-267) | [[INFRA-16](findings-cross.md#infra-16)][[VALIDATOR-13](findings-psm.md#validator-13)][[VALIDATOR-14](findings-feature.md#validator-14)][[VALIDATOR-30](findings-cross.md#validator-30)] |

Verifier corrections kept: unmapped names are legal (peptidoform.md:43); shifts inside `calculated_mz` are canonical [[SPEC-20](findings-psm.md#spec-20)]; 889 vs 755 is a fixture artefact [[MAXQUANT-23](findings-psm.md#maxquant-23)].

## Why it matters
These columns feed mass-accuracy QC, ML training, spectrum joins and identity (`scan`, `peptidoform` in `psm_id`; `observed_mz` in the OpenMS feature composite). A 0.0 m/z, a hash read as a scan, a residue mass read as a delta, or a truncated run all fail silently.

## Options considered
| Option | For | Against |
|---|---|---|
| Instrument selected ion only | One physical value | MaxQuant lacks it; contradicts [T02](t02-chimeric-spectra.md) CHI-7, DDA+ |
| **Value the identification used (mzIdentML `experimentalMassToCharge`), calibration declared** | Every tool has it | Calibration varies (declared) |
| **Typed `isotope_error`** | Raw error kept, offset explained | One column |
| Hash scans declared `surrogate` | Non-empty key | Sentinel (P-D) |
| **`spectrum_reference` beside `scan`** | Lossless | Must key rows with empty `scan` |
| Per-vendor cardinality table | Strict | Rejects OpenMS Bruker IDs ([#259](https://github.com/bigbio/qpx/issues/259) q3) |
| Map open-search masses to Unimod (scope.md:33) | Fewer strings | Invents chemistry; tolerance-dependent |

## Recommendation
The *claimed ion* is the precursor the row's identification asserts ([T02](t02-chimeric-spectra.md) CHI-7).

**SPC-1 (MUST)** `observed_mz` is the measured m/z of the claimed ion as the identification used it.
- Never computed from the row's `peptidoform` or `calculated_mz`, never a targeted extraction's theoretical m/z; null, not 0, when not measured.
- A value rebuilt from a tool's deviation starts from the tool's theoretical m/z, plus label masses that m/z omits.
- Feature rows: the quantified MS1 signal's m/z when the tool reports one (multiplets: [T05](t05-feature-rows.md) FRW-10); else the identifying PSM's ([T06](t06-identifications-per-feature.md) FID-3) unless its `isotope_error` ≠ 0; else null.
- The selected ion MAY go to cv_param `selected ion m/z` (MS:1000744).
- Validator: ≤ 0 is an error; |`mass_error_ppm`| < 0.001 on ≥ 50 % of rows is a warning.

**SPC-2 (SHOULD)** Footer `observed_mz_calibration` (psm, feature), one value per file: `raw`, `recalibrated` (by MaxQuant, MSFragger, PD Spectrum Files RC, OpenMS InternalCalibration) or `unknown`, also written by a producer that cannot tell; absent = `unknown`. ontology.md:167,197 say "selected precursor m/z".

**SPC-3 (MUST)** `calculated_mz` is the theoretical monoisotopic m/z (MS:1003053) of the written peptidoform at `charge`, counting every tag; cross-links follow psm.md:95. If a tag has no known mass it is null or the tool's value. The 20 ppm check is an error in the `dda` profile ([T13](t13-metadata-conformance.md) V6).

**SPC-4 (MUST)** `mass_error_ppm` = 1e6·(`observed_mz` − `calculated_mz`)/`calculated_mz` from the row's own columns, within ±0.2 ppm (float32 rounding); null exactly when an input is null; never copied, never isotope-corrected.

**SPC-5 (SHOULD)** New psm column `isotope_error` (int16, null): k with z·(`observed_mz` − `calculated_mz`) ≈ k·1.003355 Da.
- The engine's k (MS-GF+, Comet, MSFragger); else the k in [−1, 3] with residual ≤ max(0.1 Da, 20 ppm of the neutral mass), FragPipe's window (constants.py:116-119); else null. The validator warns on a k outside that window.
- An open-search shift is a peptidoform tag, never an isotope error, and vice versa.

**SPC-6 (MUST)** `charge` is the charge (MS:1000041) the identification assumed for the claimed ion: ≥ 1, non-null, in psm and feature ([T07](t07-unidentified-features.md) UNI-1).
- Schema `minimum: 1`: an error in basic validation (`qpxc validate`, [T07](t07-unidentified-features.md) UNI-1's basic tier) from `qpx_version` 1.4, a warning before.
- A source row without a charge is not written (counted as `skipped_rows_<view>`, [T03](t03-filtering-state.md) FLT-6), never written as 0.
- The instrument's charge stays in `mz.precursors.selected_ion_charge`. Negative-ion data is out of scope in 1.4 (scope.md:19).

**SPC-7 (MUST)** `scan` keeps every numeric component in native order (scan.md; [#259](https://github.com/bigbio/qpx/issues/259)). `scan_format`:

| Source identifier | `scan_format` | Components |
|---|---|---|
| Thermo `controllerType=0 controllerNumber=1 scan=N`; `scan=N` (MS:1000776); MGF `SCANS=` | `scan` | 1, ≥ 1 |
| `index=N` (MS:1000774); MGF position | `index` | 1, ≥ 0 |
| Other numeric native IDs: Thermo non-default, Bruker TDF (any shape), Waters, Sciex, Agilent `scanId=`, `spectrum=N` (MS:1000777) | `nativeId` | 1–4 |
| Several formats, every row with `spectrum_reference` (psm only) | `mixed` | – |
| Tool-internal ordinals (MaxQuant, MSFragger on `.d`; Mascot `query=`) or any format the producer cannot establish | `unknown` | – |

- Required on psm from 1.4, and on feature when any `scan` is non-empty ([T06](t06-identifications-per-feature.md) FID-5). It describes the non-empty scans; all-empty files write `unknown`.
- From 1.4 absent means `unknown`; files of 1.3 and earlier keep scan.md:80-81 (absent = `scan`).
- The validator checks component ranges, not a vendor table ([#259](https://github.com/bigbio/qpx/issues/259) q3).

**SPC-8 (MUST)** No hashes, sentinels or file positions under `scan`. An identifier without numeric components gives `scan = []`. A psm row with neither a non-empty `scan` nor a `spectrum_reference` is not written (counted as `skipped_rows_psm`, [T03](t03-filtering-state.md) FLT-6).

**SPC-9 (MUST)** New psm column `spectrum_reference` (string, null): the mzML nativeID or mzIdentML `spectrumID`, verbatim, never synthesised.
- Required when `scan` is empty or the format is `nativeId` or `mixed`; else SHOULD when the source has one (MaxQuant, FragPipe, PD and CDAP do not).
- mz views written from mzML have `id` = `run_file_name:nativeID` (qpx/transforms/spectra_mapping.py:665), stated in mz.md:20; the join is `mz.id = run_file_name || ':' || spectrum_reference`.
- It replaces an empty `scan` in the psm key K, for id derivation and the duplicate check (cf. `crosslinked_peptidoform`, qpx/core/data/footer_validation.py:152-153). [T01](t01-psm-rows.md) PSM-1 adopts this; it goes beyond [#259](https://github.com/bigbio/qpx/issues/259) A, which left the key alone.

**SPC-10 (MUST)** A combined spectrum (PD First ≠ Last Scan) keeps its first scan in `scan`; the range goes in cv_param `peak list scans` (MS:1000797).

**SPC-11 (MUST)** `psm.rt` is the scan start time (MS:1000016) of the identified MS2 spectrum, in seconds, on the acquisition clock (60 × `mz.scan_start_time`): never an elution or aligned time; ≥ 0; null when unknown or not placeable on that clock. Feature times: [T05](t05-feature-rows.md) FRW-7 to FRW-9. Validator: `rt` < 0 is an error; > 0.1 s from the joined mz spectrum a warning.

**SPC-12 (MUST)** `ion_mobility` (psm, feature) and `ion_mobility_start`/`_stop` (feature) are 1/K0 in V·s/cm² (MS:1002815), as `mz.inverse_ion_mobility`. A drift time goes to cv_param `ion mobility drift time` (MS:1002476, ms); a FAIMS CV to `compensation_voltage` (V, MS:1001581), which psm gains (float32, null). Absence is null, never 0. Validator: `ion_mobility` outside (0, 3] is a warning.

**SPC-13 (MUST)** `run_file_name` is `run_file_stem` of the source file name (qpx/core/files.py:16-19): basename minus one case-insensitive suffix from :12, inner dots and case kept. Every view, `grouped_runs` included, uses the identical string. Spectrum titles are parsed from the right. Two inputs with one stem stop conversion.

**SPC-14 (MUST)** Modifications (peptidoform.md rules 1–7 stay):
- (a) A modification the search defined is written as its Unimod accession when one exists, however the source spells it (name, PSI-MOD, mass).
- (b) A mass without a defined modification (open-search delta; mzIdentML `unknown modification` absent from `SearchModification`) is a signed delta with exactly four decimals, zero-padded, even when Unimod has that mass. scope.md:33 is reworded.
- (c) A name without a Unimod entry is kept only when the source gives no mass.
- (d) Unknown fixed mods or labels, or known ones the converter cannot place (MaxQuant protein-terminal), fail conversion unless the user asserts there are none; provenance records the assertion.

**SPC-15 (MUST)** In `peptidoform` only residue and terminal tags and the `?` prefix are canonical, with bodies `UNIMOD:n`, `UNIMOD:n|±mass`, `±mass` or an SPC-14(c) name.
- Labile `{}`, global `<>`, `/z`, `+` ([T02](t02-chimeric-spectra.md) CHI-1), `//`, `#g` groups, `( )` ranges and INFO are not; cross-links stay in `crosslinked_peptidoform`. scope.md:27 says "a subset of Base".
- Residues as identified: no I↔L rewriting, no J/B/Z/X; ties are separate rows ([T01](t01-psm-rows.md) PSM-6). `calculated_mz` and identity need one residue.
- The site is the final localisation (phosphoRS, PTMProphet or MaxQuant PTM score over the engine's); an unlocalised delta takes `?`; `assigned = true` candidates match the peptidoform.
- The untagged peptidoform equals `sequence`. These syntax and sequence checks (semantic_validation.py:85-88) move from `--semantic` to basic validation ([T07](t07-unidentified-features.md) UNI-1), so `''` and `'None'` fail there.

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| `observed_mz_calibration` (T11 key, [T13](t13-metadata-conformance.md) registry) | psm, feature | `raw`\|`recalibrated`\|`unknown` | SHOULD | No |
| `isotope_error` | psm | int16, null | SHOULD | No |
| `spectrum_reference`, standing in for an empty `scan` in K | psm | string, null | SPC-9 | **Yes**: ids of empty-scan and former hash rows; part of [T01](t01-psm-rows.md) PSM-1's K |
| `compensation_voltage` | psm | float32, null | optional | No |
| `scan_format` adds `mixed`, `unknown`; required; absent = `unknown` from 1.4; `nativeId` 1–4 | psm, feature, mz | string | MUST | **Yes**, versioned; mzIdentML `spectrum=N` and bare integers change declaration |
| `charge` ≥ 1 | psm, feature | `minimum: 1` | basic validation, `qpx_version` ≥ 1.4 | **Yes** for producers writing 0; older files warn |
| `observed_mz` > 0; `rt` ≥ 0; `mass_error_ppm` = formula | psm, feature | – | error | Values change: PD, CDAP mass error; CDAP `psm.rt`; MaxQuant psm `observed_mz` |
| 1/K0 documented; CV null, not 0 (psm.md:67; feature.md:72-76; feature.yaml:121-132) | psm, feature | – | doc | **Yes**: FragPipe feature ids ([T05](t05-feature-rows.md)) |
| `run_file_stem` everywhere | all | – | – | **Yes**: ids of PD `.RAW` and FragPipe dotted runs |
| SPC-14(b) masses stay masses; canonical checks (resolvable names, decimals ≠ 4, SPC-15) | psm, feature | – | `dda` | **Yes**: mzIdentML open-search peptidoforms, ids |
| FeatureFinderIdentification `observed_mz` null | feature | – | – | **Yes**: OpenMS feature ids ([T05](t05-feature-rows.md) FRW-4) |
| Docs: ontology.md:167,197; scope.md:19,27,33; scan.md:51-61,80-81; mz.md:20; peptidoform.md:107-109 [[SPEC-35](findings-cross.md#spec-35)] | docs | – | – | No |

## Producer changes
| Producer | Change needed |
|---|---|
| MaxQuant | psm `observed_mz` = m/z·(1 + ppm·10⁻⁶) + omitted label mass/z; `recalibrated`; derive `isotope_error`. Skip NULL `Modified sequence` (counted). SPC-14(d). `scan` when msms `Mass analyzer` is FTMS/ITMS or SDRF files are `.raw`, else `unknown` |
| FragPipe | `recalibrated` if `Calibrated Observed M/Z` differs from `Observed M/Z`, else `raw`. `isotope_error` from constants.py:116-142. One right-to-left `Spectrum` parser for psm, feature lookup and isobaric.py. psm `compensation_voltage`; feature CV null. `scan` for `.raw` inputs, else `unknown` |
| Proteome Discoverer | Compute the mass error; derive `isotope_error`; `run_file_stem` on both paths; SPC-10; `scan` when Spectrum File is `.raw`, else `unknown` |
| OpenMS consensusXML | Null, never 0.0. Feature `observed_mz` = run element ([T05](t05-feature-rows.md)), but null within 0.01 ppm of the hit's theoretical m/z (FeatureFinderIdentification). `[]` plus `spectrum_reference` instead of hashes; `spectrum=` is `nativeId`. Deltas from the modification's own mass, since normalize_peptidoform keeps residue masses. `original_RT` when present (OQ1) |
| OpenMS `-out_qpx` | Write its `spectrum_reference` ([#259](https://github.com/bigbio/qpx/issues/259)), all components with `scan_format`, null m/z, charge ≥ 1 (seeds to [T07](t07-unidentified-features.md)'s view) |
| mzIdentML | `qpx.core.scan.scan_from_native_id`; format from `SpectrumIDFormat`, bare integers without one `unknown`; `spectrum_reference`; MS-GF+ `IsotopeError`; SPC-14(b) |
| CDAP | psm `rt` null (half-elution to feature, [T05](t05-feature-rows.md)); missing m/z null; missing charge drops the row; compute the mass error; `scan`; phosphoRS sites in `peptidoform`, probabilities in `modification_candidates` |
| quantms/MSstats | OpenMS mass-shift fix; feature `scan_format` |
| qpx writers, merge | Write-time checks for SPC-4/7/8/15 (writers/base.py:525-533 checks psm `scan_format` only). Merge keeps an agreed `scan_format`; otherwise it writes `mixed` (psm, every row with `spectrum_reference`), else `unknown` with a warning ([T13](t13-metadata-conformance.md) META-4); carries `observed_mz_calibration` ([T13](t13-metadata-conformance.md)) |

## Reader rules
**SPC-R1 (MUST)** Read `scan_format` per file with its `qpx_version`. Violation: qpx/transforms/spectra_mapping.py:305-311,468 reads multi-file psm views as scan numbers even when every shard declares `index`.

**SPC-R2 (MUST)** Join to mz on `(run_file_name, scan)` only when both files declare one format other than `mixed`/`unknown`; otherwise on `spectrum_reference`. Under `unknown` a reader MAY join on `scan` if every pair agrees on `rt` = 60 × `scan_start_time` within 0.1 s (refines [T14](t14-reader-contract.md)'s absent-key row). qpx/core/data/mz.py:39-41 filters in minutes.

**SPC-R3 (MUST)** Mass-accuracy QC corrects by `isotope_error`, excludes rows with `±mass` tags (fitted to the measurement) and compares within one calibration value. Violation: qpx/core/semantic_validation.py:146-152 flags isotope offsets as `mass_error` at a fixed 50 ppm.

**SPC-R4 (MUST)** Compare `peptidoform` and `run_file_name` exactly: no case folding, suffix stripping or I/L merging. To match an external SDRF, apply `run_file_stem` to `comment[data file]` only. qpx MuData complies (qpx/mudata.py:183). mokume normalises the QPX run: Rust `normalize_file_key` lowercases and cuts at an inner `.d`, `.raw`, `.mzml` or `.wiff` followed by `.` or `_`, so `run.D.v1` and `run.D.v2` both become `run` (mokume:rust/crates/mokume-io/src/sdrf.rs:454-482, applied at :109); Python lowercases (mokume:python/mokume/io/feature.py:44-60, applied at :113) [[MOKUME-13](findings-cross.md#mokume-13)].

**SPC-R5 (MUST)** A null `observed_mz`, `rt`, `ion_mobility` or `compensation_voltage` means not measured. In files of 1.3 and earlier read 0 there and charge ≤ 0 as missing, and treat one-element scans of `openms-consensus` psm files without `scan_format` as possible hashes.

## Open questions for maintainers
1. OpenMS ([T05](t05-feature-rows.md) OQ1): are ProteomicsLFQ PeptideIdentification RTs aligned when "Retention time alignment" is declared and `original_RT` absent? If so, `psm.rt` becomes null.
2. Is MaxQuant `Simple mass error [ppm]` the selected-ion deviation (848 ppm on an isotope pick)? It could fill MS:1000744.
3. CDAP: is the claimed ion `QueryPrecursorMz` or `OriginalPrecursorMz` [[MZID_CDAP-18](findings-psm.md#mzid_cdap-18)]?
4. Should feature carry `spectrum_reference`, so feature files can be `mixed`?
5. SPC-9 keys empty-scan PSMs by `spectrum_reference`, touching the key [#259](https://github.com/bigbio/qpx/issues/259) A left alone: accept, or drop such rows?
6. Should negative-ion data (signed `charge`) be in scope?
