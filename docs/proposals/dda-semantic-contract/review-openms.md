# Review: OpenMS as the primary producer

[Overview](index.md) · [Strict profile](strict-profile.md) · [Consumers](consumers.md) · [OpenMS work list](openms-worklist.md)

!!! note "Review report"
    One of three reviews of the [DDA semantic contract proposal](index.md), made by reading the OpenMS default branch at `0707508` and the #378 branch at `811c466`. It is kept as written by the reviewer. It predates the restructuring around the strict profile, so some of its recommendations are already applied. Rule ids refer to the topic sections T01 to T14.

## Summary

- **The native writer was obtained and read.** It is on the OpenMS default branch (HEAD `0707508`, 2026-11-09), read from a local clone. An unmerged branch, `claude/qpx-issue-377-review-dfv5gg` (`811c466`, 2026-11-08), holds the #378 work.
- **OpenMS already holds almost everything the contract asks for.** Of the rules that apply to OpenMS, about 30 are met at HEAD (several only up to stated caveats), 11 more only on the #378 branch, and about 65 need a small export or plumbing change. Only 3 are algorithmic for a default (QT-linker) run: rank-1 ties, a precursor-level FDR, and a transfer FDR for QT-linker MBR (PIP-ECHO already has one). Four items could not be classified from source.
- **The Python consensusXML route cannot reach the contract for three things**, because ProteomicsLFQ discards them before the consensusXML is written: unaligned RTs, feature RT bounds, and the measured precursor m/z.
- **qpx currently deprecates the native route** (`qpx:qpx/converters/openms/converter.py:402`, `docs/guide/convert.md:919`), and the worklist has 30 OpenMS items against the Python converter versus one for native. That is the opposite of the stated priority.
- **Two OpenMS facts in the proposal are wrong, and one open question is now answered** (details in Part 3).

Nothing was executed. Everything below is from reading source. No file in your repos was edited; the only writes were the clone and two temporary files in the scratchpad, since deleted.

Path aliases used below (all under `src/` in the OpenMS repository):

| Alias | File |
|---|---|
| CMAE | `openms/source/FORMAT/ConsensusMapArrowExport.cpp` |
| QF | `openms/source/FORMAT/QPXFile.cpp` |
| PGAE | `openms/source/FORMAT/ProteinGroupArrowExport.cpp` |
| AIH | `openms/source/FORMAT/ArrowIOHelpers.cpp` |
| ASR | `openms/source/FORMAT/ArrowSchemaRegistry.cpp` |
| LFQ | `topp/ProteomicsLFQ.cpp` |
| ISO | `topp/IsobaricWorkflow.cpp` |
| FFId | `openms/source/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.cpp` |
| ICR | `openms/source/ANALYSIS/ID/IDConflictResolverAlgorithm.cpp` |
| FDR | `openms/source/ANALYSIS/ID/FalseDiscoveryRate.cpp` |

---

## Part 1. The native writer

**Where it is.** `QPXFile`, `QPXDatasetExport`, `QPXCollectionExport`, `QPXIdentity`, `QPXValueValidation`, `ConsensusMapArrowExport`, `ProteinGroupArrowExport`, with schemas in ASR.

**Who calls it.**
- `ProteomicsLFQ -out_qpx` (LFQ:2844-2905), `IsobaricWorkflow -out_qpx` (ISO:1358-1437) and `ProteinQuantifier -out_qpx` write exactly three files: `quantms.feature.parquet`, `quantms.psm.parquet`, `quantms.pg.parquet`.
- ProSE writes psm, mz, run, sample and provenance through `QPXDatasetExport` (`topp/ProSE.cpp:2833-2842`).
- MS1LabeledWorkflow and OpenNuXL have QPX test fixtures; I did not read those tools.

**What it stamps.** `qpx_version = "1.1"` (AIH:184), `file_type`, `creator`, `software_provider`, `creation_date`, `compression_format`, `uuid`, `primary_key`, `identity_composite`, and `scan_format` on psm and feature. None of the five proposed keys. The LFQ and isobaric paths write no run, sample or provenance view.

**State at HEAD against the contract.**

| View | What HEAD writes |
|---|---|
| psm | One row per PeptideIdentification top hit (`export_all_psms=false`, QF:1083). No `rank`, `qvalue` or `spectrum_reference` column (ASR:478-512). **The PSM q-value is not written at all**: after `FalseDiscoveryRate` the main score is the q-value, and the writer exports the main score only when it is a PEP (QF:1228-1240). The fixture `tests/topp/IsobaricWorkflow_qpx_psm.tsv` confirms it. `observed_mz` NaN becomes 0.0 (QF:1263-1265). Score names are raw OpenMS keys (`Posterior Error Probability_score`). |
| feature | Peptide, PEP, `scan` and `id_run_file_name` come from the best hit across all runs of the consensus feature (CMAE:259-298, 970-988, 1094-1110), so a row can name a donor run. Unidentified features are rows with `peptidoform = ""` (CMAE:855-862). No `identification_type`, `match_*`, `consensus_rt`, `rt_aligned`. `rt_start`/`rt_stop` are `rt ± width/2` (CMAE:1255-1281). `pg_accessions` holds every peptide evidence, `unique` means one accession, `anchor_protein` is the first accession (CMAE:1113-1191). Intensity entries are non-nullable (ASR:620-621). |
| pg | Sound base: OpenMS's own abundance per (fraction group, label), identification-only rows with null label and intensity, `pg_id` over full membership (PGAE:800-899). But evidence counts every run the feature was *quantified* in, transfers included (PGAE:271-296); `is_decoy` comes from the anchor only (PGAE:633-642); `contaminant` is always null (PGAE:644-645). |

**The #378 branch adds** a per-run identifying PSM, own-run links, `identification_type` (direct, transferred, null), `match_qvalue` (float32, from PIP-ECHO) and `consensus_rt`. It is still 1.1, still drops the q-value, and still writes `""` peptidoforms.

---

## Part 2. Per topic: what OpenMS knows, and each rule's status

Legend:
- **T** = met at HEAD
- **B** = met only on branch `811c466`
- **S** = small OpenMS change (the data is in memory; export or plumbing only)
- **A** = algorithmic change
- **–** = not applicable to OpenMS
- **?** = could not be classified from source

| Topic | What OpenMS knows (class, field, source) | Rule status |
|---|---|---|
| **T01 PSM rows** | ProteomicsLFQ keeps only the best hit per identification (`IDFilter::keepBestPeptideHits(..., false)`, LFQ:994) and reduces duplicates of (spectrum, peptidoform, charge) to one (`reduceToOnePerSpectrum`, LFQ:1125). One engine is expected (LFQ:100-107). `PeptideHit::getRank()` exists (`PeptideHit.h:299`). | PSM-1 T · PSM-2 – · PSM-3/4 S (write constant 1) · PSM-5 S (record `rank1_only`) · PSM-6 **A** (ties are discarded upstream) · PSM-7 T · PSM-8 T · PSM-9 – · PSM-10 T (`QPXIdentity.cpp:285`) |
| **T02 Chimeric** | Several identifications on one spectrum with different peptidoforms are kept as separate measurements (LFQ:105-107), each with its own `getMZ()`. IsobaricWorkflow credits a spectrum's reporters to one peptidoform and moves the others to unassigned (`creditBestPeptidoform_`, ISO:602-634). | CHI-1 T · CHI-3/5/6/10 S (group by run and spectrum reference, index by m/z) · CHI-4 T · CHI-7 T · CHI-8 T · CHI-9 S · CHI-11 T |
| **T03 Filtering** | `-psmFDR` (default 1.0, no cut), `-FDR_type` (PSM or PSM+peptide), `-proteinFDR` (0.05 for LFQ, 1.0 isobaric), `picked_proteinFDR`. Decoys are removed unless `PeptideQuantification:quantify_decoys` (LFQ:2334-2340). Unassigned identifications are always exported (LFQ:2872-2884). The decoy affix and position are known (LFQ:978-980). | FLT-1/2/3 S (caveat: cuts applied before ProteomicsLFQ, as in quantms, need a parameter) · FLT-4 S · FLT-6 S · FLT-7/8/9 S (or leave `pass_threshold` out) · FLT-10 S · FLT-11 S (new contaminant-affix parameter) · FLT-12 T · FLT-5/13/14 – |
| **T04 Confidence** | A PEP is required on input (LFQ:891-905). After FDR the score type is `q-value` (FDR:985), or `peptide q-value` under PSM+peptide (FDR:1194), and the PEP is kept as a meta value. The PSM FDR is computed once over the whole consensus map (`applyBasic(consensus, true)`, LFQ:2331; `treat_runs_separately` defaults to false, FDR:38). OpenMS has no precursor-level FDR. | CNF-1/2 S (export the q-value) · CNF-3 S (scope is always `experiment`) · CNF-4 – · CNF-5/6 S (snake_case; ProSE maps `q-value` to `global_qvalue`, `QPXDatasetExport.cpp:175`, which conflicts with CNF-5) · CNF-7 T · CNF-8 B · CNF-9 S (float64) · CNF-10 T · `precursor_qvalue`/`global_qvalue` **A**, so always null |
| **T05 Feature rows** | A `FeatureHandle` carries rt, mz, intensity, charge, width and unique id per run. The per-run `Feature` has convex hulls, isotope sub-features and meta values (`model_status`, `raw_intensity`, `masserror_ppm`, FAIMS CV, IM statistics). ProteomicsLFQ strips the meta values and clears the hulls (LFQ:1926-1966), and aligns RTs in place without storing the original (LFQ:733; `MapAlignmentTransformer.cpp:89-100`, `store_original_rt` defaults to false). | FRW-1 T · FRW-2 T (LFQ:1975; duplicate keys refused) · FRW-3/4 S (the composite includes `scan`, `QPXIdentity.cpp:284`) · FRW-5 S · FRW-6 – · FRW-7/8 S (carry the original RT per handle, or invert the transformation) · FRW-9 S (keep hull bounds) · FRW-10 ? · FRW-11 T |
| **T06 Identifications per feature** | `IDConflictResolverAlgorithm::resolve` keeps one identification per consensus feature and moves the rest to unassigned with a `feature_id` meta value holding the feature's unique id (ICR:311-357). FFId computes `n_total_ids`, `n_matching_ids` and `feature_class` (FFId:1660-1668), later stripped. | FID-1/2/4/5/6 B · FID-3 S (the branch compares the main score, which is the heavily tied q-value after FDR; PEP should come first) · FID-7/8/9 B · FID-10 S · FID-11 S (`getUniqueId()`) · FID-12 S |
| **T07 Unidentified** | Seed features from FeatureFinderMultiplex or Biosaur2, the `targeted_only` setting, handle unique ids. | UNI-1..5 S (move hitless rows to the new view) · UNI-6 S (`removed` is true today; `retained` needs identification-only rows built from unassigned identifications) · UNI-7 S |
| **T08 MBR** | **`mbr_transfer_map_indices`** on each ConsensusFeature lists the runs its identification was transferred to. It is written after FDR filtering (ICR:30-68, LFQ:2781-2784) or by PIP-ECHO (`PIPECHO/Impl.cpp:898`). **`mbr_transfer_qvalues`** holds PIP-ECHO's per-transfer q-values (`Constants.h:463-469`). Also known: `pip_echo`, `PipEcho:fdr` (0.05), `targeted_only`. Neither the HEAD writer nor qpx reads these meta values. | MBR-1/2/3/5 B (the branch infers from missing own-run evidence; it should use the annotation) · MBR-4 `match_qvalue` B, score and deltas ? · MBR-6/7 T · MBR-8/9 S (QT-linker transfers have `mbr_fdr = none`) · MBR-10 S · MBR-11 ? |
| **T09 Intensities** | LFQ intensity is the elution-model area when a model is fitted. The raw value is kept as `raw_intensity`, and failed fits are replaced by a regression estimate (`ElutionModelFitter.cpp:410-449`). ProteomicsLFQ applies no normalization (LFQ:215-229). Isobaric: isotope-correction and normalization flags (ISO:732-733), reporters below `min_reporter_intensity` set to 0 (ISO:490-493), the quantifying scan `scan_id` (ISO:518), `precursor_purity` and `quant_precursor_purity` (ISO:508-515). Protein quantity: `method`, `top:N`, `top:aggregate`, `consensus:normalize` (`PeptideAndProteinQuant.cpp:43-66`). | INT-1 T · INT-2 S · INT-4 T · INT-5/6 S (nullable entries) · INT-7 S (those zeros become null) · INT-8 T · INT-9 S (and a decision on regression-estimated fits) · INT-10 S · INT-11 S · INT-12 S |
| **T10 Protein mapping** | `PeptideEvidence` from PeptideIndexing (LFQ:964-973). `target_decoy` on hits, including `target+decoy`. One merged identification run and one inference (`mergeAllIDRuns`, LFQ:2265-2268) giving `IndistinguishableProteins`. The quantification mode (unique, strictly unique, shared with greedy resolution). | PRT-1 T · PRT-2 S · PRT-3 T · PRT-4 T (but defaults to false when the meta value is absent, QF:1244-1252) · PRT-5 S · PRT-6 S · PRT-7 S · PRT-8 ? (greedy assignment) · PRT-9/10 S · PRT-11 S · PRT-12 T once PRT-7 is fixed, and it can be MUST · PRT-13 S |
| **T11 Spectrum fields** | `getMZ()`, `getRT()`, native id, hit charge, IM and FAIMS CV meta values (FFId:602-623). **The LFQ feature m/z is theoretical** (`feat.setMZ(PrecursorMZ)`, FFId:1587); the measured deviation is `masserror_ppm`, kept only with PIP-ECHO (LFQ:1936-1937). `scan` is one integer or `[]`, never a hash (AIH:516-533). `scan_format` is `index`, `scan` or omitted (AIH:209-214). | SPC-1 S · SPC-2 S · SPC-3/4 T · SPC-5 S (T on the ProSE path) · SPC-6 S · SPC-7 S (no `nativeId`, no multi-component ids) · SPC-8 T · SPC-9 S · SPC-11 **S** (`psm.rt` is aligned) · SPC-12 S (CV column) · SPC-13 T (collisions only warn) · SPC-14/15 T |
| **T12 Protein groups** | The pg unit is the fraction group, from the experimental design. Group probability is reported as a q-value only when a protein FDR ran (PGAE:221-224, 589-592). `pg_qvalue` is always null (PGAE:594-595). | PGV-1 T (plus S for `identified_in_unit` and labels on identification-only rows) · PGV-2 S (exclude transfer runs) · PGV-3/4 S · PGV-5 T · PGV-6 S (always `dataset_inference`) · PGV-7 S · PGV-8 T · PGV-10 S · PGV-11 T · PGV-12 S (`pg_names` is `[]`) |
| **T13 Metadata** | Every value of the five new keys is a tool parameter or a constant. `provenanceView` and `paramToJSON` exist but only ProSE uses them. `QPXValueValidation` already refuses invalid tables at write time. | META-1/2/3/5 S · META-9 T · META-10 S (1.1 to 1.4) · META-11 S |
| **T14 Reader contract** | – | A strict profile lets readers drop their unknown and unassessable branches for OpenMS data. |

---

## Part 3. Where the proposal bends for other producers, and an OpenMS-first profile

### Corrections to the proposal's OpenMS facts

1. **MBR-2 and FID-1 cite `FFId_category` as OpenMS's per-row statement.** At HEAD FFId only ever sets `"internal"` (FFId:292, 367, 746), and ProteomicsLFQ strips feature meta values. The real statement is `mbr_transfer_map_indices`, with `mbr_transfer_qvalues` under PIP-ECHO.
2. **PGV-6 assigns `partial_inference` to OpenMS.** ProteomicsLFQ always merges to one identification run before inference, and the exporter requires that single run. Native output is always `dataset_inference`.
3. **The open question "are ProteomicsLFQ PeptideIdentification RTs aligned?" is answered: yes.** Feature RTs, the RTs of their identifications and of unassigned identifications are all transformed in place, and `original_RT` is not stored (`MapAlignmentTransformer.cpp:81-100`, LFQ:733).
4. **The "not verified" C++ rows in T01, T02 and T09 can now be filled in from Parts 1 and 2.**

### Complexity OpenMS and its consumers do not need

| Bend | For whom | Under an OpenMS-first priority |
|---|---|---|
| `unknown` on every footer key; null meaning unknown on `rank`, `precursor_index`, `identification_type`, `identified_in_unit`, `imputed_labels`, `group_unique`, `contaminant` (META-1, PSM-4, CHI-10, MBR-8, INT-9, PRT-9) | Converters that cannot see the source | Forbidden in the strict profile. OpenMS knows every one. |
| Multi-run unit rows, `quantification_unit_id`, feature `grouped_runs` (FRW-3/4/6, FID-1, CNF-1) | FragPipe | Move to a converter annex. |
| Grain inferred by parsing `identity_composite` into four cases (FRW-3) | MaxQuant, PD, FragPipe | One fixed composite per OpenMS mode. |
| `fdr_filter` values `applied` and `unknown`; `pass_threshold` with provenance level and FDR (FLT-7..9) | MaxQuant, Mascot, PD, mzIdentML | Decimals or `none` only. `pass_threshold` is not needed because q-values are typed. |
| Declarations moved to provenance by R1 (`decoy_policy`, `psm_qvalue_scope`, `pg_grouping`, `mbr_fdr`, `rank1_only`, and others) | Converters that cannot know them | See the conflict below. |
| Multi-engine merge, engine prefixes, `consensus_support` (PSM-2); selection of mzIdentML lists (PSM-9); alternative groups and ties (PSM-6) | mzIdentML, PD, MaxQuant | Not applicable. |
| The heuristic fallback for `transferred` (MBR-2); the 10 ppm fallback for `precursor_index` (CHI-5) | Converters | The annotation, not a heuristic. |
| "Least-processed quantity the input exports" with `ms1_trigger`, `reporter_sn`, `fragment`, `maxlfq`, `other` (INT-1/2) | FragPipe, PD, DIA-NN | The primary is always the measured signal. |
| `MULTI-MATCH-MSMS`, the Match-score rename, synthetic windows, combined spectra (SPC-10), `Quan Usage`, `mapping_set`, `skipped_runs_psm` | MaxQuant, PD, CDAP, FragPipe | Annex. |

**Conflict with revision R1, for the maintainers to decide.** R1 moved many declarations to provenance. The native LFQ and isobaric paths write no provenance view, and two reader defaults depend on such values: pivoting protein groups across units needs `pg_grouping`, and FDR re-estimation needs `decoy_policy`. For OpenMS-written files, either the profile makes a provenance view with fixed step and parameter names mandatory, or one footer claim is reinstated (for example `conformance_profile = dda-strict`), which reverses R1's removal of that key. A single claim is what lets qpxui, the portal and mokume branch once instead of per key. I would reinstate it; staying inside R1 costs a mandatory provenance view and a reader-side lookup.

### Proposed OpenMS reference profile (`dda-strict`)

Every item is known to OpenMS today. No value may be `unknown`.

**Dataset.** psm, feature, pg, run, sample and provenance always. `unidentified_feature` if and only if `targeted_only = false`.

**Footer.**

| Key | Allowed values |
|---|---|
| `qpx_version` | `1.4` |
| `acquisition_method` | `dda` |
| `fdr_filter` | `none`, or decimals at `psm`, `peptide`, `protein` |
| `unquantified_policy` | `retained` |
| `intensity_quantity` | feature: `ms1_area`, `reporter_corrected` or `reporter_uncorrected` |
| `intensity_processing` | feature: `none` or `normalized`; pg: `topn`, `sum` or `normalized` |
| `scan_format` | `scan`, `index` or `nativeId` |
| `identity_composite` | LFQ: `peptidoform, charge, run_file_name, rt`. Isobaric: `peptidoform, charge, run_file_name, scan`. |

**Fixed provenance parameters (or footer keys, per the conflict above).** `decoy_policy`, `psm_qvalue_scope = experiment`, `pg_grouping = dataset_inference`, `match_between_runs`, `mbr_fdr` (a decimal under PIP-ECHO, `none` for the QT linker), `rank1_only = true`.

**psm, never null.** `rank` (always 1), `precursor_index`, `qvalue`, `posterior_error_probability`, `is_decoy`, `observed_mz`, `rt` on the acquisition clock, `scan` with `spectrum_reference`, `protein_accessions`, `contaminant`. `feature_id` is null only for unassigned PSMs.

**feature.**
- `identification_type` is `direct`, `secondary` or `transferred`, never null.
- Direct rows: exactly one `scan`, `id_run_file_name = run_file_name`, PEP and `psm_qvalue` non-null.
- Transferred rows: none of those; `match_qvalue` non-null if and only if PIP-ECHO ran.
- `identification_conflict` and `consensus_feature_id` non-null.
- `rt`, `rt_aligned`, `consensus_rt`, `rt_start`, `rt_stop` non-null for LFQ, with `rt_start ≤ rt ≤ rt_stop`.
- `observed_mz` is measured or null, never theoretical.
- One intensity entry per label; `imputed_labels = []`.
- `pg_accessions` is the full membership of one group, with `group_unique` and the group's anchor.
- `precursor_qvalue` and `global_qvalue` are null by definition.
- Isobaric rows: `quant_scan`, `quant_ms_level`, `parent_ion_fraction`, `quant_spectrum_count = 1`.

**pg.** `identified_in_unit` non-null. `is_decoy` true only if every member is a decoy. `contaminant` non-null. `global_qvalue` non-null when a protein FDR was applied. `pg_qvalue` null by definition. `intensity` is OpenMS's own abundance.

### What to tighten, simplify or drop

1. **Add `dda-strict` as a normative profile now.** The overview defers rejecting `unknown` "until every qpx converter fills the keys"; split the profiles instead, so the lenient one serves converters.
2. **Move the FragPipe, MaxQuant, PD, mzIdentML and CDAP clauses to a per-converter annex**, so the core rules read as the OpenMS behaviour plus stated relaxations.
3. **Make PRT-12 and PGV-6 MUST for OpenMS**; drop `partial_inference` for it.
4. **Fix the OpenMS facts** listed above.
5. **Define FID-3 for OpenMS concretely**: the resolver's kept identification in its own run; in other runs the best removed identification of the same peptidoform and charge, ordered by PEP then scan.
6. **Decide two OpenMS-specific points.**
   - Whether a regression-estimated elution-model area counts as imputed.
   - Whether OpenMS should drop a transferred row when the same precursor has a direct row in that run. That would give at most one row per (peptidoform, charge, run), and readers would need no de-duplication.
7. **Re-target the worklist** from the Python converter to the C++ writer, and reverse the deprecation notice.

---

## Part 4. OpenMS-side changes in priority order, and what they retire in qpx

| # | OpenMS change | Consumer benefit | qpx code no longer needed |
|---|---|---|---|
| 1 | Merge the #378 branch: per-run identifying PSM, own-run links, `identification_type`, `match_qvalue`, `consensus_rt`. Source `transferred` from `mbr_transfer_map_indices`, order by PEP, write `match_qvalue` as float64. | qpxui and mokume can separate transfers from identifications (mokume#117). Fixes the donor run, scan and PEP on non-winner rows. | The recovery of resolver-removed identifications and the transfer heuristic in `openms_consensus/feature_adapter.py:523-634`; `record_linking.py`; `openms/psm_runs.py` |
| 2 | Export the PSM q-value to `psm.qvalue` and `feature.psm_qvalue` (and `peptide_qvalue` under PSM+peptide), with snake_case names. | `confident()` works; no native file has a typed q-value today. | `psm_adapter.py:204-246`; `feature_adapter.py:109-206` |
| 3 | Clocks and bounds: keep the original RT and the hull bounds per handle; write `rt` and `psm.rt` on the acquisition clock, plus `rt_aligned`. | XIC and mz joins in qpxui; `rt` inside its window. | Nothing to retire: this cannot be done from consensusXML at all. |
| 4 | Stamp the five footer keys and `qpx_version` 1.4. Write run, sample and provenance from the experimental design and tool parameters, reusing `QPXDatasetExport` as ProSE does. | The portal reads the filter state and intensity kind without special cases. | Most of `openms/converter.py` (label relabelling, run-name repair, the provenance stub at :544-574). Only SDRF sample metadata remains. |
| 5 | Stop writing sentinels and placeholders: unidentified features go to `unidentified_feature`; charge ≥ 1; `observed_mz` nullable (measured from `masserror_ppm`, else null; NaN no longer 0); intensity entries nullable; sub-threshold reporters and failed fits null instead of 0. | No `""` peptidoforms or fake zeros in matrices. | The migration in UNI-8; the handling of OPENMS-13, -15, -16, -M5 |
| 6 | Protein fields: feature `pg_accessions` as group membership, `group_unique`, the group's anchor, `contaminant` (new affix parameter). pg: `is_decoy` over all members, `identified_in_unit`, evidence only from own-run identifications. | mokume's group-specific filter is correct; transfer-only protein rows are visible. | `feature_adapter.py:455-490`; the evidence and decoy logic in `pg_adapter.py` |
| 7 | On the consensus path: `rank = 1`, `precursor_index`, `spectrum_reference`, `isotope_error`, `compensation_voltage`; `scan_format` `nativeId` with multi-component ids. | Chimeric counting; timsTOF and FAIMS data. | `psm_adapter.py:49-71` (hash scans) |
| 8 | Isobaric rows: `quant_scan`, `quant_ms_level`, `parent_ion_fraction`, `quant_spectrum_count = 1`. | Purity filtering; reporter de-duplication. | `quant_scan.py` |
| 9 | `identification_conflict`, `consensus_feature_id`, `spectral_count`. | Traceability across runs. | `feature_dedup.py` |
| 10 | Identification-only rows from unassigned identifications, with `unquantified_policy = retained`. | "Identified but not quantified" is distinguishable from "not identified". | The handling of OPENMS-M2 |
| 11 | Algorithmic, later: transfer score and RT/m/z deltas from PIP-ECHO; an FDR for QT-linker transfers; precursor-level q-values; keeping rank-1 ties. | Fills `match_*` and `precursor_qvalue`. | – |

**Net effect on qpx.** `qpx/converters/openms_consensus/` (about 4,100 lines) becomes a legacy path for existing consensusXML files. `qpx/converters/openms/` (about 750 lines) shrinks to validating and adding SDRF sample metadata, ontology, dataset and pepmap. `qpx/converters/quantms_msstats/` (about 710 lines) is not needed for OpenMS-produced data.

---

## Not verified from source

- Whether a protein group's `probability` after protein FDR is a group-level q-value. The writer assumes so (PGAE:221-224); I did not read that part of `FalseDiscoveryRate`.
- `FFId_category = external`: not found in the FEATUREFINDER or ANALYSIS/ID sources at HEAD. It may exist elsewhere.
- PIP-ECHO internals: whether a transfer score or RT/m/z deltas are stored, and whether decoy transfers can reach the output.
- MS1-label handling (MS1LabeledWorkflow, FeatureFinderMultiplex, `MS1LabelState`), so FRW-10.
- How greedy group resolution records which group a shared peptide is assigned to (PRT-8).
- Whether search adapters set the identification's m/z to the co-isolated ion for chimeric hits (CHI-7), and whether the search engine or ProteomicsLFQ recalibrated it (SPC-2).
- The `id_scan_id` meta value that `qpx/converters/openms_consensus/quant_scan.py:21` relies on: not found in ISO at HEAD.
- Whether seed features can carry charge 0 at HEAD, and the units of the IM meta values.
- How cuts applied before ProteomicsLFQ in the quantms pipeline are recorded.
- The ProSE writer was read only through `QPXDatasetExport.h`.
- All findings on feature-row values come from reading the code, not from a converted native fixture; qpx has none.
- The T14 reader rules were only skimmed.
