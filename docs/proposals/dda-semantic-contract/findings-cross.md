# Findings: cross-view, footer and documentation

Part of the [DDA semantic contract proposal](index.md). Each entry is one verified finding from the audit of qpx `dev` at commit `32cbbdc` (specification 1.3). A skeptic re-checked every finding against the code; "partly confirmed" entries are true only as stated in their correction. `audit-run:` evidence refers to converting the repository's own fixtures (`tests/examples`) with that commit during the audit; those outputs are not committed.

116 findings, grouped by audit area.

## OpenMS converters

### OPENMS-21: --no-include-unassigned-psms leaves feature evidence pointing at unexported PSMs {#openms-21}

*cross* · ambiguous · confirmed · checklist P10, X1

- **Now:** Removed (unassigned) same-peptide PIDs always feed the feature's scan, id_run_file_name, PEP and identification_type (record_linking.py:224-236). They reach the PSM view only when the option is on (172-182). Test lfq_noun: 25 of 68 'direct' feature rows have no linked PSM, and their scans are absent from psm.parquet.
- **Spec:** psm.md:170-173: the option only excludes unassigned IDs. Nothing says whether feature.scan or identification_type must be backed by an exported PSM.
- **Reader impact:** Joining feature.scan to psm fails for 'direct' rows. The same feature row means different things depending on the option.
- **Suggested:** State that feature identification fields reference PSMs in the psm view, or that they are independent. Keep the identifying PSM in the view regardless of the option.
- **Evidence:** `qpx/converters/openms_consensus/record_linking.py:172-196,219-245`; `docs/spec/psm.md:170-173`; `test: audit-run:openms/lfq_noun (25/68 direct features without PSM)`

### OPENMS-22: Decoys kept everywhere; is_decoy derived differently in psm/feature vs pg {#openms-22}

*cross* · ambiguous · confirmed · checklist P4, F8, G3

- **Now:** No view drops decoys. PSM and feature is_decoy come only from target_decoy meta == 'decoy' (feature_adapter.py:125-132), and are false when the meta is absent. pg falls back to accession prefixes DECOY/REV_/RANDOM_ (pg_adapter.py:154-155,394); suffix markers such as the LFQ fixture's PeptideIndexer '_rev' are not recognised. Decoy features feed the pg roll-up and peptide counts.
- **Spec:** producer-guide.md:30-34: decoys keep markers and is_decoy is true. pg.md:49: do not parse the prefix. pg.md:99-103: decoy groups may be present. The feature-view decoy policy was left open (#229).
- **Reader impact:** Without target_decoy meta, decoy PSMs and features show is_decoy=false while their pg rows show true. Feature-level decoy presence differs by producer.
- **Suggested:** One decoy derivation for all views, from meta or a declared decoy string with its position. The contract states whether decoy features and PSMs are present (footer).
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:125-132`; `qpx/converters/openms_consensus/pg_adapter.py:154-155,394`; `tests/examples/openms_topp/ProteomicsLFQ_1_out.consensusXML (PeptideIndexer:decoy_string _rev suffix)`; `docs/spec/producer-guide.md:30-34`; `docs/spec/pg.md:49,99-103`

### OPENMS-30: No footer/provenance declares FDR level, threshold, rank policy, q-value level {#openms-30}

*footer* · unspecified · partly confirmed (see correction) · checklist D8, P3, D4

- **Now:** PSM and feature footers carry only the standard keys plus identity_composite and identity_hash_version, and scan_format when the format is uniform (LFQ output: none). Provenance lists the search engines and dataProcessing actions only (provenance.py:69-180). The PSM- versus peptide-level FDR flag (Percolator:peptide_level_fdrs) is used internally but not recorded.
- **Spec:** psm.md:256-268 and feature.md:272-284 list footer keys: none for filtering, rank or q-value level. producer-guide.md:23: set scan_format.
- **Correction:** There is no footer key; confirmed. But the spec has a provenance mechanism: parameters answer 'what FDR was used?', and conversion options include qvalue_threshold ('none' when unset) and rank1_only (provenance.md:58-65). OpenMS records only include_unassigned_psms. The gap is that none of this is required, and SearchParameters flags such as Percolator:peptide_level_fdrs are not copied.
- **Reader impact:** Readers cannot tell how a file was filtered or what its q-values mean without re-reading the source.
- **Suggested:** Add normative footer keys such as `fdr_filter` (none|psm|peptide|protein with threshold), `qvalue_level` per column, `psm_ranks` and `feature_grain`, and require them in the DDA contract.
- **Evidence:** `qpx/converters/openms_consensus/provenance.py:69-180`; `qpx/converters/openms_consensus/feature_adapter.py:135-146`; `docs/spec/psm.md:256-268`; `docs/spec/feature.md:272-284`; `audit-run:openms/*/openms.*.parquet footers`

### OPENMS-31: Native -out_qpx passthrough: no semantic normalisation; two OpenMS routes differ {#openms-31}

*cross* · contradictory · partly confirmed (see correction) · checklist X1, D8

- **Now:** qpxc convert openms (deprecated) relabels intensities, normalises run names and accessions, optionally recovers PSM runs, and preserves the supplied ids and duplicate rows (converter.py:209-224,261-302; cli/convert/openms.py:16,61-69). Typed q-values are derived from whatever score names OpenMS wrote. The default feature composite is (run,peptidoform,charge,rt,scan,observed_mz) (converter.py:50-53), against 7 columns in the consensus route. The semantics of OpenMS native data could not be verified: there is no fixture.
- **Spec:** identity.md:38: 'qpxc convert openms re-derives ids'. feature.md:14-17 and convert.md:947: supplied ids are preserved (the code keeps them: override_provided_ids=False). Nothing defines the row semantics of native OpenMS files.
- **Correction:** Both routes keep supplied ids and derive only missing ones (writers/base.py:244; feature override_provided_ids=False at converter.py:216-217). The 6-column composite is only a fallback when the source footer declares none (converter.py:182-188). identity.md:38 itself says the two routes' feature_ids differ, but its 're-derives ids' contradicts feature.md:14-17 and convert.md:946-949.
- **Reader impact:** The same OpenMS run gives different feature ids, rows and confidence semantics depending on the route. The identity documentation contradicts itself.
- **Suggested:** Apply one DDA contract to both routes, or mark native files as non-conforming. Fix identity.md:38.
- **Evidence:** `qpx/converters/openms/converter.py:50-53,209-224,261-302`; `qpx/cli/convert/openms.py:16,59-69`; `docs/spec/identity.md:35-38`; `docs/spec/feature.md:14-17`; `docs/guide/convert.md:945-951`

### OPENMS-M5: OpenMS routes disagree on unidentified features {#openms-m5}

*cross* · inconsistent across producers · found by the verifier · checklist F3, X1, D1

- **Now:** openms-consensus drops features without hits (feature_adapter.py:714-716; psm.md:175). Per #229, OpenMS -out_qpx feature files contain unidentified features (PXD000279), and `qpxc convert openms` preserves rows without filtering (convert.md:946-949; openms/converter.py:261-302). Not verified: there is no native fixture.
- **Spec:** feature.md:7-10: unidentified features need a producer id. feature.md:51-52: peptidoform and sequence are required. No identification_type value covers them (feature.yaml:210-214).
- **Reader impact:** The same OpenMS analysis yields feature views with or without unidentified rows depending on the route. Their required peptidoform fields hold placeholders or nothing.
- **Suggested:** One rule for both routes: drop unidentified features, or allow them with identification_type='unidentified' and a nullable peptidoform. Declare the choice in the footer.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:714-716`; `qpx/converters/openms/converter.py:261-302`; `docs/guide/convert.md:946-949`; `GitHub bigbio/qpx #229 body ('distinct unidentified Features in PXD000279')`; `docs/spec/feature.md:7-10,51-52`

## MaxQuant converter

### MAXQUANT-06: Provenance and footer omit FDR thresholds, MBR setting and MaxQuant version {#maxquant-06}

*cross* · diverges from spec · partly confirmed (see correction) · checklist D8, P3, X1

- **Now:** The MaxQuant provenance step has only isotope_corrected=false and normalization=none, and tool_version is null (converter.py:311-336, checked in output). parameters.txt is read only for 'Fixed modifications' (constants.py:137-144), not for PSM FDR, Protein FDR, Match between runs or version. The psm, feature and pg footers have no filtering-state key.
- **Spec:** provenance.md:3 and :59: provenance answers 'what FDR was used?'. provenance.md:65: the threshold is on record even when it is 'none'. Examples use psm_fdr/protein_fdr (provenance.md:160-162).
- **Correction:** Facts hold, but the line refs are wrong. provenance.md:3 is aspirational ('should be able to answer'), with no MUST. provenance.md:65 covers qpxc's own options (qvalue_threshold), not the tool's FDR. This is a gap ('unspecified'), not a divergence. Sidecar parameters.txt is read but not recorded (MAXQUANT-M3).
- **Reader impact:** A reader cannot learn the PSM or protein FDR, whether MBR was on, or the MaxQuant version.
- **Suggested:** Parse parameters.txt/summary.txt into provenance (psm_fdr, protein_fdr, match_between_runs, version). The contract requires per-level filter state.
- **Evidence:** `qpx/converters/maxquant/converter.py:309-336`; `qpx/converters/maxquant/constants.py:137-144`; `docs/spec/provenance.md:3,59,65,160-162`; `audit-run:maxquant/simple/maxquant.provenance.parquet`

### MAXQUANT-14: Primary reporter intensity is uncorrected in feature/pg but corrected in site view {#maxquant-14}

*cross* · inconsistent across producers · partly confirmed (see correction) · checklist F6, G2, X1

- **Now:** Feature and pg use the uncorrected 'Reporter intensity N' as primary and put the corrected value in additional 'reporter_intensity_corrected' (feature_adapter.py:539-553; pg_adapter.py:516-536). Provenance says isotope_corrected=false (converter.py:320). The site view uses 'Reporter intensity corrected N' as primary (sites.py:19,111; site.md:112). The name `reporter_intensity` (pg.md:241) is never emitted.
- **Spec:** intensities.md:211-213 states the processing once per tool step, implying the same for all views. pg.md:241-242.
- **Correction:** Facts hold. The provenance key is at converter.py:122, not :319-320. This is a contradiction inside one producer (views of one MaxQuant step), not across producers. isotope_corrected=false is wrong for the site view, which is in that step's output_views.
- **Reader impact:** In one dataset, site intensities are isotope-corrected and feature/pg ones are not, while provenance says false.
- **Suggested:** One rule for primary reporter intensities in all views, or a per-view processing state.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:539-553`; `qpx/converters/maxquant/pg_adapter.py:516-536`; `qpx/converters/maxquant/sites.py:19,111`; `qpx/converters/maxquant/converter.py:319-320`; `docs/spec/site.md:112`; `docs/spec/pg.md:241`

### MAXQUANT-21: Every MaxQuant 0 becomes null, including reporter channels of quantified MS/MS {#maxquant-21}

*cross* · ambiguous · partly confirmed (see correction) · checklist D6, F6, G1

- **Now:** Zeros become null for LFQ (feature_adapter.py:526), TMT (feature_adapter.py:543), SILAC (silac.py:82) and pg (pg_adapter.py:39-42). Null pg rows are kept only with MS/MS evidence (pg_adapter.py:518,571,589). A TMT channel at 0 in a quantified MS/MS also becomes null. MaxQuant output never contains 0.
- **Spec:** intensities.md:126-127: a tool's 'not quantified' sentinel 0 becomes null, but a measured 0 stays. pg.md:137-140.
- **Correction:** Zero-to-null is real, but the spec does decide the reporter case: intensities.md:123-124 names 'an empty reporter channel of an identified plex' as a null entry. The status is therefore specified. What remains is a design concern: below-detection and not-quantified look the same.
- **Reader impact:** A reporter channel truly undetected in a quantified scan looks the same as 'not quantified', which matters for missing-value models.
- **Suggested:** The contract states per producer which zeros are sentinels. For reporter ions, decide whether 0 means below detection (keep) or null.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:526,543`; `qpx/converters/maxquant/silac.py:82`; `qpx/converters/maxquant/pg_adapter.py:39-42,516-520`; `docs/spec/intensities.md:121-129`

### MAXQUANT-22: Decoys kept everywhere; contaminant flag lost on psm/feature; decoy PSMs lack proteins {#maxquant-22}

*cross* · unspecified · partly confirmed (see correction) · checklist P4, P8, F8, G3

- **Now:** Decoys are kept in all views with is_decoy from 'Reverse' and the REV__ marker. PSM decoys have null protein_accessions (3/3 in simple; msms 'Proteins' is blank). 'Potential contaminant' is read for pg only (pg_adapter.py:460). In the simple fixture, 23 contaminant evidence rows have a target razor protein (e.g. 'P04264;CON__P04264' with razor P04264), so their features carry no contaminant signal.
- **Spec:** producer-guide.md:28-35: contaminant is flagged on pg only. pg.md:47-50: do not parse status from the prefix.
- **Correction:** Facts hold, but this is not unspecified. producer-guide.md:33 deliberately puts `contaminant` on pg only, and converter-coverage.md:281 documents the blank decoy Proteins. What remains is a design limitation (a join is required).
- **Reader impact:** Contaminant filtering at PSM and feature level needs a pg join, and prefix heuristics are wrong.
- **Suggested:** Add a contaminant bool to psm and feature, or define the join rule in the contract.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:448-449`; `qpx/converters/maxquant/psm_adapter.py:200-211`; `qpx/converters/maxquant/pg_adapter.py:459-460`; `docs/spec/producer-guide.md:28-35`; `tests/examples/maxquant/maxquant_simple/evidence.txt col 80`

### MAXQUANT-24: Canonical peptidoform not guaranteed: missing fixed mods, named unknown mods, protein-terminal fixed mods {#maxquant-24}

*cross* · diverges from spec · partly confirmed (see correction) · checklist D5, P9

- **Now:** Fixed mods come from the SDRF, else parameters.txt, else summary.txt (base_adapter.py:176-206). If none is found, the converter only warns. The simple golden fixture has none, so 1,497 Cys features lack UNIMOD:4 with a median |error| near 32,000 ppm (validate warns on 1,485 rows). Unknown mod names stay names (ptm.py:338-339), e.g. PEPT[MyCustomMod (T)]IDEK. Protein N-term fixed mods are skipped silently from parameters.txt (constants.py:129-133).
- **Spec:** peptidoform.md:17-30: all modifications inline, and a mod without Unimod entry is a signed 4-decimal mass. producer-guide.md:52. versioning.md 1.3 'Canonical peptidoform'.
- **Correction:** Keeping unknown names is allowed: peptidoform.md:43-44 says a name with no Unimod entry and no mass 'is kept as written'. The missing fixed mods (only a warning) and the skipped protein-terminal fixed mods do violate rule 1.
- **Reader impact:** The same molecule gets different peptidoform and feature_id across producers, with a wrong calculated_mz.
- **Suggested:** Fail, or flag non-canonical output in the footer, when fixed mods are unknown. Write mass shifts for unknown names.
- **Evidence:** `qpx/converters/maxquant/base_adapter.py:176-206`; `qpx/converters/maxquant/constants.py:120-134`; `qpx/converters/ptm.py:321-339`; `docs/spec/peptidoform.md:17-30`; `qpxc validate --semantic on audit-run:maxquant/simple`

### MAXQUANT-27: No scan_format footer; MaxQuant 1.5 scan header not mapped {#maxquant-27}

*cross* · diverges from spec · partly confirmed (see correction) · checklist P7, F10, D8

- **Now:** PsmWriter and FeatureWriter get no scan_format (psm_adapter.py:117; feature_adapter.py:181-185), and output footers lack it. The feature scan maps only 'MS/MS scan number' (column_mappings.yaml:159), so MQ 1.5 'MS/MS Scan Number' leaves feature.scan empty (#378). timsTOF 'Scan number' semantics not verified.
- **Spec:** producer-guide.md:23: set scan_format. scan.md:80 and versioning.md 1.3: non-scan values must be declared.
- **Correction:** A missing footer is compliant for Thermo scans: scan.md:80-81 reads a file without scan_format as scan numbers, and producer-guide.md:23 is advice. The MQ 1.5 'MS/MS Scan Number' miss is plausible (exact-case matching) but comes from open #378. timsTOF is unverified.
- **Reader impact:** The validator cannot check scans, and MQ 1.5 features lose their identifying scan.
- **Suggested:** Declare scan_format=scan for Thermo input and map case variants of the header.
- **Evidence:** `qpx/converters/maxquant/psm_adapter.py:117-125`; `qpx/converters/maxquant/feature_adapter.py:181-185`; `qpx/converters/column_mappings.yaml:159,182`; `docs/spec/producer-guide.md:23`; `docs/spec/scan.md:80`; `github #378`

### MAXQUANT-M1: NULL Modified sequence becomes literal peptidoform 'None' in PSM and feature {#maxquant-m1}

*cross* · diverges from spec · found by the verifier · checklist P9, D5, F3

- **Now:** str(row.get('Modified sequence')) turns a NULL cell into 'None', and to_proforma('None') returns 'None' (psm_adapter.py:259; feature_adapter.py:441). The PSM skip (psm_adapter.py:266-270) never fires. Verified: work/verify-maxquant/psmnull (5 rows → 5 PSMs, one 'None') and types/out (feature 'None' with a real sequence and intensity). `qpxc validate --file` reports VALID; only --conformance errors.
- **Spec:** psm.md:40 and feature.yaml:20-22: peptidoform is required ProForma. peptidoform.md:11-30 sets the canonical form.
- **Reader impact:** Bogus peptidoform 'None' rows enter identity hashes, joins and peptide counts, and basic validation passes them. Per GitHub code (SWAPS, kroganlab), timsTOF TIMS-MULTI-MATCH rows leave Modified sequence empty.
- **Suggested:** Treat NULL as missing: skip the row, or derive the peptidoform from Sequence plus fixed mods. Basic validation rejects non-ProForma peptidoforms.
- **Evidence:** `qpx/converters/maxquant/psm_adapter.py:259,265-270`; `qpx/converters/maxquant/feature_adapter.py:441`; `audit-run:verify-maxquant/psmnull/out/maxquant.psm.parquet`; `audit-run:verify-maxquant/types/out/maxquant.feature.parquet`; `gh search/code 'TIMS-MULTI-MATCH' (wilhelm-lab/SWAPS, kroganlab/MartinGordon_projects)`

### MAXQUANT-M3: Sidecar files silently change peptidoforms and labels but are not in provenance {#maxquant-m3}

*cross* · diverges from spec · found by the verifier · checklist D5, D8, P9

- **Now:** The converter reads parameters.txt, then summary.txt, from beside the input for fixed mods (base_adapter.py:192-196). SILAC also reads summary.txt (psm_adapter.py:140) and msms.txt beside evidence.txt, even without --msms-file (silac.py:100). Provenance records only CLI options with SHA-256 (cli/convert/postprocess.py:174-200). The simple run lists only msms, evidence and sdrf.
- **Spec:** provenance.md:65: an input path is recorded by file name and an input file by its SHA-256.
- **Reader impact:** One evidence.txt converted in two folders gives different peptidoforms, calculated_mz and feature_ids, and nothing in the output explains why.
- **Suggested:** Record every file read (sidecars included) with name and SHA-256, and say which source gave the fixed mods.
- **Evidence:** `qpx/converters/maxquant/base_adapter.py:192-196`; `qpx/converters/maxquant/psm_adapter.py:140`; `qpx/converters/maxquant/silac.py:100`; `qpx/cli/convert/postprocess.py:174-200`; `docs/spec/provenance.md:65`; `audit-run:maxquant/simple/maxquant.provenance.parquet`

## FragPipe converter

### FRAGPIPE-07: contam_ marker stripped from accessions in psm, feature, pepmap and pg {#fragpipe-07}

*cross* · diverges from spec · confirmed · checklist P8, P4, G3

- **Now:** protein_accessions() keeps only the decoy prefixes (rev_, REV_, DECOY_, decoy_) and strips 'contam_' through parse_uniprot_id (constants.py:28,42-58, utils.py:30-42). Verified: 'contam_sp|P02768|ALBU_HUMAN' becomes 'P02768' in psm.protein_accessions, feature anchor/pg_accessions, pepmap and pg. Only pg sets contaminant=true (pg_adapter.py:469); psm and feature have no contaminant signal.
- **Spec:** pg.md:46-48: a decoy or contaminant marker stays on the accession (CONTAM_P02768). producer-guide.md:34: contaminant -> keep the marker. pg.md:365: FragPipe contaminants (contam_) are kept.
- **Reader impact:** PSM and feature readers cannot exclude contaminants. A contaminant copy of a target protein would share its accession, which risks merged groups or pg_id collisions (not verified on real data).
- **Suggested:** Keep 'contam_' (normalised to CONTAM_) on accessions in all views. The contract should say how psm and feature readers identify contaminant evidence.
- **Evidence:** `qpx/converters/fragpipe/constants.py:28`; `qpx/converters/fragpipe/constants.py:42-58`; `qpx/converters/utils.py:30-42`; `qpx/converters/fragpipe/pg_adapter.py:469`; `docs/spec/pg.md:46-48`; `docs/spec/producer-guide.md:34`; `run: out/fp.psm.parquet CNTAMK -> [P02768]; out/fp.pg.parquet P02768 contaminant=True`

### FRAGPIPE-24: FragPipe tables in tool-mappings and converter-coverage contradict the code {#fragpipe-24}

*other* · diverges from spec · confirmed · checklist X1

- **Now:** The code builds the peptidoform from Peptide + Assigned Modifications (psm_adapter.py:203-205). run_file_name and scan come from Spectrum tokens (psm_adapter.py:229-236), observed_mz is calibrated (:239), and is_decoy is derived from the prefix (:160). On features, sequence is 'Peptide Sequence' and charge is 'Charge' (column_mappings.yaml:11,16). On pg, is_decoy and contaminant are filled (pg_adapter.py:442,469).
- **Spec:** tool-mappings.md:20 says peptidoform = 'Modified Peptide'; :28 run_file_name = 'Spectrum File' (a pepXML path); :29 scan = 'Scan' (no such column); :24,44,111-112 is_decoy and contaminant '---'; :43 feature charge '---'. converter-coverage.md:303: contaminant 'Not implemented'.
- **Reader impact:** Readers and implementers of other producers copy wrong source columns, for example taking run names from 'Spectrum File' or believing that FragPipe pg has no contaminant flag.
- **Suggested:** Regenerate the FragPipe columns of tool-mappings and converter-coverage from column_mappings.yaml and the adapters.
- **Evidence:** `docs/spec/tool-mappings.md:20`; `docs/spec/tool-mappings.md:24`; `docs/spec/tool-mappings.md:28-29`; `docs/spec/tool-mappings.md:43`; `docs/spec/tool-mappings.md:111-112`; `docs/spec/converter-coverage.md:303`; `qpx/converters/fragpipe/psm_adapter.py:203-239`; `qpx/converters/fragpipe/pg_adapter.py:442`; `qpx/converters/fragpipe/pg_adapter.py:469`

### FRAGPIPE-25: Only one psm.tsv accepted; multi-experiment studies get partial PSM coverage {#fragpipe-25}

*cross* · unspecified · confirmed · checklist P3, X1

- **Now:** The CLI takes a single --psm-file (cli/convert/fragpipe.py:25-29), which feeds both the psm view and the feature enrichment lookup (converter.py:64-76). FragPipe writes one psm.tsv per experiment folder; this layout is from tool knowledge and was not verified locally. Unless the user concatenates the files, the psm view covers one experiment, and features of other experiments get 'direct' with id_run_file_name set but no scan or PEP.
- **Spec:** Nothing. converter-coverage.md:168 measured 'three one-run experiments' without saying how the psm.tsv files were combined.
- **Reader impact:** The psm view can silently cover only part of the study. Counts and links between psm and feature then depend on the inputs the user happened to pass.
- **Suggested:** Accept several psm.tsv files (or a FragPipe results folder), and record the inputs that were covered. The contract should require the psm view's run coverage to be declared, or to equal run.parquet.
- **Evidence:** `qpx/cli/convert/fragpipe.py:25-29`; `qpx/converters/fragpipe/converter.py:64-76`; `docs/spec/converter-coverage.md:168`

### FRAGPIPE-27: compensation_voltage 0 vs null; PSM CV is a string cv_param {#fragpipe-27}

*cross* · inconsistent across producers · confirmed · checklist X1, F1

- **Now:** FragPipe writes compensation_voltage=0.0 when the column is absent or empty (feature_adapter.py:519), because CV is part of its feature identity. Other converters leave it null. On PSMs, the CV is a cv_param {'compensation_voltage': '-45.0'} string (psm_adapter.py:317-319), not a column.
- **Spec:** feature.md:76 and feature.yaml:128-131: 0 for non-FAIMS FragPipe units. There is no PSM CV field, and scores.md:126-127 reserves cv_params for non-numeric data.
- **Reader impact:** 'CV = 0' can mean non-FAIMS (FragPipe) or a real 0 V setting, while null means unknown elsewhere. A PSM CV must be parsed from strings.
- **Suggested:** Use null for non-FAIMS everywhere, and handle nulls in the identity hash. Add a typed psm.compensation_voltage.
- **Evidence:** `qpx/converters/fragpipe/feature_adapter.py:519`; `qpx/converters/fragpipe/psm_adapter.py:317-319`; `docs/spec/feature.md:76`; `qpx/core/data/schemas/feature.yaml:128-131`; `docs/spec/scores.md:126-127`

### FRAGPIPE-29: Footer and provenance omit the FragPipe pipeline: version, validator, IonQuant, MBR, protein FDR {#fragpipe-29}

*footer* · unspecified · partly confirmed (see correction) · checklist D8

- **Now:** The footer has qpx_version, creator=fragpipe, identity_composite, identity_source=derived, primary_key and hash version, but no scan_format and no filter, rank or decoy policy (verified). Provenance step 1 is 'FragPipe/MSFragger' with tool_version and parameters null. There are no steps for PeptideProphet or Percolator, Philosopher filter, ProteinProphet or IonQuant (converter.py:95-120). dataset.software_version is null (converter.py:205-211).
- **Spec:** provenance.md:3: answer 'what FDR was used?'. provenance.md:76-85: rescoring, fdr_filtering, quantification and protein_inference categories, with IonQuant and ProteinProphet named. identity.md:20-23: identity footer keys. No footer key for filtering, rank or MBR policy exists.
- **Correction:** The gist holds. However, --software-name, --software-version and --software-config let the user add a workflow provenance step, a version and a JSON config, and the format_conversion step records the inputs and their SHA-256. Nothing parses fragpipe.workflow, and no policy keys exist.
- **Reader impact:** Nothing in a FragPipe QPX dataset states the FDR levels, MBR on or off and its FDR, the validator, or the tool versions, so cross-producer comparison relies on outside knowledge.
- **Suggested:** The contract should define mandatory DDA policy keys (psm_filter_level and threshold, rank_policy, decoys_included, mbr, intensity_type). FragPipe fills them from fragpipe.workflow when given.
- **Evidence:** `qpx/converters/fragpipe/converter.py:95-120`; `qpx/converters/fragpipe/converter.py:205-211`; `docs/spec/provenance.md:3`; `docs/spec/provenance.md:76-85`; `docs/spec/identity.md:20-23`; `run: out/fp.provenance.parquet, out/fp.dataset.parquet`

### FRAGPIPE-M3: ontology.parquet source_column_name names columns the code does not read {#fragpipe-m3}

*other* · diverges from spec · found by the verifier · checklist X1, P7

- **Now:** Field provenance rows come from the YAML first-match (converter.py:194-201, core/scores.py:303-320), not from the code path. Verified rows: psm observed_mz <- 'Observed M/Z' (the code reads Calibrated, psm_adapter.py:239); psm calculated_mz <- 'Calculated M/Z' (recomputed, :165); feature calculated_mz <- 'M/Z'; pg spectral_count <- 'Combined Spectral Count' (never read; per-experiment counts are written, pg_adapter.py:39-43).
- **Spec:** ontology.md:54: source_column_name is the 'original column name in tool output'.
- **Reader impact:** Readers who audit the meaning of a field through ontology.parquet are told the uncalibrated m/z, or a combined count, was used. This hides the per-producer semantic differences the contract wants made explicit.
- **Suggested:** Emit source_column_name from the column actually used, or mark it 'derived'. The contract should require field provenance to reflect any transformation.
- **Evidence:** `qpx/converters/fragpipe/converter.py:194-201`; `qpx/core/scores.py:303-320`; `qpx/converters/column_mappings.yaml:19-39`; `docs/spec/ontology.md:54`; `run: work/verify-fragpipe/out/fp.ontology.parquet rows 3,5,8,20`

## Proteome Discoverer converter

### PD-02: Missing or unparseable 'Rank' column makes every PSM rank 1 {#pd-02}

*cross* · diverges from spec · confirmed · checklist P1, F2

- **Now:** reports.py:91 computes COALESCE(TRY_CAST(Rank),1). An export without a Rank column, or with non-numeric values, makes every hit rank 1, and no 'rank' score is written (reports.py:138-139 reads the literal 'Rank' column). feature_adapter.py:18 then aggregates rank-n hits into features, and --rank1-only filters nothing. No warning is logged.
- **Spec:** converter-coverage.md:324: features use 'The rank 1 PSMs of a precursor'. The spec has no rule on what to do when rank is unknown.
- **Reader impact:** Features and reporter sums can include alternative, wrong peptides of the same spectrum, with no trace in the output.
- **Suggested:** When Rank is absent and several hits share a scan, refuse or warn. Record in provenance whether rank was available.
- **Evidence:** `qpx/converters/proteome_discoverer/reports.py:91`; `qpx/converters/proteome_discoverer/reports.py:138-139`; `qpx/converters/proteome_discoverer/feature_adapter.py:18`; `docs/spec/converter-coverage.md:324`

### PD-03: Several rank-1 peptides on one spectrum each become 'direct' features sharing its reporters {#pd-03}

*cross* · ambiguous · confirmed · checklist P2, F2, D1

- **Now:** There is no chimeric handling. Features group by peptide and charge (feature_adapter.py:16-20), so each rank-1 PSM of one spectrum becomes its own feature with the same scan and the spectrum's full reporter ions. This covers 'Ambiguous' hits, CHIMERYS or second-search co-identifications, and other charges. Only 'Rejected' is excluded (feature_adapter.py:72), and identification_type is always 'direct' (:139). Run: two Ambiguous peptides on scan 100 gave two features, each TMT126=10.
- **Spec:** feature.md:145: 'secondary' means a co-fragmenting peptide of a chimeric spectrum. The spec has no rule separating ambiguous candidates from true chimeras, and none on whether one spectrum's quantity may feed several features.
- **Reader impact:** Protein rollups double-count one spectrum's reporter signal across peptides and proteins. Chimeric co-IDs cannot be told from ambiguity (#378: 758 of 15,582 PD LFQ spectra have several rank-1 peptides).
- **Suggested:** Contract: one spectrum's quantity feeds at most one feature unless that feature is flagged 'secondary'. Map PD Ambiguous/CHIMERYS co-IDs to secondary, or exclude them. Make psm_ambiguity a typed field.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:16-20`; `qpx/converters/proteome_discoverer/feature_adapter.py:72`; `qpx/converters/proteome_discoverer/feature_adapter.py:139`; `docs/spec/feature.md:145`; `run: audit-run:pd/cases/ambiguous_same_scan/out/pd.feature.parquet`

### PD-05: Identifying PSM picked by lowest q-value; psm and feature tie-breaks differ and are nondeterministic {#pd-05}

*cross* · ambiguous · partly confirmed (see correction) · checklist F2, F5, P6

- **Now:** Feature scan, PEP, scores, rt and precursor_qvalue come from row_number() ORDER BY _qvalue NULLS LAST, _scan (feature_adapter.py:17). The psm dedup orders by _rank, _qvalue NULLS LAST with no final key (psm_adapter.py:40). Per #378, PD's own peptide groups use the lowest-PEP PSM, and qpx's pick differs in 660 of 3,026 multi-spectrum LFQ precursors. In ~366 of 43,121 TMT features the joined psm row has a different PEP and scores, varying with thread order.
- **Spec:** feature.md:138-139 and :301: scan and id_run_file_name come from 'the best PSM', but 'best' is undefined. #378 (open) proposes the tool's own pick, else lowest PEP, ties broken by lowest scan, with the same order in both views.
- **Correction:** 'Best' is specified for PD: converter-coverage.md:324 says the identification fields come from 'the lowest q-value PSM'. Left unspecified are the tie order and the mismatch with PD's own lowest-PEP pick (#378). Ties across scans are broken by _scan. Within one scan (two engines) neither view has a final key. In case 'tie' both views kept the PEP 0.05 row although a PEP 0.001 hit existed.
- **Reader impact:** Feature PEP, q-value and scores are not PD's values and may not match the joined psm row. Re-running the conversion can change them.
- **Suggested:** Adopt #378 rule 4 for PD: lowest PEP, then q-value, then scan, as one total order used by both views.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:17`; `qpx/converters/proteome_discoverer/psm_adapter.py:40`; `docs/spec/feature.md:138-139`; `docs/spec/feature.md:301`; `bigbio/qpx#378 (Proteome Discoverer row)`

### PD-06: psm.feature_id never filled, although features are aggregated from those PSMs {#pd-06}

*cross* · unspecified · confirmed · checklist P10, P6

- **Now:** No adapter writes psm.feature_id (psm_adapter.py:48-74). The CLI summary on a converted set reports '0 features with >=1 psm'. The natural key (id_run_file_name, scan) is not unique: #378 found several PSMs per key in 440 of 5,843 LFQ and 10,661 of 43,121 TMT PD features, because of rank-n and chimeric hits.
- **Spec:** psm.md:75: feature_id is optional and authoritative when present. data-model.md:162: shared fields 'are not a guaranteed row-level link'. The spec does not say which PSMs a feature aggregates.
- **Reader impact:** Readers cannot tell which PSMs fed a feature's summed reporters. Recovering the identifying PSM needs peptidoform and charge in the key.
- **Suggested:** Contract: producers that aggregate PSMs into features fill psm.feature_id for every PSM used. The natural key includes peptidoform and charge.
- **Evidence:** `qpx/converters/proteome_discoverer/psm_adapter.py:48-74`; `docs/spec/psm.md:75`; `docs/spec/data-model.md:162`; `qpxc convert proteome-discoverer summary in audit-run:pd/cases/rank2/cli_out`; `bigbio/qpx#378`

### PD-09: No FDR filtering anywhere; the PD export's filtering state is unknown and unrecorded {#pd-09}

*cross* · unspecified · confirmed · checklist P3, D8, G1

- **Now:** The psm view has no q-value or Confidence filter (psm_adapter.py:37-40). The feature view filters only on rank, 'Rejected' and 'Excluded by…' (feature_adapter.py:18,72-76), so Low-confidence rank-1 PSMs still create 'direct' features. Filters applied in PD before export cannot be detected. The database_search provenance step has parameters=None (converter.py:91-112).
- **Spec:** #241: no filtering by default. pg.md:126-128 defines identification evidence as a PSM 'that passed the FDR threshold'. provenance.md:3 says provenance should answer 'what FDR was used?'. No footer key states the filter level or threshold.
- **Reader impact:** The feature and psm views mix confident and non-confident identifications. pg's evidence rule cannot be applied, and the same view means different things across PD exports.
- **Suggested:** Contract footer or provenance keys per view: fdr_filter_level (none|psm|peptide|protein|unknown) and fdr_threshold.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:18`; `qpx/converters/proteome_discoverer/feature_adapter.py:72-76`; `qpx/converters/proteome_discoverer/converter.py:91-112`; `docs/spec/pg.md:126-128`; `docs/spec/provenance.md:3`

### PD-11: No MBR representation: features always 'direct', while pg may include Feature Mapper quantities {#pd-11}

*cross* · unspecified · confirmed · checklist F4, D1, D2, G1

- **Now:** Features come only from PSMs. identification_type is hard-coded 'direct' (feature_adapter.py:139) and match_score/match_qvalue/match_rt_delta/match_mz_delta are never set. PeptideGroups.txt and 'Found in Sample' are not read. pg intensity is PD's protein Abundance (pg_adapter.py:111-123). In Minora with Feature Mapper workflows that abundance can include quantities mapped across runs; this PD behaviour was not verified here.
- **Spec:** feature.md:145-147 define identification_type and match_*. pg.md:131-134 allows quantity rows in units without identification (MBR). Nothing flags that a pg quantity includes transferred evidence.
- **Reader impact:** pg cannot be reproduced from features. At protein level, transferred and identified quantities are indistinguishable (mokume#117).
- **Suggested:** Contract: pg rows state their transferred share (flag or count). The PD converter could read PeptideGroups 'Found in Sample' to emit transferred features.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:139`; `qpx/converters/proteome_discoverer/pg_adapter.py:111-123`; `docs/spec/feature.md:145-147`; `docs/spec/pg.md:131-134`

### PD-23: is_decoy always false; contaminant false when the column is absent {#pd-23}

*cross* · diverges from spec · partly confirmed (see correction) · checklist P4, F8, G3

- **Now:** is_decoy is hard-coded False (psm_adapter.py:60, feature_adapter.py:123, pg_adapter.py:149), with no check for decoy rows. Contaminant is read only on pg as Contaminant=='True' (pg_adapter.py:150), which gives False when the column is missing. psm/feature contaminant status is not represented, and PD accessions carry no CONTAM_ marker.
- **Spec:** producer-guide.md:28-35: decoy and contaminant markers are kept. pg.md:95: contaminant is nullable. pg.md:99-103: decoy groups are kept when reported (#342).
- **Correction:** is_decoy is always false, and a missing Contaminant column gives false: both confirmed. But producer-guide.md:34 treats contaminants as targets (is_decoy false) and carries contaminant status only on pg. So psm/feature rows 'looking like targets' follows the spec. The real divergence is false instead of null (pg.md:95).
- **Reader impact:** Unknown contaminant status reads as 'not a contaminant', and psm/feature contaminants look like targets.
- **Suggested:** Set contaminant to null when the column is absent. Carry the PSM 'Contaminant' flag as a typed field or cv_param.
- **Evidence:** `qpx/converters/proteome_discoverer/psm_adapter.py:60`; `qpx/converters/proteome_discoverer/feature_adapter.py:123`; `qpx/converters/proteome_discoverer/pg_adapter.py:149-150`; `docs/spec/producer-guide.md:28-35`; `docs/spec/pg.md:95`

### PD-24: Two run-name normalizations; feature and pg runs disagree when Spectrum File ends in '.RAW' {#pd-24}

*cross* · diverges from spec · confirmed · checklist P7, F10

- **Now:** With 'Spectrum File', the run is a regex strip of path and a case-sensitive raw|mzML|mzml|d suffix (reports.py:99). Through InputFiles it is run_file_stem (reports.py:49; qpx/core/files.py:16-19). pg grouped_runs uses InputFiles when given (converter.py:48). Run: Spectrum File 'plex1_F1.RAW' gave psm/feature run 'plex1_F1.RAW' but pg grouped_runs ['plex1_F1'].
- **Spec:** data-model.md:22: run_file_name = raw file name without extension. data-model.md:159: feature→pg requires run_file_name ∈ grouped_runs.
- **Reader impact:** The feature→pg softlink and the run/sample joins fail silently, and the extension leaks into ids.
- **Suggested:** Use run_file_stem on both routes.
- **Evidence:** `qpx/converters/proteome_discoverer/reports.py:49`; `qpx/converters/proteome_discoverer/reports.py:99`; `qpx/converters/proteome_discoverer/converter.py:48`; `docs/spec/data-model.md:22`; `run: audit-run:pd/cases/spectrum_file_upper`

### PD-30: pepmap covers only quantified peptidoforms {#pd-30}

*pepmap* · unspecified · confirmed · checklist P8, X1

- **Now:** pepmap entries are added in the feature record after the quantity check (feature_adapter.py:99-107) and written only when a feature view exists (converter.py:53-57). Peptidoforms seen only in the psm view, and identification-only exports, get no pepmap row, although their PSMs carry Protein Accessions.
- **Spec:** producer-guide.md:16: pepmap lists every protein a peptide maps to. data-model.md:163: psm and feature join pepmap on peptidoform.
- **Reader impact:** The psm→pepmap join is partial for PD.
- **Suggested:** Build pepmap from all PSMs.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:99-107`; `qpx/converters/proteome_discoverer/converter.py:53-57`; `docs/spec/producer-guide.md:16`; `docs/spec/data-model.md:163`

### PD-31: Provenance lacks PD version and settings; nothing declares filtering, ranks or quan basis {#pd-31}

*footer* · unspecified · confirmed · checklist D8, D7

- **Now:** converter.py:91-112: the database_search step has tool_version and parameters None, and dataset software_version is None (converter.py:88). Footers carry qpx_version, the schema-default identity_composite and identity_source=derived; no scan_format (verified). Only the CLI records rank1_only (postprocess.py:204-229). No imputation, rt_aligned, identified_labels or channel_details values are written.
- **Spec:** provenance.md:3: provenance should answer 'what FDR was used?'. intensities.md:213: processing keys. The spec defines no footer key for FDR level, rank policy or intensity kind.
- **Reader impact:** Readers cannot learn the PD version (header semantics vary across versions), FDR settings, S/N vs intensity, normalization or imputation.
- **Suggested:** Contract footer/provenance keys: psm_ranks, fdr_level, fdr_threshold, intensity_kind, normalization, imputation. Converters write 'unknown' explicitly.
- **Evidence:** `qpx/converters/proteome_discoverer/converter.py:88`; `qpx/converters/proteome_discoverer/converter.py:91-112`; `qpx/cli/convert/postprocess.py:204-229`; `docs/spec/provenance.md:3`; `docs/spec/intensities.md:213`

### PD-33: Tests pin spec-divergent behaviour; no PD fixture covers chimeras, multi-engine quantification or MBR {#pd-33}

*other* · unspecified · confirmed · checklist X1

- **Now:** tests/converters/test_proteome_discoverer_converter.py pins: High→pass_threshold (115), engine dedup and rank1_only (136-148), feature exclusions and unique=True with a 2-member group (151-169), pg drop of an empty channel (197-210), and the identification-only 'LFQ' row (213-227). tests/examples has no PD fixture. Nothing tests ambiguous/chimeric spectra, two-engine feature sums, a missing Rank column, run-name normalization or SILAC pg columns.
- **Spec:** producer-guide.md:159-161: the fixtures in tests/examples are references for how a tool maps to QPX.
- **Reader impact:** Tests lock in the divergences, and there is no golden PD dataset for the contract.
- **Suggested:** Add a small real PD export fixture (multi-engine TMT plus Minora LFQ), with one test per contract rule.
- **Evidence:** `tests/converters/test_proteome_discoverer_converter.py:115`; `tests/converters/test_proteome_discoverer_converter.py:136-169`; `tests/converters/test_proteome_discoverer_converter.py:197-227`; `docs/spec/producer-guide.md:159-161`

## mzIdentML and CDAP converters

### MZID_CDAP-05: CDAP strips decoy prefix XXX_; mzIdentML keeps DECOY_ marker {#mzid_cdap-05}

*cross* · diverges from spec · confirmed · checklist P4, P8, F8, G3, X1

- **Now:** CDAP removes XXX_ from every accession (cdap/base_adapter.py:140-142) in psm.protein_accessions, feature.pg_accessions/anchor_protein and pg.pg_accessions. Fixture: 10/100 decoy PSMs (e.g. 'NP_005751.2'), 8/94 decoy features, 6/81 decoy groups all carry target-looking accessions. mzid keeps 'DECOY_P00001' (synthetic) after canonical_accession.
- **Spec:** pg.md:46-49 and producer-guide.md:33: keep the decoy marker so a decoy never shares an accession with its target; is_decoy states status.
- **Reader impact:** Decoy evidence joins target protein records (pepmap, pg, FASTA, gene maps); anyone ignoring is_decoy counts decoys as real proteins.
- **Suggested:** Keep the source decoy marker verbatim in all views; is_decoy carries the status.
- **Evidence:** `qpx/converters/cdap/base_adapter.py:140-147`; `audit-run:mzid_cdap/cdap/*.parquet decoy rows`; `docs/spec/pg.md:46-49`; `docs/spec/producer-guide.md:33`

### MZID_CDAP-09: mzIdentML footer declares a composite that does not reproduce its psm_ids {#mzid_cdap-09}

*footer* · contradictory · confirmed · checklist D8, P6

- **Now:** psm_id = derive_id([run_file, [SII ids]]) (psm_records.py:221-235); footer still says identity_composite=run_file_name,scan,peptidoform,charge with identity_source=provided (observed PXD054720 footer).
- **Spec:** identity.md:18 composite = columns ids are hashed from; identity.md:40 documents this mzid exception.
- **Reader impact:** Ids cannot be re-derived or checked; the same identification gets different psm_id from mzid and from derived producers.
- **Suggested:** Footer must state the real recipe (e.g. composite including an SII-id column) or omit identity_composite for provided ids.
- **Evidence:** `qpx/converters/mzidentml/psm_records.py:221-235`; `audit-run:mzid_cdap/pxd/PXD054720.psm.parquet footer`; `docs/spec/identity.md:18,40`

### MZID_CDAP-20: Mixed index/scan references leave scan_format undeclared, read as scan numbers {#mzid_cdap-20}

*footer* · contradictory · confirmed · checklist D8, P7

- **Now:** When records mix 'scan' and 'index' schemes (multi-run mzid with one MGF converts only the matching run, mzidentml/converter.py:163-169; or unresolved MGF lookups), declare_scan_format only warns and leaves scan_format unset (psm_records.py:38-43). From code; not exercised by fixtures.
- **Spec:** scan.md:80-81: a file without scan_format is read as holding scan numbers.
- **Reader impact:** 0-based indices of other runs are silently treated as scan numbers; spectra get attached/joined wrongly.
- **Suggested:** Per-run/per-row scan format (#259) or refuse mixed output; never leave undeclared when index values exist.
- **Evidence:** `qpx/converters/mzidentml/psm_records.py:38-43`; `qpx/converters/mzidentml/converter.py:157-195`; `docs/spec/scan.md:74-81`

### MZID_CDAP-21: No footer or provenance declares filtering state, FDR level or rank policy {#mzid_cdap-21}

*footer* · unspecified · partly confirmed (see correction) · checklist D8, P3

- **Now:** mzid psm footer: scan_format, identity keys, spectrum_content only with --include-spectra (mzidentml/converter.py:131-134). CDAP psm/feature footers have no scan_format (observed). No key for filter level/threshold, ranks kept or decoys kept. Provenance: mzid records rank1_only=False (observed) but no source thresholds; CDAP records psm_dir only (no structures) and MS-GF+ parameters=None.
- **Spec:** producer-guide.md:23 set scan_format for DDA PSM/feature; provenance.md:65 record options incl. defaults (e.g. structures); no mandated filter/rank keys.
- **Correction:** provenance.md:65 already requires recording the converter's own filter options with their defaults (qvalue_threshold 'none', rank1_only); mzid records rank1_only=False. CDAP leaves out unset defaults such as structures (cli/convert/postprocess.py _effective_parameters skips None). What is truly unspecified is the upstream/source filtering state, level and threshold, and whether decoys or ranks are present.
- **Reader impact:** From a file alone a reader cannot know if it is filtered, at which level, whether ranks>1 or decoys are present.
- **Suggested:** Mandatory footer keys, e.g. psm_filter_level, psm_filter_threshold, psm_ranks, decoys_included (values incl. 'unknown').
- **Evidence:** `qpx/converters/mzidentml/converter.py:117-136,290-316`; `qpx/converters/cdap/converter.py:184-217`; `footers of audit-run:mzid_cdap/{pxd,cdap}/*.parquet`; `docs/spec/producer-guide.md:23`; `docs/spec/provenance.md:65`

### MZID_CDAP-28: pg peptides.peptide_count means different things in mzIdentML and CDAP {#mzid_cdap-28}

*cross* · inconsistent across producers · partly confirmed (see correction) · checklist X1, G2

- **Now:** mzid: per member protein, count of PeptideHypothesis (PSM evidences) (mzidentml/pg_adapter.py:187-190). CDAP: one entry for the anchor only, count of distinct sequences in the mapping-set group (cdap/pg_adapter.py:144-146).
- **Spec:** pg.md:159-166 number of peptides per individual protein in the group.
- **Correction:** The mzid meaning also varies within mzIdentML: it is a PSM count for Mascot-style files with one PeptideEvidence per query, and a peptide-evidence count otherwise. CDAP gives distinct sequences for the anchor only.
- **Reader impact:** Min-peptide filters (e.g. >=2 peptides) keep almost every mzid group and behave differently for CDAP.
- **Suggested:** Define peptide_count = distinct peptidoforms (or sequences) passing threshold, per member protein, for every producer.
- **Evidence:** `qpx/converters/mzidentml/pg_adapter.py:187-190`; `qpx/converters/cdap/pg_adapter.py:144-146`; `docs/spec/pg.md:159-166`

### MZID_CDAP-30: Rank and PSM-count representation differ: mzid rank score, CDAP none {#mzid_cdap-30}

*cross* · inconsistent across producers · partly confirmed (see correction) · checklist P1, X1, F2

- **Now:** mzid writes a 'rank' score on every row (psm_records.py:140). CDAP writes no rank (its input holds top hits and ties). CDAP feature PSM count is a cv_param string 'psm_count' (cdap/feature_adapter.py:385-388), whereas MaxQuant/PD put spectral_count in additional_intensities (converter-coverage.md:194,324).
- **Spec:** rank only 'recommended' (psm.md:291); spectral_count is an additional_intensities name (pg.md:233).
- **Correction:** The spec does name a place for a precursor's PSM count: intensities.md:187 gives spectral_count in additional_intensities ('the precursor's PSMs', PD). CDAP's cv_param string 'psm_count' therefore diverges from the spec rather than filling a gap. The MaxQuant/PD refs are converter-coverage.md:195,324. The rank part holds.
- **Reader impact:** A rank=1 filter drops all CDAP rows (rank null); spectral counts must be read from a different place per producer.
- **Suggested:** Typed rank (null = not reported; CDAP = 1) and a single location/name for PSM or spectral counts.
- **Evidence:** `qpx/converters/mzidentml/psm_records.py:140`; `qpx/converters/cdap/feature_adapter.py:385-388`; `docs/spec/converter-coverage.md:194,324`

### MZID_CDAP-M2: mzIdentML pepmap built from every SII, incl. rank>1 and failing hits {#mzid_cdap-m2}

*pepmap* · unspecified · found by the verifier · checklist P1, P3, P8

- **Now:** _collect_pepforms adds the peptidoform and proteins of every PSM record, with no rank or pass_threshold check (mzidentml/converter.py:129,213-223). PXD054720: 307 of 2,463 pepmap peptidoforms (12.5%) have no PSM with pass_threshold=true.
- **Spec:** pepmap.md:3,9-10 describe a peptide-protein map for uniqueness and coverage of 'identified peptides'. Which identifications populate it is not stated.
- **Reader impact:** Peptide uniqueness (is_unique), protein coverage and gene maps are computed partly from rejected or alternative hits.
- **Suggested:** pepmap lists peptidoforms of PSMs passing the producer's threshold (or rank 1), or declares its population rule in the footer.
- **Evidence:** `qpx/converters/mzidentml/converter.py:129,213-223`; `verify-mzid_cdap/pxd: pepmap 2463 vs 2156 peptidoforms with a passing PSM`; `docs/spec/pepmap.md:3-10`

### MZID_CDAP-M3: CDAP contaminant entries lose their identity; pg.contaminant hard-coded false {#mzid_cdap-m3}

*cross* · diverges from spec · found by the verifier · checklist P4, P8, G3

- **Now:** CDAP cells mix RefSeq and contaminant entries, e.g. 'sp|Q29443|TRFE_BOVIN' (fixture scans 5646/5652). The accession normalizer rewrites it to bare 'Q29443' in psm.protein_accessions (log: '2 accession changes'), with no marker. pg.contaminant=False always (cdap/pg_adapter.py:167), and no contaminant parsing exists (base_adapter.py:130-147).
- **Spec:** producer-guide.md:34: keep the contaminant marker and set contaminant=true on pg. pg.md:46-49 and 95 (contaminant is bool, null). converter-coverage.md:303 documents CDAP 'Always false'.
- **Reader impact:** Bovine/keratin/trypsin groups count as sample proteins; contaminant filtering is impossible for CDAP data.
- **Suggested:** Detect contaminants from the source (non-RefSeq sp| entries or a contaminant list); keep the marker; contaminant=true, or null when unknown, never a blanket false.
- **Evidence:** `tests/examples/cdap/sample_run*.psm Protein cell at scans 5646/5652`; `qpx/converters/cdap/pg_adapter.py:167`; `verify-mzid_cdap/cdap/cdap.psm.parquet contains 'Q29443'`; `docs/spec/producer-guide.md:34`; `docs/spec/converter-coverage.md:303`

## Shared writer and validation code

### INFRA-07: q-value score names change level by view; ontology entries contradict {#infra-07}

*cross* · contradictory · confirmed · checklist D4, P5, X1

- **Now:** 'qvalue' fills psm.qvalue (PSM) and feature.precursor_qvalue (precursor) (derived.py:29,32). Its builtin entry has accession MS:1002354 'PSM-level q-value' but the description 'Run-level precursor q-value (DIA-NN Q.Value)' (builtin_scores.py:185-191). 'global_qvalue' maps to MS:1002350 'PSM-level global FDR' with the description 'DIA-NN Global.Q.Value' (:241-247). It fills feature.global_qvalue but is excluded from psm.qvalue: a PSM row with only global_qvalue gets qvalue null (writer test).
- **Spec:** psm.md:291 recommends global_qvalue as the 'experiment-level PSM q-value'. scores.md:45 defines global_qvalue as 'Precursor, per experiment'. scores.md:56 calls it 'global q-value at the experiment level' for both PSM and Feature.
- **Reader impact:** One score name means a PSM-level, precursor-level or experiment-level FDR depending on the view and producer. ontology.parquet labels DIA precursor q-values as PSM-level terms.
- **Suggested:** Use level-qualified reserved names (psm_qvalue, psm_global_qvalue, precursor_qvalue, precursor_global_qvalue, peptide_qvalue) and give each a correct PSI-MS accession. Add psm.global_qvalue, or drop the psm.md:291 recommendation.
- **Evidence:** `qpx/core/builtin_scores.py:185-191`; `qpx/core/builtin_scores.py:241-247`; `qpx/writers/derived.py:16,29-34`; `docs/spec/psm.md:291`; `docs/spec/scores.md:45,56`; `PSI-MS lookup: MS:1002354='PSM-level q-value', MS:1002350='PSM-level global FDR'`

### INFRA-08: Typed q-value columns and additional_scores can disagree; no check {#infra-08}

*cross* · ambiguous · confirmed · checklist D4

- **Now:** fill_score_columns only coalesces nulls (derived.py:58-71): a typed value set explicitly is kept even when additional_scores holds a different value of the same name and level. Writer test: precursor_qvalue=0.5 kept while additional_scores qvalue=0.014; psm qvalue=0.001 kept while 'q-value'=0.9. The score also stays in additional_scores. Even with no conflict, the value is stored twice.
- **Spec:** scores.md:38 says writers copy q-values, and scores.md:47 says the score stays in additional_scores. Nothing says what happens when they disagree. #377 (duplication) is open.
- **Reader impact:** SQL on additional_scores and SQL on the typed column can give different confidence for the same row.
- **Suggested:** The typed column is authoritative. The writer refuses a row whose same-name score differs, or moves the value out of additional_scores once copied (resolving #377).
- **Evidence:** `qpx/writers/derived.py:58-71`; `docs/spec/scores.md:36-47`; `writer test audit-run:infra/scripts/writer_exp.py rows scan 4,5 (feature) and 6 (psm)`

### INFRA-12: Write-time validation depends on the write path; one warning masks another {#infra-12}

*cross* · inconsistent across producers · confirmed · checklist X1, D6

- **Now:** Record path (write_batch): only required-non-null, PK uniqueness, pg referential and psm scan_format checks (base.py:462-479,599-617); no checks of ranges, vocabulary, score names or NaN. Table path: validate_full(strict=False) gives warnings only (base.py:414-419). _log_validation_warnings logs each check name once per file (base.py:571-576): an out-of-range peptide_qvalue warning hid the NaN/negative intensity warning in the same file. Unknown keys raise on records but are dropped with a warning on tables (base.py:298-306,385-396).
- **Spec:** producer-guide.md:60 says each writer validates on close(); serialization.md says `qpxc validate` is the strict gate. Nothing says that a written file can be invalid.
- **Reader impact:** Whether a converter's output meets the contract depends on which writer method it calls. Files with NaN intensities or out-of-range q-values are published.
- **Suggested:** Run the same value checks on every path and raise for contract violations (ranges, vocabularies, NaN). Report each issue by check and column.
- **Evidence:** `qpx/writers/base.py:414-421`; `qpx/writers/base.py:571-576`; `qpx/writers/base.py:599-617`; `audit-run:infra/scripts/writer_exp2.py: i2.feature.parquet logged only peptide_qvalue`

### INFRA-16: Canonical peptidoform is not enforced; the semantic validator accepts non-canonical tags {#infra-16}

*cross* · diverges from spec · confirmed · checklist D5, P9

- **Now:** Writers never normalize or check peptidoform. Only DIA-NN and migrate call normalize_peptidoform (diann/constants.py:75; transforms/migrate.py:85). The opt-in semantic check counts any name without ':' as canonical and accepts any mass shift (semantic_validation.py:68-81). Accepted: 'PEPTM[Oxidation]K', 'PEPTM[+15.994915]K', 'PEPTM[+16]K'. Missing fixed mods can only be detected through the calculated_mz check.
- **Spec:** peptidoform.md:13-35 rule 1: fixed mods and labels inline; rule 2: a Unimod-resolvable modification is UNIMOD:n; rule 3: otherwise a signed mass with four decimals. producer-guide.md:52.
- **Reader impact:** The same molecule can get different strings and so different feature_id/psm_id. Joins on peptidoform across producers fail silently.
- **Suggested:** The writer applies normalize_peptidoform, or rejects strings it would change. The validator flags names and masses that resolve to Unimod and masses not written with four decimals.
- **Evidence:** `qpx/core/semantic_validation.py:19-23`; `qpx/core/semantic_validation.py:68-96`; `qpx/converters/ptm.py:584-608`; `qpx/converters/diann/constants.py:75`; `docs/spec/peptidoform.md:11-46`; `audit-run:infra/scripts/writer_exp2.py output`

### INFRA-17: Reference normalize_peptidoform emits non-canonical terminal forms {#infra-17}

*cross* · diverges from spec · partly confirmed (see correction) · checklist D5

- **Now:** _add_nterm_dash runs only when the input contains '(' (ptm.py:515-548). So normalize_peptidoform('[Acetyl]PEPTMK') returns '[UNIMOD:1]PEPTMK', which the semantic validator flags as malformed. 'PEPTMK[Amidated]' becomes 'PEPTMK[UNIMOD:2]', a C-terminal mod on the side chain of K, which the validator accepts.
- **Spec:** peptidoform.md:28-29 rule 4: '[UNIMOD:1]-PEPTIDE' and 'PEPTIDE-[UNIMOD:2]'. peptidoform.md:45-46 points producers to normalize_peptidoform.
- **Correction:** Both failing inputs are themselves non-standard: '[Acetyl]PEPTMK' lacks the ProForma hyphen, and 'K[Amidated]' places the mod on K. The helper keeps the input's placement. It is a robustness gap, not wrong output from canonical ProForma. peptidoform.md:45 points readers, not producers, to the helper; producer-guide.md:52 points producers.
- **Reader impact:** Producers who follow the documented helper still write strings that are not canonical, or that change the site, so identities do not match those of other producers.
- **Suggested:** Always insert the terminal hyphens. Move modifications that are terminal-only in Unimod to the terminus. Add test vectors for terminal cases.
- **Evidence:** `qpx/converters/ptm.py:515-559`; `qpx/converters/ptm.py:584-608`; `docs/spec/peptidoform.md:28-29`; `python: normalize_peptidoform('[Acetyl]PEPTMK') -> '[UNIMOD:1]PEPTMK' (peptidoform_problem='malformed')`

### INFRA-20: Nothing declares the filtering state of psm, feature or pg files {#infra-20}

*footer* · unspecified · confirmed · checklist D8, P1, P3, P4, F4, F8

- **Now:** BaseWriter stamps qpx_version, writer_version, creator, file_type, creation_date, uuid, software_provider, compression_format, scan_format, primary_key, identity_composite, identity_hash_version and identity_source (base.py:114-135,347-353,462-463). No key gives FDR level or threshold, rank policy, decoys, transferred or unassigned rows, or upstream filtering. XicWriter has such keys (xic_selection targets|targets_and_decoys|all, xic_qvalue_threshold; xic.py:9-31). Provenance holds only converter-specific CLI options (mzIdentML rank1_only, OpenMS include_unassigned_psms). CDAP, OpenMS and MaxQuant record no upstream FDR.
- **Spec:** psm.md:256-268 and feature.md:272-284 list no such key. provenance.md:3,59 promise 'what FDR was used?', but provenance.md:65 covers only qpxc options.
- **Reader impact:** A reader cannot tell whether a psm file holds all ranks, decoys and sub-threshold PSMs or a 1%-filtered subset. Counts, FDR re-estimation and MBR-aware rollups cannot be done safely.
- **Suggested:** Add mandatory footer keys for psm, feature and pg, following the XIC precedent: rank_scope, fdr_level, fdr_threshold (or none), decoys, transferred, unassigned and upstream_filtered. Each key takes values like included/excluded/unknown.
- **Evidence:** `qpx/writers/base.py:114-135`; `qpx/writers/xic.py:9-49`; `docs/spec/xic.md:21`; `docs/spec/provenance.md:3,59,65`; `run outputs: footers of audit-run:infra/{cdap,omslfq,mzid}/*.psm.parquet; provenance parameters`

### INFRA-21: Footer keys added after the first write are invisible to Arrow readers {#infra-21}

*footer* · ambiguous · confirmed · checklist D8

- **Now:** _add_footer_metadata appends KV pairs without re-serializing ARROW:schema (base.py:359-365); only _set_scan_format re-serializes (base.py:367-383). identity_source is always added at close, after the flush (base.py:459-463). It is in the raw Parquet KV and in duckdb, but missing from pq.read_table().schema.metadata and ParquetFile.schema_arrow.metadata (checked on CDAP, OpenMS and the test files). spectrum_content is visible only when the rows did not exceed batch_size before it was set.
- **Spec:** psm.md and feature.md give reading code for footers via `parquet_file.schema_arrow.metadata` (psm.md:270-276), which misses these keys. transforms/relayout.py:139 notes 'in the footer only'.
- **Reader impact:** Readers following the spec's own example do not see identity_source or spectrum_content. Whether they see them depends on file size.
- **Suggested:** Re-serialize ARROW:schema whenever footer keys are added. Make the spec say to read raw KV metadata (ParquetFile.metadata.metadata).
- **Evidence:** `qpx/writers/base.py:359-383`; `qpx/writers/base.py:454-463`; `qpx/converters/maxquant/psm_adapter.py:116-125`; `qpx/transforms/relayout.py:139`; `python check: raw identity_source=b'derived', schema_arrow=None`

### INFRA-22: Footer key tables in the spec are incomplete or list keys never written {#infra-22}

*footer* · contradictory · partly confirmed (see correction) · checklist D8

- **Now:** Writers emit writer_version, uuid, primary_key, identity_composite, identity_hash_version, identity_source and spectrum_content (MaxQuant and mzIdentML psm). No writer emits project_accession or project_title (grep finds no writer code).
- **Spec:** psm.md:256-268 and feature.md:272-284 list seven keys and omit the identity keys, writer_version, uuid and spectrum_content (psm.md:290 mentions it only in prose). serialization.md:279-296 lists project_accession and project_title, and identity_source appears only in identity.md:21.
- **Correction:** serialization.md:283-296 does list writer_version, primary_key and identity_composite. identity_source is also explained at serialization.md:72. The registry is fragmented, not missing these keys. uuid and identity_hash_version are absent from serialization.md.
- **Reader impact:** Readers cannot know which keys to expect, or which are normative versus informational.
- **Suggested:** Keep one normative footer-key registry page (key, views, required, values) and refer to it from every view page.
- **Evidence:** `docs/spec/psm.md:256-268`; `docs/spec/feature.md:272-284`; `docs/spec/serialization.md:279-296`; `docs/spec/identity.md:14-21`; `qpx/writers/base.py:114-135`

### INFRA-30: Shared helper silently drops rows that fail to transform {#infra-30}

*other* · unspecified · confirmed · checklist P3, D8

- **Now:** BaseConverter._transform_rows_with_skip catches ValueError, TypeError, KeyError, IndexError and ArithmeticError per row and only logs the skipped count (converters/base.py:211-240). PD (psm, feature, pg), FragPipe (psm, feature, pg) and CDAP use it. The skipped counts are not written to provenance or the footer.
- **Spec:** Nothing. provenance.md:65 records options only.
- **Reader impact:** An output can be a silently filtered subset of the source. Readers cannot tell conversion losses from upstream filtering.
- **Suggested:** Record skipped-row counts per view in provenance, or fail above a threshold.
- **Evidence:** `qpx/converters/base.py:211-240`; `qpx/converters/proteome_discoverer/psm_adapter.py:43`; `qpx/converters/fragpipe/psm_adapter.py:142`; `qpx/converters/cdap/base_adapter.py:297`

### INFRA-32: Conversion filter options are recorded even when unset; rule worth keeping {#infra-32}

*other* · specified, worth keeping · partly confirmed (see correction) · checklist P3, D8

- **Now:** The CLI postprocess records every option that shapes the output, and records qvalue_threshold even when unset (cli/convert/postprocess.py:36-37). mzIdentML provenance shows rank1_only=False and OpenMS shows include_unassigned_psms=True.
- **Spec:** provenance.md:65: options are recorded with their defaults, 'qvalue_threshold (none when unset, so "not filtered" is on record)'.
- **Correction:** The rule holds only for DIA-NN and Spectronaut. No DDA converter has a qvalue_threshold option, so a DDA provenance never records 'not filtered by qpx'. The absence of the key is the only signal.
- **Reader impact:** Correct for qpx's own filtering. It covers neither the upstream tool's filtering nor the level a threshold applies at.
- **Suggested:** Keep the rule. Extend it with a standard vocabulary (fdr_level, fdr_threshold, rank_scope, decoys) that both qpxc and upstream tools fill, and mirror it in the footer (INFRA-20).
- **Evidence:** `qpx/cli/convert/postprocess.py:20-37`; `docs/spec/provenance.md:65`; `run: mzid provenance parameters rank1_only False`

### INFRA-M1: confident() helpers mix q-value levels and drop rows that lack typed q-values {#infra-m1}

*cross* · inconsistent across producers · found by the verifier · checklist D4, F4, F5, P3

- **Now:** feature.confident() = is_decoy=false AND COALESCE(precursor_qvalue, peptide_qvalue)<=t (qpx/core/data/feature.py:52-54). psm.confident() = qvalue<=t, so nulls are dropped (psm.py:23-25). Runs: MaxQuant simple 0/9,301 features and 0/1,000 PSMs pass. CDAP 86/94 pass on best-PSM q-values. OpenMS 68/68 pass on a PSM-level q-value. MBR rows without a q-value always drop.
- **Spec:** scores.md:38 calls these helpers 'the portable way to filter'. No page defines the level precedence or the null policy.
- **Reader impact:** One call gives an empty result for upstream-filtered MaxQuant/FragPipe/quantms data, a PSM-level filter for OpenMS and a precursor filter for DIA. Transferred features vanish silently.
- **Suggested:** One helper per level, no COALESCE across levels. When no typed q-value exists, raise or return the rows with an explicit 'unknown' status. Document this in scores.md.
- **Evidence:** `qpx/core/data/feature.py:52-54`; `qpx/core/data/psm.py:23-25`; `docs/spec/scores.md:38`; `run: work/verify-infra/mq (maxquant_simple): feature identification_type direct 7832 / transferred 1469, 0 pass; psm qvalue 0/1000`

## Specification text

### SPEC-10: Typed q-value fill rules copy PSM-level values into feature.precursor_qvalue {#spec-10}

*cross* · contradictory · partly confirmed (see correction) · checklist D4, P5

- **Now:** CDAP feature output (run): precursor_qvalue equals the best PSM's psm_level_qvalue (e.g. 0.00162). peptide_qvalue is filled from peptide_level_qvalue on 94/94 rows. In a writer test, additional_scores 'qvalue' ends up in precursor_qvalue.
- **Spec:** scores.md:44: feature.precursor_qvalue is filled from precursor_qvalue, qvalue, q_value, psm_level_qvalue. scores.md:42: psm.qvalue uses the same names. converter-coverage.md:181 calls CDAP peptide_qvalue 'Not implemented'. tool-mappings.md:45 shows '---'.
- **Correction:** The core holds. In CDAP output, precursor_qvalue = the best PSM's psm_level_qvalue (94/94), and peptide_qvalue is filled 94/94 although converter-coverage.md:181 says 'Not implemented'. tool-mappings.md:45 has no CDAP column, so it is not evidence.
- **Reader impact:** A 'run-level precursor q-value' can really be one PSM's q-value. The coverage docs misstate which columns are filled.
- **Suggested:** Fill precursor_qvalue only from precursor-level names, and leave PSM-level q-values in additional_scores. Fix the CDAP coverage row.
- **Evidence:** `docs/spec/scores.md:38-47`; `docs/spec/converter-coverage.md:181`; `docs/spec/tool-mappings.md:45`; `qpx/converters/cdap/feature_adapter.py:288-294`; `audit-run:spec/cdap/cdap.feature.parquet`

### SPEC-13: 'unique' means protein-unique or group-unique depending on view and producer {#spec-13}

*cross* · inconsistent across producers · partly confirmed (see correction) · checklist F9, X1

- **Now:** maxquant_full (run): 320 peptidoforms have pepmap.is_unique=true but feature.unique=false. DIA-NN: true if Protein.Group and Protein.Ids name one protein (converter-coverage.md:79). FragPipe: Protein+Mapped Proteins name one protein (:219). MaxQuant fallback: single group of any size (tool-mappings.md:62). CDAP: len(proteins)==1 (cdap/feature_adapter.py:265).
- **Spec:** feature.yaml:159-161 and versioning.md:68-69: group-based (one-protein group, no other group). pepmap.yaml:23-25 / pepmap.md:21: maps to exactly one protein. pg.md:82,85: unique to the protein group. feature.md:302 says unique indicators live in pg.
- **Correction:** pepmap.is_unique (maps to one protein) and feature.unique (one-protein group) are defined differently on purpose. The 320 MaxQuant cases (e.g. maps only to P39935, group P39935;P39936) are correct under both. Real divergences: the MaxQuant fallback (any group size) and CDAP's group of all mapped proteins. feature.md:120 vs :302 disagree on where unique lives.
- **Reader impact:** A 'unique peptide' filter returns different sets in different views and producers, conflating proteotypic with group-specific.
- **Suggested:** Use separate names (protein_unique, group_unique), one normative rule each, and correct the converter fallbacks.
- **Evidence:** `qpx/core/data/schemas/feature.yaml:159-161`; `qpx/core/data/schemas/pepmap.yaml:23-25`; `docs/spec/versioning.md:68-69`; `docs/spec/pg.md:82,85`; `docs/spec/tool-mappings.md:62`; `docs/spec/converter-coverage.md:79,219`; `audit-run:spec/mq pepmap vs feature query`

### SPEC-14: Same name pg_accessions holds group membership or evidence proteins; join undefined {#spec-14}

*cross* · ambiguous · partly confirmed (see correction) · checklist F9, P8

- **Now:** pepmap.pg_accessions = all evidence proteins (pepmap.md:61); feature.pg_accessions = assigned razor group; pg.pg_accessions = group. psm.protein_accessions = all mapped proteins, but tool-mappings.md:27 uses DIA-NN Protein.Group, and its order (razor first?) is undefined. The softlink joins list<pg_protein> to list<string> via an undefined canonical(); the pg_group_id hash bytes are unspecified.
- **Spec:** feature.yaml:151-153 (group). pepmap.yaml:15-17 ('Protein accessions with optional positions'). data-model.md:159. feature.md:103 ('Hash of the sorted pg_accessions'). psm.md:57 vs :289 (calls protein_accessions an 'optional column'; psm.yaml:99 does not).
- **Correction:** The evidence vs group split is documented (pepmap.md:61). canonical() is defined by reference to pg_id's membership set (pg.md:71, identity.md:51). Still open: pg_group_id bytes (code hashes JSON [sorted accessions] without dedup), protein_accessions order, psm.md:57 'no' vs :289 'optional', the stale DIA-NN PSM column. CDAP feature.pg_accessions = all mapped proteins.
- **Reader impact:** Joining or filtering on pg_accessions in pepmap silently mixes evidence with groups, and third-party producers cannot reproduce pg_group_id.
- **Suggested:** Rename the pepmap field (evidence_proteins). Define canonical() and the pg_group_id bytes with a test vector, and define the order of protein_accessions.
- **Evidence:** `docs/spec/pepmap.md:21,61`; `qpx/core/data/schemas/pepmap.yaml:15-17`; `docs/spec/data-model.md:93,159`; `docs/spec/feature.md:103`; `docs/spec/psm.md:57,289`; `docs/spec/tool-mappings.md:27`

### SPEC-22: Primary intensity is called 'raw' but holds MaxLFQ, normalized or summed values {#spec-22}

*cross* · contradictory · partly confirmed (see correction) · checklist F6, G2

- **Now:** pg.intensity is MaxLFQ for DIA-NN 2.x and as a FragPipe fallback (pg.md:54-57; converter-coverage.md:107,301), and normalized and scaled for PD (:325). Feature LFQ: CDAP best-PSM PrecursorArea, PD largest Precursor Abundance. TMT is summed over PSMs (cdap/feature_adapter.py:4-12; converter-coverage.md:324). Feature has no quantification_method.
- **Spec:** intensities.md:19-22: 'a data consumer can always find the raw measurement in intensities'. pg.md:110: 'Primary raw intensity'. intensities.md:213: processing keys are optional, and only in provenance.
- **Correction:** pg already names the method where it is known: cv_params quantification_method (pg.md:56-60), written only by DIA-NN and OpenMS consensus, and merge() warns on a mismatch. Feature has nothing equivalent. The contradiction between intensities.md:19-22 and pg.md:53-60,110 holds.
- **Reader impact:** The same column holds area, apex, sum, MaxLFQ or normalized values, so comparisons across producers are invalid.
- **Suggested:** Add a typed quantification_method (or footer key) on feature and pg, and drop the 'raw' guarantee or forbid derived values in intensity.
- **Evidence:** `docs/spec/intensities.md:19-22,213`; `docs/spec/pg.md:54-60,110`; `docs/spec/converter-coverage.md:107,301,324-325`; `qpx/converters/cdap/feature_adapter.py:4-12`

### SPEC-27: Decoy status parsed from accession prefix despite 'do not parse' rule {#spec-27}

*cross* · contradictory · partly confirmed (see correction) · checklist P4, G3, F8

- **Now:** DIA-NN without a Decoy column sets is_decoy when every accession has a DECOY_/decoy_/rev_/REV_ prefix (pg.md:336). Feature is_decoy for DIA-NN and MSstats comes from the accession prefix (tool-mappings.md:44).
- **Spec:** pg.md:46-49: 'is_decoy and contaminant state the status; do not parse it from the prefix'. producer-guide.md:28-35.
- **Correction:** pg.md:46-49 tells readers to rely on the flags. It does not forbid converters from inferring them. The real issue is that inference rules differ: prefix sets (mzIdentML lacks REV_, OpenMS adds RANDOM_) and any vs all (SPEC-M3). REV__ vs rev_ is not a mismatch, because startswith('REV_') matches REV__.
- **Reader impact:** Prefix vocabularies differ (REV__ vs rev_), so the same database yields different flags by producer.
- **Suggested:** Let producers derive flags from a declared prefix recorded in provenance. Readers rely on the flags only.
- **Evidence:** `docs/spec/pg.md:46-49,336`; `docs/spec/tool-mappings.md:44`; `docs/spec/producer-guide.md:28-35`

### SPEC-28: Decoy retention and contaminant flag unspecified for psm/feature {#spec-28}

*cross* · unspecified · confirmed · checklist P4, F8, P8

- **Now:** Decoy pg rows are kept (pg.md:99-103). For psm and feature nothing is said. pepmap drops decoys (converter-coverage.md:84-85). FragPipe input is pre-filtered (pg.md:365). psm and feature have no contaminant field.
- **Spec:** psm.md:45 and feature.md:60 define is_decoy only. producer-guide.md:34: the contaminant flag exists only on pg, and the marker stays on the accession (which pg.md:48 says not to parse).
- **Reader impact:** Decoy features can enter quantification, and contaminant peptides cannot be excluded without the prefix parsing the spec forbids.
- **Suggested:** State the decoy retention policy per view (footer flag). Add a contaminant bool to psm and feature, or bless a declared prefix.
- **Evidence:** `docs/spec/pg.md:48,99-103,365`; `docs/spec/psm.md:45`; `docs/spec/feature.md:60`; `docs/spec/producer-guide.md:34`; `docs/spec/converter-coverage.md:84-85`

### SPEC-29: Same quantity under different score/intensity names across producers {#spec-29}

*cross* · inconsistent across producers · confirmed · checklist G2, G3, P5

- **Now:** Protein-group PEP is named posterior_error_probability (MaxQuant, pg.md:258,322), pg_pep (DIA-NN, pg.md:356) or posterior_probability with higher_better=true (OpenMS, converter-coverage.md:109). MaxQuant LFQ maps to 'lfq' (pg.md:307) though the 'maxlfq' row lists MaxQuant (pg.md:228). CDAP writes msgf_specevalue/msgf_rawscore (run).
- **Spec:** scores.md:55: posterior_error_probability is a top-level PSM field, not an additional_scores entry. scores.md:62-63 recommend msgf_spec_evalue and msgf_raw_score. scores.md:80-81: the table lists names producers actually write.
- **Reader impact:** Readers need per-producer synonym tables to find the same quantity.
- **Suggested:** One name per quantity and level, a registry with level and direction, and a validator warning on synonyms.
- **Evidence:** `docs/spec/pg.md:228-229,258,264,307,322,356`; `docs/spec/scores.md:55,62-63,68,80-81`; `docs/spec/converter-coverage.md:109`; `audit-run:spec/cdap score names`

### SPEC-30: scores.md contradicts types.yaml and pg.yaml (cv_value null, PG cv_params, per-protein scores) {#spec-30}

*other* · contradictory · confirmed · checklist X1

- **Now:** The writer rejects the scores.md example null cv_value ('declared non-nullable', run). pg.yaml:118-120 defines pg.cv_params. score_value is a scalar float64 (types.yaml:4-8).
- **Spec:** scores.md:98,107,114: cv_value is nullable and null marks a flag. scores.md:123: PG has no cv_params. scores.md:87: per-protein scores use a 'values array index' matching pg_accessions. data-model.md:311: cv_value must not be null.
- **Reader impact:** The examples fail, and per-protein scores cannot be encoded as described.
- **Suggested:** Make cv_value nullable or fix scores.md, fix the PG row, and define a per-protein score encoding.
- **Evidence:** `docs/spec/scores.md:87,98,107,114,123`; `qpx/core/data/schemas/types.yaml:4-8,15-18`; `qpx/core/data/schemas/pg.yaml:118-120`; `docs/spec/data-model.md:311`; `audit-run:spec/ex_test.py (null cv_value REJECTED)`

### SPEC-31: Prose marks optional 1.3 columns as Required 'no' (= must be present) {#spec-31}

*other* · contradictory · confirmed · checklist D4, D7

- **Now:** The YAML marks these optional: psm.qvalue and ccs (psm.yaml:62-67,90-93); feature precursor/global_qvalue (feature.yaml:49-62), rt_aligned and irt (106-115), ccs (124), pg_group_id (162), identification_type, match_*, imputed_labels, PIF and quant_* (210-250); pg.pg_group_id and is_imputed (pg.yaml:20,113).
- **Spec:** psm.md:51,68; feature.md:58-59,66-67,73,103,145-148,154-156; pg.md:34,172 mark them 'no'. data-model.md:302-307 defines 'no' as 'the column must be present'.
- **Reader impact:** Validators following the prose reject valid files, and producers cannot tell what is mandatory.
- **Suggested:** Mark them 'optional' and generate the field tables from the YAML.
- **Evidence:** `docs/spec/data-model.md:302-307`; `docs/spec/feature.md:58-59,66-67,73,103,145-156`; `docs/spec/psm.md:51,68`; `docs/spec/pg.md:34,109,172`; `qpx/core/data/schemas/feature.yaml:49-62,106-127,162,210-250`

### SPEC-32: 1.3 changelog omits new semantic columns; matrix claims 'same columns' {#spec-32}

*other* · contradictory · confirmed · checklist D1, D4

- **Now:** Pages tag columns as 'specification 1.3': psm.qvalue, feature precursor/global_qvalue, rt_aligned, irt, channel_details, representative_label, pg_group_id (feature.md:58-67,86-87,103; psm.md:51). identification_type, match_*, imputed_labels, is_imputed, PIF and quant_scan are new as well.
- **Spec:** The versioning.md:55-97 changelog lists none of the typed q-values, identification_type or match_* columns. versioning.md:208-209: 1.2 and 1.3 have the 'Same columns'.
- **Reader impact:** Readers cannot learn from the changelog when MBR flags or typed q-values became available or what they mean.
- **Suggested:** Complete the 1.3 changelog and correct the compatibility matrix.
- **Evidence:** `docs/spec/versioning.md:55-97,208-209`; `docs/spec/feature.md:58-67,86-87,103`; `docs/spec/psm.md:51`

### SPEC-33: tool-mappings.md and ontology.md are stale against the schemas {#spec-33}

*other* · contradictory · confirmed · checklist X1

- **Now:** tool-mappings.md:17-30 maps DIA-NN PSM fields, though DIA-NN writes no PSM view (converter-coverage.md:12,130). It lists pg run_file_name and intensities (:105,109), DIA-NN feature_counts=Precursor.Quantity (:107), and FragPipe charge '---' (:43). ontology.md lists diann_qvalue/diann_global_qvalue (:290-291), pg 'intensities' (:247), and contaminant '1 = contaminant' (:243).
- **Spec:** pg.yaml:46-73 (grouped_runs, label, intensity), contaminant is a bool (pg.yaml:83). scores.md:80-81: no converter emits diann_qvalue. versioning.md:134-145 removed pg run_file_name.
- **Reader impact:** Producers and readers using these pages write or look for removed or wrong fields and terms.
- **Suggested:** Regenerate both pages from the YAML and converter-coverage, or delete the tool tables.
- **Evidence:** `docs/spec/tool-mappings.md:17-30,43,105-111`; `docs/spec/ontology.md:243,247,277,290-291`; `docs/spec/converter-coverage.md:12,130`; `docs/spec/scores.md:80-81`; `docs/spec/versioning.md:134-145`

### SPEC-34: views.md contradicts itself and data-model on analyte key and labels {#spec-34}

*other* · contradictory · confirmed · checklist X1

- **Now:** views.md has two 'Abundance matrices' sections. Protein analyte: sorted membership (:153) vs anchor_protein (:212). Unit name: first run (:157) vs sorted runs plus label (:212). Default agg: mean (:162) vs raises ValueError (:214). views.md:96-97: in a single-sample run, a feature label matching no sample (DIA-NN 'raw') is taken as that sample.
- **Spec:** data-model.md:239-246: labels must match exactly, 'This holds for a run with a single sample too'. pg.md:66: never treat anchor_protein as the analyte. versioning.md:177: 'raw' was non-conformant.
- **Reader impact:** The documented API returns different analyte keys, and silently accepts mislabelled quantities.
- **Suggested:** Delete the duplicate section, key proteins on pg_group_id or membership, and drop the label fallback.
- **Evidence:** `docs/spec/views.md:96-97,147-168,210-214`; `docs/spec/data-model.md:239-246`; `docs/spec/pg.md:66`; `docs/spec/versioning.md:177`

### SPEC-35: peptidoform.md cites a nonexistent Peptide view and a universal key {#spec-35}

*other* · contradictory · confirmed · checklist P9

- **Now:** peptidoform.md:107 lists 'Peptide (peptide_file)'. :109 says peptidoform keys with charge and run_file_name 'in all three views'. pepmap keys on peptidoform alone, and FragPipe features use quantification_unit_id (identity.md:34).
- **Spec:** versioning.md:28 and views.md:3: peptide summaries are API views, not files. psm.md:286 also mentions 'peptide views'.
- **Reader impact:** Readers look for a peptide file and assume a natural key that does not hold for every producer.
- **Suggested:** Remove the Peptide view rows and state the actual composite per view.
- **Evidence:** `docs/spec/peptidoform.md:107-109`; `docs/spec/versioning.md:28`; `docs/spec/views.md:3`; `docs/spec/identity.md:34`; `docs/spec/psm.md:286`

### SPEC-36: Footer tables incomplete; no keys for grain, rank, filter or quant method {#spec-36}

*footer* · unspecified · partly confirmed (see correction) · checklist D8, P3, P1

- **Now:** MaxQuant and CDAP feature footers have no identity_source (run), while the CDAP psm footer has one. The MaxQuant feature footer has no scan_format. spectrum_content appears only in psm.md:290.
- **Spec:** psm.md:258-268, feature.md:274-284 and serialization.md:283-296 omit spectrum_content, identity_source and identity_hash_version. identity.md:21: identity_source is absent only before 1.3. producer-guide.md:23: set scan_format.
- **Correction:** Refuted part: the MaxQuant and CDAP 1.3 feature, psm and pg footers do carry identity_source=derived and identity_hash_version=1. True part: scan_format is missing from the MaxQuant feature/psm and CDAP psm/feature footers despite producer-guide.md:23, and the footer tables omit identity_*, spectrum_content and uuid. No filter, rank or grain key exists.
- **Reader impact:** No machine-readable statement of filtering, rank policy, row grain or intensity method, and the documented identity keys are missing in real 1.3 files.
- **Suggested:** One normative footer-key table per view, with filter, rank-policy, grain and quant-method keys, checked by --conformance.
- **Evidence:** `docs/spec/psm.md:258-268,290`; `docs/spec/feature.md:274-284`; `docs/spec/serialization.md:283-296`; `docs/spec/identity.md:21`; `docs/spec/producer-guide.md:23`; `audit-run:spec/mq & cdap footers`

### SPEC-37: Canonical peptidoform rule: keep as contract core {#spec-37}

*cross* · specified, worth keeping · confirmed · checklist P9, D5

- **Now:** The rule is stated once and the same text recurs in producer-guide.md:50-52 and versioning.md:56-64. Compliance per converter was not audited here.
- **Spec:** peptidoform.md:13-41: every modification inline (fixed, isobaric, MS1 label), UNIMOD:n in upper case, signed 4-decimal mass shifts, terminal and unknown-position forms. A multiplet feature carries the lightest partner. I/L ambiguity is not addressed.
- **Reader impact:** This is what lets peptidoform serve as the cross-view join key. Without it, the same molecule gets different ids.
- **Suggested:** Keep it. Add an I/L rule and a validator check of the canonical form under --conformance.
- **Evidence:** `docs/spec/peptidoform.md:13-48`; `docs/spec/producer-guide.md:50-52`; `docs/spec/versioning.md:56-64`; `docs/spec/intensities.md:89-92`

### SPEC-M6: producer-guide says to flag 'imputation', but no such field exists {#spec-m6}

*other* · ambiguous · found by the verifier · checklist D6, D7

- **Now:** The real fields are feature.imputed_labels and pg.is_imputed (intensities.md:131-133; feature.yaml imputed_labels; pg.yaml:111-115). Only Spectronaut fills them.
- **Spec:** producer-guide.md:48: 'Imputed by the tool | The value, with `imputation` flagged'.
- **Reader impact:** Producers following the guide invent a column or cv_param, and imputed values go unmarked.
- **Suggested:** Name imputed_labels (feature) and is_imputed (pg) explicitly in producer-guide.
- **Evidence:** `docs/spec/producer-guide.md:48`; `docs/spec/intensities.md:131-133`; `qpx/core/data/schemas/pg.yaml:111-115`; `docs/spec/converter-coverage.md:201`

## Validator and read API

### VALIDATOR-01: What counts as valid depends on which entry point checks the file {#validator-01}

*cross* · ambiguous · partly confirmed (see correction) · checklist X1, P6

- **Now:** There are three tiers. (1) Writers call validate_full(strict=False) (writers/base.py:416) and only warn on duplicate PKs (base.py:557-569). Converters never run the cross-view checks (Dataset.validate is not called under qpx/converters). (2) Dataset.validate() defaults to lenient (dataset_validation.py:131). (3) The CLI always passes strict=True (cli/validate.py:136,155), and --conformance also fails on any warning (validate.py:125-126). Severity is also mixed: label and run checks are hard errors even in lenient mode, while nulls and duplicate PKs are warnings.
- **Spec:** versioning.md:120-127 makes a duplicate PK a warning on the write path and an error in qpxc validate. pg.md:75 does the same for anchor_without_membership. producer-guide.md:60 says each writer 'validates on close()'.
- **Correction:** The three tiers are specified, not ambiguous: versioning.md:120-127 and pg.md:75 define them. The writer tier is weaker than claimed. Lenient validate_full makes range and vocabulary violations warnings, so FeatureWriter persisted PEP=1.5, identification_type='MBR' and intensity -5 (tested). A null in a required column does raise in the writer. The real gap: no footer records which tier a file passed.
- **Reader impact:** A published file can be 'written and validated' and still be INVALID. Readers cannot tell which tier a file passed.
- **Suggested:** The contract should define one conformance level and record it in the footer (e.g. qpx_conformance=strict|lenient|none).
- **Evidence:** `qpx/writers/base.py:416,557-569`; `qpx/core/dataset_validation.py:131`; `qpx/cli/validate.py:125-126,136,155`; `docs/spec/versioning.md:120-127`

### VALIDATOR-08: Typed q-value columns are not checked against additional_scores or level {#validator-08}

*cross* · unspecified · confirmed · checklist D4, P5

- **Now:** The validator only checks that typed q-values lie in [0,1]. A row with psm.qvalue=0.001 and additional_scores psm_qvalue=0.9 passes --conformance, as does feature.precursor_qvalue=0.5 alongside a precursor_qvalue score of 0.001. A peptide_qvalue with no source score also passes. Nothing checks score values or higher_better (which is nullable in types.yaml:8).
- **Spec:** scores.md:38-47 says typed columns are copied from additional_scores and that the first matching name wins.
- **Reader impact:** The flat column and the nested score can disagree, so two readers that pick different fields get different FDRs.
- **Suggested:** Validate that a typed column equals its source score when both are present, or state that the typed column takes precedence and drop the duplicate (GitHub #377).
- **Evidence:** `experiment mutate.py cases p_qvalue_vs_scores_conflict, f_precursor_q_gt_global_vs_scores, f_peptide_q_set_no_source -> exit 0`; `qpx/core/data/schemas/types.yaml:4-8`

### VALIDATOR-15: psm.feature_id is checked for existence only, not for consistency {#validator-15}

*cross* · unspecified · partly confirmed (see correction) · checklist P10, F2

- **Now:** Only dangling psm.feature_id values are reported (label_validation.py:272-278). A PSM linked to a feature in another run, or with another peptidoform, passes. In the maxquant/simple output, all 95 linked PSMs point to features with a different peptidoform and run (57 also differ in charge). That dataset fails only on a dangling run_file_name.
- **Spec:** psm.md:75 calls feature_id the authoritative link. intensities.md:102-105 allows a different peptidoform only for MS1-label partners. psm.md:158-161 ranks candidates within the PSM's own run. Nothing makes these machine rules.
- **Correction:** More than existence is checked. psm.feature_id is cross-checked against feature.psm_ids in both directions (label_validation.py:288-305), and psm.label against identified_labels (:259-266). Run, peptidoform and charge consistency is not checked. maxquant/simple: 95/95 linked PSMs differ in peptidoform and run, 57 also in charge.
- **Reader impact:** Spectral counts and best-PSM lookups through link_feature_psm() can attach the wrong spectra to a feature.
- **Suggested:** Require the same run_file_name (or id_run_file_name), the same charge, and the same peptidoform unless an MS1 label differs. Report violations as errors.
- **Evidence:** `qpx/core/label_validation.py:272-278`; `experiment mutate2.py x_psm_feature_id_other_run, x_psm_feature_id_other_peptidoform -> exit 0`; `query audit-run:maxquant/simple psm JOIN feature USING(feature_id): 95 rows, 95 peptidoform mismatch, run HeLa_01 vs MFA381`

### VALIDATOR-21: Measured zero and missing values are handled differently on each read path {#validator-21}

*cross* · inconsistent across producers · partly confirmed (see correction) · checklist F6, D6, G1

- **Now:** The validator accepts both 0 and null and cannot detect zero sentinels. MuData stores an explicit 0 as a cell, but X.toarray() turns absent cells into 0 too (mudata.py:103-113). qpxc transform quantify drops measured zeros (quantify.py:269). A pg view whose intensities are all 0 gets a no_pg_quantity warning ('source carried no protein quantity', schema.py:297-313), which fails --conformance. A feature with only a null intensity yields no precursors modality at all.
- **Spec:** intensities.md:121-129: null means identified but not quantified, and 0 means measured zero. anndata.md:101 documents the sparse convention.
- **Correction:** MuData follows the documented sparse convention: explicit 0 is stored and absent means not measured (anndata.md:99-101). That is not a divergence. The real divergences: quantify.py:269 drops measured zeros, and no_pg_quantity (schema.py:306) treats an all-zero pg as 'no quantity', which fails conformance. A null-only feature giving no precursors modality matches 'not measured'.
- **Reader impact:** A true zero becomes missing in one tool and stays zero in another, so downstream statistics depend on which helper was used.
- **Suggested:** State that consumers must keep 0. Have quantify keep zeros (>=0), and do not count measured zeros as 'no quantity'.
- **Evidence:** `qpx/mudata.py:103-113`; `qpx/cli/transform/quantify.py:269`; `qpx/core/data/schema.py:297-313`; `experiment mutate.py g_measured_zero -> WARN no_pg_quantity, exit 1; MuData f_measured_zero X.data=[0.]`

### VALIDATOR-26: views.md contradicts itself and the code on design_matrix and unit naming {#validator-26}

*other* · contradictory · confirmed · checklist X1

- **Now:** In code, design_matrix with agg=None raises when a cell has several values (dataset_queries.py:524-525,538-558), and units are named by their sorted runs joined with ';' plus '|label' (sample_sql.py:19-37).
- **Spec:** views.md:158-160 says units are 'named after its first run' and that agg defaults to 'mean'. views.md:214 says the default raises. views.md:212 calls the protein feature_id the anchor_protein, while views.md:153 says it is the membership.
- **Reader impact:** Readers following the docs expect averaging or anchor keys that the API does not give.
- **Suggested:** Merge the two 'Abundance matrices' sections of views.md into one that matches the code.
- **Evidence:** `docs/spec/views.md:150-166,210-214`; `qpx/core/dataset_queries.py:524-525,538-558`

### VALIDATOR-29: Without a run view, all label and sample checks are skipped silently {#validator-29}

*cross* · unspecified · confirmed · checklist F6

- **Now:** _check_run_resolution returns at once when self.run is None (dataset_validation.py:193-194). maxquant/simple_nosdrf (feature view, no run or sample) gets only a mass_error warning. Its feature labels are never resolved to samples, and no issue reports that run.parquet is missing.
- **Spec:** data-model.md:160 requires the cross-checks only 'when the run view is present'. producer-guide.md:14-15 does not require run or sample.
- **Reader impact:** A quantitative dataset can conform while none of its quantities map to a sample.
- **Suggested:** Require run.parquet whenever feature or pg carries intensities, or emit an error-level 'missing_run_view'.
- **Evidence:** `qpx/core/dataset_validation.py:193-194`; `run: qpxc validate audit-run:maxquant/simple_nosdrf --conformance -> only mass_error warning`; `docs/spec/data-model.md:160`

### VALIDATOR-30: Physical merge drops scan_format, spectrum_content and creator from footers {#validator-30}

*footer* · diverges from spec · confirmed · checklist D8, P7

- **Now:** The merge writer is created without scan_format (dataset_merge.py:265-267), and only XIC footer keys are carried over (:252-263). Merging a psm with scan_format=scan and one with scan_format=index gives merged.psm.parquet with no scan_format, creator 'qpx' and no spectrum_content.
- **Spec:** scan.md:80-81: a file without scan_format is read as holding scan numbers. psm.md:290 defines spectrum_content.
- **Reader impact:** After a merge, zero-based spectrum indices are silently read as scan numbers, which breaks spectrum lookup.
- **Suggested:** Carry footer semantics through the merge, writing one file per distinct scan_format and spectrum_content, or refuse to merge mixed declarations.
- **Evidence:** `qpx/core/dataset_merge.py:252-267`; experiment: DatasetCollection([a(scan), b(index)]).merge -> merged.psm.parquet footer {identity_composite, creator:'qpx', identity_source:'provided'}

### VALIDATOR-31: No footer or validated declaration of FDR level, rank, decoy or MBR policy {#validator-31}

*footer* · unspecified · confirmed · checklist D8, P3, P1, F4

- **Now:** footer_issues checks qpx_version, file_type, identity_composite and primary_key (footer_validation.py:77-92); scan_format is checked for psm only. spectrum_content and creator values are never checked. No key says whether decoys or sub-threshold PSMs are kept, the FDR level and threshold, the rank policy, whether MBR rows are included, or the intensity type.
- **Spec:** serialization.md:284-298 lists the footer keys; none of these semantics is among them. provenance.md:65 puts qvalue_threshold and rank1_only in provenance parameters, which are not validated.
- **Reader impact:** Readers must infer filtering state per producer, and two files with identical schemas can mean different things.
- **Suggested:** Add validated declarations such as fdr_level, fdr_threshold, rank_policy, decoys_included and transfers_included, per view.
- **Evidence:** `qpx/core/data/footer_validation.py:77-92`; `docs/spec/serialization.md:284-298`; `docs/spec/provenance.md:65`

### VALIDATOR-33: Rules the validator does enforce today, which the contract should keep {#validator-33}

*cross* · specified, worth keeping · confirmed · checklist P4, P5, F6, D1, D5, D6, G1

- **Now:** Enforced as errors under strict: required columns and types; non-null required fields (including is_decoy and nested score_value); PEP and q-values within [0,1]; intensity finite and >=0; identification_type vocabulary; duplicate PK and declared composite; a label listed once per feature row; each (run, label) resolving to exactly one sample; run_file_name and grouped_runs present in run; pg run-disjointness; anchor implying membership; under --semantic, ProForma syntax, sequence match, known UNIMOD and calculated_mz within 20 ppm.
- **Spec:** data-model.md:239-246,300-311; intensities.md:121-129; pg.md:75,416; peptidoform.md:11-30.
- **Reader impact:** These are the only machine-guaranteed invariants a reader can rely on today.
- **Suggested:** Keep these as contract level 1, and add the missing semantic rules above as level 2.
- **Evidence:** `qpx/core/data/schema.py:558-643`; `qpx/core/data/validation.py:70-74`; `qpx/core/label_validation.py:90-129,166-199`; `qpx/core/semantic_validation.py:156-172`; `tests/conformance/test_conformance.py:67-124`

### VALIDATOR-M1: Peptidoforms missing fixed mods or labels raise only a mass_error warning {#validator-m1}

*cross* · diverges from spec · found by the verifier · checklist P9, D5, P7

- **Now:** calculated_mz is checked against the peptidoform's own mass, so if both lack a fixed mod the row passes. mass_error is hard-coded 'warning' even under strict, and runs only with --semantic (semantic_validation.py:148-152). maxquant/simple: 1,497 Cys features, none with UNIMOD:4, and 1,485 rows off by k x 57.02 Da. simple_nosdrf: no TMT, 9,118 rows off by k x 229 Da. Plain `qpxc validate` passes both.
- **Spec:** peptidoform.md:17-19: every modification is written inline, fixed mods and isobaric labels included. The 20 ppm and 50 ppm checks appear nowhere in docs/spec.
- **Reader impact:** The same peptide appears with and without carbamidomethyl or TMT across producers, which splits ids, joins and aggregation. Only --conformance fails such a file, through a warning.
- **Suggested:** Make mass_error follow strict severity, flag offsets equal to fixed-mod or label masses as a specific error, and document the thresholds.
- **Evidence:** `qpx/core/semantic_validation.py:44-65,148-152`; `query audit-run:maxquant/simple feature: round((obs-calc)*z/1.00335) = 57 (1286 rows), 114 (167 rows)`; `query audit-run:maxquant/simple_nosdrf: 457, 228, ... (TMT 229.16 Da)`; `docs/spec/peptidoform.md:17-19`

### VALIDATOR-M4: Writers persist range, vocabulary and sign violations as warnings {#validator-m4}

*cross* · diverges from spec · found by the verifier · checklist X1, D1, D6, P5

- **Now:** The write path calls validate_full(strict=False) (writers/base.py:416), and schema.py:603 then makes range, vocabulary and intensity-sign checks warnings. FeatureWriter wrote files with PEP=1.5, identification_type='MBR' and intensity -5.0, logging only a warning. A null in a required column does raise.
- **Spec:** producer-guide.md:60: each writer validates on close() and replaces the destination only when nothing failed. versioning.md:120-127 names only duplicate PK as a write-path warning. feature.md:145 defines a closed vocabulary; intensities.md:127 says negatives are invalid.
- **Reader impact:** Files produced by qpx's own writers or converters can break closed vocabularies and [0,1] ranges, so readers cannot assume those invariants without re-validating.
- **Suggested:** Make range, vocabulary and finite/non-negative checks hard errors on the write path; keep only duplicate PK as a warning.
- **Evidence:** `qpx/writers/base.py:416`; `qpx/core/data/schema.py:580,603,626`; `experiment work/verify-validator/writer: pep_1.5, idtype_MBR, neg_intensity WRITTEN; null_charge RAISED`

### VALIDATOR-M6: Charge sign and zero are unchecked, and m/z checks skip them silently {#validator-m6}

*cross* · unspecified · found by the verifier · checklist P7, F3

- **Now:** psm.yaml:33-35 and feature.yaml:33-35 declare no charge range. Mutation f_negative_charge (charge -2) passes --conformance. The calculated_mz check skips rows with charge<=0 (semantic_validation.py:50). Unidentified features are written with charge 0 (audit-run:infra/unid).
- **Spec:** scope.md:19 notes that charge=0 validates. Nothing defines negative charge or what 0 means (unknown? unidentified?).
- **Reader impact:** Charge 0 or negative values can mean unknown, unidentified or an error. Readers grouping by precursor (peptidoform, charge) and m/z checks treat them silently.
- **Suggested:** Require charge>=1 for identified PSM and feature rows, and give unidentified features a null charge or an explicit flag.
- **Evidence:** `qpx/core/semantic_validation.py:50`; `qpx/core/data/schemas/feature.yaml:33-35`; `mut/f_negative_charge -> exit 0`; `docs/spec/scope.md:19`

## mokume (consumer)

### MOKUME-05: Several mokume paths use anchor_protein as protein identity, which the spec forbids {#mokume-05}

*cross* · diverges from spec · confirmed · checklist F9

- **Now:** The low-frequency filter counts by anchor_protein but removes by the joined group name, so it does nothing for multi-member groups (feature.py:1018-1030; the Rust doc admits this at lib.rs:769-777). DEqMS peptide counts are nunique(sequence) per anchor_protein over raw rows, including decoy, MBR and null-intensity rows (stages.py:1771-1776). TissueMap keys proteins on anchor_protein (loader.py:51).
- **Spec:** data-model.md:159: do not fall back to anchor_protein. pg.md:67: leaders are not unique, and consumers MUST NOT assume a protein is in only one group. feature.md:120.
- **Reader impact:** Low-frequency filtering and DEqMS counts land on the wrong group when leaders are shared. DEqMS variance priors come from unfiltered counts.
- **Suggested:** anchor_protein is for display only. Every group join and count uses pg_group_id or pg_accessions.
- **Evidence:** `mokume/python/mokume/io/feature.py:1018-1030`; `mokume/python/mokume/pipeline/stages.py:1771-1776`; `mokume/python/mokume/tissuemap/loader.py:51`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:769-777`; `docs/spec/data-model.md:159`; `docs/spec/pg.md:67`

### MOKUME-10: How per-run features combine into a sample quantity: mokume uses MAX, the spec is silent {#mokume-10}

*cross* · unspecified · partly confirmed (see correction) · checklist F1, F6

- **Now:** At sample level mokume takes the MAX over fractions and technical replicates per peptidoform and charge (aggregation.py:263-292 and 316-320; SQL ROW_NUMBER at stages.py:128-140; Rust lib.rs:6209-6211), then sums peptidoforms and charges per sequence. According to the leads, qpx's own `ds.intensity('peptide')` sums instead; this was not re-verified on dev.
- **Spec:** Features are per run (producer-guide.md:12). How per-run quantities of fractions or replicates combine into a sample or unit quantity is not specified. pg defines units only through grouped_runs (pg.md:414).
- **Correction:** The spec is not silent. views.md:212 documents qpx's own ds.intensity('peptide'): features are summed within a unit (charge states, and fractions = runs with the same run.samples). Code: dataset_queries.py:99-134 SUM. That is API behavior, not a consumer contract, so this is an inconsistency between readers (qpx SUM vs mokume MAX), not an unspecified point. qpx's SUM also adds direct+transferred duplicates.
- **Reader impact:** For fractionated samples, mokume undercounts a peptide split across adjacent fractions, and two QPX readers disagree.
- **Suggested:** Specify fraction aggregation (sum within the quantification unit), or explicitly leave it to consumers with a recommended default.
- **Evidence:** `mokume/python/mokume/preprocessing/aggregation.py:263-292`; `mokume/python/mokume/pipeline/stages.py:128-140`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:6209-6211`; `docs/spec/pg.md:414`

### MOKUME-12: Samples come from an external SDRF, not run.samples; without an SDRF, channels collapse {#mokume-12}

*run_sample* · diverges from spec · partly confirmed (see correction) · checklist X1

- **Now:** mokume maps (run, label) to a sample through the SDRF (feature.py:755-978; sdrf.rs:108-146). It never reads run.parquet or sample.parquet outside TissueMap. Without an SDRF, sample = run_file_name (feature.py:679-686; lib.rs:8872-8881). A probe on the OpenMS TMT10 output gave 1 sample with 10 labels, and Python then keeps the max over channels (aggregation.py:316-320). grep found no TMT-without-SDRF guard in Rust.
- **Spec:** producer-guide.md:14: run.samples joins a label in a run to a sample. run.md:133-141. The validator raises an error on an unresolved (run, label) pair (label_validation.py:89-121, 192-199).
- **Correction:** The run view is optional (data-model.md:160: 'When the run view is present'). qpxc convert ships the SDRF as {prefix}.sdrf.tsv (file-naming.md:71), so SDRF-based mapping is legitimate when that file exists. The real defects: run.samples is ignored when there is no SDRF, and TMT channels silently collapse into one sample per run (probe2 re-run: 1 sample with 10 TMT labels).
- **Reader impact:** A valid, self-contained QPX dataset cannot be used as-is. TMT without an SDRF silently produces channel-collapsed output.
- **Suggested:** run.samples is the authoritative (run, label) to sample map. Consumers must use it; an SDRF is an optional override.
- **Evidence:** `mokume/python/mokume/io/feature.py:674-753,755-978`; `mokume/rust/crates/mokume-io/src/sdrf.rs:108-146`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:8872-8881`; `run: audit-run:mokume/probe2.py -> IsobaricWorkflow samples ['TMTTenPlexMethod_test'] with 10 TMT labels`; `docs/spec/producer-guide.md:14`; `qpx/core/label_validation.py:89-121,192-199`

### MOKUME-13: Labels and runs are matched case-insensitively after normalization; the spec says exact {#mokume-13}

*run_sample* · diverges from spec · partly confirmed (see correction) · checklist P7

- **Now:** mokume lowercases labels, reads the NT= value, and maps raw, lfq and label-free to 'label free sample' (feature.py:63-79; sdrf.rs:484-501). It also lowercases run names and strips path and extensions such as .raw, .mzml, .d and .zip (feature.py:44-60; sdrf.rs:454-482).
- **Spec:** intensities.md:63 and producer-guide.md:39: labels are compared exactly, case included. run.md:118: run_file_name is the file stem. The validator joins exactly.
- **Correction:** The normalization exists to match QPX against an external SDRF, and qpx's own SDRF readers also normalize spellings (intensities.md:65). The divergence only bites when two QPX runs or labels differ by case or extension, which is an edge case. Inside QPX the matching is exact (intensities.md:63; label_validation.py:89-121 raises ambiguous_label_sample when the run view is present).
- **Reader impact:** Usually harmless. But mokume can merge runs or labels that QPX treats as distinct (e.g. ITRAQ114 and iTRAQ114), and accept files the validator rejects.
- **Suggested:** Keep exact matching in the contract. Consumers should normalize no further than the spec does.
- **Evidence:** `mokume/python/mokume/io/feature.py:44-79`; `mokume/rust/crates/mokume-io/src/sdrf.rs:454-501`; `docs/spec/intensities.md:63`; `docs/spec/producer-guide.md:39`; `docs/spec/run.md:118`

### MOKUME-23: Feature-view filtering state is undeclared; mokume applies no confidence filter by default {#mokume-23}

*footer* · unspecified · confirmed · checklist P3, D8

- **Now:** By default mokume applies no PEP, q-value or FDR filter to features (filter defaults at feature.py:185-197; FDR options are opt-in at features_to_proteins.rs:75-81). It reads neither provenance nor footer filtering keys.
- **Spec:** Conversion does not filter by default (GitHub #241). provenance.md:65 records qvalue_threshold only in provenance parameters. No footer key states a view's FDR level or threshold.
- **Reader impact:** Protein results include whatever the producer kept: MaxQuant filters at 1% upstream, other producers' views may be unfiltered.
- **Suggested:** Add a footer key (e.g. `qpx.filtering`) giving the level and threshold applied to each view. Consumers check it before using the data.
- **Evidence:** `mokume/python/mokume/io/feature.py:185-197`; `mokume/rust/crates/mokume-command/src/features_to_proteins.rs:75-81`; `docs/spec/provenance.md:65`; `github_history.md: #241`

### MOKUME-28: Mixture, technical replicate and empty channels are inferred from file and sample names {#mokume-28}

*run_sample* · diverges from spec · confirmed · checklist X1

- **Now:** mixture = split_part(sample, '_', 1) (feature.py:743-746, 960-962). Rust derives the plex from sample-name tokens (lib.rs:10497-10514). The technical replicate comes from the run-name suffix (aggregation.py:226-236; feature.py:1129-1152). Empty channels are recognized by the condition string 'Empty' (aggregation.py:204; stages.py:105, 115).
- **Spec:** run.md:23-24: samples[].mixture and channel_role (reference, carrier, empty), new in 1.3. run.md:133-141: technical_replicate and biological_replicate per sample. run.fraction also exists.
- **Reader impact:** IRS and plex grouping and run normalization depend on naming conventions. The 1.3 design fields go unused.
- **Suggested:** Consumers take plex, replicates, fraction and channel role from run and run.samples.
- **Evidence:** `mokume/python/mokume/io/feature.py:743-746,1129-1152`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:10497-10514`; `mokume/python/mokume/preprocessing/aggregation.py:204,226-236`; `docs/spec/run.md:23-24,133-141`

### MOKUME-29: TissueMap loader uses anchor identity, no filters and positional channel aliases {#mokume-29}

*other* · diverges from spec · confirmed · checklist F9, X1

- **Now:** protein = anchor_protein, null is_decoy kept, intensity > 0 (loader.py:42-59). The LFQ join is on run_file_name only, taking the first sample of each run (loader.py:105-115, 360). TMT channels get positional aliases 1..N from the order of run.samples (loader.py:203-241). GIS reference falls back to the labels TMT126 and TMT131 (loader.py:254-262). No unique, MBR or FDR handling; all features are summed per protein and sample.
- **Spec:** data-model.md:159: no anchor_protein identity. The order of run.samples is unspecified (run.md:133-141). Labels are compared exactly (intensities.md:63). Reference channels are marked by channel_role (run.md:23).
- **Reader impact:** A second mokume reader gives different semantics from features2proteins on the same file.
- **Suggested:** Use channel_role for reference channels and forbid positional label aliases. All consumers share one row-selection rule.
- **Evidence:** `mokume/python/mokume/tissuemap/loader.py:42-59,105-115,136-157,203-262,360`; `docs/spec/data-model.md:159`; `docs/spec/run.md:23,133-141`

### MOKUME-M4: converter-coverage.md is stale for fields mokume depends on (psm.feature_id, peptide_qvalue) {#mokume-m4}

*other* · diverges from spec · found by the verifier · checklist P10, D4

- **Now:** MaxQuant fills psm.feature_id from Evidence ID (maxquant/psm_adapter.py:232; 95/1000 PSMs in audit-run:maxquant/simple). The writer fills CDAP feature.peptide_qvalue from peptide_level_qvalue (writers/derived.py:30; 94/94 in audit-run:mzid_cdap/cdap).
- **Spec:** converter-coverage.md:285: MaxQuant psm feature_id 'Not implemented: the Evidence ID of msms.txt is not read'. converter-coverage.md:181: CDAP peptide_qvalue 'Not implemented'.
- **Reader impact:** A consumer that reads the coverage tables (e.g. to decide whether spectral counting or peptide FDR is possible per producer) draws the wrong conclusions.
- **Suggested:** Generate the coverage tables from converter outputs or tests, and include the writer-derived typed columns.
- **Evidence:** `qpx/converters/maxquant/psm_adapter.py:232`; `qpx/writers/derived.py:28-33`; `docs/spec/converter-coverage.md:181,285`; `duckdb audit-run:maxquant/simple/maxquant.psm.parquet: 95/1000 feature_id`; `duckdb audit-run:mzid_cdap/cdap: peptide_qvalue == peptide_level_qvalue 94/94`

## External tools and standards

### EXTERNAL-03: Standards declare the FDR level and threshold per file; QPX has no such declaration {#external-03}

*footer* · unspecified · partly confirmed (see correction) · checklist P3, D8

- **Now:** mzTab 1.0 MTD false_discovery_rate lists file FDRs at PSM, peptide and/or protein level (e.g. MS:1001364 pep:global FDR 0.01). mzIdentML requires one Threshold that MUST hold a child of MS:1002703. MaxQuant filters at mqpar peptideFdr (PSM), proteinFdr and siteFdr, default 0.01. The MaxQuant QPX run records no upstream FDR in provenance, and the psm footer has no filter key.
- **Spec:** The psm footer keys (psm.md:256-268) say nothing about filtering. provenance.md:65 records only qpxc's own qvalue_threshold ('none' when unset). The upstream tool's FDR filter goes into fdr_filtering steps only when a converter writes it (provenance.md:82,154-162).
- **Correction:** mzIdentML does not require a child of MS:1002703. One Threshold is required, and it may hold a PSM-level statistic, a peptide-level statistic, a generic threshold term, or MS:1001494 'no threshold'. The mqpar peptideFdr/proteinFdr/siteFdr defaults were not checked.
- **Reader impact:** A reader cannot tell whether a psm or feature file is unfiltered, 1% PSM-filtered, or 1% peptide- or protein-filtered by the upstream tool. That decides whether decoys and low-confidence rows are present and how to re-estimate FDR.
- **Suggested:** Mandatory footer or dataset keys: filter_applied_by (tool/qpx/none), fdr_level (psm/peptide/protein), fdr_threshold, using PSI-MS terms as mzTab does. Converters fill them from mqpar, mzIdentML Threshold, PD and FragPipe settings when known, else 'unknown'.
- **Evidence:** https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/specification_document-releases/1_0-Proteomics-Release/mzTab_format_specification.pdf §6.2.18; `https://peptideatlas.org/PSI/schemas/mzIdentML/1.2/mzIdentML1.2.0_combined.html (Threshold mapping rule)`; `https://storage.jpostdb.org/JPST001663/mqpar_Normal_IgG_SPACE.xml (peptideFdr/proteinFdr/siteFdr 0.01)`; `run: audit-run:external/mq/out/maxquant.provenance.parquet`; `docs/spec/provenance.md:65`; `docs/spec/psm.md:256-268`

### EXTERNAL-09: PSI-MS has one term per FDR level; QPX global_qvalue maps to three levels {#external-09}

*cross* · contradictory · confirmed · checklist D4, P5, F5

- **Now:** PSI-MS (OLS) distinguishes MS:1002354 PSM-level q-value, MS:1002350 PSM-level global FDR, MS:1001868 distinct peptide-level q-value ('once redundant identifications … removed'), MS:1001364 peptide sequence-level global FDR, MS:1002373 protein group-level q-value and MS:1001214 protein-level global FDR. MS:1001491 percolator:Q value does not state its level; the level depends on PercolatorAdapter -peptide_level_fdrs.
- **Spec:** ontology.md:277 maps global_qvalue to MS:1002350 (PSM-level). feature.md:59 and scores.md:45 call it the experiment-level precursor q-value. ontology.md:238 maps pg global_qvalue to MS:1001214. precursor_qvalue is filled from 'qvalue', 'q_value' and 'psm_level_qvalue' (scores.md:44).
- **Reader impact:** Global_qvalue means a different statistic in each view. A PSM-level q-value of the best PSM can become feature.precursor_qvalue, and readers filter at the wrong level.
- **Suggested:** One typed column per (level, scope), each bound to exactly one PSI-MS term. Remove PSM-level names from precursor_qvalue's fill list. Engine q-values without a declared level stay in additional_scores.
- **Evidence:** `https://www.ebi.ac.uk/ols4/api/ontologies/ms/terms (MS:1002354, MS:1002350, MS:1001868, MS:1001364, MS:1002373, MS:1001214, MS:1001491)`; `docs/spec/ontology.md:238,277`; `docs/spec/feature.md:59`; `docs/spec/scores.md:44-45`; `qpx/converters/openms_consensus/feature_adapter.py:36-40`

### EXTERNAL-28: QPX claims ProForma Base compliance but passes Base constructs through {#external-28}

*other* · contradictory · confirmed · checklist P9, D5

- **Now:** ProForma 2.0 Appendix I: Base-ProForma covers Unimod and PSI-MOD, delta masses, 'N-terminal, C-terminal and labile modifications', 'ambiguity in the modification position, including support for localisation scores', 'ambiguity in the amino acid sequence' and the INFO tag. Charge, chimeric '+' and global modifications belong to the 'mass spectrum compliant' extension.
- **Spec:** QPX reads 'the Base compliance level plus unknown-position tags' (scope.md:27), yet passes through labile {…}, localisation groups #g1(0.9) and ranges (scope.md:43-46). Unknown position '?' is already Base.
- **Reader impact:** Producers relying on 'Base compliant' may write labile or grouped-localisation peptidoforms. QPX keeps them as strings with no modifications parsed, so identity and mass differ by producer.
- **Suggested:** State the true subset ('Base minus labile, groups, ranges, INFO') or implement Base. Canonical form forbids the non-parsed constructs in peptidoform and moves them to modifications/modification_candidates.
- **Evidence:** `https://arxiv.org/pdf/2109.11352 (Appendix I)`; `docs/spec/scope.md:27,41-46`; `docs/spec/peptidoform.md:13-46`

### EXTERNAL-33: DDA vs DIA is declared only via an optional SDRF-derived run field {#external-33}

*cross* · unspecified · confirmed · checklist X1, D8

- **Now:** SDRF-Proteomics: comment[proteomics data acquisition method] (e.g. 'data-dependent acquisition'; DIA variants are PRIDE:0000659 descendants). comment[label] and comment[fraction identifier] are required per data file. Tools such as CHIMERYS apply one algorithm to DDA, DIA and PRM, and MSFragger-DDA+ treats wide-window DDA specially.
- **Spec:** run.acquisition_method is optional and comes only from SDRF (run.md:125; sdrf-mapping.md:43). The psm view is 'designed primarily for DDA' (psm.md:282-283). There is no dataset- or file-level acquisition declaration.
- **Reader impact:** A DDA contract cannot be applied when a reader cannot tell DDA files from DIA ones. Converters without SDRF leave acquisition unknown.
- **Suggested:** Require acquisition_method (PSI or PRIDE term) at dataset or run level for psm/feature, filled from SDRF or the converter's knowledge (MaxQuant/FragPipe DDA, DIA-NN DIA), including 'wide-window DDA'.
- **Evidence:** `https://proteomics-sample-metadata.readthedocs.io/en/latest/sdrf.html`; `https://pmc.ncbi.nlm.nih.gov/articles/PMC8494749`; `docs/spec/run.md:125`; `docs/spec/sdrf-mapping.md:43`; `docs/spec/psm.md:282-283`

## First-principles questions

### CRITIC-07: Is posterior_error_probability a PSM-only field or also a feature field? {#critic-07}

*cross* · contradictory · confirmed · checklist P5, F5

- **Now:** psm.md says PEP is defined only in the PSM view and that features do not carry it top-level. feature.md and feature.yaml nonetheless define feature.posterior_error_probability ('PEP for the peptide match'), which the MaxQuant (evidence PEP), FragPipe (best PSM) and Spectronaut (EG.PEP) converters fill.
- **Spec:** psm.md:285-286 vs feature.md:56; feature.yaml:37-41; scores.md:55; converter-coverage.md:180
- **Reader impact:** A reader cannot tell whether feature PEP is one PSM's PEP, an evidence/precursor PEP or a peptide PEP, and each differs by producer.
- **Suggested:** Keep the column and define it as the PEP of the identifying PSM in this run, null for transferred rows. Precursor/peptide PEPs go under named scores.
- **Evidence:** `docs/spec/psm.md:285-286`; `docs/spec/feature.md:56`; `qpx/core/data/schemas/feature.yaml:37-41`; `docs/spec/converter-coverage.md:180`

### CRITIC-08: Are decoy rows present in psm/feature, and how does a reader know? {#critic-08}

*cross* · unspecified · confirmed · checklist P4, F8, G3

- **Now:** is_decoy is required on all three views. pg explicitly may contain decoys and consumers MUST filter them; psm and feature say nothing about whether decoys are retained. DIA-NN emits decoys only with --report-decoys, and pepmap drops them. No footer key exists; #229 left this open.
- **Spec:** psm.md:45; feature.md:60; pg.md:99-103; pg.md:336; converter-coverage.md:84-85
- **Reader impact:** 'No decoys found' cannot be told from 'decoys removed'. FDR re-estimation is impossible or wrong, and counts include decoys unless filtered.
- **Suggested:** Footer decoys_retained (true|false|unknown) per view. State that consumers filter is_decoy in every view.
- **Evidence:** `docs/spec/psm.md:45`; `docs/spec/feature.md:60`; `docs/spec/pg.md:99-103`; `docs/spec/converter-coverage.md:84-85`

### CRITIC-09: Does is_decoy=false mean 'target' or 'no target-decoy search / unknown'? {#critic-09}

*cross* · ambiguous · confirmed · checklist P4

- **Now:** false covers both targets and searches without target-decoy (including de novo). A peptide matching target and decoy is a target. The column is non-nullable, and no footer says whether TDA was used.
- **Spec:** psm.md:17-19,45; feature.md:60; producer-guide.md:30-35
- **Reader impact:** A de novo or non-TDA dataset looks like an FDR-controlled all-target dataset.
- **Suggested:** Footer target_decoy (true|false), or make is_decoy nullable with null meaning no TDA.
- **Evidence:** `docs/spec/psm.md:45`; `docs/spec/producer-guide.md:30-35`

### CRITIC-10: How does a reader know a PSM or feature comes from a contaminant? {#critic-10}

*cross* · unspecified · confirmed · checklist P4, P8, F9

- **Now:** Only pg has `contaminant`. Accessions keep CONTAM_ markers, but the spec says not to parse prefixes. psm and feature have no flag. Shared or unassigned features have null pg_accessions, so joining to pg cannot recover it.
- **Spec:** producer-guide.md:34; pg.md:44-50,95
- **Reader impact:** Peptide-level analyses must either parse prefixes (forbidden and producer-specific) or keep contaminants.
- **Suggested:** Add a nullable bool contaminant to psm, feature and pepmap, with a stated any/all-member rule.
- **Evidence:** `docs/spec/producer-guide.md:34`; `docs/spec/pg.md:44-50`; `docs/spec/pg.md:95`

### CRITIC-13: Is observed_mz the selected-ion m/z, a corrected monoisotopic value, or a calibrated value? {#critic-13}

*cross* · ambiguous · confirmed · checklist P7

- **Now:** Ontology maps observed_mz to MS:1002234 'selected ion m/z'. MaxQuant writes calculated + calibrated ppm error, FragPipe and Spectronaut write calibrated m/z, and PD/Comet values carry +1..+3 Da isotope errors.
- **Spec:** psm.md:47; ontology.md:171; tool-mappings.md:26,47; converter-coverage.md:182,319
- **Reader impact:** mass_error distributions and m/z matching to the mz view differ by producer. Isotope-error PSMs look like errors of thousands of ppm.
- **Suggested:** observed_mz = instrument-selected precursor m/z. Add optional isotope_error and a calibrated value; declare calibration in the footer.
- **Evidence:** `docs/spec/ontology.md:171`; `docs/spec/tool-mappings.md:26`; `docs/spec/converter-coverage.md:182`; `docs/spec/converter-coverage.md:319`

### CRITIC-14: Is mass_error_ppm corrected for isotope errors and open-search mass shifts? {#critic-14}

*cross* · ambiguous · confirmed · checklist P7

- **Now:** Only the formula is given. For FragPipe, calculated_mz includes an unassigned mass shift (so the error says nothing about the instrument), and isotope errors stay in mass_error_ppm. MaxQuant MSMS rows are null.
- **Spec:** psm.md:48; converter-coverage.md:183,230-239,275
- **Reader impact:** Mass-accuracy QC and tolerance checks are not comparable across producers.
- **Suggested:** Define it against calculated_mz of the written peptidoform, before isotope correction, and add an int isotope_error column.
- **Evidence:** `docs/spec/psm.md:48`; `docs/spec/converter-coverage.md:230-239`

### CRITIC-15: Are peptidoform and calculated_mz written identically by all producers? {#critic-15}

*cross* · specified, worth keeping · confirmed · checklist P9, D5

- **Now:** Canonical form: every modification inline (fixed mods and labels included), UNIMOD:n in upper case, signed 4-decimal mass shifts, unknown-position prefix. calculated_mz is computed from that peptidoform. A multiplet feature carries its lightest partner's peptidoform.
- **Spec:** peptidoform.md:11-41; intensities.md:89-92; producer-guide.md:50-52
- **Reader impact:** Keep this: it is what makes peptidoform a cross-view and cross-producer join key.
- **Suggested:** Keep. Add I/L and ambiguous-residue rules, and forbid the passthrough constructs (scope.md:37-46) in DDA rows.
- **Evidence:** `docs/spec/peptidoform.md:11-41`; `docs/spec/intensities.md:89-92`; `docs/spec/scope.md:37-46`

### CRITIC-16: What does charge 0 mean, and can precursor charge be unknown? {#critic-16}

*cross* · unspecified · confirmed · checklist P7

- **Now:** charge is a required non-null int16. scope.md says charge=0 validates but gives it no meaning. Whether charge is the search-assigned or instrument charge is not stated.
- **Spec:** psm.yaml:33-36; scope.md:19; feature.md:55
- **Reader impact:** Unknown-charge spectra (e.g. MGF without CHARGE) can appear as 0 and break m/z math and joins.
- **Suggested:** charge = the charge assumed by the identification, >=1; forbid 0 in psm and feature.
- **Evidence:** `qpx/core/data/schemas/psm.yaml:33-36`; `docs/spec/scope.md:19`

### CRITIC-17: In what unit and time frame is rt in each view? {#critic-17}

*cross* · specified, worth keeping · confirmed · checklist P7, F10

- **Now:** psm.rt = MS2 scan RT in seconds. feature.rt is in seconds on the run's own clock, rt_aligned is the aligned time, and irt is on the tool's scale. mz.scan_start_time is in minutes (multiply by 60 to join).
- **Spec:** psm.md:50; feature.md:65-67; mz.md:30
- **Reader impact:** Keep. The minutes-vs-seconds split with the mz view is a trap the contract must spell out.
- **Suggested:** Keep, and extend the unit and frame statement to rt_start/rt_stop (see CRITIC-37).
- **Evidence:** `docs/spec/psm.md:50`; `docs/spec/feature.md:65-67`; `docs/spec/mz.md:30`

### CRITIC-18: What quantity and unit is ion_mobility (1/K0, drift time, or CV)? {#critic-18}

*cross* · ambiguous · confirmed · checklist P7, F10

- **Now:** psm and feature describe it only as an 'ion mobility value', with no unit. The ontology page maps it to 1/K0 (Vs/cm2). The mz view has drift_time in ms, and FAIMS CV is a separate column.
- **Spec:** psm.md:67; feature.md:72-76; ontology.md:173,205; mz.md:32
- **Reader impact:** Drift-tube or TWIMS producers can write ms into what readers treat as 1/K0, so values are not comparable across instruments.
- **Suggested:** State in the view pages that ion_mobility = 1/K0 in Vs/cm2. Other mobility types go in named columns or cv_params.
- **Evidence:** `docs/spec/psm.md:67`; `docs/spec/ontology.md:173`; `docs/spec/mz.md:32`

### CRITIC-35: Does global_qvalue mean the same level in every view? {#critic-35}

*cross* · contradictory · confirmed · checklist P5, D4

- **Now:** No. On PSMs it is an 'experiment-level PSM q-value' (ontology MS:1002350, PSM-level global FDR). The feature column global_qvalue is an experiment-level precursor q-value. pg.global_qvalue is a protein-group q-value (MS:1001214).
- **Spec:** psm.md:291; scores.md:45,56; feature.md:59; pg.md:92; ontology.md:238,277
- **Reader impact:** The same name changes meaning across views, so generic tooling mis-filters.
- **Suggested:** Use level-prefixed names everywhere: psm_qvalue, precursor_global_qvalue, pg_global_qvalue.
- **Evidence:** `docs/spec/psm.md:291`; `docs/spec/feature.md:59`; `docs/spec/pg.md:92`; `docs/spec/ontology.md:277`

### CRITIC-40: Are primary intensities raw, or possibly normalized or corrected? {#critic-40}

*cross* · contradictory · partly confirmed (see correction) · checklist F6, G2

- **Now:** intensities.md promises the raw measurement is always in intensities/intensity. Yet pg.intensity may be DIA-NN 2.x PG.MaxLFQ or PD normalized/scaled abundances. Processing state exists only as optional provenance keys, which only openms-consensus sets.
- **Spec:** intensities.md:17-22,211-213; pg.md:54-57,110; converter-coverage.md:325; provenance.md:69
- **Correction:** MaxQuant also writes the processing keys isotope_corrected=false and normalization=none (qpx/converters/maxquant/converter.py:122; seen in the run's provenance), so openms-consensus is not the only one. The raw-vs-normalized contradiction holds.
- **Reader impact:** Readers renormalize already-normalized values, or compare raw values with MaxLFQ.
- **Suggested:** Mandatory per-view keys: normalization, isotope_corrected and cross_run, each allowing 'unknown'.
- **Evidence:** `docs/spec/intensities.md:17-22`; `docs/spec/pg.md:54-57`; `docs/spec/converter-coverage.md:325`; `docs/spec/provenance.md:69`

### CRITIC-41: How do null, 0 and absent differ for primary intensities? {#critic-41}

*cross* · specified, worth keeping · confirmed · checklist D6, F6, G1

- **Now:** No entry = not identified, or no sample. An entry with null intensity = identified but not quantified. 0 = measured zero; tool zero-sentinels become null; NaN and negatives are invalid; imputed values are flagged. In pg, a group with identification evidence in a unit gets a row per label (null intensity if unquantified); other combinations get no row unless a quantity was measured (e.g. MBR).
- **Spec:** intensities.md:121-133; producer-guide.md:41-48; pg.md:116-146
- **Reader impact:** Keep: this is the core of telling missingness mechanisms apart.
- **Suggested:** Keep. Extend it to additional_intensities (CRITIC-42), and require producers to name their 0 sentinels.
- **Evidence:** `docs/spec/intensities.md:121-133`; `docs/spec/producer-guide.md:41-48`; `docs/spec/pg.md:116-146`

### CRITIC-42: How do null, 0 and absent differ in additional_intensities? {#critic-42}

*cross* · unspecified · confirmed · checklist D6, G2

- **Now:** intensity_value is non-nullable, and the spec has no rule. Converters drop 0 for maxlfq/lfq/ibaq (DIA-NN, MaxQuant) but keep spectral_count 0 on MaxQuant MBR rows.
- **Spec:** types.yaml:37-40; converter-coverage.md:60,108,195,302
- **Reader impact:** An absent maxlfq may mean not computed, zero or filtered, depending on producer.
- **Suggested:** Apply the primary rule: absent pair = not reported, 0 = measured zero, sentinel 0s dropped.
- **Evidence:** `qpx/core/data/schemas/types.yaml:37-40`; `docs/spec/converter-coverage.md:195`; `docs/spec/converter-coverage.md:302`

### CRITIC-43: Do intensity_name values follow one vocabulary, case and scale? {#critic-43}

*cross* · ambiguous · confirmed · checklist G2, X1

- **Now:** MaxQuant LFQ is named `lfq` while DIA-NN and FragPipe use `maxlfq` for the same algorithm. The pg example uses 'LFQ' and 'iBAQ'. log2 values are allowed (tmt_integrator_log2_abundance). The intensities example repeats label 'LFQ' twice in one list.
- **Spec:** pg.md:206-207,228-229,244; intensities.md:158-176,181
- **Reader impact:** Name matching misses the same quantity under two names, and linear and log values mix.
- **Suggested:** A controlled snake_case vocabulary with ontology mapping, encoding scale in the name (log2_ prefix). The validator warns on unknown names.
- **Evidence:** `docs/spec/pg.md:206-207`; `docs/spec/pg.md:228-229`; `docs/spec/pg.md:244`; `docs/spec/intensities.md:158-176`

### CRITIC-55: Is feature.unique the same notion as pepmap.is_unique? {#critic-55}

*cross* · ambiguous · confirmed · checklist F9

- **Now:** No. feature.unique means the peptide's group has exactly one protein; pepmap.is_unique means the peptide maps to exactly one protein, with no group involved. Two flags with similar names have different scopes.
- **Spec:** feature.md:101; pepmap.md:21,61
- **Reader impact:** Consumers swap the two and confuse protein-level with group-level uniqueness.
- **Suggested:** Rename to proteotypic (protein level) and group_unique, and document how they relate.
- **Evidence:** `docs/spec/feature.md:101`; `docs/spec/pepmap.md:21`

### CRITIC-57: Which peptides does pepmap cover: decoys, unfiltered, all ranks? {#critic-57}

*pepmap* · unspecified · confirmed · checklist P8

- **Now:** One row per peptidoform with its evidence proteins, and no is_decoy column. DIA-NN drops decoys. Whether it covers only identified peptides, all ranks, and how I/L is handled are unspecified.
- **Spec:** pepmap.md:1-34,61; converter-coverage.md:84-85
- **Reader impact:** Uniqueness computed from pepmap may include decoy or unfiltered peptides for some producers.
- **Suggested:** pepmap = target peptidoforms present in psm or feature; add is_decoy if decoys are kept.
- **Evidence:** `docs/spec/pepmap.md:21`; `docs/spec/converter-coverage.md:84-85`

### CRITIC-58: How does a quantity resolve to a sample? {#critic-58}

*run_sample* · specified, worth keeping · confirmed · checklist X1

- **Now:** (run_file_name, or any file in grouped_runs, plus label) resolves through run.samples[] to sample_accession. Labels are canonical and matched exactly. Validation errors on a label with no sample, or with several.
- **Spec:** data-model.md:222-246; intensities.md:61-77
- **Reader impact:** Keep: this makes TMT, SILAC and LFQ sample mapping independent of the producer.
- **Suggested:** Keep, and add rules for channel_role empty/carrier channels in feature and pg.
- **Evidence:** `docs/spec/data-model.md:239-246`; `docs/spec/intensities.md:61-77`

### CRITIC-59: Does a file state DDA vs DIA, and therefore which semantics apply? {#critic-59}

*footer* · unspecified · confirmed · checklist D8

- **Now:** Only the optional run.acquisition_method exists; footers carry only creator. scan and psm semantics differ between DDA and DIA.
- **Spec:** run.md:125; producer-guide.md:21-26; scan.md:109-111
- **Reader impact:** From a feature file alone, a reader cannot tell whether scan is the identifying MS2 or the apex scan.
- **Suggested:** Footer acquisition (DDA|DIA|PRM|mixed) on psm, feature and pg.
- **Evidence:** `docs/spec/run.md:125`; `docs/spec/producer-guide.md:21-26`

### CRITIC-60: Where is the processing state of feature/pg views declared (filters, MBR, decoys, rank)? {#critic-60}

*footer* · unspecified · confirmed · checklist D8, P3, F4

- **Now:** Nowhere in the footer. Provenance parameter keys are free text (psm_fdr, match_between_runs appear only as examples), and only qpxc's own options are recorded systematically. The xic and mz views do declare their selection in the footer.
- **Spec:** serialization.md:283-296; provenance.md:17-22,58-62,65,159-162; xic.md:21; mz.md:154
- **Reader impact:** Identical columns describe different populations (MBR on/off, protein-FDR filtered) with no machine-readable signal.
- **Suggested:** A controlled footer manifest: fdr levels and thresholds, mbr, decoys_retained, rank_policy, unidentified_features, normalization.
- **Evidence:** `docs/spec/serialization.md:283-296`; `docs/spec/provenance.md:159-162`; `docs/spec/xic.md:21`

### CRITIC-61: Which footer keys exist, and how must readers read them? {#critic-61}

*footer* · ambiguous · partly confirmed (see correction) · checklist D8

- **Now:** The serialization list omits identity_source, identity_hash_version, spectrum_content and qpx.struct_lists. The docs read schema_arrow.metadata, but CDAP output has identity_source only in the Parquet key-value footer. identity_composite is described both as comma-separated and as JSON.
- **Spec:** serialization.md:283-299; identity.md:14-21; psm.md:270-278; feature.md:286-294
- **Correction:** serialization.md:298-299 describes the JSON encoding of the composite *values* for hashing, not of the footer key, which is comma-separated (:292). That sub-claim is wrong. The rest is reproduced: the Arrow schema metadata lacks identity_source while the KV footer has it.
- **Reader impact:** Readers following the documented recipe miss identity_source and other late-written keys.
- **Suggested:** One normative footer-key registry. Read keys from ParquetFile.metadata.metadata.
- **Evidence:** `docs/spec/serialization.md:283-299`; `docs/spec/psm.md:270-278`; `audit-run:critic/cdap/cdap.feature.parquet (pq.read_schema lacks identity_source; ParquetFile.metadata.metadata has it)`

### CRITIC-62: Must the 1.3 typed columns be present, and do 1.2 and 1.3 have the same columns? {#critic-62}

*other* · contradictory · confirmed · checklist D1, D7

- **Now:** The md pages mark precursor_qvalue, global_qvalue, rt_aligned, irt, identification_type, match_* and others as 'No', which means the column must be present. The yaml marks them optional. versioning says 1.2 and 1.3 have the 'Same columns', and the changelog omits these additions.
- **Spec:** feature.md:58-59,66-67,145-156; psm.md:51; data-model.md:302-307; feature.yaml:49-62,210-254; versioning.md:208-209
- **Reader impact:** An absent column could mean pre-1.3 or 'unknown'; validators and docs disagree.
- **Suggested:** Mark these columns optional in the md pages, list them all in the changelog, and define absent vs null.
- **Evidence:** `docs/spec/feature.md:58-59`; `docs/spec/data-model.md:302-307`; `qpx/core/data/schemas/feature.yaml:49-62`; `docs/spec/versioning.md:208-209`

### CRITIC-63: Can the meaning of values change without a qpx_version change? {#critic-63}

*other* · ambiguous · partly confirmed (see correction) · checklist X1

- **Now:** Yes: converter fixes change values under the same specification version (pg_id changed under 1.1). writer_version and creator are free text, and the footer records no source-tool version.
- **Spec:** versioning.md:131-133,183-198; identity.md:67-69
- **Correction:** writer_version (serialization.md:287) and provenance.tool_version (provenance.md:32,52) do record the converter and source-tool versions. What holds is that semantic changes in converter output are not signalled by qpx_version.
- **Reader impact:** Two 1.3 files can encode MBR or q-values differently with no machine-readable signal.
- **Suggested:** Footer creator_version and converter_version; any semantic change bumps the spec minor version.
- **Evidence:** `docs/spec/versioning.md:131-133`; `docs/spec/versioning.md:183-198`

### CRITIC-64: Can values be compared or joined across datasets and producers? {#critic-64}

*cross* · specified, worth keeping · confirmed · checklist X1

- **Now:** Ids and natural keys must not be joined across datasets; analysis_id tells analyses apart. Dataset-scope values (global q, MaxLFQ, coverage, groups) change when data is converted in batches. pg.intensity is comparable only when quantification_method agrees.
- **Spec:** identity.md:5-12; intensities.md:198-209; pg.md:52-60
- **Reader impact:** Keep: these explicit non-comparability rules prevent naive pooling.
- **Suggested:** Keep, and extend the rule to feature intensities and q-value columns.
- **Evidence:** `docs/spec/identity.md:5-12`; `docs/spec/intensities.md:198-209`

### CRITIC-65: Do the reference/mapping pages agree with the normative view pages? {#critic-65}

*other* · contradictory · partly confirmed (see correction) · checklist X1

- **Now:** No. tool-mappings lists pg run_file_name and pg intensities (both removed in 1.1), a DIA-NN PSM mapping (there is no DIA-NN PSM view), and feature_counts = Precursor.Quantity. ontology lists diann_qvalue and diann_global_qvalue, which scores.md says no converter emits, and still maps feature 'intensity'.
- **Spec:** tool-mappings.md:17-30,105-109; converter-coverage.md:12; pg.md:114; ontology.md:211,288-289; scores.md:80-81
- **Correction:** diann_qvalue and diann_global_qvalue are at ontology.md:290-291. The ontology maps the feature column intensities (:209), not 'intensity'. Its actually stale entry is a pg 'intensities' column (:244) that pg lacks. The tool-mappings errors (:17-30 DIA-NN PSM column, :105 run_file_name, :107, :109) are confirmed.
- **Reader impact:** Producers implementing from the mapping pages write stale semantics.
- **Suggested:** Label pages normative or informative, and generate the mapping tables from code.
- **Evidence:** `docs/spec/tool-mappings.md:105-109`; `docs/spec/ontology.md:288-289`; `docs/spec/scores.md:80-81`

### CRITIC-M5: Is the typed q-value 'portable filter' usable for every producer? {#critic-m5}

*cross* · ambiguous · found by the verifier · checklist D4, P5, F5

- **Now:** psm.confident() needs psm.qvalue; feature.confident() uses COALESCE(precursor_qvalue, peptide_qvalue), mixing levels (qpx/core/data/psm.py:23-25; feature.py:52-54). Results: MaxQuant 0/1000 PSMs and 0/9,301 features; Mascot mzIdentML 0/13,386 PSMs; CDAP 90/100 PSMs.
- **Spec:** scores.md:38 calls the typed columns 'the portable way to filter'. scores.md:47 says a column stays null without a matching score. Nothing says what PEP-only or pass_threshold-only data should do.
- **Reader impact:** The documented confidence filter silently returns nothing for PEP-only producers, and mixes precursor- and peptide-level FDR for others.
- **Suggested:** Declare the confidence basis per file (qvalue level, PEP, or pass_threshold). Make confident() level-explicit and error when the column is all-null.
- **Evidence:** `docs/spec/scores.md:38-47`; `qpx/core/data/psm.py:23-25`; `qpx/core/data/feature.py:52-54`; `work/verify-critic/{mq,mzid,cdap} confident() counts`

### CRITIC-M6: Do the coverage and mapping tables match what 1.3 converters write? {#critic-m6}

*other* · contradictory · found by the verifier · checklist X1, P10, G3

- **Now:** No, in several places: MaxQuant fills psm.feature_id (95/1000 in the run); FragPipe sets pg.contaminant from the lead accession (fragpipe/pg_adapter.py:469); Spectronaut nulls id_run_file_name for imputed rows (spectronaut/feature_adapter.py:341); MaxQuant fills rt_start/rt_stop; CDAP fills peptide_qvalue (CRITIC-33).
- **Spec:** converter-coverage.md:285 ('Not implemented'), :303 ('Not implemented'), :203 ('Filled'); tool-mappings.md:50 (MaxQuant rt_start '---'); converter-coverage.md:181.
- **Reader impact:** Readers who rely on the documented null semantics per producer misread populated or null columns.
- **Suggested:** Generate the coverage and mapping tables from fixture conversions in CI, and mark them informative rather than normative.
- **Evidence:** `docs/spec/converter-coverage.md:203`; `docs/spec/converter-coverage.md:285`; `docs/spec/converter-coverage.md:303`; `docs/spec/tool-mappings.md:50`; `qpx/converters/fragpipe/pg_adapter.py:469`; `qpx/converters/spectronaut/feature_adapter.py:341`; `work/verify-critic/mq`
