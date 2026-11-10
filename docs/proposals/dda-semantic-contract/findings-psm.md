# Findings: PSM view

Part of the annex of the [DDA semantic contract proposal](index.md). Each entry is one verified finding from the audit of qpx `dev` at commit `32cbbdc` (specification 1.3). A skeptic re-checked every finding against the code; "partly confirmed" entries are true only as stated in their correction. `audit-run:` evidence refers to converting the repository's own fixtures (`tests/examples`) with that commit during the audit; those outputs are not committed.

105 findings, grouped by audit area.

## OpenMS converters

### OPENMS-01: Every hit of every PeptideIdentification becomes a PSM; rank never read or stored {#openms-01}

*psm* · unspecified · confirmed · checklist P1, D8

- **Now:** psm_adapter.py:260-273,375-381 iterate pid.getHits() and emit one row per distinct (peptidoform,charge); hit rank is never read (streaming.py:394-400 does not parse it either). No `rank` score is written. Test fx/conflict.consensusXML: one PID on scan 24300 with a rank-1 decoy (PEP 1e-4) and a rank-2 target (PEP 0.9) gave 2 PSM rows that cannot be told apart.
- **Spec:** psm.md:12-15: the composite names one peptidoform per spectrum. psm.md:291 recommends a `rank` score. Nothing defines a rank policy (rank-1 only or all ranks) or requires `rank` when several hits per spectrum are kept.
- **Reader impact:** Counting PSMs or identified spectra, or picking the best PSM per scan, silently includes rank-2..n alternatives, which look like real co-identifications.
- **Suggested:** The contract declares a rank policy, as a footer key such as `psm_ranks=top1|all`. When more than one hit per spectrum is kept, `rank` is required in additional_scores. The OpenMS converter reads PeptideHit rank and emits it.
- **Evidence:** `qpx/converters/openms_consensus/psm_adapter.py:260-273`; `qpx/converters/openms_consensus/psm_adapter.py:375-381`; `qpx/converters/openms_consensus/streaming.py:394-400`; `docs/spec/psm.md:12-15`; `docs/spec/psm.md:291`; `test: audit-run:openms/fx/conflict.consensusXML -> out_conflict/openms.psm.parquet rows scan [24300]`

### OPENMS-02: Chimeric / multiple PIDs per spectrum emitted as plain PSMs, no marker {#openms-02}

*psm* · unspecified · confirmed · checklist P2

- **Now:** Separate PeptideIdentifications that reference one spectrum (each with a different precursor MZ) become independent PSM rows with nothing marking them as co-isolated. Fixture MS1LabeledWorkflow_6: scan 11873 gives 3 PSMs (EYNVTIDAAESTVNK, TALGYYWALTEASGK and its heavy form). There is no PSM-level counterpart of feature identification_type 'secondary', and the converter never writes 'secondary'.
- **Spec:** Nothing for the PSM view. feature.md:145 defines `secondary` (a co-fragmenting peptide of a chimeric spectrum) only for features. scope.md:41 covers only the ProForma chimeric notation.
- **Reader impact:** Readers cannot tell genuine chimeric co-identifications from rank-n alternatives (OPENMS-01) or from duplicates. FDR or spectral-count analyses treat them alike.
- **Suggested:** Add a PSM-level flag or score that marks a chimeric or secondary identification (for example `identification_type` on psm, or `chimeric_rank`), separate from the search rank.
- **Evidence:** `tests/examples/openms_topp/MS1LabeledWorkflow_6_ratios_out.consensusXML:87,160,189`; `audit-run:openms/MS1LabeledWorkflow_6_ratios_out/openms.psm.parquet`; `docs/spec/feature.md:145`; `docs/spec/scope.md:41`

### OPENMS-03: No FDR filtering, significance_threshold ignored, pass_threshold never filled {#openms-03}

*psm* · diverges from spec · partly confirmed (see correction) · checklist P3, D3, D8

- **Now:** All PSMs are written: assigned, unassigned (the default; converter.py:143-147), decoys, and every rank. PeptideIdentification significance_threshold is never read (it appears nowhere in qpx/converters/openms_consensus). The `pass_threshold` column is present but null on every row (LFQ output 0/115). Provenance records no FDR step (provenance.py:154-180; LFQ fixture: search, alignment, grouping, conversion).
- **Spec:** producer-guide.md:11: keep rows above your FDR threshold only if you say so with `pass_threshold`. psm.md:66: null when the validator does not say. psm.yaml:133-136.
- **Correction:** Code facts hold: no pass_threshold, significance_threshold unread, unassigned PIDs kept by default. But this is a gap, not a divergence. psm.md:66 and psm.yaml:133-136 allow null when the tool does not say, and every fixture PID has significance_threshold="0". Provenance does record include_unassigned_psms=True, as provenance.md:65 asks.
- **Reader impact:** Readers cannot know whether the PSM view is FDR-filtered, at which level or threshold, or which rows pass. Unassigned and conflict-removed IDs are mixed with confident ones.
- **Suggested:** Require a footer or provenance declaration of the filtering state (none/psm/peptide/protein plus threshold). Set pass_threshold from significance_threshold when it is non-zero, or from a supplied threshold, and otherwise declare 'unfiltered'.
- **Evidence:** `qpx/converters/openms_consensus/converter.py:143-147`; `qpx/converters/openms_consensus/provenance.py:154-180`; `docs/spec/producer-guide.md:11`; `docs/spec/psm.md:66`; `qpx/core/data/schemas/psm.yaml:133-136`; `audit-run:openms/ProteomicsLFQ_1_out/openms.provenance.parquet`

### OPENMS-04: Peptide-level Percolator: psm.qvalue gets peptide q-value or 1.0 placeholder {#openms-04}

*psm* · contradictory · confirmed · checklist D4, P5

- **Now:** When score_type is q-value, the primary hit score is always written as 'q-value' (psm_adapter.py:215-219), even in peptide-level mode. Only the meta path is suppressed (lines 220-224). derived.py:15-21,29 then copies 'q-value' into psm.qvalue. PSM PEP is set null in this mode (line 279). Test fx/pepfdr: the non-best ELEFANTK PSM has psm.qvalue=1.0 (the placeholder), and the other PSMs carry peptide-level values as psm.qvalue.
- **Spec:** scores.md:42: psm.qvalue is PSM-level and excludes higher-level names. feature.md:130: the OpenMS converters 'never the 1.0 placeholders'. psm.md:286: PEP is a PSM metric.
- **Reader impact:** psm.confident() and q-value filters use peptide-level q-values as PSM q-values. Placeholder 1.0 rows look like confident failures, and PSM PEP disappears in this mode.
- **Suggested:** In peptide-level mode, do not write the primary score as a PSM q-value. Name it `peptide_qvalue` and drop the 1.0 placeholders. The contract states that psm.qvalue must be PSM-level or null.
- **Evidence:** `qpx/converters/openms_consensus/psm_adapter.py:215-224,279`; `qpx/writers/derived.py:15-21,29`; `docs/spec/scores.md:42`; `docs/spec/feature.md:130`; `test: audit-run:openms/fx/pepfdr.consensusXML -> out_pepfdr/openms.psm.parquet (row scan [24218] qvalue=1.0)`

### OPENMS-16: psm.observed_mz written as 0.0 when the precursor m/z is absent {#openms-16}

*psm* · diverges from spec · confirmed · checklist P7

- **Now:** obs_mz = float(pid.getMZ()) if truthy, else 0.0 (psm_adapter.py:359). The streaming reader also defaults a missing MZ attribute to 0.0 (streaming.py:404). mass_error_ppm is then null (utils.py:117-128), but observed_mz stays 0.0. rt is the PID RT; whether the precursor m/z is instrument-reported or corrected is not recorded.
- **Spec:** psm.yaml:51-55 and psm.md:47: observed_mz is null when no measured precursor m/z is reported. versioning.md:67 says psm.observed_mz is nullable since 1.3.
- **Reader impact:** Readers filtering or plotting precursor m/z see 0.0 values as real measurements.
- **Suggested:** Write null. The contract states whether observed_mz is the instrument-reported or the corrected monoisotopic precursor.
- **Evidence:** `qpx/converters/openms_consensus/psm_adapter.py:359`; `qpx/converters/openms_consensus/streaming.py:404`; `qpx/core/data/schemas/psm.yaml:51-55`; `docs/spec/psm.md:47`; `docs/spec/versioning.md:65-67`

### OPENMS-17: Undeclared scan semantics: surrogate hashes and spectrum=/index refs read as scans {#openms-17}

*psm* · diverges from spec · partly confirmed (see correction) · checklist P7, D8

- **Now:** Native IDs with no numeric part get a blake2s int32 surrogate in scan (psm_adapter.py:49-69). scan_format is declared only when every reference has one recognised format (psm_adapter.py:385-408; scan.py:17-31). The LFQ fixture (`spectrum=N`) gets no scan_format footer, so its values are read as scan numbers. The feature view never declares scan_format (conversion.py:386-392).
- **Spec:** scan.md:80-81: a value that is not a scan number must be declared, and a file without scan_format is read as scan numbers. producer-guide.md:23-24: set scan_format. scan.md:64-70,89-92.
- **Correction:** The undeclared blake2s surrogate is confirmed. The spec contradicts itself: scan.md:64-67 lets unknown formats stay undeclared, while scan.md:80-81 and versioning.md:68-69 require a declaration. Leaving feature scan_format undeclared is documented design (scan.md:90-92). It is unknown, not shown wrong, whether `spectrum=N` is a scan number.
- **Reader impact:** Hash surrogates and spectrum indices are joined to mzML scans as if they were scan numbers, giving wrong spectra silently.
- **Suggested:** Undeclarable values must not go into scan. Leave it empty and store the native ID string (#259 spectrum_reference), or declare a format such as `surrogate`. The feature footer must state its format.
- **Evidence:** `qpx/converters/openms_consensus/psm_adapter.py:45-69,385-408`; `qpx/core/scan.py:17-44`; `qpx/converters/openms_consensus/conversion.py:386-401`; `docs/spec/scan.md:64-70,80-81,89-92`; `audit-run:openms/ProteomicsLFQ_1_out/openms.psm.parquet (no scan_format)`

### OPENMS-18: Score names not snake_case; engine raw scores dropped; collisions renamed '_2' {#openms-18}

*psm* · diverges from spec · confirmed · checklist P5

- **Now:** The primary score is written as 'q-value' or the raw score_type (e.g. 'COMET:xcorr', else 'search_score') (psm_adapter.py:218). Colliding hits are written as 'q-value_2' and so on, without naming their engine (86-99,227-233). Hit meta scores (e.g. OMSSA_score, present on all 115 fixture hits) are not exported. `qpxc validate` on the LFQ output warns that 115 entries are not snake_case ('q-value').
- **Spec:** scores.md:23-24: all score names are snake_case. scores.md:78: {tool}_{score}. psm.md:291 recommends global_qvalue and rank.
- **Reader impact:** Consumers matching 'q_value' or 'comet_xcorr' miss values. '_2' scores have no engine attribution, and the raw engine evidence is lost.
- **Suggested:** Normalise names at write time, name collision scores by engine, and export known engine meta scores.
- **Evidence:** `qpx/converters/openms_consensus/psm_adapter.py:86-99,204-245`; `tests/examples/openms_topp/ProteomicsLFQ_1_out.consensusXML:1084`; `docs/spec/scores.md:23-24,78`; `cmd: qpxc validate --dataset-path audit-run:openms/ProteomicsLFQ_1_out --semantic`

### OPENMS-19: Cross-PID duplicates of one spectrum/peptidoform silently dropped with their scores {#openms-19}

*psm* · unspecified · confirmed · checklist P6, P5

- **Now:** A shared seen-set on (peptidoform,charge,run,scan) keeps the first PID's row. Later PIDs with the same key are used only to re-rank the feature link, and their scores and engines are discarded (psm_adapter.py:375-381; record_linking.py:81-135). Inside one PID, colliding hits collapse to the lowest PEP with the others' scores appended (psm_adapter.py:293-296). The same scan with a different charge gives separate rows.
- **Spec:** psm.md:12-15: several hits of one peptidoform/charge per spectrum need a producer psm_id or an extended composite. #262 decided to keep the lowest PEP within one identification. Nothing covers PIDs from different engines or copies.
- **Reader impact:** Which engine's evidence survives depends on document order, and a reader cannot tell that a merge happened.
- **Suggested:** The contract says that one row per (spectrum, peptidoform, charge) merges all engines. Every engine's scores are kept with engine-prefixed names, and the PEP rule is the lowest PEP.
- **Evidence:** `qpx/converters/openms_consensus/psm_adapter.py:293-296,375-381`; `qpx/converters/openms_consensus/record_linking.py:81-135`; `docs/spec/psm.md:12-15`

### OPENMS-20: psm.feature_id links PSMs of other peptidoforms (even decoys) to a feature {#openms-20}

*psm* · ambiguous · confirmed · checklist P10, F2

- **Now:** Every PSM from a consensus feature's PIDs is linked to that feature's row in the PSM's run, whatever its peptidoform (record_linking.py:66-70). Same-peptidoform only ranks between copies (21-33). Test fx/conflict: the decoy IMTELLERK PSM and the rank-2 ELEFANTK PSM are both linked to the ELEFANTK feature.
- **Spec:** psm.md:75: feature_id is the producer's assignment of the PSM to a consensus feature. psm.md:158-161 ranks only among copies. intensities.md:93 allows different peptidoforms only for MS1-label partners. Nothing defines feature_id for conflicting peptides.
- **Reader impact:** Readers using psm.with_feature() as 'PSMs supporting this feature' include contradicting and decoy identifications.
- **Suggested:** feature_id links only PSMs whose peptidoform equals the feature's peptidoform (MS1 partners excepted). Conflicting ones stay null or get an explicit relation.
- **Evidence:** `qpx/converters/openms_consensus/record_linking.py:21-33,66-78`; `docs/spec/psm.md:75,147-175`; `docs/spec/intensities.md:93`; `test: audit-run:openms/fx/out_conflict/openms.psm.parquet`

### OPENMS-28: Unknown mass shifts written as absolute residue mass, 15 decimals, unsigned {#openms-28}

*psm* · diverges from spec · partly confirmed (see correction) · checklist D5, P9

- **Now:** to_proforma rewrites OpenMS toUniModString without normalisation (feature_adapter.py:240-254), and the writer does not normalise either. Verified with pyopenms 3.6: K[+54.01063456] gives 'K[182.105597942799989]', and N[+2715.9624] gives 'N[2830.005328191399713]'. Labels and fixed mods missing from the OpenMS sequence are not added: the IsobaricWorkflow fixture writes 'ELEFANTK' without TMT.
- **Spec:** peptidoform.md:17-19: fixed mods and isobaric labels inline. peptidoform.md:27-30: a signed delta mass with 4 decimals, `[+54.0106]`. producer-guide.md:52.
- **Correction:** The mass-shift bug is confirmed, including negative shifts: K[-17.0265] becomes K[111.0684...]. The suggested fix is not enough: normalize_peptidoform gives K[+182.1056], still a residue mass, so the delta must come from the modification itself. The TMT evidence comes from a synthetic fixture with unlabelled sequences; real OpenMS search output carries fixed and label mods in the AASequence, so 'labels not added' is unproven.
- **Reader impact:** Readers parse a +182 Da shift where the real shift is +54. The same molecule gets different peptidoform strings and ids across producers.
- **Suggested:** Run normalize_peptidoform (delta, sign, 4 dp) in every converter or in the writer. Check for missing labels and fixed mods against SearchParameters.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:240-254`; `qpx/converters/ptm.py (normalize_peptidoform; not called by openms converters)`; `docs/spec/peptidoform.md:17-19,27-30`; `docs/spec/producer-guide.md:52`; `cmd: python to_proforma(AASequence.fromString('PEPTIDEK[+54.01063456]R'))`

### OPENMS-M1: Primary score types mishandled: PEP-as-score lost, 'fdr' relabelled q-value {#openms-m1}

*psm* · diverges from spec · found by the verifier · checklist P5, D4

- **Now:** pep_of reads PEP only from meta keys (feature_adapter.py:31,109-114). With score_type 'Posterior Error Probability' (IDPEP output) the PEP lands in additional_scores as 'Posterior Error Probability', and psm and feature PEP are null. score_type 'fdr' is written as 'q-value' (psm_adapter.py:218,363). A Percolator meta q-value fills psm.qvalue (220-224), but feature.peptide_qvalue stays null (feature_adapter.py:158-159).
- **Spec:** psm.md:44 and scores.md:55: PEP is the top-level field, not stored in additional_scores. scores.md:23-24: names are snake_case. scores.md:42-43: typed q-values have a level.
- **Reader impact:** PEP-based filtering finds no PEP for OpenMS files scored by PEP. A raw FDR is read as a q-value. The same source gives a PSM q-value but no feature q-value.
- **Suggested:** Map the primary score by score_type: PEP goes to posterior_error_probability, a q-value goes to a level-named q-value, FDR is kept as fdr. Use the same mapping for psm and feature.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:31,109-114,149-160`; `qpx/converters/openms_consensus/psm_adapter.py:215-224,362-363`; test: audit-run:verify-openms/pepscore.consensusXML (IsobaricWorkflow, score_type=PEP, PEP meta removed) -> out_pepscore: PEP null on 4/4 psm and 4/4 feature; additional_scores 'Posterior Error Probability'; `docs/spec/psm.md:44`; `docs/spec/scores.md:23-24,55`

### OPENMS-M4: PSMs silently dropped when run unresolved or spectrum reference empty {#openms-m4}

*psm* · unspecified · found by the verifier · checklist P3, X1

- **Now:** psm_records_for_pid returns [] with no log when the run resolver gives None, for example an unassigned PID in a multi-run map with no map_index or id_merge_index (psm_adapter.py:348-350, resolver 184-199). It returns [] with only a debug log for an empty spectrum_reference (354-358). pg still counts these hits (pg_adapter.py:104-123). The fixtures had no drops (115/115).
- **Spec:** Nothing on whether the psm view must hold every source identification, or how omitted ones are reported.
- **Reader impact:** PSM counts and pg evidence can disagree, and a reader cannot tell that identifications were omitted.
- **Suggested:** The contract requires a count of omitted identifications (provenance parameter or warning), and states whether the psm view is complete.
- **Evidence:** `qpx/converters/openms_consensus/psm_adapter.py:184-199,348-358`; `qpx/converters/openms_consensus/pg_adapter.py:104-123,140-151`

## MaxQuant converter

### MAXQUANT-03: Chimeric second-peptide PSMs can duplicate psm_id; marked only in free-text cv_params {#maxquant-03}

*psm* · diverges from spec · confirmed · checklist P2, P6

- **Now:** MULTI-SECPEP rows of msms.txt become ordinary PSMs. Type survives only as cv_params {type: MULTI-SECPEP} (psm_adapter.py:348-353). The first and second peptide of one scan can share peptidoform and charge, and then both hash to the same psm_id. Reproduced: audit-run:maxquant/dup wrote 2 rows with 1 psm_id, and validate reported a duplicate PK. Real data: evidence ids 167037/167038 (full fixture) share scan 10817. #378 reports 9 duplicates on PXD003133.
- **Spec:** psm.md:9-15: the composite names one peptidoform per spectrum, and more hits need a producer psm_id or an extended composite. identity.md:29.
- **Reader impact:** The file is INVALID and joins on psm_id are ambiguous. A second (chimeric) hit cannot be told from a first hit without parsing cv_params.
- **Suggested:** Add a typed PSM field (hit_index/rank or identification_type secondary) and include it, or MaxQuant's msms 'id', in the PSM composite.
- **Evidence:** `qpx/converters/maxquant/psm_adapter.py:189-252,348-353`; `audit-run:maxquant/dup/out/maxquant.psm.parquet`; `qpxc validate: 'Primary key (psm_id) has 1 duplicate row(s)'`; `tests/examples/maxquant/maxquant_full evidence ids 167037/167038 scan 10817`; `github #378`

### MAXQUANT-04: No rank or hit index; top-hit policy and dropped alternatives undeclared {#maxquant-04}

*psm* · unspecified · partly confirmed (see correction) · checklist P1, P2

- **Now:** msms.txt holds the top hit per MS/MS plus any second peptide. Simple fixture: 998 scans have 1 row and 1 scan has 2 (MULTI-MSMS+MULTI-SECPEP). The converter writes no rank and drops 'All sequences', 'All scores', 'Score diff' and 'Localization prob' (psm_adapter.py:314-354). It skips rows with an empty Modified sequence at debug log level only (psm_adapter.py:266-270).
- **Spec:** psm.md:291 recommends a `rank` score and scores.md:58 defines it. Nothing says whether a psm view is rank-1 only. Per #229, rank joins the PK if ranks are kept, but this is not in the spec.
- **Correction:** No rank and no declared top-hit policy: true. The skip at psm_adapter.py:266-270 fires only for '' or '__'. A NULL cell becomes str(None)='None' (psm_adapter.py:259) and is written as peptidoform 'None' (see MAXQUANT-M1).
- **Reader impact:** A reader cannot tell 'one row per spectrum' from 'several co-identified peptides' or 'alternative ranks' without producer knowledge.
- **Suggested:** A footer key declares the rank policy (e.g. top1_plus_chimeric). Write rank/hit_index on every PSM.
- **Evidence:** `qpx/converters/maxquant/psm_adapter.py:314-354`; `qpx/converters/maxquant/psm_adapter.py:266-270`; `docs/spec/psm.md:291`; `docs/spec/scores.md:58`

### MAXQUANT-05: FDR-filtered MaxQuant PSMs have null pass_threshold and qvalue; confident() drops all {#maxquant-05}

*psm* · ambiguous · confirmed · checklist P3, D3, D4

- **Now:** msms.txt holds only the rows MaxQuant kept, including passing Reverse hits (3/1000 in the simple fixture). The converter leaves pass_threshold null (0/1000 filled) and psm.qvalue null (0/1000), and records no filter state. psm.confident() requires qvalue<=x (qpx/core/data/psm.py:23-25), so it returns 0 MaxQuant PSMs.
- **Spec:** producer-guide.md:11: keep rows above your threshold only if you say so with pass_threshold. psm.md:66: null means the tool does not say. scores.md:47: the column stays null without a score.
- **Reader impact:** Pre-filtered MaxQuant PSMs look like unfiltered output. The portable confidence filter discards every MaxQuant PSM.
- **Suggested:** Set pass_threshold=true on MaxQuant msms rows. Declare the filter level and threshold in the footer or provenance. confident() falls back to pass_threshold.
- **Evidence:** `qpx/converters/maxquant/psm_adapter.py:230-252`; `qpx/core/data/psm.py:23-25`; `docs/spec/producer-guide.md:11`; `docs/spec/psm.md:66`; `audit-run:maxquant/simple/maxquant.psm.parquet`

### MAXQUANT-23: psm.observed_mz is derived from the peptidoform, not measured {#maxquant-23}

*psm* · ambiguous · partly confirmed (see correction) · checklist P7

- **Now:** psm.observed_mz = calculated_mz + m/z·ppm·1e-6 (psm_adapter.py:290-295; constants.py:259-271), so peptidoform errors carry over. Simple fixture: AAAIGIDLGTTYSCVGVFQHGK z3 has msms m/z 755.716, but observed_mz is 889.483 with mass_error -0.89 ppm. The feature view uses m/z·(1+ppm) instead (feature_adapter.py:569-571), where the same errors show as about 32,000 ppm.
- **Spec:** psm.md:47: 'Experimental observed m/z'. tool-mappings.md:26 documents the derivation.
- **Correction:** The derivation is documented. The 889.48 vs 755.72 example is a fixture artifact: the TMT SDRF adds labels to LFQ msms. What remains: PSM mass_error_ppm reproduces MaxQuant's ppm by construction, so peptidoform errors stay invisible, and the PSM and feature observed_mz formulas differ.
- **Reader impact:** observed_mz is not a measurement, mass-error QC cannot catch wrong peptidoforms, and PSM and feature observed_mz are not comparable.
- **Suggested:** The contract defines observed_mz as the measured precursor m/z, or null. Derived values go elsewhere.
- **Evidence:** `qpx/converters/maxquant/psm_adapter.py:288-312`; `qpx/converters/maxquant/constants.py:259-271`; `qpx/converters/maxquant/feature_adapter.py:568-571`; `docs/spec/psm.md:47`; `docs/spec/tool-mappings.md:26,47`

### MAXQUANT-25: psm.feature_id is linked by Evidence ID without a consistency check; docs say unimplemented {#maxquant-25}

*psm* · diverges from spec · partly confirmed (see correction) · checklist P10

- **Now:** psm.feature_id is looked up from 'Evidence ID' (psm_adapter.py:229-232; feature_adapter.py:210-221) with no check of run or peptidoform. Collapsed evidence ids map to the kept row. The simple fixture's msms and evidence come from different searches: 95 PSMs were linked to MFA381 features, and 0/95 peptidoforms match.
- **Spec:** psm.md:75: feature_id is the authoritative foreign key. converter-coverage.md:285: 'Not implemented: the Evidence ID of msms.txt is not read' (stale).
- **Correction:** The doc is stale (true). The link itself is MaxQuant's own Evidence ID, the producer's assignment that psm.md:75 asks for. The 95 mismatches come from the fixture, whose msms and evidence are different searches. The missing check affects robustness only and is not a spec divergence.
- **Reader impact:** Mismatched inputs give silently wrong PSM-feature links, and the docs mislead.
- **Suggested:** Check that the link has the same run and sequence (label-aware), and update converter-coverage.
- **Evidence:** `qpx/converters/maxquant/psm_adapter.py:229-232`; `qpx/converters/maxquant/feature_adapter.py:210-221,284-289`; `docs/spec/converter-coverage.md:285`; `audit-run:maxquant/simple (join psm-feature)`

## FragPipe converter

### FRAGPIPE-01: Pre-filtered psm.tsv written without pass_threshold or any record of the FDR cut {#fragpipe-01}

*psm* · diverges from spec · partly confirmed (see correction) · checklist P3, D3, D8

- **Now:** Every psm.tsv row is written (psm_adapter.py:106-112). Philosopher has already filtered psm.tsv, usually at 1% PSM/peptide/protein FDR, but pass_threshold is null on every row (verified on a synthetic run). Provenance has only 'database_search FragPipe/MSFragger' and 'format_conversion' (converter.py:95-120): no fdr_filtering, rescoring or protein_inference step, tool_version null. No footer key either.
- **Spec:** producer-guide.md:11: keep rows above your FDR threshold only if you say so with pass_threshold. psm.md:66: pass_threshold null = validator does not say. provenance.md:3,82: the fdr_filtering category should answer 'what FDR was used?'. Nothing defines a footer key for filtering state.
- **Correction:** The code facts are correct. The format_conversion step also records the input files and their SHA-256 (postprocess.py:200-224). It is not a divergence, though: producer-guide.md:11 only requires pass_threshold when sub-threshold rows are KEPT, and FragPipe keeps none. No spec rule makes a producer declare a pre-filter, and provenance.md:3 is a statement of purpose. The status should be 'unspecified'.
- **Reader impact:** A reader cannot tell this filtered view from an unfiltered one (MaxQuant msms, mzIdentML), or at which level and threshold it was filtered. Null pass_threshold reads as 'unknown', although every row passed.
- **Suggested:** The contract should require a declared filtering state, e.g. a footer psm_filter=none|psm|peptide|protein plus threshold. FragPipe should set pass_threshold=true and add an fdr_filtering provenance step with the Philosopher filter parameters, read from fragpipe.workflow when it is available.
- **Evidence:** `qpx/converters/fragpipe/psm_adapter.py:106-112`; `qpx/converters/fragpipe/converter.py:95-120`; `docs/spec/producer-guide.md:11`; `docs/spec/psm.md:66`; `docs/spec/provenance.md:3`; `docs/spec/provenance.md:82`; `run: qpxc convert fragpipe on synthetic inputs -> audit-run:fragpipe/out/fp.psm.parquet pass_threshold all null; footer has no filter key`

### FRAGPIPE-02: No q-value anywhere; PEP is 1-Probability whatever the validator; confident() returns nothing {#fragpipe-02}

*psm* · ambiguous · confirmed · checklist P5, D4

- **Now:** PEP = 1 - 'PeptideProphet Probability' (else 'Probability'), in both the psm and feature adapters (psm_adapter.py:245-249, feature_adapter.py:401-403,471-474). FragPipe 20+ can produce that probability from Percolator rather than PeptideProphet, and this is not recorded. No q-value is read, so psm.qvalue and the feature q-value columns are always null. Ontology labels the PEP 'percolator:PEP MS:1001493' (core/scores.py:172-177).
- **Spec:** psm.md:286 and tool-mappings.md:23,33: PEP = 1 - PeptideProphet probability. scores.md:36-47: typed q-value columns. psm.confident() is 'is_decoy=false AND qvalue<=x' (core/data/psm.py:23-25), and feature.confident() is the same with COALESCE(precursor_qvalue, peptide_qvalue) (core/data/feature.py:52-54).
- **Reader impact:** psm.confident() and feature.confident() return zero FragPipe rows, although all of them passed 1% FDR. The PEP's model (PeptideProphet or Percolator) is unknown, and the ontology attribution is wrong in the PeptideProphet case.
- **Suggested:** The contract should say what a reader does when no q-value exists but the view is pre-filtered, for example confident() falling back to pass_threshold. Record the validator (PeptideProphet or Percolator) in provenance, and map the PEP ontology per converter.
- **Evidence:** `qpx/converters/fragpipe/psm_adapter.py:245-249`; `qpx/converters/fragpipe/feature_adapter.py:401-403`; `qpx/core/scores.py:172-177`; `qpx/core/data/psm.py:23-25`; `qpx/core/data/feature.py:52-54`; `qpx/writers/derived.py:15-34`; `docs/spec/psm.md:286`; `docs/spec/scores.md:36-47`; `run: audit-run:fragpipe/out/fp.ontology.parquet maps posterior_error_probability -> percolator:PEP`

### FRAGPIPE-03: Several hits per spectrum (DDA+/chimeric) written with no rank or chimeric marker {#fragpipe-03}

*psm* · unspecified · confirmed · checklist P1, P2, P6

- **Now:** Each psm.tsv row becomes one PSM. No rank or chimeric column is read, and none is derived (psm_adapter.py:144-199). Verified: scan 100 gave three rows (PEPTIDEK z2, ANTHERSK z2, PEPTIDEKR z3) with the same scan and nothing to tell them apart. The psm_id stays unique only because peptidoform or charge differ. Whether FragPipe 22 DDA+ psm.tsv carries a rank column was not verified.
- **Spec:** psm.md:12-15 covers only several hits of ONE peptidoform. scores.md:58 recommends a 'rank' score. Nothing distinguishes co-isolated (chimeric) hits from rank-2 alternatives, and #229's 'rank joins the PK' is not in the spec.
- **Reader impact:** Counting PSMs per spectrum, or assuming one peptide per scan, double-counts DDA+ spectra. A reader cannot tell whether a second peptide on a scan is a co-isolated precursor or a weaker alternative.
- **Suggested:** The contract should define a psm 'rank' (null = only the top hit kept) and a chimeric/co-isolation flag (or identification role). FragPipe writes rank 1 for all, plus a flag when several precursors share a scan in DDA+ mode.
- **Evidence:** `qpx/converters/fragpipe/psm_adapter.py:144-199`; `docs/spec/psm.md:12-15`; `docs/spec/scores.md:58`; `run: audit-run:fragpipe/out/fp.psm.parquet rows 0-2 share scan [100]`

### FRAGPIPE-04: Run names containing dots truncate run_file_name, empty scan and duplicate psm_id {#fragpipe-04}

*psm* · diverges from spec · partly confirmed (see correction) · checklist P6, P7, X1

- **Now:** The psm view takes the first dot-token of Spectrum as run_file_name and the second as scan (psm_adapter.py:45-58). 'run.D.v2.00100.00100.2' becomes run 'run' with scan []. Two such PSMs of one peptidoform share a psm_id, which the writer only logs as a warning (writers/base.py:559-569). The feature lookup (feature_adapter.py:431-442) and plex_runs (isobaric.py:249, which strips the last three tokens) parse the same id differently. Verified.
- **Spec:** psm.md:53: run_file_name is the spectrum file name without path or extension. psm.md:54: scan identifier. identity.md:34: the psm composite is run_file_name, scan, peptidoform, charge.
- **Correction:** The feature PSM lookup parses Spectrum exactly like the psm view: first dot-token, then normalize_run_name (feature_adapter.py:434-435). Only plex_runs differs: it strips the last three numeric tokens (isobaric.py:38). So the mismatch is between psm/lookup ('run') and the experiment_annotation names in feature and pg ('run.D.v2'). The rest holds: empty scan, duplicate id kept only as a warning, and the joins fail.
- **Reader impact:** PSMs from dotted raw files land on a non-existent run. They have no scan and collide on psm_id. In the same output the feature, pg and run views name the full run, so joins across views fail silently.
- **Suggested:** Parse Spectrum from the right: run = everything before the last three dot-tokens. Use one helper in all FragPipe adapters, and make a duplicate psm_id an error in conformance.
- **Evidence:** `qpx/converters/fragpipe/psm_adapter.py:45-58`; `qpx/converters/fragpipe/feature_adapter.py:431-442`; `qpx/converters/fragpipe/isobaric.py:246-252`; `qpx/writers/base.py:559-569`; `run: audit-run:fragpipe/psm_dot.parquet -> 2 rows run='run', scan=[], identical psm_id -6525761943204802137`

### FRAGPIPE-05: No scan_format declared in the psm or feature footer {#fragpipe-05}

*psm* · diverges from spec · partly confirmed (see correction) · checklist P7, D8

- **Now:** PsmWriter and FeatureWriter are opened without scan_format (psm_adapter.py:106, feature_adapter.py:178-182), and set_scan_format is never called. Scan is the MSFragger Spectrum token: Thermo scan number, or whatever MSFragger numbers for timsTOF/.d (not verified). Verified footers: identity_composite, identity_source=derived, no scan_format.
- **Spec:** producer-guide.md:23: a DDA PSM or feature holds the MS2 scan number(s); set scan_format, or the validator cannot check them. scan.md:51,80: a value that is not a scan number must be declared.
- **Correction:** scan.md:80-81 says a file without scan_format is read as holding scan numbers, so FragPipe output for Thermo data is compliant. MaxQuant, PD and CDAP also omit the key, so this is not FragPipe-specific. A real divergence exists only if MSFragger timsTOF ids are not scan numbers, which is unverified. producer-guide.md:23 ('set it') vs scan.md:51 ('can') is an inconsistency within the spec.
- **Reader impact:** A reader cannot tell a scan number from an index, especially for timsTOF FragPipe data, so a join to mzML/mz views may be wrong.
- **Suggested:** Declare scan_format=scan for Thermo-style ids. Define in the contract how MSFragger timsTOF ids are declared.
- **Evidence:** `qpx/converters/fragpipe/psm_adapter.py:106`; `qpx/converters/fragpipe/feature_adapter.py:178-182`; `docs/spec/producer-guide.md:23`; `docs/spec/scan.md:51`; `run: pq.ParquetFile(out/fp.psm.parquet).metadata -> no scan_format`

### FRAGPIPE-06: observed_mz is MSFragger-calibrated m/z; calculated_mz recomputed including open-search shift {#fragpipe-06}

*psm* · ambiguous · confirmed · checklist P7

- **Now:** observed_mz prefers 'Calibrated Observed M/Z' over 'Observed M/Z' (psm_adapter.py:239, feature_adapter.py:393). calculated_mz is recomputed from the peptidoform, including an unassigned open-search shift (psm_adapter.py:165). So mass_error_ppm is the error after calibration, or meaningless for shifted PSMs; the latter is documented in converter-coverage.md:230-236. rt = 'Retention' in seconds (psm_adapter.py:243). Charge is the search-assigned one.
- **Spec:** psm.md:47: 'Experimental observed m/z', with nothing on calibration. tool-mappings.md:26 says FragPipe observed_mz = 'Observed M/Z', and tool-mappings.md:25 says calculated_mz = 'Calculated M/Z'.
- **Reader impact:** mass_error_ppm from FragPipe (calibrated) is not comparable with producers that report the raw instrument m/z. The tool-mappings page tells readers the uncalibrated column is used.
- **Suggested:** The contract should state whether observed_mz is raw or recalibrated. Keep the other value as an additional field or cv_param, and fix tool-mappings.
- **Evidence:** `qpx/converters/fragpipe/psm_adapter.py:165`; `qpx/converters/fragpipe/psm_adapter.py:239`; `qpx/converters/fragpipe/feature_adapter.py:393`; `docs/spec/psm.md:47`; `docs/spec/tool-mappings.md:25-26`; `docs/spec/converter-coverage.md:230-236`

### FRAGPIPE-08: Decoys taken only from the Protein prefix; psm.tsv normally holds no decoys {#fragpipe-08}

*psm* · unspecified · confirmed · checklist P4, F8

- **Now:** is_decoy is true when the razor 'Protein' starts with rev_, REV_, DECOY_ or decoy_ (constants.py:28-39, psm_adapter.py:158-160). A custom decoy tag is missed. A Philosopher psm.tsv normally excludes decoys, so is_decoy is false in practice and nothing tells a reader that decoys were removed. Feature is_decoy comes from the best PSM's Protein, else the ion's prefix (feature_adapter.py:482,543-545). A test pins a target ion flagged decoy by a rev_ PSM (test_fragpipe_converter.py:296-333).
- **Spec:** psm.md:45 defines is_decoy, with target+decoy counted as target. Nothing says whether a psm or feature view must keep decoys or declare that they were removed (#229 left decoys in features open).
- **Reader impact:** An all-false is_decoy is indistinguishable from 'decoys removed', so FDR cannot be recomputed. A feature can be is_decoy=true while its anchor_protein is a target.
- **Suggested:** Footer or provenance should declare decoys_included true/false. Feature is_decoy should come from the feature's own protein mapping, consistent with anchor_protein.
- **Evidence:** `qpx/converters/fragpipe/constants.py:28-39`; `qpx/converters/fragpipe/psm_adapter.py:158-160`; `qpx/converters/fragpipe/feature_adapter.py:543-545`; `tests/converters/test_fragpipe_converter.py:296-333`; `docs/spec/psm.md:45`

### FRAGPIPE-09: psm.feature_id never filled, although the feature adapter selects a PSM per feature {#fragpipe-09}

*psm* · unspecified · confirmed · checklist P10

- **Now:** The psm adapter writes no feature_id (psm_adapter.py:178-199), so every value is null (verified, 11/11). The feature adapter picks a best PSM per feature (feature_adapter.py:444-486) but records the link only as scan plus id_run_file_name. Issue #378 lists 'fill psm.feature_id' as a FragPipe todo.
- **Spec:** psm.md:75: feature_id is the authoritative link (optional); null = unlinked. feature.md:162: psm_ids is a computed inverse of it.
- **Reader impact:** FragPipe has no PSM-feature link. Joining on (run, peptidoform, charge) fails for multi-run experiments, MBR rows and FAIMS. Dataset.link_feature_psm() returns nothing.
- **Suggested:** The contract should say which PSMs a feature links (all PSMs of the precursor in the same run, or only the identifying one). FragPipe derives feature ids first and fills psm.feature_id.
- **Evidence:** `qpx/converters/fragpipe/psm_adapter.py:178-199`; `qpx/converters/fragpipe/feature_adapter.py:444-486`; `docs/spec/psm.md:75`; `github #378 (FragPipe todo list)`; `run: out/fp.psm.parquet feature_id all null`

### FRAGPIPE-26: Peptidoform canonical: Assigned Modifications inline, UNIMOD upper case, 4-decimal shifts {#fragpipe-26}

*psm* · specified, worth keeping · confirmed · checklist P9, D5

- **Now:** Assigned Modifications, which include fixed mods and TMT labels, are resolved to UNIMOD accessions; an unknown mass is written as a signed 4-decimal shift. Verified examples: '[UNIMOD:2016]-PEPTIDEK[UNIMOD:2016]', 'PES[+123.4567]TIDEK', 'PEPTIDEK-[UNIMOD:2]'. An open-search shift is placed with '?' and modification_candidates; −1 to 3 13C isotope errors within 0.1 Da or 20 ppm are dropped (constants.py:112-173). A feature never takes a mass-shifted PSM (feature_adapter.py:464-466).
- **Spec:** peptidoform.md:13-35: canonical form (fixed mods and labels inline, UNIMOD:n, 4-decimal masses, unknown position written '[x]?'). converter-coverage.md:230-241 documents the isotope window and that real shifts within it cannot be told apart from isotope errors.
- **Reader impact:** PSM and feature peptidoforms join on identical strings, and fixed mods are explicit. Only the isotope-window caveat (for example an unlisted deamidation) remains.
- **Suggested:** Keep the rule. The contract should name the isotope-error window as producer policy and record it, e.g. as a cv_param or provenance parameter.
- **Evidence:** `qpx/converters/fragpipe/constants.py:65-173`; `qpx/converters/fragpipe/feature_adapter.py:464-466`; `docs/spec/peptidoform.md:13-35`; `docs/spec/converter-coverage.md:230-241`; `run: to_proforma checks in work/fragpipe (see transcript)`

### FRAGPIPE-28: Spec rule for psm.qvalue names omits the glycan, lib, translated and site exclusions {#fragpipe-28}

*psm* · diverges from spec · confirmed · checklist D4, P5

- **Now:** PTM-Shepherd 'Glycan q-value' is stored as additional score 'glycan_qvalue' (psm_adapter.py:295-301). derived.py:16 excludes the glycan, lib, translated and site prefixes, so psm.qvalue stays null. Hyperscore and Expectation are stored as msfragger_hyperscore and msfragger_expectation, and Purity as parent_ion_fraction (psm_adapter.py:263-292).
- **Spec:** scores.md:42 says psm.qvalue is filled from '<engine>_qvalue' except names of a higher level (peptide, peptidoform, protein, pg, gg, global, precursor). 'glycan_qvalue' literally matches '<engine>_qvalue' and is not in that exclusion list.
- **Reader impact:** A non-qpx writer following scores.md would put a glycan-assignment q-value into psm.qvalue, so psm.confident() would filter PSMs on glycan FDR.
- **Suggested:** List the code's full exclusion set in scores.md, or use an allow-list of explicit PSM-level q-value names.
- **Evidence:** `qpx/converters/fragpipe/psm_adapter.py:263-301`; `qpx/writers/derived.py:15-21`; `docs/spec/scores.md:42`

## Proteome Discoverer converter

### PD-01: PSM view keeps rank 2..n hits by default; rank exists only as nested score {#pd-01}

*psm* · unspecified · partly confirmed (see correction) · checklist P1, D8

- **Now:** psm_adapter.py:37-40 writes every rank unless --rank1-only (cli/convert/proteome_discoverer.py:46). Rank is only additional_scores 'rank' (reports.py:139). Run on synthetic export: a rank-2 'Low' hit was written next to rank 1 on scan 100 with pass_threshold=false. The policy is recorded only as provenance parameter rank1_only, and only via the CLI (cli/convert/postprocess.py:204-229). The Python API leaves it unrecorded (converter.py:91-112).
- **Spec:** psm.md:12-15 mentions several hits per spectrum (engines, cross-links) only; psm.md:291 recommends a 'rank' score. Nothing says whether a psm view is rank-1 only or how to declare it. #378 asks 'rank-1 only?' and leaves it open. #232's 'rank joins the PK' is not in the spec.
- **Correction:** The spec does document PD's choice. converter-coverage.md:323 says 'ranks stay as the rank score, and --rank1-only keeps only first-ranked hits', and provenance.md:65 names rank1_only as a recorded conversion option (CLI only). The general gap remains: there is no view-level rule or footer key for the rank policy, and the Python API records nothing.
- **Reader impact:** PSM/ID counts and FDR re-estimation include alternative hits (#378: 20,001 of 109,700 PD TMT rows are rank>1) unless the reader filters a nested score.
- **Suggested:** Contract: add a typed psm.rank column (null = unknown) and a footer key psm_ranks=all|top1. Converters write rank 1 only by default or declare the policy.
- **Evidence:** `qpx/converters/proteome_discoverer/psm_adapter.py:37-40`; `qpx/converters/proteome_discoverer/reports.py:138-139`; `qpx/cli/convert/postprocess.py:204-229`; `docs/spec/psm.md:12-15`; `docs/spec/psm.md:291`; `bigbio/qpx#378 notes (github_history.md)`; `run: audit-run:pd/cases/rank2/out/pd.psm.parquet`

### PD-07: Duplicate engine hits keep one engine's row and drop the other's scores {#pd-07}

*psm* · inconsistent across producers · confirmed · checklist P6, P5

- **Now:** psm_adapter.py:38-40 keeps one row per (run, scan, sequence, modifications, charge). The losing engine's XCorr or Ions Score is discarded; the winner's Identifying Node is kept as cv_param (psm_adapter.py:76-82). A test pins this (tests/converters/test_proteome_discoverer_converter.py:136-148). The mzIdentML converter instead keeps every hit with a provided psm_id (identity.md:40).
- **Spec:** psm.md:12-15: several hits of one peptidoform and charge per spectrum need a provided psm_id or an extended composite. #232's decision, one merged row with every engine's scores in additional_scores, is not in the spec.
- **Reader impact:** The same multi-engine search gives one row from PD and several from mzIdentML. consensus_support (scores.md:75) cannot be derived from PD output.
- **Suggested:** Contract: one row per spectrum, peptidoform and charge, with all engines' scores merged in additional_scores and consensus_support written.
- **Evidence:** `qpx/converters/proteome_discoverer/psm_adapter.py:38-40`; `qpx/converters/proteome_discoverer/psm_adapter.py:76-82`; `tests/converters/test_proteome_discoverer_converter.py:136-148`; `docs/spec/psm.md:12-15`; `docs/spec/identity.md:40`; `docs/spec/scores.md:75`

### PD-08: pass_threshold = Confidence 'High', ignoring PSM Ambiguity and the undisclosed FDR level {#pd-08}

*psm* · ambiguous · partly confirmed (see correction) · checklist D3, P3

- **Now:** psm_adapter.py:68 sets True only for 'High', False for Medium or Low, and null when Confidence is blank or absent. A 'Rejected' PSM (loser of PD's ambiguity resolution) or a rank-2 PSM can be True; ambiguity survives only as cv_param psm_ambiguity (:79). PD's strict and relaxed target FDR values are neither read nor recorded. No PSM is dropped by threshold. A test pins High→True (tests/converters/test_proteome_discoverer_converter.py:115).
- **Spec:** psm.md:66: 'search engine or validator marks the PSM as passing its threshold; null when it does not say'; producer-guide.md:11. The threshold value, its level (PSM or peptide), and its interplay with rank and ambiguity are unspecified.
- **Correction:** Mapping High to True is documented (converter-coverage.md:323) and fits psm.md:66, which makes pass_threshold the validator's own mark. Not folding in rank or ambiguity therefore follows the spec. The real gaps: no field records the threshold value or level, and without a Confidence column pass_threshold is null while every row is kept, which conflicts with producer-guide.md:11.
- **Reader impact:** Filtering on pass_threshold keeps Rejected and rank-n PSMs. Cross-producer comparisons assume an FDR that is not known.
- **Suggested:** Contract: pass_threshold applies to rank-1, non-rejected PSMs at a declared FDR. Record the threshold and its level in provenance parameters.
- **Evidence:** `qpx/converters/proteome_discoverer/psm_adapter.py:68`; `qpx/converters/proteome_discoverer/psm_adapter.py:79`; `tests/converters/test_proteome_discoverer_converter.py:115`; `docs/spec/psm.md:66`; `docs/spec/producer-guide.md:11`; `run: audit-run:pd/cases/confidence_medium`

### PD-25: scan = First Scan only; scan_format never declared {#pd-25}

*psm* · diverges from spec · partly confirmed (see correction) · checklist P7, D8

- **Now:** scan=[First Scan] (reports.py:90,116); 'Last Scan' of merged spectra is ignored. Both writers are created without scan_format (psm_adapter.py:41; feature_adapter.py:43), and the footers lack it (verified). rt = 'RT [min]' × 60, the MS2 scan time (reports.py:119).
- **Spec:** producer-guide.md:23: 'Set scan_format'. scan.md:80-81: an undeclared value is read as a scan number. psm.md:50: rt is the MS2 scan's retention time in seconds.
- **Correction:** The missing scan_format is real. But scan.md:51 makes declaring it optional ('can declare'), and scan.md:80-81 reads an undeclared value as a scan number, which PD First Scan is for Thermo data. Only producer-guide.md:23 recommends setting it. This is a recommendation gap that is inconsistent across producers, not a divergence.
- **Reader impact:** Harmless for Thermo data. For non-Thermo PD input the meaning of First Scan is undeclared, and a merged-spectrum PSM points to only one of its scans.
- **Suggested:** Declare scan_format=scan. Record Last Scan as a cv_param when it differs from First Scan.
- **Evidence:** `qpx/converters/proteome_discoverer/reports.py:90`; `qpx/converters/proteome_discoverer/reports.py:116`; `qpx/converters/proteome_discoverer/psm_adapter.py:41`; `qpx/converters/proteome_discoverer/feature_adapter.py:43`; `docs/spec/producer-guide.md:23`; `docs/spec/scan.md:80-81`

### PD-26: Modifications column, static TMT/Carbamidomethyl included, becomes a canonical UNIMOD peptidoform {#pd-26}

*psm* · specified, worth keeping · confirmed · checklist D5, P9

- **Now:** reports.py:31-42 parses tokens such as 'N-Term(Prot)(Acetyl); K9(TMTpro)' via modification_tag (qpx/converters/ptm.py:321-339). Verified: Met-loss→UNIMOD:765, GG→121, Gln->pyro-Glu→28, TMT6plex→737, Label:13C(6)15N(2)→259, C-Term(Amidated)→-[UNIMOD:2]. Unknown names are kept with a warning (ptm.py:338). An unreadable token skips the row with only a count warning (reports.py:37; base.py:228-239; 'Skipped 1 / 2 rows').
- **Spec:** peptidoform.md:17-35: canonical form with fixed modifications inline, upper-case UNIMOD and 4-decimal mass shifts. peptidoform.md:43: names without a Unimod entry are kept as written.
- **Reader impact:** Peptidoforms compare with other producers. Silently skipped rows shrink the psm and feature views with no trace in the dataset.
- **Suggested:** Keep the rule. The contract should require skipped-row counts in provenance.
- **Evidence:** `qpx/converters/proteome_discoverer/reports.py:31-42`; `qpx/converters/ptm.py:321-339`; `qpx/converters/base.py:211-240`; `docs/spec/peptidoform.md:17-35`; `script: audit-run:pd/scripts/mods.py`

### PD-27: Engine score coverage is partial, and 'delta_score' carries no engine prefix {#pd-27}

*psm* · ambiguous · confirmed · checklist P5

- **Now:** reports.py:138-147 maps only Rank, XCorr, Comet XCorr, Ions Score, Delta Score/DeltaScore (→'delta_score') and SVM Score (→percolator_svm_score). Other engines' score columns are dropped. A generic 'PEP' column is annotated percolator:PEP MS:1001493 in the ontology output. Feature additional_scores repeat the chosen PSM's scores, including rank=1 (feature_adapter.py:126).
- **Spec:** scores.md:78: tool-prefixed score names. scores.md:58: 'rank' is PSM-level. psm.md:44: PEP is a top-level field.
- **Reader impact:** 'delta_score' holds different engines' scores under one name. A feature-level 'rank' and engine scores describe one PSM, not the feature.
- **Suggested:** Prefix every engine score, keep unmapped numeric score columns under prefixed names, and drop rank from feature scores.
- **Evidence:** `qpx/converters/proteome_discoverer/reports.py:138-147`; `qpx/converters/proteome_discoverer/feature_adapter.py:126`; `docs/spec/scores.md:58`; `docs/spec/scores.md:78`; `run: audit-run:pd/cases/rank2/cli_out/proteome_discoverer.ontology.parquet`

### PD-32: psm protein_accessions = every mapped protein from 'Protein Accessions' {#pd-32}

*psm* · specified, worth keeping · confirmed · checklist P8

- **Now:** psm_adapter.py:66 splits 'Protein Accessions' (all proteins the peptide maps to, in PD's order) into protein_accessions, null when blank. Master proteins go only to feature.anchor_protein. No razor or decoy prefixes are applied.
- **Spec:** psm.md:57: protein accessions of all proteins the peptide maps to. #78: no protein inference in PSM tables.
- **Reader impact:** Consistent with the spec; the order carries no meaning.
- **Suggested:** Keep: in the contract, psm.protein_accessions lists all mapped proteins, unordered and without inference.
- **Evidence:** `qpx/converters/proteome_discoverer/psm_adapter.py:66`; `docs/spec/psm.md:57`

### PD-M3: PSM dedup on raw modification text can emit duplicate psm_id {#pd-m3}

*psm* · diverges from spec · found by the verifier · checklist P6, D5

- **Now:** psm_adapter.py:39-40 partitions on the raw _modifications text (reports.py:88), but psm_id hashes the canonical peptidoform. One scan with 'N-Term(Prot)(Acetyl)' from one engine and 'N-Term(Acetyl)' from another gave two psm rows with one psm_id. The writer only warns. Whether PD ever spells one modification two ways is unverified.
- **Spec:** psm.md:6-15: psm_id is the primary key, and the composite names one row. #378 rule 5 (open): a producer must not write two psm rows with the same (run_file_name, scan, peptidoform, charge).
- **Reader impact:** The PSM primary key is not unique, so joins on psm_id or the natural key fan out.
- **Suggested:** Deduplicate on the canonical peptidoform in both views. Contract: the PK uniqueness check fails conformance.
- **Evidence:** `qpx/converters/proteome_discoverer/psm_adapter.py:39-40`; `qpx/converters/proteome_discoverer/reports.py:88`; `docs/spec/psm.md:6-15`; `run: work/verify-pd/cases/psmcollide (warning 'Primary key (psm_id) has 1 duplicate')`

## mzIdentML and CDAP converters

### MZID_CDAP-01: mzIdentML exports every SII rank by default; rank only as additional score {#mzid_cdap-01}

*psm* · unspecified · confirmed · checklist P1, P6, D3

- **Now:** Every SpectrumIdentificationItem becomes a row (qpx/converters/mzidentml/psm_records.py:52-71); rank stored only as additional_scores 'rank' (psm_records.py:140); --rank1-only (cli/convert/mzidentml.py:40-45, off by default) drops rank>1 only (psm_records.py:63). PXD054720: 13,386 rows for 9,707 spectra, ranks 1-12, 3,679 rank>1 rows (1,869 pass_threshold=false). Same scan at another charge = another row (synthetic).
- **Spec:** psm.md:3-5 one match per row; rank is only a 'recommended additional score' (psm.md:291, scores.md:58). Nothing says whether rank>1 belongs in psm or how it joins identity (#232 'rank joins PK' not in spec).
- **Reader impact:** Unfiltered PSM/peptide counts inflated ~38% on PXD054720; filtering by rank needs a nested-score scan and is producer-dependent.
- **Suggested:** Contract: psm default rank-1 only, or typed nullable `rank` column plus footer key (e.g. psm_ranks=1|all); rank joins the identity composite when >1 kept.
- **Evidence:** `qpx/converters/mzidentml/psm_records.py:46-73,140`; `qpx/cli/convert/mzidentml.py:40-45`; `qpxc convert mzidentml on tests/examples/mzidentml/PXD054720 -> audit-run:mzid_cdap/pxd; duckdb rank/pass_threshold counts`; `docs/spec/psm.md:291`

### MZID_CDAP-02: Rank-1 ties and chimeric hits indistinguishable; CDAP flags ties only in cv_params {#mzid_cdap-02}

*psm* · inconsistent across producers · partly confirmed (see correction) · checklist P2, P1, X1

- **Now:** mzid: several rank=1 SIIs each become a row with no tie/chimeric marker (psm_records.py:62-71); synthetic SIR with PEPTIDEK and PEPTLDEK, equal scores, both rank 1 -> 2 rows. CDAP: AmbiguousMatch=1 ties (fixture scan 5646: SCHTGIGR vs SCHTGLGR, I/L) -> separate rows, only cv_param ambiguous_match='1' (cdap/psm_adapter.py:120-127). Neither represents co-isolated (chimeric) second peptides.
- **Spec:** Nothing for psm; only feature.identification_type 'secondary' for chimeric (feature.md:145).
- **Correction:** The gist holds, but the CDAP flag is set at cdap/psm_adapter.py:26-27 (_CV_PARAM_COLUMNS). Lines 120-127 are the SQL projection. In the fixture, scans 5646 and 5652 are I/L ties between targets. Scan 1901 is a tie between two decoys that differ in where the acetyl and TMT sit (a localization tie, not I/L).
- **Reader impact:** Cannot tell 'engine undecided, one is wrong' (tie) from 'two peptides co-fragmented, both right' (chimeric); spectral counts/FDR double-count ties.
- **Suggested:** Typed psm field (e.g. hit_type: unique|tie|chimeric) or rank+tie flag; CDAP AmbiguousMatch maps to tie; MaxQuant SECPEP to chimeric.
- **Evidence:** `qpx/converters/mzidentml/psm_records.py:62-71`; `qpx/converters/cdap/psm_adapter.py:120-127`; `audit-run:mzid_cdap/synth/synth.mzid (SII1_1/SII1_2) converted output`; `audit-run:mzid_cdap/cdap/cdap.psm.parquet rows with ambiguous_match=1`

### MZID_CDAP-03: mzIdentML maps 'PSM-level p-value' (MS:1002352) into posterior_error_probability {#mzid_cdap-03}

*psm* · contradictory · confirmed · checklist P5

- **Now:** _PEP_ACCESSIONS={MS:1001493, MS:1002352} (mzidentml/scores.py:18,136-148); cv_terms.py:48 comments MS:1002352 as 'global PEP', but bundled PSI-MS ontology names it 'PSM-level p-value'. Synthetic SII with p-value 0.3 -> posterior_error_probability=0.3. 'PSM-level local FDR' (MS:1002351) unused. PEP term also duplicated in additional_scores (scores.py:110-133, e.g. percolator_pep).
- **Spec:** psm.md:44 PEP = probability the PSM is incorrect. converter-coverage.md:274 endorses MS:1002352 as a 'PSI-MS PEP term' (wrong).
- **Reader impact:** p-values silently read as PEP; PEP-based filters, rollups and FDR estimates are wrong for files carrying that term.
- **Suggested:** PEP only from genuine PEP terms; p-values stay in additional_scores; fix cv_terms.py:48 and converter-coverage.md:274.
- **Evidence:** `qpx/converters/mzidentml/scores.py:18,136-148`; `qpx/core/cv_terms.py:47-48`; `qpx/core/ontology/data/psi_ms.parquet: MS:1002352='PSM-level p-value'`; `docs/spec/converter-coverage.md:274`; `synthetic SII1_1 output pep=0.3`

### MZID_CDAP-04: mzIdentML marks target+decoy peptides as decoy, contrary to spec {#mzid_cdap-04}

*psm* · diverges from spec · confirmed · checklist P4

- **Now:** is_decoy = any(PeptideEvidence isDecoy) (psm_records.py:138; crosslinks.py:157-158). Synthetic KEDITPEP with evidences DECOY_P00001 (isDecoy=true) and P00002 (false) -> is_decoy=true. CDAP uses all() (cdap/base_adapter.py:146), i.e. target.
- **Spec:** psm.md:45, feature.md:60: a peptide matching both target and decoy proteins is a target match.
- **Reader impact:** Same peptide is decoy in mzid output, target in CDAP/OpenMS output; FDR re-estimation differs by producer.
- **Suggested:** Decoy only when every evidence is decoy, for psm/feature/pg alike, stated once in the contract.
- **Evidence:** `qpx/converters/mzidentml/psm_records.py:138`; `qpx/converters/mzidentml/crosslinks.py:157-158`; `qpx/converters/cdap/base_adapter.py:146`; `synthetic SII1_3 output is_decoy=True`; `docs/spec/psm.md:45`

### MZID_CDAP-08: Multi-engine mzIdentML gives duplicate identity composites and fails conformance {#mzid_cdap-08}

*psm* · contradictory · partly confirmed (see correction) · checklist P6, P5, X1

- **Now:** Each SpectrumIdentificationList converted independently; no column names list/engine (only SII id in cv_params, psm_records.py:228-235). Synthetic: MS-GF+ and Comet hits on scan 100 PEPTIDEK z2 -> 2 rows, distinct psm_id; `qpxc validate --conformance` ERROR duplicate_composite (core/data/footer_validation.py:157-171). Comet row rt=null since rt is per SIR.
- **Spec:** psm.md:12-15 and identity.md:40: producer-supplied psm_id handles several hits of one peptidoform per spectrum. versioning.md:86-87: one composite value names one row identity, also for provided ids. #232 merge-per-spectrum not in spec.
- **Correction:** The duplicate_composite ERROR reproduces, and the problem is real. The Comet row's rt=null is an artifact of the synthetic file: SIRc1 has no 'scan start time', and RT is read per SIR (psm_adapter.py:323-339). The contradiction is inside the spec: psm.md:12-15 allows provided ids for multi-engine hits, but versioning.md:86-87 and the validator forbid them.
- **Reader impact:** Legitimate conversions fail validation; PSM counts double-count multi-engine hits; scores cannot be attributed to an engine.
- **Suggested:** Choose: one row per (spectrum, peptidoform, charge) with engine-prefixed scores (#232), or an engine/search id column joining the composite.
- **Evidence:** `qpx/converters/mzidentml/psm_records.py:221-235`; `qpx/core/data/footer_validation.py:157-171`; `audit-run:mzid_cdap/synth/synth.mzid + qpxc validate --dataset-path synth/out --conformance`; `docs/spec/psm.md:12-15`; `docs/spec/versioning.md:86-87`

### MZID_CDAP-10: mzIdentML pass_threshold is engine-defined and its threshold is not recorded {#mzid_cdap-10}

*psm* · ambiguous · partly confirmed (see correction) · checklist D3, P3, D8

- **Now:** pass_threshold = SII@passThreshold on every row incl. rank>1 (psm_records.py:146); missing attribute -> false, not null (psm_adapter.py:304); never used to filter. SIP <Threshold> (PXD054720 Mascot:SigThreshold 0.01 homology; XL fixture crosslinked PSM-level global FDR 0.05) not written to provenance (converter.py:290-316). Peptide-level MS:1002500 and MS:1003339 pass flags dropped (cv_terms.py:102-110; non-numeric).
- **Spec:** psm.md:66 'passes its threshold', null when the source does not say; producer-guide.md:11. Threshold kind/level unspecified.
- **Correction:** The mzIdentML 1.2 XSD makes passThreshold (and rank) required on SII, so 'missing attribute -> false' only affects schema-invalid files. The rest holds: <Threshold> is never read (the provenance code collects only the enzyme and modifications), and MS:1002500/MS:1003339 are dropped.
- **Reader impact:** pass_threshold=true means Mascot p<0.01 here, 5% XL FDR there, 1% PSM FDR elsewhere; readers treat it as one uniform FDR filter.
- **Suggested:** Record threshold term, value and level (provenance fdr_filtering or footer); absent attribute -> null.
- **Evidence:** `qpx/converters/mzidentml/psm_records.py:146`; `qpx/converters/mzidentml/psm_adapter.py:304`; `qpx/converters/mzidentml/converter.py:290-316`; `tests/examples/mzidentml/xl/scores_and_thresholds_1_3_0_draft.mzid.gz lines 245-255`; `audit-run:mzid_cdap/pxd/PXD054720.provenance.parquet`

### MZID_CDAP-11: CDAP input looks FDR-filtered but filtering state is recorded nowhere {#mzid_cdap-11}

*psm* · ambiguous · confirmed · checklist P3, D3, D8

- **Now:** Fixture .psm rows all have Qvalue<=0.0089 incl. 10 decoys, i.e. pre-filtered by CDAP (inferred, not documented). pass_threshold is null on all rows; provenance MS-GF+ step parameters=None (cdap/converter.py:191-201); format_conversion params only psm_dir (observed). Converter applies no filter.
- **Spec:** producer-guide.md:11 'keep rows above your FDR threshold only if you say so with pass_threshold' (ambiguous for pre-filtered inputs); provenance.md:65.
- **Reader impact:** Cannot tell CDAP data (truncated at ~1% PSM FDR) from unfiltered mzid/MaxQuant data; decoy-based FDR recomputation on truncated lists is biased.
- **Suggested:** For known pre-filtered sources set pass_threshold=true and record level/threshold (or 'unknown') in provenance/footer.
- **Evidence:** `qpx/converters/cdap/converter.py:191-217`; `awk over tests/examples/cdap/*.psm Qvalue max=0.0089`; `audit-run:mzid_cdap/cdap/cdap.provenance.parquet`

### MZID_CDAP-12: psm.qvalue picked by cvParam order when engine and rescorer q-values coexist {#mzid_cdap-12}

*psm* · ambiguous · confirmed · checklist D4, P5

- **Now:** psm.qvalue = first additional_score matching is_psm_qvalue, in cvParam order (writers/derived.py:15-21,40-55; mzidentml/scores.py:110-133). Synthetic SII with percolator:Q value 0.004 before MS-GF:QValue 0.02 -> qvalue 0.004; reversed order gives 0.02. PXD054720 (Mascot) has no q-value: qvalue null on all rows.
- **Spec:** scores.md:38-42 'first match of a row wins'.
- **Reader impact:** psm.qvalue is either the engine's target-decoy q or the rescored one depending on file order; psm.confident() differs between equivalent files.
- **Suggested:** Define precedence (post-processor/validator over engine) and record the source score name (e.g. footer psm_qvalue_source).
- **Evidence:** `qpx/writers/derived.py:15-55`; `synthetic SII2_1 output qvalue=0.004`; `docs/spec/scores.md:38-42`

### MZID_CDAP-17: CDAP rt is precursor half-elution time and 0.00 sentinel becomes rt=0 {#mzid_cdap-17}

*psm* · diverges from spec · partly confirmed (see correction) · checklist P7, F10

- **Now:** psm.rt and feature.rt = RTAtPrecursorHalfElution (cdap/psm_adapter.py:212,263; feature_adapter.py:321,365), the precursor SIC half-elution time, not the MS2 scan time. CDAP writes 0.00 when PrecursorArea=0, kept as rt=0 (fixture: 8/100 PSMs, 6/94 features). Unit assumed seconds (column_mappings.yaml:209).
- **Spec:** psm.md:50 rt = MS2 scan's retention time; feature.md:65 'precursor retention time'.
- **Correction:** The PSM rt comes from cdap/psm_adapter.py:118 and :169; the file has only 210 lines, so 212 and 263 do not exist. The behavior is confirmed: 8/100 PSMs and 6/94 features have rt=0, and every 0.00 row has PrecursorArea=0. column_mappings.yaml:209 declares rt_unit 'second' and :219 maps rt.
- **Reader impact:** RT joins to spectra and RT-prediction training are off; rt=0 rows look like early eluters.
- **Suggested:** psm.rt = MS2 scan RT or null; precursor elution time only on feature; write null for the 0.00 sentinel.
- **Evidence:** `qpx/converters/cdap/psm_adapter.py:212,263`; `qpx/converters/cdap/feature_adapter.py:321,365`; `qpx/converters/column_mappings.yaml:209,219`; `awk: all 8 RTAtPrecursorHalfElution=0.00 rows have PrecursorArea 0`

### MZID_CDAP-18: Missing observed m/z and charge become 0; corrected vs instrument m/z unstated {#mzid_cdap-18}

*psm* · diverges from spec · partly confirmed (see correction) · checklist P7, D6

- **Now:** CDAP observed_mz via _finite_float default 0.0 (cdap/psm_adapter.py:255; feature_adapter.py:183; base_adapter.py:304-318). mzid experimentalMassToCharge and chargeState default '0' (mzidentml/psm_adapter.py:299-300). CDAP uses QueryPrecursorMz/QueryCharge and drops OriginalPrecursorMz/OriginalCharge (equal in fixture). Not triggered in fixtures.
- **Spec:** psm.md:47 observed_mz null when not reported; which precursor m/z (instrument vs search-corrected) is not specified.
- **Correction:** CDAP lines are psm_adapter.py:161 and feature_adapter.py:183 (via base_adapter.py:304-318); the cited psm_adapter 207-211 and 255 do not exist. For mzid, chargeState and experimentalMassToCharge are use=required in the XSD, so the 0 defaults hit only invalid files. CDAP QueryCharge is filtered NOT NULL for features (feature_adapter.py:126).
- **Reader impact:** 0.0 m/z or charge 0 pass as measurements; mass errors computed from them are garbage.
- **Suggested:** Null for missing values; contract defines observed_mz as the m/z the search used, with the instrument value in cv_params when different.
- **Evidence:** `qpx/converters/cdap/psm_adapter.py:207-211,255`; `qpx/converters/cdap/base_adapter.py:304-318`; `qpx/converters/mzidentml/psm_adapter.py:299-300`

### MZID_CDAP-19: mzIdentML nativeID parsing drops vendor scan components {#mzid_cdap-19}

*psm* · diverges from spec · confirmed · checklist P7, D8

- **Now:** _parse_native_id keeps only the scan=/index=/spectrum= number (psm_records.py:253-269), tested: 'frame=120 scan=475'->[475], 'function=10 process=1 scan=345'->[345], 'controllerType=5 controllerNumber=1 scan=7'->[7], all declared 'scan'; Sciex 'sample=1 period=1 cycle=2740 experiment=10'->[]. XL fixture SIR cvParam 'peak list scans' 40560 ignored; scan=[26630] (index).
- **Spec:** scan.md:24-45 keep every component ([120,475], [10,1,345], [5,1,7]) with scan_format nativeId.
- **Reader impact:** Bruker/Waters PSMs from different frames/functions collide on scan; Sciex rows get empty scan; joins to the mz view break.
- **Suggested:** Parse all numeric components per scan.md; prefer MS:1000797 'peak list scans' when present.
- **Evidence:** `qpx/converters/mzidentml/psm_records.py:238-269`; `python call of _parse_native_id on vendor IDs`; `docs/spec/scan.md:24-45`; `tests/examples/mzidentml/xl/scores_and_thresholds_1_3_0_draft.mzid.gz line 344`

### MZID_CDAP-29: CDAP aggregates PSMs into features but never fills psm.feature_id {#mzid_cdap-29}

*psm* · unspecified · partly confirmed (see correction) · checklist P10, F2

- **Now:** Each CDAP feature is built from PSMs with the same raw sequence, charge and file, yet psm.feature_id is null for all rows (cdap/psm_adapter.py:274-290); qpxc summary 'feature<-psm: 0 features with >=1 psm (100 psms without a feature)'.
- **Spec:** psm.md:75 feature_id optional; converter-coverage.md:285 'Not implemented: features are aggregated from PSMs, but the link is not recorded'.
- **Correction:** The record is built at cdap/psm_adapter.py:180-196 and never sets feature_id; the file has 210 lines, so 274-290 do not exist. The behavior is confirmed and documented as not implemented (converter-coverage.md:285).
- **Reader impact:** Cannot recover which scans were summed into a feature, nor which tied PSMs were excluded.
- **Suggested:** Contract: a producer that builds features from PSMs MUST fill psm.feature_id (deterministic for CDAP).
- **Evidence:** `qpx/converters/cdap/psm_adapter.py:274-290`; `qpxc convert cdap stdout summary`; `docs/spec/converter-coverage.md:285`

### MZID_CDAP-31: CDAP 3-decimal masses padded to 4-decimal mass shifts {#mzid_cdap-31}

*psm* · ambiguous · confirmed · checklist D5, P9

- **Now:** CDAP masses (3 decimals) map to UNIMOD within tolerance (cdap/peptidoform.py:32-42); unmatched masses are padded: '+123.456' -> '[+123.4560]' (ptm.py:217-219; tested). Fixed mods and labels are inline (fixture C[UNIMOD:4], [UNIMOD:737]-). mzid uses Modification monoisotopicMassDelta via canonical_mod_tag (mzidentml/psm_adapter.py:146-161).
- **Spec:** peptidoform.md:27-30 unknown modification = signed monoisotopic mass shift with four decimals.
- **Reader impact:** The same unknown modification becomes '+123.4560' (CDAP) vs '+123.4567' elsewhere: different peptidoforms and ids across producers.
- **Suggested:** Define how source precision below 4 decimals is written (and matched), e.g. tolerance-based mass-shift comparison.
- **Evidence:** `qpx/converters/cdap/peptidoform.py:32-42`; `qpx/converters/ptm.py:217-219`; `python cdap_to_proforma('+229.163PEPTK+123.456')`; `docs/spec/peptidoform.md:27-30`

### MZID_CDAP-32: mzIdentML cv_params carry positions of the first protein only {#mzid_cdap-32}

*psm* · ambiguous · confirmed · checklist P8, X1

- **Now:** cv_params get pre/post/start/end (non-CV names) of the first PeptideEvidence only (psm_records.py:201-216); protein_accessions lists all proteins, deduplicated in file order, decoy markers kept (psm_records.py:116-131); protein_starts only for cross-links (crosslinks.py:161-174).
- **Spec:** scores.md:89-127 cv_params for CV annotations; psm.md:57 protein_accessions = all proteins; order not specified.
- **Reader impact:** Positions silently refer to protein_accessions[0]; per-protein positions for shared peptides are only in pepmap.
- **Suggested:** Positions only in pepmap/protein_starts (parallel to protein_accessions); state protein_accessions order.
- **Evidence:** `qpx/converters/mzidentml/psm_records.py:116-131,191-218`; `qpx/converters/mzidentml/crosslinks.py:161-174`

### MZID_CDAP-M1: PSM.confident() ignores pass_threshold and rank: empty for Mascot mzid {#mzid_cdap-m1}

*psm* · ambiguous · found by the verifier · checklist D3, D4, P1, P3

- **Now:** psm.py:23-25 keeps is_decoy=false AND qvalue<=x and never reads pass_threshold or rank. PXD054720 (Mascot) has qvalue null on 13,386 rows, so confident() returns 0 rows although 11,517 rows have pass_threshold=true. Rank>1 rows that carry a q-value pass (synthetic rank-2 has qvalue 0.5, so it only fails on the value). CDAP has pass_threshold null and qvalue filled.
- **Spec:** scores.md:38: typed q-values are 'the portable way to filter', used by psm.confident(). psm.md:66 defines pass_threshold. No rule combines them or says what to do when a source has no q-value.
- **Reader impact:** The 'confident PSMs' API gives nothing for Mascot/X!Tandem-style mzid and can keep non-top hits; the same call means different things per producer.
- **Suggested:** Contract confidence predicate: qvalue when present, else pass_threshold; rank=1 unless the row is declared chimeric; document the fallback.
- **Evidence:** `qpx/core/data/psm.py:23-25`; `verify-mzid_cdap/pxd: count(qvalue)=0, pass_threshold true=11,517`; `docs/spec/scores.md:38`; `docs/spec/psm.md:66`

### MZID_CDAP-M6: CDAP phosphoRS localization not used for peptidoform or candidates {#mzid_cdap-m6}

*psm* · unspecified · found by the verifier · checklist P9, D5

- **Now:** The peptidoform comes from MS-GF+ PeptideSequence (base_adapter.py:225-242). PhosphoRSPeptide and FullyLocalized are kept only as cv_param strings (psm_adapter.py:26-33), with no modification_candidates or localization scores. Features group on PeptideSequence. The fixture has no phospho columns, so this is derived from code.
- **Spec:** converter-coverage.md:251-258 (site localization table) has no CDAP row; modifications.md defines modification_candidates; nothing says which CDAP placement is canonical.
- **Reader impact:** For CPTAC phosphoproteomes the site in peptidoform may be MS-GF+'s, not the reported phosphoRS one; site-level analyses differ from CDAP's own reports.
- **Suggested:** Contract says which localization defines the peptidoform; CDAP maps phosphoRS probabilities to modification_candidates.
- **Evidence:** `qpx/converters/cdap/psm_adapter.py:26-33`; `qpx/converters/cdap/base_adapter.py:225-242`; `docs/spec/converter-coverage.md:251-258`

## Shared writer and validation code

### INFRA-05: psm.qvalue name rule is a denylist and admits higher-level q-values {#infra-05}

*psm* · ambiguous · confirmed · checklist D4, P5

- **Now:** is_psm_qvalue = name ends in q_?value and has no excluded token (peptide, peptidoform, protein, pg, gg, global, precursor, lib, translated, glycan, site) (derived.py:15-21). A probe classifies these as PSM-level: PSI-MS proteoform-level Q-value (MS:1003129), isoform-level Q-value (MS:1003133), ProSight proteoform/isoform Q-value, genes_qvalue, pep_qvalue, msgf_pep_qvalue, match_qvalue, ion_qvalue, run_qvalue. mzIdentML keeps every numeric cvParam under its normalized name (mzidentml/scores.py:110-133), so such terms reach psm.qvalue.
- **Spec:** scores.md:42 excludes only peptide, peptidoform, protein, pg, gg and global (plus precursor). The code also excludes lib, translated, glycan and site, which the spec does not list.
- **Reader impact:** A proteoform-, gene- or peptide-level q-value can show up as psm.qvalue. The spec list and the code list differ, so producers in other languages cannot reproduce the result.
- **Suggested:** Use an explicit allowlist of PSM-level names (qvalue, q_value, psm_qvalue, psm_level_q_value, plus named engine terms such as MS:1002054, MS:1001491, MS:1002354) and publish it as the normative list.
- **Evidence:** `qpx/writers/derived.py:15-21`; `qpx/converters/mzidentml/scores.py:110-133`; `docs/spec/scores.md:42`; `probe script audit-run:infra/scripts/writer_exp.py; PSI-MS term probe over bundled ontology`

### INFRA-06: OpenMS writes PID primary score as 'q-value' even under peptide-level FDR {#infra-06}

*psm* · ambiguous · confirmed · checklist D4, P5

- **Now:** _psm_additional_scores names the primary score 'q-value' whenever the score type is a q-value or FDR (openms_consensus/psm_adapter.py:216-219,363), including when peptide_level is true (:276-288). Under peptide-level Percolator FDR, non-best PSMs carry 1.0 placeholders (feature_adapter.py:166-169). The writer then copies 'q-value' into psm.qvalue, taking the first matching entry (derived.py:40-55). A score type of 'fdr' is also treated as a q-value.
- **Spec:** psm.md:51: psm.qvalue is a 'PSM-level q-value'. scores.md:38: 'the first match of a row wins'. feature.md:130 forbids the 1.0 placeholders only for features.
- **Reader impact:** In peptide-level FDR mode psm.qvalue can be a peptide-level value or a 1.0 placeholder, which would wrongly drop valid PSMs. Not verified: no peptide-level fixture.
- **Suggested:** When the identification run used peptide-level FDR, store only peptide_qvalue and leave psm.qvalue null. Never promote 'fdr' (non-monotone) to qvalue. Refuse rows with two PSM-level candidates instead of first-wins.
- **Evidence:** `qpx/converters/openms_consensus/psm_adapter.py:204-245`; `qpx/converters/openms_consensus/psm_adapter.py:276-288`; `qpx/converters/openms_consensus/psm_adapter.py:363`; `qpx/converters/openms_consensus/feature_adapter.py:163-191`; `qpx/writers/derived.py:40-55`

### INFRA-09: CV_PEP_GLOBAL is MS:1002352, a PSM-level p-value, used as PEP {#infra-09}

*psm* · contradictory · confirmed · checklist P5

- **Now:** cv_terms.py:48 defines CV_PEP_GLOBAL='MS:1002352' with the comment 'PSM-level global FDR (global PEP)'. In the bundled PSI-MS, MS:1002352 is 'PSM-level p-value'. mzidentml/scores.py:18,136-148 writes the first MS:1001493 or MS:1002352 value into posterior_error_probability. Real PEP terms are not recognized: MS-GF:PEP MS:1002056, MaxQuant:PEP MS:1001901, Byonic:PEP MS:1002265, Andromeda:PEP MS:1002995, OpenMS:ConsensusID PEP MS:1003113. Those stay only in additional_scores.
- **Spec:** converter-coverage.md:274 repeats 'A PSI-MS PEP term (MS:1001493, MS:1002352)'. psm.md:286: PEP = P(incorrect).
- **Reader impact:** A p-value can be presented as PEP, and true PEPs from non-Percolator engines leave the PEP column null.
- **Suggested:** Make PEP an explicit accession allowlist (1001493, 1002056, 1001901, 1002265, 1002995, 1003113) and remove MS:1002352. Fix cv_terms and converter-coverage.md:274.
- **Evidence:** `qpx/core/cv_terms.py:47-48`; `qpx/converters/mzidentml/scores.py:18`; `qpx/converters/mzidentml/scores.py:136-148`; `docs/spec/converter-coverage.md:274`; `PSI-MS lookup via qpx.core.ontology: MS:1002352='PSM-level p-value'`

### INFRA-23: PSM identity has no rank; equal composites get distinct ids {#infra-23}

*psm* · contradictory · confirmed · checklist P1, P6

- **Now:** The default PSM composite is run_file_name, scan, peptidoform, charge (identity.md:29). Duplicate ids give only a warning at close (base.py:557-569). mzIdentML supplies its own ids from SII ids but declares the default composite. Fixture run: 108 composite keys cover 239 rows with distinct psm_ids, the same peptidoform at ranks 1-4 of one spectrum, sometimes with different pass_threshold.
- **Spec:** serialization.md:304-305: 'rows with equal composite values carry the same id, also when the producer supplies its own ids'. identity.md:40 admits the footer does not reproduce the mzIdentML ids. #229: rank joins the PK if ranks are kept, but this is not in the spec.
- **Reader impact:** A join on (run, scan, peptidoform, charge) duplicates rows. The footer composite does not identify a row. --check-identity cannot verify these ids.
- **Suggested:** When several ranks or SIIs per spectrum are kept, the declared composite must include rank (or spectrum_identification_item_id). Duplicate ids become a write error.
- **Evidence:** `docs/spec/identity.md:29,40`; `docs/spec/serialization.md:304-305`; `qpx/writers/base.py:535-569`; `run: qpxc convert mzidentml PXD054720 -> 108 composites / 239 rows with distinct psm_id`

### INFRA-24: Rank is untyped and the rank policy differs per producer {#infra-24}

*psm* · inconsistent across producers · confirmed · checklist P1, P2

- **Now:** Rank lives only in additional_scores 'rank', and only for mzIdentML (psm_records.py:140) and PD (reports.py:139). mzIdentML keeps all ranks by default (rank up to 12; 1,758 spectra with more than one row; option --rank1-only, psm_records.py:63). PD keeps all ranks and deduplicates only identical peptidoforms per scan (psm_adapter.py:39-40). The PD feature view uses only rank 1 (feature_adapter.py:16-20). Without a rank, rank-2 alternatives cannot be told apart from chimeric co-identifications.
- **Spec:** psm.md:291 recommends a 'rank' additional score. scores.md:58 defines it. No rule says whether the psm view holds rank 1 only, and #378 asks this.
- **Reader impact:** PSM counts, spectral counts and 'identified spectra' differ by producer. Rank-2+ hits can be counted as identifications.
- **Suggested:** Add a typed psm.rank (1 = best) and a footer rank_scope (rank1 or all). Chimeric co-identifications get their own marker instead of rank.
- **Evidence:** `qpx/converters/mzidentml/psm_records.py:63,140`; `qpx/converters/proteome_discoverer/psm_adapter.py:39-40`; `qpx/converters/proteome_discoverer/feature_adapter.py:16-20`; `docs/spec/psm.md:291`; `run: mzid rank max 12; 1758 multi-row spectra`

### INFRA-25: pass_threshold means a different threshold in each producer {#infra-25}

*psm* · inconsistent across producers · confirmed · checklist D3, P3

- **Now:** mzIdentML: the SII passThreshold attribute, which here is Mascot's own threshold, true for 11,517 of 13,386 PSMs (psm_adapter.py:304, psm_records.py:146). PD: Confidence=='High', so Medium becomes False (proteome_discoverer/psm_adapter.py:68). Other converters: null. The writer checks nothing against qvalue: qvalue 0.001 with pass_threshold False was accepted.
- **Spec:** psm.md:66: 'Whether the search engine or validator marks the PSM as passing its threshold'; the threshold, its level and its value are not recorded. producer-guide.md:11: 'Keep rows above your FDR threshold only if you say so with pass_threshold' (unclear).
- **Reader impact:** Filtering on pass_threshold applies an engine score threshold for one producer, a 1% FDR for another and nothing for the rest, and the file does not say which.
- **Suggested:** Record the threshold behind pass_threshold (kind, level, value) in the footer or provenance. Rewrite producer-guide.md:11 to say rows that fail the threshold may be kept and flagged false.
- **Evidence:** `qpx/converters/mzidentml/psm_adapter.py:304`; `qpx/converters/proteome_discoverer/psm_adapter.py:68`; `docs/spec/psm.md:66`; `docs/spec/producer-guide.md:11`; `run: mzid pass_threshold true 11517/13386`

### INFRA-31: Score names are not normalized at write; OpenMS writes 'q-value' {#infra-31}

*psm* · diverges from spec · confirmed · checklist P5

- **Now:** derived.py:47 normalizes names only for matching; additional_scores keeps the raw names. OpenMS-consensus writes 'q-value' (psm_adapter.py:218,224): on the LFQ and isobaric fixtures 115/115 and 4/4 PSMs carry 'q-value'. qpxc validate gives only a warning.
- **Spec:** scores.md:23: all score names are snake_case; 'qpxc validate warns on other names'.
- **Reader impact:** SQL on score_name='q_value' or 'qvalue' misses OpenMS rows, and cross-producer score queries break.
- **Suggested:** The writer normalizes score_name (normalize_score_name) on every path, or rejects non-snake_case names.
- **Evidence:** `qpx/writers/derived.py:40-55`; `qpx/converters/openms_consensus/psm_adapter.py:216-224`; `docs/spec/scores.md:23`; `qpxc validate --file audit-run:infra/omslfq/openms.psm.parquet -> 115 non-snake_case`

### INFRA-M3: The threshold behind pass_threshold is declared in the mzIdentML source but dropped {#infra-m3}

*psm* · unspecified · found by the verifier · checklist D3, P3, D8

- **Now:** F001234.mzid has SpectrumIdentificationProtocol <Threshold> Mascot:SigThreshold=0.01 and SigThresholdType=homology. The converter reads only the SII passThreshold attribute (mzidentml/psm_adapter.py:304); grep 'Threshold' finds nothing else. The provenance Mascot step holds only enzyme/fixed_mod/variable_mod. 871 rank-2 PSMs have pass_threshold=true.
- **Spec:** psm.md:66 does not ask for the threshold. provenance.md:3 promises 'what FDR was used?'.
- **Reader impact:** pass_threshold=true here means 'Mascot homology p<0.01 per hit', not an FDR. The information to say so existed but was discarded.
- **Suggested:** Copy the protocol <Threshold> cvParams into the search step's provenance parameters (key e.g. pass_threshold_basis) and into a psm footer key.
- **Evidence:** `zcat tests/examples/mzidentml/PXD054720/F001234.mzid.gz | grep -A3 '<Threshold'`; `qpx/converters/mzidentml/psm_adapter.py:304`; `run: work/verify-infra/mzid provenance step 1 params`

## Specification text

### SPEC-11: 'global_qvalue' means PSM-level or precursor-level depending on view {#spec-11}

*psm* · contradictory · confirmed · checklist D4, P5

- **Now:** Writer test: a PSM with only global_qvalue (as psm.md recommends) gets psm.qvalue=null. On a feature, global_qvalue is copied to feature.global_qvalue (precursor, experiment). The feature.md:228 example uses it. The scope of psm.qvalue (run or experiment) is never stated.
- **Spec:** psm.md:291 recommends global_qvalue as the 'experiment-level PSM q-value'. scores.md:42 excludes 'global' from psm.qvalue. scores.md:45 sends global_qvalue to the precursor-level feature column. ontology.md:277 maps it to MS:1002350, PSM-level global FDR.
- **Reader impact:** psm.confident() ignores producers that follow psm.md, and one score name has two levels.
- **Suggested:** Use distinct names (psm_global_qvalue, precursor_global_qvalue), define the scope of psm.qvalue, and fix the psm.md recommendation.
- **Evidence:** `docs/spec/psm.md:291`; `docs/spec/scores.md:42,45,56`; `docs/spec/ontology.md:277`; `docs/spec/feature.md:226-232`; `audit-run:spec/ex/c.psm.parquet`

### SPEC-16: PSM filtering state and pass_threshold semantics are undeclared {#spec-16}

*psm* · unspecified · partly confirmed (see correction) · checklist P3, D3, D8

- **Now:** Only PD fills pass_threshold (Confidence=High, converter-coverage.md:323). The threshold's value and level are never recorded. The footer tables (psm.md:258-268, serialization.md:283-296) have no filtering key.
- **Spec:** producer-guide.md:11: 'Keep rows above your FDR threshold only if you say so with pass_threshold'. psm.md:66: null when the producer does not say. psm.md:291: pg_global_qvalue 'used for filtering'.
- **Correction:** Not only PD: mzIdentML fills pass_threshold from passThreshold, and writes false, not null, when the attribute is absent. provenance.md:65 records the converter's own qvalue_threshold ('none' = not filtered), but not the upstream filtering state. The rest holds.
- **Reader impact:** Readers cannot tell whether a PSM file is unfiltered, PSM-filtered or protein-filtered, so re-estimating FDR or counting identifications can be invalid.
- **Suggested:** Add a footer key (filter level and threshold, or none). pass_threshold = the row passes that declared threshold. Reword producer-guide.md:11.
- **Evidence:** `docs/spec/producer-guide.md:11`; `docs/spec/psm.md:66,258-268,291`; `docs/spec/serialization.md:283-296`; `docs/spec/converter-coverage.md:323`

### SPEC-17: Rank policy (rank-1 only vs all ranks) is not specified {#spec-17}

*psm* · unspecified · partly confirmed (see correction) · checklist P1, P6

- **Now:** PD keeps rank>1 hits with a 'rank' score and offers --rank1-only (converter-coverage.md:323). Per #378, 20,001 of 109,700 rows are rank 2-7. Other producers do not say. The default composite has no rank (psm.yaml:4).
- **Spec:** Only 'rank' as a recommended score name (scores.md:58; psm.md:291). The #229 decision 'if a producer retains multiple ranks, rank joins the PK' is not in the spec.
- **Correction:** mzIdentML also keeps every rank by default, with a 'rank' score and rank1_only. provenance.md:65 lists rank1_only among the recorded conversion options, so a converter's rank choice is on record in provenance, not in the footer. The spec has no default policy and no PK rule. The gist holds.
- **Reader impact:** PSM and peptide counts include alternative hypotheses in some files, and rank-2 rows look like chimeric or multi-engine rows.
- **Suggested:** Default to rank-1 only. When rank>1 rows exist, require a typed rank column in the composite and declare the policy in the footer.
- **Evidence:** `docs/spec/scores.md:58`; `docs/spec/psm.md:291`; `qpx/core/data/schemas/psm.yaml:4`; `docs/spec/converter-coverage.md:323`

### SPEC-18: Chimeric/co-isolated PSMs have no marker in the PSM view {#spec-18}

*psm* · unspecified · partly confirmed (see correction) · checklist P2

- **Now:** Only the feature view has identification_type 'secondary' (MaxQuant MULTI-SECPEP: 5,109 rows in maxquant_full, run). Chimeric ProForma ('/2+') is only passed through (scope.md:41). The PSM view has nothing.
- **Spec:** psm.md says nothing about chimeric spectra. feature.md:145 defines 'secondary' for features only.
- **Correction:** MaxQuant PSMs keep the msms.txt Type, including MULTI-SECPEP, as an untyped cv_param 'type'. In maxquant_simple, 30 PSMs carry it (ran). It is undocumented and untyped, and the spec says nothing about chimeric PSMs.
- **Reader impact:** Two PSMs on one scan cannot be told apart as chimeric co-IDs, rank alternatives or engine duplicates.
- **Suggested:** Add a PSM identification_type (primary/secondary) or a chimeric index that mirrors feature, alongside rank.
- **Evidence:** `docs/spec/feature.md:145`; `docs/spec/scope.md:41`; `docs/spec/converter-coverage.md:199`; `audit-run:spec/mq secondary rows`

### SPEC-19: Composite uniqueness rule contradicts producer-supplied duplicate-key PSM ids {#spec-19}

*psm* · contradictory · confirmed · checklist P6, X1

- **Now:** mzIdentML supplies distinct psm_ids for several hits with equal (run,scan,peptidoform,charge) and keeps the schema-default composite in the footer (identity.md:40). PD collapses multi-engine hits to the lowest rank, then the lowest q (converter-coverage.md:323).
- **Spec:** serialization.md:304-305 and versioning.md:86-87: 'rows with equal composite values carry the same id, also when the producer supplies its own ids'. psm.md:12-15 allows producer ids or an extended composite. #229 chose one merged row per spectrum and identification.
- **Reader impact:** Duplicate natural keys appear with some producers and are merged by others. The normative rule outlaws a shipped converter.
- **Suggested:** One row per (run,scan,peptidoform,charge) with engines in additional_scores (consensus_support). Otherwise require an extended composite (engine, rank) declared in the footer.
- **Evidence:** `docs/spec/serialization.md:304-305`; `docs/spec/versioning.md:86-87`; `docs/spec/psm.md:12-15`; `docs/spec/identity.md:40`; `docs/spec/converter-coverage.md:323`

### SPEC-20: observed_mz/calculated_mz/charge meaning differs by producer and is undefined {#spec-20}

*psm* · ambiguous · partly confirmed (see correction) · checklist P7

- **Now:** MaxQuant observed_mz = calculated_mz + m/z*ppm, i.e. recalibrated and monoisotopic (tool-mappings.md:26,47). FragPipe and PD keep isotope-error offsets of +1 to +3 Da (converter-coverage.md:236-238,319). FragPipe calculated_mz includes unassigned shifts (:234). Feature observed_mz comes from the PSM or from MS1 depending on the tool.
- **Spec:** psm.yaml:51-55 'Experimental m/z'. psm.md:43 'charge state of the precursor ion'. Calibration, monoisotopic correction and assigned vs instrument charge are not addressed.
- **Correction:** FragPipe calculated_mz includes the unassigned shift because the peptidoform carries it. This follows the canonical rule (peptidoform.md:47-48) and is not a divergence. The observed_mz semantics (MaxQuant recalibrated, FragPipe/PD with isotope offsets) and the charge semantics remain undefined.
- **Reader impact:** observed_mz and mass_error_ppm are not comparable across producers, e.g. for mass-accuracy QC or ML training.
- **Suggested:** observed_mz = instrument-selected precursor m/z, uncalibrated. Add an optional calibrated or monoisotopic m/z. charge = search-assigned.
- **Evidence:** `qpx/core/data/schemas/psm.yaml:51-55`; `docs/spec/tool-mappings.md:26,47`; `docs/spec/converter-coverage.md:234-238,319`

### SPEC-M5: Coverage says MaxQuant psm.feature_id is not implemented, but the converter fills it {#spec-m5}

*psm* · diverges from spec · found by the verifier · checklist P10

- **Now:** maxquant/psm_adapter.py:232 maps msms.txt Evidence ID to feature_id. Ran maxquant_simple: 95/1000 PSMs are linked and all join to feature rows. The link's meaning (one evidence row) is undocumented.
- **Spec:** converter-coverage.md:285: 'Not implemented: the Evidence ID of msms.txt is not read'.
- **Reader impact:** Readers think MaxQuant has no psm-to-feature link, and the link's meaning (the evidence row, including secondary rows) is undocumented.
- **Suggested:** Document per producer what psm.feature_id links (MaxQuant: the evidence row of the PSM) and fix the coverage row.
- **Evidence:** `qpx/converters/maxquant/psm_adapter.py:232`; `docs/spec/converter-coverage.md:285`; `audit-run:verify-spec/mqs/out (95/1000 linked)`

## Validator and read API

### VALIDATOR-02: Spec-sanctioned provided psm_ids always fail --conformance (identity_not_checked) {#validator-02}

*psm* · contradictory · confirmed · checklist P6, D8

- **Now:** stale_identity_issues emits an 'identity_not_checked' warning for every file whose identity_source is not 'derived' (footer_validation.py:112-114). --conformance fails on any warning. In practice, qpxc convert mzidentml output (PXD054720) exits 1 under --conformance and has no other issue.
- **Spec:** identity.md:40 says mzidentml supplies its own psm_ids by design. producer-guide.md:56 allows a producer to supply its own ids. producer-guide.md:157 describes --conformance as the CI gate.
- **Reader impact:** Producers who follow the spec cannot pass the gate, so the gate loses meaning for psm files with provided ids.
- **Suggested:** Make identity_not_checked an informational note that does not fail conformance, or require a re-derivable composite even for provided ids.
- **Evidence:** `qpx/core/data/footer_validation.py:112-114`; `run: qpxc validate --dataset-path audit-run:mzid_cdap/pxd --conformance -> exit 1, sole issue identity_not_checked`; `docs/spec/identity.md:40`

### VALIDATOR-03: OpenMS converter writes a non-snake_case 'q-value' score and fails conformance {#validator-03}

*psm* · diverges from spec · confirmed · checklist P5, D4

- **Now:** openms_consensus/psm_adapter.py:224 writes score_name 'q-value'. The validator only warns about this (score_name_format, validation.py:143-147), and the warning fails --conformance on all three OpenMS fixtures. The writer still fills psm.qvalue from it after normalising the name (derived.py:47). The other shipped outputs checked: maxquant/simple fails on a dangling run_file_name, maxquant/full, cdap and fragpipe pass.
- **Spec:** scores.md:24 says all score names are snake_case.
- **Reader impact:** A reader that filters additional_scores by the name 'q_value' misses OpenMS rows. Shipped converters do not meet the gate they recommend to producers.
- **Suggested:** Normalise score names on write, and make --conformance part of each converter's tests.
- **Evidence:** `qpx/converters/openms_consensus/psm_adapter.py:224`; `qpx/core/data/validation.py:143-147`; `run: --conformance on audit-run:openms/{ProteomicsLFQ_1,IsobaricWorkflow_1,MS1LabeledWorkflow_6}_out -> exit 1`

### VALIDATOR-04: Several hits of one peptidoform on one spectrum fail strict validation (duplicate_composite) {#validator-04}

*psm* · contradictory · confirmed · checklist P1, P6

- **Now:** duplicate_composite_issues groups rows by the footer composite and reports an error under strict when one composite value carries several ids (footer_validation.py:157-172). I wrote two PSMs with provided ids on the same spectrum, peptidoform and charge, with the default composite left in the footer (as mzid does). Result: an ERROR on duplicate_composite.
- **Spec:** psm.md:12-15 says several hits of one peptidoform on one spectrum need a producer-supplied psm_id. identity.md:40 says the mzid footer keeps the default composite. serialization.md:304 says rows with equal composite values carry the same id, provided ids included.
- **Reader impact:** It is unclear whether two rows with the same (run, scan, peptidoform, charge) are legal. A real mzid with such hits is rejected.
- **Suggested:** Decide on one rule. Either add rank or an engine/hit discriminator to the PSM key, or forbid duplicate hits and require producers to merge them.
- **Evidence:** `qpx/core/data/footer_validation.py:157-172`; `experiment audit-run:validator/scripts/mutate2.py case x_psm_provided_ids_two_hits_same_composite -> ERROR duplicate_composite`; `docs/spec/psm.md:12-15`; `docs/spec/identity.md:40`; `docs/spec/serialization.md:304`

### VALIDATOR-05: No rule or flag covers rank, alternative hits or chimeric PSMs {#validator-05}

*psm* · unspecified · partly confirmed (see correction) · checklist P1, P2

- **Now:** Nothing checks how many rows a spectrum has or what rank they hold. rank exists only as a recommended score name. In the outputs: mzid PXD054720 keeps ranks 2-12 (3,679 rows, 1,773 of them with pass_threshold=true). FragPipe writes 3 rows for scan [100] (two peptidoforms at charge 2, one at charge 3) with no rank. PSM.targets_only() and confident() do not look at rank.
- **Spec:** psm.md:291 recommends a 'rank' score name. producer-guide.md:11 says 'one row per match'. Rank, chimeric hits and alternative hits are not otherwise distinguished (GitHub #229 and #378 are open).
- **Correction:** Of the 3,679 rank>=2 rows, 1,810 (not 1,773) have pass_threshold=true. Chimeric hits are covered at feature level only: identification_type 'secondary' (feature.md:145). PD documents its own rank policy (converter-coverage.md:323). There is still no PSM-level rank column or chimeric/alternative flag, and no validator rule for either.
- **Reader impact:** Counting PSMs or spectra double-counts alternative hits. Rank-2 alternatives look the same as chimeric co-identifications.
- **Suggested:** Add a typed psm.rank column and a chimeric/secondary flag. Declare the rank policy in the footer and validate that rank is present when there are several rows per spectrum.
- **Evidence:** `docs/spec/psm.md:291`; `docs/spec/scores.md:58`; `query audit-run:mzid_cdap/pxd/PXD054720.psm.parquet: rank>=2 rows=3679`; `query audit-run:fragpipe/out/fp.psm.parquet: runA scan [100] -> 3 rows`

### VALIDATOR-06: FDR filtering state is neither validated nor usable from the read API {#validator-06}

*psm* · unspecified · confirmed · checklist P3, D3, D8

- **Now:** pass_threshold is optional and nothing in the validator or read side uses it (grep finds it only in converters). PSM.confident() keeps only rows with qvalue<=x (psm.py:23-25). It returns 0 rows for mzid (13,386 PSMs, no qvalue, 1,869 with pass_threshold=false), and for MaxQuant and FragPipe PSMs. No footer key records the threshold; it appears only in the provenance parameter qvalue_threshold.
- **Spec:** producer-guide.md:11 says to keep sub-threshold rows only if pass_threshold says so. provenance.md:65 records qvalue_threshold. Nothing makes either one machine-checked.
- **Reader impact:** A reader cannot tell whether the psm view is filtered, or at which level. confident() silently returns an empty set when the producer gave no q-value.
- **Suggested:** Require a footer or dataset declaration of the filter level and threshold. Make confident() fall back to pass_threshold, or raise when there is no q-value instead of returning nothing.
- **Evidence:** `qpx/core/data/psm.py:23-25`; `run: ds.psm.confident().count() on audit-run:mzid_cdap/pxd -> 0 of 13386`; `docs/spec/producer-guide.md:11`; `docs/spec/provenance.md:65`

### VALIDATOR-13: Empty PSM scans pass; the conformance fixture fails only because the id goes stale {#validator-13}

*psm* · unspecified · confirmed · checklist P7

- **Now:** The scan_format check counts only non-empty arrays (schema.py:146-157). It runs only for a single-file PSM view (data/base.py:175). The fixture psm_empty_scan (test_conformance.py:94) fails because the derived psm_id becomes stale. With the id re-derived, scan=[] passes --conformance.
- **Spec:** scan.md:83-88: non-empty arrays are checked, and unions of shards are not. producer-guide.md:23: a DDA PSM's scan is the MS2 scan number. Nothing says whether an empty scan is allowed.
- **Reader impact:** PSMs without a spectrum reference cannot be joined to mz and are indistinguishable from real ones. The test suite gives false confidence.
- **Suggested:** Forbid empty scan on PSMs, check scan_format per shard, and re-derive ids inside the conformance fixtures so each one tests the rule it names.
- **Evidence:** `qpx/core/data/schema.py:146-157`; `qpx/core/data/base.py:175`; `tests/conformance/test_conformance.py:94`; `experiment mutate.py p_empty_scan_rederived -> exit 0; mutate4.py z_psm_empty_scan_not_rederived -> only stale_identity`

### VALIDATOR-16: psm.label is not resolved to a sample or checked against the label vocabulary {#validator-16}

*psm* · unspecified · confirmed · checklist P9, D7

- **Now:** psm.label is checked only against feature.identified_labels, and only when both that column and psm.feature_id are populated (label_validation.py:259-266,441-451). Label vocabulary checks cover feature and pg only (data/label_validation.py:11-14). psm.label='TMT126' on an LFQ dataset passes --conformance.
- **Spec:** psm.md:74: label applies to MS1-labelled experiments, and is null otherwise. data-model.md:239-246 requires labels to resolve exactly to a sample.
- **Reader impact:** A PSM can name a channel that does not exist in its run.
- **Suggested:** Validate psm.label like feature labels: canonical vocabulary, exactly one run.samples entry, and null unless the run is MS1-labelled.
- **Evidence:** `qpx/core/label_validation.py:259-266,441-451`; `qpx/core/data/label_validation.py:11-14`; `experiment mutate.py p_label_on_lfq; mutate2.py x_psm_label_no_feature -> exit 0`

### VALIDATOR-32: Lengths of the parallel PSM spectral arrays are not checked {#validator-32}

*psm* · diverges from spec · confirmed · checklist P7

- **Now:** The length check runs only for XIC (validation.py:283-288). A PSM with mz_array of length 2 and intensity_array of length 1 passes --conformance.
- **Spec:** psm.md:290: mz_array, intensity_array, charge_array, ion_type_array and ion_mobility_array are parallel arrays of the same length.
- **Reader impact:** Fragment annotations can be misaligned with peaks without any warning.
- **Suggested:** Reuse the XIC length check for the PSM spectral arrays.
- **Evidence:** `qpx/core/data/validation.py:233-288`; `experiment mutate.py p_spectral_arrays_mismatch -> exit 0`; `docs/spec/psm.md:290`

### VALIDATOR-M3: ModificationView counts rank-2+ and below-threshold PSMs {#validator-m3}

*psm* · ambiguous · found by the verifier · checklist P1, P3

- **Now:** modification_view.frequency() counts every non-decoy PSM row carrying a modification, ignoring rank, pass_threshold and qvalue (views/api.py:139-151). On mzid PXD054720: Deamidated 4,323 PSMs vs 926 rank-1 passing rows; Oxidation 3,632 vs 2,375.
- **Spec:** views.md:130: 'Modification frequency across PSMs'. Nothing says which PSMs count.
- **Reader impact:** Modification statistics from a file that keeps all ranks are inflated several-fold, and are not comparable with producers that keep rank 1 only.
- **Suggested:** Count only rank-1 rows that pass the declared threshold (pass_threshold, or qvalue when present), or expose the filter and state the default.
- **Evidence:** `qpx/views/api.py:139-151`; `run: ds.modification_view.frequency() on audit-run:mzid_cdap/pxd vs duckdb filtered count`; `docs/spec/views.md:130`

## mokume (consumer)

### MOKUME-25: Spectral count counts every feature-linked PSM row, with no rank or threshold check {#mokume-25}

*psm* · unspecified · confirmed · checklist P1, P3, P10, D3

- **Now:** spectral_count.rs counts distinct psm_id values that have a feature_id, excluding decoys and short peptides (158-195, 210-233). The PSM reader reads no rank, pass_threshold or qvalue (psm.rs:14-36). Sample lookup passes no label and so fails for multiplexed runs (spectral_count.rs:197-208; sdrf.rs:123-145). feature_id is filled only by some converters (MaxQuant psm_adapter.py:232; OpenMS consensus).
- **Spec:** producer-guide.md:11 and psm.md:66: sub-threshold PSMs may be kept and flagged by pass_threshold. psm.md:291: rank is only an additional score. psm.md:75 and 146-160: feature_id is optional and partial.
- **Reader impact:** Counts depend on each producer's rank and threshold policy. Producers without feature_id cannot be counted at all.
- **Suggested:** PSMs linked by feature_id must be rank 1 and passing, or pass_threshold becomes mandatory. Declare the rank policy in the footer.
- **Evidence:** `mokume/rust/crates/mokume-pipeline/src/spectral_count.rs:158-233`; `mokume/rust/crates/mokume-io/src/psm.rs:14-36`; `qpx/converters/maxquant/psm_adapter.py:232`; `docs/spec/producer-guide.md:11`; `docs/spec/psm.md:66,75,146-160,291`

## External tools and standards

### EXTERNAL-01: PSM rank: mzIdentML requires rank with ties; QPX keeps rank only as optional score {#external-01}

*psm* · unspecified · confirmed · checklist P1, P6

- **Now:** mzIdentML 1.2 SII rank is required: '1 is the top rank. If multiple identifications have the same top score, they should all be assigned rank=1'. The QPX mzIdentML converter keeps every rank by default. In the PXD054720 run: ranks 1–12, and 1,758 spectra have >1 row. Rank is written only as additional_scores 'rank' (psm_records.py:63,140). The MaxQuant PSMs carry no rank.
- **Spec:** The composite [run_file_name,scan,peptidoform,charge] names one peptidoform per spectrum (psm.md:9-15). 'rank' is only a recommended additional score (psm.md:291; scores.md:58). No rule on whether the view holds rank-1 only, all ranks, or ties.
- **Reader impact:** A reader cannot tell whether a psm file is rank-1 only. Counting rows overcounts identifications. Ties and rank>1 rows are indistinguishable unless the nested score is unpacked.
- **Suggested:** Make rank a typed nullable int column with mzIdentML semantics (ties share a rank). Declare a footer rank policy (rank1_only/all_ranks/top_n). Default contract: rank-1 rows are the identifications.
- **Evidence:** `https://peptideatlas.org/PSI/schemas/mzIdentML/1.2/mzIdentML1.2.0_combined.html (SpectrumIdentificationItem@rank)`; `qpx/converters/mzidentml/psm_records.py:63,140`; `run: audit-run:external/mzid/out/mzidentml.psm.parquet (ranks 1-12; 1758 multi-row spectra)`; `docs/spec/psm.md:9-15,291`; `docs/spec/scores.md:58`

### EXTERNAL-02: pass_threshold copies mzIdentML passThreshold, whose meaning is engine- and threshold-specific {#external-02}

*psm* · ambiguous · confirmed · checklist D3, P3

- **Now:** mzIdentML: passThreshold is required, and 'If no such threshold has been set, value of true should be given for all results'. The fixture's <Threshold> is Mascot:SigThreshold 0.01 homology, not an FDR; its rank-2 hits have passThreshold true. QPX copies the attribute (mzidentml/psm_adapter.py:304): 871 rank-2 rows have pass_threshold=true, and no q-value. PD: pass_threshold = Confidence=='High' (proteome_discoverer/psm_adapter.py:68), whose cutoff is user-set (default 1%).
- **Spec:** pass_threshold: 'Whether the search engine or validator marks the PSM as passing its threshold. null when it does not say' (psm.md:66; psm.yaml:133). It records neither the threshold statistic nor its level or value. producer-guide.md:11 ties it to an 'FDR threshold'.
- **Reader impact:** Readers will treat pass_threshold=true as 'passes 1% FDR'. It may mean Mascot homology p<0.01, an unthresholded file (mzIdentML true-by-default), or PD High at any user cutoff.
- **Suggested:** pass_threshold is meaningful only with a declared threshold: statistic (PSI-MS term), level, value. Map MS:1001494 'no threshold' to null, not true. Record the mzIdentML <Threshold> cvParams in provenance.
- **Evidence:** `https://peptideatlas.org/PSI/schemas/mzIdentML/1.2/mzIdentML1.2.0_combined.html (passThreshold; Threshold 1..1)`; `tests/examples/mzidentml/PXD054720/F001234.mzid.gz <Threshold> Mascot:SigThreshold homology`; `qpx/converters/mzidentml/psm_adapter.py:304`; `qpx/converters/proteome_discoverer/psm_adapter.py:68`; `https://www.cores.emory.edu/eipc/_includes/documents/definitions_pd.pdf (High/Medium = 1%/5%)`; `docs/spec/psm.md:66`; `docs/spec/producer-guide.md:11`

### EXTERNAL-04: Chimeric spectra have a PSI notation ('+'); QPX psm has no chimeric marker {#external-04}

*psm* · unspecified · confirmed · checklist P2

- **Now:** ProForma 2.0 §7.2: with chimeric spectra 'multiple peptidoform sequences MUST be separated by the plus sign', e.g. EMEVEESPEK/2+ELVISLIVER/3. This is part of the 'level 2-ProForma + mass spectrum compliant' extension (Appendix I). USI 1.0 §3.6.5 allows several interpretations joined by '+'. In PSI-MS (OLS search) no CV term for chimeric or co-fragmenting identification was found.
- **Spec:** Chimeric ProForma is passed through unchanged (scope.md:41). The psm view is one peptidoform per row (psm.md:12-15) and has no field linking co-identified peptides of one spectrum. feature.identification_type 'secondary' exists only on features (feature.md:145).
- **Reader impact:** Two PSM rows on one scan may be co-fragmenting peptides (both true) or rank-1/rank-2 alternatives (one is wrong). Readers cannot tell them apart and will deduplicate or count wrongly.
- **Suggested:** Forbid '+' in psm.peptidoform: one row per co-identified peptidoform. Add psm.spectrum_role (primary/co_fragmenting/alternative) or chimeric_rank, distinct from search rank. Export '+' only when deriving a USI.
- **Evidence:** `https://arxiv.org/pdf/2109.11352 (ProForma 2.0 §7.2, Appendix I)`; `https://raw.githubusercontent.com/HUPO-PSI/usi/master/Specification/USI_FullSpecification_1.0.0_final.pdf §3.6.5`; `https://www.ebi.ac.uk/ols4/api/search?q=chimeric&ontology=ms (no hits)`; `docs/spec/scope.md:41`; `docs/spec/psm.md:12-15`; `docs/spec/feature.md:145`

### EXTERNAL-05: MaxQuant second peptides are PSM rows marked only by a cv_param {#external-05}

*psm* · inconsistent across producers · confirmed · checklist P2, D1

- **Now:** MaxQuant Type MULTI-SECPEP = 'MS1 labeling cluster identified by MS/MS as second peptide' (coxdocs). The msms.txt fixture has 30 MULTI-SECPEP rows. In the QPX run, scan 14256 has two rows (MULTI-SECPEP and MULTI-MSMS). Their only marker is cv_params {cv_name:'type'} (maxquant/psm_adapter.py:348-353), with no rank. pass_threshold is null on all 1,000 rows although msms.txt is FDR-filtered.
- **Spec:** The psm view has no identification_type or role field. The feature view has identification_type 'secondary' for MULTI-SECPEP (feature.md:145; converter-coverage.md:199). The spec defines no cv_params 'type' vocabulary.
- **Reader impact:** The chimeric status of a MaxQuant PSM is visible only through a MaxQuant-specific cv_param. Other producers' PSMs cannot be compared, and generic readers miss it.
- **Suggested:** Use the same typed role column on psm and feature (see EXTERNAL-04/32). The MaxQuant converter maps MULTI-SECPEP to it. Set pass_threshold=true for msms.txt rows (tool-filtered output) once the FDR declaration exists.
- **Evidence:** `https://cox-labs.github.io/coxdocs/output_tables.html (evidence Type)`; `tests/examples/maxquant/maxquant_simple/msms.txt (Type counts: MULTI-MSMS 935, MSMS 35, MULTI-SECPEP 30)`; `qpx/converters/maxquant/psm_adapter.py:348-353`; `run: audit-run:external/mq/out/maxquant.psm.parquet (scan [14256] two rows)`

### EXTERNAL-06: Sage reports chimeric co-identifications in the same rank column as alternatives {#external-06}

*psm* · ambiguous · confirmed · checklist P1, P2

- **Now:** Sage chimera=true takes the best candidate, subtracts its matched peaks and searches the spectrum again until report_psms candidates are found. results.sage.tsv 'rank: If report_psms > 1, then the best match will have rank = 1'. Chimeric peptides and plain alternatives share one rank column, and the file does not say which mode produced them. Sage also gives spectrum_q, peptide_q, protein_q and posterior_error.
- **Spec:** rank = 'Rank of the peptide in the search engine results (1 = best)' (scores.md:58). No converter for Sage exists (converter-coverage.md). The spec does not say that rank>1 can be a true co-identification.
- **Reader impact:** A rank filter (rank==1) silently drops valid chimeric identifications from Sage, while keeping rank>1 keeps wrong alternatives from other engines.
- **Suggested:** The contract separates 'search rank' (alternative interpretations) from chimeric co-identification. Producers record the search mode (chimera on/off, report_psms) in provenance. A rank>1 row from a chimeric search gets spectrum_role=co_fragmenting.
- **Evidence:** `https://sage-docs.vercel.app/docs/configuration/additional (chimera, report_psms)`; `https://sage-docs.vercel.app/docs/results/search (rank, spectrum_q, peptide_q, protein_q)`; `https://github.com/lazear/sage (picked-peptide & picked-protein FDR)`; `docs/spec/scores.md:58`

### EXTERNAL-07: CHIMERYS/DDA+ report several PSMs per spectrum with tool-reassigned precursors {#external-07}

*psm* · ambiguous · confirmed · checklist P2, P7, P3

- **Now:** CHIMERYS: mokapot PSM-level FDR 'specifically allowing for multiple PSMs per spectrum'; 'We specifically prevent mokapot from only retaining the top-scoring PSM per spectrum'. It has no primary/secondary distinction, and I/L isomers are 'added back with identical scores and flagged as ambiguous'. MSFragger-DDA+ searches the full isolation window and reports 'the observed precursor m/z and charge that may be different from that listed for the DDA spectrum'.
- **Spec:** observed_mz: 'Experimental observed peptide m/z' (psm.md:47). charge: 'Charge state of the precursor ion' (psm.md:43). Neither says whether the value is the instrument-selected precursor or the tool-assigned co-isolated precursor.
- **Reader impact:** For chimeric PSMs, observed_mz and charge may describe a different ion from the instrument's selected precursor. Joins on mz/charge to mzML precursor metadata or to MS1 features then fail or mismatch.
- **Suggested:** Define observed_mz and charge as those of the ion the identification claims (tool-assigned). Allow an optional selected_precursor_mz/charge (instrument) column. State that several PSMs per scan can all pass FDR.
- **Evidence:** `https://pmc.ncbi.nlm.nih.gov/articles/PMC12074992/ (CHIMERYS methods)`; `https://pmc.ncbi.nlm.nih.gov/articles/PMC11978857/ (MSFragger-DDA+)`; `docs/spec/psm.md:43,47`

### EXTERNAL-08: FDR competition unit and rank>1 q-values: PSI and tool rules not in QPX {#external-08}

*psm* · unspecified · confirmed · checklist P3, P5, P6, D4

- **Now:** Percolator TDC: 'eliminate the lower scoring target or decoy PSM(s) of each scan+expMass combination'. -U keeps all PSMs and drops peptide-level output, and peptide level keeps 'only the best scoring PSM for each unique peptide'. OpenMS PercolatorAdapter default keeps all hits per spectrum (-best_per_spectrum_only off). MS2Rescore: for ranks>1 'q-values/PEPs are computed treating each rank as an independent row … not statistically rigorous'.
- **Spec:** psm.qvalue = 'PSM, as computed by the source' (scores.md:42). The spec does not say that q-values and PEPs of rank>1, non-winning or chimeric rows may not be FDR-controlled, or what unit competed (spectrum vs spectrum+precursor mass).
- **Reader impact:** Filtering psm.qvalue<=0.01 over all rows of a multi-rank file admits rank>1 rows whose q-values were never estimated by competition, inflating identifications.
- **Suggested:** Contract: q-value/PEP are FDR-valid only for rows that won the declared competition. Footer declares competition_unit (spectrum / spectrum+precursor_mass / none). Rows outside it carry rank>1 or a null qvalue.
- **Evidence:** `https://github.com/percolator/percolator/wiki/Command-line-options`; `https://gist.github.com/dub2s/5d3ad244cadaf84c4645780dca23ceeb (Percolator 3.02 log)`; `https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_PercolatorAdapter.html`; `https://ms2rescore.readthedocs.io/en/v4.0.1/userguide/output-files/`; `docs/spec/scores.md:42`

### EXTERNAL-25: mzIdentML 1.2 peptide-level scores on PSMs have no QPX psm column {#external-25}

*psm* · unspecified · confirmed · checklist D4, P5

- **Now:** mzIdentML 1.2 added a 'peptide-level scoring' protocol term (MS:1002490). Each SII then carries the triplet MS:1002520 peptide group ID, MS:1002500 peptide passes threshold and a peptide-level score (child of MS:1002358). Percolator and Sage also give each PSM a peptide-level q-value (peptide_q).
- **Spec:** psm.qvalue excludes names containing 'peptide' (scores.md:42), and the psm view has no peptide_qvalue column (psm.md:35-57). Peptide-level q-values are typed only on features (feature.md:57).
- **Reader impact:** For a psm-only producer (mzIdentML, Sage), peptide-level confidence sits under producer-specific names in additional_scores. Readers cannot filter PSMs at a peptide-level FDR portably.
- **Suggested:** Add an optional psm.peptide_qvalue (MS:1001868 semantics) and peptide_pass_threshold, or a documented score name. The mzIdentML converter reads the MS:1002520/MS:1002500 triplet.
- **Evidence:** `https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5500760/ (mzIdentML 1.2)`; `https://sourceforge.net/p/psidev/mailman/message/35083384 (peptide-level triplet rule)`; `https://sage-docs.vercel.app/docs/results/search`; `docs/spec/scores.md:42`; `docs/spec/feature.md:57`

### EXTERNAL-26: USI index types match QPX scan_format, but scan_format is optional {#external-26}

*psm* · specified, worth keeping · confirmed · checklist P7, D8

- **Now:** USI 1.0: indexType 'scan' preferred (the Thermo scan number, not the scan index), 'index' = Nth spectrum of an MGF counting from 0, 'nativeId' = comma-separated integers in PSI-MS nativeId order (Thermo only when controllerType≠0), and 'trace'. The msRun is 'typically the root name' with or without extension. The QPX mzIdentML converter leaves scan_format undeclared on mixed schemes (mzidentml/psm_records.py:40-43).
- **Spec:** scan.md:3,20,51-61 follows USI (scan/index/nativeId; Thermo controllerType rule). scan_format is declared 'when it knows the format' (scan.md:51; psm.md:264).
- **Reader impact:** When scan_format is absent, a one-element scan cannot be told from an index (scan.md:51). USIs cannot be built and joins to mzML may pick the wrong spectrum.
- **Suggested:** Keep USI alignment. Make scan_format mandatory for psm and feature (per row if mixed), consider 'trace', and document USI construction from run_file_name+scan_format+scan+peptidoform/charge.
- **Evidence:** `https://raw.githubusercontent.com/HUPO-PSI/usi/master/Specification/USI_FullSpecification_1.0.0_final.pdf §3.6.2-3.6.4`; `https://pmc.ncbi.nlm.nih.gov/articles/PMC8405201`; `docs/spec/scan.md:3,20,51-61`; `qpx/converters/mzidentml/psm_records.py:40-43`

### EXTERNAL-27: Several searches of one spectrum: USI adds a provenance id, QPX collapses them {#external-27}

*psm* · ambiguous · partly confirmed (see correction) · checklist P6

- **Now:** USI §3.6.6: 'multiple searches of the same spectrum could yield the same interpretation but at different false discovery rates'. An optional PSM provenance identifier suffix disambiguates them. PD consensus results carry several engines (Sequest, Comet). QPX PD keeps 'the lowest rank then the lowest q-value' per identity, with 'Identifying Node' (converter-coverage.md:323).
- **Spec:** The composite names one peptidoform per spectrum. A file with several hits of one identity needs a producer psm_id or an extended composite (psm.md:12-15). It does not choose between merging engines into one row and keeping one row per engine.
- **Correction:** The USI PSM provenance identifier is a repository-assigned id, not a field producers or QPX would mint. mzTab 1.0 §5.12.1 (one row carrying several engines' scores) is the better standard anchor for the suggested decision.
- **Reader impact:** Multi-engine scores may be merged, dropped (PD keeps one) or duplicated under producer ids, depending on the converter. consensus_support (scores.md:75) is not comparable.
- **Suggested:** Contract: one row per (spectrum, peptidoform, charge); each engine's scores go to additional_scores with an engine prefix; engines and selection rule are declared in provenance. Optionally carry a search_id for USI provenance.
- **Evidence:** `https://raw.githubusercontent.com/HUPO-PSI/usi/master/Specification/USI_FullSpecification_1.0.0_final.pdf §3.6.6`; `docs/spec/converter-coverage.md:323`; `docs/spec/psm.md:12-15`; `docs/spec/scores.md:75`

### EXTERNAL-29: Isobaric sequence ambiguity (I/L, ties) has standard encodings; QPX has no rule {#external-29}

*psm* · unspecified · partly confirmed (see correction) · checklist P9, P1

- **Now:** mzTab 1.0 §5.12: ambiguous residues 'SHOULD' use B, Z, J (Leu/Ile) and X. CHIMERYS keeps the best of a group of isobaric PSMs and adds I/L isomers back 'with identical scores … flagged as ambiguous'. mzIdentML: equal top scores 'should all be assigned rank=1'. ProForma Base supports ambiguity in the amino acid sequence.
- **Spec:** Nothing on I/L or tied identifications in psm or peptidoform (peptidoform.md, psm.md). Tied peptidoforms become separate rows under the composite (psm.md:9-15). Only the DIA-NN path notes J support (converter-coverage.md:123).
- **Correction:** converter-coverage.md:123 concerns the protein molecular-weight readers (OpenMS consensusXML and others accept J in protein sequences). It is not DIA-NN peptidoform handling. The mzTab section is §5.12.4.
- **Reader impact:** An I/L-ambiguous spectrum may appear as 1 row (best only), 2 rank-1 rows (two 'identifications') or a J-sequence, depending on producer. Peptide and protein counts and pepmap differ.
- **Suggested:** Choose one rule: keep the tool's sequences, emit tied rows with the same rank plus an ambiguity flag, and never collapse to J in peptidoform. Pepmap lists the proteins of every tied sequence.
- **Evidence:** https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/specification_document-releases/1_0-Proteomics-Release/mzTab_format_specification.pdf §5.12; `https://pmc.ncbi.nlm.nih.gov/articles/PMC12074992/`; `https://peptideatlas.org/PSI/schemas/mzIdentML/1.2/mzIdentML1.2.0_combined.html`; `docs/spec/converter-coverage.md:123`

### EXTERNAL-31: Decoy representation: standards flag decoys, but whether a file retains them is undeclared {#external-31}

*psm* · unspecified · confirmed · checklist P4, F8, D8

- **Now:** mzTab 1.0 §5.12.5: decoys 'MAY be labeled' with opt_cv_MS:1002217_decoy_peptide (1/0). mzIdentML PeptideEvidence@isDecoy is optional. Sage: label -1/1. MaxQuant Reverse '+' rows appear in evidence/msms ('should be removed'). FragPipe reports are written with decoys filtered out (pg.md:365).
- **Spec:** is_decoy is required, and target+decoy counts as target (psm.md:45; feature.md:60). Producer-guide keeps decoy markers in accessions (producer-guide.md:30-35). No footer states whether decoys were retained or removed by the tool.
- **Reader impact:** A reader re-estimating FDR from is_decoy cannot tell 'no decoys in file' from 'no decoys passed'. Decoy-retaining and decoy-free files look the same.
- **Suggested:** Footer key decoys (retained_all/retained_passing/removed/none_searched) per view. Contract keeps decoys optional but declared. Map the mzTab and mzIdentML decoy flags explicitly.
- **Evidence:** https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/specification_document-releases/1_0-Proteomics-Release/mzTab_format_specification.pdf §5.12.5; `https://peptideatlas.org/PSI/schemas/mzIdentML/1.2/mzIdentML1.2.0_combined.html (PeptideEvidence@isDecoy)`; `https://sage-docs.vercel.app/docs/results/search`; `https://cox-labs.github.io/coxdocs/output_tables.html (Reverse)`; `docs/spec/psm.md:45`; `docs/spec/pg.md:365`

### EXTERNAL-32: One PSM, many proteins: mzTab and mzIdentML list all mappings, as QPX does {#external-32}

*psm* · specified, worth keeping · confirmed · checklist P8

- **Now:** mzTab 1.0 §6.5.2: 'If a PSM can be matched to multiple proteins, the same PSM should be represented on multiple rows with different accessions and the same PSM_ID'. mzIdentML SII: 'If a specific sequence can be assigned to multiple proteins … all possible PeptideEvidence elements should be referenced'. FragPipe Mapped Proteins includes I/L-substitution mappings.
- **Spec:** protein_accessions lists 'all proteins that the peptide maps to'; null when no mapping (psm.md:57). No protein inference in the psm view (psm.md:289; #78). pepmap holds every protein (producer-guide.md:16).
- **Reader impact:** Little risk, but order and I/L-derived mappings (FragPipe) are unspecified, so array equality across producers fails.
- **Suggested:** Keep 'all mapped proteins'. Specify sorted order and whether I/L-equivalent mappings are included. Never store razor/group in psm.protein_accessions.
- **Evidence:** https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/specification_document-releases/1_0-Proteomics-Release/mzTab_format_specification.pdf §6.5.2-6.5.3; `https://peptideatlas.org/PSI/schemas/mzIdentML/1.2/mzIdentML1.2.0_combined.html`; `https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html (Mapped Proteins)`; `docs/spec/psm.md:57,289`

### EXTERNAL-M1: mzIdentML 'final PSM list' flag is ignored; every SpectrumIdentificationList is converted {#external-m1}

*psm* · unspecified · found by the verifier · checklist P3, P6, P1

- **Now:** PSI-MS MS:1002439 'final PSM list' flags the SIL(s) to interpret. Among flagged lists, 'each spectrum must not be referenced from more than one SpectrumIdentificationResult'. QPX streams every SIR of every list (mzidentml/document.py:44-62), and nothing in qpx/ reads MS:1002439. The fixture has one SIL (mzid.gz:102758), so the problem is latent.
- **Spec:** psm.md:12-15 only says that files with several hits per spectrum need a producer psm_id. Nothing distinguishes intermediate from final result lists.
- **Reader impact:** A multi-list mzIdentML (per-engine intermediate lists plus a consensus/final list) mixes unvalidated and final PSMs in one psm file, duplicating spectra and inflating counts.
- **Suggested:** Contract: the psm view holds the final result list. The mzid converter keeps the SILs flagged MS:1002439 when any are present, else all lists with a warning, and records the list ids in provenance.
- **Evidence:** `https://www.ebi.ac.uk/ols4/api/ontologies/ms/terms?obo_id=MS:1002439`; `qpx/converters/mzidentml/document.py:44-62`; `grep -rn 1002439 qpx/ -> no hits`; `tests/examples/mzidentml/PXD054720/F001234.mzid.gz:102758 (single SIL_1)`; `docs/spec/psm.md:12-15`

### EXTERNAL-M2: MaxQuant packs alternative candidates in one row; QPX drops them, while mzIdentML keeps all ranks {#external-m2}

*psm* · inconsistent across producers · found by the verifier · checklist P1, P2

- **Now:** msms.txt has one row per identified MS/MS plus MULTI-SECPEP rows. Alternatives sit only in 'All sequences', 'All modified sequences' and 'All scores' (fixture: 967/1000 rows list 3 candidates). 'Delta score' = 'difference to the second best identified peptide with a different amino acid sequence' (coxdocs). QPX reads Delta score (maxquant/psm_adapter.py:328) but never 'All sequences'. mzIdentML keeps all ranks (psm_records.py:63).
- **Spec:** Rank is only an optional score (scores.md:58; psm.md:291). No footer or rule declares whether alternatives are present.
- **Reader impact:** Equivalent engine output gives rank-1(+second peptide)-only rows for MaxQuant and all-rank rows for mzIdentML. Row counts and the presence of alternatives differ silently by producer.
- **Suggested:** Every converter declares the footer rank policy (see EXTERNAL-01): MaxQuant = rank1_plus_secondary, mzIdentML = all_ranks unless --rank1-only. MaxQuant alternatives are emitted as rank>1 rows only on request.
- **Evidence:** `https://cox-labs.github.io/coxdocs/output_tables.html (msms Delta score, All scores)`; `tests/examples/maxquant/maxquant_simple/msms.txt cols 53-55 (All scores/All sequences/All modified sequences)`; `qpx/converters/maxquant/psm_adapter.py:328`; `grep 'All sequences' qpx/converters/maxquant -> no hits`; `qpx/converters/mzidentml/psm_records.py:63`

### EXTERNAL-M3: observed_mz is recalibrated for MaxQuant, raw or isotope-shifted elsewhere {#external-m3}

*psm* · ambiguous · found by the verifier · checklist P7

- **Now:** coxdocs: 'Mass error [ppm]: Mass error of the recalibrated mass-over-charge value'. QPX MaxQuant observed_mz = theoretical m/z + recalibrated ppm error (maxquant/constants.py:259-271), so it is a recalibrated value. PD observed_mz carries Comet isotope errors of +1 to +3 Da (converter-coverage.md:319). mzIdentML experimentalMassToCharge is copied as is (mzidentml/psm_adapter.py:300). DDA+ may report a non-selected precursor.
- **Spec:** observed_mz = 'Experimental observed peptide mass-to-charge ratio' (psm.md:47; feature.md:62), and mass_error_ppm is derived from it (psm.md:48). The calibration and isotope-error state are not stated.
- **Reader impact:** mass_error_ppm compares recalibrated MaxQuant values with uncalibrated or isotope-offset values from other producers. Ppm filters, calibration QC and joins to mzML precursors behave differently by producer.
- **Suggested:** Define observed_mz as the monoisotopic m/z of the claimed ion. Declare the calibration state (raw/recalibrated) in provenance. Correct isotope-error offsets or flag them (e.g. an isotope_error score).
- **Evidence:** `https://cox-labs.github.io/coxdocs/output_tables.html (msms m/z, Mass error [ppm])`; `qpx/converters/maxquant/constants.py:259-271`; `docs/spec/converter-coverage.md:319`; `qpx/converters/mzidentml/psm_adapter.py:300`; `docs/spec/psm.md:47-48`

## First-principles questions

### CRITIC-01: What is one PSM row: one per spectrum, per hit, or per ranked candidate? {#critic-01}

*psm* · ambiguous · partly confirmed (see correction) · checklist P1, P6

- **Now:** Spec says one row per match, with a PK composite (run_file_name, scan, peptidoform, charge) that allows one peptidoform per spectrum. Nothing says whether rank>1 candidates, or other charges for the same scan, are rows. PD keeps ranks 2-7 (converter-coverage.md:323). mzIdentML keeps all ranks unless --rank1-only (qpx/converters/mzidentml/psm_records.py:63).
- **Spec:** psm.md:3-15; data-model.md:61-63; producer-guide.md:11. The #229 decision 'rank joins the PK' is not in the spec.
- **Correction:** The default composite does allow several peptidoforms per spectrum, one row each. Only a repeated (peptidoform, charge) on the same spectrum needs a supplied id (psm.md:12-15). The gist holds: whether rank>1 rows exist is unspecified. The mzIdentML fixture keeps ranks 1-12 by default (13,386 rows over 9,707 spectra). The 'PD ranks 2-7' figure comes from #378, not converter-coverage.md:323.
- **Reader impact:** PSM/spectrum counts and the choice of 'the' identification per scan silently mix rank-1 hits with alternatives. Spectral counts are inflated for some producers.
- **Suggested:** Grain = one row per (spectrum, candidate). Make rank mandatory (1 = best per spectrum) and add a footer rank_policy (rank1_only|all_ranks).
- **Evidence:** `docs/spec/psm.md:3-15`; `docs/spec/data-model.md:61-63`; `docs/spec/converter-coverage.md:323`; `qpx/converters/mzidentml/psm_records.py:63`

### CRITIC-02: Where is PSM rank stored, and is rank per spectrum or per spectrum+charge? {#critic-02}

*psm* · unspecified · confirmed · checklist P1, D8

- **Now:** Rank exists only as a recommended score name `rank` in additional_scores. There is no typed column and no scope (per spectrum across charges? per engine?), and ties are undefined. No footer says whether lower ranks were dropped.
- **Spec:** psm.md:291; scores.md:58. provenance.md:65 records rank1_only only as a qpxc conversion option.
- **Reader impact:** A rank-1 filter means scanning nested scores, and it fails for producers that omit them. 'Absent because rank-1 only' cannot be told from 'unknown'.
- **Suggested:** Typed nullable int16 `rank`: per spectrum, across charges and engines; ties share a rank. Add a footer rank_policy.
- **Evidence:** `docs/spec/psm.md:291`; `docs/spec/scores.md:58`; `docs/spec/provenance.md:65`

### CRITIC-03: How are chimeric co-identifications told apart from rank-2 alternatives on PSM rows? {#critic-03}

*psm* · unspecified · partly confirmed (see correction) · checklist P2

- **Now:** PSM has no field for this; `secondary` exists only in feature.identification_type. ProForma chimeric strings (`PEPTIDE/2+SEQUENCE/3`) pass through unchanged, so one producer may write one row with a chimeric peptidoform while another writes two rows.
- **Spec:** feature.md:145; scope.md:41; psm.md:12-15
- **Correction:** The MaxQuant PSM rows keep the msms.txt Type, MULTI-SECPEP included (30/1000 rows in maxquant_simple), as the producer-specific cv_param 'type'. There is no typed or spec-defined field, so the gist holds.
- **Reader impact:** Second peptides (MaxQuant MULTI-SECPEP, Sage/MSFragger chimeric) look like low-rank alternatives or duplicate evidence, which corrupts FDR and spectral counts.
- **Suggested:** Add a PSM identification_type (primary|secondary). Forbid chimeric ProForma in peptidoform: one row per peptidoform.
- **Evidence:** `docs/spec/feature.md:145`; `docs/spec/scope.md:41`; `docs/spec/psm.md:12-15`

### CRITIC-04: Is the PSM view FDR-filtered, at which level and threshold, and by whom? {#critic-04}

*psm* · unspecified · confirmed · checklist P3, D8

- **Now:** No column or footer key states the filtering. producer-guide implies failing rows are kept only when pass_threshold flags them. Provenance records only qpxc's own qvalue_threshold ('none' by default), not upstream filtering (e.g. tool-side PSM/protein FDR). The XIC and mz views declare their selection in the footer; psm does not.
- **Spec:** producer-guide.md:11; provenance.md:65; xic.md:21 (xic_selection, xic_qvalue_threshold); mz.md:154
- **Reader impact:** Rows could be all candidates, 1% PSM-FDR survivors or protein-FDR survivors. Re-estimating FDR or comparing identification counts across producers is invalid.
- **Suggested:** Mandatory footer keys per view: filter_level (none|psm|peptide|precursor|protein), filter_threshold, filter_source (tool|qpx).
- **Evidence:** `docs/spec/producer-guide.md:11`; `docs/spec/provenance.md:65`; `docs/spec/xic.md:21`; `docs/spec/mz.md:154`

### CRITIC-05: What exactly does pass_threshold=true certify? {#critic-05}

*psm* · ambiguous · partly confirmed (see correction) · checklist D3

- **Now:** Optional bool: 'search engine or validator marks the PSM as passing its threshold'; null = does not say. The threshold value, level and tool are not recorded. PD uses Confidence=High. In the mzIdentML converter a missing passThreshold attribute becomes false, not null (qpx/converters/mzidentml/psm_adapter.py:304).
- **Spec:** psm.md:66; psm.yaml:133-136; converter-coverage.md:323; versioning.md:75
- **Correction:** passThreshold is a required attribute in the mzIdentML XSD, so 'missing -> false' (mzidentml/psm_adapter.py:304) only affects malformed files. The undefined threshold, level and tool remain real. pass_threshold is also true on 1,810 rank>=2 rows (see CRITIC-M2).
- **Reader impact:** pass_threshold=true can mean 1% PSM FDR in one file and a 5% or engine-specific cut in another, so the flags cannot be pooled.
- **Suggested:** Add a footer pass_threshold_definition (level, value, tool). Write null, never false, when the source is silent.
- **Evidence:** `docs/spec/psm.md:66`; `qpx/core/data/schemas/psm.yaml:133-136`; `docs/spec/converter-coverage.md:323`; `qpx/converters/mzidentml/psm_adapter.py:304`

### CRITIC-06: What level and scope is psm.qvalue (run vs experiment, which engine)? {#critic-06}

*psm* · ambiguous · confirmed · checklist P5, D4

- **Now:** psm.qvalue = 'PSM, as computed by the source', filled from qvalue/psm_qvalue/<engine>_qvalue but never from names containing 'global'. Yet psm.md recommends `global_qvalue` (experiment-level PSM q-value) on PSMs, and that score never reaches psm.qvalue. Run vs experiment scope is unstated, and the first matching name wins across engines.
- **Spec:** scores.md:38-47; psm.md:51,291; ontology.md:277
- **Reader impact:** psm.confident() filters on per-run, per-experiment or per-engine q-values depending on producer. Producers that follow the recommendation get a null qvalue.
- **Suggested:** psm.qvalue = experiment-wide PSM q-value of the final validator. Map the PSM `global_qvalue` into it, and declare scope in the footer.
- **Evidence:** `docs/spec/scores.md:38-47`; `docs/spec/psm.md:51`; `docs/spec/psm.md:291`; `docs/spec/ontology.md:277`

### CRITIC-11: Is psm_id truly unique: may (run, scan, peptidoform, charge) repeat? {#critic-11}

*psm* · ambiguous · confirmed · checklist P6

- **Now:** The composite names one peptidoform per spectrum. Writers only warn on duplicate ids; qpxc validate reports them as errors. For identical hits from several engines, PD keeps the lowest rank then lowest q-value. #229's 'one merged row with engine scores' is not in the spec.
- **Spec:** psm.md:9-15; versioning.md:120-127; identity.md:40; converter-coverage.md:323
- **Reader impact:** Joins and counts double-count when a file carries the same hit twice. PK uniqueness holds only if the file passed strict validation.
- **Suggested:** One row per (spectrum, peptidoform, charge) after merging engines, with engine scores side by side. Duplicates are invalid output.
- **Evidence:** `docs/spec/psm.md:9-15`; `docs/spec/versioning.md:120-127`; `docs/spec/converter-coverage.md:323`

### CRITIC-12: Which engine(s) produced each PSM in a multi-engine or consensus result? {#critic-12}

*psm* · unspecified · partly confirmed (see correction) · checklist P5, P6

- **Now:** Only a `consensus_support` count score exists; there is no engine column. PD's 'Identifying Node says which' has no defined field.
- **Spec:** scores.md:75; converter-coverage.md:323
- **Correction:** PD records the engine as cv_param identifying_node on each PSM (proteome_discoverer/psm_adapter.py:79), but that name is not registered in the spec. The gist holds: there is no defined engine field.
- **Reader impact:** Engine-specific scores appear on some rows and not others, with no way to know whose q-value or PEP was used.
- **Suggested:** Optional list<string> search_engines, or a registered cv_param name.
- **Evidence:** `docs/spec/scores.md:75`; `docs/spec/converter-coverage.md:323`

### CRITIC-19: Does protein_accessions list all mapped proteins, in what order, decoys/contaminants included? {#critic-19}

*psm* · ambiguous · confirmed · checklist P8

- **Now:** 'All proteins that the peptide maps to', null when there is no mapping. Order, razor-first, isoform completeness and decoy inclusion are unspecified. MaxQuant leaves Proteins blank for some reverse hits, and tool-mappings lists DIA-NN Protein.Group (a group, not a mapping).
- **Spec:** psm.md:57,289; tool-mappings.md:27; converter-coverage.md:281
- **Reader impact:** Shared-peptide and uniqueness analyses built from PSMs depend on the producer. The first element may or may not be the leading protein.
- **Suggested:** An unordered set of all database proteins the sequence maps to, matching pepmap's evidence proteins, with no group semantics.
- **Evidence:** `docs/spec/psm.md:57`; `docs/spec/tool-mappings.md:27`; `docs/spec/converter-coverage.md:281`

### CRITIC-20: Can a reader rely on psm.feature_id to find a PSM's quantified feature? {#critic-20}

*psm* · specified, worth keeping · partly confirmed (see correction) · checklist P10

- **Now:** feature_id is the authoritative optional FK. Null means no link was recorded, not low quality. The join is expected to be partial, and psm_ids is computed. MaxQuant, FragPipe and CDAP do not fill it.
- **Spec:** psm.md:75,147-175; data-model.md:162; converter-coverage.md:285
- **Correction:** The rule is correctly specified, but MaxQuant does fill psm.feature_id through the evidence id (95/1000 PSMs in the maxquant_simple run; maxquant/converter.py:64-85). converter-coverage.md:285 ('Not implemented') is stale.
- **Reader impact:** Keep the rule that null does not mean unidentified or unquantified. Readers should expect no links from most producers.
- **Suggested:** Keep, and add a footer feature_link (complete|partial|none) so a reader knows whether absence of a link is meaningful.
- **Evidence:** `docs/spec/psm.md:147-175`; `docs/spec/converter-coverage.md:285`

### CRITIC-M1: Can rows with equal identity_composite values carry different ids? {#critic-m1}

*psm* · contradictory · found by the verifier · checklist P1, P6

- **Now:** qpxc convert mzidentml supplies ids hashed from SII ids but still declares the footer composite run_file_name,scan,peptidoform,charge. On PXD054720 that gives 13,386 psm_ids but only 13,255 distinct composite values (the same peptidoform and charge at rank 1 and rank 2 with different Mascot scores), and qpxc validate reports VALID.
- **Spec:** serialization.md:304-305 and versioning.md:86-87: 'One identity composite value names one row identity, also for producer-supplied ids'. This conflicts with identity.md:40 and psm.md:12-15, which allow several hits of one peptidoform per spectrum.
- **Reader impact:** Readers who treat the declared composite as a natural key get duplicate rows, and re-deriving ids collides. The footer misdescribes row identity.
- **Suggested:** The declared composite must name rows uniquely (add rank or a spectrum_identification_item id). validate should flag composite duplicates.
- **Evidence:** `docs/spec/serialization.md:304-305`; `docs/spec/versioning.md:86-87`; `docs/spec/identity.md:40`; `work/verify-critic/mzid (qpxc convert mzidentml PXD054720; duckdb distinct counts; qpxc validate VALID)`

### CRITIC-M2: Does pass_threshold=true mean the spectrum's accepted identification? {#critic-m2}

*psm* · ambiguous · found by the verifier · checklist P1, D3

- **Now:** By default the mzIdentML converter keeps all ranks and copies each SII's passThreshold. On PXD054720, 11,517 rows pass over 9,707 spectra, and 1,810 of them are ranks 2-10 (qpx/converters/mzidentml/psm_records.py:62-63,140,146).
- **Spec:** psm.md:66 says only 'marks the PSM as passing its threshold'. Rank has no typed column (scores.md:58).
- **Reader impact:** Filtering on pass_threshold returns several peptides per spectrum, which inflates PSM and peptide counts and conflicts with rank-1 producers.
- **Suggested:** Define pass_threshold as acceptance of this row's peptidoform. Require a typed rank, and state that 'accepted identification' = rank 1 AND pass_threshold.
- **Evidence:** `docs/spec/psm.md:66`; `qpx/converters/mzidentml/psm_records.py:62-63`; `work/verify-critic/mzid rank x pass_threshold counts`

### CRITIC-M3: How are tied, indistinguishable hits (I/L, localisation isomers) represented? {#critic-m3}

*psm* · unspecified · found by the verifier · checklist P1, P2, F2

- **Now:** CDAP AmbiguousMatch=1 hits become separate PSM rows on the same scan, with identical scores and q-values, flagged only by the producer cv_param ambiguous_match (cdap/psm_adapter.py:27). CDAP then drops them from features (cdap/feature_adapter.py:7-9,116-119). In the fixture: 3 spectra, 6 rows (SCHTGIGR/SCHTGLGR, two K-acetyl positions).
- **Spec:** Nothing. These are neither rank alternatives nor chimeric (scope.md:41), and there is no tie or ambiguity-group field.
- **Reader impact:** Tied alternatives count as distinct PSMs and peptides, other producers keep only one, and their quantities vanish from the feature view silently.
- **Suggested:** Add an ambiguity_group or tie flag on psm, and define whether features keep one representative or drop the group.
- **Evidence:** `qpx/converters/cdap/feature_adapter.py:7-9`; `qpx/converters/cdap/feature_adapter.py:116-119`; `work/verify-critic/cdap psm rows sharing scan`
