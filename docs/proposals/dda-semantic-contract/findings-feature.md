# Findings: feature view

Part of the annex of the [DDA semantic contract proposal](index.md). Each entry is one verified finding from the audit of qpx `dev` at commit `32cbbdc` (specification 1.3). A skeptic re-checked every finding against the code; "partly confirmed" entries are true only as stated in their correction. `audit-run:` evidence refers to converting the repository's own fixtures (`tests/examples`) with that commit during the audit; those outputs are not committed.

174 findings, grouped by audit area.

## OpenMS converters

### OPENMS-05: PSM-level q-value written to feature.peptide_qvalue; precursor_qvalue never filled {#openms-05}

*feature* · contradictory · confirmed · checklist D4, F5

- **Now:** In the default (PSM-level) mode, peptide_qvalue = qvalue_of(first PID hit), which is the PSM's primary score whenever score_type is q-value/fdr (feature_adapter.py:149-160,344,745). The FDR level is not checked; only the Percolator peptide-level flag is recognised (feature_adapter.py:40,135-146). LFQ fixture: feature.peptide_qvalue equals psm.qvalue (0.00862) on 68/68 rows, and precursor_qvalue and global_qvalue are null on all rows.
- **Spec:** scores.md:43-44: peptide_qvalue is peptide level; precursor_qvalue takes qvalue/q_value. scores.md:83: peptide_qvalue is reserved for peptide-level confidence. converter-coverage.md:59 only says 'declared to be a q-value/FDR', with no level.
- **Reader impact:** feature.confident() and cross-producer comparisons mix PSM-level OpenMS values with true peptide-level values (CDAP) and precursor-level values (DIA-NN, Spectronaut). Peptide FDR is mis-stated.
- **Suggested:** Write the q-value of the identifying PSM to precursor_qvalue (or additional_scores 'psm_qvalue'). Use peptide_qvalue only when the source declares peptide-level FDR. Fix converter-coverage.md:59 to state the level.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:40,135-160,325-345,745`; `docs/spec/scores.md:43-44,83`; `docs/spec/converter-coverage.md:59,181`; `audit-run:openms/ProteomicsLFQ_1_out/openms.feature.parquet (npq=68, nprq=0)`

### OPENMS-06: Feature PEP/q-value from first PID of run, or best PSM across runs {#openms-06}

*feature* · diverges from spec · partly confirmed (see correction) · checklist F2, F5, D4

- **Now:** _confidence_by_run takes the first PID of each run in document order, not the best, with no peptide check (feature_adapter.py:325-345). In peptide-level mode it instead takes the best PSM of that peptidoform across the whole identification source, including other runs (lines 163-206). Re-attached removed IDs are ranked by primary score only, with ties in document order (lines 616-633). Test fx/pepfdr: the second ELEFANTK feature row gets PEP 0.01 from another spectrum.
- **Spec:** feature.md:138-139: scan/id_run are the best PSM's. feature.md:130: best PSM's values. feature.yaml:37-41 has feature PEP, but psm.md:286 says feature views do not carry PEP as a top-level field. #378 (open): lowest PEP, tie lowest scan, own run.
- **Correction:** PSM-level mode takes the first PID per run, not the best, with no peptide check (feature_adapter.py:325-345); confirmed. Peptide-level mode uses the source-wide best PSM, which is exactly what feature.md:130 documents, so that part conforms. The PEP contradiction is real: psm.md:286 says feature has no top-level PEP, but feature.md:56 and feature.yaml:37-41 define one.
- **Reader impact:** Feature PEP and q-value can describe a different spectrum, run, or even peptide than the row's scan. Their meaning also changes with the FDR mode.
- **Suggested:** Define one 'identifying PSM' per row: own run, tool's pick, else lowest PEP, tie lowest scan. All confidence fields come from it. Resolve the psm.md:286 vs feature.md:56 contradiction about feature PEP.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:163-206,325-345,616-633,736`; `docs/spec/feature.md:56,130,138-139`; `docs/spec/psm.md:286`; `qpx/core/data/schemas/feature.yaml:37-41`; `test: audit-run:openms/fx/out_pepfdr/openms.feature.parquet`

### OPENMS-07: Conflicting identifications: first PID defines peptide; others leak scans, no alternative row {#openms-07}

*feature* · contradictory · confirmed · checklist F2, F10

- **Now:** The peptide comes from pids[0].getHits()[0] (feature_adapter.py:724). _scan_by_run adds the spectra of every PID in the run, including other peptides (303-322). A conflicting run nulls id_run_file_name, identification_type (532-534) and the protein fields (485-487). No row is emitted for the competing peptidoform. Test fx/conflict: the ELEFANTK row has scan [24221,24300], where 24300 identified IMTELLERK; id_run, identification_type and anchor_protein are null, while PEP and q-value are still filled.
- **Spec:** feature.md:139: scan = scan of the best PSM that identified the feature. feature.md:26-28: for OpenMS, scan holds every spectrum reference for the run. Nothing covers conflicting identifications of one feature.
- **Reader impact:** feature.scan points to spectra of other peptides. Readers joining scan to psm get the wrong peptide. Quantity is attributed to an arbitrary peptidoform (document order).
- **Suggested:** The contract defines conflict handling: one row per peptidoform (marked ambiguous) or drop. scan holds only the identifying PSM's spectrum. Fix the feature.md:26-28 vs :139 contradiction.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:303-322,475-488,523-548,724`; `docs/spec/feature.md:26-28,139`; `test: audit-run:openms/fx/out_conflict/openms.feature.parquet`

### OPENMS-08: feature.scan concatenates several spectra; dedup merge destroys component structure {#openms-08}

*feature* · contradictory · partly confirmed (see correction) · checklist F10, P7

- **Now:** feature.scan is the flat concatenation of each PID's native-ID components, e.g. Waters [1,1,345]+[2,1,345] gives [1,1,345,2,1,345] (feature_adapter.py:303-322). When duplicate peaks are merged, `sorted(set(target.scan) | set(source.scan))` is applied to the flat ints (feature_dedup.py:293), giving [1,2,345]. scan is also part of the feature identity composite (feature_dedup.py:260).
- **Spec:** scan.md:89-92,124: the OpenMS converter 'can combine components from multiple supporting spectra'; do not treat the array as one native ID. feature.md:139: scan of the best PSM. producer-guide.md:23: scan number(s) with scan_format.
- **Correction:** The behaviour is real, but the line refs are wrong: feature_dedup.py has 166 lines; the composite is at :15 and the sorted(set()) merge at :48. The concatenation is documented design (scan.md:90-92,123-124; feature.md:26-28), so it contradicts only feature.md:139 and :301. The merge runs only for duplicate peaks of one run (:93-97,111-131), a rare path.
- **Reader impact:** feature.scan cannot be split back into spectra, and merged rows become meaningless. feature_id depends on how many PSMs were attached.
- **Suggested:** feature.scan = the identifying PSM's scan only, and other spectra are reached through psm.feature_id. Drop scan from the OpenMS feature composite, or use a list<list<int32>>.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:303-322`; `qpx/converters/openms_consensus/feature_dedup.py:260,278-299`; `docs/spec/scan.md:89-92,124`; `docs/spec/feature.md:139`; `docs/spec/producer-guide.md:23`

### OPENMS-09: Feature row = consensus element × run; several rows per peptidoform/charge/run {#openms-09}

*feature* · contradictory · partly confirmed (see correction) · checklist F1, F7

- **Now:** There is one row per consensus feature per run with a positive element (feature_adapter.py:348-369,734-757). For isobaric input a consensus feature is one reporter spectrum, so a row is a quant scan, not a chromatographic peak. The composite is (peptidoform,charge,run,rt,scan,observed_mz,consensus_rt) (feature_dedup.py:260). Test fx/pepfdr: two ELEFANTK/2+ rows in one run. Only byte-identical peaks are deduplicated (feature_dedup.py:338-342).
- **Spec:** data-model.md:64: 'one quantified chromatographic peak'. producer-guide.md:12: 'One row per precursor, charge and run'. feature.md:26-28 and identity.md:35 give the OpenMS composite. convert.md:1094 says 'one row per (peptidoform, charge, run, rt)'.
- **Correction:** The gist holds. Correct refs: composite at feature_dedup.py:15; peak-dedup key (run, peptidoform, charge, rt, intensities) at :93-97, not :260/:338. feature.md:11-12 lets converters declare another entity, and identity.md:35 lists the OpenMS composite, so the grain is only implied by identity_composite. The clash with data-model.md:64, producer-guide.md:12 and convert.md:1094 is real.
- **Reader impact:** Readers assuming one row per precursor and run double-count or average incorrectly. TMT rows are PSM-level quant events, so per-run peptide intensity needs aggregation the spec does not describe.
- **Suggested:** The contract states the row grain per producer, for example a footer `feature_grain=precursor_run|peak|quant_scan`. Fix producer-guide.md:12 and convert.md:1094 to match.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:348-369,734-757`; `qpx/converters/openms_consensus/feature_dedup.py:260,338-342`; `docs/spec/data-model.md:64`; `docs/spec/producer-guide.md:12`; `docs/guide/convert.md:1094`; `test: audit-run:openms/fx/out_pepfdr/openms.feature.parquet`

### OPENMS-10: Unidentified consensus features silently dropped; no 'unidentified' identification_type {#openms-10}

*feature* · unspecified · confirmed · checklist F3, D1, G2

- **Now:** feature_records_for_cf returns [] when the consensus feature has no PID or its first PID has no hits (feature_adapter.py:714-716). Fixture IsobaricWorkflow: the element with an empty PID on scan 24220 gives no feature row. The identification_type enum is validated as direct/transferred/secondary only (qpx/core/data/validation.py:71). Dropped features also never reach the pg roll-up.
- **Spec:** feature.md:7-10 allows unidentified features with a producer id. psm.md:175 documents that OpenMS omits them. #378 left open whether unidentified features belong in this view. No enum value covers them.
- **Reader impact:** Readers cannot tell 'not detected' from 'detected but unidentified'. Completeness and missing-value analyses differ by producer.
- **Suggested:** Decide explicitly: either the DDA feature view holds identified features only (declared in the footer), or unidentified rows are allowed with identification_type='unidentified' and null peptidoform.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:714-716`; `qpx/core/data/validation.py:71`; `docs/spec/feature.md:7-10`; `docs/spec/psm.md:175`; `tests/examples/openms_topp/IsobaricWorkflow_1_out.consensusXML:259-261`

### OPENMS-11: 'transferred' inferred heuristically; also set when the run's own conflicting ID was removed {#openms-11}

*feature* · ambiguous · confirmed · checklist F4, D1

- **Now:** transferred = the run has no PID on the feature and every PID's run is known (feature_adapter.py:541-546). FFId_category (internal/external), which is present on every fixture PID, is ignored. IDConflictResolver-removed PIDs of another peptide are not re-attached (604-613). Test fx/idconflict: BSA3_F1 changes from direct to transferred (scan [], PEP null) although its own MS2 identified DGDIEAEISR there; the same case without the resolver gives null (fx/conflict). 'secondary' is never written.
- **Spec:** feature.md:145: transferred = matched from another run. converter-coverage.md:55 documents the heuristic. producer-guide.md:26: transferred rows have an empty scan.
- **Reader impact:** The MBR fraction is inflated, and 'transferred' rows include runs with contradicting own-run evidence. The label depends on upstream tool choices, not on what happened.
- **Suggested:** Use the tool's own flag (FFId_category) when present. Mark a run with a conflicting own-run ID as null or 'conflict', not 'transferred'. The contract says that transferred requires the tool's statement or no own-run MS2 at all.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:523-548,604-613`; `tests/examples/openms_topp/ProteomicsLFQ_1_out.consensusXML:1088 (FFId_category)`; `docs/spec/feature.md:145`; `docs/spec/converter-coverage.md:55`; `test: audit-run:openms/fx/idconflict.consensusXML -> out_idc_mod vs out_idc_orig`

### OPENMS-12: 1.3 optional fields: match_*, imputed_labels, rt_aligned, representative_label never filled {#openms-12}

*feature* · unspecified · partly confirmed (see correction) · checklist D2, D7

- **Now:** There is no code for match_score, match_qvalue, match_rt_delta, match_mz_delta, imputed_labels, rt_aligned or representative_label (grep of qpx/converters/openms_consensus). quant_ms_level is set only when id_scan_id exists (quant_scan.py:21-22), so it stays null on IsobaricWorkflow even though scan_id equals the PID spectrum_reference. parent_ion_fraction comes from quant_precursor_purity, else precursor_purity. identified_labels and channel_details are filled for MS1 labels, with channel qvalue and evidence null.
- **Spec:** feature.md:146-148,154-156; converter-coverage.md:56. These fields are null when the tool does not say. Nothing tells 'tool did not say' apart from 'converter did not derive'.
- **Correction:** The code facts hold. But converter-coverage.md:39-43 says a null can mean not reported, ambiguous or not calculated, and :56 documents quant_ms_level as null without id_scan_id. What is missing is an in-file coverage marker, plus OpenMS entries for rt_aligned and match_* in converter-coverage.
- **Reader impact:** Readers cannot tell whether a null rt_aligned, match_* or quant_ms_level means the information is absent or was not converted. For example, consensusXML element RTs may already be aligned.
- **Suggested:** The contract lists, per producer, which optional 1.3 fields are populated (a coverage footer or converter-coverage). quant_ms_level falls back to comparing scan_id with the PID spectrum_reference.
- **Evidence:** `qpx/converters/openms_consensus/quant_scan.py:6-23`; `qpx/converters/openms_consensus/feature_adapter.py:750-756`; `docs/spec/feature.md:146-156`; `docs/spec/converter-coverage.md:56`; `audit-run:openms/IsobaricWorkflow_1_out/openms.feature.parquet (quant_ms_level null)`

### OPENMS-13: Empty reporter/MS1 channels (it=0) omitted instead of null intensity entry {#openms-13}

*feature* · diverges from spec · confirmed · checklist D6, F6

- **Now:** Sub-elements with intensity <= 0 are skipped, so the label gets no entry at all (feature_adapter.py:359-361). A channel linked twice keeps the maximum (366-368). Fixture IsobaricWorkflow: the identified MYFINGERK has 4 of 10 TMT entries, and its 6 channels with it="0.0" are absent. In MS1LabeledWorkflow_7 the YVNDLAGHEK row has only its light entry. converter-coverage.md:60 documents the omission.
- **Spec:** intensities.md:123-126: a null intensity marks an identified-but-not-quantified label, such as 'an empty reporter channel of an identified plex'; no entry means not identified or no sample. producer-guide.md:45-47.
- **Reader impact:** Readers read missing TMT and SILAC channels as 'no sample' or 'not identified' rather than 'not quantified'. Imputation and missingness models get the wrong mechanism.
- **Suggested:** For an identified plex, write every channel of the run, with null for a 0 sentinel. Keep 0 only when the tool states a measured zero.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:348-369`; `tests/examples/openms_topp/IsobaricWorkflow_1_out.consensusXML:196-205`; `audit-run:openms/IsobaricWorkflow_1_out/openms.feature.parquet (MYFINGERK n=4)`; `docs/spec/intensities.md:123-127`; `docs/spec/producer-guide.md:45-47`; `docs/spec/converter-coverage.md:60`

### OPENMS-14: Feature intensity meaning (area/apex/reporter raw vs corrected) not recorded per row {#openms-14}

*feature* · unspecified · partly confirmed (see correction) · checklist F6

- **Now:** intensities = the consensus element `it` (feature_adapter.py:359-368): the FFId feature intensity for LFQ, the reporter intensity for TMT. The feature row carries no quantification_method cv_param. IsobaricAnalyzer isotope_correction and normalization appear only in provenance parameters (provenance.py:66,107-118). pg rows do get quantification_method (pg_adapter.py:414-429).
- **Spec:** intensities.md:20-21: intensities are raw/primary. pg.md:52-58: quantification_method is on pg only. Nothing for the feature level (area vs apex, isotope-corrected or not).
- **Correction:** Isotope correction and normalization are specified: intensities.md:211-213 and provenance.md:69 put isotope_corrected and normalization keys on the tool's provenance step, and provenance.py:66,107-118 writes them. What remains unspecified is the MS1 quantity type (apex, area or isotope sum) and a feature-level quantification_method.
- **Reader impact:** Mixing producers (MaxQuant intensity, DIA-NN Precursor.Quantity, OpenMS FFId) or corrected/uncorrected TMT in one analysis without knowing it.
- **Suggested:** Require a feature-level quantification_method (footer or cv_param), including 'isotope_corrected' and 'normalized' booleans for reporter ions.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:348-369`; `qpx/converters/openms_consensus/provenance.py:66,107-118`; `docs/spec/intensities.md:20-21`; `docs/spec/pg.md:52-58`

### OPENMS-15: feature.observed_mz is consensus centroid; 0.0 sentinel; calculated_mz=observed when charge 0 {#openms-15}

*feature* · diverges from spec · confirmed · checklist F10, F7, P7

- **Now:** For LFQ, observed_mz = cf.getMZ() (the consensus centroid, the same for every run), or 0.0 when absent (feature_adapter.py:645). The run's own element m/z is used only for MS1 labels (681-691). calculated_mz = observed_mz when charge is 0 (line 648). For MS1 labels, observed_mz is the lightest measured partner while calculated_mz is the unlabelled peptide (MS1Labeled_7 SAYWNEGHIR: calc 616.79, obs 619.80).
- **Spec:** feature.md:61-62: calculated_mz null when it cannot be calculated; observed_mz is the experimental m/z, null when not reported. intensities.md:91-92: calc/obs/mass_error describe the lightest partner.
- **Reader impact:** Per-run mass accuracy and recalibration checks are wrong. A 0.0 m/z or a calculated equal to the observed m/z looks like real data. observed_mz is part of the feature identity.
- **Suggested:** Use the per-run element m/z. Write null instead of 0.0, and null calculated_mz when charge is 0. For multiplets, take observed and calculated from the same partner.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:636-666,681-691`; `docs/spec/feature.md:61-62`; `docs/spec/intensities.md:91-92`; `audit-run:openms/MS1LabeledWorkflow_7_partial_out/openms.feature.parquet`

### OPENMS-27: Shared peptides: null pg fields (no razor), unlike MaxQuant's assigned group {#openms-27}

*feature* · inconsistent across producers · confirmed · checklist F9, P8

- **Now:** The group is resolved from the first PID's protein evidence: an exact membership match, else exactly one candidate group, else null (protein_groups.py:35-60). A conflict between runs also gives null (feature_adapter.py:475-488,728-729). unique = len(group)==1 (line 468). psm.protein_accessions is the union of evidence accessions across collapsed hits, canonicalised with decoy markers kept and no order guarantee (psm_adapter.py:72-83).
- **Spec:** feature.md:99-101: for shared peptides, the group the producer assigned (MaxQuant razor), and unique is null when unresolved. versioning.md:67-69: one unique definition. psm.md:57: all mapped proteins.
- **Reader impact:** Feature-to-pg joins drop OpenMS shared peptides that MaxQuant assigns. Roll-ups over the pg link differ by producer for the same data.
- **Suggested:** The contract requires a declared shared-peptide policy (`assigned`, `unassigned_null`, `all_groups`) in the footer, plus pepmap for the full mapping.
- **Evidence:** `qpx/converters/openms_consensus/protein_groups.py:35-60`; `qpx/converters/openms_consensus/feature_adapter.py:455-488,728-729`; `qpx/converters/openms_consensus/psm_adapter.py:72-83`; `docs/spec/feature.md:99-101`; `docs/spec/versioning.md:67-69`

### OPENMS-29: Dimethyl/mTRAQ multiplet peptidoform lacks lightest label; guide contradicts spec {#openms-29}

*feature* · contradictory · confirmed · checklist D5, P9, F1

- **Now:** The feature peptidoform is hit.getSequence(), the MS1LabeledWorkflow identity with labels removed (feature_adapter.py:636-639). The PSM uses the labelled matched sequence (channels.py:156-162). convert.md:1038-1040 says the feature peptidoform is 'without labels, also for dimethyl and mTRAQ'.
- **Spec:** peptidoform.md:20-23 and intensities.md:89-90: a multiplet row carries the lightest partner's peptidoform, with 'the lightest dimethyl or mTRAQ label inline'.
- **Reader impact:** Dimethyl and mTRAQ features from OpenMS and from DIA-NN plexDIA get different peptidoforms and ids. calculated_mz and the label state disagree.
- **Suggested:** Add the lightest channel's label for dimethyl and mTRAQ, and fix convert.md:1038-1040.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:636-639`; `qpx/converters/openms_consensus/channels.py:156-162`; `docs/guide/convert.md:1038-1040`; `docs/spec/peptidoform.md:20-23`; `docs/spec/intensities.md:89-90`

### OPENMS-M2: Identified-but-unquantified precursors get no feature row (MaxQuant writes null-intensity rows) {#openms-m2}

*feature* · inconsistent across producers · found by the verifier · checklist F3, D6, X1

- **Now:** Feature rows exist only for runs with a positive consensus element (feature_adapter.py:348-369,734). Unassigned 'not mapped' PIDs reach only psm. LFQ fixture: 47/115 PSMs have no feature_id; 19 (peptidoform, charge, run) have a PSM but no feature row; no feature row has a null intensity entry.
- **Spec:** producer-guide.md:45-46 and intensities.md:123-126: identified, not quantified = an entry with null intensity. MaxQuant writes MS/MS-only evidence rows with null intensity (converter-coverage.md:194). psm.md:152-156 only notes the unlinked PSMs.
- **Reader impact:** Feature-level missingness differs by producer: OpenMS 'identified, not quantified' looks like 'not detected', so imputation and completeness metrics differ for the same data.
- **Suggested:** The contract states whether the feature view must hold a null-intensity row for every identified precursor/run. Otherwise it declares (footer) that it holds quantified precursors only.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:348-369,714-716,734`; `tests/examples/openms_topp/ProteomicsLFQ_1_out.consensusXML (47 PIDs feature_id='not mapped')`; `audit-run:openms/ProteomicsLFQ_1_out: 47 psm feature_id null; 19 precursor/run without feature row`; `docs/spec/producer-guide.md:45-46`; `docs/spec/intensities.md:123-126`; `docs/spec/converter-coverage.md:194`

### OPENMS-M3: rt vs rt_aligned undetermined although the consensusXML records an RT alignment {#openms-m3}

*feature* · ambiguous · found by the verifier · checklist F10, P7, D7

- **Now:** feature.rt is the run's first element RT (feature_adapter.py:365,742) and psm.rt is the PID RT (psm_adapter.py:313). rt_aligned is never set. The LFQ fixture declares 'Retention time alignment' (ProteomicsLFQ_1_out.consensusXML:6,16). Element RT minus PID RT has a per-run median of -1.5 to 3.8 s, so both are on one clock. Whether that clock is raw or aligned could not be verified (no mzML).
- **Spec:** feature.md:66: rt is the time on this run's own clock; rt_aligned is the time after the tool's alignment. psm.md:50: rt is the MS2 scan RT.
- **Reader impact:** If the elements are transformed, OpenMS rt holds aligned times. Joins to mzML by RT, and cross-producer RT comparisons, are then silently off.
- **Suggested:** The contract requires rt to be the raw scan clock, and rt_aligned when the producer aligned. The OpenMS converter must check for original_RT or the transformation and fill both, or state the clock in the footer.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:365,742`; `qpx/converters/openms_consensus/psm_adapter.py:313`; `tests/examples/openms_topp/ProteomicsLFQ_1_out.consensusXML:6,16`; `cmd: element-vs-PID RT medians per map 3.76/2.86/2.81/-1.54/0.72/-0.22 s`; `docs/spec/feature.md:66`; `docs/spec/psm.md:50`

## MaxQuant converter

### MAXQUANT-01: MULTI-MATCH-MSMS is 'transferred' but keeps own-run id_run_file_name, scan and PEP {#maxquant-01}

*feature* · contradictory · confirmed · checklist D1, F4, F10

- **Now:** evidence_fields.py:16 maps MULTI-MATCH-MSMS to transferred, but feature_adapter.py:656-657 nulls id_run_file_name only for exact 'MULTI-MATCH'. Scan ('MS/MS scan number') and PEP are copied as they are (feature_adapter.py:560-588). Dedup demotes only 'MULTI-MATCH' (feature_adapter.py:264). converter-coverage.md:199 (transferred includes MULTI-MATCH-MSMS) and :203 ('null for MBR rows') disagree. No fixture contains this Type.
- **Spec:** feature.md:145: id_run_file_name is null for a transferred row. producer-guide.md:26: a transferred feature has an empty scan. #378: transferred rows have empty scan/id_run_file_name and null PEP.
- **Reader impact:** Filtering MBR rows with id_run_file_name IS NULL gives a different set than identification_type='transferred'. A transferred row can carry the scan and PEP of an MS/MS that did not identify it.
- **Suggested:** Contract: transferred means null id_run_file_name, empty scan and null PEP. Derive these from identification_type rather than from the raw Type string, or declare MULTI-MATCH-MSMS direct if MaxQuant semantics say so.
- **Evidence:** `qpx/converters/maxquant/evidence_fields.py:10-17`; `qpx/converters/maxquant/feature_adapter.py:655-662`; `qpx/converters/maxquant/feature_adapter.py:260-267`; `docs/spec/feature.md:145`; `docs/spec/producer-guide.md:26`; `docs/spec/converter-coverage.md:199,203`; `github #378`

### MAXQUANT-02: Type values outside the six-key map silently look like own-run identifications {#maxquant-02}

*feature* · unspecified · confirmed · checklist D1

- **Now:** _IDENTIFICATION_TYPES has 6 keys (evidence_fields.py:10-17). Any other Type, such as PEAK or a version-specific variant, gets identification_type null with no warning. Its id_run_file_name stays run_file_name (feature_adapter.py:657) and dedup ranks it as identified. The fixtures hold only MSMS, MULTI-MSMS, MULTI-SECPEP and MULTI-MATCH, so I could not observe real PEAK rows.
- **Spec:** feature.md:145: null when the tool does not say. converter-coverage.md:199 lists 6 Types only. The enum check (qpx/core/data/validation.py:71) allows direct/transferred/secondary and has no 'unidentified' value.
- **Reader impact:** These rows look identified in their own run through id_run_file_name. WHERE identification_type='direct' silently drops them.
- **Suggested:** The contract enumerates every MaxQuant Type and its value. An unknown Type warns and leaves id_run_file_name null. Decide whether an 'unidentified' value exists.
- **Evidence:** `qpx/converters/maxquant/evidence_fields.py:10-24`; `qpx/converters/maxquant/feature_adapter.py:656-657`; `qpx/core/data/validation.py:71`; `docs/spec/converter-coverage.md:199`

### MAXQUANT-07: Many feature rows per (peptidoform, charge, run): direct, MBR and MS/MS-only mixed {#maxquant-07}

*feature* · contradictory · confirmed · checklist F1, F2, F7, X1

- **Now:** Each evidence row becomes one feature (243,822 rows in, 243,821 out). Full fixture: 14,086 (peptidoform,charge,run) keys have more than one row (61,719 rows, max 101). 12,738 keys have at least 2 quantified rows, and 9,083 keys mix direct and transferred. Composite: peptidoform,charge,run_file_name,rt_start,rt_stop,scan (feature_adapter.py:73).
- **Spec:** producer-guide.md:12: 'One row per precursor, charge and run'. identity.md:32 gives a MaxQuant composite that allows several rows. anndata.md:153 sums duplicates by default. #217: a feature is a chromatographic peak.
- **Reader impact:** Aggregating by precursor and run double-counts or mixes MS/MS and MBR quantities; MuData sums them by default.
- **Suggested:** Contract defines a feature as one chromatographic peak and fixes producer-guide. Readers choose rows by identification_type, or a precursor-level key is added.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:69-73,186-193`; `docs/spec/producer-guide.md:12`; `docs/spec/identity.md:32`; `docs/spec/anndata.md:153`; `audit-run:maxquant/full/maxquant.feature.parquet`

### MAXQUANT-08: One MS1 peak assigned to two peptidoforms yields two features with full intensity {#maxquant-08}

*feature* · unspecified · confirmed · checklist F2

- **Now:** MaxQuant can assign one isotope pattern to several peptides, and each evidence row becomes a feature with the full Intensity. Full fixture: 312 (run, m/z, charge, Intensity) groups have 2 or more Modified sequences. Example: UPS1_250amol_R2 TEAANKVEETKK and KTEAANKVEETK, m/z 674.359356, Intensity 1,586,400, both MULTI-MSMS. No flag or shared key is written.
- **Spec:** Nothing. feature.md:145 'secondary' covers chimeric MS/MS only.
- **Reader impact:** The same signal is counted twice in peptide and protein rollups, and readers cannot detect the shared peaks.
- **Suggested:** Define the rule for conflicting identifications of one peak: keep the best only, or one row each with a shared peak_id/consensus id.
- **Evidence:** `audit-run:maxquant/fullin/evidence.txt (query by Raw file,m/z,Charge,Intensity)`; `qpx/converters/maxquant/feature_adapter.py:430-486`

### MAXQUANT-09: MS/MS-only ('MSMS') rows are features with a synthetic 60 s window {#maxquant-09}

*feature* · ambiguous · confirmed · checklist F1, F3, F10, D6

- **Now:** Type MSMS rows (12,758 in full, 161 in simple) get intensities [{LFQ,null}] and spectral_count 1 (feature_adapter.py:523-531). observed_mz and mass_error are null, rt is the MS/MS RT, scan is set, and rt_start/rt_stop are a constant 60 s window (full-fixture average 60.0 s). Each MS/MS is its own row.
- **Spec:** intensities.md:123-124: a null entry means identified but not quantified, and cites MaxQuant MS/MS-only rows. feature.md:69-70 does not define the window. #217: a feature is a peak.
- **Reader impact:** Non-peaks look like features with a measured window, and the feature count is inflated.
- **Suggested:** Keep the null-intensity rule. Null rt_start/rt_stop when no peak was integrated. Consider one identified-not-quantified row per precursor and run.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:504-531`; `docs/spec/intensities.md:121-127`; `audit-run:maxquant/full/maxquant.feature.parquet`

### MAXQUANT-10: rt is on the raw clock but rt_start/rt_stop are calibrated; docs disagree {#maxquant-10}

*feature* · contradictory · confirmed · checklist F10, P7

- **Now:** rt comes from uncalibrated 'Retention time' (column_mappings.yaml:161). rt_aligned and rt_start/rt_stop come from 'Calibrated retention time' and its start/finish (feature_adapter.py:574-585). Full fixture: rt lies outside [rt_start,rt_stop] in 83,413/243,821 rows, rt_aligned lies inside in all rows, and |rt-rt_aligned| reaches 153 s.
- **Spec:** feature.md:65-66: rt uses the run's own clock and rt_aligned the aligned one. feature.md:69-70 does not say which clock the window uses. tool-mappings.md:50-51 says MaxQuant rt_start is '---', but converter-coverage.md:204 says calibrated.
- **Reader impact:** Peak width, RT containment and XIC extraction on the raw clock use the wrong windows.
- **Suggested:** The contract puts rt_start/rt_stop on the same clock as rt. The calibrated window goes to separate aligned fields. Fix tool-mappings.
- **Evidence:** `qpx/converters/column_mappings.yaml:161-164`; `qpx/converters/maxquant/feature_adapter.py:574-585`; `docs/spec/feature.md:65-70`; `docs/spec/tool-mappings.md:50-51`; `docs/spec/converter-coverage.md:204`

### MAXQUANT-11: MBR rows are well flagged, but match delta sign and donor are unspecified {#maxquant-11}

*feature* · ambiguous · confirmed · checklist F4, D2, F5

- **Now:** MULTI-MATCH rows get transferred, null id_run_file_name, empty scan, null PEP and null additional_scores (104,491 rows in full). match_score, match_rt_delta ('Match time difference'×60) and match_mz_delta are copied. match_qvalue comes from 'Match q-value', which is empty in full, so 0 rows are filled. The donor run is unknown.
- **Spec:** feature.md:145-147: deltas are the 'difference from its donor' in s and Th, with no sign convention. producer-guide.md:26.
- **Reader impact:** The MBR flag is reliable. Delta signs may differ between producers, and older MaxQuant gives no match q-value.
- **Suggested:** Keep the rule. Specify delta = this run minus donor and verify MaxQuant's sign.
- **Evidence:** `qpx/converters/maxquant/evidence_fields.py:19-30`; `qpx/converters/maxquant/feature_adapter.py:648-662`; `docs/spec/feature.md:145-147`; `audit-run:maxquant/full/maxquant.feature.parquet`

### MAXQUANT-12: TMT MBR rows get one null entry per channel ('identified, not quantified') {#maxquant-12}

*feature* · diverges from spec · partly confirmed (see correction) · checklist D6, F4, F6

- **Now:** The TMT path writes one entry per SDRF channel and nulls values <=0 (feature_adapter.py:538-544). MULTI-MATCH rows have no MS/MS, so all reporters are 0. All 1,469 MBR rows in the simple fixture have 8 null entries. The MS1 'Intensity' and spectral_count are not stored for TMT.
- **Spec:** producer-guide.md:45-46: not identified in this label or run means no entry. intensities.md:123-126: a null entry means identified but not quantified.
- **Correction:** Facts hold, but this is no divergence. A transferred row is an identification, and intensities.md:123-124 makes unquantified labels of an identified feature null entries. The gaps are elsewhere: nothing says whether isobaric transferred rows without a reporter scan should exist, and the MS1 Intensity and spectral count of EVERY isobaric feature are dropped.
- **Reader impact:** Transferred rows look like identified plexes without quantities. They carry no data yet inflate feature counts.
- **Suggested:** Isobaric transferred rows without a reporter scan get null intensities (no entries) or are dropped. Keep the MS1 Intensity as an additional_intensity.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:516-554`; `docs/spec/producer-guide.md:45-46`; `docs/spec/intensities.md:123-126`; `audit-run:maxquant/simple/maxquant.feature.parquet`

### MAXQUANT-13: TMT reporter intensities are summed over several MS/MS; quant_scan is null {#maxquant-13}

*feature* · unspecified · confirmed · checklist F6, D7

- **Now:** evidence 'Reporter intensity N' covers all MS/MS of the evidence: 'Reporter intensity count' reaches 9, and 73 simple-fixture rows are above 1. feature.scan holds only the best MS/MS. quant_scan and quant_ms_level are null. 'Reporter intensity count', 'Reporter PIF' and 'Reporter fraction' are dropped.
- **Spec:** feature.md:150-156: quant_scan is the scan the reporter intensities were read from (singular), and quant_ms_level is 2 or 3.
- **Reader impact:** Readers assume single-scan reporter quantities as with other producers. They cannot weight multi-scan sums or tell MS2 from MS3.
- **Suggested:** Allow several scans in quant_scan or add a reporter_scan_count. MaxQuant fills it from 'MS/MS IDs' or the count column.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:533-554`; `docs/spec/feature.md:150-156`; `tests/examples/maxquant/maxquant_simple/evidence.txt (Reporter intensity count 0..7)`

### MAXQUANT-15: TMT evidence without an SDRF is written as label 'LFQ' with MS1 intensity {#maxquant-15}

*feature* · diverges from spec · confirmed · checklist F6, P9, X1

- **Now:** Without SDRF channels the converter only warns and takes the LFQ path (feature_adapter.py:203-207,519-531). All 9,301/9,301 rows got label 'LFQ' (audit-run:maxquant/simple_nosdrf). The TMT label is not added to peptidoforms. CLI help says the SDRF is required for feature and pg (qpx/cli/convert/maxquant.py:71), but this is not enforced. tool-mappings.md:52 omits reporter intensities.
- **Spec:** intensities.md:53,70: 'LFQ' means label-free. producer-guide.md:39.
- **Reader impact:** An isobaric dataset is indistinguishable from label-free, and its peptidoforms and feature_ids differ from an SDRF-backed conversion.
- **Suggested:** Fail, or use a non-LFQ label, when reporter columns exist without a channel map.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:197-208,519-531`; `qpx/cli/convert/maxquant.py:69-73`; `docs/spec/intensities.md:53,70`; `docs/spec/tool-mappings.md:52`; `audit-run:maxquant/simple_nosdrf/maxquant.feature.parquet`

### MAXQUANT-16: LFQ primary may become 'LFQ intensity'; intensity meaning and method not declared {#maxquant-16}

*feature* · diverges from spec · partly confirmed (see correction) · checklist F6, G2

- **Now:** The LFQ path uses 'LFQ intensity' when it is >0, else 'Intensity' (feature_adapter.py:523-525). Provenance nonetheless says normalization=none (converter.py:319-320). The pg view writes no quantification_method cv_param (pg_adapter.py:398-415). Nothing records that MaxQuant 'Intensity' is a summed XIC over isotopes and time rather than an apex. evidence.txt normally lacks 'LFQ intensity', so the first issue is latent.
- **Spec:** intensities.md:17-22: primary intensities are raw. pg.md:52-60: write quantification_method where the converter knows it. intensities.md:211-213.
- **Correction:** Line refs are wrong: the provenance is at converter.py:121-122. evidence.txt has no 'LFQ intensity' column, so the LFQ-preference path is latent. The missing pg quantification_method, and the undeclared summed-XIC meaning, are real.
- **Reader impact:** Primary intensities could be normalised while provenance says raw. Area versus apex meaning is unknown across producers.
- **Suggested:** Primary is always 'Intensity' and LFQ goes to additional_intensities. Write quantification_method on feature and pg.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:521-531`; `qpx/converters/maxquant/converter.py:319-320`; `qpx/converters/maxquant/pg_adapter.py:398-415`; `docs/spec/pg.md:52-60`; `docs/spec/intensities.md:211-213`

### MAXQUANT-26: Feature PEP level is undefined; no q-values, so feature.confident() returns nothing {#maxquant-26}

*feature* · contradictory · confirmed · checklist F5, P5, D4

- **Now:** feature PEP comes from evidence 'PEP', which per #378 is the Best MS/MS PEP (75,493/75,493 rows on PXD003133). It is null for MULTI-MATCH. peptide_qvalue, precursor_qvalue and global_qvalue are 0% filled. feature.confident() uses COALESCE(precursor_qvalue, peptide_qvalue) (qpx/core/data/feature.py:52-54), so it returns 0 MaxQuant rows.
- **Spec:** psm.md:285-286: PEP is defined only in the PSM view. feature.md:56 lists PEP 'for the peptide match'. scores.md:38-47: typed q-values.
- **Reader impact:** Readers do not know whether feature PEP is per-PSM or aggregated, and the standard filters drop all MaxQuant data.
- **Suggested:** Define feature PEP as the best identifying PSM's PEP in this run. Define a filter fallback when q-values are absent.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:587-601`; `qpx/core/data/feature.py:52-54`; `docs/spec/psm.md:285-286`; `docs/spec/feature.md:56`; `github #378`

### MAXQUANT-28: Without proteinGroups.txt, pg_accessions holds leaders of several groups {#maxquant-28}

*feature* · contradictory · confirmed · checklist F9, P8

- **Now:** With proteinGroups.txt: pg_accessions is the razor group's 'Protein IDs', the anchor is the first 'Majority protein IDs', and unique means one group with one protein (feature_adapter.py:319-348). Without it, or when ids do not resolve: pg_accessions is 'Leading proteins' (one leader per group) and unique means one 'Proteins' entry (feature_adapter.py:350-366).
- **Spec:** feature.md:99: 'Full membership of the peptide's protein group'. tool-mappings.md:55,62 documents the fallback.
- **Reader impact:** The meaning of pg_accessions depends on the inputs supplied, and the feature-pg softlink breaks.
- **Suggested:** Fallback: null pg_accessions and unique, with all mappings in pepmap.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:319-366`; `docs/spec/feature.md:99-101`; `docs/spec/tool-mappings.md:55,62`

### MAXQUANT-29: Dedup can merge a secondary identification into a direct one {#maxquant-29}

*feature* · unspecified · confirmed · checklist F2, D1

- **Now:** Rows sharing the composite collapse. Precedence: non-'MULTI-MATCH' first, then lowest PEP, then lowest id (feature_adapter.py:230-295). Full fixture: evidence 167037 (MULTI-MSMS) and 167038 (MULTI-SECPEP), same peptidoform and scan 10817, become one 'direct' row. Secondary rows drop from 5,110 to 5,109.
- **Spec:** Nothing on collapse precedence for features. #262 (keep lowest PEP) concerns PSM hits.
- **Reader impact:** A secondary identification can be silently merged into a direct one, and collapse rules differ between producers.
- **Suggested:** The contract sets the best-identification rule: the tool's pick, else lowest PEP, ties to lowest scan (#378).
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:56-64,230-295`; `audit-run:maxquant/fullin/evidence.txt ids 167037/167038`

### MAXQUANT-30: Unidentified MS1 features are never emitted (allPeptides.txt not read) {#maxquant-30}

*feature* · unspecified · partly confirmed (see correction) · checklist F3, D1

- **Now:** Only evidence.txt is converted, so every feature row has a peptidoform from an MS/MS or MBR identification. allPeptides.txt is not supported (converter.py:237-297).
- **Spec:** feature.md:7-10 allows unidentified features with producer ids. The identification_type enum has no 'unidentified' value (validation.py:71). #378 leaves this undecided.
- **Correction:** The gist holds, but converter.py has only 220 lines, so the cited :237-297 does not exist.
- **Reader impact:** Acceptable, but readers cannot tell that MaxQuant's feature view excludes unidentified peaks.
- **Suggested:** Declare per file whether unidentified features are included (a footer key).
- **Evidence:** `qpx/converters/maxquant/converter.py:237-297`; `docs/spec/feature.md:7-10`; `qpx/core/data/validation.py:71`

### MAXQUANT-31: SILAC rule worth keeping: lightest-partner feature, labelled PSM, identified_labels {#maxquant-31}

*feature* · specified, worth keeping · confirmed · checklist F7, P9, D7

- **Now:** There is one feature per multiplet, whose peptidoform has no SILAC label (feature_adapter.py:441). Per-label intensities come from 'Intensity L/M/H', 0 becomes null, and ratios sit on the numerator (silac.py:52-86). identified_labels comes from 'Labeling State' plus msms.txt states (silac.py:89-129). The PSM carries the label mods and psm.label (psm_adapter.py:151-156,259). channel_details and representative_label are never filled.
- **Spec:** intensities.md:79-112 and peptidoform.md:17-23 define this scheme.
- **Reader impact:** Consistent. Readers expect feature and PSM peptidoforms to differ by the label.
- **Suggested:** Keep. Require channel_details where the tool reports per-label m/z and RT.
- **Evidence:** `qpx/converters/maxquant/silac.py:52-129`; `qpx/converters/maxquant/psm_adapter.py:129-156`; `docs/spec/intensities.md:79-112`

### MAXQUANT-M2: MaxQuant evidence id and raw Type not traceable on feature rows {#maxquant-m2}

*feature* · diverges from spec · found by the verifier · checklist D1, F2, X1

- **Now:** Feature cv_params is always None (feature_adapter.py:474): 0/243,821 filled in full. The evidence 'id', the ids collapsed by dedup (feature_adapter.py:284-289) and the raw 'Type' are not kept. An unmapped Type (TIMS-*, ISO-SECPEP) leaves identification_type null with no trace of the source value.
- **Spec:** feature.md:6-7: for identified features 'any producer ID remains traceable in cv_params'.
- **Reader impact:** Readers cannot go back to evidence.txt, see which evidence rows were merged, or recover the identification origin when the mapping fails.
- **Suggested:** Write cv_params maxquant_evidence_id (all merged ids) and maxquant_type. The contract requires the producer row id(s) on every feature.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:474,284-289`; `docs/spec/feature.md:6-7`; `audit-run:maxquant/full/maxquant.feature.parquet (cv_params 0 filled)`

## FragPipe converter

### FRAGPIPE-10: Feature row is a whole experiment's quantity, labelled with its first raw file {#fragpipe-10}

*feature* · contradictory · confirmed · checklist F1, X1

- **Now:** One row per (experiment, peptidoform, charge, CV) (feature_adapter.py:52-57,330-365). For a multi-run experiment the intensity aggregates all its fractions. run_file_name = runs[0] of the lexicographically sorted members (pg_adapter.py:96, feature_adapter.py:350). Verified: experiment B (F2, F10, F1) gives one row with run_file_name 'runB_F1'.
- **Spec:** data-model.md:64-65 defines a feature as 'a peptidoform+charge in one run', but data-model.md:73-74 says 'per upstream Feature entity'. feature.md:89: 'the run file name that contains the feature'. intensities.md:136: each LFQ feature row = one MS run. producer-guide.md:12: one row per precursor, charge and run.
- **Reader impact:** A per-run analysis puts whole-experiment intensity on one arbitrary fraction, and the other fractions look empty. Feature granularity differs between FragPipe and every per-run producer.
- **Suggested:** The contract should allow unit-level features explicitly: run_file_name null or a designated representative, quantification_unit_id required, and member runs stated (e.g. a feature grouped_runs). Otherwise unit-level quantities belong only in pg.
- **Evidence:** `qpx/converters/fragpipe/feature_adapter.py:52-57`; `qpx/converters/fragpipe/feature_adapter.py:330-365`; `qpx/converters/fragpipe/pg_adapter.py:96`; `docs/spec/data-model.md:64-65`; `docs/spec/data-model.md:73-74`; `docs/spec/feature.md:89`; `docs/spec/intensities.md:136`; `docs/spec/producer-guide.md:12`; `run: out/fp.feature.parquet unit B -> run_file_name runB_F1`

### FRAGPIPE-11: quantification_unit_id is a free string with no link to its runs {#fragpipe-11}

*feature* · ambiguous · partly confirmed (see correction) · checklist F1, F10

- **Now:** quantification_unit_id = the FragPipe experiment name (feature_adapter.py:351). For a TMT plex it is the user's --output-prefix (converter.py:56-59). Member runs are not on the feature; a reader recovers them only through pg.grouped_runs containing run_file_name. Two plexes converted with the same prefix would collide.
- **Spec:** feature.md:90: 'Quantification unit represented by the quantity when it may aggregate multiple raw runs'. It defines no format and no relation to pg.grouped_runs or run.parquet. identity.md:34 uses it in the FragPipe composite.
- **Correction:** The spec does define the concept: a quantification unit is one grouped_runs set (data-model.md:29-30), and views.md:212 names units by their runs joined with ';'. The real gap is that feature.quantification_unit_id is a free label with no defined mapping to that set or to pg.grouped_runs. 'Collision' is overstated: the same prefix in one folder overwrites the files.
- **Reader impact:** A reader cannot resolve a unit to its runs or samples from the feature view alone. Unit ids from different producers or conversions are not comparable.
- **Suggested:** Define the unit as the grouped_runs set (as in pg), or require its member runs in a view (run.parquet or a feature column). The id should be derived, not user-chosen.
- **Evidence:** `qpx/converters/fragpipe/feature_adapter.py:351`; `qpx/converters/fragpipe/converter.py:56-59`; `docs/spec/feature.md:90`; `docs/spec/identity.md:34`

### FRAGPIPE-12: FAIMS: CV is in the feature identity but not in the PSM lookup key {#fragpipe-12}

*feature* · diverges from spec · confirmed · checklist F2, F5, F10, P6

- **Now:** The feature composite includes compensation_voltage (feature_adapter.py:52-57), but the PSM lookup key is (run, peptidoform, charge) (feature_adapter.py:341,470). Verified: features at CV -45 and CV -60 both took the CV -60 PSM (scan 500, PEP 0.001). The PSM view stores CV only as a cv_param string, 'compensation_voltage' (psm_adapter.py:317-319); psm has no CV column.
- **Spec:** feature.md:138-139: scan and id_run_file_name are those of the best PSM that identified THE feature. feature.md:76 and identity.md:34: CV separates FragPipe features. Issue #378 todo: 'add CV to the key'.
- **Reader impact:** For FAIMS data, the scan, PEP, observed_mz and mass error of a feature can come from a spectrum acquired at another CV, which may be another ion population.
- **Suggested:** Key the PSM lookup by CV. Give psm a typed compensation_voltage column, so that the PSM and feature identities line up.
- **Evidence:** `qpx/converters/fragpipe/feature_adapter.py:52-57`; `qpx/converters/fragpipe/feature_adapter.py:341`; `qpx/converters/fragpipe/feature_adapter.py:470`; `qpx/converters/fragpipe/psm_adapter.py:317-319`; `docs/spec/feature.md:76`; `docs/spec/feature.md:138-139`; `run: out/fp.feature.parquet rows CV -45 and -60 both scan [500]`

### FRAGPIPE-13: 'Best PSM' is qpx's own lowest-PEP choice; feature PEP role is contradictory {#fragpipe-13}

*feature* · contradictory · confirmed · checklist F2, F5

- **Now:** The feature takes the PSM with the lowest PEP, then the lowest scan, among same-run PSMs of the peptidoform and charge (feature_adapter.py:72-74,476). That PSM supplies PEP, scan, m/z, ion mobility and decoy status. #378 found this differs from FragPipe's own Expectation-based pick (212 ties differ on PXD081924). Feature PEP is therefore a PSM-level PEP of one spectrum, not a peptide or precursor confidence.
- **Spec:** feature.md:56 defines feature.posterior_error_probability. psm.md:286 says the feature and peptide views do NOT carry PEP as a top-level field. feature.md:138-139,301 say 'best PSM' without a ranking rule.
- **Reader impact:** Feature PEP means different things across producers (best PSM, aggregated, or none). Rankings differ from the tool's own, so feature scan and psm rows disagree with FragPipe reports.
- **Suggested:** Adopt #378: the identifying PSM is the tool's pick, else the lowest PEP with ties broken by lowest scan; it is PSM-level. Fix psm.md:286 and state the level of feature PEP.
- **Evidence:** `qpx/converters/fragpipe/feature_adapter.py:72-74`; `qpx/converters/fragpipe/feature_adapter.py:476`; `docs/spec/feature.md:56`; `docs/spec/psm.md:286`; `docs/spec/feature.md:138-139`; `github #378`

### FRAGPIPE-14: IonQuant Match Type mapped to identification_type direct/transferred {#fragpipe-14}

*feature* · specified, worth keeping · confirmed · checklist F4, D1, D2

- **Now:** '<exp> Match Type' maps MS/MS to direct and MBR to transferred; anything else is null (feature_adapter.py:51,339,361). Transferred rows have id_run_file_name null, scan [] and PEP null (verified, experiment C). The value is validated as a closed enum of direct, transferred and secondary (core/data/validation.py:71). match_score, match_qvalue and match_*_delta are never filled.
- **Spec:** feature.md:145-147: identification_type direct|transferred|secondary, null when the tool does not say; a transferred row has null id_run_file_name. producer-guide.md:26: an MBR row has an empty scan.
- **Reader impact:** For one-run LFQ experiments, MBR rows are correctly flagged and filterable. The contract should keep this rule.
- **Suggested:** Keep the mapping. Record IonQuant's MBR FDR or ion-FDR settings in provenance, since no per-row match_qvalue exists.
- **Evidence:** `qpx/converters/fragpipe/feature_adapter.py:51`; `qpx/converters/fragpipe/feature_adapter.py:339-361`; `qpx/core/data/validation.py:71`; `docs/spec/feature.md:145-147`; `docs/spec/converter-coverage.md:199-200`; `tests/converters/test_identification_origin.py:131-144`

### FRAGPIPE-15: identification_type and id_run_file_name ambiguous for multi-run units, TMT and missing PSMs {#fragpipe-15}

*feature* · ambiguous · confirmed · checklist D1, F4, F10

- **Now:** In a multi-run experiment, Match Type MS/MS gives 'direct' while run_file_name is an arbitrary fraction and id_run null. A one-run 'direct' row gets id_run_file_name=run even when no PSM is found, leaving scan [] and PEP null (feature_adapter.py:340; verified experiment D). TMT plexes and inputs without Match Type get null. A null PEP has about six causes (MBR, multi-run, no psm.tsv, missing PSM, mass-shifted PSM, run-name mismatch).
- **Spec:** feature.md:145: 'How the tool identified the feature in THIS run'; id_run_file_name 'can be null for other reasons'. feature.md:138: id_run_file_name is the run of the best PSM. There is no rule for unit-level rows.
- **Reader impact:** A reader cannot tell whether 'direct' refers to run_file_name. id_run_file_name may name a run with no linked PSM, and null PEP or scan cannot be interpreted.
- **Suggested:** Set id_run_file_name only when a PSM is linked (#378). For unit-level rows, define 'direct' as 'MS/MS in some member run'. Consider a reason code or a documented null semantics.
- **Evidence:** `qpx/converters/fragpipe/feature_adapter.py:335-341`; `qpx/converters/fragpipe/feature_adapter.py:172-174`; `docs/spec/feature.md:138`; `docs/spec/feature.md:145`; `docs/spec/converter-coverage.md:203`; `run: out/fp.feature.parquet unit D id_run_file_name='run.D.v2', scan=[]`

### FRAGPIPE-16: Identified-but-unquantified ions are dropped instead of written with null intensity {#fragpipe-16}

*feature* · diverges from spec · confirmed · checklist D6, F6, F3

- **Now:** LFQ: '<exp> Intensity' of 0 or missing returns None and no row is written (feature_adapter.py:332-333,376-377). Verified: ANTHERSK in experiment A (Match Type MS/MS, spectral count 1, PSM in psm.tsv) has no feature. TMT: a precursor whose channels are all 0 writes no feature (isobaric.py:255-261; test_fragpipe_isobaric.py:88-90). For TMT, an individual 0 channel becomes a null entry.
- **Spec:** producer-guide.md:45-46: 'Not identified -> no entry; Identified, not quantified -> entry with null intensity'. intensities.md:123-127: a 0 sentinel is written as null, and a feature with no entry has a null intensities list.
- **Reader impact:** FragPipe evidence of 'identified in this experiment but not quantified' vanishes from the feature view, while MaxQuant keeps it. Missingness analyses (MNAR vs not identified) are biased by producer.
- **Suggested:** Emit a row with a null intensity when Match Type is MS/MS (or spectral count > 0) and the intensity is 0. Make this rule mandatory in the contract.
- **Evidence:** `qpx/converters/fragpipe/feature_adapter.py:332-333`; `qpx/converters/fragpipe/feature_adapter.py:376-377`; `qpx/converters/fragpipe/isobaric.py:255-261`; `tests/converters/test_fragpipe_isobaric.py:88-90`; `docs/spec/producer-guide.md:43-48`; `docs/spec/intensities.md:123-127`; `run: out/fp.feature.parquet has no ANTHERSK row for unit A`

### FRAGPIPE-17: Feature spectral count and RT not read; intensity meaning undefined {#fragpipe-17}

*feature* · inconsistent across producers · confirmed · checklist F6, F10, X1

- **Now:** Only '<exp> Intensity' is read (feature_adapter.py:376). '<exp> Spectral Count' is ignored, so there is no spectral_count in additional_intensities, and the per-experiment apex RT column is ignored, so rt, rt_start and rt_stop are null. That combined_ion has these columns in every FragPipe version is assumed from IonQuant output and was not verified on a real file. Whether IonQuant's ion intensity is an apex height or a peak area is not stated anywhere.
- **Spec:** intensities.md:187 lists spectral_count for MaxQuant and PD features only. feature.md:65 defines rt as the precursor retention time. converter-coverage.md:189 says FragPipe rt is 'Not implemented'. Nothing defines what the primary feature intensity measures.
- **Reader impact:** Spectral-count-based filters work for MaxQuant features but not FragPipe ones. RT-based QC or alignment is impossible, and comparing intensities across producers mixes apex and area.
- **Suggested:** Read '<exp> Spectral Count' as spectral_count and the apex RT as rt. Have the contract require a declaration of intensity type (apex, area or sum) per view or footer.
- **Evidence:** `qpx/converters/fragpipe/feature_adapter.py:353-377`; `docs/spec/intensities.md:187`; `docs/spec/feature.md:65`; `docs/spec/converter-coverage.md:189`

### FRAGPIPE-18: Feature unique comes from mapped-protein count, not group size; pg_accessions depend on --pg-file {#fragpipe-18}

*feature* · ambiguous · partly confirmed (see correction) · checklist F9

- **Now:** anchor_protein is the ion's razor 'Protein'. pg_accessions is the combined_protein group led by it, or only the lead when --pg-file is absent (feature_adapter.py:105-116,529). unique is len(Protein + Mapped Proteins) == 1, and null without the Mapped Proteins column (feature_adapter.py:506,531). Synthetic inconsistent input gave unique=true alongside 2-member pg_accessions.
- **Spec:** feature.md:99: pg_accessions = full group membership, equal to pg.pg_accessions; for a shared peptide, the assigned (razor) group. feature.md:101: unique = true when the group has exactly one protein; false for a multi-member or shared-assigned group.
- **Correction:** The code facts are correct. The conflict between the two 'unique' definitions only appears when Mapped Proteins leaves out indistinguishable members, which the synthetic input forces; in real FragPipe output those members share every peptide, so they should be listed (not verified). Lead-only pg_accessions without --pg-file is documented (converter-coverage.md:217-219); it violates feature.md:99 only against a pg view of the same groups.
- **Reader impact:** The two definitions agree only when Mapped Proteins includes every indistinguishable member. Without --pg-file, pg_accessions does not match pg rows, and the feature-pg link silently fails.
- **Suggested:** The contract should fix one definition of unique (group-size or peptide-to-protein specificity), and require pg_accessions to match pg or be null.
- **Evidence:** `qpx/converters/fragpipe/feature_adapter.py:105-116`; `qpx/converters/fragpipe/feature_adapter.py:501-531`; `docs/spec/feature.md:99-101`; `docs/spec/converter-coverage.md:213-221`; `tests/converters/test_fragpipe_converter.py:464-501`

### FRAGPIPE-19: TMT plex features: unit-level reporter sums with no PSM, PIF or quant scan {#fragpipe-19}

*feature* · unspecified · partly confirmed (see correction) · checklist D7, F6, F10

- **Now:** Per plex, intensities are the ion.tsv channel columns (a 0 becomes null). MS1 'Intensity' is stored as ms1_intensity under the first quantified channel only (isobaric.py:255-270). For a fractionated plex, no PSM enrichment happens: run_file_name is the first sorted fraction, the unit is the output prefix, and identification_type is null. parent_ion_fraction, quant_scan and quant_ms_level are null, although psm.tsv has Purity per PSM. Verified with a synthetic plex.
- **Spec:** feature.md:154-156 define parent_ion_fraction, quant_scan and quant_ms_level. intensities.md:189 documents ms1_intensity. Nothing says whether ion.tsv channel values are raw, summed or normalised reporter intensities.
- **Correction:** Placing ms1_intensity under the first quantified channel is a specified convention (intensities.md:189), and PIF only on the psm is documented (converter-coverage.md:202). The spec has a mechanism for reporter processing state, the isotope_corrected and normalization provenance keys (intensities.md:211-215); FragPipe does not fill it. A one-run plex does get PSM enrichment. Reporter aggregation semantics remain unspecified.
- **Reader impact:** A reader does not know how the reporter values were aggregated or normalised. Ratio-compression QC (PIF) is impossible at feature level, and the MS1 intensity sits under an arbitrary channel.
- **Suggested:** The contract should define reporter intensity semantics (raw, corrected or normalised; how PSMs are aggregated) and where an unlabelled MS1 quantity lives in a multiplexed row.
- **Evidence:** `qpx/converters/fragpipe/isobaric.py:246-270`; `qpx/converters/fragpipe/converter.py:56-59`; `qpx/converters/fragpipe/feature_adapter.py:340`; `docs/spec/feature.md:154-156`; `docs/spec/intensities.md:189`; `run: audit-run:fragpipe/tmt/out/plex1.feature.parquet scan=[], identification_type=None, parent_ion_fraction=None`

## Proteome Discoverer converter

### PD-04: Same spectrum from two engines summed twice in feature reporters and spectral_count {#pd-04}

*feature* · ambiguous · confirmed · checklist F2, F6, P6

- **Now:** The psm view collapses Sequest/Comet duplicates per scan and identity (psm_adapter.py:38-40). The feature SQL sums over all rank-1 rows with no per-scan dedup (feature_adapter.py:66,70). Run: one scan with TMT126=10, found by two nodes, gave feature TMT126=20 but one psm row. LFQ spectral_count counts rows, not spectra (feature_adapter.py:84). Not verified whether PD itself marks such duplicates in Quan Info.
- **Spec:** converter-coverage.md:324: 'Reporter channels are summed over them' (PSM rows). intensities.md:187: spectral_count = 'the precursor's PSMs'. Whether the unit of aggregation is the spectrum or the PSM row is unspecified.
- **Reader impact:** Multi-engine PD results inflate reporter sums and spectral counts by up to the number of engines. The feature view disagrees with the psm view of the same file.
- **Suggested:** Aggregate over distinct (run, scan) per feature. Define spectral_count as the number of distinct MS2 spectra.
- **Evidence:** `qpx/converters/proteome_discoverer/psm_adapter.py:38-40`; `qpx/converters/proteome_discoverer/feature_adapter.py:66`; `qpx/converters/proteome_discoverer/feature_adapter.py:70`; `qpx/converters/proteome_discoverer/feature_adapter.py:84`; `docs/spec/intensities.md:187`; `run: audit-run:pd/cases/two_engines_same_scan/out/pd.feature.parquet`

### PD-10: Feature precursor_qvalue holds the lowest PSM-level q-value; the spec permits this {#pd-10}

*feature* · contradictory · confirmed · checklist D4, F5, P5

- **Now:** feature_adapter.py:125 copies the chosen PSM's Percolator q-Value into precursor_qvalue. peptide_qvalue and global_qvalue stay null. psm.qvalue is written directly (psm_adapter.py:63), and no q-value appears in additional_scores (reports.py:138-147). feature.confident() filters on COALESCE(precursor_qvalue, peptide_qvalue) (qpx/core/data/feature.py:54).
- **Spec:** feature.md:58: precursor_qvalue = 'Run-level precursor q-value'. scores.md:44 lists qvalue and psm_level_qvalue as sources for it, which allows PSM-level values. scores.md:38 says the flat columns are 'filled from additional_scores'.
- **Reader impact:** The best PSM's q-value is optimistic compared with a precursor-level FDR. DIA-NN's and PD's precursor_qvalue are not comparable.
- **Suggested:** Remove PSM-level names from the precursor_qvalue sources. Leave the column null, or add a best_psm_qvalue column. Read PD PeptideGroups q-values into peptide_qvalue.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:125`; `qpx/converters/proteome_discoverer/psm_adapter.py:63`; `qpx/writers/derived.py:30-34`; `qpx/core/data/feature.py:54`; `docs/spec/feature.md:58`; `docs/spec/scores.md:38-47`

### PD-12: Identified-but-unquantified precursors get no feature row {#pd-12}

*feature* · diverges from spec · partly confirmed (see correction) · checklist D6, F3

- **Now:** A feature row is dropped in three cases. Every channel sum is null or ≤0 (feature_adapter.py:89-92,99-100). The LFQ abundance is null or 0 (:80-83). All its PSMs are Rejected or 'Excluded by Method' (WHERE clause, :18,72-76). Run: 3 identified precursors gave 1 feature row, while the psm view keeps all 3. Partially quantified channels correctly get null entries.
- **Spec:** producer-guide.md:45-46 and intensities.md:123-125: identified but not quantified is an entry with null intensity; no entry means not identified.
- **Correction:** Two drop cases diverge from the spec: precursors whose PSMs are all 'Excluded by…', and those with all channels empty or ≤0. The third case is defensible. PD's ambiguity resolution rejected the 'Rejected' PSMs, so a precursor seen only in Rejected PSMs is not 'identified'.
- **Reader impact:** A missing feature row is read as 'not identified', so missing-value models (MNAR vs MAR) get wrong input.
- **Suggested:** Emit the row with null intensities for every plex channel, and record the exclusion reason as a cv_param.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:18`; `qpx/converters/proteome_discoverer/feature_adapter.py:72-76`; `qpx/converters/proteome_discoverer/feature_adapter.py:80-92`; `qpx/converters/proteome_discoverer/feature_adapter.py:99-100`; `docs/spec/producer-guide.md:45-46`; `docs/spec/intensities.md:123-125`; `run: audit-run:pd/cases/excluded_and_empty`

### PD-18: Feature LFQ intensity mixes Minora abundance and MS2 precursor 'Intensity', taking the max over PSMs {#pd-18}

*feature* · ambiguous · confirmed · checklist F6, F1

- **Now:** Label-free: the intensity is the MAX over PSMs of 'Precursor Abundance', else 'Intensity' (column_mappings.yaml:69; feature_adapter.py:69,85), so several Minora peaks of one precursor collapse to the largest. TMT: the SUM of PSM reporter abundances (feature_adapter.py:66). Whether these are S/N or intensity, and isotope-corrected or not, depends on PD settings and is not recorded. Run: LFQ=300 taken from 'Intensity'.
- **Spec:** intensities.md:19: raw intensity. intensities.md:213: processing state via provenance keys. feature.md:84: 'primary intensity'. The spec has no vocabulary for apex vs area vs S/N.
- **Reader impact:** LFQ values from PD exports with and without Precursor Abundance differ in kind. TMT S/N sums get compared with other tools' intensities.
- **Suggested:** Contract: a footer or provenance key intensity_kind (ms1_area, ms1_apex, ms2_precursor, reporter_sn, reporter_intensity), set by PD from the column used.
- **Evidence:** `qpx/converters/column_mappings.yaml:69`; `qpx/converters/proteome_discoverer/feature_adapter.py:66`; `qpx/converters/proteome_discoverer/feature_adapter.py:69`; `qpx/converters/proteome_discoverer/feature_adapter.py:85`; `docs/spec/intensities.md:19`; `docs/spec/intensities.md:213`; `run: audit-run:pd/cases/labelfree_intensity_fallback`

### PD-19: Feature parent_ion_fraction comes from one PSM while reporters are summed over many scans {#pd-19}

*feature* · ambiguous · confirmed · checklist D7, F6

- **Now:** feature_adapter.py:114-116,140 take 1 − Isolation Interference of the chosen PSM only. quant_scan and quant_ms_level are never set. Average Reporter S/N and SPS Mass Matches are not read. Co-isolation and S/N thresholds act only through Quan Info 'Excluded by…' (feature_adapter.py:73-76); every other Quan Info value passes.
- **Spec:** feature.md:154-156: parent_ion_fraction, quant_scan and quant_ms_level describe a single quantification scan. No rule covers features summed over several PSMs.
- **Reader impact:** Interference filters on features judge a sum by one of its scans, and the SPS-MS3 origin is lost.
- **Suggested:** Contract: for multi-scan features, parent_ion_fraction is the intensity-weighted mean or null, and quant_scan is null. Record quan filters in provenance.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:73-76`; `qpx/converters/proteome_discoverer/feature_adapter.py:114-116`; `qpx/converters/proteome_discoverer/feature_adapter.py:140`; `docs/spec/feature.md:154-156`

### PD-20: Shared peptide: anchor is the first listed master; without Proteins.txt pg_accessions merges all masters {#pd-20}

*feature* · diverges from spec · partly confirmed (see correction) · checklist F9, P8

- **Now:** anchor_protein is the first entry of 'Master Protein Accessions' (feature_adapter.py:103-105), not PD's razor/quan assignment. pg_accessions is that master's group from Proteins.txt, else all masters (:112). Run: 'P55555; P12345' gave anchor P55555. Without Proteins.txt, pg_accessions=[P55555, P12345], a pseudo-group that matches no pg row.
- **Spec:** feature.md:99: for a peptide shared by several groups, the group the producer assigned it to. pg.md:70-74: the softlink needs equal membership. #322: an unresolved peptide gets null.
- **Correction:** Without Proteins.txt no pg view is written at all, so the problem is the membership required by feature.md:99, not a softlink that fails. With Proteins.txt, pg_accessions is the first master's group, which matches a pg row but is an arbitrary choice among the shared groups.
- **Reader impact:** Shared peptides attach to an arbitrary group, so the feature→pg softlink fails or points at a group PD did not quantify with.
- **Suggested:** Use PD's razor assignment, or set anchor and pg_accessions to null when unresolved (#322).
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:103-105`; `qpx/converters/proteome_discoverer/feature_adapter.py:112`; `docs/spec/feature.md:99`; `docs/spec/pg.md:70-74`; `run: audit-run:pd/cases/shared_two_masters_noprot, shared_two_masters_prot`

### PD-21: unique = peptide maps to one protein, even in a multi-member group; a test pins it {#pd-21}

*feature* · contradictory · confirmed · checklist F9

- **Now:** feature_adapter.py:136 sets unique = (len('Protein Accessions') == 1), whatever the group size. A test asserts unique=True with pg_accessions [P12345, P99999] (tests/converters/test_proteome_discoverer_converter.py:166-168). pepmap is_unique uses the same rule (qpx/converters/base.py:133).
- **Spec:** feature.md:101 and versioning.md:68-69: true only for a one-protein group that maps to no other group; false when the group has several members. pepmap.md:21 defines is_unique as 'maps to exactly one protein', a different notion.
- **Reader impact:** Unique-only quantification differs from other producers. 'unique' means protein-unique in pepmap but group-unique in feature.
- **Suggested:** Compute unique = (|pg_accessions| == 1 and a single master), and fix the test. Rename the pepmap flag or align the definitions.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:136`; `tests/converters/test_proteome_discoverer_converter.py:166-168`; `qpx/converters/base.py:133`; `docs/spec/feature.md:101`; `docs/spec/versioning.md:68-69`; `docs/spec/pepmap.md:21`

### PD-29: Feature = PSM aggregate per run, peptide and charge; rt means different things for TMT and LFQ {#pd-29}

*feature* · ambiguous · partly confirmed (see correction) · checklist F1, F7, F10

- **Now:** One row per (run, upper(Annotated Sequence), Modifications string, charge) (feature_adapter.py:16-20), not a chromatographic peak. rt is the chosen PSM's MS2 RT, replaced by 'Apex RT' when present (:109-111). rt_start and rt_stop are null (:143-144). Grouping uses raw strings while ids hash the peptidoform: 'N-Term(Prot)(Acetyl)' and 'N-Term(Acetyl)' gave two rows with the same peptidoform, charge and run (verified).
- **Spec:** #217: a feature is a chromatographic peak. feature.md:22-23: default composite includes rt. feature.md:65: rt = 'Precursor retention time'.
- **Correction:** The granularity is specified: producer-guide.md:12 says a feature is 'One row per precursor, charge and run', so PD follows the spec. #217's 'peak' is not spec text. The raw-string grouping is worse than reported. With equal rt (LFQ Apex RT) the two rows share one feature_id, and the psm dedup also uses the raw text, giving duplicate psm_ids. The writer only warns.
- **Reader impact:** TMT feature rt is an MS2 time, LFQ rt an apex, under the same composite. Duplicate precursor rows are possible.
- **Suggested:** Contract: rt = apex of the quantified signal, else null. Group on the canonical peptidoform. Declare the feature granularity in the footer.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:16-20`; `qpx/converters/proteome_discoverer/feature_adapter.py:109-111`; `qpx/converters/proteome_discoverer/feature_adapter.py:143-144`; `docs/spec/feature.md:22-23`; `run: audit-run:pd/cases/mod_spelling_variants`

### PD-M1: Feature id composite includes the picked PSM's rt though rows are per precursor/run {#pd-m1}

*feature* · diverges from spec · found by the verifier · checklist F1, D8

- **Now:** Footer identity_composite is the schema default peptidoform,charge,run_file_name,rt (verified), while rows are one per (run, raw sequence, raw mods, charge) (feature_adapter.py:16-20). rt is the winning PSM's MS2 RT or Apex RT (feature_adapter.py:108-111; reports.py:119), so feature_id depends on the best-PSM rule. Raw-text duplicates of one precursor get distinct ids (case modvar).
- **Spec:** producer-guide.md:56: a tool with one row per precursor and run declares peptidoform, charge, run_file_name, as DIA-NN does. feature.md:11-12,22-29: converters may override the default.
- **Reader impact:** feature_id changes whenever the best-PSM pick changes (e.g. #378's lowest-PEP rule), and duplicate precursor rows get distinct ids, so the PK check cannot catch them.
- **Suggested:** PD declares [peptidoform, charge, run_file_name] and groups on the canonical peptidoform. Contract: the composite must name the row's granularity, not a by-product of picking the best PSM.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:16-20`; `qpx/converters/proteome_discoverer/feature_adapter.py:108-111`; `qpx/core/data/schemas/feature.yaml:4`; `docs/spec/producer-guide.md:56`; `qpx/converters/diann/feature_adapter.py:168`; `run: work/verify-pd/cases/modvar, featcollide`

### PD-M2: feature.posterior_error_probability: spec says both 'not carried' and 'peptide match PEP' {#pd-m2}

*feature* · contradictory · found by the verifier · checklist F5, D4

- **Now:** PD writes the PEP of the lowest-q-value PSM (feature_adapter.py:17,117-118 via reports.py:120). This is neither PD's peptide-group PEP nor the lowest PEP: in case 'tie' the feature has PEP 0.05 although the same spectrum had a 0.001 hit.
- **Spec:** psm.md:285-286: PEP 'is defined only in the PSM view… The feature and peptide views do not carry PEP as a top-level field'. feature.md:56 and feature.yaml:37-41 define feature.posterior_error_probability as 'PEP for the peptide match'.
- **Reader impact:** Readers cannot tell whether feature PEP is a peptide-level PEP or one PSM's PEP. PD values come from a q-value-chosen PSM, so they are not comparable with other producers.
- **Suggested:** Define it as the identifying PSM's PEP (#378 rule 1) and drop psm.md:286's last sentence. A peptide-level PEP goes in additional_scores under a level-named score.
- **Evidence:** `docs/spec/psm.md:285-286`; `docs/spec/feature.md:56`; `qpx/core/data/schemas/feature.yaml:37-41`; `qpx/converters/proteome_discoverer/feature_adapter.py:17`; `run: work/verify-pd/cases/tie`

### PD-M4: TMT features carry no count of summed PSMs; LFQ features do {#pd-m4}

*feature* · unspecified · found by the verifier · checklist F2, F6, X1

- **Now:** _psm_count (feature_adapter.py:70) is written as spectral_count only in the label-free branch (:84). The reporter branch (:86-95) drops it, although reporter channels are summed over every used PSM (:66). A 1-PSM sum cannot be told from a 20-PSM sum.
- **Spec:** intensities.md:187 defines spectral_count for PD as 'the precursor's PSMs' without restricting it to LFQ. No rule says how a feature declares the number of spectra its quantity aggregates.
- **Reader impact:** Summed reporter intensities look like single-scan values from other tools, and readers cannot normalise or weight them by spectrum count.
- **Suggested:** Contract: any feature quantity aggregated over several spectra states the count of distinct spectra (spectral_count) for every label type.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:66`; `qpx/converters/proteome_discoverer/feature_adapter.py:70`; `qpx/converters/proteome_discoverer/feature_adapter.py:84-95`; `docs/spec/intensities.md:187`

### PD-M5: Every reporter column becomes an entry, whether or not the channel has a sample {#pd-m5}

*feature* · ambiguous · found by the verifier · checklist D6, F6

- **Now:** feature_adapter.py:86-89 writes an intensities entry for every 'Abundance <channel>' column of PSMs.txt, null when empty or ≤0. It does not consult the SDRF labels, which go only to the sample/run converter (converter.py:66-67). If PD exports an unused channel as an empty column (unverified), every feature gets an 'identified, not quantified' entry for it.
- **Spec:** intensities.md:123-126: a null entry marks a label the feature was identified in but not quantified. 'A label with no entry was not identified or has no sample in the run.'
- **Reader impact:** Missing-value models and per-channel completeness would see spurious missing values for channels with no sample.
- **Suggested:** Contract: entries only for labels that have a sample in the run (run.samples or SDRF). Converters given an SDRF limit channels to its labels.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:55-62`; `qpx/converters/proteome_discoverer/feature_adapter.py:86-89`; `qpx/converters/proteome_discoverer/converter.py:66-67`; `docs/spec/intensities.md:123-126`

## mzIdentML and CDAP converters

### MZID_CDAP-13: CDAP feature precursor/peptide q-values are best PSM's MS-GF+ values; docs disagree {#mzid_cdap-13}

*feature* · contradictory · confirmed · checklist D4, F5

- **Now:** Feature scores come from the best PSM: Qvalue as psm_level_qvalue, PepQvalue as peptide_level_qvalue (cdap/feature_adapter.py:289-299); writer copies them to precursor_qvalue and peptide_qvalue (writers/derived.py:31-32). Fixture: 94/94 features have both filled. PSM rows also carry peptide_level_qvalue (cdap/psm_adapter.py:131-137).
- **Spec:** scores.md:44 maps psm_level_qvalue to feature.precursor_qvalue ('precursor, per run'); converter-coverage.md:181 says CDAP peptide_qvalue 'Not implemented'.
- **Reader impact:** precursor_qvalue holds one PSM's PSM-level q (scope of MS-GF+ FDR unknown, not per run): level conflation sanctioned by spec; coverage doc misleads.
- **Suggested:** Copy a PSM-level q to precursor_qvalue only when the producer states precursor level; fix converter-coverage.md:181.
- **Evidence:** `qpx/converters/cdap/feature_adapter.py:289-299`; `qpx/writers/derived.py:30-34`; `audit-run:mzid_cdap/cdap/cdap.feature.parquet count(peptide_qvalue)=94`; `docs/spec/scores.md:44`; `docs/spec/converter-coverage.md:181`

### MZID_CDAP-14: CDAP feature row is a PSM rollup keyed on raw sequence, not a peak {#mzid_cdap-14}

*feature* · ambiguous · partly confirmed (see correction) · checklist F1, F7, P6

- **Now:** GROUP BY raw PeptideSequence, QueryCharge, FileName (cdap/feature_adapter.py:116-131); identity composite peptidoform,charge,run_file_name (feature_adapter.py:35); rt_start/rt_stop null (225-226). Grouping key (raw CDAP string) differs from identity key (canonical peptidoform): two raw spellings of one peptidoform (e.g. '+42.011+229.163' vs '+229.163+42.011' N-term) would give duplicate feature_id. Not seen in fixture.
- **Spec:** data-model.md:63-64 feature = one quantified chromatographic peak; identity.md:31 allows CDAP composite; producer-guide.md:12 one row per precursor, charge, run.
- **Correction:** The duplicate-id example is wrong. Raw order variants become different strings ('+42.011+229.163PEPTK' -> '[UNIMOD:1][UNIMOD:737]-PEPTK', the reverse -> '[UNIMOD:737][UNIMOD:1]-PEPTK'), so they get distinct ids, not duplicates. I found no realistic CDAP collision. The definitional ambiguity is real: data-model.md:64 says 'chromatographic peak', data-model.md:74 says 'per upstream Feature entity', and feature.md:301 says a feature 'aggregates PSMs'.
- **Reader impact:** Consumers assuming an MS1 peak (rt window, apex) get an MS2 reporter-ion sum with no peak bounds.
- **Suggested:** Contract: feature = precursor x charge x run quantity entity, with a footer key declaring its kind (peak | psm_rollup); producers group on the canonical peptidoform.
- **Evidence:** `qpx/converters/cdap/feature_adapter.py:35,116-131,225-226`; `docs/spec/data-model.md:63-64`; `docs/spec/identity.md:31`

### MZID_CDAP-15: CDAP best-PSM rule and intensity aggregation differ between TMT and LFQ {#mzid_cdap-15}

*feature* · unspecified · confirmed · checklist F2, F5, F6, F10

- **Now:** Best PSM = lowest SpecEValue, then scan, m/z, rt... (cdap/feature_adapter.py:341-372); scan, rt, observed_mz, scores come from it. TMT channels are SUMMED over all PSMs (feature_adapter.py:377-379); LFQ intensity = PrecursorArea of the best PSM only (284-286,375). PSM count only as cv_param 'psm_count' string (385-388).
- **Spec:** feature.md:138-139 'best PSM' undefined; #378 (open) proposes tool pick else lowest PEP; intensities.md says nothing about summing vs picking.
- **Reader impact:** One column means 'sum over N spectra' (TMT) or 'one PSM's MS1 area' (LFQ); scan names one of several contributing spectra.
- **Suggested:** Contract fixes the best-PSM fallback order and requires producers to declare intensity aggregation (sum/max/apex/area) per view.
- **Evidence:** `qpx/converters/cdap/feature_adapter.py:280-287,332-388`; `docs/spec/feature.md:138-139`

### MZID_CDAP-16: CDAP drops tied (AmbiguousMatch) PSMs from features but keeps them in psm {#mzid_cdap-16}

*feature* · unspecified · partly confirmed (see correction) · checklist F2, G1

- **Now:** PSMs with AmbiguousMatch=1 excluded from feature aggregation (cdap/feature_adapter.py:116-120). Fixture: 6 PSMs/3 spectra (incl. I/L pair SCHTGIGR/SCHTGLGR and a decoy) are in psm but in no feature (100 PSMs -> 94 features). A precursor identified only via tied PSMs gets no feature and no pg evidence.
- **Spec:** Nothing on features with conflicting identifications.
- **Correction:** For features the spec is silent, but not for pg. pg.md:126-131 requires a pg row when a passing PSM maps to the group, 'assigned to a feature or unassigned'. CDAP builds pg only from features, so target NP_001054.1 (I/L-tied PSMs at scans 5646/5652, Qvalue 0.0003) has no pg row: this diverges from the spec.
- **Reader impact:** psm and feature/pg disagree on what was identified; I/L-ambiguous peptides silently vanish from quantification.
- **Suggested:** Contract states a policy for tied IDs: drop with a recorded reason, one row per candidate, or merge I/L as J.
- **Evidence:** `qpx/converters/cdap/feature_adapter.py:116-120`; `audit-run:mzid_cdap/cdap outputs (100 psm vs 94 features, psm_count=1 each)`

### MZID_CDAP-22: CDAP writes identification_type 'direct' and id_run_file_name=run for all features {#mzid_cdap-22}

*feature* · specified, worth keeping · confirmed · checklist D1, D2, F4

- **Now:** identification_type='direct', id_run_file_name=run_file_name for every feature (cdap/feature_adapter.py:223-224); match_* null. Fixture 94/94 'direct'. Validator enum is direct/transferred/secondary (core/data/validation.py:71), no 'unidentified'.
- **Spec:** feature.md:145 direct = its own MS/MS; null when the tool does not say.
- **Reader impact:** Correct: CDAP has no MBR, every feature is backed by MS2 in the same run.
- **Suggested:** Keep: a tool without ID transfer writes 'direct' (not null); null only when transfer status is unknown.
- **Evidence:** `qpx/converters/cdap/feature_adapter.py:223-224`; `qpx/core/data/validation.py:71`; `docs/spec/feature.md:145`

### MZID_CDAP-23: CDAP unquantified features keep null entries; coverage doc says null list {#mzid_cdap-23}

*feature* · contradictory · confirmed · checklist D6, F6, G1

- **Now:** CDAP reporter 0 (incl. '0/?' cells) -> entry with null intensity (cdap/feature_adapter.py:280-287); a feature with no positive channel keeps 10 null entries, not a null list (fixture 6/94). Text after '/' in reporter cells and TMT10-Flags/FractionOfTotalAb/TotalAb are dropped (feature_adapter.py:325-329). pg: all-null channels -> row with null intensity (pg_adapter.py:203-207).
- **Spec:** intensities.md:123-127 null entry = identified, not quantified; converter-coverage.md:194 says CDAP intensities are 'Null when no reporter channel ... is positive'.
- **Reader impact:** Consumers testing `intensities IS NULL` (per coverage doc) miss unquantified features; meaning of dropped '/x' and Flags (interference?) is lost.
- **Suggested:** Keep the intensities.md rule, fix converter-coverage.md:194, document what CDAP '0', '/x' and Flags mean.
- **Evidence:** `qpx/converters/cdap/feature_adapter.py:280-287,325-329`; `audit-run:mzid_cdap/cdap/cdap.feature.parquet (6 all-null features)`; `docs/spec/intensities.md:123-127`; `docs/spec/converter-coverage.md:194`

### MZID_CDAP-24: CDAP reporter intensities: impurity correction and normalization state undeclared {#mzid_cdap-24}

*feature* · unspecified · confirmed · checklist F6, G2, D7

- **Now:** Reporter values taken as-is (part before '/') and summed (cdap/feature_adapter.py:325-329,377-379); provenance has no isotope_corrected/normalization keys (cdap/converter.py:191-201). I could not determine from code or fixtures whether CDAP reporter ions are impurity-corrected.
- **Spec:** intensities.md:213 state isotope_corrected and normalization on the tool's provenance step when known; no key = input does not say.
- **Reader impact:** CDAP TMT intensities get compared with MaxQuant corrected or PD normalized values without knowing the processing.
- **Suggested:** Contract: isobaric producers MUST state isotope_corrected (true/false/unknown); CDAP declares it from CDAP documentation.
- **Evidence:** `qpx/converters/cdap/feature_adapter.py:325-329,377-379`; `qpx/converters/cdap/converter.py:191-201`; `docs/spec/intensities.md:213`

## Shared writer and validation code

### INFRA-01: Writer copies PSM-level q-value names into feature.precursor_qvalue {#infra-01}

*feature* · contradictory · confirmed · checklist D4, F5

- **Now:** FEATURE_SCORE_COLUMNS maps precursor_qvalue <- {precursor_qvalue, qvalue, q_value, psm_level_qvalue} (qpx/writers/derived.py:32), coalesced at write (qpx/writers/feature.py:46-47, derived.py:58-71). CDAP puts the lowest-E-value PSM's MS-GF+ Qvalue in feature additional_scores as psm_level_qvalue (qpx/converters/cdap/feature_adapter.py:293,345-369). Ran CDAP fixture: 94/94 feature.precursor_qvalue == psm_level_qvalue of the best PSM.
- **Spec:** scores.md:44 lists qvalue, q_value, psm_level_qvalue as sources of feature.precursor_qvalue, whose scope it calls 'Precursor, per run'. scores.md:42 treats the same names as PSM-level for psm.qvalue. feature.md:58 says 'Run-level precursor q-value'.
- **Reader impact:** feature.precursor_qvalue and feature.confident() treat one PSM's PSM-level q-value as a precursor-level FDR. For DIA it is a real precursor q-value. Same column, different statistic.
- **Suggested:** Fill precursor_qvalue only from names that say they are precursor-level (precursor_qvalue). Keep a best-PSM q-value in additional_scores under a name like best_psm_qvalue. Remove qvalue, q_value and psm_level_qvalue from scores.md:44.
- **Evidence:** `qpx/writers/derived.py:28-34`; `qpx/writers/derived.py:58-71`; `qpx/converters/cdap/feature_adapter.py:289-295`; `qpx/converters/cdap/feature_adapter.py:342-369`; `docs/spec/scores.md:42-45`; `run: qpxc convert cdap on tests/examples/cdap -> audit-run:infra/cdap; duckdb: 94/94 precursor_qvalue == psm_level_qvalue`

### INFRA-02: feature.precursor_qvalue holds a different statistic for each producer {#infra-02}

*feature* · inconsistent across producers · confirmed · checklist D4, F5

- **Now:** DIA-NN and Spectronaut: a true run-level precursor q-value (diann/feature_sql.py:14; spectronaut 'qvalue' from EG.Qvalue, feature_adapter.py:39). CDAP: the best PSM's PSM q-value (INFRA-01). Proteome Discoverer: the lowest Percolator q-Value among rank-1 PSMs, set directly (proteome_discoverer/feature_adapter.py:17,125; column_mappings.yaml:58). OpenMS-consensus, MaxQuant, FragPipe, quantms-MSstats: null (verified for OpenMS and quantms runs).
- **Spec:** scores.md:44: 'Precursor, per run'. No rule says how a DDA producer turns PSM q-values into a precursor q-value. Taking the minimum over a precursor's PSMs is not forbidden.
- **Reader impact:** Filtering features at precursor_qvalue<0.01 across datasets keeps different populations. In DDA it is a best-of-N PSM statistic, which is anti-conservative.
- **Suggested:** Contract: precursor_qvalue = a precursor-level FDR computed by the tool. A best-PSM q-value is never promoted to it. Producers without a precursor-level FDR leave it null.
- **Evidence:** `qpx/converters/proteome_discoverer/feature_adapter.py:16-20`; `qpx/converters/proteome_discoverer/feature_adapter.py:125`; `qpx/converters/column_mappings.yaml:58`; `qpx/converters/diann/feature_sql.py:14`; `qpx/converters/spectronaut/feature_adapter.py:39`; `run: openms-consensus ProteomicsLFQ fixture: precursor_qvalue 0/68 filled`

### INFRA-03: OpenMS feature.peptide_qvalue is the PSM-level q-value of the first identification {#infra-03}

*feature* · diverges from spec · confirmed · checklist D4, F5

- **Now:** qvalue_of() returns hit.getScore() whenever the PID score type is 'q-value'/'qvalue'/'fdr', whatever the FDR level (openms_consensus/feature_adapter.py:149-160). _confidence_by_run takes the first PID's first hit for each run (:325-344) and it is written as peptide_qvalue (:736-745). The same score is psm.qvalue (psm_adapter.py:216-219). ProteomicsLFQ fixture: 68/68 feature.peptide_qvalue equal the linked psm.qvalue (0.00862). The consensusXML says nothing about peptide-level FDR.
- **Spec:** feature.md:126: peptide_qvalue is 'the canonical peptide-level q-value'. feature.md:132 only requires the producer to state it is a q-value, not its level. converter-coverage.md:59 says OpenMS reads any 'score explicitly declared to be a q-value/FDR'.
- **Reader impact:** Peptide-level FDR filters on OpenMS/quantms output actually filter on one PSM's PSM-level q-value. The same number is shown as PSM-level in psm and peptide-level in feature.
- **Suggested:** Write peptide_qvalue only when the producer declares peptide-level FDR (the OpenMS peptide-level FDR parameter). Otherwise keep it null and keep the value as PSM-level. Add 'level declared' to the feature.md:132 rule.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:149-160`; `qpx/converters/openms_consensus/feature_adapter.py:325-344`; `qpx/converters/openms_consensus/feature_adapter.py:736-745`; `docs/spec/feature.md:126-132`; `docs/spec/converter-coverage.md:59`; `run: qpxc convert openms-consensus ProteomicsLFQ_1_out.consensusXML; join psm.feature_id: 68/68 equal, PID score type {'q-value'}`

### INFRA-04: Coverage tables contradict write-time q-value derivation (CDAP, Spectronaut) {#infra-04}

*feature* · contradictory · partly confirmed (see correction) · checklist D4

- **Now:** The writer fills feature.peptide_qvalue from 'peptide_level_qvalue' (derived.py:31) and precursor_qvalue from 'qvalue' (derived.py:32). CDAP output: peptide_qvalue 94/94 filled, precursor_qvalue 94/94 filled. Spectronaut's EG.Qvalue, emitted as 'qvalue', therefore fills precursor_qvalue.
- **Spec:** converter-coverage.md:181: CDAP peptide_qvalue is 'Not implemented: the best PepQvalue is stored as peptide_level_qvalue in additional_scores'. Spectronaut: 'Not implemented: ... EG.Qvalue is stored as qvalue in additional_scores'. The page does not mention that the typed columns are derived from these names.
- **Correction:** The CDAP contradiction holds. The Spectronaut cell is correct for peptide_qvalue, which stays null. EG.Qvalue filling precursor_qvalue is undocumented, not contradicted: the page has no row for precursor_qvalue, global_qvalue or psm.qvalue.
- **Reader impact:** Readers who use the coverage page as the contract expect nulls, but get values whose level differs (see INFRA-01).
- **Suggested:** Build the coverage tables from converter output after the writer's derivation, or document the derived typed columns in each cell.
- **Evidence:** `docs/spec/converter-coverage.md:181`; `qpx/writers/derived.py:31-32`; `run: audit-run:infra/cdap/cdap.feature.parquet count(peptide_qvalue)=94`

### INFRA-10: Spec both forbids and defines a top-level PEP on features {#infra-10}

*feature* · contradictory · confirmed · checklist F5, P5

- **Now:** feature.yaml:37-41 defines posterior_error_probability as 'PEP for the peptide match'. Producers fill it differently: MaxQuant evidence PEP; FragPipe 1 - probability of the matching PSM; OpenMS the first or best PID's PEP for the run (openms feature_adapter.py:736-744); quantms and CDAP null.
- **Spec:** psm.md:286: 'PEP is defined only in the PSM view... The feature and peptide views do not carry PEP as a top-level field'. feature.md:56 lists feature.posterior_error_probability. converter-coverage.md:180 describes how it is filled.
- **Reader impact:** Readers cannot tell whether feature PEP is the identifying PSM's PEP, a peptide-level PEP, or an aggregate, nor whether a null means a transferred row or 'not provided'.
- **Suggested:** Define feature PEP as the PEP of the identifying PSM from this run, null for transferred rows. Delete the psm.md:286 sentence.
- **Evidence:** `qpx/core/data/schemas/feature.yaml:37-41`; `docs/spec/psm.md:285-286`; `docs/spec/feature.md:56`; `docs/spec/converter-coverage.md:180`; `qpx/converters/openms_consensus/feature_adapter.py:736-744`

### INFRA-11: identification_type is not enforced at write time and has no 'unidentified' {#infra-11}

*feature* · diverges from spec · confirmed · checklist D1, F3, F4

- **Now:** The vocabulary is direct/transferred/secondary (qpx/core/data/validation.py:70-71). The write_table path gives only a warning (base.py:414-419, non-strict). The write_batch path, used by MaxQuant, PD, FragPipe, mzIdentML and quantms, runs no value checks (base.py:599-617). Writer test: 'MBR' and 'unidentified' were written silently via write_batch. qpxc validate (strict) then reports an ERROR. Writers do no case or synonym normalization.
- **Spec:** feature.md:145: direct | transferred | secondary, null when the tool does not say. Nothing covers unidentified features (#378 undecided).
- **Reader impact:** Published files can hold values outside the vocabulary that readers do not expect. Unidentified rows cannot be marked; they are only recognizable by an empty peptidoform.
- **Suggested:** The writer raises on values outside the vocabulary on every path. Decide whether 'unidentified' is a value or such rows go to a separate view.
- **Evidence:** `qpx/core/data/validation.py:70-74`; `qpx/writers/base.py:414-421`; `qpx/writers/base.py:599-617`; `docs/spec/feature.md:145`; `writer test + qpxc validate --file audit-run:infra/fb.feature.parquet`

### INFRA-13: A tool's zero intensity is encoded three ways across producers {#infra-13}

*feature* · inconsistent across producers · confirmed · checklist D6, F6

- **Now:** No writer-level normalization exists. MaxQuant, DIA-NN and FragPipe turn a 0 into an entry with null intensity (converter-coverage.md:60,194). OpenMS-consensus drops map elements <=0, so there is no entry (converter-coverage.md:60). quantms-MSstats keeps 0.0 as written (test: Intensity 0 -> {'LFQ', 0.0}). The record path also accepts NaN, negative values and an empty list [] as well as null for 'no entries' (writer test).
- **Spec:** intensities.md:123-128: null entry = identified, not quantified; no entry = not identified; 0 = measured zero; converters turn zero sentinels into null; a feature with no entries has a null list. producer-guide.md:41-48 says the same.
- **Reader impact:** The same upstream 'not quantified' state is read as measured zero (quantms), as not identified (OpenMS), or as identified but not quantified (MaxQuant). This skews imputation and the count of missing values.
- **Suggested:** Each producer declares whether 0 is a sentinel for it. The writer normalizes to null, rejects NaN and negatives, and coerces [] to null. OpenMS emits null entries for identified-but-unquantified maps.
- **Evidence:** `docs/spec/intensities.md:121-129`; `docs/spec/producer-guide.md:41-48`; `docs/spec/converter-coverage.md:60`; `docs/spec/converter-coverage.md:194`; `qpx/converters/quantms_msstats/feature_adapter.py:385-406`; `run: audit-run:infra/scripts/msstats_exp.py ms1`

### INFRA-14: quantms-MSstats aborts on a missing intensity and keeps zeros; provenance misstates this {#infra-14}

*feature* · diverges from spec · partly confirmed (see correction) · checklist D6, G1

- **Now:** _validate_normalized_rows raises when any row has a null or non-finite intensity (quantms_msstats/feature_adapter.py:385-406), so one NA aborts the whole conversion (test ms3). Zeros pass through as measured values. The provenance config states 'missing_values: typed null' (converter.py:114), but a null intensity is never written.
- **Spec:** producer-guide.md:46: 'Identified, not quantified -> Entry with null intensity'. intensities.md:126-127: zero sentinels become null.
- **Correction:** Abort-on-NA and zero passthrough are real. The 'provenance misstates' part is weak: config 'missing_values: typed null; scan uses an empty list' is a generic statement about nullable columns, not a claim about null intensities.
- **Reader impact:** MSstats inputs with NA rows cannot be converted. Converted ones may present upstream zeros as measured zeros. The provenance text misleads.
- **Suggested:** Map NA to a null-intensity entry. Decide, with quantms, whether its 0 is a sentinel. Fix the provenance wording.
- **Evidence:** `qpx/converters/quantms_msstats/feature_adapter.py:385-410`; `qpx/converters/quantms_msstats/converter.py:110-119`; `run: msstats_exp.py ms3 -> 'rows contain invalid required values'`

### INFRA-15: imputed_labels is not checked against the labels of intensities {#infra-15}

*feature* · unspecified · partly confirmed (see correction) · checklist D7

- **Now:** Neither the writer nor qpxc validate (strict, --semantic, --check-identity) relates imputed_labels to intensities[].label. imputed_labels=['TMT999'] on an LFQ feature was written and passed validation. No code outside the converters and mudata reads it (mudata.py:213).
- **Spec:** feature.md:148: 'Labels of intensities entries whose value the tool imputed'. intensities.md:131-133: an imputed value is an ordinary entry with its label listed.
- **Correction:** The rule is not unspecified: feature.md:148 defines imputed_labels as 'Labels of intensities entries', which implies a subset. The gap is enforcement (writer and validator), so the status is diverges_from_spec, not unspecified.
- **Reader impact:** An imputed_mask built from a label that does not exist is silently wrong. Readers cannot trust that every imputed label has a value.
- **Suggested:** Normative rule: imputed_labels is a subset of the labels of non-null entries in intensities; validate it on write.
- **Evidence:** `docs/spec/feature.md:148`; `docs/spec/intensities.md:131-133`; `qpx/mudata.py:213`; `writer test row scan 3; qpxc validate --semantic --check-identity on audit-run:infra/fb.feature.parquet`

### INFRA-18: identity_source misreports namespaced ids of unidentified features {#infra-18}

*feature* · diverges from spec · confirmed · checklist F3

- **Now:** For an unidentified feature (falsy peptidoform) that has a producer id, the base class records the source 'provided' (base.py:221-226). FeatureWriter then replaces the id with derive_id([run_file_name or quantification_unit_id, provided]) (feature.py:72-86). Test: the footer says 'mixed'. The stored id is neither the provided id nor composite-derived, and 'mixed' makes --check-identity skip the whole file (footer_validation.py:111-114).
- **Spec:** identity.md:21: 'provided: every id came with its row'; 'mixed: some of each'. feature.md:8 says the producer id is namespaced but does not give the hash bytes.
- **Reader impact:** Readers cannot reproduce or verify these ids, and one unidentified row turns off identity verification for every identified row.
- **Suggested:** Add the source value 'namespaced', or count these rows as derived. Specify the namespace hash in identity.md together with test vectors.
- **Evidence:** `qpx/writers/base.py:221-236`; `qpx/writers/base.py:347-353`; `qpx/writers/feature.py:49-86`; `qpx/core/data/footer_validation.py:107-114`; `docs/spec/identity.md:14-21`; `test: audit-run:infra/unid/u.feature.parquet identity_source=b'mixed'`

### INFRA-19: Unidentified features: spec requires a producer id; writer derives one; content undefined {#infra-19}

*feature* · contradictory · confirmed · checklist F3

- **Now:** FeatureWriter treats a falsy peptidoform as unidentified (feature.py:49-51). Without a producer id it derives one from the composite and only logs a warning (feature.py:64-66,89-100). peptidoform, sequence and charge are required by the schema. The writer accepted ''/''/0. The semantic validator flags '' as 'not canonical ProForma'.
- **Spec:** feature.md:8: 'For unidentified Features, a producer-supplied ID is required'. No page says what goes into the required peptidoform, sequence and charge, or how to mark the row (#378 open).
- **Reader impact:** Unidentified rows may appear with an empty peptidoform and charge 0. Readers grouping by peptidoform get a bogus '' peptide. Conformance fails on such files.
- **Suggested:** Decide whether unidentified features belong in feature (nullable peptidoform plus identification_type='unidentified') or in a separate view. Enforce feature.md:8 or relax it.
- **Evidence:** `qpx/writers/feature.py:49-100`; `docs/spec/feature.md:8`; `qpx/core/data/schemas/feature.yaml:15-36`; `qpxc validate --semantic: "not canonical ProForma, e.g. ''"`

### INFRA-26: What one feature row is depends on the producer composite {#infra-26}

*feature* · contradictory · confirmed · checklist F1, F7

- **Now:** quantms-MSstats declares run_file_name, peptidoform, charge, rt, scan (feature_adapter.py:28). Isobaric inputs therefore give one row per MS2 scan, which is a PSM-level reporter row (test_tmt_channels_collapse_without_merging_locations: 2 scans give 2 rows). LFQ inputs give one row per precursor and run. MaxQuant adds rt_start, rt_stop and scan; OpenMS adds rt, scan, observed_mz and consensus_rt. DIA-NN, Spectronaut and CDAP use one row per precursor and run.
- **Spec:** producer-guide.md:12: 'One row per precursor, charge and run'. identity.md:27-35 lists composites that allow several rows. #217 defines a feature as a chromatographic peak.
- **Reader impact:** Summing or counting feature rows per precursor and run gives different results by producer. quantms TMT 'features' are really PSMs.
- **Suggested:** Define a feature as one quantified precursor signal per run. Footer key feature_granularity (precursor_run, peak, spectrum). Isobaric per-spectrum rows are aggregated, or flagged as spectrum-level.
- **Evidence:** `qpx/converters/quantms_msstats/feature_adapter.py:28`; `tests/converters/test_quantms_msstats_converter.py:151-189`; `docs/spec/producer-guide.md:12`; `docs/spec/identity.md:27-35`

### INFRA-27: quantms-MSstats claims own-run identification for every row {#infra-27}

*feature* · inconsistent across producers · confirmed · checklist F4, D1, F10

- **Now:** _feature_record sets id_run_file_name = run_file_name for every row and never sets identification_type, which stays null (quantms_msstats/feature_adapter.py:501-526). Requantified (transferred) ProteomicsLFQ features are therefore indistinguishable. For the same quantms data, openms-consensus writes id_run_file_name null and 'transferred' for runs without their own identification (openms_consensus/feature_adapter.py:541-547).
- **Spec:** feature.md:145: id_run_file_name is null for a transferred row. converter-coverage.md:199,203: quantms identification_type 'Not provided', id_run_file_name 'Filled'.
- **Reader impact:** Two qpx routes from one quantms analysis disagree about MBR. A reader trusting id_run_file_name counts transferred quantities as MS/MS-identified.
- **Suggested:** When the input cannot tell, write id_run_file_name null (unknown) instead of copying run_file_name. Footer transferred=unknown.
- **Evidence:** `qpx/converters/quantms_msstats/feature_adapter.py:501-526`; `qpx/converters/openms_consensus/feature_adapter.py:523-547`; `docs/spec/feature.md:145`; `docs/spec/converter-coverage.md:199,203`; `run: msstats_exp.py ms1 id_run_file_name==run_file_name, identification_type None`

### INFRA-28: quantms-MSstats: unique is dataset-wide but pg_accessions is per row {#infra-28}

*feature* · ambiguous · confirmed · checklist F9, P8, F8

- **Now:** is_unique = one distinct ProteinName per sequence across the whole file and no ';' (queries.py:167-169). pg_accessions is the union of ProteinName rows within each (run, peptidoform, charge, rt, scan) group (feature_adapter.py:463-468,501-505). Test: run2 PEPTIDEK has pg_accessions [P1] but unique=False, because run1 also listed P2. anchor_protein is the first accession after alphabetical sorting of ProteinName values. Decoys are kept when every accession starts with decoy_, rev_ or reverse_ (:80-82).
- **Spec:** tool-mappings.md (feature table): anchor_protein 'First ProteinName accession'. The provenance config describes unique (converter.py:116). The meaning of the per-row pg_accessions is not specified.
- **Reader impact:** A row's unique flag can contradict its own pg_accessions. The same peptide gets different protein sets in different runs.
- **Suggested:** Compute pg_accessions and unique from the same scope (dataset-wide per peptidoform) and take the anchor in source order.
- **Evidence:** `qpx/converters/quantms_msstats/queries.py:167-169`; `qpx/converters/quantms_msstats/feature_adapter.py:80-82`; `qpx/converters/quantms_msstats/feature_adapter.py:455-505`; `run: msstats_exp.py ms1 rows`

### INFRA-29: quantms-MSstats does not say what its intensity measures {#infra-29}

*feature* · unspecified · partly confirmed (see correction) · checklist F6, F1

- **Now:** The MSstats Intensity is stored as-is under label LFQ or the TMT channel (feature_adapter.py:350-357,518-521). Depending on the input it is an OpenMS LFQ feature intensity, a per-scan TMT reporter intensity, or a DIA-NN-derived quantity whose type is not known (not verifiable here). Provenance has no isotope_corrected or normalization keys and step 1 has parameters None (converter.py:120-130). An LFQ precursor appearing twice in one run with different intensities aborts the conversion (feature_adapter.py:425-437; test ms2).
- **Spec:** intensities.md:211-213 puts processing state in provenance, as MaxQuant (maxquant/converter.py:122) and OpenMS (openms_consensus/provenance.py:66) do. No key covers apex vs area vs reporter.
- **Correction:** The missing isotope/normalization keys are spec-compliant ('unknown'). The real gap is that no key exists for intensity kind (area/apex/reporter). The duplicate-precursor abort happens only when the rows share rt and scan, e.g. when there is no RetentionTime column; with RT they become separate rows.
- **Reader impact:** Readers cannot tell whether quantms intensities are normalized, isotope-corrected or area-based, so cross-dataset comparison is unsafe.
- **Suggested:** Require provenance or footer keys intensity_kind (ms1_area, ms1_apex, reporter) plus isotope_corrected and normalization. Use unknown when the input does not say.
- **Evidence:** `qpx/converters/quantms_msstats/feature_adapter.py:350-357`; `qpx/converters/quantms_msstats/feature_adapter.py:425-437`; `qpx/converters/quantms_msstats/converter.py:120-130`; `docs/spec/intensities.md:211-213`; `qpx/converters/maxquant/converter.py:122`

### INFRA-M2: No invariant ties identification_type to id_run_file_name, scan or PEP {#infra-m2}

*feature* · ambiguous · found by the verifier · checklist D1, F4, F10

- **Now:** The only check on identification_type is its vocabulary (validation.py:70-71). grep finds no other rule in qpx/core. quantms writes id_run_file_name=run_file_name with identification_type null (quantms_msstats/feature_adapter.py:525). MaxQuant is self-consistent (transferred rows: scan [], id_run_file_name null, PEP null).
- **Spec:** feature.md:145: id_run_file_name is null for a transferred row. producer-guide.md:26: transferred rows have empty scan. PEP or q-value on transferred rows: nothing (#378 proposes null).
- **Reader impact:** Readers cannot rely on any one field to detect MBR. A 'transferred' row may carry a scan or id_run, and a null type may hide transfers.
- **Suggested:** Make it normative and validated: direct gives id_run_file_name=run_file_name and a non-empty scan. Transferred gives id_run_file_name null, scan [], PEP null, and match_* allowed. Null type gives id_run_file_name null.
- **Evidence:** `qpx/core/data/validation.py:70-71`; `qpx/converters/quantms_msstats/feature_adapter.py:525`; `docs/spec/feature.md:145`; `docs/spec/producer-guide.md:26`; `run: work/verify-infra/mq feature grouped by identification_type/scan/id_run_file_name`

### INFRA-M4: producer-guide names a non-existent 'imputation' flag {#infra-m4}

*feature* · contradictory · found by the verifier · checklist D6, D7

- **Now:** No schema has an 'imputation' column (grep over qpx/core/data/schemas). Imputation is carried by feature.imputed_labels (feature.yaml:235) and pg.is_imputed.
- **Spec:** producer-guide.md:48: 'Imputed by the tool -> The value, with `imputation` flagged'. intensities.md:131-132 names imputed_labels / is_imputed.
- **Reader impact:** Producers following the guide look for a field that does not exist and may leave imputed values unmarked.
- **Suggested:** Replace it with feature.imputed_labels / pg.is_imputed, and link to intensities.md#missing-values.
- **Evidence:** `docs/spec/producer-guide.md:48`; `docs/spec/intensities.md:131-132`; `qpx/core/data/schemas/feature.yaml:235-239`

## Specification text

### SPEC-01: What one feature row is: 'one chromatographic peak' vs per-producer entities {#spec-01}

*feature* · contradictory · confirmed · checklist F1, F7, X1

- **Now:** MaxQuant composite adds rt_start,rt_stop,scan (identity.md:32). Running maxquant_full gave 14,086 (peptidoform,charge,run) groups with >1 row (61,719 rows): 9,639 groups mix direct and transferred rows, 2,594 have >1 quantified direct row, and 12,758 direct rows are MS/MS-only (null intensity). In FragPipe a row is one experiment (converter-coverage.md:213). CDAP and PD aggregate PSMs (cdap/feature_adapter.py:1-5, converter-coverage.md:324).
- **Spec:** data-model.md:64 'one quantified chromatographic peak: a peptidoform+charge in one run'. data-model.md:74 'per upstream Feature entity'. producer-guide.md:12 'One row per precursor, charge and run'. feature.md:301 'aggregates one or more PSMs'.
- **Reader impact:** Pivoting or summing by (peptidoform,charge,run) double counts, or mixes MS/MS-only, transferred and peak rows. Nothing says whether duplicate precursor/run rows add up or are alternatives.
- **Suggested:** Normative grain: one row per (peptidoform, charge, run). Otherwise declare the grain in the footer and say whether rows add up. Keep MS/MS-only evidence in psm, not feature.
- **Evidence:** `docs/spec/data-model.md:64-67,74`; `docs/spec/producer-guide.md:12`; `docs/spec/identity.md:29-35`; `docs/spec/converter-coverage.md:213-216,324`; `qpx/converters/cdap/feature_adapter.py:1-12`; `ran qpxc convert maxquant (maxquant_full) -> audit-run:spec/mq`

### SPEC-02: Feature top-level PEP exists in schema but psm.md says it does not {#spec-02}

*feature* · contradictory · confirmed · checklist F5, P5

- **Now:** feature.yaml:37-41 defines posterior_error_probability ('PEP for the peptide match') and feature.md:56 lists it. Converters fill it: MaxQuant PEP (null on MBR rows), FragPipe 1-prob of the matching PSM, Spectronaut EG.PEP (converter-coverage.md:180).
- **Spec:** psm.md:285-286: PEP 'is defined only in the PSM view... The feature and peptide views do not carry PEP as a top-level field'. scores.md:55: 'top-level field in the PSM view'.
- **Reader impact:** Readers cannot tell what feature PEP describes (best PSM, evidence, or precursor-level EG.PEP), or whether they should look for it at all.
- **Suggested:** Keep feature.posterior_error_probability and define its level: the identifying PSM/precursor in this run, null for transferred rows. Delete the psm.md note.
- **Evidence:** `qpx/core/data/schemas/feature.yaml:37-41`; `docs/spec/feature.md:56`; `docs/spec/psm.md:285-286`; `docs/spec/scores.md:55`; `docs/spec/converter-coverage.md:180`

### SPEC-03: Examples write feature.pg_accessions as strings; 'reserved' position fields contradicted {#spec-03}

*feature* · contradictory · confirmed · checklist F9

- **Now:** FeatureWriter rejects the feature.md example value ['P04217','P04217-2'] with 'Could not convert P... expecting tuple' (run in audit-run:spec/ex). producer-guide.md:120 fills pg_protein structs with start/end/pre/post. CDAP stores pre/post in pg_accessions (converter-coverage.md:197).
- **Spec:** feature.yaml:151-153: list<pg_protein>, 'its position fields are reserved, pg_positions holds positions'. feature.md:99 type array[struct].
- **Reader impact:** Producers who copy the examples fail. Readers do not know whether start/end/pre/post in feature.pg_accessions mean anything.
- **Suggested:** Change the examples to structs. Say that the pg_protein position fields must be null on feature (pre/post included) or are allowed. Alternatively retype it to list<string> like pg.
- **Evidence:** `docs/spec/feature.md:198,261`; `docs/spec/producer-guide.md:120`; `qpx/core/data/schemas/feature.yaml:151-153`; `docs/spec/converter-coverage.md:197`; `audit-run:spec/ex_test.py output`

### SPEC-04: id_run_file_name: 'may differ from run_file_name' vs 'null for a transferred row' {#spec-04}

*feature* · contradictory · confirmed · checklist F4, F10

- **Now:** MaxQuant, FragPipe and OpenMS write null on transferred rows (converter-coverage.md:57,203). In maxquant_full, none of the 104,491 transferred rows has an id_run_file_name, scan or PEP. FragPipe also writes null for experiments with several raw runs.
- **Spec:** feature.md:138: 'run file containing the best PSM ... (may differ from run_file_name)'. feature.yaml:207-209: 'Run file containing best PSM'. feature.md:145: 'null for a transferred row'. producer-guide.md:26: empty scan for transferred rows.
- **Reader impact:** 'May differ' invites producers to write the MBR donor run while others write null, so one column carries two meanings. The donor run cannot be recovered.
- **Suggested:** id_run_file_name = run of this row's identifying PSM (run_file_name or a fraction of the unit), null otherwise. Add explicit donor fields for transfers if needed.
- **Evidence:** `docs/spec/feature.md:138,145`; `qpx/core/data/schemas/feature.yaml:207-209`; `docs/spec/producer-guide.md:26`; `docs/spec/converter-coverage.md:57,203`; `audit-run:spec/mq (identification_type vs id_run_file_name counts)`

### SPEC-05: feature.scan has three meanings and 'best PSM' is undefined {#spec-05}

*feature* · contradictory · confirmed · checklist F2, F10, P10

- **Now:** OpenMS consensusXML concatenates every supporting spectrum (feature.md:26-27, scan.md:90-92,123). DIA uses the apex MS2 scan (scan.md:111). DDA uses the 'best PSM', chosen differently per producer: FragPipe lowest PEP then scan (converter-coverage.md:223-224), PD lowest q-value (:324), CDAP lowest E-value (cdap/feature_adapter.py:6).
- **Spec:** feature.md:139: 'Scan identifier of the best PSM'. producer-guide.md:23: 'scan number(s) of the MS2 spectrum'. 'Best' is never defined.
- **Reader impact:** Joining feature.scan to mz or psm returns wrong or unsplittable spectra (concatenated nativeId parts). Identification fields come from PSMs picked by different rules.
- **Suggested:** feature.scan = one designated identifying PSM in id_run_file_name. 'Best' = the tool's pick, else lowest PEP or q, ties to lowest scan (per #378). Multi-spectrum lists go in another column.
- **Evidence:** `docs/spec/feature.md:26-27,139`; `docs/spec/scan.md:90-92,111,123`; `docs/spec/producer-guide.md:23`; `docs/spec/converter-coverage.md:223-224,324`; `qpx/converters/cdap/feature_adapter.py:6`

### SPEC-06: Feature rows matching conflicting peptidoforms have no rule {#spec-06}

*feature* · unspecified · partly confirmed (see correction) · checklist F2

- **Now:** CDAP leaves ambiguous top-ranked PSMs out of features (cdap/feature_adapter.py:7-9). PD drops PSMs that 'PSM Ambiguity' marks Rejected (converter-coverage.md:324). OpenMS keeps the feature but nulls id_run_file_name and the protein fields on conflicts (converter-coverage.md:57,62-64).
- **Spec:** Nothing general. data-model/feature pages assume one peptidoform per feature, with no rule for a peak that several peptidoforms explain.
- **Correction:** converter-coverage.md:62-64 is about protein-group conflicts, not peptidoform conflicts. For OpenMS, a conflicting peptide assignment nulls only id_run_file_name (:57; guide/convert.md:1011-1015), and consensus features with an identical identity are deduplicated (guide/convert.md:978-980). 'secondary' is the spec's only mechanism. The gist holds.
- **Reader impact:** The same ambiguous peak is missing, emitted once, or emitted with nulls, depending on the producer, which skews quant completeness.
- **Suggested:** Normative: emit one row per (peak, peptidoform), share an id and flag the ambiguity. Alternatively emit none and record the conflict in psm.
- **Evidence:** `qpx/converters/cdap/feature_adapter.py:7-9`; `docs/spec/converter-coverage.md:57,62-64,324`

### SPEC-07: Unidentified features described but not representable; enum not in YAML {#spec-07}

*feature* · unspecified · confirmed · checklist D1, F3

- **Now:** The enum direct/transferred/secondary exists only in code (qpx/core/data/validation.py:71). A row with identification_type='unidentified' is written, then 'qpxc validate' reports an ERROR (run). The OpenMS consensus converter omits features without hits (psm.md:175).
- **Spec:** feature.md:4-10 describes unidentified Features, which need producer ids. feature.yaml:16-36 makes sequence, peptidoform and charge required. feature.yaml:210-214 types identification_type as a string with no enum and lists three values.
- **Reader impact:** Unidentified features cannot be emitted legally. The closed vocabulary cannot be discovered from the machine schema.
- **Suggested:** Either exclude unidentified features from feature (say so) or add 'unidentified' with nullable identification fields. Put the enum in feature.yaml.
- **Evidence:** `docs/spec/feature.md:4-10`; `qpx/core/data/schemas/feature.yaml:16-36,210-214`; `qpx/core/data/validation.py:71`; `docs/spec/psm.md:175`; `audit-run:spec/ex/e.feature.parquet validate ERROR`

### SPEC-08: identification_type mapping mixes MBR with imputation and MS/MS-backed transfers {#spec-08}

*feature* · inconsistent across producers · confirmed · checklist D1, F4, D6

- **Now:** Spectronaut maps EG.IsImputed=true to transferred and also lists those labels in imputed_labels (converter-coverage.md:199,201; spectronaut/feature_adapter.py:342). MaxQuant MULTI-MATCH-MSMS, which has MS/MS in the run, becomes transferred (converter-coverage.md:199). CDAP and PD hard-code direct. DIA-NN and MSstats write null.
- **Spec:** feature.md:145: direct = 'its own MS/MS', transferred = 'matched from another run'. intensities.md:131: imputed = value not measured.
- **Reader impact:** 'transferred' sometimes means imputed, and some rows with MS/MS are hidden as transferred, so filtering on identification_type depends on the producer.
- **Suggested:** Precedence: any MS/MS in this run means direct. Imputation goes only in imputed_labels/is_imputed. Make per-tool mappings normative.
- **Evidence:** `docs/spec/converter-coverage.md:199,201`; `qpx/converters/spectronaut/feature_adapter.py:342`; `docs/spec/feature.md:145`; `docs/spec/intensities.md:131-133`

### SPEC-09: Which identification feature confidence fields describe is unspecified {#spec-09}

*feature* · unspecified · confirmed · checklist F5, D2

- **Now:** MaxQuant: evidence PEP, null on transferred rows. FragPipe: PEP of the best PSM in the experiment's raw file (converter-coverage.md:223-226). PD: lowest-q PSM (:324). CDAP: best-E-value PSM. OpenMS: best PSM (feature.md:130). Only MaxQuant fills match_*.
- **Spec:** feature.md:56 'PEP for the peptide match'. feature.yaml:42-62 give q-value level names only. Nothing says whether values describe this run, the best run, or the transfer donor.
- **Reader impact:** Filtering by PEP or q drops every MBR row for some producers but not others, so the same column means different evidence.
- **Suggested:** Feature confidence = this run's identifying PSM or precursor, null for transferred rows. Transfer confidence goes only in match_*, cross-run confidence in global_qvalue.
- **Evidence:** `docs/spec/feature.md:56,130`; `qpx/core/data/schemas/feature.yaml:42-62`; `docs/spec/converter-coverage.md:180,200,223-226,324`

### SPEC-12: match_qvalue is float32 with no 0-1 range, unlike every other q-value {#spec-12}

*feature* · contradictory · confirmed · checklist D2

- **Now:** feature.yaml:220-224 types match_qvalue as float32 with no minimum or maximum. types.yaml:60 channel_detail.qvalue is float64 with no range.
- **Spec:** Other q-values are float64 with minimum 0, maximum 1 (feature.yaml:42-62, psm.yaml:62-68). versioning.md:82-85: 'YAML schemas declare minimum and maximum for q-values'.
- **Reader impact:** Precision is lost and out-of-range transfer q-values pass validation.
- **Suggested:** Make every q-value column and struct field float64 with a 0-1 range.
- **Evidence:** `qpx/core/data/schemas/feature.yaml:220-224`; `qpx/core/data/schemas/types.yaml:60`; `docs/spec/versioning.md:82-85`

### SPEC-21: Feature rt/rt_start/rt_stop: clock, unit and apex-vs-MS2 undefined {#spec-21}

*feature* · ambiguous · confirmed · checklist F10

- **Now:** MaxQuant: rt = 'Retention time' on the raw clock; rt_start/rt_stop = calibrated start/finish (maxquant/feature_adapter.py:574-583). In maxquant_full, 83,413 of 243,821 rows have rt outside [rt_start, rt_stop]; rt_aligned is always inside. Feature rt is the apex for OpenMS and PD, RTAtPrecursorHalfElution for CDAP, the MS2 time for MSstats.
- **Spec:** feature.yaml:255-260: 'RT window start/end', no unit or clock. feature.md:65: 'Precursor retention time'. feature.yaml:106-110: rt is on the run's own clock.
- **Reader impact:** RT windows contradict rt, and apex, MS2 time and midpoint are mixed. rt also feeds the identity composite.
- **Suggested:** rt = apex of the quantified peak (the MS2 time when there is no peak). rt_start/rt_stop on the same clock as rt, in seconds. Aligned values only in rt_aligned.
- **Evidence:** `qpx/core/data/schemas/feature.yaml:102-110,255-260`; `docs/spec/feature.md:65,69-70`; `qpx/converters/maxquant/feature_adapter.py:574-583`; `audit-run:spec/mq rt window query`

### SPEC-23: Missing-value encoding (null entry/no entry/no row) differs by producer {#spec-23}

*feature* · inconsistent across producers · partly confirmed (see correction) · checklist D6, F6

- **Now:** OpenMS: a map element with intensity <=0 gets no entry (converter-coverage.md:60). DIA-NN, CDAP and MaxQuant write a null entry. FragPipe TMT writes no feature for a precursor without reporter signal (:211). CDAP writes null when no channel is positive (:194).
- **Spec:** intensities.md:123-127: null entry = identified, not quantified; no entry = not identified or no sample; 'a feature with no entry at all has a null intensities list', which implies a row identified in no label. producer-guide.md:43-48.
- **Correction:** CDAP conforms: it writes a null entry per unquantified channel. The divergent producers are OpenMS consensus (drops elements <=0) and FragPipe TMT (no row).
- **Reader impact:** Not-identified and not-quantified cannot be told apart reliably, and identified but unquantified precursors are missing for some producers.
- **Suggested:** Every identified label gets an entry (null if unquantified). A feature row always has at least one entry. 'No entry' only for labels absent from the run.
- **Evidence:** `docs/spec/intensities.md:121-133`; `docs/spec/producer-guide.md:41-48`; `docs/spec/converter-coverage.md:60,194,211`

### SPEC-M1: MBR rows that have their own MS/MS keep scan, id_run_file_name and PEP {#spec-m1}

*feature* · diverges from spec · found by the verifier · checklist F4, D1, F10

- **Now:** maxquant/evidence_fields.py:16 maps MULTI-MATCH-MSMS to transferred. feature_adapter.py:655-657 nulls id_run_file_name only for 'MULTI-MATCH'. Ran 5 injected MULTI-MATCH-MSMS rows (maxquant_simple): every one became transferred with a non-null id_run_file_name, a scan and a PEP.
- **Spec:** feature.md:145: 'id_run_file_name is null for a transferred row'. producer-guide.md:26: a transferred row has an empty scan.
- **Reader impact:** Rules like 'null id_run_file_name or empty scan means MBR' fail, and transferred rows look identified by their own MS/MS.
- **Suggested:** State the precedence: MS/MS in this run means direct (or add a typed 'transferred with MS/MS' value). Keep the null rule for pure transfers and align the converter.
- **Evidence:** `qpx/converters/maxquant/evidence_fields.py:8-17`; `qpx/converters/maxquant/feature_adapter.py:655-657`; `docs/spec/feature.md:145`; `docs/spec/producer-guide.md:26`; `audit-run:verify-spec/mqmm (transferred: 1474 rows, 5 with id_run_file_name/scan/PEP)`

### SPEC-M2: feature.run_file_name is a placeholder run for multi-run units {#spec-m2}

*feature* · contradictory · found by the verifier · checklist F1, F10, X1

- **Now:** FragPipe writes run_file_name = the first raw run of the experiment or plex, plus quantification_unit_id (fragpipe/feature_adapter.py:350). The quantity can come from other fractions.
- **Spec:** feature.md:89: 'The run file name that contains the feature'. intensities.md:136: each LFQ feature row is one MS run and one sample. feature.yaml:90-94: quantification_unit_id is optional and gives no rule for run_file_name.
- **Reader impact:** Per-run QC, RT and scan lookups, and run-level joins, assign unit-level quantities to the first fraction.
- **Suggested:** run_file_name = the measured run only. Unit-level rows must set quantification_unit_id and declare the grain in the footer, or be emitted per run.
- **Evidence:** `qpx/converters/fragpipe/feature_adapter.py:350`; `docs/spec/converter-coverage.md:209,213-216`; `docs/spec/feature.md:89-90`; `docs/spec/intensities.md:136`; `qpx/core/data/schemas/feature.yaml:86-94`

## Validator and read API

### VALIDATOR-07: Feature.confident() filters on a different q-value level for each producer {#validator-07}

*feature* · inconsistent across producers · confirmed · checklist D4, F5, P5

- **Now:** Feature.confident() filters is_decoy=false AND COALESCE(precursor_qvalue, peptide_qvalue)<=x (feature.py:52-54). The writer fills precursor_qvalue from 'qvalue', 'q_value' and 'psm_level_qvalue' (derived.py:30-34). Per producer: CDAP gets its precursor_qvalue (94/94 rows) from psm_level_qvalue; OpenMS has only peptide_qvalue; MaxQuant full and FragPipe have none, so confident() returns 0 of 243,821 rows. The filter ignores global_qvalue and identification_type.
- **Spec:** scores.md:44 maps psm_level_qvalue to feature.precursor_qvalue ('precursor, per run'). scores.md:38 calls confident() 'the portable way to filter'.
- **Reader impact:** The same call keeps a PSM-level, precursor-level or peptide-level FDR, or nothing at all, depending on the producer.
- **Suggested:** Stop mapping PSM-level names into precursor_qvalue. Have confident() take an explicit level, and raise or warn when that column is all null.
- **Evidence:** `qpx/core/data/feature.py:52-54`; `qpx/writers/derived.py:30-34`; `run: ds.feature.confident().count() on audit-run:maxquant/full -> 0`; `query cdap.feature: precursor_qvalue 94/94, scores psm_level_qvalue`; `docs/spec/scores.md:38-47`

### VALIDATOR-09: Only the identification_type vocabulary is enforced; the rules tied to its values are not {#validator-09}

*feature* · unspecified · confirmed · checklist D1, D2, F4

- **Now:** The closed vocabulary direct/transferred/secondary is enforced and 'MBR' gives an error (validation.py:71). These all pass --conformance: 'transferred' with a non-empty scan and id_run_file_name equal to run_file_name; 'direct' with an empty scan; match_score or match_qvalue on a 'direct' row. There is no 'unidentified' value. The read API never reads identification_type.
- **Spec:** producer-guide.md:26 says transferred rows have an empty scan. feature.md:145 says id_run_file_name is null for a transferred row. feature.md:146-147 says match_* apply to transferred identifications.
- **Reader impact:** An MBR row can claim a spectrum in this run, and readers cannot rely on the scan or id_run_file_name invariants.
- **Suggested:** Validate that transferred implies empty scan, null id_run_file_name and null PEP; that direct implies a non-empty scan; and that match_* is non-null only on transferred rows. Decide on 'unidentified' (GitHub #378).
- **Evidence:** `qpx/core/data/validation.py:70-74`; experiment mutate.py: f_idtype_mbr_vocab -> ERROR; f_transferred_with_scan_and_same_idrun, f_direct_empty_scan, f_match_score_on_direct -> exit 0; `grep identification_type qpx (non-converter) -> validation.py:71 only`

### VALIDATOR-10: 'One row per precursor, charge and run' is not enforced, and MaxQuant breaks it {#validator-10}

*feature* · contradictory · confirmed · checklist F1, F2, P6

- **Now:** Uniqueness is checked only on the composite the footer declares (footer_validation.py:157-172). MaxQuant declares peptidoform,charge,run_file_name,rt_start,rt_stop,scan. maxquant/full has 14,086 (peptidoform, charge, run) keys with more than one row, 6,370 of them direct+transferred pairs, and still passes --conformance. If I add a second row with a different rt under the default composite, it also passes.
- **Spec:** producer-guide.md:12: 'One row per precursor, charge and run'. identity.md:27-35 lets composites include rt, scan and rt_start, which allows several rows per precursor and run.
- **Reader impact:** Readers that pivot on (peptidoform, charge, run) get collisions, which are summed or averaged silently (see VALIDATOR-19).
- **Suggested:** State whether a feature is a peak or a precursor-run aggregate. If it is a peak, drop producer-guide.md:12. If it is an aggregate, validate the natural key independently of the composite.
- **Evidence:** `qpx/core/data/footer_validation.py:157-172`; `query audit-run:maxquant/full/maxquant.feature.parquet: 14086 dup keys`; `run: qpxc validate audit-run:maxquant/full --conformance -> VALID`; `docs/spec/producer-guide.md:12`; `docs/spec/identity.md:27-35`

### VALIDATOR-11: Unidentified features are supported by the writer but rejected by --semantic {#validator-11}

*feature* · contradictory · confirmed · checklist F3, D1

- **Now:** FeatureWriter handles rows with an empty peptidoform and namespaces their producer id (writers/feature.py:61-87). sequence, peptidoform and charge are non-null, so '' is used. Under --semantic, '' fails the canonical ProForma regex and gives a peptidoform_malformed ERROR (semantic_validation.py:23,85-86). identification_type has no 'unidentified' value.
- **Spec:** feature.md:8-10 says unidentified features need a producer-supplied id. Nothing defines what sequence, peptidoform or charge they hold.
- **Reader impact:** Unidentified features cannot pass conformance, and readers would see empty-string peptidoforms mixed with real ones.
- **Suggested:** Either forbid unidentified features in the feature view, or allow null peptidoform with identification_type='unidentified' and exempt them from semantic checks.
- **Evidence:** `qpx/writers/feature.py:61-87`; `qpx/core/semantic_validation.py:23,85-86`; `experiment mutate.py f_unidentified_empty_peptidoform -> ERROR peptidoform_malformed`; `docs/spec/feature.md:8-10`

### VALIDATOR-12: Feature scan_format is never checked, and feature scan has two meanings in the spec {#validator-12}

*feature* · contradictory · confirmed · checklist P7, F10, D8

- **Now:** _declared_scan_format returns None unless the view is psm (schema.py:134-139). A feature footer declaring scan_format=scan with scan [1,2,3] passes --conformance.
- **Spec:** producer-guide.md:23 tells producers to set scan_format for a DDA PSM or feature 'or the validator cannot check them'. scan.md:90-92 says the check does not apply to features. feature.md:139 says feature.scan is the best PSM's scan, while scan.md:90 says it may combine components of several spectra.
- **Reader impact:** Readers cannot interpret feature.scan as one spectrum id, and the producer guide promises a check that does not exist.
- **Suggested:** Define feature.scan as one spectrum (best or identifying PSM) and validate it, or move combined references to a separate list column.
- **Evidence:** `qpx/core/data/schema.py:134-139`; `experiment mutate.py f_scan_three_components_scanformat_scan -> exit 0`; `docs/spec/producer-guide.md:23`; `docs/spec/scan.md:90-92`; `docs/spec/feature.md:139`

### VALIDATOR-14: Canonical peptidoform is enforced only for syntax and UNIMOD accessions, not names or decimals {#validator-14}

*feature* · diverges from spec · confirmed · checklist D5, P9

- **Now:** _tag_mass treats any tag without ':' as canonical (semantic_validation.py:78-80), and the mass-shift regex accepts any number of decimals (:25). With ids re-derived, PEPTM[Oxidation]IDEK, [+15.994915] and [+15.99] all pass --conformance. The fixture feature_noncanonical_mod (test_conformance.py:82) fails only through a stale id. Lower-case 'unimod:' and PSI-MOD accessions are caught.
- **Spec:** peptidoform.md:24-30: names that have a Unimod entry must be written as UNIMOD:n, and mass shifts are signed with four decimals.
- **Reader impact:** The same molecule can appear as several peptidoform strings, which splits identity, ids and aggregation across producers.
- **Suggested:** Reject names that resolve to a Unimod entry and mass shifts that do not have exactly four decimals. Fix the fixture so it re-derives the id.
- **Evidence:** `qpx/core/semantic_validation.py:25,68-80`; `tests/conformance/test_conformance.py:82`; experiment mutate4.py z_named_mod_rederived, z_six_decimal_shift_rederived, z_two_decimal_shift_rederived, z_psm_named_mod_rederived -> exit 0; `docs/spec/peptidoform.md:24-30`

### VALIDATOR-17: Validator requires exact label matches but the read side falls back to the run's lone sample {#validator-17}

*feature* · diverges from spec · confirmed · checklist F6, X1

- **Now:** The validator reports an error for a (run, label) pair that matches no sample (label_validation.py:90-122,180-182). Read side: the SQL gives any label of a single-sample run to that sample (sample_sql.py:79-84). MuData falls back when a run has at most one distinct label, and also matches by run-file stem (mudata_obs.py:92-106). With run.samples label 'label free sample', validation fails, intensity('peptide') returns S1 relabelled, and intensity('protein') returns nothing.
- **Spec:** data-model.md:239-246: labels MUST match exactly, including for single-sample runs. intensities.md:63 says labels are compared exactly.
- **Reader impact:** For the same invalid data, peptide-level and protein-level reads disagree, and MuData can bind a quantity to a sample the spec says is unresolvable.
- **Suggested:** Remove the read-side fallbacks, or make the fallback a spec rule and validate it as a warning only.
- **Evidence:** `qpx/core/sample_sql.py:59-85`; `qpx/core/mudata_obs.py:92-106`; `qpx/core/label_validation.py:90-122`; `experiment mutate2.py x_run_label_free_sample: validate ERROR; ds.intensity('peptide') -> S1; ds.intensity('protein') -> []`

### VALIDATOR-18: Decoy and contaminant handling differs between read helpers {#validator-18}

*feature* · inconsistent across producers · confirmed · checklist P4, F8, G3

- **Now:** is_decoy is non-null (strict null check). Feature has no targets_only() (psm.py:27 and pg.py:53 do). Feature.for_quantification() keeps decoys and null-intensity rows (feature.py:91-115): on MaxQuant full, 260 REV__ rows and 12,758 null-intensity rows. Views and intensity() filter decoys. MuData uses IS NOT TRUE (mudata.py:278-281), and quantify uses NOT coalesce(is_decoy,false). Contaminants are never filtered, and feature and psm have no contaminant flag. Decoy accession markers are not checked.
- **Spec:** pg.md:99-103: consumers MUST filter decoys themselves. producer-guide.md:33-34: decoys keep their marker and never share an accession; contaminants are flagged on pg only.
- **Reader impact:** The 'quantification-ready' export carries decoys, and contaminants enter every count and matrix.
- **Suggested:** Give every view a uniform targets_only(). Decide whether for_quantification drops decoys. Add a contaminant flag on feature/psm, or a documented rule, and check decoy markers.
- **Evidence:** `qpx/core/data/feature.py:91-115`; `qpx/mudata.py:278-281`; `qpx/cli/transform/quantify.py:305`; `run: for_quantification on audit-run:maxquant/full -> 260 REV__ rows`; `experiment mutate.py p_decoy_target_accession, f_decoy_with_target_accession, g_contaminant_decoy -> exit 0`

### VALIDATOR-19: Read helpers sum duplicate features (MS/MS plus MBR) and ignore identification_type {#validator-19}

*feature* · ambiguous · partly confirmed (see correction) · checklist F2, F4, F7

- **Now:** ds.intensity('peptide') SUMs every feature of a key within a unit, across charges, fractions and duplicate precursor rows (dataset_queries.py:124-137). build_mudata defaults to duplicates='sum' (mudata.py:84-119). On maxquant/full, 12,672 cells are summed, mostly direct+transferred pairs. qpxc transform quantify sums all charges and peptidoforms per run and label (quantify.py:261-272). No helper can exclude transferred or secondary rows.
- **Spec:** views.md:212 documents summing over charges and fractions. anndata.md:153 documents MuData sum/mean. Nothing says whether two rows for one precursor and run are additive (split peaks) or redundant.
- **Correction:** anndata.md:153 explicitly names this case (a MaxQuant MS/MS feature plus an MBR feature of one precursor in one run) and documents sum as the MuData default, so it is specified, not ambiguous. The gap is that the spec never says whether that sum double-counts. intensity('peptide') sums these upstream without counting them.
- **Reader impact:** Abundances in MBR-heavy runs can be inflated by double counting. Mokume#117 reports false positives from transfers.
- **Suggested:** Specify what several rows per precursor and run mean (additive or alternative). Add an exclude-transferred option, and default MuData duplicates to error or max until that is decided.
- **Evidence:** `qpx/core/dataset_queries.py:124-137`; `qpx/mudata.py:84-119`; `run: build_mudata(audit-run:maxquant/full) -> uns qpx_duplicates {'aggregation':'sum','cells':12672}`; `docs/spec/views.md:212`; `docs/spec/anndata.md:153`

### VALIDATOR-20: Identification counts include MBR transfers and unfiltered rows {#validator-20}

*feature* · ambiguous · confirmed · checklist F4, P3, X1

- **Now:** identification_summary, run_summary and qc_view count feature rows and distinct sequences after excluding decoys only (views/api.py:21-29,172-183). There is no q-value or identification_type filter. MaxQuant full averages per run: 9,021 features vs 5,156 non-transferred, and 6,128 vs 3,779 peptides (+62%).
- **Spec:** views.md:128-140 describes the 'Identifications per run' counts and their decoy and group handling. It does not say whether transfers count as identifications.
- **Reader impact:** 'Identifications per run' plots and counts include peptides never identified by MS/MS in that run, and are not comparable across producers.
- **Suggested:** Count identifications as direct, non-null-type rows passing a declared FDR, and report transferred rows separately.
- **Evidence:** `qpx/views/api.py:21-29,172-183`; `query audit-run:maxquant/full feature: avg n_all 9021 vs n_direct 5156 per run`

### VALIDATOR-22: Label-set fields and 1.3 quant fields are not cross-checked or range-checked {#validator-22}

*feature* · diverges from spec · confirmed · checklist D7, D6

- **Now:** All of the following pass --conformance: imputed_labels, identified_labels, representative_label or channel_details labels that are absent from intensities; additional_intensities labels that are absent from run.samples; parent_ion_fraction=3.0; quant_ms_level=7; channel_details.qvalue=5.0; additional_intensity value NaN. feature.yaml declares no ranges for these fields.
- **Spec:** feature.md:154-156 gives parent_ion_fraction as 0-1 and quant_ms_level as 2 or 3. feature.yaml:239 says imputed_labels are labels of intensities. intensities.md:93-101 defines identified_labels and channel_details per label. intensities.md:127 says NaN is not a valid intensity.
- **Reader impact:** The imputed_mask and per-channel details can point at channels that do not exist, and out-of-range values look valid.
- **Suggested:** Add schema ranges (0-1, {2,3}, q-values in [0,1], finite additional intensities) and subset checks against the intensities labels.
- **Evidence:** `qpx/core/data/schemas/feature.yaml:144-148,173-182,235-254`; experiment mutate.py f_imputed_label_not_in_intensities, f_identified_label_not_in_intensities, f_representative_label_unknown, f_channel_details_unknown_label, f_parent_ion_fraction_gt1, f_quant_ms_level_7, f_add_intensity_nan, f_add_intensity_label_not_in_run -> exit 0

### VALIDATOR-23: unique and pg_group_id are not checked, although link_feature_pg relies on pg_group_id {#validator-23}

*feature* · unspecified · confirmed · checklist F9, G1

- **Now:** unique=true on a two-member group passes. When both views have the column, link_feature_pg joins on pg_group_id (dataset_queries.py:216-220), and nothing checks it equals the hash of pg_accessions. A stale pg_group_id loses the link. A null pg_accessions with a leftover pg_group_id still links. anchor_outside_membership is a warning even in strict mode (validation.py:524).
- **Spec:** feature.yaml:161: unique means the group has exactly one protein. feature.yaml:166: pg_group_id is the hash of the sorted pg_accessions. data-model.md:159 says the link is by canonical membership.
- **Reader impact:** Files from non-qpx writers can silently mislink features to protein groups, and unique-peptide filters (mokume) give wrong results.
- **Suggested:** Validate that unique equals (len(pg_accessions)==1), that pg_group_id equals the hash of the membership, and that anchor_protein is in pg_accessions, all as strict errors.
- **Evidence:** `qpx/core/dataset_queries.py:216-220`; `qpx/core/data/validation.py:520-533`; experiment mutate.py f_unique_true_multi_member_group -> exit 0; link_feature_pg on mut/f_shared_peptide_null_group -> 1 link; mut/y_stale_pg_group_id -> 0 links

### VALIDATOR-24: Unlinked shared peptides are 'not an error' in the spec but fail --conformance {#validator-24}

*feature* · contradictory · confirmed · checklist F9, G1

- **Now:** Under --semantic, pg_link_coverage warns when fewer than 95% of features link to a pg row (label_validation.py:268-269,404-411; dataset_validation.py:352-361). Unassigned shared peptides with null pg_accessions, and decoy features with no pg row, count as unlinked. A single shared-peptide feature (0% linked) fails --conformance.
- **Spec:** feature.yaml:153,158 allow null pg_accessions and anchor for shared or unassigned peptides. pg.md:73: 'No match means no matching PG link ... (not an error)'.
- **Reader impact:** Datasets with many shared peptides are rejected by the gate, so producers are pushed to force group assignments.
- **Suggested:** Leave null-membership and decoy features out of the coverage denominator, or make coverage purely informational.
- **Evidence:** `qpx/core/label_validation.py:268-269,404-411`; `qpx/core/dataset_validation.py:352-361`; `experiment mutate3.py y_shared_null_group_and_null_groupid -> WARN pg_link_coverage, exit 1`; `docs/spec/pg.md:73`

### VALIDATOR-27: Feature PEP is defined and validated, but psm.md says PEP exists only on PSMs {#validator-27}

*feature* · contradictory · confirmed · checklist F5, P5

- **Now:** feature.yaml:37-41 defines posterior_error_probability, range-checked to [0,1]. MaxQuant fills it (null for its 104,491 transferred rows). best_id_score appears only in the docs (scores.md:122, psm.md:286); grep finds no implementation.
- **Spec:** psm.md:285-286: PEP is defined only in the PSM view, feature and peptide views do not carry it, and best_id_score should be used.
- **Reader impact:** A reader does not know which identification a feature PEP describes (best PSM, this run, or the donor), or whether it should exist at all.
- **Suggested:** Keep feature PEP, define it as the best direct PSM in this run (null for transferred rows), and remove the best_id_score text.
- **Evidence:** `qpx/core/data/schemas/feature.yaml:37-41`; `docs/spec/psm.md:285-286`; `docs/spec/scores.md:122`; `query audit-run:maxquant/full: transferred rows have null PEP`

### VALIDATOR-M2: RT window clock unspecified; rt outside [rt_start, rt_stop] passes validation {#validator-m2}

*feature* · unspecified · found by the verifier · checklist F10, P7

- **Now:** No check on rt_start<=rt_stop, rt within the window, or rt>=0. Mutations f_rt_outside_window (rt 50, window 100-10) and p_rt_negative pass --conformance. MaxQuant fills rt_start/rt_stop from calibrated RT (feature_adapter.py:580-583) and rt from the run clock. In maxquant/full, 83,413 of 243,821 rows have rt outside their window; rt_aligned falls outside in 0 rows.
- **Spec:** feature.md:65-66: rt is in seconds on the run's own clock and rt_aligned is the aligned time. feature.md:69-70 do not say which clock rt_start/rt_stop use. converter-coverage.md:204 says MaxQuant uses calibrated start/finish.
- **Reader impact:** Readers that integrate or match peaks with [rt_start, rt_stop] against rt or PSM rt mix clocks, so a third of MaxQuant features look inconsistent.
- **Suggested:** Declare rt_start/rt_stop on the same clock as rt (add aligned window fields if needed), and validate rt_start<=rt<=rt_stop and rt>=0.
- **Evidence:** `qpx/converters/maxquant/feature_adapter.py:580-583`; `query audit-run:maxquant/full: rt outside window 83413; rt_aligned outside 0`; `docs/spec/feature.md:65-70`; `docs/spec/converter-coverage.md:204`

### VALIDATOR-M5: id_run_file_name has two definitions and is never validated {#validator-m5}

*feature* · contradictory · found by the verifier · checklist F4, F10

- **Now:** The run referential check covers run_file_name only (dataset_validation.py:244-266), and no qpx/core module reads id_run_file_name except dataset_merge.py:45. Mutation f_idrun_not_in_run_view (id_run_file_name='no_such_run') passes --conformance.
- **Spec:** feature.md:138: the run of the best PSM, which 'may differ from run_file_name' (implying a donor run). feature.md:145 and converter-coverage.md:203: null for a transferred row.
- **Reader impact:** Readers that detect MBR by id_run_file_name != run_file_name (mokume) and readers that look for a null value get different answers. The donor run of a transfer cannot be recovered.
- **Suggested:** Choose one: null for transferred rows, or the donor run with identification_type governing. Validate that a non-null value exists in run.parquet.
- **Evidence:** `qpx/core/dataset_validation.py:244-266`; `docs/spec/feature.md:138,145`; `docs/spec/converter-coverage.md:203`; `mut/f_idrun_not_in_run_view -> exit 0`

## mokume (consumer)

### MOKUME-01: mokume reads `unique` as group-specific; the spec means a single-protein group {#mokume-01}

*feature* · diverges from spec · confirmed · checklist F9

- **Now:** By default mokume keeps only rows with `unique`=1, or with null `unique` and 1 accession (feature.py:224-230, 320-324; lib.rs:8460-8468). It then quantifies the whole group. As a result every peptide of a multi-member group is dropped. Running the mokume Feature reader on the dev MaxQuant full output removed 43,823 of 231,063 measured rows by `unique` alone. All 85 multi-member non-decoy groups (out of 1,063) lost every row.
- **Spec:** feature.yaml:159-161 and feature.md:101: true only when the group has exactly one protein. False for multi-member groups and for razor-assigned shared peptides. No field says whether a peptide maps to this group only.
- **Reader impact:** Multi-protein groups (isoforms, paralogs) silently disappear from mokume protein tables. Users assume that `unique` means proteotypic for its group.
- **Suggested:** Keep `unique` as the single-protein flag. Add a group-specificity flag, e.g. `shared_across_groups` or `n_groups`. Consumers that quantify groups use group-specific peptides.
- **Evidence:** `mokume/python/mokume/io/feature.py:224-230`; `mokume/python/mokume/io/feature.py:320-324`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:8455-8468`; `qpx/core/data/schemas/feature.yaml:159-161`; `docs/spec/feature.md:101`; run: audit-run:mokume/probe1.py on audit-run:maxquant/full/maxquant.feature.parquet -> unique x multi = (F,F)15125,(F,T)28698,(T,F)187240; 85/85 multi-member groups have no unique row

### MOKUME-02: MaxQuant `unique` changes meaning depending on whether proteinGroups.txt was given {#mokume-02}

*feature* · inconsistent across producers · partly confirmed (see correction) · checklist F9, X1

- **Now:** mokume applies the same `unique` filter whatever the producer or converter options (feature.py:227-230; lib.rs:8464).
- **Spec:** tool-mappings.md:62: with proteinGroups.txt it is true only for a one-protein group. Without it, it is true when the peptide maps to one protein (`Proteins`), or, without that column, when it is in a single group of any size.
- **Correction:** The dependency is real, but the data can show which meaning applies. ontology.parquet records the source column of feature pg_accessions/anchor_protein: 'Protein IDs'/'Majority protein IDs' versus 'Leading proteins'/'Leading razor protein' (tool-mappings.md:62; observed in audit-run:maxquant/full vs audit-run:maxquant/simple ontology). That is dataset-level only, not per row (row-level fallback when Protein group IDs do not resolve), and it is not recorded for `unique` itself.
- **Reader impact:** The same MaxQuant search gives different mokume protein sets depending on a converter flag. Nothing in the data tells which meaning applies.
- **Suggested:** `unique` is null whenever group size is unknown, never 'one group of any size'. Record the source column in provenance or ontology.
- **Evidence:** `docs/spec/tool-mappings.md:62`; `mokume/python/mokume/io/feature.py:227-230`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:8464`

### MOKUME-03: Unassigned shared peptides (null pg_accessions) are never quantified, even with keep-shared {#mokume-03}

*feature* · unspecified · confirmed · checklist F9, P8

- **Now:** Python drops rows with null pg_accessions (feature.py:295; stages.py:623; peptide.py:282). Rust requires a non-empty accession list (lib.rs:8467). `--keep-shared-peptides` only skips the `unique` test. So 'shared' in mokume means razor-assigned only. A null `unique` with 1 accession is treated as unique, a legacy qpx<1.1.3 fallback (feature.py:224-226; lib.rs:8462-8464).
- **Spec:** feature.yaml:154-161: anchor, pg and unique are null when the evidence spans more than one group, and 'intensities are kept either way'. pg_positions lists every evidence protein (feature.yaml:183-186). pepmap holds all proteins (producer-guide.md:16). Nothing says how consumers should use these rows.
- **Reader impact:** Unassigned shared evidence is silently lost. Consumers cannot tell razor-assigned from unassigned shared peptides without pepmap.
- **Suggested:** State whether and how consumers may use features with null pg_accessions (via pepmap or pg_positions). Name the assigned-vs-unassigned distinction explicitly.
- **Evidence:** `mokume/python/mokume/io/feature.py:292-296`; `mokume/python/mokume/pipeline/stages.py:623`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:8467`; `qpx/core/data/schemas/feature.yaml:151-161`; `docs/spec/producer-guide.md:16`

### MOKUME-04: Protein key is pg_accessions in stored order; the spec treats membership as a set {#mokume-04}

*feature* · unspecified · confirmed · checklist F9, G1

- **Now:** Python joins pg_accessions in stored order (aggregation.py:167; stages.py:85-100). Rust does the same (lib.rs:8830-8838), and the Ratio method uses the first accession (lib.rs:1876-1883, 8840-8846). The spectral-count path sorts instead (spectral_count.rs:116-118, 241-245). pg_group_id is never read. In dev outputs, 14,382 of 30,623 multi-member MaxQuant rows and 54 of 72 CDAP rows are unsorted.
- **Spec:** pg.md:414: accessions are compared as sets; order does not affect identity. feature.yaml:162-166: pg_group_id is a hash of the sorted list. No rule fixes the stored order or says that element 1 is the anchor.
- **Reader impact:** One group stored in two orders (merged datasets, two producers) becomes two proteins. Intensity and spectral-count outputs name the same group differently.
- **Suggested:** Consumers key groups on pg_group_id or the sorted set. If order carries meaning, the spec fixes it (e.g. anchor first, then sorted).
- **Evidence:** `mokume/python/mokume/preprocessing/aggregation.py:167`; `mokume/python/mokume/pipeline/stages.py:85-100`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:8830-8846`; `mokume/rust/crates/mokume-pipeline/src/spectral_count.rs:116-118,241-245`; `docs/spec/pg.md:414`; `qpx/core/data/schemas/feature.yaml:162-166`; `duckdb on audit-run:maxquant/full and audit-run:mzid_cdap/cdap feature.parquet`

### MOKUME-06: Python ignores identification_type; Rust keeps transferred features by default {#mokume-06}

*feature* · unspecified · confirmed · checklist F4, D1

- **Now:** The Python reader never reads identification_type or id_run_file_name (feature.py:568-603; MBR is mentioned only in comments at 606 and 626), so MBR quantities always enter rollups. Rust offers `--transferred-features drop` as an opt-in; the default is keep (features_to_proteins.rs:83-90; config.rs:172-180). mokume#117 reports that MBR rows inflate MaxLFQ false positives.
- **Spec:** feature.yaml:210-214 defines direct, transferred and secondary, with null meaning the tool does not say. Nothing says whether a transferred quantity is equivalent to a direct one or how consumers should treat it.
- **Reader impact:** Protein results depend on which mokume implementation runs and on an undocumented default. MBR signal is treated as identified.
- **Suggested:** State what a transferred quantity is: MS1 signal measured in this run, identity borrowed from another. Require identification_type from MBR-capable producers. Recommend that consumers expose a policy.
- **Evidence:** `mokume/python/mokume/io/feature.py:568-603`; `mokume/rust/crates/mokume-command/src/features_to_proteins.rs:83-90`; `mokume/rust/crates/mokume-core/src/config.rs:166-180`; `qpx/core/data/schemas/feature.yaml:210-214`; `github_history.md: mokume#117`

### MOKUME-07: mokume infers 'transferred' from a null id_run_file_name, but producers use null differently {#mokume-07}

*feature* · inconsistent across producers · partly confirmed (see correction) · checklist F4, F10, D1

- **Now:** In Rust an explicit type wins. Otherwise a null id_run_file_name, or one naming another run, means transferred (qpx.rs:493-524; mokume docs features2proteins.md:189-193). Producers use null differently: MaxQuant only for MULTI-MATCH; FragPipe for every multi-run experiment, even MS/MS rows; OpenMS for conflicting or spectrum-less identifications, with a null type. quantms MSstats always writes its own run and no type.
- **Spec:** feature.yaml:207-209: id_run_file_name is the 'Run file containing best PSM'. What null means is not stated. producer-guide.md:26: an MBR row has an empty scan and type transferred.
- **Correction:** The spec does say what a null means. feature.md:145: id_run_file_name 'can be null for other reasons, so filter on this field' (identification_type). converter-coverage.md:57,203 documents each producer's null rule. So mokume's fallback diverges from the spec rather than filling a gap. FragPipe fills identification_type from Match Type (fragpipe/feature_adapter.py:51,361), so its multi-run direct rows reach the heuristic only when Match Type columns are absent. The OpenMS and quantms points hold.
- **Reader impact:** Where the type is null, FragPipe fractionated direct features and OpenMS conflicting-ID features are dropped as transfers. quantms requantified features are kept as direct.
- **Suggested:** A null id_run_file_name means unknown, never transferred. Transfer is signalled only by identification_type. Define id_run_file_name for multi-run units.
- **Evidence:** `mokume/rust/crates/mokume-io/src/qpx.rs:493-524`; `mokume/docs/user-guide/features2proteins.md:189-193`; `qpx/converters/maxquant/feature_adapter.py:655-661`; `qpx/converters/fragpipe/feature_adapter.py:337-341,361`; `qpx/converters/openms_consensus/feature_adapter.py:532-547`; `qpx/converters/quantms_msstats/feature_adapter.py:525`; `docs/spec/converter-coverage.md:55`

### MOKUME-08: Feature row cardinality contradicts itself in the spec; mokume assumes duplicates are redundant {#mokume-08}

*feature* · contradictory · partly confirmed (see correction) · checklist F1, F2, F7

- **Now:** Python keeps the maximum row per (peptidoform, charge, sample) (aggregation.py:295-322; stages.py:128-140), and Rust does the same (lib.rs:2761-2772). For TMT, Rust best-row keeps the single row with the highest reporter total per (peptidoform, charge, plex) (tmt.rs:116-143, 241-250). In the dev MaxQuant full output, 14,086 (peptidoform, charge, run) triples have more than one row.
- **Spec:** producer-guide.md:12 says one row per precursor, charge and run. But the default identity composite adds rt (feature.yaml:4), and OpenMS adds scan, observed_mz and consensus_rt (feature_dedup.py:16), so several rows are legal. Whether they are alternatives or parts is not said.
- **Correction:** The contradiction is real, and data-model.md:64-67,73-74 states it more strongly: a feature is 'one quantified chromatographic peak'/'per upstream Feature entity', against producer-guide.md:12. Detail is wrong: Rust TMT best-row is opt-in; the default --tmt-row-merge is max per channel (features_to_proteins.rs:~313-316). FragPipe adds another source of extra rows: compensation_voltage is in its composite (feature.md:25).
- **Reader impact:** MAX, SUM and best-row give different quantities. A reader cannot tell whether extra rows are redundant or additive.
- **Suggested:** Define the feature row granularity per acquisition type. When several rows share (peptidoform, charge, run, label), declare whether they are alternatives or parts.
- **Evidence:** `docs/spec/producer-guide.md:12`; `qpx/core/data/schemas/feature.yaml:4`; `qpx/converters/openms_consensus/feature_dedup.py:16`; `mokume/python/mokume/preprocessing/aggregation.py:295-322`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:2761-2772`; `mokume/rust/crates/mokume-pipeline/src/tmt.rs:116-143,241-250`; `duckdb on audit-run:maxquant/full/maxquant.feature.parquet: 14086 duplicate triples`

### MOKUME-09: Direct and transferred rows coexist for one ion and run; mokume's MAX can pick the MBR value {#mokume-09}

*feature* · unspecified · confirmed · checklist F2, F4

- **Now:** In the dev MaxQuant full output, 9,272 (peptidoform, charge, run) triples have both a direct and a transferred measured row. In 799 of them the transferred intensity is larger, so mokume's MAX (aggregation.py:316-320; lib.rs:2761-2772) takes the MBR quantity even though the precursor was identified in its own run. This happens under the Rust default (keep) and always in Python.
- **Spec:** Nothing forbids a transferred row when a direct row exists for the same precursor and run, and nothing says which quantity represents the precursor (feature.yaml:210-214; producer-guide.md:12).
- **Reader impact:** The quantity of an identified precursor may come from a different, matched peak.
- **Suggested:** Allow at most one quantified row per (peptidoform, charge, run, label), or require a flag naming the representative row. Producers drop or flag transferred rows superseded by direct ones.
- **Evidence:** `duckdb on audit-run:maxquant/full/maxquant.feature.parquet: 9272 direct+transferred triples, 799 with transferred>direct`; `mokume/python/mokume/preprocessing/aggregation.py:316-320`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:2761-2772`

### MOKUME-11: quantification_unit_id is ignored; multi-run quantities land on the first raw run {#mokume-11}

*feature* · ambiguous · partly confirmed (see correction) · checklist F10, G1

- **Now:** mokume keys runs and samples on run_file_name only (feature.py:679-686; lib.rs:8872-8881) and never reads quantification_unit_id. FragPipe LFQ experiments and FragPipe TMT plexes write run_file_name as the first raw run and quantification_unit_id as the experiment or prefix (fragpipe/feature_adapter.py:337-352; converter-coverage.md:186, 206-211).
- **Spec:** feature.yaml:90-94: the unit is 'used when one quantity can aggregate multiple raw runs'. Which run_file_name such a row carries is not specified.
- **Correction:** For FragPipe the run is specified: converter-coverage.md:213-216 says run_file_name is the experiment's first raw run, a member of the pg grouped_runs. converter-coverage.md:206-210 says the same for TMT plexes. The generic rule is still missing from feature.yaml:90-94 and feature.md:90, and mokume never reads quantification_unit_id (grep finds none).
- **Reader impact:** Run-level output, run QC and run normalization assign a quantity summed over fractions to fraction 1. The other fractions look empty.
- **Suggested:** When quantification_unit_id spans runs, readers key on the unit. Define what run_file_name means on such rows (e.g. the first run, not a measurement run).
- **Evidence:** `qpx/core/data/schemas/feature.yaml:90-94`; `docs/spec/converter-coverage.md:186,206-211`; `qpx/converters/fragpipe/feature_adapter.py:337-352`; `mokume/python/mokume/io/feature.py:679-686`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:8872-8881`

### MOKUME-14: mokume accepts the legacy 'label = run file' LFQ layout, which the spec rejects {#mokume-14}

*feature* · diverges from spec · confirmed · checklist F4, F1

- **Now:** Python decides for the whole file that labels are runs when every positive label matches a run, and then re-keys sample and run on the label (feature.py:605-643, 679-686). Rust decides per row, when any label equals that row's run (qpx.rs:528-539, 548-557). Rows that list several runs block the id_run fallback (qpx.rs:510-521).
- **Spec:** intensities.md:53 and 136: an LFQ label is 'LFQ', with at most one entry per row. The validator raises the error ambiguous_label_sample (label_validation.py:192-199). Every dev converter writes 'LFQ'.
- **Reader impact:** Python and Rust can classify the same file differently, and mokume accepts files that QPX rejects.
- **Suggested:** One run per feature row; a label is never a run name. Keep any legacy reading behind an explicit flag.
- **Evidence:** `mokume/python/mokume/io/feature.py:605-643`; `mokume/rust/crates/mokume-io/src/qpx.rs:510-557`; `docs/spec/intensities.md:53,136`; `qpx/core/label_validation.py:192-199`

### MOKUME-15: Python pandas paths ignore is_decoy and only drop accessions containing 'DECOY' {#mokume-15}

*feature* · diverges from spec · confirmed · checklist P4, F8

- **Now:** TMT and iTRAQ in features2proteins (stages.py:436-438, 604-648) and all of features2peptides (peptide.py:281-304) read rows without the SQL filter (feature.py:1086-1103). They drop decoys only when an accession contains DECOY (contaminants.py:8, 38-39). A probe on the dev MaxQuant output: all 255 measured is_decoy=true (REV__) rows survive. The SQL path (feature.py:305-306) and Rust (lib.rs:8461) remove them.
- **Spec:** pg.md:46-49: is_decoy states the status; do not parse the prefix. producer-guide.md:33: REV__ and DECOY_ markers are kept on the accession.
- **Reader impact:** MaxQuant (REV__) and FragPipe (rev_) decoys enter mokume's peptide tables and TMT protein tables.
- **Suggested:** is_decoy is the only decoy signal, and consumers must filter on it.
- **Evidence:** `run: audit-run:mokume/probe3.py -> legacy pandas path 255/255 decoys kept; SQL path 0`; `mokume/python/mokume/pipeline/stages.py:436-438,604-648`; `mokume/python/mokume/normalization/peptide.py:281-304`; `mokume/python/mokume/core/contaminants.py:8,15-42`; `mokume/python/mokume/io/feature.py:303-309,1086-1103`; `docs/spec/pg.md:46-49`; `docs/spec/producer-guide.md:33`

### MOKUME-16: Null is_decoy: Python drops the row, Rust and TissueMap keep it {#mokume-16}

*feature* · specified, worth keeping · confirmed · checklist P4, F8

- **Now:** The Python SQL test `is_decoy = false` excludes NULL (feature.py:305-306). Rust's unwrap_or(false) keeps the row (lib.rs:8461; feature_mapping.rs:125-129; psm.rs:118-122). TissueMap's COALESCE(is_decoy, FALSE) keeps it (loader.py:56).
- **Spec:** feature.yaml:63-66: is_decoy is required and therefore non-nullable (loader.py:148-150). Nulls are only a warning unless validation is strict (schema.py:580, 603, 621-623).
- **Reader impact:** Files that are invalid but still accepted give results that depend on the mokume implementation.
- **Suggested:** Keep a non-null is_decoy in the contract and make a null a hard error.
- **Evidence:** `qpx/core/data/schemas/feature.yaml:63-66`; `qpx/core/data/loader.py:148-150`; `qpx/core/data/schema.py:603,621-623`; `mokume/python/mokume/io/feature.py:305-306`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:8461`; `mokume/python/mokume/tissuemap/loader.py:56`

### MOKUME-17: Contaminants are found by accession substring; the spec flags them only on pg {#mokume-17}

*feature* · unspecified · partly confirmed (see correction) · checklist P4, P8, G3

- **Now:** mokume removes groups whose members all contain CONTAMINANT, CONTAM_ or ENTRAP: case-sensitively in Python (feature.py:185-188, 292-317; contaminants.py) and case-insensitively in Rust (lib.rs:8902-8909, 8954-8990). It never reads pg.contaminant. Markers differ by producer: MaxQuant keeps CON__ verbatim (utils.py:66-91), FragPipe writes contam_ (pg.md:365) and quantms writes Cont_ (utils.py:67-68).
- **Spec:** producer-guide.md:34: feature rows have is_decoy false, and `contaminant` is true only on pg. pg.md:46-49: do not parse the prefix. The feature view has no contaminant flag.
- **Correction:** On dev, the FragPipe converter strips the contam_ marker: protein_accessions keeps only decoy prefixes (fragpipe/constants.py:28,42-58). Observed: 'contam_sp|P02768|ALBU_HUMAN' becomes 'P02768' in audit-run:fragpipe/out feature and pg. So neither Python nor Rust can remove FragPipe contaminants; the 'removed by Rust' claim is wrong. Only pg.contaminant=true survives. MaxQuant CON__ is kept verbatim and is matched by neither, as claimed.
- **Reader impact:** MaxQuant contaminants are never removed. FragPipe contaminants are removed by Rust but not by Python. Results depend on FASTA naming conventions.
- **Suggested:** Add a feature-level contaminant flag, or a normalized marker, so feature-only consumers need not parse accessions.
- **Evidence:** `mokume/python/mokume/io/feature.py:185-188,292-317`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:8902-8990`; `qpx/converters/utils.py:66-91`; `docs/spec/producer-guide.md:34`; `docs/spec/pg.md:46-49,95,365`

### MOKUME-18: mokume treats measured zeros as missing values {#mokume-18}

*feature* · diverges from spec · partly confirmed (see correction) · checklist D6, F6

- **Now:** Every path keeps only intensity > 0 (feature.py:211-212, 749; aggregation.py:201; lib.rs:8465-8466; loader.py:58). Zeros are therefore indistinguishable from missing values and can later be imputed. Null entries (identified but not quantified) are dropped as well.
- **Spec:** intensities.md:121-127 and producer-guide.md:45-47: 0 is a measured zero, null means identified but not quantified, and no entry means not identified. The spec gives consumers no guidance.
- **Correction:** The code behavior is real, but the practical impact is theoretical. Across all ~40 converter outputs in work/, primary intensity 0 occurs only in hand-made infra test files; converters write sentinel zeros as null. The example is wrong: an empty reporter channel of an identified plex is a null entry by spec (intensities.md:123-124), not 0.
- **Reader impact:** mokume merges measured-zero and missing cells, losing a distinction the spec deliberately creates (e.g. empty TMT channels).
- **Suggested:** State how consumers should treat a 0 (censored, below detection) versus null, with a recommended handling.
- **Evidence:** `mokume/python/mokume/io/feature.py:211-212,749`; `mokume/python/mokume/preprocessing/aggregation.py:201`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:8465-8466`; `docs/spec/intensities.md:121-127`; `docs/spec/producer-guide.md:45-47`

### MOKUME-19: Imputed entries are dropped via imputed_labels (keep this); producer-guide wording names no column {#mokume-19}

*feature* · specified, worth keeping · confirmed · checklist D6, D7

- **Now:** Python measured_intensity_clause (io/_qpx.py:6-10; feature.py:750) and Rust drop_imputed_entries (qpx.rs:636-652) drop entries whose label is in imputed_labels. pg.is_imputed is never read.
- **Spec:** feature.yaml:235-239 and intensities.md:131-133: an imputed value is an ordinary entry flagged in imputed_labels (is_imputed on pg). But producer-guide.md:48 says 'with `imputation` flagged', and no such column exists.
- **Reader impact:** Correct for features. The producer-guide wording could lead producers to invent an `imputation` column that mokume would never see.
- **Suggested:** Keep imputed_labels and is_imputed. Fix producer-guide.md:48 to name them.
- **Evidence:** `mokume/python/mokume/io/_qpx.py:6-10`; `mokume/rust/crates/mokume-io/src/qpx.rs:636-652`; `qpx/core/data/schemas/feature.yaml:235-239`; `docs/spec/intensities.md:131-133`; `docs/spec/producer-guide.md:48`

### MOKUME-20: 'transferred' means a matched MS1 signal (MaxQuant) or an imputed value (Spectronaut) {#mokume-20}

*feature* · inconsistent across producers · confirmed · checklist D1, F4, X1

- **Now:** Spectronaut writes the rows marked EG.IsImputed both as identification_type=transferred and in imputed_labels (spectronaut/feature_adapter.py:336-343). mokume always drops imputed entries, so Spectronaut transfers disappear even with `--transferred-features keep`, while MaxQuant MULTI-MATCH quantities are kept. Spectronaut is a DIA tool; noted for the contract.
- **Spec:** feature.yaml:214: transferred means 'matched from another run'. feature.yaml:235-239: imputed_labels marks values the tool imputed instead of measuring. The overlap between the two is not addressed (converter-coverage.md:199, 201).
- **Reader impact:** The same flag yields a measured value for one producer and an imputed value for another.
- **Suggested:** Transferred implies the value was measured in this run. Imputed values always go in imputed_labels. A row may carry both only when its value is imputed.
- **Evidence:** `qpx/converters/spectronaut/feature_adapter.py:336-343`; `docs/spec/converter-coverage.md:199,201`; `mokume/python/mokume/io/_qpx.py:6-10`; `mokume/rust/crates/mokume-io/src/qpx.rs:636-652`

### MOKUME-21: Peptide FDR filter reads only peptide_qvalue, which most DDA producers leave null {#mokume-21}

*feature* · inconsistent across producers · confirmed · checklist D4, F5, P5

- **Now:** PeptideFDRFilter and the Rust filter read only peptide_qvalue and drop rows where it is null (filters/peptide.py:387-412; lib.rs:567-574, 6709-6713). In dev outputs, MaxQuant, FragPipe and PD leave peptide_qvalue null; OpenMS and CDAP fill it. CDAP and PD fill precursor_qvalue instead.
- **Spec:** feature.yaml:42-48: optional. feature.md:124-132: peptide_qvalue only for true peptide-level q-values. The precursor and global typed columns are separate (scores.md:40-45).
- **Reader impact:** mokume's FDR filtering works for some producers only. Users may switch to another FDR level without noticing.
- **Suggested:** Declare which typed q-values each DDA view guarantees, and the FDR level at which the view was already filtered.
- **Evidence:** `mokume/python/mokume/preprocessing/filters/peptide.py:387-412`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:567-574,6709-6713`; `survey (duckdb counts) of work/{maxquant/full,fragpipe/out,pd/cases/*,openms/*,mzid_cdap/cdap}/*.feature.parquet`; `docs/spec/feature.md:124-132`; `docs/spec/scores.md:40-45`

### MOKUME-22: Score filter on precursor/global q-value uses the column only without nulls; fallback ignores aliases {#mokume-22}

*feature* · ambiguous · confirmed · checklist D4, P5

- **Now:** mokume reads the top-level precursor_qvalue or global_qvalue only when parquet stats show zero nulls in every row group (io/_qpx.py:13-21; qpx.rs:208-219). Otherwise it looks up the exact name in additional_scores (io/_qpx.py:24-46; feature.py:481-520; qpx.rs:585-611), and errors or drops the rows where that entry is missing.
- **Spec:** scores.md:38-47: typed columns are filled from aliases (qvalue, q_value and psm_level_qvalue all fill precursor_qvalue) and stay null when the source has none. derived.py:28-33 copies the value, so the original score keeps its source name.
- **Reader impact:** Partial nulls make mokume fail or filter on another source. precursor_qvalue may silently be a PSM-level q-value.
- **Suggested:** The typed column is authoritative, and null means unknown, with a stated pass/fail policy. Do not fill precursor_qvalue from PSM-level score names.
- **Evidence:** `mokume/python/mokume/io/_qpx.py:13-46`; `mokume/rust/crates/mokume-io/src/qpx.rs:208-219,585-611`; `docs/spec/scores.md:38-47`; `qpx/writers/derived.py:28-33`

### MOKUME-26: Reporter isotope-correction state is not in the data; mokume may correct twice {#mokume-26}

*feature* · ambiguous · partly confirmed (see correction) · checklist F6

- **Now:** The optional `--tmt-impurity` corrects feature.intensities on the assumption that they are raw reporter intensities (tmt.rs:88-114; features_to_proteins.rs:320-324). mokume never reads additional_intensities or provenance.
- **Spec:** intensities.md:213: correction state is recorded only as the provenance parameter isotope_corrected. MaxQuant writes raw values in intensities and corrected values in additional_intensities as reporter_intensity_corrected (converter-coverage.md:195; intensities.md:188). OpenMS records IsobaricAnalyzer's isotope_correction (openms_consensus/provenance.py:66).
- **Correction:** The spec does say what goes in `intensities`: the 'raw reporter ion signal from TMT channels' (intensities.md:219). Yet intensities.md:213 allows the primary values to be isotope-corrected and records that state only in provenance (OpenMS isotope_correction gives isotope_corrected, openms_consensus/provenance.py:66). The spec contradicts itself, and the risk of double correction in mokume holds.
- **Reader impact:** `intensities` hold raw reporter values for one producer and corrected values for another. Consumers cannot tell without reading provenance.
- **Suggested:** Fix whether feature.intensities hold uncorrected or corrected reporter values, or add a footer flag. Corrected values go in a named additional intensity.
- **Evidence:** `mokume/rust/crates/mokume-pipeline/src/tmt.rs:88-114`; `mokume/rust/crates/mokume-command/src/features_to_proteins.rs:320-324`; `docs/spec/intensities.md:188,213`; `docs/spec/converter-coverage.md:195`; `qpx/converters/openms_consensus/provenance.py:66`

### MOKUME-27: Python mokume rejects MS1 labels and lacks TMTpro-18 and deuterated TMT labels {#mokume-27}

*feature* · diverges from spec · confirmed · checklist D7, X1

- **Now:** QuantificationCategory.classify accepts only LFQ, TMT and iTRAQ and raises an error for SILAC, dimethyl and mTRAQ (labeling.py:84-117). The plex is inferred from the label count (labeling.py:91-105). Channel maps lack TMT134C, TMT135N and the *D labels (labeling.py:235-316), which become a NaN channel (aggregation.py:171). identified_labels is ignored.
- **Spec:** intensities.md:71: the vocabulary includes TMT134C, TMT135N and TMT127D, among others. intensities.md:79-110: an MS1-label multiplet is one row with several labels. feature.yaml:144-148: identified_labels. run.label_definitions (run.md:20).
- **Reader impact:** Valid spec 1.3 SILAC or TMTpro-18 files fail in Python mokume or are mis-mapped.
- **Suggested:** Take the label family and plex from run.samples and run.label_definitions, not from label counts. Consumers support the full canonical vocabulary.
- **Evidence:** `mokume/python/mokume/model/labeling.py:84-117,235-316`; `mokume/python/mokume/preprocessing/aggregation.py:171`; `docs/spec/intensities.md:71,79-110`; `qpx/core/data/schemas/feature.yaml:144-148`

### MOKUME-30: identification_type is compared as an exact string; off-vocabulary values count as identified {#mokume-30}

*feature* · ambiguous · partly confirmed (see correction) · checklist D1, F3

- **Now:** Rust treats a row as transferred only when the value equals 'transferred' (qpx.rs:501-504). Values such as 'MBR' or 'unidentified', both seen in audit-run:infra/fb.feature.parquet, count as identified. 'secondary' rows (MaxQuant MULTI-SECPEP, 5,109 rows) are kept as own-run identifications (mokume docs features2proteins.md:189-191).
- **Spec:** feature.yaml:210-214: closed set direct, transferred, secondary. It is enforced only as a warning unless validation is strict (data/validation.py:70-72; schema.py:603, 621). There is no value for unidentified features (GitHub #378, open). Whether a 'secondary' quantity is independent is not stated.
- **Correction:** 'MBR' and 'unidentified' do not come from a producer: audit-run:infra/fb.feature.parquet is a hand-made test file. It does show that the writer accepts off-vocabulary values with only a warning (validation.py:70-72; schema.py severity unless strict). The spec admits unidentified features exist (feature.md:8) but gives them no identification_type value. The rest holds.
- **Reader impact:** Off-vocabulary producers silently bypass MBR handling. Unidentified features, if ever emitted, would be treated as identified.
- **Suggested:** Make the closed enum a hard error. Add an explicit `unidentified` value or forbid unidentified rows. Define what `secondary` quantities mean.
- **Evidence:** `mokume/rust/crates/mokume-io/src/qpx.rs:501-504`; `qpx/core/data/validation.py:70-72`; `qpx/core/data/schema.py:603,621`; `qpx/core/data/schemas/feature.yaml:210-214`; `qpx/converters/maxquant/evidence_fields.py:10-17`; `github_history.md: #378`

### MOKUME-M1: OpenMS feature.peptide_qvalue is a PSM-level q-value; the peptide FDR filter mixes levels {#mokume-m1}

*feature* · diverges from spec · found by the verifier · checklist D4, F5, P5

- **Now:** Without Percolator peptide-level FDR, OpenMS writes the first identification's primary score of the run into peptide_qvalue whenever the score type is a q-value (openms_consensus/feature_adapter.py:149-160,325-345,736). In ProteomicsLFQ_1, 68/68 feature peptide_qvalue = the linked psm.qvalue, and the fixture has no peptide-level FDR. CDAP fills precursor_qvalue from psm_level_qvalue (94/94). mokume's PeptideFDRFilter reads peptide_qvalue as peptide-level.
- **Spec:** feature.md:125-132 and scores.md:43: peptide_qvalue is 'Peptide, as computed by the source', only for true peptide-level q-values (#284). converter-coverage.md:59 documents OpenMS reading any score 'declared to be a q-value/FDR', with no level check.
- **Reader impact:** mokume --filter-peptide-fdr filters OpenMS data at PSM level and silently believes it is peptide level.
- **Suggested:** A typed q-value column is filled only when the source states the level. OpenMS PSM-level q-values go to precursor_qvalue or stay in additional_scores.
- **Evidence:** `qpx/converters/openms_consensus/feature_adapter.py:149-160,325-345,736`; `docs/spec/converter-coverage.md:59`; `docs/spec/feature.md:125-132`; `docs/spec/scores.md:43-44`; `duckdb audit-run:openms/ProteomicsLFQ_1_out: 68/68 peptide_qvalue == psm.qvalue; consensusXML score_type=q-value, no peptide-level flag`; `mokume/python/mokume/preprocessing/filters/peptide.py:385-412`

### MOKUME-M2: FragPipe drops the contam_ marker from accessions; feature readers cannot see contaminants {#mokume-m2}

*feature* · diverges from spec · found by the verifier · checklist P4, G3, F9

- **Now:** fragpipe/constants.py:42-58 keeps only decoy prefixes, so 'contam_sp|P02768|ALBU_HUMAN' becomes 'P02768' in feature and pg (audit-run:fragpipe/out). pg.contaminant=true is computed from the raw value (pg_adapter.py:469). Feature rows carry no contaminant signal. mokume (feature-only, substring rules) cannot remove them.
- **Spec:** pg.md:46-49 and producer-guide.md:34: a contaminant marker stays on the accession (CONTAM_P02768), and contaminant=true on pg. pg.md:365: FragPipe keeps contam_ accessions.
- **Reader impact:** FragPipe contaminants enter every mokume protein table, and they cannot be told apart from target proteins using feature.parquet alone.
- **Suggested:** Converters keep the marker as the spec says, and the feature view gets a contaminant flag (see MOKUME-17).
- **Evidence:** `qpx/converters/fragpipe/constants.py:28,42-58`; `qpx/converters/fragpipe/pg_adapter.py:469`; duckdb audit-run:fragpipe/out fp.feature/fp.pg: anchor 'P02768', pg contaminant=true; input combined_protein.tsv 'contam_sp|P02768|ALBU_HUMAN'; `docs/spec/pg.md:46-49,365`; `docs/spec/producer-guide.md:34`

### MOKUME-M3: Any mokume confidence filter silently drops all transferred rows; spec silent on their scores {#mokume-m3}

*feature* · unspecified · found by the verifier · checklist F5, F4, D1

- **Now:** Transferred rows carry no PEP or scores: MaxQuant 0/104,491 PEP and null additional_scores; FragPipe 0/2 PEP. mokume drops rows without the selected score: SQL EXISTS (feature.py:267), Rust lib.rs:6967-6969, and pandas NaN <= t (filters/peptide.py:409). So --filter-score, --filter-peptide-fdr and the like drop all MBR rows, even with --transferred-features keep.
- **Spec:** feature.yaml:215-224 gives only match_score and match_qvalue to transferred rows. Nothing says whether PEP and typed q-values on a transferred row are null, the donor's, or the match's. converter-coverage.md:180,184 shows MaxQuant MBR rows with null PEP and scores. #378 (open) proposes null PEP.
- **Reader impact:** Users who apply a confidence filter unknowingly also apply an MBR policy, so results differ from what --transferred-features suggests.
- **Suggested:** State that identification confidence fields are null on transferred rows and that match_* fields carry the transfer confidence. Consumers should apply the MBR policy separately from score filters.
- **Evidence:** `duckdb audit-run:maxquant/full: transferred 104491 rows, 0 PEP, 0 additional_scores, 104491 match_score`; `duckdb audit-run:fragpipe/out: transferred 2 rows, 0 PEP`; `mokume/python/mokume/io/feature.py:267`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:6966-6972`; `qpx/core/data/schemas/feature.yaml:215-224`; `docs/spec/converter-coverage.md:180,184`

### MOKUME-M5: MaxQuant MULTI-MATCH-MSMS rows: type 'transferred' but id_run_file_name and scan point to own run {#mokume-m5}

*feature* · diverges from spec · found by the verifier · checklist D1, F4, F10

- **Now:** evidence_fields.py:15 maps MULTI-MATCH-MSMS to transferred. feature_adapter.py:656-657 nulls id_run_file_name only for MULTI-MATCH, so these rows keep id_run_file_name = their own run and a scan (feature_adapter.py:560-566). Seen in code only; the fixture has no such rows.
- **Spec:** feature.md:145: id_run_file_name is null for a transferred row. producer-guide.md:26: a transferred feature has an empty scan, because 'No spectrum of this run identified it'. feature.yaml:207-209: id_run_file_name = run file containing the best PSM.
- **Reader impact:** Readers that use id_run_file_name or scan to detect transfers (the mokume docs fallback, features2proteins.md:191-193) and readers that use identification_type disagree on these rows.
- **Suggested:** A transferred row always has a null id_run_file_name and an empty identifying scan. A fragmentation scan without identification goes elsewhere, e.g. quant_scan or cv_params.
- **Evidence:** `qpx/converters/maxquant/evidence_fields.py:8-17`; `qpx/converters/maxquant/feature_adapter.py:555-566,655-661`; `docs/spec/feature.md:145`; `docs/spec/producer-guide.md:26`; `docs/spec/converter-coverage.md:199,203`

## External tools and standards

### EXTERNAL-10: MaxQuant MULTI-MATCH-MSMS is undocumented; QPX says transferred but keeps id_run_file_name {#external-10}

*feature* · diverges from spec · confirmed · checklist D1, F4, F10

- **Now:** The coxdocs Type list defines MSMS, ISO-MSMS, MULTI-MSMS, MULTI-SECPEP and MULTI-MATCH; MULTI-MATCH-MSMS is not defined. QPX maps it to 'transferred' (maxquant/evidence_fields.py:8-16, with a comment that it 'was fragmented too'). id_run_file_name is nulled only for 'MULTI-MATCH' (maxquant/feature_adapter.py:655-657), so MULTI-MATCH-MSMS rows keep their own run and possibly a scan.
- **Spec:** 'id_run_file_name is null for a transferred row' (feature.md:145). producer-guide.md:26: a transferred feature has an empty scan.
- **Reader impact:** For MULTI-MATCH-MSMS rows, identification_type=transferred contradicts the non-null id_run_file_name (and scan). Readers using either rule get different MBR counts.
- **Suggested:** Decide MULTI-MATCH-MSMS semantics: probably 'transferred, MS/MS present but not identifying'. Make id_run_file_name and scan follow identification_type in every converter, and add a validator check.
- **Evidence:** `https://cox-labs.github.io/coxdocs/output_tables.html (evidence Type)`; `qpx/converters/maxquant/evidence_fields.py:8-16`; `qpx/converters/maxquant/feature_adapter.py:655-657`; `docs/spec/feature.md:145`; `docs/spec/producer-guide.md:26`

### EXTERNAL-11: MaxQuant 'Match score' is the donor's identification score, not a match-quality score {#external-11}

*feature* · ambiguous · confirmed · checklist D2, F5

- **Now:** coxdocs: 'Match score: The andromeda score of the MS/MS identification that is the source of this identification by matching between runs'. Match time and m/z difference are relative to the source feature, and 'Match q-value' is for MBR features. QPX MaxQuant run: 1,469 transferred rows, all with match_score; match_qvalue is null on all of them (the column is empty in this export), and PEP is null.
- **Spec:** match_score, match_qvalue: 'The tool's score and q-value for the match of a transferred identification' (feature.md:146; feature.yaml:215-224). The spec does not say whether the score rates the donor identification or the transfer.
- **Reader impact:** Readers may rank or filter transfers by match_score believing it measures transfer correctness. For MaxQuant it only says how good the donor PSM was.
- **Suggested:** Split the fields: donor_score / donor_pep (the identification transferred) and match_score / match_qvalue / match_pep (the transfer's own error estimate). Document the MaxQuant mapping explicitly.
- **Evidence:** `https://cox-labs.github.io/coxdocs/output_tables.html (Match score, Match q-value)`; `qpx/converters/maxquant/evidence_fields.py:20-30`; `run: audit-run:external/mq/out/maxquant.feature.parquet (transferred 1469: match_score 1469, match_qvalue 0, PEP 0)`; `docs/spec/feature.md:146`

### EXTERNAL-12: IonQuant MBR has ion/peptide/protein FDR and a probability; QPX FragPipe fills no match fields {#external-12}

*feature* · inconsistent across producers · confirmed · checklist D2, F4

- **Now:** IonQuant: --ionfdr (transferred-ion FDR, default 0.01) and --peptidefdr/--proteinfdr (default 1, i.e. off). A mixture-model posterior per transfer uses decoy transfers (m/z-shifted). The MBR ion probability threshold is always >=0.5 (1.10.12). Up to 10 donor runs (--mbrtoprun), RT tolerance 1 min. A per-experiment 'Match Type' column (MS/MS|MBR) was added in 1.9.8. QPX reads Match Type only (fragpipe/feature_adapter.py:50-51,339-340).
- **Spec:** match_score/match_qvalue/match_rt_delta/match_mz_delta are 'Not provided' for FragPipe (converter-coverage.md:200). The match-error statistic is limited to a score and a q-value (feature.md:146-147), with no probability or PEP and no level (ion vs peptide).
- **Reader impact:** Transferred rows from MaxQuant and FragPipe look alike but are controlled at different levels (ion vs none vs q-value). Readers cannot apply one MBR confidence filter.
- **Suggested:** Add a match_pep/match_probability option and record the MBR FDR level and threshold (ion/peptide/protein) in provenance (match_between_runs step, provenance.md:187-188). Read IonQuant's probability columns where they exist.
- **Evidence:** `https://github.com/Nesvilab/IonQuant (parameters)`; `https://pmc.ncbi.nlm.nih.gov/articles/PMC8131922 (IonQuant MBR FDR method)`; `https://ionquant.nesvilab.org/CHANGELOG.html`; `qpx/converters/fragpipe/feature_adapter.py:50-51,339-340`; `docs/spec/converter-coverage.md:200`

### EXTERNAL-13: MBR error control differs per tool; PIP-ECHO finds MaxQuant/FlashLFQ lists uncontrolled {#external-13}

*feature* · unspecified · confirmed · checklist F4, D2, F8

- **Now:** FlashLFQ QuantifiedPeaks has 'Peak Detection Type' (MSMS, MBR, …), 'PIP Q-Value', 'PIP PEP', 'Random RT' and 'Decoy Peptide' (decoy transfers are output). PIP-ECHO (PIP Error Control via Hybrid cOmpetition) estimated FDP at a 1% threshold: IonQuant 3.4%, MaxQuant 2.0%, FlashLFQ v1.0 6.1%. For MaxQuant and FlashLFQ v1.0 the software 'reports a list of PIP discoveries without controlling the FDR'.
- **Spec:** identification_type 'transferred' carries no error-control semantics (feature.md:145). The spec does not say whether decoy or random-RT transfer rows may appear in the feature view (is_decoy, feature.md:60).
- **Reader impact:** Equal-weight use of transferred rows (cf. mokume#117) mixes rows of unknown error rate with FDR-controlled MS/MS rows. Decoy transfers could leak in as is_decoy=false.
- **Suggested:** Contract: transferred rows MUST state their method and threshold in provenance, and match_qvalue is null when the tool does not control FDR. Decoy or random transfers are emitted only with is_decoy=true or not at all.
- **Evidence:** `https://github.com/smith-chem-wisc/FlashLFQ/wiki/Output-Formats`; `https://pmc.ncbi.nlm.nih.gov/articles/PMC12488043/ (PIP-ECHO)`; `docs/spec/feature.md:60,145`

### EXTERNAL-14: A feature with conflicting IDs: tools resolve, flag or ignore; QPX has no rule {#external-14}

*feature* · unspecified · partly confirmed (see correction) · checklist F2, X1

- **Now:** FlashLFQ 'MSMSAmbiguousPeakFinding: multiple peptides were mapped to the same peak trace and thus couldn't be quantified'. OpenMS IDMapper: an ID 'is annotated to all' overlapping features, with no conflict resolution. IDConflictResolver keeps 'only one identification … (with the best score)' and moves the rest to unassigned peptides. Dinosaur: 'no control was made for identification consistency'. mzTab-M SMF: SME_ID_REFS plus SME_ID_REF_ambiguity_code (1 ambiguous, 2 several evidence streams of the same molecule, 3 both).
- **Spec:** A feature row carries one peptidoform (feature.md:51-52). Nothing covers several peptidoforms on one chromatographic feature, an ambiguity flag, or a choice between best-only and one row per candidate.
- **Correction:** Dinosaur's 'no control was made for identification consistency' describes how the paper's evaluation matched features to IDs, not a behaviour of the tool's output. Dinosaur itself emits unidentified MS1 features.
- **Reader impact:** Depending on the producer, one MS1 peak can appear as several feature rows, each with the full intensity (double counting), or as one arbitrarily chosen peptide.
- **Suggested:** Adopt an mzTab-M-like rule: an optional feature_group_id shared by rows that come from one MS1 peak, plus id_ambiguity (none/conflicting/multiple_evidence). Consumers must not sum intensities within a group.
- **Evidence:** `https://github.com/smith-chem-wisc/FlashLFQ/wiki/Output-Formats`; `https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_IDMapper.html`; `https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_IDConflictResolver.html`; `https://pmc.ncbi.nlm.nih.gov/articles/PMC4933939 (Dinosaur)`; `https://hupo-psi.github.io/mzTab/2_0-metabolomics-release/mzTab_format_specification_2_0-M_release.html`

### EXTERNAL-15: Producers emit unidentified MS1 features; QPX requires a peptidoform and has no 'unidentified' {#external-15}

*feature* · contradictory · confirmed · checklist F3, D1

- **Now:** MaxQuant allPeptides.txt lists detected isotope patterns and multiplets with 'Identified' +/-, and mqpar has matchUnidentifiedFeatures. OpenMS ProteomicsLFQ targeted_only=false: 'include unidentified features so they can be linked to identified ones (=match between runs)'. Dinosaur and Biosaur2 detect features independently of IDs. PSI-MS has the QC term MS:4000093 for unidentified multi-run-matched MS1 feature areas.
- **Spec:** Unidentified features need a producer-supplied id (feature.md:7-10), but sequence, peptidoform and charge are required (feature.md:51-55). identification_type has no 'unidentified' value (feature.md:145; validation.py enum). The OpenMS consensus converter omits features without hits (psm.md:175).
- **Reader impact:** Producers cannot legally write unidentified features. Some drop them, others might fill placeholder peptidoforms, and a reader cannot rely on 'feature view = all detected peaks' or 'identified only'.
- **Suggested:** Decide: unidentified features go in a separate view (or are out of scope), or peptidoform/sequence become nullable with identification_type='unidentified'. Footer declares feature_scope (identified_only/all_detected).
- **Evidence:** `https://cox-labs.github.io/coxdocs/output_tables.html (allPeptides Identified)`; `https://storage.jpostdb.org/JPST001663/mqpar_Normal_IgG_SPACE.xml (matchUnidentifiedFeatures)`; `https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_ProteomicsLFQ.html`; `https://pmc.ncbi.nlm.nih.gov/articles/PMC4933939`; `docs/spec/feature.md:7-10,51-55,145`; `qpx/core/data/validation.py:71`

### EXTERNAL-16: 'Best PSM' behind feature scan/score is chosen differently by each tool {#external-16}

*feature* · inconsistent across producers · partly confirmed (see correction) · checklist F2, F5, F10

- **Now:** MaxQuant evidence: 'MS/MS scan number: … of the MS/MS with the highest peptide identification score'; 'Score: Andromeda score for the best associated MS/MS'; PEP is the evidence row's own PEP. QPX FragPipe uses the lowest-PEP psm.tsv PSM, then the lowest scan (converter-coverage.md:223-228). QPX PD uses the lowest q-value PSM (converter-coverage.md:324). OpenMS has its own ranking (psm.md:158-161).
- **Spec:** The feature's scan and id_run_file_name are those of 'the best PSM that identified the feature' (feature.md:138-139,301). 'Best' is not defined, and neither is whether PEP/q-value describe that PSM or an aggregate.
- **Correction:** psm.md:158-161 does not define OpenMS's choice of a feature's best PSM. It describes which consensus feature an OpenMS PSM's feature_id points to. The OpenMS best-PSM rule was not established.
- **Reader impact:** feature.scan, PEP and scores may come from different PSMs (highest score vs lowest PEP vs lowest q-value) depending on producer. Feature-to-PSM joins on scan and confidence comparisons are inconsistent.
- **Suggested:** Define best = the tool's own pick when it states one, else lowest PEP, then lowest q-value, then lowest scan (as in #378). Require PEP/q-value/scan to come from that one PSM, or say they are aggregates.
- **Evidence:** `https://cox-labs.github.io/coxdocs/output_tables.html (MS/MS scan number, Score)`; `docs/spec/converter-coverage.md:223-228,324`; `docs/spec/psm.md:158-161`; `docs/spec/feature.md:138-139,301`

### EXTERNAL-17: MS1 intensity is apex, area or volume depending on tool; QPX does not record which {#external-17}

*feature* · unspecified · confirmed · checklist F6, D7

- **Now:** MaxLFQ paper: isotopes 'always summed'; intensity 'full peak volume or the intensity maximum… the latter … default' (coxdocs: 'Summed up XIC'). IonQuant --intensitymode 0 apex / 1 area / 2 auto, with auto added in 1.11.18. FragPipe psm.tsv Intensity = area, ion.tsv = 'maximum intensity from all observed PSMs'. FlashLFQ: by default reports 'peak height', the sum of isotope peaks in the apex scan; integrate optional. Sage LFQ integration Sum|Apex, default Sum.
- **Spec:** intensities holds the 'Raw intensity value' (intensities.md:19,33), with no apex/area/isotope rule. Processing state is recorded only as provenance keys isotope_corrected and normalization (intensities.md:211-213).
- **Reader impact:** Intensities from different producers (or IonQuant versions) differ in kind: apex height vs integrated area vs summed isotopes. Cross-dataset comparison or merge() silently mixes scales.
- **Suggested:** Add a provenance or footer key ms1_intensity_measure (apex|area|volume|unknown) plus isotopes_summed (bool), required for MS1-quantified features. Compare across datasets only when they agree, like isotope_corrected.
- **Evidence:** `https://pmc.ncbi.nlm.nih.gov/articles/PMC4159666 (MaxLFQ)`; `https://cox-labs.github.io/coxdocs/output_tables.html (Intensity)`; `https://github.com/Nesvilab/IonQuant (--intensitymode)`; `https://ionquant.nesvilab.org/CHANGELOG.html`; `https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html`; `https://github.com/smith-chem-wisc/FlashLFQ/wiki/Output-Formats`; `https://sage-docs.vercel.app/docs/configuration/quantification`; `docs/spec/intensities.md:19,211-213`

### EXTERNAL-18: Normalized tool intensities are written as primary 'raw' intensities {#external-18}

*feature* · diverges from spec · partly confirmed (see correction) · checklist F6, G2

- **Now:** FragPipe docs: combined_ion '[experiment] Intensity: normalized ion intensities', and IonQuant --normalization defaults to 1. QPX FragPipe takes '<exp> Intensity' as the primary LFQ intensity (fragpipe/feature_adapter.py:371-373) and writes no normalization key (fragpipe/converter.py:95-119). The MaxQuant LFQ path prefers 'LFQ intensity' (MaxLFQ-normalized) over Intensity (maxquant/feature_adapter.py:521-525), yet provenance says normalization 'none' (maxquant/converter.py:121-122).
- **Spec:** intensities = raw values; normalized, LFQ and iBAQ go to additional_intensities (intensities.md:17-22,219-224; feature.md:303). The normalization key states the processing of primary intensities (intensities.md:213).
- **Correction:** The FragPipe part is confirmed. The MaxQuant part is latent: evidence.txt normally has no 'LFQ intensity' column (coxdocs and both fixtures), so the converter falls back to raw Intensity, and 'normalization none' is correct for real files. The code comment at converter.py:121 contradicts the LFQ-preferring code.
- **Reader impact:** Readers treat FragPipe (and, when present, MaxQuant LFQ) feature intensities as unnormalized raw signal and normalize them again. Provenance may assert 'none' wrongly.
- **Suggested:** Primary intensities MUST be the tool's unnormalized values when available. A tool exporting only normalized values writes them with provenance normalization=<method>. Converters must not claim 'none' without evidence.
- **Evidence:** `https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html`; `https://github.com/Nesvilab/IonQuant (--normalization default 1)`; `qpx/converters/fragpipe/feature_adapter.py:371-373`; `qpx/converters/fragpipe/converter.py:95-119`; `qpx/converters/maxquant/feature_adapter.py:521-525`; `qpx/converters/maxquant/converter.py:121-122`; `docs/spec/intensities.md:17-22,213`

### EXTERNAL-19: MaxQuant reporter column naming is version-dependent; 'Reporter intensity N' may be corrected {#external-19}

*feature* · ambiguous · partly confirmed (see correction) · checklist F6

- **Now:** The fixture evidence.txt has 'Reporter intensity corrected N' and 'Reporter intensity N'. The current coxdocs lists 'Reporter intensity'+i, 'Reporter intensity not corrected'+i and 'Reporter intensity count'+i (evidence and proteinGroups), which suggests newer versions write corrected values under 'Reporter intensity'. Definitions are not given. QPX always treats 'Reporter intensity N' as uncorrected (provenance isotope_corrected=false, maxquant/converter.py:121-122). No code reads 'not corrected'.
- **Spec:** Primary reporter intensities are raw, and isotope_corrected is stated in provenance (intensities.md:211-213). The corrected values go to reporter_intensity_corrected (intensities.md:188).
- **Correction:** The direction is reversed. In MaxQuant 1.5, plain 'Reporter intensity N' is the CORRECTED value and 'not corrected' is raw. From 1.6, plain is raw and 'corrected' is corrected. The coxdocs listing reflects the old naming. QPX mislabels MQ<=1.5 files (corrected values stored as primary with isotope_corrected=false), not newer ones.
- **Reader impact:** For MaxQuant versions using the 'not corrected' naming, QPX labels corrected intensities as uncorrected and loses the raw ones. Isotope-correction state may then be wrong.
- **Suggested:** The converter detects the naming scheme ('not corrected' present means 'Reporter intensity' is corrected) and sets isotope_corrected accordingly. The contract requires isotope_corrected for every isobaric dataset. Verify against a MaxQuant 2.x file (not verified here).
- **Evidence:** `https://cox-labs.github.io/coxdocs/output_tables.html (evidence/proteinGroups Reporter columns)`; `tests/examples/maxquant/maxquant_simple/evidence.txt header`; `qpx/converters/maxquant/converter.py:121-122`; `qpx/converters/maxquant/feature_adapter.py:539-551`; `https://archive.linux.duke.edu/cran/web/packages/PTXQC/vignettes/PTXQC-ListOfMetrics.html (corrected = isotope-corrected)`

### EXTERNAL-20: Reporter quantities may be S/N or log2 ratios; QPX has no unit declaration {#external-20}

*feature* · unspecified · confirmed · checklist F6, D7

- **Now:** In PD, the Reporter Ions Quantifier 'Reporter Abundance Based On' setting can select signal-to-noise, and PSMs below an Average Reporter S/N threshold are excluded. Sage TMT has 'sn: Use Signal/Noise instead of intensity'; its level defaults to 3 (MS3). TMT-Integrator divides by the reference channel, log2-transforms, and converts ratios to abundances using MS1 intensity by default. Its 'Quant level' is 2 or 3.
- **Spec:** intensities are 'Raw intensity value' (intensities.md:33). Provenance has isotope_corrected and normalization only (intensities.md:213). quant_ms_level is optional (feature.md:156), and only OpenMS fills it (converter-coverage.md:56).
- **Reader impact:** TMT intensities from PD (S/N) and MaxQuant (intensity), or log2 TMT-Integrator abundances, look identical. Cross-producer or cross-plex comparisons are wrong.
- **Suggested:** Provenance key reporter_quantity (intensity|signal_to_noise|area|log2_ratio) and quant_ms_level per dataset, mandatory for isobaric data. Log2 or ratio values never go into primary intensities.
- **Evidence:** `https://docs.thermofisher.com/r/Proteome-Discoverer-3.1-User-Guide/en-US1311065483v1`; `https://sage-docs.vercel.app/docs/configuration/quantification`; `https://tmt-integrator.nesvilab.org/`; `docs/spec/intensities.md:33,213`; `docs/spec/feature.md:156`

### EXTERNAL-21: Precursor purity/PIF/interference definitions differ; a PSI-MS term exists but is unmapped {#external-21}

*feature* · ambiguous · confirmed · checklist P7, D7

- **Now:** MaxQuant PIF: 'fraction the target peak makes up of the total intensity' in the window. FragPipe Purity: 'proportion of total ion abundance in the inclusion window from the precursor (including precursor isotopic peaks)'. PD Isolation Interference: % of ion current in the isolation window 'not attributed to the precursor', computed only for high-resolution precursor scans. PSI-MS MS:1003167 'MSn-1 isolation window precursor purity' fits this quantity.
- **Spec:** parent_ion_fraction covers MaxQuant PIF, FragPipe Purity and CDAP PrecursorPurity alike (scores.md:65; feature.md:154). PD's value is converted as 1-interference/100 (proteome_discoverer/reports.py:134). The ontology entry has no accession (ontology.md:318).
- **Reader impact:** Purity thresholds (e.g. ≥0.5) mean different things per producer: target peak vs whole isotope envelope, isolation vs inclusion window.
- **Suggested:** Map parent_ion_fraction to MS:1003167. State whether the precursor's isotopes are included, and keep the source column name as a cv_param.
- **Evidence:** `https://cox-labs.github.io/coxdocs/output_tables.html (PIF)`; `https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html (Purity)`; `https://docs.thermofisher.com/r/Proteome-Discoverer-3.1-User-Guide/en-US1467225483v1`; `https://www.ebi.ac.uk/ols4/api/ontologies/ms/terms (MS:1003167)`; `qpx/converters/proteome_discoverer/reports.py:132-134`; `docs/spec/ontology.md:318`

### EXTERNAL-22: The identification_type enum conflates run of origin with role in a chimeric spectrum {#external-22}

*feature* · ambiguous · confirmed · checklist D1, P2, F4

- **Now:** MaxQuant MULTI-SECPEP is an MS1 cluster 'identified by MS/MS as second peptide', so it is identified by this run's own MS/MS. Sage, CHIMERYS and MSFragger-DDA+ report co-identified peptides without any primary/secondary designation. A peptide identified as second peptide in a donor run and then transferred by MBR fits both 'transferred' and 'secondary'.
- **Spec:** identification_type is one of direct (own MS/MS), transferred (from another run) or secondary (co-fragmenting peptide of a chimeric spectrum), or null (feature.md:145; validation.py:71).
- **Reader impact:** Filtering identification_type='direct' to get 'identified in this run' drops MaxQuant second-peptide rows. Chimeric IDs from other tools appear as 'direct', so 'secondary' is MaxQuant-only.
- **Suggested:** Split the field into identification_origin (direct/transferred/unidentified) and chimeric (bool) or spectrum_role (primary/co_fragmenting). Map MULTI-SECPEP to direct+co_fragmenting.
- **Evidence:** `https://cox-labs.github.io/coxdocs/output_tables.html`; `https://sage-docs.vercel.app/docs/configuration/additional`; `https://pmc.ncbi.nlm.nih.gov/articles/PMC12074992/`; `https://pmc.ncbi.nlm.nih.gov/articles/PMC11978857/`; `docs/spec/feature.md:145`; `qpx/core/data/validation.py:71`

### EXTERNAL-23: MBR and transfer have no PSI vocabulary; QPX's terms need their own definitions {#external-23}

*feature* · unspecified · confirmed · checklist F4, D1

- **Now:** PSI-MS (OLS search) has no term for match-between-runs or identification transfer; the only hit is QC term MS:4000093 'unidentified … multiple-run-matched MS1 feature areas'. mzTab 1.0 has no MBR flag. mzTab-M models features (SMF) pointing to evidence (SME) with evidence_input_id and a rank. OpenMS FeatureFinderIdentification explicitly separates 'internal IDs (derived from the LC-MS/MS run being quantified) and external IDs (inferred from related runs)'.
- **Spec:** identification_type and match_* are QPX-specific columns (feature.md:141-148). The OpenMS converter infers 'transferred' from the absence of IDs in the run (converter-coverage.md:55) instead of reading FFId's internal/external distinction.
- **Reader impact:** No external standard anchors the meaning, so each converter infers 'transferred' by its own heuristic (absence of IDs vs the tool's flag). Results differ for the same data.
- **Suggested:** Define identification origin terms in QPX's ontology (and propose PSI-MS terms). Prefer the tool's explicit flag (FFId internal/external, MaxQuant Type, IonQuant Match Type) over inference, and record which was used.
- **Evidence:** `https://www.ebi.ac.uk/ols4/api/search?q=match%20between%20runs&ontology=ms`; `https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderIdentification.html`; `https://hupo-psi.github.io/mzTab/2_0-metabolomics-release/mzTab_format_specification_2_0-M_release.html`; `docs/spec/converter-coverage.md:55`; `docs/spec/feature.md:141-148`

### EXTERNAL-24: Same precursor in one run as both MS/MS and MBR feature: QC flags it, QPX keeps both {#external-24}

*feature* · ambiguous · partly confirmed (see correction) · checklist F1, F4

- **Now:** PTXQC MBR ID-Transfer: 'If MaxQuant only transfers peptide ID's which are not present in the target file, then each Raw file should not have any duplicates of identical peptides (incl. charge)'. It classes them as single, in-group or out-group, where out-group is 'a false ID'. FeatureFinderIdentification keeps 'at most one feature per peptide and charge state' per RT region.
- **Spec:** The default composite includes rt (feature.md:22-23), so several rows per (peptidoform, charge, run) are allowed. anndata.md:153 sums or averages 'a MaxQuant MS/MS and a match-between-runs feature' for one precursor. There is no rule for which row is authoritative.
- **Correction:** The FFId docs do not say 'per RT region'. They say at most one feature per peptide and charge state remains after filtering.
- **Reader impact:** Duplicate precursor rows (direct + transferred) are summed by default in MuData, which inflates the intensity and can propagate a false transfer.
- **Suggested:** Contract: several rows per precursor and run are allowed only as distinct MS1 peaks. When a direct row exists, a transferred duplicate is flagged and excluded from default aggregation (prefer direct, not sum).
- **Evidence:** `https://archive.linux.duke.edu/cran/web/packages/PTXQC/vignettes/PTXQC-ListOfMetrics.html`; `https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderIdentification.html`; `docs/spec/feature.md:22-23`; `docs/spec/anndata.md:153`

### EXTERNAL-30: Zero vs null is left open by PSI; QPX's stricter rule is worth keeping {#external-30}

*feature* · specified, worth keeping · confirmed · checklist D6, F6, G1

- **Now:** mzTab 1.0 §5.9: no empty cells, 'null' when unavailable, INF/NaN for ratio errors; whether zero spectral counts are 'zero or missing (null)' is 'open for debate'. mzTab-M SMF abundance_assay: 'Null or zero may be used'. Tools write 0 as a not-quantified sentinel (MaxQuant, FragPipe, DIA-NN, per converter-coverage.md:60,194).
- **Spec:** A null intensity means identified, not quantified. 0 is a measured zero. No entry means not identified or no sample. A tool's 0 sentinels become null, and NaN, infinite and negative values are invalid (intensities.md:121-133; producer-guide.md:41-48).
- **Reader impact:** Without this rule, a mzTab-derived or tool-native 0 would be read as a measured absence and bias imputation and fold changes.
- **Suggested:** Keep the rule in the contract. Require each converter to document its sentinel mapping. Add a validator warning when a source known to use 0 sentinels yields 0.0 intensities.
- **Evidence:** https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/specification_document-releases/1_0-Proteomics-Release/mzTab_format_specification.pdf §5.9; `https://hupo-psi.github.io/mzTab/2_0-metabolomics-release/mzTab_format_specification_2_0-M_release.html (abundance_assay)`; `docs/spec/intensities.md:121-133`; `docs/spec/producer-guide.md:41-48`

### EXTERNAL-M4: Spec says features carry no PEP, yet feature.posterior_error_probability exists with producer-specific meaning {#external-m4}

*feature* · contradictory · found by the verifier · checklist F5, P5

- **Now:** MaxQuant feature PEP is the evidence row's PEP (null for MBR). coxdocs says PEP 'essentially operates as a p-value'. FragPipe uses 1 - PeptideProphet probability of the best psm.tsv PSM, and Spectronaut uses EG.PEP (converter-coverage.md:180).
- **Spec:** psm.md:286: 'The feature and peptide views do not carry PEP as a top-level field'. But feature.md:56 and feature.yaml:37-41 define feature.posterior_error_probability 'for the peptide match' without naming which identification.
- **Reader impact:** Readers cannot tell whether feature PEP is the best PSM's, an evidence aggregate or precursor-level, and one spec page says the field does not exist.
- **Suggested:** Delete the psm.md:286 sentence (or the column). Define feature PEP as that of the PSM named by scan/id_run_file_name, or declare its level; null for transferred rows.
- **Evidence:** `docs/spec/psm.md:286`; `docs/spec/feature.md:56`; `qpx/core/data/schemas/feature.yaml:37-41`; `docs/spec/converter-coverage.md:180`; `https://cox-labs.github.io/coxdocs/output_tables.html (PEP)`

## First-principles questions

### CRITIC-21: What is one feature row: a chromatographic peak, a precursor per run, or a per-unit quantity? {#critic-21}

*feature* · contradictory · confirmed · checklist F1, F7

- **Now:** data-model says one quantified chromatographic peak. producer-guide says one row per precursor, charge and run. The default composite includes rt, which allows several peaks per precursor and run. A FragPipe row is an experiment (quantification unit), and quantms TMT includes scan in the composite (one row per PSM).
- **Spec:** data-model.md:64-67; producer-guide.md:12; feature.md:22-27; identity.md:27-35
- **Reader impact:** Summing feature intensities per precursor and run double-counts for some producers and not for others.
- **Suggested:** Fix the DDA grain at one row per (peptidoform, charge, run) quantity, with peak detail in xic. Declare a footer feature_grain.
- **Evidence:** `docs/spec/data-model.md:64-67`; `docs/spec/producer-guide.md:12`; `docs/spec/feature.md:22-27`; `docs/spec/identity.md:27-35`

### CRITIC-22: Can (peptidoform, charge, run_file_name) appear in several feature rows? {#critic-22}

*feature* · ambiguous · confirmed · checklist F1

- **Now:** Yes, whenever the composite includes rt, scan or rt_start (schema default, MaxQuant, OpenMS, quantms). No footer flag says the natural key is unique. The API views sum features within a unit.
- **Spec:** identity.md:27-35; views.md:212; feature.md:22-27
- **Reader impact:** PSM-to-feature joins and matrix pivots produce duplicates or silent sums.
- **Suggested:** Require (peptidoform, charge, run_file_name) to be unique, or declare natural_key_unique in the footer.
- **Evidence:** `docs/spec/identity.md:27-35`; `docs/spec/views.md:212`

### CRITIC-23: What happens when one feature has conflicting identifications (different peptidoforms)? {#critic-23}

*feature* · unspecified · confirmed · checklist F2

- **Now:** No rule. Because peptidoform is in the composite, a producer may emit one row per peptidoform sharing the same intensity, keep the best, or drop the feature. OpenMS nulls protein fields on conflicting assignments, but that rule is not in the spec.
- **Spec:** none (feature.md silent); converter-coverage.md:57,62-64
- **Reader impact:** The same MS1 signal can be counted for two peptides, or a quantified peptide can vanish silently.
- **Suggested:** One row per feature with the best peptidoform. Flag the conflict (identification_conflict bool); never duplicate an intensity.
- **Evidence:** `docs/spec/converter-coverage.md:57`; `docs/spec/converter-coverage.md:62-64`

### CRITIC-24: Are unidentified features allowed, and what fills the required peptidoform/sequence/charge? {#critic-24}

*feature* · contradictory · confirmed · checklist F3, D1

- **Now:** feature.md says unidentified features need a producer id. But sequence, peptidoform and charge are required non-null, and identification_type has no 'unidentified' value. psm.md says the OpenMS converter omits features without peptide hits. #378 is undecided.
- **Spec:** feature.md:8-10,51-55,145; psm.md:175; feature.yaml:16-23
- **Reader impact:** A producer may write empty peptidoforms, and readers cannot reliably filter to identified rows.
- **Suggested:** Disallow unidentified features in the DDA feature view, move them to a separate view, and drop feature.md:8-10.
- **Evidence:** `docs/spec/feature.md:8-10`; `docs/spec/psm.md:175`; `qpx/core/data/schemas/feature.yaml:16-23`

### CRITIC-25: How does a reader tell 'identified by MS2 in this run' from a match-between-runs transfer? {#critic-25}

*feature* · ambiguous · partly confirmed (see correction) · checklist F4, D1, D2

- **Now:** identification_type is direct/transferred/secondary, null when the tool does not say, optional, and a free string (no enum in the yaml); DIA-NN leaves it null. MaxQuant MULTI-MATCH-MSMS maps to transferred, and Spectronaut EG.IsImputed maps to transferred.
- **Spec:** feature.md:145-147; feature.yaml:210-214; converter-coverage.md:55,199; producer-guide.md:26
- **Correction:** The validator does enforce the closed set {direct, transferred, secondary}, null allowed (qpx/core/data/validation.py:71,298-304). Only the yaml lacks the enum. Optional/null and the cross-producer mapping issues hold.
- **Reader impact:** MBR rows cannot be excluded when the column is absent or null, and the same value covers different mechanisms across producers.
- **Suggested:** Required column with a validated enum {direct, transferred, secondary, unknown}. transferred = no accepted MS2 identification of this peptidoform in this run.
- **Evidence:** `docs/spec/feature.md:145`; `qpx/core/data/schemas/feature.yaml:210-214`; `docs/spec/converter-coverage.md:55`; `docs/spec/converter-coverage.md:199`

### CRITIC-26: Is a transferred (MBR) quantity the same thing as an imputed one? {#critic-26}

*feature* · contradictory · partly confirmed (see correction) · checklist F4, D1, D7

- **Now:** The spec separates them: transferred = matched from another run; imputed_labels = value imputed instead of measured. The Spectronaut converter derives both identification_type=transferred and imputed_labels (and pg.is_imputed) from the single EG.IsImputed flag.
- **Spec:** feature.md:145-148; intensities.md:131-133; converter-coverage.md:199,201
- **Correction:** The code does derive both flags from EG.IsImputed (spectronaut/feature_adapter.py:338-343), but no spec text forbids it. This is a semantic conflation across producers (Spectronaut imputed values vs MaxQuant MBR), better classed as inconsistent_across_producers than contradictory.
- **Reader impact:** Imputed background values are counted as MBR identifications, or MBR rows are dropped as imputed. MBR-aware normalization (mokume#117/118) misfires.
- **Suggested:** Keep the two flags orthogonal, and state that identification_type must never be derived from imputation flags.
- **Evidence:** `docs/spec/feature.md:145-148`; `docs/spec/converter-coverage.md:199`; `docs/spec/converter-coverage.md:201`

### CRITIC-27: For a transferred feature, is id_run_file_name the donor run or null? {#critic-27}

*feature* · contradictory · partly confirmed (see correction) · checklist F4, F10

- **Now:** The field table says 'run file containing the best PSM that identified the feature (may differ from run_file_name)', which implies the donor. The Identification Origin section says it is null for a transferred row. Spectronaut always fills it; MaxQuant and FragPipe leave it null for MBR rows.
- **Spec:** feature.md:138 vs feature.md:145; converter-coverage.md:57,203; scan.md:124
- **Correction:** feature.md:138 vs :145 is a real contradiction. But Spectronaut nulls id_run_file_name when EG.IsImputed is true (spectronaut/feature_adapter.py:341), so 'always fills' is wrong (converter-coverage.md:203 is imprecise). No converter writes a donor run. MaxQuant MULTI-MATCH-MSMS keeps its own run while typed transferred (CRITIC-M4).
- **Reader impact:** id_run_file_name != run_file_name cannot serve as an MBR detector, and donor provenance is lost.
- **Suggested:** id_run_file_name = the run of the identifying PSM. For transferred rows: the donor run if known, else null; or add a separate match_donor_run column.
- **Evidence:** `docs/spec/feature.md:138`; `docs/spec/feature.md:145`; `docs/spec/converter-coverage.md:203`

### CRITIC-28: Which run's scan numbering does feature.scan refer to? {#critic-28}

*feature* · contradictory · partly confirmed (see correction) · checklist F10

- **Now:** feature.scan is the scan of the best PSM, which may lie in id_run_file_name != run_file_name. Yet the mz-to-feature join is run_file_name+scan. producer-guide says scan is empty for transferred rows.
- **Spec:** feature.md:138-139; data-model.md:166; producer-guide.md:23,26
- **Correction:** The issue is latent only. Every converter writes id_run_file_name = run_file_name or null (maxquant feature_adapter.py:657, fragpipe :340, spectronaut :341, cdap :223, PD :138, quantms :525, openms :545), so no cross-run scan occurs today. The spec wording (feature.md:138 'may differ' with data-model.md:166) still allows one.
- **Reader impact:** A feature-to-spectrum join can fetch an unrelated spectrum from the wrong run.
- **Suggested:** feature.scan must belong to run_file_name, and is empty when the identifying spectrum is in another run.
- **Evidence:** `docs/spec/feature.md:138-139`; `docs/spec/data-model.md:166`; `docs/spec/producer-guide.md:26`

### CRITIC-29: Is feature.scan one best spectrum or all supporting spectra? {#critic-29}

*feature* · contradictory · confirmed · checklist F10

- **Now:** feature.md says the scan of the best PSM. Under the OpenMS composite, scan holds every spectrum reference in the run, as combined components rather than one native ID. The scan_format cardinality check is not applied to the feature view.
- **Spec:** feature.md:26-27,139; scan.md:90-92,124
- **Reader impact:** Decoding [a,b] as a Bruker frame/scan gives nonsense, and 'best PSM' lookups fail for OpenMS output.
- **Suggested:** feature.scan = the single identifying spectrum, with scan_format applying. All supporting spectra come via psm.feature_id.
- **Evidence:** `docs/spec/feature.md:26-27`; `docs/spec/feature.md:139`; `docs/spec/scan.md:90-92`

### CRITIC-30: Which PSM is 'best' for feature scan, id_run_file_name, PEP and scores? {#critic-30}

*feature* · unspecified · confirmed · checklist F2, F5

- **Now:** 'Best PSM' is undefined. FragPipe takes lowest PEP then lowest scan; PD takes the lowest q-value PSM; OpenMS uses its own ranking. The #378 proposal (tool's pick, else lowest PEP, ties by lowest scan) is still open.
- **Spec:** feature.md:138-139,301; converter-coverage.md:223-225,324
- **Reader impact:** Feature-level identification fields differ by producer even for identical PSM sets.
- **Suggested:** Adopt #378: the tool's designated PSM, else lowest PEP, then lowest q-value, then lowest scan, restricted to run_file_name.
- **Evidence:** `docs/spec/feature.md:138-139`; `docs/spec/converter-coverage.md:223-225`; `docs/spec/converter-coverage.md:324`

### CRITIC-31: Which identification do feature PEP, q-values and additional_scores describe? {#critic-31}

*feature* · ambiguous · confirmed · checklist F5

- **Now:** PEP is 'for the peptide match'. peptide_qvalue is the best PSM's value under the Percolator/OpenMS rule. additional_scores have no stated scope. For transferred rows MaxQuant writes null; other producers are unspecified (donor's values? the match's?).
- **Spec:** feature.md:56-59,125-132,146; converter-coverage.md:180,184
- **Reader impact:** A q-value filter may admit MBR rows carrying donor confidence, or exclude them, depending on producer.
- **Suggested:** Identification scores describe the identifying PSM in this run and are null for transferred rows. Transfer confidence goes only in match_*.
- **Evidence:** `docs/spec/feature.md:56-59`; `docs/spec/feature.md:125-132`; `docs/spec/converter-coverage.md:180`

### CRITIC-32: Can a PSM-level q-value end up in feature.precursor_qvalue? {#critic-32}

*feature* · contradictory · confirmed · checklist D4, F5

- **Now:** Yes. precursor_qvalue is filled from precursor_qvalue, qvalue, q_value or psm_level_qvalue (qpx/writers/derived.py:32). Running the CDAP fixture: all 94 features have precursor_qvalue equal to the best PSM's psm_level_qvalue.
- **Spec:** scores.md:44 vs feature.md:58 ('Run-level precursor q-value')
- **Reader impact:** A filter on the 'precursor' q-value actually uses one PSM's q-value, so levels are conflated across producers.
- **Suggested:** Fill precursor_qvalue only from explicitly precursor-level names; never promote PSM q-values.
- **Evidence:** `docs/spec/scores.md:44`; `docs/spec/feature.md:58`; `qpx/writers/derived.py:31-32`; `audit-run:critic/cdap/cdap.feature.parquet (qpxc convert cdap on tests/examples/cdap)`

### CRITIC-33: converter-coverage says CDAP does not fill feature.peptide_qvalue; does it? {#critic-33}

*feature* · contradictory · confirmed · checklist D4

- **Now:** converter-coverage says 'Not implemented: best PepQvalue stored as peptide_level_qvalue'. The writer copies peptide_level_qvalue into peptide_qvalue, and in the CDAP fixture output all 94 features have peptide_qvalue filled.
- **Spec:** converter-coverage.md:181 vs scores.md:43
- **Reader impact:** Readers who trust the coverage table assume the column is null. The docs cannot be relied on for which level a column holds.
- **Suggested:** Regenerate the coverage tables from measured 1.3 output after the typed-column fill.
- **Evidence:** `docs/spec/converter-coverage.md:181`; `docs/spec/scores.md:43`; `qpx/converters/cdap/feature_adapter.py:294`; `audit-run:critic/cdap/cdap.feature.parquet`

### CRITIC-34: Is peptide_qvalue per sequence or per peptidoform, and per run or per experiment? {#critic-34}

*feature* · ambiguous · confirmed · checklist D4, F5

- **Now:** Only 'Peptide-level q-value' is given. Neither the scope nor what 'peptide' means (stripped sequence, peptidoform or precursor) is stated. CDAP takes the best PSM's PepQvalue per run.
- **Spec:** feature.md:57,125-132; scores.md:43
- **Reader impact:** Peptide-FDR filters behave differently across producers.
- **Suggested:** peptide = peptidoform (modifications included), experiment-wide. Run-level variants take a different name.
- **Evidence:** `docs/spec/feature.md:57`; `docs/spec/feature.md:125-132`

### CRITIC-36: Is feature.rt the peak apex, the best PSM's MS2 RT, or something else? {#critic-36}

*feature* · ambiguous · confirmed · checklist F10

- **Now:** Only 'Precursor retention time' is given. Producers write different things: OpenMS the consensus apex, PD Apex RT, CDAP RTAtPrecursorHalfElution, MaxQuant evidence 'Retention time'; FragPipe does not fill it.
- **Spec:** feature.md:65; data-model.md:64; converter-coverage.md:54,189,204,324
- **Reader impact:** RT matching, alignment QC and the rt used in identity composites mean different things per producer.
- **Suggested:** feature.rt = apex of the quantified MS1 peak on the run's own clock. The MS2 RT stays in psm.
- **Evidence:** `docs/spec/feature.md:65`; `docs/spec/converter-coverage.md:54`; `docs/spec/converter-coverage.md:204`; `docs/spec/converter-coverage.md:324`

### CRITIC-37: In which unit and time frame are rt_start and rt_stop? {#critic-37}

*feature* · unspecified · confirmed · checklist F10

- **Now:** No unit or frame is stated. MaxQuant writes the calibrated start/finish while rt is the uncalibrated time (qpx/converters/maxquant/feature_adapter.py:574-583), so the window and the apex are on different clocks. DIA-NN converts minutes to seconds.
- **Spec:** feature.md:66,69-70; feature.yaml:255-260; converter-coverage.md:54,204
- **Reader impact:** Peak widths, and checks that rt lies within [rt_start, rt_stop], come out wrong.
- **Suggested:** rt_start/rt_stop in seconds on rt's clock. Aligned bounds, if needed, go in separate columns.
- **Evidence:** `docs/spec/feature.md:69-70`; `qpx/converters/maxquant/feature_adapter.py:574-583`; `docs/spec/converter-coverage.md:204`

### CRITIC-38: What physical quantity is a feature intensity (apex, area, isotope sum, reporter sum)? {#critic-38}

*feature* · unspecified · confirmed · checklist F6

- **Now:** Only 'Raw intensity value' is given, and the feature view has no quantification_method (pg has an optional cv). Producers use MaxQuant Intensity, DIA-NN Precursor.Quantity, PD's largest Precursor Abundance, FragPipe Intensity and OpenMS sub-feature intensity.
- **Spec:** intensities.md:17-22,33; feature.md:84; pg.md:52-60
- **Reader impact:** Peptide intensities are not comparable across producers.
- **Suggested:** Mandatory footer or provenance key feature_intensity_type (apex|area|volume|reporter_sum|...).
- **Evidence:** `docs/spec/intensities.md:17-22`; `docs/spec/pg.md:52-60`; `docs/spec/converter-coverage.md:324`

### CRITIC-39: For isobaric data, is a feature one PSM or a sum over PSMs, and whose PIF applies? {#critic-39}

*feature* · unspecified · confirmed · checklist F1, F6, D7

- **Now:** The quantms TMT composite includes scan (one row per PSM). PD sums reporters over rank-1 PSMs, and MaxQuant writes per evidence. parent_ion_fraction and quant_scan are single values per row, with no rule for aggregated rows.
- **Spec:** identity.md:33; converter-coverage.md:324; feature.md:154-156; intensities.md:37-49
- **Reader impact:** Reporter intensities scale with the number of PSMs summed, and PIF filtering is ambiguous.
- **Suggested:** Isobaric feature = precursor/run aggregate, with the aggregation rule and n_psms recorded. Per-PSM reporter values go on psm.
- **Evidence:** `docs/spec/identity.md:33`; `docs/spec/converter-coverage.md:324`; `docs/spec/feature.md:154-156`

### CRITIC-44: For a shared peptide, what do pg_accessions, anchor_protein and unique hold? {#critic-44}

*feature* · contradictory · confirmed · checklist F9

- **Now:** feature.md: the group the producer assigned (razor); anchor null only when no mapping. feature.yaml: anchor and pg_accessions null when the evidence spans more than one unassigned group, unique null then. MaxQuant fallback: unique true if one protein, or if 'in a single group of any size'.
- **Spec:** feature.md:99-101; feature.yaml:151-161; tool-mappings.md:62; versioning.md:68-69
- **Reader impact:** Uniqueness filters for protein quantification differ by producer, and by whether proteinGroups.txt was given.
- **Suggested:** One definition (versioning.md:68-69) in both md and yaml, null when unresolved; drop the fallback for unknown group size.
- **Evidence:** `docs/spec/feature.md:99-101`; `qpx/core/data/schemas/feature.yaml:151-161`; `docs/spec/tool-mappings.md:62`

### CRITIC-45: Is feature.pg_accessions a list of strings or of structs, and may it carry positions? {#critic-45}

*feature* · contradictory · confirmed · checklist F9

- **Now:** The md type is array[struct], but the examples show string arrays. The yaml says the position fields are reserved and positions belong in pg_positions. Yet the producer-guide example fills start/end/pre/post in it, and so does CDAP (pre/post).
- **Spec:** feature.md:99,198,261; feature.yaml:151-153; producer-guide.md:120; converter-coverage.md:197
- **Reader impact:** Positions read from pg_accessions come out null or not depending on producer.
- **Suggested:** Only the accession is meaningful in pg_accessions; positions live only in pg_positions. Fix the examples.
- **Evidence:** `docs/spec/feature.md:198`; `qpx/core/data/schemas/feature.yaml:153`; `docs/spec/producer-guide.md:120`; `docs/spec/converter-coverage.md:197`

### CRITIC-46: Is feature.run_file_name always the raw file the signal was measured in? {#critic-46}

*feature* · contradictory · confirmed · checklist F1, X1

- **Now:** feature.md says it is 'the run file name that contains the feature'. For FragPipe, run_file_name is the first raw run of the experiment when the experiment has several runs. quantification_unit_id holds the experiment, and its relation to pg.grouped_runs is undefined.
- **Spec:** feature.md:89-90; data-model.md:64; converter-coverage.md:213-216
- **Reader impact:** Per-run analyses attribute fraction-aggregated quantities to a single fraction.
- **Suggested:** run_file_name must be the acquiring run. Unit-level quantities need an explicit list of member runs.
- **Evidence:** `docs/spec/feature.md:89-90`; `docs/spec/converter-coverage.md:213-216`

### CRITIC-47: Which features were actually used to compute the protein quantity? {#critic-47}

*feature* · unspecified · confirmed · checklist G2, F9

- **Now:** No field exists ('quantified' vs 'used for quantification', #282). The unique/razor flags only approximate it, and a feature-pg link does not imply use.
- **Spec:** none; pg.md:70-74
- **Reader impact:** pg intensities cannot be reproduced or audited from the features.
- **Suggested:** Optional bool used_for_pg_quant on feature, null = unknown.
- **Evidence:** `docs/spec/pg.md:70-74`

### CRITIC-56: How does a feature link to its pg row, and when can the link fail? {#critic-56}

*feature* · specified, worth keeping · confirmed · checklist F9

- **Now:** A computed softlink: canonical pg_accessions are equal, the run is in grouped_runs, and the label matches; never anchor_protein. No match means no link, not an error. FragPipe without --pg-file writes only the lead protein, so no row matches.
- **Spec:** pg.md:70-74; data-model.md:159; converter-coverage.md:217-218
- **Reader impact:** Keep, but readers must expect silent non-links when a producer writes incomplete membership.
- **Suggested:** Keep, and require full membership or null, never a truncated list.
- **Evidence:** `docs/spec/pg.md:70-74`; `docs/spec/data-model.md:159`; `docs/spec/converter-coverage.md:217-218`

### CRITIC-M4: Does a MaxQuant MULTI-MATCH-MSMS row keep its own run while typed 'transferred'? {#critic-m4}

*feature* · diverges from spec · found by the verifier · checklist F4, D1

- **Now:** evidence_fields.py:16 maps MULTI-MATCH-MSMS to transferred, but feature_adapter.py:657 nulls id_run_file_name only for MULTI-MATCH. These rows therefore carry id_run_file_name = run_file_name, plus a scan and PEP when present. This comes from reading the code; the fixtures have no such rows.
- **Spec:** feature.md:145 says id_run_file_name is null for a transferred row. producer-guide.md:26 says scan is empty for transferred rows. converter-coverage.md:203 says it is null for MBR rows.
- **Reader impact:** MBR detection by id_run_file_name/scan disagrees with identification_type, and identification scores may leak onto transferred rows.
- **Suggested:** Null id_run_file_name, scan and identification scores on every transferred row, or classify MULTI-MATCH-MSMS explicitly with match_* filled.
- **Evidence:** `qpx/converters/maxquant/evidence_fields.py:8-17`; `qpx/converters/maxquant/feature_adapter.py:655-657`; `docs/spec/feature.md:145`; `docs/spec/producer-guide.md:26`
