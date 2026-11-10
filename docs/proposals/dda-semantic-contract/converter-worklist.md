# Converter and tooling worklist

[Overview](index.md) · [Annex overview](annex-overview.md) · [OpenMS work list](openms-worklist.md)

!!! note "Annex"
    Part of the annex of the [DDA semantic contract proposal](index.md). The OpenMS items below target the Python consensusXML converter. The work on the OpenMS native export, which has priority, is in the [OpenMS work list](openms-worklist.md).

Every verified finding that is a defect in a converter, writer, validator, document or test, rather than an open semantic question. Each line names the rule of the proposal that settles it. Lines come from the completeness pass over all findings; follow the finding link for code references and evidence. Consumer-side defects (qpx read API, MuData export, mokume) are listed in [T14](t14-reader-contract.md#reader-rules).

## MaxQuant

25 items.

- [ ] [MAXQUANT-01](findings-feature.md#maxquant-01): MULTI-MATCH-MSMS rows typed transferred keep run, scan, PEP ([MBR-3](t08-match-between-runs.md)); also [SPEC-M1](findings-feature.md#spec-m1), [MOKUME-M5](findings-feature.md#mokume-m5), [CRITIC-M4](findings-feature.md#critic-m4), [EXTERNAL-10](findings-feature.md#external-10)
- [ ] [MAXQUANT-02](findings-feature.md#maxquant-02): Types outside the six-key map (TIMS-*, ISO-SECPEP) silently null; token map + cv_param maxquant_type ([MBR-2](t08-match-between-runs.md))
- [ ] [MAXQUANT-03](findings-psm.md#maxquant-03): SECPEP twins with equal peptidoform/charge duplicate psm_id; merge ([PSM-2](t01-psm-rows.md), [CHI-6](t02-chimeric-spectra.md))
- [ ] [MAXQUANT-05](findings-psm.md#maxquant-05): FDR-cut msms written with null pass_threshold and no fdr_filter_psm ([FLT-9](t03-filtering-state.md))
- [ ] [MAXQUANT-06](findings-cross.md#maxquant-06): tool version, PSM/protein FDR and MBR setting not recorded ([META-11](t13-metadata-conformance.md), [FLT-1](t03-filtering-state.md), [MBR-8](t08-match-between-runs.md))
- [ ] [MAXQUANT-08](findings-feature.md#maxquant-08): one isotope pattern credited in full to 2-3 peptidoforms ([FID-10](t06-identifications-per-feature.md))
- [ ] [MAXQUANT-09](findings-feature.md#maxquant-09): MS/MS-only rows carry a synthetic 60 s window ([FRW-9](t05-feature-rows.md))
- [ ] [MAXQUANT-10](findings-feature.md#maxquant-10): rt on the raw clock but rt_start/rt_stop calibrated ([FRW-7](t05-feature-rows.md))
- [ ] [MAXQUANT-12](findings-feature.md#maxquant-12): MBR rows written with all-null reporter entries; drop and count ([MBR-3](t08-match-between-runs.md))
- [ ] [MAXQUANT-14](findings-cross.md#maxquant-14): reporter correction state differs between feature/pg and site and by version ([INT-1](t09-intensities-missing-values.md), [INT-2](t09-intensities-missing-values.md)); also [EXTERNAL-19](findings-feature.md#external-19)
- [ ] [MAXQUANT-15](findings-feature.md#maxquant-15): TMT evidence without SDRF written as label LFQ with MS1 intensity; must fail ([INT-5](t09-intensities-missing-values.md))
- [ ] [MAXQUANT-16](findings-feature.md#maxquant-16): latent per-row 'LFQ intensity' preference for the primary ([INT-1](t09-intensities-missing-values.md))
- [ ] [MAXQUANT-17](findings-pg.md#maxquant-17): MaxLFQ named 'lfq' instead of 'maxlfq' ([INT-12](t09-intensities-missing-values.md))
- [ ] [MAXQUANT-18](findings-pg.md#maxquant-18): experiment-level Q-value copied into unit-level pg_qvalue ([PGV-8](t12-protein-groups.md))
- [ ] [MAXQUANT-19](findings-pg.md#maxquant-19): proteinGroups PEP never read; write pg_pep ([CNF-5](t04-confidence-fields.md))
- [ ] [MAXQUANT-20](findings-pg.md#maxquant-20): 'By matching' pg rows marked only by cv_params; write identified_in_unit ([MBR-10](t08-match-between-runs.md), [PGV-1](t12-protein-groups.md))
- [ ] [MAXQUANT-22](findings-cross.md#maxquant-22): evidence 'Potential contaminant' never read; write psm/feature contaminant ([PRT-6](t10-protein-mapping.md))
- [ ] [MAXQUANT-23](findings-psm.md#maxquant-23): psm.observed_mz rebuilt from QPX calculated_mz, hiding peptidoform errors ([SPC-1](t11-spectrum-precursor-fields.md))
- [ ] [MAXQUANT-24](findings-cross.md#maxquant-24): missing fixed mods only warn; protein-terminal fixed mods skipped (SPC-14d)
- [ ] [MAXQUANT-25](findings-psm.md#maxquant-25): Evidence ID links not checked for run/charge/peptidoform ([FID-7](t06-identifications-per-feature.md))
- [ ] [MAXQUANT-27](findings-cross.md#maxquant-27): no scan_format; MQ 1.5 'MS/MS Scan Number' unmatched, empty scans ([SPC-7](t11-spectrum-precursor-fields.md), [FID-5](t06-identifications-per-feature.md))
- [ ] [MAXQUANT-28](findings-feature.md#maxquant-28): without proteinGroups.txt pg_accessions = leaders of several groups ([PRT-7](t10-protein-mapping.md))
- [ ] [MAXQUANT-29](findings-feature.md#maxquant-29): dedup ranks PEP before type and can relabel secondary/direct ([MBR-7](t08-match-between-runs.md), [CHI-9](t02-chimeric-spectra.md))
- [ ] [MAXQUANT-M1](findings-cross.md#maxquant-m1): NULL Modified sequence written as peptidoform 'None' ([UNI-1](t07-unidentified-features.md), [SPC-15](t11-spectrum-precursor-fields.md))
- [ ] [MAXQUANT-M3](findings-cross.md#maxquant-m3): parameters.txt/summary.txt sidecars change peptidoforms but are not in provenance ([META-11](t13-metadata-conformance.md))

## FragPipe

23 items.

- [ ] [EXTERNAL-12](findings-feature.md#external-12): IonQuant transfer FDR/probability not written to match_* / mbr_fdr ([MBR-4](t08-match-between-runs.md), [MBR-9](t08-match-between-runs.md))
- [ ] [FRAGPIPE-01](findings-psm.md#fragpipe-01): pre-filtered psm.tsv with null pass_threshold and no filter keys ([FLT-9](t03-filtering-state.md))
- [ ] [FRAGPIPE-02](findings-psm.md#fragpipe-02): PEP = 1-Probability labelled percolator:PEP whatever the validator ([CNF-7](t04-confidence-fields.md))
- [ ] [FRAGPIPE-04](findings-psm.md#fragpipe-04): dotted run names truncate run_file_name, empty scan, duplicate psm_id ([SPC-13](t11-spectrum-precursor-fields.md))
- [ ] [FRAGPIPE-06](findings-psm.md#fragpipe-06): calibrated observed m/z written without a calibration declaration ([SPC-2](t11-spectrum-precursor-fields.md))
- [ ] [FRAGPIPE-07](findings-cross.md#fragpipe-07): contam_ marker stripped in psm, feature, pepmap, pg ([PRT-3](t10-protein-mapping.md)); also [MOKUME-M2](findings-feature.md#mokume-m2)
- [ ] [FRAGPIPE-08](findings-psm.md#fragpipe-08): is_decoy from the razor Protein only ([PRT-4](t10-protein-mapping.md))
- [ ] [FRAGPIPE-09](findings-psm.md#fragpipe-09): psm.feature_id never filled although features pick a PSM ([FID-8](t06-identifications-per-feature.md))
- [ ] [FRAGPIPE-10](findings-feature.md#fragpipe-10): experiment quantity labelled with its first raw file, no grouped_runs ([FRW-6](t05-feature-rows.md))
- [ ] [FRAGPIPE-12](findings-feature.md#fragpipe-12): PSM lookup omits the FAIMS CV used in the feature key ([FRW-1](t05-feature-rows.md), [FID-1](t06-identifications-per-feature.md))
- [ ] [FRAGPIPE-15](findings-feature.md#fragpipe-15): one-run rows keep the run without a PSM; multi-run rows ambiguous ([FID-1](t06-identifications-per-feature.md), [FID-6](t06-identifications-per-feature.md), [FRW-6](t05-feature-rows.md))
- [ ] [FRAGPIPE-16](findings-feature.md#fragpipe-16): identified ions with 0 intensity and all-zero TMT plexes dropped ([UNI-6](t07-unidentified-features.md), [INT-6](t09-intensities-missing-values.md))
- [ ] [FRAGPIPE-17](findings-feature.md#fragpipe-17): <exp> Spectral Count and RT never read ([FID-12](t06-identifications-per-feature.md), [FRW-8](t05-feature-rows.md))
- [ ] [FRAGPIPE-18](findings-feature.md#fragpipe-18): unique from mapped-protein count; lead-only membership without --pg-file ([PRT-7](t10-protein-mapping.md), [PRT-10](t10-protein-mapping.md))
- [ ] [FRAGPIPE-20](findings-pg.md#fragpipe-20): grouped_runs re-sorted lexicographically ([PGV-5](t12-protein-groups.md))
- [ ] [FRAGPIPE-21](findings-pg.md#fragpipe-21): pg_qvalue = 1 - ProteinProphet probability; write pg_probability ([CNF-2](t04-confidence-fields.md), [PGV-8](t12-protein-groups.md))
- [ ] [FRAGPIPE-22](findings-pg.md#fragpipe-22): pg primary picked per cell among Intensity/Total Intensity/MaxLFQ ([INT-1](t09-intensities-missing-values.md))
- [ ] [FRAGPIPE-23](findings-pg.md#fragpipe-23): MBR-only pg rows not flagged ([MBR-10](t08-match-between-runs.md))
- [ ] [FRAGPIPE-27](findings-cross.md#fragpipe-27): compensation_voltage 0.0 instead of null; psm CV as string cv_param ([SPC-12](t11-spectrum-precursor-fields.md))
- [ ] [FRAGPIPE-29](findings-cross.md#fragpipe-29): no MSFragger/IonQuant versions, validator, MBR or FDR settings recorded ([META-11](t13-metadata-conformance.md), [FLT-4](t03-filtering-state.md), [MBR-9](t08-match-between-runs.md))
- [ ] [FRAGPIPE-M1](findings-pg.md#fragpipe-m1): identified channels with reporter 0 dropped ([PGV-1](t12-protein-groups.md), [INT-5](t09-intensities-missing-values.md))
- [ ] [FRAGPIPE-M2](findings-pg.md#fragpipe-m2): pg_names holds one Description instead of aligned entry names ([PGV-12](t12-protein-groups.md))
- [ ] [FRAGPIPE-M4](findings-pg.md#fragpipe-m4): ratio TMT-Integrator reports stored as log2 abundance; reject ([INT-1](t09-intensities-missing-values.md))

## Proteome Discoverer

30 items.

- [ ] [PD-02](findings-cross.md#pd-02): missing/unparseable Rank column defaults to rank 1 ([PSM-4](t01-psm-rows.md))
- [ ] [PD-03](findings-cross.md#pd-03): several rank-1 peptides of one spectrum each get its reporters as direct features ([CHI-11](t02-chimeric-spectra.md))
- [ ] [PD-04](findings-feature.md#pd-04): one spectrum from two engines summed twice in reporters and spectral_count ([PSM-2](t01-psm-rows.md), [INT-10](t09-intensities-missing-values.md), [FID-12](t06-identifications-per-feature.md))
- [ ] [PD-05](findings-cross.md#pd-05): identifying PSM by lowest q-value; psm and feature tie-breaks differ, nondeterministic ([FID-3](t06-identifications-per-feature.md), [PSM-2](t01-psm-rows.md))
- [ ] [PD-06](findings-cross.md#pd-06): psm.feature_id never filled although features aggregate those PSMs ([FID-8](t06-identifications-per-feature.md))
- [ ] [PD-07](findings-psm.md#pd-07): duplicate engine hits keep one engine's row, drop the other's scores ([PSM-2](t01-psm-rows.md))
- [ ] [PD-08](findings-psm.md#pd-08): pass_threshold = Confidence High with no declared basis ([FLT-7](t03-filtering-state.md), [FLT-8](t03-filtering-state.md))
- [ ] [PD-10](findings-feature.md#pd-10): lowest PSM q-value in feature.precursor_qvalue; move to psm_qvalue ([CNF-2](t04-confidence-fields.md))
- [ ] [PD-12](findings-feature.md#pd-12): identified-but-unquantified precursors (incl. 'Excluded by Method') dropped ([UNI-6](t07-unidentified-features.md), [UNI-7](t07-unidentified-features.md))
- [ ] [PD-13](findings-pg.md#pd-13): pg drops identified-but-unquantified channels; a test pins it ([PGV-1](t12-protein-groups.md), [META-13](t13-metadata-conformance.md))
- [ ] [PD-14](findings-pg.md#pd-14): identification-only pg as one LFQ row over all runs, sorted ([PGV-1](t12-protein-groups.md), [PGV-5](t12-protein-groups.md))
- [ ] [PD-15](findings-pg.md#pd-15): 'Exp. q-value: Combined' in unit-level pg_qvalue; move to global_qvalue ([PGV-8](t12-protein-groups.md))
- [ ] [PD-16](findings-pg.md#pd-16): experiment-wide counts repeated per file group; peptides = master only ([PGV-7](t12-protein-groups.md))
- [ ] [PD-17](findings-pg.md#pd-17): Grouped (normalized/scaled) abundance used as primary without declaration ([INT-1](t09-intensities-missing-values.md), [INT-2](t09-intensities-missing-values.md)); also [EXTERNAL-M5](findings-pg.md#external-m5)
- [ ] [PD-18](findings-feature.md#pd-18): LFQ primary = max of Precursor Abundance and Intensity over PSMs ([INT-1](t09-intensities-missing-values.md))
- [ ] [PD-19](findings-feature.md#pd-19): PIF from one PSM on reporters summed over many spectra ([INT-10](t09-intensities-missing-values.md))
- [ ] [PD-20](findings-feature.md#pd-20): shared peptide anchored on masters[0]; merged masters without Proteins.txt ([PRT-7](t10-protein-mapping.md), [PRT-8](t10-protein-mapping.md))
- [ ] [PD-21](findings-feature.md#pd-21): unique = maps to one protein even in a multi-member group; a test pins it ([PRT-10](t10-protein-mapping.md))
- [ ] [PD-22](findings-pg.md#pd-22): group members = union over every Protein Group ID of the master ([PRT-7](t10-protein-mapping.md))
- [ ] [PD-23](findings-cross.md#pd-23): is_decoy constant false; pg.contaminant false when the column is absent ([PRT-6](t10-protein-mapping.md), [FLT-10](t03-filtering-state.md))
- [ ] [PD-24](findings-cross.md#pd-24): two run-name normalisations; '.RAW' suffix splits runs ([SPC-13](t11-spectrum-precursor-fields.md))
- [ ] [PD-25](findings-psm.md#pd-25): scan = First Scan only, range lost, no scan_format ([SPC-7](t11-spectrum-precursor-fields.md), [SPC-10](t11-spectrum-precursor-fields.md))
- [ ] [PD-27](findings-psm.md#pd-27): engine score coverage partial; delta_score without engine prefix ([PSM-2](t01-psm-rows.md))
- [ ] [PD-28](findings-pg.md#pd-28): MS1-label protein abundances become duplicate LFQ pg rows ([INT-5](t09-intensities-missing-values.md), [PGV-1](t12-protein-groups.md))
- [ ] [PD-30](findings-cross.md#pd-30): pepmap built from quantified peptidoforms only ([PRT-14](t10-protein-mapping.md))
- [ ] [PD-M1](findings-feature.md#pd-m1): feature identity composite contains the picked PSM's rt ([FRW-4](t05-feature-rows.md))
- [ ] [PD-M2](findings-feature.md#pd-m2): feature PEP from the lowest-q-value PSM instead of [FID-3](t06-identifications-per-feature.md) order ([FID-3](t06-identifications-per-feature.md), [CNF-8](t04-confidence-fields.md))
- [ ] [PD-M3](findings-psm.md#pd-m3): dedup on raw modification text emits duplicate psm_id ([PSM-2](t01-psm-rows.md))
- [ ] [PD-M4](findings-feature.md#pd-m4): TMT features carry no spectral_count ([INT-11](t09-intensities-missing-values.md))
- [ ] [PD-M5](findings-feature.md#pd-m5): an entry for every reporter column whether or not the channel has a sample ([INT-5](t09-intensities-missing-values.md))

## OpenMS (consensusXML and quantms inputs)

30 items.

- [ ] [OPENMS-03](findings-psm.md#openms-03): significance_threshold unread, no filter keys written; MAY derive pass_threshold as 'other' ([FLT-1](t03-filtering-state.md), [FLT-9](t03-filtering-state.md))
- [ ] [OPENMS-04](findings-psm.md#openms-04): peptide q-values and 1.0 placeholders written to psm.qvalue; write null ([CNF-2](t04-confidence-fields.md))
- [ ] [OPENMS-05](findings-feature.md#openms-05): PSM-level q-value in feature.peptide_qvalue; move to feature.psm_qvalue ([CNF-2](t04-confidence-fields.md), [CNF-6](t04-confidence-fields.md)); also [INFRA-03](findings-feature.md#infra-03), [MOKUME-M1](findings-feature.md#mokume-m1)
- [ ] [OPENMS-06](findings-feature.md#openms-06): feature PEP/q from the run's first PID or the best PSM of another run ([FID-3](t06-identifications-per-feature.md), [FID-4](t06-identifications-per-feature.md))
- [ ] [OPENMS-07](findings-feature.md#openms-07): pids[0] defines the peptide on conflicting IDs, other scans leak in ([FID-10](t06-identifications-per-feature.md), identification_conflict)
- [ ] [OPENMS-08](findings-feature.md#openms-08): feature.scan concatenates every PID spectrum and dedup re-flattens it; one spectrum only ([FID-5](t06-identifications-per-feature.md))
- [ ] [OPENMS-10](findings-feature.md#openms-10): only the first PID is tested for hits; hitless consensus features dropped uncounted ([UNI-5](t07-unidentified-features.md), [FLT-6](t03-filtering-state.md))
- [ ] [OPENMS-11](findings-feature.md#openms-11): 'transferred' also set when IDConflictResolver removed the run's own competing PID ([MBR-2](t08-match-between-runs.md), identification_conflict)
- [ ] [OPENMS-13](findings-feature.md#openms-13): channels with it=0 omitted instead of null entries ([INT-5](t09-intensities-missing-values.md), [INT-6](t09-intensities-missing-values.md))
- [ ] [OPENMS-15](findings-feature.md#openms-15): feature observed_mz = consensus centroid or 0.0; calculated_mz = observed at charge 0 ([SPC-1](t11-spectrum-precursor-fields.md), [SPC-3](t11-spectrum-precursor-fields.md), [SPC-6](t11-spectrum-precursor-fields.md))
- [ ] [OPENMS-16](findings-psm.md#openms-16): psm.observed_mz 0.0 when absent; write null ([SPC-1](t11-spectrum-precursor-fields.md))
- [ ] [OPENMS-17](findings-psm.md#openms-17): surrogate hash scans; spectrum=/index references undeclared ([SPC-7](t11-spectrum-precursor-fields.md), [SPC-8](t11-spectrum-precursor-fields.md), [SPC-9](t11-spectrum-precursor-fields.md))
- [ ] [OPENMS-18](findings-psm.md#openms-18): non-snake_case score names, '_2' collision suffixes, engine meta scores dropped ([PSM-2](t01-psm-rows.md), [CNF-5](t04-confidence-fields.md)); also [INFRA-31](findings-psm.md#infra-31), [VALIDATOR-03](findings-psm.md#validator-03)
- [ ] [OPENMS-19](findings-psm.md#openms-19): later-PID duplicates of one spectrum/peptidoform dropped with their scores; merge ([PSM-2](t01-psm-rows.md))
- [ ] [OPENMS-20](findings-psm.md#openms-20): psm.feature_id links other peptidoforms and decoys; link supporting PSMs only ([FID-7](t06-identifications-per-feature.md))
- [ ] [OPENMS-21](findings-cross.md#openms-21): --no-include-unassigned-psms drops PSMs that features reference ([FID-8](t06-identifications-per-feature.md))
- [ ] [OPENMS-22](findings-cross.md#openms-22): is_decoy from target_decoy in psm/feature, from accession prefix in pg ([PRT-4](t10-protein-mapping.md), [PRT-5](t10-protein-mapping.md))
- [ ] [OPENMS-23](findings-pg.md#openms-23): pg evidence from any rank, decoys and the consensus feature's runs ([PGV-2](t12-protein-groups.md))
- [ ] [OPENMS-24](findings-pg.md#openms-24): QPX-computed sum written as pg primary when no abundance ([INT-4](t09-intensities-missing-values.md))
- [ ] [OPENMS-25](findings-pg.md#openms-25): pg.global_qvalue = min member-protein q-value; move to protein_global_qvalue ([PGV-8](t12-protein-groups.md))
- [ ] [OPENMS-26](findings-pg.md#openms-26): analysis-wide peptide/feature counts repeated on every unit row ([PGV-7](t12-protein-groups.md))
- [ ] [OPENMS-27](findings-feature.md#openms-27): anchor_protein without membership for ungrouped single-accession peptides ([PRT-11](t10-protein-mapping.md))
- [ ] [OPENMS-28](findings-psm.md#openms-28): unknown mass shifts written as unsigned 15-decimal residue masses (SPC-14b)
- [ ] [OPENMS-29](findings-feature.md#openms-29): dimethyl/mTRAQ multiplet peptidoform lacks the lightest label; convert.md contradicts ([FRW-10](t05-feature-rows.md))
- [ ] [OPENMS-31](findings-cross.md#openms-31): OpenMS -out_qpx upgrade passes rows through unnormalised and keeps only 4 footer keys ([META-4](t13-metadata-conformance.md))
- [ ] [OPENMS-M1](findings-psm.md#openms-m1): PEP primary score lost; score type 'fdr' relabelled q-value ([CNF-2](t04-confidence-fields.md), [CNF-7](t04-confidence-fields.md))
- [ ] [OPENMS-M2](findings-feature.md#openms-m2): precursor-runs with an accepted ID but no positive element get no identification-only row ([UNI-6](t07-unidentified-features.md), [UNI-7](t07-unidentified-features.md))
- [ ] [OPENMS-M3](findings-feature.md#openms-m3): element/PID RT clock unverified though alignment declared; aligned times go to rt_aligned ([FRW-7](t05-feature-rows.md), [SPC-11](t11-spectrum-precursor-fields.md))
- [ ] [OPENMS-M4](findings-psm.md#openms-m4): PIDs without resolvable run or spectrum reference dropped without count ([FLT-6](t03-filtering-state.md))
- [ ] [OPENMS-M5](findings-cross.md#openms-m5): OpenMS -out_qpx keeps unidentified seeds (charge 0) in feature; move to unidentified_feature ([UNI-3](t07-unidentified-features.md), [UNI-8](t07-unidentified-features.md))

## mzIdentML

13 items.

- [ ] [EXTERNAL-M1](findings-psm.md#external-m1): 'final PSM list' MS:1002439 ignored; every SIL converted ([PSM-9](t01-psm-rows.md))
- [ ] [MZID_CDAP-03](findings-psm.md#mzid_cdap-03): PSM-level p-value MS:1002352 read as PEP ([CNF-7](t04-confidence-fields.md)); also [INFRA-09](findings-psm.md#infra-09)
- [ ] [MZID_CDAP-04](findings-psm.md#mzid_cdap-04): target+decoy peptides flagged decoy (any instead of all) ([PRT-4](t10-protein-mapping.md))
- [ ] [MZID_CDAP-07](findings-pg.md#mzid_cdap-07): pg is_decoy from accession prefix while psm uses isDecoy ([PRT-5](t10-protein-mapping.md))
- [ ] [MZID_CDAP-08](findings-psm.md#mzid_cdap-08): duplicate composites; ids not reproducible from the footer composite ([PSM-2](t01-psm-rows.md), [PSM-10](t01-psm-rows.md)); also [MZID_CDAP-09](findings-cross.md#mzid_cdap-09), [CRITIC-M1](findings-psm.md#critic-m1)
- [ ] [MZID_CDAP-10](findings-psm.md#mzid_cdap-10): SIP <Threshold> never read; pass_threshold basis unrecorded ([FLT-8](t03-filtering-state.md)); also [INFRA-M3](findings-psm.md#infra-m3)
- [ ] [MZID_CDAP-12](findings-psm.md#mzid_cdap-12): psm.qvalue chosen by cvParam order between engine and rescorer ([CNF-4](t04-confidence-fields.md))
- [ ] [MZID_CDAP-19](findings-psm.md#mzid_cdap-19): nativeID parsing keeps one number, drops vendor components ([SPC-7](t11-spectrum-precursor-fields.md))
- [ ] [MZID_CDAP-20](findings-cross.md#mzid_cdap-20): mixed index/scan references leave scan_format undeclared ([SPC-7](t11-spectrum-precursor-fields.md))
- [ ] [MZID_CDAP-27](findings-pg.md#mzid_cdap-27): PAG passThreshold ignored, dataset-wide runs, inflated counts ([PGV-1](t12-protein-groups.md), [PGV-7](t12-protein-groups.md), [PGV-9](t12-protein-groups.md))
- [ ] [MZID_CDAP-28](findings-cross.md#mzid_cdap-28): peptides.peptide_count counts PeptideHypotheses vs sequences ([PGV-7](t12-protein-groups.md))
- [ ] [MZID_CDAP-32](findings-psm.md#mzid_cdap-32): positions of the first protein only, in cv_params ([PRT-1](t10-protein-mapping.md))
- [ ] [MZID_CDAP-M2](findings-cross.md#mzid_cdap-m2): pepmap built from every SII incl. rank>1 and failing hits ([PRT-14](t10-protein-mapping.md))

## CDAP

15 items.

- [ ] [MZID_CDAP-05](findings-cross.md#mzid_cdap-05): XXX_ decoy prefix stripped ([PRT-3](t10-protein-mapping.md))
- [ ] [MZID_CDAP-06](findings-pg.md#mzid_cdap-06): decoy features merged into target protein groups ([PGV-10](t12-protein-groups.md), [PRT-3](t10-protein-mapping.md))
- [ ] [MZID_CDAP-13](findings-feature.md#mzid_cdap-13): best-PSM MS-GF+ q-values in feature precursor/peptide q-value columns ([CNF-2](t04-confidence-fields.md))
- [ ] [MZID_CDAP-14](findings-feature.md#mzid_cdap-14): feature rollup keyed on raw sequence, not the canonical peptidoform ([FRW-4](t05-feature-rows.md))
- [ ] [MZID_CDAP-15](findings-feature.md#mzid_cdap-15): PSM count only as cv_param psm_count; TMT and LFQ pick rules differ ([INT-10](t09-intensities-missing-values.md), [FID-3](t06-identifications-per-feature.md), [FID-12](t06-identifications-per-feature.md))
- [ ] [MZID_CDAP-16](findings-feature.md#mzid_cdap-16): AmbiguousMatch ties dropped from features and pg evidence ([FID-10](t06-identifications-per-feature.md), [PGV-2](t12-protein-groups.md))
- [ ] [MZID_CDAP-17](findings-psm.md#mzid_cdap-17): psm.rt = precursor half-elution time; 0.00 sentinel kept ([SPC-11](t11-spectrum-precursor-fields.md), [FRW-8](t05-feature-rows.md))
- [ ] [MZID_CDAP-18](findings-psm.md#mzid_cdap-18): missing observed m/z and charge written as 0 ([SPC-1](t11-spectrum-precursor-fields.md), [SPC-6](t11-spectrum-precursor-fields.md))
- [ ] [MZID_CDAP-25](findings-pg.md#mzid_cdap-25): mapping-set groups undeclared; peptides listed for the anchor only ([PGV-6](t12-protein-groups.md), [PGV-7](t12-protein-groups.md))
- [ ] [MZID_CDAP-26](findings-pg.md#mzid_cdap-26): PSM- and peptide-level q-values written as pg scores ([PGV-8](t12-protein-groups.md))
- [ ] [MZID_CDAP-29](findings-psm.md#mzid_cdap-29): psm.feature_id never filled although features aggregate PSMs ([FID-8](t06-identifications-per-feature.md))
- [ ] [MZID_CDAP-M3](findings-cross.md#mzid_cdap-m3): contaminant entries lose identity; pg.contaminant hard-coded false ([PRT-6](t10-protein-mapping.md), [PGV-10](t12-protein-groups.md))
- [ ] [MZID_CDAP-M4](findings-pg.md#mzid_cdap-m4): QPX reporter sum written as pg primary ([INT-4](t09-intensities-missing-values.md))
- [ ] [MZID_CDAP-M5](findings-pg.md#mzid_cdap-m5): total_sequences counts features, not sequences ([PGV-7](t12-protein-groups.md))
- [ ] [MZID_CDAP-M6](findings-psm.md#mzid_cdap-m6): phosphoRS localisation ignored in peptidoform and candidates ([SPC-15](t11-spectrum-precursor-fields.md))

## quantms-MSstats

3 items.

- [ ] [INFRA-14](findings-feature.md#infra-14): aborts on NA intensity, keeps zeros, provenance misstates it ([INT-6](t09-intensities-missing-values.md), [INT-7](t09-intensities-missing-values.md))
- [ ] [INFRA-27](findings-feature.md#infra-27): id_run_file_name = run_file_name on every row, even for requantified runs ([FID-6](t06-identifications-per-feature.md), [MBR-2](t08-match-between-runs.md))
- [ ] [INFRA-28](findings-feature.md#infra-28): dataset-wide unique with per-row union pg_accessions ([PRT-7](t10-protein-mapping.md), [PRT-10](t10-protein-mapping.md))

## Writers and converter infrastructure

8 items.

- [ ] [INFRA-01](findings-feature.md#infra-01): denylist fill copies PSM-level q-value names into feature.precursor_qvalue ([CNF-6](t04-confidence-fields.md)); also [INFRA-05](findings-psm.md#infra-05)
- [ ] [INFRA-08](findings-cross.md#infra-08): typed q-value and its additional_scores copy may disagree, unchecked ([CNF-6](t04-confidence-fields.md))
- [ ] [INFRA-11](findings-feature.md#infra-11): identification_type unchecked on write_batch, no YAML enum ([MBR-12](t08-match-between-runs.md), [META-9](t13-metadata-conformance.md)); also [VALIDATOR-09](findings-feature.md#validator-09)
- [ ] [INFRA-12](findings-cross.md#infra-12): record/table/DataFrame paths validate differently; violations persisted as warnings ([META-9](t13-metadata-conformance.md)); also [VALIDATOR-01](findings-cross.md#validator-01), [VALIDATOR-M4](findings-cross.md#validator-m4)
- [ ] [INFRA-17](findings-cross.md#infra-17): normalize_peptidoform keeps non-canonical terminal placement; add hyphens and test vectors (peptidoform.md rule 4)
- [ ] [INFRA-18](findings-feature.md#infra-18): unidentified-row ids derived but identity_source says provided/mixed ([UNI-4](t07-unidentified-features.md)); also [INFRA-19](findings-feature.md#infra-19)
- [ ] [INFRA-21](findings-cross.md#infra-21): footer keys added after the first row group missing from ARROW:schema ([META-3](t13-metadata-conformance.md))
- [ ] [INFRA-30](findings-cross.md#infra-30): shared helper drops rows that fail to transform, logging only ([FLT-6](t03-filtering-state.md))

## Validator

17 items.

- [ ] [INFRA-16](findings-cross.md#infra-16): canonical peptidoform names/decimals not enforced ([SPC-15](t11-spectrum-precursor-fields.md)); also [VALIDATOR-14](findings-feature.md#validator-14)
- [ ] [VALIDATOR-02](findings-psm.md#validator-02): identity_not_checked fails --conformance for spec-sanctioned provided ids ([PSM-10](t01-psm-rows.md), [META-8](t13-metadata-conformance.md))
- [ ] [VALIDATOR-10](findings-feature.md#validator-10): precursor_run uniqueness not enforced ([FRW-12](t05-feature-rows.md))
- [ ] [VALIDATOR-11](findings-feature.md#validator-11): Writer accepts '' peptidoform rows that --semantic rejects ([UNI-1](t07-unidentified-features.md))
- [ ] [VALIDATOR-12](findings-feature.md#validator-12): feature scan_format never checked ([SPC-7](t11-spectrum-precursor-fields.md), [FID-5](t06-identifications-per-feature.md))
- [ ] [VALIDATOR-13](findings-psm.md#validator-13): empty psm scan passes ([SPC-8](t11-spectrum-precursor-fields.md))
- [ ] [VALIDATOR-15](findings-cross.md#validator-15): psm.feature_id checked for existence only ([FID-7](t06-identifications-per-feature.md))
- [ ] [VALIDATOR-22](findings-feature.md#validator-22): no ranges for PIF, quant_ms_level, channel_details.qvalue; NaN additional intensities pass ([INT-14](t09-intensities-missing-values.md), [CNF-9](t04-confidence-fields.md))
- [ ] [VALIDATOR-23](findings-feature.md#validator-23): unique, pg_group_id hash and anchor membership unchecked ([PRT-10](t10-protein-mapping.md), [PRT-13](t10-protein-mapping.md))
- [ ] [VALIDATOR-24](findings-feature.md#validator-24): allowed unlinked shared peptides fail --conformance ([META-8](t13-metadata-conformance.md))
- [ ] [VALIDATOR-28](findings-pg.md#validator-28): null-label pg row with an intensity passes ([PGV-1](t12-protein-groups.md), [INT-14](t09-intensities-missing-values.md))
- [ ] [VALIDATOR-30](findings-cross.md#validator-30): Merge drops scan_format, spectrum_content and creator from footers ([META-4](t13-metadata-conformance.md))
- [ ] [VALIDATOR-32](findings-psm.md#validator-32): lengths of the parallel PSM spectral arrays unchecked (psm.md:290)
- [ ] [VALIDATOR-M1](findings-cross.md#validator-m1): peptidoforms missing fixed mods/labels raise only a mass_error warning ([SPC-3](t11-spectrum-precursor-fields.md))
- [ ] [VALIDATOR-M2](findings-feature.md#validator-m2): rt outside [rt_start, rt_stop] and negative times pass ([FRW-9](t05-feature-rows.md), [FRW-12](t05-feature-rows.md))
- [ ] [VALIDATOR-M5](findings-feature.md#validator-m5): id_run_file_name never validated ([FID-6](t06-identifications-per-feature.md))
- [ ] [VALIDATOR-M6](findings-cross.md#validator-m6): charge sign/zero unchecked; m/z checks skip charge <= 0 ([SPC-6](t11-spectrum-precursor-fields.md))

## Documentation and tests

10 items.

- [ ] [FRAGPIPE-24](findings-cross.md#fragpipe-24): FragPipe tables in tool-mappings/converter-coverage contradict the code ([META-12](t13-metadata-conformance.md))
- [ ] [FRAGPIPE-M3](findings-cross.md#fragpipe-m3): source_column_name names columns the code does not read ([META-12](t13-metadata-conformance.md))
- [ ] [INFRA-04](findings-feature.md#infra-04): coverage tables contradict write-time q-value derivation for CDAP and Spectronaut ([META-12](t13-metadata-conformance.md)); also [CRITIC-33](findings-feature.md#critic-33)
- [ ] [INFRA-M4](findings-feature.md#infra-m4): producer-guide names a non-existent 'imputation' flag ([INT-9](t09-intensities-missing-values.md), [META-12](t13-metadata-conformance.md)); also [SPEC-M6](findings-cross.md#spec-m6)
- [ ] [MZID_CDAP-23](findings-feature.md#mzid_cdap-23): coverage says CDAP writes a null list for unquantified features; it writes null entries ([META-12](t13-metadata-conformance.md))
- [ ] [PD-33](findings-cross.md#pd-33): PD tests pin spec-divergent behaviour; no PD fixture covers chimeras, multi-engine runs or MBR ([META-13](t13-metadata-conformance.md), [META-14](t13-metadata-conformance.md))
- [ ] [SPEC-03](findings-feature.md#spec-03): feature.md examples write pg_accessions as strings that FeatureWriter rejects ([PRT-7](t10-protein-mapping.md))
- [ ] [SPEC-31](findings-cross.md#spec-31): view pages mark optional 1.3 columns Required 'no' (= must be present) ([META-12](t13-metadata-conformance.md))
- [ ] [SPEC-35](findings-cross.md#spec-35): peptidoform.md cites a nonexistent Peptide view and universal key (T11 docs list)
- [ ] [SPEC-M5](findings-psm.md#spec-m5): coverage says MaxQuant psm.feature_id is unimplemented, but it is filled ([META-12](t13-metadata-conformance.md)); also [MOKUME-M4](findings-cross.md#mokume-m4), [CRITIC-M6](findings-cross.md#critic-m6)
