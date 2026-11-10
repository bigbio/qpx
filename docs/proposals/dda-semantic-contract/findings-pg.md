# Findings: protein-group view

Part of the annex of the [DDA semantic contract proposal](index.md). Each entry is one verified finding from the audit of qpx `dev` at commit `32cbbdc` (specification 1.3). A skeptic re-checked every finding against the code; "partly confirmed" entries are true only as stated in their correction. `audit-run:` evidence refers to converting the repository's own fixtures (`tests/examples`) with that commit during the audit; those outputs are not committed.

The "Suggested" line of an entry is the auditor's suggestion at the time of the audit, not a decision. Some suggestions name footer keys or a separate view for unidentified features, which the proposal no longer uses. What was decided is in the [strict profile](strict-profile.md) and the topic sections.

48 findings, grouped by audit area.

## OpenMS converters

### OPENMS-23: pg rows/groups from any hit (any rank, decoy, unfiltered), not FDR-passing PSMs {#openms-23}

*pg* · diverges from spec · confirmed · checklist G1, G3

- **Now:** Groups are the indistinguishable groups plus a singleton for every ProteinHit (pg_adapter.py:185-206), with no protein-FDR check. A unit gets a row if any PID hit of any rank or score, assigned or unassigned and including decoys, maps to a member (104-123, 463-466), or if it has a positive quantity, which includes transferred features.
- **Spec:** pg.md:126-128: identification evidence requires at least one PSM that passed the FDR threshold. pg.md:123: passed FDR but no quantity gives a null row. #342 decided to keep decoy groups.
- **Reader impact:** Null-intensity pg rows assert identification from low-confidence, rank-n or decoy hits, so protein identification counts are inflated.
- **Suggested:** Base evidence on FDR-passing rank-1 PSMs, or declare in the footer that pg evidence is unfiltered.
- **Evidence:** `qpx/converters/openms_consensus/pg_adapter.py:104-123,185-206,432-476`; `docs/spec/pg.md:123,126-128`

### OPENMS-24: pg intensity: OpenMS abundance or qpx-computed sum; contradicts 'not computed by QPX' {#openms-24}

*pg* · contradictory · partly confirmed (see correction) · checklist G2, X1

- **Now:** intensity = OpenMS fraction_group abundance, with 0 treated as null and the method unstated beyond 'openms_protein_abundance'. Otherwise it is a qpx sum of the unique-peptide feature intensities by unmodified sequence, over all charges and modifications, including transferred and decoy-unfiltered rows (pg_adapter.py:323-340,414-429; feature_dedup.py:350-354; protein_quant.py:171). Uniqueness comes from all hits, rank-2 included. In an input that reports abundance, a group without an OpenMS value has null intensity even when the roll-up exists.
- **Spec:** intensities.md:20 and feature.md:303: additional_intensities are read from the tool, not computed by QPX. pg.md:240 lists unnormalized_unique_peptide_sum. convert.md:1096 says intensity is always the sum, contradicting convert.md:1017-1026.
- **Correction:** The behaviour is confirmed, but the cited protein_quant.py:155-180/171 and feature_dedup.py:350-354 do not exist (the files have 120 and 166 lines). The roll-up is pg_adapter.py:323-340 fed by feature_dedup.py:105-109, and zero becomes null at protein_quant.py:98,111. The spec itself sanctions the sum (pg.md:240, converter-coverage.md:107), so the contradiction is inside the spec (intensities.md:3,20; pg.md:54) and in convert.md:1096 vs :1017-1026.
- **Reader impact:** Readers cannot compare pg intensity across OpenMS files (two different methods) or with other tools. The documentation states different sources.
- **Suggested:** pg.intensity holds only tool-reported quantities. QPX-derived roll-ups go to a clearly derived name or a transform output. Fix convert.md:1096.
- **Evidence:** `qpx/converters/openms_consensus/pg_adapter.py:323-340,414-429`; `qpx/converters/openms_consensus/protein_quant.py:155-180`; `docs/spec/intensities.md:20`; `docs/spec/feature.md:303`; `docs/spec/pg.md:52-58,240`; `docs/guide/convert.md:1017-1026,1096`

### OPENMS-25: pg.global_qvalue = min member-protein q-value across all ProteinIdentifications {#openms-25}

*pg* · ambiguous · confirmed · checklist G3

- **Now:** global_qvalue is the minimum ProteinHit q-value over the members, taken over every ProteinIdentification, when the protein score_type is a q-value (pg_adapter.py:166-182,276-277,286-297). pg_qvalue is never filled. feature.pg_global_qvalue uses the same value.
- **Spec:** pg.md:92: global q-value of the protein group at the experiment level. converter-coverage.md:104 documents 'lowest member-protein q-value'.
- **Reader impact:** A protein-level q-value (possibly from several separate inferences) is read as a group-level FDR. Filters at 1% differ from MaxQuant or DIA-NN group q-values.
- **Suggested:** The contract defines global_qvalue as a group-level FDR. Put member-protein minima under a distinct name (e.g. additional_scores 'protein_qvalue_min').
- **Evidence:** `qpx/converters/openms_consensus/pg_adapter.py:166-182,260-297`; `docs/spec/pg.md:92`; `docs/spec/converter-coverage.md:104`

### OPENMS-26: pg peptide/feature counts are global over all inputs, repeated on every unit row {#openms-26}

*pg* · ambiguous · partly confirmed (see correction) · checklist G1, G2

- **Now:** Counts are computed once per group from all PID hits in every input (any rank, decoys, unassigned) and copied onto every (unit,label) row (pg_adapter.py:366-387,448,472). feature_counts counts identified (peptidoform,charge) pairs, not quantified features.
- **Spec:** pg.md:81-86: counts 'for this protein group in this file'. pg.yaml:92-97 does not give the scope.
- **Correction:** pg.md:81-86 defines the counts 'in this file', and intensities.md:206 lists pg.peptide_counts as dataset-scope (per unit only for MaxQuant), so per-file counts repeated on rows conform. feature_counts as peptide-charge pairs matches pg.md:84. What remains: the scope differs by producer, and the counts come from unfiltered, rank-n and decoy hits (overlaps OPENMS-23).
- **Reader impact:** Readers taking counts as per-unit evidence (e.g. a minimum of 2 peptides per sample) apply global counts.
- **Suggested:** Define the scope of the counts as per row (unit) and from FDR-passing evidence. Global counts go in a separate field.
- **Evidence:** `qpx/converters/openms_consensus/pg_adapter.py:366-387,447-476`; `docs/spec/pg.md:81-86`; `qpx/core/data/schemas/pg.yaml:92-97`

## MaxQuant converter

### MAXQUANT-17: MaxQuant MaxLFQ is named 'lfq'; other producers use 'maxlfq' {#maxquant-17}

*pg* · contradictory · confirmed · checklist G2, X1

- **Now:** pg_adapter.py:606 names 'LFQ intensity <exp>' as additional 'lfq'. DIA-NN and FragPipe MaxLFQ values go to 'maxlfq'.
- **Spec:** pg.md:228 lists 'maxlfq' for MaxQuant too. pg.md:229 lists 'lfq' for MaxQuant. pg.md:307 maps LFQ intensity to lfq ('MaxLFQ normalised'). pg.md:343,391 map DIA-NN and FragPipe to maxlfq.
- **Reader impact:** A query for maxlfq across producers misses MaxQuant's MaxLFQ values.
- **Suggested:** Use one name, `maxlfq`, for MaxLFQ values from any producer.
- **Evidence:** `qpx/converters/maxquant/pg_adapter.py:598-610`; `docs/spec/pg.md:228-229,307,343,391`

### MAXQUANT-18: pg_qvalue gets the experiment-level Q-value, but the spec defines it as run-level {#maxquant-18}

*pg* · contradictory · confirmed · checklist G3, D4

- **Now:** pg_qvalue and global_qvalue both get the proteinGroups 'Q-value' (pg_adapter.py:458,483-484). converter-coverage.md:300 documents this duplication.
- **Spec:** pg.md:92-93: global_qvalue is experiment-level and pg_qvalue is run-level.
- **Reader impact:** Run-level protein filtering on pg_qvalue silently applies an experiment-level FDR, for MaxQuant only.
- **Suggested:** Leave pg_qvalue null for MaxQuant. The contract forbids copying a q-value across levels.
- **Evidence:** `qpx/converters/maxquant/pg_adapter.py:456-490`; `docs/spec/pg.md:92-93`; `docs/spec/converter-coverage.md:300`

### MAXQUANT-19: Protein-group PEP in pg.md is never read {#maxquant-19}

*pg* · diverges from spec · confirmed · checklist G3

- **Now:** Only 'Score' goes into additional_scores (pg_adapter.py:465-474). A 'PEP' column is ignored. The fixture proteinGroups.txt has no PEP, so I could not verify whether current MaxQuant writes one.
- **Spec:** pg.md:322 maps `PEP` to additional_scores posterior_error_probability. pg.md:258.
- **Reader impact:** Protein PEP is lost when the source has it.
- **Suggested:** Map PEP when present, or drop it from pg.md.
- **Evidence:** `qpx/converters/maxquant/pg_adapter.py:464-474`; `docs/spec/pg.md:258,322`

### MAXQUANT-20: MBR-only protein quantities are flagged only by free-text cv_params {#maxquant-20}

*pg* · unspecified · confirmed · checklist G1, F4

- **Now:** 'Identification type <exp>' is copied as cv_params text (pg_adapter.py:403-406). Full fixture: 5,687/27,072 pg rows are 'By matching' (quantity from MBR only) and are emitted with intensity. 75 'By matching' cells with no intensity get no row (pg_adapter.py:589). is_imputed stays null.
- **Spec:** pg.md:129-134 allows MBR-only rows with a quantity. No typed field marks transfer at pg level (pg.md:291 uses cv_params).
- **Reader impact:** Readers cannot type-safely separate MS/MS-backed protein quantities from MBR-only ones (mokume#117).
- **Suggested:** Add a typed pg identification_type (direct/transferred) mirroring feature.
- **Evidence:** `qpx/converters/maxquant/pg_adapter.py:398-415,581-596`; `docs/spec/pg.md:129-134,291`; `audit-run:maxquant/full/maxquant.pg.parquet`

## FragPipe converter

### FRAGPIPE-20: grouped_runs sorted lexicographically, contrary to the fraction-order rule {#fragpipe-20}

*pg* · diverges from spec · confirmed · checklist G1, F1

- **Now:** normalize_experiment_mapping sorts each experiment's runs (pg_adapter.py:96), and plex_runs uses ORDER BY run (isobaric.py:249-251). Verified: grouped_runs=['runB_F1','runB_F10','runB_F2'] and ['plex1_F1','plex1_F10','plex1_F2']. The same sorted list chooses feature.run_file_name.
- **Spec:** pg.md:414: grouped_runs 'is stored in fraction order — sorting is not applied, because it would ... lexicographically misorder names like F1, F10, F2'.
- **Reader impact:** Fraction order is lost for FragPipe pg rows, and readers relying on grouped_runs order (fraction index) get a wrong order.
- **Suggested:** Keep experiment_annotation.tsv and psm.tsv order (first appearance), or read the fraction column. Validate the order when one is known.
- **Evidence:** `qpx/converters/fragpipe/pg_adapter.py:87-100`; `qpx/converters/fragpipe/isobaric.py:246-252`; `docs/spec/pg.md:414`; `run: out/fp.pg.parquet grouped_runs [runB_F1, runB_F10, runB_F2]`

### FRAGPIPE-21: pg_qvalue filled with 1 - ProteinProphet probability, repeated per experiment {#fragpipe-21}

*pg* · contradictory · confirmed · checklist G3, D4, P5

- **Now:** pg_qvalue is the first of 'Protein FDR', 'Protein Q-Value' or 'Combined Protein FDR', else 1 - 'Protein Probability' (pg_adapter.py:406-421, column_mappings.yaml:34). No real FragPipe file seen has an FDR column, so in practice the value is 1 - the combined probability. It is identical on every experiment row. global_qvalue is null, and the probability is also kept in additional_scores.
- **Spec:** pg.md:93: pg_qvalue is the protein group q-value at the RUN level. pg.md:405-406: FragPipe reports no PG q-value; set global_qvalue null and keep the probability in additional_scores. converter-coverage.md:300 documents 1 - probability.
- **Reader impact:** A posterior error estimate is presented as a run-level q-value. 'pg_qvalue <= 0.01' then means 'probability >= 0.99', not 1% FDR. The spec pages disagree.
- **Suggested:** Leave pg_qvalue null for FragPipe (the PEP-like value can go in additional_scores as protein_pep). Declare in provenance or footer that the input was filtered at protein FDR X.
- **Evidence:** `qpx/converters/fragpipe/pg_adapter.py:406-421`; `qpx/converters/column_mappings.yaml:31-34`; `qpx/core/data/schemas/pg.yaml:57-61`; `docs/spec/pg.md:93`; `docs/spec/pg.md:405-406`; `docs/spec/converter-coverage.md:300`; `tests/converters/test_fragpipe_converter.py:196-205`

### FRAGPIPE-22: Primary pg intensity mixes Intensity, Total Intensity and MaxLFQ per row and experiment {#fragpipe-22}

*pg* · contradictory · confirmed · checklist G2, X1

- **Now:** The primary intensity is the first positive of '<exp> Intensity', '<exp> Total Intensity' and '<exp> MaxLFQ Intensity', chosen per row and per experiment (pg_adapter.py:71). One column can hold razor or top-N intensity, a total that includes shared peptides, or MaxLFQ. Tests pin the MaxLFQ fallback (test_fragpipe_converter.py:208-231). The 'Razor Intensity', 'Razor Spectral Count' and 'Total Intensity' suffixes are recognised for experiment detection but never stored (pg_adapter.py:46-51).
- **Spec:** pg.md:394: '[exp] Intensity' is the primary razor peptide intensity, with '[exp] Total Intensity' in older reports. anndata.md:147: X = the tool's unnormalised quantity ('FragPipe razor Intensity'). pg.md:228: MaxLFQ goes in additional_intensities. converter-coverage.md:301 documents the fallback.
- **Reader impact:** Normalised MaxLFQ and raw intensities can sit side by side in pg.intensity of one file. Old reports use shared-inclusive totals while the docs claim razor.
- **Suggested:** Choose one source per file (never a per-cell fallback), declare it (e.g. footer pg_intensity_source), store the other quantities as named additional_intensities, and use 'Razor Intensity' for old reports.
- **Evidence:** `qpx/converters/fragpipe/pg_adapter.py:46-78`; `tests/converters/test_fragpipe_converter.py:208-231`; `docs/spec/pg.md:394`; `docs/spec/pg.md:228`; `docs/spec/anndata.md:147`; `docs/spec/converter-coverage.md:301`

### FRAGPIPE-23: pg rows from MBR-only quantity indistinguishable; identified-without-quantity rows kept {#fragpipe-23}

*pg* · unspecified · partly confirmed (see correction) · checklist G1, F4

- **Now:** A row is written when the intensity is > 0 or any spectral count is > 0 (pg_adapter.py:71-78). A zero intensity with counts gives a null intensity (verified contaminant row). An experiment quantified only through MBR (spectral count 0, intensity > 0) gives an ordinary row with spectral_count 0.0 (verified experiment C), with no flag. Spectral-count zeros are kept as 0.0.
- **Spec:** pg.md:118 onward and #280: a null intensity means evidence without a quantity, and no row means no evidence. The spec says nothing about a quantity that rests only on transferred identifications (andes#218: transfers never create protein evidence).
- **Correction:** The spec is not silent. pg.md:131-135 explicitly allows rows with a quantity but without identification evidence ('for example from match-between-runs') and forbids null rows there, and FragPipe complies. What is missing is only a flag that tells MBR-only rows apart; the status should be unspecified for the flag.
- **Reader impact:** A reader cannot exclude protein quantities supported only by match-between-runs, and the presence of a row is read as evidence of identification in that sample.
- **Suggested:** The contract should mark MBR-only pg rows, e.g. is_transferred, or identification evidence counts from MS/MS only (spectral_count=0 as an explicit rule).
- **Evidence:** `qpx/converters/fragpipe/pg_adapter.py:54-78`; `docs/spec/pg.md:118`; `run: out/fp.pg.parquet unit C intensity 3000, spectral_count 0.0; P02768 intensity null`

### FRAGPIPE-M1: TMT pg drops identified channels whose reporter value is 0 {#fragpipe-m1}

*pg* · diverges from spec · found by the verifier · checklist G1, D6, F6

- **Now:** _channel_records writes only channels whose intensity is non-null (pg_adapter.py:381-383), and channel_intensities turns 0 into null (isobaric.py:44-50). Verified: plex1 protein.tsv with sample_b=0 gives only a TMT126 row, with no TMT127N row, while the feature keeps a TMT127N null entry. A test pins this (test_fragpipe_isobaric.py:93-96).
- **Spec:** pg.md:109: one row per label in each unit where the group has identification evidence, including labels without a quantity. pg.md:129-131 and 143-144: a zero sentinel becomes a null intensity when the group is identified.
- **Reader impact:** In FragPipe TMT pg, a missing channel looks like 'not identified' rather than 'identified, unquantified'. The feature and pg views disagree, and missingness differs from other isobaric producers.
- **Suggested:** Write a null-intensity row for every annotated channel of each group in protein.tsv. The contract should state label completeness per unit for isobaric pg.
- **Evidence:** `qpx/converters/fragpipe/pg_adapter.py:381-383`; `qpx/converters/fragpipe/isobaric.py:44-50`; `tests/converters/test_fragpipe_isobaric.py:93-96`; `docs/spec/pg.md:109`; `docs/spec/pg.md:129-131`; `run: audit-run:fragpipe/tmt/out/plex1.pg.parquet has only the TMT126 row`

### FRAGPIPE-M2: pg_names holds one Description, not aligned entry names; spec pages disagree {#fragpipe-m2}

*pg* · contradictory · found by the verifier · checklist X1, G3

- **Now:** pg_names = [Description] (or 'Protein Description'), one element whatever the group size (pg_adapter.py:433-437). Verified: the 2-member group [P12345,Q22222] has pg_names ['desc']. Entry names, available from 'Entry Name' or the sp|..|NAME headers, are not used.
- **Spec:** pg.md:35: UniProt entry names aligned with pg_accessions (A1BG_HUMAN). pg.md:373, tool-mappings.md:103 and converter-coverage.md:297: FragPipe Description, MaxQuant 'Protein names', DIA-NN Protein.Names.
- **Reader impact:** pg_names holds entry names, descriptions or protein names depending on the producer. Index-aligned use with pg_accessions breaks for multi-member FragPipe groups.
- **Suggested:** The contract should fix pg_names as aligned entry names, null where unknown, and move descriptions to a separate field or cv_param.
- **Evidence:** `qpx/converters/fragpipe/pg_adapter.py:433-437`; `docs/spec/pg.md:35`; `docs/spec/pg.md:373`; `docs/spec/tool-mappings.md:103`; `docs/spec/converter-coverage.md:297`; `run: work/verify-fragpipe/out/fp.pg.parquet pg_names ['desc'] for 2 accessions`

### FRAGPIPE-M4: Any protein-level TMT-Integrator report is stored as log2 abundance, ratios included {#fragpipe-m4}

*pg* · diverges from spec · found by the verifier · checklist G2, D7

- **Now:** load_tmt_report only checks for 'Gene' and 'Index' columns (isobaric.py:62-73), and a ratio_protein report has them too. Verified: a ratio report with -0.75 was stored as tmt_integrator_log2_abundance and passed qpxc validate. Values are joined on Protein ID and placed only on channels whose protein.tsv value is >0.
- **Spec:** pg.md:244: tmt_integrator_log2_abundance is 'TMT-Integrator's normalized log2 abundance'. intensities.md:181: each intensity_name has one meaning. Nothing says the value is normalised across plexes rather than per row's unit.
- **Reader impact:** Log2 ratios to a reference can be read as log2 abundances. Cross-plex normalised values sit on per-plex rows without a stated scope.
- **Suggested:** Reject non-abundance reports, for example by file name or the ReferenceIntensity semantics. The contract should give each additional intensity its scale (log/linear) and normalisation scope.
- **Evidence:** `qpx/converters/fragpipe/isobaric.py:62-73`; `qpx/converters/fragpipe/pg_adapter.py:379-399`; `docs/spec/pg.md:244`; `run: work/verify-fragpipe/tmt/out/plex1.pg.parquet tmt_integrator_log2_abundance=-0.75, validate VALID`

## Proteome Discoverer converter

### PD-13: pg drops identified-but-unquantified channel rows; a test pins it {#pd-13}

*pg* · diverges from spec · partly confirmed (see correction) · checklist G1, D6

- **Now:** pg_adapter.py:113-114 skips every file-group/channel cell whose abundance is empty or ≤0. A test pins a single row when 127N is empty (tests/converters/test_proteome_discoverer_converter.py:197-210); verified by run.
- **Spec:** pg.md:120-135: in a unit with identification evidence, write one row per applicable label, with null intensity when unquantified. pg.md:143-145: zero sentinels become null.
- **Correction:** The code and test are as described. But converter-coverage.md:325 itself documents 'One row per quantified file group and channel', so the spec contradicts itself (pg.md:120-135 vs converter-coverage.md:325). The status is contradictory, not just a code divergence.
- **Reader impact:** Channel completeness and the type of missingness cannot be inferred. A plex looks unidentified in channels with blank abundance.
- **Suggested:** Write null-intensity rows for every channel of each file group where the master has PSM evidence (PSMs.txt by File ID).
- **Evidence:** `qpx/converters/proteome_discoverer/pg_adapter.py:113-114`; `tests/converters/test_proteome_discoverer_converter.py:197-210`; `docs/spec/pg.md:120-135`

### PD-14: Identification-only pg: one 'LFQ' row spanning all runs, sorted lexicographically {#pd-14}

*pg* · diverges from spec · partly confirmed (see correction) · checklist G1

- **Now:** Without abundance columns, pg_adapter.py:99-109 writes one row per master with label 'LFQ', intensity null, and grouped_runs = the sorted runs of all file groups (line 101). Rows are not split per unit, and evidence per run is not checked. A test pins label LFQ (tests/converters/test_proteome_discoverer_converter.py:213-227).
- **Spec:** pg.md:109,146: label is null for identification-only input without label context. pg.md:126-141: evidence is per unit; 'Missing evidence fields must not be treated as proof that a group was identified in every run'. pg.md:215: grouped_runs in fraction order, never sorted.
- **Correction:** The behaviour is confirmed. But converter-coverage.md:325 documents 'one LFQ row per protein with a null intensity', which contradicts pg.md:109/146 (label null without label context). The fraction-order rule is pg.md:414, not :215.
- **Reader impact:** Readers see all runs as one LFQ quantification unit, and lexicographic order breaks fraction order (F1, F10, F2).
- **Suggested:** Write one row per file group with PSM evidence, label null (or the plex labels when known), in fraction order from File ID.
- **Evidence:** `qpx/converters/proteome_discoverer/pg_adapter.py:99-109`; `tests/converters/test_proteome_discoverer_converter.py:213-227`; `docs/spec/pg.md:109`; `docs/spec/pg.md:126-141`; `docs/spec/pg.md:215`

### PD-15: Experiment-wide 'Exp. q-value: Combined' written to run-level pg_qvalue {#pd-15}

*pg* · diverges from spec · partly confirmed (see correction) · checklist G3, D4

- **Now:** column_mappings.yaml:75 maps 'Exp. q-value: Combined' and its spellings to pg_qvalue, and pg_adapter.py:147-148 sets global_qvalue to None. The same value repeats on every file-group row. Run: pg_qvalue=0.0, global_qvalue=null.
- **Spec:** pg.md:92-93: global_qvalue is experiment level, pg_qvalue run level. MaxQuant's experiment Q-value goes to global_qvalue (converter-coverage.md:299).
- **Correction:** The finding holds, but it is not PD-specific. converter-coverage.md:300 has MaxQuant write its experiment-level Q-value to pg_qvalue as well, and FragPipe writes its group q-value only to pg_qvalue (:299). It is an inconsistency across producers, and PD is missing from that table.
- **Reader impact:** A filter on global_qvalue keeps every PD group (all null), and per-run filters apply an experiment-level FDR.
- **Suggested:** Map 'Exp. q-value: Combined' to global_qvalue. Leave pg_qvalue null unless PD reports per-file q-values.
- **Evidence:** `qpx/converters/column_mappings.yaml:75`; `qpx/converters/proteome_discoverer/pg_adapter.py:147-148`; `docs/spec/pg.md:92-93`; `docs/spec/converter-coverage.md:299`; `run: audit-run:pd/cases/rank2/out/pd.pg.parquet`

### PD-16: Experiment-wide peptide counts repeated per file group; peptides lists the master only {#pd-16}

*pg* · contradictory · partly confirmed (see correction) · checklist G1, X1

- **Now:** pg_adapter.py:130-131,151-154 write the experiment-wide '# Peptides' and '# Unique Peptides' as peptide_counts on every file-group/channel row, and peptides=[{master, count}] only. A blank '# Peptides' raises 'Column peptides ... required' and aborts the whole pg view (verified with a synthetic row).
- **Spec:** pg.md:81-83: counts 'in this file'. pg.md:159: peptides per individual protein of the group. intensities.md:206 tolerates experiment-wide counts as dataset scope.
- **Correction:** pg.md:81-83 'in this file' is ambiguous (file vs run), and intensities.md:206 allows dataset-scope counts, so this is not a hard contradiction. Confirmed: peptides lists the master only, while pg.md:159 asks for one entry per member. A blank '# Peptides' raises out of convert() and aborts the whole conversion, including the ontology, provenance and dataset outputs.
- **Reader impact:** A per-unit evidence filter (e.g. ≥2 peptides in this run) silently uses global counts for PD.
- **Suggested:** Contract: state the scope of peptide_counts (unit or dataset). PD either fills per file group from PSMs or marks dataset scope.
- **Evidence:** `qpx/converters/proteome_discoverer/pg_adapter.py:130-131`; `qpx/converters/proteome_discoverer/pg_adapter.py:151-154`; `docs/spec/pg.md:81-83`; `docs/spec/pg.md:159`; `docs/spec/intensities.md:206`; `run: audit-run:pd/cases/shared_two_masters (ValueError)`

### PD-17: pg intensity is either raw per-file or grouped (normalized/scaled) abundance, unflagged {#pd-17}

*pg* · diverges from spec · partly confirmed (see correction) · checklist G2, D7, X1

- **Now:** pg_adapter.py:74-76 uses the 'Abundance: Fx' columns, else 'Abundances (Grouped)', which converter-coverage.md:325 says are normalized and scaled per PD settings. Three PD settings are not read: S/N vs intensity for TMT, summed vs Top3 protein roll-up, and imputation. There is no quantification_method cv_param or is_imputed, and provenance parameters are None (converter.py:92).
- **Spec:** pg.md:54-60: where a converter knows the method, it writes cv_params quantification_method. intensities.md:19: raw values. intensities.md:213: provenance keys isotope_corrected and normalization. pg.md:172: is_imputed.
- **Correction:** pg.md:54-60 asks for quantification_method only 'where a converter knows the method', and intensities.md:213 says a missing key means the input does not say. So a missing key is allowed. The real conflict is that converter-coverage.md:325 documents the normalized/scaled Grouped value as the primary intensity, while pg.md:110 and intensities.md:19,224 reserve that field for raw values. The status is contradictory.
- **Reader impact:** Two PD exports, or PD and MaxQuant, put differently processed quantities in intensity, and merge cannot warn.
- **Suggested:** Write quantification_method (pd_abundance or pd_grouped_abundance) and provenance normalization=unknown. Contract: intensity is raw or its processing is declared.
- **Evidence:** `qpx/converters/proteome_discoverer/pg_adapter.py:74-76`; `qpx/converters/proteome_discoverer/converter.py:92`; `docs/spec/converter-coverage.md:325`; `docs/spec/pg.md:54-60`; `docs/spec/intensities.md:213`

### PD-22: Group members = union of every Protein Group ID on the master row {#pd-22}

*pg* · ambiguous · confirmed · checklist F9, G1

- **Now:** pg_adapter.py:37-45 collects proteins per 'Protein Group IDs' token, and a master's pg_accessions is the union over every ID on its row. If a master also belongs to another PD group, that group's members join it. Not verified on real PD data. The 'Master' column name is hard-coded (pg_adapter.py:35,88).
- **Spec:** pg.md:33: pg_accessions are all proteins within this group. pg.md:67: groups may overlap.
- **Reader impact:** Groups may be over-inclusive compared with PD's. Feature and pg stay mutually consistent because both use the same function.
- **Suggested:** Take only the group that this master leads, and document the rule.
- **Evidence:** `qpx/converters/proteome_discoverer/pg_adapter.py:35-45`; `qpx/converters/proteome_discoverer/pg_adapter.py:88`; `docs/spec/pg.md:33`; `docs/spec/pg.md:67`

### PD-28: MS1-label protein abundances become duplicate 'LFQ' pg rows when PSMs lack label columns {#pd-28}

*pg* · diverges from spec · partly confirmed (see correction) · checklist X1, F7

- **Now:** The _ABUNDANCE regex (pg_adapter.py:19-21) matches 'Abundance: F1: Light, Sample' with no channel, which becomes label 'LFQ' (pg_adapter.py:115). The SILAC guard exists only in the feature adapter (feature_adapter.py:55-62). Run: Light=100 and Heavy=50 gave two rows with the same pg_id, and validation reported run_double_count.
- **Spec:** converter-coverage.md:319: 'Labelled MS1 exports stop the conversion'. pg.md:217: run-disjointness. intensities.md:84-88: MS1 labels are separate labels.
- **Correction:** The behaviour is confirmed. Run-disjointness is at pg.md:416, not :217. It occurs only when PSMs.txt has no 'Abundance' label columns, and whether real PD SILAC exports look like that is unverified.
- **Reader impact:** Duplicate keys, double-counted protein intensity, and light/heavy cannot be told apart.
- **Suggested:** Reject non-reporter channel tokens in pg columns too, or map Light/Heavy to SILAC labels.
- **Evidence:** `qpx/converters/proteome_discoverer/pg_adapter.py:19-21`; `qpx/converters/proteome_discoverer/pg_adapter.py:115`; `qpx/converters/proteome_discoverer/feature_adapter.py:55-62`; `docs/spec/converter-coverage.md:319`; `run: audit-run:pd/cases/silac_pg`

## mzIdentML and CDAP converters

### MZID_CDAP-06: CDAP pg can merge decoy features into a target protein group {#mzid_cdap-06}

*pg* · diverges from spec · confirmed · checklist G3, F8, P4

- **Now:** pg key = sorted stripped pg_accessions (cdap/pg_adapter.py:88-114); a decoy feature whose stripped membership equals a target feature's lands in the same row; is_decoy=bool_and (pg_adapter.py:103-104,152) -> false; reporter intensities summed (pg_adapter.py:187-207). No collision in the fixture; derived from code.
- **Spec:** pg.md:99-103 decoy groups are their own rows with is_decoy=true; pg.md:46-49 decoys never share accessions.
- **Reader impact:** Target protein quantities can include decoy-PSM reporter signal and the decoy group vanishes from protein-level FDR.
- **Suggested:** Fix via MZID_CDAP-05 (keep markers) and include is_decoy in any aggregation key.
- **Evidence:** `qpx/converters/cdap/pg_adapter.py:88-114,152,187-207`; `qpx/converters/cdap/base_adapter.py:140-146`

### MZID_CDAP-07: mzIdentML pg is_decoy parsed from accession prefix, PSM view uses isDecoy {#mzid_cdap-07}

*pg* · inconsistent across producers · confirmed · checklist G3, P4

- **Now:** pg.is_decoy = any accession starting with 'rev_','DECOY_','decoy_' (mzidentml/pg_adapter.py:26,231); ignores PeptideEvidence isDecoy that the PSM view uses (psm_records.py:138); misses REV_, XXX_, ###REV###. Raised in #249, unchanged on dev.
- **Spec:** pg.md:46-49 status is stated by is_decoy, 'do not parse it from the prefix'.
- **Reader impact:** A group of decoy PSMs can be is_decoy=false in pg (or vice versa) within the same dataset.
- **Suggested:** Derive pg.is_decoy from member DBSequence evidence (decoy when all members are decoy).
- **Evidence:** `qpx/converters/mzidentml/pg_adapter.py:26,231`; `qpx/converters/mzidentml/psm_records.py:138`; `issue #249`

### MZID_CDAP-25: CDAP pg groups are mapping sets, not inferred groups; overlapping, anchor-only peptides {#mzid_cdap-25}

*pg* · diverges from spec · partly confirmed (see correction) · checklist G1, G2, F9, G3

- **Now:** pg row = distinct set of mapped accessions per unit, no inference (cdap/pg_adapter.py:88-114); a protein sits in several groups (fixture: NP_958780.1 etc. in 2 of 81 groups); intensity = sum of feature reporter sums with exactly that set (187-207); peptides lists only the anchor (144-146); pg_names = accessions (156); contaminant always false (167).
- **Spec:** pg.md:64 row is a protein-group analyte; pg.md:159-166 peptides per protein in group; pg.md:95 contaminant bool, null.
- **Correction:** Overlapping groups are explicitly allowed (pg.md:67: groups are not disjoint, consumers MUST NOT assume a protein is in one group), and converter-coverage.md:303 documents contaminant 'Always false'. The divergences that hold: peptides lists only the anchor (pg_adapter.py:144-146 vs pg.md:159), and contaminant is hard-coded False although CDAP cells carry contaminant entries (see M3). Mapping-set grouping is unspecified rather than divergent.
- **Reader impact:** Shared peptides split into overlapping pseudo-groups; summing per protein double-counts; contaminant=false asserts something unknown.
- **Suggested:** Contract: pg declares its grouping method (inference vs mapping_set) and unknown contaminant is null.
- **Evidence:** `qpx/converters/cdap/pg_adapter.py:88-114,144-167,187-207`; `audit-run:mzid_cdap/cdap/cdap.pg.parquet overlap query`; `docs/spec/pg.md:64,95,159-166`

### MZID_CDAP-26: CDAP pg additional_scores carry PSM- and peptide-level q-values {#mzid_cdap-26}

*pg* · ambiguous · confirmed · checklist G3, D4

- **Now:** pg.additional_scores = min/max over features of the best-PSM scores: psm_level_qvalue, peptide_level_qvalue, msgf_specevalue, msgf_rawscore, msgf_denovo_score (cdap/pg_adapter.py:148-150,227-262; observed). global_qvalue and pg_qvalue are null.
- **Spec:** scores.md:85-87 pg scores apply to the whole group; no rule that a score name must state its level.
- **Reader impact:** The best PSM q-value on a pg row looks like protein confidence; filtering groups on it is not protein-level FDR.
- **Suggested:** Score names must carry their level; pg rows carry no lower-level scores, or only with a best_psm_ prefix.
- **Evidence:** `qpx/converters/cdap/pg_adapter.py:148-150,227-262`; `audit-run:mzid_cdap/cdap/cdap.pg.parquet additional_scores`; `docs/spec/scores.md:85-87`

### MZID_CDAP-27: mzIdentML pg ignores PDH passThreshold, uses dataset-wide runs, inflated peptide counts {#mzid_cdap-27}

*pg* · diverges from spec · partly confirmed (see correction) · checklist G1, G3, X1

- **Now:** Every ProteinAmbiguityGroup written; PDH passThreshold never read and no pass flag in pg (mzidentml/pg_adapter.py:129-143,214-254); grouped_runs = all runs feeding the list (78-127); peptides.peptide_count = number of PeptideHypothesis = PSM evidences (P0A6Y8: 2,930 in PXD054720; 187-190); every pg score higher_better=None (169,208).
- **Spec:** pg.md:125-136 evidence = PSM passing FDR in the unit; pg.md:166 'number of peptides'; scores.md:21 take direction from term/name.
- **Correction:** peptide_count counts PeptideHypothesis, i.e. PeptideEvidence refs. It equals the PSM count in PXD054720 only because Mascot writes one PeptideEvidence per query (11,332 PH = 11,332 SII refs); a producer with shared PeptideEvidences yields distinct-peptide counts. The PAG-level 'protein group passes threshold' (MS:1002415) is also dropped, and 'cluster identifier' (MS:1002407, an id) is stored as a score. All 148 PDHs pass in the fixture, so the sub-threshold claim is derived from code.
- **Reader impact:** Sub-threshold protein groups look identified; peptide_count is a PSM count; group identification asserted for every run.
- **Suggested:** Add a pg pass flag (or drop failing PDHs), count distinct peptides, set higher_better as for PSM scores.
- **Evidence:** `qpx/converters/mzidentml/pg_adapter.py:78-143,160-212`; `audit-run:mzid_cdap/pxd/PXD054720.pg.parquet peptides`; `docs/spec/pg.md:125-136,159-166`

### MZID_CDAP-M4: CDAP pg intensity is a reporter-ion sum with no quantification_method {#mzid_cdap-m4}

*pg* · diverges from spec · found by the verifier · checklist G2, F6, D7

- **Now:** pg.intensity sums each channel over the member features (pg_adapter.py:187-207), and each feature is a sum over its PSMs (feature_adapter.py:377-379). So the value is the sum of all PSM reporter ions whose protein set equals the group. pg.cv_params=None (pg_adapter.py:180): no quantification_method entry.
- **Spec:** pg.md:54-60: where a converter knows the method, cv_params carries quantification_method, and intensities are compared across datasets only when it agrees.
- **Reader impact:** CDAP protein abundances look like MaxLFQ or PD-normalized values; cross-dataset merges get no method warning.
- **Suggested:** Every pg producer that knows its method writes quantification_method (e.g. reporter_sum_all_psms); the contract lists allowed values.
- **Evidence:** `qpx/converters/cdap/pg_adapter.py:180,187-207`; `qpx/converters/cdap/feature_adapter.py:377-379`; `docs/spec/pg.md:54-60`

### MZID_CDAP-M5: CDAP peptide_counts.total_sequences counts features, not sequences {#mzid_cdap-m5}

*pg* · diverges from spec · found by the verifier · checklist G2, X1

- **Now:** total_sequences = len(list(sequence)) over features (pg_adapter.py:106,108,141-143,171), so it equals feature_counts.total_features in 81/81 groups. Fixture NP_004371.2: 1 distinct sequence, total_sequences=3. 12/81 groups have total > unique although each holds one sequence.
- **Spec:** pg.md:83: total_sequences = 'Total number of peptide sequences identified for this protein group in this file'.
- **Reader impact:** Peptide-count filters (e.g. >=2 peptides) keep single-peptide CDAP groups seen in several fractions or charges.
- **Suggested:** total_sequences = count of distinct sequences; the contract defines peptide_counts identically for every producer.
- **Evidence:** `qpx/converters/cdap/pg_adapter.py:106-108,141-143,169-172`; `verify-mzid_cdap/cdap/cdap.pg.parquet peptide_counts vs feature_counts`; `docs/spec/pg.md:81-86`

## Specification text

### SPEC-15: pg evidence = 'PSM that passed the FDR threshold', but no threshold is recorded {#spec-15}

*pg* · unspecified · confirmed · checklist G1, P3, D8

- **Now:** No column or footer key records the threshold or its level. Provenance has only free-form keys (psm_fdr etc. in the example, provenance.md:160-162; qvalue_threshold option, :65). The MaxQuant run provenance has no FDR key (run). DIA-NN never uses PG q-values as a cutoff (pg.md:326). FragPipe input is pre-filtered (pg.md:365,406).
- **Spec:** pg.md:123,126-129: a group has evidence when a PSM 'that passed the FDR threshold' maps to it. The threshold is not defined.
- **Reader impact:** A null-intensity pg row may rest on 1% PSM FDR, on protein FDR, or on no filter at all, depending on the producer.
- **Suggested:** Define normative provenance keys (psm_fdr, precursor_fdr, protein_fdr, applied_to_view) or a per-view footer filter key.
- **Evidence:** `docs/spec/pg.md:123,126-129,326,365,406`; `docs/spec/provenance.md:65,155-163,248-252`; `audit-run:spec/mq/maxquant.provenance.parquet`

### SPEC-24: pg_qvalue (run-level) holds experiment-level or 1-probability values {#spec-24}

*pg* · contradictory · confirmed · checklist G3, D4

- **Now:** MaxQuant writes the experiment-level Q-value into both global_qvalue and pg_qvalue. FragPipe writes the group q-value, or 1 - Protein Probability, into pg_qvalue (converter-coverage.md:299-300). OpenMS global_qvalue = the lowest member-protein q (:104).
- **Spec:** pg.yaml:52-61: global_qvalue = experiment level, pg_qvalue = run level. pg.md:406: for FragPipe, global_qvalue null and the probability goes to additional_scores. feature.md:132: q-value columns take only stated q-values.
- **Reader impact:** A run-level protein FDR filter actually applies experiment-level or probability-derived values.
- **Suggested:** Enforce the level per column, forbid 1-probability as a q-value, and leave MaxQuant pg_qvalue null.
- **Evidence:** `qpx/core/data/schemas/pg.yaml:52-61`; `docs/spec/pg.md:93,320,406`; `docs/spec/converter-coverage.md:104,299-300`; `docs/spec/feature.md:132`

### SPEC-25: pg rows from MBR only, identification-only labels and unit granularity are unflagged {#spec-25}

*pg* · inconsistent across producers · confirmed · checklist G1, F4

- **Now:** PD identification-only export: one LFQ row with null intensity (converter-coverage.md:325); mzIdentML: null label (pg.md:146). MaxQuant 'Identification type' only reaches cv_params (pg.md:291). views.md:214 says CDAP and Spectronaut write per-raw-file units. The code groups fractions for CDAP (cdap/pg_adapter.py:161); Spectronaut is per run (spectronaut/pg_adapter.py:244).
- **Spec:** pg.md:126-135 allows rows without identification evidence when an MBR quantity exists, with no typed flag. pg.yaml:65-68: label null only without label context. data-model.md:28-30: grouped_runs = the fractions of one sample.
- **Reader impact:** Transfer-only protein quantities look identified, and the quantification unit varies by producer, which breaks rollups.
- **Suggested:** Add a typed per-row pg identification flag, one label rule for identification-only input, and require grouped_runs to hold all fractions. Fix views.md.
- **Evidence:** `docs/spec/pg.md:126-146,291`; `docs/spec/converter-coverage.md:325`; `docs/spec/views.md:214`; `qpx/converters/cdap/pg_adapter.py:161`; `qpx/converters/spectronaut/pg_adapter.py:244`; `docs/spec/data-model.md:28-30`

### SPEC-26: pg counts scope: 'within this file' vs unit vs experiment-wide {#spec-26}

*pg* · ambiguous · partly confirmed (see correction) · checklist G1, X1

- **Now:** MaxQuant maps 'Unique peptides [exp]'/'Peptides [exp]' (pg.md:297-298), which are per experiment. intensities.md:206 says some tools report experiment-wide counts. FragPipe uses per-experiment counts.
- **Spec:** pg.md:82-86: counts 'within this file'. pg.md:12,411: 'per-file protein group quantification', although pg.md:3 says rows are per quantification unit.
- **Correction:** MaxQuant and FragPipe [exp] counts are per experiment, which is the row's unit, so they fit the row (intensities.md:206 says so). The ambiguity is the wording 'within this file' (pg.md:81-86) and intensities.md:206 allowing dataset-wide counts in the same columns with no marker.
- **Reader impact:** Counts may describe the row's unit, the experiment or the dataset, so filters like 'two peptides per unit' are unreliable.
- **Suggested:** Counts refer to the row's grouped_runs. Experiment-wide counts go under named additional_scores.
- **Evidence:** `docs/spec/pg.md:3,12,82-86,297-298,411`; `docs/spec/intensities.md:206`

### SPEC-M3: Decoy flag of groups with target and decoy members differs by producer {#spec-m3}

*pg* · inconsistent across producers · found by the verifier · checklist G3, P4

- **Now:** Any member decoy makes the group decoy in Spectronaut (pg_adapter.py:234) and mzIdentML (pg_adapter.py:231). All members must be decoys in DIA-NN (:238), OpenMS (:394), MSstats (feature_adapter.py:82) and CDAP features (:252). FragPipe uses the leader only (fragpipe/pg_adapter.py:442).
- **Spec:** psm.md:45 and feature.md:60: a peptide matching target and decoy proteins is a target. pg.md:94,99-103 gives no rule for mixed groups.
- **Reader impact:** targets_only() and protein-level FDR re-estimation keep or drop different groups for the same data.
- **Suggested:** Normative rule: a group is a decoy only if all its members are decoys, mirroring the peptide rule.
- **Evidence:** `qpx/converters/spectronaut/pg_adapter.py:234`; `qpx/converters/mzidentml/pg_adapter.py:231`; `qpx/converters/diann/pg_adapter.py:238`; `qpx/converters/openms_consensus/pg_adapter.py:394`; `qpx/converters/fragpipe/pg_adapter.py:442`; `docs/spec/pg.md:94-103`

### SPEC-M4: pg.contaminant: false vs null undefined, often parsed from the prefix {#spec-m4}

*pg* · inconsistent across producers · found by the verifier · checklist G3, P4, X1

- **Now:** CDAP hard-codes False (cdap/pg_adapter.py:167; all 810 rows). Spectronaut and mzIdentML write None. FragPipe, DIA-NN and OpenMS parse the accession prefix (is_contaminant_accession). MaxQuant and PD read the tool's column.
- **Spec:** pg.md:95: 'Contaminant flag', bool, null, with no meaning given for null. ontology.md:241: '1 = contaminant'. pg.md:46-49: the flag states the status.
- **Reader impact:** false can mean 'not checked', so contaminant filtering misses CDAP contaminants, and null and false are not comparable.
- **Suggested:** null = unknown. false only when the tool or a declared contaminant list says so. Record that list or prefix in provenance.
- **Evidence:** `qpx/converters/cdap/pg_adapter.py:167`; `qpx/converters/spectronaut/pg_adapter.py:250`; `qpx/converters/mzidentml/pg_adapter.py:246`; `qpx/converters/fragpipe/pg_adapter.py:469`; `qpx/converters/diann/pg_adapter.py:239`; `qpx/converters/openms_consensus/pg_adapter.py:407`; `docs/spec/pg.md:95`; `docs/spec/ontology.md:241`

## Validator and read API

### VALIDATOR-25: anchor_protein is still used as the protein key in several read paths {#validator-25}

*pg* · diverges from spec · partly confirmed (see correction) · checklist F9, G2

- **Now:** intensity(), the views and the MuData var now key proteins by membership, which changed from main. These still key by anchor_protein: Feature.by_protein and unique_proteins (feature.py:33-41); for_quantification, which exports anchor_protein as the protein (feature.py:91-115); qpxc transform quantify (ProteinName=anchor_protein, var_key 'anchor_protein', quantify.py:145,261-267); and site_quant occurrences (site_quant.py:24-55). views.md:212 still says the protein feature_id is the anchor_protein.
- **Spec:** pg.yaml:43 and pg.md:66: anchor_protein MUST NOT be used for feature-to-pg joins or per-protein statistics.
- **Correction:** The code follows other parts of the spec: absolute.md:67,88,198 (quantify keys by anchor_protein), site.md:113 (site quant on anchor_protein) and views.md:177-178 (by_protein matches the anchor). So this is a spec-internal contradiction with pg.yaml:43 and pg.md:66, not code diverging from spec. views.md:212 is stale: intensity() keys protein by membership (dataset_queries.py:90).
- **Reader impact:** Absolute quantification and exported tables merge distinct groups that share a leader, which double-counts shared signal.
- **Suggested:** Key every export and transform on pg_accessions or pg_group_id, keep anchor as a display column only, and fix views.md:212,214.
- **Evidence:** `qpx/core/data/feature.py:33-41,91-115`; `qpx/cli/transform/quantify.py:144-145,261-269`; `qpx/transforms/site_quant.py:24-55`; `qpx/core/dataset_queries.py:88-98`; `docs/spec/views.md:212,214`; `qpx/core/data/schemas/pg.yaml:43`

### VALIDATOR-28: A pg row with a null label can still carry an intensity {#validator-28}

*pg* · unspecified · confirmed · checklist G1

- **Now:** A pg row with label null and intensity 5.0 passes --conformance. The sample-mapping SQL skips null labels (label_validation.py:51-52). proteins_per_sample counts null-label rows of single-sample runs (views/api.py:262), while ProteinView drops them through the exact label join.
- **Spec:** pg.yaml:63-68: label and intensity are both null only when label context is absent. pg.md:117-124 is the evidence table.
- **Reader impact:** An unattributable quantity can enter some counts and not others.
- **Suggested:** Validate that a null label implies a null intensity.
- **Evidence:** `experiment mutate.py g_null_label_row_with_intensity -> exit 0`; `qpx/core/views/api.py:262 (qpx/views/api.py:262)`; `qpx/core/label_validation.py:51-52`

## mokume (consumer)

### MOKUME-24: Protein FDR uses the optional feature.pg_global_qvalue; pg.parquet is never read {#mokume-24}

*pg* · inconsistent across producers · confirmed · checklist G3, D4

- **Now:** `--filter-protein-fdr` uses the minimum feature.pg_global_qvalue per group (filters/protein.py:193-228; lib.rs:6699-6721) and errors when the column is empty. mokume reads no pg.parquet: no intensity, global_qvalue or contaminant from it. MaxQuant fills the column only with `--protein-groups-file` (converter-coverage.md:196). The FragPipe, PD and CDAP dev outputs leave it null.
- **Spec:** feature.yaml:167-172: optional, the experiment-level q-value of the group. pg holds group q-values per unit and label. How the feature value relates to the pg rows is not stated.
- **Reader impact:** Protein FDR is unavailable for many producers. Feature-only and pg-based consumers may apply different protein FDR.
- **Suggested:** Either require feature.pg_global_qvalue whenever pg has a value, or define the feature to pg join that consumers must use for protein FDR.
- **Evidence:** `mokume/python/mokume/preprocessing/filters/protein.py:193-228`; `mokume/rust/crates/mokume-pipeline/src/lib.rs:6699-6721`; `qpx/core/data/schemas/feature.yaml:167-172`; `docs/spec/converter-coverage.md:196`

## External tools and standards

### EXTERNAL-M5: PD pg primary abundances may be normalized and scaled, but no normalization key is written {#external-m5}

*pg* · diverges from spec · found by the verifier · checklist G2, F6

- **Now:** PD 3.1 has a Normalization Mode (e.g. Total Peptide Amount) and a Scaling Mode 'On All Average' (each protein's average across samples = 100). QPX PD pg intensity is 'Abundance: F1: 126, Sample' or 'Abundances (Grouped)', and that value 'is normalized and scaled if the consensus step was set so' (converter-coverage.md:325). The provenance parameters are None (proteome_discoverer/converter.py:92).
- **Spec:** Primary intensities are raw (intensities.md:19-22). The normalization key states processing, and 'No key means the input does not say' (intensities.md:213).
- **Reader impact:** PD protein intensities may be normalized or scaled to 100, yet they are stored as primary 'raw' values with no flag. Re-normalizing, or comparing with MaxQuant/OpenMS raw sums, is wrong.
- **Suggested:** A converter that cannot tell must write normalization='unknown' instead of omitting the key. PD reads the consensus-step settings when exported. Scaled values never go into primary intensities.
- **Evidence:** `https://docs.thermofisher.com/r/Proteome-Discoverer-3.1-User-Guide/en-US1454969867v1 (Normalization Mode)`; `docs/spec/converter-coverage.md:325`; `qpx/converters/proteome_discoverer/converter.py:92`; `qpx/converters/proteome_discoverer/pg_adapter.py:20-23`; `docs/spec/intensities.md:19-22,213`

### EXTERNAL-M6: anchor_protein mixes three distinct PSI concepts (anchor, leading, group representative) {#external-m6}

*pg* · ambiguous · found by the verifier · checklist F9, P8

- **Now:** PSI-MS separates three terms. MS:1001591 anchor protein is the 'representative ... from a set of sequence same-set ... proteins'. MS:1002401 leading protein marks at least one protein per group with the strongest evidence. MS:1002403 group representative is 'exactly one protein per group', a tiebreaker among the leading proteins.
- **Spec:** pg.md:39 defines anchor_protein as 'Representative protein of the group (leading protein)', while ontology.md:234 maps it to MS:1001591 'anchor protein'.
- **Reader impact:** Converters may fill anchor_protein with a same-set anchor, one of several leading proteins, or a group representative. Displays and accession-level joins then differ by producer.
- **Suggested:** Map anchor_protein to MS:1002403 (exactly one per group). Optionally add a leading_proteins list. Each converter documents which tool column fills it.
- **Evidence:** `https://www.ebi.ac.uk/ols4/api/ontologies/ms/terms?obo_id=MS:1001591`; `https://www.ebi.ac.uk/ols4/api/ontologies/ms/terms?obo_id=MS:1002401`; `https://www.ebi.ac.uk/ols4/api/ontologies/ms/terms?obo_id=MS:1002403`; `docs/spec/pg.md:39`; `docs/spec/ontology.md:234`

## First-principles questions

### CRITIC-48: Does a pg row aggregate all fractions of a sample, or one raw file? {#critic-48}

*pg* · contradictory · confirmed · checklist G1, X1

- **Now:** pg.md and data-model say grouped_runs = the fractions aggregated into one unit. views.md says the CDAP and Spectronaut converters write one pg row per raw file. Yet the CDAP fixture output groups f01 and f02 into one row.
- **Spec:** pg.md:3,40,414-416; data-model.md:29-30; views.md:214
- **Reader impact:** Protein matrices sum or average per-fraction partials, and the docs disagree with the output.
- **Suggested:** Require grouped_runs = all fractions of the unit, forbid per-fraction protein values in pg, and fix views.md.
- **Evidence:** `docs/spec/pg.md:3`; `docs/spec/views.md:214`; `audit-run:critic/cdap/cdap.pg.parquet`

### CRITIC-49: Which pg rows exist, and can a reader tell identified rows from MBR-only rows? {#critic-49}

*pg* · ambiguous · confirmed · checklist G1, F4

- **Now:** A labelled row exists when at least one PSM passed 'the FDR threshold' in the unit (null intensity allowed); MBR-only rows exist only with a quantity. The threshold is unspecified, and no flag tells evidence-backed rows from transfer-only ones (MaxQuant 'Identification type' goes only into cv_params).
- **Spec:** pg.md:116-141,291
- **Reader impact:** Counts of identified proteins per sample include MBR-only rows, though andes#218 says transfers never create evidence.
- **Suggested:** Add a bool identified_in_unit, and reference the footer filter keys for the threshold.
- **Evidence:** `docs/spec/pg.md:116-141`; `docs/spec/pg.md:291`

### CRITIC-50: What does pg.intensity measure (sum, TopN, MaxLFQ; razor or unique peptides)? {#critic-50}

*pg* · ambiguous · confirmed · checklist G2

- **Now:** pg.md calls it the 'primary raw intensity' but says it differs between tools; quantification_method is an optional cv. Producers write MaxQuant Intensity, DIA-NN PG.Quantity or MaxLFQ, FragPipe razor intensity, and OpenMS abundance or a unique-peptide sum.
- **Spec:** pg.md:52-60,110,306; converter-coverage.md:107,301
- **Reader impact:** Protein abundances from two producers are not comparable, and merge only warns.
- **Suggested:** A mandatory typed quantification_method column with a controlled vocabulary.
- **Evidence:** `docs/spec/pg.md:52-60`; `docs/spec/converter-coverage.md:107`; `docs/spec/converter-coverage.md:301`

### CRITIC-51: Is pg_qvalue run-level and global_qvalue experiment-level for every producer? {#critic-51}

*pg* · contradictory · confirmed · checklist G3, D4

- **Now:** The spec says global = experiment and pg_qvalue = run. MaxQuant writes the same experiment Q-value into both. FragPipe writes the group q-value into pg_qvalue, or 1 − Protein Probability, against the rule that only stated q-values fill q columns.
- **Spec:** pg.md:92-93,406; feature.md:132; converter-coverage.md:299-300
- **Reader impact:** A run-level protein FDR filter actually applies an experiment-level FDR or a probability.
- **Suggested:** Enforce the level semantics; probabilities go only into additional_scores.
- **Evidence:** `docs/spec/pg.md:92-93`; `docs/spec/pg.md:406`; `docs/spec/converter-coverage.md:299-300`

### CRITIC-52: What defines a protein group's membership, and is it fixed across units? {#critic-52}

*pg* · ambiguous · partly confirmed (see correction) · checklist G1, F9

- **Now:** 'All proteins within this group'; groups may overlap. Producers differ: MaxQuant Protein IDs (includes subset proteins), FragPipe lead + indistinguishable, DIA-NN Protein.Group. #51's 'groups are global, same for all runs' is not in the spec.
- **Spec:** pg.md:33,62-69,282-283; tool-mappings.md:102
- **Correction:** intensities.md:207 states that the group is computed over the dataset, which implies dataset-global membership, but there is no MUST rule. The membership definitions do differ by producer.
- **Reader impact:** The same evidence yields different groups per producer, and per-unit regrouping breaks cross-sample matrices.
- **Suggested:** Membership = the indistinguishable set, plus an optional subsumed list. Groups are dataset-global.
- **Evidence:** `docs/spec/pg.md:33`; `docs/spec/pg.md:67`; `docs/spec/pg.md:282-283`

### CRITIC-53: How are is_decoy and contaminant set for groups with mixed members? {#critic-53}

*pg* · inconsistent across producers · partly confirmed (see correction) · checklist G3

- **Now:** No rule. DIA-NN: decoy only when every accession is prefixed; contaminant when any member is CONTAM. CDAP: contaminant always false. FragPipe: not implemented (null). MaxQuant: from Reverse and Potential contaminant.
- **Spec:** pg.md:94-95,336; converter-coverage.md:103,303
- **Correction:** FragPipe does set contaminant (and is_decoy) from the lead accession marker (fragpipe/pg_adapter.py:441,469), not null; converter-coverage.md:303 is stale. CDAP's constant false is reproduced (810/810 rows). The cross-producer inconsistency holds.
- **Reader impact:** Contaminant filters remove different groups per producer, and CDAP's false is an unsupported claim.
- **Suggested:** Contaminant = any member; decoy rule fixed by the spec. Write null when not determined, never a constant false.
- **Evidence:** `docs/spec/pg.md:336`; `docs/spec/converter-coverage.md:103`; `docs/spec/converter-coverage.md:303`

### CRITIC-54: What scope do peptide_counts, feature_counts and peptides have? {#critic-54}

*pg* · ambiguous · confirmed · checklist G1

- **Now:** peptide_counts is described as 'in this file', though rows are per unit. MaxQuant counts are per experiment. `peptides` takes MaxQuant 'Peptide counts (all)', a dataset-wide count repeated on every row.
- **Spec:** pg.md:79-86,155-166,297-299; intensities.md:204-206
- **Reader impact:** A filter on '>=2 peptides in this sample' actually uses dataset-wide counts.
- **Suggested:** Counts are per the row's unit and label; dataset-wide counts get separate names.
- **Evidence:** `docs/spec/pg.md:79-86`; `docs/spec/pg.md:299`; `docs/spec/intensities.md:204-206`

### CRITIC-66: In the API views, is a protein analyte its group or its anchor_protein? {#critic-66}

*pg* · contradictory · confirmed · checklist G1

- **Now:** views.md first says 'A protein is its group: feature_id is the sorted membership', then 'At protein level... feature_id is the anchor_protein'. pg.md forbids using the anchor as the analyte.
- **Spec:** views.md:153 vs views.md:212,214; pg.md:66
- **Reader impact:** Matrices may merge distinct groups that share a leading protein.
- **Suggested:** Use the membership key everywhere and fix views.md:212.
- **Evidence:** `docs/spec/views.md:153`; `docs/spec/views.md:212`; `docs/spec/pg.md:66`
