# Converter Coverage Matrix

This page lists which QPX data views each converter produces and explains the sources and null values of selected fields.

## Views produced per converter

| Converter | PSM | Feature | PG | Pepmap | Sample | Run | Dataset | Ontology | Provenance | mz |
| ----------- | :---: | :-------: | :--: | :------: | :------: | :---: | :-------: | :--------: | :----------: | :--: |
| **MaxQuant** | Yes | Yes | Yes | If `Proteins` | If SDRF | If SDRF | Yes | If terms | Yes | No |
| **FragPipe** | Yes | Yes | Yes | Yes | If SDRF | If SDRF | Yes | If terms | Yes | No |
| **Proteome Discoverer** | Yes | If abundances | If Proteins.txt | If abundances | If SDRF | If SDRF | Yes | If terms | Yes | No |
| **DIA-NN** | No | Yes | Yes | Yes | Yes | Yes | Yes | If terms | Yes | If MS2 info |
| **Spectronaut** | No | Yes | Yes | No | If SDRF | If SDRF | Yes | Yes | Yes | No |
| **OpenMS native QPX** | Yes | Yes | Yes | No | Yes | Yes | Yes | Yes | Yes | No |
| **OpenMS consensusXML** | Yes | Yes | Yes | No | If SDRF | If SDRF | Yes | If terms | Yes | No |
| **CDAP** | Yes | Yes | Yes | No | If PDC | If PDC | Yes | Yes | Yes | No |
| **mzIdentML** | Yes | No | If in file | Yes | If SDRF | If SDRF | Yes | Yes | Yes | No |
| **QuantMS MSstats** | No | Yes | No | No | Yes | Yes | Yes | Yes | Yes | No |
| **SDRF** | No | No | No | No | Yes | Yes | No | Optional | No | No |

- **Yes** — the converter produces this view.
- **No** — the converter does not produce this view (e.g. DIA-NN has no PSM view; mzIdentML has no Feature).
- **If SDRF** — the view is produced only when an SDRF file is provided. Whether SDRF is optional or required depends on the converter; see its CLI contract below.
- **If `Proteins`** — the view is produced when `evidence.txt` has a `Proteins` column.
- **If terms** — the view is produced when the converter discovers resolvable ontology entries.
- **If in file** — the view is produced only when the input contains the source section (mzIdentML: a `ProteinDetectionList`).
- **If MS2 info** — the view is produced only when `--mzml-info-folder` contains `*_ms2_info.parquet` files (quantms-utils mzML statistics run with `--ms2_file`) and the generated feature view can be read. It holds the MS2 spectra matched to features; if no spectrum matches, no file is written.
- **If PDC** — CPTAC/PDC studies ship no SDRF, and CDAP `.psm` files carry no sample metadata. When run through `qpxc pdc2qpx` (default), the sample/run views are built from PDC GraphQL metadata, which also recovers the TMT/iTRAQ channel → biological-sample mapping. Disable with `--no-metadata`.

> The `mz` (full-spectra) view is produced by the standalone `qpxc convert mz` command, or automatically by `qpxc pdc2qpx --include-spectra`. Of the per-tool converters above, only DIA-NN emits it, and only with the feature-matched MS2 spectra (see **If MS2 info**).

## Field population: DIA-NN and OpenMS consensusXML

These tables explain selected fields whose population differs between DIA-NN
and OpenMS LFQ/TMT conversions. The OpenMS column refers to
`qpxc convert openms-consensus`. Its in-memory and streaming readers use the
same field mappings.

Source columns and annotations must be present and usable for a value to be
written. An optional field can remain null because the source does not report
it, the evidence is ambiguous, or the converter does not perform the required
calculation. **Not populated** below describes current converter behavior; it
does not imply that the quantity could never be derived.

### Feature fields

| Field | DIA-NN | OpenMS consensusXML |
| ------- | -------- | --------------------- |
| `missed_cleavages` | Calculated from the peptide sequence and the SDRF enzyme. Null when the enzyme is unavailable or unsupported. | Calculated from the peptide sequence and the SDRF enzyme, falling back to `SearchParameters.enzyme`. Null when the enzyme is unavailable or unsupported. |
| `cv_params` | Tool-specific metrics such as `proteotypic`, `irt`, `predicted_irt`, `iim` and `diann_precursor_lib_index` (`Precursor.Lib.Index`), when reported; `Quantity.Quality` is only the score `diann_quantity_quality`. `iim`, `predicted_im` and `predicted_iim` are each left out where DIA-NN reports 0. For data without ion mobility, `Predicted.IM` and `Predicted.iIM` are 0; DIA-NN 1.8 also writes `iIM` = 0, while DIA-NN 2.x reports a library `iIM`, which is kept. | The `MS1Label:*` label state of an MS1-labelled multiplet. Otherwise null: the converter has no corresponding source metrics to put in this field. These are separate from `pg.cv_params`. |
| `pg_global_qvalue` | Read from `Global.PG.Q.Value`, when present. | The resolved protein group's global q-value, using the same derivation as `pg.global_qvalue`. Null when the group is ambiguous or has no declared q-value. |
| `gg_accessions` / `gg_names` | Gene symbols from `Genes`, when present. | Gene symbols parsed from `GN=` in protein-hit descriptions, through the resolved protein group. Null when that annotation is absent or the group cannot be resolved. The converter does not read a FASTA to add missing genes. |
| `ion_mobility` | `IM` (and `IM.Start` / `IM.Stop` for `ion_mobility_start` / `ion_mobility_stop`), when present. For data acquired without ion mobility, DIA-NN reports 0 for every precursor; a 0 is written as null. | Not populated. |
| `rt_start` / `rt_stop` | Read from `RT.Start` / `RT.Stop` and converted from minutes to seconds. Null when those columns or values are absent. | Null: consensusXML feature elements provide an apex RT, not the feature's integration-window boundaries. The apex is stored in `rt`. |
| `identification_type` | Null: a report row does not say whether its identification was transferred. | `direct` for a run with a consistent identification with a spectrum reference; `transferred` for a run without identifications when every identification of the feature has a known run; otherwise null. |
| `parent_ion_fraction` / `quant_scan` / `quant_ms_level` | Not provided. | `quant_precursor_purity`, else `precursor_purity`; `scan_id`; and 2 or 3 from whether `id_scan_id` equals `scan_id`, null without `id_scan_id`. |
| `id_run_file_name` | The normalized report run name. | The resolved run of a direct identification of the exported peptide with a spectrum reference. Null for transferred features without their own identification, unresolved origins or conflicting peptide assignments. |
| `pg_positions` | Not populated from the report, which has no peptide-to-protein coordinates. Filled from a FASTA with `--fasta` or `qpxc transform protein-properties`: every one-based occurrence in each group member. | Known peptide-evidence `start` / `end` coordinates, converted from zero-based inclusive to one-based inclusive. For a resolved protein group, positions on its members; for a peptide shared across inferred groups (null `pg_accessions`), positions on every evidence protein, since each entry names its own protein. Null when no usable positions remain. A FASTA is not needed for coordinates already recorded in consensusXML. |
| `peptide_qvalue` | Null: DIA-NN `Q.Value` is precursor-level and is written to the `precursor_qvalue` column. | Read from an identification score explicitly declared to be a q-value/FDR. Null without a suitable identification in that run; raw search scores are not substituted. |
| `intensities` / `additional_intensities` | `Precursor.Quantity`; `PG.MaxLFQ`, `Precursor.Normalised` and `Ms1.Area` as `pg_maxlfq`, `normalize_intensity` and `ms1_area`. DIA-NN writes 0 for a quantity it did not measure: an unquantified precursor or channel keeps its label with a null `intensity`, and a 0 additional intensity is left out. | The consensus feature's per-map sub-feature intensities; a map element with intensity ≤ 0 has no entry. `additional_intensities` holds the ratios of an MS1-labelled multiplet; otherwise it is not populated. |

Protein-group ambiguity includes shared peptides spanning multiple inferred
groups and conflicting assignments within a run. These features retain their
quantification while the ambiguous protein-group fields remain null. See the
[OpenMS conversion guide](../guide/convert.md#openms-consensus-description) for
the evidence-selection rules, and [Feature confidence](feature.md#protein-group-fields)
for q-value semantics.

Further DIA-NN feature fields:

- **`observed_mz` / `mass_error_ppm`.** DIA-NN's `Precursor.Mz` is the theoretical
  m/z. The measured one is `Precursor.Mz` + `Ms1.Apex.Mz.Delta` (DIA-NN 2.x), and
  `mass_error_ppm` compares it with `calculated_mz`. Both are null when the offset
  is `0`, which marks a precursor without a measured MS1 apex, and for reports
  without the offset (DIA-NN 1.8 and 1.9).
- **Zeros.** A `Precursor.Quantity` of `0` is an `intensities` entry with a null
  `intensity`; an `IM`, `IM.Start`, `IM.Stop`, `iIM` or predicted ion mobility of
  `0` (data acquired without ion mobility) is null or left out of `cv_params`.
- **`unique`.** True when `Protein.Group` and `Protein.Ids` both name one protein.
  `Proteotypic` follows DIA-NN's protein-grouping level (genes by default) and
  stays in `cv_params`.
- **Pepmap.** `feature.pg_accessions` holds the peptide's protein group
  (`Protein.Group`); the pepmap view lists every protein in `Protein.Ids` for each
  converted target peptidoform, without positions. Decoy rows, which carry a
  mutated sequence under their parent's proteins, are left out.
- **Library decoys.** DIA-NN 2.x prefixes a decoy searched from the spectral
  library (a first-pass decoy kept in an empirical library) with `^` in
  `Modified.Sequence` and `Precursor.Id`; its sequence is the decoy generated
  from the same parent. The `^` is left out of `peptidoform` and the row has
  `is_library_decoy` = true. A decoy's sequence can equal a target's, so
  `is_decoy` and `is_library_decoy` are both part of the feature identity. A generated decoy shares
  `diann_precursor_lib_index` with its target; a library decoy has its own.
- **plexDIA.** A report with a `Channel` column gives one row per multiplet, with
  one `intensities` entry per channel. With an SDRF, a channel matches the
  `comment[label]` of the same name or the plexDIA label DIA-NN calls by that
  channel name (`MTRAQ4` for `4`, `SILAC heavy` for `H`), from the sdrf-pipelines
  plexDIA registry. The peptidoform is the lightest partner's: mTRAQ and dimethyl
  tags become their Unimod accessions and the SILAC tag is dropped. `observed_mz`
  and `mass_error_ppm` come from the lightest channel, and the other scalar fields
  from the channel with the best `Q.Value`, ties broken by label.
- **Modification names.** With `--original-mods`, DIA-NN reports modifications by
  name. Names Unimod knows become their accessions; with `--diann-log`, a name
  declared by `--fixed-mod` or `--var-mod` takes its declared mass.

### Protein group fields

| Field | DIA-NN | OpenMS consensusXML |
| ------- | -------- | --------------------- |
| `pg_names` | Split from the report's `Protein.Names`. | UniProt entry names taken from `db\|ACCESSION\|NAME` accessions, aligned with `pg_accessions`. Null when any member is a bare accession with no entry name, so the list never misaligns; a bare accession is not repeated as its own name. |
| `contaminant` | True when any member accession carries the `CONTAM` marker of the quantms/sdrf-pipelines contaminant database. | Same rule, shared with DIA-NN (`qpx.converters.utils.is_contaminant_accession`). |
| `global_qvalue` | Read from `Global.PG.Q.Value`, when present. | The lowest available member-protein q-value, only when the source protein score is declared to be a q-value/FDR. Null if no member has such a score. |
| `pg_qvalue` | Read from run-level `PG.Q.Value`, when present. | Null: no separate run-level protein-group q-value is reported by this path. The experiment-level value belongs in `global_qvalue`. |
| `gg_qvalue` | Read from `GG.Q.Value`, when present. | Null: this path has no gene-group q-value. Protein-level confidence is not a gene-level q-value. |
| `label` / `intensity` | The first positive value of `PG.Quantity`, `PG.MaxLFQ` and, for label-free data, the PG matrix; `quantification_method` in `cv_params` names the source. DIA-NN writes 0 for a quantity it did not compute, as QuantUMS does in `PG.MaxLFQ`: a group (or, in channel mode, a channel) without a positive quantity keeps its `label` with `intensity` = null. | OpenMS' protein abundance when the input reports it, else the interim unique-peptide sum; `quantification_method` in `cv_params` names the source. Null for a group with a PSM in the unit but no quantity. |
| `additional_intensities` | `PG.MaxLFQ` as `maxlfq` beside `PG.Quantity`, and `PG.Normalised`, `PG.TopN`, `Genes.MaxLFQ`, `Genes.MaxLFQ.Unique` and `Genes.TopN` as `normalize_intensity`, `topn`, `genes_maxlfq`, `genes_maxlfq_unique` and `genes_topn`, when reported; a 0 is left out. | With OpenMS' protein abundances: the unique-peptide sum under its `quantification_method` name and the ratios OpenMS reports. Otherwise not populated. |
| `additional_scores` | Includes `PG.Q.Value` as `pg_qvalue`, `PG.PEP` as `pg_pep`, `Lib.PG.Q.Value` as `lib_pg_qvalue`, `Protein.Q.Value` as `protein_qvalue` and the QuantUMS `PG.MaxLFQ.Quality` / `Genes.MaxLFQ.Quality` as `pg_maxlfq_quality` / `genes_maxlfq_quality`, when available; in a plexDIA report each channel's row carries that channel's values. | The existing `anchor_protein`'s `Posterior Probability_score`, stored as `posterior_probability` with `higher_better: true`. Null if absent, invalid or conflicting for that anchor. |
| `sequence_coverage` | Not populated from the report. Filled from a FASTA with `--fasta` or `qpxc transform protein-properties`: percent of the anchor covered by the dataset's target peptides whose evidence names it. | The existing anchor's recorded `ProteinHit.coverage` in percent. Unknown values, missing annotations and conflicting coverage values for that anchor remain null. Other group members are not substitutes. |
| `molecular_weight` | Not populated from the report. Filled from a FASTA with `--fasta` or `qpxc transform protein-properties`: anchor average mass in kDa. | The theoretical average molecular weight of the existing anchor's complete, unmodified `ProteinHit.sequence`, in kDa. Null for missing, ambiguous, modified or conflicting sequences; other group members are not substitutes. |
| `cv_params` | Each quantified row records whether its primary intensity came from `PG.Quantity` or the `PG.MaxLFQ` fallback. Null when the row has no quantity: DIA-NN's `0` for a group it identified but did not quantify is written as a null `intensity`. | Each quantified row records `openms_protein_abundance` when the file carries OpenMS' protein abundance, else `unnormalized_unique_peptide_sum`. Null when the row has no quantity. |

For either converter, a FASTA fills these three fields **only where they are null and only on target rows**; a value
from the producer is never overwritten, and proteins absent from the FASTA (such as DIA-NN's internal decoys) stay null.
For OpenMS consensusXML this matters when the file's `ProteinHit`s carry no sequence (MSV000085836).

The schema defines `molecular_weight` in **kDa** and allows it to be null.
OpenMS consensusXML uses the existing `anchor_protein` and average isotopic
masses, including the terminal water molecule, without post-translational
modifications. The calculation uses the recorded protein sequence and does not
load a FASTA or reconstruct a protein from identified peptides. Both readers
support the standard amino acids, U/O and the isobaric I/L code J; sequences with
B/Z/X or modification notation remain null.
Converters that receive an explicit molecular-weight value, such as MaxQuant,
can preserve it; see [Protein Group mappings](pg.md#tool-mappings).

### PSM evidence and computed links

DIA-NN does not produce a PSM view. OpenMS PSM `missed_cleavages` uses the same
sequence/enzyme calculation as the feature field above. If no exportable PSM
records remain, the consensusXML converter warns and skips the PSM file and its
output registration; other requested views can still be written.

`feature.psm_ids` and `feature.pg_ids` are computed links that QPX converters do
not materialize. OpenMS `psm.feature_id` records an explicit feature link when
both views are exported; unassigned PSMs and PSM-only exports can retain valid
identifications with a null link. A null `psm.feature_id` denotes no recorded
association. The [data-model join rules](data-model.md#entity-relationship-diagram)
describe how to compute links from the evidence that was exported.

## Field population: other converters

The tables below cover MaxQuant, FragPipe, Spectronaut, CDAP, mzIdentML and
QuantMS MSstats. Each row is a field that is null in some or all rows of at
least one converter's output; **Filled** means that converter fills it in every
row. A cell that names a source column and a condition ("when exported", "when
present") is filled from that column when the input has it. Other nulls are
marked as one of:

- **Not provided** — the tool output that the converter reads does not report the value.
- **Not implemented** — the converter does not write the field. When the input reports the value, or the value can be calculated from the input, the cell says so.
- **Out of scope** — the field does not apply to the input, for example spectrum arrays without spectra or quantities in identification-only input.

For these six converters, `modifications` is null exactly for unmodified
peptidoforms.

### How the cells were established

Each converter was run on the inputs below, and the populated fraction of every
field was measured. The converter code was then read to classify each null.
Cells about columns or options that these inputs do not exercise come from the
code alone.

| Converter | Measured on | From code only |
| ----------- | ------------- | ---------------- |
| MaxQuant | `tests/examples/maxquant`: `maxquant_full` (PXD001819 `evidence.txt` and `proteinGroups.txt`; 243,821 features, 27,072 pg rows) and `maxquant_simple` (TMT `msms.txt` and `evidence.txt`; 1,000 PSMs, 9,301 features) | `1/K0` (no timsTOF input) |
| FragPipe | PXD085647 (TMTpro 16 channels, one run; 7,395 PSMs, 4,356 features, 1,241 protein groups); a TMT `psm.tsv` from the bigbio/pmultiqc test data (1,000 PSMs); PXD081924 (FragPipe 22 label-free, three one-run experiments; 10,313 PSMs, 3,376 ions, 482 protein groups) | Multi-run experiments and match-between-runs rows |
| Spectronaut | `tests/examples/spectronaut/small`, without SDRF (168 features, 70 pg rows); PXD057632 pivot reports (three-channel precursors, 334,345 features; protein level, 44,383 pg rows) | `PG.Genes`, `PG.ProteinNames` and the SDRF-based `missed_cleavages`, which that report does not exercise; `FG.ApexIonMobility`, `EG.StartRT` and `EG.EndRT` (no input with them) |
| CDAP | `tests/examples/cdap` and CPTAC PDC000451 (12 TMT fractions; 217,102 PSMs, 151,442 features, 598,521 pg rows) | — |
| mzIdentML | `tests/examples/mzidentml`: PXD054720 with its MGF and `--include-spectra` (13,386 PSMs, 123 pg rows) and the four mzIdentML 1.3 cross-linking examples | PEP, `MS:1003044` and protein-group q-value terms, which none of these files contains |
| QuantMS MSstats | quantms `msstats_in` files of PXD018450 (OpenMS DDA LFQ; 448,395 features), PXD022992 (DIA-NN; 539,450 features) and PXD070151 (OpenMS TMT; 12,054 features) | — |

### Feature fields (other converters)

mzIdentML produces no feature view.

| Field | MaxQuant | FragPipe | Spectronaut | CDAP | QuantMS MSstats |
| ------- | ---------- | ---------- | ------------- | ------ | ----------------- |
| `posterior_error_probability` | `PEP`. Null for match-between-runs rows (`Type` `MULTI-MATCH`), which have no MS/MS. | 1 − `PeptideProphet Probability` (`Probability` in FragPipe 22) of the matching PSM (see below). | `EG.PEP`, when exported. | Not provided: `.psm` files report MS-GF+ scores and q-values, not a PEP. | Not provided: the MSstats input has no identification scores. |
| `peptide_qvalue` | Not provided: `evidence.txt` reports a PEP, not a peptide q-value. | Not implemented. | Not implemented: the precursor-level `EG.Qvalue` is stored as `qvalue` in `additional_scores`. | Not implemented: the best `PepQvalue` is stored as `peptide_level_qvalue` in `additional_scores`. | Not provided. |
| `calculated_mz` / `observed_mz` | Filled. | `calculated_mz`: `Calculated M/Z` of the matching PSM, else `M/Z`, which is also calculated. `observed_mz`: `Calibrated Observed M/Z` of the matching PSM; null without one. | `calculated_mz` from the peptidoform; `observed_mz` from `FG.PrecMzCalibrated`, when exported. `FG.PrecMz` is theoretical and is not read. | Filled. | `calculated_mz`: not implemented, although it can be calculated from the peptidoform and charge. `observed_mz`: not provided. |
| `mass_error_ppm` | From the experimental and calculated m/z; null for `Type` `MSMS` rows, where MaxQuant reports no mass error. | From the matching PSM. | Not implemented: the report's mass-accuracy columns (`F.*MassAccuracy (PPM)`) are fragment-level. | Filled. | Not provided. |
| `additional_scores` | `Score` and `Delta score`. Null for match-between-runs rows. | Not implemented. | Spectronaut q-values and C-scores, when exported. | Filled. | Not provided. |
| `predicted_rt` | Not provided. | Not implemented. | `EG.RTPredicted`, when exported. | Not provided. | Not provided. |
| `quantification_unit_id` | Not implemented. | The FragPipe experiment. | Not implemented. | Not implemented. | Not implemented. |
| `cv_params` | Not implemented. | Not implemented. | `EG.iRTPredicted`, `EG.iRTEmpirical` and `EG.IsImputed`, when exported. | Experiment type and PSM count. | Not implemented. |
| `scan` | `MS/MS scan number`. Empty for match-between-runs rows. | From the matching PSM. | Not implemented: always empty. | Filled. | The `scan=` number in `Reference`, which only isobaric inputs carry. Empty for label-free inputs, whose `Reference` is a file name, and DIA-NN inputs, which have no `Reference`. |
| `rt` | Filled. | Not implemented. | Filled. | Filled. | `RetentionTime`, which only isobaric inputs carry. Not provided by label-free and DIA-NN inputs. |
| `ion_mobility` | `1/K0`, when present. | From the matching PSM's `Ion Mobility`, when present. | `FG.ApexIonMobility`, when exported; null when 0. | Not provided. | Not provided. |
| `ccs` | `CCS`, when present. | Not provided. | Not provided. | Not provided. | Not provided. |
| `compensation_voltage` | Not implemented. | `Compensation Voltage`; 0 when the column is absent, because the field is part of FragPipe's feature identity. | Not implemented. | Not provided. | Not provided. |
| `missed_cleavages` | Filled from `Missed cleavages`. | From the matching PSM. | Calculated from the peptide sequence and the SDRF enzyme; null without `--sdrf-file`. `PEP.NrOfMissedCleavages` is not read. | Not implemented: no enzyme is declared to the converter. | Not implemented, although the required SDRF declares the enzyme. |
| `intensities` | TMT/iTRAQ: one entry per channel the SDRF assigns in the run, `Reporter intensity N` being the channel at position N (N + 1 for 0-based columns) of the plex the SDRF labels belong to. SILAC: one entry per label (`Intensity L/M/H`). The intensity is null for a channel or label MaxQuant did not quantify, for both labels of a peptide without a labelled residue, and for label-free evidence without an MS1 quantity (`Type` `MSMS`), whose `Intensity` MaxQuant leaves blank or 0. | Label-free: `<exp> Intensity`. TMT/iTRAQ: one entry per channel of `--annotation-file`, null for FragPipe's 0. An isobaric SDRF or a one-plex `ion.tsv` without the annotation stops the conversion. | `FG.Quantity`, or one entry per `EG.Channel<N>Quantity`: the report does not name its channels, so they take the SDRF's SILAC labels lightest first, or `SILAC light`/`medium`/`heavy` without an SDRF. Null for a value Spectronaut reports as `Filtered`. | Null when no reporter channel (TMT) or precursor area (label-free) is positive. | Filled. |
| `additional_intensities` | TMT corrected reporter intensities, null when no reporter intensity is positive. SILAC: `ratio_h_l` and `ratio_h_l_normalized` on the heavy label. Label-free: `MS/MS count` as `spectral_count` (0 for match-between-runs rows), also for evidence without an MS1 quantity. | TMT/iTRAQ: `Intensity` as `ms1_intensity`, once per feature under its first quantified channel. | `PG.Quantity`, when exported. Null for a report with channels, where it belongs to no single label. | Not implemented. | Not provided. |
| `pg_global_qvalue` | `Q-value` of the feature's group in `proteinGroups.txt`. Null without `--protein-groups-file`. | Not implemented. | `PG.Qvalue`, when exported. | Not provided: `.psm` files have no protein-level q-value. | Not provided. |
| `pg_positions` | Not implemented: `evidence.txt` has no coordinates, and the `Start position` / `End position` of `peptides.txt` are not read. | Not implemented. | Not implemented. | Not provided: protein cells give the flanking residues, stored as `pre` / `post` in `pg_accessions`, but no coordinates. | Not provided. |
| `gg_accessions` / `gg_names` | `Gene names`, or the group's genes from `proteinGroups.txt` (`Gene names`, or `GN=` in `Fasta headers`). Both fields hold gene symbols. | `gg_names` from `Gene`. `gg_accessions`: not implemented. | `PG.Genes`, when exported. | Not provided: RefSeq accessions without gene names. | Not provided. |
| `identification_type` | From `Type`: `direct` (`MSMS`, `MULTI-MSMS`, `ISO-MSMS`), `transferred` (`MULTI-MATCH`, `MULTI-MATCH-MSMS`), `secondary` (`MULTI-SECPEP`). | From `Match Type`: `direct` (`MS/MS`), `transferred` (`MBR`). | From `EG.IsImputed`: `transferred` when true, else `direct`. Null without the column. | `direct`. | Not provided. |
| `match_score` / `match_qvalue` / `match_rt_delta` / `match_mz_delta` | `Match score`, `Match q-value`, `Match time difference`, `Match m/z difference`. | Not provided. | Not provided. | Not provided. | Not provided. |
| `imputed_labels` | Not provided. | Not provided. | The row's labels when `EG.IsImputed` is true. A pg row's `is_imputed` is true when every precursor of the group in the run is. | Not provided. | Not provided. |
| `parent_ion_fraction` | `PIF`; null for MaxQuant's 0, which means not computed. | Not provided: on the PSM. | Not provided. | Not provided: on the PSM, from the first percentage of `PrecursorPurity`. | Not provided. |
| `id_run_file_name` | The run. Null for match-between-runs rows. | The experiment's raw run. Null for match-between-runs rows (`Match Type` `MBR`) and for an experiment of several raw runs, where the identifying fraction is unknown. | Filled. | Filled. | Filled. |
| `rt_start` / `rt_stop` | Filled from the calibrated retention time start and finish. | Not implemented. | `EG.StartRT` / `EG.EndRT`, when exported. | Not provided: `.psm` files give one retention time (`RTAtPrecursorHalfElution`) per PSM. | Not provided. |

FragPipe's `combined_*` reports hold MS1 intensities only, so the reporter
quantities of a TMT/iTRAQ study are read per plex: `--annotation-file` names the
channel columns of that plex's `ion.tsv` and `protein.tsv`, a feature takes the
plex's first raw run and the output prefix as `quantification_unit_id`, and a pg
row groups all the runs `psm.tsv` names. A precursor with no reporter quantity
writes no feature. Convert each plex with its own `--output-prefix`.

A FragPipe feature is one experiment's quantity. Its `run_file_name` is the
first raw run of the experiment (from `experiment_annotation.tsv`, else the SDRF
`source name` of the same name), which is a member of the experiment's pg
`grouped_runs`; an experiment neither resolves stops the conversion. Its
`pg_accessions` are the members of the group led by `Protein`: the lead and its
`Indistinguishable Proteins` in `combined_protein.tsv` (only the lead without
`--pg-file`). `unique` is true when `Protein` and `Mapped Proteins` name one
protein, and the pepmap view lists those proteins for each converted
peptidoform.

FragPipe features take PEP, mass error, scan, ion mobility and missed cleavages
from the best `psm.tsv` PSM (lowest PEP, then lowest scan) with the same
peptidoform and charge in the experiment's raw file. These fields are null
without `--psm-file`, for an experiment that maps to more than one raw file,
and when that raw file has no such PSM. A PSM with an unassigned mass shift is
not used, since the ion does not carry the shift.

A FragPipe PSM whose measured mass (`Calculated Peptide Mass` + `Delta Mass`)
differs from its peptidoform carries the difference as a mass shift, as in an
open or mass-offset search: on the residue when `MSFragger Localization` marks
one, otherwise at an unknown position with any marked residues as
`modification_candidates`. Its `calculated_mz` includes the shift, so
`mass_error_ppm` says nothing about the instrument for such a PSM. A difference
within 0.1 Da (or 20 ppm, the wider) of −1 to 3 <sup>13</sup>C spacings is an
isotope error of the precursor selection, not a modification: the peptidoform
leaves it out and `mass_error_ppm` shows it. Real shifts inside that window, such
as a deamidation that Assigned Modifications do not list, cannot be told apart
from it. The glycan composition is the name of the `modifications` entry of the
PSM; feature entries name the mass, and join the PSM through the peptidoform.

### Site localization

Per-site probabilities go to `modification_candidates`, one entry per scored
residue, and to a `localization_probability` score on the sites the peptidoform
occupies. A modification placed once whose probabilities sum to about 1 is one
alternative group; several placements, or independent occupancies that sum to
more, stay ungrouped.

| Converter | View | Source |
| --------- | ---- | ------ |
| MaxQuant | PSM, feature | every `<mod> Probabilities` column of `msms.txt` and `evidence.txt` |
| FragPipe | PSM | PTMProphet columns such as `STY:79.9663`; `MSFragger Localization` for a mass shift |
| DIA-NN | feature | `Site.Occupancy.Probabilities`; `PTM.Site.Confidence` and the `Peptidoform.Q.Value` columns are `additional_scores` |
| Spectronaut | feature | `EG.PTMLocalizationProbabilities`; `EG.PTMAssayProbability` is an `additional_scores` entry |
| mzIdentML | PSM | PSM-level statistics in the MS:1002505 format. Probabilities and other scores of unassigned residues need the `modification index` of a `SearchModification`; without it only the assigned residue is scored, under the term's own name. A scalar value stays in `additional_scores` |
| OpenMS | PSM, feature | the per-site values its idXML meta values hold |

FragPipe's feature view has no localization because `combined_ion.tsv` carries
none. The site reports of DIA-NN and Spectronaut are not converted to
[`site.parquet`](site.md).

`consensus_rt` (an OpenMS consensusXML field), `ion_mobility_start` and
`ion_mobility_stop` are null for all five converters. `psm_ids` and `pg_ids` are
computed links that no converter materializes.

### PSM fields (other converters)

Spectronaut and QuantMS MSstats produce no PSM view.

| Field | MaxQuant | FragPipe | CDAP | mzIdentML |
| ------- | ---------- | ---------- | ------ | ----------- |
| `posterior_error_probability` | Filled from `PEP`. | Filled: 1 − `PeptideProphet Probability` (`Probability` in FragPipe 22). | Not provided. | A PSI-MS PEP term (`MS:1001493`, `MS:1002352`) on the identification item, when present. |
| `mass_error_ppm` | From the experimental and calculated m/z; null for `Type` `MSMS` rows, where MaxQuant reports no mass error. | Filled. | Filled. | From the experimental and calculated m/z. Null when the file omits the calculated m/z. |
| `additional_scores` | Filled. | `Hyperscore`, `Expectation`, `Purity` (as `parent_ion_fraction`), and `Glycan Score` and `Glycan q-value` (as `glycan_score` and `glycan_qvalue`). | Filled. | Score terms on the identification item. Null when it has none. |
| `rt` | Filled. | Filled. | Filled. | The `scan start time` of the spectrum result. Null when the file omits it. |
| `ion_mobility` | `1/K0`, when present. | `Ion Mobility`, when present. | Not provided. | Not implemented. |
| `ccs` | `CCS`, when present. | Not provided. | Not provided. | Not provided. |
| `missed_cleavages` | Filled from `Missed cleavages`. | Filled from `Number of Missed Cleavages`. | Not implemented: no enzyme is declared to the converter. | Only the `MS:1003044` term on the identification item. Not calculated from the enzyme that the file declares. |
| `protein_accessions` | `Proteins`, which MaxQuant leaves blank for some reverse hits. | `Protein` and `Mapped Proteins`. | Filled. | Filled. |
| `cross_links`, `crosslinked_peptidoform`, `protein_starts`, `partner_protein_accessions`, `partner_protein_starts`, `crosslink_decoy_class` | Out of scope. | Out of scope. | Out of scope. | Cross-linking identifications only. |
| `ion_type_array` / `charge_array` | `Matches` of `msms.txt`, with `--spectral-data`, when it has one token per peak. | Not implemented. | Not implemented. | Not implemented. |
| `mz_array` / `intensity_array` | `Masses` / `Intensities` of `msms.txt`, with `--spectral-data`. | Not provided: `psm.tsv` has no spectra. | Not provided. | From the MGF, with `--mgf-path` and `--include-spectra`. |
| `feature_id` | Not implemented: the `Evidence ID` of `msms.txt` is not read. | Not implemented. | Not implemented: features are aggregated from PSMs, but the link is not recorded. | Out of scope: there is no feature view. |

`predicted_rt`, `charge_array`, `ion_type_array` and `ion_mobility_array` are
not implemented by any of the four converters. `qpxc convert mz` writes
spectra from mzML files as a separate `mz` view.

### Protein group fields (other converters)

QuantMS MSstats produces no pg view.

| Field | MaxQuant | FragPipe | Spectronaut | CDAP | mzIdentML |
| ------- | ---------- | ---------- | ------------- | ------ | ----------- |
| `pg_names` | `Protein names`; the accessions when the column is absent. | `Description`, the protein description. | `PG.ProteinNames`, when exported. | The accessions: `.psm` files carry no protein names. | Not implemented. |
| `gg_accessions` / `gg_names` | `Gene names`, or `GN=` in `Fasta headers`. | `Gene`. | `PG.Genes`, when exported. | Not provided. | Not implemented. |
| `global_qvalue` | `Q-value`. | Not implemented: the group q-value is written to `pg_qvalue`. | `PG.Qvalue`, when exported. | Not provided. | The protein-group q-value term `MS:1002373`, when present. |
| `pg_qvalue` | The same experiment-level `Q-value` as `global_qvalue`. | The q-value column of `combined_protein.tsv`, or 1 − `Protein Probability` when that is the only one. | `PG.QValue (Run-Wise)`, when exported. | Not provided. | Not implemented. |
| `label` / `intensity` | `intensity` is null for a group with MS/MS counts but no intensity in that experiment. | `<exp> Intensity` (`Total Intensity` in older reports, else `MaxLFQ Intensity`); null for a group with spectral counts but no intensity in that experiment. | Filled. | Null when no channel has a positive intensity. | Out of scope: identification only. |
| `additional_intensities` | `LFQ intensity`, `iBAQ` and `MS/MS count <experiment>` as `spectral_count` (label-free), corrected reporter intensities (TMT) or `Ratio X/Y [normalized] [exp]` on the numerator label (SILAC). A 0 `LFQ intensity` or `iBAQ`, which MaxQuant writes when it computed no value, is left out. Null when none is reported. | Positive `MaxLFQ Intensity`, `MaxLFQ Unique Intensity` and `Unique Intensity`, and the experiment's `Spectral Count`, `Unique Spectral Count` and `Total Spectral Count`. TMT/iTRAQ: `--tmt-report` as `tmt_integrator_log2_abundance`. | Not implemented. | Not implemented. | Out of scope. |
| `contaminant` | `Potential contaminant`. | Not implemented. | Not implemented. | Always false: not derived from the input. | Not implemented. |
| `peptide_counts` / `feature_counts` | Filled. | `peptide_counts`: `<exp> Unique Peptides` and `<exp> Total Peptides`, when both are present. `feature_counts`: not provided. | `feature_counts`: not implemented. | Filled. | Not implemented. |
| `sequence_coverage` / `molecular_weight` | `Sequence coverage [%]` / `Mol. weight [kDa]`. | `Percent Coverage` / `Protein Molecular Weight (Da)` in kDa, when present. | Not implemented. | Not provided. | Not implemented. |
| `additional_scores` | Filled. | `Protein Probability` and `Top Peptide Probability`. | Spectronaut protein-group q-values, when exported. | Filled. | Score terms on the group and its proteins. Null when there are none. |

`gg_qvalue` is null for all five converters: none reads a gene-level q-value.
`cv_params` is null except for FragPipe, which records `Protein Length` and the
`Protein Existence` level.

`qpxc transform protein-properties` and `qpxc transform gene-map` take any QPX
dataset. With the search FASTA they fill null `pg.sequence_coverage`,
`pg.molecular_weight`, `feature.pg_positions` and `gg_names`. See the
[transform guide](../guide/transform.md#protein-properties).

## Proteome Discoverer

`qpxc convert proteome-discoverer` reads the text exports of a consensus result, not the `.pdResult` file. Both header spellings of PD are accepted (`m/z [Da]` and `mz in Da`, `Abundance: 126` and `Abundance 126`). It was measured on PXD081613 (TMTpro 18-plex, Sequest and Comet with ptmRS; 109,700 PSMs of a phosphoproteome, in which Comet's isotope errors show up as +1 to +3 Da in `observed_mz`, which PD's `DeltaM` reports too), PXD073725 (TMT10, older headers, grouped abundances only; 58,948 PSMs), PXD046576 (Minora label-free; 16,359 PSMs) and PXD083739 (identification only; 29,218 PSMs). Labelled MS1 exports (`Abundance: Light`) stop the conversion.

| View | Source |
| ---- | ------ |
| `psm` | One row per run, scan, peptidoform and charge of `PSMs.txt`. Of several hits on one spectrum with the same identity, as from two search engines, the lowest rank then the lowest q-value stays, and `Identifying Node` says which. `pass_threshold` is `Confidence` = High, ranks stay as the `rank` score, and `--rank1-only` keeps only the first-ranked hits. `ptmRS <mod> Site Probabilities` give `modification_candidates`, and `Isolation Interference` gives the `parent_ion_fraction` score. |
| `feature` | The rank 1 PSMs of a precursor in a run, leaving out those `PSM Ambiguity` calls Rejected and `Quan Info` calls excluded. Reporter channels are summed over them. A label-free feature takes the largest `Precursor Abundance` (`Intensity` in an export without it), its PSM count as `spectral_count` and its `Apex RT`. The identification fields are those of the lowest q-value PSM. An export without abundance columns writes no feature view. |
| `pg` | One row per quantified file group and channel of the `IsMasterProtein` rows of `Proteins.txt`. The members of `pg_accessions` are the rows sharing its `Protein Group IDs`. `Abundance: F1: 126, Sample` is the intensity, or `Abundances (Grouped): F1, 126` when PD did not write the per-file column, and that value is normalized and scaled if the consensus step was set so. `grouped_runs` holds the fractions of the file group (`F16.1`, `F16.2`) in fraction order. An identification-only export has one `LFQ` row per protein with a null intensity. |

Reporter labels come from the column (`127N` is `TMT127N`), and must match the SDRF `comment[label]` values, which `--sdrf-file` converts to the sample and run views.

## CLI commands

| Converter | Command |
| ----------- | --------- |
| MaxQuant | `qpxc convert maxquant` |
| FragPipe | `qpxc convert fragpipe` |
| Proteome Discoverer | `qpxc convert proteome-discoverer` |
| DIA-NN | `qpxc convert diann` |
| Spectronaut | `qpxc convert spectronaut` |
| OpenMS native QPX | `qpxc convert openms` |
| OpenMS consensusXML | `qpxc convert openms-consensus` |
| CDAP | `qpxc convert cdap` |
| mzIdentML | `qpxc convert mzidentml` |
| QuantMS MSstats | `qpxc convert quantms-msstats` |
| SDRF only | `qpxc convert sdrf` |

## Input files (summary)

| Converter | Typical inputs |
| ----------- | ---------------- |
| MaxQuant | msms.txt, evidence.txt, proteinGroups.txt |
| Proteome Discoverer | PSMs.txt; optional Proteins.txt, InputFiles.txt (when PSMs.txt has no Spectrum File) and SDRF |
| FragPipe | psm.tsv, combined_ion, combined_protein; for one TMT/iTRAQ plex, psm.tsv, ion.tsv, protein.tsv, annotation.txt and optionally a TMT-Integrator abundance_protein report |
| DIA-NN | report (TSV or Parquet), pg_matrix (optional); required SDRF; optional mzML statistics folder (`--mzml-info-folder`) for feature scans and the `mz` view |
| Spectronaut | report.tsv; optional SDRF |
| OpenMS native QPX | `-out_qpx` Parquet directory; SDRF; optional companion consensusXML |
| OpenMS consensusXML | consensusXML; optional SDRF |
| CDAP | CPTAC CDAP `.psm` files in one study directory |
| mzIdentML | .mzid / .mzid.gz; optional MGF or mzML (file/folder) for spectra; optional SDRF |
| QuantMS MSstats | QuantMS-generated `*_msstats_in.csv`; required SDRF |
| SDRF | Single SDRF TSV file |

For field-level mappings from each tool’s columns to QPX, see [Tool Field Mappings](tool-mappings.md).
