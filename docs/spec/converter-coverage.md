# Converter Coverage Matrix

This page lists which QPX data views each converter produces and explains the sources and null values of selected fields.

## Views produced per converter

| Converter | PSM | Feature | PG | Pepmap | Sample | Run | Dataset | Ontology | Provenance | mz |
| ----------- | :---: | :-------: | :--: | :------: | :------: | :---: | :-------: | :--------: | :----------: | :--: |
| **MaxQuant** | Yes | Yes | Yes | No | If SDRF | If SDRF | Yes | If terms | Yes | No |
| **FragPipe** | Yes | Yes | Yes | No | If SDRF | If SDRF | Yes | If terms | Yes | No |
| **DIA-NN** | No | Yes | If PG matrix | No | Yes | Yes | Yes | If terms | Yes | If MS2 info |
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
- **If terms** — the view is produced when the converter discovers resolvable ontology entries.
- **If in file** — the view is produced only when the input contains the source section (mzIdentML: a `ProteinDetectionList`).
- **If PG matrix** — the view is produced only when `--pg-matrix-path` is provided.
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
| `cv_params` | Tool-specific metrics such as `proteotypic`, `irt`, `predicted_irt`, `iim` and `precursor_quantification_score`, when reported. | Null: the converter has no corresponding source metrics to put in this field. These are separate from `pg.cv_params`. |
| `pg_global_qvalue` | Read from `Global.PG.Q.Value`, when present. | The resolved protein group's global q-value, using the same derivation as `pg.global_qvalue`. Null when the group is ambiguous or has no declared q-value. |
| `gg_accessions` / `gg_names` | Gene symbols from `Genes`, when present. | Gene symbols parsed from `GN=` in protein-hit descriptions, through the resolved protein group. Null when that annotation is absent or the group cannot be resolved. The converter does not read a FASTA to add missing genes. |
| `rt_start` / `rt_stop` | Read from `RT.Start` / `RT.Stop` and converted from minutes to seconds. Null when those columns or values are absent. | Null: consensusXML feature elements provide an apex RT, not the feature's integration-window boundaries. The apex is stored in `rt`. |
| `id_run_file_name` | The normalized report run name. | The resolved run of a direct identification of the exported peptide with a spectrum reference. Null for transferred features without their own identification, unresolved origins or conflicting peptide assignments. |
| `pg_positions` | Not populated from the report, which has no peptide-to-protein coordinates. Filled from a FASTA with `--fasta` or `qpxc transform protein-properties`: every one-based occurrence in each group member. | Known peptide-evidence `start` / `end` coordinates, converted from zero-based inclusive to one-based inclusive. For a resolved protein group, positions on its members; for a peptide shared across inferred groups (null `pg_accessions`), positions on every evidence protein, since each entry names its own protein. Null when no usable positions remain. A FASTA is not needed for coordinates already recorded in consensusXML. |
| `peptide_qvalue` | Null: DIA-NN `Q.Value` is precursor-level and is stored as `precursor_qvalue` in `additional_scores`. | Read from an identification score explicitly declared to be a q-value/FDR. Null without a suitable identification in that run; raw search scores are not substituted. |
| `intensities` / `additional_intensities` | `Precursor.Quantity`; `PG.MaxLFQ`, `Precursor.Normalised` and `Ms1.Area` as `maxlfq`, `precursor_normalised` and `ms1_area`. DIA-NN writes 0 for a quantity it did not measure, so only positive values are kept: an unquantified channel has no entry, and a feature without any quantity keeps its row with `intensities` = null. | The consensus feature's per-map sub-feature intensities; a map element with intensity ≤ 0 has no entry. `additional_intensities` is not populated. |

Protein-group ambiguity includes shared peptides spanning multiple inferred
groups and conflicting assignments within a run. These features retain their
quantification while the ambiguous protein-group fields remain null. See the
[OpenMS conversion guide](../guide/convert.md#openms-consensus-description) for
the evidence-selection rules, and [Feature confidence](feature.md#protein-group-fields)
for q-value semantics.

### Protein group fields

| Field | DIA-NN | OpenMS consensusXML |
| ------- | -------- | --------------------- |
| `pg_names` | Split from the report's `Protein.Names`. | UniProt entry names taken from `db\|ACCESSION\|NAME` accessions, aligned with `pg_accessions`. Null when any member is a bare accession with no entry name, so the list never misaligns; a bare accession is not repeated as its own name. |
| `contaminant` | True when any member accession carries the `CONTAM` marker of the quantms/sdrf-pipelines contaminant database. | Same rule, shared with DIA-NN (`qpx.converters.utils.is_contaminant_accession`). |
| `global_qvalue` | Read from `Global.PG.Q.Value`, when present. | The lowest available member-protein q-value, only when the source protein score is declared to be a q-value/FDR. Null if no member has such a score. |
| `pg_qvalue` | Read from run-level `PG.Q.Value`, when present. | Null: no separate run-level protein-group q-value is reported by this path. The experiment-level value belongs in `global_qvalue`. |
| `gg_qvalue` | Read from `GG.Q.Value`, when present. | Null: this path has no gene-group q-value. Protein-level confidence is not a gene-level q-value. |
| `label` / `intensity` | The first positive value of `PG.Quantity`, `PG.MaxLFQ` and, for label-free data, the PG matrix; `quantification_method` in `cv_params` names the source. DIA-NN writes 0 for a quantity it did not compute, as QuantUMS does in `PG.MaxLFQ`: a group (or, in channel mode, a channel) without a positive quantity keeps its `label` with `intensity` = null. | The interim unique-peptide sum (see `cv_params`); null for a group with a PSM in the unit but no quantity. |
| `additional_intensities` | `PG.MaxLFQ` as `maxlfq` beside `PG.Quantity`; a 0 is left out. | Not populated. |
| `additional_scores` | Includes `PG.Q.Value` as `pg_qvalue` and `PG.PEP` as `pg_pep` when available. | The existing `anchor_protein`'s `Posterior Probability_score`, stored as `posterior_probability` with `higher_better: true`. Null if absent, invalid or conflicting for that anchor. |
| `sequence_coverage` | Not populated from the report. Filled from a FASTA with `--fasta` or `qpxc transform protein-properties`: percent of the anchor covered by the dataset's target peptides whose evidence names it. | The existing anchor's recorded `ProteinHit.coverage` in percent. Unknown values, missing annotations and conflicting coverage values for that anchor remain null. Other group members are not substitutes. |
| `molecular_weight` | Not populated from the report. Filled from a FASTA with `--fasta` or `qpxc transform protein-properties`: anchor average mass in kDa. | The theoretical average molecular weight of the existing anchor's complete, unmodified `ProteinHit.sequence`, in kDa. Null for missing, ambiguous, modified or conflicting sequences; other group members are not substitutes. |
| `cv_params` | Each quantified row records whether its primary intensity came from `PG.Quantity` or the `PG.MaxLFQ` fallback. | Each quantified row records `unnormalized_unique_peptide_sum`, or the corresponding top-N method selected by `--pg-top`. Null when the row has no quantity. |

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
| FragPipe | A TMT `psm.tsv` from the bigbio/pmultiqc test data (1,000 PSMs) | Feature and pg views: no `combined_ion.tsv` or `combined_protein.tsv` beyond the one-row inputs of the repository tests |
| Spectronaut | `tests/examples/spectronaut/small`, without SDRF (168 features, 70 pg rows) | `PG.Genes`, `PG.ProteinNames` and the SDRF-based `missed_cleavages`, which that report does not exercise |
| CDAP | `tests/examples/cdap` and CPTAC PDC000451 (12 TMT fractions; 217,102 PSMs, 151,442 features, 598,521 pg rows) | — |
| mzIdentML | `tests/examples/mzidentml`: PXD054720 with its MGF and `--include-spectra` (13,386 PSMs, 123 pg rows) and the four mzIdentML 1.3 cross-linking examples | PEP, `MS:1003044` and protein-group q-value terms, which none of these files contains |
| QuantMS MSstats | quantms `msstats_in` files of PXD018450 (OpenMS DDA LFQ; 448,395 features), PXD022992 (DIA-NN; 539,450 features) and PXD070151 (OpenMS TMT; 12,054 features) | — |

### Feature fields (other converters)

mzIdentML produces no feature view.

| Field | MaxQuant | FragPipe | Spectronaut | CDAP | QuantMS MSstats |
| ------- | ---------- | ---------- | ------------- | ------ | ----------------- |
| `posterior_error_probability` | `PEP`. Null for match-between-runs rows (`Type` `MULTI-MATCH`), which have no MS/MS. | 1 − `PeptideProphet Probability` of the matching PSM (see below). | `EG.PEP`, when exported. | Not provided: `.psm` files report MS-GF+ scores and q-values, not a PEP. | Not provided: the MSstats input has no identification scores. |
| `peptide_qvalue` | Not provided: `evidence.txt` reports a PEP, not a peptide q-value. | Not implemented. | Not implemented: the precursor-level `EG.Qvalue` is stored as `qvalue` in `additional_scores`. | Not implemented: the best `PepQvalue` is stored as `peptide_level_qvalue` in `additional_scores`. | Not provided. |
| `calculated_mz` / `observed_mz` | Filled. | Filled. | Filled. | Filled. | `calculated_mz`: not implemented, although it can be calculated from the peptidoform and charge. `observed_mz`: not provided. |
| `mass_error_ppm` | `Mass error [ppm]`, which MaxQuant leaves blank for `Type` `MSMS` rows. | From the matching PSM. | Not implemented: the report's mass-accuracy columns (`F.*MassAccuracy (PPM)`) are fragment-level. | Filled. | Not provided. |
| `additional_scores` | `Score` and `Delta score`. Null for match-between-runs rows. | Not implemented. | Spectronaut q-values and C-scores, when exported. | Filled. | Not provided. |
| `predicted_rt` | Not provided. | Not implemented. | `EG.RTPredicted`, when exported. | Not provided. | Not provided. |
| `quantification_unit_id` | Not implemented. | The FragPipe experiment. | Not implemented. | Not implemented. | Not implemented. |
| `cv_params` | Not implemented. | Not implemented. | `EG.iRTPredicted`, `EG.iRTEmpirical` and `EG.IsImputed`, when exported. | Experiment type and PSM count. | Not implemented. |
| `scan` | `MS/MS scan number`. Empty for match-between-runs rows. | From the matching PSM. | Not implemented: always empty. | Filled. | The `scan=` number in `Reference`, which only isobaric inputs carry. Empty for label-free inputs, whose `Reference` is a file name, and DIA-NN inputs, which have no `Reference`. |
| `rt` | Filled. | Not implemented. | Filled. | Filled. | `RetentionTime`, which only isobaric inputs carry. Not provided by label-free and DIA-NN inputs. |
| `ion_mobility` | `1/K0`, when present. | From the matching PSM's `Ion Mobility`, when present. | Not implemented. | Not provided. | Not provided. |
| `compensation_voltage` | Not implemented. | `Compensation Voltage`; 0 when the column is absent, because the field is part of FragPipe's feature identity. | Not implemented. | Not provided. | Not provided. |
| `missed_cleavages` | Filled from `Missed cleavages`. | From the matching PSM. | Calculated from the peptide sequence and the SDRF enzyme; null without `--sdrf-file`. `PEP.NrOfMissedCleavages` is not read. | Not implemented: no enzyme is declared to the converter. | Not implemented, although the required SDRF declares the enzyme. |
| `intensities` | TMT: null when every reporter intensity is zero. | Filled. | Filled. | Null when no reporter channel (TMT) or precursor area (label-free) is positive. | Filled. |
| `additional_intensities` | TMT corrected reporter intensities. Null for label-free input and when no reporter intensity is positive. | Not implemented. | `PG.Quantity` and `FG.PrecMzCalibrated`, when exported. | Not implemented. | Not provided. |
| `pg_global_qvalue` | `Q-value` of the feature's group in `proteinGroups.txt`. Null without `--protein-groups-file`. | Not implemented. | `PG.Qvalue`, when exported. | Not provided: `.psm` files have no protein-level q-value. | Not provided. |
| `pg_positions` | Not implemented: `evidence.txt` has no coordinates, and the `Start position` / `End position` of `peptides.txt` are not read. | Not implemented. | Not implemented. | Not provided: protein cells give the flanking residues, stored as `pre` / `post` in `pg_accessions`, but no coordinates. | Not provided. |
| `gg_accessions` / `gg_names` | `Gene names`, or the group's genes from `proteinGroups.txt` (`Gene names`, or `GN=` in `Fasta headers`). Both fields hold gene symbols. | `gg_names` from `Gene`. `gg_accessions`: not implemented. | `PG.Genes`, when exported. | Not provided: RefSeq accessions without gene names. | Not provided. |
| `id_run_file_name` | The run. Null for match-between-runs rows. | The FragPipe experiment. | Filled. | Filled. | Filled. |
| `rt_start` / `rt_stop` | Filled from the calibrated retention time start and finish. | Not implemented. | Not implemented. | Not provided: `.psm` files give one retention time (`RTAtPrecursorHalfElution`) per PSM. | Not provided. |

FragPipe features take PEP, mass error, scan, ion mobility and missed cleavages
from the `psm.tsv` PSM with the same peptidoform and charge in the experiment's
raw file. These fields are null without `--psm-file`, for an experiment that
maps to more than one raw file, and when that raw file has no such PSM.

`consensus_rt` (an OpenMS consensusXML field), `ion_mobility_start` and
`ion_mobility_stop` are null for all five converters. `psm_ids` and `pg_ids` are
computed links that no converter materializes.

### PSM fields (other converters)

Spectronaut and QuantMS MSstats produce no PSM view.

| Field | MaxQuant | FragPipe | CDAP | mzIdentML |
| ------- | ---------- | ---------- | ------ | ----------- |
| `posterior_error_probability` | Filled from `PEP`. | Filled: 1 − `PeptideProphet Probability`. | Not provided. | A PSI-MS PEP term (`MS:1001493`, `MS:1002352`) on the identification item, when present. |
| `mass_error_ppm` | `Mass error [ppm]`, which MaxQuant leaves blank for `Type` `MSMS` rows. | Filled. | Filled. | From the experimental and calculated m/z. Null when the file omits the calculated m/z. |
| `additional_scores` | Filled. | Filled. | Filled. | Score terms on the identification item. Null when it has none. |
| `rt` | Filled. | Filled. | Filled. | The `scan start time` of the spectrum result. Null when the file omits it. |
| `ion_mobility` | `1/K0`, when present. | `Ion Mobility`, when present. | Not provided. | Not implemented. |
| `missed_cleavages` | Filled from `Missed cleavages`. | Filled from `Number of Missed Cleavages`. | Not implemented: no enzyme is declared to the converter. | Only the `MS:1003044` term on the identification item. Not calculated from the enzyme that the file declares. |
| `protein_accessions` | `Proteins`, which MaxQuant leaves blank for some reverse hits. | Filled from `Protein`. | Filled. | Filled. |
| `cross_links` | Out of scope. | Out of scope. | Out of scope. | Cross-linking identifications only. |
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
| `label` / `intensity` | `intensity` is null for a group with MS/MS counts but no intensity in that experiment. | `intensity` is null for a group with spectral counts but no intensity in that experiment. | Filled. | Null when no channel has a positive intensity. | Out of scope: identification only. |
| `additional_intensities` | `LFQ intensity` and `iBAQ` (label-free) or corrected reporter intensities (TMT). Null when none is reported. | `MaxLFQ Intensity`, when present. | Not implemented. | Not implemented. | Out of scope. |
| `contaminant` | `Potential contaminant`. | Not implemented. | Not implemented. | Always false: not derived from the input. | Not implemented. |
| `peptide_counts` / `feature_counts` | Filled. | Filled. | `feature_counts`: not implemented. | Filled. | Not implemented. |
| `sequence_coverage` / `molecular_weight` | `Sequence coverage [%]` / `Mol. weight [kDa]`. | `Percent Coverage` / `Protein Molecular Weight (Da)` in kDa, when present. | Not implemented. | Not provided. | Not implemented. |
| `additional_scores` | Filled. | Spectral count. | Spectronaut protein-group q-values, when exported. | Filled. | Score terms on the group and its proteins. Null when there are none. |

`gg_qvalue` and `cv_params` are null for all five converters: none reads a
gene-level q-value or records the quantification method.

`qpxc transform protein-properties` and `qpxc transform gene-map` take any QPX
dataset. With the search FASTA they fill null `pg.sequence_coverage`,
`pg.molecular_weight`, `feature.pg_positions` and `gg_names`. See the
[transform guide](../guide/transform.md#protein-properties).

## CLI commands

| Converter | Command |
| ----------- | --------- |
| MaxQuant | `qpxc convert maxquant` |
| FragPipe | `qpxc convert fragpipe` |
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
| FragPipe | psm.tsv, combined_ion, combined_protein |
| DIA-NN | report (TSV or Parquet), pg_matrix (optional); required SDRF; optional mzML statistics folder (`--mzml-info-folder`) for feature scans and the `mz` view |
| Spectronaut | report.tsv; optional SDRF |
| OpenMS native QPX | `-out_qpx` Parquet directory; SDRF; optional companion consensusXML |
| OpenMS consensusXML | consensusXML; optional SDRF |
| CDAP | CPTAC CDAP `.psm` files in one study directory |
| mzIdentML | .mzid / .mzid.gz; optional MGF or mzML (file/folder) for spectra; optional SDRF |
| QuantMS MSstats | QuantMS-generated `*_msstats_in.csv`; required SDRF |
| SDRF | Single SDRF TSV file |

For field-level mappings from each tool’s columns to QPX, see [Tool Field Mappings](tool-mappings.md).
