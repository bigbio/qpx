# Tool Field Mappings

!!! info "Reference for converter developers"
    This page is a reference for converter developers. For schema definitions, see the individual view pages: [PSM](psm.md), [Feature](feature.md), [Protein Group](pg.md), [XIC](xic.md).

This page documents how fields from different proteomics search engines and tools map to the corresponding QPX fields. Each section covers a specific QPX view and lists the source column name in each supported tool.

A dash (`---`) indicates that the tool does not provide a direct mapping for that field.

For conditional fields, derived values and reasons for nulls in DIA-NN and OpenMS
consensusXML, see [Field population](converter-coverage.md#field-population-dia-nn-and-openms-consensusxml).

## PSM Field Mappings

The PSM (Peptide Spectrum Match) view captures spectrum-level identification results. The table below shows how each QPX PSM field maps to columns in common proteomics tools.

| QPX Field | MaxQuant | DIA-NN | FragPipe |
| --- | --- | --- | --- |
| `sequence` | Sequence | Stripped.Sequence | Peptide |
| `peptidoform` | Modified sequence, plus the fixed modifications MaxQuant leaves out (SDRF, else parameters.txt or summary.txt; the isobaric label from the SDRF channel labels) and the SILAC label modifications of the `Labeling state` (summary.txt `Labels<n>`) | Modified.Sequence | Modified Peptide |
| `label` | Labeling state (SILAC) | --- | --- |
| `charge` | Charge | Precursor.Charge | Charge |
| `posterior_error_probability` | PEP | PEP | 1 − PeptideProphet Probability (Probability in FragPipe 22) |
| `is_decoy` | Reverse | --- | --- |
| `calculated_mz` | From the peptidoform (`Mass` when a modification has no known mass) | --- | Calculated M/Z |
| `observed_mz` | calculated_mz + m/z × Mass error [ppm] × 10⁻⁶; msms.txt `m/z` is theoretical, and it and its mass error leave isobaric labels out | --- | Observed M/Z |
| `protein_accessions` | Proteins | Protein.Group | Protein, Mapped Proteins |
| `run_file_name` | Raw file | Run | Spectrum File |
| `scan` | Scan number | --- | Scan |
| `rt` | Retention time | RT | Retention |

!!! note "PEP direction"
    `posterior_error_probability` is the probability that the PSM is incorrect — **lower values indicate higher confidence** (lower is better). All major tools (Percolator, MaxQuant) export PEP directly as P(incorrect). FragPipe exports PeptideProphet Probability (P(correct)), so PEP must be computed as `1 - probability`.

## Feature Field Mappings

The Feature view captures quantified peptide features with intensity data. Features aggregate information across scans and are the primary view for DIA and label-free quantification workflows.

| QPX Field | MaxQuant | DIA-NN | FragPipe | QuantMS MSstats |
| --- | --- | --- | --- | --- |
| `sequence` | Sequence | Stripped.Sequence | Peptide | PeptideSequence (unmodified) |
| `peptidoform` | Modified sequence, plus the fixed modifications MaxQuant leaves out (SDRF, else parameters.txt or summary.txt; the isobaric label from the SDRF channel labels) | Modified.Sequence | Modified Peptide | PeptideSequence (ProForma-normalized) |
| `charge` | Charge | Precursor.Charge | --- | Charge / PrecursorCharge |
| `is_decoy` | Reverse | Decoy (`--report-decoys`), else accession prefix | --- | ProteinName prefix |
| `peptide_qvalue` | --- | --- | --- | --- |
| `calculated_mz` | From the peptidoform (`Mass` when a modification has no known mass) | From the peptidoform | Calculated M/Z (psm.tsv), else M/Z | --- |
| `observed_mz` | m/z × (1 + Mass error [ppm] × 10⁻⁶); evidence.txt `m/z` is theoretical without isobaric labels, and its mass error includes them | Precursor.Mz + Ms1.Apex.Mz.Delta (2.x); `Precursor.Mz` is theoretical | Calibrated Observed M/Z (psm.tsv); `M/Z` is theoretical | --- |
| `rt` | Retention time | RT | --- | RetentionTime |
| `rt_aligned` | Calibrated retention time | --- | --- | --- |
| `rt_start` | --- | RT.Start | --- | --- |
| `rt_stop` | --- | RT.Stop | --- | --- |
| `predicted_rt` | --- | Predicted.RT | --- | --- |
| `intensities` | Intensity; Intensity L/M/H (SILAC) | Precursor.Quantity (+ Channel for plexDIA) | Intensity | Intensity + Channel |
| `identified_labels` | Labeling State, and msms.txt Labeling state of its MS/MS (SILAC) | --- | --- | --- |
| `pg_accessions` | `Protein IDs` of the group in `Protein group IDs` (the razor group of a shared peptide); `Leading proteins` without proteinGroups.txt | Protein.Group (`DECOY_`-prefixed on decoy rows) | Protein + its Indistinguishable Proteins (combined_protein.tsv) | ProteinName |
| `anchor_protein` | First `Majority protein IDs` entry of that group; `Leading razor protein` without proteinGroups.txt | First Protein.Group accession | Protein | First ProteinName accession |
| `pg_positions` | --- | --- | --- | --- |
| `run_file_name` | Raw file | Run | First raw file of the experiment (experiment_annotation.tsv or SDRF) | Reference / Run resolved through SDRF |

DIA-NN `Q.Value` is precursor-level and is stored in `additional_scores` as `precursor_qvalue`; it is not mapped to the peptide-level `peptide_qvalue` field.

MaxQuant `unique` follows the group the feature takes. With proteinGroups.txt, it is true when the peptide is in a single group (`Protein group IDs` names one id) and that group has one protein; a shared peptide carries its razor group and is false. Without proteinGroups.txt, or for a row whose `Protein group IDs` do not resolve, the group falls back to `Leading proteins`, which lists one leader per group, so the group's size is unknown: `unique` is then true when the peptide maps to a single protein (`Proteins`), or, without that column, when it is in a single group of any size. In ontology.parquet the feature `pg_accessions` and `anchor_protein` name the column they were read from: `Protein IDs` and `Majority protein IDs` (proteinGroups.txt) when `Protein group IDs` resolve, otherwise `Leading proteins` and `Leading razor protein`.

## XIC Field Mappings

The XIC view stores extracted-ion chromatograms as one row per traced ion. DIA-NN writes a `*.xic.parquet` per run (enabled with `--xic`).

| QPX Field | DIA-NN | Notes |
| --- | --- | --- |
| `run_file_name` | `<run>.xic.parquet` | Strip the `.xic.parquet` suffix |
| `precursor_id` | `pr` | Preserve the modified sequence plus charge verbatim |
| `peptidoform` | `pr` | Strip the trailing charge; convert `(UniMod:N)` to ProForma `[UNIMOD:N]` |
| `charge` | `pr` | Trailing charge of the precursor identifier |
| `transition` | `feature` | `ms1` or a fragment label such as `y7^1` |
| `transition_type` | `feature` | `ms1` → `ms1`; `b`/`y` labels → `fragment` |
| `fragment_type` | `feature` | Ion series parsed from `y7^1` → `y` |
| `fragment_number` | `feature` | Parsed from `y7^1` → `7` |
| `fragment_charge` | `feature` | Parsed from `y7^1` → `1` |
| `neutral_loss` | `feature` | Suffix after the fragment charge, e.g. `-H2O` |
| `isotope` | — | Set to `0` (DIA-NN extracts monoisotopic targets) |
| `target_mz` | — | Not reported by DIA-NN XIC; null |
| `feature_ids` | — | Null unless the converter also links the trace to feature rows |
| `used_for_quantification` | — | Null when the producer does not mark it |
| `rt` | `rt` | Multiply by 60 (DIA-NN reports minutes) |
| `intensity` | `value` | Raw intensity |
| `mz` | — | Not reported by DIA-NN XIC; null |
| `label` | — | Set to `LFQ` for label-free DIA-NN output |
| `xic_id` | — | `<run_file_name>|<pr>|<feature>|<isotope>|<compensation_voltage>` |
| `is_decoy` | — | Optional, from the DIA-NN report/feature row |
| `additional_scores` | — | Null; DIA-NN XIC carries no trace-level scores |

DIA-NN also emits rows with `feature = index`, `rt = 0`, and `value = 0`; these are precursor chunk indexes, not measurements, and must be excluded.

The converter declares the trace selection in the file footer: `xic_selection=targets`, `xic_window=apex`, `xic_smoothed=false`, `xic_zero_points=kept`, `xic_rt_unit=seconds`, and `transition_format=mzPAF`, plus `xic_extraction_window_seconds` and `xic_qvalue_threshold` when provided.

## Protein Group Field Mappings

The Protein Group view captures protein-level quantification and inference results. The table below covers the three main supported tools.

| QPX Field | MaxQuant | DIA-NN | FragPipe |
| --- | --- | --- | --- |
| `pg_accessions` | Protein IDs | Protein.Group | Protein ID + Indistinguishable Proteins |
| `pg_names` | Protein names | Protein.Names | Description |
| `gg_accessions` | Gene names | Genes | Gene |
| `run_file_name` | combined | Run | --- |
| `peptide_counts` | Unique peptides [exp], Peptides [exp] | Unique.Stripped.Peptides | Unique Peptides [exp], Total Peptides [exp] |
| `feature_counts` | --- | Precursor.Quantity | --- |
| `global_qvalue` | Q-value | Global.PG.Q.Value | --- |
| `intensities` | Intensity [exp]; Intensity L/M/H [exp] (SILAC); Reporter intensity | PG.Quantity | Intensity [exp] |
| `additional_intensities` | LFQ intensity, iBAQ, MS/MS count; Ratio X/Y [normalized] [exp] (SILAC) | PG.Normalised, PG.MaxLFQ | MaxLFQ (Unique) Intensity, Unique Intensity, Spectral Count, Unique Spectral Count, Total Spectral Count [exp] |
| `is_decoy` | Reverse | Decoy (`--report-decoys`; accessions prefixed `DECOY_`) | --- |
| `contaminant` | Potential contaminant | --- | --- |

!!! note "Tool version considerations"
    Column names may vary between tool versions. The mappings above reflect the commonly used versions at the time of writing:

    - **MaxQuant**: v2.x (`evidence.txt`, `proteinGroups.txt`)
    - **DIA-NN**: v1.8+ (`report.tsv`, `pg_matrix.tsv`)
    - **FragPipe**: v20+ (combined output files)

!!! tip "Adding support for a new tool"
    To add a converter for a new tool, implement a mapping from the tool's output columns to the QPX fields listed in the tables above. See the [PSM](psm.md), [Feature](feature.md), and [Protein Group](pg.md) schema pages for field type requirements and nullability constraints.
