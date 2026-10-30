# Tool Field Mappings

!!! info "Reference for converter developers"
    This page is a reference for converter developers. For schema definitions, see the individual view pages: [PSM](psm.md), [Feature](feature.md), [Protein Group](pg.md).

This page documents how fields from different proteomics search engines and tools map to the corresponding QPX fields. Each section covers a specific QPX view and lists the source column name in each supported tool.

A dash (`---`) indicates that the tool does not provide a direct mapping for that field.

For conditional fields, derived values and reasons for nulls in DIA-NN and OpenMS
consensusXML, see [Field population](converter-coverage.md#field-population-dia-nn-and-openms-consensusxml).

## PSM Field Mappings

The PSM (Peptide Spectrum Match) view captures spectrum-level identification results. The table below shows how each QPX PSM field maps to columns in common proteomics tools.

| QPX Field | MaxQuant | DIA-NN | FragPipe |
| --- | --- | --- | --- |
| `sequence` | Sequence | Stripped.Sequence | Peptide |
| `peptidoform` | Modified sequence, plus the SILAC label modifications of the `Labeling state` (summary.txt `Labels<n>`) | Modified.Sequence | Modified Peptide |
| `label` | Labeling state (SILAC) | --- | --- |
| `charge` | Charge | Precursor.Charge | Charge |
| `posterior_error_probability` | PEP | PEP | --- |
| `is_decoy` | Reverse | --- | --- |
| `calculated_mz` | From the peptidoform (`Mass` when a modification has no known mass) | --- | Calculated M/Z |
| `observed_mz` | m/z × (1 + Mass error [ppm] × 10⁻⁶); MaxQuant's `m/z` is theoretical | --- | Observed M/Z |
| `protein_accessions` | Proteins | Protein.Group | Protein |
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
| `peptidoform` | Modified sequence | Modified.Sequence | Modified Peptide | PeptideSequence (ProForma-normalized) |
| `charge` | Charge | Precursor.Charge | --- | Charge / PrecursorCharge |
| `is_decoy` | Reverse | Decoy (`--report-decoys`), else accession prefix | --- | ProteinName prefix |
| `peptide_qvalue` | --- | --- | --- | --- |
| `calculated_mz` | From the peptidoform (`Mass` when a modification has no known mass) | --- | Calculated M/Z | --- |
| `observed_mz` | m/z × (1 + Mass error [ppm] × 10⁻⁶); MaxQuant's `m/z` is theoretical | --- | --- | --- |
| `rt` | Retention time | RT | --- | RetentionTime |
| `rt_start` | --- | RT.Start | --- | --- |
| `rt_stop` | --- | RT.Stop | --- | --- |
| `predicted_rt` | --- | Predicted.RT | --- | --- |
| `intensities` | Intensity; Intensity L/M/H (SILAC) | Precursor.Quantity | Intensity | Intensity + Channel |
| `identified_labels` | Labeling State, and msms.txt Labeling state of its MS/MS (SILAC) | --- | --- | --- |
| `pg_accessions` | Proteins | Protein.Group (`DECOY_`-prefixed on decoy rows) | --- | ProteinName |
| `anchor_protein` | --- | First Protein.Group accession | --- | First ProteinName accession |
| `pg_positions` | --- | --- | --- | --- |
| `run_file_name` | Raw file | Run | --- | Reference / Run resolved through SDRF |

DIA-NN `Q.Value` is precursor-level and is stored in `additional_scores` as `precursor_qvalue`; it is not mapped to the peptide-level `peptide_qvalue` field.

## Protein Group Field Mappings

The Protein Group view captures protein-level quantification and inference results. The table below covers the three main supported tools.

| QPX Field | MaxQuant | DIA-NN | FragPipe |
| --- | --- | --- | --- |
| `pg_accessions` | Protein IDs | Protein.Group | Group + Indistinguishable Proteins |
| `pg_names` | Protein names | Protein.Names | --- |
| `gg_accessions` | Gene names | Genes | --- |
| `run_file_name` | combined | Run | --- |
| `peptide_counts` | Unique peptides [exp], Peptides [exp] | Unique.Stripped.Peptides | Unique Peptides |
| `feature_counts` | --- | Precursor.Quantity | Precursor Ions |
| `global_qvalue` | Q-value | Global.PG.Q.Value | --- |
| `intensities` | Intensity [exp]; Intensity L/M/H [exp] (SILAC); Reporter intensity | PG.Quantity | --- |
| `additional_intensities` | LFQ intensity, iBAQ, MS/MS count; Ratio X/Y [normalized] [exp] (SILAC) | PG.Normalised, PG.MaxLFQ | --- |
| `is_decoy` | Reverse | Decoy (`--report-decoys`; accessions prefixed `DECOY_`) | --- |
| `contaminant` | Potential contaminant | --- | --- |

!!! note "Tool version considerations"
    Column names may vary between tool versions. The mappings above reflect the commonly used versions at the time of writing:

    - **MaxQuant**: v2.x (`evidence.txt`, `proteinGroups.txt`)
    - **DIA-NN**: v1.8+ (`report.tsv`, `pg_matrix.tsv`)
    - **FragPipe**: v20+ (combined output files)

!!! tip "Adding support for a new tool"
    To add a converter for a new tool, implement a mapping from the tool's output columns to the QPX fields listed in the tables above. See the [PSM](psm.md), [Feature](feature.md), and [Protein Group](pg.md) schema pages for field type requirements and nullability constraints.
