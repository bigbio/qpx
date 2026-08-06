# Protein Group View

The protein group (PG) view is a tabular Parquet file that contains the details of protein groups identified and quantified per quantification unit. A quantification unit is the group of raw files aggregated together (`grouped_runs`) — a protein quantity only exists after aggregating peptides across a sample's fractions, so the view keys on this group of raw files rather than a single raw file (for unfractionated or DIA data the group is a single file). It captures the relationship between protein groups and the grouped runs in which they were detected, including peptide counts, feature counts, quality metrics, and intensity-based quantification. The sample is resolved downstream via `(any file in grouped_runs, label) -> run.samples[].sample_accession`.

This view is analogous to outputs from tools such as MaxQuant (`proteinGroups.txt`), DIA-NN (`pg_matrix`), and FragPipe protein group reports. Each row has a mandatory opaque `pg_id`, which is the primary key. A producer-supplied ID is preserved by default; otherwise QPX derives the ID from the footer-declared `identity_composite`.

## Use cases

- Retrieve all protein groups identified or quantified in a given raw file.
- Retrieve protein group abundance by file and condition.
- Store and query FDR q-values for protein groups at both the run and experiment level.
- Support downstream statistical analysis by providing per-file protein-level quantification.

## Schema

`pg_id` is the non-null primary key. The schema-default `identity_composite` is
`[pg_accessions, grouped_runs, label]`; `label` is null only for
identification-only protein groups that carry no quantity. Fields marked with
**(nullable)** may have null values. See the full YAML schema in
[`pg.yaml`](schemas/pg.yaml).

### Row identity

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `pg_id` | Opaque primary key, supplied by the producer or derived by QPX from the footer-declared `identity_composite` | `int64` | Yes |

### Identity

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `pg_accessions` | Protein accessions of all proteins within this group | `array[string]` | Yes |
| `pg_names` | Descriptive names for the proteins in the group | `array[string]` | No |
| `gg_accessions` | Gene group accessions as a string array | `array[string]` | No |
| `gg_names` | Gene names corresponding to the proteins in the group | `array[string]` | No |
| `gg_qvalue` | Gene group q-value (e.g., DIA-NN GG.Q.Value) | `float64`, null | No |
| `anchor_protein` | Representative protein of the group (leading protein); part of the default identity composite | `string` | Yes |
| `grouped_runs` | The group of raw files aggregated into one quantification unit (fractions aggregated together; single-element for unfractionated/DIA). The sample is resolved downstream via `(any file in grouped_runs, label) -> run.samples[].sample_accession`; part of the default identity composite | `array[string]` | Yes |

### Counts

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `peptide_counts` | Peptide sequence counts for this protein group in this file | `struct` | No |
| `peptide_counts.unique_sequences` | Number of peptide sequences unique to this protein group within this file | `int` | -- |
| `peptide_counts.total_sequences` | Total number of peptide sequences identified for this protein group in this file | `int` | -- |
| `feature_counts` | Peptide feature counts (peptide-charge combinations) for this protein group | `struct` | No |
| `feature_counts.unique_features` | Number of unique peptide features specific to this protein group within this file | `int` | -- |
| `feature_counts.total_features` | Total number of peptide features identified for this protein group in this file | `int` | -- |

### Quality

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `global_qvalue` | Global q-value of the protein group at the experiment level | `float64` | No |
| `pg_qvalue` | Protein group q-value at the run level | `float64` | No |
| `is_decoy` | Whether the protein group is a decoy (`true`) or a target (`false`) | `bool` | Yes |
| `contaminant` | Contaminant flag | `bool`, null | No |
| `sequence_coverage` | Percentage of the protein sequence covered by identified peptides | `float32` | No |
| `molecular_weight` | Molecular weight of the protein in kDa | `float32` | No |

### Quantification

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `label` | Channel/label of this quantification (e.g., TMT126, LFQ); **one row per label**. Null for identification-only groups. Part of the default identity composite. | `string` | Yes (nullable) |
| `intensity` | Primary raw intensity value for this label. Null for identification-only groups. | `float32` | No |
| `additional_intensities` | Pre-computed intensity values from the upstream tool (normalized, LFQ, iBAQ, etc.) for this row's label. See [Intensities](intensities.md) | `array[struct]` | No |
| `additional_scores` | Additional scores and metrics (posterior error probability, confidence, etc.). See [Scores](scores.md) | `array[struct]` | No |

Since QPX 1.1 the protein-group quantification is **flattened**: instead of an `intensities` list, each row carries a scalar `label` + `intensity`, so there is one row per `(pg_accessions, grouped_runs, label)`. Identification-only groups (no quantity) have null `label`/`intensity`.

Each entry in `additional_intensities` contains:

| Sub-field | Description | Type |
|-----------|-------------|------|
| `label` | Label identifier (e.g., TMT126, LFQ) | `string` |
| `intensities` | Array of name-value pairs for derived intensities | `array[struct{intensity_name, intensity_value}]` |

### Peptide detail

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `peptides` | Number of peptides per individual protein in the protein group | `array[struct]` | Yes |

Each entry in `peptides` contains:

| Sub-field | Description | Type |
|-----------|-------------|------|
| `protein_name` | Protein accession | `string` |
| `peptide_count` | Number of peptides for this protein | `int` |

### CV params

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `cv_params` | Optional list of controlled vocabulary parameters for additional metadata | `array[struct{cv_name, cv_value}]` | No |

## Example

```json
{
  "pg_id": 2937018471546916399,
  "pg_accessions": ["P04217", "A0A024R4E5"],
  "pg_names": ["Alpha-1B-glycoprotein", "Alpha-1B-glycoprotein variant"],
  "gg_accessions": ["A1BG"],
  "gg_names": ["A1BG"],
  "anchor_protein": "P04217",
  "grouped_runs": ["20230101_sample_01"],
  "peptide_counts": {
    "unique_sequences": 12,
    "total_sequences": 18
  },
  "feature_counts": {
    "unique_features": 24,
    "total_features": 36
  },
  "global_qvalue": 0.001,
  "pg_qvalue": 0.005,
  "is_decoy": false,
  "contaminant": false,
  "sequence_coverage": 45.2,
  "molecular_weight": 54.3,
  "label": "TMT126",
  "intensity": 1.5e8,
  "additional_intensities": [
    {
      "label": "TMT126",
      "intensities": [
        {"intensity_name": "LFQ", "intensity_value": 1.2e8},
        {"intensity_name": "iBAQ", "intensity_value": 3.4e7}
      ]
    }
  ],
  "peptides": [
    {"protein_name": "P04217", "peptide_count": 15},
    {"protein_name": "A0A024R4E5", "peptide_count": 3}
  ],
  "additional_scores": [
    {"score_name": "posterior_error_probability", "score_value": 0.0001, "higher_better": false}
  ],
  "cv_params": null
}
```

## Recommended Intensity Names

The following names are commonly used in `additional_intensities` for the PG view. Use consistent naming so downstream tools can recognise them across datasets.

| `intensity_name` | Source tool(s) | Description |
|---|---|---|
| `maxlfq` | DIA-NN, MaxQuant, FragPipe | MaxLFQ normalised protein group quantity |
| `lfq` | MaxQuant | Label-free quantification intensity |
| `ibaq` | MaxQuant | Intensity-based absolute quantification |
| `topn` | DIA-NN | Top-N normalised protein group quantity |
| `normalize_intensity` | quantms | Median-normalised intensity |
| `spectral_count` | FragPipe, MaxQuant | Number of PSMs for razor peptides |
| `unique_spectral_count` | FragPipe | Number of PSMs for unique peptides only |
| `total_spectral_count` | FragPipe | Number of PSMs for all peptides (razor + shared) |
| `genes_maxlfq` | DIA-NN | Gene-group level MaxLFQ quantity |
| `genes_maxlfq_unique` | DIA-NN | Gene-group MaxLFQ using only proteotypic peptides |
| `ratio_h_l` | MaxQuant (SILAC) | Heavy-to-light protein ratio |
| `ratio_h_l_normalized` | MaxQuant (SILAC) | Normalised heavy-to-light ratio |
| `reporter_intensity` | MaxQuant (TMT/iTRAQ) | Raw reporter-ion intensity |
| `reporter_intensity_corrected` | MaxQuant (TMT/iTRAQ) | Isotope-purity-corrected reporter intensity |

!!! tip "Spectral counts as intensities"
    Spectral counts are quantitative measures that vary per run and per label, so they fit naturally in `additional_intensities`. Store them as float values (e.g., `3.0` instead of `3`) since the struct uses `float32`.

!!! tip "Gene-group quantities"
    DIA-NN reports protein-group level (`PG.MaxLFQ`) and gene-group level (`Genes.MaxLFQ`) quantities separately. Store gene-group quantities in `additional_intensities` with names prefixed `genes_` (e.g., `genes_maxlfq`, `genes_maxlfq_unique`). The `gg_accessions` field identifies the gene group.

## Recommended Score Names

The following score names are commonly used in `additional_scores` for the PG view. For the full list of recommended score names and naming conventions, see [Scores](scores.md).

| `score_name` | Source tool(s) | Direction | Description |
|---|---|---|---|
| `posterior_error_probability` | MaxQuant | lower is better | Protein-group PEP |
| `andromeda_score` | MaxQuant | higher is better | Andromeda protein-level score |
| `protein_probability` | FragPipe | higher is better | ProteinProphet probability |
| `top_peptide_probability` | FragPipe | higher is better | Highest PeptideProphet probability among peptides |
| `pg_maxlfq_quality` | DIA-NN | higher is better | QuantUMS quality for PG.MaxLFQ estimate |
| `genes_maxlfq_quality` | DIA-NN | higher is better | QuantUMS quality for gene-level MaxLFQ |
| `pg_pep` | DIA-NN | lower is better | Posterior error probability at the protein group level |
| `gg_qvalue` | DIA-NN | lower is better | Gene group q-value |
| `lib_pg_qvalue` | DIA-NN | lower is better | Library protein group q-value |
| `protein_qvalue` | DIA-NN | lower is better | Unique protein q-value (proteotypic evidence) |

## Tool Mappings

This section shows how output columns from common search engines and pipelines map to `pg.parquet` fields.

!!! info "Wide-to-long conversion"
    MaxQuant, DIA-NN, and FragPipe output protein groups in **wide format** (one row per protein group, one column per experiment/run). QPX uses **long format** (one row per protein group **per quantification unit**). Converters must melt wide-format columns into separate rows keyed by `grouped_runs` (the group of raw files aggregated into one quantification unit; single-element for unfractionated/DIA).

### MaxQuant (`proteinGroups.txt`)

=== "Identity & Metadata"

    | MaxQuant column | QPX field | Notes |
    |---|---|---|
    | `Protein IDs` | `pg_accessions` | Semicolon-separated → array |
    | `Majority protein IDs` | `anchor_protein` | First accession of majority IDs |
    | `Protein names` | `pg_names` | Semicolon-separated → array |
    | `Gene names` | `gg_names` | Semicolon-separated → array |
    | `Sequence coverage [%]` | `sequence_coverage` | |
    | `Mol. weight [kDa]` | `molecular_weight` | |
    | `Reverse` | `is_decoy` | `"+"` → `true`, else `false` |
    | `Potential contaminant` | `contaminant` | `"+"` → `true`, else `false` |
    | `Only identified by site` | `cv_params` | `{cv_name: "only_identified_by_site", cv_value: "true"}` |
    | `Identification type [exp]` | `cv_params` | `{cv_name: "identification_type", cv_value: "By MS/MS"}` |

=== "Counts"

    | MaxQuant column | QPX field | Notes |
    |---|---|---|
    | `Unique peptides` | `peptide_counts.unique_sequences` | Peptides exclusive to this PG |
    | `Peptides` or `Razor + unique peptides` | `peptide_counts.total_sequences` | Total assigned peptides |
    | `Peptide counts (all)` | `peptides` | Per-protein peptide counts |

=== "Quantification"

    | MaxQuant column | QPX field | Notes |
    |---|---|---|
    | `Intensity [exp]` | `intensity` | Primary raw intensity; scalar per row (one row per label) |
    | `LFQ intensity [exp]` | `additional_intensities` → `lfq` | MaxLFQ normalised |
    | `iBAQ [exp]` | `additional_intensities` → `ibaq` | Intensity-based absolute quantification |
    | `MS/MS count [exp]` | `additional_intensities` → `spectral_count` | PSM count as float |
    | `Reporter intensity [channel]` | `intensity` | For TMT/iTRAQ, one **row** per channel (`label` = channel) |
    | `Reporter intensity corrected [channel]` | `additional_intensities` → `reporter_intensity_corrected` | |
    | `Ratio H/L [exp]` | `additional_intensities` → `ratio_h_l` | SILAC ratio |
    | `Ratio H/L normalized [exp]` | `additional_intensities` → `ratio_h_l_normalized` | Normalised SILAC ratio |

=== "Quality"

    | MaxQuant column | QPX field | Notes |
    |---|---|---|
    | `Q-value` | `global_qvalue` | Protein group FDR |
    | `Score` | `additional_scores` → `andromeda_score` | |
    | `PEP` | `additional_scores` → `posterior_error_probability` | |

### DIA-NN (main report)

DIA-NN's main report is precursor-level. PG-level columns (`PG.*`, `Genes.*`) are repeated for every precursor in the group within a run. Deduplicate by `Protein.Group` + `Run` to produce one `pg.parquet` row.

=== "Identity"

    | DIA-NN column | QPX field | Notes |
    |---|---|---|
    | `Protein.Group` | `pg_accessions` | Semicolon-separated → array; also `anchor_protein` (first accession) |
    | `Protein.Names` | `pg_names` | |
    | `Genes` | `gg_names` / `gg_accessions` | |
    | `Run` | `grouped_runs` | Raw file name without path, wrapped as a single-element list |

=== "Quantification"

    | DIA-NN column | QPX field | Notes |
    |---|---|---|
    | `PG.Quantity` | `intensity` | Non-normalised protein group quantity (scalar; DIA `label` = `LFQ`) |
    | `PG.MaxLFQ` | `additional_intensities` → `maxlfq` | QuantUMS/MaxLFQ normalised |
    | `PG.Normalised` | `additional_intensities` → `normalize_intensity` | Normalised PG quantity |
    | `PG.TopN` | `additional_intensities` → `topn` | Top-N normalised quantity |
    | `Genes.MaxLFQ` | `additional_intensities` → `genes_maxlfq` | Gene-group MaxLFQ |
    | `Genes.MaxLFQ.Unique` | `additional_intensities` → `genes_maxlfq_unique` | Gene-group MaxLFQ (proteotypic only) |
    | `Genes.TopN` | `additional_intensities` → `genes_topn` | Gene-group Top-N |

=== "Quality"

    | DIA-NN column | QPX field | Notes |
    |---|---|---|
    | `Global.PG.Q.Value` | `global_qvalue` | Experiment-level PG FDR |
    | `PG.Q.Value` | `pg_qvalue` | Run-level PG FDR |
    | `PG.PEP` | `additional_scores` → `pg_pep` | |
    | `PG.MaxLFQ.Quality` | `additional_scores` → `pg_maxlfq_quality` | QuantUMS quality metric |
    | `GG.Q.Value` | `additional_scores` → `gg_qvalue` | Gene group q-value |
    | `Lib.PG.Q.Value` | `additional_scores` → `lib_pg_qvalue` | Library PG q-value |
    | `Protein.Q.Value` | `additional_scores` → `protein_qvalue` | Unique protein q-value |
    | `Genes.MaxLFQ.Quality` | `additional_scores` → `genes_maxlfq_quality` | |

### FragPipe (`combined_protein.tsv`)

FragPipe outputs are pre-filtered (no decoys or contaminants). Per-experiment columns are prefixed with the experiment name.

=== "Identity & Metadata"

    | FragPipe column | QPX field | Notes |
    |---|---|---|
    | `Protein ID` | `anchor_protein` | Leading protein accession |
    | `Protein ID` + `Indistinguishable Proteins` | `pg_accessions` | Combine into array |
    | `Description` | `pg_names` | Protein description |
    | `Gene` | `gg_names` / `gg_accessions` | |
    | `Coverage` | `sequence_coverage` | Percentage (0–100) |
    | `Protein Length` | `cv_params` | `{cv_name: "sequence_length", cv_value: "495"}` |
    | `Protein Existence` | `cv_params` | `{cv_name: "protein_existence", cv_value: "1"}` |

=== "Counts"

    | FragPipe column | QPX field | Notes |
    |---|---|---|
    | `[exp] Total Peptides` | `peptide_counts.total_sequences` | Per-experiment peptide count |

=== "Quantification"

    | FragPipe column | QPX field | Notes |
    |---|---|---|
    | `[exp] Intensity` | `intensity` | Primary razor peptide intensity (scalar per row) |
    | `[exp] MaxLFQ Intensity` | `additional_intensities` → `maxlfq` | MaxLFQ normalised |
    | `[exp] MaxLFQ Unique Intensity` | `additional_intensities` → `maxlfq_unique` | MaxLFQ using only unique peptides |
    | `[exp] Unique Intensity` | `additional_intensities` → `unique_intensity` | Intensity from unique peptides only |
    | `[exp] Spectral Count` | `additional_intensities` → `spectral_count` | Razor peptide PSMs |
    | `[exp] Unique Spectral Count` | `additional_intensities` → `unique_spectral_count` | Unique peptide PSMs |
    | `[exp] Total Spectral Count` | `additional_intensities` → `total_spectral_count` | All peptide PSMs |

=== "Quality"

    | FragPipe column | QPX field | Notes |
    |---|---|---|
    | `Protein Probability` | `additional_scores` → `protein_probability` | ProteinProphet score |
    | `Top Peptide Probability` | `additional_scores` → `top_peptide_probability` | Best PeptideProphet score |

!!! note "FragPipe q-values"
    FragPipe does not report explicit protein-group q-values. The `Protein Probability` (from ProteinProphet) is filtered at a given FDR threshold (typically 1%) before writing `combined_protein.tsv`. Set `global_qvalue` to null and store the probability in `additional_scores`.

## Notes

!!! note "Relationship to other views"
    The PG view provides per-file protein group quantification. For derived per-sample summaries (protein counts, abundances, etc.), see [API Views](views.md). For downstream absolute or differential expression results, see the [Absolute Expression](absolute.md) and [Differential Expression](differential.md) views.

!!! warning "Identity constraints"
    `pg_id` is the primary key. When QPX derives it, the default identity composite is `anchor_protein`, `grouped_runs`, and `label`. `anchor_protein` and `grouped_runs` MUST NOT be null. `label` is null **only** for identification-only protein groups that carry no quantity (e.g. mzIdentML); when a quantity exists, `label` is non-null and there is one row per label. `grouped_runs` is a **set of distinct raw files**: it MUST NOT contain duplicates, and two composites are equal when they contain the same files regardless of order (identity is compared set-wise). The list is **stored in fraction order** — sorting is not applied, because it would destroy fraction ordering (and lexicographically misorder names like `F1, F10, F2`). Each record represents a single protein group quantified in one quantification unit (the group of raw files/fractions aggregated together) for one label. A protein quantity only exists after aggregating peptides across a sample's fractions, so the PG view identifies this group of raw files rather than a single raw file; for unfractionated or DIA data the list has a single element.

    Within one `(anchor_protein, label)`, the `grouped_runs` sets across rows MUST be **disjoint** — every raw file contributes to at most one row, so no measurement is counted twice (the *run-disjointness* invariant enforced by validation).
