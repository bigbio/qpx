# Feature View

The feature view captures quantified peptide information, including intensity
across labels and protein-group mappings. Each row has a mandatory opaque
`feature_id`, which is the primary key. For identified Features, QPX derives the
ID from the upstream properties recorded in the Parquet footer's
`identity_composite`; any producer ID remains traceable in `cv_params`. For
unidentified Features, a producer-supplied ID is required and QPX namespaces it
by its run or quantification unit; the original ID remains traceable in
`cv_params`.
Converters may override the schema default when their producer represents a
different Feature-level entity.

When enriching an existing native OpenMS QPX bundle, QPX preserves supplied
opaque IDs and their cross-view references. Run references are normalized to
the SDRF basename convention before any missing IDs are derived; supplied IDs
are not recalculated from the normalized run names.

For a de novo workflow without a database search, set `is_decoy` to `false` and
record a `de_novo_peptide_sequencing` step in `provenance.parquet`. Protein-mapping
fields may be null; `anchor_protein` is an annotation and is not part of a
feature's identity. The schema default composite is `[peptidoform, charge,
run_file_name, rt]`. FragPipe `combined_ion`, which has no per-feature RT and may
separate the same precursor by FAIMS voltage, instead declares
`[quantification_unit_id, peptidoform, charge, compensation_voltage]`. Identity
OpenMS consensusXML declares `[peptidoform, charge, run_file_name, rt, scan,
observed_mz, consensus_rt]`, where `scan` contains every spectrum reference for
the run. Identity is meaningful within a file only; distinct QPX files must not
be joined on `feature_id` alone ([Row Identity](identity.md)).

## Use Cases

- **Quantified peptide information**: Stores peptide intensities linked to sample metadata, enabling downstream quantitative analysis and integration with SDRF annotations.
- **Peptide-level statistics**: Enables algorithms that operate at the peptide level (e.g., peptide-to-protein rollup, normalization, missing value imputation).
- **Integration with sample metadata**: Each feature carries label information, connecting quantification data to experimental design described in the SDRF.

## Schema

### Identity

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `feature_id` | Opaque primary key, supplied by the producer or derived by QPX from the footer-declared `identity_composite` | int64 | yes |

### Core Identification Fields

These fields are shared with the PSM view and describe the peptide identification associated with the feature.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `sequence` | Unmodified peptide amino acid sequence | string | yes |
| `peptidoform` | Peptide sequence with modifications in ProForma notation | string | yes |
| `modifications` | Structured list of modifications with name, accession, position, and localization scores | array[struct], null | no |
| `modification_candidates` | Competing modification interpretations per position, combined over the feature's identifications by the producer (see [Competing interpretations](modifications.md#competing-interpretations)); optional column | array[struct], null | optional |
| `charge` | Charge of the quantified analyte | int16 | yes |
| `posterior_error_probability` | Posterior error probability (PEP) for the peptide match | float64, null | no |
| `peptide_qvalue` | Peptide-level q-value; null when the search engine does not provide one | float64, null | optional |
| `precursor_qvalue` | Run-level precursor q-value, filled from `additional_scores` ([typed q-values](scores.md#typed-q-value-columns)); specification 1.3 | float64, null | No |
| `global_qvalue` | Experiment-level precursor q-value, filled from `additional_scores`; specification 1.3 | float64, null | No |
| `is_decoy` | Whether the peptide is a decoy match (`true`) or a target match (`false`); use `false` when no target-decoy search was used. A peptide that matches both target and decoy proteins (OpenMS `target_decoy = target+decoy`) is a target match | bool | yes |
| `calculated_mz` | Theoretical peptide mass-to-charge ratio based on identified sequence and modifications; null when it cannot be calculated | float32, null | yes |
| `observed_mz` | Experimental observed peptide mass-to-charge ratio; null when the producer reports no measured precursor m/z | float32, null | yes |
| `mass_error_ppm` | Mass error in ppm: 1e6 × (observed_mz − calculated_mz) / calculated_mz | float32, null | no |
| `missed_cleavages` | Number of missed enzymatic cleavages | int16, null | no |
| `rt` | Precursor retention time (in seconds) | float32, null | no |
| `rt_aligned` | Retention time after the tool's alignment between runs, in seconds; `rt` is the time on this run's own clock. MaxQuant `Calibrated retention time`; specification 1.3 | float32, null | No |
| `irt` | Indexed retention time on the tool's scale, filled from the `irt` or `empirical_irt` entry of `cv_params`; specification 1.3 | float32, null | No |
| `consensus_rt` | Parent OpenMS `ConsensusFeature` retention time | float32, null | optional |
| `rt_start` | Start of the retention time window for the feature | float32, null | no |
| `rt_stop` | End of the retention time window for the feature | float32, null | no |
| `predicted_rt` | Predicted retention time of the peptide (in seconds) | float32, null | no |
| `ion_mobility` | Ion mobility value for the precursor ion | float32, null | no |
| `ion_mobility_start` | Start ion mobility value for the precursor ion | float32, null | no |
| `ion_mobility_stop` | Stop ion mobility value for the precursor ion | float32, null | no |
| `compensation_voltage` | FAIMS compensation voltage; `0` for non-FAIMS FragPipe units | float32, null | optional |
| `additional_scores` | List of score structures with name, value, and direction indicator | array[struct], null | no |
| `cv_params` | Optional list of controlled vocabulary parameters for additional metadata | array[struct], null | no |

### Quantification Fields

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `intensities` | Primary intensity-based abundance of the feature across labels | array[struct], null | no |
| `additional_intensities` | Pre-computed intensity values from the upstream tool (e.g., normalized, LFQ, iBAQ) as named key-value pairs per label | array[struct], null | no |
| `channel_details` | MS1-label multiplets: `{label, mz, rt, qvalue, evidence}` of each partner ([MS1 labels](intensities.md#ms1-labels-silac-dimethyl-mtraq)); specification 1.3 | array[struct], null | optional |
| `representative_label` | The channel whose report row gave the row's scalars, where a tool reports one row per channel; specification 1.3 | string, null | optional |
| `identified_labels` | MS1-label multiplets: the labels whose partner was identified by MS/MS ([MS1 labels](intensities.md#ms1-labels-silac-dimethyl-mtraq)) | array[string], null | optional |
| `run_file_name` | The run file name that contains the feature | string | yes |
| `quantification_unit_id` | Quantification unit represented by the quantity when it may aggregate multiple raw runs | string, null | optional |

!!! tip "Intensity structure"
    For details on the `intensities` and `additional_intensities` data structures, including examples for LFQ and TMT experiments, see [Intensities](intensities.md).

### Protein Group Fields

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `pg_accessions` | Full membership of the peptide's protein group; for a peptide shared by several groups, the group the producer assigned it to (e.g. MaxQuant's razor group) | array[struct], null | no |
| `anchor_protein` | One protein accession that represents the protein group; null when no protein mapping was performed (converters normalize blank/whitespace to null) | string, null | yes |
| `unique` | True when the peptide's protein group has exactly one protein and the peptide maps to no other group; false when the group has several members or the peptide is shared with other groups; null when the group is unresolved | bool, null | no |
| `pg_positions` | Peptide start and end positions within each protein in the protein group | array[struct], null | optional |
| `pg_group_id` | Hash of the sorted `pg_accessions`, equal to `pg.pg_group_id`, to group and join without comparing lists; specification 1.3 | int64, null | No |
| `pg_global_qvalue` | Global q-value of the protein group at the experiment level | float64, null | optional |
| `gg_accessions` | Gene accessions associated with the protein group | array[string], null | optional |
| `gg_names` | Gene names associated with the protein group | array[string], null | optional |

Each entry in `pg_positions` contains:

| Sub-field | Description | Type |
|-----------|-------------|------|
| `protein_accession` | Protein accession within the protein group | `string` |
| `start` | 1-based start position of the peptide in the protein sequence | `int` |
| `end` | 1-based end position of the peptide in the protein sequence (inclusive) | `int` |

!!! note "Gene and protein inference data"
    Gene accessions, gene names, and unique peptide indicators are optionally included in the feature file for convenience. Protein-level scores are stored in the [Protein Group View](pg.md). When `pg_ids` is populated, join it directly to `pg.pg_id`. Otherwise join on the group membership (`pg_accessions` compared as a set with `pg.pg_accessions`), require `feature.run_file_name` to be a member of `pg.grouped_runs`, and match `pg.label` to a label in `feature.intensities` (`Dataset.link_feature_pg()`).

!!! info "Optional vs nullable"
    `pg_global_qvalue` and `peptide_qvalue` are **optional** — the column may be absent from the file entirely if the search engine does not provide that q-value. When present, individual values may be null.

!!! warning "Which q-value to read"
    `peptide_qvalue` is the canonical peptide-level q-value and is the field to read when filtering or reporting peptide-level FDR. DIA-NN does not report a peptide-level q-value: its precursor-level `Q.Value` remains in `additional_scores` under `precursor_qvalue`, and the DIA-NN converter leaves `peptide_qvalue` null.

    Other producer-specific q-values may also appear in `additional_scores`; use the documented score name and level rather than treating every q-value as peptide-level confidence.

    Percolator with peptide-level FDR scores each peptide once: its best PSM carries the peptide's q-value and PEP, every other PSM gets 1.0. The OpenMS converters report the best PSM's values and never the 1.0 placeholders; `peptide_qvalue` is null when the best PSM is not in the input.

    A value belongs in `peptide_qvalue` **only** when the producer states it is a q-value. Writing a raw search-engine score (a Percolator SVM score, an E-value) into this column makes it indistinguishable from an FDR downstream — leave the field null instead.

### Spectra Reference

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `id_run_file_name` | The run file containing the best PSM that identified the feature (may differ from `run_file_name`) | string, null | no |
| `scan` | Scan identifier of the best PSM that identified the feature, as an array of integer components (e.g., `[43920]` for single-scan instruments, `[10, 1, 345]` for Waters function/process/scan) | array[int32] | yes |

### Identification Origin

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `identification_type` | How the tool identified the feature in this run: `direct` (its own MS/MS), `transferred` (matched from another run), `secondary` (a co-fragmenting peptide of a chimeric spectrum). Null when the tool does not say; `id_run_file_name` is null for a transferred row and can be null for other reasons, so filter on this field. | string, null | No |
| `match_score`, `match_qvalue` | The tool's score and q-value for the match of a transferred identification | float, null | No |
| `match_rt_delta`, `match_mz_delta` | Difference of a transferred identification from its donor, in seconds and Th | float, null | No |
| `imputed_labels` | Labels of `intensities` entries whose value the tool imputed instead of measuring. Null when none are or the tool does not say. | array[string], null | No |

### Isobaric Quantification Scan

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `parent_ion_fraction` | Fraction of the isolation window's signal from the target precursor, 0-1. Isolation interference is one minus this. | float, null | No |
| `quant_scan` | Scan the reporter intensities were read from, as the integer components `scan` uses. Differs from `scan` under SPS-MS3. | array[int32], null | No |
| `quant_ms_level` | MS level of that scan, 2 or 3, when the producer records it | int32, null | No |

### Optional Cross-References

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `psm_ids` | The PSM rows this feature maps to (through `psm.psm_id`). A **computed inverse softlink**, not materialized by qpx: the authoritative direction is `psm.feature_id`, and `Dataset.link_feature_psm()` recovers this list by grouping `(feature_id, psm_id)` on `feature_id` ([bigbio/qpx#267](https://github.com/bigbio/qpx/issues/267)). Remains an optional producer hardlink a producer MAY populate; null otherwise. | array[int64], null | optional |
| `pg_ids` | References to protein-group rows through `pg.pg_id`; null when no explicit links are available | array[int64], null | optional |

## Shared Fields

Several fields in the feature view use structures shared across other QPX views:

- For details on the `modifications` field structure, see [Modifications](modifications.md).
- For details on `intensities` and `additional_intensities`, see [Intensities](intensities.md).
- For details on `additional_scores` and score semantics, see [Scores](scores.md).
- For details on `cv_params` usage and recommended terms, see [Scores & CV Terms](scores.md).

## Example

### Feature with TMT Intensities

```json
{
  "feature_id": 8049980460766804816,
  "sequence": "AADLLTSFLGHK",
  "peptidoform": "AADLLTSFLGHK",
  "modifications": null,
  "charge": 3,
  "is_decoy": false,
  "calculated_mz": 424.2345,
  "observed_mz": 424.2350,
  "rt": 2345.67,
  "rt_start": 2340.12,
  "rt_stop": 2351.23,
  "predicted_rt": 2348.00,
  "ion_mobility": null,
  "ion_mobility_start": null,
  "ion_mobility_stop": null,
  "run_file_name": "20200101_TMT_fraction01",
  "scan": [15234],
  "id_run_file_name": "20200101_TMT_fraction01",
  "pg_accessions": ["P04217", "P04217-2"],
  "anchor_protein": "P04217",
  "pg_positions": [
    {"protein_accession": "P04217", "start": 25, "end": 36},
    {"protein_accession": "P04217-2", "start": 25, "end": 36}
  ],
  "intensities": [
    {"label": "TMT126", "intensity": 15234.5},
    {"label": "TMT127N", "intensity": 18456.7},
    {"label": "TMT127C", "intensity": 12890.3},
    {"label": "TMT128N", "intensity": 21045.8}
  ],
  "additional_intensities": [
    {
      "label": "TMT126",
      "intensities": [
        {"intensity_name": "normalize_intensity", "intensity_value": 0.1523},
        {"intensity_name": "ibaq", "intensity_value": 4567.8}
      ]
    },
    {
      "label": "TMT127N",
      "intensities": [
        {"intensity_name": "normalize_intensity", "intensity_value": 0.1846},
        {"intensity_name": "ibaq", "intensity_value": 5432.1}
      ]
    }
  ],
  "additional_scores": [
    {
      "score_name": "global_qvalue",
      "score_value": 0.0012,
      "higher_better": false
    }
  ],
  "cv_params": null
}
```

### Feature with LFQ Intensity

```json
{
  "feature_id": -5948337580339590069,
  "sequence": "VLHPLEGAVVIIFK",
  "peptidoform": "[UNIMOD:1]-VLHPLEGAVVIIFK",
  "modifications": [
    {
      "name": "Acetyl",
      "accession": "UNIMOD:1",
      "positions": [
        {"position": 0, "amino_acid": null, "scores": null}
      ]
    }
  ],
  "charge": 2,
  "is_decoy": false,
  "calculated_mz": 782.4721,
  "observed_mz": 782.4725,
  "rt": 3567.89,
  "run_file_name": "20200101_LFQ_rep1",
  "scan": [28901],
  "anchor_protein": "P68871",
  "pg_accessions": ["P68871"],
  "pg_positions": [
    {"protein_accession": "P68871", "start": 2, "end": 15}
  ],
  "intensities": [
    {"label": "LFQ", "intensity": 98765.4}
  ],
  "additional_intensities": null
}
```

## File Metadata

Feature Parquet files store file-level metadata as key-value pairs in the Parquet footer. The following metadata fields are defined:

| Field | Description |
|-------|-------------|
| `qpx_version` | Version of the QPX format used to generate the file |
| `software_provider` | Software that wrote the file (`qpx` for the qpx writers) |
| `scan_format` | Format of scan identifiers: `scan`, `index`, or `nativeId` |
| `creator` | Tool whose output the data was converted from (e.g. `maxquant`, `diann`) |
| `file_type` | Type of the file (value: `feature_file`) |
| `creation_date` | Date when the file was created |
| `compression_format` | Compression algorithm used: `zstd` (default), `snappy`, `gzip`, `lzo`, or `none` |

!!! tip "Reading file metadata in Python"
    ```python
    import pyarrow.parquet as pq

    parquet_file = pq.ParquetFile("experiment.feature.parquet")
    metadata = parquet_file.schema_arrow.metadata
    for key, value in metadata.items():
        print(f"{key.decode()}: {value.decode()}")
    ```

## Notes

!!! note "Works for both DDA and DIA"
    Unlike the PSM view (which is DDA-specific), the feature view supports both DDA and DIA workflows. For DIA experiments, the feature view is the primary peptide-level output.

- **Relationship to PSM view**: A feature aggregates one or more PSMs into a single quantified peptide entry. The `scan` and `id_run_file_name` fields link back to the best PSM that identified the feature. In DDA-LFQ workflows, both views may exist; in DIA workflows, only the feature view is typically generated.
- **Relationship to protein group view**: The feature view contains protein group mappings (`pg_accessions`, `anchor_protein`) and peptide positions (`pg_positions`) that connect each peptide to its inferred protein groups. Gene annotations, unique peptide indicators, and protein-level scores are stored in the [Protein Group View](pg.md). The PG view provides the complete protein-level perspective with aggregated intensities and peptide counts.
- **Intensities vs additional_intensities**: Use `intensities` for raw/primary measurements across experimental labels (TMT/iTRAQ tags or LFQ). Use `additional_intensities` for values pre-computed by the upstream tool (e.g., normalized intensities, LFQ, iBAQ). QPX reads these from the tool's output — it does not compute them. This separation keeps experimental design aspects distinct from data processing aspects.
- **Protein group accessions**: The `pg_accessions` field should contain all proteins within a protein group. The `anchor_protein` is the representative protein selected by the search engine or inference algorithm to represent the group. The `pg_positions` field maps the peptide's start and end coordinates within each protein in the group.
