# PSM View

The PSM (Peptide Spectrum Match) view captures spectrum-level identification
results. Each row represents a single match between a mass spectrum and a
peptide sequence, including the identification scores, optional spectral data,
and protein mappings. Its mandatory opaque `psm_id` is the primary key. A
producer-supplied ID is preserved by default; otherwise QPX derives the ID from
the footer-declared `identity_composite`, whose schema default is
`[peptidoform, charge, run_file_name, scan]`.

De novo sequencing workflows set `is_decoy` to `false` and may omit
`protein_accessions`. Record the workflow with a `de_novo_peptide_sequencing`
step in `provenance.parquet`.

## Use Cases

- **AI/ML training**: Provides peptide-spectrum pairs with optional spectral arrays (m/z, intensity, ion types) for training intensity prediction, de novo sequencing, and clustering models.
- **Spectrum-level analysis**: Enables detailed inspection of individual identifications, including retention time, charge state, and search engine scores.
- **DDA identification results**: Designed primarily for data-dependent acquisition (DDA) workflows where each spectrum yields one or more peptide identifications.

## Schema

### Identity

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `psm_id` | Opaque primary key, supplied by the producer or derived by QPX from the footer-declared `identity_composite` | int64 | yes |

### Core Identification Fields

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `sequence` | Unmodified peptide amino acid sequence | string | yes |
| `peptidoform` | Peptide sequence with modifications in ProForma notation | string | yes |
| `modifications` | Structured list of modifications with name, accession, position, and localization scores | array[struct], null | no |
| `charge` | Charge state of the precursor ion | int16 | yes |
| `posterior_error_probability` | Posterior error probability (PEP) for the peptide-spectrum match — the probability that the PSM is incorrect. **Lower values indicate higher confidence** (lower is better). Ranges from 0.0 (confident) to 1.0 (likely incorrect) | float64, null | no |
| `is_decoy` | Whether the PSM is a decoy match (`true`) or a target match (`false`); use `false` when no target-decoy search was used | bool | yes |
| `calculated_mz` | Theoretical peptide mass-to-charge ratio based on identified sequence and modifications; null when it cannot be calculated | float32, null | yes |
| `observed_mz` | Experimental observed peptide mass-to-charge ratio | float32 | yes |
| `mass_error_ppm` | Mass error in ppm: 1e6 × (observed_mz − calculated_mz) / calculated_mz | float32, null | no |
| `missed_cleavages` | Number of missed enzymatic cleavages | int16, null | no |
| `rt` | MS2 scan's retention time (in seconds) | float32, null | no |
| `predicted_rt` | Predicted retention time of the peptide (in seconds) | float32, null | no |
| `run_file_name` | Spectrum file name without path or extension | string | yes |
| `scan` | Scan identifier as an array of integer components (e.g., `[43920]` for single-scan instruments, `[10, 1, 345]` for Waters function/process/scan) | array[int32] | yes |
| `additional_scores` | List of score structures with name, value, and direction indicator | array[struct], null | no |
| `cv_params` | Optional list of controlled vocabulary parameters for additional metadata | array[struct], null | no |

### Optional Fields

These fields are optional and may not exist in the file at all. They are included based on conversion settings or user preference.

| Field | Description                                                                                                                                             | Type | Required |
|-------|---------------------------------------------------------------------------------------------------------------------------------------------------------|------|----------|
| `protein_accessions` | Protein accessions of all proteins that the peptide maps to. Optional because protein mapping can be recovered from the feature and protein group views | array[string], null | no |
| `cross_links` | Cross-link information for XL-MS experiments. Each entry describes one cross-link site. `null` for non-cross-linked PSMs                                | array[struct], null | no |
| `ion_mobility` | Ion mobility value for the precursor ion                                                                                                                | float32, null | no |
| `mz_array` | Array of m/z values for the spectrum                                                                                                                    | array[float32], null | no |
| `intensity_array` | Array of intensity values for the spectrum                                                                                                              | array[float32], null | no |
| `charge_array` | Array of fragment ion charge values                                                                                                                     | array[int32], null | no |
| `ion_type_array` | Array of fragment ion type annotations (e.g., b1, y2, a2)                                                                                               | array[string], null | no |
| `ion_mobility_array` | Array of fragment ion mobility values                                                                                                                   | array[float32], null | no |
| `feature_id` | **Authoritative** (optional) foreign key to the linked Feature row through `feature.feature_id`; the producer's assignment of this PSM to a consensus feature (null when unlinked). It is NOT recomputable from shared identification fields, so it is materialized. Its inverse, `feature.psm_ids`, is a **computed softlink** — `Dataset.link_feature_psm()` groups `(feature_id, psm_id)` by `feature_id`; qpx does not materialize it ([bigbio/qpx#267](https://github.com/bigbio/qpx/issues/267)). | int64, null | no |

!!! note "Nullable vs Optional"
    Core fields marked as "not required" are **nullable** -- the column always exists in the file but individual values may be null. Optional fields (protein accessions, spectral data) may be **absent from the file entirely**, depending on conversion settings. Protein mappings can be recovered by joining with the feature and protein group views.

### Cross-Linking Fields

The `cross_links` field supports XL-MS (cross-linking mass spectrometry) experiments following the mzIdentML 1.3 specification. Each entry in the array represents one cross-link site.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `xl_type` | Type of cross-link: `"inter"` (between two peptides), `"intra"` (within same peptide), or `"dead-end"` (one reactive end) | string | yes |
| `partner_sequence` | Plain amino acid sequence of the beta (partner) peptide. `null` for dead-end and intra-peptide links | string, null | no |
| `partner_peptidoform` | Beta peptide in ProForma notation. `null` for dead-end and intra-peptide links | string, null | no |
| `donor_position` | Cross-link attachment position on the alpha peptide (1-indexed) | int32 | yes |
| `acceptor_position` | Attachment position on the beta peptide (1-indexed). `null` for dead-end links | int32, null | no |
| `linker_name` | Name of the cross-linker reagent (e.g., `"DSS"`, `"BS3"`, `"DSSO"`) | string | yes |
| `linker_accession` | XLMOD controlled vocabulary accession for the linker (e.g., `"XLMOD:02001"`) | string, null | no |
| `linker_mass` | Cross-linker mass in Daltons | float64 | yes |

!!! example "Cross-linked PSM"
    ```json
    {
      "sequence": "AAKPEPTIDER",
      "cross_links": [
        {
          "xl_type": "inter",
          "partner_sequence": "LKSEQUENCER",
          "partner_peptidoform": "LKSEQUENCER",
          "donor_position": 3,
          "acceptor_position": 2,
          "linker_name": "DSS",
          "linker_accession": "XLMOD:02001",
          "linker_mass": 138.0681
        }
      ]
    }
    ```

!!! note "Non-cross-linked PSMs"
    For regular (non-XL) PSMs, the `cross_links` field is `null`. This keeps the schema clean -- the struct-based design adds no overhead for non-cross-linked datasets.

## Shared Fields

Several fields in the PSM view use structures shared across other QPX views:

- For details on the `modifications` field structure, see [Modifications](modifications.md).
- For details on `additional_scores` and score semantics, see [Scores](scores.md).
- For details on `cv_params` usage and recommended terms, see [Scores & CV Terms](scores.md).

!!! note "A null `feature_id` means no linked feature"
    `feature_id` records a PSM's feature link in the exported dataset. A null
    value means no link is recorded; it does not establish identification
    quality or whether the precursor was quantified.

    OpenMS unassigned PeptideIdentifications are retained by default to preserve
    identification evidence. Some repeat evidence for a precursor linked through
    another spectrum, while others have no detected or mapped feature. Also,
    when only the PSM view is exported, **all `feature_id` values are null**, even
    for identifications assigned to a consensus feature in the source.

    The API filters the recorded links:

    ```python
    psm.with_feature()     # feature_id IS NOT NULL
    psm.without_feature()  # feature_id IS NULL
    ```

    **A `psm` -> `feature` join is expected to be partial.** When both views are
    present, join from `psm.with_feature()` on `feature_id`. The conversion option
    `--no-include-unassigned-psms` excludes source unassigned identifications;
    it does not populate feature links when the feature view is not exported.

    The OpenMS consensus converter omits features without peptide hits.

## Example

### Basic PSM Record

```json
{
  "psm_id": -4279600662514899531,
  "sequence": "AAAAAAAAAAGAAGGR",
  "peptidoform": "[Acetyl]-AAAAAAAAAAGAAGGR",
  "charge": 2,
  "scan": [42164],
  "rt": 5140.98,
  "calculated_mz": 635.3311,
  "observed_mz": 635.3315,
  "is_decoy": false,
  "posterior_error_probability": 5.58e-20,
  "predicted_rt": null,
  "run_file_name": "20200101_sample_A",
  "protein_accessions": ["Q86U42-2", "Q86U42"],
  "modifications": [
    {
      "name": "Acetyl",
      "accession": "UniMod:1",
      "positions": [
        {
          "position": 0,
          "amino_acid": null,
          "scores": []
        }
      ]
    }
  ],
  "additional_scores": [
    {
      "score_name": "andromeda_score",
      "score_value": 175.73,
      "higher_better": true
    },
    {
      "score_name": "andromeda_delta_score",
      "score_value": 160.47,
      "higher_better": true
    },
    {
      "score_name": "parent_ion_fraction",
      "score_value": 0.0,
      "higher_better": true
    }
  ],
  "cv_params": [
    {"cv_name": "dissociation method", "cv_value": "HCD"},
    {"cv_name": "normalized collision energy", "cv_value": "28"}
  ]
}
```

### PSM with Spectral Data

When spectral arrays are included, the record also contains peak-level data:

```json
{
  "psm_id": -4279600662514899531,
  "sequence": "AAAAAAAAAAGAAGGR",
  "peptidoform": "[Acetyl]-AAAAAAAAAAGAAGGR",
  "charge": 2,
  "scan": [42164],
  "rt": 5140.98,
  "calculated_mz": 635.3311,
  "observed_mz": 635.3315,
  "is_decoy": false,
  "run_file_name": "20200101_sample_A",
  "mz_array": [175.119, 289.163, 360.200, 431.236, 488.258],
  "intensity_array": [1234.5, 5678.9, 3456.7, 2345.6, 1234.5],
  "charge_array": [1, 1, 1, 1, 1],
  "ion_type_array": ["y1", "y2", "b3", "b4", "b5"]
}
```

## File Metadata

PSM Parquet files store file-level metadata as key-value pairs in the Parquet footer. The following metadata fields are defined:

| Field | Description |
|-------|-------------|
| `qpx_version` | Version of the QPX format used to generate the file |
| `software_provider` | Name and version of the software that generated the data |
| `scan_format` | Format of scan identifiers: `scan`, `index`, or `nativeId` |
| `creator` | Name of the tool or person who created the file |
| `file_type` | Type of the file (value: `psm_file`) |
| `creation_date` | Date when the file was created |
| `compression_format` | Compression algorithm used: `zstd` (default), `snappy`, `gzip`, `lzo`, or `none` |

!!! tip "Reading file metadata in Python"
    ```python
    import pyarrow.parquet as pq

    parquet_file = pq.ParquetFile("experiment.psm.parquet")
    metadata = parquet_file.schema_arrow.metadata
    for key, value in metadata.items():
        print(f"{key.decode()}: {value.decode()}")
    ```

## Notes

!!! warning "DDA-specific view"
    The PSM view is designed primarily for **DDA** (data-dependent acquisition) methods. It is **not recommended** for DIA experiments, where the feature view should be used instead. Generating a PSM file for DIA data would produce duplicated information relative to the feature view.

!!! note "PEP is a PSM-level metric"
    `posterior_error_probability` is defined only in the PSM view. It represents the probability that a *specific* peptide-spectrum match is incorrect (lower is better). All major tools export PEP as P(incorrect): Percolator (`posterior_error_prob`), MaxQuant (`PEP`). FragPipe exports PeptideProphet Probability (P(correct)), so converters must compute `PEP = 1 - probability`. The feature and peptide views do not carry PEP as a top-level field; use `additional_scores` or `best_id_score` if a derived PEP value is needed at those levels.

- **Relationship to feature view**: The PSM view captures individual spectrum matches, while the [Feature View](feature.md) aggregates these into quantified peptide features with intensity data. A single feature may correspond to multiple PSMs across different scans.
- **Protein inference**: Protein inference results should not be the primary focus of the PSM view. `protein_accessions` is an **optional** column — protein mappings can be recovered by joining through the feature and protein group views. When included, it is useful for peptide filtering and protein-level browsing. For full protein group information, use the protein group (PG) view.
- **Spectral arrays**: The `mz_array` and `intensity_array` are parallel arrays of the same length. For large-scale spectral storage, the dedicated mass spectra (mz) view is recommended.
- **Recommended additional scores**: `global_qvalue` (experiment-level PSM q-value), `rank` (peptide rank in search results), and `pg_global_qvalue` (protein group q-value used for filtering).
