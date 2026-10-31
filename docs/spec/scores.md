# Scores & CV Terms

QPX provides two structured fields for attaching metadata to records beyond the core columns: `additional_scores` for numeric quality metrics and `cv_params` for controlled vocabulary annotations. Both are stored as arrays of structs, making them extensible without schema changes.

## Additional scores

The `additional_scores` field captures search engine scores, quality metrics, and other numeric values associated with a PSM, feature, peptide, or protein group. Each entry is self-describing: it carries its name, value, and an indication of whether higher values are better.

### Struct definition

```text
additional_scores: array[struct{
    score_name:    string,       -- Score identifier in snake_case (e.g. "comet_xcorr")
    score_value:   float,        -- Numeric value
    higher_better: bool, null    -- true = higher is better; false = lower is better; null = unknown
}]
```

The `higher_better` field makes the data self-describing so that downstream consumers can interpret scores without looking up external ontology definitions.

A producer writes `null` when it does not know the direction; it never defaults to `true`. Most PSI-MS score terms declare no order, so QPX converters take the direction from the term when it has one, treat e-values, q-values, p-values, PEP and FDR as lower-is-better and engine match scores and probabilities as higher-is-better, and leave anything else `null`. A consumer ranking PSMs by a score whose `higher_better` is `null` must decide the direction itself.

!!! important "snake_case naming convention"
    All score names in QPX use `snake_case` -- no colons, dots, spaces, or mixed case. This ensures score names are valid identifiers in SQL, Python, R, and any query language without quoting. The [Ontology Mapping](ontology.md) maps each snake_case name to its proper ontology term and accession. `qpxc validate` warns on other `additional_scores` names.

### Example

```json
[
  {"score_name": "comet_xcorr",    "score_value": 3.42,   "higher_better": true},
  {"score_name": "global_qvalue",  "score_value": 0.0012, "higher_better": false},
  {"score_name": "rank",           "score_value": 1.0,    "higher_better": false}
]
```

### Typed q-value columns

Filtering on a nested score reads about 100 times the data of a flat column, and a filter on one name misses other producers. Writers therefore copy q-values from `additional_scores` to flat columns (specification 1.3), which is the portable way to filter: `psm.confident()` and `feature.confident()` use them. Score names are compared in `snake_case`, and the first match of a row wins.

| Column | Level and scope | Filled from |
| ------ | --------------- | ----------- |
| `psm.qvalue` | PSM, as computed by the source | `qvalue`, `q_value`, `psm_qvalue`, `psm_level_qvalue`, `<engine>_qvalue`, `<engine>_q_value`; not names of a higher level (`peptide`, `peptidoform`, `protein`, `pg`, `gg`, `global`, `precursor`) |
| `feature.peptide_qvalue` | Peptide, as computed by the source | `peptide_qvalue`, `peptide_level_qvalue` |
| `feature.precursor_qvalue` | Precursor, per run | `precursor_qvalue`, `qvalue`, `q_value`, `psm_level_qvalue` |
| `feature.global_qvalue` | Precursor, per experiment | `global_qvalue` |

A column stays null when the source has no such score, and the score stays in `additional_scores`.

### Recommended score names

The following score names are commonly used across QPX views. All names are `snake_case`. For formal ontology names and accessions, see the [Ontology Mapping](ontology.md#scores).

| Score name | Typical view(s) | Direction | Description |
| ---------- | ---------------- | --------- | ----------- |
| `posterior_error_probability` | PSM | lower is better | Posterior error probability — probability that the PSM is incorrect. Ranges 0.0–1.0. This is a **top-level field** in the PSM view, not stored in `additional_scores` |
| `global_qvalue` | PSM, Feature | lower is better | Global q-value at the experiment level |
| `pg_global_qvalue` | PSM, Feature | lower is better | Protein group global q-value used to filter at the protein group level |
| `rank` | PSM | lower is better | Rank of the peptide in the search engine results (1 = best) |
| `comet_xcorr` | PSM | higher is better | Cross-correlation score from the Comet search engine |
| `comet_deltacn` | PSM | higher is better | Delta CN score from Comet |
| `comet_expect` | PSM | lower is better | Expectation value from Comet |
| `msgf_raw_score` | PSM | higher is better | Raw score from MS-GF+ |
| `msgf_spec_evalue` | PSM | lower is better | Spectral E-value from MS-GF+ |
| `sage_hyperscore` | PSM | higher is better | Hyperscore from Sage search engine |
| `precursor_qvalue` | Feature | lower is better | Run-level precursor q-value (`Q.Value`), as emitted by the DIA-NN converter |
| `pg_qvalue` | Feature | lower is better | Run-level protein-group q-value, as emitted by the DIA-NN converter |
| `pg_pep` | PG | lower is better | Protein-group posterior error probability (`PG.PEP`), as emitted by the DIA-NN converter |
| `diann_lib_qvalue` | Feature | lower is better | Library q-value from DIA-NN |
| `diann_cscore` | Feature | higher is better | Confidence score from DIA-NN |
| `diann_evidence` | Feature | higher is better | Evidence score from DIA-NN |
| `diann_ms1_profile_corr` | Feature | higher is better | MS1 profile correlation from DIA-NN |
| `diann_mass_evidence` | Feature | higher is better | Mass evidence score from DIA-NN |
| `diann_quantity_quality` | Feature | higher is better | Quantity quality score from DIA-NN |
| `consensus_support` | PSM | higher is better | Number of search engines supporting the identification |

!!! tip "Naming convention for new scores"
    When adding tool-specific scores, use `snake_case` with the tool name as prefix: `{tool}_{score}`. For example, `comet_xcorr`, `diann_cscore`, `sage_hyperscore`. Register the mapping to the proper ontology term in `ontology.parquet` (see [Ontology Mapping](ontology.md)).

!!! warning "Names in this table are the ones producers actually write"
    Earlier revisions listed `diann_qvalue` and `diann_global_qvalue`, which no converter emits — the DIA-NN converter writes `precursor_qvalue`, `global_qvalue` and `pg_qvalue`. The table now records the emitted names so a consumer can match on them.

    `feature.peptide_qvalue` is reserved for peptide-level confidence. DIA-NN `Q.Value` is precursor-level, so it remains in `additional_scores` under `precursor_qvalue` and must not be copied into `feature.peptide_qvalue`.

### Protein-level additional scores

At the protein group level, `additional_scores` follows the same struct definition. The `score_value` is a float that applies to the entire protein group entry. For per-protein scores within a group, the values array index corresponds to the `pg_accessions` array index.

## Controlled vocabulary terms

The `cv_params` field stores key-value annotations drawn from controlled vocabularies. Unlike `additional_scores` (which are always numeric), CV parameters can carry string values or be value-less (presence-only).

### Struct definition

```text
cv_params: array[struct{
    cv_name:  string,        -- Term name from a controlled vocabulary
    cv_value: string, null   -- Term value; null if the term is a flag (presence-only)
}]
```

### Example

```json
[
  {"cv_name": "ms level",          "cv_value": "2"},
  {"cv_name": "deconvoluted data", "cv_value": null},
  {"cv_name": "prot:FDR threshold", "cv_value": "0.01"},
  {"cv_name": "number of unmatched peaks", "cv_value": "3"}
]
```

!!! note
    The `cv_name` is always required. The `cv_value` is optional -- a `null` value indicates a boolean/flag-style term where mere presence is meaningful (e.g., "deconvoluted data").

## Where these fields are used

| View | `additional_scores` | `cv_params` |
| ---- | :-----------------: | :---------: |
| PSM (`psm_file`) | Yes | Yes |
| Feature (`feature_file`) | Yes | Yes |
| API Views | Yes (as `best_id_score`) | -- |
| Protein Group (`pg_file`) | Yes | -- |
| MZ (`mz_file`) | -- | Yes |

!!! warning
    Do not store search engine scores in `cv_params`. Numeric scores belong in `additional_scores`, which provides the `higher_better` flag and enforces float typing. Reserve `cv_params` for non-numeric annotations and metadata.

## Further reading

- [Modifications](modifications.md) -- localization scores use a similar score struct pattern
- [Intensities](intensities.md) -- additional processed values for features and protein groups
- [QPX Format Overview](index.md) -- full list of views and concepts
