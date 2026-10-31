# Processing Provenance

`provenance.parquet` captures the processing chain that generated a QPX dataset -- which tools ran, in what order, with what parameters. Its purpose is **understanding**, not reproduction: a user looking at the data should be able to answer questions like "what FDR was used?", "what search tolerance?", or "why is my protein missing?" without reading pipeline logs.

See the full YAML schema in [`provenance.yaml`](schemas/provenance.yaml).

!!! info "Understanding, not reproduction"
    `provenance.parquet` records the key decisions and parameters that shape results. It is not a workflow execution log or a Nextflow `-resume` file. For full reproducibility, refer to the pipeline's own configuration and execution artifacts.

---

## PyArrow Schema

```python
import pyarrow as pa

PARAMETER = pa.struct(
    [
        pa.field("key", pa.string()),  # parameter name (snake_case)
        pa.field("value", pa.string()),  # parameter value (as string)
    ]
)

provenance_schema = pa.schema(
    [
        # Step identity
        pa.field("step_order", pa.int32()),  # execution order (1, 2, 3...)
        pa.field("step_category", pa.string()),  # broad category
        pa.field("step_name", pa.string()),  # specific step name
        # Tool
        pa.field("tool_name", pa.string()),  # tool name (e.g. "comet", "percolator")
        pa.field("tool_version", pa.string(), nullable=True),  # tool version
        pa.field("tool_uri", pa.string(), nullable=True),  # container URI or tool URL
        # Parameters -- curated summary of key settings
        pa.field("parameters", pa.list_(PARAMETER), nullable=True),
        # Full configuration -- raw JSON for complete reference
        pa.field("config", pa.string(), nullable=True),
        # What this step produced
        pa.field("output_views", pa.list_(pa.string()), nullable=True),
    ]
)
```

## Field Reference

| Field | Description | Type | Required |
|---|---|---|---|
| `step_order` | Execution order (1, 2, 3...) | int32 | yes |
| `step_category` | Broad category of the processing step | string | yes |
| `step_name` | Specific name for the step | string | yes |
| `tool_name` | Name of the tool or software | string | yes |
| `tool_version` | Version of the tool | string | no |
| `tool_uri` | Container URI (e.g., Docker/Singularity image) or tool URL | string | no |
| `parameters` | Key-value pairs of important settings (curated summary) | list[struct{key, value}] | no |
| `config` | Full tool/step configuration as a JSON string | string | no |
| `output_views` | Which QPX views this step generated (e.g., `["psm"]`, `["feature", "pg"]`) | list[string] | no |

!!! tip "`parameters` vs `config`"
    **`parameters`** is the curated summary -- a short list of the most important settings, queryable with SQL `UNNEST`. Use it for quick lookups like "what FDR was used?".

    **`config`** is the complete reference -- the full tool configuration as a JSON string. Use it when you need every setting, not just the highlights. For a quantms pipeline, this would be the Nextflow `params.json`; for a Comet search, the full `comet.params` converted to JSON.

## Recorded conversion options

`qpxc convert` adds the options a conversion ran with to the `parameters` of its `format_conversion` step, defaults included: every option that shapes the output, such as `qvalue_threshold` (`none` when unset, so "not filtered" is on record), `structures` or `rank1_only`. An input path is recorded by file name, and an input file also by its SHA-256 as `<option>_sha256`. Output locations and resource limits are not recorded.

## Step Categories

Each processing step belongs to a category that groups steps by function.

| `step_category` | Description | Typical tools |
|---|---|---|
| `workflow` | Pipeline orchestrator | quantms, nf-core, Nextflow |
| `raw_conversion` | Raw file conversion to open format | ThermoRawFileParser, msconvert |
| `database_search` | Sequence database search | Comet, MSGF+, MaxQuant, Sage, DIA-NN |
| `de_novo_peptide_sequencing` | Peptide sequence prediction without a sequence database search | Casanovo, InstaNovo, π-PrimeNovo |
| `rescoring` | PSM rescoring and validation | Percolator, PeptideProphet, mokapot |
| `fdr_filtering` | FDR threshold application | OpenMS, custom scripts |
| `quantification` | Peptide/protein quantification | FlashLFQ, DIA-NN, IonQuant |
| `protein_inference` | Protein group inference | Epifany, ProteinProphet, IDPicker |
| `normalization` | Intensity normalization | quantms, MSstats |
| `statistical_analysis` | Differential expression / statistics | MSstats, limma, DEqMS, proDA |
| `expression` | Absolute expression computation | quantms (iBAQ) |

## Example

### quantms LFQ pipeline

```json
[
  {
    "step_order": 1,
    "step_category": "workflow",
    "step_name": "pipeline_execution",
    "tool_name": "quantms",
    "tool_version": "1.3.0",
    "tool_uri": null,
    "parameters": [
      {"key": "nextflow_version", "value": "23.10.0"},
      {"key": "profile", "value": "docker"},
      {"key": "acquisition_method", "value": "lfq"}
    ],
    "config": "{\"input\":\"samplesheet.csv\",\"database\":\"UP000005640_9606.fasta\",\"acquisition_method\":\"lfq\",\"search_engines\":\"comet\",\"psm_fdr\":0.01,\"protein_fdr\":0.01,\"enable_mod_localization\":true,\"precursor_mass_tolerance\":10,\"fragment_mass_tolerance\":0.02,\"variable_mods\":\"Oxidation (M)\",\"fixed_mods\":\"Carbamidomethyl (C)\",\"max_mods\":3,\"match_between_runs\":true}",
    "output_views": null
  },
  {
    "step_order": 2,
    "step_category": "raw_conversion",
    "step_name": "raw_to_mzml",
    "tool_name": "thermorawfileparser",
    "tool_version": "1.4.2",
    "tool_uri": "docker://ghcr.io/bigbio/thermorawfileparser:1.4.2",
    "parameters": [
      {"key": "output_format", "value": "mzML"},
      {"key": "peak_picking", "value": "true"}
    ],
    "output_views": ["mz"]
  },
  {
    "step_order": 3,
    "step_category": "database_search",
    "step_name": "sequence_search",
    "tool_name": "comet",
    "tool_version": "2024.01.0",
    "tool_uri": "docker://ghcr.io/bigbio/quantms-comet:2024.01.0",
    "parameters": [
      {"key": "precursor_mass_tolerance", "value": "10ppm"},
      {"key": "fragment_mass_tolerance", "value": "0.02Da"},
      {"key": "isotope_error", "value": "0/1/2"},
      {"key": "fasta_database", "value": "UP000005640_9606.fasta"},
      {"key": "fasta_checksum", "value": "sha256:a1b2c3d4..."}
    ],
    "output_views": ["psm"]
  },
  {
    "step_order": 4,
    "step_category": "rescoring",
    "step_name": "psm_rescoring",
    "tool_name": "percolator",
    "tool_version": "3.06.1",
    "tool_uri": "docker://ghcr.io/bigbio/quantms-percolator:3.06.1",
    "parameters": [
      {"key": "train_fdr", "value": "0.01"},
      {"key": "test_fdr", "value": "0.01"}
    ],
    "output_views": ["psm"]
  },
  {
    "step_order": 5,
    "step_category": "fdr_filtering",
    "step_name": "multi_level_fdr",
    "tool_name": "openms",
    "tool_version": "3.1.0",
    "tool_uri": null,
    "parameters": [
      {"key": "psm_fdr", "value": "0.01"},
      {"key": "peptide_fdr", "value": "0.01"},
      {"key": "protein_fdr", "value": "0.01"}
    ],
    "output_views": null
  },
  {
    "step_order": 6,
    "step_category": "protein_inference",
    "step_name": "protein_grouping",
    "tool_name": "epifany",
    "tool_version": "3.1.0",
    "tool_uri": null,
    "parameters": [
      {"key": "algorithm", "value": "bayesian"},
      {"key": "protein_fdr", "value": "0.01"}
    ],
    "output_views": ["pg"]
  },
  {
    "step_order": 7,
    "step_category": "quantification",
    "step_name": "peptide_quantification",
    "tool_name": "flashlfq",
    "tool_version": "1.2.0",
    "tool_uri": null,
    "parameters": [
      {"key": "match_between_runs", "value": "true"},
      {"key": "mbr_rt_tolerance", "value": "2.5min"}
    ],
    "output_views": ["feature"]
  },
  {
    "step_order": 8,
    "step_category": "normalization",
    "step_name": "intensity_normalization",
    "tool_name": "quantms",
    "tool_version": "1.3.0",
    "tool_uri": null,
    "parameters": [
      {"key": "method", "value": "median"},
      {"key": "log_transform", "value": "true"}
    ],
    "output_views": null
  },
  {
    "step_order": 9,
    "step_category": "statistical_analysis",
    "step_name": "group_comparison",
    "tool_name": "msstats",
    "tool_version": "4.10.0",
    "tool_uri": null,
    "parameters": [
      {"key": "normalization", "value": "equalizeMedians"},
      {"key": "imputation", "value": "MinProb"},
      {"key": "fdr_threshold", "value": "0.05"}
    ],
    "output_views": ["de"]
  },
  {
    "step_order": 10,
    "step_category": "expression",
    "step_name": "ibaq_calculation",
    "tool_name": "quantms",
    "tool_version": "1.3.0",
    "tool_uri": null,
    "parameters": [
      {"key": "method", "value": "ibaq"},
      {"key": "log_transform", "value": "true"}
    ],
    "output_views": ["ae"]
  }
]
```

## Common Queries

```sql
-- Full processing chain in order
SELECT step_order, step_category, tool_name, tool_version
FROM 'PXD014414.provenance.parquet'
ORDER BY step_order;

-- What search engine and tolerances were used?
SELECT tool_name, tool_version, parameters
FROM 'PXD014414.provenance.parquet'
WHERE step_category = 'database_search';

-- What FDR thresholds were applied anywhere in the pipeline?
SELECT s.step_name, p.key, p.value
FROM 'PXD014414.provenance.parquet' s,
     UNNEST(s.parameters) AS p
WHERE p.key LIKE '%fdr%';

-- Which tools produced the feature view?
SELECT step_name, tool_name, tool_version
FROM 'PXD014414.provenance.parquet'
WHERE list_contains(output_views, 'feature');

-- All container images used
SELECT step_name, tool_name, tool_uri
FROM 'PXD014414.provenance.parquet'
WHERE tool_uri IS NOT NULL;
```

## User Story: "Why is my protein missing?"

A researcher notices that a protein they expected is absent from the results. They can trace through the provenance to understand why:

| Check | Query | What to look for |
|---|---|---|
| Was it in the database? | `database_search` → `fasta_database` | Protein must be in the FASTA file |
| Were tolerances too tight? | `database_search` → `precursor_mass_tolerance` | Tight tolerances may miss modified peptides |
| Was FDR too strict? | `fdr_filtering` → `psm_fdr`, `protein_fdr` | Strict thresholds remove borderline identifications |
| Was it collapsed into another group? | `protein_inference` → `algorithm` | Check the PG view for related protein groups |
| Was match-between-runs off? | `quantification` → `match_between_runs` | MBR increases identifications across runs |

## Relationship to Existing Metadata

`provenance.parquet` complements, not replaces, the existing metadata files:

| Existing | Purpose | Provenance adds |
|---|---|---|
| `dataset.parquet` → `software_name`, `software_version` | Quick project-level lookup | Full tool chain with all versions |
| `run.parquet` → `enzymes`, `modification_parameters` | Per-run search config | Pipeline-level search engine settings |
| `run.parquet` → `instrument`, `dissociation_method` | Instrument info | Stays in run.parquet (per-run, not per-pipeline) |

---

## Related Pages

- [Dataset Metadata](dataset.md) -- project-level metadata
- [Run Metadata](run.md) -- per-run instrument and search parameters
- [QPX Format Overview](index.md)
