# Absolute Expression

Absolute expression (AE) quantification determines the baseline amount of a target protein in a sample. In proteomics, the primary computational method is the intensity-based [absolute quantification (iBAQ)](https://www.nature.com/articles/nature11848) method.

## Use cases

- Store and retrieve baseline protein abundance (iBAQ values) per sample.
- Understand expression profiles of a protein across different conditions, tissues, and organisms.
- Provide a proxy for protein copy number estimation and concentration in biological samples.
- Enable fast visualization of absolute expression results.
- Integrate with scverse ecosystem tools (scanpy, scvi-tools, muon) for multi-omics analysis.

## Format

The absolute expression view uses **AnnData** (`.h5ad`) as its primary format. AnnData is a matrix-format standard from the [scverse](https://scverse.org/) ecosystem that naturally represents the samples-by-proteins structure of absolute expression data.

For general AnnData concepts and conventions shared with the differential expression view, see [AnnData Concepts](anndata.md).

### AnnData structure

```text
AnnData (n_obs x n_vars = samples x proteins)
    obs:    sample metadata (organism, tissue, disease, ...)
    var:    protein metadata (gene_name, ...)
    X:      ibaq_log (log2 of ibaq_raw, the primary quantification matrix)
    layers: ibaq_raw, ibaq_ppb, copies_per_cell, concentration_nm
    uns:    file-level metadata (qpx_version, project_accession, ...)
```

### Slots

| AnnData slot | Content | Description |
|-------------|---------|-------------|
| `X` | `ibaq_log` | **Primary data matrix** -- the logarithm of `ibaq_raw` in the base of `uns["log_base"]`, normally 2 (samples x proteins); NaN for a protein not detected in a sample |
| `obs` | Sample metadata | One row per sample with experimental annotations |
| `var` | Protein metadata | One row per protein with gene names |
| `layers["ibaq_raw"]` | Raw iBAQ | Untransformed iBAQ values |
| `layers["ibaq_ppb"]` | Relative iBAQ | Parts-per-billion normalized iBAQ |
| `layers["copies_per_cell"]` | Copy number | Estimated protein copies per cell |
| `layers["concentration_nm"]` | Concentration | Estimated protein concentration in nM |
| `uns` | File metadata | QPX version, project accession, factor values, etc. |

!!! note "Convention for primary matrix"
    `ibaq_log` is chosen as the primary quantification (maps to `X`). This parallels the scRNA-seq convention where `X` holds log-normalized counts and `layers["counts"]` holds raw counts.

### obs (sample metadata)

Each row in `obs` represents one sample. The index is the sample accession.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `sample_accession` | Sample accession from SDRF (used as index) | `string` | Yes |
| `organism` | Organism (promoted from SDRF) | `string` | No |
| `organism_part` | Tissue or organ (promoted from SDRF) | `string` | No |
| `disease` | Disease condition (promoted from SDRF) | `string` | No |
| `cell_line` | Cell line (promoted from SDRF) | `string` | No |
| `biological_replicate` | Biological replicate number | `int32` | No |
| `technical_replicate` | Technical replicate number | `int32` | No |

Additional SDRF factor values may be included as extra columns in `obs`.

### var (protein metadata)

Each row in `var` is one protein group, and the index is unique. `uns["var_key"]` states what it holds:

- `pg_accessions`: the group's sorted members joined by `;`, as in the `proteins` modality of the [MuData view](anndata.md). This is the key to use, since it names one group.
- `anchor_protein`: the group's leading protein. Groups that share a leader are one row, with their peptides pooled. `qpxc transform quantify` writes this, because it computes iBAQ per protein.

`pg_id` is not a key here: it depends on the run and label ([Row Identity](identity.md)).

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `protein` | The key above, used as index | `string` | Yes |
| `gene_name` | Gene symbol | `string` | No |

### uns (file-level metadata)

| Key | Description | Type |
|-----|-------------|------|
| `qpx_version` | Version of the QPX format | `string` |
| `file_type` | Always `"absolute_expression"` | `string` |
| `project_accession` | Project accession in PRIDE Archive | `string` |
| `project_title` | Project title | `string` |
| `factor_value` | Factor value from SDRF | `string` |
| `creation_date` | Date when the file was created | `string` |
| `creator` | Name of the tool that created the file | `string` |
| `log_base` | Base of the logarithm in `X` | `int` |
| `var_key` | What the `var` index holds: `pg_accessions` or `anchor_protein` | `string` |
| `aggregation` | How a sample's value combines its runs and fractions, as the tool states it; `qpxc transform quantify` writes `mokume_peptides_to_protein_per_sample` | `string` |
| `quantification_parameters` | JSON object of the method and its settings, such as enzyme, peptide length, ploidy and cell protein concentration | `string` |
| `standards` | Spike-in standards behind absolute amounts: a table with one row per standard and the columns `protein`, `amount` and `unit`, for example `fmol` | table |

## Producing and checking the file

`qpxc transform quantify --method ibaq --fasta ... -o PXD000000.ae.h5ad` writes this layout, and `qpxc validate --dataset-path` checks an `.ae.h5ad` file it finds in the dataset: unique indexes, no negative layer, `X` as the declared logarithm of `ibaq_raw`, the `standards` columns, and a warning when `log_base`, `var_key` or `aggregation` is missing.

## Example

### Creating an AE AnnData

```python
import anndata as ad
import numpy as np
import pandas as pd

# Sample metadata (obs)
obs = pd.DataFrame(
    {
        "organism": ["Homo sapiens", "Homo sapiens"],
        "organism_part": ["heart", "liver"],
        "disease": ["normal", "normal"],
        "biological_replicate": [1, 1],
    },
    index=["PXD000000-Sample-1", "PXD000000-Sample-2"],
)

# Protein metadata (var)
var = pd.DataFrame(
    {
        "gene_name": ["A1BG", "HBB", "TP53"],
    },
    index=["P04217", "P68871", "P04637"],
)

# Primary matrix: ibaq_log (samples x proteins)
X = np.array(
    [
        [12.4714, 7.8734, 11.1957],  # Sample-1
        [10.2697, 12.7290, 11.7552],  # Sample-2
    ]
)

# Create AnnData
adata = ad.AnnData(X=X, obs=obs, var=var)

# Add alternative quantifications as layers
adata.layers["ibaq_raw"] = np.array(
    [
        [5678.9, 234.5, 2345.6],
        [1234.5, 6789.0, 3456.7],
    ]
)

# Add file-level metadata
adata.uns["qpx_version"] = "1.3"
adata.uns["file_type"] = "absolute_expression"
adata.uns["project_accession"] = "PXD000000"
adata.uns["log_base"] = 2
adata.uns["var_key"] = "pg_accessions"
adata.uns["aggregation"] = "mokume_peptides_to_protein_per_sample"

# Save
adata.write("PXD000000.ae.h5ad")
```

### Querying AE data

```python
import anndata as ad

adata = ad.read_h5ad("PXD000000.ae.h5ad")

# Expression of a specific gene across all samples
tp53_idx = adata.var["gene_name"] == "TP53"
print(adata[:, tp53_idx].X)

# All heart tissue samples
heart = adata[adata.obs["organism_part"] == "heart"]
print(heart.X)

# Most abundant proteins in a sample
sample = adata[0, :]
top_proteins = sample.var.index[np.argsort(sample.X.flatten())[::-1][:10]]
```

## File naming

AE AnnData files follow the QPX naming convention:

```text
{PREFIX}.ae.h5ad
```

Example: `PXD000000.ae.h5ad`

## Notes

!!! tip "Missing values"
    Proteins not detected in a sample result in `NaN` values in the AnnData matrix. This is consistent with how scRNA-seq handles dropout events.

!!! tip "Multi-omics integration"
    AE AnnData can be concatenated with scRNA-seq AnnData for joint analysis (e.g., CITE-seq + bulk proteomics comparison). The scverse ecosystem provides tools like `muon` for multi-modal data integration.

!!! note "Relationship to other views"
    The AE view derives from the [Protein Group View](pg.md), which provides per-file protein quantification. AE aggregates across files and computes iBAQ-based absolute quantities. For differential comparisons between conditions, see the [Differential Expression View](differential.md).

!!! warning "Protein group encoding"
    With `var_key` `anchor_protein`, the `var` index is a single representative protein accession. The full protein group membership is available in the [Protein Group View](pg.md).
