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
    X:      ibaq_log (primary quantification matrix)
    layers: ibaq_raw, ibaq_ppb, copies_per_cell, concentration_nm
    uns:    file-level metadata (qpx_version, project_accession, ...)
```

### Slots

| AnnData slot | Content | Description |
|-------------|---------|-------------|
| `X` | `ibaq_log` | **Primary data matrix** -- log-transformed iBAQ values (samples x proteins) |
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

Each row in `var` represents one protein. The index is the protein accession.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `protein` | Protein accession (UniProt), used as index | `string` | Yes |
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
        [8.48, 6.23, 7.91],  # Sample-1
        [7.12, 9.45, 8.03],  # Sample-2
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
    Protein groups in the `var` index are written as a single representative protein accession. The full protein group membership is available in the [Protein Group View](pg.md).
