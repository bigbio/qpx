# Differential Expression

The differential expression (DE) view stores statistical results comparing protein abundance between experimental conditions. It contains fold changes, p-values, and associated metrics for each protein across contrasts. This is a key output of the proteomics data analysis workflow, enabling identification of differentially expressed proteins.

## Use cases

- Store differentially expressed proteins between contrasts with fold changes and p-values.
- Enable visualization using [Volcano Plots](https://en.wikipedia.org/wiki/Volcano_plot_(statistics)) and other statistical plots.
- Enable integration with other omics data for multi-omics analysis.
- Store metadata about the statistical method and contrast definitions alongside the results.

## Format

The differential expression view uses **AnnData** (`.h5ad`) as its primary format. DE results are stored in AnnData's `uns` (unstructured metadata) slot, with the field names of scanpy's `rank_genes_groups`, but not its layout. The AnnData object can optionally include an expression matrix in `X` (from the AE view) alongside the DE results.

For general AnnData concepts and conventions shared with the absolute expression view, see [AnnData Concepts](anndata.md).

### AnnData structure

```text
AnnData (n_obs x n_vars = samples x proteins)
    obs:    sample metadata (organism, tissue, disease, ...)
    var:    protein metadata (gene_name, ...)
    X:      protein quantification matrix (e.g. MaxLFQ, TMT intensity)
    uns:    DE results + file-level metadata
```

### DE results in `uns`

DE results are stored in `uns["de_results"]` as a dictionary keyed by contrast identifier: one entry per contrast, each a set of arrays parallel to `names`. A contrast is one test, whatever it compares, so an F test or an interaction is an entry too.

```python
adata.uns["de_results"] = {
    "cancer_vs_normal": {
        "names": np.array(["P04217", "P68871", ...]),  # protein accessions
        "gene_names": np.array(["A1BG", "HBB", ...]),  # gene symbols
        "logfoldchanges": np.array([-1.542, 0.891, ...]),  # log2 fold change
        "scores": np.array([-8.567, 4.321, ...]),  # t-value / test statistic
        "pvals": np.array([0.0001, 0.002, ...]),  # raw p-values
        "pvals_adj": np.array([0.005, 0.045, ...]),  # adjusted p-values (BH)
        "se": np.array([0.18, 0.206, ...]),  # standard error
        "df": np.array([29.05, 35.0, ...]),  # degrees of freedom
        "is_significant": np.array([True, True, ...]),  # pre-computed significance
        "issue": np.array([None, None, ...]),  # inference issues
        "condition_test": "squamous cell carcinoma",  # test condition
        "condition_reference": "normal",  # reference condition
    }
}
```

### Fields per contrast

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `names` | Tested entities, of the type `uns["entity_type"]` states (sorted by p-value) | `array[string]` | Yes |
| `gene_names` | Gene symbols corresponding to each protein | `array[string]` | No |
| `logfoldchanges` | Log2 fold change (test / reference). Not required for an F or likelihood-ratio test, which has no single fold change | `array[float64]` | Yes |
| `scores` | Test statistic, of the type `statistic_type` states | `array[float64]` | No |
| `statistic_type` | `t`, `z`, `F` or `LRT` | `string` | No |
| `pvals` | Raw p-values | `array[float64]` | Yes |
| `pvals_adj` | Adjusted p-values (multiple testing corrected) | `array[float64]` | Yes |
| `se` | Standard error of the log2 fold change | `array[float64]` | No |
| `df` | Degrees of freedom; moderated tests give fractions, such as limma's 29.05 | `array[float64]` | No |
| `ave_expr` | Mean log2 abundance over the samples of the contrast | `array[float64]` | No |
| `n_test`, `n_reference` | Samples with a value in the test and the reference condition | `array[int32]` | No |
| `missing_fraction` | Fraction of the contrast's samples without a value, 0-1 | `array[float64]` | No |
| `is_significant` | Pre-computed significance at the stored FDR threshold | `array[bool]` | No |
| `issue` | Issue with protein quantification, if any | `array[string]` | No |
| `condition_test` | Test condition in the contrast. Required without `contrast_definition` | `string` | Yes |
| `condition_reference` | Reference/control condition in the contrast. Required without `contrast_definition` | `string` | Yes |
| `contrast_definition` | The tested hypothesis as the tool states it, for what is not a pair of conditions: `(A_treated-A_control)-(B_treated-B_control)` for an interaction, or the dropped terms of an F test | `string` | No |

### uns (file-level metadata)

| Key | Description | Type |
|-----|-------------|------|
| `qpx_version` | Version of the QPX format | `string` |
| `file_type` | Always `"differential_expression"` | `string` |
| `project_accession` | Project accession in PRIDE Archive | `string` |
| `entity_type` | What `names` holds: `protein_group` (the default), `gene_group`, `peptidoform`, `site` or `proteoform` | `string` |
| `method_details` | JSON object of the model: tool, version, design formula, covariates, normalization, imputation | `string` |
| `statistical_method` | Statistical method used (e.g., `msstats_group_comparison`, `limma`, `deqms`) | `string` |
| `correction_method` | Multiple testing correction method (e.g., `BH`) | `string` |
| `fdr_threshold` | FDR threshold used for `is_significant` | `string` |
| `factor_names` | JSON array of factor names from SDRF | `string` |
| `contrasts` | JSON array of contrast identifiers | `string` |
| `creation_date` | Date when the file was created | `string` |
| `creator` | Name of the tool that created the file | `string` |

### Tool mapping

QPX's DE fields map to common statistical tool outputs:

| QPX field | scanpy | MSstats | PyDESeq2 | limma (R) | edgeR (R) |
|-----------|--------|---------|----------|-----------|-----------|
| `names` | `names` | `Protein` | (index) | (rownames) | (rownames) |
| `logfoldchanges` | `logfoldchanges` | `log2FC` | `log2FoldChange` | `logFC` | `logFC` |
| `scores` | `scores` | `Tvalue` | `stat` | `t` | -- |
| `pvals` | `pvals` | `pvalue` | `pvalue` | `P.Value` | `PValue` |
| `pvals_adj` | `pvals_adj` | `adj.pvalue` | `padj` | `adj.P.Val` | `FDR` |
| `se` | -- | `SE` | `lfcSE` | -- | -- |
| `df` | -- | `DF` | -- | `df.residual` | -- |

## Reading and checking the file

`qpxc validate --dataset-path` checks a `.de.h5ad` file it finds in the dataset: unique indexes, the required fields, arrays of one length, p-values in [0, 1] and `entity_type`. The MuData view reshapes `uns["de_results"]` into a contrasts x entities matrix, where an entity a contrast did not test is NaN, as is the fold change of an F test. It also reads a long `pandas.DataFrame` in `uns["de_results"]` with the columns `contrast`, `protein` and `log2FC`, and optionally `pvals_adj`, `scores`, `pvals` and `se`.

## Example

### Creating a DE AnnData

```python
import anndata as ad
import numpy as np
import pandas as pd

# Protein metadata
var = pd.DataFrame(
    {
        "gene_name": ["A1BG", "HBB", "TP53"],
    },
    index=["P04217", "P68871", "P04637"],
)

# Sample metadata (optional — can be empty if only storing DE results)
obs = pd.DataFrame(index=["Sample-1", "Sample-2", "Sample-3", "Sample-4"])

# Create AnnData (X can be empty or from AE)
adata = ad.AnnData(
    X=np.zeros((len(obs), len(var))),
    obs=obs,
    var=var,
)

# Add DE results
adata.uns["de_results"] = {
    "cancer_vs_normal": {
        "names": np.array(["P04217", "P68871", "P04637"]),
        "gene_names": np.array(["A1BG", "HBB", "TP53"]),
        "logfoldchanges": np.array([-1.542, 0.891, 2.345]),
        "scores": np.array([-8.567, 4.321, 12.456]),
        "pvals": np.array([0.0001, 0.002, 0.00001]),
        "pvals_adj": np.array([0.005, 0.045, 0.001]),
        "se": np.array([0.18, 0.206, 0.188]),
        "df": np.array([37.0, 35.0, 37.0]),
        "is_significant": np.array([True, True, True]),
        "issue": np.array([None, None, None]),
        "condition_test": "squamous cell carcinoma",
        "condition_reference": "normal",
    }
}

# Add file-level metadata
adata.uns["qpx_version"] = "1.3"
adata.uns["file_type"] = "differential_expression"
adata.uns["statistical_method"] = "msstats_group_comparison"
adata.uns["correction_method"] = "BH"
adata.uns["fdr_threshold"] = "0.05"
adata.uns["contrasts"] = '["cancer_vs_normal"]'

# Save
adata.write("PXD000000.de.h5ad")
```

### Querying DE results

```python
import pandas as pd

adata = ad.read_h5ad("PXD000000.de.h5ad")

# Get significant up-regulated proteins in a contrast
contrast = adata.uns["de_results"]["cancer_vs_normal"]
de_df = pd.DataFrame(
    {
        "protein": contrast["names"],
        "gene_name": contrast["gene_names"],
        "log2fc": contrast["logfoldchanges"],
        "pvalue": contrast["pvals"],
        "adj_pvalue": contrast["pvals_adj"],
        "significant": contrast["is_significant"],
    }
)

# Filter: significant and up-regulated
up_regulated = de_df[(de_df["significant"]) & (de_df["log2fc"] > 1.0)]
print(up_regulated.sort_values("adj_pvalue"))
```

## File naming

DE AnnData files follow the QPX naming convention:

```text
{PREFIX}.de.h5ad
```

Example: `PXD000000.de.h5ad`

## Notes

!!! note "Bundling AE and DE"
    AE and DE results can be stored in the same AnnData file. The expression matrix (`X` and `layers`) comes from the AE view, and DE results are added to `uns["de_results"]`. This enables a single file that contains both the expression data and the statistical comparisons.

!!! note "scanpy"
    `de_results` shares field names with scanpy's `rank_genes_groups` but not its layout, which has one structured array per field with a column per group, so `sc.pl.rank_genes_groups` does not read it. Build `uns["rank_genes_groups"]` from it to use scanpy's plots.

!!! note "Relationship to other views"
    The DE view takes as input the protein quantification from the [Protein Group View](pg.md) or the [Absolute Expression View](absolute.md). For absolute expression values per sample, see the [Absolute Expression View](absolute.md).

!!! warning "Protein group encoding"
    Protein groups in the DE results use a single representative protein accession. The full protein group membership is available in the [Protein Group View](pg.md).
