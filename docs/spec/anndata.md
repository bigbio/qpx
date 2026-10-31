# AnnData Concepts

QPX uses [AnnData](https://anndata.readthedocs.io/) (`.h5ad`) as the primary format for expression views -- both [Absolute Expression](absolute.md) and [Differential Expression](differential.md). This page describes the general AnnData conventions shared across expression views.

## Why AnnData?

AnnData is the standard interchange format for the [scverse](https://scverse.org/) ecosystem (scanpy, scvi-tools, muon). It provides:

- **Matrix-form representation** (samples x proteins) natural for expression data
- **Multi-layer storage** for alternative quantifications (raw, normalized, log-transformed)
- **Structured metadata** for samples (`obs`), proteins (`var`), and unstructured data (`uns`)
- **HDF5 backend** for efficient storage and partial reads
- **Ecosystem interop** with scanpy plotting, scvi-tools, muon multi-omics, and Lamin.ai

## AnnData structure overview

```text
AnnData
  X            (n_obs x n_vars)    Primary data matrix
  obs          (n_obs x ...)       Observation (sample) metadata
  var          (n_vars x ...)      Variable (protein) metadata
  layers       (n_obs x n_vars)    Alternative quantification matrices
  uns          (dict)              Unstructured metadata (DE results, file metadata)
  obsm/varm    (optional)          Embeddings and dimensionality reductions
```

![AnnData structure diagram](../images/anndata.png)

## Conventions

### obs index = sample accession

The `obs` index uses sample accessions from the SDRF. Sample-level metadata (organism, tissue, disease, replicate info) is stored as columns in `obs`.

### var index = protein accession

The `var` index uses UniProt protein accessions. Protein-level metadata (gene names) is stored as columns in `var`.

### X = primary quantification

The primary quantification matrix occupies `X`. For absolute expression, this is `ibaq_log`. The choice of log-transformed values as the primary matrix follows the scRNA-seq convention where `X` holds log-normalized counts.

### layers = alternative quantifications

Additional quantification types are stored as AnnData layers. Each layer has the same dimensions as `X` (samples x proteins).

### uns = file metadata + DE results

File-level metadata (QPX version, project accession, creation date) and differential expression results are stored in `uns`.

## File naming

Expression AnnData files follow the QPX naming convention:

| View | Extension | Example |
|------|-----------|---------|
| Absolute Expression | `.ae.h5ad` | `PXD000000.ae.h5ad` |
| Differential Expression | `.de.h5ad` | `PXD000000.de.h5ad` |

## Reading and writing

### Reading

```python
import anndata as ad

# Read AE data
adata = ad.read_h5ad("PXD000000.ae.h5ad")
print(adata)
# AnnData object with n_obs x n_vars = 120 x 5432
#     obs: 'organism', 'organism_part', 'disease', ...
#     var: 'gene_name'
#     layers: 'ibaq_raw', 'ibaq_ppb', 'copies_per_cell', 'concentration_nm'
```

### Writing

```python
adata.write("PXD000000.ae.h5ad")
```

### Using the QPX library API (proposed)

```python
# AE
ae = project.ae()
adata = ae.to_anndata(x_column="ibaq_log")

# DE
de = project.de()
adata.uns["de_results"] = de.to_de_results()

# Save combined
adata.write("PXD000000.ae.h5ad")
```

## Missing values

When a protein is not detected in a sample, the corresponding cell in `X` and `layers` contains `NaN`. This is consistent with how scRNA-seq handles dropout events and is the standard behavior when pivoting from long-form to matrix-form.

## Multi-omics integration

AnnData enables integration of proteomics data with other omics modalities:

- **CITE-seq + bulk proteomics**: Concatenate protein-level AnnData objects from different technologies.
- **muon**: Use `muon.MuData` to combine proteomics and transcriptomics AnnData objects.
- **Lamin.ai**: Register QPX AnnData output as a Lamin Artifact with schema validation and ontology-backed labels.

## MuData (`.h5mu`) {#mudata-h5mu}

Converters that produce expression output also emit a single combined **MuData** container, `<prefix>.h5mu`, that bundles the per-precursor and per-protein matrices (and the standalone expression views, when present) into one multi-modal object for the scverse / `muon` ecosystem. It is written by `qpx/mudata.py` (`build_mudata`).

### Modalities

| Modality (`mdata.mod[...]`) | `obs` (rows) | `var` (columns) | `X` |
|---|---|---|---|
| `precursors` | run (`run_file_name`, joined to sample metadata) | precursor identity `peptidoform\|charge` | sparse precursor intensity matrix |
| `proteins` | quantification unit (`unit`, the first of its runs, listed in `runs`; joined to sample metadata) | protein group, its sorted `pg_accessions` joined by `;` (with `anchor_protein` and `gene_name` in `var`) | sparse protein-group intensity matrix |
| `expression` | Absolute Expression observations | AE variables | Absolute Expression matrix (loaded from the AE view when present) |
| `differential` | Differential Expression observations | DE variables | Differential Expression matrix (loaded from the DE view when present) |

A modality is included only when its source data is available, so a given `.h5mu` may carry a subset of these four.

Converters write every label as its own observation (`run|label`, or `unit|label` for proteins). `build_mudata()` keeps label-free data at one observation per run (or unit) unless `all_intensity_labels=True`.

- **Units.** A protein group is quantified once for the runs in `grouped_runs` (the fractions of one sample or plex), so the `proteins` modality has one observation per unit, not a copy of its value per fraction. Run-level columns (`run_accession`, `fraction`) are left out of its `obs`.
- **Sample metadata.** `obs` carries the sample of each run and label from `run.samples` and its characteristics from `sample.parquet` (organism, disease, ...), when present.
- **Cells.** `X` is sparse: a cell that was not quantified is not stored. Decoys are left out. The features of one precursor in one run and label are summed into one cell, as MaxQuant sums the evidence rows of a peptide in a run.
- **Layers.** `X` holds the primary `intensity`, the tool's unnormalised quantity (MaxQuant `Intensity`, DIA-NN `PG.Quantity`, FragPipe razor `Intensity`; see the tool mappings in [Protein Group](pg.md#recommended-intensity-names) and [Feature](feature.md)). Each additional intensity (`additional_intensities`, such as `maxlfq`, `lfq` or `ibaq`) is a layer named after its `intensity_name`, on the observations and variables of `X`. A cell the tool imputed is 1 in the layer `imputed_mask`, present when any is (`feature.imputed_labels`, `pg.is_imputed`).

### Cross-modality mapping

`mdata.varp["feature_mapping"]` is a sparse precursor↔protein adjacency matrix over the concatenated `var` axis: each non-zero entry links a `precursors` column to the `proteins` column(s) of the groups its features belong to, derived from `feature.pg_accessions`.

### Identity note

The `precursors` modality's `var_names` is the cross-run precursor identity `peptidoform|charge`. This intentionally collapses the per-run `feature_id` rows of the row-level `feature` view (the same precursor measured across multiple runs) into a single matrix column. `feature_id` remains the primary key of the `feature` view; the MuData `var` axis is a coarser, analysis-oriented identity. The `qpx_version` is stored under `mdata.uns` (see [Versioning](versioning.md)).

## Further reading

- [Absolute Expression](absolute.md) -- AE-specific AnnData schema and examples
- [Differential Expression](differential.md) -- DE-specific AnnData schema and examples
- [AnnData documentation](https://anndata.readthedocs.io/) -- Official AnnData reference
- [scanpy](https://scanpy.readthedocs.io/) -- Analysis framework for AnnData
