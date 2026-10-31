# File Extensions & Naming

QPX defines a consistent file naming convention that encodes the project identity, data view, and serialization format directly in the filename. This enables both humans and automated tools to identify file contents without reading metadata.

## Convention

All QPX files follow a single naming pattern:

```text
{PREFIX}.{view}.{format}
```

Where:

- **`{PREFIX}`** -- A short, unique identifier for the dataset. See [Choosing a prefix](#choosing-a-prefix) below.
- **`{view}`** -- One of the QPX spec-defined views (e.g., `psm`, `feature`, `pg`, `mz`, `ae`, `de`, `dataset`, `sample`, `run`, `ontology`, `provenance`).
- **`{format}`** -- The serialization format extension (`parquet`, `h5ad`, `h5mu`, `tsv`).

### Choosing a prefix

The prefix identifies your dataset. Choose one based on your context:

| Context | Recommended prefix | Example |
|---------|-------------------|---------|
| PRIDE / ProteomeXchange submission | ProteomeXchange accession | `PXD014414` |
| Local research project | Short project or experiment name | `my_phospho_study` |
| Lab notebook / internal | Lab identifier or experiment code | `EXP2024_001` |
| Multi-condition study | Descriptive short name | `heart_proteome_aging` |

The only requirement is that the prefix is consistent across all files in the same dataset. Avoid spaces and special characters — use underscores or hyphens.

## Examples

### PRIDE / ProteomeXchange dataset

| File Name | View | Format |
|-----------|------|--------|
| `PXD014414.psm.parquet` | PSM | Parquet |
| `PXD014414.feature.parquet` | Feature | Parquet |
| `PXD014414.pg.parquet` | Protein Group | Parquet |
| `PXD014414.mz.parquet` | Mass Spectra | Parquet |
| `PXD014414.ae.h5ad` | Absolute Expression | AnnData |
| `PXD014414.de.h5ad` | Differential Expression | AnnData |
| `PXD014414.h5mu` | Combined expression (MuData) | MuData |
| `PXD014414.dataset.parquet` | Dataset Metadata | Parquet |
| `PXD014414.sample.parquet` | Sample Metadata | Parquet |
| `PXD014414.run.parquet` | Run Metadata | Parquet |
| `PXD014414.ontology.parquet` | Ontology Mapping | Parquet |
| `PXD014414.provenance.parquet` | Processing Provenance | Parquet |
| `PXD014414.sdrf.tsv` | Original SDRF | TSV |

### Local research project

| File Name | View | Format |
|-----------|------|--------|
| `my_phospho_study.psm.parquet` | PSM | Parquet |
| `my_phospho_study.feature.parquet` | Feature | Parquet |
| `my_phospho_study.pg.parquet` | Protein Group | Parquet |
| `my_phospho_study.mz.parquet` | Mass Spectra | Parquet |
| `my_phospho_study.ae.h5ad` | Absolute Expression | AnnData |
| `my_phospho_study.de.h5ad` | Differential Expression | AnnData |
| `my_phospho_study.dataset.parquet` | Dataset Metadata | Parquet |
| `my_phospho_study.sample.parquet` | Sample Metadata | Parquet |
| `my_phospho_study.run.parquet` | Run Metadata | Parquet |
| `my_phospho_study.ontology.parquet` | Ontology Mapping | Parquet |
| `my_phospho_study.provenance.parquet` | Processing Provenance | Parquet |

!!! info "Why preserve the original SDRF?"
    `qpxc convert` copies the SDRF it was given to `{prefix}.sdrf.tsv`. The `.sdrf.tsv` file is kept alongside `sample.parquet` and `run.parquet` to maintain a direct link to the original experimental metadata as submitted to ProteomeXchange. For local projects without an SDRF, this file can be omitted -- the `sample.parquet` and `run.parquet` views contain all the necessary metadata.

## Complete Project Layout

A typical QPX project directory contains:

=== "PRIDE dataset"

    ```
    PXD014414/
      # Metadata
      PXD014414.dataset.parquet          # Project-level metadata
      PXD014414.sample.parquet           # Biological samples
      PXD014414.run.parquet              # Data acquisition runs
      PXD014414.ontology.parquet         # Field-to-ontology mapping
      PXD014414.provenance.parquet       # Processing chain & parameters
      PXD014414.sdrf.tsv                 # Original SDRF (provenance)

      # Data views
      PXD014414.psm.parquet              # Peptide spectrum matches
      PXD014414.feature.parquet          # Quantified peptide features
      PXD014414.pg.parquet               # Protein groups
      PXD014414.mz.parquet               # Mass spectra

      # Expression views
      PXD014414.ae.h5ad                  # Absolute expression (AnnData)
      PXD014414.de.h5ad                  # Differential expression (AnnData)
    ```

=== "Local project"

    ```
    my_phospho_study/
      # Metadata
      my_phospho_study.dataset.parquet   # Project-level metadata
      my_phospho_study.sample.parquet    # Biological samples
      my_phospho_study.run.parquet       # Data acquisition runs
      my_phospho_study.ontology.parquet  # Field-to-ontology mapping
      my_phospho_study.provenance.parquet # Processing chain & parameters

      # Data views
      my_phospho_study.psm.parquet       # Peptide spectrum matches
      my_phospho_study.feature.parquet   # Quantified peptide features
      my_phospho_study.pg.parquet        # Protein groups
      my_phospho_study.mz.parquet        # Mass spectra

      # Expression views
      my_phospho_study.ae.h5ad           # Absolute expression (AnnData)
      my_phospho_study.de.h5ad           # Differential expression (AnnData)
    ```

!!! tip "Versioning multiple analysis runs"
    If a project is reprocessed with different tools or parameters, use subdirectories or an external artifact management system (e.g., lamindb) to distinguish analysis runs -- not filename conventions. Each project directory should contain a single set of QPX files representing one analysis.

## See Also

- [Serialization](serialization.md) -- serialization formats and Parquet details
- [QPX Format Overview](index.md) -- overview of all QPX views and their relationships
