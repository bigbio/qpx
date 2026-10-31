# QPX Format Overview

The **QPX format** is a modern, scalable data format designed specifically for proteomics data analysis. It addresses the limitations of existing formats -- XML-based HUPO-PSI standards (mzML, mzIdentML) and tab-delimited formats like mzTab -- which struggle with large-scale datasets and advanced analytical use cases such as AI/ML model development.

QPX organizes proteomics results into multiple **views**, each capturing a different aspect of the data (identifications, quantifications, spectra, metadata). Data views use Apache Parquet as the primary serialization format. Expression views (absolute and differential expression) use AnnData (`.h5ad`) for native interoperability with the scverse ecosystem.

!!! note
    QPX does not aim to replace mzTab or individual tool output formats. Its goal is to provide a unified, performance-oriented format that enables AI-related use cases and easy integration of proteomics results across tools and platforms.

## Data model

The diagram below shows the QPX views and how they relate to each other. Arrows indicate data flow from raw spectra through identification and quantification to final expression results. For the full conceptual model — every term defined, entity-relationship and quantification-unit diagrams, and how a measurement flows from raw file to protein quantity — see the [Data Model](data-model.md).

```mermaid
graph LR
    MZ[mz<br/>Mass Spectra] --> PSM[psm<br/>Peptide Spectrum Matches]
    PSM --> FEAT[feature<br/>Peptide Features]
    FEAT --> PG[pg<br/>Protein Groups]
    FEAT --> API["API Views<br/>(on demand)"]
    PG --> API
    PG --> AE[ae<br/>Absolute Expression]
    PG --> DE[de<br/>Differential Expression]
    PSM -.-> PPM[pepmap<br/>Peptide-Protein Map]
    PSM --> SITE[site<br/>Modification Sites]
    PG -.-> SITE
    PSM --> PF[proteoform<br/>Proteoforms]
    PG -.-> PF
    FEAT -.-> PPM
    SAMPLE[sample<br/>Sample Metadata] -.-> FEAT
    SAMPLE -.-> PG
    SAMPLE -.-> AE
    RUN[run<br/>Run Metadata] -.-> FEAT
    RUN -.-> PG
    RUN -.-> API
    DATASET[dataset<br/>Dataset Metadata] -.-> MZ
    DATASET -.-> PSM
    DATASET -.-> FEAT

    style MZ fill:#e1f5fe
    style PSM fill:#e8f5e9
    style FEAT fill:#e8f5e9
    style PG fill:#fff3e0
    style PPM fill:#e8f5e9
    style SITE fill:#fff3e0
    style PF fill:#fff3e0
    style API fill:#e8f5e9,stroke-dasharray: 5 5
    style AE fill:#fce4ec
    style DE fill:#fce4ec
    style SAMPLE fill:#f3e5f5
    style RUN fill:#f3e5f5
    style DATASET fill:#f3e5f5
```

## Views at a glance

### Data Views

| View           | File pattern                    | Format    | Description                                                                 |
| -------------- | ------------------------------- | --------- | --------------------------------------------------------------------------- |
| `mz`           | `{PREFIX}.mz.parquet`           | Parquet   | Raw and processed mass spectra (MS1 and MS2)                                |
| `psm`          | `{PREFIX}.psm.parquet`          | Parquet   | Peptide spectrum matches from database search engines (primarily DDA)       |
| `feature`      | `{PREFIX}.feature.parquet`      | Parquet   | Quantified peptide features per MS run, with sample and channel intensities |
| `pg`           | `{PREFIX}.pg.parquet`           | Parquet   | Protein groups with quantification across runs and channels                 |
| `pepmap` | `{PREFIX}.pepmap.parquet` | Parquet | Deduplicated peptide-to-protein mapping with gene names and uniqueness flags |
| `site` | `{PREFIX}.site.parquet` | Parquet | Protein modification-site claims with a q-value over the claims, the probability that the residue is modified and the evidence fraction of each candidate modification |
| `protein` | `{PREFIX}.protein.parquet` | Parquet | Name, organism and sequence of each protein the dataset names |
| `proteoform` | `{PREFIX}.proteoform.parquet` | Parquet | Canonical, isoform, variant and chain members of protein families with member-level q-values (sequence entries, not modification combinations) |

### API Views

These are **not** stored as standalone files. They are programmable, on-demand views that the API computes by joining and aggregating the primary data views. Users can request different fields and aggregations depending on their needs (e.g., protein counts, peptide abundances, identification summaries). See [API Views](views.md) for details.

### Expression Views

| View           | File pattern                    | Format    | Description                                                                 |
| -------------- | ------------------------------- | --------- | --------------------------------------------------------------------------- |
| `ae`           | `{PREFIX}.ae.h5ad`              | AnnData   | Absolute expression matrix (iBAQ values per protein per sample)             |
| `de`           | `{PREFIX}.de.h5ad`              | AnnData   | Differential expression results (fold changes, p-values between contrasts)  |

### Metadata

| View           | File pattern                       | Format    | Description                                                              |
| -------------- | ---------------------------------- | --------- | ------------------------------------------------------------------------ |
| `dataset`      | `{PREFIX}.dataset.parquet`         | Parquet   | Project-level metadata, search parameters, and software provenance       |
| `sample`       | `{PREFIX}.sample.parquet`          | Parquet   | Biological sample metadata (one row per sample)                          |
| `run`          | `{PREFIX}.run.parquet`             | Parquet   | Data acquisition run metadata (one row per run)                          |
| `ontology`     | `{PREFIX}.ontology.parquet`        | Parquet   | Field-to-ontology mapping (makes the dataset self-describing)            |
| `provenance`   | `{PREFIX}.provenance.parquet`      | Parquet   | Processing chain: tools, versions, parameters, and FDR thresholds        |
| `sdrf`         | `{PREFIX}.sdrf.tsv`                | TSV       | Original SDRF file, preserved for provenance                             |

!!! tip
    File extensions follow the pattern `{PREFIX}.{view}.{format}` -- for example, `PXD014414.feature.parquet` or `PXD014414.sample.parquet`. See [File Extensions & Naming](file-naming.md) for the full convention.

## Which view do I need?

Use the decision guide below to find the right starting point for your use case.

| Use case                                          | Recommended view(s)                  |
| ------------------------------------------------- | ------------------------------------ |
| Train an ML model on PSM-level spectral features  | `psm` (with spectral arrays) + `mz` |
| Retrieve peptide intensities per sample            | `feature`                            |
| Build a volcano plot                               | `de`                                 |
| Compare protein abundance across conditions        | `ae` or `pg`                         |
| Look up which proteins a peptide maps to           | [`pepmap`](pepmap.md), `feature`, or `psm` |
| Get a quick protein or peptide report for a web portal | [API Views](views.md)            |
| Understand the experimental design                 | `sample` + `run`                     |
| Get project-level metadata and search parameters   | `dataset`                            |
| Integrate with scanpy / scverse ecosystem          | `.ae.h5ad` + `.de.h5ad`             |
| List all files in a QPX project                    | `dataset` or glob `{PREFIX}-*.parquet` |

!!! note
    Some views are better suited for specific acquisition methods. The `psm` view is primarily intended for DDA experiments. The `feature` view covers both DDA and DIA workflows.

## Shared concepts

Several data structures are reused across multiple views. Each concept has its own page with full definitions and examples.

| Concept                                | Used in views                      | Description                                                        |
| -------------------------------------- | ---------------------------------- | ------------------------------------------------------------------ |
| [Peptidoform](peptidoform.md)          | PSM, Feature                       | Peptide sequence with modifications in ProForma notation           |
| [Modifications](modifications.md)      | PSM, Feature                       | Structured modification representation with localization scores    |
| [Intensities](intensities.md)          | Feature, Protein Group             | Primary and additional intensity measurements across channels      |
| [Scores & CV Terms](scores.md)         | PSM, Feature, Protein Group        | Search engine scores and controlled vocabulary parameters          |
| [Scan Numbers](scan.md)               | PSM, Feature, MZ                   | Instrument-specific scan identifiers and format conventions        |

## Current status

The QPX format is at **version 1.3**, primarily implemented in the [quantms workflow](https://github.com/bigbio/quantms). The format is open and can be adopted by any software tool. During the pre-2.0 stabilisation period, minor versions may include documented breaking changes; readers must check each file's `qpx_version` and the [versioning changelog](versioning.md). Every view file includes a `qpx_version` metadata field to identify which specification version produced it.
