# Transform Commands

Transform and process data within the quantms.io ecosystem.

## Overview

The `transform` command group provides tools for processing and transforming quantms.io data into various downstream formats. These commands enable absolute and differential expression analysis, metadata mapping, and data format conversions.

## Available Commands

- [ae](#ae) - Convert iBAQ to absolute expression format
- [differential](#differential) - Convert MSstats differential expression data
- [gene](#gene) - Map gene information to proteins
- [ibaq](#ibaq) - Process iBAQ quantification files
- [spectra](#spectra) - Map spectrum information
- [uniprot](#uniprot) - Map latest UniProt annotations
- [anndata](#anndata) - Merge AE files into AnnData format

---

## ae

Convert iBAQ absolute expression data to quantms.io format.

### Description {#ae-description}

This command transforms iBAQ (intensity-Based Absolute Quantification) data into the standardized quantms.io absolute expression format. It integrates protein quantification with sample metadata from SDRF files.

### Parameters {#ae-parameters}

| Parameter           | Type   | Required | Default | Description                             |
| ------------------- | ------ | -------- | ------- | --------------------------------------- |
| `--ibaq-file`       | Path   | Yes      | -       | iBAQ file path                          |
| `--sdrf-file`       | Path   | Yes      | -       | SDRF file path for metadata             |
| `--protein-file`    | Path   | No       | -       | Protein file with specific requirements |
| `--project-file`    | Path   | No       | -       | quantms.io project file                 |
| `--output-folder`   | Path   | Yes      | -       | Output directory for generated files    |
| `--output-prefix`   | String | No       | -       | Prefix for output files                 |
| `--delete-existing` | Flag   | No       | False   | Delete existing files in output folder  |

### Usage Examples {#ae-examples}

#### Basic Example {#ae-example-basic}

```bash
quantmsioc transform ae \
    --ibaq-file tests/examples/AE/PXD016999.1-ibaq.tsv \
    --sdrf-file tests/examples/AE/PXD016999-first-instrument.sdrf.tsv \
    --output-folder ./output
```

#### With Project Metadata {#ae-example-metadata}

```bash
quantmsioc transform ae \
    --ibaq-file tests/examples/AE/PXD016999.1-ibaq.tsv \
    --sdrf-file tests/examples/AE/PXD016999-first-instrument.sdrf.tsv \
    --project-file tests/examples/AE/project.json \
    --output-folder ./output \
    --output-prefix ae_with_metadata \
    --delete-existing
```

#### Filter Specific Proteins {#ae-example-filter}

```bash
quantmsioc transform ae \
    --ibaq-file tests/examples/AE/PXD016999.1-ibaq.tsv \
    --sdrf-file tests/examples/AE/PXD016999-first-instrument.sdrf.tsv \
    --protein-file tests/examples/fasta/Homo-sapiens.fasta \
    --output-folder ./output
```

### Input File Formats {#ae-input-formats}

**iBAQ File**: Tab-separated file with protein accessions and iBAQ intensities

```
ProteinID    Sample1    Sample2    Sample3
P12345       1000000    950000     1050000
Q67890       500000     480000     520000
```

**SDRF File**: Standard PRIDE SDRF format with sample metadata

### Output Files {#ae-output}

- **Output**: `{output-prefix}-{uuid}.absolute.parquet`
- **Format**: Parquet file containing absolute expression quantification
- **Schema**: Conforms to quantms.io absolute expression specification

### Common Issues {#ae-issues}

**Issue**: Mismatched sample names between iBAQ and SDRF

- **Solution**: Ensure column names in iBAQ file match sample identifiers in SDRF

**Issue**: Missing protein accessions

- **Solution**: Provide `--protein-file` to filter and validate protein IDs

### Best Practices {#ae-best-practices}

- Always provide project metadata file when available for better data provenance
- Use `--delete-existing` flag carefully to avoid accidental data loss
- Validate SDRF file format before processing
- Check sample name consistency across input files

---

## differential

Convert MSstats differential expression data to quantms.io format.

### Description {#differential-description}

Transforms differential expression analysis results from MSstats into the standardized quantms.io differential expression format. Supports FDR-based filtering and protein-specific subsetting.

### Parameters {#differential-parameters}

| Parameter           | Type   | Required | Default | Description                             |
| ------------------- | ------ | -------- | ------- | --------------------------------------- |
| `--msstats-file`    | Path   | Yes      | -       | MSstats differential expression file    |
| `--sdrf-file`       | Path   | Yes      | -       | SDRF file for metadata extraction       |
| `--project-file`    | Path   | No       | -       | quantms.io project file                 |
| `--protein-file`    | Path   | No       | -       | Protein file with specific requirements |
| `--fdr-threshold`   | Float  | No       | 0.05    | FDR threshold to filter results         |
| `--output-folder`   | Path   | Yes      | -       | Output directory for generated files    |
| `--output-prefix`   | String | No       | -       | Prefix for output files                 |
| `--delete-existing` | Flag   | No       | False   | Delete existing files in output folder  |
| `--verbose`         | Flag   | No       | False   | Enable verbose logging                  |

### Usage Examples {#differential-examples}

#### Basic Example {#differential-example-basic}

```bash
quantmsioc transform differential \
    --msstats-file tests/examples/DE/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv \
    --sdrf-file tests/examples/DE/PXD033169.sdrf.tsv \
    --output-folder ./output
```

#### With Custom FDR Threshold {#differential-example-fdr}

```bash
quantmsioc transform differential \
    --msstats-file tests/examples/DE/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv \
    --sdrf-file tests/examples/DE/PXD033169.sdrf.tsv \
    --fdr-threshold 0.01 \
    --output-folder ./output \
    --output-prefix de_stringent \
    --verbose
```

#### With Project Metadata {#differential-example-metadata}

```bash
quantmsioc transform differential \
    --msstats-file tests/examples/DE/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv \
    --sdrf-file tests/examples/DE/PXD033169.sdrf.tsv \
    --project-file tests/examples/DE/project.json \
    --fdr-threshold 0.05 \
    --output-folder ./output \
    --delete-existing
```

### Input File Format {#differential-input-format}

**MSstats File**: CSV file with comparison results

```
Protein,Label,log2FC,SE,Tvalue,DF,pvalue,adj.pvalue
P12345,Condition2-Condition1,2.5,0.3,8.33,10,0.0001,0.001
Q67890,Condition2-Condition1,-1.8,0.4,-4.5,10,0.002,0.01
```

### Output Files {#differential-output}

- **Output**: `{output-prefix}-{uuid}.differential.parquet`
- **Format**: Parquet file containing differential expression results
- **Schema**: Conforms to quantms.io differential expression specification

### Common Issues {#differential-issues}

**Issue**: No significant results after FDR filtering

- **Solution**: Increase `--fdr-threshold` or check input data quality

**Issue**: Memory errors with large comparison files

- **Solution**: Process comparisons in batches or increase available memory

### Best Practices {#differential-best-practices}

- Use FDR threshold of 0.05 or lower for publication-quality results
- Enable verbose mode to monitor filtering statistics
- Validate comparison group names match SDRF metadata
- Include project file for complete data provenance

---

## gene

Map gene information from FASTA to parquet format.

### Description {#gene-description}

Maps gene names and information from a FASTA file to protein identifications in quantms.io PSM or feature files. This command enriches protein data with gene-level metadata extracted from FASTA headers.

### Parameters {#gene-parameters}

| Parameter         | Type    | Required | Default | Description                                  |
| ----------------- | ------- | -------- | ------- | -------------------------------------------- |
| `--parquet-path`  | Path    | Yes      | -       | PSM or feature parquet file path             |
| `--fasta`         | Path    | Yes      | -       | FASTA file path                              |
| `--output-folder` | Path    | Yes      | -       | Output directory for generated files         |
| `--file-num`      | Integer | No       | 10      | Number of rows to read in each batch         |
| `--partitions`    | String  | No       | -       | Fields for splitting files (comma-separated) |
| `--species`       | String  | No       | `human` | Species name                                 |

### Usage Examples {#gene-examples}

#### Basic Example {#gene-example-basic}

```bash
quantmsioc transform gene \
    --parquet-path ./output/psm.parquet \
    --fasta tests/examples/fasta/Homo-sapiens.fasta \
    --output-folder ./output
```

#### With Partitioning {#gene-example-partition}

```bash
quantmsioc transform gene \
    --parquet-path ./output/feature.parquet \
    --fasta tests/examples/fasta/Homo-sapiens.fasta \
    --output-folder ./output \
    --file-num 20 \
    --partitions reference_file_name \
    --species human
```

### Output Files {#gene-output}

- **Output**: Enhanced parquet file(s) with gene information
- **Format**: Parquet file in output folder
- **Added Fields**: Gene names and metadata from FASTA headers

### Best Practices {#gene-best-practices}

- Use species-specific FASTA files for accurate gene annotation
- Adjust `--file-num` based on available memory for large files
- Use partitioning for better file organization in large datasets

---

## ibaq

Convert feature data to iBAQ format.

### Description {#ibaq-description}

Transforms feature-level quantification data into iBAQ (intensity-Based Absolute Quantification) format. This command integrates feature quantification with sample metadata from SDRF files to generate iBAQ values.

### Parameters {#ibaq-parameters}

| Parameter         | Type   | Required | Default | Description                          |
| ----------------- | ------ | -------- | ------- | ------------------------------------ |
| `--feature-file`  | Path   | Yes      | -       | Feature file path                    |
| `--sdrf-file`     | Path   | Yes      | -       | SDRF file for metadata extraction    |
| `--output-folder` | Path   | Yes      | -       | Output directory for generated files |
| `--output-prefix` | String | No       | `ibaq`  | Prefix for output files              |

### Usage Examples {#ibaq-examples}

#### Basic Example {#ibaq-example-basic}

```bash
quantmsioc transform ibaq \
    --feature-file ./output/feature.parquet \
    --sdrf-file ./metadata.sdrf.tsv \
    --output-folder ./output
```

#### With Custom Prefix {#ibaq-example-prefix}

```bash
quantmsioc transform ibaq \
    --feature-file ./output/feature.parquet \
    --sdrf-file ./metadata.sdrf.tsv \
    --output-folder ./output \
    --output-prefix ibaq_quantification
```

### Output Files {#ibaq-output}

- **Output**: `{output-prefix}-{uuid}.ibaq.parquet`
- **Format**: Parquet file containing iBAQ quantification values
- **Content**: Protein-level iBAQ values per sample

### Best Practices {#ibaq-best-practices}

- Ensure feature file contains all necessary quantification data
- Verify SDRF metadata matches sample identifiers in feature file
- Use iBAQ output for absolute protein quantification analysis

---

## spectra

Map spectrum information from mzML to parquet format.

### Description {#spectra-description}

Enriches PSM or feature data with additional spectral information extracted from mzML files. This command maps spectrum metadata and peak information to the corresponding peptide-spectrum matches.

### Parameters {#spectra-parameters}

| Parameter          | Type    | Required | Default | Description                                  |
| ------------------ | ------- | -------- | ------- | -------------------------------------------- |
| `--parquet-path`   | Path    | Yes      | -       | PSM or feature parquet file path             |
| `--mzml-directory` | Path    | Yes      | -       | Directory containing mzML files              |
| `--output-folder`  | Path    | Yes      | -       | Output directory for generated files         |
| `--file-num`       | Integer | No       | 10      | Number of rows to read in each batch         |
| `--partitions`     | String  | No       | -       | Fields for splitting files (comma-separated) |

### Usage Examples {#spectra-examples}

#### Basic Example {#spectra-example-basic}

```bash
quantmsioc transform spectra \
    --parquet-path ./output/psm.parquet \
    --mzml-directory ./mzml_files \
    --output-folder ./output
```

#### With Batch Processing and Partitioning {#spectra-example-batch}

```bash
quantmsioc transform spectra \
    --parquet-path ./output/psm.parquet \
    --mzml-directory ./mzml_files \
    --output-folder ./output \
    --file-num 20 \
    --partitions reference_file_name
```

### Output Files {#spectra-output}

- **Output**: Enhanced PSM/feature parquet file(s) with spectral information
- **Format**: Parquet file in output folder
- **Added Fields**: Spectrum metadata and peak information from mzML files

### Best Practices {#spectra-best-practices}

- Ensure mzML files are in the specified directory with correct naming
- Adjust `--file-num` based on available memory and file size
- Use partitioning for organized output when processing large datasets

---

## uniprot

Map feature data to latest UniProt version.

### Description {#uniprot-description}

Maps peptides and features to the latest UniProt protein database using a FASTA file. This command updates protein identifications to match current UniProt accessions and annotations.

### Parameters {#uniprot-parameters}

| Parameter         | Type   | Required | Default | Description                          |
| ----------------- | ------ | -------- | ------- | ------------------------------------ |
| `--feature-file`  | Path   | Yes      | -       | Feature file path                    |
| `--fasta`         | Path   | Yes      | -       | UniProt FASTA file path              |
| `--output-folder` | Path   | Yes      | -       | Output directory for generated files |
| `--output-prefix` | String | No       | -       | Prefix for output files              |

### Usage Examples {#uniprot-examples}

#### Basic Mapping {#uniprot-example-basic}

```bash
quantmsioc transform uniprot \
    --feature-file ./output/feature.parquet \
    --fasta tests/examples/fasta/Homo-sapiens.fasta \
    --output-folder ./output
```

#### With Custom Prefix {#uniprot-example-prefix}

```bash
quantmsioc transform uniprot \
    --feature-file ./output/feature.parquet \
    --fasta ./uniprot_human_2024.fasta \
    --output-folder ./output \
    --output-prefix feature_updated
```

### Output Files {#uniprot-output}

- **Output**: `{output-prefix}-{uuid}.feature.parquet`
- **Format**: Parquet file with updated UniProt mappings
- **Content**: Feature data mapped to latest UniProt protein identifications

### Best Practices {#uniprot-best-practices}

- Use the latest UniProt FASTA file for most current annotations
- Run this command when updating to a new UniProt release
- Verify FASTA file matches the organism of your study

---

## anndata

Merge multiple AE files into a file in AnnData format.

### Description {#anndata-description}

Combines multiple absolute expression (AE) files from a directory into a single AnnData object (H5AD format). This command is useful for integrating data from multiple experiments for downstream analysis with scanpy or other Python-based tools.

### Parameters {#anndata-parameters}

| Parameter         | Type   | Required | Default | Description                          |
| ----------------- | ------ | -------- | ------- | ------------------------------------ |
| `--directory`     | Path   | Yes      | -       | Directory for storing AE files       |
| `--output-folder` | Path   | Yes      | -       | Output directory for generated files |
| `--output-prefix` | String | No       | -       | Prefix for output files              |

### Usage Examples {#anndata-examples}

#### Basic Example {#anndata-example-basic}

```bash
quantmsioc transform anndata \
    --directory ./ae_files \
    --output-folder ./output
```

#### With Custom Prefix {#anndata-example-prefix}

```bash
quantmsioc transform anndata \
    --directory ./ae_files \
    --output-folder ./output \
    --output-prefix merged_ae
```

### Output Files {#anndata-output}

- **Output**: `{output-prefix}-{uuid}.h5ad`
- **Format**: AnnData H5AD file (HDF5-based format)
- **Structure**:
  - `X`: Protein expression matrix
  - `obs`: Sample metadata
  - `var`: Protein metadata

### Use Cases {#anndata-use-cases}

- Integration with scanpy for dimensionality reduction and clustering
- Compatibility with Python machine learning libraries
- Multi-experiment data integration
- Cross-study meta-analysis

### Best Practices {#anndata-best-practices}

- Ensure all AE files in the directory have consistent format
- Validate sample metadata consistency before merging
- Use the output with scanpy or other AnnData-compatible tools
- Consider memory requirements for large datasets

---

## Related Commands

- [Convert Commands](cli-convert.md) - Convert raw data to quantms.io format
- [Visualization Commands](cli-visualize.md) - Visualize transformed data
- [Statistics Commands](cli-stats.md) - Analyze transformed data
