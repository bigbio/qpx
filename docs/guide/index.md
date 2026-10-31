# CLI Reference

The qpx command-line tool provides a comprehensive set of commands for converting, transforming, querying, and validating mass spectrometry proteomics data.

## Overview

The qpx CLI is organized into five main command groups:

### [Convert Commands](convert.md)

Convert various mass spectrometry data formats to the QPX standard format:

- **OpenMS Conversion**: Enrich native `-out_qpx` output or convert consensusXML to QPX data files
- **DIA-NN Conversion**: Convert DIA-NN reports to feature and protein group formats
- **MaxQuant Conversion**: Convert MaxQuant PSM, feature, and protein group data
- **FragPipe Conversion**: Convert FragPipe PSM data to QPX format
- **mzIdentML Conversion**: Convert mzIdentML format PSM data
- **SDRF Conversion**: Convert SDRF metadata to QPX sample/run format

[View detailed documentation →](convert.md)

### [Transform Commands](transform.md)

Transform and process data within the QPX ecosystem:

- **Gene Mapping**: Map gene information from FASTA to protein data
- **Protein Quantification**: Compute protein-level quantification via mokume (DirectLFQ, MaxLFQ, iBAQ, TopN, etc.)

[View detailed documentation →](transform.md)

### Query Commands

Query QPX datasets using SQL or structured commands:

- Execute SQL queries on QPX Parquet files
- Filter and aggregate data
- Join multiple datasets
- Export query results

### Info Commands

Inspect QPX datasets and metadata:

- Display dataset schema and statistics
- Show sample and run information
- List available columns and data types
- Validate data integrity

### Validate Commands

Validate QPX data against schemas:

- Check data conformance to QPX specifications
- Validate required fields and data types
- Verify referential integrity
- Generate validation reports

## Python API

Core data operations are available through the Python API:

```python
import qpx

# Load a dataset
with qpx.open_dataset("path/to/dataset") as ds:
    # Access data views
    psm_df = ds.psm.to_df()
    feature_df = ds.feature.to_df()

    # Filter and query
    targets = ds.psm.targets_only().to_df()
    count = ds.feature.count()

    # Validate against canonical schema
    results = ds.validate()
    for name, result in results.items():
        print(f"{name}: {result.summary}")
        for issue in result.issues:
            print(f"  [{issue.severity}] {issue.message}")

    # Run SQL queries
    df = ds.sql("SELECT anchor_protein, COUNT(*) AS n FROM feature GROUP BY 1")
```

## Quick Start

### Installation

```bash
pip install qpx
```

### Basic Usage

View all available commands:

```bash
qpxc --help
```

View help for a specific command group:

```bash
qpxc convert --help
```

View detailed help for a specific command:

```bash
qpxc convert diann --help
```

## Common Options

Most commands support the following common options:

- `--verbose`: Enable verbose logging for debugging
- `--output-folder`: Specify the output directory
- `--output-prefix`: Specify the output file prefix

## Example Data Processing Workflow

A typical data processing workflow:

```bash
# 1. Convert raw data from MaxQuant
qpxc convert maxquant \
    --evidence-file evidence.txt \
    --msms-file msms.txt \
    --output-folder ./output

# 2. Protein quantification (DirectLFQ)
qpxc transform quantify \
    --feature-path ./output/diann.feature.parquet \
    --method directlfq \
    -o ./output/proteins.parquet

# 3. Query the data
qpxc query sql \
    --dataset-path ./output \
    --sql "SELECT anchor_protein, AVG(intensity) FROM feature GROUP BY anchor_protein" \
    --output results.csv

# 4. Validate the dataset
qpxc validate \
    --dataset-path ./output

# 5. Inspect dataset information
qpxc info --dataset-path ./output
qpxc info schema --dataset-path ./output --structure feature
```

For further analysis, use the Python API:

```python
import qpx

# Load and inspect the dataset
with qpx.open_dataset("./output") as ds:
    # Validate the dataset
    results = ds.validate()
    for name, result in results.items():
        print(f"{name}: {result.summary}")

    # Query with SQL
    top_proteins = ds.sql("SELECT anchor_protein, COUNT(*) AS n FROM feature GROUP BY 1 ORDER BY n DESC LIMIT 10")
    print(top_proteins)

    # Access data views as DataFrames
    feature_df = ds.feature.to_df()
    print(f"Features: {len(feature_df)} rows")
```

## Getting Help

- Each command provides detailed help information using the `--help` parameter
- See [Format Specification](../spec/index.md) for output file formats and detailed schema information
- Check the [Troubleshooting Guide](../troubleshooting.md) for common issues
- Visit the [GitHub Repository](https://github.com/bigbio/qpx) to report issues
