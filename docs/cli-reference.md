# CLI Reference

The qpx command-line tool provides a comprehensive set of commands for converting, transforming, visualizing, and analyzing mass spectrometry proteomics data.

## Overview

The qpx CLI is organized into five main command groups:

### [Convert Commands](cli-convert.md)

Convert various mass spectrometry data formats to the QPX standard format:

- **DIA-NN Conversion**: Convert DIA-NN reports to feature and protein group formats
- **MaxQuant Conversion**: Convert MaxQuant PSM, feature, and protein group data
- **FragPipe Conversion**: Convert FragPipe PSM data to QPX format
- **QuantMS Conversion**: Create QPX data files from mzTab format
- **idXML Conversion**: Convert OpenMS idXML format PSM data

[View detailed documentation →](cli-convert.md)

### [Transform Commands](cli-transform.md)

Transform and process data within the QPX ecosystem:

- **Absolute Expression (AE)**: Convert iBAQ absolute expression data to QPX format ([format specification](https://io.quantms.org/format-specification/#absolute))
- **Differential Expression (DE)**: Convert MSstats differential expression analysis results
- **Gene Mapping**: Map gene information to protein data
- **iBAQ Transformation**: Process iBAQ quantification files
- **Spectra Mapping**: Map spectra information
- **UniProt Mapping**: Update UniProt annotation information
- **AnnData Merging**: Merge multiple AE files into AnnData format

[View detailed documentation →](cli-transform.md)

### [Visualization Commands](cli-visualize.md)

Create various data visualization plots:

- PSM/Peptide distribution plots
- iBAQ distribution plots
- Intensity distribution plots (KDE, box plots)
- Peptide-protein distribution plots

[View detailed documentation →](cli-visualize.md)

### [Statistics Commands](cli-stats.md)

Perform statistical analysis on QPX data:

- Project AE data statistics
- PSM data statistics
- Generate comprehensive statistical reports

[View detailed documentation →](cli-stats.md)

### [Project Management Commands](cli-project.md)

Manage project metadata and files:

- Create project.json from PRIDE projects
- Attach files to project metadata

[View detailed documentation →](cli-project.md)

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
# 1. Convert raw data
qpxc convert maxquant-psm \
    --msms-file msms.txt \
    --output-folder ./output

# 2. Transform to absolute expression data
qpxc transform ae \
    --ibaq-file ibaq.tsv \
    --sdrf-file metadata.sdrf.tsv \
    --output-folder ./output

# 3. Generate visualization
qpxc visualize plot ibaq-distribution \
    --ibaq-path ./output/ae.parquet \
    --save-path ./plots/distribution.svg

# 4. Generate statistical report
qpxc stats analyze psm \
    --parquet-path ./output/psm.parquet \
    --save-path ./stats/report.txt
```

## Getting Help

- Each command provides detailed help information using the `--help` parameter
- See [Format Specification](format-specification.md) for output file formats
- View the [online format specification](https://io.quantms.org/format-specification/) for detailed schema information
- Visit the [GitHub Repository](https://github.com/bigbio/qpx) to report issues

