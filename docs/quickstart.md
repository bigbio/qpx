# Quick Start

## Quick Start Flow

![Quick Start Flow](images/workflow2.png)

Get started with quantmsio in a few simple steps: from installation to data conversion, transformation, and visualization.

## Installation

```bash
pip install quantmsio
```

## Basic Usage

```bash
# View available commands
quantmsioc --help

# Convert MaxQuant data
quantmsioc convert maxquant-psm \
    --msms-file msms.txt \
    --output-folder ./output

# Transform to absolute expression
quantmsioc transform ae \
    --ibaq-file ibaq.tsv \
    --sdrf-file metadata.sdrf.tsv \
    --output-folder ./output

# Generate visualizations
quantmsioc visualize plot ibaq-distribution \
    --ibaq-path ./output/ae.parquet \
    --save-path ./plots/distribution.svg
```

[View complete CLI documentation →](cli-reference.md)
