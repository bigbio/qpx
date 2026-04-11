# QPX

A standardized format and toolkit for mass spectrometry proteomics data

[![PyPI version](https://badge.fury.io/py/qpx.svg)](https://badge.fury.io/py/qpx)
[![Python application](https://github.com/bigbio/qpx/actions/workflows/python-app.yml/badge.svg?branch=dev)](https://github.com/bigbio/qpx/actions/workflows/python-app.yml)

---

## What is QPX?

QPX is a comprehensive ecosystem for proteomics data that provides:

- **Standardized Data Format**: Parquet-based format for efficient storage and processing of proteomics data
- **Universal Converter**: Convert data from MaxQuant, DIA-NN, FragPipe, quantms, mzIdentML, and more
- **Complete Toolkit**: Process, analyze, visualize, and share your proteomics results
- **Python API & CLI**: Flexible tools for both programmatic and command-line usage

---

## Architecture Overview

![QPX Architecture](images/qpx-architecture.svg)

---

## Quick Start

### Installation

=== "pip (Recommended)"

    ```bash
    pip install qpx
    ```

=== "conda"

    ```bash
    conda create -n qpx python=3.10
    conda activate qpx
    pip install qpx
    ```

=== "From Source"

    ```bash
    git clone https://github.com/bigbio/qpx.git
    cd qpx
    pip install .
    ```

### Verify Installation

```bash
qpxc --version
qpxc --help
```

### Your First Conversion

```bash
# Convert MaxQuant PSM data to QPX parquet
qpxc convert maxquant \
    --msms-file msms.txt \
    --output-folder ./output \
    --structures psm

# Convert DIA-NN data
qpxc convert diann \
    --report-path report.tsv \
    --sdrf-file data.sdrf.tsv \
    --output-folder ./output

# Query the output
qpxc query sql \
    --dataset-path ./output \
    --sql "SELECT * FROM psm LIMIT 10"
```

### Inspect with Python

```python
import pyarrow.parquet as pq

table = pq.read_table("output/psm.parquet")
df = table.to_pandas()
print(f"Total PSMs: {len(df)}")
print(df.head())
```

---

## Supported Input Formats

| Software      | PSM              | Feature              | Protein Group        |
| ------------- | ---------------- | -------------------- | -------------------- |
| **MaxQuant**  | msms.txt         | evidence.txt         | proteinGroups.txt    |
| **DIA-NN**    | -                | report.tsv           | pg_matrix.tsv        |
| **FragPipe**  | psm.tsv          | combined_peptide.tsv | combined_protein.tsv |
| **quantms**   | mzTab PSM section| mzTab + MSstats      | mzTab PRT section    |
| **mzIdentML** | .mzid / .mzid.gz | -                    | -                    |
| **SDRF**      | -                | -                    | - (sample + run)     |

All conversions produce standardized Parquet and AnnData files following the [QPX specification](spec/index.md).

---

## Ecosystem

QPX is part of the [quantms ecosystem](https://quantms.org):

| Tool | Description |
|------|-------------|
| [quantms](https://github.com/bigbio/quantms) | DDA proteomics Nextflow pipeline |
| [mokume](https://mokume.quantms.org) | Protein quantification library |
| [pmultiqc](https://pmultiqc.quantms.org) | Interactive QC reporting |
| [portal.quantms.org](https://portal.quantms.org) | Browse reanalyzed datasets |
