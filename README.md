# QPX

[![Python application](https://github.com/bigbio/qpx/actions/workflows/python-app.yml/badge.svg?branch=dev)](https://github.com/bigbio/qpx/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/qpx/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/qpx/actions/workflows/python-publish.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/qpx/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/qpx/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_Coverage)
[![PyPI version](https://badge.fury.io/py/qpx.svg)](https://badge.fury.io/py/qpx)

A Python package for working with mass spectrometry data in the QPX format.

![QPX Architecture](docs/images/qpx-architecture.svg)

## Features

- **Convert** data from OpenMS (native QPX and consensusXML), DIA-NN, MaxQuant, Spectronaut, FragPipe, CPTAC CDAP (`.psm`), mzIdentML, and SDRF to QPX Parquet format
- **Transform** QPX data: gene mapping, protein quantification (DirectLFQ, MaxLFQ, iBAQ, TopN, …), accession normalization, metadata updates
- **Query** datasets with SQL, filter rows, or preview with `head`
- **Inspect** dataset summaries, Arrow schemas, and Parquet metadata
- **Validate** datasets against the canonical QPX schema
- **Ontology** management for PSI-MS and PRIDE CV terms

### MuData Export (quantification results only)

QPX datasets can be exported to [MuData](https://mudata.readthedocs.io/) — the multi-modal container from the scverse ecosystem. This export is available **only for quantification results** (precursor/protein intensities and, optionally, protein expression and differential expression results).

![QPX MuData Structure](docs/images/mudata-architecture.svg)

```python
from qpx import Dataset
from qpx.mudata import build_mudata

ds = Dataset("path/to/PXD000000/")
mdata = build_mudata(ds)  # auto-detects label & available modalities
mdata.write("PXD000000.h5mu")  # serialize to HDF5
```

Use `build_mudata(ds, all_intensity_labels=True)` to represent every TMT/iTRAQ
channel as a separate `run|label` observation.

> MuData support is included in the default QPX installation.

### Performance

![QPX Benchmark](docs/images/qpx-benchmark.svg)

## Installation

### Install from PyPI

```bash
pip install qpx                # includes MuData export (mudata, anndata, scipy)

# Optional extras
pip install "qpx[quantify]"    # protein quantification via mokume[directlfq]
pip install "qpx[transforms]"  # gene mapping / BioPython helpers
pip install "qpx[plotting]"    # plotting dependencies
pip install "qpx[mzidentml]"   # mzIdentML conversion (lxml)
pip install "qpx[pdc]"         # PDC/CPTAC download (pridepy)
pip install "qpx[all]"         # all optional extras above
```

### Install from GitHub (latest dev)

```bash
pip install git+https://github.com/bigbio/qpx.git
```

### Install from Source

```bash
# Clone the repository
git clone https://github.com/bigbio/qpx.git
cd qpx

# Install the package locally
pip install .
```

### Install and build with uv

[uv](https://docs.astral.sh/uv/) is a fast Python package installer and resolver. The project supports PEP 621 and can be installed, built, and published with uv.

**Prerequisites:** [Install uv](https://docs.astral.sh/uv/getting-started/installation/) (e.g. `curl -LsSf https://astral.sh/uv/install.sh | sh` or `pip install uv`).

```bash
# Install from GitHub
uv pip install "qpx @ git+https://github.com/bigbio/qpx.git"

# With optional extras
uv pip install "qpx[quantify,transforms,plotting] @ git+https://github.com/bigbio/qpx.git"
```

**From a local clone:**

```bash
git clone https://github.com/bigbio/qpx.git
cd qpx

# Create a venv, install the project and its dependencies (recommended)
uv sync

# Or install in editable mode with optional dev dependencies
uv sync --extra dev

# Run the CLI without installing globally
uv run qpxc --help
```

**Build distributable packages** (sdist and wheel in `dist/`):

```bash
uv build
```

**Publish to PyPI** (after configuring credentials or trusted publishing):

```bash
uv build
uv publish
```

The `pyproject.toml` uses PEP 621 metadata with Hatchling as the build backend.

### Development Installation

For development with all dependencies:

```bash
# Using uv (recommended for fast installs)
uv sync --extra dev

# Or using pip
pip install -e ".[dev]"
```

### System Dependencies

QPX depends on pyOpenMS, which requires certain system libraries. If you encounter errors related to missing shared libraries (e.g., `libglib-2.0.so.0`), install the required system dependencies:

**Ubuntu/Debian:**

```bash
sudo apt-get update
sudo apt-get install -y libglib2.0-0
```

**macOS:**

```bash
brew install glib
```

**Using Conda/Mamba (Recommended for pyOpenMS):**

Using mamba (faster dependency resolution):

```bash
mamba env create -f environment.yml
conda activate qpx
pip install git+https://github.com/bigbio/qpx.git
```

Or with conda:

```bash
conda env create -f environment.yml
conda activate qpx
pip install git+https://github.com/bigbio/qpx.git
```

## Usage

The package provides a command-line interface (`qpxc`) with the following command groups:

```bash
qpxc [OPTIONS] COMMAND [ARGS]...

Commands:
  convert    Convert external tool outputs to QPX format.
  transform  Transform QPX data into derived representations.
  query      Query and inspect QPX datasets.
  info       Show information about a QPX dataset.
  validate   Validate a QPX dataset or structure against the canonical schema.
  ontology   Manage CV ontology data (PSI-MS, PRIDE CV).
  pdc2qpx    Download a PDC/CPTAC study and convert it to a QPX dataset.
```

### pdc2qpx

One-shot download + conversion of PDC/CPTAC studies (requires `pip install qpx[pdc]`). `-a/--accession` takes a single ID, comma-separated IDs, or a CSV with a `pdc_study_id`/`pdc_id` column:

```bash
# Entire QPX including full spectra (for quantms reanalysis)
qpxc pdc2qpx -a PDC000109 \
    --download-dir ./downloads \
    --output-folder ./qpx/PDC000109 \
    --include-spectra --max-cpus 24

# Many studies from a CSV (each -> ./qpx/<study>/)
qpxc pdc2qpx -a cptac_lfq.csv --download-dir ./downloads --output-folder ./qpx
```

### Convert

```bash
qpxc convert [openms | openms-consensus | diann | maxquant | spectronaut | fragpipe | mzidentml | cdap | mz | sdrf] [OPTIONS]
```

### Transform

```bash
qpxc transform [gene-map | quantify | normalize-accessions | update-metadata] [OPTIONS]
```

### Query

```bash
# Run SQL against a dataset
qpxc query sql --dataset-path ./PXD014414 --sql "SELECT anchor_protein, COUNT(*) FROM feature GROUP BY 1"

# Filter rows
qpxc query filter --dataset-path ./PXD014414 --structure feature --condition "charge >= 3"

# Preview first N rows
qpxc query head --dataset-path ./PXD014414 --structure feature -n 20
```

### Info & Validate

```bash
# Dataset summary
qpxc info --dataset-path ./PXD014414

# Validate against canonical schema
qpxc validate --dataset-path ./PXD014414
```

## Configuration

Most commands support a `--verbose` flag that enables more detailed logging to stdout. The CLI uses standard logging configuration and does not require environment variables.

## Development

### Project Structure

```text
qpx/
├── cli/                    # Click CLI (entry point: qpx.cli.main:main)
│   ├── main.py             # Top-level CLI group
│   ├── pdc2qpx.py          # pdc2qpx command (PDC download + convert)
│   └── convert.py          # convert subcommands (openms, openms-consensus, maxquant, diann, spectronaut, fragpipe, mzidentml, cdap, mz, sdrf)
├── pipeline/               # High-level orchestration (pdc2qpx: download + CDAP + mz)
├── converters/             # Tool-specific converters
│   ├── openms/             # OpenMS native QPX enrichment
│   ├── openms_consensus/   # OpenMS consensusXML converter
│   ├── cdap/               # CPTAC CDAP (.psm) converter
│   ├── diann/              # DIA-NN converter
│   ├── maxquant/           # MaxQuant converter
│   ├── spectronaut/        # Spectronaut converter
│   ├── fragpipe/           # FragPipe converter
│   ├── mzidentml/          # mzIdentML converter
│   └── sdrf.py             # Shared SDRF converter
├── core/                   # Core logic & formats
│   ├── data/               # Schema definitions (YAML + Python)
│   │   └── schemas/        # YAML schema files for all structures
│   ├── engine.py           # DuckDB engine wrapper
│   ├── scores.py           # Score normalization & ontology
│   └── ontology/           # OBO ontology registry
├── writers/                # Parquet writers (one per structure)
├── views/                  # Analytical views (protein, peptide, QC)
└── dataset.py              # Main Dataset class entry point
```

### Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests
5. Submit a pull request

## License

This project is licensed under the Apache-2.0 License - see the LICENSE file for details.

## Core contributors and collaborators

The project is run by different groups:

- Yasset Perez-Riverol (PRIDE Team, European Bioinformatics Institute - EMBL-EBI, U.K.)
- Ping Zheng (Chongqing Key Laboratory of Big Data for Bio Intelligence, Chongqing University of Posts and Telecommunications, Chongqing, China)

IMPORTANT: If you contribute with the following specification, please make sure to add your name to the list of contributors.

## Code of Conduct

As part of our efforts toward delivering open and inclusive science, we follow the [Contributor Covenant Code of Conduct for Open Source Projects](https://www.contributor-covenant.org/version/2/0/code_of_conduct/).

## How to cite

## Copyright notice

```text
Copyright 2025 BigBio

Licensed under the Apache License, Version 2.0.
See the LICENSE file for details.
```
