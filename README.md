# QPX

[![Python application](https://github.com/bigbio/qpx/actions/workflows/python-app.yml/badge.svg?branch=dev)](https://github.com/bigbio/qpx/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/qpx/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/qpx/actions/workflows/python-publish.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/qpx/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/qpx/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_Coverage)
[![PyPI version](https://badge.fury.io/py/qpx.svg)](https://badge.fury.io/py/qpx)

A Python package for converting, querying and validating quantitative proteomics data in the QPX format. QPX stores linked records as Apache Parquet views and exports quantification matrices to MuData for the scverse ecosystem.

![QPX Architecture](docs/images/qpx-architecture.svg)

## Features

- **Convert** data from OpenMS/quantms (native QPX, consensusXML and MSstats input), DIA-NN, MaxQuant, Spectronaut, FragPipe, Proteome Discoverer, CPTAC CDAP (`.psm`), mzIdentML, mzML and SDRF to QPX
- **Transform** QPX data: gene mapping, protein quantification (DirectLFQ, MaxLFQ, iBAQ, TopN, Top3 and sum), accession normalization, protein properties, site quantification and metadata updates
- **Query** datasets with SQL or the Python API, filter rows, or preview with `head`; export results to pandas, Arrow or Polars
- **Inspect** dataset summaries, Arrow schemas, and Parquet metadata
- **Validate** schemas, values and links between views against the canonical QPX specification
- **Prepare for publication** with `transform relayout` (sorting, ID remapping and numeric precision reduction), or upgrade older datasets with `migrate`
- **Ontology** management for PSI-MS and PRIDE CV terms

Available views depend on the input tool. The [data model](docs/spec/data-model.md) covers PSMs, precursor features, protein groups, spectra, fragments, peptide-to-protein mappings, proteins, modification sites, proteoforms and metadata. A dataset keeps these views in one folder as `{prefix}.{view}.parquet` files; optional analysis matrices use HDF5.

### MuData Export (quantification results only)

QPX datasets can be exported to [MuData](https://mudata.readthedocs.io/) — the multi-modal container from the scverse ecosystem. This export is available **only for quantification results** (precursor/protein intensities and, optionally, protein expression and differential expression results).

![QPX MuData Structure](docs/images/mudata-architecture.svg)

```python
from qpx import Dataset
from qpx.mudata import build_mudata

ds = Dataset("path/to/PXD000000/")
mdata = build_mudata(ds)  # detects available modalities and expands multiplexed labels
mdata.write("PXD000000.h5mu")  # serialize to HDF5
```

By default, multiplexed data (TMT, iTRAQ and plexDIA) expand over every channel:
`run|label` for precursor observations and `grouped_runs|label` for protein-group
observations. Label-free data keep one observation per run or quantification
unit; `all_intensity_labels=True` forces the expanded form. Use
`intensity_label="TMT126"` to select one channel.

Protein observations are quantification units: the fractions in one
`grouped_runs` set share one protein-group quantity per label. Variables are
`peptidoform|charge` for precursors and sorted `pg_accessions` joined by `;` for
protein groups. Additional intensities become layers. When both modalities are
present, `varp["feature_mapping"]` links the precursor and protein-group variables.
Decoys are left out by default. Repeated values in one cell are summed, or
averaged with `duplicates="mean"`. Null intensities are omitted from the sparse
matrices; their implicit cells read as zero. See [MuData](docs/spec/anndata.md#mudata-h5mu)
for `include_decoys=True`, imputation masks and optional expression modalities.

> MuData support is included in the default QPX installation.

### Performance

![QPX Benchmark](docs/images/qpx-benchmark.svg)

Local measurements use the LFQ and TMT outputs of [PXD007683](https://www.ebi.ac.uk/pride/archive/projects/PXD007683). The LFQ query comparison uses a gzip mzTab PSM section read with pandas and the corresponding QPX PSM view read with DuckDB. Each query ran in five fresh processes with a warm OS page cache; the table shows medians. Imports, conversion and `relayout` are outside the timed interval.

| LFQ PSM query | gzip mzTab + pandas | QPX + DuckDB | Matching result |
| --- | ---: | ---: | --- |
| Count all PSMs | 2.361 s | 1.71 ms | 542,403 |
| Count PSMs with `charge = 2` | 2.332 s | 2.53 ms | 382,341 |
| Read `sequence` and `charge` | 2.320 s | 26.43 ms | 542,403 rows; identical multisets |

This compares two reading workflows and engines. The text workflow parses the full PSM section for each query; DuckDB can read only the needed Parquet columns and answer an unfiltered count from the footer. The figure's memory values are increases in process peak RSS above the imported process, not total process memory. Measurements used an Intel Core Ultra 7 270K Plus, 24 logical CPUs, 122.9 GiB RAM, Python 3.13.9, DuckDB 1.5.2, pandas 2.3.3 and PyArrow 24, with 24 DuckDB threads. Results depend on the data, query and machine.

Storage compares related outputs after `qpxc transform relayout`: LFQ gzip mzTab + MSstats occupy **59.59 MiB**, versus **23.74 MiB** for six QPX Parquet views; TMT gzip MSstats occupies **18.17 MiB**, versus **5.13 MiB** for the QPX feature view alone. LFQ feature rows differ between MSstats and QPX (518,248 and 573,203). TMT stores 1,362,592 channel entries in 123,872 feature rows with 11 intensities each. These are different representations and contents. `relayout` also sorts records, remaps IDs and reduces the precision of selected measured values, so these sizes are not a lossless compression ratio.

## Installation

### Install from PyPI

```bash
pip install qpx                # includes MuData export (mudata, anndata, scipy)

# Optional extras
pip install "qpx[quantify]"    # protein quantification via mokume[directlfq]
pip install "qpx[plotting]"    # plotting dependencies
pip install "qpx[mzidentml]"   # mzIdentML conversion (lxml)
pip install "qpx[pdc]"         # PDC/CPTAC download (pridepy)
pip install "qpx[all]"         # all optional extras above
```

### Install from GitHub (latest dev)

```bash
pip install git+https://github.com/bigbio/qpx.git@dev
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
uv pip install "qpx[quantify,plotting] @ git+https://github.com/bigbio/qpx.git"
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
  migrate    Upgrade a dataset to the current QPX specification.
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
qpxc convert COMMAND [OPTIONS]
# openms, openms-consensus, quantms-msstats, diann, diann-xic, maxquant, spectronaut,
# fragpipe, proteome-discoverer, mzidentml, cdap, mz, sdrf
```

Use `qpxc convert <command> --help` for the required input files and available output views. The `openms` command enriches native OpenMS QPX output and is deprecated in favor of `openms-consensus` for consensusXML input.

### Transform

```bash
qpxc transform COMMAND [OPTIONS]
# gene-map, quantify, normalize-accessions, protein-properties,
# site, update-metadata, relayout

# Write a sorted copy for selective reads and publication
qpxc transform relayout --dataset ./PXD014414 --output-folder ./PXD014414_sorted

# Upgrade an older dataset
qpxc migrate --help
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
│   ├── convert/            # Tool commands, shared options and postprocessing
│   ├── transform/          # Transform commands
│   ├── query.py            # SQL, filter and head commands
│   └── migrate.py          # Dataset migration command
├── pipeline/               # High-level orchestration (pdc2qpx: download + CDAP + mz)
├── converters/             # Tool-specific converters
│   ├── openms/             # OpenMS native QPX enrichment
│   ├── openms_consensus/   # OpenMS consensusXML converter
│   ├── cdap/               # CPTAC CDAP (.psm) converter
│   ├── diann/              # DIA-NN converter
│   ├── maxquant/           # MaxQuant converter
│   ├── spectronaut/        # Spectronaut converter
│   ├── fragpipe/           # FragPipe converter
│   ├── proteome_discoverer/ # Proteome Discoverer exports converter
│   ├── mzidentml/          # mzIdentML converter
│   ├── quantms_msstats/    # QuantMS *_msstats_in.csv + SDRF converter
│   └── sdrf.py             # Shared SDRF converter
├── core/                   # Core logic & formats
│   ├── data/               # Schema definitions (YAML + Python)
│   │   └── schemas/        # YAML schema files for all structures
│   ├── engine.py           # DuckDB engine wrapper
│   ├── scores.py           # Score normalization & ontology
│   └── ontology/           # OBO ontology registry
├── writers/                # Parquet writers (one per structure)
├── views/                  # Analytical views (protein, peptide, QC)
├── transforms/             # Gene mapping, quantification, relayout and metadata
├── collection.py           # Multi-dataset queries
├── mudata.py               # Quantification matrices and MuData export
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
