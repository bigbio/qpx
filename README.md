# QPX

[![Python application](https://github.com/bigbio/QPX/actions/workflows/python-app.yml/badge.svg?branch=dev)](https://github.com/bigbio/QPX/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/QPX/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/QPX/actions/workflows/python-publish.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/QPX/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/QPX/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_Coverage)
[![PyPI version](https://badge.fury.io/py/qpx.svg)](https://badge.fury.io/py/qpx)

A Python package for working with mass spectrometry data in the QPX format.

## Features

- Convert data from various mass spectrometry formats to QPX format
- Analyze and process QPX data
- Visualize results
- Manage project metadata
- Transform data between different formats

## Installation

### Install from PyPI

```bash
# To install the stable release from PyPI:
pip install qpx
```

### Install from Source (Without PyPI)

```bash
# Fork the repository on GitHub
# Clone the repository
git clone https://github.com/your-username/QPX.git
cd QPX

# Install the package locally
pip install .
```

### Development Installation

For development with all dependencies:

```bash
# Using Poetry (recommended)
poetry install

# Or using pip
pip install -r requirements.txt
pip install -e .
```

## Usage

The package provides a command-line interface (CLI) with several command groups:

### Main CLI

```bash
Usage: cli [OPTIONS] COMMAND [ARGS]...

  qpx - A tool for converting and analyzing mass spectrometry proteomics
  data

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  convert    Convert external formats to QPX format.
  project    Project management commands.
  stats      Statistical analysis of QPX data.
  transform  Transform QPX data into different representations.
  visualize  Visualize QPX data.
```

### Convert Commands

Convert data from various external formats to QPX:

```bash
Usage: convert [OPTIONS] COMMAND [ARGS]...

  Convert external formats to QPX format.

Options:
  --help  Show this message and exit.

Commands:
  diann             Convert DIA-NN report to QPX format
  diann-pg          Convert DIA-NN report to protein group format
  fragpipe          Convert FragPipe PSMs from psm.tsv to parquet file in
                    QPX
  idxml             Convert IdXML to PSM parquet file in QPX
  idxml-batch       Convert multiple IdXML files to a single merged PSM parquet
                    file
  maxquant-feature  Convert feature data from MaxQuant evidence.txt to parquet
                    format
  maxquant-pg       Convert MaxQuant proteinGroups.txt to QPX protein
                    group format
  maxquant-psm      Convert PSM data from MaxQuant msms.txt to parquet format
  quantms-feature   Convert feature data from mzTab to QPX format.
  quantms-pg        Convert protein groups from mzTab quantms TMT and LFQ...
  quantms-psm       Convert PSM data from mzTab to QPX format.
```

### Transform Commands

Transform data within the QPX ecosystem:

```bash
Usage: transform [OPTIONS] COMMAND [ARGS]...

  Transform QPX data into different representations.

Options:
  --help  Show this message and exit.

Commands:
  ae            Convert IBAQ absolute file into QPX format
  anndata       Merge multiple AE files into a file in AnnData format.
  differential  Convert a MSstats differential file into a QPX file
                format
  gene          Map gene information from FASTA to parquet format
  ibaq          Convert feature data to IBAQ format
  spectra       Map spectrum information from mzML to parquet format
  uniprot       Map feature data to latest UniProt version
```

### Visualization Commands

Visualize QPX data:

```bash
Usage: visualize [OPTIONS] COMMAND [ARGS]...

  Visualize QPX data.

Options:
  --help  Show this message and exit.

Commands:
  plot  Visualization commands for QPX data
```

### Statistics Commands

Analyze QPX data:

```bash
Usage: stats [OPTIONS] COMMAND [ARGS]...

  Statistical analysis of QPX data.

Options:
  --help  Show this message and exit.

Commands:
  analyze  Statistical analysis commands for QPX data
```

### Project Management Commands

Manage project metadata:

```bash
Usage: project [OPTIONS] COMMAND [ARGS]...

  Project management commands.

Options:
  --help  Show this message and exit.

Commands:
  attach  Register the file to project.json.
  create  Generate a project file from original PRIDE accession
```

## Configuration

Most commands support a `--verbose` flag that enables more detailed logging to stdout. The CLI uses standard logging configuration and does not require environment variables.

## Development

### Project Structure

```
qpx/
├── __init__.py
├── qpxc.py                 # CLI entry point (poetry script: qpxc)
├── commands/               # CLI command groups
│   ├── convert/            # Converters: quantms, maxquant, diann, idxml, fragpipe
│   ├── transform/          # Transforms: ibaq, ae, gene, spectra, anndata, differential, uniprot
│   └── utils/              # Utility CLIs: project(create/attach), stats(analyze), plot
├── core/                   # Core logic & formats
│   ├── quantms/            # quantms feature/psm/pg, mztab helpers
│   ├── diann/, maxquant/, fragpipe/, idxml_utils/ ...
│   └── project.py, duckdb.py, format.py, common.py
├── operate/                # High-level operations (stats, plotting, tools)
│   ├── plots.py, query.py, statistics.py, tools.py
│   └── ...
└── utils/                  # Utilities
    ├── logger.py           # Basic logger getter
    ├── file_utils.py       # File helpers (e.g., AE file discovery)
    ├── pride_utils.py      # PRIDE archive helpers
    ├── mztab_utils.py      # mzTab helpers
    ├── system.py           # System utilities
    └── constants.py        # Constants and configurations
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

    Copyright 2025 BigBio
    
    Licensed under the Apache License, Version 2.0.
    See the LICENSE file for details.