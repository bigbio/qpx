# quantms.io

[![Python application](https://github.com/bigbio/quantms.io/actions/workflows/python-app.yml/badge.svg?branch=dev)](https://github.com/bigbio/quantms.io/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/quantms.io/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/quantms.io/actions/workflows/python-publish.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/quantms.io/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/quantms.io/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_Coverage)
[![PyPI version](https://badge.fury.io/py/quantmsio.svg)](https://badge.fury.io/py/quantmsio)

A Python package for working with mass spectrometry data in the quantms.io format.

## Features

- Convert data from various mass spectrometry formats to quantms.io format
- Analyze and process quantms.io data
- Visualize results
- Manage project metadata
- Transform data between different formats

## Installation

### Install from PyPI

```bash
# To install the stable release from PyPI:
pip install quantmsio
```

### Install from Source (Without PyPI)

```bash
# Fork the repository on GitHub
# Clone the repository
git clone https://github.com/your-username/quantms.io.git
cd quantms.io

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

  quantmsio - A tool for converting and analyzing mass spectrometry proteomics
  data

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  convert    Convert external formats to quantms.io format.
  project    Project management commands.
  stats      Statistical analysis of quantms.io data.
  transform  Transform quantms.io data into different representations.
  visualize  Visualize quantms.io data.
```

### Convert Commands

Convert data from various external formats to quantms.io:

```bash
Usage: convert [OPTIONS] COMMAND [ARGS]...

  Convert external formats to quantms.io format.

Options:
  --help  Show this message and exit.

Commands:
  diann             Convert DIA-NN report to quantms.io format
  diann-pg          Convert DIA-NN report to protein group format
  fragpipe          Convert FragPipe PSMs from psm.tsv to parquet file in
                    quantms.io
  idxml             Convert IdXML to PSM parquet file in quantms io
  idxml-batch       Convert multiple IdXML files to a single merged PSM parquet
                    file
  maxquant-feature  Convert feature data from MaxQuant evidence.txt to parquet
                    format
  maxquant-pg       Convert MaxQuant proteinGroups.txt to quantms.io protein
                    group format
  maxquant-psm      Convert PSM data from MaxQuant msms.txt to parquet format
  quantms-feature   Convert feature data from mzTab to quantms.io format.
  quantms-pg        Convert protein groups from mzTab quantms TMT and LFQ...
  quantms-psm       Convert PSM data from mzTab to quantms.io format.
```

### Transform Commands

Transform data within the quantms.io ecosystem:

```bash
Usage: transform [OPTIONS] COMMAND [ARGS]...

  Transform quantms.io data into different representations.

Options:
  --help  Show this message and exit.

Commands:
  ae            Convert IBAQ absolute file into quantms.io format
  anndata       Merge multiple AE files into a file in AnnData format.
  differential  Convert a MSstats differential file into a quantms.io file
                format
  gene          Map gene information from FASTA to parquet format
  ibaq          Convert feature data to IBAQ format
  spectra       Map spectrum information from mzML to parquet format
  uniprot       Map feature data to latest UniProt version
```

### Visualization Commands

Visualize quantms.io data:

```bash
Usage: visualize [OPTIONS] COMMAND [ARGS]...

  Visualize quantms.io data.

Options:
  --help  Show this message and exit.

Commands:
  plot  Visualization commands for quantms.io data
```

### Statistics Commands

Analyze quantms.io data:

```bash
Usage: stats [OPTIONS] COMMAND [ARGS]...

  Statistical analysis of quantms.io data.

Options:
  --help  Show this message and exit.

Commands:
  analyze  Statistical analysis commands for quantms.io data
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
quantmsio/
├── __init__.py
├── quantmsioc.py           # CLI entry point (poetry script: quantmsioc)
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

    This information is free; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This information is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this work; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
