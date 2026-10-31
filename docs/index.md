# QPX

A standardized format and toolkit for mass spectrometry proteomics data

[![Python application](https://github.com/bigbio/qpx/actions/workflows/python-app.yml/badge.svg?branch=dev)](https://github.com/bigbio/qpx/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/qpx/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/qpx/actions/workflows/python-publish.yml)
[![PyPI version](https://badge.fury.io/py/qpx.svg)](https://badge.fury.io/py/qpx)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/qpx/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)

---

## What is QPX?

QPX is a comprehensive ecosystem for proteomics data that provides:

- **Standardized Data Format**: Parquet-based format for efficient storage and processing of proteomics data
- **Universal Converter**: Convert data from OpenMS, MaxQuant, DIA-NN, Spectronaut, FragPipe, Proteome Discoverer, CPTAC CDAP, mzIdentML, and more
- **Complete Toolkit**: Process, analyze, visualize, and share your proteomics results
- **Python API & CLI**: Flexible tools for both programmatic and command-line usage

---

## Architecture Overview

![QPX Architecture](images/qpx-architecture.svg)

QPX provides a comprehensive proteomics data processing architecture with core modules for data conversion, transformation, visualization, statistical analysis, and project management.

### Performance

![QPX Benchmark](images/qpx-benchmark.svg)

---

## Supported Input Formats

| Software                | PSM              | Feature              | Protein Group            |
| ----------------------- | ---------------- | -------------------- | ------------------------ |
| **MaxQuant**            | msms.txt         | evidence.txt         | proteinGroups.txt        |
| **DIA-NN**              | -                | report.tsv           | pg_matrix.tsv            |
| **Spectronaut**         | -                | report.tsv           | report.tsv (PG.Quantity) |
| **FragPipe**            | psm.tsv          | combined_ion.tsv     | combined_protein.tsv     |
| **Proteome Discoverer** | PSMs.txt         | PSMs.txt             | Proteins.txt         |
| **OpenMS native QPX**   | psm.parquet      | feature.parquet      | pg.parquet               |
| **OpenMS consensusXML** | consensusXML     | consensusXML         | consensusXML             |
| **CPTAC CDAP**          | .psm             | .psm reporter ions   | feature-derived          |
| **mzIdentML**           | .mzid / .mzid.gz | -                    | -                        |
| **SDRF**                | -                | -                    | - (sample + run)         |

All conversions produce standardized Parquet and AnnData files following the [QPX specification](spec/index.md).

---

## Documentation

| Section                                             | Description                                     |
| --------------------------------------------------- | ----------------------------------------------- |
| **[Quick Start](quickstart.md)**                    | Installation and basic usage                    |
| **[Format Specification](spec/index.md)**           | Data format schemas and specifications          |
| **[User Guide](guide/index.md)**                    | CLI reference and command documentation         |
| **[Examples & Tutorials](examples/index.md)**       | Usage examples and integrations                 |
| **[Troubleshooting](troubleshooting.md)**           | Common issues and solutions                     |
| **[Community & Support](community.md)**             | Get help, contribute, license & acknowledgments |

---

**Ready to get started?** Install qpx and [check out the examples](examples/index.md)!
