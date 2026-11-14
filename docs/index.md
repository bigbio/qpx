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

`qpxc convert` turns tool outputs into a QPX dataset: one folder of `{prefix}.{view}.parquet` files (PSMs, precursor features, protein groups, spectra, fragments, peptide-to-protein mappings, proteins, modification sites, proteoforms and metadata) and optional HDF5 analysis matrices (`.h5mu`, `.ae.h5ad`, `.de.h5ad`). Available views depend on the input tool. Datasets are read with SQL (DuckDB), the Python API, scverse tools or any Parquet reader, and maintained with `qpxc transform` and `qpxc migrate`.

### Performance

![QPX Benchmark](images/qpx-benchmark.svg)

The figure uses the LFQ and TMT outputs of [PXD007683](https://www.ebi.ac.uk/pride/archive/projects/PXD007683). Queries compare a gzip mzTab PSM section read with pandas against the QPX PSM view read with DuckDB; each query ran in five fresh processes with a warm OS page cache, and the figure shows medians, with imports, conversion and `relayout` outside the timed interval. The text workflow parses the whole PSM section for every query, while DuckDB reads only the needed Parquet columns and answers an unfiltered count from the footer. Memory is the increase in peak RSS above the imported process. Storage compares related outputs after `qpxc transform relayout`, which sorts records, remaps IDs and reduces the precision of selected measured values, so the sizes are not a lossless compression ratio. Measurements used an Intel Core Ultra 7 270K Plus (24 logical CPUs, 122.9 GiB RAM), Python 3.13.9, DuckDB 1.5.2, pandas 2.3.3 and PyArrow 24; results depend on the data, query and machine.

---

## Supported Input Formats

| Software | PSM | Feature | Protein Group |
| --- | --- | --- | --- |
| **MaxQuant** | msms.txt | evidence.txt | proteinGroups.txt |
| **DIA-NN** | - | report.tsv | pg_matrix.tsv |
| **Spectronaut** | - | report.tsv | report.tsv (PG.Quantity) |
| **FragPipe** | psm.tsv | combined_ion.tsv | combined_protein.tsv |
| **Proteome Discoverer** | PSMs.txt | PSMs.txt | Proteins.txt |
| **OpenMS native QPX** | psm.parquet | feature.parquet | pg.parquet |
| **OpenMS consensusXML** | consensusXML | consensusXML | consensusXML |
| **CPTAC CDAP** | .psm | .psm reporter ions | feature-derived |
| **mzIdentML** | .mzid / .mzid.gz | - | - |
| **QuantMS MSstats** | - | `*_msstats_in.csv` | - |
| **SDRF** | - | - | - (sample + run) |

All conversions produce standardized Parquet and AnnData files following the [QPX specification](spec/index.md). `qpxc convert mz` reads spectra from mzML files into the `mz` view.

---

## Documentation

| Section | Description |
| --- | --- |
| **[Quick Start](quickstart.md)** | Installation and basic usage |
| **[Format Specification](spec/index.md)** | Data format schemas and specifications |
| **[User Guide](guide/index.md)** | CLI reference and command documentation |
| **[Examples & Tutorials](examples/index.md)** | Usage examples and integrations |
| **[Troubleshooting](troubleshooting.md)** | Common issues and solutions |
| **[Community & Support](community.md)** | Get help, contribute, license & acknowledgments |

---

**Ready to get started?** Install qpx and [check out the examples](examples/index.md)!
