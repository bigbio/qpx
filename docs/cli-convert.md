# Convert Commands

Convert various mass spectrometry data formats to the quantms.io standard format.

## Overview

The `convert` command group provides converters for multiple proteomics software outputs, enabling standardization of data formats for downstream analysis. All commands generate parquet-format output files following the quantms.io specification.

## Available Commands

- [diann](#diann) - Convert DIA-NN report to feature format
- [diann-pg](#diann-pg) - Convert DIA-NN report to protein group format
- [maxquant-psm](#maxquant-psm) - Convert MaxQuant PSM data
- [maxquant-feature](#maxquant-feature) - Convert MaxQuant feature data
- [maxquant-pg](#maxquant-pg) - Convert MaxQuant protein groups
- [fragpipe](#fragpipe) - Convert FragPipe PSM data
- [quantms-psm](#quantms-psm) - Convert mzTab to PSM format
- [quantms-feature](#quantms-feature) - Convert mzTab to feature format
- [quantms-pg](#quantms-pg) - Convert mzTab to protein group format
- [idxml](#idxml) - Convert single idXML file to PSM format
- [idxml-batch](#idxml-batch) - Convert multiple idXML files to merged PSM format

---

## diann

Convert DIA-NN report files to quantms.io feature format.

### Description {#diann-description}

This command converts DIA-NN (Data-Independent Acquisition by Neural Networks) report files to the standardized quantms.io parquet format. It processes peptide-level quantification data and can optionally partition the output based on specified fields.

### Parameters {#diann-parameters}

| Parameter             | Type    | Required | Default   | Description                                      |
| --------------------- | ------- | -------- | --------- | ------------------------------------------------ |
| `--report-path`       | Path    | Yes      | -         | DIA-NN report file path (typically `report.tsv`) |
| `--qvalue-threshold`  | Float   | Yes      | 0.05      | Q-value threshold for filtering peptides         |
| `--mzml-info-folder`  | Path    | Yes      | -         | Folder containing mzML info files                |
| `--sdrf-path`         | Path    | Yes      | -         | SDRF file path for metadata                      |
| `--output-folder`     | Path    | Yes      | -         | Output directory for generated files             |
| `--protein-file`      | Path    | No       | -         | Protein file with specific requirements          |
| `--output-prefix`     | String  | No       | `feature` | Prefix for output files                          |
| `--partitions`        | String  | No       | -         | Field(s) for splitting files (comma-separated)   |
| `--duckdb-max-memory` | String  | No       | -         | Maximum memory for DuckDB (e.g., "4GB")          |
| `--duckdb-threads`    | Integer | No       | -         | Number of threads for DuckDB                     |
| `--batch-size`        | Integer | No       | 100       | Number of files to process simultaneously        |
| `--verbose`           | Flag    | No       | False     | Enable verbose logging                           |

### Usage Examples {#diann-examples}

#### Basic Example {#diann-example-basic}

Convert a DIA-NN report with default settings:

```bash
quantmsioc convert diann \
    --report-path tests/examples/diann/small/diann_report.tsv \
    --qvalue-threshold 0.05 \
    --mzml-info-folder tests/examples/diann/small/mzml \
    --sdrf-path tests/examples/diann/small/PXD019909-DIA.sdrf.tsv \
    --output-folder ./output
```

#### Advanced Example with Partitioning {#diann-example-advanced}

Convert with file partitioning based on reference_file_name:

```bash
quantmsioc convert diann \
    --report-path tests/examples/diann/full/diann_report.tsv.gz \
    --qvalue-threshold 0.01 \
    --mzml-info-folder tests/examples/diann/full/mzml \
    --sdrf-path tests/examples/diann/full/PXD036609.sdrf.tsv \
    --output-folder ./output \
    --partitions reference_file_name \
    --duckdb-max-memory 8GB \
    --duckdb-threads 4 \
    --verbose
```

### Output Files {#diann-output}

- **Output**: `{output-prefix}-{uuid}.feature.parquet`
- **Format**: Parquet file containing feature-level quantification data
- **Schema**: Conforms to quantms.io feature specification

### Common Issues {#diann-issues}

**Issue**: Out of memory errors with large files

- **Solution**: Increase `--duckdb-max-memory` parameter (e.g., `8GB`, `16GB`)

**Issue**: Slow processing

- **Solution**: Increase `--duckdb-threads` to utilize more CPU cores

**Issue**: Missing mzML info files

- **Solution**: Ensure all mzML info TSV files are in the specified folder with correct naming

### Best Practices {#diann-best-practices}

- Use Q-value threshold of 0.05 or lower for high-confidence results
- Enable partitioning for large datasets to improve memory usage
- Use verbose mode during initial testing to diagnose issues
- Ensure SDRF file correctly matches sample names in DIA-NN report

---

## diann-pg

Convert DIA-NN report files to quantms.io protein group format.

### Description {#diann-pg-description}

This command converts DIA-NN protein quantification matrices to the standardized quantms.io protein group format. It combines information from the main report file and the protein group matrix.

### Parameters {#diann-pg-parameters}

| Parameter             | Type    | Required | Default | Description                               |
| --------------------- | ------- | -------- | ------- | ----------------------------------------- |
| `--report-path`       | Path    | Yes      | -       | DIA-NN report file path                   |
| `--pg-matrix-path`    | Path    | Yes      | -       | DIA-NN protein quantities table file path |
| `--sdrf-path`         | Path    | Yes      | -       | SDRF file path for metadata               |
| `--output-folder`     | Path    | Yes      | -       | Output directory for generated files      |
| `--output-prefix`     | String  | No       | `pg`    | Prefix for output files                   |
| `--duckdb-max-memory` | String  | No       | -       | Maximum memory for DuckDB (e.g., "4GB")   |
| `--duckdb-threads`    | Integer | No       | -       | Number of threads for DuckDB              |
| `--batch-size`        | Integer | No       | 100     | Number of files to process simultaneously |
| `--verbose`           | Flag    | No       | False   | Enable verbose logging                    |

### Usage Examples {#diann-pg-examples}

#### Basic Example {#diann-pg-example-basic}

```bash
quantmsioc convert diann-pg \
    --report-path tests/examples/diann/full/diann_report.tsv.gz \
    --pg-matrix-path tests/examples/diann/full/diann_report.pg_matrix.tsv \
    --sdrf-path tests/examples/diann/full/PXD036609.sdrf.tsv \
    --output-folder ./output
```

#### High-Performance Example {#diann-pg-example-performance}

```bash
quantmsioc convert diann-pg \
    --report-path tests/examples/diann/full/diann_report.tsv.gz \
    --pg-matrix-path tests/examples/diann/full/diann_report.pg_matrix.tsv \
    --sdrf-path tests/examples/diann/full/PXD036609.sdrf.tsv \
    --output-folder ./output \
    --duckdb-max-memory 16GB \
    --duckdb-threads 8 \
    --output-prefix protein_groups \
    --verbose
```

### Output Files {#diann-pg-output}

- **Output**: `{output-prefix}-{uuid}.pg.parquet`
- **Format**: Parquet file containing protein group quantification data
- **Schema**: Conforms to quantms.io protein group specification

### Best Practices {#diann-pg-best-practices}

- Ensure both report and pg_matrix files are from the same DIA-NN run
- Use adequate memory allocation for large datasets
- Validate SDRF metadata matches the sample columns in the matrix file

---

## maxquant-psm

Convert MaxQuant PSM data from `msms.txt` to quantms.io parquet format.

### Description {#maxquant-psm-description}

Converts MaxQuant Peptide-Spectrum Match (PSM) data to the standardized quantms.io format. This command processes the `msms.txt` file generated by MaxQuant.

### Parameters {#maxquant-psm-parameters}

| Parameter         | Type    | Required | Default | Description                          |
| ----------------- | ------- | -------- | ------- | ------------------------------------ |
| `--msms-file`     | Path    | Yes      | -       | MaxQuant msms.txt file path          |
| `--output-folder` | Path    | Yes      | -       | Output directory for generated files |
| `--batch-size`    | Integer | No       | 1000000 | Read batch size for processing       |
| `--output-prefix` | String  | No       | `psm`   | Prefix for output files              |
| `--spectral-data` | Flag    | No       | False   | Include spectral data fields         |
| `--verbose`       | Flag    | No       | False   | Enable verbose logging               |

### Usage Examples {#maxquant-psm-examples}

#### Basic Example {#maxquant-psm-example-basic}

```bash
quantmsioc convert maxquant-psm \
    --msms-file tests/examples/maxquant/maxquant_simple/msms.txt \
    --output-folder ./output
```

#### With Spectral Data {#maxquant-psm-example-spectral}

```bash
quantmsioc convert maxquant-psm \
    --msms-file tests/examples/maxquant/maxquant_simple/msms.txt \
    --output-folder ./output \
    --spectral-data \
    --batch-size 500000 \
    --output-prefix psm_with_spectra \
    --verbose
```

### Output Files {#maxquant-psm-output}

- **Output**: `{output-prefix}-{uuid}.psm.parquet`
- **Format**: Parquet file containing PSM-level data
- **Schema**: Conforms to quantms.io PSM specification

### Best Practices {#maxquant-psm-best-practices}

- Adjust `--batch-size` based on available memory
- Use `--spectral-data` flag if downstream analysis requires spectral information
- Ensure sufficient disk space for large msms.txt files

---

## maxquant-feature

Convert MaxQuant feature data from `evidence.txt` to quantms.io parquet format.

### Description {#maxquant-feature-description}

Converts MaxQuant feature-level quantification data to the standardized quantms.io format. This command processes the `evidence.txt` file, which contains peptide feature intensities across samples.

### Parameters {#maxquant-feature-parameters}

| Parameter               | Type    | Required | Default   | Description                                    |
| ----------------------- | ------- | -------- | --------- | ---------------------------------------------- |
| `--evidence-file`       | Path    | Yes      | -         | MaxQuant evidence.txt file path                |
| `--sdrf-file`           | Path    | Yes      | -         | SDRF file for metadata extraction              |
| `--output-folder`       | Path    | Yes      | -         | Output directory for generated files           |
| `--protein-file`        | Path    | No       | -         | Protein file with specific requirements        |
| `--protein-groups-file` | Path    | No       | -         | MaxQuant proteinGroups.txt for Q-value mapping |
| `--partitions`          | String  | No       | -         | Field(s) for splitting files (comma-separated) |
| `--batch-size`          | Integer | No       | 1000000   | Read batch size                                |
| `--output-prefix`       | String  | No       | `feature` | Prefix for output files                        |
| `--verbose`             | Flag    | No       | False     | Enable verbose logging                         |

### Usage Examples {#maxquant-feature-examples}

#### Basic Example {#maxquant-feature-example-basic}

```bash
quantmsioc convert maxquant-feature \
    --evidence-file tests/examples/maxquant/maxquant_full/evidence.txt.gz \
    --sdrf-file tests/examples/maxquant/maxquant_full/PXD001819.sdrf.tsv \
    --output-folder ./output
```

#### With Protein Groups Q-value Mapping {#maxquant-feature-example-qvalue}

```bash
quantmsioc convert maxquant-feature \
    --evidence-file tests/examples/maxquant/maxquant_full/evidence.txt.gz \
    --sdrf-file tests/examples/maxquant/maxquant_full/PXD001819.sdrf.tsv \
    --protein-groups-file tests/examples/maxquant/maxquant_full/proteinGroups.txt \
    --output-folder ./output \
    --batch-size 500000 \
    --verbose
```

### Output Files {#maxquant-feature-output}

- **Output**: `{output-prefix}-{uuid}.feature.parquet`
- **Format**: Parquet file containing feature-level quantification
- **Schema**: Conforms to quantms.io feature specification

### Common Issues {#maxquant-feature-issues}

**Issue**: Memory errors with compressed evidence files

- **Solution**: Reduce `--batch-size` or increase available RAM

**Issue**: Missing Q-value information

- **Solution**: Provide `--protein-groups-file` for accurate Q-value mapping

### Best Practices {#maxquant-feature-best-practices}

- Always provide `--protein-groups-file` when available for better data quality
- Ensure SDRF sample names match MaxQuant experiment names
- Use compressed evidence files (.gz) to save disk space

---

## maxquant-pg

Convert MaxQuant protein groups from `proteinGroups.txt` to quantms.io format.

### Description {#maxquant-pg-description}

Converts MaxQuant protein group quantification to the standardized quantms.io protein group format.

### Parameters {#maxquant-pg-parameters}

| Parameter               | Type    | Required | Default | Description                          |
| ----------------------- | ------- | -------- | ------- | ------------------------------------ |
| `--protein-groups-file` | Path    | Yes      | -       | MaxQuant proteinGroups.txt file      |
| `--sdrf-file`           | Path    | Yes      | -       | SDRF file for metadata extraction    |
| `--output-folder`       | Path    | Yes      | -       | Output directory for generated files |
| `--batch-size`          | Integer | No       | 1000000 | Batch size (for logging purposes)    |
| `--output-prefix`       | String  | No       | `pg`    | Prefix for output files              |
| `--verbose`             | Flag    | No       | False   | Enable verbose logging               |

### Usage Examples {#maxquant-pg-examples}

#### Basic Example {#maxquant-pg-example-basic}

```bash
quantmsioc convert maxquant-pg \
    --protein-groups-file tests/examples/maxquant/maxquant_full/proteinGroups.txt \
    --sdrf-file tests/examples/maxquant/maxquant_full/PXD001819.sdrf.tsv \
    --output-folder ./output
```

### Output Files {#maxquant-pg-output}

- **Output**: `{output-prefix}-{uuid}.pg.parquet`
- **Format**: Parquet file containing protein group data
- **Schema**: Conforms to quantms.io protein group specification

---

## fragpipe

Convert FragPipe PSM data to quantms.io parquet format.

### Description {#fragpipe-description}

Converts FragPipe PSM (Peptide-Spectrum Match) data from `psm.tsv` files to the standardized quantms.io format.

### Parameters {#fragpipe-parameters}

| Parameter         | Type    | Required | Default | Description                          |
| ----------------- | ------- | -------- | ------- | ------------------------------------ |
| `--msms-file`     | Path    | Yes      | -       | FragPipe psm.tsv file path           |
| `--output-folder` | Path    | Yes      | -       | Output directory for generated files |
| `--batch-size`    | Integer | No       | 1000000 | Read batch size                      |
| `--output-prefix` | String  | No       | -       | Prefix for output files              |

### Usage Examples {#fragpipe-examples}

#### Basic Example {#fragpipe-example-basic}

```bash
quantmsioc convert fragpipe \
    --msms-file /path/to/psm.tsv \
    --output-folder ./output
```

#### With Custom Settings {#fragpipe-example-custom}

```bash
quantmsioc convert fragpipe \
    --msms-file /path/to/psm.tsv \
    --output-folder ./output \
    --batch-size 500000 \
    --output-prefix fragpipe_psm
```

### Output Files {#fragpipe-output}

- **Output**: `{output-prefix}-{uuid}.psm.parquet`
- **Format**: Parquet file containing PSM data
- **Schema**: Conforms to quantms.io PSM specification

---

## quantms-psm

Convert mzTab PSM data to quantms.io parquet format.

### Description {#quantms-psm-description}

Converts PSM data from mzTab format to the quantms.io standardized parquet format. Can work with existing DuckDB indexes or create new ones from mzTab files.

### Parameters {#quantms-psm-parameters}

| Parameter         | Type   | Required    | Default | Description                                              |
| ----------------- | ------ | ----------- | ------- | -------------------------------------------------------- |
| `--mztab-path`    | Path   | Conditional | -       | Input mzTab file path (required if creating new indexer) |
| `--database-path` | Path   | Conditional | -       | DuckDB database file path (existing or to be created)    |
| `--output-folder` | Path   | Yes         | -       | Output directory for generated files                     |
| `--output-prefix` | String | No          | `psm`   | Prefix for output files                                  |
| `--spectral-data` | Flag   | No          | False   | Include spectral data fields                             |
| `--verbose`       | Flag   | No          | False   | Enable verbose logging                                   |

### Usage Examples {#quantms-psm-examples}

#### Create from mzTab {#quantms-psm-example-create}

```bash
quantmsioc convert quantms-psm \
    --mztab-path tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz \
    --output-folder ./output \
    --verbose
```

#### Use Existing Database {#quantms-psm-example-existing}

```bash
quantmsioc convert quantms-psm \
    --database-path ./existing_database.duckdb \
    --output-folder ./output \
    --spectral-data
```

### Output Files {#quantms-psm-output}

- **Output**: `{output-prefix}-{uuid}.psm.parquet`
- **Format**: Parquet file containing PSM data
- **Schema**: Conforms to quantms.io PSM specification

### Best Practices {#quantms-psm-best-practices}

- Reuse database files when processing multiple outputs from the same mzTab
- Use `--spectral-data` flag when spectral information is needed for downstream analysis

---

## quantms-feature

Convert mzTab feature data to quantms.io parquet format.

### Description {#quantms-feature-description}

Converts feature-level quantification data from mzTab format to the quantms.io standardized format, including MSstats quantification data.

### Parameters {#quantms-feature-parameters}

| Parameter         | Type   | Required    | Default   | Description                          |
| ----------------- | ------ | ----------- | --------- | ------------------------------------ |
| `--mztab-path`    | Path   | Conditional | -         | Input mzTab file path                |
| `--database-path` | Path   | Conditional | -         | DuckDB database file path            |
| `--output-folder` | Path   | Yes         | -         | Output directory for generated files |
| `--output-prefix` | String | No          | `feature` | Prefix for output files              |
| `--sdrf-file`     | Path   | Yes         | -         | SDRF file path for metadata          |
| `--msstats-file`  | Path   | Yes         | -         | MSstats input file path              |
| `--verbose`       | Flag   | No          | False     | Enable verbose logging               |

### Usage Examples {#quantms-feature-examples}

#### Basic Example {#quantms-feature-example-basic}

```bash
quantmsioc convert quantms-feature \
    --mztab-path tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz \
    --sdrf-file tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv \
    --msstats-file tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz \
    --output-folder ./output
```

### Output Files {#quantms-feature-output}

- **Output**: `{output-prefix}-{uuid}.feature.parquet`
- **Format**: Parquet file containing feature quantification
- **Schema**: Conforms to quantms.io feature specification

---

## quantms-pg

Convert mzTab protein group data to quantms.io parquet format.

### Description {#quantms-pg-description}

Converts protein group quantification from mzTab format to the quantms.io standardized format. Supports both TMT and LFQ data, with optional TopN and iBAQ intensity calculations.

### Parameters {#quantms-pg-parameters}

| Parameter         | Type    | Required    | Default | Description                           |
| ----------------- | ------- | ----------- | ------- | ------------------------------------- |
| `--mztab-path`    | Path    | Conditional | -       | Input mzTab file path                 |
| `--database-path` | Path    | Conditional | -       | DuckDB database file path             |
| `--msstats-file`  | Path    | Yes         | -       | MSstats input file for quantification |
| `--sdrf-file`     | Path    | Yes         | -       | SDRF file path for metadata           |
| `--output-folder` | Path    | Yes         | -       | Output directory for generated files  |
| `--output-prefix` | String  | No          | `pg`    | Prefix for output files               |
| `--compute-topn`  | Flag    | No          | True    | Compute TopN intensity                |
| `--compute-ibaq`  | Flag    | No          | True    | Compute iBAQ intensity                |
| `--topn`          | Integer | No          | 3       | Number of peptides for TopN intensity |
| `--verbose`       | Flag    | No          | False   | Enable verbose logging                |

### Usage Examples {#quantms-pg-examples}

#### LFQ Data with All Intensities {#quantms-pg-example-lfq}

```bash
quantmsioc convert quantms-pg \
    --mztab-path tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz \
    --msstats-file tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz \
    --sdrf-file tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv \
    --output-folder ./output \
    --compute-topn \
    --compute-ibaq \
    --topn 3
```

#### TMT Data (Skip iBAQ) {#quantms-pg-example-tmt}

```bash
quantmsioc convert quantms-pg \
    --mztab-path tests/examples/quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_openms.mzTab.gz \
    --msstats-file tests/examples/quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz \
    --sdrf-file tests/examples/quantms/dda-plex-full/PXD007683-TMT.sdrf.tsv \
    --output-folder ./output \
    --no-compute-ibaq
```

### Output Files {#quantms-pg-output}

- **Output**: `{output-prefix}-{uuid}.pg.parquet`
- **Format**: Parquet file containing protein group quantification
- **Schema**: Conforms to quantms.io protein group specification

### Best Practices {#quantms-pg-best-practices}

- Use `--no-compute-ibaq` for TMT/iTRAQ labeled data
- Adjust `--topn` value based on dataset characteristics (typically 3-5)
- Enable verbose mode for large datasets to monitor progress

---

## idxml

Convert a single OpenMS idXML file to quantms.io PSM format.

### Description {#idxml-description}

Converts PSM data from OpenMS idXML format to the quantms.io standardized parquet format. Can optionally attach spectral information from corresponding mzML files.

### Parameters {#idxml-parameters}

| Parameter              | Type   | Required | Default | Description                                  |
| ---------------------- | ------ | -------- | ------- | -------------------------------------------- |
| `--idxml-file`         | Path   | Yes      | -       | idXML file containing identifications        |
| `--output-folder`      | Path   | Yes      | -       | Output directory for generated files         |
| `--mzml-file`          | Path   | No       | -       | Optional mzML file to attach spectra by scan |
| `--output-prefix-file` | String | No       | `psm`   | Prefix for output files                      |
| `--spectral-data`      | Flag   | No       | False   | Include spectral data fields                 |

### Usage Examples {#idxml-examples}

#### Basic Example {#idxml-example-basic}

```bash
quantmsioc convert idxml \
    --idxml-file tests/examples/idxml/SF_200217_pPeptideLibrary_pool1_HCDnlETcaD_OT_rep2_consensus_fdr_pep_luciphor.idXML \
    --output-folder ./output
```

#### With Spectral Data {#idxml-example-spectral}

```bash
quantmsioc convert idxml \
    --idxml-file tests/examples/idxml/SF_200217_pPeptideLibrary_pool1_HCDnlETcaD_OT_rep2_consensus_fdr_pep_luciphor.idXML \
    --mzml-file tests/examples/idxml/SF_200217_pPeptideLibrary_pool1_HCDnlETcaD_OT_rep1.mzML \
    --output-folder ./output \
    --spectral-data \
    --output-prefix-file idxml_psm_with_spectra
```

### Output Files {#idxml-output}

- **Output**: `{output-prefix-file}-{uuid}.psm.parquet`
- **Format**: Parquet file containing PSM data
- **Schema**: Conforms to quantms.io PSM specification

---

## idxml-batch

Convert multiple OpenMS idXML files to a single merged PSM parquet file.

### Description {#idxml-batch-description}

Batch converts multiple idXML files and merges them into a single quantms.io PSM parquet file. Supports both folder-based and file-list-based input, with flexible mzML matching strategies.

### Parameters {#idxml-batch-parameters}

| Parameter              | Type   | Required    | Default      | Description                                                           |
| ---------------------- | ------ | ----------- | ------------ | --------------------------------------------------------------------- |
| `--idxml-folder`       | Path   | Conditional | -            | Folder containing idXML files (mutually exclusive with --idxml-files) |
| `--idxml-files`        | String | Conditional | -            | Comma-separated list of idXML file paths                              |
| `--output-folder`      | Path   | Yes         | -            | Output directory for merged parquet file                              |
| `--output-prefix-file` | String | No          | `merged-psm` | Prefix for output files                                               |
| `--mzml-folder`        | Path   | No          | -            | Folder containing mzML files (basename matching)                      |
| `--mzml-files`         | String | No          | -            | Comma-separated list of mzML file paths                               |
| `--verbose`            | Flag   | No          | False        | Enable verbose logging                                                |

### Usage Examples {#idxml-batch-examples}

#### Folder-Based Conversion {#idxml-batch-example-folder}

```bash
quantmsioc convert idxml-batch \
    --idxml-folder ./idxml_files \
    --output-folder ./output \
    --output-prefix-file batch_psm
```

#### File List with Index Matching {#idxml-batch-example-list}

```bash
quantmsioc convert idxml-batch \
    --idxml-files file1.idXML,file2.idXML,file3.idXML \
    --mzml-files file1.mzML,file2.mzML,file3.mzML \
    --output-folder ./output \
    --verbose
```

#### Folder with Basename Matching {#idxml-batch-example-basename}

```bash
quantmsioc convert idxml-batch \
    --idxml-folder ./idxml_files \
    --mzml-folder ./mzml_files \
    --output-folder ./output \
    --output-prefix-file merged_with_spectra \
    --verbose
```

### Matching Strategies {#idxml-batch-matching}

The command supports three mzML matching strategies:

1. **Folder-Folder**: Matches files by basename (filename without extension)
2. **List-List**: Matches files by position in the list (index-based)
3. **Folder-List**: Matches folder files by basename with list files

### Output Files {#idxml-batch-output}

- **Output**: `{output-prefix-file}-{uuid}.psm.parquet`
- **Format**: Single merged parquet file containing PSM data from all inputs
- **Schema**: Conforms to quantms.io PSM specification

### Best Practices {#idxml-batch-best-practices}

- Use verbose mode to monitor matching and conversion progress
- Ensure consistent naming when using basename matching
- Verify file order when using index-based matching
- Check temporary directory has sufficient space for large batches

---

## Related Commands

- [Transform Commands](cli-transform.md) - Further process converted data
- [Visualization Commands](cli-visualize.md) - Create plots from converted data
- [Statistics Commands](cli-stats.md) - Generate statistics from converted data
