# Troubleshooting

Common issues and solutions when using qpx.

## Installation Issues

### ModuleNotFoundError: No module named 'qpx'

**Problem**: Python cannot find the qpx package after installation.

**Solutions**:

1. Ensure you're using the correct Python environment:

   ```bash
   which python
   which qpxc
   ```

2. Reinstall in the active environment:

   ```bash
   pip install --force-reinstall qpx
   ```

3. If using conda, ensure the environment is activated:

   ```bash
   conda activate qpx
   ```

### Python version incompatibility

**Problem**: Installation fails with Python version errors.

**Solution**: qpx requires Python 3.10 or higher. Check your version:

```bash
python --version
```

If needed, install a compatible Python version:

```bash
# Using conda
conda create -n qpx python=3.10
conda activate qpx
pip install qpx
```

### Missing dependencies

**Problem**: Import errors for packages like `venn` or `pyopenms`.

**Solution**: Install the missing package (`anndata` / `mudata` / `scipy` ship with core `qpx`):

```bash
pip install venn pyopenms
```

## Conversion Issues

### File not found errors

**Problem**: `FileNotFoundError` when running convert commands.

**Solutions**:

1. Use absolute paths:

   ```bash
   qpxc convert maxquant \
       --msms-file /full/path/to/msms.txt \
       --output-folder /full/path/to/output
   ```

2. Verify the file exists:

   ```bash
   ls -la path/to/file.txt
   ```

### Memory errors with large files

**Problem**: `MemoryError` or system becomes unresponsive with large datasets.

**Solutions**:

1. Process files in batches if supported by the command

2. Increase available memory or use a machine with more RAM

3. For DIA-NN reports, use the `--qvalue-threshold` to filter data:

   ```bash
   qpxc convert diann \
       --report-path report.tsv \
       --sdrf-path metadata.sdrf.tsv \
       --mzml-info-folder ./mzml_info \
       --qvalue-threshold 0.01 \
       --output-folder ./output
   ```

4. For MaxQuant, use the `--batch-size` and `--max-memory` options:

   ```bash
   qpxc convert maxquant \
       --msms-file msms.txt \
       --output-folder ./output \
       --batch-size 50000 \
       --max-memory 8GB
   ```

### Invalid file format errors

**Problem**: `ValueError` or parsing errors when reading input files.

**Solutions**:

1. Verify the file format matches the expected format for the converter

2. Check for file corruption:

   ```bash
   head -20 input_file.txt
   ```

3. Ensure the file encoding is UTF-8:

   ```bash
   file input_file.txt
   ```

## Output Issues

### Empty output files

**Problem**: Parquet files are created but contain no data.

**Solutions**:

1. Check if input data passes quality filters (q-value, PEP thresholds)

2. Verify column names match expected format for the software

3. Use `--verbose` flag to see processing details:

   ```bash
   qpxc convert maxquant \
       --msms-file msms.txt \
       --output-folder ./ \
       --verbose
   ```

### Missing columns in output

**Problem**: Expected columns are not present in the output Parquet file.

**Solutions**:

1. Check if the input file contains the required source columns

2. For spectral data, ensure `--spectral-data` flag is used:

   ```bash
   qpxc convert maxquant \
       --msms-file msms.txt \
       --output-folder ./ \
       --spectral-data
   ```

3. Review the [Format Specification](spec/index.md) for required vs optional fields

## Validation Issues

### Validating a dataset

Use the `validate` command to check a dataset against the canonical schemas:

```bash
# Validate all structures in a dataset
qpxc validate --dataset-path ./PXD014414

# Validate a specific structure
qpxc validate --dataset-path ./PXD014414 --structure feature

# Validate a single Parquet file
qpxc validate --file ./data.feature.parquet
```

### Common validation errors

**Missing required column**: A required column is absent from the Parquet file. Check the [schema reference](spec/schemas.md) for required fields.

**Type mismatch**: A column has a different Arrow type than the schema expects. This usually means the data was written with an older version of qpx or a different tool.

**Null values in non-nullable columns**: Required columns should not contain null values. Check your input data and conversion pipeline.

**Duplicate primary key**: Rows with identical primary key values exist. This may indicate duplicate entries in the source data.

### Programmatic validation

You can also validate from Python:

```python
import qpx

with qpx.open_dataset("./PXD014414") as ds:
    results = ds.validate()
    for name, result in results.items():
        print(result.summary)
        for issue in result.issues:
            print(f"  [{issue.severity}] {issue.message}")
```

## SDRF Issues

### Sample name mismatches

**Problem**: Samples in data files don't match SDRF sample names.

**Solutions**:

1. Ensure `source name` column in SDRF matches file names (without extension)

2. Check for whitespace or case sensitivity issues:

   ```python
   import pandas as pd

   sdrf = pd.read_csv("experiment.sdrf.tsv", sep="\t")
   print(sdrf["source name"].unique())
   ```

### Missing factor values

**Problem**: Factor values are not extracted from SDRF.

**Solution**: Ensure factor columns follow the format `factor value[factor_name]`:

```text
source name    factor value[disease]    factor value[organism part]
sample1        healthy                  liver
sample2        cancer                   liver
```

### Converting SDRF to QPX metadata

Use the dedicated SDRF converter to produce `sample.parquet` and `run.parquet`:

```bash
qpxc convert sdrf \
    --sdrf-file metadata.sdrf.tsv \
    --output-folder ./output
```

## Query Issues

### SQL query errors

**Problem**: SQL queries fail with column or table not found errors.

**Solutions**:

1. Check available structures in the dataset:

   ```bash
   qpxc info --dataset-path ./PXD014414
   ```

2. Check the schema for a specific structure:

   ```bash
   qpxc info schema --dataset-path ./PXD014414 --structure feature
   ```

3. Table names in SQL match the QPX structure names: `psm`, `feature`, `pg`, `mz`, `sample`, `run`, `dataset`, `ontology`, `provenance`.

### Memory issues with large queries

**Solutions**:

1. Use `--duckdb-memory` to increase DuckDB memory (query commands):

   ```bash
   qpxc query sql \
       --dataset-path ./PXD014414 \
       --sql "SELECT * FROM feature" \
       --duckdb-memory 32GB
   ```

   For convert commands, use `--max-memory` instead:

   ```bash
   qpxc convert diann \
       --report-path report.tsv \
       --sdrf-file data.sdrf.tsv \
       --output-folder ./output \
       --max-memory 32GB
   ```

2. Use `--limit` or SQL `LIMIT` to restrict results

3. Export large results directly to Parquet:

   ```bash
   qpxc query sql \
       --dataset-path ./PXD014414 \
       --sql "SELECT * FROM feature" \
       --output results.parquet \
       --output-format parquet
   ```

## Performance Issues

### Slow processing

**Solutions**:

1. Use SSD storage for input/output files

2. Increase available RAM

3. For large datasets, consider processing samples in parallel

4. Use compressed input files (.gz) to reduce I/O

### High memory usage

**Solutions**:

1. Close other applications to free memory

2. Process smaller batches of data

3. Use streaming/chunked processing where available

## Getting More Help

If your issue isn't listed here:

1. **Search existing issues**: [GitHub Issues](https://github.com/bigbio/qpx/issues)

2. **Enable verbose logging**: Add `--verbose` to any command for detailed output

3. **Create a new issue**: Include:
   - qpx version (`qpxc --version`)
   - Python version (`python --version`)
   - Operating system
   - Complete error message
   - Minimal reproducible example
