# Convert Commands

Convert various mass spectrometry data formats to the quantms.io standard format.

```python exec="1" session="doc_utils" result="ansi"
import click
import textwrap

def get_click_type_display(param):
    param_type = param.type
    type_str = str(param_type)
    if 'Path' in type_str:
        if hasattr(param_type, 'dir_okay') and not param_type.dir_okay:
            return 'FILE'
        elif hasattr(param_type, 'file_okay') and not param_type.file_okay:
            return 'DIRECTORY'
        else:
            return 'PATH'
    elif isinstance(param_type, click.types.FloatParamType):
        return 'FLOAT'
    elif isinstance(param_type, click.types.IntParamType):
        return 'INTEGER'
    elif param.is_flag:
        return 'FLAG'
    else:
        return 'TEXT'

def generate_params_table(command):
    table = '<table>\n<thead>\n<tr>\n'
    table += '<th>Parameter</th><th>Type</th><th>Required</th><th>Default</th><th>Description</th>\n'
    table += '</tr>\n</thead>\n<tbody>\n'
    for param in command.params:
        if isinstance(param, click.Option) and param.name not in ['help']:
            param_names = param.opts
            param_name = param_names[0] if param_names else f"--{param.name}"
            param_type = get_click_type_display(param)
            required = 'Yes' if param.required else 'No'

            # Extract default value from help text if it contains "(default: ...)"
            description = param.help or ''
            default_from_help = None
            if '(default:' in description.lower():
                import re
                match = re.search(r'\(default:\s*([^)]+)\)', description, re.IGNORECASE)
                if match:
                    default_from_help = match.group(1).strip()

            # Determine default value display
            if default_from_help:
                default = default_from_help
            elif param.default is not None and str(param.default) != 'Sentinel.UNSET':
                if param.is_flag:
                    default = '-'
                elif isinstance(param.default, (int, float)):
                    default = str(param.default)
                elif isinstance(param.default, str):
                    default = f'<code>{param.default}</code>'
                else:
                    default = str(param.default)
            else:
                default = '-'

            table += f'<tr>\n<td><code>{param_name}</code></td>\n<td>{param_type}</td>\n<td>{required}</td>\n<td>{default}</td>\n<td>{description}</td>\n</tr>\n'
    table += '</tbody>\n</table>'
    return table

def generate_description(command):
    if command.help:
        help_text = command.help
        if 'Example' in help_text:
            description = help_text.split('Example')[0].strip()
        else:
            description = help_text.strip()
        lines = description.split('\n')
        if len(lines) > 1:
            description = '\n'.join(lines[1:]).strip()
            return f'<p>{description}</p>'
    return ''

def generate_example(command, default_text=''):
    if command.help and 'Example' in command.help:
        example_section = command.help.split('Example')[1]
        if ':' in example_section:
            example_section = example_section.split(':', 1)[1]
        example_section = textwrap.dedent(example_section).strip()
        output = ''
        if default_text:
            output += f'<p>{default_text}</p>\n'
        output += f'<pre><code class="language-bash">{example_section}</code></pre>'
        return output
    return ''
```

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

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.diann import convert_diann_cmd
print(generate_description(convert_diann_cmd))
```

### Parameters {#diann-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.diann import convert_diann_cmd
print(generate_params_table(convert_diann_cmd))
```

### Usage Examples {#diann-examples}

#### Basic Example {#diann-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.diann import convert_diann_cmd
print(generate_example(convert_diann_cmd, 'Convert a DIA-NN report with default settings:'))
```

#### Advanced Example {#diann-example-advanced}

Convert with file partitioning based on reference_file_name:

```bash
quantmsioc convert diann \
    --report-path tests/examples/diann/small/diann_report.tsv \
    --qvalue-threshold 0.05 \
    --mzml-info-folder tests/examples/diann/small/mzml \
    --sdrf-path tests/examples/diann/small/PXD019909-DIA.sdrf.tsv \
    --output-folder ./output
```

#### Advanced Example with Partitioning {#diann-example-partitioning}

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

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.diann import convert_diann_pg_cmd
print(generate_description(convert_diann_pg_cmd))
```

### Parameters {#diann-pg-parameters}

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.diann import convert_diann_pg_cmd
print(generate_params_table(convert_diann_pg_cmd))
```

### Usage Examples {#diann-pg-examples}

#### Basic Example {#diann-pg-example-basic}

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.diann import convert_diann_pg_cmd
print(generate_example(convert_diann_pg_cmd, 'Convert DIA-NN protein groups with default settings:'))
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

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.maxquant import convert_maxquant_psm_cmd
print(generate_description(convert_maxquant_psm_cmd))
```

### Parameters {#maxquant-psm-parameters}

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.maxquant import convert_maxquant_psm_cmd
print(generate_params_table(convert_maxquant_psm_cmd))
```

### Usage Examples {#maxquant-psm-examples}

#### Basic Example {#maxquant-psm-example-basic}

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.maxquant import convert_maxquant_psm_cmd
print(generate_example(convert_maxquant_psm_cmd, 'Convert MaxQuant PSM data with default settings:'))
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

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.maxquant import convert_maxquant_feature_cmd
print(generate_description(convert_maxquant_feature_cmd))
```

### Parameters {#maxquant-feature-parameters}

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.maxquant import convert_maxquant_feature_cmd
print(generate_params_table(convert_maxquant_feature_cmd))
```

### Usage Examples {#maxquant-feature-examples}

#### Basic Example {#maxquant-feature-example-basic}

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.maxquant import convert_maxquant_feature_cmd
print(generate_example(convert_maxquant_feature_cmd, 'Convert MaxQuant feature data with default settings:'))
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

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.maxquant import convert_maxquant_pg_cmd
print(generate_description(convert_maxquant_pg_cmd))
```

### Parameters {#maxquant-pg-parameters}

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.maxquant import convert_maxquant_pg_cmd
print(generate_params_table(convert_maxquant_pg_cmd))
```

### Usage Examples {#maxquant-pg-examples}

#### Basic Example {#maxquant-pg-example-basic}

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.maxquant import convert_maxquant_pg_cmd
print(generate_example(convert_maxquant_pg_cmd, 'Convert MaxQuant protein groups with default settings:'))
```

### Output Files {#maxquant-pg-output}

- **Output**: `{output-prefix}-{uuid}.pg.parquet`
- **Format**: Parquet file containing protein group data
- **Schema**: Conforms to quantms.io protein group specification

---

## fragpipe

Convert FragPipe PSM data to quantms.io parquet format.

### Description {#fragpipe-description}

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.fragpipe import convert_fragpipe_psm_cmd
print(generate_description(convert_fragpipe_psm_cmd))
```

### Parameters {#fragpipe-parameters}

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.fragpipe import convert_fragpipe_psm_cmd
print(generate_params_table(convert_fragpipe_psm_cmd))
```

### Usage Examples {#fragpipe-examples}

#### Basic Example {#fragpipe-example-basic}

```python exec="1" session="doc_utils" html="1"
from quantmsio.commands.convert.fragpipe import convert_fragpipe_psm_cmd
print(generate_example(convert_fragpipe_psm_cmd, 'Convert FragPipe PSM data with default settings:'))
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

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.quantms import convert_quantms_psm_cmd
print(generate_description(convert_quantms_psm_cmd))
```

### Parameters {#quantms-psm-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.quantms import convert_quantms_psm_cmd
print(generate_params_table(convert_quantms_psm_cmd))
```

### Usage Examples {#quantms-psm-examples}

#### Basic Example {#quantms-psm-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.quantms import convert_quantms_psm_cmd
print(generate_example(convert_quantms_psm_cmd, 'Convert PSM data with default settings:'))
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

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.quantms import convert_quantms_feature_cmd
print(generate_description(convert_quantms_feature_cmd))
```

### Parameters {#quantms-feature-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.quantms import convert_quantms_feature_cmd
print(generate_params_table(convert_quantms_feature_cmd))
```

### Usage Examples {#quantms-feature-examples}

#### Basic Example {#quantms-feature-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.quantms import convert_quantms_feature_cmd
print(generate_example(convert_quantms_feature_cmd, 'Convert feature data with default settings:'))
```

### Output Files {#quantms-feature-output}

- **Output**: `{output-prefix}-{uuid}.feature.parquet`
- **Format**: Parquet file containing feature quantification
- **Schema**: Conforms to quantms.io feature specification

---

## quantms-pg

Convert mzTab protein group data to quantms.io parquet format.

### Description {#quantms-pg-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.quantms import convert_quantms_pg_cmd
print(generate_description(convert_quantms_pg_cmd))
```

### Parameters {#quantms-pg-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.quantms import convert_quantms_pg_cmd
print(generate_params_table(convert_quantms_pg_cmd))
```

### Usage Examples {#quantms-pg-examples}

#### Basic Example {#quantms-pg-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.quantms import convert_quantms_pg_cmd
print(generate_example(convert_quantms_pg_cmd, 'Convert protein groups with default settings:'))
```

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

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.idxml import convert_idxml_file
print(generate_description(convert_idxml_file))
```

### Parameters {#idxml-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.idxml import convert_idxml_file
print(generate_params_table(convert_idxml_file))
```

### Usage Examples {#idxml-examples}

#### Basic Example {#idxml-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.idxml import convert_idxml_file
print(generate_example(convert_idxml_file, 'Convert a single idXML file:'))
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

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.idxml import convert_idxml_batch
print(generate_description(convert_idxml_batch))
```

### Parameters {#idxml-batch-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.idxml import convert_idxml_batch
print(generate_params_table(convert_idxml_batch))
```

### Usage Examples {#idxml-batch-examples}

#### Basic Example {#idxml-batch-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.convert.idxml import convert_idxml_batch
print(generate_example(convert_idxml_batch, 'Convert multiple idXML files:'))
```

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
