# Transform Commands

Transform and process data within the quantms.io ecosystem.

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

The `transform` command group provides tools for processing and transforming quantms.io data into various downstream formats. These commands enable absolute and differential expression analysis, metadata mapping, and data format conversions.

## Available Commands

- [ae](#ae) - Convert iBAQ to absolute expression format
- [differential](#differential) - Convert MSstats differential expression data
- [gene](#gene) - Map gene information to proteins
- [ibaq](#ibaq) - Process iBAQ quantification files
- [spectra](#spectra) - Map spectrum information
- [uniprot](#uniprot) - Map latest UniProt annotations
- [anndata](#anndata) - Merge AE files into AnnData format

---

## ae

Convert iBAQ absolute expression data to quantms.io format.

### Description {#ae-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.ae import convert_ibaq_absolute_cmd
print(generate_description(convert_ibaq_absolute_cmd))
```

**Format Specification**: For details about the AE format structure and fields, see the [Absolute Expression Format Specification](https://io.quantms.org/format-specification/#absolute).

### Parameters {#ae-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.ae import convert_ibaq_absolute_cmd
print(generate_params_table(convert_ibaq_absolute_cmd))
```

### Usage Examples {#ae-examples}

#### Basic Example {#ae-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.ae import convert_ibaq_absolute_cmd
print(generate_example(convert_ibaq_absolute_cmd, 'Convert iBAQ data with default settings:'))
```

#### With Project Metadata {#ae-example-metadata}

```bash
quantmsioc transform ae \
    --ibaq-file tests/examples/AE/PXD016999.1-ibaq.tsv \
    --sdrf-file tests/examples/AE/PXD016999-first-instrument.sdrf.tsv \
    --project-file tests/examples/AE/project.json \
    --output-folder ./output \
    --output-prefix ae_with_metadata \
    --delete-existing
```

#### Filter Specific Proteins {#ae-example-filter}

```bash
quantmsioc transform ae \
    --ibaq-file tests/examples/AE/PXD016999.1-ibaq.tsv \
    --sdrf-file tests/examples/AE/PXD016999-first-instrument.sdrf.tsv \
    --protein-file tests/examples/fasta/Homo-sapiens.fasta \
    --output-folder ./output
```

### Input File Formats {#ae-input-formats}

**iBAQ File**: Tab-separated file with protein accessions and iBAQ intensities

```
ProteinID    Sample1    Sample2    Sample3
P12345       1000000    950000     1050000
Q67890       500000     480000     520000
```

**SDRF File**: Standard PRIDE SDRF format with sample metadata

### Output Files {#ae-output}

- **Output**: `{output-prefix}-{uuid}.absolute.parquet`
- **Format**: Parquet file containing absolute expression quantification
- **Schema**: Conforms to [quantms.io absolute expression specification](https://io.quantms.org/format-specification/#absolute)

### Common Issues {#ae-issues}

**Issue**: Mismatched sample names between iBAQ and SDRF

- **Solution**: Ensure column names in iBAQ file match sample identifiers in SDRF

**Issue**: Missing protein accessions

- **Solution**: Provide `--protein-file` to filter and validate protein IDs

### Best Practices {#ae-best-practices}

- Always provide project metadata file when available for better data provenance
- Use `--delete-existing` flag carefully to avoid accidental data loss
- Validate SDRF file format before processing
- Check sample name consistency across input files

---

## differential

Convert MSstats differential expression data to quantms.io format.

### Description {#differential-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.de import convert_msstats_differential_cmd
print(generate_description(convert_msstats_differential_cmd))
```

### Parameters {#differential-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.de import convert_msstats_differential_cmd
print(generate_params_table(convert_msstats_differential_cmd))
```

### Usage Examples {#differential-examples}

#### Basic Example {#differential-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.de import convert_msstats_differential_cmd
print(generate_example(convert_msstats_differential_cmd, 'Convert MSstats differential expression data:'))
```

#### With Custom FDR Threshold {#differential-example-fdr}

```bash
quantmsioc transform differential \
    --msstats-file tests/examples/DE/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv \
    --sdrf-file tests/examples/DE/PXD033169.sdrf.tsv \
    --fdr-threshold 0.01 \
    --output-folder ./output \
    --output-prefix de_stringent \
    --verbose
```

#### With Project Metadata {#differential-example-metadata}

```bash
quantmsioc transform differential \
    --msstats-file tests/examples/DE/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv \
    --sdrf-file tests/examples/DE/PXD033169.sdrf.tsv \
    --project-file tests/examples/DE/project.json \
    --fdr-threshold 0.05 \
    --output-folder ./output \
    --delete-existing
```

### Input File Format {#differential-input-format}

**MSstats File**: CSV file with comparison results

```
Protein,Label,log2FC,SE,Tvalue,DF,pvalue,adj.pvalue
P12345,Condition2-Condition1,2.5,0.3,8.33,10,0.0001,0.001
Q67890,Condition2-Condition1,-1.8,0.4,-4.5,10,0.002,0.01
```

### Output Files {#differential-output}

- **Output**: `{output-prefix}-{uuid}.differential.parquet`
- **Format**: Parquet file containing differential expression results
- **Schema**: Conforms to quantms.io differential expression specification

### Common Issues {#differential-issues}

**Issue**: No significant results after FDR filtering

- **Solution**: Increase `--fdr-threshold` or check input data quality

**Issue**: Memory errors with large comparison files

- **Solution**: Process comparisons in batches or increase available memory

### Best Practices {#differential-best-practices}

- Use FDR threshold of 0.05 or lower for publication-quality results
- Enable verbose mode to monitor filtering statistics
- Validate comparison group names match SDRF metadata
- Include project file for complete data provenance

---

## gene

Map gene information from FASTA to parquet format.

### Description {#gene-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.gene import map_gene_message_cmd
print(generate_description(map_gene_message_cmd))
```

### Parameters {#gene-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.gene import map_gene_message_cmd
print(generate_params_table(map_gene_message_cmd))
```

### Usage Examples {#gene-examples}

#### Basic Example {#gene-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.gene import map_gene_message_cmd
print(generate_example(map_gene_message_cmd, 'Map gene information to parquet file:'))
```

#### With Partitioning {#gene-example-partition}

```bash
quantmsioc transform gene \
    --parquet-path ./output/feature.parquet \
    --fasta tests/examples/fasta/Homo-sapiens.fasta \
    --output-folder ./output \
    --file-num 20 \
    --partitions reference_file_name \
    --species human
```

### Output Files {#gene-output}

- **Output**: Enhanced parquet file(s) with gene information
- **Format**: Parquet file in output folder
- **Added Fields**: Gene names and metadata from FASTA headers

### Best Practices {#gene-best-practices}

- Use species-specific FASTA files for accurate gene annotation
- Adjust `--file-num` based on available memory for large files
- Use partitioning for better file organization in large datasets

---

## ibaq

Convert feature data to iBAQ format.

### Description {#ibaq-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.ibaq import convert_ibaq_file_cmd
print(generate_description(convert_ibaq_file_cmd))
```

### Parameters {#ibaq-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.ibaq import convert_ibaq_file_cmd
print(generate_params_table(convert_ibaq_file_cmd))
```

### Usage Examples {#ibaq-examples}

#### Basic Example {#ibaq-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.ibaq import convert_ibaq_file_cmd
print(generate_example(convert_ibaq_file_cmd, 'Convert feature data to iBAQ format:'))
```

#### With Custom Prefix {#ibaq-example-prefix}

```bash
quantmsioc transform ibaq \
    --feature-file ./output/feature.parquet \
    --sdrf-file ./metadata.sdrf.tsv \
    --output-folder ./output \
    --output-prefix ibaq_quantification
```

### Output Files {#ibaq-output}

- **Output**: `{output-prefix}-{uuid}.ibaq.parquet`
- **Format**: Parquet file containing iBAQ quantification values
- **Content**: Protein-level iBAQ values per sample

### Best Practices {#ibaq-best-practices}

- Ensure feature file contains all necessary quantification data
- Verify SDRF metadata matches sample identifiers in feature file
- Use iBAQ output for absolute protein quantification analysis

---

## spectra

Map spectrum information from mzML to parquet format.

### Description {#spectra-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.spectra import map_spectrum_message_cmd
print(generate_description(map_spectrum_message_cmd))
```

### Parameters {#spectra-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.spectra import map_spectrum_message_cmd
print(generate_params_table(map_spectrum_message_cmd))
```

### Usage Examples {#spectra-examples}

#### Basic Example {#spectra-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.spectra import map_spectrum_message_cmd
print(generate_example(map_spectrum_message_cmd, 'Map spectrum information to parquet:'))
```

#### With Batch Processing and Partitioning {#spectra-example-batch}

```bash
quantmsioc transform spectra \
    --parquet-path ./output/psm.parquet \
    --mzml-directory ./mzml_files \
    --output-folder ./output \
    --file-num 20 \
    --partitions reference_file_name
```

### Output Files {#spectra-output}

- **Output**: Enhanced PSM/feature parquet file(s) with spectral information
- **Format**: Parquet file in output folder
- **Added Fields**: Spectrum metadata and peak information from mzML files

### Best Practices {#spectra-best-practices}

- Ensure mzML files are in the specified directory with correct naming
- Adjust `--file-num` based on available memory and file size
- Use partitioning for organized output when processing large datasets

---

## uniprot

Map feature data to latest UniProt version.

### Description {#uniprot-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.uniprot import map_latest_uniprot_cmd
print(generate_description(map_latest_uniprot_cmd))
```

### Parameters {#uniprot-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.uniprot import map_latest_uniprot_cmd
print(generate_params_table(map_latest_uniprot_cmd))
```

### Usage Examples {#uniprot-examples}

#### Basic Mapping {#uniprot-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.uniprot import map_latest_uniprot_cmd
print(generate_example(map_latest_uniprot_cmd, 'Map features to latest UniProt:'))
```

#### With Custom Prefix {#uniprot-example-prefix}

```bash
quantmsioc transform uniprot \
    --feature-file ./output/feature.parquet \
    --fasta ./uniprot_human_2024.fasta \
    --output-folder ./output \
    --output-prefix feature_updated
```

### Output Files {#uniprot-output}

- **Output**: `{output-prefix}-{uuid}.feature.parquet`
- **Format**: Parquet file with updated UniProt mappings
- **Content**: Feature data mapped to latest UniProt protein identifications

### Best Practices {#uniprot-best-practices}

- Use the latest UniProt FASTA file for most current annotations
- Run this command when updating to a new UniProt release
- Verify FASTA file matches the organism of your study

---

## anndata

Merge multiple AE files into a file in AnnData format.

### Description {#anndata-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.anndata import merge_ae_files_cmd
print(generate_description(merge_ae_files_cmd))
```

### Parameters {#anndata-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.anndata import merge_ae_files_cmd
print(generate_params_table(merge_ae_files_cmd))
```

### Usage Examples {#anndata-examples}

#### Basic Example {#anndata-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.transform.anndata import merge_ae_files_cmd
print(generate_example(merge_ae_files_cmd, 'Merge AE files into AnnData format:'))
```

#### With Custom Prefix {#anndata-example-prefix}

```bash
quantmsioc transform anndata \
    --directory ./ae_files \
    --output-folder ./output \
    --output-prefix merged_ae
```

### Output Files {#anndata-output}

- **Output**: `{output-prefix}-{uuid}.h5ad`
- **Format**: AnnData H5AD file (HDF5-based format)
- **Structure**:
  - `X`: Protein expression matrix
  - `obs`: Sample metadata
  - `var`: Protein metadata

### Use Cases {#anndata-use-cases}

- Integration with scanpy for dimensionality reduction and clustering
- Compatibility with Python machine learning libraries
- Multi-experiment data integration
- Cross-study meta-analysis

### Best Practices {#anndata-best-practices}

- Ensure all AE files in the directory have consistent format
- Validate sample metadata consistency before merging
- Use the output with scanpy or other AnnData-compatible tools
- Consider memory requirements for large datasets

---

## Related Commands

- [Convert Commands](cli-convert.md) - Convert raw data to quantms.io format
- [Visualization Commands](cli-visualize.md) - Visualize transformed data
- [Statistics Commands](cli-stats.md) - Analyze transformed data
