# Statistics Commands

Perform statistical analysis on quantms.io data.

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

The `stats` command group provides tools for generating comprehensive statistical summaries of quantms.io data files. These commands help assess data quality, completeness, and provide key metrics for experimental reports.

## Available Commands

All statistics commands are accessed through the `analyze` subcommand:

- [project-ae](#project-ae) - Generate statistics for absolute expression data
- [psm](#psm) - Generate statistics for PSM data

---

## project-ae

Generate comprehensive statistics for a project's absolute expression data.

### Description {#project-ae-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.stats import project_ae_statistics_cmd
print(generate_description(project_ae_statistics_cmd))
```

### Parameters {#project-ae-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.stats import project_ae_statistics_cmd
print(generate_params_table(project_ae_statistics_cmd))
```

### Usage Examples {#project-ae-examples}

#### Print to Console {#project-ae-example-console}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.stats import project_ae_statistics_cmd
print(generate_example(project_ae_statistics_cmd, 'Generate project statistics:'))
```

#### Save to File {#project-ae-example-file}

```bash
quantmsioc stats analyze project-ae \
    --absolute-path ./output/ae.parquet \
    --parquet-path ./output/psm.parquet \
    --save-path ./reports/project_statistics.txt
```

### Output Format {#project-ae-output}

The command generates a text report with the following metrics:

```
Number of proteins: 2,547
Number of peptides: 12,384
Number of samples: 24
Number of peptidoforms: 15,921
Number of msruns: 24
iBAQ Number of proteins: 2,547
iBAQ Number of samples: 24
```

### Metrics Explained {#project-ae-metrics}

| Metric                      | Description                         |
| --------------------------- | ----------------------------------- |
| **Number of proteins**      | Unique protein identifications      |
| **Number of peptides**      | Unique peptide sequences identified |
| **Number of samples**       | Biological samples in the dataset   |
| **Number of peptidoforms**  | Unique modified peptide forms       |
| **Number of msruns**        | Total MS runs performed             |
| **iBAQ Number of proteins** | Proteins with iBAQ quantification   |
| **iBAQ Number of samples**  | Samples with iBAQ data              |

### Use Cases {#project-ae-use-cases}

- **Quality Control**: Verify expected number of identifications
- **Publication Reporting**: Generate summary statistics for methods sections
- **Data Completeness**: Assess coverage across samples
- **Comparative Analysis**: Compare statistics across different processing pipelines

### Best Practices {#project-ae-best-practices}

- Run statistics after data processing to verify completeness
- Compare statistics with expected values based on sample type and instrument
- Use for QC before downstream analysis
- Include in supplementary materials for publications

---

## psm

Generate statistics for PSM (Peptide-Spectrum Match) data.

### Description {#psm-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.stats import psm_statistics_cmd
print(generate_description(psm_statistics_cmd))
```

### Parameters {#psm-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.stats import psm_statistics_cmd
print(generate_params_table(psm_statistics_cmd))
```

### Usage Examples {#psm-examples}

#### Print to Console {#psm-example-console}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.stats import psm_statistics_cmd
print(generate_example(psm_statistics_cmd, 'Generate PSM statistics:'))
```

#### Save to File {#psm-example-file}

```bash
quantmsioc stats analyze psm \
    --parquet-path ./output/psm.parquet \
    --save-path ./reports/psm_statistics.txt
```

### Output Format {#psm-output}

```
Number of proteins: 1,823
Number of peptides: 8,642
Number of peptidoforms: 11,205
Number of PSMs: 45,892
Number of msruns: 12
```

### Metrics Explained {#psm-metrics}

| Metric                     | Description                                      |
| -------------------------- | ------------------------------------------------ |
| **Number of proteins**     | Unique proteins with at least one PSM            |
| **Number of peptides**     | Unique peptide sequences (without modifications) |
| **Number of peptidoforms** | Unique peptide sequences with modifications      |
| **Number of PSMs**         | Total peptide-spectrum matches                   |
| **Number of msruns**       | Number of MS runs contributing data              |

### Understanding the Metrics {#psm-understanding}

#### Peptide vs Peptidoform

- **Peptide**: Amino acid sequence (e.g., `PEPTIDE`)
- **Peptidoform**: Peptide + modifications (e.g., `PEPTIDE[+16]`)

#### Expected Ratios

Typical ratios for quality data:

- **PSMs per peptide**: 2-5 (varies with replicates)
- **Peptidoforms per peptide**: 1-3 (depends on PTM analysis)
- **Peptides per protein**: 3-20 (depends on protein abundance and coverage)

### Use Cases {#psm-use-cases}

- **Quality Assessment**: Verify identification rates
- **Method Optimization**: Compare different search parameters
- **Replication Analysis**: Assess consistency across runs
- **FDR Validation**: Ensure sufficient identifications after filtering

### Best Practices {#psm-best-practices}

- Compare PSM counts before and after FDR filtering
- Monitor peptidoform counts to assess modification analysis quality
- Track msrun numbers to verify all files processed correctly
- Use as input for sample size calculations in future experiments

---

## Statistical Interpretation Guide

### Data Quality Indicators

#### Good Quality Signs

- Consistent protein counts across replicates
- Expected PSM-to-peptide ratios
- Complete data across all msruns
- Reasonable peptidoform diversity

#### Potential Issues

- **Very low PSM counts**: Search parameter issues, poor sample quality
- **High peptidoform-to-peptide ratio**: Over-prediction of modifications
- **Missing msruns**: File processing errors
- **Extreme variation**: Batch effects, contamination

### Comparative Analysis

When comparing datasets:

1. **Normalize by sample amount**: Account for loading differences
2. **Consider instrument type**: Different platforms yield different numbers
3. **Account for search space**: Database size affects identification rates
4. **Match FDR thresholds**: Ensure fair comparison

### Reporting Guidelines

For publications, report:

- Total unique proteins, peptides, and PSMs
- Number of biological replicates
- Number of technical replicates (msruns)
- FDR thresholds applied
- Protein and peptide identification rates

---

## Automation and Integration

### Batch Processing

Process multiple files:

```bash
#!/bin/bash
for file in ./output/*.psm.parquet; do
    filename=$(basename "$file" .parquet)
    quantmsioc stats analyze psm \
        --parquet-path "$file" \
        --save-path "./reports/${filename}_stats.txt"
done
```

### Integration with Reports

Use output in automated reports:

```bash
quantmsioc stats analyze psm \
    --parquet-path ./output/psm.parquet \
    --save-path ./reports/stats.txt

# Extract specific metrics
grep "Number of proteins" ./reports/stats.txt >> ./summary_report.md
```

### Combining with Visualization

Generate both statistics and plots:

```bash
# Generate statistics
quantmsioc stats analyze psm \
    --parquet-path ./output/psm.parquet \
    --save-path ./reports/stats.txt

# Generate visualizations
quantmsioc visualize plot peptide-distribution \
    --feature-path ./output/feature.parquet \
    --save-path ./plots/peptide_dist.svg
```

---

## Advanced Usage

### Quality Control Workflow

Complete QC workflow example:

```bash
#!/bin/bash

# 1. Generate PSM statistics
quantmsioc stats analyze psm \
    --parquet-path ./output/psm.parquet \
    --save-path ./qc/psm_stats.txt

# 2. Generate AE statistics (if available)
if [ -f ./output/ae.parquet ]; then
    quantmsioc stats analyze project-ae \
        --absolute-path ./output/ae.parquet \
        --parquet-path ./output/psm.parquet \
        --save-path ./qc/ae_stats.txt
fi

# 3. Create visualizations
quantmsioc visualize plot box-intensity \
    --feature-path ./output/feature.parquet \
    --save-path ./qc/intensity_boxplot.svg

quantmsioc visualize plot peptide-distribution \
    --feature-path ./output/feature.parquet \
    --save-path ./qc/peptide_distribution.svg

echo "QC report generated in ./qc/"
```

### Custom Thresholds

Define quality thresholds:

```bash
#!/bin/bash

# Generate statistics
quantmsioc stats analyze psm \
    --parquet-path ./output/psm.parquet \
    --save-path ./stats.txt

# Check thresholds
proteins=$(grep "Number of proteins" ./stats.txt | awk '{print $4}' | tr -d ',')
peptides=$(grep "Number of peptides" ./stats.txt | awk '{print $4}' | tr -d ',')

if [ "$proteins" -lt 1000 ]; then
    echo "WARNING: Low protein count ($proteins)"
fi

if [ "$peptides" -lt 5000 ]; then
    echo "WARNING: Low peptide count ($peptides)"
fi
```

---

## Related Commands

- [Convert Commands](cli-convert.md) - Generate data files for analysis
- [Transform Commands](cli-transform.md) - Process data before statistics
- [Visualization Commands](cli-visualize.md) - Create visual representations of statistics
