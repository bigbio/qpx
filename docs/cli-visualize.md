# Visualization Commands

Create various data visualization plots from quantms.io data.

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

The `visualize` command group provides tools for creating publication-quality visualizations from quantms.io parquet files. All plots are saved in vector formats (SVG, PDF) for high-resolution output.

## Available Commands

All visualization commands are accessed through the `plot` subcommand:

- [psm-peptides](#psm-peptides) - Plot peptides by condition in LFQ
- [ibaq-distribution](#ibaq-distribution) - Plot iBAQ distribution
- [kde-intensity](#kde-intensity) - Plot KDE intensity distribution
- [peptide-distribution](#peptide-distribution) - Plot peptide distribution across proteins
- [box-intensity](#box-intensity) - Plot intensity box plots

---

## psm-peptides

Plot peptides by condition in label-free quantification (LFQ) experiments.

### Description {#psm-peptides-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_peptides_cmd
print(generate_description(plot_peptides_cmd))
```

### Parameters {#psm-peptides-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_peptides_cmd
print(generate_params_table(plot_peptides_cmd))
```

### Usage Examples {#psm-peptides-examples}

#### Basic Example {#psm-peptides-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_peptides_cmd
print(generate_example(plot_peptides_cmd, 'Plot peptides by condition:'))
```

#### Generate PDF Output {#psm-peptides-example-pdf}

```bash
quantmsioc visualize plot psm-peptides \
    --psm-parquet-path ./output/psm.parquet \
    --sdrf-path ./metadata.sdrf.tsv \
    --save-path ./plots/peptides_by_condition.pdf
```

### Output {#psm-peptides-output}

- **Format**: SVG or PDF (based on file extension in `--save-path`)
- **Content**: Bar plot showing peptide counts per condition
- **Axes**:
  - X-axis: Experimental conditions
  - Y-axis: Number of identified peptides

### Interpretation {#psm-peptides-interpretation}

- **High variation**: May indicate batch effects or quality issues
- **Low counts**: May suggest technical problems with specific samples
- **Consistent counts**: Indicates good data quality and reproducibility

### Best Practices {#psm-peptides-best-practices}

- Use SVG format for publications (scalable vector graphics)
- Verify SDRF metadata correctly defines experimental conditions
- Compare with expected peptide yields for your sample type

---

## ibaq-distribution

Plot the distribution of iBAQ (intensity-Based Absolute Quantification) values.

### Description {#ibaq-distribution-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_ibaq_distribution_cmd
print(generate_description(plot_ibaq_distribution_cmd))
```

### Parameters {#ibaq-distribution-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_ibaq_distribution_cmd
print(generate_params_table(plot_ibaq_distribution_cmd))
```

### Usage Examples {#ibaq-distribution-examples}

#### Plot All Samples {#ibaq-distribution-example-all}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_ibaq_distribution_cmd
print(generate_example(plot_ibaq_distribution_cmd, 'Plot iBAQ distribution:'))
```

#### Plot Specific Sample {#ibaq-distribution-example-specific}

```bash
quantmsioc visualize plot ibaq-distribution \
    --ibaq-path tests/examples/AE/PXD016999.1-ibaq.tsv \
    --select-column Sample_001 \
    --save-path ./plots/ibaq_sample001.svg
```

### Output {#ibaq-distribution-output}

- **Format**: SVG or PDF
- **Content**: Distribution plot (histogram + kernel density estimate)
- **Axes**:
  - X-axis: log10(iBAQ intensity)
  - Y-axis: Density or frequency

### Interpretation {#ibaq-distribution-interpretation}

- **Bimodal distribution**: May indicate distinct protein abundance classes
- **Long tail**: High-abundance proteins (housekeeping, abundant structural proteins)
- **Narrow range**: Limited dynamic range, possible detection issues

### Best Practices {#ibaq-distribution-best-practices}

- Log-transform intensities for better visualization
- Compare distributions across samples to identify outliers
- Check for batch effects if distributions vary significantly

---

## kde-intensity

Plot Kernel Density Estimation (KDE) of intensity distributions across samples.

### Description {#kde-intensity-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_kde_intensity_distribution_cmd
print(generate_description(plot_kde_intensity_distribution_cmd))
```

### Parameters {#kde-intensity-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_kde_intensity_distribution_cmd
print(generate_params_table(plot_kde_intensity_distribution_cmd))
```

### Usage Examples {#kde-intensity-examples}

#### Plot Default Samples {#kde-intensity-example-default}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_kde_intensity_distribution_cmd
print(generate_example(plot_kde_intensity_distribution_cmd, 'Plot KDE intensity distribution:'))
```

#### Plot More Samples {#kde-intensity-example-more}

```bash
quantmsioc visualize plot kde-intensity \
    --feature-path ./output/feature.parquet \
    --save-path ./plots/intensity_kde_all.svg \
    --num-samples 20
```

### Output {#kde-intensity-output}

- **Format**: SVG or PDF
- **Content**: Overlaid KDE curves for each sample
- **Axes**:
  - X-axis: log10(intensity)
  - Y-axis: Density
- **Legend**: Sample identifiers

### Interpretation {#kde-intensity-interpretation}

- **Overlapping curves**: Good sample-to-sample consistency
- **Shifted curves**: Potential batch effects or normalization issues
- **Different shapes**: Sample-specific technical issues

### Best Practices {#kde-intensity-best-practices}

- Limit to 10-20 samples for readability
- Use this plot to identify samples requiring normalization
- Compare before and after normalization

---

## peptide-distribution

Plot the distribution of peptides across proteins.

### Description {#peptide-distribution-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_peptide_distribution_cmd
print(generate_description(plot_peptide_distribution_cmd))
```

### Parameters {#peptide-distribution-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_peptide_distribution_cmd
print(generate_params_table(plot_peptide_distribution_cmd))
```

### Usage Examples {#peptide-distribution-examples}

#### Basic Example {#peptide-distribution-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_peptide_distribution_cmd
print(generate_example(plot_peptide_distribution_cmd, 'Plot peptide distribution:'))
```

#### Show Top 50 Proteins {#peptide-distribution-example-top50}

```bash
quantmsioc visualize plot peptide-distribution \
    --feature-path ./output/feature.parquet \
    --save-path ./plots/peptide_per_protein_top50.svg \
    --num-samples 50
```

### Output {#peptide-distribution-output}

- **Format**: SVG or PDF
- **Content**: Bar plot showing peptide counts per protein
- **Axes**:
  - X-axis: Protein identifiers (top N by peptide count)
  - Y-axis: Number of identified peptides

### Interpretation {#peptide-distribution-interpretation}

- **High peptide counts**: Abundant proteins with good coverage
- **Single peptide proteins**: May be less confident identifications
- **Distribution shape**: Reflects proteome complexity

### Best Practices {#peptide-distribution-best-practices}

- Focus on top proteins for initial quality assessment
- Filter single-peptide identifications for high-confidence datasets
- Compare expected vs. observed peptide counts for key proteins

---

## box-intensity

Plot box plots of intensity distributions across samples.

### Description {#box-intensity-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_box_intensity_distribution_cmd
print(generate_description(plot_box_intensity_distribution_cmd))
```

### Parameters {#box-intensity-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_box_intensity_distribution_cmd
print(generate_params_table(plot_box_intensity_distribution_cmd))
```

### Usage Examples {#box-intensity-examples}

#### Basic Example {#box-intensity-example-basic}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.plot import plot_box_intensity_distribution_cmd
print(generate_example(plot_box_intensity_distribution_cmd, 'Plot box intensity distribution:'))
```

#### Plot All Samples {#box-intensity-example-all}

```bash
quantmsioc visualize plot box-intensity \
    --feature-path ./output/feature.parquet \
    --save-path ./plots/intensity_boxplot_all.svg \
    --num-samples 50
```

### Output {#box-intensity-output}

- **Format**: SVG or PDF
- **Content**: Box plots for each sample
- **Axes**:
  - X-axis: Sample identifiers
  - Y-axis: log10(intensity)
- **Elements**: Box (IQR), whiskers (1.5×IQR), outliers (points)

### Interpretation {#box-intensity-interpretation}

- **Aligned medians**: Good normalization
- **Similar IQR**: Consistent quantification across samples
- **Many outliers**: May indicate contamination or technical issues
- **Different ranges**: Batch effects or loading differences

### Best Practices {#box-intensity-best-practices}

- Use this plot to identify samples requiring normalization
- Check for systematic differences between batches or conditions
- Compare before and after normalization to verify effectiveness
- Flag samples with unusual distributions for further investigation

---

## General Plotting Tips

### Output Formats

- **SVG**: Recommended for publications (scalable, editable)
- **PDF**: Alternative vector format (portable)
- **PNG**: Raster format (not recommended for publications)

### Color Considerations

- Plots use color-blind friendly palettes by default
- Ensure sufficient contrast for grayscale printing

### Size and Resolution

- Vector formats (SVG/PDF) scale without quality loss
- Suitable for both presentations and manuscripts

### Customization

For advanced customization beyond these commands, consider:

1. Exporting data to CSV and using custom plotting scripts
2. Using the quantmsio Python API for programmatic access
3. Importing SVG files into vector graphics editors

---

## Related Commands

- [Convert Commands](cli-convert.md) - Prepare data for visualization
- [Transform Commands](cli-transform.md) - Process data before plotting
- [Statistics Commands](cli-stats.md) - Generate numeric summaries
