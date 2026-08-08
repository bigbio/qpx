# Transform Commands

Transform and process data within the QPX ecosystem.

```python exec="1" session="doc_utils" result="ansi"
import click
import textwrap


def get_click_type_display(param):
    param_type = param.type
    type_str = str(param_type)
    if "Path" in type_str:
        if hasattr(param_type, "dir_okay") and not param_type.dir_okay:
            return "FILE"
        elif hasattr(param_type, "file_okay") and not param_type.file_okay:
            return "DIRECTORY"
        else:
            return "PATH"
    elif isinstance(param_type, click.types.FloatParamType):
        return "FLOAT"
    elif isinstance(param_type, click.types.IntParamType):
        return "INTEGER"
    elif param.is_flag:
        return "FLAG"
    else:
        return "TEXT"


def generate_params_table(command):
    table = "<table>\n<thead>\n<tr>\n"
    table += "<th>Parameter</th><th>Type</th><th>Required</th><th>Default</th><th>Description</th>\n"
    table += "</tr>\n</thead>\n<tbody>\n"
    for param in command.params:
        if isinstance(param, click.Option) and param.name not in ["help"]:
            param_names = param.opts
            param_name = param_names[0] if param_names else f"--{param.name}"
            param_type = get_click_type_display(param)
            required = "Yes" if param.required else "No"

            # Extract default value from help text if it contains "(default: ...)"
            description = param.help or ""
            default_from_help = None
            if "(default:" in description.lower():
                import re

                match = re.search(r"\(default:\s*([^)]+)\)", description, re.IGNORECASE)
                if match:
                    default_from_help = match.group(1).strip()

            # Determine default value display
            if default_from_help:
                default = default_from_help
            elif param.default is not None and str(param.default) != "Sentinel.UNSET":
                if param.is_flag:
                    default = "-"
                elif isinstance(param.default, (int, float)):
                    default = str(param.default)
                elif isinstance(param.default, str):
                    default = f"<code>{param.default}</code>"
                else:
                    default = str(param.default)
            else:
                default = "-"

            table += f"<tr>\n<td><code>{param_name}</code></td>\n<td>{param_type}</td>\n<td>{required}</td>\n<td>{default}</td>\n<td>{description}</td>\n</tr>\n"
    table += "</tbody>\n</table>"
    return table


def generate_description(command):
    if command.help:
        help_text = command.help
        if "Example" in help_text:
            description = help_text.split("Example")[0].strip()
        else:
            description = help_text.strip()
        lines = description.split("\n")
        if len(lines) > 1:
            description = "\n".join(lines[1:]).strip()
            return f"<p>{description}</p>"
    return ""


def generate_example(command, default_text=""):
    if command.help and "Example" in command.help:
        example_section = command.help.split("Example")[1]
        if ":" in example_section:
            example_section = example_section.split(":", 1)[1]
        example_section = textwrap.dedent(example_section).strip()
        output = ""
        if default_text:
            output += f"<p>{default_text}</p>\n"
        output += f'<pre><code class="language-bash">{example_section}</code></pre>'
        return output
    return ""
```

## Overview

The `transform` command group provides tools for processing and transforming QPX data into various downstream formats. These commands enable gene annotation and protein-level quantification from feature data.

## Available Commands

- [gene-map](#gene-map) - Map genes from FASTA
- [normalize-accessions](#normalize-accessions) - Normalize protein accession formats (full ↔ bare)
- [update-metadata](#update-metadata) - Update sample/run metadata from a revised SDRF
- [quantify](#quantify) - Protein quantification via mokume (DirectLFQ, MaxLFQ, iBAQ, TopN, etc.)

---

## gene-map

Map gene information from FASTA to parquet format.

### Description {#gene-map-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_gene_map_cmd

print(generate_description(transform_gene_map_cmd))
```

### Parameters {#gene-map-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_gene_map_cmd

print(generate_params_table(transform_gene_map_cmd))
```

### Usage Examples {#gene-map-examples}

#### Basic Example {#gene-map-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_gene_map_cmd

print(generate_example(transform_gene_map_cmd, "Map gene information to parquet file:"))
```

#### With Species Parameter {#gene-map-example-species}

```bash
qpxc transform gene-map \
    --parquet-path ./output/feature.parquet \
    --fasta tests/examples/fasta/Homo-sapiens.fasta \
    --output-folder ./output \
    --species human
```

### Output Files {#gene-map-output}

- **Output**: Enhanced parquet file(s) with gene information
- **Format**: Parquet file in output folder
- **Added Fields**: Gene names and metadata from FASTA headers

### Best Practices {#gene-map-best-practices}

- Use species-specific FASTA files for accurate gene annotation
- Enable verbose mode for debugging

---

## normalize-accessions

Normalize protein accession formats between full UniProt form (`sp|ACC|NAME`) and bare form (`ACC`).

### Description {#normalize-accessions-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_normalize_accessions_cmd

print(generate_description(transform_normalize_accessions_cmd))
```

### Parameters {#normalize-accessions-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_normalize_accessions_cmd

print(generate_params_table(transform_normalize_accessions_cmd))
```

### Usage Examples {#normalize-accessions-examples}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_normalize_accessions_cmd

print(generate_example(transform_normalize_accessions_cmd, "Normalize protein accession formats:"))
```

---

## update-metadata

Update `sample.parquet` and `run.parquet` metadata from a revised SDRF file, with safety checks on protected fields.

### Description {#update-metadata-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_update_metadata_cmd

print(generate_description(transform_update_metadata_cmd))
```

### Parameters {#update-metadata-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_update_metadata_cmd

print(generate_params_table(transform_update_metadata_cmd))
```

### Usage Examples {#update-metadata-examples}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_update_metadata_cmd

print(generate_example(transform_update_metadata_cmd, "Update dataset metadata from SDRF:"))
```

---

## quantify

Compute protein-level quantification from QPX feature data using [mokume](https://github.com/bigbio/mokume).

### Description {#quantify-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_quantify_cmd

print(generate_description(transform_quantify_cmd))
```

### Parameters {#quantify-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_quantify_cmd

print(generate_params_table(transform_quantify_cmd))
```

### Supported Methods

| Method | Description | Extra Requirements |
|--------|-------------|-------------------|
| `directlfq` | DirectLFQ intensity traces (default) | `pip install "qpx[quantify]"` |
| `maxlfq` | MaxLFQ delayed normalization | `pip install "qpx[quantify]"` |
| `topn` | Average of N most intense peptides | `pip install "qpx[quantify]"` |
| `top3` | Average of 3 most intense peptides | `pip install "qpx[quantify]"` |
| `ibaq` | Intensity-Based Absolute Quantification | `pip install "qpx[quantify]"` + `--fasta` |
| `sum` | Sum of all peptide intensities | `pip install "qpx[quantify]"` |

### Usage Examples {#quantify-examples}

#### DirectLFQ (default) {#quantify-example-directlfq}

```bash
qpxc transform quantify \
    --feature-path ./qpx_output/feature.parquet \
    --method directlfq \
    -o proteins_directlfq.parquet
```

#### iBAQ (requires FASTA) {#quantify-example-ibaq}

```bash
qpxc transform quantify \
    --feature-path ./qpx_output/feature.parquet \
    --method ibaq --fasta proteome.fasta \
    -o proteins_ibaq.tsv
```

#### MaxLFQ with 8 threads {#quantify-example-maxlfq}

```bash
qpxc transform quantify \
    --feature-path ./qpx_output/feature.parquet \
    --method maxlfq --threads 8 \
    -o proteins_maxlfq.parquet
```

#### TopN with normalization {#quantify-example-topn}

```bash
qpxc transform quantify \
    --feature-path ./qpx_output/feature.parquet \
    --method topn --topn-n 5 --normalize \
    -o proteins_top5.parquet
```

### Output Files {#quantify-output}

- **Parquet**: `.parquet` files with protein-level quantification
- **TSV**: `.tsv` files (tab-separated) — determined by output file extension
- **Content**: Protein accessions, sample IDs, and quantified intensities

### Common Issues {#quantify-issues}

**Issue**: `mokume is not installed` / `DirectLFQ is not installed`

- **Solution**: Install with `pip install "qpx[quantify]"` (pulls `mokume[directlfq]>=0.1.0`)

**Issue**: `--fasta option is required for the ibaq method`

- **Solution**: Provide a FASTA database file with `--fasta`

### Best Practices {#quantify-best-practices}

- Ensure QPX feature.parquet contains valid `anchor_protein`, `intensities`, and `run_file_name` fields
- If using older QPX datasets, ensure fields have been migrated from legacy names (`precursor_charge` → `charge`, `id_scan` → `scan`)
- Decoy entries (`is_decoy=true`) and zero-intensity rows are automatically filtered
- Use `--normalize` for cross-sample normalization
- Use `--threads` to control parallelism for MaxLFQ

---

## Related Commands

- [Convert Commands](convert.md) - Convert raw data to QPX format
- [Visualization Commands](visualize.md) - Visualize transformed data
- [Statistics Commands](stats.md) - Analyze transformed data
