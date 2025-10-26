# Project Management Commands

Manage project metadata and file attachments for quantms.io projects.

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

The `project` command group provides tools for creating and managing project-level metadata, including integration with PRIDE Archive, SDRF file handling, and file attachment management.

## Available Commands

- [create](#create) - Create project.json from PRIDE accession
- [attach](#attach) - Attach files to project metadata

---

## create

Generate a project file from a PRIDE project accession and SDRF metadata.

### Description {#create-description}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.project import generate_pride_project_json_cmd
print(generate_description(generate_pride_project_json_cmd))
```

### Parameters {#create-parameters}

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.project import generate_pride_project_json_cmd
print(generate_params_table(generate_pride_project_json_cmd))
```

### Usage Examples {#create-examples}

#### Basic Example

```python exec="1" html="1" session="doc_utils"
from quantmsio.commands.utils.project import generate_pride_project_json_cmd
print(generate_example(generate_pride_project_json_cmd, 'Create project metadata with full parameters:'))
```

#### With Software Information

```bash
quantmsioc project create \
    --project-accession PXD016999 \
    --sdrf-file tests/examples/AE/PXD016999-first-instrument.sdrf.tsv \
    --output-folder ./project_metadata \
    --software-name MaxQuant \
    --software-version 2.0.3.0 \
    --delete-existing
```

#### Complete Workflow Example

```bash
# Create project metadata
quantmsioc project create \
    --project-accession PXD033169 \
    --sdrf-file ./PXD033169.sdrf.tsv \
    --output-folder ./project \
    --software-name OpenMS \
    --software-version 2.8.0

# The SDRF file is automatically attached to the project
echo "Project file created: ./project/project.json"
```

### Output Files

The command generates the following files in the output folder:

1. **project.json**: Main project metadata file
2. **{project_accession}.sdrf.tsv**: Copy of the SDRF file attached to the project

### project.json Structure

The generated `project.json` file contains:

```json
{
  "accession": "PXD001234",
  "title": "Project title from PRIDE",
  "description": "Project description from PRIDE",
  "organism": ["Homo sapiens"],
  "instrument": ["Q Exactive HF"],
  "quantification_method": "label free",
  "publication": {
    "title": "Publication title",
    "doi": "10.1234/journal.1234567",
    "pubmed_id": "12345678"
  },
  "samples": [
    {
      "sample_id": "Sample_001",
      "condition": "Control",
      "biological_replicate": "1"
    }
  ],
  "software": {
    "name": "MaxQuant",
    "version": "2.0.3.0"
  },
  "quantmsio_version": "1.0.0",
  "files": [
    {
      "name": "PXD001234.sdrf.tsv",
      "type": "sdrf",
      "checksum": "abc123..."
    }
  ]
}
```

### Metadata Sources

The command integrates metadata from multiple sources:

| Source                 | Information Retrieved                                                      |
| ---------------------- | -------------------------------------------------------------------------- |
| **PRIDE Archive**      | Project title, description, organism, publication info, instrument details |
| **SDRF File**          | Sample metadata, experimental design, conditions, replicates               |
| **Command Parameters** | Software information, quantmsio version                                    |

### Common Issues

**Issue**: Project not found in PRIDE

- **Solution**: Verify the project accession is correct and the project is public in PRIDE

**Issue**: Network timeout when fetching PRIDE metadata

- **Solution**: Check internet connection or try again later

**Issue**: SDRF sample names don't match data files

- **Solution**: Ensure SDRF file is correctly formatted and sample names match your data

### Best Practices {#create-best-practices}

- Run this command at the beginning of data processing to establish provenance
- Include software name and version for reproducibility
- Verify SDRF file format before running (use PRIDE SDRF validator)
- Keep the project.json file with your processed data
- Use `--delete-existing` flag carefully to avoid accidental data loss

### Validation {#create-validation}

After creating the project file, validate it:

```bash
# Check the project file was created
ls -lh ./project_metadata/project.json

# View the project metadata
cat ./project_metadata/project.json | python -m json.tool

# Verify SDRF was attached
ls -lh ./project_metadata/*.sdrf.tsv
```

---

## attach

Attach additional files to an existing project metadata file.

### Description {#attach-description}

Adds references to data files in the project.json metadata. This command is useful for tracking all files associated with a project.

### Parameters {#attach-parameters}

| Parameter            | Type   | Required | Default | Description                                                                        |
| -------------------- | ------ | -------- | ------- | ---------------------------------------------------------------------------------- |
| `--project-file`     | Path   | Yes      | -       | Existing project.json file path                                                    |
| `--attach-file`      | Path   | Yes      | -       | File to attach to the project                                                      |
| `--category`         | Choice | Yes      | -       | File category: sdrf-file, psm-file, feature-file, absolute-file, differential-file |
| `--is-folder`        | Flag   | No       | False   | Indicates if the file is a folder                                                  |
| `--partitions`       | String | No       | -       | Fields used for splitting files, separated by comma                                |
| `--replace-existing` | Flag   | No       | False   | Whether to delete old files                                                        |

### Usage Examples {#attach-examples}

#### Attach PSM File

```bash
quantmsioc project attach \
    --project-file ./project/project.json \
    --attach-file ./output/psm-abc123.psm.parquet \
    --category psm-file
```

#### Attach Multiple Files

```bash
# Attach PSM file
quantmsioc project attach \
    --project-file ./project/project.json \
    --attach-file ./output/psm.parquet \
    --category psm-file

# Attach feature file
quantmsioc project attach \
    --project-file ./project/project.json \
    --attach-file ./output/feature.parquet \
    --category feature-file

# Attach absolute expression file
quantmsioc project attach \
    --project-file ./project/project.json \
    --attach-file ./output/ae.parquet \
    --category absolute-file
```

#### Complete Project Assembly

```bash
#!/bin/bash

PROJECT_FILE="./project/project.json"

# Create project metadata
quantmsioc project create \
    --project-accession PXD001234 \
    --sdrf-file ./metadata.sdrf.tsv \
    --output-folder ./project \
    --software-name MaxQuant \
    --software-version 2.0.3.0

# Attach all processed files
for file in ./output/*.parquet; do
    # Determine file category from filename
    if [[ $file == *"psm"* ]]; then
        category="psm-file"
    elif [[ $file == *"feature"* ]]; then
        category="feature-file"
    elif [[ $file == *"absolute"* ]]; then
        category="absolute-file"
    elif [[ $file == *"differential"* ]]; then
        category="differential-file"
    else
        continue
    fi

    quantmsioc project attach \
        --project-file "$PROJECT_FILE" \
        --attach-file "$file" \
        --category "$category"
done

echo "All files attached to project"
```

### File Categories

Supported file category values:

| Category            | Description                         |
| ------------------- | ----------------------------------- |
| `psm-file`          | Peptide-spectrum match data         |
| `feature-file`      | Feature-level quantification        |
| `absolute-file`     | Absolute expression data            |
| `differential-file` | Differential expression results     |
| `sdrf-file`         | Sample and data relationship format |

### Output

The command updates the project.json file by adding a file entry to the appropriate section based on the category.

### Best Practices {#attach-best-practices}

- Attach files immediately after creating them to maintain accurate file tracking
- Use correct file categories for proper organization
- Keep project.json file backed up as it tracks all project data

---

## Project Metadata Best Practices

### Complete Project Setup Workflow

```bash
#!/bin/bash

# Define variables
PROJECT_ID="PXD001234"
SDRF_FILE="./metadata/${PROJECT_ID}.sdrf.tsv"
OUTPUT_DIR="./processed"
PROJECT_DIR="./project"

# Step 1: Create project metadata
echo "Creating project metadata..."
quantmsioc project create \
    --project-accession "$PROJECT_ID" \
    --sdrf-file "$SDRF_FILE" \
    --output-folder "$PROJECT_DIR" \
    --software-name MaxQuant \
    --software-version 2.0.3.0

# Step 2: Process data (example with MaxQuant)
echo "Processing data..."
quantmsioc convert maxquant-psm \
    --msms-file ./raw/msms.txt \
    --output-folder "$OUTPUT_DIR"

quantmsioc convert maxquant-feature \
    --evidence-file ./raw/evidence.txt \
    --sdrf-file "$SDRF_FILE" \
    --output-folder "$OUTPUT_DIR"

quantmsioc convert maxquant-pg \
    --protein-groups-file ./raw/proteinGroups.txt \
    --sdrf-file "$SDRF_FILE" \
    --output-folder "$OUTPUT_DIR"

# Step 3: Attach all processed files
echo "Attaching processed files to project..."
for file in "$OUTPUT_DIR"/*.parquet; do
    filename=$(basename "$file")

    # Determine category
    if [[ $filename == *"psm"* ]]; then
        category="psm-file"
    elif [[ $filename == *"feature"* ]]; then
        category="feature-file"
    fi

    quantmsioc project attach \
        --project-file "$PROJECT_DIR/project.json" \
        --attach-file "$file" \
        --category "$category"
done

echo "Project setup complete!"
echo "Project metadata: $PROJECT_DIR/project.json"
```

### Version Control Integration

Track project metadata with git:

```bash
# Initialize git repository for project
cd ./project
git init
git add project.json *.sdrf.tsv
git commit -m "Initial project metadata for $PROJECT_ID"

# After attaching files
git add project.json
git commit -m "Attached processed data files"

# Tag releases
git tag -a v1.0 -m "Initial data release"
```

### Data Sharing

Prepare project for sharing:

```bash
#!/bin/bash

PROJECT_DIR="./project"
ARCHIVE_NAME="project_data_$(date +%Y%m%d).tar.gz"

# Create archive with project metadata and data
tar -czf "$ARCHIVE_NAME" \
    "$PROJECT_DIR/project.json" \
    "$PROJECT_DIR"/*.sdrf.tsv \
    ./output/*.parquet

echo "Project archive created: $ARCHIVE_NAME"
echo "SHA256: $(sha256sum $ARCHIVE_NAME)"
```

---

## Related Commands

- [Convert Commands](cli-convert.md) - Generate data files to attach to projects
- [Transform Commands](cli-transform.md) - Process data for project workflows
- [Statistics Commands](cli-stats.md) - Generate project statistics
