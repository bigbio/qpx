# Validate Command

Validate QPX datasets and individual structures against their canonical schemas.

```python exec="1" session="validate_utils" result="ansi"
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
    elif isinstance(param_type, click.Choice):
        return "CHOICE"
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
            description = param.help or ""
            default_from_help = None
            if "(default:" in description.lower():
                import re

                match = re.search(r"\(default:\s*([^)]+)\)", description, re.IGNORECASE)
                if match:
                    default_from_help = match.group(1).strip()
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

The `validate` command checks QPX datasets and individual Parquet files against their canonical schemas. It verifies column presence, type matching, null values in required columns, and primary key uniqueness.

## Parameters {#validate-parameters}

```python exec="1" html="1" session="validate_utils"
from qpx.cli.validate import validate_cmd

print(generate_params_table(validate_cmd))
```

## Description {#validate-description}

```python exec="1" html="1" session="validate_utils"
from qpx.cli.validate import validate_cmd

print(generate_description(validate_cmd))
```

## Usage Examples {#validate-examples}

```python exec="1" html="1" session="validate_utils"
from qpx.cli.validate import validate_cmd

print(generate_example(validate_cmd, "Validate QPX datasets and files:"))
```

## Validation Checks

The validator performs the following checks on each structure:

| Check | Description |
|-------|-------------|
| **Column presence** | All required columns defined in the canonical schema must exist |
| **Type matching** | Column Arrow types must match the schema definition |
| **Null values** | Non-nullable columns must not contain null values |
| **Primary key uniqueness** | Primary key columns must have unique values |

When a dataset holds several views, it also checks the references between them:

| Check | Severity | Description |
|-------|----------|-------------|
| **Run references** | error | Every `run_file_name` of feature, psm and mz, and every member of `pg.grouped_runs`, is a run in `run.parquet` |
| **Label → sample** | error | Each `(run, label)` of `feature.intensities` and each `(grouped_runs, label)` of pg resolves to exactly one `run.samples` entry, and a feature row lists each label once. A feature label that matches no sample of a single-sample run is only a warning, because the query helpers fall back to that sample |
| **PSM label** | warning, error with `--strict` | `psm.label` is one of the `identified_labels` of the feature it links to |
| **ID references** | warning, error with `--strict` | `psm.feature_id`, `feature.psm_ids`, `feature.pg_ids` and the site/proteoform `pg_id`/`best_psm_id` links name existing rows |

A cross-view check whose query fails is reported as skipped, with the reason (I/O error, out of memory or possible schema drift), and fails strict validation.

## Programmatic Validation

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

## Exit Codes

| Code | Meaning |
|------|---------|
| `0` | All validated structures are valid |
| `1` | One or more structures have validation errors |
