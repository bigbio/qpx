# Query Commands

Query and inspect QPX datasets using SQL, filters, or quick previews.

```python exec="1" session="query_utils" result="ansi"
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

The `query` command group provides tools for querying QPX datasets. You can run arbitrary SQL queries, filter data structures by conditions, or quickly preview the first rows of any structure.

## Available Commands

- [sql](#sql) - Run arbitrary SQL queries against a dataset
- [filter](#filter) - Filter a data structure by a SQL condition
- [head](#head) - Show the first N rows of a structure

---

## sql

Run an arbitrary SQL query against a QPX dataset.

### Description {#sql-description}

```python exec="1" html="1" session="query_utils"
from qpx.cli.query import query_sql_cmd

print(generate_description(query_sql_cmd))
```

### Parameters {#sql-parameters}

```python exec="1" html="1" session="query_utils"
from qpx.cli.query import query_sql_cmd

print(generate_params_table(query_sql_cmd))
```

### Usage Examples {#sql-examples}

```python exec="1" html="1" session="query_utils"
from qpx.cli.query import query_sql_cmd

print(generate_example(query_sql_cmd, "Run SQL queries against QPX datasets:"))
```

---

## filter

Filter a QPX data structure by a SQL condition.

### Description {#filter-description}

```python exec="1" html="1" session="query_utils"
from qpx.cli.query import query_filter_cmd

print(generate_description(query_filter_cmd))
```

### Parameters {#filter-parameters}

```python exec="1" html="1" session="query_utils"
from qpx.cli.query import query_filter_cmd

print(generate_params_table(query_filter_cmd))
```

### Usage Examples {#filter-examples}

```python exec="1" html="1" session="query_utils"
from qpx.cli.query import query_filter_cmd

print(generate_example(query_filter_cmd, "Filter data structures by conditions:"))
```

---

## head

Show the first N rows of a QPX data structure.

### Description {#head-description}

```python exec="1" html="1" session="query_utils"
from qpx.cli.query import query_head_cmd

print(generate_description(query_head_cmd))
```

### Parameters {#head-parameters}

```python exec="1" html="1" session="query_utils"
from qpx.cli.query import query_head_cmd

print(generate_params_table(query_head_cmd))
```

### Usage Examples {#head-examples}

```python exec="1" html="1" session="query_utils"
from qpx.cli.query import query_head_cmd

print(generate_example(query_head_cmd, "Quick preview of data structure contents:"))
```

---

## Best Practices

- Use `--duckdb-threads` and `--duckdb-memory` to control resource usage for large datasets, and `--temp-dir` to move DuckDB's spill files off a small or read-only working directory
- Use `query head` to preview data before running complex SQL queries
- Use `query filter` for simple conditions; use `query sql` for joins and aggregations
- Export large results to Parquet format with `--output-format parquet` to preserve types
- Query relayouted files through `qpxc query`, which reads `additional_scores`, `cv_params` and `additional_intensities` as lists where DuckDB opening a file sees a struct ([Struct lists](../spec/serialization.md#struct-lists))
