# Convert Commands

Convert various mass spectrometry data formats to the QPX standard format.

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

The `convert` command group provides converters for multiple proteomics software outputs, enabling standardization of data formats for downstream analysis. All commands generate parquet-format output files following the QPX specification.

**Shared metadata options.** The commands `diann`, `maxquant`, `fragpipe`, `mzidentml`, `openms-consensus`, `quantms-msstats`, `spectronaut` and `cdap` accept five options: `--software-name`, `--software-version`, `--software-config`, `--enrich-pride` and `--compute-integrity`. `--software-version` and `--software-config` require `--software-name`. `--enrich-pride` needs `--project-accession`. Without it, or when PRIDE cannot be reached, the command prints a warning and continues. `--compute-integrity` covers the Parquet and h5ad files, not the `.h5mu` file. The steps run in this order: conversion, `--fasta`, `--software-name`, `--enrich-pride`, MuData, `--compute-integrity`. The `--fasta` option and MuData exist only on `diann` and `openms-consensus`.

**Modifications.** The converters write `peptidoform` in ProForma with Unimod accessions (`PEPTM[UNIMOD:35]K`). Modification names (`maxquant`, `spectronaut`) and delta masses (`fragpipe`, `cdap`) are resolved with the Unimod release that ships with pyopenms (`ModificationsDB`). A name is looked up among Unimod names and synonyms, with or without the vendor site suffix (`Phospho (STY)`, `Deamidation (NQ)`, `Pyro-glu from E`), plus a few vendor names Unimod lacks (`GlyGly`, `TMT10plex`, `TMT16plex`); a name Unimod does not know stays a name in the peptidoform, with square brackets replaced by parentheses, and is logged. A name for a terminus (`Acetyl (Protein N-term)`, `Amidated (Protein C-term)`) is placed on that terminus. A mass gets the closest Unimod entry within 0.01 Da that Unimod defines at that residue or terminus. Entries of equal mass at the same site resolve only to the modification search engines are configured with (iTRAQ4plex `UNIMOD:214`, not mTRAQ:13C(3)15N(1) `UNIMOD:889`; Butyryl `UNIMOD:1289`, not Crotonaldehyde `UNIMOD:253`). Any other tie (Isovaleryl, 2-methylbutyryl and Pentanoyl on K) and a mass without an entry at that site, such as TMT6plex on Y, are kept as the signed mass (`[+999.1234]`). OpenMS' own `ModificationsDB` entries (record ids from 99900) are not used. `calculated_mz` uses the same Unimod masses.

## Available Commands

- [openms](#openms) - Enrich native OpenMS `-out_qpx` output into a full QPX dataset
- [openms-consensus](#openms-consensus) - Convert an OpenMS consensusXML (+ SDRF) to QPX (interim; pg intensity = OpenMS abundance or unnormalized unique-peptide sum)
- [diann](#diann) - Convert DIA-NN report to QPX format
- [diann-xic](#diann-xic) - Merge DIA-NN XIC exports into one QPX xic view
- [spectronaut](#spectronaut) - Convert Spectronaut report to QPX format
- [maxquant](#maxquant) - Convert MaxQuant output to QPX format
- [fragpipe](#fragpipe) - Convert FragPipe output to QPX format
- [proteome-discoverer](#proteome-discoverer) - Convert Proteome Discoverer text exports to QPX format
- [mzidentml](#mzidentml) - Convert mzIdentML file to PSM format
- [cdap](#cdap) - Convert CPTAC CDAP `.psm` files to QPX format
- [mz](#mz) - Convert an mzML spectra directory to QPX `mz.parquet` (full spectra)
- [sdrf](#sdrf) - Convert SDRF to sample and run parquet files
- [quantms-msstats](#quantms-msstats) - Convert QuantMS `*_msstats_in.csv` plus SDRF to QPX

---

## Memory {#memory}

`--max-memory` limits DuckDB's buffer pool. Above that limit DuckDB spills to `--temp-dir`, so a *lower* limit usually gives a lower peak, while the Python side of the converter comes on top in every case.

| Converter | Options | What sets the peak | If it runs out of memory |
|---|---|---|---|
| `diann` | `--max-memory` | DuckDB, plus Python columns of about one million report rows at a time. A 10× Parquet report peaked at 3.8 GB by default or with `8GB`, and at 2.4 GB with `2GB` | Lower `--max-memory` and put `--temp-dir` on a large disk. Convert a Parquet report rather than a TSV one |
| `spectronaut`, `cdap`, `quantms-msstats` | `--max-memory` | DuckDB plus one batch of Python rows | Lower `--max-memory` and put `--temp-dir` on a large disk |
| `proteome-discoverer` | `--max-memory` | DuckDB holds the PSM export twice, once per view; `Proteins.txt` is read whole in pandas | Lower `--max-memory`, then use a machine with more RAM |
| `fragpipe` | `--max-memory` | DuckDB, batches in pandas and the protein table | Lower `--max-memory`, then use a machine with more RAM |
| `maxquant` | `--max-memory`, with no real effect | Python batches of 100,000 rows and the deduplicated tables | Use a machine with more RAM; a larger `--max-memory` does not help |
| `mzidentml` | None | The parsed document and the spectra index | Use a machine with more RAM |
| `openms-consensus` | `--streaming` | The in-memory reader peaks near 7× the file size. The streaming reader is the default above 1 GB | Force `--streaming` |

Writing the output as Hive partitions does not lower the conversion's peak: partition afterwards ([Hive partitioning](../spec/serialization.md#hive-partitioning)).

## diann

Convert DIA-NN report files to QPX format.

### Description {#diann-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_diann_cmd

print(generate_description(convert_diann_cmd))
```

### Parameters {#diann-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_diann_cmd

print(generate_params_table(convert_diann_cmd))
```

### Usage Examples {#diann-examples}

#### Basic Example - Feature Data {#diann-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_diann_cmd

print(generate_example(convert_diann_cmd, "Convert a DIA-NN report with default settings:"))
```

#### Advanced Example {#diann-example-advanced}

Convert with a q-value cutoff, mzML info and resource limits:

```bash
qpxc convert diann \
    --report-path tests/examples/diann/full/diann_report.tsv.gz \
    --qvalue-threshold 0.01 \
    --mzml-info-folder tests/examples/diann/full/mzml \
    --sdrf-file tests/examples/diann/full/PXD036609.sdrf.tsv \
    --output-folder ./output \
    --max-memory 8GB \
    --max-cpus 4 \
    --verbose
```

#### Protein Groups with the PG Matrix {#diann-example-pg}

The pg view comes from the report. The pg_matrix file is optional: it fills the
MaxLFQ quantity of a label-free group the report has none for.

```bash
qpxc convert diann \
    --report-path tests/examples/diann/full/diann_report.tsv.gz \
    --pg-matrix-path tests/examples/diann/full/diann_report.pg_matrix.tsv \
    --sdrf-file tests/examples/diann/full/PXD036609.sdrf.tsv \
    --output-folder ./output \
    --max-memory 16GB \
    --max-cpus 8 \
    --verbose
```

#### XIC Traces from DIA-NN `--xic` Export {#diann-example-xic}

Merge the per-run `*.xic.parquet` files written by DIA-NN's `--xic` option into the dataset's single `xic` view:

```bash
qpxc convert diann \
    --report-path tests/examples/diann/full/diann_report.tsv.gz \
    --sdrf-file tests/examples/diann/full/PXD036609.sdrf.tsv \
    --xic-folder ./diann_xic \
    --output-folder ./output
```

### Output Files {#diann-output}

Output files generated:

- **Feature**: `{output-prefix}.feature.parquet` (always produced)
- **Protein Group**: `{output-prefix}.pg.parquet` (produced with `--pg-matrix-path`)
- **Peptide-protein map**: `{output-prefix}.pepmap.parquet` (every protein in `Protein.Ids` for each converted peptidoform)
- **Mass spectra**: `{output-prefix}.mz.parquet` (produced when `--mzml-info-folder` contains `*_ms2_info.parquet`; feature-matched MS2 spectra only)
- **XIC**: `{output-prefix}.xic.parquet` (produced when `--xic-folder` contains DIA-NN `*.xic.parquet` exports)
- **MuData**: `{output-prefix}.h5mu` (written by default; `--no-mudata` skips it)

### Common Issues {#diann-issues}

**Issue**: Out of memory errors with large files

- **Solution**: Lower `--max-memory`, do not raise it: see [Memory](#memory). A Parquet report is read in place; a TSV report is first converted once to a temporary Parquet file in the temp directory (`--temp-dir`).

**Issue**: Slow processing

- **Solution**: Increase `--max-cpus` to utilize more CPU cores

**Issue**: Missing mzML info files

- **Solution**: Ensure all mzML info TSV files are in the specified folder with correct naming

### Best Practices {#diann-best-practices}

- Use Q-value threshold of 0.05 or lower for high-confidence results
- See [Memory](#memory) for large reports, and [Hive partitioning](../spec/serialization.md#hive-partitioning) for partitioning the output afterwards
- Use verbose mode during initial testing to diagnose issues
- Ensure SDRF file correctly matches sample names in DIA-NN report
- When passing a pg_matrix, ensure it comes from the same DIA-NN run as the report

---

## diann-xic

Merge DIA-NN per-run XIC exports into one dataset-level QPX `xic.parquet`.

### Description {#diann-xic-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_diann_xic_cmd

print(generate_description(convert_diann_xic_cmd))
```

### Parameters {#diann-xic-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_diann_xic_cmd

print(generate_params_table(convert_diann_xic_cmd))
```

### Usage Example {#diann-xic-example}

```bash
qpxc convert diann-xic \
    --xic-folder ./diann_xic \
    --output-folder ./output \
    --output-prefix PXD014414 \
    --xic-selection targets \
    --xic-window apex \
    --xic-extraction-window-seconds 10 \
    --xic-qvalue-threshold 0.01
```

`--xic-folder` must contain the `*.xic.parquet` files written by DIA-NN's
`--xic` option. All runs are merged into `PXD014414.xic.parquet`; `run_file_name`
distinguishes them. The file footer records `xic_selection`, `xic_window`,
`xic_extraction_window_seconds`, `xic_qvalue_threshold`, `xic_smoothed`,
`xic_zero_points`, `xic_rt_unit` and `transition_format`. Use
`qpxc convert diann --xic-folder ...` instead when you are also converting the
DIA-NN report and want the XIC view in the same output folder.

### Output Files {#diann-xic-output}

- **XIC**: `{output-prefix}.xic.parquet`

---

## spectronaut

Convert Spectronaut report files to QPX format.

### Description {#spectronaut-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_spectronaut_cmd

print(generate_description(convert_spectronaut_cmd))
```

### Parameters {#spectronaut-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_spectronaut_cmd

print(generate_params_table(convert_spectronaut_cmd))
```

### Usage Examples {#spectronaut-examples}

#### Basic Example {#spectronaut-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_spectronaut_cmd

print(generate_example(convert_spectronaut_cmd, "Convert Spectronaut data with default settings:"))
```

#### With SDRF Metadata {#spectronaut-example-sdrf}

```bash
qpxc convert spectronaut \
    --report-path report.tsv \
    --sdrf-file data.sdrf.tsv \
    --output-folder ./qpx_output \
    --qvalue-threshold 0.01 \
    --verbose
```

#### Large Dataset with Custom DuckDB Settings {#spectronaut-example-large}

```bash
qpxc convert spectronaut \
    --report-path large_report.tsv \
    --sdrf-file data.sdrf.tsv \
    --output-folder ./qpx_output \
    --max-memory 16GB \
    --max-cpus 8 \
    --verbose
```

### Report Layouts {#spectronaut-layouts}

Both Spectronaut layouts are read: the long report and the pivot report, whose
columns are one per run and quantity (`[1] run.raw.EG.Channel1Quantity`). A
pivot report is melted to one row per precursor and run, and a run where all of
a precursor's quantities are `Filtered` has no feature row. Pivot reports
usually carry no scores, so `--qvalue-threshold` warns and applies no filter
when `EG.Qvalue` is absent.

`EG.Channel<N>Quantity` columns become one `intensities` entry per label.
Spectronaut does not name its channels: they are assumed lightest first, and
take the SDRF's SILAC labels in that order, or `SILAC light`, `medium` and
`heavy` without an SDRF. If channel 1 of your workflow is not the lightest
label, the labels are wrong. Check them against the workflow's labeling setup.

A report converts to the views it has columns for, and warns about the others.
A report with channels writes no `pg` view, because `PG.Quantity` is not per
channel. A protein-level report, without precursor columns, writes no `feature`
view.

### Output Files {#spectronaut-output}

Output files generated:

- **Feature**: `{output-prefix}.feature.parquet` (always produced)
- **Protein Group**: `{output-prefix}.pg.parquet` (always produced)
- **Ontology**: `{output-prefix}.ontology.parquet` (score CV terms)
- **Provenance**: `{output-prefix}.provenance.parquet` (conversion metadata)

### Supported Spectronaut Columns {#spectronaut-columns}

The converter auto-detects column names from the Spectronaut report, supporting both standard and custom column naming. Key mapped columns include:

| Spectronaut Column       | QPX Field                  | Notes                              |
| ------------------------ | -------------------------- | ---------------------------------- |
| `R.FileName`             | `run_file_name`            | File extensions auto-stripped      |
| `PEP.StrippedSequence`   | `sequence`                 | Bare amino acid sequence           |
| `EG.ModifiedSequence`    | `peptidoform`              | Converted to ProForma notation     |
| `FG.Charge`              | `charge`                   | Precursor charge state             |
| `FG.Quantity`            | `intensities[0].intensity` | Raw precursor intensity            |
| `EG.Qvalue`              | `additional_scores`        | Mapped to PSI-MS CV terms          |
| `PG.ProteinGroups`       | `pg_accessions`            | Semicolon-delimited protein groups |
| `PG.Quantity`            | PG `intensities`           | Protein group quantity             |

### Common Issues {#spectronaut-issues}

**Issue**: European decimal separator (comma instead of period)

- **Solution**: The converter auto-detects decimal separators by inspecting the `EG.Qvalue` column

**Issue**: Out of memory with large reports

- **Solution**: Lower `--max-memory` so DuckDB spills to disk, as [Memory](#memory) describes

**Issue**: Missing modifications in output

- **Solution**: Ensure `EG.ModifiedSequence` column is present in the Spectronaut export

### Best Practices {#spectronaut-best-practices}

- Export all required columns from Spectronaut (at minimum: `R.FileName`, `PEP.StrippedSequence`, `EG.ModifiedSequence`, `FG.Charge`, `FG.Quantity`, `PG.ProteinGroups`)
- Provide SDRF file for complete sample and run metadata
- Use Q-value threshold of 0.01 or 0.05 for high-confidence results

---

## maxquant

Convert MaxQuant output to QPX format.

### Description {#maxquant-description}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_maxquant_cmd

print(generate_description(convert_maxquant_cmd))
```

### Parameters {#maxquant-parameters}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_maxquant_cmd

print(generate_params_table(convert_maxquant_cmd))
```

### Usage Examples {#maxquant-examples}

#### Basic Example {#maxquant-example-basic}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_maxquant_cmd

print(generate_example(convert_maxquant_cmd, "Convert MaxQuant data with default settings:"))
```

#### PSM Data Only {#maxquant-example-psm}

```bash
qpxc convert maxquant \
    --msms-file tests/examples/maxquant/maxquant_simple/msms.txt \
    --output-folder ./output \
    --structures psm \
    --spectral-data \
    --output-prefix maxquant_psm
```

#### Feature Data with Protein Groups {#maxquant-example-feature}

```bash
qpxc convert maxquant \
    --evidence-file tests/examples/maxquant/maxquant_full/evidence.txt.gz \
    --protein-groups-file tests/examples/maxquant/maxquant_full/proteinGroups.txt \
    --sdrf-file tests/examples/maxquant/maxquant_full/PXD001819.sdrf.tsv \
    --output-folder ./output \
    --structures feature \
    --verbose
```

#### All Structures {#maxquant-example-all}

```bash
qpxc convert maxquant \
    --msms-file tests/examples/maxquant/maxquant_full/msms.txt.gz \
    --evidence-file tests/examples/maxquant/maxquant_full/evidence.txt.gz \
    --protein-groups-file tests/examples/maxquant/maxquant_full/proteinGroups.txt \
    --sdrf-file tests/examples/maxquant/maxquant_full/PXD001819.sdrf.tsv \
    --output-folder ./output \
    --structures psm,feature,pg \
    --verbose
```

### Modification sites {#maxquant-sites}

`--sites-file` converts a `<modification>Sites.txt` table into the [site view](../spec/site.md#quantified-sites). Repeat it for several modifications. Each experiment's runs come from `--evidence-file` or the SDRF.

```bash
qpxc convert maxquant \
    --evidence-file evidence.txt \
    --protein-groups-file proteinGroups.txt \
    --sites-file "Phospho (STY)Sites.txt" \
    --output-folder ./qpx_output
```

### Output Files {#maxquant-output}

Depending on `--structures` parameter:

- **PSM**: `{output-prefix}.psm.parquet`
- **Feature**: `{output-prefix}.feature.parquet`
- **Protein Group**: `{output-prefix}.pg.parquet`
- **Peptide-protein map**: `{output-prefix}.pepmap.parquet`, from `Proteins` of `evidence.txt`
- **Protein**: `{output-prefix}.protein.parquet`, from `proteinGroups.txt`
- **Site**: `{output-prefix}.site.parquet`, with `--sites-file`

### Common Issues {#maxquant-issues}

**Issue**: Memory errors with large evidence files

- **Solution**: Use a machine with more RAM: `--max-memory` does not change MaxQuant's peak ([Memory](#memory))

**Issue**: Missing Q-value information

- **Solution**: Provide `--protein-groups-file` for accurate Q-value mapping

**Issue**: `The SDRF TMT labels [...] match no TMT plex`

- **Solution**: Each `Reporter intensity N` column is read as the channel at that position of the TMT or iTRAQ plex (TMT6 to TMTpro 18, iTRAQ 4 and 8; partial plexes too) that holds all SDRF `comment[label]` values, so labels from different plexes cannot be mapped. Every conversion given `--sdrf-file`, a PSM-only one too, resolves the plex first and stops before writing any file.

### Best Practices {#maxquant-best-practices}

- Use `--structures` to control which output files are generated
- Always provide `--protein-groups-file` when available: features then take the protein group their `Protein group IDs` name, the membership the pg view carries
- Ensure SDRF sample names match MaxQuant experiment names
- Provide `--sdrf-file` with the PSM and feature conversions: MaxQuant leaves fixed modifications out of `Modified sequence`, so the SDRF's `MT=fixed` modifications, else those of the `parameters.txt` or `summary.txt` next to the input, are added to the peptidoforms to match `calculated_mz`. TMT, iTRAQ and SILAC labels among the SDRF's are skipped: the isobaric label comes from the SDRF channel labels, a SILAC label from `Labeling state`. A warning is logged for each fixed modification that cannot be placed, and when no source declares one
- Use compressed files (.gz) to save disk space
- Use `--spectral-data` flag if downstream analysis requires spectral information

---

## fragpipe

Convert FragPipe output to QPX format.

### Description {#fragpipe-description}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_fragpipe_cmd

print(generate_description(convert_fragpipe_cmd))
```

### Parameters {#fragpipe-parameters}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_fragpipe_cmd

print(generate_params_table(convert_fragpipe_cmd))
```

### Usage Examples {#fragpipe-examples}

#### Basic Example {#fragpipe-example-basic}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_fragpipe_cmd

print(generate_example(convert_fragpipe_cmd, "Convert FragPipe PSM data with default settings:"))
```

#### With Custom Settings {#fragpipe-example-custom}

```bash
qpxc convert fragpipe \
    --psm-file /path/to/psm.tsv \
    --output-folder ./output \
    --output-prefix fragpipe_psm
```

#### TMT and iTRAQ {#fragpipe-example-tmt}

`combined_ion.tsv` and `combined_protein.tsv` hold MS1 intensities only, so convert each plex from its own reports. Its `annotation.txt` names the channel columns, and its `psm.tsv` the raw runs:

```bash
qpxc convert fragpipe \
    --psm-file plex1/psm.tsv \
    --ion-file plex1/ion.tsv \
    --pg-file plex1/protein.tsv \
    --annotation-file plex1/annotation.txt \
    --tmt-report tmt-report/abundance_protein_MD.tsv \
    --sdrf-file study.sdrf.tsv \
    --output-folder ./output \
    --output-prefix plex1
```

The labels come from the channel of each annotation line (`127N` is `TMT127N`), and they must match the SDRF `comment[label]` values for `qpxc validate` and `design_matrix()` to resolve the samples.

### Output Files {#fragpipe-output}

- **Output**: `{output-prefix}.psm.parquet`
- **Format**: Parquet file containing PSM data
- **Schema**: Conforms to QPX PSM specification

---

## proteome-discoverer

Convert the text exports of a Proteome Discoverer consensus result to QPX format. Export `PSMs`, `Proteins` and `Input Files` from the result; the `.pdResult` file itself is not read.

### Parameters {#proteome-discoverer-parameters}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_proteome_discoverer_cmd

print(generate_params_table(convert_proteome_discoverer_cmd))
```

### Usage Examples {#proteome-discoverer-examples}

```bash
qpxc convert proteome-discoverer \
    --psm-file study_PSMs.txt \
    --proteins-file study_Proteins.txt \
    --input-files study_InputFiles.txt \
    --sdrf-file study.sdrf.tsv \
    --output-folder ./output
```

What each view takes from the exports is in [Converter coverage](../spec/converter-coverage.md#proteome-discoverer).

---

## mzidentml

Convert mzIdentML (.mzid) files to QPX PSM parquet format.

### Description {#mzidentml-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_mzidentml_cmd

print(generate_description(convert_mzidentml_cmd))
```

### Parameters {#mzidentml-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_mzidentml_cmd

print(generate_params_table(convert_mzidentml_cmd))
```

### Usage Examples {#mzidentml-examples}

#### Basic Example {#mzidentml-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_mzidentml_cmd

print(generate_example(convert_mzidentml_cmd, "Convert an mzIdentML file with default settings:"))
```

#### With Spectral Data from an MGF {#mzidentml-example-spectral}

```bash
qpxc convert mzidentml \
    --mzid-path /path/to/results.mzid \
    --mgf-path /path/to/spectra.mgf \
    --include-spectra \
    --output-folder ./output \
    --output-prefix psm_with_spectra
```

`--mgf-path` alone fills scan numbers and m/z; `--include-spectra` also attaches the whole spectrum of each PSM as `mz_array` and `intensity_array`.

### Supported Native ID Formats {#mzidentml-native-id}

The converter supports multiple native ID formats for scan number extraction:

| Format | Vendor/Source | Example |
|--------|---------------|---------|
| `scan=XXX` | Thermo | `controllerType=0 controllerNumber=1 scan=12345` |
| `cycle=XXX` | Waters/Agilent | `sample=1 period=1 cycle=1055 experiment=4` |
| `index=XXX` | Generic | `index=500` |
| `spectrum=XXX` | Various | `spectrum=999` |

### Output Files {#mzidentml-output}

- **Output**: `{output-prefix}.psm.parquet`
- **Format**: Parquet file containing PSM-level data
- **Schema**: Conforms to QPX PSM specification

### Supported mzIdentML Features {#mzidentml-features}

- **Compressed files**: Supports both `.mzid` and `.mzid.gz` formats
- **Modifications**: Full support for UNIMOD and custom modifications
- **Scores**: Extracts all CV-term scores with `higher_better` flag annotation
- **Decoy detection**: Automatic detection via `isDecoy` attribute
- **Multi-file support**: Handles mzIdentML referencing multiple spectra files

### Best Practices {#mzidentml-best-practices}

- Use compressed `.mzid.gz` files to save disk space
- Enable `--include-spectra` only when spectral arrays are needed for downstream analysis

### Common Issues {#mzidentml-issues}

**Issue**: No spectra attached from the MGF

- **Solution**: Pass `--include-spectra` with `--mgf-path`, and check that the MGF's `index=` values match the spectrum references of the mzIdentML

**Issue**: Scan numbers not extracted correctly

- **Solution**: Check if your native ID format is supported; the converter auto-detects common formats

---

## cdap

Convert CPTAC CDAP `.psm` files to QPX format.

### Description {#cdap-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_cdap_cmd

print(generate_description(convert_cdap_cmd))
```

### Parameters {#cdap-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_cdap_cmd

print(generate_params_table(convert_cdap_cmd))
```

### Usage Examples {#cdap-examples}

#### Basic Example {#cdap-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_cdap_cmd

print(generate_example(convert_cdap_cmd, "Convert one CPTAC CDAP study directory with default settings:"))
```

#### Select Output Structures {#cdap-example-structures}

```bash
qpxc convert cdap \
    --psm-dir /data/CPTAC/PDC000440 \
    --output-folder ./qpx_output/PDC000440 \
    --project-accession PDC000440 \
    --structures psm,feature
```

### Output Files {#cdap-output}

Depending on `--structures`:

- **PSM**: `{output-prefix}.psm.parquet`
- **Feature**: `{output-prefix}.feature.parquet`, one row per peptidoform, charge and fraction run, with one `intensities` entry per channel (null where CDAP reports no signal)
- **Protein Group**: `{output-prefix}.pg.parquet`, one row per group, channel and quantification unit. `grouped_runs` is the plex: from `--sdrf-file`, else the runs that differ only in a `_fNN` suffix. `pdc2qpx` takes it from PDC metadata
- **Ontology**: `{output-prefix}.ontology.parquet`
- **Provenance**: `{output-prefix}.provenance.parquet`
- **Dataset**: `{output-prefix}.dataset.parquet`

### Best Practices {#cdap-best-practices}

- Point `--psm-dir` to one CPTAC study directory containing `*.psm` files
- Use `--structures` to skip views that are not needed
- See [Memory](#memory) for full-scale CPTAC studies

---

## mz

Convert a directory of mzML spectra to a QPX `mz.parquet` (full spectra).

### Description {#mz-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_mz_cmd

print(generate_description(convert_mz_cmd))
```

### Parameters {#mz-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_mz_cmd

print(generate_params_table(convert_mz_cmd))
```

### Usage Examples {#mz-examples}

#### All MS Levels (full spectra) {#mz-example-all}

```bash
qpxc convert mz \
    --mzml-dir /data/CPTAC/PDC000109/mzml \
    --output ./qpx_output/PDC000109.mz.parquet
```

#### MS2 Only {#mz-example-ms2}

```bash
qpxc convert mz \
    --mzml-dir /data/CPTAC/PDC000109/mzml \
    --output ./PDC000109.mz.parquet \
    --ms-levels 2
```

### Output Files {#mz-output}

A single `mz.parquet` with one row per spectrum:

- `id`, `run_file_name`, `scan` — spectrum identity and linkage to PSM / feature
- `ms_level`, `scan_start_time`, `total_ion_current`, `centroid`
- `precursors` — precursor m/z, charge, isolation window (MS2+)
- `mz`, `intensity` — full peak arrays

### Best Practices {#mz-best-practices}

- Reads `.mzML` and `.mzML.gz` directly (no manual decompression needed)
- Omit `--ms-levels` (or use `1,2`) for full spectra; precursor-level LFQ
  reanalysis (e.g. quantms) needs MS1 in addition to MS2
- Each spectrum carries `run_file_name` + `scan` so it joins back to PSM / feature

---

## sdrf

Convert SDRF metadata files to QPX sample and run parquet format.

`run.acquisition_start` keeps the first valid `comment[acquisition date]` per run. The optional `--mzml-dir` reads mzML `run/@startTimeStamp` headers for runs without an SDRF date. SDRF takes precedence; text, fractional seconds and timezone suffixes are preserved without parsing, timezone conversion or added time components. This mzML source is currently available only with `qpxc convert sdrf --mzml-dir`.

### Description {#sdrf-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_sdrf_cmd

print(generate_description(convert_sdrf_cmd))
```

### Parameters {#sdrf-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_sdrf_cmd

print(generate_params_table(convert_sdrf_cmd))
```

### Usage Examples {#sdrf-examples}

#### Basic Example {#sdrf-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_sdrf_cmd

print(generate_example(convert_sdrf_cmd, "Convert SDRF metadata with default settings:"))
```

### Output Files {#sdrf-output}

- **Sample**: `{output-prefix}.sample.parquet`
- **Run**: `{output-prefix}.run.parquet`
- **Format**: Parquet files containing sample and run metadata
- **Schema**: Conforms to QPX sample and run specifications

### Best Practices {#sdrf-best-practices}

- Ensure SDRF file follows the PRIDE SDRF specifications
- Use verbose mode to diagnose parsing issues
- The converter automatically maps SDRF characteristics to QPX ontology terms

---

## quantms-msstats

Convert a QuantMS-generated `*_msstats_in.csv` table and its SDRF metadata to
QPX. The SDRF is required and is authoritative for run, sample, and label
mapping.

### Description {#quantms-msstats-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_quantms_msstats_cmd

print(generate_description(convert_quantms_msstats_cmd))
```

### Parameters {#quantms-msstats-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_quantms_msstats_cmd

print(generate_params_table(convert_quantms_msstats_cmd))
```

### Usage Example {#quantms-msstats-example}

```bash
qpxc convert quantms-msstats \
    --msstats-file PXD007683.sdrf_openms_design_msstats_in.csv \
    --sdrf-file PXD007683.sdrf.tsv \
    --output-folder ./qpx_output \
    --project-accession PXD007683
```

### Output and validation {#quantms-msstats-output}

- Produces `feature`, `sample`, `run`, `dataset`, `ontology`, and `provenance` Parquet files.
- Does not produce `psm.parquet` because MSstats does not contain complete spectrum-identification evidence.
- Does not produce `pg.parquet`; protein aggregation remains a downstream analysis step.
- Supports QuantMS `ProteinName`, `PeptideSequence`, `Charge`/`PrecursorCharge`, `Intensity`, `Run`/`Reference`, optional `RetentionTime`, and optional `Channel` columns.
- Resolves `Run` and `Reference` through SDRF aliases and rejects missing, ambiguous, or conflicting mappings.
- Collapses TMT/iTRAQ channel rows into one Feature while retaining RT and scan location when present.
- Rejects conflicting intensity values for the same Feature and channel instead of selecting or summing them.
- Ignores MSstats `IsotopeLabelType`. A labelled table needs a `Channel` column holding the SDRF labels: `H` in MSstats means a reference standard, not a heavy channel of the experiment, so it cannot be aliased to a label.

---

## openms

Enrich native OpenMS `-out_qpx` Parquet output into a full QPX dataset.

!!! warning "Deprecated"

    Use [`openms-consensus`](#openms-consensus), which reads the consensusXML and is the supported OpenMS route. What OpenMS calls a "QPX collection" is the files of one dataset, not a [collection](../spec/collection.md). OpenMS quantification tools write no run, sample or provenance view, so without this enrichment `pg` rows cannot be mapped to samples and `design_matrix` raises.

### Description {#openms-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_openms_cmd

print(generate_description(convert_openms_cmd))
```

### Parameters {#openms-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_openms_cmd

print(generate_params_table(convert_openms_cmd))
```

### Usage Examples {#openms-examples}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_openms_cmd

print(generate_example(convert_openms_cmd))
```

### Output Files {#openms-output}

- Upgrades and validates the OpenMS `psm`, `feature`, and `pg` Parquet files found in `--qpx-dir`, normalizing run references while preserving supplied IDs.
- Recovers PSM runs from exact spectrum evidence when the companion `--consensusxml` contains peptide identifications. Missing or ambiguous spectrum matches stop conversion before existing core outputs are replaced.
- Preserves records with duplicate IDs and reports warnings without deduplicating rows or inventing replacement IDs. `qpxc validate` reports these duplicates as errors.
- Generates `run`, `sample`, `ontology`, `provenance`, and `dataset` views.
- Generates a MuData file when the converted quantification views are sufficient for export.

---

## openms-consensus

### Description {#openms-consensus-description}

Convert an OpenMS `.consensusXML` (optionally with an SDRF) directly to QPX. This is the
**interim quantms production path** while OpenMS `-out_qpx` is not yet emitting
QPX format 1.1. The consensusXML carries per-run peptide-feature intensities,
PSMs, and the protein-inference graph; the SDRF supplies sample/label/fraction
metadata and the `grouped_runs` quantification units.

Several consensusXML files (for example one per sample group, as nf-core/mhcquant
writes them) can be converted into one dataset by giving `--consensusxml` a
comma-separated list.
Each file is read with its own column-to-run mapping and identification
metadata; feature, PSM and pg rows go into the same views, pg spans all inputs,
and provenance lists every input file. No run may appear in more than one input.

Runs are taken from the consensusXML column headers. A single-run map promoted
by FileConverter has no column filename; when the identifications record
exactly one primary MS run, that run's name is used. When the same group-merged
identifications were copied into every run's map (FeatureFinderIdentification
given a merged idXML, then linked), each copy is attributed to the run its
spectrum came from (`id_merge_index` into its own ProteinIdentification's
`spectra_data`), so every spectrum yields one PSM. Consensus features that are
identical in every feature identity column (two isobaric targets linked to the
same peak) yield one feature row, from the higher-quality consensus feature.

If PSM output is requested but no exportable PSM records remain (for example,
the identifications lack spectrum references), the converter logs a warning and
does not create `psm.parquet` or register it in the returned outputs or provenance.
Other requested views are still exported. A PSM-only request with no exportable
records completes with a warning and no output files.

Protein properties are taken from the group's existing `anchor_protein`:
`ProteinHit.coverage` populates `pg.sequence_coverage` (percent), and
`Posterior Probability_score` is retained in `pg.additional_scores` as
`posterior_probability` with `higher_better: true`. These describe the
representative protein, not an aggregate probability or coverage of the group.
Missing or unknown values remain null. If repeated records for that accession
contain conflicting values, the conflicting field remains null; values from
other group members are not substituted. Both consensusXML readers use this
mapping.

`pg.molecular_weight` is the theoretical average molecular weight of that same
anchor's complete, unmodified `ProteinHit.sequence`, in **kDa**. It remains null
when the sequence is absent, contains ambiguous residues such as B/Z/X or
modification notation, or conflicts across records for the anchor. U/O and the
isobaric I/L code J are supported. Other group members and identified peptides
are not used to reconstruct the sequence; this path does not load a FASTA.

`feature.pg_positions` retains the assigned peptide's known positions within
the resolved protein group. OpenMS's zero-based inclusive coordinates are
converted to QPX's one-based inclusive coordinates. Repeated evidence is
deduplicated, while distinct positions for the same protein are preserved.
Unknown positions and proteins without direct peptide evidence are omitted.
Positions come from the run's own identifications when available; otherwise
matching identifications on the consensus feature can supply sequence positions.
`feature.id_run_file_name` is populated only when the run has a resolved direct
identification of the exported peptide with a spectrum reference. Transferred
features without their own identification, unresolved source runs and conflicting
peptide assignments keep this field null.

!!! warning "Protein intensity is OpenMS' abundance or an unnormalized rollup"
    ProteomicsLFQ, IsobaricWorkflow and MS1LabeledWorkflow record their protein
    abundance in the consensusXML header (`fraction_group_level_abundance`), and
    `intensity` is that value (`quantification_method` = `openms_protein_abundance`);
    a zero is "not quantified" and becomes null. Files without one (quantms runs
    ProteinQuantifier separately) get a stopgap instead: the **unnormalized sum of
    the group's unique peptides** for that channel (the quantms `unique_peptides`
    policy, no normalization), marked `unnormalized_unique_peptide_sum`. Where
    OpenMS' abundance is present, the sum stays in `additional_intensities` under
    that name. The sum does not reproduce OpenMS' aggregation, which depends on its
    own settings.

!!! note "Quantification units"
    `grouped_runs` comes from the SDRF, else from the `fraction_group` of the map
    list. Without either, each run is its own unit, with a warning: fractions of one
    sample then stay separate, which an SDRF fixes.

!!! note "MS1 labels (SILAC, dimethyl, mTRAQ)"
    One feature row per multiplet per run, one `intensities` entry per channel. A
    two-plex calls its non-light channel `SILAC heavy` and a three-plex is
    light/medium/heavy, as the native OpenMS writer does; labels that cannot be
    resolved stop the conversion. The feature peptidoform is the peptide identity
    MS1LabeledWorkflow reports, without labels, also for dimethyl and mTRAQ.
    `observed_mz` and `rt` are the lightest measured partner's; when that is not the
    run's lightest channel, `mass_error_ppm` is null. A PSM is the partner the
    spectrum was matched with (`MS1Label:labeled_sequence`), with `psm.label`, and
    the label state is kept in `cv_params`. Ratios go to the numerator label as
    `ratio_h_l` and `ratio_h_l_normalized` (feature: `MS1Label:evidence_ratio`; pg:
    `fraction_group_level_ratio`), relative to the run's or fraction group's lowest
    channel.

### Parameters {#openms-consensus-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_openms_consensus_cmd

print(generate_params_table(convert_openms_consensus_cmd))
```

### Usage Examples {#openms-consensus-examples}

```bash
qpxc convert openms-consensus \
  --consensusxml results.consensusXML \
  --sdrf-file experiment.sdrf.tsv \
  --output-folder ./qpx_output \
  --output-prefix PXD001819
```

Several consensusXML files into one dataset:

```bash
qpxc convert openms-consensus \
  --consensusxml PBMC007_1_resolved.consensusXML,PBMC009_1_resolved.consensusXML \
  --sdrf-file PXD011628.sdrf.tsv \
  --structures feature,psm,run,sample \
  --output-folder ./qpx_output \
  --output-prefix PXD011628
```

Record the workflow, its parameters and PRIDE metadata:

```bash
qpxc convert openms-consensus \
  --consensusxml results.consensusXML \
  --sdrf-file experiment.sdrf.tsv \
  --output-folder ./qpx_output \
  --output-prefix PXD001819 \
  --project-accession PXD001819 \
  --software-name nf-core/mhcquant \
  --software-version 3.1.0 \
  --software-config params.json \
  --enrich-pride
```

### Output Files {#openms-consensus-output}

- `<prefix>.feature.parquet` — one row per `(peptidoform, charge, run, rt)` with per-run/channel intensities.
- `<prefix>.psm.parquet` — one row per spectrum match (scan, PEP, q-value, decoy).
- `<prefix>.pg.parquet` — protein groups (`pg_accessions`, `grouped_runs`, peptide/feature counts, `global_qvalue`, decoy, genes); one row per channel with a populated `label` and an `intensity` = unnormalized sum of the group's unique peptides (stamped with a `quantification_method` cv_param; null where a group has no unique-peptide signal).
- `<prefix>.run.parquet`, `<prefix>.sample.parquet` — from the SDRF (when provided).
- `<prefix>.h5mu` — the dataset's MuData view, written by default; disable with `--no-mudata`. Best-effort: a view that cannot be built is reported and the conversion still succeeds, because the Parquet views are the dataset's source of truth.

---

## Related Commands

- [Transform Commands](transform.md) - Further process converted data
- [Visualization Commands](visualize.md) - Create plots from converted data
- [Statistics Commands](stats.md) - Generate statistics from converted data
