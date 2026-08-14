# Converter Coverage Matrix

This page lists which QPX data views each converter produces. Use it to see at a glance what output to expect from each tool.

## Views produced per converter

| Converter | PSM | Feature | PG | Pepmap | Sample | Run | Dataset | Ontology | Provenance | mz |
|-----------|:---:|:-------:|:--:|:------:|:------:|:---:|:-------:|:--------:|:----------:|:--:|
| **MaxQuant** | Yes | Yes | Yes | No | If SDRF | If SDRF | Yes | If SDRF | No | No |
| **FragPipe** | Yes | Yes | Yes | No | If SDRF | If SDRF | Yes | If SDRF | No | No |
| **DIA-NN** | No | Yes | Yes | No | Yes | Yes | Yes | If terms | Yes | No |
| **Spectronaut** | No | Yes | Yes | No | If SDRF | If SDRF | Yes | Yes | Yes | No |
| **OpenMS native QPX** | Yes | Yes | Yes | No | Yes | Yes | Yes | Yes | Yes | No |
| **OpenMS consensusXML** | Yes | Yes | Yes | No | If SDRF | If SDRF | Yes | If terms | Yes | No |
| **CDAP** | Yes | Yes | Yes | No | If PDC | If PDC | Yes | Yes | Yes | No |
| **mzIdentML** | Yes | No | No | Yes | If SDRF | If SDRF | Yes | Yes | Yes | No |
| **QuantMS MSstats** | No | Yes | No | No | Yes | Yes | Yes | Yes | Yes | No |
| **SDRF** | No | No | No | No | Yes | Yes | No | Optional | No | No |

- **Yes** — the converter produces this view.
- **No** — the converter does not produce this view (e.g. DIA-NN has no PSM view; mzIdentML has no Feature/PG).
- **If SDRF** — the view is produced only when an SDRF file is provided. Whether SDRF is optional or required depends on the converter; see its CLI contract below.
- **If terms** — the view is produced when the converter discovers resolvable ontology entries.
- **If PDC** — CPTAC/PDC studies ship no SDRF, and CDAP `.psm` files carry no sample metadata. When run through `qpxc pdc2qpx` (default), the sample/run views are built from PDC GraphQL metadata, which also recovers the TMT/iTRAQ channel → biological-sample mapping. Disable with `--no-metadata`.

> The `mz` (full-spectra) view is produced by the standalone `qpxc convert mz` command, or automatically by `qpxc pdc2qpx --include-spectra`; it is not emitted by the per-tool converters above.

## CLI commands

| Converter | Command |
|-----------|---------|
| MaxQuant | `qpxc convert maxquant` |
| FragPipe | `qpxc convert fragpipe` |
| DIA-NN | `qpxc convert diann` |
| Spectronaut | `qpxc convert spectronaut` |
| OpenMS native QPX | `qpxc convert openms` |
| OpenMS consensusXML | `qpxc convert openms-consensus` |
| CDAP | `qpxc convert cdap` |
| mzIdentML | `qpxc convert mzidentml` |
| QuantMS MSstats | `qpxc convert quantms-msstats` |
| SDRF only | `qpxc convert sdrf` |

## Input files (summary)

| Converter | Typical inputs |
|-----------|----------------|
| MaxQuant | msms.txt, evidence.txt, proteinGroups.txt |
| FragPipe | psm.tsv, combined_ion, combined_protein |
| DIA-NN | report (TSV or Parquet), pg_matrix (optional); required SDRF |
| Spectronaut | report.tsv; optional SDRF |
| OpenMS native QPX | `-out_qpx` Parquet directory; SDRF; optional companion consensusXML |
| OpenMS consensusXML | consensusXML; optional SDRF |
| CDAP | CPTAC CDAP `.psm` files in one study directory |
| mzIdentML | .mzid / .mzid.gz; optional MGF or mzML (file/folder) for spectra; optional SDRF |
| QuantMS MSstats | QuantMS-generated `*_msstats_in.csv`; required SDRF |
| SDRF | Single SDRF TSV file |

For field-level mappings from each tool’s columns to QPX, see [Tool Field Mappings](tool-mappings.md).
