# pdc2qpx

One-shot download + conversion of PDC/CPTAC studies into QPX datasets.

`pdc2qpx` orchestrates these steps in a single command:

1. Download each study's CDAP `.psm` files (and, with `--include-spectra`, the mzML spectra) from PDC via [pridepy](https://github.com/bigbio/pridepy).
2. Convert the `.psm` files to QPX psm/feature/pg/dataset/ontology/provenance views (`CdapConverter`).
3. Build `sample`/`run` views from PDC GraphQL metadata (on by default) — this recovers the TMT/iTRAQ channel → biological-sample mapping that `.psm` files lack. Disable with `--no-metadata`.
4. With `--include-spectra`, convert the mzML files to a full-spectra `mz.parquet`, linked to PSM/feature by `run_file_name` + `scan`.

The download layer requires the optional `pdc` extra:

```bash
pip install qpx[pdc]
```

## Usage

```bash
qpxc pdc2qpx -a PDC000109 --download-dir ./downloads --output-folder ./qpx/PDC000109 [OPTIONS]
```

`-a/--accession` accepts a **single** ID, **comma-separated** IDs, or a **CSV** with a `pdc_study_id`/`pdc_id` column (the same input forms as `pridepy download-pdc-files -a`). With more than one study, each is written to `<output-folder>/<study>/` and a failed study (e.g. one with no CDAP `.psm`) is reported without aborting the rest — unless `--stop-on-error` is set.

## Parameters

| Option | Description | Default |
| ------ | ----------- | ------- |
| `-a`, `--accession` | PDC study accession(s): single ID, comma-separated, or a CSV with `pdc_study_id`/`pdc_id` | required |
| `--download-dir` | Directory for downloaded files (written under `<download-dir>/<study>/`) | required |
| `--output-folder` | QPX output dir. One study writes here; multiple studies write to `<output-folder>/<study>/` | required |
| `--include-spectra` | Also download mzML and produce a full-spectra `<study>.mz.parquet` | off |
| `--no-metadata` | Skip building `sample`/`run` views from PDC metadata | off (views built) |
| `--ms-levels` | MS levels for the mz view (e.g. `2` or `1,2`) | all |
| `--max-cpus` | Threads for the CDAP conversion | one per CPU (DuckDB default) |
| `--max-memory` | Memory limit for the CDAP conversion | 80% of RAM (DuckDB default) |
| `--skip-download` | Use files already present under `<download-dir>/<study>/` | off |
| `--stop-on-error` | With multiple studies, abort on the first failure | off (continue) |

## Examples

### One study, entire QPX including full spectra

For quantms-style reanalysis you need MS1 (precursor LFQ quantification) as well as MS2 (identification), so include all MS levels (the default when `--ms-levels` is omitted):

```bash
qpxc pdc2qpx -a PDC000109 \
    --download-dir ./downloads \
    --output-folder ./qpx/PDC000109 \
    --include-spectra --max-cpus 24
```

### Many studies from a CSV

Each study is written to `./qpx/<study>/`; studies without CDAP `.psm` are skipped and reported at the end.

```bash
qpxc pdc2qpx -a cptac_lfq.csv \
    --download-dir ./downloads \
    --output-folder ./qpx --max-cpus 24
```

### Re-run on already-downloaded files

```bash
qpxc pdc2qpx -a PDC000109 \
    --download-dir ./downloads \
    --output-folder ./qpx/PDC000109 \
    --include-spectra --skip-download
```

## Output

Per study, prefixed with the study accession (under `--output-folder` for a single study, or `<output-folder>/<study>/` for many):

- `<study>.psm.parquet`, `<study>.feature.parquet`, `<study>.pg.parquet`
- `<study>.sample.parquet`, `<study>.run.parquet` (unless `--no-metadata`) — biological sample metadata and the channel → sample mapping from PDC
- `<study>.dataset.parquet`, `<study>.ontology.parquet`, `<study>.provenance.parquet`
- `<study>.mz.parquet` (only with `--include-spectra`) — full spectra, joinable to PSM/feature by `run_file_name` + `scan`

## Notes

- Only PDC studies that publish CDAP `.psm` files can produce the base QPX views; `--include-spectra` additionally requires mzML to be available on PDC for that study.
- `sample`/`run` views come from PDC GraphQL metadata (`studyExperimentalDesign` + `biospecimenPerStudy` + `fileMetadata`), not from the `.psm` files. A metadata fetch failure is logged and skipped without breaking the rest of the conversion.
- Downloads go through pridepy's PDC layer (signed CloudFront URLs, md5 validation, automatic 403 URL re-signing).
- See [Converters](convert.md) for the underlying `qpxc convert cdap` and `qpxc convert mz` commands.
