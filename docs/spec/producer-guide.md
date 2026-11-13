# Producer Guide

How to write QPX files from your own tool. The [views](data-model.md) define the columns; this page says how to decide what goes into them, walks through the writers, and tells you how to check the result.

## Decisions

### Which files

| You have | Write | Notes |
|---|---|---|
| Peptide-spectrum matches | `psm` | One row per match. Keep rows above your FDR threshold only if you say so with `pass_threshold` |
| Precursor quantities per run | `feature` | One row per precursor, charge and run. A multiplet (SILAC, dimethyl, mTRAQ) is one row with several `intensities` entries |
| Protein group quantities | `pg` | One row per group, quantification unit and label |
| Sample metadata | `sample`, `run` | From an SDRF, or by hand; `run.samples` joins a label in a run to a sample |
| Provenance | `dataset`, `provenance`, `ontology` | `dataset` is required for a dataset folder |
| Which proteins a peptide maps to | `pepmap` | Every protein, where `feature.pg_accessions` has only the assigned group |
| Spectra | `mz` | Optional |

### `scan`

| Row | Value | Details |
|---|---|---|
| DDA PSM or feature | The scan number(s) of the MS2 spectrum | Set `scan_format`, or the validator cannot check them |
| Scan id that is not a scan number | Its digits, with `scan_format` declaring how to read them | [Scan Numbers](scan.md#the-scan_format-metadata-field) |
| DIA feature | The MS2 scan at the precursor's apex, else empty | [DIA](scan.md#dia) |
| Feature transferred by match-between-runs | Empty, with `identification_type` `transferred` | No spectrum of this run identified it |

### Decoys and contaminants

| Case | `is_decoy` | Accession |
|---|---|---|
| Target | `false` | As in the search database |
| Decoy | `true` | Keep the marker (`DECOY_P04114`, `REV__P04114`): a decoy never shares an accession with its target |
| Contaminant | `false`, and `contaminant` `true` on `pg` | Keep the marker (`CONTAM_P02768`) |
| No database (de novo) | `false` | Null `anchor_protein` |

### Labels

`intensities[].label` and `run.samples[].label` are the [canonical names](intensities.md#label-vocabulary), compared exactly: `TMT126`, `ITRAQ114`, `SILAC heavy`, `LFQ`. `normalize_label` in `qpx.converters.channel_labels` turns common spellings into these.

### Missing quantities

| What you know | Write | Details |
|---|---|---|
| Not identified in this label or run | No entry | |
| Identified, not quantified | Entry with null `intensity` | [Missing values](intensities.md#missing-values) |
| Measured as zero | `0.0` | Never as a stand-in for not measured |
| Imputed by the tool | The value, with `imputation` flagged | Never NaN |

### Peptidoforms

Write the [canonical form](peptidoform.md#canonical-form): fixed modifications and labels inline, `UNIMOD:n` in upper case, signed four-decimal mass shifts. `normalize_peptidoform` in `qpx.converters.ptm` rewrites other spellings. [Scope](scope.md#proforma-subset) lists the subset QPX parses.

### Ids

Let the writer derive ids. Choose a composite that names one row: the schema default fits a tool with a retention time per precursor, and a tool with one row per precursor and run declares `peptidoform, charge, run_file_name`, as DIA-NN does. Every derived id must match the [test vectors](identity.md#test-vectors), and [Row Identity](identity.md) specifies the bytes, for producers in other languages. To supply your own id, fill the column: `identity_source` then says `provided`.

## Writing a Dataset

Each writer validates on `close()` and replaces the destination only when nothing failed.

```python
from pathlib import Path

from qpx.converters.ptm import compute_precursor_mz
from qpx.writers import DatasetWriter, FeatureWriter, PgWriter, PsmWriter, RunWriter, SampleWriter

out = Path("mytool_output")
out.mkdir(exist_ok=True)

PEPTIDOFORM = "PEPTM[UNIMOD:35]IDEK"
SEQUENCE = "PEPTMIDEK"
CHARGE = 2
MZ = compute_precursor_mz(PEPTIDOFORM, CHARGE)
RUN = "run_01"
PROTEIN = "P12345"
OXIDATION = [
    {
        "name": "Oxidation",
        "accession": "UNIMOD:35",
        "positions": [{"position": 5, "amino_acid": "M", "scores": None}],
    }
]
PRECURSOR = {
    "sequence": SEQUENCE,
    "peptidoform": PEPTIDOFORM,
    "modifications": OXIDATION,
    "charge": CHARGE,
    "is_decoy": False,
    "calculated_mz": MZ,
    "observed_mz": MZ + 0.001,
    "run_file_name": RUN,
    "scan": [1234],
}

with SampleWriter(out / "mytool.sample.parquet") as writer:
    writer.write_batch([{"sample_accession": "S1", "organism": "Homo sapiens", "organism_part": "liver"}])

with RunWriter(out / "mytool.run.parquet") as writer:
    writer.write_batch(
        [
            {
                "run_accession": RUN,
                "run_file_name": RUN,
                "samples": [{"sample_accession": "S1", "label": "LFQ", "biological_replicate": 1}],
            }
        ]
    )

with PsmWriter(out / "mytool.psm.parquet", scan_format="scan") as writer:
    writer.write_batch([{**PRECURSOR, "qvalue": 0.001}])

composite = ("peptidoform", "charge", "run_file_name")
with FeatureWriter(out / "mytool.feature.parquet", identity_composite=composite, scan_format="scan") as writer:
    writer.write_batch(
        [
            {
                **PRECURSOR,
                "anchor_protein": PROTEIN,
                "pg_accessions": [{"accession": PROTEIN, "start": 10, "end": 18, "pre": "K", "post": "R"}],
                "precursor_qvalue": 0.001,
                "intensities": [{"label": "LFQ", "intensity": 1000.0}],
            }
        ]
    )

with PgWriter(out / "mytool.pg.parquet") as writer:
    writer.write_batch(
        [
            {
                "pg_accessions": [PROTEIN],
                "anchor_protein": PROTEIN,
                "grouped_runs": [RUN],
                "is_decoy": False,
                "label": "LFQ",
                "intensity": 1000.0,
                "peptides": [{"protein_name": PROTEIN, "peptide_count": 1}],
            }
        ]
    )

with DatasetWriter(out / "mytool.dataset.parquet") as writer:
    writer.write_batch([{"project_accession": "PXD000000", "creation_date": "2026-10-31"}])
```

- A record key outside the schema raises, and a value the schema cannot hold raises naming the column and row. [Serialization](serialization.md#writing-from-python) has the details.
- `write_table()` and `write_dataframe()` take Arrow tables and DataFrames, which is much faster than records.
- Use `x_<namespace>_<name>` for columns of your own ([Extensions](extensions.md)).
- To fill `psm.feature_id` first, get the feature ids from `FeatureWriter.derive_ids(table, identity_composite=composite)`.

## Checking the Result

```bash
qpxc validate --dataset-path mytool_output --conformance
```

`--conformance` runs every check, `--semantic` and `--check-identity` included, and fails on warnings as well as errors; informational notices that provided ids were not re-derived do not count as failures. It is meant as the gate for a producer's CI.

### Golden datasets

The converters' fixtures in [`tests/examples`](https://github.com/bigbio/qpx/tree/main/tests/examples) are references for how a real tool's output maps to QPX: convert one and read the output next to your own. The fixtures of `tests/conformance` start from the dataset above and change one thing per case, such as a dropped column, a duplicate id or an unknown label, and each one must fail `--conformance`. Use them to check a validator of your own in the same way.
