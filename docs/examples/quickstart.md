# Quick Start Examples

## Example 1: Convert MaxQuant Data

Convert MaxQuant output files to QPX format.

```bash
# Convert PSM data only
qpxc convert maxquant \
    --msms-file tests/examples/maxquant/maxquant_simple/msms.txt \
    --output-folder ./output \
    --structures psm

# Convert all structures (PSM + feature + PG)
qpxc convert maxquant \
    --evidence-file tests/examples/maxquant/maxquant_full/evidence.txt.gz \
    --sdrf-file tests/examples/maxquant/maxquant_full/PXD001819.sdrf.tsv \
    --protein-groups-file tests/examples/maxquant/maxquant_full/proteinGroups.txt \
    --output-folder ./output
```

**Expected Output**:

- `output/maxquant.psm.parquet` (first command)
- `output/maxquant.feature.parquet`, `output/maxquant.pg.parquet` and the metadata tables (second command)

## Example 2: Convert DIA-NN Data

Process DIA-NN report files for Data-Independent Acquisition (DIA) data.

```bash
# Convert to feature format
qpxc convert diann \
    --report-path tests/examples/diann/small/diann_report.tsv \
    --qvalue-threshold 0.05 \
    --mzml-info-folder tests/examples/diann/small/mzml \
    --sdrf-file tests/examples/diann/small/PXD019909-DIA.sdrf.tsv \
    --output-folder ./output

# Convert features + protein groups
qpxc convert diann \
    --report-path tests/examples/diann/full/diann_report.tsv.gz \
    --pg-matrix-path tests/examples/diann/full/diann_report.pg_matrix.tsv \
    --sdrf-file tests/examples/diann/full/PXD036609.sdrf.tsv \
    --output-folder ./output
```

**Expected Output**:

- `output/diann.feature.parquet`
- `output/diann.pg.parquet`

## Example 3: Protein Quantification

Compute protein-level quantification from QPX feature data using mokume.

```bash
# DirectLFQ quantification (default)
qpxc transform quantify \
    --feature-path ./output/maxquant.feature.parquet \
    --method directlfq \
    -o proteins_directlfq.parquet

# iBAQ quantification (requires FASTA)
qpxc transform quantify \
    --feature-path ./output/maxquant.feature.parquet \
    --method ibaq \
    --fasta proteome.fasta \
    -o proteins_ibaq.tsv
```

**Expected Output**:

- `proteins_directlfq.parquet` or `proteins_ibaq.tsv`

---

[← Back to Examples Overview](index.md) | [View Workflows →](workflows.md)
