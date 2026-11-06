# Workflow Examples

Complete processing workflows for different analysis scenarios.

## Complete Label-Free Quantification (LFQ) Workflow

Process a complete LFQ experiment from raw MaxQuant output to publication-ready results.

```bash
#!/bin/bash

# Define paths
RAW_DIR="data/raw"
OUTPUT_DIR="output"
PLOTS_DIR="plots"
REPORTS_DIR="reports"
PROJECT_DIR="project"

# Create directories
mkdir -p $OUTPUT_DIR $PLOTS_DIR $REPORTS_DIR $PROJECT_DIR

# Step 1: Convert MaxQuant data (all structures in one command)
echo "Converting MaxQuant data..."
qpxc convert maxquant \
    --msms-file $RAW_DIR/msms.txt \
    --evidence-file $RAW_DIR/evidence.txt.gz \
    --sdrf-file $RAW_DIR/experiment.sdrf.tsv \
    --protein-groups-file $RAW_DIR/proteinGroups.txt \
    --output-folder $OUTPUT_DIR \
    --verbose

# Step 2: Protein quantification (DirectLFQ)
echo "Running protein quantification..."
qpxc transform quantify \
    --feature-path $OUTPUT_DIR/diann.feature.parquet \
    --method directlfq \
    -o $OUTPUT_DIR/proteins_directlfq.parquet

# Step 3: Validate the QPX output
echo "Validating dataset..."
qpxc validate --dataset-path $OUTPUT_DIR

# Step 4: Inspect dataset info
echo "Dataset info..."
qpxc info --dataset-path $OUTPUT_DIR

echo "Workflow complete!"
echo "Results:"
echo "  - Converted data: $OUTPUT_DIR/"
```

## Protein Quantification Workflow

Compute protein-level quantification from QPX feature data.

```bash
#!/bin/bash

# Both TSV and Parquet outputs are supported — the file extension
# determines the format (.tsv for tabular inspection, .parquet for
# efficient columnar storage).

# iBAQ quantification (requires FASTA database)
qpxc transform quantify \
    --feature-path ./output/diann.feature.parquet \
    --method ibaq \
    --fasta proteome.fasta \
    -o ./output/proteins_ibaq.tsv

# MaxLFQ quantification with 8 threads
qpxc transform quantify \
    --feature-path ./output/diann.feature.parquet \
    --method maxlfq \
    --threads 8 \
    -o ./output/proteins_maxlfq.parquet

echo "Protein quantification complete!"
```

---

[← Back to Examples Overview](index.md) | [View Integration Examples →](integration.md)
