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

# Step 1: Convert MaxQuant data
echo "Converting MaxQuant data..."
quantmsioc convert maxquant-psm \
    --msms-file $RAW_DIR/msms.txt \
    --output-folder $OUTPUT_DIR \
    --verbose

quantmsioc convert maxquant-feature \
    --evidence-file $RAW_DIR/evidence.txt.gz \
    --sdrf-file $RAW_DIR/experiment.sdrf.tsv \
    --protein-groups-file $RAW_DIR/proteinGroups.txt \
    --output-folder $OUTPUT_DIR \
    --verbose

quantmsioc convert maxquant-pg \
    --protein-groups-file $RAW_DIR/proteinGroups.txt \
    --sdrf-file $RAW_DIR/experiment.sdrf.tsv \
    --output-folder $OUTPUT_DIR \
    --verbose

# Step 2: Calculate absolute expression
echo "Calculating absolute expression..."
quantmsioc transform ae \
    --ibaq-file $RAW_DIR/ibaq.tsv \
    --sdrf-file $RAW_DIR/experiment.sdrf.tsv \
    --output-folder $OUTPUT_DIR

# Step 3: Generate statistics
echo "Generating statistics..."
quantmsioc stats analyze psm \
    --parquet-path $OUTPUT_DIR/psm-*.psm.parquet \
    --save-path $REPORTS_DIR/psm_statistics.txt

quantmsioc stats analyze project-ae \
    --absolute-path $OUTPUT_DIR/ae-*.absolute.parquet \
    --parquet-path $OUTPUT_DIR/psm-*.psm.parquet \
    --save-path $REPORTS_DIR/project_statistics.txt

# Step 4: Create visualizations
echo "Creating visualizations..."
quantmsioc visualize plot box-intensity \
    --feature-path $OUTPUT_DIR/feature-*.feature.parquet \
    --save-path $PLOTS_DIR/intensity_boxplot.svg \
    --num-samples 20

quantmsioc visualize plot peptide-distribution \
    --feature-path $OUTPUT_DIR/feature-*.feature.parquet \
    --save-path $PLOTS_DIR/peptide_distribution.svg

quantmsioc visualize plot ibaq-distribution \
    --ibaq-path $OUTPUT_DIR/ae-*.absolute.parquet \
    --save-path $PLOTS_DIR/ibaq_distribution.svg

# Step 5: Create project metadata
echo "Creating project metadata..."
quantmsioc project create \
    --project-accession PXD001234 \
    --sdrf-file $RAW_DIR/experiment.sdrf.tsv \
    --output-folder $PROJECT_DIR \
    --software-name MaxQuant \
    --software-version 2.0.3.0

echo "Workflow complete!"
echo "Results:"
echo "  - Converted data: $OUTPUT_DIR/"
echo "  - Statistics: $REPORTS_DIR/"
echo "  - Plots: $PLOTS_DIR/"
echo "  - Project metadata: $PROJECT_DIR/project.json"
```

## Differential Expression Analysis Workflow

Analyze differential expression from MSstats output.

```bash
#!/bin/bash

# Convert differential expression data
quantmsioc transform differential \
    --msstats-file tests/examples/DE/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv \
    --sdrf-file tests/examples/DE/PXD033169.sdrf.tsv \
    --project-file tests/examples/DE/project.json \
    --fdr-threshold 0.05 \
    --output-folder ./output \
    --verbose

echo "Differential expression analysis complete!"
```

---

[← Back to Examples Overview](examples-overview.md) | [View Integration Examples →](examples-integration.md)
