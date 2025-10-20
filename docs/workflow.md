# Example Workflow

## Complete Workflow Overview

![Complete Workflow](images/workflow.png)

The complete quantmsio workflow includes data conversion, statistical analysis, transformation, and visualization - from raw data to publication-ready plots.

## Detailed Example

A complete example of processing proteomics data with quantmsio:

```bash
# 1. Convert MaxQuant output
quantmsioc convert maxquant-psm \
    --msms-file data/msms.txt \
    --output-folder ./output

quantmsioc convert maxquant-feature \
    --evidence-file data/evidence.txt.gz \
    --sdrf-file data/experiment.sdrf.tsv \
    --protein-groups-file data/proteinGroups.txt \
    --output-folder ./output

# 2. Calculate absolute expression
quantmsioc transform ae \
    --ibaq-file data/ibaq.tsv \
    --sdrf-file data/experiment.sdrf.tsv \
    --output-folder ./output

# 3. Generate statistics
quantmsioc stats analyze psm \
    --parquet-path ./output/psm.parquet \
    --save-path ./reports/statistics.txt

# 4. Create visualizations
quantmsioc visualize plot box-intensity \
    --feature-path ./output/feature.parquet \
    --save-path ./plots/intensity_distribution.svg

# 5. Create project metadata
quantmsioc project create \
    --project-accession PXD001234 \
    --sdrf-file data/experiment.sdrf.tsv \
    --output-folder ./project
```

[View more examples →](examples-overview.md)
