# Integration Examples

Examples of integrating qpx with other tools and frameworks.

## Nextflow Pipeline Integration

Real production Nextflow pipeline from the QPX project for processing MaxQuant PSM data.

### Pipeline Features

- Converts Thermo RAW files to mzML format using ThermoRawFileParser
- Processes MaxQuant PSM data to QPX parquet format
- Extracts spectral information from mzML files
- Supports containerization (Docker/Singularity)
- Configurable resource allocation for HPC environments
- Automatic error handling and completion reporting

### Run the Pipeline

```bash
# Using Docker
nextflow run nf-mq-psm.nf -profile docker

# Using Singularity
nextflow run nf-mq-psm.nf -profile singularity

# With custom parameters
nextflow run nf-mq-psm.nf \
    --raw_dir ./raw_files \
    --mzml_dir ./mzml_output \
    --output_dir ./results \
    --msms_file ./maxquant/msms.txt \
    --chunksize 2000000 \
    -profile docker
```

**Full pipeline source**: [nextflow/nf-mq-psm/](https://github.com/bigbio/qpx/tree/main/nextflow/nf-mq-psm)

## Python API Usage

Use qpx programmatically in Python scripts.

```python
#!/usr/bin/env python3
"""
Complete proteomics data processing using qpx Python API
"""

import logging
from pathlib import Path
from qpx.core.maxquant.maxquant import MaxQuant
from qpx.core.ae import AbsoluteExpressionHander
from qpx.operate.statistics import ParquetStatistics, IbaqStatistics
from qpx.operate.plots import (
    plot_intensity_box_of_samples,
    plot_peptide_distribution_of_protein
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    # Define paths
    raw_dir = Path("data/raw")
    output_dir = Path("output")
    plots_dir = Path("plots")

    # Create directories
    output_dir.mkdir(exist_ok=True)
    plots_dir.mkdir(exist_ok=True)

    # Step 1: Convert MaxQuant PSM data
    logger.info("Converting MaxQuant PSM data...")
    mq = MaxQuant(spectral_data=False)
    psm_output = output_dir / "psm.parquet"
    mq.process_psm_file(
        msms_path=str(raw_dir / "msms.txt"),
        output_path=str(psm_output),
        chunksize=1000000
    )

    # Step 2: Convert feature data
    logger.info("Converting MaxQuant feature data...")
    feature_output = output_dir / "feature.parquet"
    mq.process_feature_file(
        evidence_path=str(raw_dir / "evidence.txt"),
        output_path=str(feature_output),
        sdrf_path=str(raw_dir / "experiment.sdrf.tsv"),
        chunksize=1000000
    )

    # Step 3: Calculate absolute expression
    logger.info("Calculating absolute expression...")
    ae_handler = AbsoluteExpressionHander()
    ae_handler.load_ibaq_file(str(raw_dir / "ibaq.tsv"))
    ae_handler.load_sdrf_file(str(raw_dir / "experiment.sdrf.tsv"))
    ae_handler.convert_ibaq_to_quantms(
        output_folder=str(output_dir),
        output_file_prefix="ae"
    )

    # Step 4: Generate statistics
    logger.info("Generating statistics...")
    psm_stats = ParquetStatistics(str(psm_output))
    print(f"Number of proteins: {psm_stats.get_number_of_proteins()}")
    print(f"Number of peptides: {psm_stats.get_number_of_peptides()}")
    print(f"Number of PSMs: {psm_stats.get_number_of_psms()}")

    # Step 5: Create visualizations
    logger.info("Creating visualizations...")
    plot_intensity_box_of_samples(
        str(feature_output),
        str(plots_dir / "intensity_box.svg"),
        num_samples=20
    )

    plot_peptide_distribution_of_protein(
        str(feature_output),
        str(plots_dir / "peptide_dist.svg"),
        num_samples=50
    )

    logger.info("Processing complete!")
    logger.info(f"Output files: {output_dir}/")
    logger.info(f"Plots: {plots_dir}/")

if __name__ == "__main__":
    main()
```

---

[← Back to Examples Overview](examples-overview.md) | [View QC Examples →](examples-qc.md)

