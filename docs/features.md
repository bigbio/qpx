# Key Features

## Data Conversion

![Data Conversion Flow](images/data_conversion.png)

Convert data from leading proteomics software to standardized quantms.io format:

1. **MaxQuant** result files:

   - `msms.txt`: MS/MS scan information (PSM data)
   - `evidence.txt`: Peptide evidence (Feature data)
   - `proteinGroups.txt`: Protein identification results
   - `*.sdrf.tsv`: SDRF-Proteomics metadata (optional)

2. **DIA-NN** result files:

   - `report.tsv` or `report.parquet`: DIA-NN main report
   - `pg_matrix.tsv`: Protein group matrix
   - `*.sdrf.tsv`: SDRF-Proteomics metadata (optional)
   - `*ms_info.parquet`: mzML statistics folder (optional)

3. **FragPipe** result files:

   - `psm.tsv`: PSM identifications
   - `*.sdrf.tsv`: SDRF-Proteomics metadata (optional)

4. **OpenMS/idXML** files:

   - `*.idXML`: Identification results
   - `*.mzML`: Corresponding spectra files (optional, for batch processing)

5. **mzTab** (quantms pipeline) files:

   - `*.mzTab`: Standard proteomics format (PSM/PEP/PRT sections)
   - `*msstats*.csv`: MSstats/MSstatsTMT input files (optional)

**Output format**: All conversions produce standardized parquet files:

- PSM: `*-psm.parquet`
- Feature: `*-feature.parquet`
- Protein Group: `*-pg.parquet`

[View all converters →](cli-convert.md)

## Data Transformation

Transform and enrich your data:

- **Absolute Expression (AE)**: iBAQ-based absolute quantification
- **Differential Expression (DE)**: MSstats statistical analysis
- **Gene Mapping**: Add gene annotations
- **UniProt Integration**: Latest protein annotations
- **AnnData Export**: Integration with single-cell tools

[View all transformations →](cli-transform.md)

## Visualization

Generate publication-quality plots:

- Intensity distributions (KDE, box plots)
- Peptide distribution analysis
- iBAQ distribution plots
- Quality control visualizations

[View visualization tools →](cli-visualize.md)

## Statistical Analysis

Comprehensive data statistics:

- PSM data statistics (proteins, peptides, PSMs count)
- Project-level absolute expression statistics
- Sample-level quantitative metrics
- Automated report generation

[View statistical tools →](cli-stats.md)

## Project Management

Track and share your work:

- PRIDE Archive integration
- Project metadata management
- File registration and tracking

[View project tools →](cli-project.md)
