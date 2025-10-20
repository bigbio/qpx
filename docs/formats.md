# Supported Formats

## Format Conversion Flow

![Format Conversion Flow](images/format_conversion.png)

## Input Formats

| Software     | PSM         | Feature      | Protein Group     |
| ------------ | ----------- | ------------ | ----------------- |
| **MaxQuant** | msms.txt    | evidence.txt | proteinGroups.txt |
| **DIA-NN**   | report.tsv  | report.tsv   | pg_matrix.tsv     |
| **FragPipe** | psm.tsv     | -            | -                 |
| **OpenMS**   | idXML       | -            | -                 |
| **mzTab**    | PSM section | PEP section  | PRT section       |

## Output Format

All conversions produce standardized **parquet files** following the quantms.io specification:

- **PSM files**: `*-psm.parquet` - Peptide-spectrum matches
- **Feature files**: `*-feature.parquet` - Quantified peptide features
- **Protein Group files**: `*-pg.parquet` - Protein quantification
- **AE files**: `*-absolute.parquet` - Absolute expression (iBAQ)
- **DE files**: `*-differential.parquet` - Differential expression

[View format specification →](format-specification.md)
