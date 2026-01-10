# QPX Data Model Review: DIA-NN Coverage Analysis

**Date**: 2026-01-10
**Version**: QPX 1.0
**Reviewer**: Claude Code

## Executive Summary

This document reviews the QPX data model's capability to represent all DIA-NN original file format columns. The analysis covers the current mappings in `qpx/core/common.py` and `qpx/core/diann/diann.py` against the complete DIA-NN output specification.

### Overall Assessment

**Status: Partially Complete - Ready for Core Use Cases, Gaps Exist for Advanced Features**

The QPX data model covers the **essential** DIA-NN columns required for standard proteomics workflows (identification, quantification, filtering). However, there are notable gaps in:
1. Ion mobility data fields
2. Quality/scoring metrics used for advanced analysis
3. MS1-based quantification columns
4. Library-referenced values

---

## 1. Currently Mapped DIA-NN Columns

### 1.1. Feature-Level Mappings (DIANN_MAP)

| DIA-NN Column | QPX Field | Status |
|---------------|-----------|--------|
| `Precursor.Quantity` | `intensity` | Mapped |
| `RT.Start` | `rt_start` | Mapped |
| `RT.Stop` | `rt_stop` | Mapped |
| `RT` | `rt` | Mapped |
| `Predicted.RT` | `predicted_rt` | Mapped |
| `Protein.Group` | `pg_accessions` | Mapped |
| `Protein.Ids` | `mp_accessions` | Mapped |
| `PEP` | `posterior_error_probability` | Mapped |
| `Global.Q.Value` | `global_qvalue` | Mapped |
| `Global.PG.Q.Value` | `pg_global_qvalue` | Mapped |
| `Q.Value` | `qvalue` | Mapped (via additional_scores) |
| `PG.Q.Value` | `pg_qvalue` | Mapped (via additional_scores) |
| `Precursor.Normalised` | `normalize_intensity` | Mapped |
| `PG.MaxLFQ` | `lfq` | Mapped |
| `Quantity.Quality` | `precursor_quantification_score` | Mapped (via cv_params) |
| `Precursor.Charge` | `precursor_charge` | Mapped |
| `Stripped.Sequence` | `sequence` | Mapped |
| `Modified.Sequence` | `peptidoform` | Mapped |
| `Genes` | `gg_names` | Mapped |
| `Run` | `reference_file_name` | Mapped |

### 1.2. Protein Group Mappings (DIANN_PG_MAP)

| DIA-NN Column | QPX Field | Status |
|---------------|-----------|--------|
| `Protein.Group` | `pg_accessions` | Mapped |
| `Protein.Names` | `pg_names` | Mapped |
| `Genes` | `gg_accessions` | Mapped |
| `Run` | `reference_file_name` | Mapped |
| `Global.PG.Q.Value` | `global_qvalue` | Mapped |
| `PG.MaxLFQ` | `lfq` | Mapped |
| `PG.Q.Value` | `qvalue` | Mapped |
| `Proteotypic` | `proteotypic` | Mapped (internal use) |
| `Stripped.Sequence` | `stripped_sequence` | Mapped (internal use) |
| `Precursor.Id` | `precursor_id` | Mapped (internal use) |

---

## 2. Missing DIA-NN Columns

### 2.1. Ion Mobility Fields (HIGH PRIORITY for TIMS-TOF data)

| DIA-NN Column | Description | QPX Support | Recommendation |
|---------------|-------------|-------------|----------------|
| `IM` | Measured ion mobility value | **NOT MAPPED** | Add to feature schema |
| `iIM` | Reference ion mobility from spectral library | **NOT MAPPED** | Add to cv_params or new field |
| `Predicted.IM` | Predicted ion mobility based on iIM | **NOT MAPPED** | Add to feature schema |
| `Predicted.iIM` | iIM predicted based on measured IM | **NOT MAPPED** | Add to cv_params |

**Note**: The QPX schema has `ion_mobility`, `start_ion_mobility`, and `stop_ion_mobility` fields defined, but the DIA-NN converter does not map the actual DIA-NN ion mobility columns (`IM`, `iIM`, `Predicted.IM`, `Predicted.iIM`). Currently hardcoded to `None` in `diann.py:634`.

### 2.2. Quality/Scoring Metrics (MEDIUM PRIORITY)

| DIA-NN Column | Description | QPX Support | Recommendation |
|---------------|-------------|-------------|----------------|
| `CScore` | Combined score for identification | **NOT MAPPED** | Add to additional_scores |
| `Evidence` | Evidence score component | **NOT MAPPED** | Add to additional_scores |
| `Spectrum.Similarity` | Similarity to reference spectrum | **NOT MAPPED** | Add to additional_scores |
| `Averagine` | MS1 isotopic profile quality | **NOT MAPPED** | Add to additional_scores |
| `Mass.Evidence` | Mass accuracy evidence | **NOT MAPPED** | Add to additional_scores |
| `Decoy.Evidence` | Decoy-specific evidence score | **NOT MAPPED** | Add to additional_scores |
| `Decoy.Score` | Decoy scoring metric | **NOT MAPPED** | Add to additional_scores |

### 2.3. MS1 Quantification (MEDIUM PRIORITY)

| DIA-NN Column | Description | QPX Support | Recommendation |
|---------------|-------------|-------------|----------------|
| `Ms1.Area` | Non-normalized MS1 peak area | **NOT MAPPED** | Add to additional_intensities |
| `Ms1.Translated` | Translated MS1 quantity | **NOT MAPPED** | Add to additional_intensities |

### 2.4. Library Reference Values (LOW PRIORITY)

| DIA-NN Column | Description | QPX Support | Recommendation |
|---------------|-------------|-------------|----------------|
| `Lib.Q.Value` | Q-value from spectral library | **NOT MAPPED** | Add to additional_scores |
| `Lib.PG.Q.Value` | Protein group q-value from library | **NOT MAPPED** | Add to additional_scores |
| `Lib.PEP` | PEP from spectral library | **NOT MAPPED** | Add to additional_scores |
| `Lib.Index` | Index in spectral library | **NOT MAPPED** | Optional - internal use |

### 2.5. Fragment-Level Data (REMOVED in DIA-NN 2.0)

| DIA-NN Column | Description | QPX Support | Notes |
|---------------|-------------|-------------|-------|
| `Fragment.Quant.Raw` | Raw fragment quantities | N/A | Removed in DIA-NN 2.0 |
| `Fragment.Correlations` | Fragment correlation scores | N/A | Removed in DIA-NN 2.0 |

### 2.6. Identification Metadata (LOW PRIORITY)

| DIA-NN Column | Description | QPX Support | Recommendation |
|---------------|-------------|-------------|----------------|
| `Precursor.Id` | Unique precursor identifier | **NOT MAPPED** | Consider for feature schema |
| `First.Protein.Description` | Protein description text | **NOT MAPPED** | Optional metadata |
| `Decoy` | Decoy indicator (DIA-NN 2.0+) | **NOT MAPPED** | Map to is_decoy field |

### 2.7. DIA-NN 2.0 PG Matrix Additions

| DIA-NN Column | Description | QPX Support | Recommendation |
|---------------|-------------|-------------|----------------|
| `N.Sequences` | Number of sequences per protein | **NOT MAPPED** | Map to peptide_counts |
| `N.Proteotypic.Sequences` | Number of proteotypic sequences | **NOT MAPPED** | Map to peptide_counts |

---

## 3. QPX Schema Capabilities Not Utilized for DIA-NN

The QPX schema defines fields that could accommodate DIA-NN data but are not currently populated:

### 3.1. In Feature Schema (format.py)

```python
# Defined but set to None for DIA-NN:
- ion_mobility          # Could map IM
- start_ion_mobility    # Could derive from IM range
- stop_ion_mobility     # Could derive from IM range
- gg_accessions         # Set to None, could parse from Genes
- scan_reference_file_name  # Set to None
```

### 3.2. In Additional Scores Structure

The `additional_scores` field is a flexible array that could accommodate all DIA-NN quality metrics:

```python
additional_scores: [
    {"score_name": "qvalue", "score_value": 0.01},
    {"score_name": "pg_qvalue", "score_value": 0.02},
    {"score_name": "global_qvalue", "score_value": 0.005},
    # Missing but could add:
    {"score_name": "cscore", "score_value": 2.5},
    {"score_name": "evidence", "score_value": 1.8},
    {"score_name": "spectrum_similarity", "score_value": 0.95},
    {"score_name": "averagine", "score_value": 0.85},
    {"score_name": "mass_evidence", "score_value": 1.2},
]
```

### 3.3. In Additional Intensities Structure

```python
additional_intensities: [
    {
        "sample_accession": "Sample-1",
        "channel": "LFQ",
        "intensities": [
            {"intensity_name": "lfq", "intensity_value": 1234.5},
            # Missing but could add:
            {"intensity_name": "ms1_area", "intensity_value": 5678.9},
            {"intensity_name": "ms1_translated", "intensity_value": 4567.8},
        ]
    }
]
```

---

## 4. Recommendations

### 4.1. High Priority (Essential for TIMS-TOF workflows)

1. **Map Ion Mobility Fields**: Update `DIANN_MAP` to include:
   ```python
   "IM": "ion_mobility",
   ```

2. **Update DIA-NN Converter**: Remove the hardcoded `None` assignments for ion mobility in `diann.py:634-636` and properly extract these values.

### 4.2. Medium Priority (Enhanced quality metrics)

3. **Expand DIANN_MAP for Quality Scores**:
   ```python
   # Add to DIANN_USECOLS for reading
   DIANN_QUALITY_COLS = [
       "CScore", "Evidence", "Spectrum.Similarity",
       "Averagine", "Mass.Evidence", "Ms1.Area"
   ]
   ```

4. **Populate additional_scores** with DIA-NN quality metrics during conversion.

5. **Add MS1 Quantification**: Map `Ms1.Area` to `additional_intensities` with intensity_name `"ms1_area"`.

### 4.3. Low Priority (Nice-to-have)

6. **Library Reference Values**: Add `Lib.Q.Value`, `Lib.PG.Q.Value` to additional_scores for users doing library quality assessment.

7. **DIA-NN 2.0 Decoy Field**: Map the new `Decoy` column to `is_decoy` field.

---

## 5. DIA-NN Version Compatibility

### DIA-NN 1.x vs 2.0 Changes

| Change | Impact on QPX |
|--------|---------------|
| `File.Name` replaced by `Run` | Already handled (uses `Run`) |
| `PG.Quantity` and `PG.Normalised` removed | N/A - using MaxLFQ instead |
| `Fragment.Quant.Raw` removed | N/A - not previously mapped |
| `Fragment.Correlations` removed | N/A - not previously mapped |
| `Decoy` column added | Should map to `is_decoy` |
| Report format changed to Parquet | Already supported via DuckDB |
| `N.Sequences` added to pg_matrix | Should map to peptide_counts |

---

## 6. Gap Analysis Summary

### Coverage Statistics

| Category | Mapped | Missing | Coverage |
|----------|--------|---------|----------|
| Core Identification | 10/10 | 0 | 100% |
| Quantification | 4/6 | 2 | 67% |
| Quality Scores | 4/12 | 8 | 33% |
| Ion Mobility | 0/4 | 4 | 0% |
| Library Reference | 0/4 | 4 | 0% |
| **Overall** | **18/36** | **18** | **50%** |

### Priority Matrix

| Priority | Fields | Effort | Impact |
|----------|--------|--------|--------|
| **High** | Ion Mobility (4 fields) | Medium | High for TIMS-TOF |
| **Medium** | Quality Scores (8 fields) | Low | Medium for QC |
| **Medium** | MS1 Quantification (2 fields) | Low | Medium for MS1 users |
| **Low** | Library Reference (4 fields) | Low | Low - specialized use |

---

## 7. Conclusion

The QPX data model is **ready for standard DIA-NN workflows** involving:
- Peptide/protein identification
- Label-free quantification (MS2-based)
- FDR filtering and quality control
- Integration with sample metadata via SDRF

**Key Gaps** exist for:
1. **TIMS-TOF/Ion Mobility data** - Critical for Bruker instrument users
2. **Quality scoring metrics** - Useful for advanced QC and troubleshooting
3. **MS1-based quantification** - For users preferring MS1 over MS2 quantities

The QPX schema architecture (additional_scores, additional_intensities, cv_params) provides sufficient flexibility to accommodate all missing DIA-NN columns without schema changes. The primary work required is in the DIA-NN converter (`diann.py`) to read and map these additional columns.

---

## Appendix A: Complete DIA-NN Column Reference

Based on documentation and GitHub issues, the known DIA-NN report columns include:

### Identification
- `Precursor.Id`, `Modified.Sequence`, `Stripped.Sequence`
- `Protein.Group`, `Protein.Ids`, `Protein.Names`
- `Genes`, `Proteotypic`
- `Precursor.Charge`

### Quality/Filtering
- `Q.Value`, `PG.Q.Value`, `Global.Q.Value`, `Global.PG.Q.Value`
- `PEP`, `CScore`, `Evidence`
- `Spectrum.Similarity`, `Averagine`, `Mass.Evidence`
- `Decoy.Evidence`, `Decoy.Score`, `Decoy` (2.0+)
- `Lib.Q.Value`, `Lib.PG.Q.Value`, `Lib.PEP`, `Lib.Index`

### Quantification
- `Precursor.Quantity`, `Precursor.Normalised`
- `PG.Quantity` (removed 2.0), `PG.Normalised` (removed 2.0), `PG.MaxLFQ`
- `Ms1.Area`, `Ms1.Translated`
- `Quantity.Quality`
- `Fragment.Quant.Raw` (removed 2.0), `Fragment.Correlations` (removed 2.0)

### Retention Time
- `RT`, `RT.Start`, `RT.Stop` (or `RT.End`), `Predicted.RT`

### Ion Mobility
- `IM`, `iIM`, `Predicted.IM`, `Predicted.iIM`

### File Reference
- `Run` (replaces `File.Name` in 2.0)

### Protein Group Matrix (report.pg_matrix.tsv)
- `Protein.Group`, `Protein.Names`, `Genes`
- `N.Sequences` (2.0+), `N.Proteotypic.Sequences` (2.0+)
- Sample columns with MaxLFQ values

---

## Appendix B: Sources

- [DIA-NN GitHub Repository](https://github.com/vdemichev/DiaNN)
- [DIA-NN Documentation](https://vdemichev.github.io/DiaNN/)
- [DIA-NN Issue #1068 - Column Documentation](https://github.com/vdemichev/DiaNN/issues/1068)
- [DIA-NN Report Files Reference](https://deepwiki.com/vdemichev/DiaNN/6.1-report-files)
- QPX Source Code: `qpx/core/common.py`, `qpx/core/diann/diann.py`, `qpx/core/format.py`
