# PSM Schema Documentation and Examples

## Summary

This PR provides refined PSM (Peptide Spectrum Match) schema documentation and examples for multiple configurations, enabling team discussion and feedback before release.

## What's Included

### 1. Schema Review

- Verified consistency between `docs/psm.avsc` (Avro schema) and `qpx/core/format.py` (PyArrow schema)
- All 22 PSM fields have correct types, nullable settings, and descriptions
- Note: Minor field order difference in `protein_accessions` position (does not affect functionality)

### 2. Documentation

The PSM section in `docs/format-specification.md` (L700-851) is complete with:

- Complete field description tables
- Tool mapping (MaxQuant, DIA-NN, FragPipe, mzTab)
- Use case descriptions (DDA/DIA, with/without spectral data)
- File metadata documentation

### 3. Example Generation (Demonstration)

The following MaxQuant PSM examples can be generated to demonstrate the conversion workflow:

| Configuration | Description           | Command                                                                             |
| ------------- | --------------------- | ----------------------------------------------------------------------------------- |
| Basic PSM     | Without spectral data | `qpxc convert maxquant-psm --msms-file msms.txt --output-folder ./`                 |
| PSM + Spectra | With spectral arrays  | `qpxc convert maxquant-psm --msms-file msms.txt --output-folder ./ --spectral-data` |

See "How to Generate Examples" section below for complete commands.

### 4. MaxQuant to PSM Conversion Verification

Verified the conversion from MaxQuant `msms.txt` to PSM parquet format:

- **Field Mapping**: All 22 fields correctly mapped via `MAXQUANT_PSM_MAP`
- **Modifications**: Correctly parsed with name, accession, positions structure
- **Additional Scores**: Collected andromeda_score, andromeda_delta_score, parent_ion_fraction
- **Spectral Arrays**: mz_array, intensity_array, charge_array, ion_type_array all correctly processed

## PSM Schema Fields (22 total)

### Core Identification Fields (15)

| Field                       | Type         | Description                                    |
| --------------------------- | ------------ | ---------------------------------------------- |
| sequence                    | string       | Unmodified peptide sequence                    |
| peptidoform                 | string       | Peptide sequence with modifications (ProForma) |
| modifications               | list[struct] | Modification details with positions and scores |
| precursor_charge            | int32        | Precursor ion charge                           |
| posterior_error_probability | float32      | PEP value                                      |
| is_decoy                    | int32        | Decoy indicator (1=decoy, 0=target)            |
| calculated_mz               | float32      | Theoretical m/z                                |
| observed_mz                 | float32      | Experimental m/z                               |
| additional_scores           | list[struct] | Search engine scores                           |
| predicted_rt                | float32      | Predicted retention time (seconds)             |
| reference_file_name         | string       | Reference file name                            |
| cv_params                   | list[struct] | CV parameters                                  |
| scan                        | string       | Scan number                                    |
| rt                          | float32      | MS2 retention time (seconds)                   |
| ion_mobility                | float32      | Ion mobility value                             |

### PSM Unique Fields (7)

| Field              | Type          | Description                          |
| ------------------ | ------------- | ------------------------------------ |
| protein_accessions | list[string]  | Protein accessions                   |
| number_peaks       | int32         | Number of peaks                      |
| mz_array           | list[float32] | m/z values array                     |
| intensity_array    | list[float32] | Intensity values array               |
| charge_array       | list[int32]   | Fragment ion charge array            |
| ion_type_array     | list[string]  | Ion type annotations (b, y, a, etc.) |
| ion_mobility_array | list[float32] | Fragment ion mobility array          |

## How to Generate Examples

```bash
# Basic PSM (without spectral data)
qpxc convert maxquant-psm \
    --msms-file tests/examples/maxquant/maxquant_simple/msms.txt \
    --output-folder examples/psm/ \
    --output-prefix maxquant-dda \
    --verbose

# PSM with spectral data
qpxc convert maxquant-psm \
    --msms-file tests/examples/maxquant/maxquant_simple/msms.txt \
    --output-folder examples/psm/ \
    --output-prefix maxquant-dda-spectra \
    --spectral-data \
    --verbose
```

## How to Read Examples

```python
import pyarrow.parquet as pq

# Read PSM file
table = pq.read_table("examples/psm/maxquant-dda-*.psm.parquet")
df = table.to_pandas()

# View sample data
print(df[['sequence', 'peptidoform', 'precursor_charge', 'scan']].head())

# View modifications
print(df['modifications'].head())

# View additional scores
print(df['additional_scores'].head())
```

## Related Files

```
docs/
|-- psm.avsc                    # Avro schema definition
+-- format-specification.md     # Format specification (PSM: L700-851)

qpx/core/
|-- format.py                   # PyArrow schema (PEPTIDE_FIELDS, PSM_UNIQUE_FIELDS)
|-- common.py                   # PSM_SCHEMA, MAXQUANT_PSM_MAP
+-- maxquant/maxquant.py        # MaxQuant converter (process_psm_file)

tests/examples/maxquant/
|-- maxquant_simple/msms.txt    # Test data for PSM conversion
+-- maxquant_full/              # Full example with spectral data
```

## Feedback Welcome

This PR is intended to facilitate discussion on the PSM schema. Please review:

1. Are the field definitions clear and complete?
2. Do the examples adequately demonstrate different use cases?
3. Is the MaxQuant field mapping appropriate?
4. Any suggestions for additional configurations or examples?
