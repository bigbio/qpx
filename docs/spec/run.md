# Run Metadata

`run.parquet` stores one row per data acquisition run, capturing the technical and instrument properties of each MS run. Each row links to one or more biological samples through a `samples` list, which supports TMT/iTRAQ multiplexing where a single run contains multiple labeled channels. The file is derived from the community [SDRF](https://github.com/bigbio/proteomics-sample-metadata) (Sample and Data Relationship Format) standard.

See the full YAML schema in [`run.yaml`](schemas/run.yaml).

---

## Instrument, Enzyme, and Dissociation Fields

Instrument, enzyme, and dissociation method fields store **plain strings** (human-readable names). CV accession mappings are stored separately in `ontology.parquet`, keeping `run.parquet` simple and consistent with `sample.parquet`.

---

## Per-Run Search and Channel Fields {#per-run-fields}

`modification_parameters` lists the modifications of the SDRF rows of that run, so that a phospho run and a global run of one study keep their own. These fields are in specification 1.3:

| Field | Source | Meaning |
| ----- | ------ | ------- |
| `label_definitions` | `comment[modification parameters]` of each label's rows | For each SILAC, dimethyl or mTRAQ label of the run, its label modifications with their masses. `SILAC heavy` alone does not say whether a plex used K4R6 or K8R10 |
| `gradient_length` | `comment[gradient length]` | LC gradient in seconds. `120 min`, `2 h` and `7200 s` are read; a bare number is minutes |
| `acquisition_start` | `comment[acquisition date]`, else mzML `run/@startTimeStamp` | Run acquisition start, preserved as source text |
| `samples[].channel_role` | `comment[channel role]`, or `comment[reference channel]` | `reference` for the channel `comment[reference channel]` names, `sample` for the others; otherwise the value of `comment[channel role]`, such as `carrier` or `empty` |
| `samples[].mixture` | `comment[mixture]` | Identifier of the labelled mixture, to separate TMT plexes |

`acquisition_start` is an optional, nullable string. The first valid `comment[acquisition date]` value of each run takes precedence; SDRF missing values such as `not available` and `not applicable` are ignored. If no valid value exists, `qpxc convert sdrf --mzml-dir` reads the `startTimeStamp` attribute from the run's `.mzML` or `.mzML.gz` header. The mzML source is currently exposed only by this command. Missing files or attributes leave the value null. Text is preserved verbatim, including date-only values (`2011-05-31`), fractional seconds (`2016-10-06T22:36:56.803000`) and timezone suffixes (`2015-02-19T20:08:48Z`); QPX does not parse dates, convert timezones or fill in missing time components.

A sample may occupy several labels of a run, as in pulsed SILAC; helpers keep them apart by `(sample, label, unit)`.

---

## MODIFICATION Structure {#modification-structure}

Modification parameters use a dedicated struct with explicit fields, since modifications have well-known properties (fixed/variable, position, target residue) that benefit from a typed schema:

```python
MODIFICATION = pa.struct(
    [
        pa.field("accession", pa.string()),  # UNIMOD accession (required), e.g. "UNIMOD:21"
        pa.field("name", pa.string(), nullable=True),  # human-readable name, e.g. "Phospho"
        pa.field("fixed", pa.bool_(), nullable=True),  # true = fixed, false = variable
        pa.field("position", pa.string(), nullable=True),  # e.g. "Anywhere", "Protein N-term"
        pa.field("target_amino_acid", pa.string(), nullable=True),  # e.g. "S,T,Y", "C"
        pa.field("mass_delta", pa.float64(), nullable=True),  # SDRF MM=, else the Unimod mass
    ]
)

# Generic ontology term (accession + optional name/properties) used by additional_terms
ONTOLOGY_TERM = pa.struct(
    [
        pa.field("accession", pa.string()),  # CV accession, e.g. "MS:1003221" (DIA)
        pa.field("name", pa.string(), nullable=True),  # human-readable term name
        pa.field(
            "properties",
            pa.list_(
                pa.struct(
                    [
                        pa.field("key", pa.string()),
                        pa.field("value", pa.string()),
                    ]
                )
            ),
            nullable=True,
        ),
    ]
)
```

---

## PyArrow Schema

```python
import pyarrow as pa

# (ONTOLOGY_TERM and MODIFICATION defined above)

run_schema = pa.schema(
    [
        # Identity
        pa.field("run_accession", pa.string()),  # run file stem
        pa.field("run_file_name", pa.string()),  # PK: raw data file name (without extension)
        pa.field("file_name", pa.string(), nullable=True),  # original file name with extension (e.g. "S1_Frontal_1.raw")
        # Sample-channel mapping (one run can contain multiple samples via TMT/iTRAQ)
        pa.field(
            "samples",
            pa.list_(
                pa.struct(
                    [
                        pa.field("sample_accession", pa.string()),  # FK to sample.parquet
                        pa.field("label", pa.string()),  # channel label, e.g. "TMT126", "LFQ"
                        pa.field("biological_replicate", pa.int32(), nullable=True),
                        pa.field("technical_replicate", pa.int32(), nullable=True),
                    ]
                )
            ),
        ),
        # Run-level properties
        pa.field("fraction", pa.string(), nullable=True),
        # Instrument and method (plain strings; CV mappings in ontology.parquet)
        pa.field("instrument", pa.string(), nullable=True),
        pa.field("enzymes", pa.list_(pa.string()), nullable=True),
        pa.field("dissociation_method", pa.string(), nullable=True),
        pa.field("acquisition_method", pa.string(), nullable=True),
        pa.field("acquisition_start", pa.string(), nullable=True),
        # Other ontology-backed run terms (acquisition method as CV, labeling, etc.)
        pa.field("additional_terms", pa.list_(ONTOLOGY_TERM), nullable=True),
        # Modification parameters -- list of typed modification structs
        pa.field("modification_parameters", pa.list_(MODIFICATION), nullable=True),
    ]
)
```

## Field Reference

*Required* follows [Required, nullable and optional columns](data-model.md#required-nullable-and-optional-columns): an *optional* column may be left out of the file.

| Field | Description | Type | Required |
|---|---|---|---|
| `run_accession` | Run identifier. Converters write the run file stem, whatever the SDRF `assay name` says, so it equals `run_file_name`; not the primary key | string | yes |
| `run_file_name` | Raw data file name (without extension); the primary key | string | yes |
| `file_name` | Original file name with extension (e.g. `S1_Frontal_1.raw`) | string | optional |
| `samples` | List of sample-channel mappings with replicate info (FK to sample.parquet) | list[sample_channel] | yes |
| `fraction` | Fraction identifier | string | optional |
| `instrument` | Mass spectrometer name | string | optional |
| `enzymes` | Proteolytic enzyme name(s) | list[string] | optional |
| `dissociation_method` | Fragmentation method name (e.g. HCD, CID, ETD) | string | optional |
| `acquisition_method` | Data acquisition method name (e.g. DDA, DIA); from `comment[proteomics data acquisition method]` | string | optional |
| `additional_terms` | Other ontology-backed run terms (acquisition method, labeling, etc.) | list[ontology_term] | optional |
| `modification_parameters` | Modifications configured in the database search | list[modification_param] | optional |
| `label_definitions` | For each SILAC, dimethyl or mTRAQ label of the run, its label modifications; see above | list[label_definition] | optional |
| `gradient_length` | LC gradient in seconds; see above | float32 | optional |
| `acquisition_start` | Run acquisition start as source text; see above | string | optional |

### Samples list detail

Each element in the `samples` list contains:

| Subfield | Description | Type |
|---|---|---|
| `sample_accession` | FK to `sample.parquet` | string |
| `label` | Channel label (e.g. `TMT126`, `TMT127N`, `LFQ`) | string |
| `biological_replicate` | Biological replicate number | int32 (nullable) |
| `technical_replicate` | Technical replicate number | int32 (nullable) |

!!! note "Why replicates live in the samples struct"
    `biological_replicate` and `technical_replicate` are stored per sample-run pair rather than at the run or sample level. This makes the `samples` list a self-contained design matrix: each entry fully describes one sample's role in one run, without requiring joins to additional tables.

!!! example "Label-free vs TMT"
    In a **label-free** experiment, each run maps to one sample with label `"LFQ"`:
    ```json
    {"samples": [
        {"sample_accession": "Sample_01", "label": "LFQ",
         "biological_replicate": 1, "technical_replicate": 1}
    ]}
    ```

    In a **TMT-10plex** experiment, each run maps to 10 samples:
    ```json
    {"samples": [
        {"sample_accession": "Sample_01", "label": "TMT126",
         "biological_replicate": 1, "technical_replicate": 1},
        {"sample_accession": "Sample_02", "label": "TMT127N",
         "biological_replicate": 2, "technical_replicate": 1},
        {"sample_accession": "Sample_03", "label": "TMT127C",
         "biological_replicate": 3, "technical_replicate": 1},
        ...
    ]}
    ```

### MODIFICATION detail

Each MODIFICATION struct contains:

| Subfield | Description | Type | Required |
|---|---|---|---|
| `accession` | UNIMOD accession (e.g. `UNIMOD:21`, `UNIMOD:4`) | string | yes |
| `name` | Human-readable name (e.g. `Phospho`, `Carbamidomethyl`) | string | no |
| `fixed` | `true` = fixed modification, `false` = variable modification | bool | no |
| `position` | Where the modification can occur (e.g. `Anywhere`, `Protein N-term`) | string | no |
| `target_amino_acid` | Target residue(s), comma-separated (e.g. `S,T,Y`, `C`) | string | no |

## Example Data

**Modification parameters:**

```json
[
    {
        "accession": "UNIMOD:4",
        "name": "Carbamidomethyl",
        "fixed": true,
        "position": "Anywhere",
        "target_amino_acid": "C"
    },
    {
        "accession": "UNIMOD:21",
        "name": "Phospho",
        "fixed": false,
        "position": "Anywhere",
        "target_amino_acid": "S,T,Y"
    }
]
```

**Enzymes:**

```json
["Trypsin"]
```

**Instrument:**

```json
"Q Exactive HF"
```

**Dissociation method:**

```json
"HCD"
```

!!! note "CV accession resolution"
    CV accessions for instrument, enzyme, and dissociation method names are stored in `ontology.parquet`, not in `run.parquet`. This keeps `run.parquet` simple and consistent with `sample.parquet`.

!!! tip "Relationship to per-PSM modifications"
    The `modification_parameters` list in `run.parquet` captures what was configured in the search engine. For per-peptide modification evidence with localization scores, see the [Modifications](modifications.md) specification.

!!! tip "Reading run.parquet"
    ```python
    import pyarrow.compute as pc
    import pyarrow.parquet as pq

    runs = pq.read_table("PXD014414.run.parquet")
    print(runs.column_names)
    # ['run_accession', 'run_file_name', 'samples', 'fraction',
    #  'instrument', 'enzymes', ...]

    # Get all sample-channel mappings for a specific run
    run_row = runs.filter(pc.equal(runs.column("run_accession"), "Run_01")).to_pydict()
    print(run_row["samples"])
    ```

!!! tip "Querying with DuckDB"
    ```sql
    -- All runs with their instruments
    SELECT run_accession, instrument
    FROM 'PXD014414.run.parquet';

    -- Unnest samples to get the full design matrix
    SELECT r.run_accession, s.sample_accession, s.label,
           s.biological_replicate, s.technical_replicate
    FROM 'PXD014414.run.parquet' r, UNNEST(r.samples) AS s;

    -- Join runs with samples to get full metadata
    SELECT r.run_accession, r.run_file_name,
           rs.biological_replicate, rs.technical_replicate,
           s.organism, s.disease
    FROM 'PXD014414.run.parquet' r,
         UNNEST(r.samples) AS rs,
         'PXD014414.sample.parquet' s
    WHERE rs.sample_accession = s.sample_accession;
    ```

---

## Related Pages

- [Sample Metadata](sample.md) -- `sample.parquet`
- [Dataset Metadata](dataset.md) -- `dataset.parquet`
- [SDRF Conversion](sdrf-mapping.md) -- Column mapping and ingestion process
- [Modifications](modifications.md) -- Per-peptide modification evidence
- [QPX Format Overview](index.md)
