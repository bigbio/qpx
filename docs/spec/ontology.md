# Field Ontology Mapping

QPX fields are mapped to controlled vocabulary (CV) and ontology terms where applicable. This mapping provides formal semantic definitions for each field, enabling interoperability with proteomics standards, biological databases, and data management platforms like lamindb.

!!! important "snake_case everywhere"
    All QPX field names and score names use **snake_case** -- no colons, dots, spaces, or mixed case. This ensures every name is a valid identifier in SQL, Python, R, and any query language without quoting. The proper ontology term name (e.g., `Comet:xcorr` for `comet_xcorr`) and accession (e.g., `MS:1002252`) are stored in this mapping.

QPX uses three categories of ontologies:

1. **PSI-MS** ([MS](https://www.ebi.ac.uk/ols4/ontologies/ms)) -- mass spectrometry terms for instrument, acquisition, and identification fields.
2. **Biological ontologies** -- for sample metadata fields: NCBI Taxonomy, UBERON, EFO/MONDO, Cell Ontology, etc.
3. **UniMod** ([UNIMOD](https://www.unimod.org/)) -- for post-translational modification definitions.

!!! info "Fields without ontology mappings"
    Fields with blank ontology entries are QPX-specific fields without standardized CV terms. These fields are defined entirely by the QPX specification.

## `ontology.parquet` -- machine-readable mapping {#ontology-parquet}

Each QPX dataset includes an `ontology.parquet` file that stores the field-to-ontology mapping in a queryable format. This makes the dataset self-describing -- any consumer can look up the proper ontology term for any field or score name without external documentation.

See the full YAML schema in [`ontology.yaml`](schemas/ontology.yaml).

!!! note "Primary key"
    Rows are unique on `(field_name, view)`. Converters collapse duplicates
    (first-wins) when the same field is accumulated both as a discovered score
    and as a mapped field, so validation does not report duplicate ontology PKs.

### Schema

```python
import pyarrow as pa

ontology_schema = pa.schema(
    [
        pa.field("field_name", pa.string()),  # snake_case QPX field name (e.g. "comet_xcorr")
        pa.field("ontology_name", pa.string(), nullable=True),  # proper ontology term (e.g. "Comet:xcorr")
        pa.field("ontology_accession", pa.string(), nullable=True),  # CV accession (e.g. "MS:1002252")
        pa.field("ontology_source", pa.string(), nullable=True),  # ontology prefix (e.g. "MS", "UBERON", "UNIMOD")
        pa.field("ontology_version", pa.string(), nullable=True),  # version of the ontology (e.g. "4.1.235")
        pa.field("view", pa.string()),  # which view (e.g. "psm", "feature", "sample")
        pa.field("description", pa.string(), nullable=True),  # human-readable description
        pa.field("source_column_name", pa.string(), nullable=True),  # original column name in tool output
        pa.field("source_tool", pa.string(), nullable=True),  # tool that produced this field
    ]
)
```

### Example rows

| field_name | ontology_name | ontology_accession | ontology_source | view | description |
|---|---|---|---|---|---|
| `comet_xcorr` | Comet:xcorr | MS:1002252 | MS | psm | Cross-correlation score from Comet |
| `organism` | organism | NCBITaxon | NCBITaxon | sample | Species |
| `organism_part` | organism part | UBERON | UBERON | sample | Anatomical structure |
| `rt` | retention time | MS:1000894 | MS | psm | Retention time (seconds) |
| `anchor_protein` | anchor protein | MS:1001591 | MS | feature | Representative protein |

### Querying

```sql
-- Look up the ontology term for a field
SELECT ontology_name, ontology_accession
FROM 'PXD014414.ontology.parquet'
WHERE field_name = 'comet_xcorr';

-- All fields with PSI-MS mappings
SELECT field_name, ontology_name, ontology_accession
FROM 'PXD014414.ontology.parquet'
WHERE ontology_source = 'MS';

-- All sample metadata ontology mappings
SELECT field_name, ontology_name, ontology_source
FROM 'PXD014414.ontology.parquet'
WHERE view = 'sample';
```

---

## Sample Metadata (`sample.parquet`) {#sample}

Sample metadata fields map to biological ontologies. These mappings are critical for integration with lamindb/bionty and for validating annotations against community-standard registries.

| Field | Ontology | Ontology Accession | Description |
|---|---|---|---|
| `sample_accession` | --- | --- | Unique sample identifier (= SDRF `source name`) |
| `organism` | NCBI Taxonomy | [NCBITaxon](https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon) | Species (e.g., `Homo sapiens` = NCBITaxon:9606) |
| `organism_part` | UBERON | [UBERON](https://www.ebi.ac.uk/ols4/ontologies/uberon) | Anatomical structure (e.g., `heart` = UBERON:0000948) |
| `disease` | EFO / MONDO | [EFO](https://www.ebi.ac.uk/ols4/ontologies/efo), [MONDO](https://www.ebi.ac.uk/ols4/ontologies/mondo) | Disease state (e.g., `breast cancer` = MONDO:0007254) |
| `cell_line` | Cellosaurus | [Cellosaurus](https://www.cellosaurus.org/) | Cell line (e.g., `HeLa` = CVCL_0030) |
| `cell_type` | Cell Ontology | [CL](https://www.ebi.ac.uk/ols4/ontologies/cl) | Cell type (e.g., `T cell` = CL:0000084) |
| `sex` | PATO | [PATO](https://www.ebi.ac.uk/ols4/ontologies/pato) | Biological sex (e.g., `female` = PATO:0000383) |
| `age` | --- | --- | Age of the specimen (free text) |
| `developmental_stage` | HsapDv / MmusDv | [HsapDv](https://www.ebi.ac.uk/ols4/ontologies/hsapdv) | Developmental stage |
| `ancestry` | HANCESTRO | [HANCESTRO](https://www.ebi.ac.uk/ols4/ontologies/hancestro) | Ancestry category |
| `individual` | --- | --- | Individual/patient identifier |
| `sample_description` | --- | --- | Free-text sample description |

!!! tip "lamindb / bionty compatibility"
    Each sample metadata field maps to a bionty registry:

    | QPX field | bionty registry | Example validation |
    |---|---|---|
    | `organism` | `bt.Organism` | `bt.Organism.from_source(name="Homo sapiens")` |
    | `organism_part` | `bt.Tissue` | `bt.Tissue.from_source(name="heart")` |
    | `disease` | `bt.Disease` | `bt.Disease.from_source(name="breast cancer")` |
    | `cell_line` | `bt.CellLine` | `bt.CellLine.from_source(name="HeLa")` |
    | `cell_type` | `bt.CellType` | `bt.CellType.from_source(name="T cell")` |

    See [Compatibility with lamindb and AnnData](#compatibility) for details.

---

## Run Metadata (`run.parquet`) {#run}

| Field | Ontology | Ontology Accession | Description |
|---|---|---|---|
| `run_accession` | --- | --- | Unique run identifier (= SDRF `assay name`) |
| `run_file_name` | --- | --- | Raw data file name without path or extension |
| `samples` | --- | --- | Sample-channel mapping (FK to `sample.parquet`) |
| `fraction` | --- | --- | Fraction identifier |
| `instrument` | PSI-MS | [MS](https://www.ebi.ac.uk/ols4/ontologies/ms) | Mass spectrometer (e.g., `Q Exactive HF` = MS:1002523) |
| `enzymes` | PSI-MS | [MS](https://www.ebi.ac.uk/ols4/ontologies/ms) | Proteolytic enzymes (e.g., `Trypsin` = MS:1001251) |
| `dissociation_method` | PSI-MS | [MS](https://www.ebi.ac.uk/ols4/ontologies/ms) | Fragmentation method (e.g., `HCD` = MS:1000422) |
| `modification_parameters` | UniMod | [UNIMOD](https://www.unimod.org/) | Search modifications (e.g., `Phospho` = UNIMOD:21) |
| `additional_terms` | PSI-MS | [MS](https://www.ebi.ac.uk/ols4/ontologies/ms) | Other ontology-backed terms (acquisition method, labeling, etc.) |

---

## Dataset Metadata (`dataset.parquet`) {#dataset}

| Field | Ontology | Ontology Accession | Description |
|---|---|---|---|
| `project_accession` | --- | --- | ProteomeXchange accession or local project identifier |
| `project_title` | --- | --- | Title of the project |
| `project_description` | --- | --- | Description of the project |
| `pubmed_id` | --- | --- | PubMed ID of the associated publication |
| `qpx_version` | --- | --- | Version of the QPX format specification |
| `software_name` | --- | --- | Software that generated the data |
| `software_version` | --- | --- | Version of the software |
| `creation_date` | --- | --- | ISO 8601 creation date |

---

## PSM View (`psm.parquet`) {#psm}

| Field | Ontology Name | Ontology Accession | Description |
|---|---|---|---|
| `sequence` | AA sequence | `MS:1001344` | Amino acid sequence |
| `peptidoform` | peptidoform | `MS:1003049` | Peptide with modifications in ProForma notation |
| `modifications` | protein modifications | `MS:1000933` | Structured modification list with UNIMOD accessions |
| `charge` | charge state | `MS:1000041` | Charge state of the precursor ion |
| `posterior_error_probability` | posterior error probability | `MS:1001493` | PEP score: probability the PSM is incorrect (lower is better) |
| `is_decoy` | decoy peptide | `MS:1002217` | Whether the PSM is from a decoy sequence |
| `calculated_mz` | theoretical monoisotopic m/z | `MS:1003053` | Theoretical m/z from molecular composition |
| `observed_mz` | selected ion m/z | `MS:1002234` | Experimentally observed m/z |
| `rt` | retention time | `MS:1000894` | Retention time of the MS2 scan (seconds) |
| `predicted_rt` | predicted retention time | `MS:1000897` | Predicted retention time (seconds) |
| `run_file_name` | --- | --- | Spectrum file name without path or extension |
| `scan` | scan number | `MS:1003057` | Scan identifier as array of integer components |
| `protein_accessions` | --- | --- | Protein accessions the peptide maps to (optional column) |
| `ion_mobility` | ion mobility drift time | `MS:1002476` | Ion mobility value for the precursor |
| `mz_array` | m/z array | `MS:1000514` | Array of fragment m/z values |
| `intensity_array` | intensity array | `MS:1000515` | Array of fragment intensity values |
| `charge_array` | --- | --- | Array of fragment ion charge values |
| `ion_type_array` | --- | --- | Array of fragment ion type annotations (b, y, a, etc.) |
| `ion_mobility_array` | --- | --- | Array of fragment ion mobility values |
| `mass_error_ppm` | --- | --- | Mass error in ppm |
| `missed_cleavages` | missed cleavages | `MS:1003044` | Number of missed enzymatic cleavages |
| `additional_scores` | --- | --- | Score array with name, value, and direction (see [Scores](scores.md)) |
| `cv_params` | --- | --- | Controlled vocabulary parameters |

---

## Feature View (`feature.parquet`) {#feature}

| Field | Ontology Name | Ontology Accession | Description |
|---|---|---|---|
| `sequence` | AA sequence | `MS:1001344` | Amino acid sequence |
| `peptidoform` | peptidoform | `MS:1003049` | Peptide with modifications in ProForma notation |
| `modifications` | protein modifications | `MS:1000933` | Structured modification list |
| `charge` | charge state | `MS:1000041` | Charge of the quantified analyte |
| `is_decoy` | decoy peptide | `MS:1002217` | Whether the peptide is from a decoy sequence |
| `calculated_mz` | theoretical monoisotopic m/z | `MS:1003053` | Theoretical m/z |
| `observed_mz` | selected ion m/z | `MS:1002234` | Experimentally observed m/z |
| `rt` | retention time | `MS:1000894` | Precursor retention time (seconds) |
| `rt_start` | --- | --- | Start of the retention time window |
| `rt_stop` | --- | --- | End of the retention time window |
| `predicted_rt` | predicted retention time | `MS:1000897` | Predicted retention time (seconds) |
| `ion_mobility` | ion mobility drift time | `MS:1002476` | Ion mobility value for the precursor |
| `ion_mobility_start` | --- | --- | Start ion mobility value |
| `ion_mobility_stop` | --- | --- | End ion mobility value |
| `run_file_name` | --- | --- | Run file containing the feature |
| `id_run_file_name` | --- | --- | Run file containing the best PSM for this feature |
| `scan` | scan number | `MS:1003057` | Scan identifier of the best PSM |
| `intensities` | --- | --- | Primary intensity across labels (see [Intensities](intensities.md)) |
| `additional_intensities` | --- | --- | Tool-provided intensities (normalized, LFQ, iBAQ) |
| `pg_accessions` | --- | --- | Protein group accessions |
| `anchor_protein` | anchor protein | `MS:1001591` | Representative protein of the protein group |
| `pg_positions` | --- | --- | Peptide start/end positions in each protein |
| `pg_global_qvalue` | protein-level global FDR | `MS:1001214` | Global q-value of the protein group |
| `unique` | --- | --- | Unique peptide indicator |
| `gg_accessions` | --- | --- | Gene group identifiers |
| `gg_names` | --- | --- | Gene group names |
| `mass_error_ppm` | --- | --- | Mass error in ppm |
| `missed_cleavages` | missed cleavages | `MS:1003044` | Number of missed enzymatic cleavages |
| `posterior_error_probability` | posterior error probability | `MS:1001493` | PEP for the peptide match |
| `additional_scores` | --- | --- | Score array (see [Scores](scores.md)) |
| `cv_params` | --- | --- | Controlled vocabulary parameters |

---

## Protein Group View (`pg.parquet`) {#pg}

| Field | Ontology Name | Ontology Accession | Description |
|---|---|---|---|
| `pg_accessions` | --- | --- | Protein accessions within this group |
| `pg_names` | --- | --- | Descriptive names for the proteins |
| `gg_accessions` | --- | --- | Gene group accessions |
| `gg_names` | --- | --- | Gene names |
| `anchor_protein` | anchor protein | `MS:1001591` | Representative protein of the group |
| `grouped_runs` | --- | --- | Raw files aggregated into this protein-group quantification unit |
| `peptide_counts` | --- | --- | Peptide sequence counts (unique + total) |
| `feature_counts` | --- | --- | Feature counts (unique + total) |
| `global_qvalue` | protein-level global FDR | `MS:1001214` | Global q-value at the experiment level |
| `pg_qvalue` | --- | --- | Run-level protein group q-value |
| `is_decoy` | decoy peptide | `MS:1002217` | Whether the protein group is a decoy |
| `contaminant` | --- | --- | Contaminant indicator (1 = contaminant) |
| `sequence_coverage` | sequence coverage | `MS:1001093` | Percentage of protein sequence covered by peptides |
| `molecular_weight` | molecular mass | `MS:1000224` | Molecular weight of the protein (kDa) |
| `intensities` | --- | --- | Primary intensity across labels |
| `additional_intensities` | --- | --- | Tool-provided intensities |
| `peptides` | --- | --- | Peptide counts per individual protein in the group |
| `additional_scores` | --- | --- | Additional scores and metrics |
| `cv_params` | --- | --- | Controlled vocabulary parameters |

---

## Mass Spectra View (`mz.parquet`) {#mz}

| Field | Ontology Name | Ontology Accession | Description |
|---|---|---|---|
| `id` | spectrum identifier nativeID format | `MS:1000767` | Unique identifier for the scan |
| `ms_level` | ms level | `MS:1000511` | MS level (1 = MS1, 2 = MS2, etc.) |
| `centroid` | centroid spectrum | `MS:1000127` | Whether data is centroided |
| `scan_start_time` | scan start time | `MS:1000016` | Start time of the scan (minutes) |
| `inverse_ion_mobility` | inverse reduced ion mobility | `MS:1002815` | Inverse ion mobility (1/K0) for TIMS data |
| `ion_injection_time` | ion injection time | `MS:1000927` | Ion injection time (milliseconds) |
| `total_ion_current` | total ion current | `MS:1000285` | Total ion current for the scan |
| `precursors` | precursor | `MS:1000456` | Precursor ions for this MS/MS scan |
| `mz` | m/z array | `MS:1000514` | Array of m/z values |
| `intensity` | intensity array | `MS:1000515` | Array of intensity values |
| `cv_params` | --- | --- | Additional CV parameters |

---

## Score Fields {#scores}

Score names used in `additional_scores` that have PSI-MS ontology mappings. All QPX score names use **snake_case** -- the ontology name column shows the proper CV term.

| QPX Name (snake_case) | Ontology Name | Ontology Accession | Direction |
|---|---|---|---|
| `posterior_error_probability` | posterior error probability | `MS:1001493` | lower is better |
| `global_qvalue` | PSM-level global FDR | `MS:1002350` | lower is better |
| `pg_global_qvalue` | protein-level global FDR | `MS:1001214` | lower is better |
| `precursor_qvalue` | DIA-NN:Q.Value | --- | lower is better |
| `andromeda_score` | Andromeda:score | `MS:1002338` | higher is better |
| `andromeda_delta_score` | Andromeda:delta score | `MS:1003433` | higher is better |
| `percolator_score` | percolator:score | `MS:1001492` | higher is better |
| `openms_score` | OpenMS:Best PSM Score | `MS:1003114` | higher is better |
| `msgf_raw_score` | MS-GF:RawScore | `MS:1002049` | higher is better |
| `comet_xcorr` | Comet:xcorr | `MS:1002252` | higher is better |
| `comet_deltacn` | Comet:deltacn | `MS:1002253` | higher is better |
| `comet_expect` | Comet:expectation value | `MS:1002257` | lower is better |
| `msgf_spec_evalue` | MS-GF:SpecEValue | `MS:1002052` | lower is better |
| `sage_hyperscore` | Sage:hyperscore | --- | higher is better |
| `diann_qvalue` | --- | --- | lower is better |
| `diann_global_qvalue` | --- | --- | lower is better |
| `diann_cscore` | --- | --- | higher is better |
| `consensus_support` | --- | --- | higher is better |

---

## QPX-specific Fields {#qpx-specific}

Fields defined by QPX without standardized ontology terms.

| Field | View(s) | Description |
|---|---|---|
| `run_file_name` | PSM, Feature | Spectrum/run file name without path or extension |
| `grouped_runs` | PG | Raw files aggregated into a protein-group quantification unit |
| `pg_accessions` | Feature, PG | Protein group accessions |
| `pg_names` | PG | Descriptive protein names |
| `gg_accessions` | PG | Gene group accessions |
| `gg_names` | PG | Gene names |
| `rt_start` / `rt_stop` | Feature | Retention time window boundaries |
| `ion_mobility_start` / `ion_mobility_stop` | Feature | Ion mobility window boundaries |
| `intensities` / `additional_intensities` | Feature, PG | Intensity structures (see [Intensities](intensities.md)) |
| `peptide_counts` / `feature_counts` | PG | Count structures for peptides and features |
| `contaminant` | PG | Contaminant indicator |
| `pg_positions` | Feature | Peptide positions within each protein |
| `id_run_file_name` / `id_scan` | Feature | Reference to best PSM for the feature |
| `protein_accessions` | PSM | Protein accessions (optional column) |
| `charge_array` / `ion_type_array` / `ion_mobility_array` | PSM | Fragment ion annotation arrays |
| `parent_ion_fraction` | Score | Fraction of target peak intensity in the isolation window |
| `precursor_quantification_score` | Score | QuantUMS-derived quantification confidence: 1.0 / (1.0 + SD) |

---

## Compatibility with lamindb and AnnData {#compatibility}

### AnnData

The QPX expression views (`ae.h5ad`, `de.h5ad`) use AnnData natively. The ontology mapping provides semantic definitions for the `.obs` (sample) and `.var` (protein/gene) slots:

| AnnData slot | QPX source | Ontology-backed fields |
|---|---|---|
| `.obs` | `sample.parquet` | `organism` (NCBITaxon), `organism_part` (UBERON), `disease` (EFO/MONDO), `cell_type` (CL), `cell_line` (Cellosaurus), `sex` (PATO) |
| `.var` | `pg.parquet` | `protein` (UniProt), `gene_name` (HGNC/Ensembl) |
| `.uns` | `dataset.parquet` | `project_accession`, `qpx_version` |

When building AnnData objects from QPX Parquet files, the `.obs` columns inherit the ontology definitions from `sample.parquet`. This means any AnnData created from QPX data carries structured, ontology-backed annotations that tools like scanpy and lamindb can validate.

### lamindb / bionty

[lamindb](https://lamin.ai/) uses [bionty](https://lamin.ai/docs/bionty) registries to validate biological entities. QPX's sample metadata fields map directly to bionty registries:

| QPX field | bionty registry | Ontology source | Example |
|---|---|---|---|
| `organism` | `bt.Organism` | NCBI Taxonomy | `Homo sapiens` → NCBITaxon:9606 |
| `organism_part` | `bt.Tissue` | UBERON | `heart` → UBERON:0000948 |
| `disease` | `bt.Disease` | MONDO | `breast cancer` → MONDO:0007254 |
| `cell_type` | `bt.CellType` | Cell Ontology | `T cell` → CL:0000084 |
| `cell_line` | `bt.CellLine` | Cellosaurus | `HeLa` → CVCL_0030 |
| `sex` | --- | PATO | `female` → PATO:0000383 |
| `developmental_stage` | `bt.DevelopmentalStage` | HsapDv | --- |

This means QPX datasets can be registered in lamindb with validated, ontology-backed annotations:

```python
import lamindb as ln
import bionty as bt
import anndata as ad

# Read QPX AnnData
adata = ad.read_h5ad("PXD014414.ae.h5ad")

# Validate against bionty registries
bt.Organism.from_source(name="Homo sapiens")
bt.Tissue.from_source(name="heart")
bt.Disease.from_source(name="breast cancer")

# Register as a lamindb artifact with validated metadata
artifact = ln.Artifact.from_anndata(adata, description="PXD014414 absolute expression")
artifact.save()
```

### Parquet views and lamindb

The Parquet data views (`psm`, `feature`, `pg`, `mz`) can also be registered as lamindb artifacts. The ontology mapping allows lamindb to understand what each field represents:

```python
import lamindb as ln

# Register Parquet files as artifacts
artifact = ln.Artifact("PXD014414.feature.parquet", description="PXD014414 features")
artifact.save()

# Link to biological context via validated sample metadata
artifact.organisms.add(bt.Organism.from_source(name="Homo sapiens"))
artifact.tissues.add(bt.Tissue.from_source(name="heart"))
```

The key benefit: because QPX defines clear ontology mappings for sample metadata, lamindb can automatically validate and link QPX datasets to the correct biological entities in its registry. No manual curation of field semantics is needed.

---

## Ontology references

| Ontology | Prefix | URL | Used for |
|---|---|---|---|
| PSI-MS | `MS:` | [OLS4](https://www.ebi.ac.uk/ols4/ontologies/ms) | MS instrument, acquisition, identification fields |
| UniMod | `UNIMOD:` | [UniMod](https://www.unimod.org/) | Post-translational modifications |
| NCBI Taxonomy | `NCBITaxon:` | [OLS4](https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon) | Organism / species |
| UBERON | `UBERON:` | [OLS4](https://www.ebi.ac.uk/ols4/ontologies/uberon) | Anatomical structures / tissues |
| EFO | `EFO:` | [OLS4](https://www.ebi.ac.uk/ols4/ontologies/efo) | Experimental factors, diseases |
| MONDO | `MONDO:` | [OLS4](https://www.ebi.ac.uk/ols4/ontologies/mondo) | Disease ontology |
| Cell Ontology | `CL:` | [OLS4](https://www.ebi.ac.uk/ols4/ontologies/cl) | Cell types |
| Cellosaurus | `CVCL_` | [Cellosaurus](https://www.cellosaurus.org/) | Cell lines |
| PATO | `PATO:` | [OLS4](https://www.ebi.ac.uk/ols4/ontologies/pato) | Phenotypic qualities (sex, etc.) |
| HANCESTRO | `HANCESTRO:` | [OLS4](https://www.ebi.ac.uk/ols4/ontologies/hancestro) | Human ancestry |
| HsapDv | `HsapDv:` | [OLS4](https://www.ebi.ac.uk/ols4/ontologies/hsapdv) | Human developmental stages |
