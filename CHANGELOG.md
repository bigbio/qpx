# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Spectronaut converter**: `qpxc convert spectronaut` — full support for Spectronaut report TSV files, producing feature.parquet and pg.parquet with DuckDB-accelerated batch processing
- **CPTAC CDAP converter**: `qpxc convert cdap` — convert CPTAC CDAP `.psm` study directories to QPX psm/feature/pg/dataset/ontology/provenance views
- **Full-spectra mz converter**: `qpxc convert mz` — convert a directory of mzML / `.mzML.gz` files to a single `mz.parquet`; each spectrum carries `run_file_name` + `scan` for linking back to PSM/feature
- **pdc2qpx pipeline**: `qpxc pdc2qpx` — one-shot PDC/CPTAC download (via pridepy, `qpx[pdc]` extra) + CDAP + full-spectra conversion into an entire QPX dataset
- **Shared channel-label resolution**: `qpx/converters/channel_labels.py` — single source for canonical TMT/iTRAQ/LFQ labels via sdrf-pipelines `channel_map`, used by the OpenMS `consensusXML` and `-out_qpx` paths
- **openms-consensus interim protein intensity**: the `openms-consensus` converter now fills `pg.intensity` with an interim, **unnormalized sum of each group's unique peptides** per `(protein group, grouped_runs, label)` (the quantms `unique_peptides` policy) instead of leaving it null, until OpenMS `-out_qpx` ships the authoritative quant. Every quantified row is stamped with a `quantification_method` cv_param; `--pg-top N` bounds the peptides used (`0` = all; `3` mirrors the ProteomicsLFQ/IsobaricWorkflow default)
- **openms-consensus channel/SDRF consistency check**: when an SDRF is provided, the converter compares the isobaric channels read from the consensusXML maps against the SDRF `comment[label]` set and logs a warning for any channel present in one but not the other (e.g. a mis-declared plex or the wrong SDRF)
- **Mandatory identity ids**: the `feature`, `psm` and `pg` views each carry a single required `int64` identity column — `feature_id` / `psm_id` / `pg_id` — that is the primary key of the view. The id is derived by the writer as an opaque hash of a footer-declared `identity_composite` of existing columns (feature `[peptidoform, charge, run_file_name, rt]`, psm `[peptidoform, charge, run_file_name, scan]`, pg `[anchor_protein, grouped_runs, label]`), or accepted under the view's producer-ID rule. Each file self-describes its `primary_key` + `identity_composite` in the footer, and primary-key uniqueness validation catches any hash collision. See bigbio/qpx#229.
- **Cross-view reference columns (optional)**: `feature.psm_ids` (`list<int64>`), `feature.pg_ids` (`list<int64>`) and `psm.feature_id` (`int64`, nullable) let the views reference each other by id (feature ↔ PSM ↔ protein group); populated where a converter can resolve the mapping, null otherwise
- **Writer identity enforcement**: `FeatureWriter` re-derives identified Feature ids by default, stashing an overridden producer id as a `provided_feature_id` cv_param. Unidentified Features without a producer id may derive one from the declared composite only when the resulting full-file primary key is unique; successful fallback is reported as a warning. Producer ids are namespaced by run. Existing-file transformations can preserve the stored id and its footer-declared composite. Identity-bearing writers validate the complete output file on close and raise on null or duplicate primary keys.
- **MuData (`.h5mu`) export**: OpenMS (`-out_qpx`), `openms-consensus`, and DIA-NN conversions now also emit a `<prefix>.h5mu` MuData container (`qpx/mudata.py` + `orchestrator._write_mudata`), giving parity with quantmsdiann. It bundles up to four modalities — `precursors` (peptidoform×charge feature intensities), `proteins` (protein-group intensities), `expression` (Absolute Expression), `differential` (Differential Expression) — plus a `varp["feature_mapping"]` precursor↔protein adjacency. See the MuData section of the serialization spec
- **De-novo / no-database-search support**: features and PSMs from de-novo or search-free workflows (no protein mapping, no target-decoy) are now accepted

### Changed

- **Feature/PSM/PG primary key is now the derived id** (`feature_id` / `psm_id` / `pg_id`) instead of the composite tuple; the previous composite is retained as the footer-declared `identity_composite` the id is derived from
- **Deterministic ids across write paths**: the id is now derived from the **persisted (Arrow-cast) values** in a single place, so `write_batch`, `write_table`, and `write_dataframe` produce identical ids for the same row (previously a float32 field such as `rt` could hash differently on the record vs table path). The `write_table`/`write_dataframe` paths derive the id before the non-nullable cast, so a frame/table without a precomputed id is accepted
- **Injective identity encoding**: `canonical()` uses JSON (quoted/escaped) instead of a delimiter join, so composites whose list elements or fields contain the delimiters (e.g. run names with commas) can no longer alias to the same id. Ordered list fields such as multi-component `scan` identifiers retain their order; only set-valued `grouped_runs` is canonicalized without regard to order
- **Referential validation** also flags reciprocal feature↔psm desync (an edge present one way but contradicted by the other), as a warning, only where the opposite direction is populated
- **OpenMS consensus identity** retains every spectrum reference associated with a run and the parent `consensus_rt`, preventing distinct ConsensusFeatures from collapsing onto the same derived `feature_id`
- **BREAKING — pg view flattened**: the `pg` view changes from one row per `(anchor_protein, grouped_runs)` carrying a nested `intensities: list<{label, intensity}>` column to **one row per label** with scalar `label` + `intensity` columns (identification-only groups become a single row with null `label`/`intensity`). Converters still build the natural per-group `intensities` list; the writer materializes the flat layout. Regenerate pg files
- **De-novo relaxation**: `feature.is_decoy`, `feature.anchor_protein`, and `psm.is_decoy` are now **nullable** (previously required) so de-novo / no-inference records are representable — consumers must no longer assume these are non-null
- **BREAKING — Feature/PSM identity composites** (format 1.1, issue #217): the
  schema-default Feature composite changes from `[sequence, charge,
  run_file_name, anchor_protein]` to `[peptidoform, charge, run_file_name, rt]`;
  the PSM composite changes from `[sequence, charge, run_file_name, scan]` to
  `[peptidoform, charge, run_file_name, scan]`. Measured across ~13M real rows,
  `anchor_protein` is functionally redundant and apex `rt` distinguishes
  repeated peaks where the producer reports it. A converter may declare a
  different producer-specific Feature composite when its upstream entity uses
  other measured properties. The actual primary keys are the mandatory opaque
  ID columns described above. Regenerate Feature and PSM files.
- **Parquet output size**: writers now apply `BYTE_STREAM_SPLIT` encoding to high-entropy float columns (rt, rt_start/stop, predicted_rt, calculated/observed m/z, intensity arrays) and raise the ZSTD level to 9. Encoding-only and fully lossless — no schema change; output reads unchanged with pyarrow and DuckDB. Measured ~16% smaller on a 14 GB feature.parquet.

### Removed

- **BREAKING — QuantMS/mzTab converter**: removed `qpxc convert quantms`, `QuantMSConverter`, its adapters, and the mzTab+MSstats loader. Use `qpxc convert openms-consensus` (consensusXML) — see the deprecation of `qpxc convert openms` below.

### Deprecated

- **`qpxc convert openms` (the OpenMS `-out_qpx` parquet path)**: deprecated in favour of `qpxc convert openms-consensus`. OpenMS `-out_qpx` mis-assigns **every** PSM's `run_file_name` to the first run (the per-PSM `map_index` is dropped; OpenMS#9872) and emits duplicate PSMs (one feature-attached + one unassigned; OpenMS#9871) — so on any multi-run / multi-fraction / multi-replicate design ~half the PSMs carry the wrong run. `qpxc convert openms-consensus` reads the consensusXML directly and resolves the correct run per PSM (via `map_index`/`id_merge_index`/feature-element runs), so it does not have these defects. The `openms` command now emits a deprecation warning; it will be reconsidered once OpenMS ships an `-out_qpx` that carries the correct per-PSM run.

### Fixed

- **DIA-NN q-value handling — feature/pg consistency, filtering is opt-in** (bigbio/qpx#241): the DIA-NN feature and pg views are now consistent and **neither filters by default**. `--qvalue-threshold` is opt-in (default unset): with no threshold the converter emits every feature and protein group the report contains, since DIA-NN already FDR-filters its main report and all per-row q-value columns (feature: `Q.Value`, `Global.Q.Value`, `PEP`, `Lib.Q.Value`, `Translated.Q.Value`; pg: `PG.Q.Value`, `Global.PG.Q.Value`, `Protein.Q.Value`, `Lib.PG.Q.Value`, `GG.Q.Value`) are carried through unconditionally for downstream filtering at any level. This resolves the earlier feature/pg inconsistency and count inflation by making both views as-reported rather than by forcing a pg filter. When `--qvalue-threshold` **is** given, each view filters at its own level — the feature view on precursor `Q.Value` and the pg view on `Global.PG.Q.Value` (falling back to `PG.Q.Value`). The Spectronaut converter's `--qvalue-threshold` is likewise opt-in.
- **Deterministic DIA-NN merged-group PG quantity**: when annotation noise merges protein-group rows that disagree on `PG.Quantity` / `PG.MaxLFQ`, the pg view now emits the max of the finite values (order-independent) instead of a row-order-dependent "first", and warns when the values actually disagree.
- **DIA-NN `observed_mz` NULL instead of 0.0 sentinel** (bigbio/qpx#244): a feature with no `Precursor.Mz` and no `ms_info` source now carries a truthful `NULL` observed m/z rather than a `0.0` that corrupts downstream mass-error / PPM math; report runs with no matching `*_ms_info.parquet` are warned about before being dropped.
- **Controlled unidentified Feature identity fallback**: an unidentified Feature (null `peptidoform`) with no producer id derives a candidate `feature_id` from the declared composite. The fallback is accepted with a warning only after full-file primary-key uniqueness succeeds; a collision remains a hard error. A producer-supplied id remains namespaced by run and preserved as a `provided_feature_id` cv_param.
- **openms-consensus isobaric detection**: real quantms `IsobaricWorkflow` output stamps TMT/iTRAQ consensusXML with `experiment_type="label-free"` while the maps still carry `tmt6plex_*`/`itraq*plex_*` labels. The converter now detects channels from the **map label** (not `experiment_type`), so real quantms TMT no longer collapses all reporter channels into a single `LFQ` label. Verified on real cluster output (PXD000001 TMT → TMT126–131; BSA/PXD002395 LFQ unchanged)
- **RT unit conversion**: DIA-NN and MaxQuant converters now correctly convert retention time from minutes to seconds in feature and PSM parquet output
- **Code quality**: Spectronaut converter refactored to reduce cyclomatic complexity, fix logging f-string interpolation, remove unused arguments, and eliminate duplicate code
- **CDAP label-free intensity label**: label-free `PrecursorArea` intensities are now emitted with the `"LFQ"` label (aligned with the FragPipe/MaxQuant converters) so downstream label-free consumers (mokume) recognize them as primary intensities
- **QPX TMT/iTRAQ channel labels**: OpenMS `-out_qpx` enrichment relabels feature/pg intensities from filenames/bare indices to plex-aware canonical reporter names; `run.samples[].label` is normalized to match — all from the shared sdrf-pipelines `channel_map` vocabulary

## [1.0.0] - 2025-03-22

### Added

- **CLI command groups**: `convert`, `transform`, `query`, `info`, `validate`, `ontology`
- **Converters**: DIA-NN, MaxQuant, QuantMS (mzTab LFQ & TMT), FragPipe, mzIdentML, SDRF
- **Transform commands**:
  - `gene-map` — map gene names from FASTA to QPX parquet data
  - `quantify` — protein-level quantification via mokume (DirectLFQ, MaxLFQ, iBAQ, TopN, Sum)
  - `normalize-accessions` — normalize protein accession formats (full ↔ bare UniProt)
  - `update-metadata` — update sample/run metadata from a revised SDRF
- **Query commands**: `sql`, `filter`, `head` for interactive dataset exploration
- **Info commands**: dataset summary, Arrow schema display, Parquet metadata inspection
- **Validate command**: schema validation with column presence, type matching, null checks, and PK uniqueness
- **Ontology management**: `info`, `update`, `build`, `search` for PSI-MS and PRIDE CV terms
- **Python API**: `qpx.open_dataset()`, `qpx.read_feature()`, `qpx.read_psm()`, `qpx.read_pg()`, etc.
- **Dataset class**: unified access to all QPX structures with DuckDB-backed SQL queries
- **DatasetCollection**: work with multiple datasets simultaneously
- **Writers**: Parquet writers for all QPX structures (Feature, PSM, PG, Sample, Run, MzSpectra, PepMap, Ontology, Provenance)
- **Views**: analytical views for protein, peptide, QC, and sample summaries
- **Schema validation**: YAML-defined canonical schemas for all data structures
- **Score normalization**: automatic PSI-MS ontology mapping for search engine scores
- **PRIDE API integration**: `--enrich-pride` flag for automatic project metadata enrichment
- **S3 support**: read QPX datasets from S3-compatible storage
- **MkDocs documentation**: full CLI reference with auto-generated parameter tables and usage examples

### Changed

- Replaced per-converter `_safe_float_val` with shared `safe_float` utility
- Optimized `_table_exists` to use parameterized `information_schema` query
- Increased DIA-NN `file_num` default from 50 to 100 for larger cohort support
- Defensive handling of None/NaN in gene mapping (`_resolve_gene_names`, `annotate_dataframe`)

### Fixed

- `.gitignore` duplicate entry removal
- SQL injection safety in `_table_exists` (parameterized query)
- Gene mapping crash on None/NaN protein accessions
- Gene mapping `list(dict)` producing keys instead of preserving struct

## [0.0.4] - 2025-03-01

### Added

- Initial PyPI release
- Basic converter framework for QuantMS, DIA-NN, MaxQuant
- QPX Parquet schema definition
- CLI entry point (`qpxc`)
