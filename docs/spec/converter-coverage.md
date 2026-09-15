# Converter Coverage Matrix

This page lists which QPX data views each converter produces and explains the sources and null values of selected fields.

## Views produced per converter

| Converter | PSM | Feature | PG | Pepmap | Sample | Run | Dataset | Ontology | Provenance | mz |
| ----------- | :---: | :-------: | :--: | :------: | :------: | :---: | :-------: | :--------: | :----------: | :--: |
| **MaxQuant** | Yes | Yes | Yes | No | If SDRF | If SDRF | Yes | If SDRF | No | No |
| **FragPipe** | Yes | Yes | Yes | No | If SDRF | If SDRF | Yes | If SDRF | No | No |
| **DIA-NN** | No | Yes | Yes | No | Yes | Yes | Yes | If terms | Yes | No |
| **Spectronaut** | No | Yes | Yes | No | If SDRF | If SDRF | Yes | Yes | Yes | No |
| **OpenMS native QPX** | Yes | Yes | Yes | No | Yes | Yes | Yes | Yes | Yes | No |
| **OpenMS consensusXML** | Yes | Yes | Yes | No | If SDRF | If SDRF | Yes | If terms | Yes | No |
| **CDAP** | Yes | Yes | Yes | No | If PDC | If PDC | Yes | Yes | Yes | No |
| **mzIdentML** | Yes | No | No | Yes | If SDRF | If SDRF | Yes | Yes | Yes | No |
| **QuantMS MSstats** | No | Yes | No | No | Yes | Yes | Yes | Yes | Yes | No |
| **SDRF** | No | No | No | No | Yes | Yes | No | Optional | No | No |

- **Yes** — the converter produces this view.
- **No** — the converter does not produce this view (e.g. DIA-NN has no PSM view; mzIdentML has no Feature/PG).
- **If SDRF** — the view is produced only when an SDRF file is provided. Whether SDRF is optional or required depends on the converter; see its CLI contract below.
- **If terms** — the view is produced when the converter discovers resolvable ontology entries.
- **If PDC** — CPTAC/PDC studies ship no SDRF, and CDAP `.psm` files carry no sample metadata. When run through `qpxc pdc2qpx` (default), the sample/run views are built from PDC GraphQL metadata, which also recovers the TMT/iTRAQ channel → biological-sample mapping. Disable with `--no-metadata`.

> The `mz` (full-spectra) view is produced by the standalone `qpxc convert mz` command, or automatically by `qpxc pdc2qpx --include-spectra`; it is not emitted by the per-tool converters above.

## Field population: DIA-NN and OpenMS consensusXML

These tables explain selected fields whose population differs between DIA-NN
and OpenMS LFQ/TMT conversions. The OpenMS column refers to
`qpxc convert openms-consensus`. Its in-memory and streaming readers use the
same field mappings.

Source columns and annotations must be present and usable for a value to be
written. An optional field can remain null because the source does not report
it, the evidence is ambiguous, or the converter does not perform the required
calculation. **Not populated** below describes current converter behavior; it
does not imply that the quantity could never be derived.

### Feature fields

| Field | DIA-NN | OpenMS consensusXML |
| ------- | -------- | --------------------- |
| `missed_cleavages` | Calculated from the peptide sequence and the SDRF enzyme. Null when the enzyme is unavailable or unsupported. | Calculated from the peptide sequence and the SDRF enzyme, falling back to `SearchParameters.enzyme`. Null when the enzyme is unavailable or unsupported. |
| `cv_params` | Tool-specific metrics such as `proteotypic`, `irt`, `predicted_irt`, `iim` and `precursor_quantification_score`, when reported. | Null: the converter has no corresponding source metrics to put in this field. These are separate from `pg.cv_params`. |
| `pg_global_qvalue` | Read from `Global.PG.Q.Value`, when present. | The resolved protein group's global q-value, using the same derivation as `pg.global_qvalue`. Null when the group is ambiguous or has no declared q-value. |
| `gg_accessions` / `gg_names` | Gene symbols from `Genes`, when present. | Gene symbols parsed from `GN=` in protein-hit descriptions, through the resolved protein group. Null when that annotation is absent or the group cannot be resolved. The converter does not read a FASTA to add missing genes. |
| `rt_start` / `rt_stop` | Read from `RT.Start` / `RT.Stop` and converted from minutes to seconds. Null when those columns or values are absent. | Null: consensusXML feature elements provide an apex RT, not the feature's integration-window boundaries. The apex is stored in `rt`. |
| `id_run_file_name` | The normalized report run name. | The resolved run of a direct identification of the exported peptide with a spectrum reference. Null for transferred features without their own identification, unresolved origins or conflicting peptide assignments. |
| `pg_positions` | Not populated from the report, which has no peptide-to-protein coordinates. Filled from a FASTA with `--fasta` or `qpxc transform protein-properties`: every one-based occurrence in each group member. | Known peptide-evidence `start` / `end` coordinates, converted from zero-based inclusive to one-based inclusive. For a resolved protein group, positions on its members; for a peptide shared across inferred groups (null `pg_accessions`), positions on every evidence protein, since each entry names its own protein. Null when no usable positions remain. A FASTA is not needed for coordinates already recorded in consensusXML. |
| `peptide_qvalue` | Null: DIA-NN `Q.Value` is precursor-level and is stored as `precursor_qvalue` in `additional_scores`. | Read from an identification score explicitly declared to be a q-value/FDR. Null without a suitable identification in that run; raw search scores are not substituted. |

Protein-group ambiguity includes shared peptides spanning multiple inferred
groups and conflicting assignments within a run. These features retain their
quantification while the ambiguous protein-group fields remain null. See the
[OpenMS conversion guide](../guide/convert.md#openms-consensus-description) for
the evidence-selection rules, and [Feature confidence](feature.md#protein-group-fields)
for q-value semantics.

### Protein group fields

| Field | DIA-NN | OpenMS consensusXML |
| ------- | -------- | --------------------- |
| `pg_names` | Split from the report's `Protein.Names`. | UniProt entry names taken from `db\|ACCESSION\|NAME` accessions, aligned with `pg_accessions`. Null when any member is a bare accession with no entry name, so the list never misaligns; a bare accession is not repeated as its own name. |
| `contaminant` | True when any member accession carries the `CONTAM` marker of the quantms/sdrf-pipelines contaminant database. | Same rule, shared with DIA-NN (`qpx.converters.utils.is_contaminant_accession`). |
| `global_qvalue` | Read from `Global.PG.Q.Value`, when present. | The lowest available member-protein q-value, only when the source protein score is declared to be a q-value/FDR. Null if no member has such a score. |
| `pg_qvalue` | Read from run-level `PG.Q.Value`, when present. | Null: no separate run-level protein-group q-value is reported by this path. The experiment-level value belongs in `global_qvalue`. |
| `gg_qvalue` | Read from `GG.Q.Value`, when present. | Null: this path has no gene-group q-value. Protein-level confidence is not a gene-level q-value. |
| `additional_scores` | Includes `PG.Q.Value` as `pg_qvalue` when available. | The existing `anchor_protein`'s `Posterior Probability_score`, stored as `posterior_probability` with `higher_better: true`. Null if absent, invalid or conflicting for that anchor. |
| `sequence_coverage` | Not populated from the report. Filled from a FASTA with `--fasta` or `qpxc transform protein-properties`: percent of the anchor covered by the dataset's target peptides whose evidence names it. | The existing anchor's recorded `ProteinHit.coverage` in percent. Unknown values, missing annotations and conflicting coverage values for that anchor remain null. Other group members are not substitutes. |
| `molecular_weight` | Not populated from the report. Filled from a FASTA with `--fasta` or `qpxc transform protein-properties`: anchor average mass in kDa. | The theoretical average molecular weight of the existing anchor's complete, unmodified `ProteinHit.sequence`, in kDa. Null for missing, ambiguous, modified or conflicting sequences; other group members are not substitutes. |
| `cv_params` | Each quantified row records whether its primary intensity came from `PG.Quantity` or the `PG.MaxLFQ` fallback. | Each quantified row records `unnormalized_unique_peptide_sum`, or the corresponding top-N method selected by `--pg-top`. Null when the row has no quantity. |

For either converter, a FASTA fills these three fields **only where they are null and only on target rows**; a value
from the producer is never overwritten, and proteins absent from the FASTA (such as DIA-NN's internal decoys) stay null.
For OpenMS consensusXML this matters when the file's `ProteinHit`s carry no sequence (MSV000085836).

The schema defines `molecular_weight` in **kDa** and allows it to be null.
OpenMS consensusXML uses the existing `anchor_protein` and average isotopic
masses, including the terminal water molecule, without post-translational
modifications. The calculation uses the recorded protein sequence and does not
load a FASTA or reconstruct a protein from identified peptides. Both readers
support the standard amino acids, U/O and the isobaric I/L code J; sequences with
B/Z/X or modification notation remain null.
Converters that receive an explicit molecular-weight value, such as MaxQuant,
can preserve it; see [Protein Group mappings](pg.md#tool-mappings).

### PSM evidence and computed links

DIA-NN does not produce a PSM view. OpenMS PSM `missed_cleavages` uses the same
sequence/enzyme calculation as the feature field above. If no exportable PSM
records remain, the consensusXML converter warns and skips the PSM file and its
output registration; other requested views can still be written.

`feature.psm_ids` and `feature.pg_ids` are computed links that QPX converters do
not materialize. OpenMS `psm.feature_id` records an explicit feature link when
both views are exported; unassigned PSMs and PSM-only exports can retain valid
identifications with a null link. A null `psm.feature_id` denotes no recorded
association. The [data-model join rules](data-model.md#entity-relationship-diagram)
describe how to compute links from the evidence that was exported.

## CLI commands

| Converter | Command |
| ----------- | --------- |
| MaxQuant | `qpxc convert maxquant` |
| FragPipe | `qpxc convert fragpipe` |
| DIA-NN | `qpxc convert diann` |
| Spectronaut | `qpxc convert spectronaut` |
| OpenMS native QPX | `qpxc convert openms` |
| OpenMS consensusXML | `qpxc convert openms-consensus` |
| CDAP | `qpxc convert cdap` |
| mzIdentML | `qpxc convert mzidentml` |
| QuantMS MSstats | `qpxc convert quantms-msstats` |
| SDRF only | `qpxc convert sdrf` |

## Input files (summary)

| Converter | Typical inputs |
| ----------- | ---------------- |
| MaxQuant | msms.txt, evidence.txt, proteinGroups.txt |
| FragPipe | psm.tsv, combined_ion, combined_protein |
| DIA-NN | report (TSV or Parquet), pg_matrix (optional); required SDRF |
| Spectronaut | report.tsv; optional SDRF |
| OpenMS native QPX | `-out_qpx` Parquet directory; SDRF; optional companion consensusXML |
| OpenMS consensusXML | consensusXML; optional SDRF |
| CDAP | CPTAC CDAP `.psm` files in one study directory |
| mzIdentML | .mzid / .mzid.gz; optional MGF or mzML (file/folder) for spectra; optional SDRF |
| QuantMS MSstats | QuantMS-generated `*_msstats_in.csv`; required SDRF |
| SDRF | Single SDRF TSV file |

For field-level mappings from each tool’s columns to QPX, see [Tool Field Mappings](tool-mappings.md).
