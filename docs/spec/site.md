# Modification Sites

`site.parquet` holds **protein modification-site claims**: one row per claim that a residue of a protein carries a given modification, in the runs listed in `grouped_runs`. Three columns describe a claim, and they answer different questions:

- **`site_qvalue`** — error control over the claims of this view. A claim is **false** when the residue does not carry this modification in the sample, whatever the reason: the peptide or the protein is wrong, the modification is on another residue, or it is another modification. The **counting unit** is the distinct row (protein, position, modification, runs), not the PSM: the q-value is the estimated proportion of false claims among the rows at least as well supported. It is a set-level statistic, not the probability that one row is correct. A producer that has no estimate at this unit leaves it null, and records the method and the population it estimated over in [provenance](provenance.md).
- **`modified_probability`** — the probability that the residue carries *some* modification, from the spectra of the peptides that cover it. It is **not stoichiometry** (the fraction of the protein's copies that are modified), which this view does not report.
- **`evidence_fraction`** — the share of the residue's modification evidence that supports *this* modification. Within one spectrum the candidate explanations are mutually exclusive and get posterior probabilities; the probability mass placing a modification on the residue is summed over spectra, and the evidence fraction is the part of it placing this modification. With the `alternatives` of the row it sums to 1; a row without alternatives holds the residue's whole evidence, so its fraction is 1. Mass spectrometry often cannot tell near-isobaric modifications apart (acetyl and trimethyl differ by 0.036 Da), and the fractions show how the evidence splits between them.

!!! warning "An evidence fraction is not a presence probability"
    Two modifications that both occur at a residue — on different molecules, or in different runs — are two rows, each with its own `site_qvalue`. Their evidence fractions can both be near 0.5 although both are confidently present. Read presence from `site_qvalue`, and use the evidence fraction to see whether the data tell competing modifications apart.

The YAML schema definition lives at `qpx/core/data/schemas/site.yaml` in the source repository.

## Use Cases

- **Site-level PTM reporting**: report a modification at a protein position with error control over the reported sites, independent of the PSM q-value.
- **Ambiguous modifications**: keep both readings of a residue (acetyl 0.71, trimethyl 0.29) instead of a forced choice.
- **Prior annotations**: tell sites that confirm a known annotation (`prior_evidence`) from new ones.

## Schema

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `site_id` | Opaque int64 identity, derived from the identity composite (`protein_accession`, `position`, `modification_name`, `grouped_runs`, `label`) | int64 | yes |
| `protein_accession` | Protein the position refers to | string | yes |
| `family` | Protein family (the canonical entry of an isoform, variant or chain) | string, null | no |
| `position` | Position in the protein: 1-based for a residue, 0 for the protein N-terminus, length + 1 for the C-terminus (see below) | int32 | yes |
| `amino_acid` | Residue at the position; null for a protein terminus | string, null | yes |
| `modification_name` | Modification name (e.g. `Phospho`); part of the identity | string | yes |
| `modification_accession` | Modification accession (e.g. `UNIMOD:21`) | string, null | no |
| `mass_delta` | Monoisotopic mass shift (Da); not part of the identity, so modifications of equal mass (dimethyl, ethyl) are distinct rows | float64 | yes |
| `grouped_runs` | Runs whose PSMs support the claim | list\<string\> | yes |
| `pg_id` | Protein group row holding the protein (see below) | int64, null | no |
| `prior_evidence` | Evidence class of the claim before the search (e.g. `curated`, `catalogue`, `sample`) | string, null | no |
| `site_qvalue` | q-value over the site-modification claims (above) | float64, null | no |
| `modified_probability` | Probability that the residue carries some modification | float64, null | no |
| `evidence_fraction` | Share of the residue's modification evidence supporting this modification | float64, null | no |
| `alternatives` | The other modifications the same evidence places on the residue, with their evidence fractions, largest first | list\<site_alternative\>, null | no |
| `psm_count` | PSMs at the PSM cutoff that place the modification there | int32, null | no |
| `best_psm_id` | Best of those PSMs (`psm.psm_id`) | int64, null | no |
| `is_decoy` | Decoy or knockoff claim | bool | yes |
| `label` | Channel of this quantification, null for a claim without one; part of the identity; specification 1.3 | string, null | no |
| `intensity` | Quantity of the modified site in this label and these runs; specification 1.3 | float32, null | no |
| `stoichiometry` | Fraction of the protein's copies at this residue that carry the modification; specification 1.3 | float64, null | no |
| `additional_scores` | Additional scores, e.g. a localization probability | list\<score\> | no |
| `cv_params` | CV parameters | list\<cv_param\> | no |

`qpxc transform relayout` can store `additional_scores` and `cv_params` as a struct with one field per name ([Struct lists](serialization.md#struct-lists)).

The `site_alternative` struct:

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Modification name |
| `accession` | string, null | Modification accession |
| `mass_delta` | float64 | Monoisotopic mass shift (Da) |
| `evidence_fraction` | float64 | Its share of the residue's modification evidence |

!!! note "Linking to pg"
    `pg_id` names one pg row, and a pg row has one label. A site has no label, so it links to the identification pg row (null `label`) with the same `grouped_runs`. When the site spans several pg rows — several labels, or runs grouped differently in pg — `pg_id` is null, and the protein is found through `pg_accessions` membership.

!!! note "Positions and termini"
    Positions follow the [Modifications](modifications.md#position-format-rules) convention, applied to the protein instead of the peptide: a residue is 1-based and the modification is on its side chain, the protein N-terminus is `0` and the C-terminus is the protein length + 1, both with a null `amino_acid`. Protein N-terminal acetylation of Met1 is position `0`; acetylation of the side chain of a Lys at residue 1 is position `1`.

!!! note "Aggregated over spectra, unlike psm and feature"
    A site row summarises every spectrum that places a modification on the residue. The `psm` view reports one identification, and the `feature` view carries the modifications of its best PSM; neither aggregates. The competing explanations of one spectrum (acetyl or trimethyl on the same peptide) are separate PSM rows, each with its own scores; the site view is where their evidence is combined.

!!! note "Not localization"
    None of these is a localization score. Which residue of a peptide carries a modification is reported per PSM in `psm.modifications[].positions[].scores`, under the [localization score names](modifications.md#localization-scores) (e.g. `localization_probability`, `luciphor_site_probability`, see [Scores](scores.md)). A localization-only error rate (false localization rate) is a different quantity from `site_qvalue`, which also counts wrong peptides, proteins and modifications; a producer that estimates one reports it in `additional_scores`.

## Validation

Under `qpxc validate` (strict):

- `site_qvalue`, `modified_probability`, `evidence_fraction` and the alternatives' fractions are finite and in [0, 1];
- a row's evidence fraction and its alternatives' sum to 1, and a row without alternatives has a fraction of 1 (or null);
- alternatives are listed in decreasing evidence fraction;
- a populated `pg_id` resolves to a pg row with a null label and the same `grouped_runs`;
- that pg row's `pg_accessions` hold the site's `protein_accession`;
- a populated `best_psm_id` resolves to a `psm.psm_id` when the dataset has a psm view.

`prior_evidence` is an open vocabulary: producers name their own evidence classes.

## Python API

```python
import qpx

ds = qpx.open_dataset("PXD014414/")

# Sites on one protein
tp53 = ds.site.by_protein("P04637")

# Acetylation sites at site q <= 0.05 on residues modified with probability >= 0.5
acetyl = ds.site.targets_only().by_modification("Acetyl").confident(0.05, min_modified_probability=0.5)

# ... where the evidence tells acetyl from its alternatives
resolved = acetyl.confident(None, min_evidence_fraction=0.95)
print(resolved.to_df())
```

## Relationship to Other Views

- `pg_id` links a site to its protein group in [pg](pg.md).
- `best_psm_id` links to the [psm](psm.md) row with the best evidence.
- A [proteoform](proteoform.md) of the same `family` can carry the residue in a variant or chain. Site rows do not establish that several modifications occur on the same molecule.

## Quantified sites

A quantified row holds one site, label and set of runs, and `intensity` with `stoichiometry` say how much of it there is. Do not read `modified_probability` as stoichiometry, which is what these columns are for. A pg link (`pg_id`) stays null on a labelled row, because validation expects the link to name a null-label pg row.

| Producer | What it writes |
|----------|----------------|
| `qpxc convert maxquant --sites-file` | Each `<modification>Sites.txt` per experiment: `Intensity [L/M/H] <experiment>`, or `Reporter intensity corrected N` with an isobaric SDRF; `Occupancy` as `stoichiometry`; `Localization prob` as `modified_probability`; scores and `PEP` in `additional_scores`; the other proteins and positions in `cv_params`. `grouped_runs` holds the experiment's raw files, from the evidence file or SDRF. `psm_count` covers all experiments. |
| `qpxc transform site --modification` | For any converter whose feature view has `pg_positions`: per run and label, the summed intensity of the peptides that place the modification, on their `anchor_protein`. `stoichiometry` is their share of the intensity of all peptides covering the residue, an apparent value since peptides ionize differently. `modified_probability` is the best probability score of the placement. |

DIA-NN and Spectronaut site tables are not read: use `qpxc transform site`, which reads the localization both converters keep on the feature view.
