# Crosslinks (proposed)

!!! warning "Proposal, not part of the specification"
    This page proposes how QPX holds crosslinking (XL-MS) results above the PSM level. Nothing on it is implemented: no schema, writer or validator check exists yet. Discussion is in issue #394; field names and placement may change.

A crosslinking result is reported at three levels, each with its own error control:

| Level | Claim | View |
|-------|-------|------|
| CSM | This spectrum is two linked peptides (or one looped peptide) | [`psm`](psm.md), with the columns below |
| Residue pair | These two residues are linked in the sample | `residue_pair` (new) |
| Protein pair | These two proteins are linked (interact, or a homo-oligomer for one protein) | `protein_pair` (new) |

A peptide-pair level is optional: a producer that estimates it gives the CSM row a `peptide_pair_qvalue` score and records the estimate in `crosslink_fdr` (`MS:1003338`).

## Link classes

Three columns describe a link, and they answer different questions:

- **`cross_links[].xl_type`** (unchanged): what joins the peptides. `inter` is a link between two peptides, `intra` a link within one peptide (mzIdentML's looplink), `dead-end` a linker with one reacted end, `noncovalent` two peptides associated without a link.
- **`xl_protein_scope`** (new): where the two peptides come from, the split of `MS:1003343` (FDR applied separately to self and heteromeric crosslinks), tested in this order, the first that holds wins:
    1. `self`: both peptides come from one protein sequence, the CV's self crosslink.
    2. `self_proteoform_family`: the two peptides come from distinct sequences of one proteoform family, as in [proteoforms](proteoform.md) (Shortreed & Smith: a gene's canonical entry, its isoforms, variants and chains); homologs are in different proteoform families. `MS:1003343` defines heteromeric as "crosslinks between distinct protein sequences", so the CV would call these heteromeric. QPX counts them as self, because their peptides usually cannot tell the sequences apart; a producer that follows the CV text puts them with `heteromeric`, and `crosslink_fdr` says which population they were estimated in.
    3. `heteromeric`: otherwise.
- **`xl_peptides_overlap`** (new, bool): in a `self` link, the two peptides are the same or overlap in the sequence. One chain cannot hold the same stretch twice (unless the sequence repeats), so this is the case that indicates a homo-oligomer (xiNET's homomultimeric link); other self links cannot be told from a link within one chain.
- **`xl_partner_neighbor`** (new, bool): a peptide of the CSM has a neighbor (Lin, Plubell, Keich & Noble 2021): a different peptide of near-equal precursor mass and similar fragment spectrum, in a protein the CSM does not list. Its partner may then be a homolog of the true one, an error reversed or shuffled decoys cannot model. A producer that estimates such CSMs apart says so in `crosslink_fdr`; on a proteome-sized database they are most of the false heteromeric CSMs a rescorer accepts. The producer states its neighbor criterion in `crosslink_fdr.neighbor_rule`. An I/L-identical copy is not a neighbor: it is ordinary protein ambiguity, and its protein is listed in the peptide's `protein_accessions`.

## Decoys

One rule throughout:

- A **peptide** is a decoy when every protein it maps to is a decoy (`xl_peptides[].is_decoy`); a peptide that maps to a target and a decoy is a target.
- A **CSM** has `crosslink_decoy_class` `target_target`, `target_decoy` or `decoy_decoy` from its peptides, and `is_decoy` is true for the latter two.
- A **residue pair** or **protein pair** takes the `decoy_class` of its best CSM, the CSM whose score ranks it.

The mzIdentML converter marks a CSM a decoy when any peptide evidence is a decoy (`crosslinks.py`, `_is_decoy`). Under this rule it would change for a peptide that maps to both a target and a decoy accession.

## psm: crosslink columns

All nullable, and null on rows that are not crosslink spectrum matches. Crosslink PSMs carry provided ids, because the schema composite names only the alpha peptide: the file footer states `identity_source` `provided` or `mixed` ([Row Identity](identity.md)).

| Field | Description | Type |
|-------|-------------|------|
| `xl_protein_scope` | `self`, `self_proteoform_family` or `heteromeric` (above) | string, null |
| `xl_peptides_overlap` | In a `self` link, the two peptides overlap: a candidate homo-oligomer link (above) | bool, null |
| `xl_partner_neighbor` | A peptide has a neighbor in a protein the CSM does not list (above) | bool, null |
| `xl_peptides` | Per-peptide evidence, alpha first as in `peptidoform`; one entry for a looplink | list\<xl_peptide\>, null |
| `xl_child_spectra` | Other spectra supporting a peptide: MS3 scans of a cleavable linker, ETD or EThcD spectra of the same precursor | list\<xl_child_spectrum\>, null |
| `residue_pair_ids` | `residue_pair.residue_pair_id` of the residue pairs this CSM supports (`MS:1003344`); a soft link | list\<int64\>, null |
| `fdr_estimate_id` | `crosslink_fdr` row of the CSM q-value | int32, null |

The CSM q-value is the PSM's q-value column. The PEP of the combined identification (`MS:1003336`, shared by the parent and child spectra in mzIdentML) is the PSM's `posterior_error_probability`. The pair score the search ranked by and the score the CSM was finally estimated by are named entries of `additional_scores`; `crosslink_fdr.score_name` says which one a q-value used.

The link positions stay where they are: `cross_links[]` with `donor_position` and `acceptor_position` (a looplink has both on one peptide), `crosslinked_peptidoform` and `partner_sequence` / `partner_peptidoform`. `xl_peptides` holds only the evidence those columns lack, so it does not repeat a sequence or a position.

The `xl_peptide` struct:

| Field | Type | Description |
|-------|------|-------------|
| `protein_accessions` | list\<string\> | Proteins the peptide maps to (for the alpha peptide, equal to `protein_accessions`) |
| `protein_starts` | list\<int32\> | 1-based start of the peptide in each of those proteins |
| `is_decoy` | bool, null | Decoy peptide (rule above); null when the producer does not say |
| `stub` | string, null | Linker stub observed on this peptide (`MS:1003346`, `MS:1003390`), or `intact` |
| `score` | float64, null | The peptide's own score |
| `stub_scores` | list\<score\>, null | Per-stub scores of a cleavable linker |
| `n_candidate_sites` | int16, null | Link sites the peptide offered; 1 means unambiguous |
| `linker_end` | string, null | End of the reagent that reacted, for heterobifunctional linkers |

With `xl_peptides[1].protein_accessions` and `protein_starts`, `partner_protein_accessions` and `partner_protein_starts` repeat a value; one of the two would be deprecated.

The `xl_child_spectrum` struct (mzIdentML child SIIs of `MS:1003332`):

| Field | Type | Description |
|-------|------|-------------|
| `peptide_index` | int16 | Index into `xl_peptides` |
| `run_file_name` | string, null | Null for the CSM's own run |
| `scan` | int32 | Scan number, as `psm.scan` |
| `ms_level` | int16 | 2 or 3 |
| `fragmentation` | string, null | e.g. `CID`, `HCD`, `ETD`, `EThcD` |
| `peptidoform` | string, null | The peptide identified in the child spectrum, with its stub |
| `score` | float64, null | Score of the child identification |
| `rank` | int16, null | Rank of the child identification |

## residue_pair

One row per **ambiguity set of linked residues**: every (protein, position) site the evidence cannot tell apart on each side. The set is the counting unit of `residue_pair_qvalue`: the q-value is the estimated proportion of false sets among those at least as well supported. A set is false when no pair of its sites is linked in the sample. Two overlapping sets (`A:12 — B:5` and `A:12, C:12 — B:5`) are two rows.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `residue_pair_id` | Opaque int64 identity, hashed from the canonical key below | int64 | yes |
| `protein_a`, `position_a`, `amino_acid_a` | Anchor site of side A: the first of `sites_a` | string, int32, string null | yes |
| `protein_b`, `position_b`, `amino_acid_b` | Anchor site of side B | string, int32, string null | yes |
| `sites_a`, `sites_b` | Every site of each side (`link_site`), sorted by (`protein_accession`, `position`); the side whose first site sorts first is A | list\<link_site\> | yes |
| `is_ambiguous` | A side has more than one site | bool | yes |
| `xl_protein_scope`, `xl_peptides_overlap`, `xl_partner_neighbor` | Of the best CSM | string, bool, bool | yes |
| `family_a`, `family_b` | Proteoform families of the anchors, as in [proteoforms](proteoform.md) | string, null | no |
| `pg_id_a`, `pg_id_b` | Identification pg row holding the anchor, as `site.pg_id` | int64, null | no |
| `linkers` | Reagents seen in the supporting CSMs; not part of the identity | list\<string\> | yes |
| `grouped_runs` | Runs whose CSMs support the set | list\<string\> | yes |
| `csm_count` | Scored CSMs of the set | int32 | yes |
| `peptide_pair_count` | Distinct peptidoform pairs among them | int32, null | no |
| `best_psm_id` | Best CSM (`psm.psm_id`) | int64, null | no |
| `is_decoy`, `decoy_class` | Of the best CSM (rule above); null when the producer exports targets only | bool null, string null | no |
| `score` | Score the q-value ranked by, full precision: the best CSM's | float64, null | no |
| `residue_pair_qvalue` | q-value over the sets of its population (`MS:1002677`) | float64, null | no |
| `residue_pair_qvalue_min_csms` | q-value over the sets of at least `min_csms` CSMs only; its own population | float64, null | no |
| `fdr_estimate_id`, `fdr_estimate_id_min_csms` | `crosslink_fdr` rows of the two q-values | int32, null | no |
| `pass_threshold` | Passes the producer's threshold (`MS:1003340`) | bool, null | no |
| `additional_scores`, `cv_params` | As elsewhere | list | no |

The `link_site` struct: `protein_accession` (string), `position` (int32), `amino_acid` (string, null), `site_kind` (string, null), `is_decoy` (bool, null). `position` follows [site](site.md): 1-based residue, 0 for the protein N-terminus, length + 1 for the C-terminus. `site_kind` is new: `residue`, `protein_n_term` or `protein_c_term`. A link at a peptide N-terminus that is not a protein N-terminus (crosslinking after digestion) has `residue` and the position of that first residue.

The **canonical key** of `residue_pair_id` is `sites_a` and `sites_b` as the sorted lists of (`protein_accession`, `position`) pairs, then `grouped_runs` sorted, serialized as in [Row Identity](identity.md).

## protein_pair

One row per claimed protein pair. A `heteromeric` pair is two proteins; a `self` link gives a pair of one protein with itself, a candidate homo-oligomer (a strong one when `xl_peptides_overlap`), as in mzIdentML 1.3, which reports self pairs under `MS:1002676`. A producer that estimates protein-pair FDR for heteromeric pairs only leaves self pairs out.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `protein_pair_id` | Opaque int64 identity, hashed from the protein sets and `grouped_runs` | int64 | yes |
| `protein_a`, `protein_b` | Anchor proteins, the first of each set; equal for a self pair | string | yes |
| `proteins_a`, `proteins_b` | Proteins each side cannot be told apart from (shared peptides, near-copies) | list\<string\> | yes |
| `self_pair` | Both sides are one protein | bool | yes |
| `family_a`, `family_b`, `pg_id_a`, `pg_id_b` | As in `residue_pair` | string null, int64 null | no |
| `grouped_runs` | Runs whose CSMs support the pair | list\<string\> | yes |
| `residue_pair_count`, `csm_count` | Supporting residue pairs and CSMs | int32 | yes |
| `best_residue_pair_id`, `best_psm_id` | Best supporting residue pair and CSM | int64, null | no |
| `is_decoy`, `decoy_class` | Of the best CSM | bool null, string null | no |
| `score` | Score the q-value ranked by | float64, null | no |
| `protein_pair_qvalue` | q-value over the protein pairs of its population (`MS:1002676`) | float64, null | no |
| `pass_threshold` | `MS:1003341` | bool, null | no |
| `fdr_estimate_id` | `crosslink_fdr` row | int32, null | no |
| `additional_scores`, `cv_params` | As elsewhere | list | no |

## crosslink_fdr: how each q-value was estimated

One row per population a q-value was estimated over: the CSM, peptide-pair, residue-pair and protein-pair q-values each name theirs (`fdr_estimate_id`), and the min-CSM variant is its own population with its own counts.

| Field | Description | Type |
|-------|-------------|------|
| `fdr_estimate_id` | Identity of the population | int32 |
| `analysis` | The analysis the population belongs to: the runs and settings searched together (one per producer run, or one per cohort group) | string |
| `level` | `csm`, `peptide_pair`, `residue_pair` or `protein_pair` | string |
| `xl_protein_scope` | The class estimated, or null when the classes were estimated together | string, null |
| `xl_partner_neighbor` | The population is the CSMs with a neighbor partner (true), the others (false), or both (null) | bool, null |
| `neighbor_rule` | The producer's neighbor criterion, e.g. `same length, <= 2 substitutions, I = L` | string, null |
| `min_csms` | Only sets with at least this many CSMs; null for none | int32, null |
| `score_name` | Score ranked by, as named in `additional_scores` | string |
| `score_training` | For a learned score, the population it was trained on | string, null |
| `estimator` | `td_minus_dd` (false targets = TD − DD) or `td_plus_dd` (false targets = TD + DD, for classes whose decoys are only TD and DD) | string |
| `pseudo_count` | Added to the false-target count | int32 |
| `min_false` | Floor of the false-target count | int32 |
| `n_target_target`, `n_target_decoy`, `n_decoy_decoy` | Rows of each class in the population | int32 |
| `dropped` | Rows left out before estimating, and why (e.g. no score) | string, null |
| `decoy_model` | How decoys were made, as CV terms (`MS:1001194`, `MS:1001195` / `MS:1001196`, `MS:1001197`) | list\<cv_param\> |
| `runs` | Runs of the population | list\<string\> |

The FDR at a rank is `min(1, max(false + pseudo_count, min_false) / TT)`, with `false` from `estimator` over the rows at or above the rank, and FDR 1 when TT is 0; tied scores take the value at the end of their group, and the q-value is the running minimum from the lowest score up. andes uses `td_minus_dd`, `pseudo_count` 1 and `min_false` 2 (`td_plus_dd` for links of one sequence). A q-value can be recomputed from `score`, `decoy_class` and these fields only when the decoy rows are kept and the score is stored at full precision; a validator check would be opt-in. Looplinks and noncovalent pairs get a null q-value unless a producer records an estimator for them.

When one FDR is estimated over several projects, `crosslink_group_fdr` gives the group FDR (separate FDR, Lin et al. 2021; subgroup FDR, Fu & Qian) of each project at the threshold of the combined estimate, to show whether combining moved errors into a weaker project. It is descriptive, not a q-value, so it has its own formula: `max(TD − DD, 0) / TT` (`TD + DD` for links of one sequence), without pseudo-count or floor.

| Field | Description | Type |
|-------|-------------|------|
| `fdr_estimate_id` | The combined CSM population | int32 |
| `project` | Project label, as the producer names its inputs | string |
| `q_cut` | CSM q cut of the combined estimate | float64 |
| `n_target_target`, `n_target_decoy`, `n_decoy_decoy` | That project's CSMs at or below the cut | int32 |
| `group_fdr` | The project's own estimate by the formula above | float64 |

## run: linkers

`run.linkers` lists the reagents searched. The SDRF `comment[cross-linker]` gives only the name and the XLMOD accession; the other fields come from XLMOD or the search settings, and are null when neither has them.

| Field | Type | Description |
|-------|------|-------------|
| `linker_name` | string | Reagent name, as `cross_links[].linker_name` |
| `linker_accession` | string, null | XLMOD accession, as `cross_links[].linker_accession` |
| `linker_mass` | float64, null | Mass of the intact link (Da), as `cross_links[].linker_mass`; mzIdentML donor `massDelta` |
| `cleavable` | bool, null | MS-cleavable |
| `end_1`, `end_2` | string, null | Reactive sites of each end (`K, S, T, Y, Protein N-term`) |
| `stubs` | list\<{name, mass, pairs_with}\>, null | Stub fragments of a cleavable linker (`MS:1003390`) |
| `label` | string, null | Isotope label of a light/heavy reagent pair |

## Relation to mzIdentML 1.3

| QPX | mzIdentML 1.3 |
|-----|---------------|
| CSM (`psm` row) | Two SIIs sharing `MS:1002511`; peptides with donor/acceptor modifications `MS:1002509` / `MS:1002510` |
| `xl_child_spectra` | Child SIIs grouped by `MS:1003332`; the group's PEP `MS:1003336` is `psm.posterior_error_probability` |
| Looplink (`xl_type` `intra`) | One SII with `MS:1003329` |
| CSM q-value (the PSM's q-value) | `MS:1003337` (crosslinked PSM-level global FDR) |
| Peptide-pair q | `MS:1003338`; passes threshold `MS:1003339` |
| `residue_pair_qvalue`, `pass_threshold` | `MS:1002677`, `MS:1003340`: encoded strings (`ID.a\|b:position:score:pass`) on the ProteinDetectionHypothesis of ordinary protein groups, the two sides paired by ID |
| `protein_pair_qvalue`, `pass_threshold` | `MS:1002676`, `MS:1003341`, encoded the same way; self pairs included |
| `xl_protein_scope` | `MS:1003343` splits self (`self`; QPX adds `self_proteoform_family`) from heteromeric (`heteromeric`) |
| `run.linkers` | `SearchModification` with donor/acceptor terms, XLMOD accession, `MS:1003390` |

In the CV, "global" FDR is the counterpart of local FDR (`MS:1002350` against `MS:1002351`), not FDR over all classes together; `MS:1003343` is the term that says the classes were estimated apart. q-values are a third family (`MS:1002354`), and there is no residue-pair or protein-pair q-value term: mapping `residue_pair_qvalue` to `MS:1002677` and `protein_pair_qvalue` to `MS:1002676` puts a q-value under a global-FDR term. Requesting residue-pair and protein-pair q-value terms from PSI would close that gap. mzIdentML's donor/acceptor ordering rule (longer peptide, then mass, then alphabetical) is a fallback for unmarked or symmetric reagents; QPX keeps the alpha/beta order of `peptidoform`, and records the reagent end per peptide.

The mzIdentML converter would keep these protein groups in `pg`, and in addition parse the encoded residue-pair and protein-pair values into `residue_pair` and `protein_pair`, pairing the `.a` and `.b` sides by their ID; today `pg_adapter._float_value` returns null for them and they are lost. It would also apply the decoy rule above.

## Possible follow-up views

Like `feature` and `pg` hold quantities and `site` holds positions, two views would carry crosslinks further. Neither is part of this proposal; they are listed so the identity and keys above leave room for them.

- **Quantified crosslinks** (`crosslink_feature`): one row per residue pair (or CSM-level species) per run and label, with its intensity, as `feature` does for peptidoforms. Quantitative XL-MS compares conformations between conditions, with light/heavy reagents or label-free intensities.
- **Crosslinks on structures** (`crosslink_structure`): one row per residue pair per structure model (PDB entry, AlphaFold or ModelArchive model, with version and chains), with the residues' positions in the model, the Cα–Cα (or side-chain) distance, the linker's maximum span and whether the link fits it. A structure-centric view lets a reader validate links against structures and compare models.

## Open questions

1. Crosslink columns in `psm`, or a separate `csm` view?
2. `crosslink_fdr` and `crosslink_group_fdr`: tables, or provenance parameters?
3. Ambiguity as lists (`sites_a`, `sites_b`) or one row per alternative with a group id, as `modification_candidates` does?
4. Per-project group FDR in the format, or only in producer reports?
5. Calibration columns (entrapment or truth-set error) in production files, or only in benchmarks?
6. Which of `partner_protein_*` and `xl_peptides[].protein_*` to keep?
7. Schema, `file_type`, collection and merge registration, versioning, ontology and converter-coverage entries: one PR per view or one for all?
8. Are the follow-up views worth a proposal of their own?
9. `cross_links[].xl_type` `inter` / `intra`: most XL-MS software reads these at protein level (an intralink within one protein, a looplink within one peptide), which is why `xl_protein_scope` is a separate column. Renaming the values to `crosslink` / `looplink` / `dead-end` / `noncovalent` at the next schema version would remove the clash.

## References

- psi-ms.obo 4.2.17 and the mzIdentML 1.3 crosslinking extension.
- Leitner et al., Structure 2020, https://doi.org/10.1016/j.str.2020.09.011
- Lenz et al., Nat Commun 2021, https://doi.org/10.1038/s41467-021-23666-z
- Fischer & Rappsilber, Anal Chem 2017, https://pmc.ncbi.nlm.nih.gov/articles/PMC5423704
- Lin, Plubell, Keich & Noble, J Proteome Res 2021, https://pmc.ncbi.nlm.nih.gov/articles/PMC8489664
- xiFDR, FDR groups, https://xifdr-py.readthedocs.io/en/stable/how_it_works.html
- Shortreed & Smith, proteoform families, https://pmc.ncbi.nlm.nih.gov/articles/PMC9676067/
