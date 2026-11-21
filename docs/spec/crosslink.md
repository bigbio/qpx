# Cross-Links (proposed)

!!! warning "Proposal, not part of the specification"
    This page proposes how QPX holds cross-linking (XL-MS) results above the PSM level. Nothing on it is implemented: no schema, writer or validator check exists yet. Discussion is in the linked issue; field names and placement may change.

A cross-linking result is reported at three levels, each with its own error control:

| Level | Claim | View |
|-------|-------|------|
| CSM | This spectrum is two linked peptides (or one looped peptide) | [`psm`](psm.md), with the columns below |
| Residue pair | These two residues are linked in the sample | `crosslink` (new) |
| Protein pair | These two proteins interact | `interaction` (new) |

A fourth, the peptide pair, is optional: a producer that estimates it records it in `crosslink_fdr` and gives its q-value as a named score of the psm row (`peptide_pair_qvalue`, `MS:1003338`).

## Link vocabulary

Two columns describe a link, and they answer different questions:

- **`cross_links[].xl_type`** (unchanged): what joins the peptides. `inter` is a link between two peptides, `intra` a link within one peptide (mzIdentML's looplink), `dead-end` a linker with one reacted end, `noncovalent` two peptides associated without a link.
- **`xl_protein_scope`** (new): where the two peptides come from. `same_family` when both map to one protein family (homologs included), `same_sequence` when both peptides have one bare sequence (a self-link of one sequence cannot be told from a homodimer), `different_family` otherwise. Producers that estimate FDR separately by class do so per `xl_protein_scope`; mzIdentML has only the boolean `MS:1003343` (self vs heteromeric).

## psm: cross-link columns

All nullable, and null on rows that are not cross-link spectrum matches. A cross-link row carries a provided `psm_id` (`identity_source` `provided`), because the schema composite names only the alpha peptide ([Row Identity](identity.md)).

| Field | Description | Type |
|-------|-------------|------|
| `xl_protein_scope` | `same_family`, `same_sequence` or `different_family` (above) | string, null |
| `xl_peptides` | The linked peptides, alpha first as in `peptidoform`; one entry for a looplink | list\<xl_peptide\>, null |
| `xl_child_spectra` | Other spectra supporting a peptide of the CSM: MS3 scans of a cleavable linker, ETD or EThcD spectra of the same precursor | list\<xl_child_spectrum\>, null |
| `crosslink_ids` | `crosslink.crosslink_id` of the residue pairs this CSM supports (`MS:1003344`); a soft link | list\<int64\>, null |

The pair score the search ranked by and the score the CSM was finally estimated by (they differ after rescoring) are both named entries of `additional_scores`; `crosslink_fdr.score_name` says which one a q-value used.

`xl_peptides` becomes the one source for the partner peptide; `partner_protein_accessions`, `partner_protein_starts` and `protein_starts` would be deprecated in a later version. Until then they must agree with it.

The `xl_peptide` struct:

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | string | Bare sequence |
| `peptidoform` | string | ProForma of this peptide alone, link not tagged |
| `link_position` | int32, null | 1-based residue, 0 for the peptide N-terminus, length + 1 for the C-terminus |
| `link_residue` | string, null | Residue at the link; null for a terminus |
| `site_kind` | string, null | `residue`, `peptide_n_term`, `protein_n_term` or `protein_c_term` |
| `reagent_end` | string, null | End of the reagent that reacted here, for heterobifunctional linkers |
| `protein_accessions` | list\<string\> | Proteins the peptide maps to |
| `protein_starts` | list\<int32\> | 1-based start of the peptide in each of those proteins |
| `is_decoy` | bool, null | Decoy peptide; null when the producer does not say |
| `form` | string, null | Linker form observed on this peptide: `intact` or a stub name |
| `score` | float64, null | The peptide's own score |
| `form_scores` | list\<score\>, null | Per-form scores of a cleavable linker |
| `n_candidate_sites` | int16, null | Link sites the peptide offered; 1 means unambiguous |

The `xl_child_spectrum` struct (mzIdentML parent/child SIIs, `MS:1003332`):

| Field | Type | Description |
|-------|------|-------------|
| `peptide_index` | int16 | Index into `xl_peptides` |
| `spectrum_id` | string | Native id of the child spectrum |
| `ms_level` | int16 | 2 or 3 |
| `fragmentation` | string, null | e.g. `CID`, `HCD`, `ETD`, `EThcD` |
| `peptidoform` | string, null | The peptide identified in the child spectrum, with its stub |
| `score` | float64, null | Score of the child identification |
| `posterior_error_probability` | float64, null | `MS:1003336` |
| `rank` | int16, null | Rank of the child identification |

## crosslink: residue pairs

One row per **ambiguity set of linked residues**: every (protein, position) end the evidence cannot tell apart on each side. The set is the counting unit of `crosslink_qvalue`: the q-value is the estimated proportion of false sets among those at least as well supported. A set is false when no pair of its ends is linked in the sample. Two overlapping sets (`A:12 — B:5` and `A:12, C:12 — B:5`) are two rows.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `crosslink_id` | Opaque int64 identity, derived from the ends and `grouped_runs` | int64 | yes |
| `protein_a`, `position_a`, `amino_acid_a` | Anchor end of side A: the first of `ends_a`. Position as in [site](site.md): 1-based residue, 0 protein N-terminus, length + 1 C-terminus | string, int32, string null | yes |
| `protein_b`, `position_b`, `amino_acid_b` | Anchor end of side B | string, int32, string null | yes |
| `ends_a`, `ends_b` | Every end of each side, sorted; sides ordered so that `ends_a` sorts first | list\<link_end\> | yes |
| `is_ambiguous` | A side has more than one end | bool | yes |
| `xl_protein_scope` | Of the best CSM | string | yes |
| `family_a`, `family_b` | Protein families of the anchors | string, null | no |
| `pg_id_a`, `pg_id_b` | Identification pg row holding the anchor, as `site.pg_id` | int64, null | no |
| `crosslinkers` | Reagents seen in the supporting CSMs; not part of the identity | list\<string\> | yes |
| `grouped_runs` | Runs whose CSMs support the set | list\<string\> | yes |
| `csm_count` | Scored CSMs of the set | int32 | yes |
| `peptide_pair_count` | Distinct peptidoform pairs among them | int32, null | no |
| `best_psm_id` | Best CSM (`psm.psm_id`) | int64, null | no |
| `is_decoy` | The best CSM is a decoy (`target_decoy` or `decoy_decoy`); null when the producer exports targets only | bool, null | no |
| `decoy_class` | `target_target`, `target_decoy` or `decoy_decoy` of the best CSM | string, null | no |
| `score` | Score the q-value ranked by, full precision: the best CSM's | float64, null | no |
| `crosslink_qvalue` | q-value over the sets of this `xl_protein_scope` (`MS:1002677`) | float64, null | no |
| `crosslink_qvalue_min_csms` | Same estimate with only sets of at least `crosslink_fdr.min_csms` CSMs; null for a set with fewer | float64, null | no |
| `pass_threshold` | Passes the producer's threshold (`MS:1003340`) | bool, null | no |
| `fdr_estimate_id` | Row of `crosslink_fdr` this q-value comes from | int32, null | no |
| `additional_scores` | Additional scores | list\<score\> | no |
| `cv_params` | CV parameters | list\<cv_param\> | no |

The `link_end` struct: `protein_accession` (string), `position` (int32), `amino_acid` (string, null), `site_kind` (string, null), `is_decoy` (bool, null). A side is a target when it has a target end.

## interaction: protein pairs

One row per protein pair claimed from links between different families.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `interaction_id` | Opaque int64 identity, derived from the protein sets and `grouped_runs` | int64 | yes |
| `protein_a`, `protein_b` | Anchor proteins, the first of each set | string | yes |
| `proteins_a`, `proteins_b` | Proteins each side cannot be told apart from (homologs, shared peptides) | list\<string\> | yes |
| `family_a`, `family_b`, `pg_id_a`, `pg_id_b` | As in `crosslink` | string null, int64 null | no |
| `grouped_runs` | Runs whose CSMs support the pair | list\<string\> | yes |
| `residue_pair_count` | `crosslink` rows supporting the pair | int32 | yes |
| `csm_count` | Scored CSMs supporting the pair | int32 | yes |
| `best_crosslink_id`, `best_psm_id` | Best supporting residue pair and CSM | int64, null | no |
| `is_decoy`, `decoy_class` | As in `crosslink` | bool null, string null | no |
| `score` | Score the q-value ranked by | float64, null | no |
| `interaction_qvalue` | q-value over the protein pairs (`MS:1002676`) | float64, null | no |
| `pass_threshold` | `MS:1003341` | bool, null | no |
| `fdr_estimate_id` | Row of `crosslink_fdr` | int32, null | no |
| `additional_scores`, `cv_params` | As elsewhere | list | no |

## crosslink_fdr: how each q-value was estimated

One row per population a q-value was estimated over. Without it, a reader of a q-value cannot tell which CSMs competed, which decoys were counted, or what the formula was.

| Field | Description | Type |
|-------|-------------|------|
| `fdr_estimate_id` | Identity of the population | int32 |
| `level` | `csm`, `peptide_pair`, `residue_pair` or `interaction` | string |
| `xl_protein_scope` | The class estimated, or null when classes were pooled | string, null |
| `score_name` | Score ranked by, as named in `additional_scores` | string |
| `score_training` | For a learned score (e.g. Percolator), the population it was trained on | string, null |
| `estimator` | The formula, written out: e.g. `fdr = min(1, max(TD + 1 - DD, 2) / TT); fdr = 1 if TT = 0; ties take the value at the end of their group; q = running minimum from the lowest score` | string |
| `min_csms` | Eligibility of `crosslink_qvalue_min_csms` | int32, null |
| `decoy_model` | How decoys were made, as CV terms (`MS:1001194`, `MS:1001195` / `MS:1001196`, `MS:1001197`) | list\<cv_param\> |
| `n_target_target`, `n_target_decoy`, `n_decoy_decoy` | Rows of each class in the population | int32 |
| `dropped` | Rows left out before estimating, and why (e.g. no score) | string, null |
| `runs` | Runs of the population | list\<string\> |

A q-value can be recomputed from `score`, `decoy_class` and `estimator` only when the decoy rows are kept and the score is stored at full precision; a validator check would be opt-in. Looplinks and noncovalent pairs get a null q-value unless a producer defines and records an estimator for them.

When one FDR is estimated over several projects, `crosslink_fdr_shares` gives each project's part of the CSMs passing a pooled cut, to show whether pooling moved errors into a weaker project:

| Field | Description | Type |
|-------|-------------|------|
| `fdr_estimate_id` | The pooled population | int32 |
| `project` | Project label, as the producer names its inputs | string |
| `q_cut` | Pooled CSM q cut | float64 |
| `n_target_target`, `n_target_decoy`, `n_decoy_decoy` | That project's CSMs at or below the cut | int32 |
| `share_estimate`, `share_estimator` | The project's own decoy estimate and its formula (e.g. `max(TD - DD, 0) / TT`) | float64, string |

## run: cross-linkers

`run.crosslinkers` lists the reagents searched; `comment[cross-linker]` in the SDRF fills it where present.

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Reagent name |
| `accession` | string, null | XLMOD accession |
| `intact_mass` | float64 | Mass of the intact link (Da); mzIdentML donor `massDelta` |
| `cleavable` | bool | MS-cleavable |
| `end_1`, `end_2` | string | Reactive sites of each end (`K, S, T, Y, Protein N-term`) |
| `stubs` | list\<{name, mass, pairs_with}\>, null | Stub fragments of a cleavable linker (`MS:1003390`) |
| `label` | string, null | Isotope label of a light/heavy reagent pair |

## Relation to mzIdentML 1.3

| QPX | mzIdentML 1.3 |
|-----|---------------|
| CSM (`psm` row, `xl_peptides`) | Two SIIs sharing `MS:1002511`; peptides with donor/acceptor modifications `MS:1002509` / `MS:1002510` |
| `xl_child_spectra` | Child SIIs grouped by `MS:1003332`, PEP `MS:1003336` |
| Looplink (`xl_type` `intra`) | One SII with `MS:1003329` |
| Peptide-pair q | `MS:1003338` |
| `crosslink_qvalue`, `pass_threshold` | `MS:1002677`, `MS:1003340` in a ProteinAmbiguityGroup |
| `interaction_qvalue`, `pass_threshold` | `MS:1002676`, `MS:1003341` |
| `xl_protein_scope`, `crosslink_fdr` | Only `MS:1003343` (FDR separate for self and heteromeric links) |
| `run.crosslinkers` | `SearchModification` with donor/acceptor terms, XLMOD accession, `MS:1003390` |

mzIdentML's FDR terms are global; QPX q-values are per population, as `crosslink_fdr` states. mzIdentML's donor/acceptor ordering rule (longer peptide, then mass, then alphabetical) is a fallback for unmarked or symmetric reagents; QPX keeps the alpha/beta order of `peptidoform` and records the reagent end per peptide.

The mzIdentML converter would need two changes: read residue-pair and protein-pair ProteinAmbiguityGroups into these views instead of `pg`, and apply the side rule above (a side is a target when it has a target accession) for decoys in ambiguity sets.

## Open questions

1. Cross-link columns in `psm`, or a separate `csm` view?
2. Names: `crosslink`, `interaction`, `crosslink_fdr`, `crosslink_fdr_shares`; tables or provenance parameters?
3. Ambiguity as lists (`ends_a`, `ends_b`) or one row per alternative with a group id, as `modification_candidates` does?
4. Per-project shares in the format, or only in producer reports?
5. Calibration columns (entrapment or truth-set error) in production files, or only in benchmarks?
6. Schema, `file_type`, collection and merge registration, versioning, ontology and converter-coverage entries: one PR per view or one for all?
