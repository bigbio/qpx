# Crosslinks (proposed)

!!! warning "Proposal, not part of the specification"
    This page proposes how QPX holds crosslinking (XL-MS) results above the PSM level. Nothing on it is implemented: no schema, writer or validator check exists yet. Discussion is in issue #394; field names and placement may change.

A crosslinking result is reported at three levels, the levels most XL-MS software reports (xiFDR, Scout, MaxLynx, XlinkX, pLink), each with its own q-value:

| Level | Claim | View |
|-------|-------|------|
| CSM (crosslink spectrum match) | This spectrum is two linked peptides, or one looped peptide | [`psm`](psm.md), with the columns below |
| Residue pair | These two residues are linked in the sample | `residue_pair` (new) |
| Protein pair (PPI) | These two proteins are linked: they interact, or form a homo-oligomer when both sides are one protein | `protein_pair` (new) |

A producer that also controls the error of peptide pairs gives the CSM row a `peptide_pair_qvalue` score (`MS:1003338`).

The page keeps what a reader needs to use and judge the links: who is linked to whom, at which residues, how well supported, and at what error. Search-engine internals (per-ion scores, candidate lists, encoded mzIdentML values) stay in the producer's own files; the [Relation to mzIdentML 1.3](#relation-to-mzidentml-13) table maps the fields that have a counterpart there.

## Link types

Two columns describe a link, and they answer different questions:

- **`cross_links[].xl_type`** (unchanged): what joins the peptides. `inter` is a crosslink between two peptides, `intra` a looplink within one peptide, `dead-end` a monolink (one reacted end), `noncovalent` two peptides associated without a link. Open question 9 proposes the plainer values `crosslink`, `looplink`, `monolink` and `noncovalent`.
- **`link_scope`** (new): whether the two peptides come from one protein or two. `intra_protein` or `inter_protein`, the words of most XL-MS software (pLink `Protein_Type`, XlinkX `Crosslink Type`, Scout, MaxLynx). The CV names the same split self and heteromeric (`MS:1003343`). A looplink is always `intra_protein`; a monolink has no partner and gets null; a noncovalent pair is classified like a crosslink. Two peptides from distinct sequences of one proteoform family (a gene's canonical entry, its isoforms, variants and chains, as in [proteoforms](proteoform.md)) count as `intra_protein`, because their peptides usually cannot tell those sequences apart. The CV text would call them heteromeric; a producer that follows it says so in `crosslink_fdr`.

Two flags mark links that need care:

- **`possible_homomer`** (new, bool): an `intra_protein` crosslink (`xl_type` `inter`) whose two peptides are the same or overlap in the sequence. One chain cannot hold the same stretch twice unless the sequence repeats, so such a link points to a homo-oligomer (xiNET's homomultimeric link). Null for looplinks, monolinks and `inter_protein` links.
- **`partner_neighbor`** (new, bool): a peptide of the CSM has a neighbor (Lin, Plubell, Keich & Noble 2021): a different peptide of near-equal mass and similar fragment spectrum, in a protein the CSM does not list. The partner may then be a homolog of the true one, an error that reversed or shuffled decoys do not model; on a proteome-sized database such CSMs make up most of the false `inter_protein` CSMs a rescorer accepts. An I/L-identical copy is not a neighbor: its protein is listed in the peptide's `protein_accessions`.

## Decoys

- A **peptide** is a decoy when every protein it maps to is a decoy; a peptide that maps to a target and a decoy is a target.
- A **CSM** has `crosslink_decoy_class` `target_target`, `target_decoy` or `decoy_decoy` from its two peptides, and `is_decoy` is true for the latter two.
- A **residue pair** or **protein pair** takes the class of its best CSM, the CSM whose score ranks it.

The mzIdentML converter marks a CSM a decoy when any peptide evidence is a decoy (`crosslinks.py`, `_is_decoy`); under this rule it would change for a peptide that maps to both a target and a decoy.

## psm: crosslink columns

All nullable, and null on rows that are not crosslink spectrum matches. Crosslink PSMs carry provided ids, because the schema composite names only the alpha peptide: the file footer states `identity_source` `provided` or `mixed` ([Row Identity](identity.md)).

| Field | Description | Type |
|-------|-------------|------|
| `link_scope` | `intra_protein` or `inter_protein` (above) | string, null |
| `possible_homomer` | Overlapping peptides of one protein (above) | bool, null |
| `partner_neighbor` | A peptide has a neighbor in a protein the CSM does not list (above) | bool, null |
| `xl_peptides` | Per-peptide evidence, alpha first as in `peptidoform`; one entry for a looplink | array[xl_peptide], null |
| `supporting_spectra` | Further spectra of the CSM's peptides: MS3 scans of a cleavable linker's stubs, or ETD/EThcD spectra of the same precursor | array[xl_spectrum], null |
| `residue_pair_ids` | `residue_pair.residue_pair_id` of the residue pairs this CSM supports (`MS:1003344`) | array[int64], null |
| `fdr_estimate_id` | `crosslink_fdr` row of the CSM q-value | int32, null |

The CSM q-value is the PSM's `qvalue`, and the PEP of the whole identification (parent and further spectra together, `MS:1003336`) its `posterior_error_probability`. Scores go in `additional_scores` by name; `crosslink_fdr.score_name` says which one a q-value ranked by.

The link positions stay where they are: `cross_links[]` with `donor_position` and `acceptor_position` (a looplink has both on one peptide), `partner_sequence` and `partner_peptidoform`. `xl_peptides` holds only what those columns lack:

| `xl_peptide` field | Type | Description |
|--------------------|------|-------------|
| `protein_accessions` | array[string] | Proteins the peptide maps to (for the alpha peptide, equal to `protein_accessions`) |
| `protein_starts` | array[int32] | 1-based start of the peptide in each of those proteins |
| `is_decoy` | bool, null | Decoy peptide (rule above); null when the producer does not say |
| `stub` | string, null | Linker stub observed on this peptide (`MS:1003346`), or `intact` |
| `score` | float64, null | The peptide's own score |

With `xl_peptides[1].protein_accessions` and `protein_starts`, `partner_protein_accessions` and `partner_protein_starts` repeat a value; one of the two would be deprecated.

| `xl_spectrum` field | Type | Description |
|---------------------|------|-------------|
| `peptide_index` | int16 | Index into `xl_peptides` of the peptide the spectrum identifies |
| `scan` | int32 | Its scan number, one value (`psm.scan` is a list because a PSM may merge scans) |
| `ms_level` | int16 | 2 or 3 |
| `peptidoform` | string, null | The peptide identified in the spectrum, with its stub |
| `score` | float64, null | Score of that identification |

A spectrum from another run is not supported; mzIdentML allows it, but no producer we know writes one.

## residue_pair

One row per **linked residue pair**, written as protein and position on each side. When the evidence cannot tell several sites apart (a peptide shared by proteins, or several linkable residues), the row is the set of those sites, and that set is the unit the q-value counts: a set is false when no pair of its sites is linked in the sample. Two overlapping sets (`A:12 — B:5` and `A:12, C:12 — B:5`) are two rows.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `residue_pair_id` | Opaque int64 identity, hashed from the key below | int64 | yes |
| `protein_a`, `position_a`, `residue_a` | Side A: protein, position and amino acid of its first site | string, int32, string | yes |
| `protein_b`, `position_b`, `residue_b` | Side B, likewise | string, int32, string | yes |
| `gene_a`, `gene_b` | Gene names of the two proteins | string, null | no |
| `sites_a`, `sites_b` | Every site of each side, as (`protein`, `position`), sorted; the side whose first site sorts first is A | array[struct] | yes |
| `is_ambiguous` | A side has more than one site | bool | yes |
| `link_scope`, `possible_homomer`, `partner_neighbor` | Of the best CSM | string, bool, bool | yes |
| `linkers` | Reagents of the supporting CSMs | array[string] | yes |
| `grouped_runs` | Runs whose CSMs support the pair | array[string] | yes |
| `csm_count` | Supporting CSMs | int32 | yes |
| `best_psm_id` | Best CSM (`psm.psm_id`) | int64, null | no |
| `is_decoy`, `decoy_class` | Of the best CSM; null when the producer exports targets only | bool, string, null | no |
| `score` | Score the q-value ranked by (the best CSM's), full precision | float64, null | no |
| `qvalue` | Residue-pair q-value (`MS:1002677`) | float64, null | no |
| `fdr_estimate_id` | `crosslink_fdr` row of the q-value | int32, null | no |
| `pg_id_a`, `pg_id_b` | Identification `pg` row holding each side, as `site.pg_id` | int64, null | no |
| `additional_scores`, `cv_params` | As elsewhere | list | no |

`position` follows [site](site.md): 1-based residue, 0 for the protein N-terminus, length + 1 for the C-terminus, so a terminal link needs no extra field; `residue` is then the terminal residue. A link at a peptide N-terminus that is not a protein N-terminus (crosslinking after digestion) has the position of that first residue.

The **key** of `residue_pair_id` is `sites_a` and `sites_b`, then `grouped_runs` sorted, serialized as in [Row Identity](identity.md). The linker and the isotope label are not part of it: the light and heavy forms of one residue pair, and links made by two reagents, are one row, and `linkers` lists the reagents. Their separate intensities belong in a quantitative view (see [Possible follow-up views](#possible-follow-up-views)).

## protein_pair

One row per protein pair. An `inter_protein` pair is two proteins; an `intra_protein` link gives a pair of one protein with itself, a candidate homo-oligomer (a strong one when `possible_homomer`). mzIdentML 1.3 also reports such self pairs (`MS:1002676`). A producer that controls the error of `inter_protein` pairs only leaves the self pairs out.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `protein_pair_id` | Opaque int64 identity, hashed from `proteins_a`, `proteins_b` and `grouped_runs` | int64 | yes |
| `protein_a`, `protein_b` | The first protein of each side; equal for a self pair | string | yes |
| `gene_a`, `gene_b` | Gene names | string, null | no |
| `proteins_a`, `proteins_b` | Proteins each side cannot be told apart from (shared peptides) | array[string] | yes |
| `link_scope` | `intra_protein` for a self pair, else `inter_protein` | string | yes |
| `grouped_runs` | Runs whose CSMs support the pair | array[string] | yes |
| `residue_pair_count`, `csm_count` | Supporting residue pairs and CSMs | int32 | yes |
| `best_psm_id` | Best supporting CSM | int64, null | no |
| `is_decoy`, `decoy_class` | Of the best CSM | bool, string, null | no |
| `score` | Score the q-value ranked by | float64, null | no |
| `qvalue` | Protein-pair q-value (`MS:1002676`) | float64, null | no |
| `fdr_estimate_id` | `crosslink_fdr` row | int32, null | no |
| `pg_id_a`, `pg_id_b` | As in `residue_pair` | int64, null | no |
| `additional_scores`, `cv_params` | As elsewhere | list | no |

## crosslink_fdr: how each q-value was estimated

The FDR of crosslinks is usually estimated apart for `intra_protein` and `inter_protein` links, and sometimes apart for links with a neighbor partner, so a q-value means little without saying over which rows it was estimated. One row per such population; the CSM, residue-pair and protein-pair q-values each point to theirs with `fdr_estimate_id`.

| Field | Description | Type |
|-------|-------------|------|
| `fdr_estimate_id` | Identity of the population | int32 |
| `level` | `csm`, `peptide_pair`, `residue_pair` or `protein_pair` | string |
| `link_scope` | The link type estimated, or null when both were estimated together | string, null |
| `partner_neighbor` | The population is the CSMs with a neighbor partner (true), those without (false), or both (null) | bool, null |
| `score_name` | Score ranked by, as named in `additional_scores` | string |
| `false_count` | How false targets are counted: `td_minus_dd` (TD − DD) or `td_plus_dd` (TD + DD) | string |
| `n_target_target`, `n_target_decoy`, `n_decoy_decoy` | Rows of each class in the population | int32 |
| `runs` | Runs of the population | array[string] |

The FDR at a rank is the false-target count over the rows at that rank or better, divided by their TT count, capped at 1; the q-value is the lowest FDR at that rank or any worse rank, so it never falls as the score gets worse, whichever way the score runs. `td_minus_dd` is the usual count. `td_plus_dd` is for links whose two peptides have one sequence (homotypic links): the decoy of such a link is a decoy–decoy pair of one sequence, so there the DD rows are the decoys that stand for false targets, not a correction for TD rows counted twice, and they are added rather than subtracted. It is not a looplink estimator: looplinks and noncovalent pairs get a null q-value unless a producer records a population for them. andes adds a pseudo-count of 1 and a floor of 2 to the false count; a producer that does likewise says so in its documentation.

When one FDR is estimated over several projects together, `crosslink_project_fdr` gives, for each project, the FDR of its own CSMs at the q-value cut of the combined estimate, to show whether combining moved errors into a weaker project. It is descriptive, not a q-value: `max(false, 0) / TT` over that project's CSMs at or below the cut.

| Field | Description | Type |
|-------|-------------|------|
| `fdr_estimate_id` | The combined CSM population | int32 |
| `project` | Project, as the producer names its inputs | string |
| `q_cut` | CSM q-value cut of the combined estimate | float64 |
| `n_target_target`, `n_target_decoy`, `n_decoy_decoy` | That project's CSMs at or below the cut | int32 |
| `project_fdr` | The project's own FDR at the cut | float64 |

## run: linkers

`run.linkers` lists the reagents searched. The SDRF `comment[cross-linker]` gives the name and XLMOD accession; the other fields come from XLMOD or the search settings, and are null when neither has them.

| Field | Type | Description |
|-------|------|-------------|
| `linker_name` | string | Reagent name, as `cross_links[].linker_name` |
| `linker_accession` | string, null | XLMOD accession |
| `linker_mass` | float64, null | Mass of the intact link (Da) |
| `cleavable` | bool, null | MS-cleavable |
| `sites` | string, null | Residues each end reacts with (`K, S, T, Y, Protein N-term`) |
| `stubs` | array[{name, mass}], null | Stub masses of a cleavable linker (`MS:1003390`) |
| `label` | string, null | Isotope label of a light/heavy reagent pair |

## Relation to mzIdentML 1.3

mzIdentML 1.3 is the reference for what a crosslinking result contains; QPX keeps the subset a reader uses, in plain columns instead of encoded strings.

| QPX | mzIdentML 1.3 |
|-----|---------------|
| CSM (`psm` row) | Two SIIs sharing `MS:1002511`; donor/acceptor modifications `MS:1002509` / `MS:1002510` |
| Looplink (`xl_type` `intra`) | One SII with `MS:1003329` |
| `supporting_spectra` | Child SIIs grouped by `MS:1003332`; the group's PEP `MS:1003336` is `psm.posterior_error_probability` |
| CSM q-value (`psm.qvalue`) | `MS:1003337` |
| `peptide_pair_qvalue` | `MS:1003338`; passes threshold `MS:1003339` |
| `psm.residue_pair_ids` | `MS:1003344` (residue-pair ref) |
| `residue_pair.qvalue` | `MS:1002677`, passes threshold `MS:1003340`: encoded strings (`ID.a\|b:position:score:pass`) on ordinary protein groups |
| `protein_pair.qvalue` | `MS:1002676`, passes threshold `MS:1003341`, encoded the same way; self pairs included |
| `link_scope` | `MS:1003343`: self (`intra_protein`) and heteromeric (`inter_protein`) estimated apart |
| `run.linkers` | `SearchModification` with donor/acceptor terms, XLMOD accession, `MS:1003390` |

The CV has no residue-pair or protein-pair q-value term: `MS:1002677` and `MS:1002676` are global-FDR terms (global as opposed to local FDR, `MS:1002350` / `MS:1002351`), and q-values are their own family (`MS:1002354`). Requesting q-value terms for the two levels from PSI would close that gap. The passes-threshold flags are left out: a reader applies a cut to `qvalue`. mzIdentML's donor/acceptor ordering rule (longer peptide, then mass, then alphabetical) is a fallback for symmetric reagents; QPX keeps the alpha/beta order of `peptidoform`.

The mzIdentML converter would parse the encoded residue-pair and protein-pair values into `residue_pair` and `protein_pair`, pairing the `.a` and `.b` sides by their ID; today `pg_adapter._float_value` returns null for them and they are lost.

## Possible follow-up views

Neither is part of this proposal; they are listed so the identity and keys above leave room for them.

- **Quantified crosslinks** (`crosslink_feature`): one row per residue pair per run and label, with its intensity, as `feature` does for peptidoforms. Quantitative XL-MS compares conformations between conditions, with light/heavy reagents or label-free intensities.
- **Crosslinks on structures** (`crosslink_structure`): one row per residue pair per structure model (PDB entry, AlphaFold or ModelArchive model, with version and chains), with the Cα–Cα distance, the linker's maximum span and whether the link fits it.

## Open questions

1. Crosslink columns in `psm`, or a separate `csm` view?
2. `crosslink_fdr` and `crosslink_project_fdr`: tables, or provenance parameters?
3. Ambiguity as lists (`sites_a`, `sites_b`) or one row per alternative with a group id, as `modification_candidates` does?
4. Per-project FDR in the format, or only in producer reports?
5. Calibration columns (entrapment or truth-set error) in production files, or only in benchmarks?
6. Which of `partner_protein_*` and `xl_peptides[].protein_*` to keep?
7. Schema, `file_type`, collection and merge registration, versioning, ontology and converter-coverage entries: one PR per view or one for all?
8. Are the follow-up views worth a proposal of their own?
9. `cross_links[].xl_type`: rename `inter` / `intra` / `dead-end` to `crosslink` / `looplink` / `monolink` at the next schema version? Most XL-MS software reads `inter` / `intra` at protein level, which is why `link_scope` is a separate column.
10. Heterobifunctional reagents (which end reacted on which peptide) and the reacting end per peptide: worth a field, or left to `peptidoform`?

## References

- psi-ms.obo 4.2.17 and the mzIdentML 1.3 crosslinking extension 1.0.0.
- Leitner et al., Structure 2020, https://doi.org/10.1016/j.str.2020.09.011
- Lenz et al., Nat Commun 2021, https://doi.org/10.1038/s41467-021-23666-z
- Fischer & Rappsilber, Anal Chem 2017, https://pmc.ncbi.nlm.nih.gov/articles/PMC5423704
- Lin, Plubell, Keich & Noble, J Proteome Res 2021 (neighbor peptides), https://pmc.ncbi.nlm.nih.gov/articles/PMC8489664
- xiFDR, https://xifdr-py.readthedocs.io/en/stable/how_it_works.html
- Scout, https://github.com/diogobor/Scout
- MaxLynx, Yılmaz et al., Anal Chem 2022, https://pmc.ncbi.nlm.nih.gov/articles/PMC8792900
- pLink 2 output, https://github.com/pFindStudio/pLink2/wiki/CSV-result
- Shortreed & Smith, proteoform families, https://pmc.ncbi.nlm.nih.gov/articles/PMC9676067/
