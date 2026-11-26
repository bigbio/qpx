# Crosslinks (proposed)

!!! warning "Proposal, not part of the specification"
    This page proposes how QPX holds crosslinking (XL-MS) results above the PSM level. Nothing on it is implemented: no schema, writer or validator check exists yet. Discussion is in issue #394; field names and placement may change.

A crosslinking result is reported at four levels, the levels of the literature and the CV (Fischer & Rappsilber 2017; Leitner et al. 2020), each with its own q-value:

| Level | Claim | Where | q-value |
|-------|-------|-------|---------|
| `csm` (crosslink spectrum match) | This spectrum is two linked peptides, or one looped peptide | [`psm`](psm.md) row | `psm.qvalue` |
| `peptide_pair` | These two peptides are linked | `psm` row | `peptide_pair_qvalue` in `additional_scores` |
| `residue_pair` | These two residues are linked in the sample | `residue_pair` (new) | `residue_pair_qvalue` |
| `protein_pair` | These two proteins are linked: two proteins that interact, or one protein with itself | `protein_pair` (new) | `protein_pair_qvalue` |

The page adds as few names as it can. The two views follow the layout of [site](site.md) (`protein_accession`, `position`, `amino_acid`, `family`, `grouped_runs`, `label`, `psm_count`, `best_psm_id`, `pg_id`); the q-values follow the typed q-value names (`site_qvalue`, `pg_qvalue`, `proteoform_qvalue`); the peptides keep the existing `psm` columns; the reagents go in `run.modification_parameters` and the FDR settings in [provenance](provenance.md). Search-engine internals (per-ion scores, candidate lists, encoded mzIdentML values) stay in the producer's own files.

## Link classes

Two columns describe a link, and they answer different questions:

- **`cross_links[].xl_type`** (unchanged): what joins the peptides. `inter` is a crosslink between two peptides, `intra` a looplink within one peptide, `dead-end` a monolink, `noncovalent` two peptides associated without a link. Open question 6 proposes plainer values.
- **`crosslink_protein_class`** (new, named like the existing `crosslink_decoy_class`): `self` when both peptides come from one protein sequence, `heteromeric` when they come from distinct sequences, the split of `MS:1003343` (Leitner et al. 2020; Lenz et al. 2021). A looplink is `self`; a monolink has no partner and gets null; a noncovalent pair is classified like a crosslink.

Two peptides from distinct entries of one protein `family` (the canonical entry, its isoforms, variants and chains, as in [site](site.md) and [proteoforms](proteoform.md)) are `self` in QPX, because their peptides usually cannot tell those entries apart. This deviates from `MS:1003343`, whose text calls a link between distinct sequences heteromeric; a producer that follows the CV text records it in [provenance](#how-the-q-values-were-estimated).

**Overlapping peptides.** A `self` crosslink whose two peptides are the same or overlap in the protein cannot come from one chain unless the sequence repeats, so it points to a homo-oligomer (xiNET's homomultimeric link). QPX stores no flag for it: it follows from columns the PSM already has. The query compares the first protein and start of each side, so it misses an overlap in a shared protein listed later, or at a repeat other than the first (open question 9).

```sql
SELECT *
FROM (SELECT *, unnest(cross_links) AS link FROM psm)
WHERE link.xl_type = 'inter'
  AND protein_accessions[1] = partner_protein_accessions[1]
  AND protein_starts[1] < partner_protein_starts[1] + length(link.partner_sequence)
  AND partner_protein_starts[1] < protein_starts[1] + length(sequence)
```

**Neighbor partners.** A peptide may have a neighbor (Lin, Plubell, Keich & Noble 2021): a different peptide of near-equal mass and similar fragment spectrum, in a protein the CSM does not list. Its partner may then be a homolog of the true one, an error that reversed or shuffled decoys do not model. A producer that marks such CSMs, or estimates their FDR apart, says so in `cv_params` and in provenance; QPX adds no column. An I/L-identical copy is not a neighbor but ordinary protein ambiguity: its protein belongs in `protein_accessions`.

## Decoys

The PSM rules apply unchanged: a peptide that maps to a target and a decoy is a target (`psm.is_decoy`), and a CSM has `crosslink_decoy_class` `target_target`, `target_decoy` or `decoy_decoy`, with `is_decoy` true for the latter two. A residue pair or protein pair takes `is_decoy` and `crosslink_decoy_class` from its best CSM, the CSM whose score ranks it.

A decoy side names a decoy accession: a decoy never shares an accession with its target ([Producer Guide](producer-guide.md)), and the mzIdentML converter adds `DECOY_` where a producer gives it one (xiFDR does). That is what keeps a decoy pair's identity apart from its target pair's, since the identity composite holds no decoy class: with the targets' accessions, a target and a decoy protein pair share a key in 7 of the 10 xiFDR files we converted from PRIDE.

The mzIdentML converter marks a CSM a decoy when any peptide evidence is a decoy (`crosslinks.py`, `_is_decoy`); under the PSM rule it would change for a peptide that maps to both a target and a decoy. In the PRIDE files we converted, that happens only where decoys are known by accession prefix (Kojak: the alpha peptide of 168 PSMs and the beta peptide of 713 CSMs).

## psm: new columns

The peptides of a CSM are already described: `cross_links[]` holds `xl_type`, `partner_sequence`, `partner_peptidoform`, `donor_position`, `acceptor_position` and the `linker_*` fields; `protein_accessions` and `protein_starts` place the alpha peptide, `partner_protein_accessions` and `partner_protein_starts` the beta peptide; `crosslinked_peptidoform` names the linked species. The CSM q-value is `qvalue`, and the PEP of the whole identification (parent and child spectra together, `MS:1003336`) is `posterior_error_probability`. Per-peptide scores go in `additional_scores` as `alpha_score` and `beta_score`.

Three columns are new, all nullable and null on rows that are not crosslink spectrum matches:

| Field | Description | Type |
|-------|-------------|------|
| `crosslink_protein_class` | `self` or `heteromeric` (above) | string, null |
| `child_spectra` | Further spectra of the CSM's peptides: MS3 scans of a cleavable linker's stubs, or ETD/EThcD spectra of the same precursor (`MS:1003332`) | array[child_spectrum], null |
| `residue_pair_ids` | `residue_pair.residue_pair_id` of the residue pairs this CSM supports (`MS:1003344`) | array[int64], null |

| `child_spectrum` field | Type | Description |
|------------------------|------|-------------|
| `scan` | array[int32] | Scan identifier, as `psm.scan` |
| `ms_level` | int16 | 2 or 3 |
| `peptide` | string | `alpha` or `beta`: the peptide the spectrum identifies |
| `peptidoform` | string, null | That peptide as identified in the spectrum, with its linker stub as a modification (`MS:1003346`) |

Crosslink PSMs carry provided ids, because the schema composite names only the alpha peptide: the file footer states `identity_source` `provided` or `mixed` ([Row Identity](identity.md)). A child spectrum from another run is not supported; mzIdentML allows it, but no producer we know writes one.

## residue_pair

One row per **linked residue pair** and label, laid out as a [site](site.md) row with an A and a B side. When the evidence cannot tell several sites apart (a peptide shared by proteins, or several linkable residues), the other sites are the side's `alternatives` (other sites, unlike `site.alternatives`, which lists other modifications of one residue), and the set is the unit the q-value counts: a set is false when no pair of its sites is linked in the sample. Two overlapping sets (`A:12 — B:5` and `A:12, C:12 — B:5`) are two rows. A peptide that occurs more than once in one protein (polyubiquitin, filaggrin) has one `psm.protein_starts` entry per accession today, so its other sites cannot be derived from the psm view (open question 9). A noncovalent pair links no residues, so it has no `residue_pair` row.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `residue_pair_id` | Opaque int64 identity, derived from the identity composite (both sides with their alternatives, `grouped_runs`, `label`) | int64 | yes |
| `protein_accession_a`, `position_a`, `amino_acid_a` | Side A: protein, position and residue of its first site | string, int32, string null | yes |
| `protein_accession_b`, `position_b`, `amino_acid_b` | Side B, likewise; the side whose site sorts first is A, and between equal sites the side whose sorted `alternatives` sort first (element by element, null first) | string, int32, string null | yes |
| `family_a`, `family_b` | Protein family of each side, as `site.family` | string, null | no |
| `alternatives_a`, `alternatives_b` | The other sites each side may be, as (`protein_accession`, `position`), sorted | array[struct], null | no |
| `crosslink_protein_class` | `self` or `heteromeric` | string | yes |
| `linker_names` | Reagents of the supporting CSMs, as `cross_links[].linker_name` | array[string] | yes |
| `grouped_runs` | Runs whose CSMs support the pair | array[string] | yes |
| `label` | Channel of the reagent or sample, null without one; part of the identity, as in `site` | string, null | no |
| `psm_count` | CSMs at the CSM cutoff that support the pair: the `csm_fdr` of provenance, or `psm.pass_threshold` where the producer reports no CSM q-value | int32, null | no |
| `best_psm_id` | Best of those CSMs (`psm.psm_id`) | int64, null | no |
| `residue_pair_qvalue` | q-value over the residue pairs (`MS:1002677`, see below) | float64, null | no |
| `is_decoy`, `crosslink_decoy_class` | Of the best CSM | bool, string null | yes |
| `pg_id_a`, `pg_id_b` | Protein group row holding each side, as `site.pg_id` | int64, null | no |
| `additional_scores`, `cv_params` | As elsewhere; the score the q-value ranked by is here under its name | list | no |

`position` follows [site](site.md): 1-based residue, 0 for the protein N-terminus, length + 1 for the C-terminus, so a terminal link needs no extra field. A link at a peptide N-terminus that is not a protein N-terminus (crosslinking after digestion) has the position of that first residue. Gene names come from the protein group rows through `pg_id_a` and `pg_id_b`.

A quantified residue pair adds `intensity`, as a quantified site does; it is not part of this proposal.

## protein_pair

One row per protein pair and label. A `heteromeric` pair is two proteins; a `self` pair is one protein, or one protein `family`, with itself. A self link may join two sites of one molecule or of two copies; a self pair points to a homo-oligomer only when its peptides overlap (above). mzIdentML 1.3 also reports such self pairs (`MS:1002676`). A producer that controls the error of heteromeric pairs only leaves the self pairs out.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `protein_pair_id` | Opaque int64 identity, derived from the identity composite (both sides with their alternatives, `grouped_runs`, `label`) | int64 | yes |
| `protein_accession_a`, `protein_accession_b` | The first protein of each side; for a self pair equal, or two entries of one `family`. The side whose protein sorts first is A, and between equal proteins the side whose sorted `alternatives` sort first, as in `residue_pair` | string | yes |
| `family_a`, `family_b` | Protein family of each side | string, null | no |
| `alternatives_a`, `alternatives_b` | The other proteins each side cannot be told apart from (shared peptides) | array[string], null | no |
| `crosslink_protein_class` | `self` or `heteromeric` | string | yes |
| `grouped_runs` | As in `residue_pair` | array[string] | yes |
| `label` | As in `residue_pair` | string, null | no |
| `residue_pair_count`, `psm_count` | Supporting residue pairs, and CSMs at the CSM cutoff as in `residue_pair` | int32, null | no |
| `best_psm_id` | Best supporting CSM | int64, null | no |
| `protein_pair_qvalue` | q-value over the protein pairs (`MS:1002676`, see below) | float64, null | no |
| `is_decoy`, `crosslink_decoy_class` | Of the best CSM | bool, string null | yes |
| `pg_id_a`, `pg_id_b` | As in `residue_pair` | int64, null | no |
| `additional_scores`, `cv_params` | As elsewhere | list | no |

## How the q-values were estimated

A crosslink q-value means little without saying how it was estimated: the FDR is usually estimated apart for self and heteromeric links (`MS:1003343`), and sometimes apart for other subsets. QPX records this where it records every FDR setting, in the `parameters` of the [provenance](provenance.md) step that estimated it, next to the existing `psm_fdr`, `peptide_fdr` and `protein_fdr`:

| Key | Value |
|-----|-------|
| `csm_fdr`, `peptide_pair_fdr`, `residue_pair_fdr`, `protein_pair_fdr` | The threshold applied at each level, if any |
| `fdr_score` | Score ranked by, as named in `additional_scores` |
| `fdr_false_count` | `td_minus_dd` (TD − DD) or `td_plus_dd` (TD + DD) |
| `fdr_separate_self_heteromeric` | `true` when the two classes were estimated apart (`MS:1003343`) |
| `fdr_subsets` | Any further subsets estimated apart, e.g. CSMs with a neighbor partner |

The FDR at a rank is the false-target count over the rows at that rank or better, divided by their target-target count, bounded to [0, 1] (TD − DD can be negative), and 1 where no target-target row is at that rank or better; the q-value is the lowest FDR at that rank or any worse rank, whichever way the score runs. `td_minus_dd` is the usual count. `td_plus_dd` is for links whose two peptides have one sequence: the decoy of such a link is a decoy–decoy pair of one sequence, so its DD rows stand for false targets and are added, not subtracted. It is not a looplink estimator; looplinks and noncovalent pairs get a null q-value unless a producer states an estimate for them. A producer that adjusts the count (andes never counts fewer false targets than DD rows, and adds a pseudo-count and a floor) describes it in the step's `config`.

When one FDR is estimated over several projects together, the error of each project's part at the combined cut is a producer report, not a QPX table.

## run: reagents

A crosslinker is a search modification, as in mzIdentML (`SearchModification` with its XLMOD accession), so `run.modification_parameters` holds it, filled from the SDRF `comment[cross-linker]`:

| `modification_parameters` field | For a crosslinker |
|---------------------------------|-------------------|
| `accession` | XLMOD accession (`XLMOD:02001`) |
| `name` | Reagent name, as `cross_links[].linker_name` |
| `fixed` | `false` |
| `position`, `target_amino_acid` | Where each end reacts (`K,S,T,Y`, `Protein N-term`) |
| `mass_delta` | Mass of the intact link, as `cross_links[].linker_mass` |

Stub masses of a cleavable reagent (`MS:1003390`) come from XLMOD by accession. A light/heavy reagent pair is described in `label_definitions`, as a SILAC label is. The `accession` field today names UNIMOD accessions; it would also accept XLMOD ones.

## Relation to mzIdentML 1.3

mzIdentML 1.3 is the reference for what a crosslinking result contains; QPX keeps the subset a reader uses, in plain columns instead of encoded strings.

| QPX | mzIdentML 1.3 |
|-----|---------------|
| CSM (`psm` row) | Two SIIs sharing `MS:1002511`; donor/acceptor modifications `MS:1002509` / `MS:1002510` |
| Looplink (`xl_type` `intra`) | One SII with `MS:1003329`; in mzIdentML 1.2, a donor and an acceptor on one Peptide |
| Noncovalent pair (`xl_type` `noncovalent`) | Two SIIs sharing `MS:1003331`; older xiFDR versions write a crosslink pair at location -1 instead |
| `child_spectra` | Child SIIs grouped by `MS:1003332`; the group's PEP `MS:1003336` is `psm.posterior_error_probability` |
| `psm.qvalue` of a CSM | `MS:1003337` (an FDR, see below) |
| `peptide_pair_qvalue` | `MS:1003338` (an FDR, see below); passes threshold `MS:1003339` |
| `psm.residue_pair_ids` | `MS:1003344` (residue-pair ref) |
| `residue_pair_qvalue` | `MS:1002677` (an FDR, see below), passes threshold `MS:1003340`: encoded strings (`ID.a\|b:position:score:pass`) on ordinary protein groups |
| `protein_pair_qvalue` | `MS:1002676` (an FDR, see below), passes threshold `MS:1003341`, encoded the same way; self pairs included |
| `crosslink_protein_class` | `MS:1003343`: self and heteromeric estimated apart |
| `run.modification_parameters` | `SearchModification` with donor/acceptor terms, XLMOD accession, `MS:1003390` |

The CV has no crosslink q-value term at any level: `MS:1003337`, `MS:1003338`, `MS:1002677` and `MS:1002676` are global-FDR terms (global as opposed to local FDR, `MS:1002350` / `MS:1002351`), and q-values are their own family (`MS:1002354`, for PSMs). A value under one of them is an FDR, which equals the q-value only when the FDR does not fall at worse ranks: FDRs of 0.10 and then 0.05 give the first row a q-value of 0.05. A converter therefore fills `psm.qvalue`, `peptide_pair_qvalue`, `residue_pair_qvalue` and `protein_pair_qvalue` only with a q-value, one the producer states or one derived from the ranked FDRs; otherwise it keeps the FDR in `additional_scores` under its own name and leaves the q-value null. Requesting q-value terms for the levels from PSI would close that gap. The passes-threshold flags of the pair levels are left out: a reader applies a cut to the q-value. `psm.pass_threshold` stays, since the CSM is often the only level with a cut: no file we converted from PRIDE (mzIdentML 1.1 and 1.2) has a CSM q-value. mzIdentML's donor/acceptor ordering rule (longer peptide, then mass, then alphabetical) is a fallback for symmetric reagents; QPX keeps the alpha/beta order of `peptidoform`.

The mzIdentML converter would parse the encoded residue-pair and protein-pair values into `residue_pair` and `protein_pair`, pairing the `.a` and `.b` sides by their ID; today `pg_adapter._float_value` returns null for them and they are lost. What those values hold, and how the converter should read them, is open question 11. xiFDR also writes some SIIs twice in one result, and the protein groups reference both copies, so `psm_count` and `best_psm_id` come from the psm rows, not from those references.

## Possible follow-up

**Crosslinks on structures** (`crosslink_structure`): one row per residue pair per structure model (PDB entry, AlphaFold or ModelArchive model, with version and chains), with the Cα–Cα distance, the linker's maximum span and whether the link fits it. Not part of this proposal; listed so the residue-pair identity leaves room for it.

## Open questions

1. Crosslink columns in `psm`, or a separate `csm` view?
2. Ambiguity as `alternatives_a` / `alternatives_b`, or one row per alternative with a group id, as `modification_candidates` does?
3. Calibration columns (entrapment or truth-set error) in production files, or only in benchmarks?
4. Schema, `file_type`, collection and merge registration, versioning, ontology and converter-coverage entries: one PR per view or one for all?
5. Is the structure view worth a proposal of its own?
6. `cross_links[].xl_type`: rename `inter` / `intra` / `dead-end` to `crosslink` / `looplink` / `monolink` at the next schema version? Most XL-MS software reads `inter` / `intra` at protein level (Mascot: intralinks are within one protein, looplinks within one peptide), which is why `crosslink_protein_class` is a separate column.
7. Heterobifunctional reagents (which end reacted on which peptide): worth a field, or left to `peptidoform`?
8. `run.modification_parameters.accession` accepting XLMOD accessions: enough, or does a reagent need fields a modification lacks?
9. `psm.protein_starts` holds one start per accession. Should a peptide that occurs more than once in a protein list each start, repeating the accession, so its residue pairs can be derived? In the 10 xiFDR files from PRIDE, QPX rebuilds xiFDR's residue pairs from the CSMs (`protein_starts` + position − 1), except 173, and 168 of those are in such proteins.
10. Monolinks: in the PRIDE files they are Unimod "… quenched monolink" modifications on the peptide (Kojak `UNIMOD:1878` and `UNIMOD:1879` on 20,772 PSMs, OpenPepXL `UNIMOD:1020` on 642), and the mzIdentML converter keeps them in `peptidoform` with no `cross_links` entry. Should converters write them as `xl_type` `dead-end`, or should monolinks stay modifications? Not every `Xlink:` Unimod entry is a monolink: only those named "… quenched monolink" are.
11. The encoded residue-pair and protein-pair values: xiFDR writes its FDR there, -1 when it estimated none, and other negative values too, and one protein-pair id can carry different values on different groups (9 in PXD062462). Should the converter derive q-values from them, which needs every ranked pair, failing pairs and decoys included, and map -1 to null? Which value does a conflicting id take?
12. Noncovalent pairs: does a pair of proteins seen only as a noncovalent association count as a `protein_pair` row, and with what q-value?

## References

- psi-ms.obo 4.2.17 and the mzIdentML 1.3 crosslinking extension 1.0.0.
- Leitner et al., Structure 2020, https://doi.org/10.1016/j.str.2020.09.011
- Lenz et al., Nat Commun 2021, https://doi.org/10.1038/s41467-021-23666-z
- Fischer & Rappsilber, Anal Chem 2017, https://pmc.ncbi.nlm.nih.gov/articles/PMC5423704
- Lin, Plubell, Keich & Noble, J Proteome Res 2021 (neighbor peptides), https://pmc.ncbi.nlm.nih.gov/articles/PMC8489664
- xiFDR, https://xifdr-py.readthedocs.io/en/stable/how_it_works.html
- Shortreed & Smith, proteoform families, https://pmc.ncbi.nlm.nih.gov/articles/PMC9676067/
