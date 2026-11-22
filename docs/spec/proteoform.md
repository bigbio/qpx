# Proteoforms

`proteoform.parquet` reports, for each protein **family**, the canonical entry and each **isoform**, **sequence variant** and **processed chain** the evidence supports. A family groups a canonical protein (e.g. a UniProt entry) with its alternative sequence entries, so alternative entries do not split the protein's evidence; each member is a separate claim.

!!! warning "Scope: sequence-entry members"
    A member is a sequence entry. The view does not represent proteoforms as combinations of modifications on one molecule, and [site](site.md) rows of a family do not establish that their modifications co-occur.

The YAML schema definition lives at `qpx/core/data/schemas/proteoform.yaml` in the source repository.

## Use Cases

- **Proteogenomics**: report which sequence variants of a protein are observed, with error control over the reported members.
- **Isoform and chain evidence**: separate the evidence specific to an isoform or a processed chain from the evidence shared with the canonical form.
- **Proteoform quantification**: quantify a member from the peptides that identify it, per run or channel.

## Error control: `proteoform_qvalue`

A member claim is **false** when the member is not present in the sample. `proteoform_qvalue` is the estimated proportion of false members among the members of the view at least as well supported, estimated **at the member level**: each member counted once, for example by target-decoy or knockoff competition between members. The method and the population it was estimated over belong in [provenance](provenance.md).

A PSM-level q-value does not control error over members, even the q-value of a member's best PSM: many true PSMs collapse onto few true members while false PSMs each make a false member. A producer that has only PSM-level estimates leaves `proteoform_qvalue` null and keeps them in `additional_scores` (e.g. `best_psm_qvalue`). `best_psm_id` names a supporting PSM, not the source of the q-value.

## Evidence

- `specific_peptides` counts peptides that no other member of the **family** explains. Specificity is within the family only: such a peptide can also map to another family. Whether it counts for this family is decided by the protein inference behind `pg_id`, and the member's `intensity` uses only specific peptides that the inference assigns to the member's protein group.
- `evidence` is `specific` when the member has peptides specific to it within the family, and `shared` when only peptides shared within the family support it. A `shared` row establishes support for the **family**, not the presence of this member — a canonical row from shared-only evidence does not show that the canonical sequence itself is present.

## Schema

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `proteoform_id` | Opaque int64 identity, derived from the identity composite (`accession`, `grouped_runs`, `label`) | int64 | yes |
| `family` | Canonical entry shared by all members of the family | string | yes |
| `accession` | This member's accession | string | yes |
| `kind` | `canonical`, `isoform`, `variant` or `chain` | string | yes |
| `description` | What distinguishes the member (e.g. `p.Arg175His`, an isoform identifier, a chain range) | string, null | no |
| `edits` | For a variant, the changes from the canonical sequence (see [Sequence Edits](sequence-edits.md)); optional column | list\<sequence_edit\>, null | no |
| `pg_id` | Protein group row holding the family, for this row's label and runs (see below) | int64, null | no |
| `grouped_runs` | Runs grouped into this quantification unit | list\<string\> | yes |
| `source` | Database or annotation the member comes from | string, null | no |
| `prior_evidence` | Evidence class of the member before the search (e.g. `curated`, `catalogue`, `sample`) | string, null | no |
| `specific_peptides` | Peptides no other member of the family explains | int32, null | no |
| `shared_peptides` | Peptides shared with other members of the family | int32, null | no |
| `psm_count` | PSMs at their cutoff that support the member | int32, null | no |
| `evidence` | `specific` or `shared` (above) | string, null | no |
| `proteoform_qvalue` | Member-level q-value (above); null without a member-level estimate | float64, null | no |
| `best_psm_id` | A supporting PSM, the best one (`psm.psm_id`) | int64, null | no |
| `is_decoy` | Decoy or knockoff member | bool | yes |
| `label` | Channel/label of the quantity; null for identification only | string, null | no |
| `intensity` | Quantity from the member's peptides for this label | float32, null | no |
| `additional_intensities` | Other quantities for this label | list\<additional_intensity\> | no |
| `additional_scores` | Additional scores | list\<score\> | no |
| `cv_params` | CV parameters | list\<cv_param\> | no |

!!! note "One row per label, one pg row per label"
    Like [pg](pg.md), the view is flattened: a quantified member has one row per label, and an identification-only member has one row with `label` and `intensity` null. Writers accept the natural `intensities: list<{label, intensity}>` per member and unroll it. `pg_id` names one pg row, so each label row links to the pg row of its own label: give `pg_id` on each `intensities` entry. `ProteoformWriter` refuses a record-level `pg_id` on a member with several labels.

## Validation

Under `qpxc validate` (strict):

- `proteoform_qvalue` is finite and in [0, 1];
- `kind` is one of `canonical`, `isoform`, `variant`, `chain`, and `evidence` one of `specific`, `shared` (or null);
- a populated `pg_id` resolves to a pg row with the same `label` and `grouped_runs`;
- that pg row's `pg_accessions` hold the row's `family`;
- a populated `best_psm_id` resolves to a `psm.psm_id` when the dataset has a psm view.

`source` and `prior_evidence` are open vocabularies: producers name their own sources and evidence classes.

## Python API

```python
import qpx

ds = qpx.open_dataset("PXD014414/")

# Variants at a member-level q-value of 0.05 (members without one are left out)
variants = ds.proteoform.targets_only().by_kind("variant").confident(0.05)
print(variants.to_df()[["family", "accession", "description", "proteoform_qvalue"]])

# All members of one family
print(ds.proteoform.by_family("P04637").to_df())
```

## Relationship to Other Views

- `pg_id` links the family to its protein group in [pg](pg.md); the pg q-value controls the protein group, the proteoform q-value the member's own claim.
- `best_psm_id` links to a supporting [psm](psm.md) row.
- [site](site.md) rows of the same `family` report modification sites.
