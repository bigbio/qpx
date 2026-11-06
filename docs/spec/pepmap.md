# Peptide-Protein Mapping

`pepmap.parquet` is a **deduplicated mapping** between peptidoforms and the proteins they belong to. Each row represents a unique peptidoform, with a nested list of protein accessions (including optional positional info).

The YAML schema definition lives at `qpx/core/data/schemas/pepmap.yaml` in the source repository.

## Use Cases

- **Peptide uniqueness analysis**: Identify peptides that map to a single protein vs. shared peptides.
- **Protein coverage**: Determine which regions of a protein are covered by identified peptides via `start`/`end` positions in the `pg_accessions` struct.
- **Space-efficient storage**: At scale (10M+ PSMs), storing protein mappings as a separate table avoids duplicating protein accession strings across millions of PSM rows.

## Schema

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `sequence` | Unmodified peptide amino acid sequence | string | yes |
| `peptidoform` | ProForma notation | string | yes |
| `pg_accessions` | Protein accessions with optional positions (start, end, pre, post) | list\<pg_protein\> | no |
| `sequence_edits` | Sequence changes the peptide carries relative to a reference protein it maps to (see [Sequence Edits](sequence-edits.md)); optional column | list\<peptide_sequence_edit\>, null | optional |
| `is_unique` | `true` if peptide maps to exactly one protein | bool | no |

The `pg_protein` struct contains:

| Field | Type | Description |
|-------|------|-------------|
| `accession` | string | UniProt protein accession |
| `start` | int32 | Start position in protein sequence (1-indexed) |
| `end` | int32 | End position in protein sequence (1-indexed) |
| `pre` | string | Flanking amino acid before the peptide |
| `post` | string | Flanking amino acid after the peptide |

!!! note "Primary key"
    The primary key is `peptidoform`. Each unique peptidoform appears exactly once.

## Python API

```python
import qpx

ds = qpx.open_dataset("PXD014414/")

# Access the mapping
mapping = ds.pepmap

# Filter by protein
brca1_peptides = mapping.by_protein("P38398")
print(brca1_peptides.to_df())

# Filter by peptide
proteins = mapping.by_peptide("PEPTIDEK")
print(proteins.to_df())

# Get only unique peptides
unique = mapping.unique_peptides()
print(f"Unique peptides: {unique.count()}")
```

## Relationship to Other Views

The `pepmap` complements the `pg_accessions` field on Feature views. `pepmap` holds the **evidence proteins**, all the proteins the peptide maps to; `pg_accessions` holds the **assigned group**, which for a razor assignment names fewer proteins than the evidence:

- **Feature `pg_accessions`**: A nested list of `pg_protein` structs on each row. Simple but duplicates strings at scale.
- **`pepmap.parquet`**: A deduplicated lookup table with positional info and uniqueness flags. More space-efficient for large datasets.

Both share the same `pg_protein` struct type from `types.yaml`.

## Writing a Peptide-Protein Map

```python
from qpx.writers import PepMapWriter

records = [
    {
        "sequence": "PEPTIDEK",
        "peptidoform": "PEPTIDEK",
        "pg_accessions": [
            {"accession": "P12345", "start": 100, "end": 107, "pre": "K", "post": "R"},
        ],
        "is_unique": True,
    },
]

with PepMapWriter("exp.pepmap.parquet") as w:
    w.write_batch(records)
```

## Related Pages

- [PSM View](psm.md) -- carries optional `protein_accessions` list
- [Feature View](feature.md) -- `pg_accessions` and `anchor_protein` fields
- [Protein Group View](pg.md) -- `pg_accessions` and protein group inference
