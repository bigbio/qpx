# Protein View

The protein view holds what the sequence database says of each protein a dataset names: its name, organism and sequence. Every other view keeps accessions only, so join this one on them. It is optional, and specification 1.3 added it.

File: `{PREFIX}.protein.parquet`. One row per accession.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `accession` | Protein accession as the other views name it | string | Yes |
| `entry_name` | Entry name of the database, such as `ALBU_HUMAN` | string, null | No |
| `description` | Protein name from the FASTA header | string, null | No |
| `gene_name` | Gene symbol from the header's `GN=` | string, null | No |
| `organism` | Organism from the header's `OS=` | string, null | No |
| `taxonomy_id` | NCBI taxonomy id from the header's `OX=` | int32, null | No |
| `length` | Sequence length in residues | int32, null | No |
| `sequence` | Amino acid sequence | string, null | No |

## Producers

- `qpxc transform protein-properties --fasta` and the `--fasta` of `qpxc convert` write it for the accessions of the `pg` and `feature` views that the FASTA has, with sequences. A value a producer recorded stays.
- `qpxc convert maxquant` writes it from `proteinGroups.txt`: the `Fasta headers` and `Sequence lengths` of each `Protein IDs` member, where each cell has one token per member. It has no sequences.

A header not in UniProt form yields only its description. Decoys stay out, since a FASTA's decoy entries are skipped and MaxQuant's `REV__` rows are dropped.

## The FASTA itself

The provenance view's `protein_properties_from_fasta` step records the file's `fasta`, `fasta_sha256` and, with `--fasta-uri`, the `fasta_uri` others can fetch it from. Use `sequence` when the file is not at hand.

## Peptide positions

`feature.pg_positions` holds where a peptide sits in its proteins, and `pepmap.pg_accessions` also has per-protein `start`, `end`, `pre` and `post`. The `start`, `end`, `pre` and `post` of `feature.pg_accessions` are reserved: no converter fills them, so do not read positions there.
