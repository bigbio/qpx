# Row Identity

`psm_id`, `feature_id` and `pg_id` are opaque signed 64-bit primary keys. A producer may supply them; otherwise the writer hashes them from the file's `identity_composite`. This page states what an id means and how to reproduce one.

## Scope

An id is meaningful **within one file** of one analysis, and only for the rows of that view. It says nothing about the analysis it came from:

- Two analyses of the same raw files share ids. DIA-NN 1.8.1 and 2.5.1 on one project give 97.6 % equal `feature_id` values, with intensities a factor of 13 apart.
- The same input converted under two accessions gets identical ids.

So do not join across datasets on an id or on a natural key (`peptidoform`, `charge`, `run_file_name`) alone. `dataset.analysis_id` names the analysis: converters and `merge()` give each one a new value. `DatasetCollection` and `merge()` keep datasets apart with `source_dataset`, and `merge()` namespaces the ids of each source.

## Footer keys

| Key | Meaning |
|-----|---------|
| `identity_composite` | The columns the ids are hashed from, in hashing order |
| `primary_key` | The id column |
| `identity_hash_version` | Version of the bytes and the hash below; currently `1` |
| `identity_source` | `derived`: the writer hashed every id from the composite. `provided`: every id came with its row, from a producer or an earlier file. `mixed`: some of each. Absent before specification 1.3 |

`qpxc validate --check-identity` re-derives the ids of a `derived` file and reports those that differ. Files of one view that declare different composites or hash versions get a `mixed_identity` warning, because equal rows do not have equal ids.

## Composites

| View | Producer | `identity_composite` |
|------|----------|----------------------|
| `psm` | schema default | `run_file_name, scan, peptidoform, charge` |
| `feature` | schema default | `peptidoform, charge, run_file_name, rt` |
| `feature` | DIA-NN, Spectronaut, CDAP | `peptidoform, charge, run_file_name` |
| `feature` | MaxQuant | `peptidoform, charge, run_file_name, rt_start, rt_stop, scan` |
| `feature` | quantms MSstats | `run_file_name, peptidoform, charge, rt, scan` |
| `feature` | FragPipe `combined_ion` | `quantification_unit_id, peptidoform, charge, compensation_voltage` |
| `feature` | `qpxc convert openms-consensus` | `peptidoform, charge, run_file_name, rt, scan, observed_mz, consensus_rt` |
| `pg` | schema default | `pg_accessions, grouped_runs, label` |

OpenMS `-out_qpx` files declare their own composite in the footer, and `qpxc convert openms` re-derives ids with that composite. The two OpenMS routes declare different feature composites, so their `feature_id` values differ for the same data.

`qpxc convert mzidentml` supplies its `psm_id` values (`identity_source` `provided`): the hash of `[run_file_name, spectrum_identification_item_ids]`, because a spectrum can hold several hits of one peptidoform. Its footer still names the schema default, which does not reproduce them. On a cross-linked PSM the validator reads `crosslinked_peptidoform` in place of `peptidoform`, which names only the alpha peptide.

## JSON and URLs

Ids are signed 64-bit integers, and almost all of them lie beyond ±2^53, where a JavaScript number can no longer hold them exactly. Rounded ids do not collide; lookups simply fail. In JSON and in URLs an id is therefore a decimal string. `qpxc query --output-format json` writes every `*_id` and `*_ids` column this way; an alias in your own SQL (`SELECT feature_id AS f`) is an ordinary integer column, so cast it with `CAST(feature_id AS VARCHAR)`.

## How an id is derived

1. Take the composite's values in composite order, **as persisted**: a `float32` column contributes the value it reads back, so `1234.56` hashes as `1234.56005859375`. Hashing the raw double gives another id.
2. In `pg_accessions` and `grouped_runs`, drop duplicates and sort by each element's JSON encoding. Other lists, such as `scan`, keep their order.
3. Encode the values as one JSON array in UTF-8: no whitespace, non-ASCII characters as `\uXXXX`, strings escaped as JSON does, null as `null`, and a float as its shortest round-trip decimal of the 64-bit value, written like Python's `repr` (`1234.56005859375`, `1e-05`).
4. Compute BLAKE2b with a digest size of 8 bytes. The digest size is a parameter of the hash, so truncating a longer BLAKE2b digest gives a different result.
5. Read those 8 bytes as a big-endian signed integer.

In Python this is `qpx.core.data.identity.derive_composite_id(row, composite)`. The JDK, Node's `crypto` and DuckDB SQL have no BLAKE2b with this digest size; libsodium's generic hash and Bouncy Castle's `Blake2bDigest(64)` do.

### Test vectors

| `identity_composite` | Bytes | Id |
|----------------------|-------|----|
| `run_file_name, scan, peptidoform, charge` | `["run_1",[1234],"PEPTIDEK",2]` | `-5127882623363935776` |
| `peptidoform, charge, run_file_name, rt` | `["[UNIMOD:1]-AC[UNIMOD:4]DEK",3,"a,b",1234.56005859375]` | `6417640111734626350` |
| `pg_accessions, grouped_runs, label` | `[["P04114","Q9Y6K9"],["f1","f2"],"TMT126"]` | `-5459222103841755848` |
| `pg_accessions, grouped_runs, label` | `[["P04114"],["r"],null]` | `-4507348376872481825` |

## Changes to the derivation

Any change to how ids are derived moves the specification version, so that `qpx_version` tells such files apart, and a change to the bytes or the hash also moves `identity_hash_version`. Specification 1.3 changes the PSM composite to the order agreed in [#229](https://github.com/bigbio/qpx/issues/229); see [Versioning](versioning.md).
