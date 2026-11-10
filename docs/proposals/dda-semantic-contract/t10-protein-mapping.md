# T10. Protein mapping on PSMs and features: accessions, groups, uniqueness, decoys and contaminants

[← T09](t09-intensities-missing-values.md) · [Overview](index.md) · [Strict profile](strict-profile.md) · [Annex overview](annex-overview.md) · [T11 →](t11-spectrum-precursor-fields.md)

!!! note "Annex: producer-neutral analysis"
    Part of the annex of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text, and they cover every producer. For OpenMS-written files the [strict profile](strict-profile.md) is the contract and wins where the two differ. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** `psm.protein_accessions` is the set of mapped entries (no inference, no decoy entries on target rows); `feature.pg_accessions` is one group's full membership or null; `anchor_protein` is its representative (MS:1002403), never an identity; a new `group_unique` marks group-specific peptides and `unique` derives from it; `is_decoy` follows one all-entries rule; markers stay verbatim and a nullable `contaminant` on psm and feature replaces parsing.

## Questions a reader has
- `psm.protein_accessions`: every mapped protein or a razor one? Ordered? Decoy and contaminant entries?
- `feature.pg_accessions`: membership or evidence proteins, strings or structs, a set, matching a pg row? What does a shared peptide carry, assigned or not?
- `unique`: single-protein group, group-specific or proteotypic? `anchor_protein`: leader, same-set anchor or representative; a key?
- `is_decoy` for target+decoy peptides and mixed groups? Contaminant PSMs and features? Does `pg.contaminant = false` mean "checked"?

## What happens today
**Spec.**
- `psm.md:57`: "all proteins that the peptide maps to"; order, decoys, duplicates open. `psm.md:289` points to feature and pg, which hold only the assigned group (`pepmap.md:61`). [#78](https://github.com/bigbio/qpx/issues/78): no inference.
- `pg_accessions` is `array[struct]` (`feature.md:99`) with "reserved" positions (`feature.yaml:151-153`); examples write strings that FeatureWriter rejects (`feature.md:198,261`) [[SPEC-03](findings-feature.md#spec-03)]; `producer-guide.md:120` fills positions.
- `feature.md:100` nulls `anchor_protein` only "when no protein mapping was performed". `feature.yaml:154-158` also nulls it for unassigned shared peptides ([#322](https://github.com/bigbio/qpx/issues/322)).
- `unique`: "a one-protein group that maps to no other group" (`feature.md:101`, `versioning.md:68-69`); `pepmap.md:21` `is_unique`: one protein; `pg.md:82`: unique to the group. `feature.md:117` plans the name `is_unique`.
- `pg.md:39` calls the anchor the "leading protein". `ontology.md:68,212,234` map it to MS:1001591 (same-set anchor). `pg.yaml:43`, `pg.md:66` and `data-model.md:159` forbid keying on it, yet `absolute.md:67`, `site.md:113` and `views.md:177-178,212` do.
- Target+decoy is target (`psm.md:45`, `feature.md:60`). Markers stay, "do not parse" (`pg.md:46-49`). No mixed-group or null-`contaminant` rule (`pg.md:94-95`); contaminants flagged on pg only (`producer-guide.md:34`).
- `pg_group_id` hashes "the sorted pg_accessions" (`feature.yaml:166`). `identity.md:51` dedups for `pg_id`; `qpx/writers/derived.py:101-112` does not.

| Producer | Behaviour | Findings |
|---|---|---|
| MaxQuant | psm: `Proteins`, blank on reverse hits (`qpx/converters/maxquant/psm_adapter.py:205-211`). With proteinGroups.txt: the `Protein group IDs` group (the razor group if shared), `unique` = one id and one member (`feature_adapter.py:319-348`). Without it: `Leading proteins` (one leader per group), and `unique` = one `Proteins` entry, i.e. proteotypic (`:350-366`). Evidence `Potential contaminant` is never read | [[MAXQUANT-22](findings-cross.md#maxquant-22)] [[MAXQUANT-28](findings-feature.md#maxquant-28)] [[MOKUME-02](findings-feature.md#mokume-02)] |
| FragPipe | `contam_` is stripped in all views (`qpx/converters/fragpipe/constants.py:42-58`); `canonical_accession` leaves `contam_sp\|P02768\|ALBU_HUMAN` whole (`qpx/converters/utils.py:69-73`). psm `is_decoy` comes from the razor `Protein` only (`psm_adapter.py:160`), the feature's from its best PSM (`feature_adapter.py:482`). `unique = len(Protein+Mapped Proteins)==1` (`:531`). Without `--pg-file`: lead only (`:529`) | [[FRAGPIPE-07](findings-cross.md#fragpipe-07)] [[FRAGPIPE-08](findings-psm.md#fragpipe-08)] [[FRAGPIPE-18](findings-feature.md#fragpipe-18)] [[MOKUME-M2](findings-feature.md#mokume-m2)] |
| Proteome Discoverer | psm: `Protein Accessions` (`psm_adapter.py:66`). Anchor = `masters[0]`, also when shared (`feature_adapter.py:105`); no Proteins.txt: merged masters (`:112`). `unique = len(mapped)==1` (`:136`), pinned on a 2-member group (`tests/converters/test_proteome_discoverer_converter.py:166`). `is_decoy` is always false. `pg.contaminant` is false if the column is absent (`pg_adapter.py:150`). Members are a union over `Protein Group IDs` (`pg_adapter.py:37-45`) | [[PD-20](findings-feature.md#pd-20)] [[PD-21](findings-feature.md#pd-21)] [[PD-22](findings-pg.md#pd-22)] [[PD-23](findings-cross.md#pd-23)] [[PD-32](findings-psm.md#pd-32)] |
| OpenMS consensusXML | psm/feature `is_decoy` only from `target_decoy`, false if absent (`feature_adapter.py:125-132`); pg falls back to DECOY/REV_/RANDOM_ prefixes (`pg_adapter.py:154,394`), missing `_rev`. Shared peptides: null group. Ungrouped single-accession peptide: anchor without membership (`feature_adapter.py:462`) | [[OPENMS-22](findings-cross.md#openms-22)] [[OPENMS-27](findings-feature.md#openms-27)] |
| OpenMS `-out_qpx` | Protein fields come from OpenMS's writer; `qpx/converters/openms/` passes them through. Verified at HEAD: feature `pg_accessions` holds every peptide evidence, `unique` means one accession, and `anchor_protein` is the first accession; pg `is_decoy` comes from the anchor only, and `contaminant` is always null. See the [OpenMS worklist](openms-worklist.md). | — |
| mzIdentML | `is_decoy = any(isDecoy)` (`psm_records.py:138`, `crosslinks.py:157-158`). pg uses prefixes (`pg_adapter.py:26,231`). Anchor = the first PDH flagged MS:1002403 or MS:1002401 (`pg_adapter.py:198-200`). cv_params hold the first protein's positions (`psm_records.py:201-216`) | [[MZID_CDAP-04](findings-psm.md#mzid_cdap-04)] [[MZID_CDAP-07](findings-pg.md#mzid_cdap-07)] [[MZID_CDAP-32](findings-psm.md#mzid_cdap-32)] |
| CDAP | `XXX_` is stripped (`cdap/base_adapter.py:140-142`), so 10/100 decoy PSMs carry target accessions. Feature group = all mapped proteins, pre/post on 94/94 rows. `pg.contaminant` is always False (`pg_adapter.py:167`) despite `sp\|Q29443\|TRFE_BOVIN` entries | [[MZID_CDAP-05](findings-cross.md#mzid_cdap-05)] [[MZID_CDAP-06](findings-pg.md#mzid_cdap-06)] [[MZID_CDAP-25](findings-pg.md#mzid_cdap-25)] [[MZID_CDAP-M3](findings-cross.md#mzid_cdap-m3)] |
| quantms-MSstats | `unique` is dataset-wide (`queries.py:19-21`), but `pg_accessions` is a per-row union: `[P1]` with `unique=false` | [[INFRA-28](findings-feature.md#infra-28)] |
| Mixed pg groups | Decoy if any member is (Spectronaut, mzIdentML), if all are (DIA-NN, OpenMS), or if the leader is (`fragpipe/pg_adapter.py:442`) | [[SPEC-M3](findings-pg.md#spec-m3)] [[SPEC-M4](findings-pg.md#spec-m4)] |

**Measured.** MaxQuant full: 14,382 of 30,623 multi-member lists are unsorted, all with the anchor first. maxquant_simple: `Potential contaminant` flags 113 of 9,301 evidence rows. It flags 22 rows that have `CON__` in `Proteins` and a target razor protein, misses 4 such rows, and flags one reverse hit (`REV__CON__Q3SZ57`, blank `Proteins`). The tool's flag is therefore not an "any entry" rule.

## Why it matters
- mokume reads `unique` as group-specific: on MaxQuant full it drops 43,823 measured rows, all rows of the 85 multi-member groups [[MOKUME-01](findings-feature.md#mokume-01)], and depends on whether proteinGroups.txt was given [[MOKUME-02](findings-feature.md#mokume-02)].
- Keys on stored order split groups [[MOKUME-04](findings-feature.md#mokume-04)]. Keys on the anchor merge them [[VALIDATOR-25](findings-pg.md#validator-25)] [[MOKUME-05](findings-cross.md#mokume-05)].
- Decoy status differs by producer [[MZID_CDAP-04](findings-psm.md#mzid_cdap-04)]. Stripped markers merge decoys into target records [[MZID_CDAP-05](findings-cross.md#mzid_cdap-05)] [[MZID_CDAP-06](findings-pg.md#mzid_cdap-06)]. Readers matching substrings keep 255 `REV__` rows [[MOKUME-15](findings-feature.md#mokume-15)].
- FragPipe and CDAP contaminants cannot be found even by parsing, and mokume misses `CON__` [[MOKUME-17](findings-feature.md#mokume-17)] [[CRITIC-10](findings-cross.md#critic-10)].

## Options considered
| Question | Option | For | Against |
|---|---|---|---|
| Uniqueness | Redefine `unique` as group-specific | Matches mokume and `pg.md:82` | Silent meaning change (P-A) |
| | **Add `group_unique`; derive `unique`** | Additive; separates assigned from unassigned shared peptides | One column |
| Shared peptide | Always null ([#322](https://github.com/bigbio/qpx/issues/322)) | Simple | Loses razor assignments |
| | **Group only if the source names it** | Matches the dev yaml | Readers handle nulls |
| Contaminant on psm/feature | Join to pg | No column | Fails for null membership and psm-only data |
| | **Typed nullable `contaminant`** | P-C | New column |
| Unknown membership | Anchor alone | Keeps the razor protein | Breaks `pg.md:75` |
| | **Both null; `razor_protein` cv_param** | Keeps `pg.md:75` | Lookup in cv_params |
| Decoy entries on target rows | Keep them (as mzIdentML and mzTab do) | Lossless | Readers must parse |
| | **Leave them out** | No parsing | Decoy mappings of target peptides are not in psm |
| `pg_accessions` type | **Structs with null positions; `list<string>` in 2.0** | Non-breaking now | Two steps |

## Recommendation
**PRT-1 (MUST)** `psm.protein_accessions` holds the distinct accessions of every entry the source maps the peptide to, including contaminant entries. Sources are mzIdentML PeptideEvidence, MaxQuant `Proteins`, FragPipe `Protein` + `Mapped Proteins`, PD `Protein Accessions`, OpenMS PeptideEvidence and the CDAP cell, with the engine's own matching (I/L). It holds no razor protein, group or inference ([#78](https://github.com/bigbio/qpx/issues/78)). Null means the source lists no entries (MaxQuant reverse hits); the value is never an empty list. Order carries no meaning: producers SHOULD sort by code point and permute `protein_starts` to match. Positions go to `protein_starts`, `pg_positions` or pepmap, never to cv_params.

**PRT-2 (MUST)** Rows with `is_decoy = false` hold no decoy entry in `protein_accessions`, `pg_accessions` or `pg_positions`, nor does a target peptidoform's pepmap row. This narrows `psm.md:57` for target+decoy peptides, so that target counts never include decoy proteins.

**PRT-3 (MUST)** Accessions follow `pg.md:33,44-50`: the bare accession, with the database's decoy or contaminant marker kept verbatim as a prefix or suffix (`REV__P04114`, `rev_P12345`, `contam_P02768`, `CON__P04264`, `XXX_NP_005751.2`, `P12345_rev`). An entry is the same string in every view. Producers never strip a marker. They add one only to separate a decoy that the source merged with its target (DIA-NN `DECOY_`).

**PRT-4 (MUST)** In psm and feature, a peptide is a decoy if and only if every entry it maps to is a decoy.
- Per-entry flags (mzIdentML `isDecoy`, OpenMS protein `target_decoy`) are combined with all().
- Row-level tool flags (MaxQuant `Reverse`, OpenMS peptide `target_decoy`, DIA-NN `Decoy`) are used as they are, including on rows with no entries.
- Otherwise the database marker is applied to every mapped entry, never to the razor protein alone.

psm and feature rows of one peptidoform have equal, non-null `is_decoy`; without target-decoy search it is `false`, with provenance `decoy_policy = not_searched` ([T03](t03-filtering-state.md) FLT-10). A marker used SHOULD be recorded as provenance `decoy_marker` and `decoy_marker_position` (`prefix` | `suffix`) of `format_conversion`.

**PRT-5 (MUST)** pg `is_decoy` is true if and only if every member is a decoy ([T12](t12-protein-groups.md) PGV-10). A feature and its linked pg row have equal `is_decoy`.

**PRT-6 (MUST)** psm and feature get `contaminant` (bool, nullable), with pg's meaning: the identification or group is a potential contaminant. The value comes from:
1. the tool's flag: MaxQuant `Potential contaminant` (msms rows through `Evidence ID`) or PD `Contaminant`;
2. otherwise the marker or contaminant list recorded as provenance `contaminant_marker`: true when any mapped entry (psm, feature) or member (pg) matches it, false when none does;
3. otherwise null.

A constant `false` is never written. `contaminant` never changes `is_decoy`, so `REV__CON__` rows are both.

**PRT-7 (MUST)** `feature.pg_accessions` is the full membership of the one group, in the producer's grouping ([T12](t12-protein-groups.md) PGV-6), that the feature is attributed to. If the dataset has a pg view, some pg row has the same membership set. Unknown membership is null. The column never holds:
- leaders of several groups (MaxQuant `Leading proteins`);
- a lead-only list (FragPipe without `--pg-file`);
- a union of several masters or `ProteinName` values (PD, quantms-MSstats);
- a mapping set, unless provenance records `pg_grouping = mapping_set` (CDAP).

Entries are distinct, with null positions, and SHOULD be sorted by code point. FeatureWriter accepts `list<string>`, which makes the `feature.md` examples valid.

**PRT-8 (MUST)** A peptide is *shared* when the grouping puts its target entries in more than one group; entries in no group are ignored. It carries a group only when the source names the assignment (MS:1003015 razor peptide): MaxQuant `Protein group IDs` + `Leading razor protein`, or FragPipe `Protein`. Choosing by position or sort order (PD `masters[0]`) is forbidden. An unassigned shared peptide keeps its intensities but has null `pg_accessions`, `anchor_protein` and `pg_group_id`. Its evidence proteins are in psm and pepmap.

**PRT-9 (MUST)** New feature column `group_unique` (bool, nullable):
- true: the peptide is evidence for exactly one group, the one in `pg_accessions`;
- false: it is evidence for several groups, whether assigned or not;
- null: the producer has no grouping for the peptide (MaxQuant without proteinGroups.txt), or the peptide is in no group.

Under `mapping_set`, grouped rows are true. A producer that writes any non-null `pg_accessions` SHOULD fill `group_unique` on every row (MUST under the [T13](t13-metadata-conformance.md) `dda` profile).

**PRT-10 (MUST)** `unique` is null when `pg_accessions` is null, false when the group has more than one member, and otherwise equals `group_unique`. This is `versioning.md:68-69` unchanged. `unique` is deprecated in 1.4, and 2.0 removes it. This reverses the `is_unique` rename in `feature.md:117`, because `pepmap.is_unique` already uses that name with another meaning.

**PRT-11 (MUST)** `anchor_protein` (pg, feature) is the group representative, MS:1002403, replacing MS:1001591 in `ontology.md` and "leading protein" in `pg.md:39`: mzIdentML the MS:1002403 PDH, else the first MS:1002401 PDH; MaxQuant the first `Majority protein IDs`; FragPipe `Protein`; PD the master; DIA-NN the first `Protein.Group`; otherwise the first member in source order.

The anchor is a member of `pg_accessions`. On feature it is non-null exactly when `pg_accessions` is, and it equals the anchor of pg rows with the same membership. With unknown membership the razor protein MAY go to cv_param `razor_protein`.

**PRT-12 (MUST)** In one file, psm rows with the same `sequence` have the same `protein_accessions` and `contaminant`. Feature rows with the same `sequence` have the same `pg_accessions`, `anchor_protein`, `group_unique`, `unique` and `contaminant` across runs, charges, labels (TMT, SILAC) and modifications ([#51](https://github.com/bigbio/qpx/issues/51)): MUST when provenance records `pg_grouping = dataset_inference` ([T12](t12-protein-groups.md) PGV-6), which native OpenMS output always is; SHOULD otherwise (merged files; consensusXML of other OpenMS workflows with per-run inference).

**PRT-13 (MUST)** `pg_group_id` = `derive_id([members deduplicated and sorted per identity.md:51])`, with a test vector.

**PRT-14 (SHOULD)** pepmap has one row per target peptidoform of the rank-1 psm rows and the feature rows ([T01](t01-psm-rows.md) PSM-8). Decoys are left out, because pepmap has no `is_decoy`. Entries follow PRT-1. `is_unique` means exactly one entry (MS:1001363).

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| `contaminant` | psm, feature | bool, null = unknown | optional | No |
| `group_unique` | feature | bool, null | SHOULD with `pg_accessions`; MUST under `dda` | No |
| `unique` derived, deprecated | feature | bool, null | — | Removed in 2.0 |
| `pg_accessions` positions null; writer takes strings | feature | `list<pg_protein>`; `list<string>` in 2.0 | — | Retype in 2.0 |
| No decoy entries on target rows | psm | — | — | Meaning narrowed (1.3 OpenMS target+decoy rows hold them) |
| Anchor null if and only if membership null; corrected `unique` | feature | — | — | Data change: MaxQuant without proteinGroups.txt, FragPipe without `--pg-file`, PD, OpenMS, CDAP, quantms |
| `contaminant` null = unknown; all-members decoy rule | pg | bool | — | No |
| `anchor_protein` → MS:1002403 | ontology, pg.md | — | — | No |
| `decoy_marker`, `decoy_marker_position`, `contaminant_marker` | provenance | string | when markers are used | No |
| cv_param `razor_protein` | feature | accession | MAY | No |
| `pg_group_id` dedup + test vector | feature, pg | int64 | — | Only for lists with duplicates |
| Footer | — | none new; `decoy_policy`, `contaminant_policy` ([T03](t03-filtering-state.md) FLT-10/11) and `pg_grouping` ([T12](t12-protein-groups.md) PGV-6) are provenance parameters | — | — |
| `dda` validator errors: anchor outside membership (today a warning, `qpx/core/data/validation.py:524`); membership without anchor; duplicate members; non-null positions; `unique` ≠ PRT-10; `pg_group_id` ≠ hash; `is_decoy` differing between feature and pg or psm; PRT-12. `pg_link_coverage` skips null-membership and decoy rows (`qpx/core/label_validation.py:404-411`) | all | — | — | Profile only |

## Producer changes
| Producer | Change needed |
|---|---|
| MaxQuant | Unresolved `Protein group IDs`: null membership, anchor, `unique`, `group_unique`; `razor_protein`; warn that proteinGroups.txt is needed. `group_unique` = one id. `contaminant` from evidence `Potential contaminant` (msms via `Evidence ID`), else `CON__`. Sort members |
| FragPipe | Keep `contam_` in `protein_accessions()`; `canonical_accession` accepts lower-case markers. `is_decoy` over `Protein` + `Mapped Proteins` of the row's own peptidoform. Without `--pg-file`: null membership. `group_unique` from combined_protein |
| Proteome Discoverer | `group_unique` = one master. Several masters, or no Proteins.txt: null group. Members from the master's own group. Derive `unique` and fix the test. `contaminant`: the PSM `Contaminant` column (to confirm), else Proteins.txt entries with `Contaminant` = True; pg null when absent. `is_decoy` false with provenance `decoy_policy = removed` ([T03](t03-filtering-state.md)) |
| OpenMS consensusXML | One decoy rule: per-entry `target_decoy`, else PeptideIndexer `decoy_string` and position from the header parameters. Drop the pg prefix fallback. No anchor without a group. `group_unique` from candidates. `contaminant` by marker |
| OpenMS `-out_qpx` | Upstream change to OpenMS's writer: `pg_accessions` = the group's membership (PRT-7; today every peptide evidence), the group's anchor (PRT-11; today the first accession), `unique` per PRT-10 (today one accession), pg `is_decoy` over all members (PRT-5; today the anchor only), plus `group_unique` and `contaminant`. Until then `group_unique` and `contaminant` are absent or null, i.e. unknown. See the [OpenMS worklist](openms-worklist.md). |
| mzIdentML | `is_decoy = all(isDecoy)` (psm, crosslinks). pg decoy from evidence. Anchor per PRT-11. No position cv_params. pepmap from rank-1 rows. `contaminant` null unless a marker is given |
| CDAP | Keep `XXX_`. Pre/post go to pepmap. `contaminant` from a supplied list, else null. Provenance `pg_grouping = mapping_set` |
| quantms-MSstats | Membership = one `ProteinName` split on `;`. Several values per sequence: null group, `group_unique = false`. `unique` per PRT-10. Anchor = first in source order |

## Reader rules
**PRT-R1 (MUST)** Treat `protein_accessions` and `pg_accessions` as sets. Key groups on `pg_group_id` or the sorted distinct membership, never on stored order, the first element or `anchor_protein`. An anchor-keyed output (iBAQ, sites) declares `var_key = anchor_protein` and is not a group quantity. Violations:
- mokume stored order and Ratio's first accession (`mokume:python/mokume/preprocessing/aggregation.py:167`, `mokume:rust/crates/mokume-pipeline/src/lib.rs:8830-8846`);
- anchor keys in mokume (`mokume:python/mokume/io/feature.py:1018-1030`, `mokume:python/mokume/pipeline/stages.py:1771-1776`, `mokume:python/mokume/tissuemap/loader.py:51`) and qpx (`qpx/core/data/feature.py:33-41,91-115`, `qpx/cli/transform/quantify.py:145,262`); `views.md:212` is stale.

**PRT-R2 (MUST)** Exclude decoys by `is_decoy = true`, never by accession text. A null `is_decoy` makes the file invalid and is never read as target. Violations:
- substring matching (`mokume:python/mokume/core/contaminants.py:8,15-42`);
- null read as target (`mokume:rust/crates/mokume-pipeline/src/lib.rs:8461`, `mokume:python/mokume/tissuemap/loader.py:56`, `qpx/mudata.py:278-281`, `qpx/cli/transform/quantify.py:305`);
- `Feature` has no `targets_only()`, and `for_quantification` keeps decoys.

**PRT-R3 (MUST)** Remove contaminants by the typed `contaminant`; null is unknown, never false. For 1.3 psm and feature files use the linked pg row, else unknown. Violations: mokume substrings (`mokume:python/mokume/io/feature.py:185-188,292-317`; `lib.rs:8902-8909`); qpx views and MuData ignore `contaminant`.

**PRT-R4 (MUST)** Group-specific = `group_unique = true`; single-protein group = `unique = true`; proteotypic = `pepmap.is_unique`. Never substitute one for another. The [#347](https://github.com/bigbio/qpx/issues/347) fallback (null `unique` with one member counts as unique) applies only below `qpx_version` 1.1.3. Violation: mokume filters on `unique` and applies the fallback to every version (`mokume:python/mokume/io/feature.py:224-230`; `lib.rs:8462-8464`).

**PRT-R5 (MUST)** Null `pg_accessions` means unattributed: with `group_unique = false` an unassigned shared peptide, with null no known grouping. Never attribute it by anchor or first evidence protein; an attribution through pepmap is a declared inference. Violation: mokume drops these rows even with `--keep-shared-peptides` (`mokume:python/mokume/io/feature.py:295`; `lib.rs:8467`).

## Open questions for maintainers
1. Retype `feature.pg_accessions` to `list<string>` in 2.0, when `unique` is removed?
2. Does any PD export name a shared peptide's razor or quan group? If so, PRT-8 uses it.
3. Should the tool's contaminant flag outrank the "any entry" rule? MaxQuant differs on 5 fixture rows, and [T12](t12-protein-groups.md) PGV-10 assumes the flag outranks it.
4. Should decoy mappings of target peptides get a separate column (PRT-2)?
5. Rename pepmap `pg_accessions` (evidence proteins) to `protein_accessions` in 2.0?
