# T04. Confidence fields: q-value levels, PEP, score names and the confident() helpers

[← T03](t03-filtering-state.md) · [Overview](index.md) · [Strict profile](strict-profile.md) · [Annex overview](annex-overview.md) · [T05 →](t05-feature-rows.md)

!!! note "Annex: producer-neutral analysis"
    Part of the annex of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text, and they cover every producer. For OpenMS-written files the [strict profile](strict-profile.md) is the contract and wins where the two differ. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** The column name fixes each typed confidence column's level; its scope is fixed or stated by the ontology property `scope` on the column's `ontology.parquet` row. Writers fill typed columns only from level-explicit names; in 1.4 the score entry stays beside the typed column with the same value, and is removed in the next minor version. The identifying PSM's q-value gets `feature.psm_qvalue` and is never promoted to precursor or peptide level. In DDA, feature PEP is the identifying PSM's. `confident()` takes an explicit level, settles nulls only through [T03](t03-filtering-state.md) (provenance `<level>_fdr`, FLT-R4), and raises instead of returning an empty set.

## Questions a reader has
- What level and scope does each column hold: `psm.qvalue`, PSM PEP, feature `peptide_qvalue` / `precursor_qvalue` / `global_qvalue` / PEP / `match_qvalue`, and pg `pg_qvalue` / `global_qvalue`? Is `global_qvalue` the same in every view?
- Can a best-PSM q-value become a precursor or peptide q-value? What did 1.3 files put there?
- Whose PEP does a feature carry? What does a transferred row hold?
- Which q-value wins when there are several? Can 1.0 placeholders or `1 − probability` appear?
- Typed column or `additional_scores` ([#377](https://github.com/bigbio/qpx/issues/377)): which is right?
- Why does `confident()` return nothing for MaxQuant, FragPipe and Mascot?

## What happens today
**Spec.**
- `qvalue`, `q_value` and `psm_level_qvalue` fill both `psm.qvalue` and `feature.precursor_qvalue` ("Precursor, per run") (docs/spec/scores.md:40-45).
- `psm.qvalue` uses a denylist of name tokens. The code's list (qpx/writers/derived.py:15-21) is longer than the spec's (scores.md:42).
- `global_qvalue` is a PSM value (psm.md:291), a precursor value (scores.md:45) and a protein-group value (pg.md:92). The ontology maps it to the list statistics MS:1002350 and MS:1001214 (ontology.md:238,277).
- psm.md:286 says features carry no PEP and points to `best_id_score`, which exists nowhere (it also appears at scores.md:122). feature.yaml:37-41 and feature.md:56 define a feature PEP "for the peptide match".
- `match_qvalue` is float32 with no range (feature.yaml:220-224). `channel_detail.qvalue` has no range (types.yaml:60).

| Producer | Behaviour | Findings |
|---|---|---|
| OpenMS consensus | **PSM-level mode.** The q-value (primary score or Percolator meta value) is written as `q-value` and becomes `psm.qvalue`. The PSM q-value of each run's first PID becomes `feature.peptide_qvalue` (qpx/converters/openms_consensus/feature_adapter.py:325-345,736-745). **Peptide-level mode.** `psm.qvalue` holds peptide q-values or 1.0 placeholders, and PSM PEP is null (psm_adapter.py:215-224,279). Feature PEP and `peptide_qvalue` are the peptide's best values over the whole identification source, which may be another run (feature_adapter.py:163-206). Score type `fdr` becomes `q-value` (psm_adapter.py:363). A PEP primary score stays in `additional_scores`, and PEP is null. | [[OPENMS-04](findings-psm.md#openms-04)] [[OPENMS-05](findings-feature.md#openms-05)] [[OPENMS-06](findings-feature.md#openms-06)] [[OPENMS-M1](findings-psm.md#openms-m1)] [[INFRA-03](findings-feature.md#infra-03)] [[INFRA-06](findings-psm.md#infra-06)] |
| MaxQuant | No q-values. Feature PEP comes from the evidence and is null on 104,491 MBR rows. The experiment `Q-value` goes to both pg `global_qvalue` and `pg_qvalue` (maxquant/pg_adapter.py:458,483-484). | [[MAXQUANT-05](findings-psm.md#maxquant-05)] [[MAXQUANT-18](findings-pg.md#maxquant-18)] [[MAXQUANT-19](findings-pg.md#maxquant-19)] [[MAXQUANT-26](findings-feature.md#maxquant-26)] |
| FragPipe | PEP = 1 − `Probability` (fragpipe/psm_adapter.py:245-249), labelled `percolator:PEP`. No q-values. pg `pg_qvalue` is in practice 1 − probability (fragpipe/pg_adapter.py:406-421). | [[FRAGPIPE-02](findings-psm.md#fragpipe-02)] [[FRAGPIPE-21](findings-pg.md#fragpipe-21)] |
| Proteome Discoverer | `precursor_qvalue` is the PSM q-value of the lowest-q PSM (proteome_discoverer/feature_adapter.py:17,125). `Exp. q-value: Combined` goes to `pg_qvalue` (column_mappings.yaml:75; pg_adapter.py:147-148). | [[PD-10](findings-feature.md#pd-10)] [[PD-15](findings-pg.md#pd-15)] [[PD-M2](findings-feature.md#pd-m2)] |
| mzIdentML | `psm.qvalue` is the first matching cvParam, whether it is the engine's or the rescorer's. PEP is also read from MS:1002352, a PSM-level **p-value** (mzidentml/scores.py:18,136-148; core/cv_terms.py:48). | [[MZID_CDAP-03](findings-psm.md#mzid_cdap-03)] [[MZID_CDAP-12](findings-psm.md#mzid_cdap-12)] [[INFRA-05](findings-psm.md#infra-05)] [[INFRA-09](findings-psm.md#infra-09)] |
| CDAP | The best PSM's `Qvalue` (`psm_level_qvalue`, cdap/feature_adapter.py:289-295) becomes `precursor_qvalue` (94/94). `PepQvalue` fills `peptide_qvalue`, contrary to converter-coverage.md:181. pg scores carry PSM- and peptide-level q-values. | [[INFRA-01](findings-feature.md#infra-01)] [[MZID_CDAP-13](findings-feature.md#mzid_cdap-13)] [[MZID_CDAP-26](findings-pg.md#mzid_cdap-26)] |
| quantms/MSstats | No confidence fields in the input; all null. | – |
| DIA-NN, Spectronaut (reference) | DIA-NN: true run- and experiment-level precursor q-values (diann/feature_sql.py:14-16) and a run-level precursor PEP. Spectronaut: `EG.Qvalue` is `qvalue`; the experiment-wide `PG.Qvalue` is `pg_qvalue` on features (spectronaut/feature_adapter.py:39-40) but `global_qvalue` on pg (converter-coverage.md:299). | [[INFRA-02](findings-feature.md#infra-02)] [[INFRA-07](findings-cross.md#infra-07)] |
| Shared code | The writer fills nulls only and keeps an unchecked duplicate (derived.py:58-71). | [[INFRA-08](findings-cross.md#infra-08)] [[INFRA-M1](findings-cross.md#infra-m1)] [[VALIDATOR-07](findings-feature.md#validator-07)] [[VALIDATOR-08](findings-cross.md#validator-08)] [[CRITIC-M5](findings-cross.md#critic-m5)] |

## Why it matters
- `feature.precursor_qvalue` is a run-level precursor FDR for DIA, an anti-conservative best-of-N PSM q-value for PD and CDAP, and empty for MaxQuant, FragPipe and OpenMS.
- `feature.confident()` returns 0 of 243,821 MaxQuant features. `psm.confident()` returns 0 of 13,386 Mascot PSMs, although 11,517 have `pass_threshold = true`.
- mokume's `--filter-peptide-fdr` filters OpenMS data on PSM-level q-values [[MOKUME-M1](findings-feature.md#mokume-m1)].
- Run-level `pg_qvalue` filters read experiment values or 1 − probability [[SPEC-24](findings-pg.md#spec-24)]; score filters drop all transfers [[MOKUME-M3](findings-feature.md#mokume-m3)].

## Options considered
| Question | Option | For | Against |
|---|---|---|---|
| Best-PSM q-value on features | Keep in `precursor_qvalue` | No change | Mixes levels; anti-conservative |
| | Only as a score | No schema change | Reads about 100× the data (scores.md:38) |
| | **Typed `feature.psm_qvalue`** | Pairs with PEP, `scan` and `id_run_file_name` ([#378](https://github.com/bigbio/qpx/issues/378)) | One more optional column |
| Scope | One fixed scope per column | Simple | Scope depends on configuration |
| | **Ontology property where scope varies** | Additive; `unknown` allowed; read with SQL | Readers must check `ontology.parquet` |
| [#377](https://github.com/bigbio/qpx/issues/377) duplication | **Copy, plus an equality check (1.4)** | Compatible with 1.3; consumers look q-values up by name in `additional_scores` today | Two stores until the entries are removed |
| | Move to the typed column | One store; smaller files | Breaks lookups by name; deferred to the next minor version, announced in the changelog |
| Fill names | Denylist | Catches new names | Admits gene, match q-values [[INFRA-05](findings-psm.md#infra-05)] |
| | **Allowlist** | Testable | Converters name the level |
| Feature PEP | Drop it | Matches psm.md:286 | Loses data most converters fill |
| | **The identifying PSM's PEP** | What MaxQuant, FragPipe and PD write | Nulls OpenMS peptide-level PEP |
| `confident()` | Silent fallback to `pass_threshold` or PEP | Non-empty results | Mixes bases silently |
| | **Explicit level; nulls settled only by [T03](t03-filtering-state.md)'s `<level>_fdr` parameters and the declared `pass_threshold` basis; raise when nothing is assessable** | Portable; matches [T03](t03-filtering-state.md) FLT-R3, [T14](t14-reader-contract.md) RDR-3 | Callers must choose |

## Recommendation
**CNF-1 (MUST)** Columns have the levels and scopes below. "global" in a QPX name means experiment scope at the row's own level. It never means a PSI-MS "global FDR" (MS:1002350, MS:1001364, MS:1001214, MS:1002369), which are statistics of a whole list ([T03](t03-filtering-state.md)).

| Column | Level | Scope | PSI-MS |
|---|---|---|---|
| `psm.qvalue` | PSM | ontology property `scope` (CNF-3) | MS:1002354 |
| `psm.posterior_error_probability` | PSM | – | CNF-7 |
| `feature.psm_qvalue` (new) | identifying PSM ([T06](t06-identifications-per-feature.md)) | ontology property `scope` (CNF-3) | MS:1002354 |
| `feature.posterior_error_probability` | identifying PSM (DIA: precursor, CNF-8) | – | CNF-7 |
| `feature.precursor_qvalue` | precursor (peptidoform, charge) | run | none |
| `feature.global_qvalue` | precursor | experiment | none |
| `feature.peptide_qvalue` | peptide | ontology property `scope` (CNF-3) | MS:1001868 (per stripped sequence; open question 2) |
| `feature.channel_details.qvalue` | the label's precursor | run | none |
| `pg.pg_qvalue` / `pg.global_qvalue` | protein group | the unit (`grouped_runs`) / experiment | MS:1002373 |
| `feature.pg_global_qvalue` | protein group | experiment; equals `pg.global_qvalue` of the row's group | MS:1002373 |
| `feature.match_qvalue` | transfer ([T08](t08-match-between-runs.md)) | – | none |

`run` means estimated within the row's run, or its unit when `identity_composite` starts with `quantification_unit_id` ([T05](t05-feature-rows.md) FRW-3). `experiment` means one estimate over all runs of the analysis that produced the file.

**CNF-2 (MUST)** A value enters a typed column only when the source states its level. For a fixed-scope column, the source must also state that scope. For `psm.qvalue`, `psm_qvalue` and `peptide_qvalue`, the scope goes in the CNF-3 ontology property. Forbidden:
- any PSM's PSM-level q-value, the best one included, in `precursor_qvalue` or `peptide_qvalue`;
- an experiment-scope value in a run-scope column, or the reverse;
- `1 − probability` in a q-value column;
- a non-monotone FDR (OpenMS `fdr`) as a q-value.

Tool placeholders (Percolator's peptide-level 1.0) are written as null.

A typed q-value column is filled exactly when the source estimated FDR at that level and scope. A level that was not estimated is null on every row of the file. This is legitimate for every producer, and it is not a defect and not `unknown`: OpenMS always estimates PSM-level FDR, and it sometimes stops there or at peptide level (for example in immunopeptidomics) without a protein-level FDR.
- `feature.peptide_qvalue` is filled when peptide-level FDR was estimated.
- `feature.psm_qvalue` and `psm.qvalue` are filled when PSM-level FDR was estimated. OpenMS always estimates it, so in OpenMS-written files `psm.qvalue` is filled on every row and `feature.psm_qvalue` on every direct and secondary row ([strict profile](strict-profile.md)).
- `pg.global_qvalue` is nullable for every producer, OpenMS included ([T12](t12-protein-groups.md) PGV-8).

**CNF-3 (MUST)** The scope is stated by the ontology property `scope` ([T13](t13-metadata-conformance.md) META-1) on the `ontology.parquet` rows of `psm.qvalue`, `feature.psm_qvalue` (when filled) and `feature.peptide_qvalue`: `run`, `experiment` or `unknown`. It describes the column, so it is neither a footer key nor a provenance parameter; the scope of a filter cut stays [T03](t03-filtering-state.md) FLT-4's provenance `<level>_fdr_scope`.
- A writer that fills the column writes the property. An absent property or row means unknown. `psm.qvalue` and `feature.psm_qvalue` describe the same PSMs and carry the same value.
- A merge writes `run` only if every input states `run`; otherwise it writes `unknown` and warns ([T13](t13-metadata-conformance.md) META-4).

**CNF-4 (MUST)** With several PSM-level q-values, `psm.qvalue` takes, in order:
1. the value the source marks as final (mzIdentML MS:1002354, or the tool's documented final column);
2. else a post-processor's (Percolator MS:1001491, mokapot, MS2Rescore, ConsensusID);
3. else the engine's (MS-GF:QValue MS:1002054).

The others stay as `<tool>_qvalue`; PSM PEP follows the same order. File and column order never decide. mzIdentML reads only the lists [T01](t01-psm-rows.md) PSM-9 selects. ontology.parquet records the winner's `source_column_name` and `source_tool`.

**CNF-5 (MUST)** Every q-value, PEP, probability or FDR entry in `additional_scores` has a level, stated by the property `level` in the `properties` list of its ontology.parquet row ([T13](t13-metadata-conformance.md) META-1): `psm`, `precursor`, `peptide`, `protein`, `protein_group`, `gene_group`, `transfer` or `other`. A reserved name has a fixed level, which its property repeats; a `<tool>_` name gets its level from the property alone.
- Reserved names:
  - scope from CNF-3: `psm_qvalue`, `peptide_qvalue`;
  - run scope: `precursor_qvalue`, `protein_qvalue`, `pg_qvalue`;
  - experiment scope: `protein_global_qvalue`, `pg_global_qvalue`, and `global_qvalue` on feature and pg rows only;
  - no scope guarantee: `peptide_pep`, `precursor_pep`, `pg_pep`, `pg_probability`.
- `higher_better` is false except for probabilities.
- The bare names `qvalue`, `q_value`, `q-value` and `fdr` are deprecated. During 1.4 they only raise a validator warning, because CNF-6 keeps the entries for existing readers. They become an error when the entries are removed.
- Feature rows carry no `rank`. pg rows carry no reserved PSM-, precursor- or peptide-level names; a lower-level summary on pg uses a non-reserved name whose ontology property `level` states its level ([T12](t12-protein-groups.md) PGV-8).

**CNF-6 (MUST)** The writer fills typed columns only from these names:

| Column | Filled from |
|---|---|
| `psm.qvalue` | `psm_qvalue`, `psm_level_qvalue`; deprecated `qvalue`, `q_value` |
| `feature.psm_qvalue` | `psm_qvalue`, `psm_level_qvalue` |
| `feature.peptide_qvalue` | `peptide_qvalue`, `peptide_level_qvalue` |
| `feature.precursor_qvalue` | `precursor_qvalue` |
| `feature.global_qvalue` | `global_qvalue` |
| `psm.rank` ([T01](t01-psm-rows.md) PSM-3) | `rank` |

- In 1.4 the writer fills the typed column and keeps the score entry, with the same value (dual-write), because consumers look q-values and `rank` up by name in `additional_scores` today. The validator checks that the two are equal.
- Removing the entries is deferred to the next minor version and announced in the changelog. This settles the q-value part of [#377](https://github.com/bigbio/qpx/issues/377) (open).
- Two entries that would fill one column on the same row are an error.
- Where the two differ, the typed column wins (1.3 and 1.4 files), and the validator reports the inequality.

**CNF-7 (MUST)** `psm.posterior_error_probability` is the probability that the PSM is incorrect.
- Sources: MS:1001493, MS:1002056, MS:1001901, MS:1002995, MS:1002265, MS:1003113 or MS:1002351 (PSM-level local FDR); the OpenMS score type `Posterior Error Probability`; or 1 − a PSM-level probability (MS:1002357, PeptideProphet, iProphet, FragPipe `Probability`).
- ontology.parquet names the source term, and uses MS:1002351 when the term is unknown.
- Never a p-value (MS:1002352), an e-value or a peptide-level PEP; a peptide-level PEP goes to `peptide_pep`.
- PEP stays float64, which decides [#377](https://github.com/bigbio/qpx/issues/377) B2 (open). `BYTE_STREAM_SPLIT` is allowed; −log10 is not.

**CNF-8 (MUST)** On a feature row with an identifying PSM ([T06](t06-identifications-per-feature.md) FID-1), PEP, `psm_qvalue` and the `psm`-level `additional_scores` describe that PSM. They equal its psm row when that row is present.
- Rows without an identifying PSM (transferred rows, rows of a multi-run unit; [T06](t06-identifications-per-feature.md) FID-1) have null PEP and `psm_qvalue` and no `psm`-level entries. Unidentified rows ([T07](t07-unidentified-features.md) UNI-3) have a null PEP and null q-values at every level. In runs whose `run.acquisition_method` is `dia`, PEP is the tool's run-level precursor PEP.
- Transferred rows also have null `precursor_qvalue`, and null `peptide_qvalue` unless its scope is `experiment`; they MAY carry `global_qvalue` and `pg_global_qvalue`. Transfer confidence goes only in `match_*` ([T08](t08-match-between-runs.md)).
- Delete the last sentence of psm.md:286, scores.md:122 and the PSM `global_qvalue` recommendation at psm.md:291.

**CNF-9 (MUST)** Every q-value column and struct field, `match_qvalue` and `channel_detail.qvalue` included, is float64 in [0, 1].

**CNF-10 (SHOULD)** Rows with rank > 1 and `secondary` rows keep the source's `qvalue` and PEP. FDR control holds only over rank-1 rows ([T01](t01-psm-rows.md) PSM-R1), and over secondary rows only as far as [T02](t02-chimeric-spectra.md) CHI-8 allows.

**CNF-11 (MUST)** ontology.parquet and ontology.md follow CNF-1 and CNF-7 (qpx/core/builtin_scores.py:185-191,241-254; ontology.md:164,220,238,276-278). This covers the properties: the row of every typed confidence column carries the `level` of CNF-1 (`psm`, `precursor`, `peptide`, `protein_group` or `transfer`), the rows CNF-3 names carry `scope`, and score rows carry `level` per CNF-5. The properties sit in the `properties` list ([T13](t13-metadata-conformance.md) META-1); the ontology gets no separate `level` column.

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| `psm_qvalue` column | feature | float64 [0,1], nullable, optional | no | No |
| Property `scope` on the rows of `psm.qvalue`, `feature.psm_qvalue` and `feature.peptide_qvalue` (no footer key) | ontology | `run` \| `experiment` \| `unknown` | when the column is filled | No |
| Property `level` in the new `properties` list ([T13](t13-metadata-conformance.md) META-1) | ontology | CNF-5 values | on every confidence column and score name (CNF-5, CNF-11) | No |
| `match_qvalue` → float64 [0,1]; `channel_detail.qvalue` range | feature | float64 | – | **Yes** (physical type) |
| Allowlist fill; in 1.4 the score entry is kept beside the typed column (same value, equality checked); the entries are removed in the next minor version | psm, feature | CNF-6 | – | **Yes** in 1.4 only in that `<engine>_qvalue` no longer fills `psm.qvalue`; name lookups keep working. **Yes** in the next minor version: name lookups fail once the entries are removed (announced in the changelog) |
| Level corrections: PD and CDAP `precursor_qvalue` and OpenMS PSM-level `peptide_qvalue` move to `psm_qvalue`; the peptide-level feature PEP of the OpenMS consensusXML converter becomes null | feature | – | – | **Yes** (meaning) |
| Deprecate bare q-value names | all | warning | – | No |

## Producer changes
| Producer | Change needed |
|---|---|
| OpenMS consensus (and upstream `-out_qpx`) | **PSM-level:** emit `psm_qvalue`; `feature.psm_qvalue` from the identifying PSM; `peptide_qvalue` null. **Peptide-level** (consensusXML converter only, when the file holds no PSM-level values): `psm.qvalue`, PSM and feature PEP null; write `peptide_qvalue`, `peptide_pep`. The native writer always exports the PSM-level q-value and PEP, because OpenMS always estimates PSM-level FDR, and adds `peptide_qvalue` when peptide-level FDR was estimated ([strict profile](strict-profile.md)). **Both:** `fdr` → `openms_fdr`; a PEP primary score fills PEP. **Scope (CNF-3):** ProteomicsLFQ computes the PSM FDR once over the whole consensus map, so the PSM q-value rows get `scope = experiment`; for consensusXML of other workflows, from the runs each ProteinIdentification covers (one → `run`, all → `experiment`, else `unknown`). OpenMS has no precursor-level FDR, so `precursor_qvalue` and feature `global_qvalue` are null. **Native writer (verified at HEAD):** no `qvalue` column, and the PSM q-value is not written at all (the main score is exported only when it is a PEP); score names are raw OpenMS keys. The unmerged [#378](https://github.com/bigbio/qpx/issues/378) branch writes `match_qvalue` as float32 (CNF-9: float64). See the [OpenMS worklist](openms-worklist.md). |
| MaxQuant | pg `pg_qvalue` null. Map proteinGroups `PEP` to `pg_pep` when present. |
| FragPipe | pg `pg_qvalue` null. The probability goes to `pg_probability`. Record the PSM validator in provenance and fix the PEP ontology row. |
| Proteome Discoverer | `precursor_qvalue` becomes `psm_qvalue` (`unknown`). `Exp. q-value: Combined` goes to `global_qvalue`. MAY read PeptideGroups q-values (`experiment`). |
| mzIdentML | MS:1002354, MS:1002054 and MS:1001491 (PSM-level) go to `psm_qvalue` per CNF-4, and MS:1002055 to `peptide_qvalue`. Drop MS:1002352 from PEP. Scope is `run` when each list covers one SpectraData. |
| CDAP | `Qvalue` becomes `psm_qvalue`. `PepQvalue` stays `peptide_qvalue` (`unknown`). Drop the pg q-values. |
| Spectronaut | `qvalue` → `precursor_qvalue`. Feature `pg_qvalue` → `pg_global_qvalue`. `pg_qvalue_runwise` → `pg_qvalue`. |
| DIA-NN, quantms/MSstats | None. |
| Writer and docs | Implement CNF-6 (allowlist, dual-write and its equality check; qpx/writers/derived.py) and the ontology properties of CNF-3, CNF-5 and CNF-11. Fix converter-coverage.md:59,181,274,299-300 and scores.md:40-56. |

## Reader rules
- **CNF-R1 (MUST)** The typed column is authoritative. A null makes the row unassessable unless [T03](t03-filtering-state.md) settles it: FLT-R3 from the provenance `<level>_fdr` of that level, or, for `psm.qvalue`, a `pass_threshold = true` whose ontology row states a PSM-level FDR ≤ t (properties `level` and `fdr`; FLT-R4). Helpers never count unassessable rows as passing. They report the count, and raise an error when every row is unassessable ([T14](t14-reader-contract.md) RDR-3). A q-value column that is null on every row means FDR was not estimated at that level (CNF-2), so filtering at that level is not available: the reader says so and filters at a level that exists. It never reads a null as 1.0 and never drops every row silently.
- **CNF-R2 (MUST)** Never `COALESCE` across levels or scopes. Never use a PSM-level value for peptide or protein FDR.
- **CNF-R3 (MUST)** When the CNF-3 `scope` property is `unknown` or absent, warn and do not report an experiment-wide FDR.
- **CNF-R4 (MUST)** `PSM.confident(max_qvalue)` keeps target rows with `qvalue ≤ t` that are rank 1 where rank is known ([T01](t01-psm-rows.md) PSM-R1, PSM-R3); null `qvalue` rows follow CNF-R1. It is a q-value filter: when `qvalue` is non-null, `pass_threshold` is not consulted. A separate `PSM.accepted()` implements [T03](t03-filtering-state.md) FLT-R4 (source acceptance), and identification counts use `accepted()`. `Feature.confident(max_qvalue, column=…)` has no default; `column` ∈ {`psm_qvalue`, `precursor_qvalue`, `peptide_qvalue`, `global_qvalue`}. Breaking for callers ([T14](t14-reader-contract.md)).
- **CNF-R5 (MUST)** On a transferred row, a column that CNF-8 leaves null is not unassessable: identification filters pass the row through, report how many they passed, and leave its fate to the reader's transfer policy ([T08](t08-match-between-runs.md) MBR-R5, [T14](t14-reader-contract.md) RDR-3). A filled experiment-scope column is filtered as on other rows.
- **CNF-R6 (MUST)** In files with `qpx_version` < 1.4, the level of feature `precursor_qvalue` and `peptide_qvalue` from PD, CDAP or OpenMS is unknown. Readers warn.

Current violations:
- qpx: `psm.confident()` ignores rank and `pass_threshold` and drops null `qvalue` rows (qpx/core/data/psm.py:23-25); `feature.confident()` uses `COALESCE(precursor_qvalue, peptide_qvalue)` (qpx/core/data/feature.py:52-54); `ProteinView.intensity` passes null `global_qvalue` rows (qpx/views/api.py:45-49).
- mokume:
  - It uses the typed column only without nulls, else a name lookup in `additional_scores` (mokume:python/mokume/io/_qpx.py:13-46; mokume:rust/crates/mokume-io/src/qpx.rs:208-219,585-611).
  - It treats `protein_qvalue` as `pg_global_qvalue` (mokume:python/mokume/io/feature.py:601-603; qpx.rs:30).
  - Its peptide FDR filter ignores scope and the declared filtering state (mokume:python/mokume/preprocessing/filters/peptide.py:387-413; mokume:rust/crates/mokume-pipeline/src/lib.rs:567-574).
  - Its score filters drop transferred rows (mokume:python/mokume/io/feature.py:263-270; lib.rs:6966-6972).
- MuData export filters only on `is_decoy` (qpx/mudata.py:279-280), which conforms.

## Open questions for maintainers
1. Should psm get a typed `peptide_qvalue` for producers that write only PSMs (mzIdentML, Sage) [[EXTERNAL-25](findings-psm.md#external-25)]?
2. Is `peptide_qvalue` per peptidoform (PD peptide groups, DIA-NN `Peptidoform.Q.Value`) or per stripped sequence (as MS:1001868 defines it)? Fix one, or record it in provenance?
3. *Resolved in the overview and by the maintainer:* the scope of a column (CNF-3, an ontology property) stays separate from the scope of a filter cut ([T03](t03-filtering-state.md) FLT-4, the provenance parameter `<level>_fdr_scope`).
4. Should `feature.global_qvalue` be renamed `precursor_global_qvalue` in 2.0?
