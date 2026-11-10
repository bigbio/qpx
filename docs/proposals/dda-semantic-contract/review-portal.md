# Review: the quantms.io portal as a consumer

[Overview](index.md) · [Strict profile](strict-profile.md) · [Consumers](consumers.md) · [OpenMS work list](openms-worklist.md)

!!! note "Review report"
    One of three reviews of the [DDA semantic contract proposal](index.md), made by reading bigbio/quantms-portal at `728f5e4` and its backend at `c562424`. It is kept as written by the reviewer. It predates the restructuring around the strict profile and the decisions taken afterwards, so parts of it are superseded: it still speaks of footer keys (the proposal now adds none; see [where a dataset says what it holds](strict-profile.md#where-the-dataset-says-what-it-holds)) and of a separate view for unidentified features (they stay in the feature view as marked rows). It also lists as open some points decided since; see [Decisions taken](index.md#decisions-taken). Where it differs from the [strict profile](strict-profile.md), the strict profile holds. Rule ids refer to the topic sections T01 to T14.

Path prefixes used below:
- `BE/` = `quantms-portal-backend:src/quantms_portal_backend/`
- `BR/` = `quantms-portal-backend:`
- `FE/` = `quantms-portal:frontend/src/`
- Proposal sections = the topic sections T01 to T14

## Summary

- **The portal reads five views (psm, feature, pg, run, sample) with plain DuckDB/pyarrow and filters only on `is_decoy`.** It reads no Parquet footer key, no QPX provenance, no id column, and none of `unique`, `rank`, `pass_threshold`, `identification_type`, `is_imputed`, `additional_intensities` or `cv_params`.
- **Four proposed changes would blank portal output for its main producers unless a typed replacement ships with them:** PGV-8, INT-4, CNF-2 and CNF-6 (details in 3b).
- **About half the proposal is irrelevant to this consumer** (chimeric spectra, PSM–feature links, unidentified-feature view, spectrum fields, ids, MaxQuant/FragPipe/PD/CDAP fixes).
- **What the portal most needs is not in the proposal:** an accepted-row predicate usable in SQL without provenance, a strict "release" profile, a stored dataset summary, and an always-filled protein quantity.
- **Nothing here was measured on live portal data.** "Wrong today" statements come from portal code plus the proposal's findings.

## 1. How the portal consumes QPX

### 1a. What is read

| Consumer | Views | Columns and fields |
|---|---|---|
| Dataset KPIs at publish and restat: `BE/qpx_stats.py:103-295, 457-572, 668-739` | psm, feature, pg, run, sample | `sequence`, `is_decoy`, `protein_accessions`, `pg_accessions`; row counts of run and sample; `sample.organism`, `run.instrument` |
| Explorer: `BE/apps/qpx_explore/store.py`, `queries.py`, `qc.py` | all five, as `SELECT *` views (`store.py:396`) | feature: `sequence`, `peptidoform`, `charge`, `run_file_name`, `peptide_qvalue`, `posterior_error_probability`, `rt`, `mass_error_ppm`, `observed_mz`, `calculated_mz`, `missed_cleavages`, `intensities[].intensity`, `modifications[].name`, `pg_accessions[].accession`, `anchor_protein`, `scan`. pg: `pg_accessions`, `anchor_protein`, `pg_names`, `global_qvalue`, `intensity`, `label`, `grouped_runs`, `peptide_counts.{unique,total}_sequences`, `sequence_coverage`, `molecular_weight`, `contaminant`. run: `samples[].{sample_accession,label}`, `fraction`, `instrument`, `acquisition_method`. sample: `organism`, `organism_part`, `disease`, `cell_line`, `sex`, `pooled_sample` |
| Frontend SQL sent to the explorer's `/sql`: `FE/qpxAnalysis.js:53-58, 139-149`; `FE/components/QxOverviewStats.vue:202-258`; `FE/views/QpxExplore.vue:441-479` | pg, run, feature, psm | pg `anchor_protein`, `grouped_runs[1]`, `label`, `intensity`, `is_decoy`, `contaminant`, `pg_names`, `gg_names`; `run.samples` |
| Peptide-search index: `BE/apps/peptide_search/index_builder.py:345-352, 536-590, 756-770`; source choice at `build_full_index.py:181-188` | psm (MSNet) or feature (QPX datasets); `msstats_in.csv` when present | `sequence`, `peptidoform`, `protein_accessions` or `pg_accessions[].accession`, `is_decoy`, `posterior_error_probability`, `charge` |
| Evaluation gate: `BE/evaluation/view.py`, `checks/*.py` | psm, feature, pg, run | PEP, typed `global_qvalue`/`qvalue`, `additional_scores` by name, `observed_mz`, `calculated_mz`, `charge`, `run.samples`, accession lists |
| Co-expression: `BE/apps/coexpression/builder.py:202-283` | pg plus SDRF | `anchor_protein`, `gg_names`, `pg_names`, `run_file_name`, `global_qvalue`, `is_decoy`, `contaminant`, `intensities[].{label,intensity}` (an older pg layout) |

- **Score names looked up in `additional_scores`:** `qvalue`, `q-value`, `q_value`, `pep`, `posterior_error_probability`, `precursor_qvalue` (`BE/evaluation/checks/ms_confidence.py:29, 72`). The planned atlas FDR also reads `'q-value'` and `'precursor_qvalue'` by name (`BR/data_processing/algorithms/atlas-fdr/ALGORITHM.md:22-23, 55-58`).
- **Footer keys read:** none. Only `num_rows` (`BE/storage.py:190`).
- **Publishable kinds:** psm, feature, pg, run, sample, dataset, ontology, provenance, precursor (`BE/qpx_names.py:11, 27`). pepmap, mz, xic and site are never published.
- **Producers ingested:** quantms (OpenMS-based DDA) and quantmsdiann (DIA-NN) outputs from the PRIDE FTP (`BE/ingest_workflow.py:36-64`; `BR/docs/backend-architecture.md:176`), and MSNet PSM-only bundles (`BE/ingest_msnet.py:30-33`). No MaxQuant, FragPipe or PD.
- **Downstream apps** (statistics, compass, protein-metadata, admin metrics) read the peptide-search index, not QPX. The cell-line atlas reads a build made in another repository.

### 1b. What it assumes

| Id | Assumption | Code |
|---|---|---|
| A1 | A non-decoy row is an accepted identification; no q-value, rank or threshold filter | `BE/qpx_stats.py:163, 179, 205, 260`; `BE/apps/qpx_explore/queries.py:46, 63, 80`; `BE/evaluation/checks/ms_yield.py:12`; `index_builder.py:553, 756` |
| A2 | `is_decoy` is complete. Nulls are dropped in some places and kept in others | dropped: `queries.py:46`, `store.py:471`, `FE/qpxAnalysis.js:57`; kept: `qpx_stats.py:163`, `BE/evaluation/view.py:102`, `builder.py:221` |
| A3 | Contaminants are recognisable from accession text, with four different patterns; a null `pg.contaminant` is read as false | `qpx_stats.py:22, 33-35`; `BE/evaluation/checks/bio.py:11, 140`; `qc.py:54`; `builder.py:59, 223`; `FE/qpxAnalysis.js:57` |
| A4 | One feature row is one quantified precursor in one run | `qpx_stats.py:242-246, 260`; `qc.py:47`; `store.py:462`; `index_builder.py:583, 766` (`n_obs`) |
| A5 | Every feature row is an own-run identification (transfers never excluded) | same sites as A4; `store.py:1010-1016` |
| A6 | One psm row is one PSM | `qpx_stats.py:179`; `store.py:699-703`; `queries.py:42-50` |
| A7 | `feature.peptide_qvalue` is the feature's q-value; DIA precursor q and DDA PSM q are found by score name; `pg.global_qvalue` is a group-level experiment q; the applied cut is inferred from `max(q)` | `store.py:1139, 1166`; `queries.py:116`; `FE/views/QpxExplore.vue:212-214, 789`; `ms_confidence.py:29, 59, 72`; `builder.py:219` |
| A8 | PEP is a PSM-level error probability | `ms_confidence.py:43`; `index_builder.py:584-585, 765`; `store.py:560` |
| A9 | A pg row means the group is identified. The group key varies: raw list, sorted list, or anchor | `store.py:471`; `qpx_stats.py:228`, `ms_yield.py:15`; `store.py:951, 997`, `FE/qpxAnalysis.js:55`, `builder.py:228` |
| A10 | `feature.pg_accessions` is a list of structs and holds the peptide's proteins | `store.py:981, 1016`; `index_builder.py:764` |
| A11 | `pg.intensity` is filled, linear, not imputed, and summable over runs and channels; 0 means missing | `FE/qpxAnalysis.js:9, 55-57`; `qc.py:181-182`; `store.py:676`; `builder.py:256, 271` |
| A12 | `(grouped_runs[1], label)` resolves to one sample through `run.samples`; co-expression instead uses the SDRF with positional channel indices and a run-only fallback | `FE/qpxAnalysis.js:54-56`; `qc.py:117-123, 181`; `builder.py:36-53, 192, 297` |
| A13 | DDA means a psm table exists; DIA means only a feature table | `BE/evaluation/view.py:37` |
| A14 | `rt` is in seconds; m/z > 0 is real; `pg_names` are gene names | `store.py:557`, `qc.py:50`; `BE/evaluation/checks/ms_integrity.py:98-102`; `store.py:673`, `queries.py:95-98`, `builder.py:135` |
| A15 | Across datasets, a protein is a taxon-scoped gene from the portal's own accession-to-gene map; intensities are never pooled across datasets | `qpx_stats.py:33, 51, 734-737`; `BE/apps/peptide_search/shard.py:137-141`; `BE/apps/coexpression/README.md` |
| A16 | Layout differences are handled by sniffing columns; `qpx_version` is never read | `store.py:144-149, 597-598, 661-670` |

## 2. Do QPX 1.3 producers satisfy these?

| Id | 1.3 reality for the portal's producers | What the user sees |
|---|---|---|
| A1, A6 | OpenMS-consensus psm holds every hit of every identification with no rank (OPENMS-01, OPENMS-03). The DIA feature table goes up to 5% precursor q (`BR/docs/status/2026-11-03-peptide-search-review.md`, finding B4) | PSM and peptide KPIs are upper bounds. DIA peptide counts are at 5%, DDA at about 1% |
| A2 | OpenMS `is_decoy` is false when the `target_decoy` annotation is absent, while pg uses prefixes (OPENMS-22). DIA-NN reports no decoys (`ALGORITHM.md:30`) | The "Decoy %" tile is understated, or shows 0 for DIA |
| A3 | psm and feature have no contaminant column; `pg.contaminant` may be null (SPEC-M4). `bio.py:11` anchors `^CONTAM_`, but published accessions can be `sp\|CONTAM_…` (`qpx_stats.py:27-32`) | Contaminant fraction differs per module; the evaluation check can read 0 |
| A4, A5 | OpenMS rows are consensus element × run, several per key; TMT rows are reporter spectra (OPENMS-09). Transfers are set by a heuristic (OPENMS-11). DIA-NN type is null | "Features" is a row count. Peptide-search observation counts and per-run protein abundance include transfers |
| A7 | OpenMS `peptide_qvalue` is the PSM q-value of the first identification (OPENMS-05). For DIA-NN it is null | The explorer's q-value column and filter are mislabelled for DDA and empty for DIA |
| A7 | OpenMS `pg.global_qvalue` is the minimum member-protein q over all inferences (OPENMS-25) | Co-expression's 1% protein cut and `ms.protein_fdr` use a non-group statistic |
| A8 | Feature PEP is from the first identification, or the best across runs (OPENMS-06) | Minor |
| A9 | OpenMS pg rows come from any hit: any rank, decoy, unfiltered (OPENMS-23). Counts are dataset-wide and repeated per row (OPENMS-26) | The protein KPI is inflated. Three code paths give three group counts |
| A10 | OpenMS shared peptides have null group and anchor (OPENMS-27) | Shared peptides lose protein attribution in peptide search and the protein pages |
| A11 | OpenMS `pg.intensity` is either the OpenMS abundance or a QPX unique-peptide sum that includes transferred rows (OPENMS-24). Nothing declares the kind | Sample map, group comparison and heatmap run on an undeclared quantity. "Never imputed" is not checked |
| A12 | Published TMT QPX can lack run and sample views (`qpx_stats.py:366, 617`) | Sample analyses are unavailable for those datasets; co-expression relies on channel position |
| A13 | A feature-only DDA dataset is classed DIA | The wrong evaluation checks run |
| A14 | The spec gives no unit; `pg_names` are entry names or FASTA descriptions | The "genes" column shows entry names |

I also checked two things directly with DuckDB 2.0:
- `qpxc transform relayout` stores `additional_scores` as a struct (`qpx:docs/spec/serialization.md:80`). On that layout the lookup at `ms_confidence.py:71-72` fails with a binder error. I confirmed the error; I could not check whether published files are relaid out.
- `count(DISTINCT pg_accessions)` treats `[A,B]` and `[B,A]` as two groups, so `store.py:471` is order-sensitive.

## 3. The proposal from the portal's side

### 3a. Rules that directly fix a portal problem

| Portal problem | Rules |
|---|---|
| The cut is unknown and inferred from `max(q)`; DDA and DIA differ (B4) | FLT-1, FLT-2, FLT-5, PGV-3 (`fdr_filter`); CNF-1, CNF-2 (one level per typed column) |
| DDA or DIA guessed from files | META-5 (`acquisition_method`) |
| PSM counts include lower ranks and failing rows | PSM-3, PSM-4 (`rank`); FLT-7, FLT-9; FLT-R4 |
| Protein count inflated and inconsistent | PGV-1, PGV-2, MBR-10 (`identified_in_unit`), PRT-13 (`pg_group_id`), PGV-R1 |
| Decoy handling | PRT-4, PRT-5, PGV-10, PRT-R2 |
| Contaminants by regex | PRT-6 (typed on psm and feature), PGV-10, FLT-R6 |
| Transfers counted as evidence or summed | MBR-1, MBR-2, MBR-3, MBR-7; MBR-R3, MBR-R4; RDR-5 |
| Feature row meaning and counts | FRW-2, FRW-3, FRW-4; INT-5, INT-6; UNI-1; FRW-R6; RDR-9 |
| Undeclared quantity, imputation, zero versus null | INT-1, INT-2, INT-7, INT-9; INT-R1, INT-R3 |
| Sample resolution; missing run view | RDR-8; PGV-5; META-6 (the `dda` profile requires a run view) |
| Anchor used as key; `pg_names` shown as genes | PRT-11, PRT-R1, RDR-7; PGV-12 |
| Units | FRW-7 (seconds); SPC-4; SPC-R5 |
| "Not measured is not detected" (`BR/docs/decisions/DR-009-measured-vs-detected.md`) | UNI-6, PGV-4 (`unquantified_policy`) |

### 3b. Changes that break the portal or force work

| Change | Effect | Code sites |
|---|---|---|
| **PGV-8:** the OpenMS member-protein minimum leaves `pg.global_qvalue`, which becomes null | Co-expression drops every protein of quantms DDA datasets, because null is read as 1.0. `ms.protein_fdr` warns. Explorer q-values go blank | `builder.py:219`; `ms_confidence.py:57-63`; `store.py:675, 944`; `queries.py:99` |
| **INT-4:** a QPX-computed pg sum leaves `pg.intensity`; the file declares `intensity_quantity = none` | If quantms pg files use the QPX sum, the sample map, comparison, heatmap, co-expression, per-sample QC and protein ranking all go empty. The value moves into nested `additional_intensities`, which the portal never reads | `FE/qpxAnalysis.js:55-57`; `qc.py:172-182`; `store.py:676, 997`; `queries.py:100-105`; `FE/views/QpxExplore.vue:443` |
| **CNF-2:** OpenMS PSM-level `feature.peptide_qvalue` moves to the new `psm_qvalue` | The explorer's feature q-value column, filter, sort and histogram go null for DDA too | `store.py:973-982, 1139, 1166, 1177`; `queries.py:109-124`; `BE/apps/qpx_explore/api.py:216`; `FE/views/QpxExplore.vue:212-214, 338, 789`; `FE/qpxFormat.js:8` |
| **CNF-5, CNF-6:** typed q-values leave `additional_scores`; bare `q-value` and `qvalue` are deprecated. This applies to DIA-NN files too | `ms.dia_precursor_q` warns on every 1.4 DIA dataset. `fdr_verifiable` loses its fallback, and its typed list lacks `precursor_qvalue`, `peptide_qvalue` and `psm_qvalue`. The atlas-FDR inputs break | `ms_confidence.py:20, 24-32, 66-76`; `ALGORITHM.md:22-23, 55-58`; `BR/tests/evaluation/test_ms_confidence.py:53-84` |
| **PGV-7:** `peptide_counts` becomes unit-scope | `max(peptide_counts…)` no longer gives the dataset total; the portal must use `peptides[]` or the feature table | `store.py:597-603, 665-666, 945` |
| **PGV-1, UNI-6, INT-5:** identification-only rows and null entries are added | "Quantified features" and group counts by row existence rise. Filters on non-null intensity and `identified_in_unit` are needed | `qpx_stats.py:226-235, 260`; `store.py:462, 471`; `ms_yield.py:15`; `qc.py:47` |
| **PGV-12:** descriptions leave `pg_names` | The `GN=` fallback fails; gene names depend on `gg_names` being filled | `builder.py:125-138`; `store.py:673, 942`; `FE/qpxAnalysis.js:139, 148` |
| **PRT-8:** shared peptides get null group and anchor | Protein attribution needs psm or pepmap, and pepmap is not publishable | `index_builder.py:542-544, 764`; `qpx_names.py:11` |
| **UNI-3:** new `unidentified_feature` view | Not publishable. If such a file ever reached a dataset folder, the explorer's suffix `feature.parquet` (no leading dot) would union it into `feature` | `qpx_names.py:27`; `store.py:46, 368` |
| **Footer keys** | New portal code is needed. The explorer's SQL sandbox blocks table functions, so frontend SQL cannot read footers; the backend must expose them in `/overview`. Multi-file `psm/` folders need the META-4 combine rule | `store.py:81, 452-489, 858-861`; `BE/evaluation/view.py:13`; `qpx_stats.py:144-147` |
| **Retyping `feature.pg_accessions` to `list<string>`** (2.0 open question) | The `x -> x.accession` lambdas break | `store.py:981, 1016, 1025`; `index_builder.py:764` |

These do not break the portal, because nothing reads them: the id changes, `unique` deprecation, `lfq` → `maxlfq`, `match_qvalue` type, nullable `psm.pass_threshold`, and nullable `pg.peptides`. One side effect remains: any re-conversion changes file bytes, so every dataset gets a new content-addressed reference and all indexes rebuild.

### 3c. Irrelevant or over-engineered for this consumer

- **T02** (`precursor_index`, CHI-*): the portal never counts per precursor.
- **T06** (identifying PSM, links, `identification_conflict`, `consensus_feature_id`): the portal never joins psm to feature.
- **T07** unidentified-feature view and its ids.
- **T11** `spectrum_reference`, the `scan_format` extension, `isotope_error`, ion-mobility units: `scan` is only displayed, and no USI is built.
- **T08** `match_*` definitions and MBR FDR parameters.
- **T13** validator severities, registry mechanics, ARROW-schema parity.
- All MaxQuant, FragPipe, PD, mzIdentML and CDAP producer rules.
- PSM-2 multi-engine merging; INT-10 to INT-13; PRT-9 and PRT-10.

Two rules run against how the portal has to work on mixed legacy data:
- SPC-R4 (exact run-name comparison) versus the fuzzy match at `ms_integrity.py:42-49`.
- FLT-R6 (typed contaminant flag only) versus the accession regexes.

They become usable only after `qpxc migrate` has been run over the published datasets.

### 3d. What the portal needs and the proposal does not give

1. **An accepted-row predicate in SQL without provenance.** FLT-R4 needs the provenance parameters `pass_threshold_level` and `pass_threshold_fdr`; FLT-R2, PGV-R5 and CNF-R3 also read provenance. The portal's templates, the ad-hoc SQL editor and its MCP tools cannot do that join. Feature rows have no acceptance column at all.
2. **A strict release profile.** The `dda` profile allows `unknown` and all ranks. The portal wants files where `fdr_filter` is a decimal, rows are rank-1 and accepted, and a run view exists, so that `WHERE NOT is_decoy` stays correct.
3. **`decoy_policy` where a reader can see it.** The planned atlas FDR pools decoys (`ALGORITHM.md:15-16, 34, 161-162`). That is a default behaviour depending on the value, so by the revision's own test it belongs in the footer.
4. **A stored dataset summary.** The portal recomputes KPIs with several scans at publish (`qpx_stats.py:103-295`). `qpx:qpx/converters/summary.py` only prints. The RDR-9 counts could be written to the dataset view or a sidecar.
5. **An always-filled, declared protein quantity** (the opposite of INT-4), and a typed group-level acceptance when `global_qvalue` is null (PGV-9 is only SHOULD).
6. **Cross-dataset protein identity.** `pg_group_id` is stable only for identical accession strings. There is no rule for filling `gg_names`, the accession namespace or database version, or per-accession taxon.
7. **Comparable quantities.** INT-R1 only forbids pooling. The portal avoids pooling anyway, so the two declaration keys are enough.
8. **The same contract for DIA.** The portal applies one SQL to both, and half its corpus is DIA.
9. **The portal as a tested reader.** T14's violation tables and RDR-14 fixtures cover only qpx and mokume.

## 4. Consumer-first subset, in priority order

1. `fdr_filter` (FLT-1, FLT-2, PGV-3), with decimals required under a release profile.
2. Typed q-values with one level each (CNF-1, CNF-2), kept filled: `psm.qvalue`, `feature.precursor_qvalue`, `feature.psm_qvalue`, `pg.global_qvalue`.
3. `is_decoy` never null, under one rule; typed `contaminant` on psm, feature and pg (PRT-4, PRT-5, PRT-6, PGV-10).
4. `acquisition_method` (META-5).
5. `psm.rank` and `pass_threshold` (PSM-3, PSM-4, FLT-7, FLT-9).
6. `identification_type` filled, and `pg.identified_in_unit` (MBR-1, MBR-2, MBR-3, MBR-10, PGV-1, PGV-2).
7. `intensity_quantity` and `intensity_processing`, null versus 0, imputation flags (INT-2, INT-6, INT-7, INT-9).
8. Protein key: `pg_group_id`, anchor as a member, bare accessions (PRT-3, PRT-11, PRT-13), plus filled `gg_names`.
9. Samples: a run view always present and exact `run.samples` resolution (META-6, RDR-8, PGV-5).
10. Feature row meaning and seconds (FRW-2, FRW-3, FRW-4, FRW-7), with the guide recommending `count(DISTINCT (peptidoform, charge, run_file_name))`.
11. Count definitions (RDR-9), stored as a summary.

What I would simplify:
- **Publish one-line SQL predicates per view** in the consumer guide, built only on typed columns and footer keys.
- **Do not remove a value consumers filter on without a typed replacement.** For PGV-8, write a group acceptance (`pass_threshold` plus the `protein` entry of `fdr_filter`). For INT-4, keep the sum in `pg.intensity` declared as `intensity_processing = sum`, with provenance saying QPX computed it.
- **Dual-write for one minor version** under CNF-6 (typed column plus the score entry), and keep `peptide_qvalue` readable until readers move to `psm_qvalue`.
- **Use four flat keys or a JSON value for `fdr_filter`** rather than a `;`-separated mini-language; existing `qpx.*` keys are already JSON.
- **Add a `release` profile** with no `unknown`, rank-1 rows and a declared cut; the portal's evaluation gate can run it at publish.
- **Defer** T02, T06 links, T07, most of T11 and the id changes until after this subset.

## Not verified

- Which route produced the portal's quantms files (`qpxc convert openms-consensus` or the OpenMS native `-out_qpx` writer, which the proposal did not audit), and their `qpx_version`. `BE/apps/coexpression/README.md` says at least two pg layouts are live.
- Whether live quantms `pg.intensity` is the OpenMS abundance or the QPX sum, and whether `pg.global_qvalue` is filled. These decide the real impact of INT-4 and PGV-8. One query on the pg `cv_params` value `quantification_method` would settle the first.
- Whether published files were relaid out into struct lists.
- MSNet tables' conformance; they are not produced by a qpx converter.
- The cell-line atlas build code, which lives in another repository.
- I did not run the portal's tests.
