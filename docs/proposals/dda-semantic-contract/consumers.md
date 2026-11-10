# Consumers: qpxui and the quantms.io portal

[Overview](index.md) · [Strict profile](strict-profile.md) · [OpenMS work list](openms-worklist.md)

!!! warning "Status: proposal, not normative"
    Part of the [DDA semantic contract proposal](index.md). This page checks the contract against how two consumers read QPX today.

The contract only helps if consumers can use it. This page is based on reading the source of qpxui (`bigbio/qpxui` at `ddeda63`) and of the quantms.io portal (`bigbio/quantms-portal` at `728f5e4` and its backend at `c562424`). Nothing was run against live data. The full reviews are in [review-qpxui.md](review-qpxui.md) and [review-portal.md](review-portal.md).

## How they read QPX

- **Both read typed columns with DuckDB and almost nothing else.** qpxui reads one footer key (`qpx.sort_order`). The portal reads none. Neither reads provenance, `qpx_version`, `rank`, `pass_threshold`, `identification_type` or `imputed_labels`. The portal's frontend SQL cannot read footers at all.
- **Their only row filter is "not a decoy".** Every non-decoy row is treated as an accepted identification. No q-value threshold is applied anywhere.
- **The portal ingests quantms (OpenMS) and quantmsdiann (DIA-NN) output and MSNet PSM bundles.** It ingests no MaxQuant, FragPipe or Proteome Discoverer files.

So a rule reaches these consumers only if it reduces to a `WHERE` clause on typed columns. That is why the [strict profile](strict-profile.md) fills every column and gives one predicate per question.

## What is wrong today

| Consumer assumption | What files actually hold | What the user sees |
|---|---|---|
| Every non-decoy psm row is an identification | OpenMS consensusXML conversions hold every hit of every identification, with no rank and no cut | PSM and peptide counts are upper bounds |
| Every feature row is identified in its run | Rows include match-between-runs transfers | Per-run peptide and protein counts are inflated |
| Rows of one precursor in one run can be summed | A direct and a transferred row can coexist | Quantities are too high |
| `feature.peptide_qvalue` is a peptide-level q-value | For OpenMS it is the PSM-level q-value of the first identification; for DIA-NN it is null | The q-value column is mislabelled or empty |
| A pg row means the group was identified | Rows can come from any hit, or rest only on transfers | Protein counts are inflated |
| `anchor_protein` identifies a protein | Different groups can share an anchor; shared peptides have no anchor | Groups are merged; shared peptides vanish from protein pages |
| Contaminants can be found in the accession text | Markers differ by producer and are sometimes stripped; psm and feature have no flag | Contaminant fractions differ between modules. The portal uses four different patterns |
| `is_decoy` is never null | It can be null; each code path treats that differently | Decoy fractions are understated |
| `pg.intensity` is filled and of a known kind | It is an OpenMS abundance or a QPX sum, undeclared | Heatmaps and comparisons run on an unknown quantity |
| `rt` and the peak window share a clock | Not for every producer | Retention charts draw windows that miss their point |

Two defects are in the consumers themselves and are independent of the contract:

- **qpxui reads the pg view in the 1.0 layout** (a `run_file_name` column and an `intensities` list). On pg files of 1.1 or later it silently loses protein quantities, per-run protein rows and group counts (`src/lib/queries/quant.ts:16-20`, `src/lib/queries/proteins.ts:21`).
- **The portal uses three different protein-group keys** (the raw accession list, the sorted list and the anchor), so one dataset shows three different group counts (`apps/qpx_explore/store.py:471`, `qpx_stats.py:228`, `apps/qpx_explore/store.py:951`).

## What consumers change

With files that follow the strict profile, each consumer needs about six predicates. The exact SQL is in the [strict profile](strict-profile.md#reading-these-files).

| Change | qpxui | Portal |
|---|---|---|
| Count PSMs and peptides with `NOT is_decoy AND qvalue <= :t` | `src/lib/queries/header.ts`, `overview.ts`, `ptms.ts` | `qpx_stats.py`, `apps/qpx_explore/queries.py` |
| Count identifications only from `direct` and `secondary` feature rows | `src/lib/queries/qc.ts` | `qpx_stats.py`, `apps/peptide_search/index_builder.py` |
| Key protein groups on `pg_group_id`; show `anchor_protein` | `src/lib/report.ts:222-228`, `queries/quant.ts` | `apps/qpx_explore/store.py`, `apps/coexpression/builder.py` |
| Use the typed `contaminant` column | `src/lib/queries/quant.ts:26-28` | `qpx_stats.py:22`, `evaluation/checks/bio.py:11`, `apps/qpx_explore/qc.py:54` |
| Count identified groups with `identified_in_unit` | `src/lib/queries/qc.ts:47` | `qpx_stats.py:226-235` |
| Read the q-value of a feature's PSM from `psm_qvalue` | `src/lib/queries/overview.ts:108-110` | `apps/qpx_explore/store.py:1139` |
| Resolve samples by exact `(run, label)` match | `src/lib/report.ts:337-346` | `apps/coexpression/builder.py:36-53` |
| Read pg in the current layout | `src/lib/queries/quant.ts`, `proteins.ts` | – |

## Changes that would have emptied consumer output

Four rules of the first draft removed a value that these consumers filter on. They were changed or need a replacement:

| Draft rule | Effect on consumers | Resolution |
|---|---|---|
| q-values leave `additional_scores` once typed columns exist | The portal's evaluation checks and its planned atlas FDR look scores up by name | Both are written during 1.4 ([CNF-6](t04-confidence-fields.md)) |
| A QPX-computed sum may not fill `pg.intensity` | Protein quantities go empty in both consumers for consensusXML conversions | Allowed when the tool gives none, declared as `intensity_processing = sum` ([INT-4](t09-intensities-missing-values.md)) |
| OpenMS's PSM-level value leaves `feature.peptide_qvalue` | q-value columns, filters and charts go empty | Consumers switch to `psm_qvalue`. The old column must not keep a value of the wrong level |
| OpenMS's member-protein minimum leaves `pg.global_qvalue` | The portal's co-expression reads null as 1.0 and drops every protein | The strict profile requires a group-level `global_qvalue` from OpenMS |

Which route produced the portal's current quantms files, and whether their `pg.intensity` is an OpenMS abundance or a QPX sum, could not be established. One query on those files would settle how urgent the second and fourth rows are.

## What the contract does not give consumers yet

1. **Stored dataset counts.** qpxui defers distinct counts on large files, and the portal recomputes its own with different rules, so one dataset shows different numbers in different tools. The writer could store the counts of spectra, PSMs, peptidoforms and protein groups in the dataset view.
2. **The sample on pg rows.** Resolving `(grouped_runs, label)` to a sample is awkward SQL. An optional `sample_accession` column on pg rows would make the protein matrix one `SELECT`.
3. **One encoding for name-value lists.** `cv_params`, `additional_scores` and `additional_intensities` exist in three physical layouts, which breaks lookups by name.
4. **The same contract for DIA.** Half the portal's datasets are DIA, and it applies one SQL to both.
5. **Protein identity across datasets.** `pg_group_id` is stable only for identical accession strings. Filling gene names and stating the database version would help the portal's cross-dataset views.
6. **Tested readers.** qpxui and the portal should run the same reader fixtures as qpx ([RDR-14](t14-reader-contract.md)).

## What these consumers do not need

Neither consumer uses the following, so they are deferred in the [order of work](index.md#order-of-work): the chimeric precursor index beyond counting spectra, the links from PSMs to features, conflict flags, the `unidentified_feature` view, the `match_*` details, most spectrum-field rules, and every rule that exists only for MaxQuant, FragPipe, Proteome Discoverer, mzIdentML or CDAP.
