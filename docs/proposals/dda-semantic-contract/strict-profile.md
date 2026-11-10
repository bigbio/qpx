# OpenMS strict profile

[Overview](index.md) · [Consumers](consumers.md) · [OpenMS work list](openms-worklist.md)

!!! warning "Status: proposal, not normative"
    This is the core of the [DDA semantic contract proposal](index.md), targeting specification 1.4.

This profile states what a DDA dataset **written by OpenMS** contains: ProteomicsLFQ and IsobaricWorkflow with `-out_qpx`, as used by the quantms pipeline. It is strict on purpose. Every column listed is filled as stated, and no footer value is `unknown`, so that qpxui, the quantms.io portal and mokume can read these files with plain SQL on typed columns, without special cases.

OpenMS can meet this profile because it holds the information internally; the [OpenMS work list](openms-worklist.md) names the export changes needed. Files from other producers (MaxQuant, FragPipe, Proteome Discoverer, mzIdentML, CDAP, and the Python consensusXML converter) follow the relaxed rules in the [annex](index.md#annex). Where the annex and this page differ, this page wins for OpenMS-written files.

## Dataset

| Requirement | Detail |
|---|---|
| Views | `psm`, `feature`, `pg`, `run`, `sample` and `provenance` are always written. |
| Runs and samples | Every `run_file_name` used in psm, feature and pg has a row in `run`. Every `(run_file_name, label)` of an intensity entry resolves to exactly one sample through `run.samples`, by exact string equality. |
| Unidentified features | Optional. OpenMS may write MS1 features without an identification as rows of the feature view, marked `identification_type = 'unidentified'`. Readers ignore them by default. See [feature rows](#feature-rows). |
| Version | `qpx_version = 1.4`. |

## Footer

| Key | Views | Allowed values in this profile |
|---|---|---|
| `acquisition_method` | psm, feature, pg | `dda` |
| `fdr_filter` | psm, feature, pg | `none`, or thresholds at `psm`, `peptide` and `protein`, for example `psm:0.01;protein:0.01`. Never `applied` or `unknown`. |
| `unquantified_policy` | feature, pg | `retained` |
| `intensity_quantity` | feature, pg | `ms1_area`, `reporter_corrected` or `reporter_uncorrected` |
| `intensity_processing` | feature, pg | feature: `none` or `normalized`. pg: `topn`, `sum` or `normalized`. |
| `scan_format` | psm, feature | `scan`, `index` or `nativeId` |
| `identity_composite` | feature | Label-free: `peptidoform, charge, run_file_name`. Isobaric: `peptidoform, charge, run_file_name, scan`. |

The footer has no key that names the profile. A dataset claims it in the `provenance` view instead: the step that wrote the QPX files carries the parameter `conformance_profile = dda-strict`. Provenance is an ordinary table, so consumers read the claim with a plain `SELECT`. A merge keeps the claim only if every input dataset has it. A reader that needs a guarantee runs `qpxc validate --profile dda-strict`.

## What the profile fixes, so no declaration is needed

| Fact | Value for OpenMS-written files |
|---|---|
| PSM rank | Only the best hit per precursor is exported. `rank` is always 1. |
| FDR levels | OpenMS always estimates PSM-level FDR, so PSM q-values are always present. Peptide-level FDR (common for immunopeptidomics) and protein-level FDR are estimated when the analysis asks for them. A peptide or protein q-value column is filled exactly when FDR was estimated at its level, and is null on every row otherwise. |
| Scope of PSM q-values | One estimate over the whole experiment. |
| Protein grouping | One inference over the whole dataset, so a group is the same in every run. |
| Primary intensity | The measured signal: the feature area for label-free data, the reporter-ion intensity for isobaric data. Protein quantities are OpenMS's own. |
| Imputation | OpenMS imputes nothing. A value that OpenMS's own algorithm estimates, such as the area it derives when the elution-model fit fails, counts as measured. `imputed_labels` is `[]` and pg `is_imputed` is false. A tool that fills in values afterwards, such as mokume, flags them as imputed. |
| Precursor-level FDR | OpenMS has none. `feature.precursor_qvalue` and `feature.global_qvalue` are null by definition. |

Settings that vary between analyses go to the `provenance` view under fixed parameter names: `decoy_policy`, the FDR settings, whether match-between-runs ran and with which method, and `mbr_fdr` (a threshold under PIP-ECHO, `none` otherwise).

## psm rows

One row is one identified precursor of one MS/MS spectrum.

| Column | Rule |
|---|---|
| `run_file_name`, `scan`, `peptidoform`, `charge` | Unique together. `charge ≥ 1`. `scan` holds one spectrum, and `spectrum_reference` holds the native id. |
| `rank` | Always 1. |
| `precursor_index` | Never null. 1 for the spectrum's first precursor; 2, 3, … for co-isolated precursors identified in the same spectrum ([CHI-3](t02-chimeric-spectra.md)). |
| `qvalue` | Never null. The PSM-level q-value; OpenMS always estimates PSM-level FDR. |
| `posterior_error_probability` | Never null. |
| `is_decoy` | Never null. True only if every protein the peptide maps to is a decoy. |
| `contaminant` | Never null. |
| `observed_mz` | The measured precursor m/z, or null. Never 0 and never a theoretical value. |
| `rt` | Seconds on the run's acquisition clock, not the aligned clock. |
| `protein_accessions` | Every protein the peptide maps to. No decoy entries on target rows. |
| `feature_id` | The feature this PSM supports. Null only when it supports none. |

## feature rows

One row is one quantity of one precursor (peptidoform and charge) in one run.

- **Label-free:** at most one quantified row per `(peptidoform, charge, run_file_name)`. OpenMS does not write a transferred row beside a quantified direct row of the same precursor and run, and records in provenance how many it left out.
- **Isobaric:** one row per quantified reporter spectrum. Readers add the rows of a precursor in a run.

The rules below apply to identified rows, that is every row except the optional unidentified ones described after the table.

| Column | Rule |
|---|---|
| `peptidoform`, `sequence`, `charge` | Never null or empty. Canonical peptidoform. `charge ≥ 1`. |
| `identification_type` | Never null: `direct`, `secondary` or `transferred` ([MBR-1](t08-match-between-runs.md)). The transfer status comes from OpenMS's own annotation of transferred runs. |
| `scan`, `id_run_file_name`, `posterior_error_probability`, `psm_qvalue` | On `direct` and `secondary` rows: exactly one spectrum, `id_run_file_name = run_file_name`, and the PEP and `psm_qvalue` of that PSM, all non-null. On `transferred` rows: `scan` is empty and the others are null. |
| `peptide_qvalue` | Filled whenever peptide-level FDR was estimated. |
| `match_qvalue` | On `transferred` rows, non-null exactly when PIP-ECHO ran. Null on other rows. |
| `intensities` | One entry per label that has a sample in the run. A null intensity means not quantified; 0 is never used for missing. A `transferred` row has at least one measured value. |
| `rt`, `rt_start`, `rt_stop` | Seconds on the acquisition clock, with `rt_start ≤ rt ≤ rt_stop` (label-free). `rt_aligned` holds the aligned time. |
| `observed_mz` | The measured m/z of the quantified signal, or null. Never the theoretical value. |
| `is_decoy`, `contaminant` | Never null. |
| `pg_accessions`, `anchor_protein`, `pg_group_id` | The full membership of the one protein group the peptide is attributed to, its representative, and its id. All three are null together when the peptide is attributed to no group. |
| `group_unique` | Never null. True when the peptide is evidence for exactly one group. |
| `consensus_feature_id`, `identification_conflict` | Never null. Rows of one cross-run feature share the id. |
| Isobaric rows | `quant_scan`, `quant_ms_level` and `parent_ion_fraction` are filled, and `quant_spectrum_count = 1`. |

### Unidentified features (optional)

OpenMS may also write MS1 features that have no identification. They are rows of the same feature view ([T07](t07-unidentified-features.md)):

- `identification_type = 'unidentified'`.
- `peptidoform` and `sequence` are null, never empty strings. `charge` is the feature's charge, or null when none was assigned, never 0.
- `is_decoy` is false. The scan, PEP, q-values, `match_*`, `contaminant` and all protein columns are null.
- Intensities, retention times, `observed_mz` and `consensus_feature_id` follow the rules above.
- No PSM links to such a row, and it never counts as evidence for a peptide or protein.

Readers ignore these rows by default. The test is `peptidoform IS NOT NULL`.

## pg rows

One row is one protein group in one quantification unit (the fractions of one sample) for one label.

| Column | Rule |
|---|---|
| `pg_accessions`, `pg_group_id`, `anchor_protein` | The group's full membership, its id, and a representative that is a member. |
| `grouped_runs`, `label` | Never null. The runs of the unit, and the label. |
| `identified_in_unit` | Never null. True when a PSM from a run of this unit, not a transfer, supports the group. |
| `intensity` | OpenMS's own protein abundance. Null means identified but not quantified. |
| `global_qvalue` | The group's q-value over the experiment. Filled whenever protein-level FDR was estimated; null on every row otherwise. |
| `pg_qvalue` | Null by definition; OpenMS has no per-unit protein q-value. |
| `is_decoy` | Never null. True only if every member is a decoy. |
| `contaminant` | Never null. |

## Reading these files

Each question a consumer asks is one `WHERE` clause on typed columns. `:t` is the reader's threshold, usually 0.01.

| Question | Predicate |
|---|---|
| Is this PSM an accepted identification? | `NOT is_decoy AND qvalue <= :t` |
| Was this precursor identified in this run? | `NOT is_decoy AND identification_type IN ('direct', 'secondary') AND psm_qvalue <= :t` |
| Does this feature row carry a usable quantity? | `NOT is_decoy AND peptidoform IS NOT NULL`, and an entry with a non-null intensity. |
| Was this protein group identified in this unit? | `NOT is_decoy AND identified_in_unit AND global_qvalue <= :t` |
| Is this row a contaminant? | `contaminant` |

Three points go with these predicates:

- **A q-value level that was not estimated cannot be filtered on.** If a q-value column is null on every row, the reader filters at a level that exists and says so. For an analysis without protein-level FDR, that means no q-value condition on protein groups; the PSM predicates always work. A reader never treats null as 1.0 and never drops every row silently.
- **Transferred quantities are kept by default.** They improve quantification. Leaving them out is the reader's choice: add `identification_type <> 'transferred'`, or under PIP-ECHO keep those with `match_qvalue <= :t`.
- **Protein groups are keyed on `pg_group_id`.** `anchor_protein` is for display only.

### Reference queries (DuckDB)

These were run on small synthetic tables that follow this profile.

```sql
-- Identified spectra, PSMs and peptidoforms at 1% PSM FDR
SELECT count(DISTINCT (run_file_name, scan))                  AS identified_spectra,
       count(DISTINCT (run_file_name, scan, precursor_index)) AS psms,
       count(DISTINCT peptidoform)                            AS peptidoforms
FROM psm
WHERE NOT is_decoy AND qvalue <= 0.01;

-- Precursors identified per run (transfers are not identifications)
SELECT run_file_name,
       count(DISTINCT (peptidoform, charge)) AS identified_precursors
FROM feature
WHERE NOT is_decoy
  AND identification_type IN ('direct', 'secondary')
  AND psm_qvalue <= 0.01
GROUP BY run_file_name;

-- Precursor quantities per run, with the sample attached
SELECT f.peptidoform, f.charge, f.run_file_name, s.sample_accession,
       sum(i.intensity)                                AS intensity,
       bool_and(f.identification_type = 'transferred') AS transferred_only
FROM feature f, unnest(f.intensities) AS t(i),
     run r, unnest(r.samples) AS u(s)
WHERE NOT f.is_decoy
  AND f.peptidoform IS NOT NULL
  AND i.intensity IS NOT NULL
  AND r.run_file_name = f.run_file_name
  AND s.label = i.label
GROUP BY ALL;

-- Protein groups identified per unit at 1% protein FDR
-- (drop the q-value condition when the file has no protein q-values)
SELECT grouped_runs, count(DISTINCT pg_group_id) AS identified_groups
FROM pg
WHERE NOT is_decoy AND identified_in_unit AND global_qvalue <= 0.01
GROUP BY grouped_runs;

-- Protein quantities per sample
SELECT p.pg_group_id, any_value(p.anchor_protein) AS anchor_protein,
       s.sample_accession, sum(p.intensity) AS intensity
FROM pg p, run r, unnest(r.samples) AS u(s)
WHERE NOT p.is_decoy AND NOT p.contaminant
  AND p.intensity IS NOT NULL
  AND r.run_file_name = p.grouped_runs[1]
  AND s.label = p.label
GROUP BY p.pg_group_id, s.sample_accession;
```

For precursor quantities per sample, add the runs that are fractions of one sample. Technical replicates of one sample are combined by a method the reader names, not by summing.

## Transition rules

Two annex rules were changed so that existing consumers keep working:

- **Typed columns and score entries are both written in 1.4.** q-values and `rank` fill their typed columns and also stay in `additional_scores` with the same value. The entries are removed in a later minor version ([CNF-6](t04-confidence-fields.md)).
- **A protein quantity that QPX computes may fill `pg.intensity`**, but only when the tool provides none in the whole file, and then with `intensity_processing = sum` ([INT-4](t09-intensities-missing-values.md)). This does not arise for OpenMS-written files, which carry OpenMS's own abundance.

Two changes need consumer work:

- OpenMS's PSM-level q-value moves from `feature.peptide_qvalue` to `feature.psm_qvalue`.
- `peptidoform`, `sequence` and `charge` can be null in the feature view, on unidentified rows only.

See [Consumers](consumers.md).

## Not in this profile yet

These need algorithm work in OpenMS, or are not needed by the current consumers, and are deferred:

- keeping rank-1 ties as separate PSM rows;
- a precursor-level FDR;
- an FDR for transfers made without PIP-ECHO;
- transfer scores and the retention-time and m/z differences to the donor.

## Decisions taken

The maintainers decided these points in [#383](https://github.com/bigbio/qpx/issues/383):

| Point | Decision |
|---|---|
| A transferred row beside a quantified direct row of the same precursor and run | OpenMS leaves the transferred row out. |
| Unidentified features | They may be stored in the feature view, optionally, marked `identification_type = 'unidentified'`. Readers ignore them. There is no separate view. |
| FDR levels | PSM-level FDR is always estimated, so PSM q-values are always present. Protein-level FDR is not: OpenMS sometimes stops at peptide-level FDR, for example for immunopeptidomics. A peptide or protein q-value column is null when its level was not estimated. This holds for every producer. |
| Transferred quantities | Kept by default. Filtering them out is the reader's choice. |
| Estimated values | A value the quantification tool's own algorithm estimates counts as measured; that is the tool's decision. A value added afterwards by another tool counts as imputed. |
| Profile claim | A provenance parameter, not a footer key. |

## To verify in OpenMS

- What OpenMS's protein quantification does with a direct and a transferred feature of the same precursor in one run, so that the feature rows still add up to the protein quantity after the transferred row is left out.
- That ProteomicsLFQ never writes two direct features for one precursor in one run.
