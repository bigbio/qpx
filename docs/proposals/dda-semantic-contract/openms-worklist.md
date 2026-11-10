# OpenMS work list

[Overview](index.md) · [Strict profile](strict-profile.md) · [Consumers](consumers.md)

!!! warning "Status: proposal, not normative"
    Part of the [DDA semantic contract proposal](index.md). This page lists what the OpenMS QPX export needs to change to meet the [strict profile](strict-profile.md).

OpenMS is the producer the maintainers control, so the contract is written for it first. This page is based on reading the OpenMS C++ source: the default branch at `0707508` (2026-11-09) and the unmerged branch `claude/qpx-issue-377-review-dfv5gg` at `811c466`, which holds the work of [#378](https://github.com/bigbio/qpx/issues/378). Nothing was executed; no natively written fixture exists in qpx yet. The full review is in [review-openms.md](review-openms.md).

The work is tracked in the OpenMS repository as [OpenMS/OpenMS#10869](https://github.com/OpenMS/OpenMS/issues/10869).

## Summary

- **OpenMS already holds almost everything the contract asks for.**
    - About 30 of the rules that apply to OpenMS are met today, several with caveats.
    - 11 more are met on the #378 branch.
    - About 65 need only a small export change, because the data is in memory.
    - 4 need algorithm work: keeping rank-1 ties, a precursor-level FDR, an FDR for transfers made without PIP-ECHO, and transfer scores with the differences to the donor.
- **The native writer is the route to invest in.** qpx currently marks it as deprecated in favour of the Python consensusXML converter (`qpx/converters/openms/converter.py:402`, `docs/guide/convert.md:919`). The deprecation notice cites two exporter defects, PSMs assigned to the first run and conflicting identities (OpenMS#9872, OpenMS#9871). Fixing them belongs to item 1 below, after which the deprecation should be reversed.
- **The consensusXML route cannot meet the contract on three points**, because ProteomicsLFQ discards the data before the consensusXML is written:
    - retention times on the acquisition clock (they are aligned in place, and the original is not kept);
    - the bounds of a feature's elution peak;
    - the measured precursor m/z of a label-free feature (the stored value is theoretical).

## What the native writer does today

Source files are under `src/openms/source/FORMAT/` in OpenMS (`QPXFile.cpp`, `ConsensusMapArrowExport.cpp`, `ProteinGroupArrowExport.cpp`, `ArrowIOHelpers.cpp`, `ArrowSchemaRegistry.cpp`). `ProteomicsLFQ -out_qpx` and `IsobaricWorkflow -out_qpx` write three files: feature, psm and pg. Files are stamped `qpx_version = 1.1`.

| View | State at the default branch |
|---|---|
| psm | One row per best hit. There is no `rank`, `qvalue` or `spectrum_reference` column. **The PSM q-value is not written at all**: after FDR estimation the main score is the q-value, and the writer exports the main score only when it is a PEP. A missing `observed_mz` is written as 0.0. Score names are raw OpenMS keys. |
| feature | Peptide, PEP, `scan` and `id_run_file_name` come from the best hit across all runs of the consensus feature, so a row can name another run. Unidentified features are rows with an empty peptidoform. There is no `identification_type`, `match_*`, `consensus_rt` or `rt_aligned`. `rt_start` and `rt_stop` are `rt ± width/2`. `pg_accessions` holds every protein the peptide maps to, `unique` means one accession, and `anchor_protein` is the first accession. Intensity entries cannot be null. |
| pg | A sound base: OpenMS's own abundance per fraction group and label, identification-only rows, and an id over the full membership. But evidence counts every run in which the feature was quantified, transfers included; `is_decoy` comes from the anchor only; `contaminant` is always null. |
| Other | These two tools write no run, sample, provenance or ontology view, so the dataset says nothing about its FDR cuts or intensity kind. |

The #378 branch adds a per-run identifying PSM, links from PSMs to features of their own run, `identification_type`, `match_qvalue` and `consensus_rt`. It still drops the q-value and still writes empty peptidoforms.

OpenMS already has the code for the companion views: the search engine ProSE writes run, sample and provenance views through `QPXDatasetExport` (`src/topp/ProSE.cpp:2833-2842`). ProteomicsLFQ and IsobaricWorkflow do not call it yet.

Where OpenMS states a transfer: the consensus feature annotation `mbr_transfer_map_indices` lists the runs an identification was transferred to, and `mbr_transfer_qvalues` holds PIP-ECHO's q-value per transfer. Neither the native writer nor qpx reads them today.

## Changes, in priority order

| # | Change in the OpenMS export | Benefit for consumers | qpx code this retires |
|---|---|---|---|
| 1 | Merge the #378 branch. Take the transfer status from `mbr_transfer_map_indices`, order candidate PSMs by PEP, and write `match_qvalue` as float64. | Transfers can be told from identifications. Fixes rows that name another run's scan and PEP. | The transfer heuristic and recovery of removed identifications in `openms_consensus/feature_adapter.py:523-634`; `record_linking.py`; `openms/psm_runs.py` |
| 2 | Always export the PSM q-value to `psm.qvalue` and `feature.psm_qvalue`; PSM-level FDR is always estimated. Export the peptide q-value to `feature.peptide_qvalue` and the protein-group q-value to `pg.global_qvalue` when FDR was estimated at that level; otherwise the column stays null. Use snake_case score names. | Confidence filters work. No native file has a typed q-value today. | `psm_adapter.py:204-246`; `feature_adapter.py:109-206` |
| 3 | Keep the original retention time and the peak bounds per feature. Write `rt` and `psm.rt` on the acquisition clock, and `rt_aligned` beside them. | Chromatogram and spectrum lookups in qpxui; `rt` lies inside its window. | Nothing: this cannot be done from consensusXML. |
| 4 | Write the run, sample, provenance and ontology views from the experimental design and the tool parameters, reusing `QPXDatasetExport` as ProSE does, and stamp `qpx_version` 1.4. The run view carries `acquisition_method`. Provenance carries the FDR cuts per level, `unquantified_rows`, `decoy_policy`, `contaminant_policy`, `match_between_runs` (`true` or `false`) with `mbr_fdr`, `pg_grouping = dataset_inference` and `conformance_profile = dda-strict`. The ontology view carries the properties of `intensity` (quantity, processing), of each q-value column (level, scope) and of `observed_mz` (calibration). No new footer keys. | Consumers read the filter state and intensity kind with a plain `SELECT`. Samples resolve without an SDRF. | Most of `openms/converter.py`. Only SDRF sample metadata remains. |
| 5 | Stop writing sentinels: unidentified features, when written at all, are rows with `identification_type = 'unidentified'` and a null peptidoform, sequence and (if unassigned) charge; identified rows have `charge ≥ 1`; `observed_mz` is measured or null (for label-free features the measured value can be rebuilt from the theoretical m/z and the `masserror_ppm` deviation, which today survives ProteomicsLFQ only with PIP-ECHO); intensity entries can be null; reporters below the minimum become null instead of 0. | No empty peptidoforms or false zeros in tables and matrices. | The handling of four converter defects |
| 6 | Protein fields. Feature: `pg_accessions` as the group's membership, `group_unique`, the group's anchor, `contaminant` (needs a contaminant-marker parameter). pg: `is_decoy` over all members, `identified_in_unit`, evidence only from identifications of the unit's own runs. | Group-specific peptide filters are correct. Protein rows that rest only on transfers are visible. | `feature_adapter.py:455-490`; the evidence and decoy logic in `pg_adapter.py` |
| 7 | psm: `rank = 1`, `precursor_index`, `spectrum_reference`, `isotope_error`, `compensation_voltage`; `scan_format = nativeId` with multi-part ids. | Correct counts for chimeric spectra; timsTOF and FAIMS data. | The hashed scans in `psm_adapter.py:49-71` |
| 8 | Isobaric rows: `quant_scan`, `quant_ms_level`, `parent_ion_fraction`, `quant_spectrum_count = 1`. | Purity filtering. | `quant_scan.py` |
| 9 | `identification_conflict`, `consensus_feature_id`, `spectral_count`. | Traceability across runs. | `feature_dedup.py` |
| 10 | Guarantee one feature row per precursor and run for label-free data. ProteomicsLFQ already produces one feature per precursor in a run. If the internal data ever holds more than one candidate, the export decides which one to write, consistently with the protein quantification. | Readers that sum or count rows are correct. | – |
| 11 | Identification-only feature rows from identifications that support no quantified feature, with the provenance parameter `unquantified_rows = retained`. | "Identified but not quantified" can be told from "not identified". | – |
| 12 | MS1-label multiplets (SILAC, dimethyl, mTRAQ). Write the label state to the typed columns `feature.identified_labels` and `psm.label` instead of `cv_params`. In a labelled run `identified_labels` is never null: `[]` when no partner was identified by MS/MS, as on every transferred row. For dimethyl and mTRAQ, write the lightest label inline in the feature peptidoform; the export currently writes the label-free peptide identity there, so `calculated_mz` lacks the label mass and the mass error is left null. | Ratios and per-channel identification state without parsing `cv_params`; calculated and measured m/z describe the same ion. | – |
| 13 | Later, algorithm work: a transfer score and the differences to the donor from PIP-ECHO; an FDR for transfers made by the linker; precursor-level q-values; rank-1 ties. | Fills `match_*` and `precursor_qvalue`. | – |

**Effect on qpx.** `qpx/converters/openms_consensus/` (about 4,100 lines) becomes a legacy path for existing consensusXML files. `qpx/converters/openms/` (about 750 lines) shrinks to validation and adding SDRF sample metadata. `qpx/converters/quantms_msstats/` (about 710 lines) is not needed for OpenMS-produced data.

The maintainers take ProteomicsLFQ to produce one feature per precursor in a run. The strict profile relies on that; it was not checked against the source here.

## Not verified from source

- Whether a protein group's probability after protein FDR is a group-level q-value, as the writer assumes.
- PIP-ECHO internals: whether a transfer score or differences to the donor are stored, and whether decoy transfers can reach the output.
- MS1-label handling inside MS1LabeledWorkflow and FeatureFinderMultiplex. The export code was read: it writes the label-free peptide identity on feature rows, the matched peptidoform on psm rows, and the label state as `cv_params`.
- How OpenMS records which group a shared peptide is assigned to.
- Whether search adapters set a chimeric hit's m/z to the co-isolated ion, and whether that m/z was recalibrated.
- How FDR cuts applied before ProteomicsLFQ in the quantms pipeline are recorded.
