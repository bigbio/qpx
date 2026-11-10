# DDA semantic contract for QPX

**Proposal for specification 1.4**

!!! warning "Status: proposal, not normative"
    Nothing here changes the specification until the maintainers adopt it. The analysis is of qpx `dev` at commit `32cbbdc` (specification 1.3). Discussion: [#383](https://github.com/bigbio/qpx/issues/383).

QPX fixes column names, types and nullability. It does not yet fix what the rows and values *mean*. Two DDA files can pass validation with the same schema and still hold different things, and a reader cannot tell which from the file. This proposal defines a semantic contract so that everyone reading a DDA QPX file knows what its columns and entries mean.

## Priorities

The contract is built around two priorities set by the maintainers:

1. **Consumers must work well with QPX.** qpxui and the quantms.io portal read typed columns with SQL and little else. A rule counts only if it lets them get correct numbers that way.
2. **OpenMS is the producer that matters most.** The maintainers control its export, so it can be made exact. Support for MaxQuant, mzIdentML and other converters comes second.

## How the proposal is organised

| Part | Page | Content |
|---|---|---|
| Core | [OpenMS strict profile](strict-profile.md) | What an OpenMS-written DDA dataset contains: views, footer, the rules for psm, feature and pg rows, and the SQL predicates for reading them. Nothing is `unknown`. |
| Core | [Consumers](consumers.md) | How qpxui and the portal read QPX, what is wrong today, and what they change. |
| Core | [OpenMS work list](openms-worklist.md) | What the OpenMS export changes, in priority order. |
| Annex | [Annex overview](annex-overview.md) | The producer-neutral analysis: all footer and column proposals, breaking changes, producer status and open questions. |
| Annex | Topic sections [T01](t01-psm-rows.md) to [T14](t14-reader-contract.md) | Detailed analysis and numbered rules per topic, including the relaxations other producers need. |
| Annex | Findings: [PSM](findings-psm.md), [feature](findings-feature.md), [protein group](findings-pg.md), [cross-view](findings-cross.md) | 443 verified findings with evidence. |
| Annex | [Converter worklist](converter-worklist.md) | 174 converter and tooling defects, each mapped to a rule. |
| Annex | Reviews: [qpxui](review-qpxui.md), [portal](review-portal.md), [OpenMS](review-openms.md) | The three reviews behind the core pages. |

For OpenMS-written files the strict profile is the contract, and it wins where it differs from the annex.

## What the audit found

Each number was measured on the repository's own fixtures with qpx `dev`, unless marked.

- **Readers return nothing or the wrong thing.** `feature.confident()` keeps 0 of 243,821 MaxQuant features and `psm.confident()` keeps 0 of 13,386 Mascot PSMs, because no typed q-value was written.
- **Rows are not what consumers assume.** A Mascot file gives 13,386 PSM rows for 9,707 spectra. 43% of MaxQuant feature rows are match-between-runs transfers, and in 9,272 target precursor-run pairs a measured direct row and a measured transferred row coexist.
- **One column, several meanings.** `feature.precursor_qvalue` is a precursor-level FDR for DIA-NN, the best PSM's q-value for Proteome Discoverer and CDAP, and empty elsewhere. For OpenMS, `feature.peptide_qvalue` holds a PSM-level value.
- **Files do not say how they were filtered.** No footer states the FDR cut, or whether decoys or transfers were kept.
- **The OpenMS native export has gaps** (read from the OpenMS source, not measured): the PSM q-value is not written at all, a feature row can name another run's scan and PEP, unidentified features are rows with an empty peptidoform, and retention times are aligned without the original being kept.
- **Consumers filter only on "not a decoy"** (read from their source): qpxui and the portal treat every non-decoy row as an accepted identification and apply no q-value threshold.

## The four guiding questions

These are the answers under the [strict profile](strict-profile.md). The annex has the general rules.

**How many identifications per feature?**
One feature row is one precursor (peptidoform and charge) in one run. A `direct` or `secondary` row has exactly one identifying PSM from its own run, and that PSM supplies the row's `scan`, PEP and `psm_qvalue`. Other PSMs of the same precursor in that run link to the row through `psm.feature_id`. A `transferred` row has no PSM.

**What if the spectrum is chimeric?**
Each co-isolated precursor identified in a spectrum is its own psm row. The rows share `run_file_name` and `scan` and differ in `precursor_index` (1, 2, …). Each has rank 1. A chimeric spectrum's reporter ions feed at most one feature.

**What if a feature has no identification?**
By default it is not written. A producer may store such MS1 features in the feature view, marked `identification_type = 'unidentified'` with a null peptidoform, and readers ignore those rows. An identified precursor without a quantity is a feature row whose intensity is null.

**How are match-between-runs and identification in the same run told apart?**
By `identification_type`, which is never null: `direct` or `secondary` means identified by MS/MS in this run; `transferred` means the signal was measured here but the identity came from another run. Transferred rows have no scan, PEP or `psm_qvalue`; their own confidence is `match_qvalue`. They never count as identifications or as protein evidence. On protein rows, `identified_in_unit` says whether the group has MS/MS evidence in that unit.

## Proposed footer keys

A key is in the footer only if a reader's default behaviour depends on its value. That leaves five new keys.

| Key | Views | Values in the strict profile | Why readers need it |
|---|---|---|---|
| `acquisition_method` | psm, feature, pg | `dda` | Decides whether the DDA or DIA rules apply. |
| `fdr_filter` | psm, feature, pg | `none`, or thresholds per level such as `psm:0.01;protein:0.01` | Says which rows were already cut. |
| `unquantified_policy` | feature, pg | `retained` | A missing row then means "not identified". |
| `intensity_quantity` | feature, pg | `ms1_area`, `reporter_corrected`, `reporter_uncorrected` | Decides whether files can be pooled, and whether reporters may be corrected again. |
| `intensity_processing` | feature, pg | `none`, `normalized`, `topn`, `sum` | Same reason. |

No footer key names the profile. A dataset claims it with the provenance parameter `conformance_profile = dda-strict`, which lets a reader skip every "unknown" branch. Other producers may also write `unknown` and further values; see the [annex overview](annex-overview.md#proposed-footer-keys), which also lists what replaced the roughly 25 keys of the first draft.

## Order of work

1. **OpenMS export.** Merge the [#378](https://github.com/bigbio/qpx/issues/378) branch, then export the PSM q-value. These two give consumers correct identification counts and a working confidence filter. The remaining steps are in the [OpenMS work list](openms-worklist.md), tracked as [OpenMS/OpenMS#10869](https://github.com/OpenMS/OpenMS/issues/10869).
2. **Consumer fixes that do not depend on new files.**
    - qpxui: read the pg view in its current layout.
    - Portal: key protein groups one way.
    - Both: stop treating a null `is_decoy` as a target.
3. **Footer plumbing in qpx.** The five keys, kept through merge and migrate ([T13](t13-metadata-conformance.md)).
4. **Strict-profile validation.** `qpxc validate --profile dda-strict` and the provenance claim, with fixtures written by OpenMS, and the reference queries as shared reader tests for qpx, qpxui, the portal and mokume.
5. **Consumer predicates.** Each consumer adopts the predicates of the [strict profile](strict-profile.md#reading-these-files).
6. **Other producers.** The converter work in the annex, by producer.

Deferred until a consumer needs them: writing unidentified features, conflict flags, the `match_*` details beyond `match_qvalue`, most spectrum-field rules, and the rules that exist only for other producers.

## Annex

The annex keeps the full producer-neutral analysis. It matters for three reasons:

- **Other producers.** Files converted from MaxQuant, FragPipe, Proteome Discoverer, mzIdentML, CDAP or consensusXML cannot always fill every column. The annex says what they write instead: null, or `unknown` in the footer.
- **Evidence.** Every rule in the strict profile traces back to findings there.
- **Readers of mixed data.** The reader rules in [T14](t14-reader-contract.md) say how to handle files that are not strict.

Start with the [annex overview](annex-overview.md).

## Decisions taken

The maintainers decided these points in [#383](https://github.com/bigbio/qpx/issues/383). They are worked into the strict profile and the annex.

| Point | Decision |
|---|---|
| A transferred row beside a quantified direct row of the same precursor and run | Not written. The contract requires one row per precursor and run; how OpenMS resolves such a pair internally is decided in the OpenMS export. |
| Unidentified features | They may be stored in the feature view, optionally, marked `identification_type = 'unidentified'`. Readers ignore them. There is no separate view. This keeps [#229](https://github.com/bigbio/qpx/issues/229)/[#232](https://github.com/bigbio/qpx/issues/232) and settles [#378](https://github.com/bigbio/qpx/issues/378). |
| FDR levels | PSM-level FDR is always estimated, so PSM q-values are always present. Protein-level FDR is not: OpenMS sometimes stops at peptide-level FDR, for example for immunopeptidomics. A peptide or protein q-value column is null when its level was not estimated. This holds for every producer. |
| Transferred quantities | Kept by default, because they improve quantification. Filtering them out is the reader's choice. |
| Estimated values | A value the quantification tool's own algorithm estimates counts as measured; that is the tool's decision. A value added afterwards by another tool, such as mokume, counts as imputed. |
| Features per precursor | ProteomicsLFQ is taken to produce one feature per precursor in a run, so label-free files have one feature row per precursor and run. |
| MS1-label multiplets | One feature row per multiplet with the light peptidoform and one intensity entry per channel. PSMs keep their labelled peptidoform. |
| Profile claim | A provenance parameter, not a footer key. |
| PSI alignment | Not a goal. PSI terms are cited only where they pin down a meaning. |
| Footer | Five new keys. Everything else is a typed column, a fixed rule or a provenance parameter. |

The open points that remain concern other producers; they are in the [annex overview](annex-overview.md#open-questions-for-maintainers).

## How this was produced

1. **Audit.** Eleven audits of the specification, the DDA converters, the writer and validator, the qpx read API, mokume, the GitHub history and external tools. A second pass tried to refute every finding: 443 survived, 119 of them with a correction.
2. **Design.** Fourteen topic sections, each drafted, reviewed adversarially and harmonised.
3. **Footer reduction.** From about 30 proposed keys to five.
4. **Consumer and OpenMS review.** qpxui, the portal and the OpenMS C++ export were read against the proposal, which led to the strict profile and this structure.

Limits:

- The consumer and OpenMS reviews are from reading source. Nothing was run against live portal data, and no OpenMS-written fixture exists in qpx.
- It is not established which route produced the portal's current quantms files.
- Proteome Discoverer and FragPipe behaviour comes from converter code and synthetic cases.
