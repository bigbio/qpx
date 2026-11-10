# T01. PSM rows: grain, ranks, alternatives, ties, multi-engine hits and psm_id uniqueness

[Overview](index.md) · [Strict profile](strict-profile.md) · [Annex overview](annex-overview.md) · [T02 →](t02-chimeric-spectra.md)

!!! note "Annex: producer-neutral analysis"
    Part of the annex of the [DDA semantic contract proposal](index.md) for QPX. Rules here are proposals for specification 1.4, not normative text, and they cover every producer. For OpenMS-written files the [strict profile](strict-profile.md) is the contract and wins where the two differ. Bracketed ids such as [OPENMS-01](findings-psm.md#openms-01) link to the verified audit findings; code references are to `dev` at commit `32cbbdc`.


> **Decision in one line:** One psm row per hypothesis per spectrum, unique on K = `(run_file_name, scan, peptidoform, charge)` with engines merged; a typed `rank` (the source's rank per precursor ion, ties share 1; null = unknown) and no footer key; rank stays out of the identity; only rank-1 rows are identifications.

## Questions a reader has
- How many rows per spectrum, and what tells them apart? Is a scan with two charges two rows?
- Are rank-2+ alternatives present? Where is rank, and where is that choice declared?
- Is rank per spectrum, per precursor or per charge query? Are decoys ranked? How are rank-1 ties (I/L, localisation isomers) shown?
- Several engines on one spectrum: one row or one per engine? Which PEP? What does `consensus_support` count?
- Is `psm_id` unique, and is K? Does rank join the identity ([#229](https://github.com/bigbio/qpx/issues/229))?
- Which mzIdentML lists are converted? Which rows may carry `feature_id`?
- What is an identified spectrum, a PSM, a peptide?

## What happens today
**Spec.**
- `docs/spec/psm.md:3-5`, `producer-guide.md:11`: "one row per match".
- `psm.md:12-15` gives several hits of one peptidoform and charge producer ids or an extended composite, against `serialization.md:304-305` and `versioning.md:86-87` (one composite value, one row identity, provided ids included). `identity.md:40`: the mzIdentML footer does not reproduce its ids.
- Rank is only a recommended score, with no scope or tie rule (`psm.md:291`, `scores.md:58`). `consensus_support` = "Number of search engines supporting the identification" (`scores.md:75`, `qpx/core/builtin_scores.py:77`).
- The rank choice appears only as the CLI provenance option `rank1_only` (`provenance.md:65`).

**Issues.** [#229](https://github.com/bigbio/qpx/issues/229): one merged row per spectrum and identification; its per-engine variant and "rank joins the PK" note never reached the spec. [#262](https://github.com/bigbio/qpx/issues/262): lowest PEP among colliding hits. [#241](https://github.com/bigbio/qpx/issues/241): no filtering by default. [#378](https://github.com/bigbio/qpx/issues/378) (open) rule 5: no two rows share K.

| Producer | Behaviour | Findings |
|---|---|---|
| mzIdentML | Every SII of every SIL (`qpx/converters/mzidentml/document.py:44-62`, `psm_records.py:52-71`); MS:1002439 never read; rank only a score (`:140`), `--rank1-only` off; `psm_id` hashed from SII ids (`:221-235`) under the default composite. PXD054720 (Mascot, one SIL, 485 cross-link rows): 13,386 rows, 9,707 spectra, ranks 1-12, 3,679 rank>1 rows, 1,810 of them `pass_threshold=true` | [[MZID_CDAP-01](findings-psm.md#mzid_cdap-01)] [[MZID_CDAP-08](findings-psm.md#mzid_cdap-08)] [[MZID_CDAP-09](findings-cross.md#mzid_cdap-09)] [[EXTERNAL-M1](findings-psm.md#external-m1)] [[CRITIC-M2](findings-psm.md#critic-m2)] |
| Proteome Discoverer | All ranks unless `--rank1-only`; dedup on raw modification text, rank then q-value, no final key (`proteome_discoverer/psm_adapter.py:37-40`, `reports.py:88`); losing engine's scores dropped; missing `Rank` becomes 1 (`reports.py:91`). [#378](https://github.com/bigbio/qpx/issues/378): 20,001/109,700 TMT rows are rank 2-7 | [[PD-01](findings-psm.md#pd-01)] [[PD-02](findings-cross.md#pd-02)] [[PD-05](findings-cross.md#pd-05)] [[PD-07](findings-psm.md#pd-07)] [[PD-M3](findings-psm.md#pd-m3)] |
| OpenMS consensusXML | Every hit of every PID, rank never read (`openms_consensus/psm_adapter.py:260-273`); in-PID collisions keep the lowest PEP (`:296`), others' scores renamed `_2`, `_3` (`:86-98,227-233`); later-PID collisions dropped (`:375-381`). ConsensusID's fraction (n−1)/(N−1) (pyopenms: 2 of 3 engines → 0.5) fills the count-defined `consensus_support` (`:234-240`). Synthetic scan 24300: a rank-1 decoy and a rank-2 target both link to the ELEFANTK feature (`record_linking.py:21-33`) | [[OPENMS-01](findings-psm.md#openms-01)] [[OPENMS-19](findings-psm.md#openms-19)] |
| MaxQuant | Top hit plus MULTI-SECPEP, no rank; `Type` only a cv_param (`maxquant/psm_adapter.py:348-353`); `All sequences` dropped. Rows repeating K share a `psm_id`: 35 keys in maxquant_full evidence, 34 of them MULTI-SECPEP pairs, one MULTI-MSMS + MULTI-SECPEP (reverse, ids 167037/167038, scan 10817); 9 on PXD003133 ([#378](https://github.com/bigbio/qpx/issues/378)) | [[MAXQUANT-03](findings-psm.md#maxquant-03)] [[MAXQUANT-04](findings-psm.md#maxquant-04)] [[EXTERNAL-M2](findings-psm.md#external-m2)] |
| FragPipe | One row per psm.tsv row, no rank (`fragpipe/psm_adapter.py:144-199`); fragpipe.workflow not read. Synthetic scan 100: 3 rows at 2 charges, roles indistinguishable. Dotted run names give empty `scan` and duplicate ids | [[FRAGPIPE-03](findings-psm.md#fragpipe-03)] [[FRAGPIPE-04](findings-psm.md#fragpipe-04)] |
| CDAP | Top hits plus ties; `AmbiguousMatch` only a cv_param (`cdap/psm_adapter.py:26-27`); ties left out of features (`cdap/feature_adapter.py:116-120`) | [[MZID_CDAP-02](findings-psm.md#mzid_cdap-02)] [[MZID_CDAP-16](findings-feature.md#mzid_cdap-16)] [[CRITIC-M3](findings-psm.md#critic-m3)] |
| quantms MSstats | Writes no psm view | – |

Infrastructure:
- Writers only warn on a duplicate `psm_id` (`qpx/writers/base.py:557-569`).
- `duplicate_composite` already checks K over all files of a view, for every `identity_source`, reading `crosslinked_peptidoform` for cross-links; error under strict validation (`qpx/core/data/footer_validation.py:148-172`).
- Every provided-id file gets `identity_not_checked`, which fails `--conformance` (`:112-114`) [[VALIDATOR-02](findings-psm.md#validator-02)] [[VALIDATOR-04](findings-psm.md#validator-04)].

**Correction.** The 108 PXD054720 composite values carrying 239 ids [[INFRA-23](findings-psm.md#infra-23)] [[CRITIC-M1](findings-psm.md#critic-m1)] are cross-link species differing only in `crosslinked_peptidoform`; with the documented substitution there is no duplicate K (re-checked).

## Why it matters
- **Counts.** PXD054720: 13,386 rows for 9,707 identified spectra (+38%). `ModificationView` counts 4,323 Deamidated PSMs where 926 are passing rank-1 rows [[VALIDATOR-M3](findings-psm.md#validator-m3)]. Mascot passes 1,810 rank>1 rows, so "passing" is not "the spectrum's answer".
- **Hidden rank.** Once decoys are dropped, a rank-2 target that lost to a decoy looks like the identification (OpenMS scan 24300).
- **Joins.** Feature joins need a unique natural key ([#378](https://github.com/bigbio/qpx/issues/378)); MaxQuant and Proteome Discoverer break it.
- **Multi-engine results.** One search gives one row with scores dropped (PD), one merged row (OpenMS) or several rows (mzIdentML); `consensus_support` is a count or a fraction.

## Options considered
| Option | For | Against |
|---|---|---|
| Rank-1 only is mandatory | Simplest for readers | Breaks [#241](https://github.com/bigbio/qpx/issues/241); rescoring needs rank>1 |
| **All ranks allowed; typed `rank` (null = unknown)** | Lossless; one filter (`rank = 1`) for every producer | Readers must apply the filter |
| Rank stays a score | No schema change | About 100x the read cost (`scores.md:38`); absent means "rank 1" or "unknown" depending on producer |
| Rank per spectrum (literal mzIdentML SIR) | Same as mzIdentML | Sage, MaxQuant and DDA+ co-identifications become rank 2+, so a rank filter drops them [[EXTERNAL-06](findings-psm.md#external-06)] |
| **Rank per precursor ion** | A co-identification is rank 1 for its own ion ([T02](t02-chimeric-spectra.md) CHI-4); equals mzIdentML for single-precursor spectra | Charges ranked in separate source queries stay rank 1 each |
| **Merged row per hypothesis ([#229](https://github.com/bigbio/qpx/issues/229) default)** | Unique key; matches mzTab 1.0 §5.12.1 | The top-level PEP comes from one engine or the source's combination |
| Per-engine rows, `search_engine` in the composite ([#229](https://github.com/bigbio/qpx/issues/229) variant) | Engine-separate rows | New column; every count and join must deduplicate; conflicts with [#378](https://github.com/bigbio/qpx/issues/378) rule 5; never implemented |
| Rank joins the composite ([#229](https://github.com/bigbio/qpx/issues/229) note) | Allows the same hypothesis twice | Legalises duplicates; rescoring changes ids; unneeded once equal K merge |
| Ties collapsed to one row or `J` | Fewer rows | Loses the source's interpretation |
| **Ties as separate rows sharing rank 1** | mzIdentML rule; no new column | Readers must count spectra, not rows |

## Recommendation
**PSM-1 (MUST)**
- A psm row is one identification hypothesis for one spectrum. A spectrum with no hypothesis has no row.
- K = `(run_file_name, scan, peptidoform, charge)`, with `crosslinked_peptidoform` in place of `peptidoform` when non-null. K is unique within a psm view, across its files.
- Different peptidoforms or charges on one scan are separate rows.
- Where `scan` is empty, `spectrum_reference` stands in for it in K ([T11](t11-spectrum-precursor-fields.md) SPC-9); a row with neither is not written ([T11](t11-spectrum-precursor-fields.md) SPC-8). An empty `scan` from a parsable identifier (FragPipe dotted run names) is a converter defect ([T11](t11-spectrum-precursor-fields.md) SPC-13), never a reason for supplied ids.

**PSM-2 (MUST)** Hits with the same K are merged into one row: hits from several engines, SILs or PIDs, repeated source rows, and raw spellings of one canonical peptidoform.
- Each engine's scores stay in `additional_scores` under `<engine>_<score>` names (`scores.md` convention), never `_2` suffixes.
- `posterior_error_probability`, `qvalue` and `pass_threshold` take the source's combined values (ConsensusID, Percolator over merged results) when it reports them. Otherwise all three come from one merged hit: lowest PEP ([#262](https://github.com/bigbio/qpx/issues/262)), then lowest q-value, then engine name, nulls last (the [T06](t06-identifications-per-feature.md) order).
- `consensus_support` is the number of engines that reported this hypothesis, written only when the producer knows it (the engines it merged, or a count the source gives). A support measure of another kind keeps its own name: OpenMS ConsensusID's fraction becomes `consensusid_support`.
- No file or view holds one row per engine. This replaces the per-engine variant of [#229](https://github.com/bigbio/qpx/issues/229), which conflicts with [#378](https://github.com/bigbio/qpx/issues/378) rule 5 and was never implemented.

**PSM-3 (MUST)** `rank` (int16, ≥ 1, nullable) is the source's rank of the row's hypothesis among the hypotheses for one precursor ion: the isolation target, or one co-isolated precursor (MaxQuant second peptide, DDA+, CHIMERYS, Sage `chimera=true`; [T02](t02-chimeric-spectra.md) CHI-4).
- It is copied, never recomputed after rows are dropped: a target that lost to a decoy keeps rank 2. Ties share a rank (mzIdentML 1.2 `SpectrumIdentificationItem@rank`).
- Hits without a rank attribute (OpenMS PeptideIdentification) are ranked by score in the declared direction; equal scores share a rank.
- A source consensus rank across engines wins; otherwise a merged row takes the lowest engine rank.
- Where the source ranked the charge states of one precursor in separate queries (Comet, MSFragger `file.scan.scan.charge`), each query's ranks are kept, so several rank-1 rows of one ion are alternatives (PSM-6).
- A source rank of another scope (Sage co-identifications at rank 2+) goes to `<engine>_rank`. mzIdentML rank 0 (PMF) gives null.

**PSM-4 (MUST)** No footer key declares the rank policy; the `rank` column alone carries it.
- `rank` is non-null whenever the source states a rank. A source that reports only top hits gives 1; rank-1 ties appear as far as the source reports them.
- Otherwise `rank` is null, meaning unknown. A producer that does not know the rank never defaults it to 1.
- Removal of lower ranks is recorded only by the provenance option `rank1_only` (`provenance.md:65`).

**PSM-5 (SHOULD)** Converters drop no ranks by default ([#241](https://github.com/bigbio/qpx/issues/241)); `--rank1-only` is recorded as `rank1_only`. The writer fills `rank` on every entry point, the Python API included.

**PSM-6 (MUST)** Rank-1 alternatives stay one row each, with the source's sequence: ties (I/L, localisation isomers), charge states from separate queries, disagreeing engines, SILAC partner hypotheses. They are never collapsed to one row, to `J`, or to one isomer.
- An *alternative group* is two or more rank-1 rows sharing `(run_file_name, scan, precursor_index)` ([T02](t02-chimeric-spectra.md) CHI-3); rows with a null index form one group per `(run_file_name, scan)`. Rows of different indices are co-identifications, never alternatives.
- No tie column is added.

**PSM-7 (MUST)** A row with `rank > 1` has a null `feature_id` (also [T06](t06-identifications-per-feature.md) FID-7).

**PSM-8 (SHOULD)** pepmap and pg evidence draw on rank-1 rows only ([T10](t10-protein-mapping.md) PRT-14, [T12](t12-protein-groups.md) PGV-2).

**PSM-9 (MUST)** mzIdentML: when any SIL carries MS:1002439 (final PSM list), only flagged lists are converted; otherwise all lists are converted and merged per PSM-2. The provenance parameter `spectrum_identification_lists` records the SIL ids used.

**PSM-10 (MUST)** The psm `identity_composite` stays `run_file_name, scan, peptidoform, charge`, with the PSM-1 substitutions (`crosslinked_peptidoform`, `spectrum_reference`). Rank and engine never join it, which reverses the [#229](https://github.com/bigbio/qpx/issues/229) note: rank is an attribute of a hypothesis, and after PSM-2 two ranks never share K.
- Producer-supplied ids are allowed only when they map 1:1 to K.
- With K unique, `identity_not_checked` becomes informational ([T13](t13-metadata-conformance.md)) and no longer fails `--conformance` on its own.

## Schema and footer changes
| Change | View | Type / values | Required | Breaking? |
|---|---|---|---|---|
| `rank` column | psm | int16, nullable, ≥ 1; optional column; null = unknown | Non-null whenever the source states a rank (PSM-4) | No: additive (a 1.3 `qpxc validate` warns `undeclared_column`, `footer_validation.py:139-145`) |
| `rank` score | psm | Writer fills the typed column from it and, in 1.4, keeps the entry with the same value (the validator checks equality), exactly as [T04](t04-confidence-fields.md) CNF-6 does for q-values. Removing the entry is deferred to the next minor version and announced in the changelog (settles the rank part of [#377](https://github.com/bigbio/qpx/issues/377)). In 1.3 files the score is read as the rank. No rank on feature rows ([T04](t04-confidence-fields.md)) | – | No in 1.4: lookups of the score `rank` keep working. **Yes** in the next minor version, when the entry is removed |
| K unique for every `identity_source`; `psm.md:12-15` allowance removed; `identity.md:40` rewritten | psm | Spec text; the validator already checks K | Yes | **Yes** for producers: multi-SIL mzIdentML, MaxQuant K pairs, PD spelling variants |
| `consensus_support` = engine count (OpenMS fraction → `consensusid_support`); engine prefixes replace `_2` suffixes | psm | Score names | – | **Yes**: OpenMS score names and meaning |
| mzIdentML `psm_id` derived from K | psm | `identity_source=derived` | – | **Yes**: id values change; re-convert |

## Producer changes
| Producer | Change needed |
|---|---|
| mzIdentML | `rank` from SII@rank (0 → null); PSM-9 needs a pre-pass (in the 1.2 XSD a SIL's cvParams follow its results; `document.py:44-62` streams results first); merge across SILs; derived `psm_id`, SII ids in cv_params; all ranks kept unless `--rank1-only` |
| Proteome Discoverer | `rank` from `Rank`; no `Rank` gives null and a warning; partition on the canonical peptidoform; merge engine-prefixed scores and write `consensus_support`; pick per PSM-2; drop `rank` from feature scores (`feature_adapter.py:126`) |
| OpenMS consensusXML | `rank` per PSM-3; merge cross-PID collisions (`psm_adapter.py:375-381`); engine-prefixed names; ConsensusID fraction → `consensusid_support`; link only same-peptidoform, same-charge, rank-1 target rows (`record_linking.py:30-33`) |
| OpenMS `-out_qpx` / `qpxc convert openms` | C++ writer, verified at HEAD: one row per PeptideIdentification top hit (ProteomicsLFQ keeps only the best hit, and one row per spectrum, peptidoform and charge), so K is unique; no `rank` column; score names are raw OpenMS keys. Needed: `rank = 1` on every row (PSM-4) and `scores.md` score names. The Python upgrade (`qpx/converters/openms/converter.py`) keeps the source's `rank` or leaves it null. See the [OpenMS worklist](openms-worklist.md). |
| MaxQuant | `rank = 1`; merge every pair of rows repeating K (lowest PEP), SECPEP pairs included; [T02](t02-chimeric-spectra.md) marks SECPEP rows |
| FragPipe | psm.tsv has no rank. Read fragpipe.workflow (new input): `msfragger.output_report_topN = 1` gives `rank = 1`, with DDA+ extra rows marked by [T02](t02-chimeric-spectra.md); without the file, null |
| CDAP | `rank = 1`; AmbiguousMatch rows form an alternative group |
| quantms MSstats | None (no psm view) |
| Writer and validator | Fill `rank`; check `rank > 1 ⇒ feature_id IS NULL`; warn when `rank` is null and a spectrum has several rows |

## Reader rules
**PSM-R1 (MUST)** An identification row has `rank = 1` (1.3 files: score `rank` = 1). Rows with `rank > 1` are never counted as identifications or used as peptide, protein or modification evidence, whatever their `pass_threshold` or `qvalue`; their q-values and PEPs are not FDR-controlled ([T04](t04-confidence-fields.md) CNF-10).

**PSM-R2 (MUST)** Over accepted rows ([T03](t03-filtering-state.md) FLT-R4: target, rank 1, passing the declared filter):
- identified spectra = `COUNT(DISTINCT (run_file_name, scan))`, with `spectrum_reference` in place of an empty `scan` ([T11](t11-spectrum-precursor-fields.md) SPC-9);
- PSMs (identified precursors) = `COUNT(DISTINCT (run_file_name, scan, precursor_index))`, with the same substitution ([T02](t02-chimeric-spectra.md) CHI-3, CHI-R2): an alternative group counts once, each co-isolated precursor counts. Rows with a null index count once per spectrum, and the reader reports how many spectra that covers ([T02](t02-chimeric-spectra.md) CHI-R7). Such a count is a lower bound on identified precursors;
- peptidoforms = `COUNT(DISTINCT peptidoform)`; each tied peptidoform counts, as the source ranked it 1.

**PSM-R3 (MUST)** With a null or absent `rank` a reader does not assume top-1: counts are upper bounds, with a warning when a spectrum has several rows.

**PSM-R4 (MUST)** A reader joins on K only in a view that passed strict validation; otherwise it deduplicates K first.

**Current violations:**
- `qpx/core/data/psm.py:23-25`: `confident()` ignores rank.
- `qpx/converters/summary.py:109,121`: `n_psms = COUNT(*)`; without a feature view, `n_peptidoforms` counts peptidoforms of all psm rows (any rank, decoys).
- `qpx/views/api.py:139-151`: `ModificationView.frequency()` counts every non-decoy row.
- `qpx/core/dataset_queries.py:308-313`: `link_feature_psm()` returns links unchecked, so the OpenMS rank-2 link reads as feature evidence.
- mokume reads no rank, `pass_threshold` or q-value (`mokume:rust/crates/mokume-io/src/psm.rs:14-25`) and counts linked `psm_id` rows, not spectra (`mokume:rust/crates/mokume-pipeline/src/spectral_count.rs:158-195,235-249`) [[MOKUME-25](findings-psm.md#mokume-25)]. Its hard error on a duplicate `psm_id` (`:222-230`) fits PSM-1.
- The MuData export (`qpx/mudata.py`) does not read psm.

## Open questions for maintainers
1. Should `qpxc validate --profile dda` require `rank = 1` on every row, or accept all ranks ([#241](https://github.com/bigbio/qpx/issues/241)) and rely on PSM-R1 ([#378](https://github.com/bigbio/qpx/issues/378) "where do you stop")?
2. Should a typed `search_engines list<string>` replace parsing prefixes and PD's `identifying_node` [[CRITIC-12](findings-psm.md#critic-12)]?
3. FragPipe: is psm.tsv always top-1 per query, DDA+ aside? Then fragpipe.workflow is not needed.
4. PD `PSM Ambiguity` (Selected/Rejected) resolves ties: cv_param, or mark Rejected rows?
5. MaxQuant K pairs point to two evidence rows, but the merged PSM can link only one (`feature_id` is single-valued). Accept, or merge the evidence rows too ([T06](t06-identifications-per-feature.md))?
6. Should a duplicate K be fatal at write time, ending the tolerance of `versioning.md:120-127`?
7. Confirm reversing the [#229](https://github.com/bigbio/qpx/issues/229) "rank joins the PK" note, dropping its per-engine variant (PSM-2, PSM-10), and the OpenMS `consensus_support` rename.
