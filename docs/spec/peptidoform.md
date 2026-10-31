# Peptidoform

A **peptidoform** is a specific molecular form of a peptide: the amino acid sequence together with all of its modifications. In QPX, peptidoforms are encoded using the [ProForma notation](https://github.com/HUPO-PSI/ProForma), which provides a standardized, human-readable string representation for modified peptide sequences.

The `peptidoform` field appears alongside the bare `sequence` field (amino acids only) in every identification- and quantification-level view. While `sequence` gives the unmodified backbone, `peptidoform` captures the complete molecular identity needed to distinguish different modified forms of the same peptide.

## ProForma notation

ProForma is a HUPO-PSI standard for encoding peptide sequences with modifications in a single string. Modifications are placed in square brackets immediately after the residue they modify (or before the sequence for N-terminal modifications).

### Canonical form

Every QPX producer writes the peptidoform in one canonical form, so the same
molecule gets the same string, and therefore the same identity, whichever
tool produced it:

1. **Every modification is written inline**, including fixed modifications
   (carbamidomethyl cysteine), isobaric labels (TMT, iTRAQ) and the MS1 label
   (SILAC, dimethyl, mTRAQ) of an identified PSM, even where the search engine
   leaves them out of its own sequence column. A feature that holds a whole
   MS1-label multiplet carries its lightest partner's peptidoform, since the
   partners differ only in the label
   ([MS1 labels](intensities.md#ms1-labels-silac-dimethyl-mtraq)).
2. **A modification with a Unimod entry is written as `UNIMOD:n`** in upper
   case. Names (`Oxidation`, `GlyGly (K)`) and PSI-MOD accessions are mapped to
   their Unimod accession.
3. **A modification without a Unimod entry is written as its signed
   monoisotopic mass shift in Daltons with four decimals**, `[+79.9663]` or
   `[-17.0265]`. `CHEMMOD:` is not used. An accession QPX cannot check the mass
   of keeps the reported mass after a pipe, `[UNIMOD:1842|+54.0106]`.
4. N-terminal modifications precede the sequence and a hyphen,
   `[UNIMOD:1]-PEPTIDE`; C-terminal ones follow a hyphen, `PEPTIDE-[UNIMOD:2]`.
5. A modification whose position is unknown, such as the mass shift of an open
   search that is not localized to one residue, comes first and ends with a
   question mark: `[+129.0804]?[UNIMOD:1]-PEPTIDE`.
6. A glycan follows rule 2: its Unimod accession when one matches its mass and
   residue (`N[UNIMOD:1531]`), otherwise its mass (`N[+2715.9624]`), so every
   reader can parse it. Its composition, when the producer reports one, is the
   name of its entry in `modifications` (`Glycan:HexNAc2Hex9`).
7. `calculated_mz` is the theoretical m/z of the peptidoform at the reported
   charge whenever every modification in it has a known mass.

A name that resolves to no Unimod entry and comes with no mass is kept as
written (`K[Custom (K)]`); parentheses inside a bracketed tag belong to the
tag. Readers that accept other ProForma spellings can bring them to the
canonical form with `qpx.converters.ptm.normalize_peptidoform`.

### Examples

**Variable modification**

```text
PEPTIDM[UNIMOD:35]K
```

Methionine at position 7 carries an oxidation (Unimod 35).

**N-terminal modification**

```text
[UNIMOD:1]-PEPTIDMK
```

Acetylation on the N-terminus of the peptide.

**Fixed and variable modifications together**

```text
[UNIMOD:1]-PEPT[UNIMOD:21]IDC[UNIMOD:4]M[UNIMOD:35]K
```

N-terminal acetylation, a phosphorylated threonine at position 4, the fixed
carbamidomethylation of the cysteine and an oxidised methionine.

**TMT-labeled peptide**

```text
[UNIMOD:737]-PEPTIDM[UNIMOD:35]K[UNIMOD:737]
```

TMT6plex (Unimod 737) on the N-terminus and on the C-terminal lysine, with an
additional oxidation on methionine. TMTpro is Unimod 2016, iTRAQ4plex 214 and
iTRAQ8plex 730.

**Modification without a Unimod entry**

```text
PEPTIDK[+54.0106]R
```

**Mass shift of unknown position**

```text
[+129.0804]?PEPTIDEK
```

## Where peptidoform is used

The `peptidoform` field is present in the following QPX views:

| View                  | Field name     | Description                                                    |
| --------------------- | -------------- | -------------------------------------------------------------- |
| PSM (`psm_file`)      | `peptidoform`  | The identified peptidoform for this spectrum match              |
| Feature (`feature_file`) | `peptidoform` | The quantified peptidoform for this feature                   |
| Peptide (`peptide_file`) | `peptidoform` | The peptidoform summarized across samples                     |

In all three views, the `peptidoform` string serves as a key component of the record identity alongside `charge` and `run_file_name`.

## Relationship to the modifications field

The `peptidoform` string encodes modifications inline, which is compact and human-readable. However, some use cases require structured access to modification details -- for example, querying all phosphorylation sites with their localization probabilities.

For this purpose, QPX provides a separate `modifications` field that represents the same information as a structured array of records. See the [Modifications](modifications.md) page for the full struct definition and examples.

!!! warning
    The `peptidoform` string and the `modifications` struct should be consistent with each other. If both are present, consumers should treat the `peptidoform` as the canonical representation and the `modifications` field as supplementary detail.

## Further reading

- [ProForma specification (HUPO-PSI)](https://github.com/HUPO-PSI/ProForma)
- [Modifications](modifications.md) -- structured modification representation with localization scores
- [QPX Format Overview](index.md) -- full list of views and concepts
