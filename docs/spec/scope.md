# Scope

QPX describes **bottom-up proteomics**: peptides identified from MS/MS or DIA spectra, their quantities per run, and protein groups inferred from them. This page says what the format models, what it deliberately leaves out, and how much of ProForma it reads.

## What QPX models

| Area | Covered |
|------|---------|
| Acquisition | DDA, DIA, PRM and other targeted methods (see [Extensions](extensions.md) for targeted transitions and curve fits) |
| Quantification | Label-free, MS1 labels (SILAC, dimethyl, mTRAQ), isobaric labels (TMT, iTRAQ) |
| Identification | One row per peptide-spectrum match, including cross-links ([PSM](psm.md)) |
| Quantities | Precursors ([Feature](feature.md)), protein groups ([PG](pg.md)), sites ([Site](site.md)), fragments ([Fragment](fragment.md)), proteoforms ([Proteoform](proteoform.md)) |
| Metadata | [Sample](sample.md) and [run](run.md) tables from an SDRF, [provenance](provenance.md) |

## Not modelled

These data may pass validation, but no column carries what makes them different, so a reader cannot trust them:

- **Top-down proteomics.** Intact proteoforms with charges such as 22 fit `charge`, and `charge = 0` validates, yet no field describes intact-mass deconvolution, and sequence-level fields assume a digest.
- **Spectral libraries.** A library entry with an empty `run_file_name` validates, but a row in `feature` or `psm` states a measurement in a run. Keep libraries in their own formats.
- **Metabolomics, lipidomics and glycomics without a peptide.** Every quantity joins to a peptide or protein accession.
- **Raw signal.** The [mz](mz.md) view holds centroided peaks and acquisition parameters. Profile data, ion-mobility frames and XIC traces stay in the instrument or mzML files.
- **Statistical conclusions from the format itself.** QPX stores what a tool reported. Differential results are stored in the [differential](differential.md) file with the method named; QPX does not make them comparable across tools.

## ProForma subset

QPX reads and writes this subset of [ProForma 2.0](https://github.com/HUPO-PSI/ProForma), which is the *Base* compliance level plus unknown-position tags:

| Construct | Example | Support |
|-----------|---------|---------|
| Residue tag by accession | `M[UNIMOD:35]` | Written and parsed |
| Residue tag by name | `M[Oxidation]` | Parsed; canonical form uses the accession |
| Mass shift | `M[+15.9949]` | Written when no Unimod entry matches; parsed |
| Terminal tags | `[UNIMOD:1]-PEPTIDE`, `PEPTIDE-[UNIMOD:2]` | Written and parsed |
| Unknown position | `[Phospho]?PEPTIDE` | Written and parsed; lists no positions |

Everything else is **passed through unchanged**: the string survives in `peptidoform`, and `from_proforma` reports no modifications for it rather than guess them.

| Construct | Example |
|-----------|---------|
| Cross-link and chimeric species | `EMEVTK[XLMOD:02001#XL1]SESPEK//EMEVTK[#XL1]SESPEK`, `PEPTIDE/2+SEQUENCE/3` |
| Charge suffix | `PEPTIDE/2` |
| Global modification | `<[UNIMOD:4]@C>ACDEK`, `<13C>PEPTIDEK` |
| Labile modification | `{Glycan:Hex}PEPTIDE` |
| Localisation group | `PEPT[Phospho#g1(0.9)]IDES[#g1(0.1)]K` |
| Range | `PRT(ES)[+79.97]EM` |

Cross-linked species have their own column, `crosslinked_peptidoform`, documented in [PSM](psm.md). Localisation probabilities go in `modifications[].positions[].scores`, not in the string.
