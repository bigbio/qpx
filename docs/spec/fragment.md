# Fragment View

The fragment view holds what a DIA tool measured of each fragment of a precursor: its quantity, peak shape and, where a tool has it, the extracted ion chromatogram. It is optional, can be large, and is written only on request (`--fragments`). Specification 1.3 added it.

File: `{PREFIX}.fragment.parquet`. One row per fragment of a feature, channel and run. Identity composite: `feature_id, label, ion_type, fragment_charge, neutral_loss, fragment_index`.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `fragment_id` | Opaque identity, derived from the composite | int64 | Yes |
| `feature_id` | The feature (`feature.feature_id`) the fragment belongs to; `qpxc validate` checks it | int64 | Yes |
| `run_file_name`, `peptidoform`, `charge` | The precursor and run, repeated for reading without a join | string, string, int16 | Yes |
| `label` | Channel | string, null | No |
| `fragment_index` | 1-based position in the tool's list, where it names fragments by position | int32, null | No |
| `ion_type` | Ion series and number, such as `y8`; null where the tool does not name the fragment | string, null | No |
| `fragment_charge` | Charge of the fragment | int32, null | No |
| `neutral_loss` | Neutral loss such as `H2O`; null for none | string, null | No |
| `mz` / `observed_mz` | Theoretical and measured m/z | float64, null | No |
| `intensity` | Raw quantity; null for not measured | float32, null | No |
| `corrected_intensity` | Quantity after the tool's correction | float32, null | No |
| `peak_height` | Peak height | float32, null | No |
| `correlation` | Correlation of the fragment's elution with the precursor's other fragments | float32, null | No |
| `excluded` | Whether the tool left the fragment out of the precursor's quantity | bool, null | No |
| `rt_array`, `intensity_array` | Extracted ion chromatogram, seconds and intensities. No converter fills them yet: they are for a producer that has them, as OpenMS' XIC files | list\<float32\>, null | No |
| `additional_scores` | Other per-fragment values | list\<score\>, null | No |

`qpxc transform relayout` can store `additional_scores` as a struct with one field per name ([Struct lists](serialization.md#struct-lists)).

## Producers

| | DIA-NN `--fragments` | Spectronaut `--fragments` |
|---|---|---|
| `fragment_index` | Position in `Fragment.Quant.Raw` | --- |
| `ion_type`, `fragment_charge`, `neutral_loss`, `mz` | `Fragment.Info`, where the report has it; else null, and the position is all that names the fragment | `F.FrgIon`, `F.Charge`, `F.FrgLossType`, `F.TheoreticalMz` |
| `intensity` | `Fragment.Quant.Raw` | `F.PeakArea` |
| `corrected_intensity` | `Fragment.Quant.Corrected` | `F.NormalizedPeakArea` |
| `correlation` | `Fragment.Correlations` | --- |
| `observed_mz`, `peak_height`, `excluded` | --- | `F.MeasuredMz`, `F.PeakHeight`, `F.ExcludedFromQuantification` |
| `additional_scores` | --- | `F.InterferenceScore`, `F.Log10SignalToNoise`, `F.MeasuredRelativeIntensity`, `F.PredictedRelativeIntensity`, `F.CalibratedMassAccuracy (PPM)` |

DIA-NN lists every library fragment, with 0 for one it did not measure. A fragment with neither quantity has no row, so `fragment_index` can skip values.

## Targeted runs

A transition of a PRM or SRM run is a fragment row. See [Extensions](extensions.md#targeted-runs-prm-srm-and-heavy-standards).
