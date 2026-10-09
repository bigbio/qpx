# Mass Spectra View

The mass spectra (mz) view is a Parquet file that stores raw mass spectrometry spectral data in a columnar format. This view is based on the [mz_parquet](https://github.com/lazear/mz_parquet) format developed by Michael Lazear, adapted for the QPX ecosystem with additional metadata support through controlled vocabulary (CV) parameters.

## Use cases

- Retrieve precursor mass, retention time, and intensity values for all spectra in a file.
- Enable visualization and interactive scanning at the mass spectra level.
- Provide structured spectral data for AI/ML training and prediction tasks.
- Store spectra in a columnar format that is significantly more efficient than mzML for large-scale retrieval.

## Schema

See the full YAML schema in [`mz.yaml`](schemas/mz.yaml).

### Scan identification

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `id` | Unique identifier for the scan or spectrum | `string` | Yes |
| `run_file_name` | Run file containing this spectrum; links to PSM/feature `run_file_name` | `string` | No |
| `scan` | Scan identifier parsed from the native ID as an array of integer components (see [Scan Numbers](scan.md)); links to PSM/feature `scan` | `array[int32]` | No |
| `ms_level` | The MS level (1 for MS1, 2 for MS2, etc.) | `int` | Yes |
| `centroid` | Whether the data is centroided (`true`) or profile mode (`false`) | `boolean` | Yes |

### Timing and mobility

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `scan_start_time` | Start time of the scan in minutes | `float` | Yes |
| `inverse_ion_mobility` | Inverse ion mobility, used for TIMS data (1/K0) | `float` | No |
| `ion_injection_time` | Ion injection time in milliseconds | `float` | Yes |
| `total_ion_current` | Total ion current (TIC) for the scan | `float` | Yes |

### Precursors

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `precursors` | List of precursor ions for this MS/MS scan | `array[struct]` | No |

Each precursor struct contains:

| Sub-field | Description | Type | Required |
|-----------|-------------|------|----------|
| `selected_ion_mz` | m/z value of the selected precursor ion | `float` | Yes |
| `selected_ion_charge` | Charge state of the selected precursor ion | `int` | No |
| `selected_ion_intensity` | Intensity of the selected precursor ion | `float` | No |
| `isolation_window_target` | Target m/z for the isolation window | `float` | No |
| `isolation_window_lower` | Lower bound of the isolation window | `float` | No |
| `isolation_window_upper` | Upper bound of the isolation window | `float` | No |
| `spectrum_ref` | Reference to another spectrum (for linking to external datasets) | `string` | No |

### Spectral data

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `mz` | Array of m/z values for the scan | `array[float]` | Yes |
| `intensity` | Array of intensity values corresponding to the m/z values | `array[float]` | Yes |

### CV params

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `cv_params` | Optional list of controlled vocabulary parameters for additional metadata | `array[struct{name, value}]` | No |

Each CV param struct contains:

| Sub-field | Description | Type |
|-----------|-------------|------|
| `name` | Name of the CV term (e.g., from PSI-MS or other ontologies) | `string` |
| `value` | Value associated with the CV term | `string` |

## Example

An MS2 spectrum with a single precursor ion:

```json
{
  "id": "controllerType=0 controllerNumber=1 scan=15236",
  "ms_level": 2,
  "centroid": true,
  "scan_start_time": 42.156,
  "inverse_ion_mobility": null,
  "ion_injection_time": 35.0,
  "total_ion_current": 2.45e6,
  "precursors": [
    {
      "selected_ion_mz": 547.2891,
      "selected_ion_charge": 2,
      "selected_ion_intensity": 8.3e5,
      "isolation_window_target": 547.29,
      "isolation_window_lower": 0.7,
      "isolation_window_upper": 0.7,
      "spectrum_ref": null
    }
  ],
  "mz": [
    110.0713, 120.0808, 136.0757, 147.1128, 175.1190,
    234.1448, 262.1397, 349.1718, 462.2559, 575.3399,
    688.4240, 801.5080, 914.5921, 1013.6605
  ],
  "intensity": [
    5200.0, 12300.0, 8900.0, 45000.0, 23000.0,
    67000.0, 34000.0, 89000.0, 120000.0, 95000.0,
    78000.0, 56000.0, 34000.0, 12000.0
  ],
  "cv_params": [
    {"name": "MS:1000016", "value": "42.156"},
    {"name": "MS:1000512", "value": "FTMS + c NSI d Full ms2 547.29@hcd30.00"}
  ]
}
```

## Notes

!!! note "Preferred over embedding spectra in PSM view"
    While the PSM view can optionally carry `mz_array` and `intensity_array` fields, the mz view is the preferred format for storing and retrieving spectral data. Separating spectra from identifications allows independent access patterns: spectra can be queried by retention time or m/z range without loading identification data, and vice versa.

!!! tip "Querying with DuckDB"
    The columnar Parquet layout enables efficient SQL-based access to spectra:

    ```sql
    -- Retrieve all MS2 spectra in a retention time window
    SELECT id, scan_start_time, precursors, mz, intensity
    FROM 'sample.mz.parquet'
    WHERE ms_level = 2 AND scan_start_time BETWEEN 30.0 AND 45.0

    -- Count spectra per MS level
    SELECT ms_level, COUNT(*) as num_spectra
    FROM 'sample.mz.parquet'
    GROUP BY ms_level
    ```

!!! warning "Array alignment"
    The `mz` and `intensity` arrays MUST have the same length within each record. Each element at index `i` in `intensity` corresponds to the peak at index `i` in `mz`.
