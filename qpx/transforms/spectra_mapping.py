"""Spectra mapping transform — mzML spectral data extraction for PSM annotation.

This transform reads mzML files and extracts spectral arrays (m/z and intensity)
for identified scans, enabling PSM-level spectral annotation. It uses pyopenms
for mzML parsing with support for multiple native ID formats from different
instrument vendors.

The transform can:
1. Extract spectra for individual scans by scan number
2. Batch-extract spectra for all PSMs in a Dataset
3. Write annotated PSMs with spectral data to MZ Parquet files

Usage:
    mapper = SpectraMappingTransform(mzml_directory="/data/mzml/")

    # Extract a single spectrum
    spectrum = mapper.get_spectrum("run1", scan_number=12345)

    # Annotate PSMs in a Dataset with spectral arrays
    annotated_df = mapper.annotate_dataset_psms(dataset)

    # Write spectral data to .mz.parquet
    mapper.write_mz_parquet(dataset, "output.mz.parquet")
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from qpx.core.scan import scan_format_from_native_id, scan_from_native_id

logger = logging.getLogger(__name__)


# --- Native ID patterns for different instrument vendors ---
NATIVE_ID_PATTERNS = [
    "scan=(?<SCAN>\\d+)",  # Standard Thermo format
    "cycle=(?<SCAN>\\d+)",  # Waters/Agilent format
    "index=(?<SCAN>\\d+)",  # Generic index format
    "spectrum=(?<SCAN>\\d+)",  # Alternative format
]

# PSI-MS "ion injection time"; pyopenms keeps it on the scan's acquisition.
_ISOLATION_TARGET = "isolation window target m/z"
_FILTER_STRING = "filter string"
_ION_INJECTION_TIME = "MS:1000927"


# Peaks per written batch: bounds memory whatever the file size, while keeping row
# groups large enough to compress well.
_MZ_BATCH_PEAKS = 10_000_000


class _MzRecordWriter:
    """pyopenms mzML consumer that writes QPX mz records in peak-bounded batches."""

    def __init__(self, writer, build_record):
        self.acquisition_start: str | None = None
        self._writer = writer
        self._build_record = build_record
        self._records: list[dict] = []
        self._peaks = 0
        self._index = 0
        self.written = 0

    def consume_spectrum(self, spectrum) -> None:
        """Buffer one spectrum's record; write the batch once it holds enough peaks."""
        index = self._index
        self._index += 1
        record = self._build_record(spectrum, index)
        if record is None:
            return
        self._records.append(record)
        self._peaks += len(record["mz"])
        if self._peaks >= _MZ_BATCH_PEAKS:
            self.flush()

    def flush(self) -> None:
        """Write the buffered records."""
        if self._records:
            self._writer.write_table(pa.Table.from_pylist(self._records, schema=self._writer.arrow_schema))
            self.written += len(self._records)
            self._records, self._peaks = [], 0

    def experimental_settings(self, settings) -> None:
        """Keep the mzML ``startTimeStamp``."""
        self.acquisition_start = settings.getDateTime().get() or None

    def ignore(self, *_args) -> None:
        """Expected sizes and chromatograms are not written."""

    # pyopenms MzMLFile.transform calls camelCase consumer hooks; route them here.
    _HOOKS = {
        "consumeSpectrum": "consume_spectrum",
        "consumeChromatogram": "ignore",
        "setExperimentalSettings": "experimental_settings",
        "setExpectedSize": "ignore",
    }

    def __getattr__(self, name):
        if name in self._HOOKS:
            return getattr(self, self._HOOKS[name])
        raise AttributeError(name)


def _total_ion_current(spectrum, intensity) -> float:
    """The scan's reported TIC (MS:1000285); the summed peak intensity when the mzML has none."""
    if spectrum.metaValueExists("total ion current"):
        return float(spectrum.getMetaValue("total ion current"))
    return float(intensity.sum()) if len(intensity) else 0.0


def _drift_value(spectrum, unit) -> float | None:
    """The spectrum's drift value when it is in *unit*; OpenMS holds 1/K0, drift time and FAIMS CV in one slot."""
    return float(spectrum.getDriftTime()) if spectrum.getDriftTimeUnit() == unit else None


def _peak_ion_mobility(spectrum, oms) -> dict:
    """The per-peak ion mobility array and its unit, for spectra that carry one."""
    if not spectrum.containsIMData():
        return {"ion_mobility_array": None, "ion_mobility_array_unit": None}
    index, unit = spectrum.getIMData()
    units = {oms.DriftTimeUnit.VSSC: "Vs/cm^2", oms.DriftTimeUnit.MILLISECOND: "ms"}
    # get_data() is a view into the array object, so the object has to outlive the copy.
    peak_array = spectrum.getFloatDataArrays()[index]
    return {"ion_mobility_array": peak_array.get_data().copy(), "ion_mobility_array_unit": units.get(unit)}


def _ion_injection_time(spectrum) -> Optional[float]:
    """Ion injection time (ms) of the scan, None when the mzML does not report it."""
    for acquisition in spectrum.getAcquisitionInfo():
        if acquisition.metaValueExists(_ION_INJECTION_TIME):
            return float(acquisition.getMetaValue(_ION_INJECTION_TIME))
    return None


def _isolation_target(precursor) -> float:
    """The isolation window target, which differs from the selected ion in a DIA or wide-window scan."""
    if precursor.metaValueExists(_ISOLATION_TARGET):
        return float(precursor.getMetaValue(_ISOLATION_TARGET))
    return float(precursor.getMZ())


def _spectrum_cv_params(spectrum) -> list[dict]:
    """Activation, collision energy, filter string and analyzer of a spectrum, as far as the mzML states them."""
    params: list[dict] = []
    if spectrum.metaValueExists(_FILTER_STRING):
        filter_string = str(spectrum.getMetaValue(_FILTER_STRING))
        params.append({"cv_name": "filter_string", "cv_value": filter_string})
        analyzer = filter_string.split(" ", 1)[0]
        if analyzer.endswith("MS"):
            params.append({"cv_name": "analyzer", "cv_value": analyzer})
    for precursor in spectrum.getPrecursors()[:1]:
        methods = [
            method.decode() if isinstance(method, bytes) else method for method in precursor.getActivationMethodsAsString()
        ]
        if methods:
            params.append({"cv_name": "activation_method", "cv_value": ";".join(methods)})
        if precursor.getActivationEnergy() > 0:
            params.append({"cv_name": "collision_energy", "cv_value": str(precursor.getActivationEnergy())})
    return params


class _MzMLCache:
    """
    Cache for loaded mzML experiments and spectrum lookup objects.

    Avoids reloading the same mzML file for multiple scan lookups.
    This is an internal implementation detail of SpectraMappingTransform.
    """

    def __init__(self):
        self._experiments: dict[str, object] = {}
        self._lookups: dict[str, object] = {}
        self._patterns: dict[str, str] = {}
        self._component_indices: dict[str, dict[tuple[int, ...], Optional[int]]] = {}

    def _component_index(self, mzml_path: str, components: tuple[int, ...]) -> Optional[int]:
        """Find a unique complete native ID; never fall back to one component."""
        if mzml_path not in self._component_indices:
            indices: dict[tuple[int, ...], Optional[int]] = {}
            for index, spectrum in enumerate(self._experiments[mzml_path]):
                key = tuple(scan_from_native_id(spectrum.getNativeID()))
                if key:
                    indices[key] = None if key in indices else index
            self._component_indices[mzml_path] = indices
        indices = self._component_indices[mzml_path]
        index = indices.get(components)
        if index is None:
            logger.warning("No unique spectrum for native ID components %s in %s", components, mzml_path)
        return index

    def get_spectrum(
        self,
        mzml_path: str,
        scan_number: Union[int, tuple[int, ...]],
        native_id_pattern: Optional[str] = None,
    ) -> tuple[int, np.ndarray, np.ndarray]:
        """
        Get a spectrum from a cached mzML file.

        Args:
            mzml_path: Full path to the mzML file.
            scan_number: Scan number or complete native ID components to retrieve.
            native_id_pattern: Optional custom regex for native ID parsing.

        Returns:
            Tuple of (num_peaks, mz_array, intensity_array).
            Returns (0, empty_array, empty_array) if scan not found.
        """
        try:
            import pyopenms as oms
            from pyopenms import SpectrumLookup
        except ImportError:
            raise ImportError("pyopenms is required for mzML parsing. Install it with: pip install pyopenms")

        # Load experiment if not cached
        if mzml_path not in self._experiments:
            logger.debug(f"Loading mzML file: {mzml_path}")
            exp = oms.MSExperiment()
            oms.MzMLFile().load(mzml_path, exp)
            self._experiments[mzml_path] = exp
            # Reset lookup for new file
            self._lookups.pop(mzml_path, None)
            self._patterns.pop(mzml_path, None)

        exp = self._experiments[mzml_path]

        if isinstance(scan_number, tuple):
            index = self._component_index(mzml_path, scan_number)
            if index is None:
                return 0, np.array([], dtype=np.float32), np.array([], dtype=np.float32)
            return self._spectrum_arrays(exp.getSpectrum(index))

        # Build spectrum lookup if not cached
        if mzml_path not in self._lookups:
            patterns = [native_id_pattern] if native_id_pattern else NATIVE_ID_PATTERNS

            for pattern in patterns:
                try:
                    lookup = SpectrumLookup()
                    lookup.readSpectra(exp, pattern)
                    # Verify the pattern works by looking up the requested scan
                    lookup.findByScanNumber(scan_number)
                    self._lookups[mzml_path] = lookup
                    self._patterns[mzml_path] = pattern
                    logger.debug(f"Using native ID pattern: {pattern}")
                    break
                except Exception:
                    logger.debug("Native ID pattern %s failed for %s", pattern, mzml_path, exc_info=True)
                    continue

            if mzml_path not in self._lookups:
                logger.warning(f"Could not find working native ID pattern for {mzml_path}")
                return 0, np.array([], dtype=np.float32), np.array([], dtype=np.float32)

        lookup = self._lookups[mzml_path]

        try:
            index = lookup.findByScanNumber(scan_number)
        except (IndexError, RuntimeError):
            return 0, np.array([], dtype=np.float32), np.array([], dtype=np.float32)

        return self._spectrum_arrays(exp.getSpectrum(index))

    def get_spectrum_at(self, mzml_path: str, position: int) -> tuple[int, np.ndarray, np.ndarray]:
        """The spectrum at a 0-based *position* in an mzML file (``scan_format`` index)."""
        import pyopenms as oms

        if mzml_path not in self._experiments:
            exp = oms.MSExperiment()
            oms.MzMLFile().load(mzml_path, exp)
            self._experiments[mzml_path] = exp
        exp = self._experiments[mzml_path]
        if isinstance(position, int) and 0 <= position < exp.getNrSpectra():
            return self._spectrum_arrays(exp.getSpectrum(position))
        return 0, np.array([], dtype=np.float32), np.array([], dtype=np.float32)

    @staticmethod
    def _spectrum_arrays(spectrum) -> tuple[int, np.ndarray, np.ndarray]:
        """Return peak arrays with the transform's stable output dtypes."""
        mz_array, intensity_array = spectrum.get_peaks()
        return (
            len(mz_array),
            np.array(mz_array, dtype=np.float32),
            np.array(intensity_array, dtype=np.float32),
        )

    def clear(self):
        """Clear all cached experiments and lookups."""
        self._experiments.clear()
        self._lookups.clear()
        self._patterns.clear()
        self._component_indices.clear()


def _psm_scan_format(psm) -> str | None:
    """The scan_format a single PSM file declares; several files may declare different ones."""
    paths = psm.file_paths
    if len(paths) != 1 or not Path(paths[0]).is_file():
        return None
    declared = (pq.read_metadata(paths[0]).metadata or {}).get(b"scan_format")
    return declared.decode() if declared else None


class SpectraMappingTransform:
    """
    Extract spectral arrays from mzML files and annotate PSM data.

    This transform maps scan numbers from PSM identifications to their
    corresponding m/z and intensity arrays in mzML files. It uses pyopenms
    for efficient mzML parsing with caching to avoid reloading files.

    The transform supports multiple instrument native ID formats (Thermo,
    Waters, Agilent, etc.) and automatically detects the correct format.

    Usage:
        mapper = SpectraMappingTransform(mzml_directory="/data/mzml/")

        # Get a single spectrum
        n_peaks, mz, intensity = mapper.get_spectrum("run1", 12345)

        # Annotate all PSMs in a dataset
        psm_df = mapper.annotate_dataset_psms(dataset)

        # Write MZ spectral data to Parquet
        mapper.write_mz_parquet(dataset, "output.mz.parquet")
    """

    # Common mzML file extensions to try
    MZML_EXTENSIONS = [".mzML", ".mzml", ".MZML", ".mzML.gz", ".mzml.gz", ".MZML.gz"]

    def __init__(
        self,
        mzml_directory: Union[str, Path],
        native_id_pattern: Optional[str] = None,
    ):
        """
        Initialize the spectra mapping transform.

        Args:
            mzml_directory: Directory containing mzML files.
            native_id_pattern: Optional custom regex for native ID parsing.
                If not provided, common patterns are tried automatically.
        """
        self._mzml_dir = Path(mzml_directory)
        if not self._mzml_dir.is_dir():
            raise NotADirectoryError(f"mzML directory not found: {mzml_directory}")

        self._native_id_pattern = native_id_pattern
        self._cache = _MzMLCache()
        self._acquisition_starts: dict[str, str] = {}

    def _resolve_mzml_path(self, run_file_name: str) -> Optional[Path]:
        """
        Resolve a run file name to an mzML file path.

        Tries the run file name with common mzML extensions in the configured
        directory.

        Args:
            run_file_name: Run file name (without extension).

        Returns:
            Path to the mzML file, or None if not found.
        """
        for ext in self.MZML_EXTENSIONS:
            candidate = self._mzml_dir / f"{run_file_name}{ext}"
            if candidate.exists():
                return candidate

        # Also try the name as-is (might already have extension)
        candidate = self._mzml_dir / run_file_name
        if candidate.exists():
            return candidate

        return None

    def get_spectrum(
        self,
        run_file_name: str,
        scan_number: Union[int, tuple[int, ...]],
    ) -> tuple[int, np.ndarray, np.ndarray]:
        """
        Get a single spectrum by run file name and scan number.

        Args:
            run_file_name: Name of the MS run (without mzML extension).
            scan_number: Scan number or complete native ID components to retrieve.

        Returns:
            Tuple of (num_peaks, mz_array, intensity_array).
            Returns (0, empty_array, empty_array) if scan or file not found.
        """
        mzml_path = self._resolve_mzml_path(run_file_name)
        if mzml_path is None:
            logger.warning(f"mzML file not found for run: {run_file_name}")
            return 0, np.array([], dtype=np.float32), np.array([], dtype=np.float32)

        return self._cache.get_spectrum(
            str(mzml_path),
            scan_number,
            native_id_pattern=self._native_id_pattern,
        )

    def get_spectra_batch(
        self,
        run_file_name: str,
        scan_numbers: list[Union[int, tuple[int, ...]]],
        position: bool = False,
    ) -> list[tuple[int, np.ndarray, np.ndarray]]:
        """
        Get multiple spectra from the same run file.

        More efficient than calling get_spectrum() repeatedly because the mzML
        file is loaded only once and kept in cache.

        Args:
            run_file_name: Name of the MS run.
            scan_numbers: Scan numbers or complete native ID component tuples.
            position: The values are 0-based file positions (``scan_format`` index).

        Returns:
            List of (num_peaks, mz_array, intensity_array) tuples.
        """
        mzml_path = self._resolve_mzml_path(run_file_name)
        if mzml_path is None:
            logger.warning(f"mzML file not found for run: {run_file_name}")
            empty = np.array([], dtype=np.float32)
            return [(0, empty.copy(), empty.copy()) for _ in scan_numbers]

        mzml_str = str(mzml_path)
        if position:
            return [self._cache.get_spectrum_at(mzml_str, value) for value in scan_numbers]
        return [self._cache.get_spectrum(mzml_str, scan, native_id_pattern=self._native_id_pattern) for scan in scan_numbers]

    def annotate_dataset_psms(
        self,
        dataset,
        batch_size: int = 20,
    ) -> pd.DataFrame:
        """
        Annotate PSM data from a Dataset with spectral arrays.

        Retrieves m/z and intensity arrays for each PSM's scan number from
        the corresponding mzML file. Uses batched processing by run file
        to optimize mzML loading.

        Args:
            dataset: A qpx.Dataset with PSM data.
            batch_size: Number of run files to process in each batch.

        Returns:
            DataFrame with PSM data plus num_peaks, mz_array, intensity_array columns.
        """
        if dataset.psm is None:
            raise ValueError("Dataset does not contain PSM data.")

        annotated_chunks = []
        position = _psm_scan_format(dataset.psm) == "index"

        for run_files, batch_df in dataset.psm.iter_batches(
            partition_by="run_file_name",
            batch_size=batch_size,
        ):
            for run_file in run_files:
                run_mask = batch_df["run_file_name"] == run_file
                run_df = batch_df[run_mask].copy()

                if run_df.empty:
                    continue

                # PSM.scan identifies one spectrum; preserve every component.
                scan_numbers = []
                for scan_val in run_df["scan"]:
                    if isinstance(scan_val, (list, tuple, np.ndarray)):
                        components = tuple(int(value) for value in scan_val)
                        scan_numbers.append(components[0] if len(components) == 1 else components)
                    elif isinstance(scan_val, (int, np.integer)):
                        scan_numbers.append(int(scan_val))
                    else:
                        scan_numbers.append(())

                # Batch-extract spectra for this run
                spectra = self.get_spectra_batch(run_file, scan_numbers, position=position)

                # Add spectral data to DataFrame
                run_df["num_peaks"] = [s[0] for s in spectra]
                run_df["mz_array"] = [s[1].tolist() for s in spectra]
                run_df["intensity_array"] = [s[2].tolist() for s in spectra]

                annotated_chunks.append(run_df)

        if not annotated_chunks:
            logger.warning("No PSMs were annotated with spectral data.")
            return pd.DataFrame()

        result = pd.concat(annotated_chunks, ignore_index=True)
        n_with_spectra = (result["num_peaks"] > 0).sum()
        logger.info(
            f"Annotated {n_with_spectra}/{len(result)} PSMs with spectral data ({n_with_spectra / len(result) * 100:.1f}%)"
        )
        return result

    def write_mz_parquet(
        self,
        dataset,
        output_path: Union[str, Path],
        ms_levels: Optional[list[int]] = None,
    ) -> Path:
        """
        Write spectral data from mzML files to .mz.parquet format.

        Reads mzML files for all runs in the Dataset and writes spectral data
        (m/z arrays, intensity arrays, precursor info, etc.) in the QPX MZ
        schema format.

        Args:
            dataset: A qpx.Dataset with run data (to identify mzML files).
            output_path: Path for the output .mz.parquet file.
            ms_levels: Optional list of MS levels to include (e.g., [1, 2]).
                If None, all MS levels are included.

        Returns:
            Path to the written file.
        """
        from qpx.writers.mz import MzWriter

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if dataset.run is None:
            raise ValueError("Dataset does not contain run data.")

        # Get run file names from the dataset
        run_df = dataset.run.to_df()
        run_file_names = run_df["run_file_name"].tolist()

        total_spectra = 0
        with MzWriter(output_path) as writer:
            for run_name in run_file_names:
                mzml_path = self._resolve_mzml_path(run_name)
                if mzml_path is None:
                    logger.warning("Skipping run %s: mzML file not found", run_name)
                    continue
                total_spectra += self._write_mzml_spectra(mzml_path, run_name, writer, ms_levels)
            self._declare_mz_contents(writer, ms_levels)

        logger.info("Wrote %d spectra from %d runs to %s", total_spectra, len(run_file_names), output_path)
        return output_path

    def write_mz_parquet_from_dir(
        self,
        output_path: Union[str, Path],
        ms_levels: Optional[list[int]] = None,
    ) -> Path:
        """Write spectral data from every mzML in the configured directory.

        Unlike :meth:`write_mz_parquet`, this does not require a Dataset with a
        run view — it scans ``mzml_directory`` directly. Each ``run_file_name``
        is the mzML filename without its (optionally gzipped) extension. This is
        the entry point for converting downloaded raw spectra (e.g. CPTAC PDC
        ``.mzML.gz`` files) into a standalone QPX ``mz.parquet``.

        Args:
            output_path: Path for the output ``.mz.parquet`` file.
            ms_levels: Optional list of MS levels to include (e.g. ``[1, 2]``).
                If ``None``, all MS levels are written.

        Returns:
            Path to the written file.
        """
        from qpx.writers.mz import MzWriter

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        mzml_paths = self._discover_mzml_files()
        if not mzml_paths:
            raise FileNotFoundError(f"No mzML files found in {self._mzml_dir}")

        total_spectra = 0
        with MzWriter(output_path) as writer:
            formats: set[str | None] = set()
            for mzml_path in mzml_paths:
                run_name = self._mzml_run_name(mzml_path)
                total_spectra += self._write_mzml_spectra(mzml_path, run_name, writer, ms_levels, formats)
            if len(formats) == 1 and None not in formats:
                writer.set_scan_format(formats.pop())
            self._declare_mz_contents(writer, ms_levels)

        logger.info("Wrote %d spectra from %d mzML files to %s", total_spectra, len(mzml_paths), output_path)
        return output_path

    def _discover_mzml_files(self) -> list[Path]:
        """List all mzML files (including ``.gz``) in the configured directory."""
        found: set[Path] = set()
        for ext in self.MZML_EXTENSIONS:
            found.update(self._mzml_dir.glob(f"*{ext}"))
        return sorted(found)

    @staticmethod
    def _mzml_run_name(mzml_path: Path) -> str:
        """Run name = filename with its mzML(.gz) extension stripped."""
        name = mzml_path.name
        for ext in (".mzML.gz", ".mzml.gz", ".MZML.gz", ".mzML", ".mzml", ".MZML"):
            if name.endswith(ext):
                return name[: -len(ext)]
        return mzml_path.stem

    def _write_mzml_spectra(self, mzml_path, run_name, writer, ms_levels, formats: set | None = None) -> int:
        """Stream one mzML and write its spectra (optionally filtered by MS level).

        Returns the number of spectra written; *formats* collects their scan formats.
        """
        try:
            import pyopenms as oms
        except ImportError:
            raise ImportError("pyopenms is required for mzML parsing. Install it with: pip install pyopenms")

        logger.info("Processing mzML: %s", mzml_path)

        # Stream spectra rather than loading the whole experiment, so memory stays
        # bounded by one peak-limited batch instead of growing with the file.
        def build_record(spectrum, index):
            if ms_levels is not None and spectrum.getMSLevel() not in ms_levels:
                return None
            if formats is not None:
                formats.add(scan_format_from_native_id(spectrum.getNativeID()))
            return self._build_mz_record(spectrum, run_name, index, oms)

        consumer = _MzRecordWriter(writer, build_record)
        oms.MzMLFile().transform(str(mzml_path), consumer)
        consumer.flush()
        if consumer.acquisition_start:
            self._acquisition_starts[run_name] = consumer.acquisition_start
        return consumer.written

    def _declare_mz_contents(self, writer, ms_levels) -> None:
        """Say in the footer what the view holds and when its runs were acquired."""
        writer.add_footer_metadata(
            {
                b"spectrum_content": b"all_spectra" if ms_levels is None else b"selected_ms_levels",
                b"acquisition_start": json.dumps(self._acquisition_starts, sort_keys=True).encode(),
            }
        )

    @staticmethod
    def _build_mz_record(spectrum, run_name, index, oms) -> dict:
        """Convert one pyopenms spectrum into a QPX mz record.

        ``run_file_name`` and ``scan`` (parsed from the native ID) are emitted
        so the spectrum can be linked back to PSM/feature records.
        """
        mz_array, intensity_array = spectrum.get_peaks()
        native_id = spectrum.getNativeID()
        scan_id = f"{run_name}:{native_id}" if native_id else f"{run_name}:index={index}"

        ms_level = spectrum.getMSLevel()
        precursors = None
        if ms_level >= 2 and spectrum.getPrecursors():
            precursors = []
            for precursor in spectrum.getPrecursors():
                precursors.append(
                    {
                        "selected_ion_mz": float(precursor.getMZ()),
                        "selected_ion_charge": (int(precursor.getCharge()) if precursor.getCharge() > 0 else None),
                        "selected_ion_intensity": (float(precursor.getIntensity()) if precursor.getIntensity() > 0 else None),
                        "isolation_window_target": _isolation_target(precursor),
                        "isolation_window_lower": float(precursor.getIsolationWindowLowerOffset()),
                        "isolation_window_upper": float(precursor.getIsolationWindowUpperOffset()),
                        "spectrum_ref": None,
                    }
                )

        return {
            "id": scan_id,
            "run_file_name": run_name,
            "scan": scan_from_native_id(native_id) or None,
            "ms_level": ms_level,
            "centroid": spectrum.getType() == oms.SpectrumSettings.SpectrumType.CENTROID,
            "scan_start_time": float(spectrum.getRT() / 60.0),  # Convert seconds to minutes
            # timsTOF PASEF MS2 carries its precursor's 1/K0; frames and DIA windows have none.
            "inverse_ion_mobility": _drift_value(spectrum, oms.DriftTimeUnit.VSSC),
            "drift_time": _drift_value(spectrum, oms.DriftTimeUnit.MILLISECOND),
            "compensation_voltage": _drift_value(spectrum, oms.DriftTimeUnit.FAIMS_COMPENSATION_VOLTAGE),
            "ion_injection_time": _ion_injection_time(spectrum),
            "total_ion_current": _total_ion_current(spectrum, intensity_array),
            "precursors": precursors,
            "mz": mz_array,
            "intensity": intensity_array,
            **_peak_ion_mobility(spectrum, oms),
            "cv_params": _spectrum_cv_params(spectrum) or None,
        }

    def close(self):
        """Release cached mzML data to free memory."""
        self._cache.clear()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
