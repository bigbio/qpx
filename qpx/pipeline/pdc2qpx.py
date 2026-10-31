"""pdc2qpx — one-shot PDC/CPTAC download + QPX conversion.

Given a PDC study accession, download its CDAP ``.psm`` files (and optionally the
mzML spectra), then convert them into a QPX dataset:

- ``.psm`` files  -> psm / feature / pg / dataset / ontology / provenance views
  (via :class:`qpx.converters.cdap.CdapConverter`)
- mzML files      -> a full-spectra ``mz.parquet``
  (via :meth:`qpx.transforms.spectra_mapping.SpectraMappingTransform.write_mz_parquet_from_dir`)

The download layer is :mod:`pridepy` (an optional dependency, ``qpx[pdc]``).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def _download(study: str, file_type: str, download_dir: Path, download_threads: int) -> None:
    """Download one PDC file type into ``download_dir/<study>/`` via pridepy."""
    try:
        from pridepy.pdc.downloader import download_pdc_files
    except ImportError:
        raise ImportError("pridepy is required for pdc2qpx. Install it with: pip install qpx[pdc]") from None

    logger.info("Downloading PDC %s file_type=%s -> %s", study, file_type, download_dir)
    download_pdc_files(
        accession=study,
        file_type=file_type,
        output_folder=str(download_dir),
        skip_if_downloaded_already=True,
        checksum_check=True,
        download_threads=download_threads,
        retry=True,
    )


# Core CDAP columns that identify a CPTAC/PDC CDAP .psm header.
_CDAP_SIGNATURE_COLUMNS = ("FileName", "ScanNum", "PeptideSequence", "Protein")


def _assert_cdap_psm(study_dir: Path) -> None:
    """Pre-flight check that *study_dir* holds CDAP-format ``.psm`` files.

    Fails early with a clear message -- before the DuckDB conversion -- when the
    study provides no ``.psm`` files or when they are not CDAP output, instead of
    letting :class:`~qpx.converters.cdap.CdapConverter` fail later with an opaque
    "column not found" SQL error.
    """
    psm_files = sorted(study_dir.glob("*.psm")) if study_dir.is_dir() else []
    if not psm_files:
        raise FileNotFoundError(
            f"No .psm files found in {study_dir}. pdc2qpx needs CDAP .psm output; this PDC study may not provide harmonized PSMs."
        )
    with open(psm_files[0], encoding="utf-8") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
    missing = [col for col in _CDAP_SIGNATURE_COLUMNS if col not in header]
    if missing:
        raise ValueError(
            f"{psm_files[0].name} is not a CDAP-format .psm (missing columns: "
            f"{missing}); pdc2qpx only supports CPTAC/PDC CDAP output."
        )


def _fetch_run_to_plex(study: str, download_threads: int) -> Optional[dict[str, dict]]:
    """The study's run -> plex map; a failure only logs a warning, and fraction names group the runs instead."""
    try:
        from qpx.pdc.metadata import map_runs_to_plex

        return map_runs_to_plex(study, threads=download_threads)
    except (OSError, RuntimeError, ValueError, KeyError, TypeError) as exc:
        logger.warning("PDC run->plex map unavailable for %s: %s", study, exc)
        return None


def _add_pdc_metadata(
    study: str, output_folder: Path, download_threads: int, outputs: dict[str, Path], run_to_plex: Optional[dict[str, dict]]
) -> None:
    """Build the sample/run views from PDC metadata into ``outputs``; a failure only logs a warning."""
    try:
        from qpx.converters.cdap.pdc_sample_run import build_sample_run_from_pdc

        built = build_sample_run_from_pdc(study, output_folder, prefix=study, threads=download_threads, run_to_plex=run_to_plex)
        if built:
            outputs.update(built)
    except (OSError, RuntimeError, ValueError, KeyError, TypeError) as exc:
        # Best-effort: network (OSError/URLError), GraphQL (RuntimeError), or
        # data-shape failures must not break the main conversion. The builder
        # writes via temp files + atomic rename, so a failure here never
        # clobbers a previously-good sample/run parquet.
        logger.warning("PDC metadata sample/run build skipped for %s: %s", study, exc)


def run_pdc2qpx(
    study: str,
    download_dir: Path,
    output_folder: Path,
    *,
    include_spectra: bool = False,
    include_metadata: bool = True,
    ms_levels: Optional[list[int]] = None,
    max_cpus: int | None = None,
    max_memory: str | None = None,
    download_threads: int = 24,
    skip_download: bool = False,
) -> dict[str, Path]:
    """Download a PDC study and convert it to a QPX dataset.

    Args:
        study: PDC study accession (e.g. ``"PDC000109"``).
        download_dir: Directory for downloaded files; PDC writes them under
            ``download_dir/<study>/``.
        output_folder: Directory for generated QPX parquet files.
        include_spectra: Also download mzML and produce a full-spectra
            ``<study>.mz.parquet``.
        include_metadata: Build ``<study>.sample.parquet`` and
            ``<study>.run.parquet`` from PDC GraphQL metadata (channel ->
            biological-sample mapping). On by default; metadata failures are
            logged and skipped without breaking the conversion.
        ms_levels: MS levels for the mz view (``None`` = all). Precursor-level
            LFQ reanalysis needs MS1 as well as MS2.
        max_cpus: Threads for the CDAP conversion.
        max_memory: Memory limit for the CDAP conversion (e.g. ``"16GB"``).
        download_threads: Parallel HTTP Range threads per file for downloads.
        skip_download: Skip downloading and use files already present under
            ``download_dir/<study>/`` (useful for re-runs and tests).

    Returns:
        Mapping of produced outputs: ``{"base": output_folder}`` and, when
        *include_spectra* is set, ``{"mz": <study>.mz.parquet}``.
    """
    download_dir = Path(download_dir)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    # PDC downloads land under <download_dir>/<study>/. The same directory holds
    # both .psm and .mzML.gz files; the CDAP and mz converters glob their own.
    study_dir = download_dir / study

    # 1. Base QPX from CDAP .psm files.
    if not skip_download:
        _download(study, "psm", download_dir, download_threads)
    _assert_cdap_psm(study_dir)

    from qpx.converters.cdap import CdapConverter
    from qpx.converters.cdap.pdc_sample_run import run_units_from_plexes

    run_to_plex = _fetch_run_to_plex(study, download_threads) if include_metadata else None
    logger.info("Converting CDAP .psm -> QPX base views (%s)", study)
    converter = CdapConverter(max_memory=max_memory, max_cpus=max_cpus)
    converter.convert(
        psm_dir=study_dir,
        output_folder=output_folder,
        output_prefix=study,
        project_accession=study,
        run_units=run_units_from_plexes(run_to_plex) if run_to_plex else None,
    )
    outputs: dict[str, Path] = {"base": output_folder}

    # 2. Optional sample/run views from PDC metadata (channel -> sample mapping).
    if include_metadata:
        _add_pdc_metadata(study, output_folder, download_threads, outputs, run_to_plex)

    # 3. Optional full-spectra mz.parquet from mzML files.
    if include_spectra:
        if not skip_download:
            _download(study, "mzml", download_dir, download_threads)
        if not any(p.name.lower().endswith((".mzml", ".mzml.gz")) for p in study_dir.iterdir()):
            raise FileNotFoundError(f"--include-spectra set but no mzML files found in {study_dir}")

        from qpx.transforms.spectra_mapping import SpectraMappingTransform

        mz_out = output_folder / f"{study}.mz.parquet"
        logger.info("Converting mzML -> %s", mz_out)
        SpectraMappingTransform(mzml_directory=study_dir).write_mz_parquet_from_dir(mz_out, ms_levels=ms_levels)
        outputs["mz"] = mz_out

    logger.info("pdc2qpx complete for %s -> %s", study, output_folder)
    return outputs


def parse_accessions(accession: str) -> list[str]:
    """Resolve PDC study accessions from a single ``-a`` argument.

    Accepts a single ID, comma/newline-separated IDs, or a path to a CSV with a
    ``pdc_study_id``/``pdc_id`` column. Mirrors the input forms of pridepy's
    ``parse_accessions`` so ``pdc2qpx`` and ``pridepy download-pdc-files`` take
    the same ``-a`` argument -- but keeps parsing dependency-free (pridepy is an
    optional extra). Blank lines and ``#`` comments are ignored; order is
    preserved and duplicates dropped.
    """
    source = Path(accession).expanduser()
    if source.is_file() and source.suffix.lower() == ".csv":
        import csv

        with source.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            column = next((c for c in ("pdc_study_id", "pdc_id") if c in (reader.fieldnames or [])), None)
            if column is None:
                raise ValueError(f"CSV must contain a pdc_study_id or pdc_id column: {source}")
            raw = [str(row.get(column) or "") for row in reader]
    else:
        raw = accession.replace(",", "\n").split()
    cleaned = [item.strip() for item in raw if item.strip() and not item.strip().startswith("#")]
    if not cleaned:
        raise ValueError(f"No PDC study accession found in: {accession!r}")
    return list(dict.fromkeys(cleaned))


def run_pdc2qpx_batch(
    studies: list[str],
    download_dir: Path,
    output_root: Path,
    *,
    include_spectra: bool = False,
    include_metadata: bool = True,
    ms_levels: Optional[list[int]] = None,
    max_cpus: int | None = None,
    max_memory: str | None = None,
    download_threads: int = 24,
    skip_download: bool = False,
    continue_on_error: bool = True,
) -> dict[str, dict]:
    """Run :func:`run_pdc2qpx` for each study into ``output_root/<study>/``.

    Returns ``{study: {"status": "ok"|"failed", "outputs": dict|None,
    "error": str|None}}``. With *continue_on_error* (default) a study's expected
    failure (missing/non-CDAP ``.psm``, network) is recorded and the batch
    continues; set it ``False`` to abort on the first failure.
    """
    output_root = Path(output_root)
    results: dict[str, dict] = {}
    for index, study in enumerate(studies, start=1):
        logger.info("pdc2qpx %d/%d: %s", index, len(studies), study)
        try:
            outputs = run_pdc2qpx(
                study,
                download_dir,
                output_root / study,
                include_spectra=include_spectra,
                include_metadata=include_metadata,
                ms_levels=ms_levels,
                max_cpus=max_cpus,
                max_memory=max_memory,
                download_threads=download_threads,
                skip_download=skip_download,
            )
            results[study] = {"status": "ok", "outputs": outputs, "error": None}
        except (OSError, RuntimeError, ValueError, KeyError, TypeError) as exc:
            if not continue_on_error:
                raise
            logger.warning("pdc2qpx: study %s failed: %s", study, exc)
            results[study] = {"status": "failed", "outputs": None, "error": str(exc)}

    n_ok = sum(1 for record in results.values() if record["status"] == "ok")
    logger.info("pdc2qpx batch complete: %d/%d studies ok", n_ok, len(studies))
    return results
