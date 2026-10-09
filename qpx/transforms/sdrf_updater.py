"""Update QPX dataset metadata from a revised SDRF file.

Re-generates sample.parquet (and optionally run.parquet metadata fields)
from a new SDRF, while protecting fields that affect identification or
quantification.

Protected fields (changes blocked with warning):
    - instrument, enzymes, modification_parameters, fraction
    - label / sample-run channel mapping
    - comment[data file] (run identity)

Allowed fields (sample metadata updates):
    - disease, organism_part, cell_type, cell_line, sex, age,
      developmental_stage, ancestry, individual, sample_description
    - Any additional characteristics[*] columns
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import pyarrow.parquet as pq

from qpx.transforms.utils import discover_qpx_file_prefix

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)


def _read_sdrf_normalized(path: str | Path) -> "pd.DataFrame":
    """Read SDRF and normalize column names to lowercase."""
    import pandas as pd

    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lower().strip() for c in df.columns]
    return df


def _check_data_files(old: "pd.DataFrame", new: "pd.DataFrame") -> list[str]:
    """Check that data file references haven't changed."""
    col = "comment[data file]"
    if col not in old.columns or col not in new.columns:
        return []
    warnings: list[str] = []
    old_files = set(old[col].dropna().unique())
    new_files = set(new[col].dropna().unique())
    if old_files == new_files:
        return warnings
    removed = old_files - new_files
    added = new_files - old_files
    if removed:
        warnings.append(
            f"BLOCKED: {len(removed)} data files removed from SDRF — "
            f"this would break run references. Removed: {sorted(removed)[:5]}"
        )
    if added:
        warnings.append(f"WARNING: {len(added)} new data files in SDRF — these won't have corresponding feature/PG data.")
    return warnings


def _check_label_mapping(old: "pd.DataFrame", new: "pd.DataFrame") -> list[str]:
    """Check that sample-run-label mapping is unchanged."""
    col = "comment[label]"
    required = ["source name", "comment[data file]", col]
    if not all(k in old.columns and k in new.columns for k in required):
        return []
    old_map = set(zip(old["source name"], old["comment[data file]"], old[col]))
    new_map = set(zip(new["source name"], new["comment[data file]"], new[col]))
    if old_map != new_map:
        return ["BLOCKED: sample-run-label mapping changed — this would invalidate intensity assignments."]
    return []


def _check_run_fields(old: "pd.DataFrame", new: "pd.DataFrame") -> list[str]:
    """Check protected run-level fields (instrument, enzymes, etc.)."""
    warnings: list[str] = []
    protected_cols = [
        "comment[instrument]",
        "comment[cleavage agent details]",
        "comment[fraction identifier]",
        "comment[dissociation method]",
    ]
    for col in protected_cols:
        if col in old.columns and col in new.columns:
            old_vals = set(old[col].dropna().unique())
            new_vals = set(new[col].dropna().unique())
            if old_vals != new_vals:
                warnings.append(
                    f"BLOCKED: '{col}' changed from {old_vals} to {new_vals} — this affects identification/quantification."
                )
    return warnings


def _check_modifications(old: "pd.DataFrame", new: "pd.DataFrame") -> list[str]:
    """Check that modification parameters are unchanged."""
    mod_cols_old = [c for c in old.columns if "modification parameters" in c]
    mod_cols_new = [c for c in new.columns if "modification parameters" in c]
    if not mod_cols_old and not mod_cols_new:
        return []
    old_mods = {v for c in mod_cols_old for v in old[c].dropna().unique()}
    new_mods = {v for c in mod_cols_new for v in new[c].dropna().unique()}
    if old_mods != new_mods:
        return ["BLOCKED: modification parameters changed — this would invalidate PSM/feature identifications."]
    return []


def _check_protected_changes(
    old_sdrf_path: str | Path | None,
    new_sdrf_path: str | Path,
) -> list[str]:
    """Compare old vs new SDRF for protected field changes.

    Returns list of warning messages. Empty list means safe to proceed.
    If old_sdrf_path is None, skip comparison (first-time update).
    """
    if old_sdrf_path is None:
        return []

    old = _read_sdrf_normalized(old_sdrf_path)
    new = _read_sdrf_normalized(new_sdrf_path)

    warnings: list[str] = []
    warnings.extend(_check_data_files(old, new))
    warnings.extend(_check_label_mapping(old, new))
    warnings.extend(_check_run_fields(old, new))
    warnings.extend(_check_modifications(old, new))
    return warnings


def update_metadata(
    dataset_path: Path,
    new_sdrf_path: Path,
    old_sdrf_path: Path | None = None,
    force: bool = False,
) -> dict:
    """Update sample.parquet and run.parquet metadata from a new SDRF.

    Args:
        dataset_path: QPX dataset directory.
        new_sdrf_path: Path to the updated SDRF TSV file.
        old_sdrf_path: Path to the original SDRF (for safety checks).
            If None, skips protected-field validation.
        force: If True, apply changes even with protected-field warnings.

    Returns:
        Dict with stats: {"samples_updated", "runs_updated", "warnings"}.
    """
    from qpx.converters.sdrf import SdrfConverter

    # Step 1: Check for protected field changes
    check_warnings = _check_protected_changes(old_sdrf_path, new_sdrf_path)
    blocked = [w for w in check_warnings if w.startswith("BLOCKED")]

    if blocked and not force:
        for w in check_warnings:
            logger.warning(w)
        raise ValueError(
            f"{len(blocked)} protected field(s) changed in the new SDRF. Use --force to override (data integrity not guaranteed)."
        )

    for w in check_warnings:
        if w.startswith("WARNING"):
            logger.warning(w)

    # Step 2: Re-generate sample.parquet and run.parquet from new SDRF
    file_prefix = discover_qpx_file_prefix(dataset_path)
    sample_path = dataset_path / f"{file_prefix}.sample.parquet"
    run_path = dataset_path / f"{file_prefix}.run.parquet"

    # Back up originals
    if sample_path.exists():
        backup = sample_path.with_name(f"{sample_path.name}.bak")
        if not backup.exists():
            import shutil

            shutil.copy2(sample_path, backup)
            logger.info("Backed up original sample.parquet → %s", backup.name)

    if run_path.exists():
        backup = run_path.with_name(f"{run_path.name}.bak")
        if not backup.exists():
            import shutil

            shutil.copy2(run_path, backup)
            logger.info("Backed up original run.parquet → %s", backup.name)

    # Step 3: Convert new SDRF
    converter = SdrfConverter()
    converter.convert(
        sdrf_path=str(new_sdrf_path),
        sample_output=str(sample_path),
        run_output=str(run_path),
    )

    # Count results
    n_samples = len(pq.read_table(str(sample_path)))
    n_runs = len(pq.read_table(str(run_path)))

    logger.info(
        "Updated metadata: %d samples, %d runs from %s",
        n_samples,
        n_runs,
        new_sdrf_path.name,
    )

    return {
        "samples_updated": n_samples,
        "runs_updated": n_runs,
        "warnings": check_warnings,
    }
