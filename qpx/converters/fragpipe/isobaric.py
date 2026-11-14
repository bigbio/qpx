"""FragPipe reporter-ion channels: annotation.txt and the TMT-Integrator protein report."""

from __future__ import annotations

from pathlib import Path

import duckdb
import pandas as pd

from qpx.converters.channel_labels import normalize_label
from qpx.converters.utils import safe_float
from qpx.core.sql import escape_path

TMT_REPORT_INTENSITY = "tmt_integrator_log2_abundance"
_MS1_INTENSITY = "ms1_intensity"


def load_channel_annotation(path: str | Path) -> dict[str, str]:
    """Map each channel column named by a FragPipe annotation.txt to its canonical label."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    pairs = [line.split(maxsplit=1) for line in lines if line.strip()]
    channels = {pair[1].strip(): normalize_label(pair[0]) for pair in pairs if len(pair) == 2}
    if not channels:
        raise ValueError(f"{path} names no channels; each line is '<channel> <sample>', as '126 sample_1'")
    return channels


def require_channel_columns(channels: dict[str, str], columns: set[str], path: str) -> None:
    """Stop when the report does not hold the columns the annotation names."""
    missing = sorted(set(channels) - columns)
    if missing:
        raise ValueError(f"{path} has no channel column for annotated sample(s): {', '.join(missing)}")


def plex_runs(psm_path: str) -> list[str]:
    """The raw runs of one plex, which FragPipe names in the psm.tsv spectrum ids."""
    rows = duckdb.sql(
        f"""SELECT DISTINCT regexp_replace("Spectrum", '(\\.\\d+){{2}}\\.\\d+$', '') AS run
        FROM read_csv('{escape_path(psm_path)}', delim='\t', header=true, all_varchar=true) ORDER BY run"""
    ).fetchall()
    return [run for (run,) in rows if run]


def channel_intensities(row, channels: dict[str, str]) -> list[dict] | None:
    """One entry per channel, null for FragPipe's zero; ``None`` when no channel was quantified."""
    entries = []
    for column, label in channels.items():
        value = safe_float(row.get(column))
        entries.append({"label": label, "intensity": value if value is not None and value > 0 else None})
    return entries if any(entry["intensity"] is not None for entry in entries) else None


def ms1_additional_intensity(row, intensities: list[dict]) -> list[dict] | None:
    """The precursor's MS1 intensity, once per feature under its first quantified channel."""
    ms1 = safe_float(row.get("Intensity"))
    if ms1 is None or ms1 <= 0:
        return None
    label = next(entry["label"] for entry in intensities if entry["intensity"] is not None)
    return [{"label": label, "intensities": [{"intensity_name": _MS1_INTENSITY, "intensity_value": ms1}]}]


def load_tmt_report(path: str | Path, channels: dict[str, str]) -> dict[str, dict[str, float]]:
    """Read a TMT-Integrator ``abundance_protein_*.tsv`` as protein -> label -> log2 abundance."""
    report = pd.read_csv(path, sep="\t", low_memory=False)
    if "Gene" not in report.columns or "Index" not in report.columns:
        raise ValueError(f"{path} is not a protein-level TMT-Integrator report (abundance_protein_*.tsv)")
    columns = {column: label for column, label in channels.items() if column in report.columns}
    if not columns:
        raise ValueError(f"{path} has no column for the samples of the annotation")
    return {
        protein: {columns[column]: float(value) for column, value in values.items() if pd.notna(value)}
        for protein, values in report.set_index("Index")[list(columns)].to_dict("index").items()
    }
