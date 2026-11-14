"""DIA-NN ``*.xic.parquet`` exports -> one QPX ``xic.parquet`` view."""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable
from dataclasses import replace
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from qpx.converters.diann.constants import to_proforma
from qpx.core.data import XicSchema
from qpx.core.engine import create_converter_connection
from qpx.writers import XicExtraction, XicWriter

logger = logging.getLogger(__name__)

_FRAGMENT = re.compile(r"^([by])(\d+)\^(\d+)(.*)$")
_CHARGE = re.compile(r"^(.*?)(\d+)$")
_REQUIRED_COLUMNS = {"pr", "feature", "rt", "value"}


def _run_name(path: Path) -> str:
    """The run name encoded in a DIA-NN XIC file name."""
    for suffix in (".xic.parquet", ".parquet"):
        if path.name.endswith(suffix):
            return path.name[: -len(suffix)]
    return path.stem


def _parse_precursor(precursor_id: str) -> tuple[str, int]:
    """Split a DIA-NN ``pr`` value into a ProForma peptidoform and charge."""
    match = _CHARGE.match(precursor_id)
    if match is None:
        raise ValueError(f"Cannot parse the precursor charge from DIA-NN pr value {precursor_id!r}")
    modified_sequence, charge = match.groups()
    return to_proforma(modified_sequence), int(charge)


def _parse_transition(label: str) -> dict:
    """Split a DIA-NN ``feature`` value into QPX transition fields.

    DIA-NN writes the charge before a neutral loss (``y9^1-H2O``) while mzPAF
    writes the loss before the charge (``y9-H2O^1``); fragment labels are
    normalised to the mzPAF order.
    """
    if label == "ms1":
        return {
            "transition": "ms1",
            "transition_type": "ms1",
            "fragment_type": None,
            "fragment_number": None,
            "fragment_charge": None,
            "neutral_loss": None,
        }
    match = _FRAGMENT.match(label)
    if match is None:
        return {
            "transition": label,
            "transition_type": "fragment",
            "fragment_type": None,
            "fragment_number": None,
            "fragment_charge": None,
            "neutral_loss": None,
        }
    fragment_type, fragment_number, fragment_charge, neutral_loss = match.groups()
    transition = f"{fragment_type}{fragment_number}{neutral_loss or ''}^{fragment_charge}"
    return {
        "transition": transition,
        "transition_type": "fragment",
        "fragment_type": fragment_type,
        "fragment_number": int(fragment_number),
        "fragment_charge": int(fragment_charge),
        "neutral_loss": neutral_loss or None,
    }


def _trace_query(path: Path) -> str:
    """Aggregate one DIA-NN XIC file into one row per (precursor, transition)."""
    location = path.as_posix().replace("'", "''")
    return f"""
        SELECT pr, feature,
               list(CAST(rt * 60.0 AS REAL) ORDER BY rt) AS rt,
               list(CAST(value AS REAL) ORDER BY rt) AS intensity
        FROM read_parquet('{location}')
        WHERE feature <> 'index'
        GROUP BY pr, feature
    """


def _export_runs(paths: list[Path]) -> list[str]:
    """The run of each DIA-NN XIC export, which must be distinct and carry the trace columns."""
    if not paths:
        raise ValueError("No DIA-NN XIC files were provided")
    runs = [_run_name(path) for path in paths]
    duplicates = sorted({run for run in runs if runs.count(run) > 1})
    if duplicates:
        raise ValueError(f"Duplicate run names across XIC files: {', '.join(duplicates)}")
    for path in paths:
        missing = _REQUIRED_COLUMNS - set(pq.read_schema(path).names)
        if missing:
            raise ValueError(f"{path} is missing DIA-NN XIC column(s): {', '.join(sorted(missing))}")
    return runs


# XIC fields a DIA-NN export does not carry.
_NOT_EXPORTED = dict.fromkeys(
    (
        "feature_ids",
        "feature_scores",
        "compensation_voltage",
        "target_mz",
        "target_rt",
        "rt_start",
        "rt_stop",
        "target_ion_mobility",
        "ion_mobility_start",
        "ion_mobility_stop",
        "library_intensity",
        "used_for_quantification",
        "mz",
        "is_decoy",
        "additional_scores",
        "cv_params",
    )
)


def _run_table(con, path: Path, run: str, label: str, parsed: tuple[dict, dict]) -> pa.Table:
    """The traces of one run's export as a validated XIC table.

    *parsed* holds the precursors and transitions already parsed, keyed by their
    DIA-NN values, and is shared across the runs.
    """
    precursors, transitions = parsed
    rows = []
    for precursor_id, transition, rt, intensity in con.execute(_trace_query(path)).fetchall():
        if precursor_id not in precursors:
            precursors[precursor_id] = _parse_precursor(precursor_id)
        if transition not in transitions:
            transitions[transition] = _parse_transition(transition)
        peptidoform, charge = precursors[precursor_id]
        annotation = transitions[transition]
        rows.append(
            {
                "xic_id": f"{run}|{precursor_id}|{annotation['transition']}",
                "run_file_name": run,
                "label": label,
                "precursor_id": precursor_id,
                "peptidoform": peptidoform,
                "charge": charge,
                **annotation,
                "isotope": 0,
                "rt": rt,
                "intensity": intensity,
                **_NOT_EXPORTED,
            }
        )
    table = pa.Table.from_pylist(rows, schema=XicSchema.get_arrow_schema())
    result = XicSchema.validate_full(table, strict=True)
    if result.errors:
        raise ValueError(f"XIC schema validation failed for {path}:\n" + "\n".join(issue.message for issue in result.errors))
    return table


# DIA-NN exports the traces around the apex of the target precursors by default.
_DEFAULT_EXTRACTION = XicExtraction(selection="targets", window="apex")


def convert_diann_xic(
    xic_paths: Iterable[str | Path],
    output_path: str | Path,
    *,
    label: str = "LFQ",
    extraction: XicExtraction | None = None,
    duckdb_memory: str | None = None,
    duckdb_threads: int | None = None,
) -> int:
    """Merge DIA-NN per-run XIC exports into one dataset-level QPX XIC file.

    Parameters
    ----------
    xic_paths:
        DIA-NN ``*.xic.parquet`` files, one per run.
    output_path:
        Destination of the QPX ``{PREFIX}.xic.parquet`` file.
    label:
        Channel label written for every trace. DIA-NN XIC exports are
        label-free, so the default is ``"LFQ"``.
    extraction:
        What the export contains (``selection``: ``targets``,
        ``targets_and_decoys`` or ``all``), whether its traces cover the apex
        region or the full extraction window (``window``), and, when known, the
        extraction-window width in seconds and the q-value threshold that
        selected the precursors. Written to the file footer; defaults to the
        targets around the apex. The transitions are always declared as mzPAF,
        unsmoothed and with their zero points kept.
    duckdb_memory, duckdb_threads:
        DuckDB memory limit (e.g. ``"16GB"``) and thread count for reading the exports.

    Returns
    -------
    int
        Number of chromatogram traces written.
    """
    paths = [Path(path) for path in xic_paths]
    runs = _export_runs(paths)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    extraction = replace(extraction or _DEFAULT_EXTRACTION, transition_format="mzPAF", smoothed=False, zero_points="kept")

    parsed: tuple[dict, dict] = ({}, {})
    total = 0
    with (
        create_converter_connection(memory_limit=duckdb_memory, threads=duckdb_threads) as con,
        XicWriter(output_path, extraction=extraction) as writer,
    ):
        for path, run in zip(paths, runs):
            table = _run_table(con, path, run, label, parsed)
            writer.write_table(table)
            total += table.num_rows
            logger.info("%s: %d traces", path.name, table.num_rows)

    logger.info("Wrote %s (%d traces across %d runs)", output_path, total, len(paths))
    return total
