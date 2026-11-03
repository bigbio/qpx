"""DIA-NN ``*.xic.parquet`` exports -> one QPX ``xic.parquet`` view."""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from qpx.converters.diann.constants import to_proforma
from qpx.core.data import XicSchema
from qpx.writers import XicWriter

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
    """Split a DIA-NN ``feature`` value into QPX transition fields."""
    if label == "ms1":
        return {
            "transition_type": "ms1",
            "fragment_type": None,
            "fragment_number": None,
            "fragment_charge": None,
            "neutral_loss": None,
        }
    match = _FRAGMENT.match(label)
    if match is None:
        return {
            "transition_type": "fragment",
            "fragment_type": None,
            "fragment_number": None,
            "fragment_charge": None,
            "neutral_loss": None,
        }
    fragment_type, fragment_number, fragment_charge, neutral_loss = match.groups()
    return {
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


def convert_diann_xic(
    xic_paths: Iterable[str | Path],
    output_path: str | Path,
    *,
    label: str = "LFQ",
    compression: str = "zstd",
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
    compression:
        Parquet compression codec passed to :class:`~qpx.writers.xic.XicWriter`.

    Returns
    -------
    int
        Number of chromatogram traces written.
    """
    paths = [Path(path) for path in xic_paths]
    if not paths:
        raise ValueError("No DIA-NN XIC files were provided")

    runs = [_run_name(path) for path in paths]
    duplicates = sorted({run for run in runs if runs.count(run) > 1})
    if duplicates:
        raise ValueError(f"Duplicate run names across XIC files: {', '.join(duplicates)}")

    for path in paths:
        columns = set(pq.read_schema(path).names)
        missing = _REQUIRED_COLUMNS - columns
        if missing:
            raise ValueError(f"{path} is missing DIA-NN XIC column(s): {', '.join(sorted(missing))}")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()
    precursor_cache: dict[str, tuple[str, int]] = {}
    transition_cache: dict[str, dict] = {}
    total = 0
    try:
        with XicWriter(output_path, compression=compression) as writer:
            for path, run in zip(paths, runs):
                traces = con.execute(_trace_query(path)).fetchall()
                rows = []
                for precursor_id, transition, rt, intensity in traces:
                    precursor = precursor_cache.get(precursor_id)
                    if precursor is None:
                        precursor = _parse_precursor(precursor_id)
                        precursor_cache[precursor_id] = precursor
                    annotation = transition_cache.get(transition)
                    if annotation is None:
                        annotation = _parse_transition(transition)
                        transition_cache[transition] = annotation
                    peptidoform, charge = precursor
                    rows.append(
                        {
                            "xic_id": f"{run}|{precursor_id}|{transition}",
                            "run_file_name": run,
                            "label": label,
                            "feature_id": None,
                            "precursor_id": precursor_id,
                            "peptidoform": peptidoform,
                            "charge": charge,
                            "transition": transition,
                            **annotation,
                            "product_mz": None,
                            "library_intensity": None,
                            "rt": rt,
                            "intensity": intensity,
                            "is_decoy": None,
                            "cv_params": None,
                        }
                    )

                table = pa.Table.from_pylist(rows, schema=XicSchema.get_arrow_schema())
                result = XicSchema.validate_full(table, strict=True)
                if result.errors:
                    raise ValueError(
                        f"XIC schema validation failed for {path}:\n" + "\n".join(issue.message for issue in result.errors)
                    )
                writer.write_table(table)
                total += table.num_rows
                logger.info("%s: %d traces", path.name, table.num_rows)
    finally:
        con.close()

    logger.info("Wrote %s (%d traces across %d runs)", output_path, total, len(paths))
    return total
