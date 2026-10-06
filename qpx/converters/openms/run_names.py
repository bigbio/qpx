"""Normalize native OpenMS run references to the SDRF basename convention."""

from collections.abc import Callable
from functools import cache, partial
from pathlib import Path

import pyarrow as pa

from qpx.converters.channel_labels import experiment_runs_from_sdrf
from qpx.converters.openms.psm_runs import PsmRunResolver, build_psm_run_lookup
from qpx.core.constants import FEATURE, PG, PSM
from qpx.core.files import run_file_stem

RunNormalizer = Callable[[pa.Table, str], pa.Table]

_RUN_COLUMNS = {
    FEATURE: ("run_file_name", "id_run_file_name"),
    PSM: ("run_file_name",),
    PG: ("grouped_runs", "run_file_name"),
}


def normalize_run_names(
    table: pa.Table,
    view: str,
    known_runs: frozenset[str] = frozenset(),
    psm_run_resolver: PsmRunResolver | None = None,
) -> pa.Table:
    """Normalize stored run references before deriving missing identities.

    Only acquisition-file suffixes are removed, so meaningful dots in run names
    survive. Missing references and grouped-run order are retained; this does
    not change supplied IDs and cross-references. Declared canonical run names
    are preserved; optional PSM recovery runs last because its result is already
    canonical.
    """

    @cache
    def normalize(run):
        return run if run in known_runs else run_file_stem(run)

    for name in _RUN_COLUMNS.get(view, ()):
        if name not in table.column_names:
            continue
        values = table.column(name).to_pylist()
        if name == "grouped_runs":
            values = [
                [normalize(run) if run is not None else None for run in group] if group is not None else None for group in values
            ]
        else:
            values = [normalize(run) if run is not None else None for run in values]
        field = table.schema.field(name)
        table = table.set_column(table.schema.get_field_index(name), field, pa.array(values, type=field.type))
    if view == PSM and psm_run_resolver is not None:
        table = psm_run_resolver.apply(table)
    return table


def build_run_normalizer(
    sdrf_path: str | None,
    consensusxml_path: str | Path | None,
    maplist: dict[int, dict[str, str]],
    has_psms: bool,
) -> RunNormalizer:
    """Prepare declared run names and optional companion evidence once per import."""
    experiment_runs = experiment_runs_from_sdrf(sdrf_path) or {}
    known_runs = {run for runs in experiment_runs.values() for run in runs}
    known_runs.update(run_file_stem(entry["name"]) for entry in maplist.values() if entry.get("name"))
    resolver = build_psm_run_lookup(consensusxml_path) if consensusxml_path and has_psms else None
    return partial(normalize_run_names, known_runs=frozenset(known_runs), psm_run_resolver=resolver)
