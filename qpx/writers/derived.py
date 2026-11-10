"""Flat columns derived at write time, so every converter fills them the same way."""

from __future__ import annotations

import re
from collections.abc import Callable

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc

from qpx.core.data.identity import derive_id
from qpx.core.scores import normalize_score_name

_QVALUE = re.compile(r"(?:^|_)q_?value$")
_NOT_PSM_LEVEL = re.compile(r"(?:^|_)(?:peptide|peptidoform|protein|pg|gg|global|precursor|lib|translated|glycan|site)(?:_|$)")


def is_psm_qvalue(name: str) -> bool:
    """``qvalue``, ``q_value``, ``psm_level_qvalue``, ``percolator_q_value``: a q-value of no higher level."""
    return bool(_QVALUE.search(name)) and not _NOT_PSM_LEVEL.search(name)


def _in_names(*names: str) -> Callable[[str], bool]:
    return frozenset(names).__contains__


# Typed column -> additional score names that fill it, normalized to snake_case (scores.md).
PSM_SCORE_COLUMNS: dict[str, Callable[[str], bool]] = {"qvalue": is_psm_qvalue}
FEATURE_SCORE_COLUMNS: dict[str, Callable[[str], bool]] = {
    "peptide_qvalue": _in_names("peptide_qvalue", "peptide_level_qvalue"),
    "precursor_qvalue": _in_names("precursor_qvalue", "qvalue", "q_value", "psm_level_qvalue"),
    "global_qvalue": _in_names("global_qvalue"),
}


_IRT_PARAMS = ["irt", "empirical_irt"]


def _score_values(scores: pa.Array, wanted: Callable[[str], bool]) -> np.ndarray:
    """The first wanted score of each row of a ``list<score>`` array, NaN where it has none."""
    values = np.full(len(scores), np.nan)
    entries = scores.flatten()
    if len(entries) == 0:
        return values
    encoded = pc.call_function("dictionary_encode", [entries.field("score_name")])
    keep = np.array([wanted(normalize_score_name(name)) for name in encoded.dictionary.to_pylist()], dtype=bool)
    indices = encoded.indices.fill_null(-1).to_numpy(zero_copy_only=False)
    chosen = np.flatnonzero(np.append(keep, False)[indices])
    if len(chosen):
        parents = pc.call_function("list_parent_indices", [scores]).to_numpy()
        numbers = entries.field("score_value").to_numpy(zero_copy_only=False).astype(np.float64)
        # Reversed so that the first wanted score of a row is the one left in place.
        values[parents[chosen[::-1]]] = numbers[chosen[::-1]]
    return values


def fill_score_columns(table: pa.Table, columns: dict[str, Callable[[str], bool]]) -> pa.Table:
    """Fill each null of the typed *columns* from ``additional_scores``."""
    if "additional_scores" not in table.schema.names or not table.num_rows:
        return table
    scores = pa.concat_arrays(table.column("additional_scores").chunks)
    for name, wanted in columns.items():
        index = table.schema.get_field_index(name)
        if index < 0:
            continue
        derived = _score_values(scores, wanted)
        typed = pa.array(derived, mask=np.isnan(derived), type=table.schema.field(index).type)
        filled = pc.call_function("coalesce", [table.column(index), typed])
        table = table.set_column(index, table.schema.field(index), filled)
    return table


def fill_irt(table: pa.Table) -> pa.Table:
    """Fill a null ``irt`` from the ``irt`` or ``empirical_irt`` entry of ``cv_params``."""
    index = table.schema.get_field_index("irt")
    if index < 0 or "cv_params" not in table.schema.names or not table.num_rows:
        return table
    params = pa.concat_arrays(table.column("cv_params").chunks)
    entries = params.flatten()
    if len(entries) == 0:
        return table
    chosen = pc.call_function("is_in", [entries.field("cv_name")], pc.SetLookupOptions(pa.array(_IRT_PARAMS)))
    rows = np.flatnonzero(chosen.fill_null(False).to_numpy(zero_copy_only=False))
    if not len(rows) > 0:
        return table
    values = np.full(len(params), np.nan)
    parents = pc.call_function("list_parent_indices", [params]).to_numpy()
    numbers = pc.call_function("cast", [entries.field("cv_value").take(pa.array(rows))], pc.CastOptions(pa.float64()))
    values[parents[rows]] = numbers.to_numpy(zero_copy_only=False)
    typed = pa.array(values, mask=np.isnan(values), type=table.schema.field(index).type)
    return table.set_column(index, table.schema.field(index), pc.call_function("coalesce", [table.column(index), typed]))


def _members(entry) -> tuple[str, ...]:
    """The accessions of a ``pg_accessions`` value: strings on pg, structs on feature."""
    return tuple(member["accession"] if isinstance(member, dict) else member for member in entry)


def fill_group_ids(table: pa.Table) -> pa.Table:
    """Fill ``pg_group_id`` with the hash of the group's sorted membership; null without a membership."""
    index = table.schema.get_field_index("pg_group_id")
    if index < 0 or "pg_accessions" not in table.schema.names:
        return table
    ids: dict[tuple[str, ...], int] = {}
    values = []
    for entry in table.column("pg_accessions").to_pylist():
        members = _members(entry) if entry else ()
        if members and members not in ids:
            ids[members] = derive_id([sorted(members)])
        values.append(ids[members] if members else None)
    return table.set_column(index, table.schema.field(index), pa.array(values, type=pa.int64()))
