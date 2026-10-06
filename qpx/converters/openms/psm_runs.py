"""Recover native OpenMS PSM runs from exact companion spectrum evidence."""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pyarrow as pa

from qpx.converters.openms_consensus.feature_adapter import feature_map_info, to_proforma
from qpx.converters.openms_consensus.psm_adapter import _cf_element_runs, _run_resolver
from qpx.converters.openms_consensus.streaming import StreamingConsensusMap
from qpx.core.scan import scan_from_native_id

_SIGNATURE_COLUMNS = ("peptidoform", "charge", "scan", "rt", "observed_mz")


def _spectrum_key(row: dict, position: int) -> tuple:
    """Require complete numerical evidence before matching a source spectrum."""
    rt, mz = row["rt"], row["observed_mz"]
    valid_numbers = all(value is not None and math.isfinite(value) for value in (rt, mz))
    if not row["scan"] or not valid_numbers or mz <= 0:
        raise ValueError(f"Cannot recover PSM run at row {position}: incomplete numerical spectrum signature")
    return row["peptidoform"], row["charge"], tuple(row["scan"]), rt, mz


@dataclass
class PsmRunResolver:
    """Match complete spectrum evidence without trusting the native run label."""

    _references: list[tuple]
    _indexes: dict[tuple, dict] = field(default_factory=dict, init=False, repr=False)

    def _typed_index(self, rt_type: pa.DataType, mz_type: pa.DataType) -> dict:
        """Compare XML measurements at the native Parquet column precision."""
        types = (rt_type, mz_type)
        if types not in self._indexes:
            retention_times = pa.array([row[1] for row in self._references], type=pa.float64()).cast(rt_type)
            masses = pa.array([row[2] for row in self._references], type=pa.float64()).cast(mz_type)
            index: dict[tuple, set[str | None]] = defaultdict(set)
            for row, rt, mz in zip(self._references, retention_times.to_pylist(), masses.to_pylist()):
                index[(*row[0], rt, mz)].add(row[3])
            self._indexes[types] = index
        return self._indexes[types]

    def apply(self, table: pa.Table) -> pa.Table:
        """Replace only run_file_name; reject absent or ambiguous XML matches."""
        required = (*_SIGNATURE_COLUMNS, "run_file_name")
        missing = [column for column in required if column not in table.column_names]
        if missing:
            raise ValueError(f"Cannot recover PSM runs: missing signature columns {missing}")
        types = tuple(table.schema.field(column).type for column in ("rt", "observed_mz"))
        if any(dtype not in (pa.float32(), pa.float64()) for dtype in types):
            raise ValueError("Cannot recover PSM runs: rt and observed_mz must be float32 or float64")
        index = self._typed_index(*types)
        restored = []
        for position, row in enumerate(table.select(_SIGNATURE_COLUMNS).to_pylist()):
            key = _spectrum_key(row, position)
            candidates = index.get(key, set())
            if not candidates:
                raise ValueError(f"Cannot recover PSM run at row {position}: no exact companion XML match")
            if len(candidates) != 1 or None in candidates:
                raise ValueError(f"Cannot recover PSM run at row {position}: ambiguous or unresolved companion XML run")
            restored.append(next(iter(candidates)))
        run_field = table.schema.field("run_file_name")
        return table.set_column(table.schema.get_field_index(run_field.name), run_field, pa.array(restored, type=run_field.type))


def _pid_references(pid, run: str | None) -> list[tuple]:
    """Keep every numeric native-ID component and each hit's peptidoform."""
    reference = pid.getSpectrumReference()
    if not reference and pid.metaValueExists("spectrum_reference"):
        reference = pid.getMetaValue("spectrum_reference")
    scan = scan_from_native_id(str(reference or ""))
    if not scan:
        return []
    return [
        ((to_proforma(hit.getSequence()), int(hit.getCharge()), tuple(scan)), pid.getRT(), pid.getMZ(), run)
        for hit in pid.getHits()
    ]


def build_psm_run_lookup(consensusxml_path: str | Path) -> PsmRunResolver | None:
    """Read PID evidence once; header-only companions retain label-only usage."""
    consensus = StreamingConsensusMap(str(consensusxml_path))
    resolve_run = _run_resolver(consensus)
    map_info = feature_map_info(consensus)
    references = []
    has_identifications = False
    for kind, item in consensus.iter_all():
        identifications = item.getPeptideIdentifications() if kind == "element" else [item]
        runs = _cf_element_runs(item, map_info) if kind == "element" else None
        for pid in identifications:
            has_identifications = True
            references.extend(_pid_references(pid, resolve_run(pid, runs)))
    return PsmRunResolver(references) if has_identifications else None
