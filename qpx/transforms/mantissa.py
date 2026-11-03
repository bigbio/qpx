"""Round measured values to the mantissa bits their measurement supports."""

from __future__ import annotations

import json

import numpy as np
import pyarrow as pa

MANTISSA_BITS_FOOTER_KEY = b"qpx.mantissa_bits"

# Measured values keep 12 (intensities, mass errors) or 16 (retention times, ion
# mobility) mantissa bits: a relative error of at most 2^-13 or 2^-17, without bias.
# On MSV000085836 the log2 ratio between repeated PSMs of one precursor varies by
# 0.22 (MAD), 3,000 times the error 12 bits add to it (RMS 7.1e-5). Scores,
# q-values, PEP and m/z keep every bit, so the rows that pass a threshold stay the
# same. A path names a column, then struct fields; list levels are left out, and a
# path that ends above the floats covers all of them. The cv_params paths are the
# iRT and ion mobility terms of DIA-NN and Spectronaut, which are floats once
# relayout stores the column as a struct (qpx.core.struct_lists); in the list
# their values are strings and stay as they are.
MANTISSA_BITS = {
    "feature": {
        "intensities.intensity": 12,
        "additional_intensities": 12,
        "mass_error_ppm": 12,
        "rt": 16,
        "rt_start": 16,
        "rt_stop": 16,
        "rt_aligned": 16,
        "consensus_rt": 16,
        "predicted_rt": 16,
        "irt": 16,
        "channel_details.rt": 16,
        "ion_mobility": 16,
        "ion_mobility_start": 16,
        "ion_mobility_stop": 16,
        "ccs": 16,
        "cv_params.irt": 16,
        "cv_params.empirical_irt": 16,
        "cv_params.predicted_irt": 16,
        "cv_params.iim": 16,
        "cv_params.predicted_im": 16,
        "cv_params.predicted_iim": 16,
    },
    "psm": {
        "intensity_array": 12,
        "mass_error_ppm": 12,
        "rt": 16,
        "predicted_rt": 16,
        "ion_mobility": 16,
        "ion_mobility_array": 16,
        "ccs": 16,
    },
    "pg": {"intensity": 12, "additional_intensities": 12},
    "mz": {
        "intensity": 12,
        "precursors.selected_ion_intensity": 12,
        "total_ion_current": 12,
        "scan_start_time": 16,
        "inverse_ion_mobility": 16,
        "drift_time": 16,
        "ion_mobility_array": 16,
    },
}

_MANTISSA = {32: (23, np.uint32), 64: (52, np.uint64)}


def _holds_floats(kind: pa.DataType, parts: list[str]) -> bool:
    """Whether the struct path *parts* of a value of *kind* leads to floats."""
    if pa.types.is_floating(kind):
        return not parts
    if pa.types.is_list(kind):
        return _holds_floats(kind.value_type, parts)
    if pa.types.is_struct(kind):
        return any(_holds_floats(field.type, parts[1:]) for field in kind if not parts or field.name == parts[0])
    return False


def rounding_plan(bits: dict[str, int], schema: pa.Schema) -> dict[str, int]:
    """The paths of *bits* that *schema* holds floats at, leaving out the columns of its identity composite.

    A rounded composite value would no longer name its row.
    """
    composite = set((schema.metadata or {}).get(b"identity_composite", b"").decode().split(","))
    plan = {}
    for path, keep in bits.items():
        column, *parts = path.split(".")
        if column in schema.names and column not in composite and _holds_floats(schema.field(column).type, parts):
            plan[path] = keep
    return plan


def mantissa_footer(plan: dict[str, int]) -> bytes:
    """The ``qpx.mantissa_bits`` footer value: the mantissa bits each rounded path keeps."""
    return json.dumps(plan, separators=(",", ":")).encode()


def _bitround(array: pa.Array, keep: int) -> pa.Array:
    """*array* rounded to *keep* mantissa bits, to nearest with ties to even, as numcodecs' BitRound does."""
    mantissa, unsigned = _MANTISSA[array.type.bit_width]
    drop = mantissa - keep
    values = array.to_numpy(zero_copy_only=False)
    integers = values.view(unsigned)
    half = unsigned((1 << (drop - 1)) - 1)
    rounded = ((integers + ((integers >> unsigned(drop)) & unsigned(1)) + half) & ~unsigned((1 << drop) - 1)).view(values.dtype)
    # A NaN or an infinity keeps its bits: rounding can carry a NaN's mantissa into the sign.
    rounded = np.where(np.isfinite(values), rounded, values)
    mask = array.is_null().to_numpy(zero_copy_only=False) if array.null_count else None
    return pa.array(rounded, type=array.type, mask=mask)


def _round(array: pa.Array, parts: list[str], keep: int) -> pa.Array:
    """*array* with the floats at the struct path *parts*, or all of them, rounded."""
    kind = array.type
    if pa.types.is_floating(kind):
        return _bitround(array, keep)
    if pa.types.is_list(kind):
        # ListArray.from_arrays takes no null mask with the offsets of a slice.
        array = pa.concat_arrays([array]) if array.offset else array
        mask = array.is_null() if array.null_count else None
        return pa.ListArray.from_arrays(array.offsets, _round(array.values, parts, keep), type=kind, mask=mask)
    if pa.types.is_struct(kind):
        children = [
            _round(array.field(index), parts[1:], keep) if not parts or field.name == parts[0] else array.field(index)
            for index, field in enumerate(kind)
        ]
        mask = array.is_null() if array.null_count else None
        return pa.StructArray.from_arrays(children, fields=list(kind), mask=mask)
    return array


def round_column(column: pa.ChunkedArray, path: str, keep: int) -> pa.ChunkedArray:
    """*column* with the floats at *path* (its name, then struct fields) rounded to *keep* mantissa bits."""
    parts = path.split(".")[1:]
    return pa.chunked_array([_round(chunk, parts, keep) for chunk in column.chunks], type=column.type)


def round_table(table: pa.Table, plan: dict[str, int]) -> pa.Table:
    """*table* with the floats at each path of *plan* rounded to its mantissa bits."""
    for path, keep in plan.items():
        index = table.schema.get_field_index(path.split(".")[0])
        table = table.set_column(index, table.field(index), round_column(table.column(index), path, keep))
    return table
