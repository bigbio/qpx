"""Measured values are rounded to the mantissa bits their measurement supports."""

import numpy as np
import pyarrow as pa
import pytest

from qpx.transforms.mantissa import round_column, rounding_plan

# The step between floats of 12 mantissa bits at 1.0.
STEP = 2.0**-12


def _rounded(values, kind, path="x", keep=12):
    return round_column(pa.chunked_array([pa.array(values, kind)]), path, keep).to_pylist()


def test_values_round_to_nearest_with_ties_to_even():
    halfway_to_odd, halfway_to_even, above_half = 1 + STEP / 2, 1 + 1.5 * STEP, 1 + STEP / 2 + 2.0**-20

    rounded = _rounded([halfway_to_odd, halfway_to_even, above_half, -above_half], pa.float32())

    assert rounded == [1.0, 1 + 2 * STEP, 1 + STEP, -(1 + STEP)]


@pytest.mark.parametrize("keep", [12, 16])
def test_the_relative_error_is_at_most_half_a_step(keep):
    values = np.random.default_rng(0).lognormal(10, 3, 10_000)

    rounded = np.array(_rounded(values, pa.float64(), keep=keep))

    assert np.abs(rounded / values - 1).max() <= 2.0 ** -(keep + 1)


def test_nulls_and_values_that_are_not_finite_keep_their_bits():
    rounded = _rounded([None, float("inf"), float("-inf"), float("nan"), -0.0], pa.float64())

    assert rounded[:3] == [None, float("inf"), float("-inf")]
    assert np.isnan(rounded[3]) and np.signbit(rounded[4])


def test_a_struct_path_rounds_its_floats_in_lists_and_slices():
    kind = pa.list_(pa.struct([("label", pa.string()), ("intensity", pa.float32())]))
    rows = pa.array([[{"label": "a", "intensity": 1 + STEP / 2 + 2.0**-20}], None, [{"label": "b", "intensity": None}]], kind)
    column = pa.chunked_array([rows.slice(1, 2), rows.slice(0, 1)])

    rounded = round_column(column, "intensities.intensity", 12)

    assert rounded.to_pylist() == [None, [{"label": "b", "intensity": None}], [{"label": "a", "intensity": 1 + STEP}]]
    assert round_column(column, "intensities.label", 12).equals(column)


def test_the_plan_leaves_out_the_identity_composite_and_columns_without_floats():
    schema = pa.schema([("rt", pa.float32()), ("consensus_rt", pa.float32()), ("charge", pa.int16())])
    schema = schema.with_metadata({b"identity_composite": b"peptidoform,consensus_rt"})

    plan = rounding_plan({"rt": 16, "consensus_rt": 16, "charge": 12, "ccs": 16}, schema)

    assert plan == {"rt": 16}
