"""Unit tests for the lossless spectra codec (``qpx.core.spectra_codec``).

The codec is bit-exact: decoding ``encode_mz`` / ``encode_intensity`` must return
the original ``float32`` arrays unchanged. m/z is delta-encoded (requires a sorted
array) and intensity is bitshuffled; both are lossless.
"""

import numpy as np
import pytest

from qpx.core.spectra_codec import decode_intensity, decode_mz, encode_intensity, encode_mz


def _sorted_mz(n: int, seed: int = 0) -> np.ndarray:
    """A strictly-increasing ``float32`` m/z array of length *n*."""
    rng = np.random.default_rng(seed)
    return np.ascontiguousarray(np.cumsum(rng.uniform(0.01, 2.0, size=n)).astype(np.float32))


def _intensity(n: int, seed: int = 1) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.ascontiguousarray(rng.uniform(0.0, 1e6, size=n).astype(np.float32))


@pytest.mark.parametrize("n", [0, 1, 2, 3, 50, 1000])
def test_mz_roundtrip_bit_exact(n):
    mz = _sorted_mz(n)
    decoded = decode_mz(encode_mz(mz), n)
    assert decoded.dtype == np.float32
    np.testing.assert_array_equal(decoded, mz)


@pytest.mark.parametrize("n", [0, 1, 2, 3, 50, 1000])
def test_intensity_roundtrip_bit_exact(n):
    intensity = _intensity(n)
    decoded = decode_intensity(encode_intensity(intensity), n)
    assert decoded.dtype == np.float32
    np.testing.assert_array_equal(decoded, intensity)


def test_mz_single_peak_blob_is_just_first_value():
    mz = np.array([175.1190], dtype=np.float32)
    blob = encode_mz(mz)
    assert len(blob) == 4  # 4-byte first uint32, no varint deltas
    np.testing.assert_array_equal(decode_mz(blob, 1), mz)


def test_empty_arrays_roundtrip():
    assert encode_mz(np.array([], dtype=np.float32)) == b""
    assert encode_intensity(np.array([], dtype=np.float32)) == b""
    assert decode_mz(b"", 0).size == 0
    assert decode_intensity(b"", 0).size == 0


def test_non_contiguous_input():
    mz = _sorted_mz(10)[::2]  # non-contiguous slice
    np.testing.assert_array_equal(decode_mz(encode_mz(mz), len(mz)), mz)

    intensity = _intensity(10)[::2]
    np.testing.assert_array_equal(decode_intensity(encode_intensity(intensity), len(intensity)), intensity)


def test_mz_preserves_full_float32_range():
    """Large and small m/z values both survive the uint32 reinterpretation."""
    mz = np.array([50.0, 150.0, 500.0, 2000.0, 4000.0], dtype=np.float32)
    np.testing.assert_array_equal(decode_mz(encode_mz(mz), len(mz)), mz)


def test_intensity_handles_zeros_and_negatives():
    """Zeros and (unusual) non-positive intensities round-trip bit-exactly."""
    intensity = np.array([0.0, 1e-6, 12345.0, -1.0], dtype=np.float32)
    np.testing.assert_array_equal(decode_intensity(encode_intensity(intensity), len(intensity)), intensity)
