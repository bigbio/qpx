"""Lossless codec for fragment spectra (m/z + intensity arrays).

The storage design separates feature metadata (Parquet) from fragment spectra,
deduplicated by ``(run_file_name, scan)`` and stored in the ``mz`` view whose
``mz`` / ``intensity`` columns carry the structurally-transformed bytes
below. Parquet's own zstd compression then applies entropy coding on top.

Both encodings are **bit-exact**: decoding returns the original ``float32`` arrays.

m/z
    ``float32`` -> bit-reinterpret ``uint32`` (positive floats are monotonic under
    IEEE-754 bit order) -> per-spectrum delta -> little-endian varint. Stored as
    ``[4-byte first value][varint deltas...]``. Requires the array to be sorted.

intensity
    ``float32`` -> bit-reinterpret ``uint32`` -> 4-byte-plane transpose
    (bitshuffle). This groups the highly-redundant exponent/sign bytes together so
    the downstream compressor is effective. The 23-bit mantissa is near-random and
    is the lossless entropy floor.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "encode_mz",
    "decode_mz",
    "encode_intensity",
    "decode_intensity",
]


def _varint_encode(values: np.ndarray) -> bytes:
    """Encode non-negative integers as little-endian base-128 varints."""
    out = bytearray()
    for v in values:
        v = int(v)
        while v >= 0x80:
            out.append((v & 0x7F) | 0x80)
            v >>= 7
        out.append(v)
    return bytes(out)


def _varint_decode(data: bytes, n: int) -> np.ndarray:
    """Decode ``n`` little-endian base-128 varints."""
    out = np.zeros(n, dtype=np.uint64)
    j = 0
    val = 0
    shift = 0
    for b in data:
        val |= (b & 0x7F) << shift
        if not b & 0x80:
            out[j] = val
            j += 1
            val = 0
            shift = 0
            if j >= n:
                break
        else:
            shift += 7
    return out


def encode_mz(mz) -> bytes:
    """Losslessly encode a sorted ``float32`` m/z array.

    Parameters
    ----------
    mz:
        Array-like of ``float32`` m/z values, sorted ascending.

    Returns
    -------
    bytes
        ``[4-byte first uint32][varint deltas]``.
    """
    arr = np.ascontiguousarray(mz, dtype=np.float32)
    if arr.size == 0:
        return b""
    u32 = arr.view(np.uint32)
    first = int(u32[0])
    deltas = np.diff(u32.astype(np.uint64)) if arr.size > 1 else np.array([], dtype=np.uint64)
    return first.to_bytes(4, "little") + _varint_encode(deltas)


def decode_mz(blob: bytes, n_peaks: int) -> np.ndarray:
    """Decode ``blob`` produced by :func:`encode_mz` into ``float32`` m/z values."""
    if n_peaks == 0:
        return np.array([], dtype=np.float32)
    first = int.from_bytes(blob[:4], "little")
    out = np.zeros(n_peaks, dtype=np.uint32)
    out[0] = first
    if n_peaks > 1:
        deltas = _varint_decode(blob[4:], n_peaks - 1)
        out[1:] = first + np.cumsum(deltas, dtype=np.uint64).astype(np.uint32)
    return out.view(np.float32)


def encode_intensity(intensity) -> bytes:
    """Losslessly encode a ``float32`` intensity array via bitshuffle (byte-plane transpose)."""
    arr = np.ascontiguousarray(intensity, dtype=np.float32)
    if arr.size == 0:
        return b""
    planes = arr.view(np.uint32).view(np.uint8).reshape(-1, 4)
    return planes.T.reshape(-1).tobytes()


def decode_intensity(blob: bytes, n_peaks: int) -> np.ndarray:
    """Decode ``blob`` produced by :func:`encode_intensity` into ``float32`` intensities."""
    if n_peaks == 0:
        return np.array([], dtype=np.float32)
    planes = np.frombuffer(blob, dtype=np.uint8).reshape(4, n_peaks)
    return planes.T.reshape(-1).view(np.uint32).view(np.float32)
