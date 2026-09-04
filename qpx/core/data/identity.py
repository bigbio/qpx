"""Opaque identity hashing for QPX view primary keys.

The ``feature``, ``psm`` and ``pg`` views each carry a single mandatory
identity column (``feature_id`` / ``psm_id`` / ``pg_id``) that is the primary
key of the view. When no producer ID is supplied, the *writer* derives an
**opaque** signed 64-bit integer deterministically from a footer-declared
composite of existing columns (the view's ``identity_composite``). It is not
meant to be parsed or reversed — it is only an identity/equality token.

Because the id is a fixed-width hash of the composite, distinct composites can
in principle collide onto the same id. That is not silently tolerated: the id
is the primary key, and the primary-key **uniqueness** validation
(``duplicate_pk``) catches any collision, so a real clash surfaces as a
validation issue rather than a silent data loss.
"""

from __future__ import annotations

import hashlib
import json


def _normalize(value):
    """Recursively normalize a composite value into a JSON-serializable form."""
    if isinstance(value, (list, tuple)):
        return [_normalize(v) for v in value]
    return value


def canonical(values: list, *, unordered_list_indices: tuple[int, ...] = ()) -> bytes:
    """Encode composite *values* into a canonical byte string.

    Uses canonical JSON (fixed separators, sorted keys). JSON quotes and escapes
    strings, so — unlike a plain delimiter join — distinct composites can never
    alias to the same bytes even when list elements or fields contain the
    delimiters themselves (e.g. run names with commas/brackets): ``["a,b", "c"]``
    and ``["a", "b,c"]`` encode differently. ``None`` becomes JSON ``null``.
    List order and duplicate elements are preserved unless its top-level
    position is explicitly listed in *unordered_list_indices*. QPX uses that
    exception for set-valued fields such as ``grouped_runs`` and
    ``pg_accessions``; ordered identifiers such as a multi-component ``scan``
    keep their component order.
    """
    normalized = [_normalize(value) for value in values]
    for index in unordered_list_indices:
        value = normalized[index]
        if isinstance(value, list):
            by_encoding = {json.dumps(item, separators=(",", ":"), sort_keys=True, ensure_ascii=True): item for item in value}
            normalized[index] = [by_encoding[key] for key in sorted(by_encoding)]
    return json.dumps(normalized, separators=(",", ":"), sort_keys=True, ensure_ascii=True).encode("utf-8")


def derive_id(values: list, *, unordered_list_indices: tuple[int, ...] = ()) -> int:
    """Derive a signed 64-bit opaque identity from composite *values*.

    Uses BLAKE2b truncated to 8 bytes. The identity is opaque; uniqueness is
    guaranteed not by construction but by the primary-key uniqueness validation.
    """
    encoded = canonical(values, unordered_list_indices=unordered_list_indices)
    return int.from_bytes(hashlib.blake2b(encoded, digest_size=8).digest(), "big", signed=True)


# Sentinel for values that cannot be used as a cache key (e.g. dicts).
_UNHASHABLE = object()


def _cache_key(value):
    """A hashable stand-in for *value*, or ``_UNHASHABLE`` when there is none.

    The type is part of the key. Python hashes ``1``, ``1.0`` and ``True`` as the
    same dict key, but they encode to ``1``, ``1.0`` and ``true`` — caching on the
    bare value would hand one row another row's encoding and silently change its id.
    """
    if isinstance(value, (str, int, float, bool, type(None))):
        return (type(value), value)
    if isinstance(value, (list, tuple)):
        parts = []
        for item in value:
            key = _cache_key(item)
            if key is _UNHASHABLE:
                return _UNHASHABLE
            parts.append(key)
        return (type(value), tuple(parts))
    return _UNHASHABLE


def _encode_element(value, *, unordered: bool) -> str:
    """Encode one composite element exactly as :func:`canonical` would.

    ``json.dumps`` with ``separators=(",", ":")`` writes a list as
    ``[<e0>,<e1>,...]`` with no padding, so encoding each element independently
    and joining with ``,`` reproduces the whole-list encoding byte for byte.
    """
    normalized = _normalize(value)
    if unordered and isinstance(normalized, list):
        by_encoding = {json.dumps(item, separators=(",", ":"), sort_keys=True, ensure_ascii=True): item for item in normalized}
        normalized = [by_encoding[key] for key in sorted(by_encoding)]
    return json.dumps(normalized, separators=(",", ":"), sort_keys=True, ensure_ascii=True)


def derive_ids(columns: list[list], *, unordered_list_indices: tuple[int, ...] = ()) -> list[int]:
    """Derive identities for many rows at once, column-wise.

    Equivalent to calling :func:`derive_id` per row — verified byte for byte by
    the test suite — but encodes each *distinct* value in a column only once.

    The per-row path re-encoded every value on every row, which dominated large
    conversions: a 231M-row DIA-NN report has only ~5,800 distinct run names and
    a handful of charges, yet each row paid a full ``json.dumps`` for both. On a
    real profile this was 48% of the feature conversion (bigbio/qpx#291).

    *columns* is one list per composite field, all the same length.
    """
    if not columns:
        return []
    n = len(columns[0])
    if any(len(col) != n for col in columns):
        raise ValueError("all composite columns must have the same length")

    encoded_columns: list[list[str]] = []
    for index, column in enumerate(columns):
        unordered = index in unordered_list_indices
        cache: dict = {}
        encoded: list[str] = []
        for value in column:
            key = _cache_key(value)
            if key is _UNHASHABLE:
                encoded.append(_encode_element(value, unordered=unordered))
                continue
            fragment = cache.get(key)
            if fragment is None:
                fragment = _encode_element(value, unordered=unordered)
                cache[key] = fragment
            encoded.append(fragment)
        encoded_columns.append(encoded)

    blake2b = hashlib.blake2b
    from_bytes = int.from_bytes
    return [
        from_bytes(blake2b(("[" + ",".join(row) + "]").encode("utf-8"), digest_size=8).digest(), "big", signed=True)
        for row in zip(*encoded_columns)
    ]
