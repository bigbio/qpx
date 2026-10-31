"""Parquet encoding options shared by the QPX writers."""

from __future__ import annotations

import pyarrow as pa

# High-entropy float leaves that compress poorly under the default
# dictionary/PLAIN encoding (≈1.1x with ZSTD on continuous values). Parquet's
# BYTE_STREAM_SPLIT encoding reorders the bytes of IEEE floats so that ZSTD can
# exploit the per-byte regularity, measuring 15-28% smaller per column on real
# data. Paths use pyarrow's dotted leaf convention (``col``, ``a.list.element``,
# nested struct fields joined by ``.``). Only leaves present in a given schema
# are touched, so this list can name columns from any QPX view.
#
# NOTE: this is deliberately NOT applied to low-cardinality / dictionary-friendly
# float columns (e.g. ``additional_scores.score_value``,
# ``posterior_error_probability``), where BYTE_STREAM_SPLIT measured *larger*
# than the default dictionary encoding. Those keep their default encoding.
#
# PEP was a candidate for float32 + BYTE_STREAM_SPLIT (audit projected ~52%),
# but empirical validation over the real qpx corpus found genuine PEP values as
# small as ~1e-300..1e-318 (DIA/quantms feature+psm files): these underflow to
# zero in float32 and would corrupt any downstream -log10(PEP)/q-value use. The
# float32 narrowing was therefore held (PEP stays float64), so it is not split.
#
# The intensity leaves (``intensities.list.element.intensity`` and the nested
# ``additional_intensities…intensity_value``) are intentionally NOT split.
# BYTE_STREAM_SPLIT helps continuous LFQ intensities (~-20%) but slightly *hurts*
# TMT reporter-ion columns (+4-7%), which dominate TMT feature files (~64% of
# the file). Rather than probe the data at write time, we simply omit them here:
# they fall back to dictionary encoding, which is safe everywhere (never worse
# than PLAIN) and avoids the TMT regression. RT/mz leaves (a proven win) stay.
#
# The mz view's peak arrays are split as well: on Thermo DIA, timsTOF DDA-PASEF and
# diaPASEF spectra, m/z measured ~36% smaller than dictionary encoding; intensity
# 10% smaller (Thermo, DDA-PASEF) or 2% larger (diaPASEF), ~25% smaller files overall.
_BYTE_STREAM_SPLIT_LEAVES: tuple[str, ...] = (
    "predicted_rt",
    "rt",
    "rt_start",
    "rt_stop",
    "calculated_mz",
    "observed_mz",
    "mz.list.element",
    "intensity.list.element",
)

# ZSTD level 9 (vs the pyarrow default of 3) buys a few extra percent at modest
# write-cost; only codecs that accept a level get one.
_COMPRESSION_LEVELS: dict[str, int] = {"zstd": 9, "gzip": 9}


def _schema_leaf_paths(schema: pa.Schema) -> list[str]:
    """Return every leaf column path in pyarrow's dotted convention."""

    def _walk(field: pa.Field, prefix: str) -> list[str]:
        dtype = field.type
        if pa.types.is_struct(dtype):
            out: list[str] = []
            for i in range(dtype.num_fields):
                child = dtype.field(i)
                out.extend(_walk(child, f"{prefix}.{child.name}"))
            return out
        if pa.types.is_list(dtype) or pa.types.is_large_list(dtype):
            return _walk(dtype.value_field, f"{prefix}.list.element")
        return [prefix]

    leaves: list[str] = []
    for top in schema:
        leaves.extend(_walk(top, top.name))
    return leaves


def parquet_write_options(schema: pa.Schema, compression: str) -> dict:
    """Build QPX ParquetWriter options for *schema* and *compression*."""
    options: dict = {"compression": compression}
    level = _COMPRESSION_LEVELS.get(str(compression).lower())
    if level is not None:
        options["compression_level"] = level

    leaves = _schema_leaf_paths(schema)
    bss = {leaf: "BYTE_STREAM_SPLIT" for leaf in _BYTE_STREAM_SPLIT_LEAVES if leaf in leaves}
    # Always stamp the explicit v2.6 format and dictionary-encode non-BSS leaves so
    # every artifact (feature/psm and pg alike) uses the same encoding, rather than
    # falling back to pyarrow defaults when no byte-stream-split leaf is present.
    options["version"] = "2.6"
    options["use_dictionary"] = [leaf for leaf in leaves if leaf not in bss]
    if bss:
        options["column_encoding"] = bss
    return options
