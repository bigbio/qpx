"""Evidence-table fields recording how MaxQuant identified and isolated a feature, and msms.txt peak annotations."""

import re

from qpx.converters.utils import safe_float

_SECONDS_PER_MINUTE = 60.0
# MaxQuant's ``Type`` column. A MULTI-MATCH-MSMS row was fragmented too, but its
# identification still came from another run.
_IDENTIFICATION_TYPES = {
    "MSMS": "direct",
    "MULTI-MSMS": "direct",
    "ISO-MSMS": "direct",
    "MULTI-SECPEP": "secondary",
    "MULTI-MATCH": "transferred",
    "MULTI-MATCH-MSMS": "transferred",
}


def identification_fields(row: dict, evidence_type: str) -> dict:
    """The identification type, match metadata and precursor purity of one evidence row."""
    rt_delta = safe_float(row.get("Match time difference"))
    return {
        "identification_type": _IDENTIFICATION_TYPES.get(evidence_type),
        "match_score": safe_float(row.get("Match score")),
        "match_qvalue": safe_float(row.get("Match q-value")),
        "match_rt_delta": None if rt_delta is None else rt_delta * _SECONDS_PER_MINUTE,
        "match_mz_delta": safe_float(row.get("Match m/z difference")),
        # MaxQuant writes 0 when it did not compute the fraction: an identified precursor has some signal.
        "parent_ion_fraction": safe_float(row.get("PIF")) or None,
    }


_FRAGMENT_CHARGE = re.compile(r"\((\d+)\+\)")


def peak_annotations(matches, peak_count: int) -> tuple[list[str], list[int]] | tuple[None, None]:
    """Ion types and charges parallel to the peaks from a ``Matches`` cell: ``y3;b4(2+);y5-H2O``.

    Without one token per peak the two cannot be aligned, so there are none.
    """
    tokens = str(matches).split(";") if matches else []
    if not tokens or len(tokens) != peak_count:
        return None, None
    charges = [(found := _FRAGMENT_CHARGE.search(token)) and int(found.group(1)) or 1 for token in tokens]
    return [_FRAGMENT_CHARGE.sub("", token) for token in tokens], charges
