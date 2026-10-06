"""Encode numeric native spectrum identifiers without discarding components."""

import re

_NATIVE_COMPONENT = re.compile(r"([A-Za-z][A-Za-z0-9_]*)=(\d+)(?=\s|$)")
_NATIVE_FORMAT_KEYS = {
    ("controllertype", "controllernumber", "scan"),
    ("frame", "scan"),
    ("frame", "scan", "precursor"),
    ("frame", "windowgroup", "scan"),
    ("merged", "frame", "scanstart", "scanend"),
    ("function", "process", "scan"),
    ("sample", "period", "cycle", "experiment"),
}


def scan_format_from_native_id(native_id: str) -> str | None:
    """Classify explicit native IDs; ambiguous ordinals and unknown formats stay unset."""
    components = _NATIVE_COMPONENT.findall(native_id)
    if " ".join(f"{key}={value}" for key, value in components) != " ".join(native_id.split()):
        return None
    keys = tuple(key.lower() for key, _ in components)
    if keys == ("scan",):
        return "scan"
    if keys == ("index",):
        return "index"
    if keys not in _NATIVE_FORMAT_KEYS:
        return None
    if keys == ("controllertype", "controllernumber", "scan") and [int(value) for _, value in components[:2]] == [0, 1]:
        return "scan"
    return "nativeId"


def scan_from_native_id(native_id: str) -> list[int]:
    """Keep numeric components in native order, including repeated values.

    Default Thermo controller components are omitted according to the QPX scan
    convention. An identifier without numeric components returns an empty list.
    """
    components = [(key.lower(), int(value)) for key, value in _NATIVE_COMPONENT.findall(native_id)]
    fields = dict(components)
    if len(components) == 3 and fields.get("controllertype") == 0 and fields.get("controllernumber") == 1 and "scan" in fields:
        return [fields["scan"]]
    return [value for _, value in components]
