"""YAML schema loader — reads spec YAML files and resolves them to PyArrow types.

This is the single engine that converts the canonical YAML definitions in
``qpx/core/data/schemas/`` into :class:`ViewSchema` instances used by the rest
of the library.
"""

from __future__ import annotations

import re
from pathlib import Path

import pyarrow as pa
import yaml

from qpx.core.data.schema import FieldDef, ViewSchema

SPECS_DIR = Path(__file__).parent / "schemas"

# Mapping from YAML type names to PyArrow primitive types
PRIMITIVE_MAP = {
    "string": pa.string(),
    "int16": pa.int16(),
    "int32": pa.int32(),
    "int64": pa.int64(),
    "float32": pa.float32(),
    "float64": pa.float64(),
    "bool": pa.bool_(),
    "binary": pa.binary(),
}

# Module-level cache for custom types loaded from types.yaml
_custom_types: dict[str, pa.DataType] | None = None


# --- Type resolution ---


def _load_custom_types() -> dict[str, pa.DataType]:
    """Load and cache all custom struct types from types.yaml."""
    global _custom_types
    if _custom_types is not None:
        return _custom_types

    with open(SPECS_DIR / "types.yaml", encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    _custom_types = {}
    for type_name, type_def in raw.items():
        _custom_types[type_name] = _build_struct_type(type_def["fields"], _custom_types)
    return _custom_types


def _build_struct_type(fields_def: dict, custom_types: dict[str, pa.DataType]) -> pa.StructType:
    """Build a PyArrow struct type from a YAML fields definition."""
    arrow_fields = []
    for fname, fdef in fields_def.items():
        nullable = fdef.get("nullable", False)
        arrow_type = _resolve_type(fdef["type"], custom_types)
        arrow_fields.append(pa.field(fname, arrow_type, nullable=nullable))
    return pa.struct(arrow_fields)


def _resolve_type(type_str: str, custom_types: dict[str, pa.DataType]) -> pa.DataType:
    """Resolve a YAML type string to a PyArrow DataType.

    Handles primitives, ``list<X>``, ``map<K, V>``, and custom type references.
    """
    # Primitives
    if type_str in PRIMITIVE_MAP:
        return PRIMITIVE_MAP[type_str]

    # list<inner>
    m = re.match(r"^list<(.+)>$", type_str)
    if m:
        inner = _resolve_type(m.group(1), custom_types)
        return pa.list_(inner)

    # map<key, value>  (value may itself contain angle brackets)
    if type_str.startswith("map<") and type_str.endswith(">"):
        inner = type_str[4:-1]
        key_str, val_str = _split_generic_args(inner)
        return pa.map_(
            _resolve_type(key_str.strip(), custom_types),
            _resolve_type(val_str.strip(), custom_types),
        )

    # Custom type reference (defined in types.yaml)
    if type_str in custom_types:
        return custom_types[type_str]

    raise ValueError(f"Unknown type in YAML schema: {type_str!r}")


def _split_generic_args(s: str) -> tuple[str, str]:
    """Split ``K, V`` where V may contain nested angle brackets."""
    depth = 0
    for i, c in enumerate(s):
        if c == "<":
            depth += 1
        elif c == ">":
            depth -= 1
        elif c == "," and depth == 0:
            return s[:i], s[i + 1 :]
    raise ValueError(f"Cannot split generic args: {s!r}")


# --- Schema loading ---


def load_schema(name: str) -> ViewSchema:
    """Load a schema YAML spec and return a :class:`ViewSchema` instance.

    Parameters
    ----------
    name : str
        Schema name matching the YAML filename (without extension),
        e.g. ``"feature"``, ``"psm"``, ``"pg"``.
    """
    custom_types = _load_custom_types()

    with open(SPECS_DIR / f"{name}.yaml", encoding="utf-8") as f:
        spec = yaml.safe_load(f)

    fields: dict[str, FieldDef] = {}
    for fname, fdef in spec["fields"].items():
        arrow_type = _resolve_type(fdef["type"], custom_types)
        required = fdef.get("required", False)
        optional = fdef.get("optional", False)
        nullable = fdef.get("nullable", not required)
        doc = fdef.get("doc", "")
        fields[fname] = FieldDef(fname, arrow_type, nullable=nullable, optional=optional, doc=doc)

    return ViewSchema(
        view_name=spec["name"],
        file_type=spec["file_type"],
        primary_key=spec["primary_key"],
        fields=fields,
        doc=spec.get("doc", ""),
        extra_columns=spec.get("extra_columns", False),
        identity_composite=spec.get("identity_composite"),
    )
