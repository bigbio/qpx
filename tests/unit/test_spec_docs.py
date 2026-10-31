"""The field tables of the spec pages agree with the YAML schemas."""

import re
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
SCHEMAS = ROOT / "qpx" / "core" / "data" / "schemas"
SPEC = ROOT / "docs" / "spec"
VIEWS = [
    "feature",
    "psm",
    "pg",
    "dataset",
    "run",
    "sample",
    "mz",
    "pepmap",
    "site",
    "proteoform",
    "provenance",
    "protein",
    "fragment",
]
TYPES = yaml.safe_load((SCHEMAS / "types.yaml").read_text(encoding="utf-8"))
TYPES = TYPES.get("types", TYPES)
SPELLINGS = {"float": "float32", "int": "int32", "boolean": "bool"}
_CELL_SEPARATOR = re.compile(r"(?<!\\)\|")


def _schema_fields(view):
    return yaml.safe_load((SCHEMAS / f"{view}.yaml").read_text(encoding="utf-8"))["fields"]


def _documented_fields(view):
    """Field name -> Type cell of every table with a Field and a Type column."""
    fields = {}
    header = None
    for line in (SPEC / f"{view}.md").read_text(encoding="utf-8").splitlines():
        if not line.startswith("|"):
            header = None
            continue
        cells = [cell.strip() for cell in _CELL_SEPARATOR.split(line)[1:-1]]
        if cells[0].lower() == "field":
            header = [cell.lower() for cell in cells]
        elif header and "type" in header and len(cells) == len(header):
            for name in re.findall(r"`([a-z_0-9]+)`", cells[0]):
                fields[name] = cells[header.index("type")]
    return fields


def _struct_members():
    members = set()
    for definition in TYPES.values():
        if isinstance(definition, dict):
            members.update(definition.get("fields", definition))
    return members


def _normal(type_text):
    text = type_text.replace("`", "").replace("\\", "").replace("&lt;", "<").replace("&gt;", ">").lower()
    text = re.sub(r",\s*null$", "", text).replace(" ", "")
    text = re.sub(r"(?:array|list)\[(.*)\]$", r"list<\1>", text)
    head = re.fullmatch(r"list<(\w+)>", text)
    if head:
        return f"list<{SPELLINGS.get(head.group(1), head.group(1))}>"
    return SPELLINGS.get(text, text)


@pytest.mark.parametrize("view", VIEWS)
def test_every_schema_field_is_documented(view):
    assert sorted(set(_schema_fields(view)) - set(_documented_fields(view))) == []


@pytest.mark.parametrize("view", VIEWS)
def test_every_documented_field_exists(view):
    extra = set(_documented_fields(view)) - set(_schema_fields(view)) - _struct_members()
    assert sorted(extra) == []


@pytest.mark.parametrize("view", VIEWS)
def test_documented_types_match_the_schema(view):
    schema = _schema_fields(view)
    wrong = []
    for name, type_text in _documented_fields(view).items():
        if name not in schema:
            continue
        expected = str(schema[name]["type"]).replace(" ", "")
        if expected.removeprefix("list<").removesuffix(">") in TYPES or expected in TYPES:
            # A nested type: the page may say list<struct> or name the struct, never a scalar list
            if not re.search(r"struct|" + re.escape(expected.removeprefix("list<").removesuffix(">")), type_text.lower()):
                wrong.append((name, type_text, expected))
        elif "," not in type_text.replace(", null", "") and _normal(type_text) != expected:
            wrong.append((name, type_text, expected))
    assert wrong == []


@pytest.mark.parametrize("view", VIEWS)
def test_only_nullable_fields_are_documented_as_null(view):
    schema = _schema_fields(view)
    wrong = [
        name
        for name, type_text in _documented_fields(view).items()
        if name in schema and "null" in type_text.lower() and schema[name].get("nullable") is False
    ]
    assert wrong == []
