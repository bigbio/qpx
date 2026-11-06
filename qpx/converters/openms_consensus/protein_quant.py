"""Protein abundances OpenMS records in the consensusXML header.

ProteomicsLFQ, IsobaricWorkflow and MS1LabeledWorkflow store their protein
quantification as ``indistinguishable_proteins_<N>_fraction_group_level_*`` user
parameters next to the group itself. pyopenms does not expose them, so both
readers take them from the header, which precedes the bulk of the file.
"""

from __future__ import annotations

import re
from typing import NamedTuple

from qpx.converters.openms_consensus.channels import channels_by_column, list_value, ratio_name
from qpx.converters.openms_consensus.streaming import _iterparse, _localname, _user_params

_GROUP_PARAM = re.compile(r"^indistinguishable_proteins_(\d+)$")
_PREFIX = "indistinguishable_proteins_{}_fraction_group_level_"
_HEADER_END = ("UnassignedPeptideIdentification", "mapList", "consensusElementList")


class ProteinQuant(NamedTuple):
    """Fraction-group level quantification of one protein group."""

    accessions: frozenset[str]
    abundances: dict[tuple[int, int], float]
    ratios: dict[tuple[int, int], dict[str, float]]


def _column(meta: dict[str, str], prefix: str, name: str) -> list[str]:
    return list_value(meta.get(f"{prefix}{name}", ""))


def _group_quant(number: str, accessions: frozenset[str], meta: dict[str, str]) -> ProteinQuant | None:
    prefix = _PREFIX.format(number)
    if f"{prefix}abundance" not in meta:
        return None
    abundances = {
        (int(group), int(label)): float(value)
        for group, label, value in zip(
            _column(meta, prefix, "fraction_group"), _column(meta, prefix, "label"), _column(meta, prefix, "abundance")
        )
    }
    ratios: dict[tuple[int, int], dict[str, float]] = {}
    keys = list(zip(_column(meta, prefix, "ratio_fraction_group"), _column(meta, prefix, "ratio_label")))
    for name in ("ratio", "ratio_normalized"):
        for (group, label), value in zip(keys, _column(meta, prefix, name)):
            ratios.setdefault((int(group), int(label)), {})[name] = float(value)
    return ProteinQuant(accessions, abundances, ratios)


def read_protein_quant(path: str) -> list[ProteinQuant]:
    """Fraction-group abundances and ratios per protein group; empty when the file has none."""
    protein_hits: dict[str, str] = {}
    quants: list[ProteinQuant] = []
    for event, el in _iterparse(path, events=("start", "end")):
        tag = _localname(el.tag)
        if event == "start" and tag in _HEADER_END:
            break
        if event == "end" and tag == "ProteinHit":
            protein_hits[el.attrib.get("id", "")] = el.attrib.get("accession", "")
            el.clear()
        elif event == "end" and tag == "ProteinIdentification":
            meta = _user_params(el)
            for name, value in meta.items():
                match = _GROUP_PARAM.match(name)
                if match:
                    accessions = frozenset(protein_hits[ref] for ref in value.split(",")[1:] if ref in protein_hits)
                    quant = _group_quant(match.group(1), accessions, meta)
                    if quant is not None:
                        quants.append(quant)
            el.clear()
    return quants


def fraction_group_units(cm, runs: dict[int, str]) -> dict[int, tuple[str, ...]]:
    """Fraction group -> its runs, from the map list; empty unless every map declares one."""
    headers = cm.getColumnHeaders()
    if not headers or not all(headers[idx].metaValueExists("fraction_group") for idx in headers):
        return {}
    units: dict[int, set[str]] = {}
    for idx in headers:
        units.setdefault(int(headers[idx].getMetaValue("fraction_group")), set()).add(runs[idx])
    return {group: tuple(sorted(members)) for group, members in units.items()}


def _reference_labels(quant: ProteinQuant, label_of) -> dict[int, str]:
    """Per fraction group, the reference label: the lowest channel with a ratio."""
    references: dict[int, int] = {}
    for group, channel in quant.ratios:
        references[group] = min(references.get(group, channel), channel)
    return {group: label_of[(group, channel)] for group, channel in references.items()}


def quant_by_unit(cm, runs: dict[int, str], labels: dict[int, str], path: str) -> dict[tuple, dict]:
    """``(group accessions, unit, label)`` -> ``{"abundance", "ratios"}`` of one input.

    A zero abundance is OpenMS' "not quantified" and stays absent.
    """
    units = fraction_group_units(cm, runs)
    if not units:
        return {}
    headers = cm.getColumnHeaders()
    channels = channels_by_column(cm)
    label_of = {(int(headers[idx].getMetaValue("fraction_group")), channels[idx]): labels[idx] for idx in headers}
    result: dict[tuple, dict] = {}
    for quant in read_protein_quant(path):
        references = _reference_labels(quant, label_of)
        for (group, channel), abundance in quant.abundances.items():
            label = label_of[(group, channel)]
            result[(quant.accessions, units[group], label)] = {"abundance": abundance or None, "ratios": {}}
        for (group, channel), values in quant.ratios.items():
            label = label_of[(group, channel)]
            if label == references[group]:
                continue
            ratios = result[(quant.accessions, units[group], label)]["ratios"]
            ratios[ratio_name(label, references[group])] = values.get("ratio")
            if "ratio_normalized" in values:
                ratios[ratio_name(label, references[group], normalized=True)] = values["ratio_normalized"]
    return result
