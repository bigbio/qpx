"""Streaming (low-memory) consensusXML reader.

``pyopenms.ConsensusXMLFile().load()`` pulls the whole ``ConsensusMap`` into
memory, which does not scale to 100+ GB files. This module
provides :class:`StreamingConsensusMap`, which mimics the exact pyopenms API
surface the ``openms_consensus`` adapters use, backed by a defused ``iterparse``:

* the small identification header (protein hits + indistinguishable groups) and
  the column map list are parsed once and held in memory;
* the large ``<consensusElement>`` list is streamed on each iteration, one element
  at a time (cleared as we go), and the unassigned peptide identifications stream
  the same way.

Because the shim exposes the same methods the adapters call, the existing,
validated adapters run **unchanged** on streamed data — so the rows match the
pyopenms path, with a few float values (PEP, mass error) differing in the last
bit; peak memory is about half and run time about double. ``getSequence()``
returns a real ``pyopenms.AASequence`` (built per-peptide, no map), so ProForma /
modification handling is byte-identical.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pyopenms as oms
from defusedxml.ElementTree import iterparse


def _is_gzip(path: str) -> bool:
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _iterparse(path: str, events: tuple[str, ...]):
    """``iterparse`` over a consensusXML file, decompressing gzip transparently."""
    with gzip.open(path, "rb") if _is_gzip(path) else open(path, "rb") as source:
        yield from iterparse(source, events=events)


def exceeds_xml_size(path: str, limit: int) -> bool:
    """Whether the consensusXML, decompressed, is larger than *limit* bytes."""
    if Path(path).stat().st_size > limit:
        return True
    if not _is_gzip(path):
        return False
    # gzip stores the decompressed size modulo 4 GiB, so count it up to the limit.
    decompressed = 0
    with gzip.open(path, "rb") as stream:
        while chunk := stream.read(64 * 1024**2):
            decompressed += len(chunk)
            if decompressed > limit:
                return True
    return False


def _localname(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _f32(value: str | float) -> float:
    """Parse to float32 then back to float, matching pyopenms's ``float`` storage of an intensity.

    OpenMS stores sub-feature intensities as C++ ``float`` (32-bit); reading the
    XML string as a 64-bit Python float would keep digits pyopenms rounds away, so
    quantize to float32 to stay byte-identical with the pyopenms path.
    """
    return float(np.float32(value))


def _user_params(element) -> dict[str, str]:
    """Collect the direct-child ``<UserParam>`` name->value pairs of an element."""
    params: dict[str, str] = {}
    for child in element:
        if _localname(child.tag) == "UserParam":
            params[child.attrib.get("name", "")] = child.attrib.get("value", "")
    return params


def _parse_spectra_data(prot_meta: dict[str, str]) -> list[str]:
    """Parse the ProteinIdentification ``spectra_data`` StringList into a path list.

    OpenMS serialises the merged run's primary MS run paths as a stringList
    UserParam ``value="[runA.mzML,runB.mzML]"``. Returns them in order (the
    ``id_merge_index`` ordering); empty when absent.
    """
    raw = prot_meta.get("spectra_data")
    if not raw:
        return []
    inner = raw.strip()
    if inner.startswith("[") and inner.endswith("]"):
        inner = inner[1:-1]
    return [p.strip() for p in inner.split(",") if p.strip()]


def _se_engines(search_meta: dict[str, str]) -> list[tuple[str, str]]:
    """Engines recorded as ``SE:<engine>`` search parameters, as (engine, value)."""
    return [(key[3:], value) for key, value in search_meta.items() if key.startswith("SE:")]


# Parameters that say how a tool was run on this machine, not what it computed.
_EXECUTION_PARAMETERS = {"in", "out", "log", "debug", "threads", "no_progress", "test", "force", "ini", "write_ini"}


def _tool_parameters(element) -> tuple[tuple[str, str], ...]:
    """The ``parameter: <name>`` user params of a ``dataProcessing`` element, without the execution ones."""
    prefix = "parameter: "
    return tuple(
        (name.removeprefix(prefix), value)
        for name, value in _user_params(element).items()
        if name.startswith(prefix) and name.removeprefix(prefix) not in _EXECUTION_PARAMETERS
    )


def _software(element) -> tuple[str, str, str, tuple[tuple[str, str], ...]] | None:
    """A ``dataProcessing`` element's software as (name, version, processing action, parameters); None without a name."""
    software = next((child for child in element if _localname(child.tag) == "software"), None)
    if software is None or not software.attrib.get("name"):
        return None
    action = next((child for child in element if _localname(child.tag) == "processingAction"), None)
    action_name = action.attrib.get("name", "") if action is not None else ""
    return software.attrib["name"], software.attrib.get("version", ""), action_name, _tool_parameters(element)


class ProvenanceHeader(NamedTuple):
    """The consensusXML header fields that provenance records."""

    search_engines: list[tuple[str, str]]
    original_search_engines: list[tuple[str, str]]
    data_processing: list[tuple[str, str, str, tuple[tuple[str, str], ...]]]


def read_provenance_header(path: str) -> ProvenanceHeader:
    """Read the search engines and data processing of a consensusXML.

    OpenMS writes ``dataProcessing`` and ``IdentificationRun`` before the unassigned
    identifications and the map list, so parsing stops at the first of those
    instead of walking the whole header.
    """
    header = ProvenanceHeader([], [], [])
    for event, el in _iterparse(path, events=("start", "end")):
        tag = _localname(el.tag)
        if event == "start" and tag in ("UnassignedPeptideIdentification", "mapList", "consensusElementList"):
            break
        if event == "start" and tag == "IdentificationRun" and el.attrib.get("search_engine"):
            header.search_engines.append((el.attrib["search_engine"], el.attrib.get("search_engine_version", "")))
        elif event == "end" and tag == "SearchParameters":
            header.original_search_engines.extend(_se_engines(_user_params(el)))
            el.clear()
        elif event == "end" and tag == "dataProcessing":
            software = _software(el)
            if software is not None:
                header.data_processing.append(software)
            el.clear()
        elif event == "end" and tag in ("ProteinHit", "ProteinIdentification"):
            el.clear()
    return header


# ---------------------------------------------------------------------------
# Lightweight shims mimicking the pyopenms objects the adapters read
# ---------------------------------------------------------------------------


class _MetaMixin:
    """metaValueExists / getMetaValue backed by a ``_meta`` dict of strings."""

    _meta: dict[str, str]

    def metaValueExists(self, key: str) -> bool:
        return key in self._meta

    def getMetaValue(self, key: str):
        return self._meta.get(key)


class _ColHeader(_MetaMixin):
    def __init__(self, filename: str, label: str, meta: dict[str, str]):
        self.filename = filename
        self.label = label
        self._meta = meta


class _PrecursorPosition(NamedTuple):
    """The m/z and retention time of a precursor or consensus feature."""

    mz: float
    rt: float


class _PeptideIdentification(_MetaMixin):
    __slots__ = ("_hits", "_mz", "_rt", "_ref", "_score_type", "_higher", "_meta")

    def __init__(self, hits, position: _PrecursorPosition, ref, score_type, higher, meta):
        self._hits = hits
        self._mz = position.mz
        self._rt = position.rt
        self._ref = ref
        self._score_type = score_type
        self._higher = higher
        self._meta = meta

    def getHits(self):
        return self._hits

    @property
    def identifier(self):
        """Source IdentificationRun ID, matching the streamed protein header."""
        return self._meta.get("identification_run_ref", "")

    def getMZ(self):
        return self._mz

    def getRT(self):
        return self._rt

    def getSpectrumReference(self):
        return self._ref

    def getScoreType(self):
        return self._score_type

    def isHigherScoreBetter(self) -> bool:
        return self._higher


class _SubFeature:
    __slots__ = ("_map", "_rt", "_it", "_mz")

    def __init__(self, map_index: int, rt: float, intensity: float, mz: float):
        self._map = map_index
        self._rt = rt
        self._it = intensity
        self._mz = mz

    def getMapIndex(self) -> int:
        return self._map

    def getRT(self) -> float:
        return self._rt

    def getIntensity(self) -> float:
        return self._it

    def getMZ(self) -> float:
        """Element m/z (pyopenms accessor name)."""
        return self._mz


class _ConsensusFeature(_MetaMixin):
    __slots__ = ("_charge", "_position", "_subs", "_pids", "_quality", "_uid", "_meta")

    def __init__(self, charge, position: _PrecursorPosition, subs, pids, quality=0.0, uid=0, meta=None):
        self._meta = meta or {}
        self._uid = uid
        self._charge = charge
        self._position = position
        self._subs = subs
        self._pids = pids
        self._quality = quality

    def getCharge(self):
        return self._charge

    def getQuality(self):  # pylint: disable=invalid-name
        """Consensus feature quality (pyopenms accessor name)."""
        return self._quality

    def getUniqueId(self) -> int:  # pylint: disable=invalid-name
        """Consensus feature unique id (pyopenms accessor name)."""
        return self._uid

    def getMZ(self):
        return self._position.mz

    @property
    def consensus_rt(self):
        """Parent ConsensusFeature retention time."""
        return self._position.rt

    def getFeatureList(self):
        return self._subs

    def getPeptideIdentifications(self):
        return self._pids


class _Group:
    __slots__ = ("accessions",)

    def __init__(self, accessions: list[str]):
        self.accessions = accessions


class _SearchParameters:
    """The SearchParameters UserParams, with the pyopenms meta value accessors."""

    __slots__ = ("_meta",)

    def __init__(self, meta):
        self._meta = meta or {}

    def metaValueExists(self, name):  # pylint: disable=invalid-name
        """Whether the search parameters carry ``name`` (pyopenms accessor name)."""
        return name in self._meta

    def getMetaValue(self, name):  # pylint: disable=invalid-name
        """Value of search parameter ``name`` (pyopenms accessor name)."""
        return self._meta.get(name)


class _ProteinIdentification:
    __slots__ = ("_hits", "_groups", "_score_type", "_run_paths", "_identifier", "_search_meta")

    def __init__(self, hits, groups, score_type, run_paths=None, identifier="", search_meta=None):
        self._hits = hits
        self._groups = groups
        self._score_type = score_type
        self._run_paths = run_paths or []
        self._identifier = identifier
        self._search_meta = search_meta or {}

    @property
    def identifier(self):
        """Source IdentificationRun ID referenced by peptide identifications."""
        return self._identifier

    def getHits(self):
        return self._hits

    def getIndistinguishableProteins(self):
        return self._groups

    def getScoreType(self):
        return self._score_type

    def getSearchParameters(self):  # pylint: disable=invalid-name
        """SearchParameters UserParams (pyopenms accessor name)."""
        return _SearchParameters(self._search_meta)

    def getPrimaryMSRunPath(self, output):
        # Mirror pyopenms: append the merged run's spectra_data (as bytes) so the
        # id_merge_index -> run mapping matches between both read paths.
        output.extend(p.encode() if isinstance(p, str) else p for p in self._run_paths)


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------


def _parse_hit_score(attrs: dict[str, str]) -> float:
    """Use NaN for absent scores; adapters normalize it to None, not zero.

    The hit stores it as pyopenms does, so both readers agree on any version.
    """
    return float(attrs.get("score") or "nan")


def _parse_protein_hit(element) -> oms.ProteinHit:
    attrs = element.attrib
    hit = oms.ProteinHit()
    hit.setAccession(attrs.get("accession", ""))
    hit.setScore(_parse_hit_score(attrs))
    hit.setDescription(attrs.get("description", ""))
    hit.setCoverage(float(attrs.get("coverage", "-1")))
    hit.setSequence(attrs.get("sequence", ""))
    for key, value in _user_params(element).items():
        hit.setMetaValue(key, value)
    return hit


def _parse_peptide_evidences(hit_el, ph_to_acc: dict[str, str]) -> list[oms.PeptideEvidence]:
    """Keep the parallel protein-reference and coordinate lists aligned."""
    refs = hit_el.attrib.get("protein_refs", "").split()
    starts = [int(value) for value in hit_el.attrib.get("start", "").split()]
    ends = [int(value) for value in hit_el.attrib.get("end", "").split()]
    evidences = []
    for index, ref in enumerate(refs):
        if ref not in ph_to_acc:
            continue
        evidence = oms.PeptideEvidence()
        evidence.setProteinAccession(ph_to_acc[ref])
        evidence.setStart(starts[index] if index < len(starts) else -1)
        evidence.setEnd(ends[index] if index < len(ends) else -1)
        evidences.append(evidence)
    return evidences


def _parse_peptide_hit(hit_el, ph_to_acc: dict[str, str]) -> oms.PeptideHit:
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString(hit_el.attrib.get("sequence", "")))
    hit.setCharge(int(hit_el.attrib.get("charge") or 0))
    hit.setScore(_parse_hit_score(hit_el.attrib))
    hit.setPeptideEvidences(_parse_peptide_evidences(hit_el, ph_to_acc))
    for key, value in _user_params(hit_el).items():
        hit.setMetaValue(key, value)
    return hit


def _parse_peptide_id(pid_el, ph_to_acc: dict[str, str]) -> _PeptideIdentification:
    hits = [_parse_peptide_hit(h, ph_to_acc) for h in pid_el if _localname(h.tag) == "PeptideHit"]
    a = pid_el.attrib
    mz = float(a["MZ"]) if a.get("MZ") not in (None, "") else 0.0
    rt = float(a["RT"]) if a.get("RT") not in (None, "") else 0.0
    higher = a.get("higher_score_better", "false").lower() == "true"
    meta = _user_params(pid_el)
    meta["identification_run_ref"] = a.get("identification_run_ref", "")
    return _PeptideIdentification(
        hits, _PrecursorPosition(mz, rt), a.get("spectrum_reference", ""), a.get("score_type", ""), higher, meta
    )


def _parse_consensus_element(el, ph_to_acc: dict[str, str]) -> _ConsensusFeature:
    charge = int(el.attrib.get("charge") or 0)
    mz = 0.0
    rt = 0.0
    subs: list[_SubFeature] = []
    pids: list[_PeptideIdentification] = []
    for child in el:
        tag = _localname(child.tag)
        if tag == "centroid":
            mz = float(child.attrib.get("mz") or 0.0)
            rt = float(child.attrib.get("rt") or 0.0)
        elif tag == "groupedElementList":
            for e in child:
                if _localname(e.tag) == "element":
                    subs.append(
                        _SubFeature(
                            int(e.attrib["map"]),
                            float(e.attrib.get("rt") or 0.0),
                            _f32(e.attrib.get("it") or 0.0),  # intensity is float32 in OpenMS
                            float(e.attrib.get("mz") or 0.0),
                        )
                    )
        elif tag == "PeptideIdentification":
            pids.append(_parse_peptide_id(child, ph_to_acc))
    # quality is a C++ float in OpenMS, like the sub-feature intensities
    uid = el.attrib.get("id", "").removeprefix("e_")
    quality = _f32(el.attrib.get("quality") or 0.0)
    return _ConsensusFeature(
        charge, _PrecursorPosition(mz, rt), subs, pids, quality, int(uid) if uid.isdigit() else 0, _user_params(el)
    )


# ---------------------------------------------------------------------------
# The streaming map
# ---------------------------------------------------------------------------


@dataclass
class _IdentificationHeader:
    """Fields collected while reading one protein-identification header."""

    hits: list = field(default_factory=list)
    meta: dict = field(default_factory=dict)
    score_type: str = ""
    identifier: str = ""
    search_meta: dict = field(default_factory=dict)


_PROTEIN_HEADER_EVENTS = {
    ("start", "IdentificationRun"),
    ("start", "ProteinIdentification"),
    ("end", "ProteinHit"),
    ("end", "ProteinIdentification"),
}


class StreamingConsensusMap:
    """A read-only, low-memory stand-in for ``pyopenms.ConsensusMap``.

    Only the methods the ``openms_consensus`` adapters call are implemented. The
    identification header (protein hits + groups) and column maps are parsed once;
    ``__iter__`` and ``getUnassignedPeptideIdentifications`` re-stream the file.
    """

    def __init__(self, path: str):
        self._path = str(path)
        self._headers: dict[int, _ColHeader] = {}
        self._prots: list[_ProteinIdentification] = []
        self._enzyme: str = ""
        self._ph_to_acc: dict[str, str] = {}
        self._parse_header()

    # -- header (proteins + maps), parsed once --------------------------------

    def _parse_header(self) -> None:
        header = _IdentificationHeader()
        for event, el in _iterparse(self._path, events=("start", "end")):
            tag = _localname(el.tag)
            if (event, tag) == ("end", "map"):
                self._headers[int(el.attrib["id"])] = _ColHeader(
                    el.attrib.get("name", ""), el.attrib.get("label", ""), _user_params(el)
                )
                el.clear()
            elif (event, tag) == ("end", "SearchParameters"):
                # The enzyme used for the search; needed to count missed cleavages.
                # Kept here because the streaming reader never revisits the header.
                enzyme = el.attrib.get("enzyme", "")
                if enzyme and not self._enzyme:
                    self._enzyme = enzyme
                header.search_meta = _user_params(el)
                el.clear()
            elif (event, tag) in _PROTEIN_HEADER_EVENTS:
                self._parse_protein_header(event, el, header)
            elif (event, tag) == ("end", "dataProcessing"):
                el.clear()
            elif (event, tag) == ("start", "consensusElementList"):
                break  # header done; the bulk element list follows
            if event == "end" and tag in ("consensusElement", "UnassignedPeptideIdentification"):
                el.clear()  # never accumulate the bulk while scanning the header
        if not self._prots:
            self._prots.append(
                _ProteinIdentification(
                    header.hits, self._build_groups(header.meta), header.score_type, _parse_spectra_data(header.meta)
                )
            )

    def _parse_protein_header(self, event, el, header) -> None:
        """Update identification-run and protein fields for one header event."""
        tag = _localname(el.tag)
        if (event, tag) == ("start", "IdentificationRun"):
            header.identifier = el.attrib.get("id", "")
            header.search_meta = {}
        elif (event, tag) == ("start", "ProteinIdentification"):
            header.hits = []
            header.meta = {}
            header.score_type = el.attrib.get("score_type", "")
        elif (event, tag) == ("end", "ProteinHit"):
            acc = el.attrib.get("accession", "")
            phid = el.attrib.get("id", "")
            self._ph_to_acc[phid] = acc
            header.hits.append(_parse_protein_hit(el))
            el.clear()
        elif (event, tag) == ("end", "ProteinIdentification"):
            # The indistinguishable groups are ProteinIdentification-level
            # UserParams: name "indistinguishable_proteins_N", value
            # "score,PH_a,PH_b,...". Collect them before the element clears.
            header.meta = _user_params(el)
            self._prots.append(
                _ProteinIdentification(
                    header.hits,
                    self._build_groups(header.meta),
                    header.score_type,
                    _parse_spectra_data(header.meta),
                    header.identifier,
                    header.search_meta,
                )
            )
            el.clear()
            # do NOT stop: the <mapList> comes after the IdentificationRun.

    def _build_groups(self, prot_meta: dict[str, str]) -> list[_Group]:
        groups: list[_Group] = []
        for name, value in prot_meta.items():
            if name.startswith("indistinguishable_proteins_"):
                parts = value.split(",")
                # first token is the group score; the rest are PH ids
                accs = [self._ph_to_acc[p] for p in parts[1:] if p in self._ph_to_acc]
                if accs:
                    groups.append(_Group(accs))
        return groups

    # -- pyopenms-compatible surface ------------------------------------------

    def getColumnHeaders(self) -> dict[int, _ColHeader]:
        return self._headers

    def getProteinIdentifications(self) -> list[_ProteinIdentification]:
        return list(self._prots)

    def getUnassignedPeptideIdentifications(self) -> list[_PeptideIdentification]:
        return [pid for _, pid in self.iter_identifications(include_assigned=False)]

    def iter_identifications(self, *, include_assigned: bool = True):
        """Stream identification context without constructing feature/sub-feature objects."""
        root = elements = None
        for event, el in _iterparse(self._path, events=("start", "end")):
            tag = _localname(el.tag)
            if event == "start":
                root = el if root is None else root
                elements = el if tag == "consensusElementList" else elements
            elif tag == "UnassignedPeptideIdentification":
                yield "unassigned", _parse_peptide_id(el, self._ph_to_acc)
                root.clear()
            elif tag == "PeptideIdentification" and include_assigned:
                yield "assigned", _parse_peptide_id(el, self._ph_to_acc)
                el.clear()
            elif tag == "consensusElement":
                elements.clear()

    def __iter__(self):
        for kind, obj in self.iter_all():
            if kind == "element":
                yield obj

    def iter_all(self):
        """Single pass yielding ``("element", ConsensusFeature)`` then, in file
        order, ``("unassigned", PeptideIdentification)``. The converter uses
        one traversal for identification context and another for all output views."""
        root = elements = None
        for event, el in _iterparse(self._path, events=("start", "end")):
            tag = _localname(el.tag)
            if event == "start":
                root = el if root is None else root
                elements = el if tag == "consensusElementList" else elements
            elif tag == "consensusElementList":
                elements = None
            elif tag == "consensusElement" and elements is not None:
                yield "element", _parse_consensus_element(el, self._ph_to_acc)
                # Detach the element: a cleared child still costs memory per row.
                elements.clear()
            elif tag == "UnassignedPeptideIdentification":
                yield "unassigned", _parse_peptide_id(el, self._ph_to_acc)
                root.clear()
