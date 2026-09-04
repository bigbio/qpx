"""Streaming (low-memory) consensusXML reader.

``pyopenms.ConsensusXMLFile().load()`` pulls the whole ``ConsensusMap`` into
memory (~0.8x the file size), which does not scale to 100+ GB files. This module
provides :class:`StreamingConsensusMap`, which mimics the exact pyopenms API
surface the ``openms_consensus`` adapters use, backed by a defused ``iterparse``:

* the small identification header (protein hits + indistinguishable groups) and
  the column map list are parsed once and held in memory;
* the large ``<consensusElement>`` list is streamed on each iteration, one element
  at a time (cleared as we go), and the unassigned peptide identifications stream
  the same way.

Because the shim exposes the same methods the adapters call, the existing,
validated adapters run **unchanged** on streamed data — so the output is
identical to the pyopenms path; only peak memory changes. ``getSequence()``
returns a real ``pyopenms.AASequence`` (built per-peptide, no map), so ProForma /
modification handling is byte-identical.
"""

from __future__ import annotations

import numpy as np
from defusedxml.ElementTree import iterparse


def _localname(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _f32(value: str | float) -> float:
    """Parse to float32 then back to float, matching pyopenms's ``float`` storage.

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


# ---------------------------------------------------------------------------
# Lightweight shims mimicking the pyopenms objects the adapters read
# ---------------------------------------------------------------------------


class _ColHeader:
    __slots__ = ("filename", "label")

    def __init__(self, filename: str, label: str):
        self.filename = filename
        self.label = label


class _MetaMixin:
    """metaValueExists / getMetaValue backed by a ``_meta`` dict of strings."""

    _meta: dict[str, str]

    def metaValueExists(self, key: str) -> bool:
        return key in self._meta

    def getMetaValue(self, key: str):
        return self._meta.get(key)


class _Evidence:
    __slots__ = ("_acc",)

    def __init__(self, acc: str):
        self._acc = acc

    def getProteinAccession(self) -> str:
        return self._acc


class _PeptideHit(_MetaMixin):
    __slots__ = ("_seq", "_charge", "_accs", "_meta", "_score")

    def __init__(self, seq: str, charge: int, accs: list[str], meta: dict[str, str], score: float):
        self._seq = seq
        self._charge = charge
        self._accs = accs
        self._meta = meta
        self._score = score

    def getSequence(self):
        import pyopenms as oms

        return oms.AASequence.fromString(self._seq)

    def getCharge(self) -> int:
        return self._charge

    def getScore(self):
        return self._score

    def getPeptideEvidences(self) -> list[_Evidence]:
        return [_Evidence(a) for a in self._accs]


class _PeptideIdentification(_MetaMixin):
    __slots__ = ("_hits", "_mz", "_rt", "_ref", "_score_type", "_higher", "_meta")

    def __init__(self, hits, mz, rt, ref, score_type, higher, meta):
        self._hits = hits
        self._mz = mz
        self._rt = rt
        self._ref = ref
        self._score_type = score_type
        self._higher = higher
        self._meta = meta

    def getHits(self):
        return self._hits

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
    __slots__ = ("_map", "_rt", "_it")

    def __init__(self, map_index: int, rt: float, intensity: float):
        self._map = map_index
        self._rt = rt
        self._it = intensity

    def getMapIndex(self) -> int:
        return self._map

    def getRT(self) -> float:
        return self._rt

    def getIntensity(self) -> float:
        return self._it


class _ConsensusFeature:
    __slots__ = ("_charge", "_mz", "_rt", "_subs", "_pids")

    def __init__(self, charge, mz, rt, subs, pids):
        self._charge = charge
        self._mz = mz
        self._rt = rt
        self._subs = subs
        self._pids = pids

    def getCharge(self):
        return self._charge

    def getMZ(self):
        return self._mz

    @property
    def consensus_rt(self):
        """Parent ConsensusFeature retention time."""
        return self._rt

    def getFeatureList(self):
        return self._subs

    def getPeptideIdentifications(self):
        return self._pids


class _ProteinHit(_MetaMixin):
    __slots__ = ("_acc", "_score", "_desc", "_meta")

    def __init__(self, acc, score, desc, meta):
        self._acc = acc
        self._score = score
        self._desc = desc
        self._meta = meta

    def getAccession(self):
        return self._acc

    def getScore(self):
        return self._score

    def getDescription(self):
        # pyopenms exposes the FASTA header via the "Description" UserParam.
        return self._meta.get("Description", self._desc)


class _Group:
    __slots__ = ("accessions",)

    def __init__(self, accessions: list[str]):
        self.accessions = accessions


class _ProteinIdentification:
    __slots__ = ("_hits", "_groups", "_score_type", "_run_paths")

    def __init__(self, hits, groups, score_type, run_paths=None):
        self._hits = hits
        self._groups = groups
        self._score_type = score_type
        self._run_paths = run_paths or []

    def getHits(self):
        return self._hits

    def getIndistinguishableProteins(self):
        return self._groups

    def getScoreType(self):
        return self._score_type

    def getPrimaryMSRunPath(self, output):
        # Mirror pyopenms: append the merged run's spectra_data (as bytes) so the
        # id_merge_index -> run mapping matches between both read paths.
        output.extend(p.encode() if isinstance(p, str) else p for p in self._run_paths)


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------


def _parse_peptide_hit(hit_el, ph_to_acc: dict[str, str]) -> _PeptideHit:
    meta = _user_params(hit_el)
    refs = hit_el.attrib.get("protein_refs", "").split()
    accs = [ph_to_acc[r] for r in refs if r in ph_to_acc]
    charge = int(hit_el.attrib.get("charge") or 0)
    score = _f32(hit_el.attrib["score"]) if hit_el.attrib.get("score") not in (None, "") else None
    return _PeptideHit(hit_el.attrib.get("sequence", ""), charge, accs, meta, score)


def _parse_peptide_id(pid_el, ph_to_acc: dict[str, str]) -> _PeptideIdentification:
    hits = [_parse_peptide_hit(h, ph_to_acc) for h in pid_el if _localname(h.tag) == "PeptideHit"]
    a = pid_el.attrib
    mz = float(a["MZ"]) if a.get("MZ") not in (None, "") else 0.0
    rt = float(a["RT"]) if a.get("RT") not in (None, "") else 0.0
    higher = a.get("higher_score_better", "false").lower() == "true"
    # Only the map_index / id_merge_index metas are read off a PID by the adapters.
    meta = _user_params(pid_el)
    return _PeptideIdentification(hits, mz, rt, a.get("spectrum_reference", ""), a.get("score_type", ""), higher, meta)


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
                        )
                    )
        elif tag == "PeptideIdentification":
            pids.append(_parse_peptide_id(child, ph_to_acc))
    return _ConsensusFeature(charge, mz, rt, subs, pids)


# ---------------------------------------------------------------------------
# The streaming map
# ---------------------------------------------------------------------------


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
        ph_hits: list[_ProteinHit] = []
        ph_meta: dict[str, str] = {}  # ProteinIdentification-level UserParams (groups live here)
        prot_score_type = ""
        for event, el in iterparse(self._path, events=("start", "end")):
            tag = _localname(el.tag)
            if event == "end" and tag == "map":
                self._headers[int(el.attrib["id"])] = _ColHeader(el.attrib.get("name", ""), el.attrib.get("label", ""))
                el.clear()
            elif event == "end" and tag == "SearchParameters":
                # The enzyme used for the search; needed to count missed cleavages.
                # Kept here because the streaming reader never revisits the header.
                enzyme = el.attrib.get("enzyme", "")
                if enzyme and not self._enzyme:
                    self._enzyme = enzyme
                el.clear()
            elif event == "start" and tag == "ProteinIdentification":
                ph_hits = []
                ph_meta = {}
                prot_score_type = el.attrib.get("score_type", "")
            elif event == "end" and tag == "ProteinHit":
                acc = el.attrib.get("accession", "")
                phid = el.attrib.get("id", "")
                self._ph_to_acc[phid] = acc
                score = _f32(el.attrib["score"]) if el.attrib.get("score") not in (None, "") else None
                ph_hits.append(_ProteinHit(acc, score, el.attrib.get("description", ""), _user_params(el)))
                el.clear()
            elif event == "end" and tag == "ProteinIdentification":
                # The indistinguishable groups are ProteinIdentification-level
                # UserParams: name "indistinguishable_proteins_N", value
                # "score,PH_a,PH_b,...". Collect them before the element clears.
                ph_meta = _user_params(el)
                self._prots.append(
                    _ProteinIdentification(ph_hits, self._build_groups(ph_meta), prot_score_type, _parse_spectra_data(ph_meta))
                )
                el.clear()
                # do NOT stop: the <mapList> comes after the IdentificationRun.
            elif event == "start" and tag == "consensusElementList":
                break  # header done; the bulk element list follows
            if event == "end" and tag in ("consensusElement", "UnassignedPeptideIdentification"):
                el.clear()  # never accumulate the bulk while scanning the header
        if not self._prots:
            self._prots.append(
                _ProteinIdentification(ph_hits, self._build_groups(ph_meta), prot_score_type, _parse_spectra_data(ph_meta))
            )

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
        # Streamed fresh; unassigned IDs are few relative to the element list.
        out: list[_PeptideIdentification] = []
        for _, el in iterparse(self._path, events=("end",)):
            tag = _localname(el.tag)
            if tag == "UnassignedPeptideIdentification":
                out.append(_parse_peptide_id(el, self._ph_to_acc))
                el.clear()
            elif tag == "consensusElement":
                el.clear()
        return out

    def __iter__(self):
        for kind, obj in self.iter_all():
            if kind == "element":
                yield obj

    def iter_all(self):
        """Single pass yielding ``("element", ConsensusFeature)`` then, in file
        order, ``("unassigned", PeptideIdentification)`` — for the single-pass
        streaming converter (one parse builds feature/psm/pg)."""
        depth_in_list = False
        for event, el in iterparse(self._path, events=("start", "end")):
            tag = _localname(el.tag)
            if event == "start" and tag == "consensusElementList":
                depth_in_list = True
            elif event == "end" and tag == "consensusElementList":
                depth_in_list = False
            elif event == "end" and tag == "consensusElement" and depth_in_list:
                yield "element", _parse_consensus_element(el, self._ph_to_acc)
                el.clear()
            elif event == "end" and tag == "UnassignedPeptideIdentification":
                yield "unassigned", _parse_peptide_id(el, self._ph_to_acc)
                el.clear()
