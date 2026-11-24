"""Where a consensusXML identification was measured: its run and its scan."""

from __future__ import annotations

import hashlib
from collections import Counter

from qpx.converters.openms_consensus.channels import column_labels
from qpx.converters.openms_consensus.protein_groups import identification_identifier
from qpx.core.files import run_file_stem as _run_stem
from qpx.core.scan import scan_from_native_id

# scan is a list<int32>; keep any surrogate within the signed 32-bit range.
_INT32_MASK = 0x7FFFFFFF


def _surrogate_scan(spectrum_ref: str) -> int:
    """Deterministic int32 surrogate ordinal for a nativeID with no scan token.

    Derived only from the spectrum reference string, so the same spectrum maps to
    the same value in both the pyopenms and streaming paths. Used as a last resort
    (instead of dropping the PSM) for nativeID schemes that expose no recognizable
    ordinal at all.
    """
    digest = hashlib.blake2s(str(spectrum_ref).encode("utf-8"), digest_size=4).digest()
    return int.from_bytes(digest, "big") & _INT32_MASK


def _scan_of(spectrum_ref: str) -> list[int]:
    """Parse numeric native ID components with the shared QPX scan convention."""
    ref = str(spectrum_ref or "")
    if not ref:
        return []
    scans = scan_from_native_id(ref)
    if scans:
        return scans
    return [_surrogate_scan(ref)]


def primary_run_stems(protein_identification) -> list[str]:
    """Run stems of one ProteinIdentification's primary MS run paths (``spectra_data``), in order."""
    get_paths = protein_identification.getPrimaryMSRunPath
    paths: list = []
    try:
        get_paths(paths)  # pyopenms < 3.6 fills an output list
    except TypeError:
        paths = get_paths()  # pyopenms >= 3.6 returns the paths
    return [_run_stem(p.decode() if isinstance(p, (bytes, bytearray)) else str(p)) for p in paths]


def column_runs(cm) -> dict[int, str]:
    """Map index -> run_file_name from the consensusXML column headers.

    The single source of the column -> run mapping for every adapter. A lone column
    without a filename (a single-run map promoted by FileConverter) takes the one
    primary MS run all ProteinIdentifications record, when there is exactly one.
    """
    headers = cm.getColumnHeaders()
    runs = {idx: _run_stem(headers[idx].filename) for idx in headers}
    if len(runs) == 1 and not next(iter(runs.values())):
        id_runs = {run for prot in cm.getProteinIdentifications() for run in primary_run_stems(prot)}
        if len(id_runs) == 1:
            runs = dict.fromkeys(runs, id_runs.pop())
    return runs


def _identifier_of(identification) -> str:
    """Identification run identifier, or ``""`` for objects that carry none."""
    try:
        return identification_identifier(identification) or ""
    except AttributeError:
        return ""


def _merge_index_runs(cm) -> dict[str, list[str]]:
    """Run stems indexed by ``id_merge_index``, per ProteinIdentification identifier.

    When OpenMS merges several identification runs, each moved PeptideIdentification
    is tagged with ``id_merge_index`` = the position of its ORIGINAL MS run in the
    primary MS run paths (the ``spectra_data`` StringList) of the ProteinIdentification
    it references — NOT a map-column index. Resolving through the PID's own
    ProteinIdentification (matched by identifier) mirrors OpenMS'
    ``IdentifierMSRunMapper``. Identifiers without recorded run paths are omitted.
    """
    runs: dict[str, list[str]] = {}
    for prot in cm.getProteinIdentifications():
        stems = primary_run_stems(prot)
        if stems:
            runs.setdefault(_identifier_of(prot), []).extend(stems)
    return runs


def _run_resolver(cm):
    """Build a callable mapping a PeptideIdentification to its run_file_name.

    ``map_index`` is a global map-column index (label-free: one map per run, so it
    is authoritative) — except when merged identifications were copied into every
    run's map (FeatureFinderIdentification fed a group-merged idXML, then linked).
    Then several ProteinIdentifications record the same multi-run ``spectra_data``
    (OpenMS' ``IdentifierMSRunMapper`` rejects that layout, so quantms output never
    has it), ``map_index`` names the map the copy landed in, and ``id_merge_index``
    names the spectrum's run. That run wins when it is one of the label-free map
    runs; when the two agree nothing changes.

    ``id_merge_index`` is otherwise only used without a ``map_index``: through the
    PID's own ProteinIdentification (:func:`_merge_index_runs`), falling back to
    all run paths in document order when the own list cannot hold the index
    (e.g. one ProteinIdentification per run). Callers with a consensus feature
    pass ``cf_runs`` (the runs its positive-intensity elements map to); with
    exactly one such run it is authoritative for that feature's PIDs.
    """
    map_run = column_runs(cm)
    map_runs = set(map_run.values())
    label_free_maps = {idx for idx, label in column_labels(cm).items() if label == "LFQ"}
    sole_run = next(iter(map_runs)) if len(map_runs) == 1 else None
    runs_by_identifier = _merge_index_runs(cm)
    run_lists = [primary_run_stems(prot) for prot in cm.getProteinIdentifications()]
    all_runs = [run for runs in run_lists for run in runs]
    list_counts = Counter(tuple(runs) for runs in run_lists if len(runs) > 1)
    copied_lists = {runs for runs, count in list_counts.items() if count > 1}

    def merged_run(pid) -> tuple[str | None, bool]:
        """The id_merge_index run, and whether it comes from a copied multi-run list."""
        if not pid.metaValueExists("id_merge_index"):
            return None, False
        idx = int(pid.getMetaValue("id_merge_index"))
        own = runs_by_identifier.get(_identifier_of(pid), [])
        if 0 <= idx < len(own):
            return own[idx], tuple(own) in copied_lists
        return (all_runs[idx] if 0 <= idx < len(all_runs) else None), False

    def resolve(pid, cf_runs=None) -> str | None:
        origin, copied_ids = merged_run(pid)
        if pid.metaValueExists("map_index"):
            idx = int(pid.getMetaValue("map_index"))
            run = map_run.get(idx)
            if run:
                moved = copied_ids and idx in label_free_maps and origin in map_runs
                return origin if moved else run
        # Assigned PID: fall back to the consensus feature's single element run.
        if cf_runs and len(cf_runs) == 1:
            return next(iter(cf_runs))
        # Unassigned PID in a merged multi-run consensusXML: id_merge_index selects
        # the original MS run, NOT a map column.
        if origin:
            return origin
        return sole_run

    return resolve
