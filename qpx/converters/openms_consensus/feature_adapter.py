"""consensusXML -> QPX feature records (core extraction).

Uses pyopenms to read the ConsensusMap. Each consensus feature is a peptide
linked across the map columns; a map column is one *run* (label-free) or one
*channel of a run* (isobaric, ``experiment_type == labeled_MS2``). We emit one
feature record per ``(peptidoform, charge, run_file_name, rt)`` carrying the
per-channel intensities as ``intensities: list<{label, intensity}>`` — matching
the QPX feature schema. Protein-level quantity is not produced here.
"""

from __future__ import annotations

import ast
import re
from typing import Optional

from qpx.converters.channel_labels import normalize_label
from qpx.core.files import run_file_stem as _run_stem

# OpenMS isobaric channel labels look like ``tmt6plex_126`` / ``itraq4plex_114``.
_CHANNEL_RE = re.compile(r"(tmt|itraq)\d*plex?_(\d+[NC]?)", re.IGNORECASE)


def _canonical_channel(label: Optional[str]) -> str:
    """Canonicalize an OpenMS map label to the QPX channel label.

    ``tmt6plex_126`` -> ``TMT126``; ``itraq4plex_114`` -> ``ITRAQ114``. Falls back
    to :func:`normalize_label`, then the raw string, so unknown labels are kept
    rather than dropped.
    """
    if not label:
        return "LFQ"
    m = _CHANNEL_RE.search(str(label))
    if m:
        family = "TMT" if m.group(1).lower() == "tmt" else "ITRAQ"
        return normalize_label(f"{family}{m.group(2).upper()}")
    return normalize_label(str(label))


def _map_label(raw_label) -> str:
    """Resolve a consensusXML map label to a QPX channel, or ``LFQ``.

    Detection is by the label itself, NOT ``experiment_type``: quantms/OpenMS
    write ``experiment_type="label-free"`` even for isobaric runs (the maps still
    carry ``tmt6plex_126`` / ``itraq4plex_114`` labels; label-free maps carry the
    literal ``"label-free"``). Only a real isobaric channel pattern becomes a
    channel; everything else (``"label-free"``, empty, unknown) is ``LFQ``.
    """
    return _canonical_channel(raw_label) if _CHANNEL_RE.search(str(raw_label or "")) else "LFQ"


def consensus_channels(cm) -> set[str]:
    """Canonical isobaric channels present in the consensusXML maps (excludes LFQ)."""
    headers = cm.getColumnHeaders()
    return {lab for i in headers if (lab := _map_label(headers[i].label)) != "LFQ"}


def check_channels_vs_sdrf(cm, sdrf_path) -> list[str]:
    """Return warnings where the consensusXML channels disagree with the SDRF.

    The consensusXML map labels are the channel identity the converter uses; the
    SDRF ``comment[label]`` set is the experiment's declared ground truth. A
    channel present in one but not the other is a metadata inconsistency worth
    surfacing (e.g. a mis-declared plex, or the wrong SDRF). Returns human-readable
    messages — empty when the two agree, or when the SDRF has no label column.
    """
    from qpx.converters.channel_labels import read_sdrf_labels

    sdrf_raw = read_sdrf_labels(sdrf_path)
    if not sdrf_raw:
        return []
    sdrf_channels = {normalize_label(v) for v in sdrf_raw} - {"LFQ"}
    cm_channels = consensus_channels(cm)
    messages: list[str] = []
    only_cm = sorted(cm_channels - sdrf_channels)
    only_sdrf = sorted(sdrf_channels - cm_channels)
    if only_cm:
        messages.append(f"channel(s) in the consensusXML but not the SDRF comment[label]: {only_cm}")
    if only_sdrf:
        messages.append(f"channel(s) in the SDRF comment[label] but not the consensusXML: {only_sdrf}")
    return messages


def to_proforma(aa_sequence) -> str:
    """Convert an OpenMS ``AASequence`` to a ProForma peptidoform string.

    OpenMS ``toUniModString()`` uses ``(UniMod:N)`` with a leading/trailing ``.``
    for terminal mods (e.g. ``.(UniMod:737)THSQEEM(UniMod:35)QHMQR``); ProForma
    uses ``[UNIMOD:N]`` and ``[..]-`` / ``-[..]`` for N-/C-terminal mods
    (``[UNIMOD:737]-THSQEEM[UNIMOD:35]QHMQR``).
    """
    s = aa_sequence.toUniModString().replace("UniMod:", "UNIMOD:")
    s = re.sub(r"^\.\(([^)]+)\)", r"[\1]-", s)  # N-terminal
    s = re.sub(r"\.\(([^)]+)\)$", r"-[\1]", s)  # C-terminal
    s = re.sub(r"\(([^)]+)\)", r"[\1]", s)  # internal residue mods
    return s


def _parse_site_scores(raw) -> dict[int, float]:
    """Parse a Luciphor ``Luciphor_site_scores`` UserParam into {position: score}.

    The value is a stringified dict (e.g. ``{1: -31.36, 9: 31.36}``) keyed by
    1-based amino-acid index. Positional ``None`` values are dropped so a site
    without a localisation score contributes no score entry.
    """
    if not raw:
        return {}
    try:
        parsed = ast.literal_eval(str(raw))
        if not isinstance(parsed, dict):
            return {}
        out = {}
        for k, v in parsed.items():
            try:
                pos = int(k)
            except (TypeError, ValueError):
                continue
            if v is not None:
                try:
                    out[pos] = float(v)
                except (TypeError, ValueError):
                    pass
        return out
    except (ValueError, SyntaxError, TypeError):
        # Malformed literal (ast.literal_eval) or a non-mapping value -> no scores.
        return {}


def to_modifications(aa_sequence, site_scores: dict[int, list[dict]] | None = None) -> list[dict]:
    """Convert a pyopenms ``AASequence`` to the QPX ``modifications`` structure.

    Each modification is ``{name, accession, positions: [{position, amino_acid,
    scores}]}`` — the residue position (1-based) plus, when a site-localisation
    score dict is supplied, the per-site localisation score(s) (as QPX ``score``
    structs). N- and C-terminal modifications are reported at position 1 /
    ``len(sequence)`` with an ``amino_acid`` of ``None`` (no single residue
    carries them), mirroring ``to_proforma``'s ``[UNIMOD:N]-`` / ``-[UNIMOD:N]``
    rendering.
    """
    mods: list[dict] = []
    seen: dict[tuple[str, str | None], dict] = {}

    def _add_mod(name: str, accession: str | None, position: int, amino_acid: str | None, scores: list[dict] | None) -> None:
        key = (name, accession)
        entry = seen.get(key)
        if entry is None:
            entry = {"name": name, "accession": accession, "positions": []}
            seen[key] = entry
            mods.append(entry)
        entry["positions"].append({"position": position, "amino_acid": amino_acid, "scores": scores})

    def _acc(accession) -> str | None:
        # Normalise pyopenms "UniMod:N" to QPX's canonical "UNIMOD:N".
        if not accession:
            return None
        return str(accession).replace("UniMod:", "UNIMOD:")

    for i in range(aa_sequence.size()):
        mod = aa_sequence[i].getModification()
        if not mod:
            continue
        _add_mod(
            mod.getId() or "",
            _acc(mod.getUniModAccession()),
            i + 1,
            aa_sequence[i].getOneLetterCode(),
            site_scores.get(i + 1) if site_scores else None,
        )

    nterm = aa_sequence.getNTerminalModification() if hasattr(aa_sequence, "getNTerminalModification") else None
    if nterm:
        _add_mod(nterm.getId() or "", _acc(nterm.getUniModAccession()), 1, None, None)
    cterm = aa_sequence.getCTerminalModification() if hasattr(aa_sequence, "getCTerminalModification") else None
    if cterm:
        _add_mod(cterm.getId() or "", _acc(cterm.getUniModAccession()), aa_sequence.size(), None, None)
    return mods or None


# Site-localisation algorithms whose UserParams OpenMS may stamp onto peptide
# hits: ``<name>_pep_score`` (peptide-level score) and ``<name>_site_scores``
# (a ``{position: score}`` dict keyed by 1-based residue index). Multiple
# algorithms may be present in one file (e.g. a run localised by Ascore and
# another by PhosphoRS), so all matches are collected rather than a single one.
_LOCALIZATION_ALGORITHMS = ("Luciphor", "Ascore", "PhosphoRS", "AScore")


def localization_scores(hit) -> tuple[list[dict] | None, dict[int, list[dict]]]:
    """Return ``(additional_scores, site_scores)`` from a hit's site-localisation UserParams.

    For each known localisation algorithm present on the hit, ``<Algo>_pep_score``
    and ``<Algo>_global_flr`` (when present) become ``additional_scores`` entries
    (the QPX ``score`` struct); ``<Algo>_site_scores`` is parsed into
    ``{position: [score, ...]}`` so :func:`to_modifications` can attach the
    per-site localisation score(s) to the matching modification position.
    """
    additional: list[dict] = []
    site_scores: dict[int, list[dict]] = {}
    for algo in _LOCALIZATION_ALGORITHMS:
        pep_key = f"{algo}_pep_score"
        if hit.metaValueExists(pep_key):
            additional.append(
                {"score_name": pep_key.lower(), "score_value": float(hit.getMetaValue(pep_key)), "higher_better": None}
            )
        flr_key = f"{algo}_global_flr"
        if hit.metaValueExists(flr_key):
            additional.append(
                {"score_name": flr_key.lower(), "score_value": float(hit.getMetaValue(flr_key)), "higher_better": None}
            )
        site_key = f"{algo}_site_scores"
        if hit.metaValueExists(site_key):
            parsed = _parse_site_scores(hit.getMetaValue(site_key))
            for pos, score in parsed.items():
                site_scores.setdefault(pos, []).append(
                    {"score_name": f"{algo.lower()}_site_score", "score_value": score, "higher_better": None}
                )
    return (additional or None), site_scores


def load_consensus_map(consensusxml_path: str):
    """Parse a consensusXML into a ``ConsensusMap`` — the single parse point.

    The feature/psm/pg adapters each only *read* the map, so one loaded instance
    is shared across all three (see ``OpenMSConsensusConverter.convert``) instead
    of re-parsing the (possibly large) XML three times.
    """
    import pyopenms as oms

    cm = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(str(consensusxml_path), cm)
    return cm


def _pid_scans(pid) -> list[int]:
    """Scan numbers parsed from one identification's spectrum reference.

    Uses the shared parser in :mod:`psm_adapter` (imported lazily to avoid a
    circular import — psm_adapter imports this module) so feature.scan and
    psm.scan are produced by the same logic, including Sciex ``cycle=`` ordinals
    and the deterministic surrogate fallback.
    """
    ref = pid.getSpectrumReference() if hasattr(pid, "getSpectrumReference") else ""
    if not ref and pid.metaValueExists("spectrum_reference"):
        ref = pid.getMetaValue("spectrum_reference")
    from qpx.converters.openms_consensus.psm_adapter import _scan_of

    return _scan_of(ref)


def _scan_by_run(pids, map_info: dict[int, tuple[str, str]], cf_runs: Optional[set[str]] = None) -> dict[str, list[int]]:
    """Attribute each identification's scan(s) to its own run.

    A consensus feature links spectra from several runs, so scans are resolved
    per ID via its ``map_index`` (falling back to the sole run only when every
    map is the same physical run, e.g. isobaric channels) rather than copying
    one ID's scan onto every run's record. In a multi-run isobaric consensusXML
    the PID carries only a local ``id_merge_index`` (no global ``map_index``); the
    caller's ``cf_runs`` (the feature's positive-intensity element runs) then
    attributes the scan — each such feature lives in a single run.
    """
    runs_in_map = {info[0] for info in map_info.values()}
    sole_run = next(iter(runs_in_map)) if len(runs_in_map) == 1 else None
    scan_by_run: dict[str, list[int]] = {}
    for pid in pids:
        scans = _pid_scans(pid)
        if not scans:
            continue
        pid_run = None
        if pid.metaValueExists("map_index"):
            pid_run = map_info.get(int(pid.getMetaValue("map_index")), (None, None))[0]
        if pid_run is None and cf_runs and len(cf_runs) == 1:
            pid_run = next(iter(cf_runs))
        pid_run = pid_run or sole_run
        if pid_run is not None:
            run_scans = scan_by_run.setdefault(pid_run, [])
            run_scans.extend(scan for scan in scans if scan not in run_scans)
    return scan_by_run


def _group_subfeatures_by_run(cf, map_info: dict[int, tuple[str, str]]) -> dict[str, dict]:
    """Group a consensus feature's sub-feature intensities by run.

    For isobaric, one run carries several channel maps (same rt); for label-free,
    one map == one run. One entry per label — if a channel is linked twice, keep
    the max rather than emit it twice (which would double-count downstream).
    """
    by_run: dict[str, dict] = {}
    for sub in cf.getFeatureList():
        intensity = float(sub.getIntensity())
        if intensity <= 0:
            continue
        run, label = map_info.get(sub.getMapIndex(), (None, None))
        if run is None:
            continue
        entry = by_run.setdefault(run, {"rt": float(sub.getRT()), "labels": {}})
        entry["labels"][label] = max(entry["labels"].get(label, 0.0), intensity)
    return by_run


def consensus_features_to_records(consensusxml_path: str | None = None, cm=None, anchor_map=None) -> list[dict]:
    """Return QPX feature record dicts extracted from a consensusXML.

    Pass either ``consensusxml_path`` (loaded here) or an already-loaded ``cm``.
    ``anchor_map`` (accession -> protein-group leader, from
    ``pg_adapter.accession_to_anchor``) makes each feature's ``anchor_protein``
    the same group leader the pg view uses, so the feature->pg join is reliable;
    without it the anchor falls back to the peptide's first protein evidence.
    """
    cm = cm if cm is not None else load_consensus_map(consensusxml_path)
    map_info = feature_map_info(cm)
    records: list[dict] = []
    for cf in cm:
        records.extend(feature_records_for_cf(cf, map_info, anchor_map))
    return records


def feature_map_info(cm) -> dict[int, tuple[str, str]]:
    """Map index -> (run_file_name, label) from the consensusXML column headers.

    Label comes from the map label itself (isobaric channel -> canonical, else
    LFQ) — experiment_type is unreliable (quantms writes "label-free" for TMT).
    """
    headers = cm.getColumnHeaders()
    return {idx: (_run_stem(headers[idx].filename), _map_label(headers[idx].label)) for idx in headers}


def feature_records_for_cf(cf, map_info: dict[int, tuple[str, str]], anchor_map=None) -> list[dict]:
    """Feature records for one consensus feature (one per run, channels as intensities)."""
    pids = cf.getPeptideIdentifications()
    if not pids or not pids[0].getHits():
        return []
    by_run = _group_subfeatures_by_run(cf, map_info)
    cf_runs = set(by_run)
    scan_by_run = _scan_by_run(pids, map_info, cf_runs=cf_runs)
    hit = pids[0].getHits()[0]
    seq_obj = hit.getSequence()
    peptidoform = to_proforma(seq_obj)
    sequence = seq_obj.toUnmodifiedString()
    additional_scores, site_scores = localization_scores(hit)
    modifications = to_modifications(seq_obj, site_scores)
    charge = int(cf.getCharge() or hit.getCharge() or 0)
    is_decoy = False
    if hit.metaValueExists("target_decoy"):
        is_decoy = "decoy" in str(hit.getMetaValue("target_decoy")).lower()
    observed_mz = float(cf.getMZ()) if cf.getMZ() else 0.0
    streamed_consensus_rt = getattr(cf, "consensus_rt", None)
    consensus_rt = float(cf.getRT() if streamed_consensus_rt is None else streamed_consensus_rt)
    calculated_mz = float(seq_obj.getMZ(charge)) if charge else observed_mz
    evidences = hit.getPeptideEvidences()
    anchor_protein = evidences[0].getProteinAccession() if evidences else None
    if isinstance(anchor_protein, bytes):
        anchor_protein = anchor_protein.decode()
    # Resolve to the protein-group leader so feature.anchor_protein matches the
    # pg view (peptide evidence order alone does not identify the leader).
    if anchor_map and anchor_protein is not None:
        anchor_protein = anchor_map.get(anchor_protein, anchor_protein)

    records: list[dict] = []
    for run, entry in _group_subfeatures_by_run(cf, map_info).items():
        intensities = [{"label": label, "intensity": inten} for label, inten in entry["labels"].items()]
        records.append(
            {
                "sequence": sequence,
                "peptidoform": peptidoform,
                "modifications": modifications,
                "charge": charge,
                "run_file_name": run,
                "scan": scan_by_run.get(run, []),
                "rt": entry["rt"],
                "intensities": intensities,
                "is_decoy": is_decoy,
                "calculated_mz": calculated_mz,
                "observed_mz": observed_mz,
                "consensus_rt": consensus_rt,
                "anchor_protein": anchor_protein,
                "additional_scores": additional_scores,
            }
        )
    return records
