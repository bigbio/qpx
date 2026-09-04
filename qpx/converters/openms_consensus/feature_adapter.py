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
import logging
import re
from typing import Optional

from qpx.converters.channel_labels import normalize_label
from qpx.core.cleavage import count_missed_cleavages
from qpx.core.files import run_file_stem as _run_stem

# OpenMS isobaric channel labels look like ``tmt6plex_126`` / ``itraq4plex_114``.
_log = logging.getLogger(__name__)

# OpenMS isobaric channel labels look like ``tmt6plex_126`` / ``itraq4plex_114``.
_CHANNEL_RE = re.compile(r"(tmt|itraq)\d*plex?_(\d+[NC]?)", re.IGNORECASE)

# Meta values OpenMS uses for the posterior error probability, best first.
_PEP_META_KEYS = ("Posterior Error Probability_score", "PEP", "pep")

# Score types that mean the PeptideIdentification's primary score IS a q-value.
_QVALUE_SCORE_TYPES = ("q-value", "qvalue", "fdr")


def mass_error_ppm(calculated_mz, observed_mz) -> float | None:
    """Relative precursor mass error in ppm, or ``None`` when it is not measurable.

    ``None`` is reserved for values that are absent or not measurable. Both m/z
    values must be positive: ``feature_records_for_cf`` substitutes 0.0 when the
    ConsensusFeature carries no m/z, and treating that as a measurement would
    report roughly -1e6 ppm instead of "not measured".

    An exact match returns ``0.0``, not ``None``. The OpenMS adapters are this
    helper's only callers and their ``observed_mz`` is a real measurement, so a
    zero error is a genuine result and must not be suppressed. (A producer that
    echoed the theoretical m/z back would report a constant 0.0 ppm — but no such
    producer reaches this function, and the guard cost real measurements.)
    """
    if calculated_mz is None or observed_mz is None:
        return None
    if calculated_mz <= 0 or observed_mz <= 0:
        return None
    return float((observed_mz - calculated_mz) / calculated_mz * 1e6)


def sdrf_enzyme(sdrf_path) -> str | None:
    """First enzyme declared in the SDRF's ``comment[cleavage agent details]``.

    The SDRF is the experiment's declared ground truth and is what the DIA-NN and
    Spectronaut adapters already use, so preferring it here keeps missed-cleavage
    counts consistent across converters rather than depending on which producer
    wrote the file.
    """
    if not sdrf_path:
        return None
    try:
        from qpx.core.sdrf import SDRFHandler

        enzymes = SDRFHandler(str(sdrf_path)).get_enzymes()
        if enzymes:
            return str(enzymes[0])
    except (OSError, KeyError, TypeError, ValueError):
        _log.debug("Could not load enzyme from SDRF %s", sdrf_path)
    return None


def resolve_enzyme(cm, sdrf_path=None) -> str | None:
    """Enzyme to count missed cleavages against: SDRF first, consensusXML second.

    The SDRF wins because it is the declared experimental design and the other
    converters already key on it. The consensusXML ``SearchParameters`` is the
    fallback for a conversion run without an SDRF. A disagreement between the two
    is a metadata inconsistency worth surfacing — the same treatment
    :func:`check_channels_vs_sdrf` gives channel mismatches — but it is not fatal.
    """
    from_sdrf = sdrf_enzyme(sdrf_path)
    from_xml = consensus_enzyme(cm)
    if from_sdrf and from_xml and from_sdrf.strip().lower() != from_xml.strip().lower():
        _log.warning(
            "enzyme mismatch: SDRF declares %r, consensusXML SearchParameters says %r; using the SDRF for missed-cleavage counts",
            from_sdrf,
            from_xml,
        )
    return from_sdrf or from_xml


def consensus_enzyme(cm) -> str | None:
    """Digestion enzyme used for the search, or ``None`` when it is not recorded.

    Read from the consensusXML ``<SearchParameters enzyme="...">``. Needed to count
    missed cleavages, which the OpenMS path otherwise leaves null while the DIA-NN
    path populates it (bigbio/qpx#300).

    Works for both readers: the streaming map captures the attribute while parsing
    its header, and a pyopenms ConsensusMap exposes it through the protein
    identification's search parameters. ``unknown_enzyme`` counts as absent.
    """
    enzyme = getattr(cm, "_enzyme", None)
    if not enzyme:
        try:
            prots = cm.getProteinIdentifications()
            if prots:
                params = prots[0].getSearchParameters()
                digestion = getattr(params, "digestion_enzyme", None)
                if digestion is not None:
                    enzyme = digestion.getName()
        except (AttributeError, IndexError, RuntimeError):
            return None
    if not enzyme:
        return None
    enzyme = str(enzyme).strip()
    if not enzyme or enzyme.lower() in ("unknown_enzyme", "unknown", "no enzyme"):
        return None
    return enzyme


def pep_of(hit) -> float | None:
    """Posterior error probability for a PeptideHit, or ``None`` when absent."""
    for mv in _PEP_META_KEYS:
        if hit.metaValueExists(mv):
            return float(hit.getMetaValue(mv))
    return None


def qvalue_of(hit, score_type: str) -> float | None:
    """Peptide-level q-value for a PeptideHit, or ``None`` when unavailable.

    quantms runs FDR before export, so the PeptideIdentification's primary score
    is usually the q-value itself — but only when ``score_type`` says so. A raw
    search-engine score must never be written to ``peptide_qvalue``: a consumer
    cannot tell the two apart once it is in the column, and a Percolator SVM
    score silently read as an FDR is worse than a null.
    """
    if str(score_type or "").lower() not in _QVALUE_SCORE_TYPES:
        return None
    score = hit.getScore()
    return float(score) if score is not None else None


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
        # N-terminal position is always 0 (amino_acid null); see the QPX
        # modifications spec.
        _add_mod(nterm.getId() or "", _acc(nterm.getUniModAccession()), 0, None, None)
    cterm = aa_sequence.getCTerminalModification() if hasattr(aa_sequence, "getCTerminalModification") else None
    if cterm:
        # C-terminal position is always length + 1 (amino_acid null).
        _add_mod(cterm.getId() or "", _acc(cterm.getUniModAccession()), aa_sequence.size() + 1, None, None)
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


def _pid_run(pid, map_info: dict[int, tuple[str, str]], cf_runs: Optional[set[str]] = None) -> str | None:
    """Resolve one identification to the run that produced it."""
    pid_run = None
    if pid.metaValueExists("map_index"):
        pid_run = map_info.get(int(pid.getMetaValue("map_index")), (None, None))[0]
    if pid_run is None and cf_runs and len(cf_runs) == 1:
        pid_run = next(iter(cf_runs))
    if pid_run is None:
        runs_in_map = {info[0] for info in map_info.values()}
        if len(runs_in_map) == 1:
            pid_run = next(iter(runs_in_map))
    return pid_run


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
    scan_by_run: dict[str, list[int]] = {}
    for pid in pids:
        scans = _pid_scans(pid)
        if not scans:
            continue
        pid_run = _pid_run(pid, map_info, cf_runs)
        if pid_run is not None:
            run_scans = scan_by_run.setdefault(pid_run, [])
            run_scans.extend(scan for scan in scans if scan not in run_scans)
    return scan_by_run


def _confidence_by_run(
    pids,
    map_info: dict[int, tuple[str, str]],
    cf_runs: Optional[set[str]] = None,
) -> dict[str, tuple[float | None, float | None]]:
    """Return the first identification's PEP and q-value for each run."""
    confidence_by_run: dict[str, tuple[float | None, float | None]] = {}
    for pid in pids:
        hits = pid.getHits()
        pid_run = _pid_run(pid, map_info, cf_runs)
        if not hits or pid_run is None or pid_run in confidence_by_run:
            continue
        hit = hits[0]
        score_type = str(pid.getScoreType() or "")
        confidence_by_run[pid_run] = (pep_of(hit), qvalue_of(hit, score_type))
    return confidence_by_run


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


def consensus_features_to_records(consensusxml_path: str | None = None, cm=None, group_map=None, sdrf_path=None) -> list[dict]:
    """Return QPX feature record dicts extracted from a consensusXML.

    Pass either ``consensusxml_path`` (loaded here) or an already-loaded ``cm``.
    ``sdrf_path`` is optional and only supplies the digestion enzyme for
    ``missed_cleavages``; without it the enzyme falls back to the consensusXML
    ``SearchParameters``.

    ``group_map`` (accession -> full protein-group membership, from
    ``pg_adapter.accession_to_group``) makes each feature stamp BOTH the same
    group leader the pg view uses as ``anchor_protein`` AND the full
    ``pg_accessions`` membership, so the feature->pg join is unambiguous even
    when two distinct groups share a leader; without it the anchor falls back to
    the peptide's first protein evidence and ``pg_accessions`` is null.
    """
    cm = cm if cm is not None else load_consensus_map(consensusxml_path)
    map_info = feature_map_info(cm)
    # Resolve the enzyme once, as the converter does, so this public entry point
    # reports missed_cleavages too. Without it a caller of the library API got a
    # null column while the same data through the CLI got a value.
    enzyme = resolve_enzyme(cm, sdrf_path)
    records: list[dict] = []
    for cf in cm:
        records.extend(feature_records_for_cf(cf, map_info, group_map, enzyme=enzyme))
    return records


def feature_map_info(cm) -> dict[int, tuple[str, str]]:
    """Map index -> (run_file_name, label) from the consensusXML column headers.

    Label comes from the map label itself (isobaric channel -> canonical, else
    LFQ) — experiment_type is unreliable (quantms writes "label-free" for TMT).
    """
    headers = cm.getColumnHeaders()
    return {idx: (_run_stem(headers[idx].filename), _map_label(headers[idx].label)) for idx in headers}


def feature_records_for_cf(cf, map_info: dict[int, tuple[str, str]], group_map=None, enzyme=None) -> list[dict]:
    """Feature records for one consensus feature (one per run, channels as intensities).

    ``pg_accessions`` carries the full protein-group membership; the feature->pg
    association is computed on read via the ``Dataset`` softlink (bigbio/qpx#269),
    so ``pg_ids`` is left unset here.
    """
    pids = cf.getPeptideIdentifications()
    if not pids or not pids[0].getHits():
        return []
    by_run = _group_subfeatures_by_run(cf, map_info)
    cf_runs = set(by_run)
    scan_by_run = _scan_by_run(pids, map_info, cf_runs=cf_runs)
    confidence_by_run = _confidence_by_run(pids, map_info, cf_runs=cf_runs)
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
    orig = evidences[0].getProteinAccession() if evidences else None
    if isinstance(orig, bytes):
        orig = orig.decode()
    # Resolve to the full protein-group membership so feature.anchor_protein AND
    # feature.pg_accessions match the pg view (peptide evidence order alone does
    # not identify the leader, nor the group when leaders are shared). Keep the
    # group in the order _merge_protein_ids/_build_groups produced so pg_accessions
    # lines up row-for-row with the pg view's pg_accessions for that group.
    group = group_map.get(orig) if (group_map and orig is not None) else None
    anchor_protein = group[0] if group else orig
    pg_accessions = [{"accession": a, "start": None, "end": None, "pre": None, "post": None} for a in group] if group else None
    # A peptide is unique when its resolved group holds exactly one protein.
    # With no resolved group the answer is unknown: a lone peptide evidence is
    # not proof of uniqueness, and claiming True there would invent information.
    unique = (len(group) == 1) if group else None
    error_ppm = mass_error_ppm(calculated_mz, observed_mz) if charge > 0 else None
    # Missed cleavages are a property of the peptide and the search enzyme, both
    # of which are in hand; the DIA-NN path already reports them (bigbio/qpx#300).
    missed = count_missed_cleavages(sequence, enzyme) if enzyme else None

    records: list[dict] = []
    for run, entry in by_run.items():
        intensities = [{"label": label, "intensity": inten} for label, inten in entry["labels"].items()]
        posterior_error_probability, peptide_qvalue = confidence_by_run.get(run, (None, None))
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
                "posterior_error_probability": posterior_error_probability,
                "peptide_qvalue": peptide_qvalue,
                "is_decoy": is_decoy,
                "calculated_mz": calculated_mz,
                "observed_mz": observed_mz,
                "mass_error_ppm": error_ppm,
                "missed_cleavages": missed,
                "unique": unique,
                "consensus_rt": consensus_rt,
                "anchor_protein": anchor_protein,
                "pg_accessions": pg_accessions,
                "additional_scores": additional_scores,
            }
        )
    return records
