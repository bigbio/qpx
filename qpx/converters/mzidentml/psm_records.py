"""QPX PSM records and identities from parsed mzIdentML data."""

from __future__ import annotations

import copy
import logging
import re
from typing import NamedTuple, Optional

from qpx.converters.mzidentml.crosslinks import _build_interlink, _build_looplink, _group_siis
from qpx.converters.mzidentml.scores import _attach_site_localization_scores, _extract_pep, _extract_scores
from qpx.core.cv_terms import MULTI_SPECTRA_ACCESSIONS
from qpx.core.data.identity import derive_id
from qpx.core.scores import normalize_score_name

logger = logging.getLogger(__name__)


class _SpectrumContext(NamedTuple):
    """Spectrum reference and retention time shared by its SII groups."""

    scan: list[int]
    run_file: str
    rt: float | None


def scan_schemes(records: list[dict]) -> set[str | None]:
    """The formats of the spectrum references of *records*."""
    return {record.get("_spectrum_id_scheme") for record in records}


def declare_scan_format(writer, schemes: set[str | None]) -> None:
    """Declare the scan format of a PSM file when all its spectrum references have the same one."""
    if len(schemes) == 1 and None not in schemes:
        writer.set_scan_format(next(iter(schemes)))
    else:
        logger.warning("Spectrum references have formats %s; scan_format is left undeclared", sorted(map(str, schemes)))


def _build_psm_records(parsed: dict, rank1_only: bool = False) -> list[dict]:
    """Convert parsed mzIdentML data into QPX PSM record dicts, optionally only the first-ranked hits."""
    spectra_data = parsed["spectra_data"]
    spectrum_results = parsed["spectrum_results"]

    records = []
    for sir in spectrum_results:
        scan, scan_scheme = _parse_native_id(sir["spectrum_id"])
        run_file = spectra_data.get(sir["spectra_data_ref"], "unknown")
        rt = sir.get("rt")
        siis = sir["siis"]
        context = _SpectrumContext(scan, run_file, rt)

        # Classify SIIs
        groups = _group_siis(siis)

        for group in groups:
            if rank1_only and group["alpha"]["rank"] > 1:
                continue
            record = _build_one_psm(group, context, parsed)
            if record is not None:
                # Transient hint (dropped by the writer) telling
                # ``_attach_spectra`` whether ``scan`` is an acquisition scan
                # number or a 0-based file index, so it fetches by the right key.
                record["_spectrum_id_scheme"] = scan_scheme
                records.append(record)

    return records


def _build_one_psm(group: dict, context: _SpectrumContext, parsed: dict) -> Optional[dict]:
    """Build a single QPX PSM dict from a classified SII group."""
    xl_type = group["xl_type"]
    alpha_sii = group["alpha"]
    beta_sii = group.get("beta")
    peptides = parsed["peptides"]
    linker_info = parsed["linker_info"]
    alpha_pep = peptides.get(alpha_sii["peptide_ref"], {})
    sequence = alpha_pep.get("sequence", "")
    if not sequence:
        return None

    confidence_fields = _psm_confidence_fields(alpha_sii, parsed)
    cross_links = None
    if xl_type == "intra":
        cross_links = _build_looplink(alpha_pep, linker_info)
    elif xl_type in ("inter", "noncovalent"):
        cross_links = _build_interlink(alpha_pep, beta_sii, peptides, linker_info, xl_type)
    peptidoform = alpha_pep["peptidoform"]

    # Peptide modifications are shared across PSMs; copy before adding per-SII scores.
    modifications, candidates = _attach_site_localization_scores(
        copy.deepcopy(alpha_pep.get("modifications")), sequence, alpha_sii["cv_params"], parsed["modification_tags"]
    )
    spectrum_fields = _psm_spectrum_fields(alpha_sii, context)
    source_fields = _psm_source_fields(context.run_file, alpha_sii, beta_sii, parsed["peptide_evidence"])
    return {
        "sequence": sequence,
        "peptidoform": peptidoform,
        "modifications": modifications,
        "modification_candidates": candidates,
        **confidence_fields,
        **spectrum_fields,
        **source_fields,
        "cross_links": cross_links,
    }


def _resolve_proteins(
    evidence_refs: list[str],
    pep_evidence: dict,
    db_sequences: dict,
) -> list[str]:
    """Resolve PeptideEvidenceRefs to protein accessions."""
    seen = set()
    accessions = []
    for ref in evidence_refs:
        pe = pep_evidence.get(ref, {})
        db_ref = pe.get("db_ref", "")
        acc = db_sequences.get(db_ref, "")
        if acc and acc not in seen:
            seen.add(acc)
            accessions.append(acc)
    return accessions


def _psm_confidence_fields(alpha_sii, parsed) -> dict:
    """Resolve accessions, decoy status and spectrum-level confidence."""
    pep_evidence = parsed["peptide_evidence"]
    proteins = _resolve_proteins(alpha_sii["peptide_evidence_refs"], pep_evidence, parsed["db_sequences"])
    is_decoy = any(pep_evidence.get(ref, {}).get("is_decoy", False) for ref in alpha_sii["peptide_evidence_refs"])
    additional_scores = _extract_scores(alpha_sii["cv_params"])
    additional_scores.append({"score_name": "rank", "score_value": float(alpha_sii["rank"]), "higher_better": False})
    pep_value = _extract_pep(alpha_sii["cv_params"])
    return {
        "posterior_error_probability": pep_value,
        "is_decoy": is_decoy,
        "additional_scores": additional_scores,
        "pass_threshold": alpha_sii["pass_threshold"],
        "protein_accessions": proteins or None,
    }


def _psm_spectrum_fields(alpha_sii, context) -> dict:
    """Precursor m/z and mass error, spectrum location and the empty peak-array fields of one PSM."""
    _calc_mz = alpha_sii["calc_mz"]
    _obs_mz = alpha_sii["exp_mz"]
    mass_error_ppm = 1e6 * (_obs_mz - _calc_mz) / _calc_mz if _calc_mz and _obs_mz else None
    missed_cleavages = _missed_cleavages(alpha_sii["cv_params"])
    return {
        "charge": alpha_sii["charge"],
        "calculated_mz": _calc_mz,
        "observed_mz": _obs_mz,
        "mass_error_ppm": mass_error_ppm,
        "predicted_rt": None,
        "run_file_name": context.run_file,
        "scan": context.scan,
        "rt": context.rt,
        "ion_mobility": None,
        "missed_cleavages": missed_cleavages,
        "mz_array": None,
        "intensity_array": None,
        "charge_array": None,
        "ion_type_array": None,
        "ion_mobility_array": None,
    }


def _missed_cleavages(cv_params) -> int | None:
    """Missed cleavages from the first MS:1003044 cvParam; None when it is absent or not an integer."""
    # Extract missed_cleavages from MS:1003044 cvParam on the SII
    missed_cleavages = None
    for cv in cv_params:
        if cv["accession"] == "MS:1003044":
            try:
                missed_cleavages = int(cv["value"])
            except (ValueError, TypeError):
                pass
            break

    return missed_cleavages


def _psm_source_fields(run_file, alpha_sii, beta_sii, pep_evidence) -> dict:
    """Build source identity and the first evidence's positional fields."""
    psm_id, cv_params = _source_sii_identity(run_file, alpha_sii, beta_sii)
    # Which SIIs make up one identification based on several spectra (MS2 and MS3 of a cleavable link).
    cv_params += [
        {"cv_name": normalize_score_name(cv["name"]), "cv_value": cv["value"]}
        for cv in alpha_sii["cv_params"]
        if cv["accession"] in MULTI_SPECTRA_ACCESSIONS
    ]

    # Extract positional info (pre/post/start/end) from first PeptideEvidence
    for ref in alpha_sii["peptide_evidence_refs"]:
        pe = pep_evidence.get(ref, {})
        pe_pre = pe.get("pre")
        pe_post = pe.get("post")
        pe_start = pe.get("start")
        pe_end = pe.get("end")
        if pe_pre is not None:
            cv_params.append({"cv_name": "pre", "cv_value": str(pe_pre)})
        if pe_post is not None:
            cv_params.append({"cv_name": "post", "cv_value": str(pe_post)})
        if pe_start is not None:
            cv_params.append({"cv_name": "start", "cv_value": str(pe_start)})
        if pe_end is not None:
            cv_params.append({"cv_name": "end", "cv_value": str(pe_end)})
        break  # Use first PeptideEvidence for positional info

    return {"psm_id": psm_id, "cv_params": cv_params or None}


def _source_sii_identity(run_file: str, alpha_sii: dict, beta_sii: dict | None) -> tuple[int, list[dict]]:
    """Map one mzIdentML SII group to a stable QPX PSM id and trace fields."""
    source_ids = [alpha_sii["id"]]
    if beta_sii is not None:
        source_ids.append(beta_sii["id"])
    if any(not source_id for source_id in source_ids):
        raise ValueError("SpectrumIdentificationItem is missing its required id")
    cv_params = [
        {
            "cv_name": "mzidentml_spectrum_identification_item_id",
            "cv_value": source_id,
        }
        for source_id in source_ids
    ]
    return derive_id([run_file, source_ids]), cv_params


def _parse_native_id(spectrum_id: str) -> tuple[list[int], str | None]:
    """Extract the numeric spectrum reference and its nativeID scheme.

    Returns ``(values, scheme)`` where ``scheme`` is one of:

    * ``"scan"``  — the value is an acquisition scan number
      (``scan=N``, ``spectrum=N``, or a bare integer);
    * ``"index"`` — the value is a 0-based position in the spectrum file
      (``index=N``). This is NOT a scan number: it must be looked up by
      position, never by scan number, or peaks/RT from an unrelated spectrum
      whose ``SCANS=`` happens to equal the index get attached;
    * ``None``    — no numeric reference could be parsed.

    The scheme is what lets spectra attachment fetch the CORRECT spectrum.
    """
    m = re.search(r"scan=(\d+)", spectrum_id)
    if m:
        return [int(m.group(1))], "scan"
    m = re.search(r"index=(\d+)", spectrum_id)
    if m:
        return [int(m.group(1))], "index"
    m = re.search(r"spectrum=(\d+)", spectrum_id)
    if m:
        return [int(m.group(1))], "scan"
    # Bare integer spectrumID — historically treated as a scan number.
    m = re.match(r"^(\d+)$", spectrum_id.strip())
    if m:
        return [int(m.group(1))], "scan"
    # Unparseable: return an empty marker (NOT ``[0]``) so a spectrum with no
    # parseable id does not alias onto a real scan-0 spectrum, and two distinct
    # unparseable spectra are not silently collapsed onto the same sentinel.
    return [], None
