"""mzIdentML PSM adapter — parses mzIdentML XML (including 1.3 XL-MS) to QPX PSMs.

Uses lxml for XML parsing because pyopenms cannot handle mzIdentML 1.3
crosslinking modifications (``crosslink donor``, ``crosslink acceptor``).

//TODO: Use pyopenms in the futue when solving the issues with cross-linking modifications, as it provides a more robust and tested mzIdentML parser. For now, this custom parser is necessary to handle the specific features of mzIdentML 1.3 used in XL-MS data.
"""

from __future__ import annotations

import copy
import gzip
import logging
import re
from typing import Any, Optional

try:
    from lxml import etree
except ImportError:
    etree: Any = None

from qpx.converters.base import BaseConverter
from qpx.converters.ptm import build_proforma
from qpx.core.cv_terms import (
    CV_CROSSLINK_SII,
    CV_LOOPLINK_SII,
    CV_NONCOVALENT_SII,
    CV_PEP,
    CV_PEP_GLOBAL,
    CV_SCAN_START_TIME,
    CV_XL_ACCEPTOR,
    CV_XL_DONOR,
    SITE_LOCALIZATION_ACCESSIONS,
    SKIP_SCORE_ACCESSIONS,
)
from qpx.core.data.identity import derive_id
from qpx.core.files import run_file_stem
from qpx.core.scores import is_higher_better, normalize_score_name
from qpx.writers.psm import PsmWriter

logger = logging.getLogger(__name__)

# Alias for backward compatibility with local names
_SKIP_SCORE_ACCESSIONS = SKIP_SCORE_ACCESSIONS


class MzIdentMLPsmAdapter(BaseConverter):
    """Convert mzIdentML files (including XL-MS 1.3) to QPX PSM Parquet."""

    def convert(
        self,
        mzid_path: str,
        output_path: str,
        creator: str = "mzidentml",
    ) -> None:
        """Parse *mzid_path* and write QPX PSM records to *output_path*."""
        parsed = self._parse_mzidentml(mzid_path)
        records = self._build_psm_records(parsed)

        if not records:
            logger.warning("No PSM records produced from %s", mzid_path)
            return

        self._track_scores(records)

        with PsmWriter(output_path, creator=creator, compression=self._compression) as writer:
            writer.write_batch(records)

        logger.info(
            "Wrote %d PSMs (%d with cross-links) to %s",
            len(records),
            sum(1 for r in records if r.get("cross_links")),
            output_path,
        )

    # ------------------------------------------------------------------
    # XML parsing
    # ------------------------------------------------------------------

    def _parse_mzidentml(self, path: str) -> dict:
        """Parse mzIdentML XML into intermediate data structures.

        Transparently handles gzip-compressed files (``*.gz``).
        """
        if etree is None:
            raise ImportError("lxml is required for mzIdentML conversion. Install it with: pip install qpx[mzidentml]")
        if str(path).endswith(".gz"):
            with gzip.open(path, "rb") as fh:
                tree = etree.parse(fh)
        else:
            tree = etree.parse(path)
        root = tree.getroot()
        ns = root.nsmap.get(None, "")
        self._ns = ns

        return {
            "spectra_data": self._parse_spectra_data(root, ns),
            "db_sequences": self._parse_db_sequences(root, ns),
            "peptides": self._parse_peptides(root, ns),
            "peptide_evidence": self._parse_peptide_evidence(root, ns),
            "linker_info": self._parse_linker_info(root, ns),
            "spectrum_results": self._parse_spectrum_results(root, ns),
        }

    def _parse_spectra_data(self, root, ns: str) -> dict[str, str]:
        """Map SpectraData id → file stem."""
        result = {}
        for sd in root.iter(f"{{{ns}}}SpectraData"):
            sd_id = sd.get("id")
            location = sd.get("location", "")
            stem = run_file_stem(location)
            result[sd_id] = stem
        return result

    def _parse_db_sequences(self, root, ns: str) -> dict[str, str]:
        """Map DBSequence id → protein accession."""
        result = {}
        for db in root.iter(f"{{{ns}}}DBSequence"):
            result[db.get("id")] = db.get("accession", "")
        return result

    def _parse_peptides(self, root, ns: str) -> dict[str, dict]:
        """Parse Peptide elements into {id: {sequence, mods, xl_mods}}.

        ``xl_mods`` contains crosslink donor/acceptor modifications keyed
        by their ``value`` attribute (the cross-link pair ID).
        """
        peptides = {}
        for pep in root.iter(f"{{{ns}}}Peptide"):
            pep_id = pep.get("id")
            seq_el = pep.find(f"{{{ns}}}PeptideSequence")
            sequence = seq_el.text if seq_el is not None else ""

            mods = []
            xl_mods: dict[str, list[dict]] = {}  # keyed by xl pair value

            for mod in pep.findall(f"{{{ns}}}Modification"):
                location = int(mod.get("location", "0"))
                mass_delta = float(mod.get("monoisotopicMassDelta", "0"))
                residues = mod.get("residues", "")

                # Check for XL cvParams
                is_xl = False
                for cv in mod.findall(f"{{{ns}}}cvParam"):
                    acc = cv.get("accession", "")
                    if acc in (CV_XL_DONOR, CV_XL_ACCEPTOR):
                        is_xl = True
                        xl_value = cv.get("value", "")
                        role = "donor" if acc == CV_XL_DONOR else "acceptor"
                        xl_mods.setdefault(xl_value, []).append(
                            {
                                "role": role,
                                "position": location,
                                "mass_delta": mass_delta,
                                "residues": residues,
                            }
                        )
                        break

                if not is_xl:
                    # Regular modification
                    mod_name = ""
                    mod_accession = None
                    for cv in mod.findall(f"{{{ns}}}cvParam"):
                        mod_name = cv.get("name", "")
                        mod_accession = cv.get("accession")
                        break  # Take first cvParam as the mod identity
                    mods.append(
                        {
                            "name": mod_name or f"CHEMMOD:+{mass_delta}",
                            "accession": mod_accession,
                            "positions": [
                                {
                                    "position": location,
                                    "amino_acid": residues or None,
                                    "scores": None,
                                }
                            ],
                        }
                    )

            peptides[pep_id] = {
                "sequence": sequence,
                "modifications": mods or None,
                "xl_mods": xl_mods,
            }
        return peptides

    def _parse_peptide_evidence(self, root, ns: str) -> dict[str, dict]:
        """Map PeptideEvidence id → {peptide_ref, db_ref, is_decoy, start, end, pre, post}."""
        result = {}
        for pe in root.iter(f"{{{ns}}}PeptideEvidence"):
            start_val = pe.get("start")
            end_val = pe.get("end")
            result[pe.get("id")] = {
                "peptide_ref": pe.get("peptide_ref"),
                "db_ref": pe.get("dBSequence_ref"),
                "is_decoy": pe.get("isDecoy", "false").lower() == "true",
                "start": int(start_val) if start_val is not None else None,
                "end": int(end_val) if end_val is not None else None,
                "pre": pe.get("pre", None),
                "post": pe.get("post", None),
            }
        return result

    def _parse_linker_info(self, root, ns: str) -> dict:
        """Extract cross-linker name, accession and mass from SearchModification elements."""
        linker_info: dict[str, dict] = {}  # keyed by XL pair group number

        for sm in root.iter(f"{{{ns}}}SearchModification"):
            mass = float(sm.get("massDelta", "0"))
            cvs = sm.findall(f"{{{ns}}}cvParam")

            # Look for crosslink donor/acceptor cvParam
            xl_group = None
            xl_role = None
            linker_name = None
            linker_accession = None

            for cv in cvs:
                acc = cv.get("accession", "")
                if acc == CV_XL_DONOR:
                    xl_group = cv.get("value", "")
                    xl_role = "donor"
                elif acc == CV_XL_ACCEPTOR:
                    xl_group = cv.get("value", "")
                    xl_role = "acceptor"

            if xl_group is None:
                continue

            # Find the linker name from other cvParams in same SearchModification
            for cv in cvs:
                name = cv.get("name", "")
                acc = cv.get("accession", "")
                if name.startswith("Xlink:"):
                    linker_name = name.replace("Xlink:", "")
                    linker_accession = acc
                    break

            if xl_role == "donor" and linker_name:
                linker_info[xl_group] = {
                    "name": linker_name,
                    "accession": linker_accession,
                    "mass": mass,
                }

        # If only one linker found, use it as default
        if len(linker_info) == 1:
            default = next(iter(linker_info.values()))
            linker_info["_default"] = default
        elif linker_info:
            # Use any linker as default (they're usually the same linker)
            linker_info["_default"] = next(iter(linker_info.values()))

        return linker_info

    def _parse_spectrum_results(self, root, ns: str) -> list[dict]:
        """Parse SpectrumIdentificationResult elements."""
        results = []
        for sir in root.iter(f"{{{ns}}}SpectrumIdentificationResult"):
            spectrum_id = sir.get("spectrumID", "")
            spectra_data_ref = sir.get("spectraData_ref", "")

            # Extract retention time from SIR cvParams
            rt_seconds = None
            for cv in sir.findall(f"{{{ns}}}cvParam"):
                if cv.get("accession") == CV_SCAN_START_TIME:
                    try:
                        rt_value = float(cv.get("value", ""))
                        unit_name = cv.get("unitName", "")
                        # Convert minutes to seconds if needed
                        if "minute" in unit_name.lower():
                            rt_value *= 60.0
                        rt_seconds = rt_value
                    except (ValueError, TypeError):
                        pass
                    break

            siis = []
            for sii in sir.findall(f"{{{ns}}}SpectrumIdentificationItem"):
                sii_data = {
                    "id": sii.get("id", ""),
                    "charge": int(sii.get("chargeState", "0")),
                    "exp_mz": float(sii.get("experimentalMassToCharge", "0")),
                    "calc_mz": float(sii.get("calculatedMassToCharge", "0")),
                    "rank": int(sii.get("rank", "1")),
                    "pass_threshold": sii.get("passThreshold", "false").lower() == "true",
                    "peptide_ref": sii.get("peptide_ref", ""),
                    "cv_params": [],
                    "peptide_evidence_refs": [],
                }

                for cv in sii.findall(f"{{{ns}}}cvParam"):
                    sii_data["cv_params"].append(
                        {
                            "accession": cv.get("accession", ""),
                            "name": cv.get("name", ""),
                            "value": cv.get("value", ""),
                        }
                    )

                for per in sii.findall(f"{{{ns}}}PeptideEvidenceRef"):
                    sii_data["peptide_evidence_refs"].append(per.get("peptideEvidence_ref", ""))

                siis.append(sii_data)

            results.append(
                {
                    "spectrum_id": spectrum_id,
                    "spectra_data_ref": spectra_data_ref,
                    "rt": rt_seconds,
                    "siis": siis,
                }
            )
        return results

    # ------------------------------------------------------------------
    # PSM record building
    # ------------------------------------------------------------------

    def _build_psm_records(self, parsed: dict) -> list[dict]:
        """Convert parsed mzIdentML data into QPX PSM record dicts."""
        spectra_data = parsed["spectra_data"]
        db_sequences = parsed["db_sequences"]
        peptides = parsed["peptides"]
        pep_evidence = parsed["peptide_evidence"]
        linker_info = parsed["linker_info"]
        spectrum_results = parsed["spectrum_results"]

        records = []
        for sir in spectrum_results:
            scan, scan_scheme = _parse_native_id(sir["spectrum_id"])
            run_file = spectra_data.get(sir["spectra_data_ref"], "unknown")
            rt = sir.get("rt")
            siis = sir["siis"]

            # Classify SIIs
            groups = _group_siis(siis)

            for group in groups:
                record = self._build_one_psm(
                    group,
                    scan,
                    run_file,
                    peptides,
                    pep_evidence,
                    db_sequences,
                    linker_info,
                    rt=rt,
                )
                if record is not None:
                    # Transient hint (dropped by the writer) telling
                    # ``_attach_spectra`` whether ``scan`` is an acquisition scan
                    # number or a 0-based file index, so it fetches by the right key.
                    record["_spectrum_id_scheme"] = scan_scheme
                    records.append(record)

        return records

    def _build_one_psm(
        self,
        group: dict,
        scan: list[int],
        run_file: str,
        peptides: dict,
        pep_evidence: dict,
        db_sequences: dict,
        linker_info: dict,
        rt: float | None = None,
    ) -> Optional[dict]:
        """Build a single QPX PSM dict from a classified SII group."""
        xl_type = group["xl_type"]
        alpha_sii = group["alpha"]
        beta_sii = group.get("beta")

        alpha_pep = peptides.get(alpha_sii["peptide_ref"], {})
        sequence = alpha_pep.get("sequence", "")
        if not sequence:
            return None

        # Protein accessions from PeptideEvidence
        proteins = self._resolve_proteins(
            alpha_sii["peptide_evidence_refs"],
            pep_evidence,
            db_sequences,
        )

        # Is decoy?
        is_decoy = any(pep_evidence.get(ref, {}).get("is_decoy", False) for ref in alpha_sii["peptide_evidence_refs"])

        # Scores
        additional_scores = _extract_scores(alpha_sii["cv_params"])

        # Extract PEP if present
        pep_value = _extract_pep(alpha_sii["cv_params"])

        # Build cross_links
        cross_links = None
        if xl_type == "intra":
            cross_links = self._build_looplink(
                alpha_pep,
                linker_info,
            )
        elif xl_type in ("inter", "noncovalent"):
            cross_links = self._build_interlink(
                alpha_pep,
                beta_sii,
                peptides,
                linker_info,
                xl_type,
            )

        # Build peptidoform from sequence and parsed modifications
        peptidoform = _build_peptidoform(sequence, alpha_pep.get("modifications"))

        # Enrich modifications with site localization scores from SII cvParams.
        # ``alpha_pep["modifications"]`` is shared by every PSM that references
        # this Peptide id, and site-localization scoring mutates positions in
        # place — deep-copy first so PSMs sharing a peptide_ref do not accumulate
        # each other's phosphoRS/Ascore values.
        modifications = _attach_site_localization_scores(
            copy.deepcopy(alpha_pep.get("modifications")),
            alpha_sii["cv_params"],
        )

        _calc_mz = alpha_sii["calc_mz"]
        _obs_mz = alpha_sii["exp_mz"]
        mass_error_ppm = 1e6 * (_obs_mz - _calc_mz) / _calc_mz if _calc_mz and _obs_mz else None

        # Extract missed_cleavages from MS:1003044 cvParam on the SII
        missed_cleavages = None
        for cv in alpha_sii["cv_params"]:
            if cv["accession"] == "MS:1003044":
                try:
                    missed_cleavages = int(cv["value"])
                except (ValueError, TypeError):
                    pass
                break

        # Build cv_params from rank, pass_threshold, and PeptideEvidence positional info
        psm_id, cv_params = _source_sii_identity(run_file, alpha_sii, beta_sii)
        rank = alpha_sii.get("rank", 1)
        pass_threshold = alpha_sii.get("pass_threshold", True)
        if rank > 1:
            cv_params.append({"cv_name": "rank", "cv_value": str(rank)})
        if not pass_threshold:
            cv_params.append({"cv_name": "pass_threshold", "cv_value": "false"})

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

        record = {
            "psm_id": psm_id,
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "charge": alpha_sii["charge"],
            "posterior_error_probability": pep_value,
            "is_decoy": is_decoy,
            "calculated_mz": _calc_mz,
            "observed_mz": _obs_mz,
            "mass_error_ppm": mass_error_ppm,
            "additional_scores": additional_scores or None,
            "predicted_rt": None,
            "run_file_name": run_file,
            "cv_params": cv_params or None,
            "scan": scan,
            "rt": rt,
            "ion_mobility": None,
            "missed_cleavages": missed_cleavages,
            "protein_accessions": proteins or None,
            "cross_links": cross_links,
            "mz_array": None,
            "intensity_array": None,
            "charge_array": None,
            "ion_type_array": None,
            "ion_mobility_array": None,
        }
        return record

    def _resolve_proteins(
        self,
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

    def _build_looplink(
        self,
        alpha_pep: dict,
        linker_info: dict,
    ) -> Optional[list[dict]]:
        """Build cross_links for a looplink (intra-peptide) PSM."""
        xl_mods = alpha_pep.get("xl_mods", {})
        if not xl_mods:
            return None

        links = []
        for xl_value, mods in xl_mods.items():
            donor = next((m for m in mods if m["role"] == "donor"), None)
            acceptor = next((m for m in mods if m["role"] == "acceptor"), None)
            if not donor or not acceptor:
                continue

            linker = linker_info.get(xl_value, linker_info.get("_default", {}))
            links.append(
                {
                    "xl_type": "intra",
                    "partner_sequence": None,
                    "partner_peptidoform": None,
                    "donor_position": donor["position"],
                    "acceptor_position": acceptor["position"],
                    "linker_name": linker.get("name", "unknown"),
                    "linker_accession": linker.get("accession"),
                    "linker_mass": donor["mass_delta"],
                }
            )

        return links or None

    def _build_interlink(
        self,
        alpha_pep: dict,
        beta_sii: Optional[dict],
        peptides: dict,
        linker_info: dict,
        xl_type: str,
    ) -> Optional[list[dict]]:
        """Build cross_links for an inter-peptide or noncovalent PSM."""
        if beta_sii is None:
            return None

        beta_pep = peptides.get(beta_sii["peptide_ref"], {})
        beta_sequence = beta_pep.get("sequence", "")

        # Find donor/acceptor positions from xl_mods
        alpha_xl = alpha_pep.get("xl_mods", {})
        beta_xl = beta_pep.get("xl_mods", {})

        # Match by shared xl_value — handle both orientations:
        # alpha=donor/beta=acceptor OR alpha=acceptor/beta=donor
        links = []
        all_xl_values = set(alpha_xl.keys()) | set(beta_xl.keys())
        for xl_value in all_xl_values:
            a_mods = alpha_xl.get(xl_value, [])
            b_mods = beta_xl.get(xl_value, [])

            donor = next((m for m in a_mods if m["role"] == "donor"), None)
            acceptor = next((m for m in b_mods if m["role"] == "acceptor"), None)

            # If donor not in alpha, check if donor is in beta (swapped orientation)
            if not donor:
                donor = next((m for m in b_mods if m["role"] == "donor"), None)
                acceptor = next((m for m in a_mods if m["role"] == "acceptor"), None)

            if not donor:
                continue

            linker = linker_info.get(xl_value, linker_info.get("_default", {}))

            links.append(
                {
                    "xl_type": xl_type,
                    "partner_sequence": beta_sequence,
                    "partner_peptidoform": beta_sequence,
                    "donor_position": donor["position"],
                    "acceptor_position": acceptor["position"] if acceptor else None,
                    "linker_name": linker.get("name", "unknown"),
                    "linker_accession": linker.get("accession"),
                    "linker_mass": donor["mass_delta"],
                }
            )

        # If no xl_mods found (e.g., noncovalent without mods), create a minimal link
        if not links and xl_type == "noncovalent":
            links.append(
                {
                    "xl_type": "noncovalent",
                    "partner_sequence": beta_sequence,
                    "partner_peptidoform": beta_sequence,
                    "donor_position": 0,
                    "acceptor_position": None,
                    "linker_name": "noncovalent",
                    "linker_accession": None,
                    "linker_mass": 0.0,
                }
            )

        return links or None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_peptidoform(sequence: str, modifications: list[dict] | None) -> str:
    """Build a ProForma string from mzIdentML parsed modifications.

    Extracts (position, accession) pairs from the modification dicts and
    delegates to ``build_proforma()``.  Accessions are normalised: if the
    accession already looks like ``UNIMOD:N`` it is kept as-is; otherwise
    ``MOD:NNNNN`` (PSI-MOD) accessions and plain names are used directly
    as the tag.
    """
    if not sequence:
        return ""
    if not modifications:
        return sequence

    mods: list[tuple[int, str]] = []
    for mod in modifications:
        accession = mod.get("accession") or mod.get("name", "")
        # Normalise UNIMOD prefix casing
        tag = re.sub(r"(?i)^unimod:", "UNIMOD:", accession) if accession else mod.get("name", "")
        for pos_info in mod.get("positions", []):
            position = pos_info.get("position", 0)
            mods.append((position, tag))

    return build_proforma(sequence, mods)


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


def _parse_scan(spectrum_id: str) -> list[int]:
    """Return the scan/index value(s) from a spectrumID, without the scheme.

    Prefer :func:`_parse_native_id` when the nativeID scheme matters (e.g. to
    attach the correct spectrum). Returns an empty list when nothing parses.
    """
    return _parse_native_id(spectrum_id)[0]


def _group_siis(siis: list[dict]) -> list[dict]:
    """Group SII elements into classified cross-link groups.

    Returns list of dicts: ``{xl_type, alpha, beta?}``
    where ``xl_type`` is one of: ``"none"``, ``"inter"``, ``"intra"``, ``"noncovalent"``.
    """
    groups = []

    # Index SIIs by their crosslink group value
    xl_groups: dict[str, list[dict]] = {}  # group_value → [sii, ...]
    looplink_siis = []
    noncovalent_groups: dict[str, list[dict]] = {}
    regular_siis = []

    for sii in siis:
        xl_group = None
        is_looplink = False
        is_noncovalent = False

        for cv in sii["cv_params"]:
            if cv["accession"] == CV_CROSSLINK_SII:
                xl_group = cv["value"]
            elif cv["accession"] == CV_LOOPLINK_SII:
                is_looplink = True
            elif cv["accession"] == CV_NONCOVALENT_SII:
                is_noncovalent = True
                xl_group = cv["value"]

        if is_looplink:
            looplink_siis.append(sii)
        elif is_noncovalent and xl_group:
            noncovalent_groups.setdefault(xl_group, []).append(sii)
        elif xl_group:
            xl_groups.setdefault(xl_group, []).append(sii)
        else:
            regular_siis.append(sii)

    # Regular PSMs
    for sii in regular_siis:
        groups.append({"xl_type": "none", "alpha": sii})

    # Looplinks (intra-peptide)
    for sii in looplink_siis:
        groups.append({"xl_type": "intra", "alpha": sii})

    # Inter-peptide crosslinks (paired SIIs)
    for _group_val, paired in xl_groups.items():
        if len(paired) >= 2:
            # First SII is alpha (donor), second is beta (acceptor)
            groups.append(
                {
                    "xl_type": "inter",
                    "alpha": paired[0],
                    "beta": paired[1],
                }
            )
        else:
            # Orphan SII — treat as regular
            groups.append({"xl_type": "none", "alpha": paired[0]})

    # Noncovalent associations
    for _group_val, paired in noncovalent_groups.items():
        if len(paired) >= 2:
            groups.append(
                {
                    "xl_type": "noncovalent",
                    "alpha": paired[0],
                    "beta": paired[1],
                }
            )
        else:
            groups.append({"xl_type": "none", "alpha": paired[0]})

    return groups


def _attach_site_localization_scores(
    modifications: list[dict] | None,
    cv_params: list[dict],
) -> list[dict] | None:
    """Enrich modification positions with site localization scores from SII cvParams.

    In mzIdentML, site localization scores (phosphoRS, ptmRS, Ascore) appear as
    cvParams on the SpectrumIdentificationItem.  Their ``value`` field typically
    encodes per-position probabilities, e.g. ``"S3: 95.3, T5: 2.1, S7: 2.6"``.
    When a single float is present it is applied to the first modification position.
    """
    if not modifications:
        return modifications

    # Collect site-localization cvParams
    site_cvs = [cv for cv in cv_params if cv["accession"] in SITE_LOCALIZATION_ACCESSIONS]
    if not site_cvs:
        return modifications

    for cv in site_cvs:
        score_name = normalize_score_name(cv["name"])
        value_str = cv.get("value", "")
        if not value_str:
            continue

        # Try to parse position-keyed format: "S3: 95.3, T5: 2.1"
        pos_scores: dict[int, float] = {}
        for token in re.split(r"[,;]\s*", value_str):
            m = re.match(r"([A-Za-z])(\d+)\s*:\s*([\d.eE+-]+)", token.strip())
            if m:
                pos = int(m.group(2))
                try:
                    pos_scores[pos] = float(m.group(3))
                except ValueError:
                    pass

        if pos_scores:
            # Attach scores to matching modification positions
            for mod in modifications:
                for pos_info in mod.get("positions", []):
                    p = pos_info.get("position", 0)
                    if p in pos_scores:
                        scores = pos_info.get("scores") or []
                        scores.append(
                            {
                                "score_name": score_name,
                                "score_value": pos_scores[p],
                                "higher_better": is_higher_better(score_name),
                            }
                        )
                        pos_info["scores"] = scores
        else:
            # Single float value — apply to first modification position
            try:
                val = float(value_str)
            except (ValueError, TypeError):
                continue
            for mod in modifications:
                positions = mod.get("positions", [])
                if positions:
                    scores = positions[0].get("scores") or []
                    scores.append(
                        {
                            "score_name": score_name,
                            "score_value": val,
                            "higher_better": is_higher_better(score_name),
                        }
                    )
                    positions[0]["scores"] = scores
                    break

    return modifications


def _extract_scores(cv_params: list[dict]) -> list[dict]:
    """Extract numeric scores from SII cvParams."""
    scores = []
    for cv in cv_params:
        acc = cv["accession"]
        if acc in _SKIP_SCORE_ACCESSIONS:
            continue
        value_str = cv.get("value", "")
        if not value_str:
            continue
        # Try to parse as float
        try:
            value = float(value_str)
        except (ValueError, TypeError):
            continue
        normalized = normalize_score_name(cv["name"])
        scores.append(
            {
                "score_name": normalized,
                "score_value": value,
                "higher_better": is_higher_better(normalized),
            }
        )
    return scores


_PEP_ACCESSIONS = frozenset({CV_PEP, CV_PEP_GLOBAL})


def _extract_pep(cv_params: list[dict]) -> float | None:
    """Extract posterior error probability from SII cvParams.

    Looks for MS:1001493 (PEP) or MS:1002352 (global PEP).
    Returns float value or None.
    """
    for cv in cv_params:
        if cv["accession"] in _PEP_ACCESSIONS:
            try:
                return float(cv.get("value", ""))
            except (ValueError, TypeError):
                pass
    return None
