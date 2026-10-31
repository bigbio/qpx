"""mzIdentML PSM adapter — parses mzIdentML XML (including 1.3 XL-MS) to QPX PSMs.

Uses lxml for XML parsing because pyopenms cannot handle mzIdentML 1.3
crosslinking modifications (``crosslink donor``, ``crosslink acceptor``).

//TODO: Use pyopenms in the futue when solving the issues with cross-linking modifications, as it provides a more robust and tested mzIdentML parser. For now, this custom parser is necessary to handle the specific features of mzIdentML 1.3 used in XL-MS data.
"""

from __future__ import annotations

import logging
from itertools import chain
from typing import Iterator

from qpx.converters.base import BaseConverter
from qpx.converters.mzidentml.document import MzidDocument, stream_results
from qpx.converters.mzidentml.psm_records import _build_psm_records, declare_scan_format, scan_schemes
from qpx.converters.ptm import build_proforma, canonical_mod_tag
from qpx.core.cv_terms import (
    CV_MODIFICATION_INDEX,
    CV_SCAN_START_TIME,
    CV_XL_ACCEPTOR,
    CV_XL_DONOR,
    MODIFICATION_NON_IDENTITY_ACCESSIONS,
)
from qpx.core.files import run_file_stem
from qpx.writers.psm import PsmWriter

logger = logging.getLogger(__name__)

# Results of a file converted to records at a time.
_RESULT_BATCH = 20_000


def modification_identity_param(element, ns: str):
    """The cvParam that names a Modification or SearchModification; ``None`` when none does."""
    for cv in element.findall(f"{{{ns}}}cvParam"):
        if cv.get("accession") not in MODIFICATION_NON_IDENTITY_ACCESSIONS:
            return cv
    return None


class MzIdentMLPsmAdapter(BaseConverter):
    """Convert mzIdentML files (including XL-MS 1.3) to QPX PSM Parquet."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        # Sections of the last file that its records refer to; set with its first result.
        self.parsed: dict = {}

    def convert(
        self,
        mzid_path: str,
        output_path: str,
        creator: str = "mzidentml",
    ) -> None:
        """Parse *mzid_path* and write QPX PSM records to *output_path*, declaring the scan format."""
        batches = (records for records in self.iter_psm_batches(mzid_path, MzidDocument()) if records)
        first = next(batches, None)
        if first is None:
            logger.warning("No PSM records produced from %s", mzid_path)
            return
        total = 0
        with_links = 0
        formats: set = set()
        with PsmWriter(output_path, creator=creator, compression=self._compression) as writer:
            for records in chain([first], batches):
                self._track_scores(records)
                formats |= scan_schemes(records)
                writer.write_batch(records)
                total += len(records)
                with_links += sum(1 for r in records if r.get("cross_links"))
            declare_scan_format(writer, formats)
        logger.info("Wrote %d PSMs (%d with cross-links) to %s", total, with_links, output_path)

    # ------------------------------------------------------------------
    # XML parsing
    # ------------------------------------------------------------------

    def iter_psm_batches(self, path, document: MzidDocument, rank1_only: bool = False) -> Iterator[list[dict]]:
        """Yield QPX PSM records in batches, parsing *path* once; *document* keeps the rest of the file."""
        results: list[dict] = []
        for result in stream_results(path, document):
            if not self.parsed:
                self.parsed = self._parse_sections(document.root, document.ns)
                document.release_sequences()
            results.append(self._parse_spectrum_result(result, document.ns))
            if len(results) == _RESULT_BATCH:
                yield _build_psm_records({**self.parsed, "spectrum_results": results}, rank1_only)
                results = []
        if results:
            yield _build_psm_records({**self.parsed, "spectrum_results": results}, rank1_only)

    def _parse_sections(self, root, ns: str) -> dict:
        """Parse the sections that come before the first result into intermediate data structures."""
        return {
            "spectra_data": self._parse_spectra_data(root, ns),
            "db_sequences": self._parse_db_sequences(root, ns),
            "peptides": self._parse_peptides(root, ns),
            "peptide_evidence": self._parse_peptide_evidence(root, ns),
            "linker_info": self._parse_linker_info(root, ns),
            "modification_tags": self._parse_modification_tags(root, ns),
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

            mods, xl_mods, tags = self._parse_peptide_modifications(pep, ns, sequence or "")

            peptides[pep_id] = {
                "sequence": sequence,
                "peptidoform": build_proforma(sequence, tags) if sequence else "",
                "modifications": mods or None,
                "xl_mods": xl_mods,
            }
        return peptides

    @staticmethod
    def _regular_modification(mod, ns: str, location: int, mass_text, sequence: str) -> tuple[dict, str]:
        """The QPX modification and canonical ProForma tag of a non-cross-link Modification."""
        identity = modification_identity_param(mod, ns)
        mod_name = identity.get("name", "") if identity is not None else ""
        mod_accession = identity.get("accession") if identity is not None else None
        residue = sequence[location - 1] if 0 < location <= len(sequence) else None
        mass = float(mass_text) if mass_text else None
        terminus = "N" if location == 0 else "C" if location > len(sequence) else None
        tag = canonical_mod_tag(mod_accession, mod_name, mass, residue, terminus)
        head = tag.split("|", 1)[0]
        is_unimod = head.startswith("UNIMOD:")
        return {
            "name": mod_name if is_unimod and mod_name else head,
            "accession": head if is_unimod else None,
            "positions": [{"position": location, "amino_acid": mod.get("residues", "") or None, "scores": None}],
        }, tag

    @staticmethod
    def _parse_peptide_modifications(pep, ns: str, sequence: str) -> tuple[list, dict, list[tuple[int, str]]]:
        """Separate regular modifications from crosslink donor/acceptor sites.

        Also returns the ``(position, tag)`` pairs of the regular modifications
        in canonical ProForma form, checked against the reported mass.
        """
        mods = []
        tags: list[tuple[int, str]] = []
        xl_mods: dict[str, list[dict]] = {}  # keyed by xl pair value

        for mod in pep.findall(f"{{{ns}}}Modification"):
            location = int(mod.get("location", "0"))
            mass_text = mod.get("monoisotopicMassDelta")
            mass_delta = float(mass_text) if mass_text else 0.0
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
                mod_dict, tag = MzIdentMLPsmAdapter._regular_modification(mod, ns, location, mass_text, sequence)
                mods.append(mod_dict)
                tags.append((location, tag))

        return mods, xl_mods, tags

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

    @staticmethod
    def _parse_modification_tags(root, ns: str) -> dict[int, str]:
        """Canonical tag of each SearchModification with a ``modification index``, the key of MS:1002505 scores."""
        tags: dict[int, str] = {}
        for search_mod in root.iter(f"{{{ns}}}SearchModification"):
            index = search_mod.find(f"{{{ns}}}cvParam[@accession='{CV_MODIFICATION_INDEX}']")
            identity = modification_identity_param(search_mod, ns)
            if index is None or identity is None:
                continue
            residue = search_mod.get("residues", ".")[0]
            tag = canonical_mod_tag(
                identity.get("accession"),
                identity.get("name", ""),
                float(search_mod.get("massDelta", "0")),
                None if residue == "." else residue,
            )
            tags[int(index.get("value"))] = tag.split("|", 1)[0]
        return tags

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

    def _parse_spectrum_result(self, sir, ns: str) -> dict:
        """Parse one SpectrumIdentificationResult element."""
        siis = []
        for sii in sir.findall(f"{{{ns}}}SpectrumIdentificationItem"):
            sii_data = {
                "id": sii.get("id", ""),
                "charge": int(sii.get("chargeState", "0")),
                "exp_mz": float(sii.get("experimentalMassToCharge", "0")),
                # Optional in mzIdentML; cross-link SIIs commonly leave it out.
                "calc_mz": float(sii.get("calculatedMassToCharge")) if sii.get("calculatedMassToCharge") else None,
                "rank": int(sii.get("rank", "1")),
                "pass_threshold": sii.get("passThreshold", "false").lower() == "true",
                "peptide_ref": sii.get("peptide_ref", ""),
                "cv_params": [
                    {"accession": cv.get("accession", ""), "name": cv.get("name", ""), "value": cv.get("value", "")}
                    for cv in sii.findall(f"{{{ns}}}cvParam")
                ],
                "peptide_evidence_refs": [
                    per.get("peptideEvidence_ref", "") for per in sii.findall(f"{{{ns}}}PeptideEvidenceRef")
                ],
            }
            siis.append(sii_data)
        return {
            "spectrum_id": sir.get("spectrumID", ""),
            "spectra_data_ref": sir.get("spectraData_ref", ""),
            "rt": self._spectrum_rt(sir, ns),
            "siis": siis,
        }

    @staticmethod
    def _spectrum_rt(sir, ns: str) -> float | None:
        """Read scan start time in seconds from the result cvParams."""
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

        return rt_seconds
