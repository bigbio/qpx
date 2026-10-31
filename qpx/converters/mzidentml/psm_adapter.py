"""mzIdentML PSM adapter — parses mzIdentML XML (including 1.3 XL-MS) to QPX PSMs.

Uses lxml for XML parsing because pyopenms cannot handle mzIdentML 1.3
crosslinking modifications (``crosslink donor``, ``crosslink acceptor``).

//TODO: Use pyopenms in the futue when solving the issues with cross-linking modifications, as it provides a more robust and tested mzIdentML parser. For now, this custom parser is necessary to handle the specific features of mzIdentML 1.3 used in XL-MS data.
"""

from __future__ import annotations

import gzip
import logging
from typing import Any

try:
    from lxml import etree
except ImportError:
    etree: Any = None

from qpx.converters.base import BaseConverter
from qpx.converters.mzidentml.psm_records import _build_psm_records
from qpx.converters.ptm import build_proforma, canonical_mod_tag
from qpx.core.cv_terms import CV_SCAN_START_TIME, CV_XL_ACCEPTOR, CV_XL_DONOR
from qpx.core.files import run_file_stem
from qpx.writers.psm import PsmWriter

logger = logging.getLogger(__name__)


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
        records = _build_psm_records(parsed)

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
        mod_name = ""
        mod_accession = None
        for cv in mod.findall(f"{{{ns}}}cvParam"):
            mod_name = cv.get("name", "")
            mod_accession = cv.get("accession")
            break  # Take first cvParam as the mod identity
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
            rt_seconds = self._spectrum_rt(sir, ns)

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
