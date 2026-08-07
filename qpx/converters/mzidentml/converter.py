"""MzIdentML full-dataset converter.

Orchestrates production of a complete QPX dataset directory from an
mzIdentML file and optional MGF spectra file.  Composes the existing
:class:`MzIdentMLPsmAdapter` for PSM extraction and adds pepmap,
provenance, ontology, and dataset metadata structures.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Any

try:
    from lxml import etree
except ImportError:
    etree: Any = None

from qpx._version import __version__
from qpx.converters.mzidentml.pg_adapter import MzIdentMLPgAdapter
from qpx.converters.mzidentml.psm_adapter import MzIdentMLPsmAdapter
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.core.scores import field_ontology_entries, score_ontology_entries
from qpx.writers.pepmap import PepMapWriter
from qpx.writers.psm import PsmWriter

logger = logging.getLogger(__name__)


class MzIdentMLConverter(BaseOrchestrator):
    """Convert mzIdentML (+ optional MGF) to a complete QPX dataset."""

    def __init__(self, compression: str = "zstd"):
        self._compression = compression

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def convert(
        self,
        mzid_path: str | Path,
        output_folder: str | Path,
        output_prefix: str = "mzidentml",
        mgf_path: str | Path | None = None,
        include_spectra: bool = False,
        project_accession: str | None = None,
    ) -> Path:
        """Convert mzIdentML to QPX dataset.

        Args:
            mzid_path: Path to ``.mzid`` or ``.mzid.gz`` file.
            output_folder: Directory for output files.
            output_prefix: File prefix (e.g. ``"PXD054720"``).
            mgf_path: Optional MGF file for spectra attachment.
            include_spectra: If ``True`` and *mgf_path* given, attach
                spectra arrays to PSM records.
            project_accession: PRIDE / ProteomeXchange accession.

        Returns:
            Path to *output_folder*.
        """
        if etree is None:
            raise ImportError("lxml is required for mzIdentML conversion. Install it with: pip install qpx[mzidentml]")
        mzid_path = Path(mzid_path)
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)

        # 1. Parse mzIdentML and build PSM records via the existing adapter
        adapter = MzIdentMLPsmAdapter()
        try:
            parsed = adapter._parse_mzidentml(str(mzid_path))
            records = adapter._build_psm_records(parsed)
        finally:
            adapter.close()

        if not records:
            logger.warning("No PSM records produced from %s", mzid_path)
            return output_folder

        # 2. Optionally enrich PSMs with spectra from MGF
        if include_spectra and mgf_path is not None:
            records = self._attach_spectra(records, mgf_path)

        # 3. Track scores for ontology generation
        adapter._track_scores(records)
        discovered_scores = adapter.get_discovered_scores()

        # 4. Write PSM parquet
        psm_path = output_folder / f"{output_prefix}.psm.parquet"
        with PsmWriter(psm_path, creator="mzidentml", compression=self._compression) as writer:
            writer.write_batch(records)
        logger.info("Wrote %d PSMs to %s", len(records), psm_path)

        # 5. Write pg (protein groups) parquet
        pg_path = output_folder / f"{output_prefix}.pg.parquet"
        with MzIdentMLPgAdapter(compression=self._compression) as pg_adapter:
            pg_adapter.convert(mzid_path=mzid_path, output_path=pg_path, creator="mzidentml")

        # 6. Write pepmap parquet
        pepmap_path = output_folder / f"{output_prefix}.pepmap.parquet"
        pepmap_records = self._build_pepmap(
            records,
            parsed["peptide_evidence"],
            parsed["peptides"],
            parsed["db_sequences"],
        )
        if pepmap_records:
            with PepMapWriter(pepmap_path, creator="mzidentml", compression=self._compression) as writer:
                writer.write_batch(pepmap_records)
            logger.info("Wrote %d pepmap entries to %s", len(pepmap_records), pepmap_path)

        # 7. Write provenance parquet
        provenance_records = self._build_provenance(mzid_path)
        self._write_provenance(output_folder, output_prefix, provenance_records)

        # 8. Write ontology parquet
        ontology_entries = score_ontology_entries(discovered_scores, view="psm")
        ontology_entries.extend(field_ontology_entries(view="psm"))
        self._write_ontology(output_folder, output_prefix, ontology_entries)

        # 9. Write dataset metadata parquet
        sw_name = None
        sw_version = None
        if provenance_records:
            sw_name = provenance_records[0].get("tool_name")
            sw_version = provenance_records[0].get("tool_version")
        self._write_dataset(
            output_folder,
            output_prefix,
            project_accession,
            software_name=sw_name or "qpx",
            software_version=sw_version or __version__,
            provenance_records=provenance_records,
        )

        return output_folder

    # ------------------------------------------------------------------
    # Spectra attachment
    # ------------------------------------------------------------------

    @staticmethod
    def _attach_spectra(records: list[dict], mgf_path: str | Path) -> list[dict]:
        """Attach mz_array and intensity_array from MGF to matching PSM records."""
        from qpx.converters.mzidentml.mgf_parser import MgfSpectraIndex

        mgf_index = MgfSpectraIndex(mgf_path)
        matched_by_scan = 0
        matched_by_index = 0
        for rec in records:
            scans = rec.get("scan") or []
            if not scans:
                continue
            value = scans[0]
            # Honour the nativeID scheme: an ``index=`` value is a 0-based file
            # position, NOT a scan number. Looking it up by scan number would
            # match an unrelated spectrum whose ``SCANS=`` equals the index and
            # attach the wrong peaks/RT. A ``scan=`` value must never be treated
            # as an index either — use the lookup that matches the scheme only.
            scheme = rec.get("_spectrum_id_scheme")
            if scheme == "index":
                spectrum = mgf_index.get_spectrum_by_index(value)
                if spectrum is not None:
                    matched_by_index += 1
            else:
                # "scan" (or a legacy record without a recorded scheme)
                spectrum = mgf_index.get_spectrum(value)
                if spectrum is not None:
                    matched_by_scan += 1

            if spectrum is not None:
                rec["mz_array"] = spectrum["mz_array"]
                rec["intensity_array"] = spectrum["intensity_array"]
                # Fill RT from MGF only if not already set from mzIdentML
                if rec.get("rt") is None and spectrum.get("rt") is not None:
                    rec["rt"] = spectrum["rt"]

        total_matched = matched_by_scan + matched_by_index
        logger.info(
            "Attached spectra to %d / %d PSMs from MGF index (%d spectra, %d by scan, %d by index)",
            total_matched,
            len(records),
            len(mgf_index),
            matched_by_scan,
            matched_by_index,
        )
        return records

    # ------------------------------------------------------------------
    # Pepmap construction
    # ------------------------------------------------------------------

    @staticmethod
    def _build_pepmap(
        records: list[dict],
        pep_evidence: dict[str, dict],
        peptides: dict[str, dict],
        db_sequences: dict[str, str],
    ) -> list[dict]:
        """Build deduplicated peptide-to-protein map grouped by peptidoform.

        Uses PSM records for peptidoform strings and PeptideEvidence XML
        data for per-protein positional info (start, end, pre, post).
        """
        # Step 1: Build (sequence, protein) → position lookup from PeptideEvidence
        position_lookup: dict[tuple[str, str], dict] = {}
        for pe in pep_evidence.values():
            pep_data = peptides.get(pe.get("peptide_ref", ""), {})
            seq = pep_data.get("sequence", "")
            acc = db_sequences.get(pe.get("db_ref", ""), "")
            if seq and acc and (seq, acc) not in position_lookup:
                position_lookup[(seq, acc)] = {
                    "start": pe.get("start"),
                    "end": pe.get("end"),
                    "pre": pe.get("pre"),
                    "post": pe.get("post"),
                }

        # Step 2: Group PSM records by peptidoform
        pepform_data: dict[str, dict] = {}
        for psm in records:
            seq = psm.get("sequence", "")
            pf = psm.get("peptidoform", seq)
            if not seq:
                continue
            if pf not in pepform_data:
                pepform_data[pf] = {"sequence": seq, "proteins": set()}
            for prot in psm.get("protein_accessions") or []:
                pepform_data[pf]["proteins"].add(prot)

        # Step 3: Build pepmap records with pg_accessions structs
        pepmap_records: list[dict] = []
        for pf, data in pepform_data.items():
            seq = data["sequence"]
            proteins = data["proteins"]
            pg_accessions = []
            for prot in sorted(proteins):
                pos = position_lookup.get((seq, prot), {})
                pg_accessions.append(
                    {
                        "accession": prot,
                        "start": pos.get("start"),
                        "end": pos.get("end"),
                        "pre": pos.get("pre"),
                        "post": pos.get("post"),
                    }
                )
            pepmap_records.append(
                {
                    "sequence": seq,
                    "peptidoform": pf,
                    "pg_accessions": pg_accessions or None,
                    "is_unique": len(proteins) == 1,
                }
            )

        return pepmap_records

    # ------------------------------------------------------------------
    # Provenance extraction
    # ------------------------------------------------------------------

    @staticmethod
    def _build_provenance(mzid_path: str | Path) -> list[dict]:
        """Extract processing provenance from mzIdentML AnalysisSoftware elements."""
        mzid_path = Path(mzid_path)

        if str(mzid_path).endswith(".gz"):
            with gzip.open(mzid_path, "rb") as fh:
                tree = etree.parse(fh)
        else:
            tree = etree.parse(str(mzid_path))

        root = tree.getroot()
        ns = root.nsmap.get(None, "")

        provenance_records: list[dict] = []
        step_order = 0

        # Collect parameters from AnalysisProtocolCollection (enzyme + modifications)
        params_by_sip: dict[str, list[dict]] = {}
        for sip in root.iter(f"{{{ns}}}SpectrumIdentificationProtocol"):
            sip_id = sip.get("id", "")
            params: list[dict] = []

            # Enzyme(s)
            for enzyme in sip.iter(f"{{{ns}}}Enzyme"):
                enzyme_name_el = enzyme.find(f"{{{ns}}}EnzymeName")
                if enzyme_name_el is not None:
                    cv = enzyme_name_el.find(f"{{{ns}}}cvParam")
                    name = cv.get("name", "") if cv is not None else ""
                else:
                    name = enzyme.get("name", "")
                if name:
                    params.append({"key": "enzyme", "value": name})

            # Modifications (fixed and variable)
            for mod in sip.iter(f"{{{ns}}}SearchModification"):
                fixed = mod.get("fixedMod", "false").lower() == "true"
                residues = mod.get("residues", "").strip()
                cv = mod.find(f"{{{ns}}}cvParam")
                mod_name = cv.get("name", "") if cv is not None else ""
                if mod_name:
                    label = "fixed_mod" if fixed else "variable_mod"
                    value = f"{mod_name} ({residues})" if residues else mod_name
                    params.append({"key": label, "value": value})

            if params:
                params_by_sip[sip_id] = params

        # Build a mapping from AnalysisSoftware id → SIP parameters
        # (SpectrumIdentificationProtocol references analysisSoftware_ref)
        sw_params: dict[str, list[dict]] = {}
        for sip in root.iter(f"{{{ns}}}SpectrumIdentificationProtocol"):
            sw_ref = sip.get("analysisSoftware_ref", "")
            sip_id = sip.get("id", "")
            if sw_ref and sip_id in params_by_sip:
                sw_params[sw_ref] = params_by_sip[sip_id]

        for sw in root.iter(f"{{{ns}}}AnalysisSoftware"):
            step_order += 1

            name_attr = sw.get("name", "")
            version = sw.get("version")
            sw_id = sw.get("id", "")

            # Try to get the proper software name from SoftwareName/cvParam
            tool_name = name_attr
            software_name_el = sw.find(f"{{{ns}}}SoftwareName")
            if software_name_el is not None:
                cv = software_name_el.find(f"{{{ns}}}cvParam")
                if cv is not None:
                    tool_name = cv.get("name", name_attr)

            parameters = sw_params.get(sw_id) or None

            provenance_records.append(
                {
                    "step_order": step_order,
                    "step_category": "database_search",
                    "step_name": "spectrum_identification",
                    "tool_name": tool_name,
                    "tool_version": version,
                    "tool_uri": sw.get("uri"),
                    "parameters": parameters,
                    "config": None,
                    "output_views": ["psm"],
                }
            )

        return provenance_records
