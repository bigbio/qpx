"""mzIdentML PG adapter — converts ProteinDetectionList to pg.parquet.

Parses ``ProteinAmbiguityGroup`` elements from a mzIdentML file and
emits one ``pg.parquet`` row per protein group, using ``PgWriter``.
"""

from __future__ import annotations

import logging
from pathlib import Path

from qpx.converters.base import BaseConverter
from qpx.converters.mzidentml.document import MzidDocument, stream_results
from qpx.converters.utils import canonical_accession, uniprot_entry_name
from qpx.core.files import run_file_stem
from qpx.core.scores import normalize_score_name
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)

# CV accession constants
_CV_GROUP_REPRESENTATIVE = "MS:1002403"
_CV_LEADING_PROTEIN = "MS:1002401"
_CV_PROTEIN_GROUP_QVALUE = "MS:1002373"

_DECOY_PREFIXES = ("rev_", "DECOY_", "decoy_")


class MzIdentMLPgAdapter(BaseConverter):
    """Convert mzIdentML ProteinDetectionList to QPX pg.parquet."""

    def convert(
        self,
        mzid_path: str | Path,
        output_path: str | Path,
        creator: str = "mzidentml",
        document: MzidDocument | None = None,
    ) -> None:
        """Parse *mzid_path* and write QPX PG records to *output_path*.

        Args:
            mzid_path: Path to ``.mzid`` or ``.mzid.gz`` file.
            output_path: Destination ``pg.parquet`` path.
            creator: Creator tag written to Parquet metadata.
            document: The file's tree when its results were already streamed.
        """
        mzid_path = Path(mzid_path)
        fallback_run = mzid_path.name.replace(".mzid.gz", "").replace(".mzid", "")

        if document is None:
            document = MzidDocument()
            for _ in stream_results(mzid_path, document):
                document.release_sequences()
        root, ns = document.root, document.ns
        grouped_runs_by_list, fallback_runs = self._build_grouped_runs(root, ns, fallback_run, document.spectra_refs)

        records = list(self._iter_protein_groups(root, ns, document.db_accessions, grouped_runs_by_list, fallback_runs))

        if not records:
            logger.warning("No protein group records produced from %s", mzid_path)
            return

        with PgWriter(output_path, creator=creator) as writer:
            writer.write_batch(records)

        logger.info("Wrote %d protein groups to %s", len(records), output_path)

    # ------------------------------------------------------------------
    # XML parsing helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _ordered_unique(values: list[str]) -> list[str]:
        """Return non-empty values once, preserving document order."""
        return list(dict.fromkeys(value for value in values if value))

    @classmethod
    def _build_grouped_runs(
        cls,
        root,
        ns: str,
        fallback_run: str,
        spectra_refs: dict[str, list[str]],
    ) -> tuple[dict[str, list[str]], list[str]]:
        """Resolve each ProteinDetectionList to its contributing spectra runs."""
        spectra_data = {
            element.get("id"): run_file_stem(element.get("location", ""))
            for element in root.iter(f"{{{ns}}}SpectraData")
            if element.get("id")
        }
        fallback_runs = cls._ordered_unique(list(spectra_data.values())) or [fallback_run]

        runs_by_identification_list: dict[str, list[str]] = {}
        for identification in root.iter(f"{{{ns}}}SpectrumIdentification"):
            list_ref = identification.get("spectrumIdentificationList_ref")
            if not list_ref:
                continue
            runs = list(runs_by_identification_list.get(list_ref, []))
            runs.extend(
                spectra_data.get(input_spectra.get("spectraData_ref"), "")
                for input_spectra in identification.findall(f"{{{ns}}}InputSpectra")
            )
            runs_by_identification_list[list_ref] = cls._ordered_unique(runs)

        # Some producers omit AnalysisCollection/InputSpectra.  The result-level
        # spectraData_ref is the authoritative fallback for those files.
        for identification_list in root.iter(f"{{{ns}}}SpectrumIdentificationList"):
            list_id = identification_list.get("id")
            if not list_id:
                continue
            runs = list(runs_by_identification_list.get(list_id, []))
            runs.extend(spectra_data.get(ref, "") for ref in spectra_refs.get(list_id, []))
            runs_by_identification_list[list_id] = cls._ordered_unique(runs)

        runs_by_protein_list: dict[str, list[str]] = {}
        for detection in root.iter(f"{{{ns}}}ProteinDetection"):
            protein_list_ref = detection.get("proteinDetectionList_ref")
            if not protein_list_ref:
                continue
            runs: list[str] = []
            for input_ids in detection.findall(f"{{{ns}}}InputSpectrumIdentifications"):
                runs.extend(runs_by_identification_list.get(input_ids.get("spectrumIdentificationList_ref", ""), []))
            resolved = cls._ordered_unique(runs)
            if resolved:
                runs_by_protein_list[protein_list_ref] = resolved

        return runs_by_protein_list, fallback_runs

    def _iter_protein_groups(
        self,
        root,
        ns: str,
        dbseq_map: dict[str, str],
        grouped_runs_by_list: dict[str, list[str]],
        fallback_runs: list[str],
    ):
        """Yield one pg record dict per ProteinAmbiguityGroup."""
        for protein_list in root.iter(f"{{{ns}}}ProteinDetectionList"):
            grouped_runs = grouped_runs_by_list.get(protein_list.get("id", ""), fallback_runs)
            for pag in protein_list.findall(f"{{{ns}}}ProteinAmbiguityGroup"):
                record = self._build_pag_record(pag, ns, dbseq_map, grouped_runs)
                if record is not None:
                    yield record

    @staticmethod
    def _float_value(value: str | None) -> float | None:
        try:
            return float(value) if value is not None else None
        except (TypeError, ValueError):
            return None

    @classmethod
    def _extract_pag_scores(cls, pag, ns: str) -> tuple[float | None, list[dict]]:
        """Separate group-level q-value from other PAG scores."""
        global_qvalue: float | None = None
        additional_scores: list[dict] = []
        for cv in pag.findall(f"{{{ns}}}cvParam"):
            cv_acc = cv.get("accession", "")
            score_val = cls._float_value(cv.get("value"))
            if score_val is None:
                continue
            if cv_acc == _CV_PROTEIN_GROUP_QVALUE:
                global_qvalue = score_val
                continue
            additional_scores.append(
                {
                    "score_name": normalize_score_name(cv.get("name", "") or cv_acc),
                    "score_value": score_val,
                    "higher_better": None,
                }
            )
        return global_qvalue, additional_scores

    def _pag_members(self, pdhs, ns: str, dbseq_map: dict, additional_scores: list) -> tuple:
        """Collect member accessions, peptide counts, representative and scores."""
        pg_accessions: list[str] = []
        member_names: list[str | None] = []
        anchor_protein: str | None = None
        peptides: list[dict] = []
        for pdh in pdhs:
            source_id = dbseq_map.get(pdh.get("dBSequence_ref", ""), pdh.get("dBSequence_ref", ""))
            acc = canonical_accession(source_id)
            if acc:
                pg_accessions.append(acc)
                member_names.append(uniprot_entry_name(source_id))

            # Count PeptideHypothesis elements as peptide count for this protein
            pep_count = len(pdh.findall(f"{{{ns}}}PeptideHypothesis"))
            if acc:
                peptides.append({"protein_name": acc, "peptide_count": pep_count})

            # Extract cvParams
            for cv in pdh.findall(f"{{{ns}}}cvParam"):
                cv_acc = cv.get("accession", "")
                cv_name = cv.get("name", "")
                cv_val = cv.get("value")

                if cv_acc in (_CV_GROUP_REPRESENTATIVE, _CV_LEADING_PROTEIN):
                    if anchor_protein is None:
                        anchor_protein = acc
                else:
                    score_val = self._float_value(cv_val)
                    if score_val is not None:
                        additional_scores.append(
                            {
                                "score_name": normalize_score_name(cv_name or cv_acc),
                                "score_value": score_val,
                                "higher_better": None,
                            }
                        )

        return pg_accessions, member_names if all(member_names) else None, anchor_protein, peptides

    def _build_pag_record(self, pag, ns: str, dbseq_map: dict[str, str], grouped_runs: list[str]) -> dict | None:
        """Build a single pg record from a ProteinAmbiguityGroup element."""
        pdhs = pag.findall(f"{{{ns}}}ProteinDetectionHypothesis")
        if not pdhs:
            return None

        global_qvalue, additional_scores = self._extract_pag_scores(pag, ns)

        pg_accessions, pg_names, anchor_protein, peptides = self._pag_members(pdhs, ns, dbseq_map, additional_scores)

        if not pg_accessions:
            return None

        # Fallback anchor: first accession
        if anchor_protein is None:
            anchor_protein = pg_accessions[0]

        is_decoy = any(acc.startswith(_DECOY_PREFIXES) for acc in pg_accessions)

        return {
            "pg_accessions": pg_accessions,
            "pg_names": pg_names,
            "gg_accessions": None,
            "gg_names": None,
            "gg_qvalue": None,
            "anchor_protein": anchor_protein,
            "grouped_runs": list(grouped_runs),
            "global_qvalue": global_qvalue,
            "pg_qvalue": None,
            "intensities": None,
            "additional_intensities": None,
            "is_decoy": is_decoy,
            "contaminant": None,
            "peptides": peptides,
            "peptide_counts": None,
            "feature_counts": None,
            "sequence_coverage": None,
            "molecular_weight": None,
            "additional_scores": additional_scores if additional_scores else None,
            "cv_params": None,
        }
