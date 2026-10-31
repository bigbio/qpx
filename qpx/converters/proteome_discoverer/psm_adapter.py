"""Proteome Discoverer ``PSMs.txt`` to the psm view."""

from __future__ import annotations

import re

from qpx.converters.localization import localize_modifications
from qpx.converters.proteome_discoverer.reports import PdPsmTable, peptide_identity, split_accessions
from qpx.converters.ptm import compute_precursor_mz, modification_tag
from qpx.writers.psm import PsmWriter

_SITE_PROBABILITIES = re.compile(r"^ptmRS:? (.+) Site Probabilities$")
_SITE = re.compile(r"\((\d+)\):\s*([\d.]+)")


def _ptmrs_sites(row: dict, columns: dict[str, str]) -> dict[str, dict[int, float]]:
    """Per-site probabilities of ptmRS, ``S(2): 0.1; T(4): 99.3``, by canonical tag."""
    sites = {}
    for column, tag in columns.items():
        text = row.get(column)
        if text:
            sites[tag] = {int(position): float(value) / 100 for position, value in _SITE.findall(text)}
    return sites


class ProteomeDiscovererPsmAdapter(PdPsmTable):
    """Write one psm row per run, scan, peptidoform and charge."""

    _site_columns: dict[str, str] = {}

    def convert(self, output_path: str, *, rank1_only: bool = False, chunksize: int = 100_000) -> None:
        """Write the loaded PSMs; of the hits an engine repeats on one spectrum, the best one stays."""
        columns = self.get_table_columns("pd_psms")
        self._site_columns = {
            column: modification_tag(match.group(1)) for column in columns if (match := _SITE_PROBABILITIES.match(column))
        }
        rank_filter = "WHERE _rank = 1" if rank1_only else ""
        query = f"""SELECT * FROM pd_psms {rank_filter}
            QUALIFY row_number() OVER (
                PARTITION BY _run, _scan, _sequence, _modifications, _charge ORDER BY _rank, _qvalue NULLS LAST) = 1"""
        with PsmWriter(output_path, creator="proteome_discoverer", compression=self._compression) as writer:
            for batch in self._query_batched(query, chunksize):
                records = self._transform_rows_with_skip(batch.to_pylist(), self._record, "Proteome Discoverer PSM")
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

    def _record(self, row: dict) -> dict:
        sequence, charge = row["_sequence"], row["_charge"]
        peptidoform, modifications = peptide_identity(sequence, row["_modifications"])
        modifications, candidates = localize_modifications(modifications, sequence, _ptmrs_sites(row, self._site_columns))
        confidence = row.get(self._resolved.get("confidence"))
        return {
            **self.spectrum_fields(row),
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "modification_candidates": candidates,
            "charge": charge,
            "is_decoy": False,
            "calculated_mz": compute_precursor_mz(peptidoform, charge),
            "additional_scores": self.search_scores(row),
            "qvalue": row["_qvalue"],
            "predicted_rt": None,
            "cv_params": self._cv_params(row),
            "protein_accessions": split_accessions(row.get(self._resolved.get("protein_accessions"))) or None,
            # PD's High confidence is the strict target FDR of its validator node.
            "pass_threshold": confidence == "High" if confidence else None,
            "mz_array": None,
            "intensity_array": None,
            "charge_array": None,
            "ion_type_array": None,
            "ion_mobility_array": None,
        }

    def _cv_params(self, row: dict) -> list[dict] | None:
        params = [
            {"cv_name": name, "cv_value": str(value)}
            for field, name in (("identifying_node", "identifying_node"), ("psm_ambiguity", "psm_ambiguity"))
            if (value := row.get(self._resolved.get(field)))
        ]
        return params or None
