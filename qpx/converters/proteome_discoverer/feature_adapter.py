"""Proteome Discoverer ``PSMs.txt`` to the feature view: the PSMs of one precursor in one run."""

from __future__ import annotations

import re

from qpx.converters.channel_labels import normalize_label
from qpx.converters.proteome_discoverer.reports import PdPsmTable, peptide_identity, split_accessions
from qpx.converters.ptm import compute_precursor_mz
from qpx.converters.utils import safe_float
from qpx.core.sql import sql_build
from qpx.writers.feature import FeatureWriter

_REPORTER_COLUMN = re.compile(r"^Abundance:? (\d{3}[NCD]?)$")
_ABUNDANCE_COLUMN = re.compile(r"^Abundance:? ")
_FEATURE_SQL = """SELECT * EXCLUDE (_best) FROM (
    SELECT *, $aggregates, row_number() OVER (g ORDER BY _qvalue NULLS LAST, _scan) AS _best
    FROM pd_psms WHERE _rank = 1 AND $used
    WINDOW g AS (PARTITION BY _run, _sequence, _modifications, _charge)
) WHERE _best = 1"""


class ProteomeDiscovererFeatureAdapter(PdPsmTable):
    """Sum the reporter abundances of a precursor's PSMs, or take its largest MS1 abundance."""

    _pepmap: dict[str, tuple[str, list[str]]] = {}
    _protein_groups: dict[str, list[str]] = {}
    _genes: dict[str, str] = {}
    _channels: dict[str, str] = {}

    def load_protein_groups(self, groups: dict[str, list[str]], genes: dict[str, str]) -> None:
        """The members and gene of each master protein, from the Proteins export."""
        self._protein_groups, self._genes = groups, genes

    def convert(self, output_path: str, *, chunksize: int = 100_000) -> bool:
        """Write the features; returns whether the export held any quantity."""
        self._pepmap = {}
        self._channels = self._reporter_columns()
        ms1 = self._resolved.get("ms1_abundance")
        if not self._channels and ms1 is None:
            self.logger.warning("This PSM export has no abundance columns, so no feature view is written")
            return False
        with FeatureWriter(output_path, creator="proteome_discoverer") as writer:
            for batch in self._query_batched(self._query(ms1), chunksize):
                records = self._transform_rows_with_skip(batch.to_pylist(), self._record, "Proteome Discoverer feature")
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)
        return True

    def write_pepmap(self, output_path: str) -> bool:
        """Write the proteins of each converted peptidoform."""
        return self._write_pepmap(output_path, self._pepmap, "proteome_discoverer")

    def _reporter_columns(self) -> dict[str, str]:
        """Reporter column -> canonical label; other abundance columns (SILAC) stop the conversion."""
        columns = sorted(self.get_table_columns("pd_psms"))
        channels = {column: normalize_label(match.group(1)) for column in columns if (match := _REPORTER_COLUMN.match(column))}
        other = [column for column in columns if _ABUNDANCE_COLUMN.match(column) and column not in channels]
        if other:
            raise ValueError(f"Proteome Discoverer MS1-label abundance columns are not converted: {', '.join(other)}")
        return channels

    def _query(self, ms1: str | None) -> str:
        sums = [
            f'SUM(TRY_CAST("{column}" AS DOUBLE)) OVER g AS "_channel_{position}"'
            for position, column in enumerate(self._channels)
        ]
        ms1_sql = f'MAX(TRY_CAST("{ms1}" AS DOUBLE)) OVER g' if ms1 else "CAST(NULL AS DOUBLE)"
        aggregates = [*sums, f"{ms1_sql} AS _ms1", "COUNT(*) OVER g AS _psm_count"]
        ambiguity = self._resolved.get("psm_ambiguity")
        used = [f"COALESCE(\"{ambiguity}\", '') <> 'Rejected'"] if ambiguity else []
        quan_info = self._resolved.get("quan_info")
        if quan_info:
            # The PSMs PD left out of its own quantification: "Excluded by Method"
            used.append(f"lower(replace(COALESCE(\"{quan_info}\", ''), ' ', '')) NOT LIKE 'excludedby%'")
        return sql_build(_FEATURE_SQL, aggregates=", ".join(aggregates), used=" AND ".join(used) or "true")

    def _quantities(self, row: dict) -> tuple[list[dict] | None, list[dict] | None]:
        ms1 = row["_ms1"] if row["_ms1"] and row["_ms1"] > 0 else None
        if not self._channels:
            if ms1 is None:
                return None, None
            additional = [{"intensity_name": "spectral_count", "intensity_value": float(row["_psm_count"])}]
            return [{"label": "LFQ", "intensity": ms1}], [{"label": "LFQ", "intensities": additional}]
        intensities = []
        for position, label in enumerate(self._channels.values()):
            value = row[f"_channel_{position}"]
            intensities.append({"label": label, "intensity": value if value and value > 0 else None})
        quantified = next((entry["label"] for entry in intensities if entry["intensity"] is not None), None)
        if quantified is None:
            return None, None
        if ms1 is None:
            return intensities, None
        return intensities, [{"label": quantified, "intensities": [{"intensity_name": "ms1_intensity", "intensity_value": ms1}]}]

    def _record(self, row: dict) -> dict | None:
        intensities, additional = self._quantities(row)
        if intensities is None:
            return None
        sequence, charge = row["_sequence"], row["_charge"]
        peptidoform, modifications = peptide_identity(sequence, row["_modifications"])
        masters = split_accessions(row.get(self._resolved.get("anchor_protein")))
        mapped = split_accessions(row.get(self._resolved.get("protein_accessions")))
        anchor = masters[0] if masters else None
        if mapped:
            self._pepmap.setdefault(peptidoform, (sequence, mapped))
        spectrum = self.spectrum_fields(row)
        apex_minutes = safe_float(row.get(self._resolved.get("apex_rt")))
        if apex_minutes is not None:
            spectrum["rt"] = apex_minutes * 60
        members = self._protein_groups.get(anchor) or masters
        gene = self._genes.get(anchor)
        interference = next(
            (s["score_value"] for s in self.search_scores(row) or [] if s["score_name"] == "parent_ion_fraction"), None
        )
        return {
            **spectrum,
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "charge": charge,
            "is_decoy": False,
            "calculated_mz": compute_precursor_mz(peptidoform, charge),
            "precursor_qvalue": row["_qvalue"],
            "additional_scores": self.search_scores(row),
            "predicted_rt": None,
            "cv_params": None,
            "intensities": intensities,
            "additional_intensities": additional,
            "pg_accessions": [
                {"accession": accession, "start": None, "end": None, "pre": None, "post": None} for accession in members
            ]
            or None,
            "anchor_protein": anchor,
            "unique": len(mapped) == 1 if mapped else None,
            "gg_names": [gene] if gene else None,
            "id_run_file_name": row["_run"],
            "identification_type": "direct",
            "parent_ion_fraction": interference,
            "ion_mobility_start": None,
            "ion_mobility_stop": None,
            "rt_start": None,
            "rt_stop": None,
        }
