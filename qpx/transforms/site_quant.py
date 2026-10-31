"""Derive site.parquet quantities from the feature view.

A modified peptide's intensity is the quantity of the site it places, per run and
label. Stoichiometry is the modified peptides' share of the intensity of all peptides
of the protein covering the residue in that run and label: an apparent value, since
peptides ionize differently.
"""

from __future__ import annotations

from pathlib import Path

import pyarrow.parquet as pq

from qpx.converters.ptm import name_to_unimod, unimod_mass, unimod_title
from qpx.core.engine import create_connection
from qpx.core.sql import escape_path, sql_build
from qpx.writers import SiteWriter

_BATCH_ROWS = 100_000
_SITE_SQL = """
    WITH occurrences AS (
        SELECT sequence, run_file_name, modifications, entry.label AS label, entry.intensity AS intensity,
               anchor_protein, unnest(list_filter(pg_positions, lambda p: p.protein_accession = anchor_protein)) AS occurrence
        FROM (SELECT *, unnest(intensities) AS entry FROM read_parquet('$feature') WHERE NOT coalesce(is_decoy, false))
        WHERE entry.intensity > 0
    ),
    placed AS (
        SELECT o.*, unnest(modifications) AS modification FROM occurrences o WHERE modifications IS NOT NULL
    ),
    positioned AS (
        SELECT anchor_protein, run_file_name, label, intensity, occurrence, modification, unnest(modification.positions) AS spot,
               length(sequence) AS peptide_length
        FROM placed
        WHERE modification.accession IN ($accessions)
    ),
    sites AS (
        SELECT anchor_protein AS protein_accession,
               CASE WHEN spot.position = 0 THEN IF(occurrence.start <= 2, 0, occurrence.start)
                    WHEN spot.position > peptide_length THEN occurrence."end"
                    ELSE occurrence.start + spot.position - 1 END AS position,
               spot.amino_acid AS amino_acid, modification.accession AS modification_accession,
               run_file_name, label, sum(intensity) AS intensity,
               max(list_max(list_transform(
                   list_filter(spot.scores, lambda s: s.score_name LIKE '%probabilit%'), lambda s: s.score_value
               ))) AS modified_probability
        FROM positioned
        GROUP BY ALL
    ),
    covering AS (
        SELECT s.protein_accession, s.position, s.modification_accession, s.run_file_name, s.label,
               sum(o.intensity) AS total
        FROM sites s
        JOIN (SELECT DISTINCT sequence, run_file_name, label, intensity, anchor_protein, occurrence FROM occurrences) o
          ON o.anchor_protein = s.protein_accession AND o.run_file_name = s.run_file_name AND o.label = s.label
         AND o.occurrence.start <= GREATEST(s.position, 1) AND o.occurrence."end" >= GREATEST(s.position, 1)
        GROUP BY ALL
    )
    SELECT s.*, least(s.intensity / c.total, 1.0) AS stoichiometry
    FROM sites s
    LEFT JOIN covering c USING (protein_accession, position, modification_accession, run_file_name, label)
    ORDER BY protein_accession, position, modification_accession, run_file_name, label
"""


def _modification_fields(modifications: list[str]) -> dict[str, dict]:
    """The requested modifications by accession; each is a Unimod name or accession."""
    fields = {}
    for modification in modifications:
        uid = name_to_unimod(modification)
        if uid is None:
            raise ValueError(f"Unknown modification: {modification}")
        fields[f"UNIMOD:{uid}"] = {"modification_name": unimod_title(uid), "mass_delta": unimod_mass(uid)}
    return fields


def derive_sites(feature_path: Path, modifications: list[str], output_path: Path) -> int:
    """Write the site view of the requested modifications; returns its rows."""
    columns = set(pq.read_schema(feature_path).names)
    needed = {"pg_positions", "modifications", "anchor_protein", "intensities"}
    if not needed <= columns or not pq.read_table(feature_path, columns=["pg_positions"]).column(0).drop_null():
        raise ValueError("The feature view has no pg_positions: run 'qpxc transform protein-properties --fasta' first")
    fields = _modification_fields(modifications)
    accessions = ", ".join(f"'{accession}'" for accession in fields)
    sql = sql_build(_SITE_SQL, feature=escape_path(str(feature_path)), accessions=accessions)
    con = create_connection()
    rows = 0
    try:
        with SiteWriter(output_path) as writer:
            for batch in con.sql(sql).to_arrow_reader(_BATCH_ROWS):
                records = [
                    {
                        **{key: value for key, value in row.items() if key != "run_file_name"},
                        **fields[row["modification_accession"]],
                        "grouped_runs": [row["run_file_name"]],
                        "is_decoy": False,
                    }
                    for row in batch.to_pylist()
                ]
                if records:
                    writer.write_batch(records)
                rows += len(records)
    finally:
        con.close()
    return rows
