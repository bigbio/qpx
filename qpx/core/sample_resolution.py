"""SQL resolving feature intensities and pg rows to their run.samples entries."""

from __future__ import annotations

from qpx.core.sql import sql_build

# run.samples flattened to one row per run and sample, so that a sample is
# resolved by an equi-join on (run_file_name, label) instead of an UNNEST per
# joined row, which plans a product of the joined rows and every sample.
#
# ``unit`` names the quantification unit of the run: the runs carrying the same
# samples, labels and technical replicates are the fractions of one unit (one
# multiplex mixture, or one injection series of a label-free sample), and the
# first of them by name is the unit's name. A technical replicate or another
# plex holding the same bridge sample is another unit.
RUN_SAMPLES_SQL = """
    WITH run_units AS (
        SELECT run_file_name, samples,
               min(run_file_name) OVER (
                   PARTITION BY list_sort(list_transform(
                       samples, lambda s: concat_ws('|', s.sample_accession, s.label, s.technical_replicate)))
               ) AS unit
        FROM run
    )
    SELECT r.run_file_name, r.unit, s.sample_accession, s.label, s.biological_replicate,
           len(r.samples) = 1 AS single_sample
    FROM run_units r
    CROSS JOIN UNNEST(r.samples) AS _s(s)
"""

# Feature intensities with their samples. The label of an intensity names its
# sample in the run. A feature with a single intensity in a single-sample run
# resolves to that sample whatever its label, since label-free converters write
# a generic label (DIA-NN's 'raw') where run.samples carries 'LFQ'; several
# labels in such a run (a heavy standard, pulsed SILAC) stay unresolved rather
# than adding up. The exact match and the fallback are two equi-joins. The
# intensities are unnested in the select list: a lateral CROSS JOIN UNNEST runs
# ~10x slower once a filter (is_decoy) is pushed into it.
_FEATURE_SAMPLES_SQL = """
    WITH run_samples AS ($run_samples),
    feature_intensities AS (
        SELECT * EXCLUDE (i), i.label AS intensity_label, i.intensity
        FROM (
            SELECT $columns, f.run_file_name, UNNEST(f.intensities) AS i, len(f.intensities) = 1 AS single_intensity
            FROM $source f
        )
    )
    SELECT fi.* EXCLUDE (intensity_label, single_intensity),
           rs.unit, rs.sample_accession, rs.label, rs.biological_replicate
    FROM feature_intensities fi
    JOIN run_samples rs ON rs.run_file_name = fi.run_file_name AND rs.label = fi.intensity_label
    UNION ALL
    SELECT fi.* EXCLUDE (intensity_label, single_intensity),
           rs.unit, rs.sample_accession, rs.label, rs.biological_replicate
    FROM feature_intensities fi
    JOIN run_samples rs
      ON rs.run_file_name = fi.run_file_name AND rs.single_sample AND rs.label IS DISTINCT FROM fi.intensity_label
    WHERE fi.single_intensity
"""


def feature_samples_sql(source: str, columns: tuple[str, ...]) -> str:
    """One row per intensity of the features in *source* and its sample.

    Each row carries the feature *columns*, ``run_file_name``, ``intensity``,
    the run's quantification ``unit``, ``sample_accession``, the sample's
    ``label`` and ``biological_replicate``.
    """
    return sql_build(
        _FEATURE_SAMPLES_SQL,
        run_samples=RUN_SAMPLES_SQL,
        columns=", ".join(f"f.{column}" for column in columns),
        source=source,
    )
