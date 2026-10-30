"""SQL resolving feature intensities and pg rows to their run.samples entries."""

from __future__ import annotations

from qpx.core.sql import sql_build

# run.samples flattened to one row per run and sample, so that a sample is
# resolved by an equi-join on (run_file_name, label) instead of an UNNEST per
# joined row, which plans a product of the joined rows and every sample.
RUN_SAMPLES_SQL = """
    SELECT r.run_file_name, s.sample_accession, s.label, s.biological_replicate,
           len(r.samples) = 1 AS single_sample
    FROM run r
    CROSS JOIN UNNEST(r.samples) AS _s(s)
"""

# Feature intensities with their samples. The label of an intensity names its
# sample in the run; in a single-sample run every label resolves to that sample,
# since label-free converters write a generic label (DIA-NN's 'raw') where
# run.samples carries 'LFQ'. The exact match and the fallback are two equi-joins.
_FEATURE_SAMPLES_SQL = """
    WITH run_samples AS ($run_samples),
    feature_intensities AS (
        SELECT $columns, f.run_file_name, i.label AS intensity_label, i.intensity
        FROM $source f
        CROSS JOIN UNNEST(f.intensities) AS _i(i)
    )
    SELECT fi.* EXCLUDE (intensity_label), rs.sample_accession, rs.label, rs.biological_replicate
    FROM feature_intensities fi
    JOIN run_samples rs ON rs.run_file_name = fi.run_file_name AND rs.label = fi.intensity_label
    UNION ALL
    SELECT fi.* EXCLUDE (intensity_label), rs.sample_accession, rs.label, rs.biological_replicate
    FROM feature_intensities fi
    JOIN run_samples rs
      ON rs.run_file_name = fi.run_file_name AND rs.single_sample AND rs.label IS DISTINCT FROM fi.intensity_label
"""


def feature_samples_sql(source: str, columns: tuple[str, ...]) -> str:
    """One row per intensity of the features in *source* and its sample.

    Each row carries the feature *columns*, ``run_file_name``, ``intensity``,
    ``sample_accession``, the sample's ``label`` and ``biological_replicate``.
    """
    return sql_build(
        _FEATURE_SAMPLES_SQL,
        run_samples=RUN_SAMPLES_SQL,
        columns=", ".join(f"f.{column}" for column in columns),
        source=source,
    )
