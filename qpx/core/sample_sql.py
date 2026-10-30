"""SQL that resolves feature and protein-group quantities to their samples.

A quantity reaches its sample through ``run.samples[]``: the run it was measured
in (``feature.run_file_name``, or any member of ``pg.grouped_runs``) and its
channel ``label``. These common table expressions do that resolution for the
Dataset helpers and the API views.

Both lists are unnested on their own, before any join, so every join is an
equi-join on keys. Unnesting them inside one ``FROM`` list next to ``feature``
and ``run`` made DuckDB plan a cross product of features and runs, which on 405
runs took 32.5 s instead of 0.27 s and on fractionated TMT filled the temp
directory.
"""

from __future__ import annotations

from qpx.core.sql import sql_build

# A quantification unit is named by its distinct member runs, sorted and joined
# by this separator (a single-run unit keeps its run_file_name), then ``|`` and
# the channel label. The MuData proteins modality names its observations the
# same way.
UNIT_SEPARATOR = ";"

# The runs of one quantification unit are the fractions of one sample context:
# they carry the same run.samples entries (sample, label, biological and
# technical replicate) and differ only in fraction, the rule the converters use
# to build pg.grouped_runs from an SDRF (fraction_groups_from_sdrf). Technical
# replicates and separate plexes stay apart.
_RUN_UNITS_CTE = """
    run_units AS (
        SELECT run_file_name,
               array_to_string(
                   list_sort(list(run_file_name) OVER (PARTITION BY list_sort(samples))), '$separator'
               ) AS unit_runs
        FROM run
    )"""

# One row per (run, sample channel). ``n_samples`` drives the label-free
# fallback below. ``$run_extra_inner`` / ``$run_extra`` carry optional run
# columns (e.g. ``fraction``) through the unnest.
_RUN_SAMPLES_CTE = """
    run_samples AS (
        SELECT run_file_name,
               n_samples,
               unit_runs,
               s.sample_accession AS sample_accession,
               s.label AS sample_label,
               s.biological_replicate AS biological_replicate,
               s.technical_replicate AS technical_replicate$run_extra
        FROM (
            SELECT r.run_file_name, len(r.samples) AS n_samples, u.unit_runs,
                   UNNEST(r.samples) AS s$run_extra_inner
            FROM run r
            JOIN run_units u ON u.run_file_name = r.run_file_name
        )
    )"""

# The label carried by a feature intensity is not always the SDRF channel label
# in run.samples: label-free converters emit a generic feature label (DIA-NN's
# ``raw``, ``LFQ``) while the run carries ``label free sample`` or ``LFQ``. A
# label-free run has exactly one sample, so a single-sample run takes every
# intensity whatever its label; a multiplexed run needs the exact label. The
# two cases are two equi-joins, so neither multiplies rows by the channel count.
_FEATURE_SAMPLES_CTES = """
    feature_intensities AS (
        SELECT $feature_columns, i.label AS intensity_label, i.intensity AS intensity
        FROM (
            SELECT $feature_columns_inner, UNNEST(f.intensities) AS i
            FROM $source f$feature_where
        )
    ),
    feature_samples AS (
        SELECT fi.*, rs.* EXCLUDE (run_file_name)
        FROM feature_intensities fi
        JOIN run_samples rs
          ON rs.run_file_name = fi.run_file_name
         AND rs.sample_label = fi.intensity_label
        UNION ALL
        SELECT fi.*, rs.* EXCLUDE (run_file_name)
        FROM feature_intensities fi
        JOIN run_samples rs ON rs.run_file_name = fi.run_file_name
        WHERE rs.n_samples = 1
          AND (rs.sample_label = fi.intensity_label) IS NOT TRUE
    )"""

# Each distinct grouped_runs set (a quantification unit) resolved to the samples
# of its member runs. Since QPX 1.1 pg is flattened (scalar label / intensity);
# a pg row then joins its unit on (grouped_runs, label). Resolving each distinct
# set once, instead of every pg row against every member run, keeps the unnest
# to a few hundred rows. ``DISTINCT`` returns a quantity once per sample when
# several member runs are fractions of that sample. ``unit_runs`` names the set.
_PG_UNIT_SAMPLES_CTE = """
    pg_unit_samples AS (
        SELECT DISTINCT u.grouped_runs, u.unit_runs, rs.sample_accession, rs.sample_label
        FROM (
            SELECT grouped_runs,
                   array_to_string(list_sort(list_distinct(grouped_runs)), '$separator') AS unit_runs,
                   UNNEST(list_distinct(grouped_runs)) AS run_file_name
            FROM (SELECT DISTINCT grouped_runs FROM pg)
        ) u
        JOIN run_samples rs USING (run_file_name)
    )"""


def _feature_column_lists(columns: tuple[str, ...]) -> tuple[str, str]:
    """Return the outer and the ``f.``-qualified select lists of feature columns."""
    if "run_file_name" not in columns:
        columns = (*columns, "run_file_name")
    return ", ".join(columns), ", ".join(f"f.{name}" for name in columns)


def feature_sample_ctes(
    feature_columns: tuple[str, ...],
    *,
    source: str = "feature",
    feature_where: str = "",
    run_extra: tuple[str, str] = ("", ""),
) -> str:
    """Return the ``WITH`` clause that defines ``feature_samples``.

    ``feature_samples`` has one row per (feature, intensity label, sample) the
    label resolves to: the requested feature columns plus ``run_file_name``,
    ``intensity_label``, ``intensity``, ``n_samples``, ``unit_runs`` (the runs of
    the run's quantification unit), ``sample_accession``, ``sample_label``,
    ``biological_replicate``, ``technical_replicate`` and any optional run
    columns.

    Parameters
    ----------
    feature_columns : tuple of str
        Fixed feature column names to carry through (never user input).
    source : str
        The feature relation (a view name or a parenthesised subquery).
    feature_where : str
        Optional fixed ``WHERE`` predicate on the feature alias ``f``.
    run_extra : tuple of (str, str)
        Optional run columns as ``(", r.fraction", ", fraction")``: the select
        list inside the unnest and the column names outside it.
    """
    outer, inner = _feature_column_lists(feature_columns)
    run_extra_inner, run_extra_outer = run_extra
    where = f"\n            WHERE {feature_where}" if feature_where else ""
    return "WITH" + sql_build(
        _RUN_UNITS_CTE + "," + _RUN_SAMPLES_CTE + "," + _FEATURE_SAMPLES_CTES,
        separator=UNIT_SEPARATOR,
        run_extra=run_extra_outer,
        run_extra_inner=run_extra_inner,
        feature_columns=outer,
        feature_columns_inner=inner,
        source=source,
        feature_where=where,
    )


def pg_sample_ctes() -> str:
    """Return the ``WITH`` clause that defines ``pg_unit_samples``.

    ``pg_unit_samples`` maps each distinct ``grouped_runs`` set, named by
    ``unit_runs``, to the ``(sample_accession, sample_label)`` pairs of its
    member runs; join a pg row to it on ``grouped_runs`` and
    ``label = sample_label``.
    """
    return "WITH" + sql_build(
        _RUN_UNITS_CTE + "," + _RUN_SAMPLES_CTE + "," + _PG_UNIT_SAMPLES_CTE,
        separator=UNIT_SEPARATOR,
        run_extra="",
        run_extra_inner="",
    )
