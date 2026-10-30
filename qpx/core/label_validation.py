"""Cross-view run and label checks: the runs a view names exist, and each quantity label resolves to one sample."""

from __future__ import annotations

import duckdb

from qpx.core.data.schema import ValidationIssue
from qpx.core.sql import sql_build, validate_identifier

# Examples quoted in an issue that covers many (run, label) pairs.
_EXAMPLES = 5

# Each protein-group unit (its runs and label) is resolved once, not once per pg
# row, so fractionated multi-plex pg stays small (bigbio/qpx#H4).
GROUPED_RUN_SAMPLE_MAPPING_SQL = """
                WITH units AS (
                    -- pg is flattened since 1.1: scalar label; ID-only rows (null
                    -- label) have no sample to resolve.
                    SELECT list_sort(list_distinct(grouped_runs)) AS grouped_runs,
                           label,
                           count(*) AS pg_rows
                    FROM pg
                    WHERE label IS NOT NULL
                    GROUP BY ALL
                ),
                sample_matches AS (
                    SELECT DISTINCT u.grouped_runs, u.label, rs.sample_accession
                    FROM units u
                    CROSS JOIN UNNEST(u.grouped_runs) AS _g(run_file_name)
                    JOIN run r USING (run_file_name)
                    CROSS JOIN UNNEST(r.samples) AS _s(rs)
                    WHERE rs.label = u.label
                )
                SELECT u.grouped_runs,
                       u.label,
                       any_value(u.pg_rows) AS pg_rows,
                       COUNT(sm.sample_accession) AS sample_count,
                       LIST(sm.sample_accession ORDER BY sm.sample_accession)
                           FILTER (WHERE sm.sample_accession IS NOT NULL) AS samples
                FROM units u
                LEFT JOIN sample_matches sm
                  ON sm.grouped_runs = u.grouped_runs AND sm.label = u.label
                GROUP BY u.grouped_runs, u.label
                HAVING COUNT(sm.sample_accession) != 1
                ORDER BY u.label, u.grouped_runs
                """

_DANGLING_RUN_SQL = """
                SELECT v.run_file_name, count(*) AS row_count
                FROM $view v
                WHERE v.run_file_name IS NOT NULL
                  AND NOT EXISTS (SELECT 1 FROM run r WHERE r.run_file_name = v.run_file_name)
                GROUP BY v.run_file_name
                ORDER BY v.run_file_name
                """

_DUPLICATE_LABEL_SQL = """
                SELECT count(*)
                FROM feature
                WHERE len(intensities) <> len(list_distinct(list_transform(intensities, entry -> entry.label)))
                """

# Each (run, label) of the feature intensities that does not match exactly one
# sample, with the run's sample count; a run missing from run.parquet is the
# dangling-run check's.
_FEATURE_LABEL_SAMPLE_SQL = """
                WITH labels AS (
                    SELECT DISTINCT f.run_file_name, entry.label
                    FROM feature f, UNNEST(f.intensities) AS _i(entry)
                ),
                samples AS (
                    SELECT r.run_file_name, rs.label, rs.sample_accession
                    FROM run r, UNNEST(r.samples) AS _s(rs)
                )
                SELECT l.run_file_name,
                       l.label,
                       COUNT(DISTINCT s.sample_accession) AS sample_count,
                       max(len(r.samples)) AS run_sample_count
                FROM labels l
                JOIN run r USING (run_file_name)
                LEFT JOIN samples s ON s.run_file_name = l.run_file_name AND s.label = l.label
                GROUP BY l.run_file_name, l.label
                HAVING COUNT(DISTINCT s.sample_accession) != 1
                ORDER BY l.run_file_name, l.label
                """

# Multi-run units whose grouped_runs are not stored in ascending numeric
# fraction order; units with a run that has no numeric fraction are skipped.
_FRACTION_ORDER_SQL = """
                WITH members AS (
                    SELECT grouped_runs,
                           unnest(grouped_runs) AS run_file_name,
                           unnest(range(len(grouped_runs))) AS position
                    FROM (SELECT DISTINCT grouped_runs FROM pg WHERE len(grouped_runs) > 1)
                ),
                fractions AS (
                    SELECT m.grouped_runs,
                           list(TRY_CAST(r.fraction AS DOUBLE) ORDER BY m.position) AS fractions
                    FROM members m
                    JOIN run r USING (run_file_name)
                    GROUP BY m.grouped_runs
                )
                SELECT count(*)
                FROM fractions
                WHERE list_count(fractions) = len(fractions)
                  AND len(fractions) > 1
                  AND fractions <> list_sort(fractions)
                """

_PSM_LABEL_SQL = """
                SELECT count(*)
                FROM psm p
                JOIN feature f ON f.feature_id = p.feature_id
                WHERE p.label IS NOT NULL
                  AND f.identified_labels IS NOT NULL
                  AND NOT list_contains(f.identified_labels, p.label)
                """


def skipped_check_reason(exc: duckdb.Error) -> str:
    """Why a cross-structure check could not run, for its skipped issue.

    An I/O failure says nothing about the data: DuckDB spills large queries to
    its temp directory (``.tmp`` under the working directory by default), so a
    read-only working directory fails the query (bigbio/qpx#300). Running out of
    memory is a resource limit, not a data problem either.
    """
    if isinstance(exc, duckdb.OutOfMemoryException):
        return "out of memory; raise DuckDB's memory limit or give it a temp directory with room to spill"
    if isinstance(exc, duckdb.IOException):
        return "I/O error; if DuckDB could not create its temp directory, run from a writable directory"
    return "possible schema drift"


def _issue(structure: str, check: str, severity: str, column: str | None, message: str) -> ValidationIssue:
    return ValidationIssue(structure=structure, check=check, severity=severity, column=column, message=message)


def _examples(pairs: list[tuple]) -> str:
    shown = ", ".join(f"({run!r}, {label!r})" for run, label in pairs[:_EXAMPLES])
    return shown + (", ..." if len(pairs) > _EXAMPLES else "")


def dangling_run_issues(engine, view: str) -> list[ValidationIssue]:
    """Runs that *view* rows name but run.parquet lacks; those rows resolve to no sample."""
    rows = engine.execute(sql_build(_DANGLING_RUN_SQL, view=validate_identifier(view))).fetchall()
    return [
        _issue(
            view,
            "dangling_run",
            "error",
            "run_file_name",
            f"{count} {view} row(s) name run {run!r}, which is not a run_file_name in run.parquet",
        )
        for run, count in rows
    ]


def feature_label_issues(engine) -> list[ValidationIssue]:
    """Repeated labels in one feature row, and (run, label) pairs that do not resolve to exactly one sample.

    A label that matches no sample of a run with a single sample is a warning:
    the quantity queries fall back to that sample.
    """
    issues = []
    duplicates = engine.execute(_DUPLICATE_LABEL_SQL).fetchall()[0][0]
    if duplicates:
        message = f"{duplicates} feature row(s) repeat a label in intensities; each label must appear once per row"
        issues.append(_issue("feature", "duplicate_label", "error", "intensities", message))
    rows = engine.execute(_FEATURE_LABEL_SAMPLE_SQL).fetchall()
    fallback = [(run, label) for run, label, count, run_count in rows if count == 0 and run_count == 1]
    unresolved = [(run, label) for run, label, count, run_count in rows if count == 0 and run_count != 1]
    ambiguous = [(run, label) for run, label, count, _ in rows if count > 1]
    if fallback:
        message = (
            f"{len(fallback)} (run, label) pair(s) match no run.samples label and resolve only because "
            f"the run has a single sample: {_examples(fallback)}"
        )
        issues.append(_issue("feature", "label_single_sample_fallback", "warning", "intensities", message))
    if unresolved:
        message = f"{len(unresolved)} (run, label) pair(s) resolve to no sample in run.samples: {_examples(unresolved)}"
        issues.append(_issue("feature", "unresolved_label", "error", "intensities", message))
    if ambiguous:
        message = f"{len(ambiguous)} (run, label) pair(s) resolve to several samples in run.samples: {_examples(ambiguous)}"
        issues.append(_issue("feature", "ambiguous_label_sample", "error", "intensities", message))
    return issues


def psm_label_issues(engine, severity: str) -> list[ValidationIssue]:
    """PSMs whose label is not among the identified labels of the feature they link to."""
    psm_columns = {row[0] for row in engine.execute('DESCRIBE "psm"').fetchall()}
    feature_columns = {row[0] for row in engine.execute('DESCRIBE "feature"').fetchall()}
    if not {"label", "feature_id"} <= psm_columns or "identified_labels" not in feature_columns:
        return []
    count = engine.execute(_PSM_LABEL_SQL).fetchall()[0][0]
    if not count:
        return []
    message = f"{count} psm row(s) have a label that is not in identified_labels of the feature they link to"
    return [_issue("psm", "label_not_identified", severity, "label", message)]


def fraction_order_issues(engine) -> list[ValidationIssue]:
    """Protein-group units whose grouped_runs are not in fraction order."""
    if "fraction" not in {row[0] for row in engine.execute('DESCRIBE "run"').fetchall()}:
        return []
    count = engine.execute(_FRACTION_ORDER_SQL).fetchall()[0][0]
    if not count:
        return []
    message = f"{count} grouped_runs list(s) are not stored in ascending run.fraction order"
    return [_issue("pg", "grouped_runs_fraction_order", "warning", "grouped_runs", message)]


def grouped_run_sample_issues(engine) -> list[ValidationIssue]:
    """Protein-group units (grouped runs and label) that do not resolve to exactly one sample."""
    return [
        _issue(
            "pg",
            "ambiguous_grouped_run_sample",
            "error",
            "grouped_runs",
            f"{pg_rows} pg row(s) of label {label!r} across grouped_runs={grouped_runs!r} resolve to "
            f"{sample_count} samples ({samples or []!r}); exactly one sample is required",
        )
        for grouped_runs, label, pg_rows, sample_count, samples in engine.execute(GROUPED_RUN_SAMPLE_MAPPING_SQL).fetchall()
    ]
