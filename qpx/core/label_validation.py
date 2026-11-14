"""Cross-view checks: the runs and samples a view names exist, each quantity label resolves to one sample, and links resolve."""

from __future__ import annotations

import logging

import duckdb

from qpx.core.data.schema import ValidationIssue
from qpx.core.data.validation import _issue
from qpx.core.sql import sql_build

_log = logging.getLogger("qpx.dataset")

# The distinct pg.grouped_runs elements that run.parquet does not list. pg is
# reduced to its distinct runs before the anti-join: a lateral UNNEST per pg row
# with a correlated lookup each ran out of memory on a fractionated TMT pg.
_DANGLING_GROUPED_RUN_SQL = """
                WITH pg_runs AS (
                    SELECT DISTINCT grouped_run
                    FROM (SELECT UNNEST(grouped_runs) AS grouped_run FROM pg)
                    WHERE grouped_run IS NOT NULL
                )
                SELECT g.grouped_run
                FROM pg_runs g
                WHERE NOT EXISTS (SELECT 1 FROM run r WHERE r.run_file_name = g.grouped_run)
                ORDER BY g.grouped_run
                """

# run.samples as one row per (run, channel). It is unnested in the select list of
# its own subquery: an UNNEST(r.samples) placed after a join to run is evaluated
# for every joined row, and under DuckDB 2.0 that made the pg sample-mapping check
# spill gigabytes on a fractionated TMT dataset of a few MB.
_RUN_CHANNELS_SQL = """
                    SELECT run_file_name, sc.label AS label, sc.sample_accession AS sample_accession
                    FROM (SELECT run_file_name, UNNEST(samples) AS sc FROM run)
"""

# pg is flattened since 1.1 (scalar label). Every protein group of a plex repeats
# the same (grouped runs, label) unit, so each distinct unit is resolved once
# instead of once per pg row. ID-only rows (null label) have no sample to resolve.
_GROUPED_RUN_SAMPLE_MAPPING_SQL = sql_build(
    """
                WITH units AS MATERIALIZED (
                    SELECT ROW_NUMBER() OVER () AS unit_id, *
                    FROM (
                        SELECT list_sort(list_distinct(grouped_runs)) AS unit_runs,
                               label,
                               COUNT(*) AS pg_rows,
                               MIN(anchor_protein) AS example_protein
                        FROM pg
                        WHERE label IS NOT NULL
                        GROUP BY ALL
                    )
                ),
                unit_members AS (
                    SELECT unit_id, label, UNNEST(unit_runs) AS run_file_name
                    FROM units
                ),
                run_channels AS ($run_channels),
                unit_samples AS (
                    SELECT um.unit_id, rc.sample_accession
                    FROM unit_members um
                    JOIN run_channels rc
                      ON rc.run_file_name = um.run_file_name
                     AND rc.label = um.label
                )
                SELECT u.example_protein,
                       u.unit_runs,
                       u.label,
                       u.pg_rows,
                       COUNT(DISTINCT us.sample_accession) AS sample_count,
                       list_sort(LIST(DISTINCT us.sample_accession)
                           FILTER (WHERE us.sample_accession IS NOT NULL)) AS samples
                FROM units u
                LEFT JOIN unit_samples us USING (unit_id)
                GROUP BY u.unit_id, u.example_protein, u.unit_runs, u.label, u.pg_rows
                HAVING COUNT(DISTINCT us.sample_accession) != 1
                ORDER BY u.example_protein, u.label, u.unit_runs
                """,
    run_channels=_RUN_CHANNELS_SQL,
)


# Each distinct (run, label) of the feature intensities, with the number of
# intensities that carry it, resolved once against the run's channels. A pair
# must match exactly one channel entry: the sample-joined views join a quantity
# to every entry it matches, so a sample channel the run lists twice counts it
# twice. Runs that are not in the run view are left to the run_file_name check.
_FEATURE_LABEL_SAMPLE_SQL = sql_build(
    """
                WITH units AS MATERIALIZED (
                    SELECT ROW_NUMBER() OVER () AS unit_id, *
                    FROM (
                        SELECT run_file_name, i.label AS label, COUNT(*) AS n_values
                        FROM (SELECT run_file_name, UNNEST(intensities) AS i FROM feature)
                        GROUP BY ALL
                    )
                    WHERE run_file_name IN (SELECT run_file_name FROM run)
                ),
                run_channels AS ($run_channels),
                unit_samples AS (
                    SELECT u.unit_id, rc.sample_accession
                    FROM units u
                    JOIN run_channels rc
                      ON rc.run_file_name = u.run_file_name
                     AND rc.label = u.label
                )
                SELECT u.run_file_name,
                       u.label,
                       u.n_values,
                       COUNT(DISTINCT us.sample_accession) AS sample_count,
                       list_sort(LIST(DISTINCT us.sample_accession)
                           FILTER (WHERE us.sample_accession IS NOT NULL)) AS samples
                FROM units u
                LEFT JOIN unit_samples us USING (unit_id)
                GROUP BY u.unit_id, u.run_file_name, u.label, u.n_values
                HAVING COUNT(us.sample_accession) != 1
                ORDER BY u.run_file_name, u.label
                """,
    run_channels=_RUN_CHANNELS_SQL,
)

# Feature rows that list a (non-null) label more than once in intensities.
_DUPLICATE_FEATURE_LABEL_SQL = """
                SELECT COUNT(*)
                FROM (SELECT list_transform(intensities, lambda x: x.label) AS labels FROM feature)
                WHERE len(list_distinct(labels)) < len(list_filter(labels, lambda x: x IS NOT NULL))
                """

# Views whose rows name the acquisition run they were measured in.
_RUN_KEYED_VIEWS = ("feature", "psm", "mz")

# Per view: the distinct run_file_name values that run.parquet does not list,
# with the number of rows that carry each. The view is reduced to its distinct
# values before the anti-join, so memory follows the number of runs.
_DANGLING_RUN_FILE_NAME_SQL = {
    view: sql_build(
        """
                WITH view_runs AS (
                    SELECT run_file_name, COUNT(*) AS n_rows
                    FROM $view
                    WHERE run_file_name IS NOT NULL
                    GROUP BY run_file_name
                )
                SELECT v.run_file_name, v.n_rows
                FROM view_runs v
                WHERE NOT EXISTS (SELECT 1 FROM run r WHERE r.run_file_name = v.run_file_name)
                ORDER BY v.run_file_name
                """,
        view=view,
    )
    for view in _RUN_KEYED_VIEWS
}

# How many example values an aggregated issue message lists.
_MAX_EXAMPLES = 5


def _examples(items: list[str]) -> str:
    """The first few of the formatted *items*, for an issue message that counts them all."""
    shown = ", ".join(items[:_MAX_EXAMPLES])
    return f"{shown}, ..." if len(items) > _MAX_EXAMPLES else shown


def _feature_label_issues(duplicated: int, units: list[tuple]) -> list[ValidationIssue]:
    """Issues for feature labels that list a sample twice or name no single sample.

    *units* are the (run_file_name, label, n_values, sample_count, samples)
    rows of the (run, label) pairs that match no channel entry of their run or
    several: of several samples, or of one sample the run lists more than once.
    """
    issues = []
    if duplicated:
        message = (
            f"{duplicated} feature row(s) list a label more than once in intensities; a label names one sample, "
            "which would receive several quantities from one row"
        )
        issues.append(_issue("feature", "duplicate_intensity_label", "error", "intensities", message))
    unresolved = [unit for unit in units if unit[3] == 0]
    if unresolved:
        issues.append(_label_pair_issue(unresolved, "match no sample_channel label"))
    ambiguous = [unit for unit in units if unit[3] > 1]
    if ambiguous:
        issues.append(_label_pair_issue(ambiguous, "match the labels of several samples"))
    repeated = [unit for unit in units if unit[3] == 1]
    if repeated:
        issues.append(_label_pair_issue(repeated, "match a sample channel listed more than once in the samples"))
    return issues


def _label_pair_issue(units: list[tuple], what: str) -> ValidationIssue:
    """One ambiguous_label_sample issue for the (run, label, n_values, sample_count, samples) *units*."""
    examples = [f"({run!r}, {label!r})" + (f" -> {samples!r}" if count else "") for run, label, _, count, samples in units]
    message = (
        f"{sum(unit[2] for unit in units)} feature intensities in {len(units)} (run_file_name, label) pair(s) "
        f"{what} of their run ({_examples(examples)}); each label must name exactly one sample in run.samples"
    )
    return _issue("feature", "ambiguous_label_sample", "error", "intensities", message)


def _view_columns(engine, view: str) -> set[str]:
    """Column names of a registered view."""
    return {row[0] for row in engine.execute(sql_build('DESCRIBE "$view"', view=view)).fetchall()}


def _skipped_issue(structure: str, check: str, what: str, exc: duckdb.Error, *, strict: bool) -> ValidationIssue:
    """The issue for a cross-structure check that could not run: an error under strict validation."""
    reason = _skipped_check_reason(exc)
    _log.warning("%s could not run (%s): %s", what, reason, exc)
    return _issue(structure, check, "error" if strict else "warning", None, f"{what} could not run ({reason}): {exc}")


def _skipped_check_reason(exc: duckdb.Error) -> str:
    """Why a cross-structure check could not run, for its skipped issue.

    Only a binder or catalog error (a view or column the check reads is missing
    or has another type) points at the data. DuckDB spills large queries to its
    temp directory (``.tmp`` under the working directory by default): a
    read-only working directory fails the query with an I/O error
    (bigbio/qpx#300), and a query that outgrows the memory limit and the temp
    space fails with an out-of-memory error.
    """
    if isinstance(exc, duckdb.OutOfMemoryException):
        return (
            "DuckDB ran out of memory or temp space, which says nothing about the data; "
            "raise the memory limit (Dataset(duckdb_memory=...)) or free temp space"
        )
    if isinstance(exc, duckdb.IOException):
        return "I/O error; if DuckDB could not create its temp directory, run from a writable directory"
    if isinstance(exc, (duckdb.BinderException, duckdb.CatalogException)):
        return "possible schema drift"
    return "query error"


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

# Below this share of linked features, feature->pg queries silently lose rows.
_MIN_PG_LINK_COVERAGE = 0.95


_DANGLING_FEATURE_IDS_SQL = """
                    SELECT DISTINCT p.feature_id AS feature_id
                    FROM psm p
                    WHERE p.feature_id IS NOT NULL
                      AND NOT EXISTS (SELECT 1 FROM feature f WHERE f.feature_id = p.feature_id)
                    ORDER BY feature_id
                    """

_DANGLING_PSM_IDS_SQL = """
                    SELECT DISTINCT pid AS psm_id
                    FROM feature f, UNNEST(f.psm_ids) AS _u(pid)
                    WHERE pid IS NOT NULL
                      AND NOT EXISTS (SELECT 1 FROM psm p WHERE p.psm_id = pid)
                    ORDER BY psm_id
                    """

_FEATURE_PSM_DESYNC_SQL = """
                    SELECT DISTINCT p.psm_id AS psm_id, p.feature_id AS feature_id
                    FROM psm p JOIN feature f ON f.feature_id = p.feature_id
                    WHERE p.feature_id IS NOT NULL
                      AND f.psm_ids IS NOT NULL AND len(f.psm_ids) > 0
                      AND NOT list_contains(f.psm_ids, p.psm_id)
                    ORDER BY psm_id
                    """

_FEATURE_PSM_INVERSE_DESYNC_SQL = """
                    SELECT DISTINCT f.feature_id, p.psm_id, p.feature_id AS actual_feature_id
                    FROM feature f
                    CROSS JOIN UNNEST(f.psm_ids) AS _u(referenced_psm_id)
                    JOIN psm p ON p.psm_id = referenced_psm_id
                    WHERE p.feature_id IS NOT NULL
                      AND p.feature_id <> f.feature_id
                    ORDER BY f.feature_id, p.psm_id
                    """


def _feature_psm_references(engine):
    """Query dangling and contradictory feature/PSM references."""
    feature_columns = {row[0] for row in engine.execute('DESCRIBE "feature"').fetchall()}
    psm_columns = {row[0] for row in engine.execute('DESCRIBE "psm"').fetchall()}
    has_psm_feature_id = "feature_id" in psm_columns
    # feature.psm_ids is an OPTIONAL producer hardlink (qpx computes the
    # inverse on read instead). Only run its referential/desync checks when
    # a producer actually populated it — an absent/all-null column cannot
    # desync from the authoritative psm.feature_id.
    has_feature_psm_ids = "psm_ids" in feature_columns
    feature_psm_ids_populated = False
    if has_feature_psm_ids:
        (feature_psm_ids_populated,) = engine.execute(
            "SELECT EXISTS (SELECT 1 FROM feature WHERE psm_ids IS NOT NULL AND len(psm_ids) > 0)"
        ).fetchall()[0]

    dangling_feature_ids = []
    if has_psm_feature_id:
        dangling_feature_ids = engine.execute(_DANGLING_FEATURE_IDS_SQL).fetchall()

    dangling_psm_ids = []
    if feature_psm_ids_populated:
        dangling_psm_ids = engine.execute(_DANGLING_PSM_IDS_SQL).fetchall()

    # Reciprocal desync: a psm points at a feature whose (non-empty) psm_ids
    # does not list it back — i.e. the two directions describe different edges.
    desync = []
    inverse_desync = []
    if feature_psm_ids_populated and has_psm_feature_id:
        desync = engine.execute(_FEATURE_PSM_DESYNC_SQL).fetchall()
        # The inverse contradiction: a feature lists a psm that explicitly
        # points at a different feature. A null psm.feature_id means that
        # direction is unpopulated and is therefore not a contradiction.
        inverse_desync = engine.execute(_FEATURE_PSM_INVERSE_DESYNC_SQL).fetchall()
    return dangling_feature_ids, dangling_psm_ids, desync, inverse_desync


def psm_reference_issues(dangling_feature_ids, desync, severity: str) -> list[ValidationIssue]:
    """Build psm reference issues in reporting order."""
    issues = []
    for (feature_id,) in dangling_feature_ids:
        issues.append(
            ValidationIssue(
                structure="psm",
                check="dangling_feature_id",
                severity=severity,
                column="feature_id",
                message=(f"psm.feature_id {feature_id!r} does not resolve to a feature.feature_id in feature.parquet"),
            )
        )
    for psm_id, feature_id in desync:
        issues.append(
            ValidationIssue(
                structure="psm",
                check="feature_psm_desync",
                severity=severity,
                column="feature_id",
                message=(
                    f"psm.psm_id {psm_id!r} points to feature {feature_id!r}, but that "
                    f"feature.psm_ids does not list this psm back (reciprocal cross-ref desync)"
                ),
            )
        )
    return issues


def _feature_reference_issues(dangling_psm_ids, inverse_desync, severity: str) -> list[ValidationIssue]:
    """Build feature reference issues in reporting order."""
    issues = []
    for (psm_id,) in dangling_psm_ids:
        issues.append(
            ValidationIssue(
                structure="feature",
                check="dangling_psm_id",
                severity=severity,
                column="psm_ids",
                message=(f"feature.psm_ids contains {psm_id!r}, which does not resolve to a psm.psm_id in psm.parquet"),
            )
        )
    for feature_id, psm_id, actual_feature_id in inverse_desync:
        issues.append(
            ValidationIssue(
                structure="feature",
                check="feature_psm_desync",
                severity=severity,
                column="psm_ids",
                message=(
                    f"feature.feature_id {feature_id!r} lists psm {psm_id!r}, but that "
                    f"psm.feature_id points to feature {actual_feature_id!r} "
                    "(reciprocal cross-ref desync)"
                ),
            )
        )
    return issues


def pg_link_coverage_issues(diagnostics: dict) -> list[ValidationIssue]:
    """Warn when few features link to a pg row, which usually means two producers named groups differently."""
    linked, unlinked = diagnostics["n_features_linked"], diagnostics["n_features_without_pg"]
    if not linked and not unlinked or linked >= _MIN_PG_LINK_COVERAGE * (linked + unlinked):
        return []
    share = linked / (linked + unlinked)
    message = f"Only {share:.1%} of features link to a pg row by membership, run and label ({unlinked} do not)"
    return [_issue("feature", "pg_link_coverage", "warning", "pg_accessions", message)]


def run_sample_issues(engine) -> list[ValidationIssue]:
    """Samples run.samples names that the sample view lacks; their quantities resolve to no metadata."""
    rows = engine.execute(
        """
        SELECT DISTINCT entry.sample_accession
        FROM run, UNNEST(run.samples) AS t(entry)
        WHERE entry.sample_accession NOT IN (SELECT sample_accession FROM sample)
        ORDER BY 1
        """
    ).fetchall()
    return [
        _issue("run", "dangling_sample", "error", "samples", f"run.samples names {sample!r}, which is not in the sample view")
        for (sample,) in rows
    ]


def fragment_link_issues(engine, severity: str) -> list[ValidationIssue]:
    """Fragments that name a feature the feature view does not have."""
    (count,) = engine.execute(
        'SELECT count(*) FROM "fragment" v WHERE NOT EXISTS (SELECT 1 FROM feature f WHERE f.feature_id = v.feature_id)'
    ).fetchall()[0]
    if not count:
        return []
    message = f"{count} fragment row(s) have a feature_id that is not in the feature view"
    return [_issue("fragment", "dangling_feature_id", severity, "feature_id", message)]


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
