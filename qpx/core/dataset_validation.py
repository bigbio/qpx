"""Dataset schema, cross-structure and file integrity checks."""

from __future__ import annotations

import logging
from pathlib import Path

import duckdb

from qpx.core.data.schema import ValidationIssue, ValidationResult
from qpx.core.data.validation import _LINK_CONTEXT_MESSAGE, _LINK_MEMBERSHIP_MESSAGE, _LINK_SQL, _issue
from qpx.core.files import sha256_file
from qpx.core.sql import sql_build

_log = logging.getLogger("qpx.dataset")

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
# intensities that carry it, resolved once against the run's channels. Runs that
# are not in the run view are left to the run_file_name check.
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
                HAVING COUNT(DISTINCT us.sample_accession) != 1
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
    rows of the (run, label) pairs that resolve to no sample or to several.
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
    return issues


def _label_pair_issue(units: list[tuple], what: str) -> ValidationIssue:
    """One ambiguous_label_sample issue for the (run, label, n_values, sample_count, samples) *units*."""
    examples = [f"({run!r}, {label!r})" + (f" -> {samples!r}" if count > 1 else "") for run, label, _, count, samples in units]
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


def _psm_reference_issues(dangling_feature_ids, desync, severity: str) -> list[ValidationIssue]:
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


def _verify_one_file(
    fpath: Path,
    name: str,
    expected_sha: str,
    expected_size: int | None,
    expected_rows: int | None,
) -> tuple[list[str], list[str]]:
    """Check one recorded file's size, checksum and row count.

    Returns ``(errors, warnings)``.
    """
    errors: list[str] = []
    warnings: list[str] = []

    if expected_size is not None:
        actual_size = fpath.stat().st_size
        if actual_size != expected_size:
            # Cheap and decisive: a size change is a content change, and
            # reporting it by name beats a bare checksum mismatch.
            errors.append(f"Size mismatch: {name} ({actual_size} != {expected_size})")

    if sha256_file(fpath) != expected_sha:
        errors.append(f"Checksum mismatch: {name}")

    if expected_rows is None:
        return errors, warnings
    if expected_rows < 0:
        # The -1 sentinel compute_integrity stores for a file whose metadata it
        # could not read. Never silently treat it as a verified count.
        warnings.append(f"Row count was not recorded for {name}; it could not be read at packaging time")
        return errors, warnings

    actual_rows = _row_count(fpath)
    if actual_rows is None:
        warnings.append(f"Could not read a row count for {name}")
    elif actual_rows != expected_rows:
        errors.append(f"Row count mismatch: {name} ({actual_rows} != {expected_rows})")
    return errors, warnings


def _as_int(value) -> int | None:
    """Coerce a stored integrity number to int, or None if it is not one.

    Integrity dicts come back through pandas, so the values are often numpy
    integers and ``isinstance(value, int)`` is False for them — type-testing
    here silently skipped every comparison.
    """
    if value is None or isinstance(value, bool):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _row_count(path: Path) -> int | None:
    """Row count of a Parquet or AnnData file, or None when it cannot be read."""
    try:
        if path.suffix == ".h5ad":
            import anndata

            return int(anndata.read_h5ad(path, backed="r").n_obs)
        import pyarrow.parquet as pq

        return int(pq.read_metadata(path).num_rows)
    # pylint: disable-next=broad-exception-caught
    except Exception as exc:  # noqa: BLE001 - verification must not raise on a bad file
        _log.debug("Could not read a row count for %s: %s", path, exc)
        return None


def _unrecorded_file_warnings(path: Path, dataset_meta_path: Path, verified: set[Path]) -> list[str]:
    """Report data files missing from the stored integrity record."""
    warnings = []
    # Files present but absent from the record. Verification that only walks
    # the stored keys can answer "are the recorded files intact?" but never
    # "is anything here that should not be?", so an added file passed clean.
    for actual in sorted(path.rglob("*.parquet")) + sorted(path.rglob("*.h5ad")):
        resolved = actual.resolve()
        if resolved == dataset_meta_path or resolved in verified:
            continue
        warnings.append(f"File not covered by the integrity record: {actual.relative_to(path).as_posix()}")
    return warnings


class _DatasetValidationMixin:
    """Validate registered structures, cross-references and file integrity."""

    # --- Validation ---
    def validate(self, structures: list[str] | None = None, *, strict: bool = False) -> dict[str, ValidationResult]:
        """Validate the dataset or specific structures against their schemas.

        Parameters
        ----------
        structures : list[str] | None
            Structure names to validate.  If *None*, validates all available.

        Returns
        -------
        dict[str, ValidationResult]
            Mapping of structure name to its validation result.
        """
        targets = structures or self.available_structures
        results = {name: self._validate_structure(name, strict=strict) for name in targets}
        present = {name for name in results if self._structures.get(name) is not None}
        self._check_run_resolution(results, present, strict=strict)
        self._check_cross_references(results, present, strict=strict)
        return results

    def _validate_structure(self, name: str, *, strict: bool) -> ValidationResult:
        """Validate one registered structure, or report it missing."""
        struct = self._structures.get(name)
        if struct is not None:
            return struct.validate(strict=strict)
        result = ValidationResult(structure=name)
        result.issues.append(
            ValidationIssue(
                structure=name,
                check="missing_structure",
                severity="error",
                column=None,
                message=f"Structure '{name}' not found in dataset at {self.path}",
            )
        )
        return result

    def _check_run_resolution(self, results: dict[str, ValidationResult], present: set[str], *, strict: bool) -> None:
        """Cross-structure checks that resolve rows through the run view.

        A row that names no acquisition run, or a label that names no single
        sample of its runs, is dropped (or counted twice) by every sample-joined
        view, so these findings are errors.
        """
        if self.run is None:
            return
        # Every pg.grouped_runs element must be a real run.run_file_name, and each
        # pg label must name exactly one sample across its grouped runs.
        if "pg" in present:
            self._check_grouped_runs_referential(results["pg"], strict=strict)
            self._check_grouped_runs_sample_mapping(results["pg"], strict=strict)
        # feature, psm and mz rows each name the run they were measured in, and
        # each feature label must name exactly one sample of that run.
        for view in _RUN_KEYED_VIEWS:
            if view in present:
                self._check_run_file_names(view, results[view], strict=strict)
        if "feature" in present:
            self._check_feature_label_samples(results["feature"], strict=strict)

    def _check_cross_references(self, results: dict[str, ValidationResult], present: set[str], *, strict: bool) -> None:
        """Cross-structure checks of the id links between views."""
        # Cross-structure invariant: the optional feature<->psm cross-references
        # must resolve. Every non-null psm.feature_id must name a real
        # feature.feature_id, and every feature.psm_ids element a real psm.psm_id.
        # Unresolved references are warnings during normal use and errors under
        # strict validation, but only when both views are present.
        if {"feature", "psm"} <= present:
            self._check_feature_psm_referential(results["feature"], results["psm"], strict=strict)

        # feature.pg_ids is the canonical feature -> protein-group reference.
        # When both views are present, every populated id must resolve to pg.pg_id.
        if {"feature", "pg"} <= present:
            self._check_feature_pg_referential(results["feature"], strict=strict)

        # site/proteoform links: pg_id names one pg row of the row's run (and, for a
        # proteoform, label) context; best_psm_id names a psm row.
        for view in ("site", "proteoform"):
            if view in present:
                self._check_pg_link(view, results[view], strict=strict)

    def _check_run_file_names(self, view: str, result: ValidationResult, *, strict: bool) -> None:
        """Flag ``run_file_name`` values of *view* that are not in run.run_file_name.

        Each distinct value is looked up once. The check is skipped when the view
        has no ``run_file_name`` column: it is optional in mz, and a missing
        required column is already reported by the schema check.
        """
        try:
            if "run_file_name" not in _view_columns(self._engine, view):
                return
            rows = self._engine.execute(_DANGLING_RUN_FILE_NAME_SQL[view]).fetchall()
        except duckdb.Error as exc:
            what = f"{view} run_file_name referential check"
            result.issues.append(_skipped_issue(view, "referential_check_skipped", what, exc, strict=strict))
            return
        if not rows:
            return
        result.issues.append(
            _issue(
                view,
                "dangling_run_file_name",
                "error",
                "run_file_name",
                f"{sum(n for _, n in rows)} {view} row(s) name {len(rows)} run_file_name value(s) that are not "
                f"in run.parquet ({_examples([repr(name) for name, _ in rows])}); these rows resolve to no run or sample",
            )
        )

    def _check_feature_label_samples(self, result: ValidationResult, *, strict: bool) -> None:
        """Require each feature label to name exactly one sample of its run, once per row.

        ``run.samples[].label`` is the canonical channel name and must match
        ``feature.intensities[].label`` exactly (see the data-model spec). Each
        distinct (run_file_name, label) pair is resolved once; rows of a run that
        is not in the run view are left to :meth:`_check_run_file_names`.
        """
        try:
            duplicated = self._engine.execute(_DUPLICATE_FEATURE_LABEL_SQL).fetchall()[0][0]
            units = self._engine.execute(_FEATURE_LABEL_SAMPLE_SQL).fetchall()
        except duckdb.Error as exc:
            what = "feature label sample-mapping check"
            result.issues.append(_skipped_issue("feature", "sample_mapping_check_skipped", what, exc, strict=strict))
            return
        result.issues.extend(_feature_label_issues(duplicated, units))

    def _check_pg_link(self, view: str, result: ValidationResult, *, strict: bool = False) -> None:
        """Flag site/proteoform links that do not resolve, or resolve to another context.

        ``pg_id`` must name a pg row: for a proteoform the row of its own label and
        grouped runs, for a site (no label) the identification row (null label) of
        its grouped runs. ``best_psm_id`` must name a psm row when the psm view is
        present.
        The pg row must hold the site's protein or the proteoform's family in ``pg_accessions``.
        """
        severity = "error" if strict else "warning"
        checks = []
        if self.pg is not None:
            checks.append(
                ("dangling_pg_id", "pg_id", _LINK_SQL[view]["dangling_pg"], "have a pg_id that does not resolve to a pg.pg_id")
            )
            checks.append(("pg_link_context", "pg_id", _LINK_SQL[view]["pg_context"], _LINK_CONTEXT_MESSAGE[view]))
            checks.append(("pg_link_membership", "pg_id", _LINK_SQL[view]["pg_membership"], _LINK_MEMBERSHIP_MESSAGE[view]))
        if self.psm is not None:
            checks.append(
                (
                    "dangling_psm_id",
                    "best_psm_id",
                    _LINK_SQL[view]["dangling_psm"],
                    "have a best_psm_id that does not resolve to a psm.psm_id",
                )
            )
        for check, column, sql, what in checks:
            try:
                count = self._engine.execute(sql).fetchall()[0][0]
            except duckdb.Error as exc:
                reason = _skipped_check_reason(exc)
                _log.warning("%s link check could not run (%s): %s", view, reason, exc)
                result.issues.append(
                    ValidationIssue(
                        structure=view,
                        check="referential_check_skipped",
                        severity=severity,
                        column=column,
                        message=f"{view} link check could not run ({reason}): {exc}",
                    )
                )
                continue
            if count:
                result.issues.append(
                    ValidationIssue(
                        structure=view, check=check, severity=severity, column=column, message=f"{count} {view} row(s) {what}"
                    )
                )

    def _check_feature_pg_referential(self, feature_result: ValidationResult, *, strict: bool = False) -> None:
        """Flag feature.pg_ids elements that do not resolve to pg.pg_id."""
        severity = "error" if strict else "warning"
        try:
            feature_columns = {row[0] for row in self._engine.execute('DESCRIBE "feature"').fetchall()}
            if "pg_ids" not in feature_columns:
                return
            dangling_pg_ids = self._engine.execute(
                """
                SELECT DISTINCT referenced_pg_id
                FROM feature f, UNNEST(f.pg_ids) AS _u(referenced_pg_id)
                WHERE referenced_pg_id IS NOT NULL
                  AND NOT EXISTS (SELECT 1 FROM pg p WHERE p.pg_id = referenced_pg_id)
                ORDER BY referenced_pg_id
                """
            ).fetchall()
        except duckdb.Error as exc:
            reason = _skipped_check_reason(exc)
            _log.warning("feature->pg referential check could not run (%s): %s", reason, exc)
            feature_result.issues.append(
                ValidationIssue(
                    structure="feature",
                    check="referential_check_skipped",
                    severity=severity,
                    column="pg_ids",
                    message=f"feature->pg referential check could not run ({reason}): {exc}",
                )
            )
            return

        for (pg_id,) in dangling_pg_ids:
            feature_result.issues.append(
                ValidationIssue(
                    structure="feature",
                    check="dangling_pg_id",
                    severity=severity,
                    column="pg_ids",
                    message=f"feature.pg_ids contains {pg_id!r}, which does not resolve to a pg.pg_id in pg.parquet",
                )
            )

    def _check_feature_psm_referential(
        self, feature_result: ValidationResult, psm_result: ValidationResult, *, strict: bool = False
    ) -> None:
        """Flag feature<->psm cross-references that do not resolve.

        ``psm.feature_id`` is the single source of truth for the feature<->psm
        association (bigbio/qpx#267); qpx does not materialize its inverse,
        ``feature.psm_ids`` (computed on read via
        :meth:`~qpx.dataset.Dataset.link_feature_psm`). This check therefore always
        validates the authoritative direction — appending an issue to *psm_result*
        for each non-null ``psm.feature_id`` that is not a ``feature.feature_id``.
        The ``feature.psm_ids`` checks (dangling ``psm_ids`` elements and
        **reciprocal desync**) run ONLY when a producer actually supplies a
        populated ``feature.psm_ids`` (an optional producer hardlink); when the
        column is absent or all-null they are skipped — there is nothing to desync,
        since the inverse is computed, not stored. Reciprocal desync is also only
        flagged where the opposite direction is populated, so an unpopulated
        cross-ref is never a false positive. A query failure (e.g. the optional
        cross-ref columns are absent on older files) is surfaced as
        ``referential_check_skipped`` rather than masking schema validation.
        Issues are warnings normally and errors under strict validation.
        """
        severity = "error" if strict else "warning"
        try:
            dangling_feature_ids, dangling_psm_ids, desync, inverse_desync = _feature_psm_references(self._engine)
        except duckdb.Error as exc:
            reason = _skipped_check_reason(exc)
            _log.warning("feature<->psm referential check could not run (%s): %s", reason, exc)
            skipped = ValidationIssue(
                structure="psm",
                check="referential_check_skipped",
                severity="error" if strict else "warning",
                column=None,
                message=f"feature<->psm referential check could not run ({reason}): {exc}",
            )
            psm_result.issues.append(skipped)
            return

        psm_result.issues.extend(_psm_reference_issues(dangling_feature_ids, desync, severity))
        feature_result.issues.extend(_feature_reference_issues(dangling_psm_ids, inverse_desync, severity))

    def _check_grouped_runs_referential(self, pg_result: ValidationResult, *, strict: bool = False) -> None:
        """Flag pg.grouped_runs values that are not present in run.run_file_name.

        Appends an ``error`` issue to *pg_result* for each dangling grouped-run
        token; each distinct token is looked up once. A query failure (e.g. a
        missing column) does not mask schema validation: it is reported as
        ``referential_check_skipped``.
        """
        try:
            rows = self._engine.execute(_DANGLING_GROUPED_RUN_SQL).fetchall()
        except duckdb.Error as exc:
            # Do not fail silently: a query error here can mean genuine schema
            # drift (a referenced pg/run column was removed/renamed), which is how
            # the last release's bug hid. Surface it in the validation result.
            reason = _skipped_check_reason(exc)
            _log.warning("pg grouped_runs referential check could not run (%s): %s", reason, exc)
            pg_result.issues.append(
                ValidationIssue(
                    structure="pg",
                    check="referential_check_skipped",
                    severity="error" if strict else "warning",
                    column=None,
                    message=f"grouped_runs referential check could not run ({reason}): {exc}",
                )
            )
            return

        for (grouped_run,) in rows:
            pg_result.issues.append(
                ValidationIssue(
                    structure="pg",
                    check="dangling_grouped_run",
                    severity="error",
                    column="grouped_runs",
                    message=(
                        f"grouped_runs contains {grouped_run!r}, which is not a "
                        "run_file_name in run.parquet; PG rows keyed by it are "
                        "dropped from all sample-joined views"
                    ),
                )
            )

    def _check_grouped_runs_sample_mapping(self, pg_result: ValidationResult, *, strict: bool) -> None:
        """Require each PG label to resolve to exactly one sample across its grouped runs.

        Each distinct (grouped-run set, label) unit is resolved once, and a unit
        that resolves to no sample or to several is reported once, with the
        number of pg rows that share it.
        """
        try:
            rows = self._engine.execute(_GROUPED_RUN_SAMPLE_MAPPING_SQL).fetchall()
        except duckdb.Error as exc:
            # Do not fail silently on drift (see _check_grouped_runs_referential).
            reason = _skipped_check_reason(exc)
            _log.warning("pg grouped_runs sample-mapping check could not run (%s): %s", reason, exc)
            pg_result.issues.append(
                ValidationIssue(
                    structure="pg",
                    check="sample_mapping_check_skipped",
                    severity="error" if strict else "warning",
                    column=None,
                    message=f"grouped_runs sample-mapping check could not run ({reason}): {exc}",
                )
            )
            return

        for example_protein, grouped_runs, label, pg_rows, sample_count, samples in rows:
            pg_result.issues.append(
                ValidationIssue(
                    structure="pg",
                    check="ambiguous_grouped_run_sample",
                    severity="error",
                    column="grouped_runs",
                    message=(
                        f"{pg_rows} pg row(s) with label {label!r} across "
                        f"grouped_runs={grouped_runs!r} (e.g. PG {example_protein!r}) "
                        f"resolve to {sample_count} samples ({samples or []!r}); "
                        "exactly one sample is required"
                    ),
                )
            )

    # --- Integrity ---
    def _require_local(self, operation: str) -> None:
        """Raise if the dataset is S3-backed; integrity and save require local paths."""
        if self._is_s3:
            raise NotImplementedError(
                f"{operation} is only supported for local datasets. This dataset is S3-backed (path={self.path!r})."
            )

    def compute_integrity(self) -> dict:
        """Compute checksums, row counts, and file sizes for all Parquet files.

        Only supported for local datasets. For S3-backed datasets, raises
        NotImplementedError.

        Returns
        -------
        dict
            Integrity fields suitable for writing to dataset.parquet (file_checksums,
            file_row_counts, file_sizes_bytes, total_structures, packaged_at).
        """
        self._require_local("compute_integrity")
        from datetime import datetime, timezone

        import pyarrow.parquet as pq

        path = Path(self.path)
        checksums, row_counts, sizes = {}, {}, {}

        # Parquet files. rglob, not glob: Hive-partitioned structures live in
        # subdirectories (feature/run_file_name=.../part-0.parquet) and a
        # non-recursive scan silently left them unchecksummed, so the integrity
        # record claimed to cover data it had never read (bigbio/qpx#287).
        for f in sorted(path.rglob("*.parquet")):
            name = f.relative_to(path).as_posix()
            sizes[name] = f.stat().st_size
            checksums[name] = sha256_file(f)
            try:
                row_counts[name] = pq.read_metadata(f).num_rows
            # pylint: disable-next=broad-exception-caught
            except Exception as exc:  # noqa: BLE001 - one bad file must not abort the whole record
                # -1 is the "could not be read" sentinel. Record it, but say so:
                # storing it silently produced an integrity record that looked
                # complete for a file whose metadata was unreadable.
                _log.warning("Could not read Parquet metadata for %s: %s", name, exc)
                row_counts[name] = -1

        # H5AD files (AnnData from downstream tools)
        for f in sorted(path.rglob("*.h5ad")):
            name = f.relative_to(path).as_posix()
            sizes[name] = f.stat().st_size
            checksums[name] = sha256_file(f)
            try:
                import anndata

                row_counts[name] = anndata.read_h5ad(f, backed="r").n_obs
            # pylint: disable-next=broad-exception-caught
            except Exception as exc:  # noqa: BLE001 - one bad file must not abort the whole record
                _log.warning("Could not read AnnData obs count for %s: %s", name, exc)
                row_counts[name] = -1

        return {
            "file_checksums": checksums,
            "file_row_counts": row_counts,
            "file_sizes_bytes": sizes,
            "total_structures": len(checksums),
            "packaged_at": datetime.now(timezone.utc).isoformat(),
        }

    def verify_integrity(self) -> dict[str, list[str]]:
        """Verify dataset files against stored integrity data.

        Only supported for local datasets. For S3-backed datasets, returns
        a result with a single warning and no verification performed.

        Returns
        -------
        dict[str, list[str]]
            Dict with 'errors' and 'warnings' lists.
        """
        errors, warnings = [], []
        verified: set[Path] = set()
        if self._is_s3:
            warnings.append("Integrity verification is not supported for S3-backed datasets.")
            return {"errors": errors, "warnings": warnings}
        if self.dataset_meta is None:
            errors.append("No dataset.parquet found — cannot verify integrity")
            return {"errors": errors, "warnings": warnings}

        meta_df = self.dataset_meta.to_df()
        if meta_df.empty or "file_checksums" not in meta_df.columns:
            warnings.append("No integrity data stored in dataset.parquet")
            return {"errors": errors, "warnings": warnings}

        stored_checksums = meta_df["file_checksums"].iloc[0]
        if not isinstance(stored_checksums, dict):
            warnings.append("file_checksums is null")
            return {"errors": errors, "warnings": warnings}

        def _stored(field: str) -> dict:
            value = meta_df[field].iloc[0] if field in meta_df.columns else None
            return value if isinstance(value, dict) else {}

        stored_row_counts = _stored("file_row_counts")
        stored_sizes = _stored("file_sizes_bytes")

        # Skip only the metadata file carrying this integrity record: nested
        # *.dataset.parquet files are ordinary recorded inputs and must verify.
        path = Path(self.path).resolve()
        dataset_meta_path = Path(self.dataset_meta.file_paths[0]).resolve()
        for name, expected_sha in stored_checksums.items():
            # Keys are dataset-root-relative POSIX paths. Resolve and confirm the
            # result stays inside the dataset before reading it, so a crafted or
            # corrupted key cannot walk outside the directory being verified.
            fpath = (path / name).resolve()
            try:
                fpath.relative_to(path)
            except ValueError:
                errors.append(f"Integrity entry escapes the dataset directory: {name}")
                continue
            if fpath == dataset_meta_path:
                continue
            if not fpath.exists():
                errors.append(f"Missing file: {name}")
                continue
            verified.add(fpath)

            file_errors, file_warnings = _verify_one_file(
                fpath,
                name,
                expected_sha,
                _as_int(stored_sizes.get(name)),
                _as_int(stored_row_counts.get(name)),
            )
            errors.extend(file_errors)
            warnings.extend(file_warnings)

        warnings.extend(_unrecorded_file_warnings(path, dataset_meta_path, verified))

        return {"errors": errors, "warnings": warnings}
