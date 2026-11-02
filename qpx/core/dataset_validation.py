"""Dataset schema, cross-structure and file integrity checks."""

from __future__ import annotations

import logging
from functools import partial
from pathlib import Path

import duckdb

from qpx.core.anndata_validation import validate_anndata_view
from qpx.core.data.schema import ValidationIssue, ValidationResult
from qpx.core.data.validation import _LINK_CONTEXT_MESSAGE, _LINK_MEMBERSHIP_MESSAGE, _LINK_SQL, _issue
from qpx.core.files import sha256_file
from qpx.core.label_validation import (
    _DANGLING_GROUPED_RUN_SQL,
    _DANGLING_RUN_FILE_NAME_SQL,
    _DUPLICATE_FEATURE_LABEL_SQL,
    _FEATURE_LABEL_SAMPLE_SQL,
    _GROUPED_RUN_SAMPLE_MAPPING_SQL,
    _RUN_KEYED_VIEWS,
    _examples,
    _feature_label_issues,
    _feature_psm_references,
    _feature_reference_issues,
    _skipped_check_reason,
    _skipped_issue,
    _view_columns,
    fraction_order_issues,
    fragment_link_issues,
    pg_link_coverage_issues,
    psm_label_issues,
    psm_reference_issues,
    run_sample_issues,
)
from qpx.core.ontology_fields import ontology_field_issues

_log = logging.getLogger("qpx.dataset")


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
    def validate(
        self, structures: list[str] | None = None, *, strict: bool = False, semantic: bool = False, check_identity: bool = False
    ) -> dict[str, ValidationResult]:
        """Validate the dataset or specific structures against their schemas.

        Parameters
        ----------
        structures : list[str] | None
            Structure names to validate.  If *None*, validates all available, and reports
            those whose files could not be read as errors.
        strict : bool
            Report the issues that are warnings by default as errors.
        semantic : bool
            Also check PSM and feature peptidoforms, modification positions and m/z values.
        check_identity : bool
            Also re-derive each id from the composite its footer declares.

        Returns
        -------
        dict[str, ValidationResult]
            Mapping of structure name to its validation result.
        """
        targets = structures or [name for name in self.structure_names() if name in self._structures or name in self._skipped]
        if not targets:
            message = f"No QPX views found in {self.path}"
            issue = ValidationIssue(structure="dataset", check="no_structures", severity="error", column=None, message=message)
            return {"dataset": ValidationResult(structure="dataset", issues=[issue])}
        results = {
            name: self._validate_structure(name, strict=strict, semantic=semantic, check_identity=check_identity)
            for name in targets
        }
        present = {name for name in results if self._structures.get(name) is not None}
        if "dataset" in present and (rows := self.dataset_meta.count()) != 1:
            message = f"The dataset view describes one dataset, so it has one row; it has {rows}"
            results["dataset"].issues.append(
                ValidationIssue(structure="dataset", check="single_row", severity="error", column=None, message=message)
            )
        self._check_run_resolution(results, present, strict=strict)
        self._check_cross_references(results, present, semantic=semantic, strict=strict)
        if structures is None and not self._is_s3:
            for view in ("ae", "de"):
                for path in sorted(Path(self.path).glob(f"{self._file_prefix or '*'}.{view}.h5ad")):
                    results[view] = validate_anndata_view(path, view)
        return results

    def _validate_structure(self, name: str, *, strict: bool, semantic: bool, check_identity: bool) -> ValidationResult:
        """Validate one registered structure, or report it missing or unreadable."""
        struct = self._structures.get(name)
        if struct is not None:
            return struct.validate(strict=strict, semantic=semantic, check_identity=check_identity)
        if name in self._skipped:
            check, message = "unreadable_structure", f"Structure '{name}' could not be read: {self._skipped[name]}"
        else:
            check, message = "missing_structure", f"Structure '{name}' not found in dataset at {self.path}"
        return ValidationResult(structure=name, issues=[_issue(name, check, "error", None, message)])

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
            self._cross_view_issues(results["pg"], "fraction order", partial(fraction_order_issues, self._engine), strict)
        # feature, psm and mz rows each name the run they were measured in, and
        # each feature label must name exactly one sample of that run.
        for view in _RUN_KEYED_VIEWS:
            if view in present:
                self._check_run_file_names(view, results[view], strict=strict)
        if "feature" in present:
            self._check_feature_label_samples(results["feature"], strict=strict)
        # Every sample that run.samples names is in the sample view.
        if {"run", "sample"} <= present:
            self._cross_view_issues(results["run"], "run-sample", partial(run_sample_issues, self._engine), strict)

    def _check_cross_references(
        self, results: dict[str, ValidationResult], present: set[str], *, semantic: bool, strict: bool
    ) -> None:
        """Cross-structure checks of the id links between views."""
        severity = "error" if strict else "warning"
        # Cross-structure invariant: the optional feature<->psm cross-references
        # must resolve. Every non-null psm.feature_id must name a real
        # feature.feature_id, and every feature.psm_ids element a real psm.psm_id.
        # Unresolved references are warnings during normal use and errors under
        # strict validation, but only when both views are present. A PSM's label
        # must be one its feature was identified in.
        if {"feature", "psm"} <= present:
            self._check_feature_psm_referential(results["feature"], results["psm"], strict=strict)
            self._cross_view_issues(results["psm"], "psm label", partial(psm_label_issues, self._engine, severity), strict)

        # feature.pg_ids is the canonical feature -> protein-group reference.
        # When both views are present, every populated id must resolve to pg.pg_id.
        if {"feature", "pg"} <= present:
            self._check_feature_pg_links(results["feature"], semantic=semantic, strict=strict)

        if {"feature", "fragment"} <= present:
            self._cross_view_issues(
                results["fragment"], "fragment link", partial(fragment_link_issues, self._engine, severity), strict
            )

        if "ontology" in present:
            self._cross_view_issues(
                results["ontology"], "ontology fields", partial(ontology_field_issues, self._engine, present), strict
            )

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

    @staticmethod
    def _cross_view_issues(result: ValidationResult, check: str, issues, strict: bool) -> None:
        """Add the issues of the *issues* callable, or a skipped-check issue when its query fails."""
        try:
            result.issues.extend(issues())
        except duckdb.Error as exc:
            what = f"{result.structure} {check} check"
            result.issues.append(_skipped_issue(result.structure, "referential_check_skipped", what, exc, strict=strict))

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

    def _check_feature_pg_links(self, feature_result: ValidationResult, *, semantic: bool, strict: bool) -> None:
        """Flag dangling feature.pg_ids and, for a semantic check, a low share of linked features."""
        self._check_feature_pg_referential(feature_result, strict=strict)
        if semantic:
            self._cross_view_issues(
                feature_result,
                "pg link coverage",
                lambda: pg_link_coverage_issues(self.feature_link_diagnostics()),
                strict,
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

        psm_result.issues.extend(psm_reference_issues(dangling_feature_ids, desync, severity))
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
