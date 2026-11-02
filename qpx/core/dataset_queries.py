"""Dataset softlinks, abundance queries and analysis views."""

from __future__ import annotations

import logging
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

from qpx.core.convert import QueryResult
from qpx.core.data.identity import derive_id
from qpx.core.query import _escape_sql_string
from qpx.core.sample_resolution import RUN_SAMPLES_SQL, feature_samples_sql
from qpx.core.sql import escape_path, sql_build

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)

# The analyte a peptide-level abundance is keyed by.
_PEPTIDE_KEYS = {
    "peptidoform": "peptidoform",
    "precursor": "peptidoform || '/' || charge",
    "sequence": "sequence",
}
# How design_matrix(by="sample") combines the units of a sample, as DuckDB aggregates.
_UNIT_AGGREGATES = {"mean": "avg", "median": "median", "sum": "sum", "min": "min", "max": "max"}
# The row keys of design_matrix.
_MATRIX_ROWS = {"sample": ("sample_accession",), "unit": ("sample_accession", "label", "unit")}
# pg and feature with an additional intensity ($name) in place of the primary one.
# A pg row carries the additional intensities of its own label.
_PG_ADDITIONAL_SQL = """(
    SELECT p.* REPLACE (x.intensity_value AS intensity)
    FROM pg p
    CROSS JOIN UNNEST(p.additional_intensities) AS _a(a)
    CROSS JOIN UNNEST(a.intensities) AS _x(x)
    WHERE a.label = p.label AND x.intensity_name = '$name'
)"""
_FEATURE_ADDITIONAL_SQL = """(
    SELECT * REPLACE (
        list_filter(
            list_transform(additional_intensities, lambda a: {
                'label': a.label,
                'intensity': list_filter(a.intensities, lambda x: x.intensity_name = '$name')[1].intensity_value
            }),
            lambda e: e.intensity IS NOT NULL
        ) AS intensities
    )
    FROM feature
)"""
_ADDITIONAL_NAMES_SQL = """
    SELECT DISTINCT x.intensity_name
    FROM (SELECT UNNEST(a.intensities) AS x FROM (SELECT UNNEST(additional_intensities) AS a FROM $table))
    ORDER BY 1
"""


class _DatasetQueryMixin:
    """Query links, intensities and analysis views across registered structures."""

    # --- Cross-structure queries ---
    def sql(self, query: str) -> QueryResult:
        """Execute arbitrary SQL across registered structures."""
        return QueryResult(self._engine.connection, query)

    # --- Feature <-> protein-group softlink (bigbio/qpx#269) ---
    #
    # The feature->pg association is a *computed softlink*: instead of persisting
    # feature.pg_ids on write, it is derived on read by joining the registered
    # ``feature`` and ``pg`` views. The join is label-aware so a feature links
    # only to the pg rows for the labels it actually carries (LFQ: one row; TMT /
    # plexDIA: one row per channel the feature has), never over-linking to a
    # same-membership pg row that exists only for a channel the feature lacks.
    #
    # Match rule (all three must hold):
    #   1. membership  — canonical(feature.pg_accessions) == canonical(pg.pg_accessions),
    #                    set-wise (sorted accession strings), the same canonical
    #                    membership that keys the pg identity;
    #   2. run         — feature.run_file_name IN pg.grouped_runs;
    #   3. label       — a label from feature.intensities IS NOT DISTINCT FROM pg.label.
    # ``pg.pg_id`` is read straight from the matched pg row (never re-derived).
    # A feature whose group has no matching pg row — or which carries no
    # quantified label — simply produces no rows: identified but not quantified.
    # ``grouped_runs`` is unnested so the run test is an equality: all three
    # conditions are then hash-join keys instead of a filter applied to every
    # same-membership pair. The view is unordered; link_feature_pg() sorts.
    # Files with ``pg_group_id`` (spec 1.3) compare that flat key instead of lists.
    _LINK_FEATURE_PG_SQL = """
        WITH feat_labels AS (
            SELECT
                f.feature_id AS feature_id,
                f.run_file_name AS run_file_name,
                $feature_membership AS membership,
                UNNEST(list_transform(f.intensities, lambda i: i.label)) AS feature_label
            FROM feature f
            $feature_filter
        ),
        pg_runs AS (
            SELECT
                p.pg_id AS pg_id,
                p.label AS label,
                $pg_membership AS membership,
                UNNEST(p.grouped_runs) AS run_file_name
            FROM pg p
            $pg_filter
        )
        SELECT DISTINCT
            fl.feature_id AS feature_id,
            pr.pg_id AS pg_id,
            pr.label AS label
        FROM feat_labels fl
        JOIN pg_runs pr
            ON pr.membership = fl.membership
           AND pr.run_file_name = fl.run_file_name
           AND pr.label IS NOT DISTINCT FROM fl.feature_label
    """

    # Fixed identifiers for the DuckDB temp views that back the softlinks. The
    # link diagnostics reference these views by these literal names inside
    # constant-literal queries (no f-string / concat SQL), so there is no
    # formatted-SQL construction to flag.
    _FEATURE_PG_LINK_VIEW = "qpx_feature_pg_link"
    _FEATURE_PSM_LINK_VIEW = "qpx_feature_psm_link"

    def _feature_pg_link_sql(self, group: list[str] | None = None) -> str:
        """The softlink query, on ``pg_group_id`` when both views have it, optionally of one protein group."""
        by_id = all("pg_group_id" in self._engine.view_columns(view) for view in ("feature", "pg"))
        feature_membership = "f.pg_group_id" if by_id else "list_sort(list_transform(f.pg_accessions, lambda x: x.accession))"
        pg_membership = "p.pg_group_id" if by_id else "list_sort(p.pg_accessions)"
        feature_filter = pg_filter = ""
        if group is not None:
            members = "[" + ", ".join(f"'{_escape_sql_string(accession)}'" for accession in sorted(group)) + "]"
            key = str(derive_id([sorted(group)])) if by_id else members
            feature_filter = f"WHERE {feature_membership} = {key}"
            pg_filter = f"WHERE {pg_membership} = {key}"
        return sql_build(
            self._LINK_FEATURE_PG_SQL,
            feature_membership=feature_membership,
            pg_membership=pg_membership,
            feature_filter=feature_filter,
            pg_filter=pg_filter,
        )

    def _register_feature_pg_link_view(self) -> None:
        """Register the feature->pg softlink as the ``qpx_feature_pg_link`` view."""
        self._engine.register_view(self._FEATURE_PG_LINK_VIEW, self._feature_pg_link_sql())

    def _register_feature_psm_link_view(self) -> None:
        """Register the feature<->psm softlink as the ``qpx_feature_psm_link`` view."""
        self._engine.register_view(self._FEATURE_PSM_LINK_VIEW, self._LINK_FEATURE_PSM_SQL)

    def link_feature_pg(self, group: list[str] | None = None) -> QueryResult:
        """Compute the feature->pg softlink as ``(feature_id, pg_id, label)`` rows.

        ``group`` limits it to the protein group with exactly these members, which
        is filtered before the join and takes a fraction of the whole link's time.

        A label-aware read-side join of the registered ``feature`` and ``pg``
        views (bigbio/qpx#269): a feature links to a pg row only when their
        canonical ``pg_accessions`` membership matches set-wise, the feature's
        ``run_file_name`` is in the pg row's ``grouped_runs``, and the pg row's
        ``label`` is one of the labels in the feature's ``intensities``. ``pg_id``
        is taken directly from the matched pg row (not re-derived). Works the same
        on the local and S3 read paths (it is a DuckDB query over the views).

        Returns
        -------
        QueryResult
            Lazy result of ``(feature_id, pg_id, label)`` — one row per matched
            (feature, channel). Features identified but not quantified in a given
            channel/fraction produce no rows (this is correct, not an error).
        """
        self._require_feature_pg("link_feature_pg")
        if group is not None:
            return QueryResult(self._engine.connection, f"{self._feature_pg_link_sql(group)} ORDER BY feature_id, pg_id")
        self._register_feature_pg_link_view()
        return QueryResult(
            self._engine.connection, "SELECT feature_id, pg_id, label FROM qpx_feature_pg_link ORDER BY feature_id, pg_id"
        )

    def features_without_pg_link(self) -> QueryResult:
        """Return ``feature_id`` for features that resolve to no pg row.

        The features whose computed softlink (:meth:`link_feature_pg`) is
        empty — no pg row shares their canonical membership, run and a carried
        label. A linked feature can still lack a quantity, since the matched pg
        row may have a null ``intensity``. Minimal by design; a fuller
        conversion summary is built on top of the softlink.
        """
        self._require_feature_pg("features_without_pg_link")
        self._register_feature_pg_link_view()
        return QueryResult(
            self._engine.connection,
            "SELECT feature_id FROM feature "
            "WHERE feature_id NOT IN (SELECT feature_id FROM qpx_feature_pg_link) "
            "ORDER BY feature_id",
        )

    def _require_feature_pg(self, operation: str) -> None:
        """Raise if the feature or pg view is not available for the softlink."""
        if self.feature is None or self.pg is None:
            raise ValueError(
                f"{operation}() requires both the 'feature' and 'pg' structures; available: {self.available_structures}"
            )

    # --- PSM <-> Feature softlink (bigbio/qpx#267) ---
    #
    # The association is *asymmetric*. ``psm.feature_id`` is the authoritative,
    # producer-supplied optional foreign key: which consensus feature a PSM
    # belongs to (openms_consensus stamps it; null when a producer does not
    # supply it). It is NOT recomputable from shared identification fields
    # ((peptidoform, charge, run) is many-to-one and ambiguous), so it stays
    # materialized. Its inverse, ``feature.psm_ids``, is redundant and is NOT
    # materialized by qpx — it is computed on read here by grouping the persisted
    # ``psm.feature_id`` by feature: the (feature_id, psm_id) pairs below ARE the
    # feature<->psm association, and grouping by feature_id gives feature.psm_ids.
    _LINK_FEATURE_PSM_SQL = """
        SELECT feature_id, psm_id
        FROM psm
        WHERE feature_id IS NOT NULL
        ORDER BY feature_id, psm_id
    """

    def link_feature_psm(self) -> QueryResult:
        """Compute the feature<->psm softlink as ``(feature_id, psm_id)`` rows.

        The inverse of the authoritative ``psm.feature_id`` foreign key
        (bigbio/qpx#267), read from the persisted psm view: every PSM with a
        non-null ``feature_id`` yields one ``(feature_id, psm_id)`` pair. Grouping
        the result by ``feature_id`` recovers each feature's ``psm_ids`` (which
        qpx no longer materializes). Works the same on the local and S3 read paths
        (it is a DuckDB query over the psm view).

        Returns
        -------
        QueryResult
            Lazy result of ``(feature_id, psm_id)`` — one row per assigned PSM.
            PSMs with no feature assignment (null ``feature_id``) produce no rows.
        """
        self._require_feature_psm("link_feature_psm")
        return QueryResult(self._engine.connection, self._LINK_FEATURE_PSM_SQL)

    def psms_without_feature(self) -> QueryResult:
        """Return ``psm_id`` for PSMs not assigned to any feature.

        Convenience for the complement of :meth:`link_feature_psm`: the PSMs whose
        authoritative ``psm.feature_id`` is null — identified but not tied to a
        quantified consensus feature.
        """
        self._require_feature_psm("psms_without_feature")
        return QueryResult(self._engine.connection, "SELECT psm_id FROM psm WHERE feature_id IS NULL ORDER BY psm_id")

    def _require_feature_psm(self, operation: str) -> None:
        """Raise if the feature or psm view is not available for the softlink."""
        if self.feature is None or self.psm is None:
            raise ValueError(
                f"{operation}() requires both the 'feature' and 'psm' structures; available: {self.available_structures}"
            )

    def _link_scalar(self, sql: str) -> int:
        """Run a constant-literal scalar COUNT."""
        return int(self._engine.execute(sql).fetchall()[0][0])

    def feature_link_diagnostics(self) -> dict:
        """Return feature<->pg and feature<->psm softlink counts for reporting.

        Public diagnostics API consumed by the conversion summary. Each count is
        computed with a constant-literal ``COUNT`` / ``COUNT(DISTINCT ...)`` over
        the registered softlink views (:attr:`_FEATURE_PG_LINK_VIEW`,
        :attr:`_FEATURE_PSM_LINK_VIEW`), which are (re)registered here on demand.

        Returns
        -------
        dict
            Keys ``n_feature_pg_links``, ``n_features_linked``,
            ``n_features_without_pg``, ``n_features_with_psm`` and
            ``n_psms_without_feature``. The feature<->pg group is ``None`` unless
            both the ``feature`` and ``pg`` views exist; the feature<->psm group
            is ``None`` unless both the ``feature`` and ``psm`` views exist.
        """
        diagnostics: dict = {
            "n_feature_pg_links": None,
            "n_features_linked": None,
            "n_features_without_pg": None,
            "n_features_with_psm": None,
            "n_psms_without_feature": None,
        }

        if self.feature is not None and self.pg is not None:
            self._register_feature_pg_link_view()
            # One pass over the join: each count read from the view re-ran it.
            links, linked, unlinked = self._engine.execute(
                """
                WITH link AS MATERIALIZED (SELECT feature_id FROM qpx_feature_pg_link),
                linked AS (SELECT DISTINCT feature_id FROM link)
                SELECT
                    (SELECT COUNT(*) FROM link),
                    (SELECT COUNT(*) FROM linked),
                    (SELECT COUNT(*) FROM feature WHERE feature_id NOT IN (SELECT feature_id FROM linked))
                """
            ).fetchall()[0]
            diagnostics["n_feature_pg_links"] = int(links)
            diagnostics["n_features_linked"] = int(linked)
            diagnostics["n_features_without_pg"] = int(unlinked)

        if self.feature is not None and self.psm is not None:
            self._register_feature_psm_link_view()
            diagnostics["n_features_with_psm"] = self._link_scalar("SELECT COUNT(DISTINCT feature_id) FROM qpx_feature_psm_link")
            diagnostics["n_psms_without_feature"] = self._link_scalar("SELECT COUNT(*) FROM psm WHERE feature_id IS NULL")

        return diagnostics

    # --- Analysis helpers ---

    def _abundance_sql(self, level: str, key: str | None = None, value_col: str = "intensity") -> str:
        """Long-form abundance SQL: sample_accession, label, unit, feature_id, intensity.

        A row holds the quantity of one feature in one sample and label of one
        quantification unit; the fractions of the unit are summed. *value_col*
        names the quantity: the primary ``intensity`` or an additional intensity.
        """
        name = None if value_col == "intensity" else _escape_sql_string(value_col)
        if level == "protein":
            if key is not None:
                raise ValueError("key applies to level='peptide'; proteins are keyed by their group membership")
            if self.pg is None or self.run is None:
                raise ValueError("level='protein' requires pg and run structures.")
            # Since QPX 1.1 the pg view is flattened: one row per label with scalar
            # p.label / p.intensity. Any member run of grouped_runs resolves the
            # label to the unit's sample (pg.md), so the first one is joined. ID-only
            # rows (null label / intensity) never match a sample label. The group,
            # not its anchor protein, is the analyte: groups can share a leader.
            return sql_build(
                """
            WITH run_samples AS ($run_samples)
            SELECT rs.sample_accession, rs.label, rs.unit,
                   array_to_string(list_sort(p.pg_accessions), ';') AS feature_id,
                   sum(p.intensity) AS intensity
            FROM $source p
            JOIN run_samples rs ON rs.run_file_name = p.grouped_runs[1] AND rs.label = p.label
            WHERE p.is_decoy = false
            GROUP BY ALL
            """,
                run_samples=RUN_SAMPLES_SQL,
                source="pg" if name is None else sql_build(_PG_ADDITIONAL_SQL, name=name),
            )
        if level == "peptide":
            expression = _PEPTIDE_KEYS.get(key or "peptidoform")
            if expression is None:
                raise ValueError(f"key must be one of {sorted(_PEPTIDE_KEYS)}, got {key!r}")
            if self.feature is None or self.run is None:
                raise ValueError("level='peptide' requires feature and run structures.")
            return sql_build(
                """
            SELECT sample_accession, label, unit, $key AS feature_id, sum(intensity) AS intensity
            FROM ($feature_samples)
            WHERE is_decoy = false
            GROUP BY ALL
            """,
                key=expression,
                feature_samples=feature_samples_sql(
                    "feature" if name is None else sql_build(_FEATURE_ADDITIONAL_SQL, name=name),
                    ("sequence", "peptidoform", "charge", "is_decoy"),
                ),
            )
        raise ValueError(f"level must be 'protein' or 'peptide', got '{level}'")

    def intensity(self, level: str = "protein", key: str | None = None) -> QueryResult:
        """Lazy long-form intensity query — scalable for large datasets.

        Returns a QueryResult that stays lazy until you materialize it.
        DuckDB handles memory management, so this works on datasets that
        don't fit in memory.

        Parameters
        ----------
        level : str
            "protein" (uses pg + run) or "peptide" (uses feature + run).
        key : str | None
            The peptide analyte: ``"peptidoform"`` (default), ``"precursor"``
            (peptidoform and charge) or ``"sequence"``. Proteins are keyed by their
            group membership, the sorted ``pg_accessions`` joined by ``;``.

        Returns
        -------
        QueryResult
            Long-form table with columns: sample_accession, label, unit,
            feature_id, intensity. ``unit`` names the quantification unit (the first
            run of the fractions that hold the same samples), so technical
            replicates and a bridge sample in several plexes stay separate rows.
            Call .to_df(), .to_arrow(), .to_polars(), or iterate row-by-row.

        Examples
        --------
        Materialize when data fits in memory:
            df = ds.intensity("protein").to_df()

        Arrow is more compact than pandas:
            table = ds.intensity("protein").to_arrow()
        """
        sql = self._abundance_sql(level, key)
        return QueryResult(self._engine.connection, sql)

    def _additional_intensity_names(self, level: str) -> list[str]:
        """Names of the additional intensities of the level's table."""
        table = "pg" if level == "protein" else "feature"
        rows = self._engine.execute(sql_build(_ADDITIONAL_NAMES_SQL, table=table)).fetchall()
        return [row[0] for row in rows]

    def _warn_multilabel_samples(self, agg: str) -> None:
        """Warn when a sample carries several labels in one unit, which by='sample' combines."""
        count, example = self._engine.execute(
            sql_build(
                """
            WITH run_samples AS ($run_samples)
            SELECT count(*), any_value(sample_accession)
            FROM (SELECT unit, sample_accession FROM run_samples GROUP BY ALL HAVING count(DISTINCT label) > 1)
            """,
                run_samples=RUN_SAMPLES_SQL,
            )
        ).fetchall()[0]
        if count:
            logger.warning(
                "%d sample(s) carry several labels within one unit (e.g. %s); by='sample' combines them with %s, "
                "by='unit' keeps them apart",
                count,
                example,
                agg,
            )

    def design_matrix(
        self,
        level: str = "protein",
        value_col: str = "intensity",
        fillna: float | None = None,
        output_path: str | Path | None = None,
        *,
        key: str | None = None,
        by: str = "sample",
        agg: str = "mean",
    ) -> "pd.DataFrame | Path":
        """Pivot intensity data into a matrix of samples (or units) by features.

        Rows are samples. A sample measured in several quantification units —
        technical replicates, or a bridge channel in several plexes — gets *agg*
        of its unit values, after the fractions of each unit are summed.
        ``by="unit"`` keeps one row per (sample_accession, label, unit) instead,
        which also keeps apart the labels of a sample that one unit carries twice
        (pulsed SILAC, a heavy standard).

        For small-to-medium datasets, returns a pandas DataFrame.
        For large datasets, pass *output_path* to write the pivot directly
        to Parquet via DuckDB (out-of-core, no full pandas materialization).

        For streaming access to large data without pivoting, use
        :meth:`intensity` instead.

        Parameters
        ----------
        level : str
            "protein" (uses pg + run) or "peptide" (uses feature + run).
        value_col : str
            The quantity to pivot: ``"intensity"`` (the primary quantity) or the
            ``intensity_name`` of an additional intensity, such as ``"lfq"``.
        fillna : float | None
            Value to fill missing cells; *None* (the default) keeps NaN, so a
            quantity that was not measured is not read as zero.
            Ignored when *output_path* is set.
        output_path : str | Path | None
            If set, writes the pivoted matrix directly to this Parquet path
            using DuckDB PIVOT (out-of-core). Returns the Path instead of a
            DataFrame.
        key : str | None
            The peptide analyte (see :meth:`intensity`).
        by : str
            ``"sample"`` (default) or ``"unit"``.
        agg : str
            How ``by="sample"`` combines the units of a sample: ``"mean"``
            (default), ``"median"``, ``"sum"``, ``"min"`` or ``"max"``.

        Returns
        -------
        pd.DataFrame | Path
            DataFrame (in-memory) or Path (when output_path is set).
        """
        import pandas as pd

        rows = _MATRIX_ROWS.get(by)
        if rows is None:
            raise ValueError(f"by must be one of {sorted(_MATRIX_ROWS)}, got {by!r}")
        if agg not in _UNIT_AGGREGATES:
            raise ValueError(f"agg must be one of {sorted(_UNIT_AGGREGATES)}, got {agg!r}")
        sql = self._abundance_sql(level, key, value_col)
        if value_col != "intensity":
            names = self._additional_intensity_names(level)
            if value_col not in names:
                raise ValueError(f"value_col must be 'intensity' or an additional intensity {names}, got {value_col!r}")
        if by == "sample":
            self._warn_multilabel_samples(agg)

        if output_path is not None:
            # Out-of-core: DuckDB PIVOT → Parquet, never touches pandas
            output_path = Path(output_path)
            # DuckDB 1.6 development builds crash on PIVOT of a subquery holding CTEs, but not of a view.
            self._engine.execute(f"CREATE OR REPLACE TEMP VIEW _matrix_source AS {sql}")
            pivot_sql = sql_build(
                """COPY (
                PIVOT _matrix_source ON feature_id USING $agg(intensity) GROUP BY $rows
            ) TO '$out_path' (FORMAT PARQUET)""",
                agg=_UNIT_AGGREGATES[agg],
                rows=", ".join(rows),
                out_path=escape_path(str(output_path)),
            )
            self._engine.execute(pivot_sql)
            self._engine.execute("DROP VIEW _matrix_source")
            return output_path

        # In-memory path for small/medium datasets
        df = self._engine.execute(sql).fetchdf()

        if df.empty:
            return pd.DataFrame()

        grouped = df.groupby([*rows, "feature_id"])["intensity"]
        values = grouped.sum(min_count=1) if agg == "sum" else grouped.agg(agg)
        matrix = values.unstack("feature_id")
        matrix.columns.name = None

        if fillna is not None:
            matrix = matrix.fillna(fillna)

        return matrix

    # --- View accessors (cached) ---
    @cached_property
    def protein_view(self):
        """Protein intensity per sample (pg joined with run)."""
        from qpx.views.api import ProteinView

        return ProteinView(self)

    @cached_property
    def peptide_view(self):
        """Peptide intensity per sample (feature joined with run)."""
        from qpx.views.api import PeptideView

        return PeptideView(self)

    @cached_property
    def identification_summary(self):
        """Per-run counts of proteins, peptides and features."""
        from qpx.views.api import IdentificationSummaryView

        return IdentificationSummaryView(self)

    @cached_property
    def run_summary(self):
        """Per-run statistics with the intensity distribution."""
        from qpx.views.api import RunSummaryView

        return RunSummaryView(self)

    @cached_property
    def modification_view(self):
        """Modification frequency across PSMs."""
        from qpx.views.api import ModificationView

        return ModificationView(self)

    @cached_property
    def qc_view(self):
        """Dataset-wide QC metrics."""
        from qpx.views.api import QualityControlView

        return QualityControlView(self)

    @cached_property
    def sample_summary(self):
        """Per-sample and per-run peptide and protein counts."""
        from qpx.views.api import SampleSummaryView

        return SampleSummaryView(self)

    @cached_property
    def ae_view(self):
        """AnnData view of the absolute expression file."""
        from qpx.views.api import AbsoluteExpressionView

        return AbsoluteExpressionView(self)
