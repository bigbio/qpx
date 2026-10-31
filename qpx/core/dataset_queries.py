"""Dataset softlinks, abundance queries and analysis views."""

from __future__ import annotations

from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

from qpx.core.convert import QueryResult
from qpx.core.sample_sql import feature_sample_ctes, pg_sample_ctes
from qpx.core.sql import escape_path, sql_build

if TYPE_CHECKING:
    import pandas as pd


# design_matrix aggregations and the DuckDB functions of its out-of-core pivot.
_AGGREGATIONS = ("mean", "median", "sum", "min", "max")
_SQL_AGGREGATES = {"mean": "avg", "median": "median", "sum": "sum", "min": "min", "max": "max"}

# Matrix cells that would receive more than one value, the number of samples
# measured in more than one quantification unit, and one of them for the error.
_REPEATED_CELLS_SQL = """
    WITH base AS MATERIALIZED ($base_sql),
    repeated_samples AS (
        SELECT sample_accession FROM base
        GROUP BY sample_accession HAVING count(DISTINCT quantification_unit) > 1
    )
    SELECT (SELECT count(*) FROM (SELECT 1 FROM base GROUP BY $keys, feature_id HAVING count(*) > 1)),
           (SELECT count(*) FROM repeated_samples),
           (SELECT min(sample_accession) FROM repeated_samples)
"""

# Why a sample has several quantification units, for the design_matrix error.
# Converters that write a pg row per raw file (grouped_runs holds one run: CDAP,
# Spectronaut) leave each fraction of a sample a protein unit of its own that
# holds part of the sample's quantity; peptide units join such fractions.
_REPEATED_UNIT_CAUSES = (
    "technical replicates",
    "a reference channel repeated across plexes",
    "several labels of one sample",
)
_PER_FRACTION_PG_CAUSE = (
    "fractions whose pg rows name only their own run in grouped_runs "
    "(CDAP, Spectronaut; agg='sum' adds up their partial quantities)"
)


# Since QPX 1.1 the pg view is flattened: one row per label with scalar
# p.label / p.intensity (no intensities list to UNNEST). Each pg row matches the
# sample carrying its label in its quantification unit (pg_unit_samples). ID-only
# rows (null label / intensity) never match a sample label and drop out.
_PROTEIN_ABUNDANCE_SQL = """
    SELECT us.sample_accession,
           p.anchor_protein AS feature_id,
           p.intensity,
           us.unit_runs || '|' || p.label AS quantification_unit
    FROM pg p
    JOIN pg_unit_samples us
      ON us.grouped_runs = p.grouped_runs
     AND us.sample_label = p.label
    WHERE p.is_decoy = false
"""

# A peptidoform's features are summed within its quantification unit (charge
# states, fractions of one sample), never across units, and modified forms (a
# heavy standard, a phosphosite) stay apart. The unit's channel is the sample
# label, except in a single-sample unit whose features carry several labels
# (SILAC or plexDIA L/H without SDRF labels): there each label is its own unit,
# not a sum. That is decided per unit, over all its runs, not per feature or per
# run, so a peptide quantified only in H, or a fraction that saw only L, still
# lands in the unit's H or L channel. Label-free units carry one feature label
# (DIA-NN's "raw", "LFQ") and keep the sample label. Only single-sample runs are
# counted: a multiplexed run matches each label exactly. The distinct labels of
# each run are taken first, so the join to the units is on a few rows per run.
_UNIT_LABELS_CTE = """
    unit_labels AS (
        SELECT u.unit_runs, count(DISTINCT rl.intensity_label) AS n_unit_labels
        FROM (
            SELECT DISTINCT run_file_name, intensity_label
            FROM feature_intensities
            WHERE run_file_name IN (SELECT run_file_name FROM run_samples WHERE n_samples = 1)
        ) rl
        JOIN run_units u USING (run_file_name)
        GROUP BY u.unit_runs
    )"""

_PEPTIDE_ABUNDANCE_SQL = """
    SELECT sample_accession,
           peptidoform AS feature_id,
           SUM(intensity) AS intensity,
           unit_runs || '|' || unit_label AS quantification_unit
    FROM (
        SELECT fs.sample_accession, fs.peptidoform, fs.intensity, fs.unit_runs,
               CASE WHEN ul.n_unit_labels > 1 THEN coalesce(fs.intensity_label, fs.sample_label)
                    ELSE fs.sample_label END AS unit_label
        FROM feature_samples fs
        LEFT JOIN unit_labels ul USING (unit_runs)
    )
    GROUP BY sample_accession, unit_runs, unit_label, peptidoform
"""


def _design_matrix_keys(value_col: str, agg: str | None, by: str) -> list[str]:
    """Validate the design_matrix options and return its row keys."""
    if value_col != "intensity":
        raise ValueError(f"value_col must be 'intensity', got {value_col!r}")
    if agg is not None and agg not in _AGGREGATIONS:
        raise ValueError(f"agg must be None or one of {_AGGREGATIONS}, got {agg!r}")
    if by not in ("sample", "unit"):
        raise ValueError(f"by must be 'sample' or 'unit', got {by!r}")
    return ["sample_accession"] if by == "sample" else ["sample_accession", "quantification_unit"]


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
    _LINK_FEATURE_PG_SQL = """
        WITH feat_labels AS (
            SELECT
                f.feature_id AS feature_id,
                f.run_file_name AS run_file_name,
                list_sort(list_transform(f.pg_accessions, lambda x: x.accession)) AS membership,
                UNNEST(list_transform(f.intensities, lambda i: i.label)) AS feature_label
            FROM feature f
        ),
        pg_runs AS (
            SELECT
                p.pg_id AS pg_id,
                p.label AS label,
                list_sort(p.pg_accessions) AS membership,
                UNNEST(p.grouped_runs) AS run_file_name
            FROM pg p
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

    def _register_feature_pg_link_view(self) -> None:
        """Register the feature->pg softlink as the ``qpx_feature_pg_link`` view."""
        self._engine.register_view(self._FEATURE_PG_LINK_VIEW, self._LINK_FEATURE_PG_SQL)

    def _register_feature_psm_link_view(self) -> None:
        """Register the feature<->psm softlink as the ``qpx_feature_psm_link`` view."""
        self._engine.register_view(self._FEATURE_PSM_LINK_VIEW, self._LINK_FEATURE_PSM_SQL)

    def link_feature_pg(self) -> QueryResult:
        """Compute the feature->pg softlink as ``(feature_id, pg_id, label)`` rows.

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

    def _abundance_sql(self, level: str) -> str:
        """Build the long-form abundance SQL for a given level.

        Returns SQL that produces columns: sample_accession, feature_id,
        intensity, quantification_unit.
        """
        if level == "protein":
            if self.pg is None or self.run is None:
                raise ValueError("level='protein' requires pg and run structures.")
            return pg_sample_ctes() + _PROTEIN_ABUNDANCE_SQL
        if level == "peptide":
            if self.feature is None or self.run is None:
                raise ValueError("level='peptide' requires feature and run structures.")
            ctes = feature_sample_ctes(("peptidoform",), feature_where="f.is_decoy = false")
            return ctes + "," + _UNIT_LABELS_CTE + _PEPTIDE_ABUNDANCE_SQL
        raise ValueError(f"level must be 'protein' or 'peptide', got '{level}'")

    def intensity(self, level: str = "protein") -> QueryResult:
        """Lazy long-form intensity query — scalable for large datasets.

        Returns a QueryResult that stays lazy until you materialize it.
        DuckDB handles memory management, so this works on datasets that
        don't fit in memory.

        Parameters
        ----------
        level : str
            "protein" (uses pg + run) or "peptide" (uses feature + run).

        Returns
        -------
        QueryResult
            Long-form table with columns: sample_accession, feature_id,
            intensity, quantification_unit. A sample measured in several
            quantification units (technical replicates, a reference channel in
            several plexes) has one row per unit. The unit is named by its runs,
            sorted and joined by ``;``, then ``|`` and the channel label
            (``F1;F2|TMT126``). Call .to_df(), .to_arrow(), .to_polars(), or
            iterate row-by-row.

        Examples
        --------
        Lazy iteration (constant memory):
            for row in ds.intensity("protein"):
                sample, protein, intensity, unit = row

        Materialize when data fits in memory:
            df = ds.intensity("protein").to_df()

        Write directly to Parquet (out-of-core):
            ds.intensity("protein").to_arrow()  # Arrow is more compact than pandas
        """
        sql = self._abundance_sql(level)
        return QueryResult(self._engine.connection, sql)

    def design_matrix(
        self,
        level: str = "protein",
        value_col: str = "intensity",
        fillna: float | None = None,
        output_path: str | Path | None = None,
        agg: str | None = None,
        by: str = "sample",
    ) -> "pd.DataFrame | Path":
        """Pivot intensity data into a samples-by-features matrix (DataFrame or Parquet).

        Each cell holds one quantification unit's value. A sample measured in
        several units (technical replicates, a reference channel repeated across
        plexes) has several values per feature; with the default ``agg=None``
        that raises ``ValueError`` instead of summing them.

        Parameters
        ----------
        level : str
            "protein" (uses pg + run) or "peptide" (uses feature + run).
        value_col : str
            Only ``"intensity"`` (the primary quantity) is supported.
        fillna : float | None
            Value for missing cells; None (default) keeps NaN. Not used with *output_path*.
        output_path : str | Path | None
            Pivot out-of-core with DuckDB into this Parquet file and return its path.
        agg : {None, "mean", "median", "sum", "min", "max"}
            How to combine several values of one cell.
        by : {"sample", "unit"}
            Rows are ``sample_accession``, or ``(sample_accession, quantification_unit)``.
        """
        import pandas as pd

        keys = _design_matrix_keys(value_col, agg, by)
        sql = self._abundance_sql(level)
        if agg is None:
            self._require_one_value_per_cell(sql, keys, level)
        if output_path is not None:
            return self._pivot_to_parquet(sql, keys, _SQL_AGGREGATES[agg or "sum"], Path(output_path))

        df = self._engine.execute(sql).fetchdf()
        if df.empty:
            return pd.DataFrame()
        grouped = df.groupby([*keys, "feature_id"])["intensity"]
        cells = grouped.sum(min_count=1) if agg in (None, "sum") else grouped.agg(agg)
        matrix = cells.unstack("feature_id")
        matrix.columns.name = None
        return matrix if fillna is None else matrix.fillna(fillna)

    def _require_one_value_per_cell(self, sql: str, keys: list[str], level: str) -> None:
        """Raise when a matrix cell would receive more than one value."""
        cells, samples, example = self._engine.execute(
            sql_build(_REPEATED_CELLS_SQL, base_sql=sql, keys=", ".join(keys))
        ).fetchall()[0]
        if not cells:
            return
        if keys == ["sample_accession"] and samples:
            causes = _REPEATED_UNIT_CAUSES + ((_PER_FRACTION_PG_CAUSE,) if level == "protein" else ())
            cause = (
                f"{samples} samples have more than one quantification unit (e.g. {example!r}): "
                f"{', '.join(causes[:-1])}, or {causes[-1]}"
            )
            remedy = ", or by='unit' for one row per quantification unit"
        else:
            cause = f"{cells} cells hold several values within one quantification unit (e.g. groups sharing a leader)"
            remedy = ""
        raise ValueError(
            f"design_matrix(level={level!r}): {cause}. Pass agg='mean' (or 'median', 'sum', 'min', 'max') "
            f"to combine them{remedy}."
        )

    def _pivot_to_parquet(self, sql: str, keys: list[str], aggregate: str, output_path: Path) -> Path:
        """Pivot the long form into *output_path* out-of-core with DuckDB."""
        pivot_sql = sql_build(
            """COPY (
                PIVOT (SELECT $keys, feature_id, intensity FROM ($base_sql))
                ON feature_id USING $aggregate(intensity) GROUP BY $keys
            ) TO '$out_path' (FORMAT PARQUET)""",
            base_sql=sql,
            keys=", ".join(keys),
            aggregate=aggregate,
            out_path=escape_path(str(output_path)),
        )
        self._engine.execute(pivot_sql)
        return output_path

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
