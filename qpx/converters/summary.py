"""Conversion summary — a concise report of a produced QPX dataset.

Emitted at the end of every ``qpxc convert <tool>`` run so the operator sees, at
a glance, what the conversion produced: row counts for the PSM / feature / pg
views, the number of distinct peptidoforms (and, when cheap, proteins and
genes), plus feature->pg link diagnostics built on the computed softlink
(:meth:`qpx.dataset.Dataset.link_feature_pg`).

The link diagnostics count features with no matching protein-group link
(:meth:`~qpx.dataset.Dataset.features_without_pg_link`). Link presence does not
indicate whether a feature or its protein group has a quantified intensity.

Everything is computed with DuckDB ``COUNT`` / ``COUNT(DISTINCT ...)`` over the
registered views (fast, and identical on the local and S3 read paths); whole
tables are never pulled into Python. Each metric is guarded for the view being
absent, since a conversion may produce only some structures.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from qpx.dataset import Dataset

_log = logging.getLogger(__name__)


def _scalar(dataset: Dataset, sql: str) -> int | None:
    """Run a scalar aggregate query, returning ``None`` on any failure.

    Used so a single missing/renamed column (e.g. ``gg_accessions`` on an older
    file) degrades that one metric to ``None`` instead of aborting the summary.
    """
    try:
        row = dataset.sql(sql).fetchone()
    except Exception as exc:  # noqa: BLE001 - a summary must never raise  # pylint: disable=broad-exception-caught
        _log.debug("conversion summary metric failed for %r: %s", sql, exc)
        return None
    if not row or row[0] is None:
        return 0
    return int(row[0])


def conversion_summary(dataset: Dataset | str | Path) -> dict:
    """Compute a conversion summary for a QPX dataset.

    Parameters
    ----------
    dataset:
        An open :class:`~qpx.dataset.Dataset`, or a path (local folder or
        ``s3://`` prefix) that is opened for the duration of the summary.

    Returns
    -------
    dict
        Keys (any may be ``None`` when the backing view/column is absent):
        ``accession``, ``n_psms``, ``n_features``, ``n_protein_groups``,
        ``n_peptidoforms``, ``n_proteins``, ``n_genes``; only when BOTH the
        ``feature`` and ``pg`` views exist — ``n_feature_pg_links``,
        ``n_features_linked`` and ``n_features_without_pg``; and only when BOTH
        the ``feature`` and ``psm`` views exist — ``n_features_with_psm`` and
        ``n_psms_without_feature``.
    """
    ds, opened = _as_dataset(dataset)
    try:
        return _collect(ds)
    finally:
        if opened:
            ds.close()


def _as_dataset(dataset: Dataset | str | Path):
    """Return ``(dataset, opened)`` — opening a path if one was given."""
    from qpx.dataset import Dataset as _Dataset

    if isinstance(dataset, (str, Path)):
        return _Dataset(dataset), True
    return dataset, False


def _collect(ds: Dataset) -> dict:
    summary: dict = {
        "accession": None,
        "n_psms": None,
        "n_features": None,
        "n_protein_groups": None,
        "n_peptidoforms": None,
        "n_proteins": None,
        "n_genes": None,
        "n_feature_pg_links": None,
        "n_features_linked": None,
        "n_features_without_pg": None,
        "n_features_with_psm": None,
        "n_psms_without_feature": None,
    }

    if ds.dataset_meta is not None:
        row = None
        try:
            row = ds.sql("SELECT project_accession FROM dataset LIMIT 1").fetchone()
        except Exception as exc:  # noqa: BLE001  # pylint: disable=broad-exception-caught
            _log.debug("could not read project_accession: %s", exc)
        if row and row[0]:
            summary["accession"] = str(row[0])

    if ds.psm is not None:
        summary["n_psms"] = _scalar(ds, "SELECT COUNT(*) FROM psm")
    if ds.feature is not None:
        summary["n_features"] = _scalar(ds, "SELECT COUNT(*) FROM feature")
    if ds.pg is not None:
        summary["n_protein_groups"] = _scalar(ds, "SELECT COUNT(*) FROM pg")

    # Peptidoforms: distinct across feature, falling back to psm.
    if ds.feature is not None:
        summary["n_peptidoforms"] = _scalar(ds, "SELECT COUNT(DISTINCT peptidoform) FROM feature")
    elif ds.psm is not None:
        summary["n_peptidoforms"] = _scalar(ds, "SELECT COUNT(DISTINCT peptidoform) FROM psm")

    # Proteins / genes: distinct accessions unnested from the pg membership lists.
    if ds.pg is not None:
        summary["n_proteins"] = _scalar(
            ds,
            "SELECT COUNT(DISTINCT acc) FROM (SELECT UNNEST(pg_accessions) AS acc FROM pg)",
        )
        summary["n_genes"] = _scalar(
            ds,
            "SELECT COUNT(DISTINCT gg) FROM (SELECT UNNEST(gg_accessions) AS gg FROM pg)",
        )

    # Link diagnostics — feature->pg and feature<->psm softlink counts, computed
    # by the public Dataset diagnostics API (each group is guarded to None when the
    # backing views are absent, so copying every key is safe).
    _collect_link_diagnostics(ds, summary)

    return summary


def _collect_link_diagnostics(ds: Dataset, summary: dict) -> None:
    """Populate feature<->pg and feature<->psm diagnostics via the public API.

    Delegates to :meth:`qpx.dataset.Dataset.feature_link_diagnostics`, which
    registers the softlink views and computes the counts with constant-literal
    queries. A failure degrades the diagnostics to ``None`` without aborting the
    summary.
    """
    try:
        diagnostics = ds.feature_link_diagnostics()
    except Exception as exc:  # noqa: BLE001 - a summary must never raise  # pylint: disable=broad-exception-caught
        _log.debug("conversion summary link diagnostics failed: %s", exc)
        return
    for key in (
        "n_feature_pg_links",
        "n_features_linked",
        "n_features_without_pg",
        "n_features_with_psm",
        "n_psms_without_feature",
    ):
        summary[key] = diagnostics.get(key)


def _fmt(value: int | None) -> str:
    return f"{value:,}" if value is not None else "n/a"


def format_conversion_summary(summary: dict) -> str:
    """Render a :func:`conversion_summary` dict as a compact, aligned block."""
    accession = summary.get("accession")
    header = f"QPX conversion summary ({accession}):" if accession else "QPX conversion summary:"

    lines = [header]
    label_width = len("feature->pg links")

    def add(label: str, value: int | None) -> None:
        if value is None:
            return
        lines.append(f"  {label:<{label_width}}: {_fmt(value)}")

    add("psms", summary.get("n_psms"))
    add("features", summary.get("n_features"))
    add("protein groups", summary.get("n_protein_groups"))
    add("peptidoforms", summary.get("n_peptidoforms"))
    add("proteins", summary.get("n_proteins"))
    add("genes", summary.get("n_genes"))

    n_links = summary.get("n_feature_pg_links")
    if n_links is not None:
        linked = summary.get("n_features_linked")
        without = summary.get("n_features_without_pg")
        detail = []
        if linked is not None:
            detail.append(f"{_fmt(linked)} features linked")
        if without is not None:
            detail.append(f"{_fmt(without)} features without a pg link")
        suffix = f" ({'; '.join(detail)})" if detail else ""
        lines.append(f"  {'feature->pg links':<{label_width}}: {_fmt(n_links)}{suffix}")

    with_psm = summary.get("n_features_with_psm")
    if with_psm is not None:
        without_feature = summary.get("n_psms_without_feature")
        suffix = f" ({_fmt(without_feature)} psms without a feature)" if without_feature is not None else ""
        lines.append(f"  {'feature<-psm':<{label_width}}: {_fmt(with_psm)} features with >=1 psm{suffix}")

    return "\n".join(lines)


def log_conversion_summary(output_folder: Path | str, logger: logging.Logger | None = None) -> None:
    """Open ``output_folder`` and log its conversion summary at INFO.

    Wrapped so a summary failure NEVER fails an otherwise-successful conversion:
    any exception is logged as a warning and swallowed. Intended to be called by
    each ``qpxc convert`` subcommand after it writes its output folder.
    """
    log = logger or _log
    try:
        summary = conversion_summary(output_folder)
        log.info("%s", format_conversion_summary(summary))
    except Exception as exc:  # noqa: BLE001 - never break a good conversion  # pylint: disable=broad-exception-caught
        log.warning("Could not build conversion summary for %s: %s", output_folder, exc)
