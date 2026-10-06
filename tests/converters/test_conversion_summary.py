"""Tests for the conversion summary report (qpx.converters.summary).

Runs a real DIA-NN conversion on the shipped small example, then checks the
summary counts and the feature->pg link diagnostics built on the computed
softlink. Also checks that the diagnostics degrade gracefully when the ``pg``
view is absent (``--structures feature`` style).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import qpx
from qpx.converters.summary import (
    conversion_summary,
    format_conversion_summary,
    log_conversion_summary,
)

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples" / "diann" / "small"

_REPORT = EXAMPLES_DIR / "diann_report.tsv"
_SDRF = EXAMPLES_DIR / "PXD019909-DIA.sdrf.tsv"
_PG_MATRIX = EXAMPLES_DIR / "diann_report.pg_matrix.tsv"
_MZML_TSV_DIR = EXAMPLES_DIR / "mzml"


def _parse_scan_number(spectrum_id: str) -> int:
    m = re.search(r"scan=(\d+)", str(spectrum_id))
    return int(m.group(1)) if m else 0


def _prepare_ms_info_parquet(source_dir: Path, dest_dir: Path) -> None:
    """Convert the shipped *_mzml_info.tsv files into *_ms_info.parquet files.

    The feature adapter expects ``rt``, ``scan`` and ``precursor_mz`` columns;
    the shipped TSVs use different names, so bridge them here.
    """
    import pandas as pd

    for tsv_path in source_dir.glob("*_mzml_info.tsv"):
        df = pd.read_csv(tsv_path, sep="\t")
        out = pd.DataFrame(
            {
                "rt": df["Retention_Time"].astype(float),
                "scan": df["SpectrumID"].apply(_parse_scan_number).astype(int),
                "precursor_mz": pd.to_numeric(df["Exp_Mass_To_Charge"], errors="coerce"),
            }
        )
        stem = tsv_path.stem.replace("_mzml_info", "")
        out.to_parquet(dest_dir / f"{stem}_ms_info.parquet", index=False)


@pytest.fixture(scope="module")
def converted_output(tmp_path_factory):
    """Run a DIA-NN feature + pg conversion once for this module."""
    from qpx.converters.diann.converter import DiaNNConverter

    if not _REPORT.exists():
        pytest.skip(f"Test data not found: {_REPORT}")

    output = tmp_path_factory.mktemp("summary_diann_output")
    mzml_parquet_dir = tmp_path_factory.mktemp("summary_mzml_info")
    _prepare_ms_info_parquet(_MZML_TSV_DIR, mzml_parquet_dir)

    converter = DiaNNConverter(report_path=str(_REPORT), sdrf_path=str(_SDRF))
    converter.convert_features(
        mzml_info_folder=str(mzml_parquet_dir),
        qvalue_threshold=0.05,
        output_folder=str(output),
        output_prefix="diann_test",
    )
    converter.convert_pg(
        pg_matrix_path=str(_PG_MATRIX),
        output_folder=str(output),
        output_prefix="diann_test",
    )
    return output


def test_counts_are_present_and_positive(converted_output):
    ds = qpx.Dataset(str(converted_output))
    try:
        summary = conversion_summary(ds)
    finally:
        ds.close()

    assert summary["n_features"] is not None and summary["n_features"] > 0
    assert summary["n_protein_groups"] is not None and summary["n_protein_groups"] > 0
    assert summary["n_peptidoforms"] is not None and summary["n_peptidoforms"] > 0
    # Proteins/genes are cheap distinct-unnest counts over the pg membership.
    assert summary["n_proteins"] is not None and summary["n_proteins"] > 0


def test_link_diagnostics_are_consistent(converted_output):
    ds = qpx.Dataset(str(converted_output))
    try:
        summary = conversion_summary(ds)
        # Cross-check the distinct linked feature count against the raw softlink.
        link_feature_ids = {fid for fid, _pid, _label in ds.link_feature_pg().fetchall()}
    finally:
        ds.close()

    n_features = summary["n_features"]
    linked = summary["n_features_linked"]
    without = summary["n_features_without_pg"]
    total_links = summary["n_feature_pg_links"]

    assert linked is not None and without is not None and total_links is not None
    # Every linked feature is a real feature; the partition never exceeds n_features.
    assert linked + without <= n_features
    # feature_id is the feature PK, so the partition is exact here.
    assert linked + without == n_features
    # Distinct linked features match the softlink; a link exists per linked feature.
    assert linked == len(link_feature_ids)
    assert total_links >= linked


def test_conversion_summary_accepts_a_path(converted_output):
    """Passing a folder path opens and closes its own Dataset."""
    summary = conversion_summary(str(converted_output))
    assert summary["n_features"] > 0
    assert summary["n_feature_pg_links"] is not None


def test_format_block_renders_expected_lines(converted_output):
    summary = conversion_summary(str(converted_output))
    block = format_conversion_summary(summary)

    assert block.startswith("QPX conversion summary")
    assert "features" in block
    assert "protein groups" in block
    assert "peptidoforms" in block
    assert "feature->pg links" in block
    assert "features without a pg link" in block


def test_diagnostics_degrade_without_pg(converted_output):
    """With only the feature view, link diagnostics are absent (not an error)."""
    ds = qpx.Dataset(str(converted_output), structures=["feature"])
    try:
        summary = conversion_summary(ds)
    finally:
        ds.close()

    assert summary["n_features"] is not None and summary["n_features"] > 0
    assert summary["n_peptidoforms"] is not None and summary["n_peptidoforms"] > 0
    # No pg view -> no protein/gene counts and no link diagnostics.
    assert summary["n_protein_groups"] is None
    assert summary["n_proteins"] is None
    assert summary["n_feature_pg_links"] is None
    assert summary["n_features_linked"] is None
    assert summary["n_features_without_pg"] is None

    block = format_conversion_summary(summary)
    assert "feature->pg links" not in block
    assert "protein groups" not in block


def test_log_summary_never_raises_on_bad_path(caplog):
    """A summary failure is swallowed with a warning, never propagated."""
    import logging

    with caplog.at_level(logging.WARNING):
        log_conversion_summary("/nonexistent/path/does/not/exist")
    # No exception raised; the folder simply has no structures -> empty summary,
    # or a warning is logged. Either way the call must return cleanly.
