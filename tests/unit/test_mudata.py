"""Regression tests for qpx.mudata._attach_uns_metadata.

Previously _attach_uns_metadata only filtered None and float NaN; it missed
pd.NA (from pandas nullable dtypes when reading parquet). That caused
mdata.write() to fail on real data with:
    IORegistryError: No method registered for writing <class 'pandas._libs.missing.NAType'>

These tests ensure all scalar NA flavors are skipped while normal values
(including non-scalar dicts/lists) are preserved.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

mudata = pytest.importorskip("mudata")
anndata = pytest.importorskip("anndata")

from qpx.converters.orchestrator import BaseOrchestrator  # noqa: E402
from qpx.dataset import Dataset  # noqa: E402
from qpx.mudata import (  # noqa: E402
    _LABEL_FIELD_QUERIES,
    _attach_uns_metadata,
    _detect_intensity_label,
    _detect_label_field,
    build_mudata,
)
from qpx.version import QPX_SPEC_VERSION  # noqa: E402
from qpx.writers import FeatureWriter, PgWriter, RunWriter  # noqa: E402
from tests.conftest import make_feature_record, make_pg_record, make_run_record  # noqa: E402


class _FakeResult:
    def __init__(self, df: pd.DataFrame) -> None:
        self._df = df

    def fetchdf(self) -> pd.DataFrame:
        return self._df


class _FakeEngine:
    """Minimal stand-in for qpx.core.engine.DuckDBEngine used by _attach_uns_metadata.

    _attach_uns_metadata only calls engine.execute(sql).fetchdf(), so this shim
    is sufficient and avoids spinning up a real DuckDB + parquet fixture.
    """

    def __init__(self, df: pd.DataFrame) -> None:
        self._df = df

    def execute(self, sql: str) -> _FakeResult:  # noqa: ARG002
        return _FakeResult(self._df)


class _FakeRowResult:
    """Returns a single-row fetchone() result."""

    def __init__(self, value) -> None:
        self._value = value

    def fetchone(self):
        return (self._value,) if self._value is not None else None

    def fetchdf(self) -> pd.DataFrame:
        return pd.DataFrame()


class _TableAwareEngine:
    """Fake engine that dispatches based on SQL table references."""

    def __init__(self, table_labels: dict[str, str]) -> None:
        self._table_labels = table_labels

    def execute(self, sql: str, params=None) -> _FakeRowResult:  # noqa: ARG002
        for table, label in self._table_labels.items():
            if table in sql:
                if "typeof" in sql.lower():
                    return _FakeRowResult('STRUCT("label" VARCHAR, intensity FLOAT)[]')
                return _FakeRowResult(label)
        return _FakeRowResult(None)


def _make_empty_mdata() -> "mudata.MuData":
    return mudata.MuData({"precursors": anndata.AnnData(X=np.zeros((1, 1)))})


def _patch_export_dataset(monkeypatch, *, feature, pg):
    """Install a minimal Dataset double and return its close tracker."""
    closed = []
    dataset = SimpleNamespace(
        psm=object(),
        feature=feature,
        pg=pg,
        close=lambda: closed.append(True),
    )
    monkeypatch.setattr("qpx.dataset.Dataset", lambda *_args, **_kwargs: dataset)
    return closed


def _write_mudata(tmp_path, prefix):
    """Invoke the automatic MuData export under test."""
    return BaseOrchestrator()._write_mudata(tmp_path, prefix)


def test_attach_uns_metadata_skips_all_na_flavors():
    """All scalar NA types must be filtered out of mdata.uns.

    Regression for IORegistryError on pd.NA when writing .h5mu.
    """
    df = pd.DataFrame(
        {
            "software_name": ["qpx"],
            "software_version": ["0.1.0"],
            "file_checksums": pd.array([pd.NA], dtype="object"),
            "file_row_counts": [None],
            "file_sizes_bytes": [np.nan],
            "creation_date": [pd.NaT],
            "total_structures": [5],
        }
    )
    mdata = _make_empty_mdata()

    _attach_uns_metadata(_FakeEngine(df), mdata)

    assert "file_checksums" not in mdata.uns
    assert "file_row_counts" not in mdata.uns
    assert "file_sizes_bytes" not in mdata.uns
    assert "creation_date" not in mdata.uns

    assert mdata.uns["software_name"] == "qpx"
    assert mdata.uns["software_version"] == "0.1.0"
    assert mdata.uns["total_structures"] == 5


def test_attach_uns_metadata_preserves_non_scalar_values():
    """Dicts and lists must not be misclassified as NA by pd.isna()."""
    df = pd.DataFrame(
        {
            "tags": [["a", "b"]],
            "config": [{"key": "value"}],
        }
    )
    mdata = _make_empty_mdata()

    _attach_uns_metadata(_FakeEngine(df), mdata)

    assert mdata.uns["tags"] == ["a", "b"]
    assert mdata.uns["config"] == {"key": "value"}


def test_attach_uns_metadata_always_stamps_qpx_versions():
    """Version identity is present even when dataset.parquet has no row."""
    mdata = _make_empty_mdata()

    _attach_uns_metadata(_FakeEngine(pd.DataFrame()), mdata)

    assert mdata.uns["qpx_version"] == QPX_SPEC_VERSION
    assert isinstance(mdata.uns["writer_version"], str)
    assert mdata.uns["writer_version"]


def test_attach_uns_metadata_allows_hdf5_write(tmp_path):
    """End-to-end: mdata.write() must succeed when dataset row contains pd.NA.

    This is the symptom we originally observed on real DIA-NN conversion output.
    """
    df = pd.DataFrame(
        {
            "software_name": ["qpx"],
            "file_checksums": pd.array([pd.NA], dtype="object"),
        }
    )
    mdata = _make_empty_mdata()
    _attach_uns_metadata(_FakeEngine(df), mdata)

    out = tmp_path / "roundtrip.h5mu"
    mdata.write(out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_label_field_query_generates_table_specific_sql():
    """_LABEL_FIELD_QUERIES must embed the table name, not hardcode 'feature'."""
    sql_feat = _LABEL_FIELD_QUERIES[("label", "feature")]
    sql_pg = _LABEL_FIELD_QUERIES[("label", "pg")]
    assert "feature" in sql_feat and "pg" not in sql_feat
    assert "pg" in sql_pg and "feature" not in sql_pg


def test_detect_intensity_label_per_table():
    """DIA-NN feature uses 'raw', pg uses 'LFQ' — each must be detected independently.

    Regression: previously _detect_intensity_label only inspected the feature table,
    causing _build_protein_adata to query pg WHERE label='raw' and get zero rows.
    """
    engine = _TableAwareEngine({"feature": "raw", "pg": "LFQ"})
    assert _detect_intensity_label(engine, "feature") == "raw"
    assert _detect_intensity_label(engine, "pg") == "LFQ"


def _write_quant_bundle(tmp_path, prefix, feature_intensities, protein_intensities):
    with FeatureWriter(tmp_path / f"{prefix}.feature.parquet") as writer:
        writer.write_batch([make_feature_record(intensities=feature_intensities)])
    with PgWriter(tmp_path / f"{prefix}.pg.parquet") as writer:
        writer.write_batch([make_pg_record(intensities=protein_intensities)])

    run = make_run_record()
    run["samples"] = [
        {
            "sample_accession": f"{prefix}_{entry['label']}",
            "label": entry["label"],
            "biological_replicate": index,
            "technical_replicate": 1,
        }
        for index, entry in enumerate(feature_intensities, start=1)
    ]
    with RunWriter(tmp_path / f"{prefix}.run.parquet") as writer:
        writer.write_batch([run])


def test_orchestrator_mudata_exports_all_channels_from_requested_prefix(tmp_path):
    """Automatic MuData export must include every channel and never mix prefixes."""
    _write_quant_bundle(
        tmp_path,
        "a",
        [{"label": "TMT126", "intensity": 900.0}],
        [{"label": "TMT126", "intensity": 9000.0}],
    )
    _write_quant_bundle(
        tmp_path,
        "z",
        [
            {"label": "TMT126", "intensity": 100.0},
            {"label": "TMT127N", "intensity": 200.0},
        ],
        [
            {"label": "TMT126", "intensity": 1000.0},
            {"label": "TMT127N", "intensity": 2000.0},
        ],
    )

    output = _write_mudata(tmp_path, "z")
    assert output == tmp_path / "z.h5mu"

    mdata = mudata.read_h5mu(output)
    precursor = mdata.mod["precursors"]
    protein = mdata.mod["proteins"]
    assert mdata.uns["qpx_version"] == QPX_SPEC_VERSION
    assert mdata.uns["writer_version"]
    assert list(precursor.obs_names) == ["run_01|TMT126", "run_01|TMT127N"]
    assert list(precursor.obs["sample_accession"]) == ["z_TMT126", "z_TMT127N"]
    np.testing.assert_allclose(precursor.X.toarray()[:, 0], [100.0, 200.0])
    np.testing.assert_allclose(protein.X.toarray()[:, 0], [1000.0, 2000.0])
    assert not (tmp_path / "z.tmp.h5mu").exists()


def test_orchestrator_mudata_expected_failure_is_best_effort(tmp_path, monkeypatch):
    """Expected assembly failures must not fail the primary conversion."""
    closed = _patch_export_dataset(monkeypatch, feature=object(), pg=None)

    def fail_build(dataset, modalities=None, all_intensity_labels=False):
        _ = dataset, modalities, all_intensity_labels
        raise ValueError("invalid MuData input")

    monkeypatch.setattr("qpx.mudata.build_mudata", fail_build)

    output = _write_mudata(tmp_path, "broken")

    assert output is None
    assert closed == [True]


def test_orchestrator_mudata_write_failure_leaves_no_artifact(tmp_path, monkeypatch):
    """An interrupted HDF5 write must leave neither a final nor temporary file."""
    closed = _patch_export_dataset(monkeypatch, feature=object(), pg=None)

    def partial_write(path):
        """Write a partial file before simulating an HDF5 failure."""
        Path(path).write_bytes(b"partial")
        raise OSError("simulated write failure")

    partial_mdata = SimpleNamespace(
        mod={"precursors": object()},
        write=partial_write,
    )
    monkeypatch.setattr("qpx.mudata.build_mudata", lambda *_args, **_kwargs: partial_mdata)
    (tmp_path / "broken.h5mu").write_bytes(b"stale")

    output = _write_mudata(tmp_path, "broken")

    assert output is None
    assert not (tmp_path / "broken.h5mu").exists()
    assert not (tmp_path / "broken.tmp.h5mu").exists()
    assert closed == [True]


def test_orchestrator_mudata_skips_psm_only_dataset(tmp_path, monkeypatch):
    """Automatic export must not create an empty MuData for PSM-only output."""
    closed = _patch_export_dataset(monkeypatch, feature=None, pg=None)

    def fail_if_called(*_args, **_kwargs):
        pytest.fail("build_mudata must not run without feature or pg input")

    monkeypatch.setattr("qpx.mudata.build_mudata", fail_if_called)

    output = _write_mudata(tmp_path, "psm-only")

    assert output is None
    assert not (tmp_path / "psm-only.h5mu").exists()
    assert closed == [True]


def test_orchestrator_mudata_rejects_missing_required_modality(tmp_path, monkeypatch):
    """A failed core modality must not produce a partial MuData artifact."""
    closed = _patch_export_dataset(monkeypatch, feature=object(), pg=object())

    incomplete_mdata = SimpleNamespace(mod={"proteins": object()})
    monkeypatch.setattr("qpx.mudata.build_mudata", lambda *_args, **_kwargs: incomplete_mdata)

    output = _write_mudata(tmp_path, "incomplete")

    assert output is None
    assert not (tmp_path / "incomplete.h5mu").exists()
    assert not (tmp_path / "incomplete.tmp.h5mu").exists()
    assert closed == [True]


def test_build_mudata_explicit_label_keeps_run_level_contract(tmp_path):
    """Explicit single-label export keeps run observations and matches sample metadata."""
    _write_quant_bundle(
        tmp_path,
        "z",
        [
            {"label": "TMT126", "intensity": 100.0},
            {"label": "TMT127N", "intensity": 200.0},
        ],
        [
            {"label": "TMT126", "intensity": 1000.0},
            {"label": "TMT127N", "intensity": 2000.0},
        ],
    )

    dataset = Dataset(tmp_path, file_prefix="z")
    try:
        mdata = build_mudata(
            dataset,
            intensity_label="TMT127N",
            modalities=["precursors"],
        )
    finally:
        dataset.close()

    precursor = mdata.mod["precursors"]
    assert list(precursor.obs_names) == ["run_01"]
    assert list(precursor.obs["sample_accession"]) == ["z_TMT127N"]
    np.testing.assert_allclose(precursor.X.toarray()[:, 0], [200.0])


def _write_multi_run_bundle(tmp_path, prefix, runs, sample_label=None):
    """Write feature/pg/run parquet for arbitrary run x channel layouts.

    Parameters
    ----------
    runs : list of (run_file_name, [(label, feature_intensity, protein_intensity), ...])
        One entry per MS run, each carrying one or more channels.
    sample_label : str, optional
        Force the ``run.samples[].label`` to this value (used to reproduce the
        label-free case where the intensity label "LFQ" does not equal the SDRF
        sample label "label free sample"). Defaults to the channel label.
    """
    features, pgs, run_records = [], [], []
    for run_name, channels in runs:
        feat_int = [{"label": lab, "intensity": fi} for lab, fi, _ in channels]
        prot_int = [{"label": lab, "intensity": pi} for lab, _, pi in channels]
        features.append(make_feature_record(run_file_name=run_name, intensities=feat_int))
        pgs.append(make_pg_record(run_file_name=run_name, intensities=prot_int))
        samples = [
            {
                "sample_accession": f"{prefix}-{run_name}-{lab}",
                "label": sample_label if sample_label is not None else lab,
                "biological_replicate": idx,
                "technical_replicate": 1,
            }
            for idx, (lab, _, _) in enumerate(channels, start=1)
        ]
        record = make_run_record(run_accession=f"assay_{run_name}", run_file_name=run_name)
        record["samples"] = samples
        run_records.append(record)

    with FeatureWriter(tmp_path / f"{prefix}.feature.parquet") as writer:
        writer.write_batch(features)
    with PgWriter(tmp_path / f"{prefix}.pg.parquet") as writer:
        writer.write_batch(pgs)
    with RunWriter(tmp_path / f"{prefix}.run.parquet") as writer:
        writer.write_batch(run_records)


def test_build_mudata_tmt_default_expands_all_channels(tmp_path):
    """Multiplexed (TMT) data must expand obs over run x channel by default.

    Regression: build_mudata auto-detected ONE label (e.g. TMT126), collapsing
    an 11-plex run to a single channel with every sample mapped to Sample-1.
    """
    _write_multi_run_bundle(
        tmp_path,
        "tmt",
        runs=[
            ("run_01", [("TMT126", 100.0, 1000.0), ("TMT127N", 200.0, 2000.0), ("TMT128N", 300.0, 3000.0)]),
            ("run_02", [("TMT126", 400.0, 4000.0), ("TMT127N", 500.0, 5000.0), ("TMT128N", 600.0, 6000.0)]),
        ],
    )

    dataset = Dataset(tmp_path, file_prefix="tmt")
    try:
        mdata = build_mudata(dataset, modalities=["precursors", "proteins"])
    finally:
        dataset.close()

    precursor = mdata.mod["precursors"]
    protein = mdata.mod["proteins"]

    # obs = runs x channels = 2 x 3 = 6, not 2 collapsed runs.
    assert precursor.n_obs == 6
    assert protein.n_obs == 6
    assert set(precursor.obs["intensity_label"]) == {"TMT126", "TMT127N", "TMT128N"}
    assert set(protein.obs["intensity_label"]) == {"TMT126", "TMT127N", "TMT128N"}

    # Each (run, channel) maps to its own real sample_accession (not all Sample-1).
    assert precursor.obs["sample_accession"].nunique() == 6
    expected = {
        ("run_01", "TMT126"): "tmt-run_01-TMT126",
        ("run_01", "TMT127N"): "tmt-run_01-TMT127N",
        ("run_02", "TMT128N"): "tmt-run_02-TMT128N",
    }
    lookup = {
        (r, lab): acc
        for r, lab, acc in zip(
            precursor.obs["run_file_name"],
            precursor.obs["intensity_label"],
            precursor.obs["sample_accession"],
        )
    }
    for key, acc in expected.items():
        assert lookup[key] == acc


def test_build_mudata_lfq_default_stays_run_level(tmp_path):
    """Label-free data must stay obs = runs with a single LFQ label.

    The intensity label "LFQ" intentionally differs from the SDRF sample label
    "label free sample"; sample_accession must still resolve per run.
    """
    _write_multi_run_bundle(
        tmp_path,
        "lfq",
        runs=[
            ("run_01", [("LFQ", 100.0, 1000.0)]),
            ("run_02", [("LFQ", 200.0, 2000.0)]),
        ],
        sample_label="label free sample",
    )

    dataset = Dataset(tmp_path, file_prefix="lfq")
    try:
        mdata = build_mudata(dataset, modalities=["precursors", "proteins"])
    finally:
        dataset.close()

    precursor = mdata.mod["precursors"]
    assert precursor.n_obs == 2
    assert list(precursor.obs_names) == ["run_01", "run_02"]
    assert list(precursor.obs["sample_accession"]) == ["lfq-run_01-LFQ", "lfq-run_02-LFQ"]


def test_build_mudata_lfq_all_labels_keeps_sample_accession(tmp_path):
    """all_intensity_labels=True (orchestrator path) must not drop LFQ metadata.

    Regression: the intensity label "LFQ" never matches the sample label
    "label free sample", so the all-labels branch produced obs without any
    sample_accession/run_accession columns.
    """
    _write_multi_run_bundle(
        tmp_path,
        "lfq",
        runs=[
            ("run_01", [("LFQ", 100.0, 1000.0)]),
            ("run_02", [("LFQ", 200.0, 2000.0)]),
        ],
        sample_label="label free sample",
    )

    dataset = Dataset(tmp_path, file_prefix="lfq")
    try:
        mdata = build_mudata(dataset, all_intensity_labels=True, modalities=["precursors"])
    finally:
        dataset.close()

    precursor = mdata.mod["precursors"]
    assert "sample_accession" in precursor.obs.columns
    assert list(precursor.obs["sample_accession"]) == ["lfq-run_01-LFQ", "lfq-run_02-LFQ"]


def test_detect_label_field_uses_label_for_current_schema(tmp_path):
    """The current intensities struct is {label, intensity}; queries must use i.label.

    Guards the schema drift where some paths hardcoded i.channel.
    """
    _write_multi_run_bundle(
        tmp_path,
        "tmt",
        runs=[("run_01", [("TMT126", 100.0, 1000.0), ("TMT127N", 200.0, 2000.0)])],
    )

    dataset = Dataset(tmp_path, file_prefix="tmt")
    try:
        assert _detect_label_field(dataset._engine, "feature") == "label"
        assert _detect_label_field(dataset._engine, "pg") == "label"
        # The label queries must actually match rows (non-empty intensity matrix).
        mdata = build_mudata(dataset, modalities=["precursors"])
        assert mdata.mod["precursors"].X.nnz > 0
    finally:
        dataset.close()


def test_build_mudata_feature_mapping_survives_third_modality(tmp_path):
    """bigbio/qpx#252: ``varp['feature_mapping']`` was sized against only
    precursors+proteins, so with a 3rd modality present the assignment raised
    'incorrect shape' and was silently swallowed. The mapping must be built
    against the full global var axis and be present + correct."""
    import scipy.sparse as _sp

    _write_quant_bundle(
        tmp_path,
        "m3",
        [{"label": "LFQ", "intensity": 100.0}],
        [{"label": "LFQ", "intensity": 1000.0}],
    )
    # 3rd modality: an absolute-expression AnnData (<prefix>.pe.h5ad) with its
    # own var axis, disjoint names from the protein accessions.
    expr = anndata.AnnData(
        X=np.array([[1.0, 2.0, 3.0, 4.0]], dtype=float),
        var=pd.DataFrame(index=["G1", "G2", "G3", "G4"]),
    )
    expr.write_h5ad(tmp_path / "m3.pe.h5ad")

    dataset = Dataset(tmp_path, file_prefix="m3")
    try:
        mdata = build_mudata(dataset, modalities=["precursors", "proteins", "expression"])
    finally:
        dataset.close()

    assert "expression" in mdata.mod, "3rd modality must be built for this regression"
    # The axis is genuinely larger than precursors+proteins alone.
    assert mdata.n_vars > mdata.mod["precursors"].n_vars + mdata.mod["proteins"].n_vars

    assert "feature_mapping" in mdata.varp, "feature_mapping must NOT be silently dropped"
    mapping = mdata.varp["feature_mapping"]
    assert mapping.shape == (mdata.n_vars, mdata.n_vars)

    # The PEPTIDEK|2 <-> P12345 precursor/protein link is set at the correct
    # global positions (per-modality offsets on the global axis).
    offsets = {}
    running = 0
    for name, adata in mdata.mod.items():
        offsets[name] = running
        running += adata.n_vars
    prec_pos = offsets["precursors"] + list(mdata.mod["precursors"].var_names).index("PEPTIDEK|2")
    prot_pos = offsets["proteins"] + list(mdata.mod["proteins"].var_names).index("P12345")
    csr = _sp.csr_matrix(mapping)
    assert csr[prec_pos, prot_pos]
    assert csr[prot_pos, prec_pos]
    assert csr.nnz == 2


class TestArrowPivotEquivalence:
    """The Arrow pivot path must produce exactly what the pandas path produced.

    The Arrow path exists only for speed: on a large experiment the pandas path
    hashed one Python string per row (hundreds of millions of them) and copied the
    whole frame to build observation_id. Dictionary-encoding does the same work
    once per distinct value, so these must stay bit-identical.
    """

    @staticmethod
    def _frame():
        return pd.DataFrame(
            {
                "observation_id": ["r2|L", "r1|L", "r1|L", "r3|L", "r2|L"],
                "precursor_id": ["PEP|2", "PEP|2", "PEP|3", "OTHER|2", "OTHER|2"],
                "intensity": [10.0, 20.0, 30.0, 40.0, 50.0],
            }
        )

    def test_matches_pandas_pivot(self):
        """The Arrow pivot produces the same sparse matrix as pandas."""
        import pyarrow as pa

        from qpx.mudata import _pivot_arrow_to_sparse, _pivot_to_sparse

        df = self._frame()
        rows = pd.Index(sorted(df["observation_id"].unique()), name="observation_id")
        cols = pd.Index(sorted(df["precursor_id"].unique()), name="precursor_id")

        expected = _pivot_to_sparse(df, "observation_id", "precursor_id", "intensity", rows, cols)
        actual = _pivot_arrow_to_sparse(pa.Table.from_pandas(df), "observation_id", "precursor_id", "intensity", rows, cols)

        assert actual.shape == expected.shape
        assert (actual != expected).nnz == 0

    def test_unmatched_keys_are_dropped_not_misplaced(self):
        """A value absent from the index must be dropped, exactly as get_indexer does."""
        import pyarrow as pa

        from qpx.mudata import _pivot_arrow_to_sparse, _pivot_to_sparse

        df = self._frame()
        rows = pd.Index(["r1|L", "r2|L"], name="observation_id")  # r3|L deliberately absent
        cols = pd.Index(sorted(df["precursor_id"].unique()), name="precursor_id")

        expected = _pivot_to_sparse(df, "observation_id", "precursor_id", "intensity", rows, cols)
        actual = _pivot_arrow_to_sparse(pa.Table.from_pandas(df), "observation_id", "precursor_id", "intensity", rows, cols)

        assert (actual != expected).nnz == 0
        assert actual.sum() == pytest.approx(110.0)  # 40.0 from the dropped r3 row is gone

    def test_null_keys_are_dropped_before_dictionary_lookup(self):
        """A null Arrow dictionary code is unmatched, not a NumPy index."""
        import pyarrow as pa

        from qpx.mudata import _pivot_arrow_to_sparse, _pivot_to_sparse

        df = pd.DataFrame(
            {
                "observation_id": ["r1|L", None],
                "precursor_id": ["PEP|2", "PEP|2"],
                "intensity": [10.0, 20.0],
            }
        )
        rows = pd.Index(["r1|L"], name="observation_id")
        cols = pd.Index(["PEP|2"], name="precursor_id")

        expected = _pivot_to_sparse(df, "observation_id", "precursor_id", "intensity", rows, cols)
        actual = _pivot_arrow_to_sparse(pa.Table.from_pandas(df), "observation_id", "precursor_id", "intensity", rows, cols)

        assert (actual != expected).nnz == 0
        assert actual.sum() == pytest.approx(10.0)

    def test_sorted_unique_matches_pandas(self):
        """Arrow unique values match pandas sorting and ordering."""
        import pyarrow as pa

        from qpx.mudata import _sorted_unique

        df = self._frame()
        table = pa.Table.from_pandas(df)
        expected = pd.Index(sorted(df["precursor_id"].unique()), name="precursor_id")
        assert list(_sorted_unique(table, "precursor_id")) == list(expected)

    def test_prepare_observations_arrow_matches_pandas(self):
        """Arrow observation preparation matches the pandas implementation."""
        import pyarrow as pa

        from qpx.mudata import _prepare_observations, _prepare_observations_arrow

        df = pd.DataFrame(
            {
                "run_file_name": ["r2", "r1", "r1"],
                "intensity_label": ["L", "L", "M"],
                "precursor_id": ["PEP|2", "PEP|2", "PEP|3"],
                "intensity": [1.0, 2.0, 3.0],
            }
        )
        _, pandas_col, pandas_index, pandas_keys = _prepare_observations(df, None)
        _, arrow_col, arrow_index, arrow_keys = _prepare_observations_arrow(pa.Table.from_pandas(df), None)

        assert arrow_col == pandas_col
        assert list(arrow_index) == list(pandas_index)
        assert list(arrow_keys["run_file_name"]) == list(pandas_keys["run_file_name"])
        assert list(arrow_keys["intensity_label"]) == list(pandas_keys["intensity_label"])
