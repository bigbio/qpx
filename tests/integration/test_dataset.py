"""Tests for QPX Dataset class: auto-discovery, SQL, partial loading, register_external."""

import pandas as pd
import pyarrow as pa
import pytest

from qpx.core.convert import QueryResult
from qpx.core.data.dataset import DatasetMeta
from qpx.core.data.feature import Feature
from qpx.core.data.pg import PG
from qpx.core.data.run import Run
from qpx.core.data.sample import Sample
from qpx.dataset import Dataset
from qpx.writers import FeatureWriter
from tests.conftest import make_feature_record

# ---------------------------------------------------------------------------
# Auto-discovery tests
# ---------------------------------------------------------------------------


class TestDatasetDiscovery:
    """Test that Dataset auto-discovers Parquet files by suffix."""

    def test_discovers_all_structures(self, dataset_dir):
        """Dataset finds all Parquet files matching QPX suffixes."""
        with Dataset(dataset_dir) as ds:
            available = ds.available_structures
            assert "feature" in available
            assert "pg" in available
            assert "sample" in available
            assert "run" in available
            assert "dataset" in available

    def test_feature_accessor(self, dataset_dir):
        """Dataset.feature returns a Feature data structure."""
        with Dataset(dataset_dir) as ds:
            assert ds.feature is not None
            assert isinstance(ds.feature, Feature)

    def test_pg_accessor(self, dataset_dir):
        """Dataset.pg returns a PG data structure."""
        with Dataset(dataset_dir) as ds:
            assert ds.pg is not None
            assert isinstance(ds.pg, PG)

    def test_sample_accessor(self, dataset_dir):
        """Dataset.sample returns a Sample data structure."""
        with Dataset(dataset_dir) as ds:
            assert ds.sample is not None
            assert isinstance(ds.sample, Sample)

    def test_run_accessor(self, dataset_dir):
        """Dataset.run returns a Run data structure."""
        with Dataset(dataset_dir) as ds:
            assert ds.run is not None
            assert isinstance(ds.run, Run)

    def test_dataset_meta_accessor(self, dataset_dir):
        """Dataset.dataset_meta returns a DatasetMeta data structure."""
        with Dataset(dataset_dir) as ds:
            assert ds.dataset_meta is not None
            assert isinstance(ds.dataset_meta, DatasetMeta)

    def test_missing_structure_returns_none(self, dataset_dir):
        """Accessor for a missing structure returns None."""
        with Dataset(dataset_dir) as ds:
            # PSM was not created in dataset_dir
            assert ds.psm is None
            assert ds.mz is None
            assert ds.ontology is None
            assert ds.provenance is None

    def test_feature_count(self, dataset_dir):
        """Feature structure has correct row count."""
        with Dataset(dataset_dir) as ds:
            assert ds.feature.count() == 3

    def test_pg_count(self, dataset_dir):
        """PG structure has correct row count."""
        with Dataset(dataset_dir) as ds:
            assert ds.pg.count() == 2


# ---------------------------------------------------------------------------
# Partial loading tests
# ---------------------------------------------------------------------------


class TestDatasetPartialLoading:
    """Test Dataset with structures= parameter for selective loading."""

    def test_load_only_feature(self, dataset_dir):
        """Loading only 'feature' ignores other structures."""
        with Dataset(dataset_dir, structures=["feature"]) as ds:
            assert ds.feature is not None
            assert ds.pg is None
            assert ds.run is None
            assert ds.sample is None
            assert ds.available_structures == ["feature"]

    def test_load_feature_and_run(self, dataset_dir):
        """Loading 'feature' and 'run' loads both."""
        with Dataset(dataset_dir, structures=["feature", "run"]) as ds:
            assert ds.feature is not None
            assert ds.run is not None
            assert ds.pg is None
            assert set(ds.available_structures) == {"feature", "run"}

    def test_load_nonexistent_structure(self, dataset_dir):
        """Requesting a structure that doesn't exist in the directory is fine (just None)."""
        with Dataset(dataset_dir, structures=["psm"]) as ds:
            assert ds.psm is None
            assert ds.available_structures == []


# ---------------------------------------------------------------------------
# Cross-structure SQL tests
# ---------------------------------------------------------------------------


class TestDatasetSQL:
    """Test Dataset.sql for cross-structure SQL queries."""

    def test_simple_sql(self, dataset_dir):
        """Dataset.sql can query a single structure."""
        with Dataset(dataset_dir) as ds:
            result = ds.sql("SELECT COUNT(*) AS n FROM feature")
            assert isinstance(result, QueryResult)
            row = result.fetchone()
            assert row[0] == 3

    def test_cross_structure_sql(self, dataset_dir):
        """Dataset.sql can join across structures."""
        with Dataset(dataset_dir) as ds:
            result = ds.sql("""
                SELECT f.sequence, r.run_accession
                FROM feature f
                JOIN run r ON f.run_file_name = r.run_file_name
            """)
            df = result.to_df()
            assert "sequence" in df.columns
            assert "run_accession" in df.columns
            assert len(df) == 3

    def test_sql_returns_queryresult(self, dataset_dir):
        """Dataset.sql returns a QueryResult with all materialization options."""
        with Dataset(dataset_dir) as ds:
            result = ds.sql("SELECT * FROM feature LIMIT 1")
            # to_df
            df = result.to_df()
            assert isinstance(df, pd.DataFrame)

    def test_sql_aggregation(self, dataset_dir):
        """Dataset.sql supports aggregation queries."""
        with Dataset(dataset_dir) as ds:
            result = ds.sql("""
                SELECT run_file_name, COUNT(*) AS n_features
                FROM feature
                GROUP BY run_file_name
                ORDER BY run_file_name
            """)
            df = result.to_df()
            assert len(df) == 2  # Two runs
            assert "n_features" in df.columns

    def test_sql_with_pg(self, dataset_dir):
        """Dataset.sql can query PG structure."""
        with Dataset(dataset_dir) as ds:
            result = ds.sql("SELECT anchor_protein FROM pg ORDER BY anchor_protein")
            df = result.to_df()
            assert len(df) == 2
            assert "anchor_protein" in df.columns


# ---------------------------------------------------------------------------
# available_structures tests
# ---------------------------------------------------------------------------


class TestDatasetAvailableStructures:
    """Test Dataset.available_structures property."""

    def test_returns_list(self, dataset_dir):
        """available_structures returns a list."""
        with Dataset(dataset_dir) as ds:
            assert isinstance(ds.available_structures, list)

    def test_contains_discovered_structures(self, dataset_dir):
        """available_structures includes all discovered structures."""
        with Dataset(dataset_dir) as ds:
            available = ds.available_structures
            assert "feature" in available
            assert "pg" in available


# ---------------------------------------------------------------------------
# register_external tests
# ---------------------------------------------------------------------------


class TestDatasetRegisterExternal:
    """Test Dataset.register_external for adding external views."""

    def test_register_external_parquet(self, dataset_dir, tmp_path):
        """register_external adds a queryable view from an external file."""
        # Create an external Parquet file (e.g., mokume output)
        external_path = tmp_path / "external.feature.parquet"
        with FeatureWriter(external_path) as w:
            w.write_batch(
                [
                    make_feature_record(sequence="EXTERNAL_PEP", anchor_protein="PEXT1"),
                ]
            )

        with Dataset(dataset_dir) as ds:
            ds.register_external("mokume_output", external_path)
            result = ds.sql("SELECT COUNT(*) FROM mokume_output")
            count = result.fetchone()[0]
            assert count == 1

    def test_register_external_queryable_with_sql(self, dataset_dir, tmp_path):
        """External views can be joined with internal structures via SQL."""
        external_path = tmp_path / "extra.feature.parquet"
        with FeatureWriter(external_path) as w:
            w.write_batch(
                [
                    make_feature_record(
                        sequence="PEPTIDEK",
                        anchor_protein="P12345",
                        run_file_name="run_01",
                    ),
                ]
            )

        with Dataset(dataset_dir) as ds:
            ds.register_external("extra_features", external_path)
            result = ds.sql("""
                SELECT e.sequence, r.run_accession
                FROM extra_features e
                JOIN run r ON e.run_file_name = r.run_file_name
            """)
            df = result.to_df()
            assert len(df) == 1
            assert df["run_accession"].iloc[0] == "assay_01"


# ---------------------------------------------------------------------------
# Dataset context manager tests
# ---------------------------------------------------------------------------


class TestDatasetContextManager:
    """Test Dataset context manager lifecycle."""

    def test_context_manager_opens_and_closes(self, dataset_dir):
        """Dataset works as a context manager."""
        with Dataset(dataset_dir) as ds:
            assert ds.feature is not None
            assert ds.feature.count() == 3
        # After exit, engine should be closed

    def test_repr(self, dataset_dir):
        """Dataset repr includes path and structure names."""
        with Dataset(dataset_dir) as ds:
            r = repr(ds)
            assert "Dataset" in r
            assert "feature" in r

    def test_custom_duckdb_settings(self, dataset_dir):
        """Dataset accepts custom DuckDB memory/thread settings."""
        with Dataset(dataset_dir, duckdb_memory="1GB", duckdb_threads=1) as ds:
            assert ds.feature.count() == 3


# ---------------------------------------------------------------------------
# Dataset data structure operations through Dataset
# ---------------------------------------------------------------------------


class TestDatasetOperations:
    """Test data structure operations when accessed through Dataset."""

    def test_feature_filter_through_dataset(self, dataset_dir):
        """Feature.filter works when accessed through Dataset."""
        with Dataset(dataset_dir) as ds:
            filtered = ds.feature.filter("charge = 2")
            assert filtered.count() == 2

    def test_feature_to_df_through_dataset(self, dataset_dir):
        """Feature.to_df works when accessed through Dataset."""
        with Dataset(dataset_dir) as ds:
            df = ds.feature.to_df()
            assert isinstance(df, pd.DataFrame)
            assert len(df) == 3

    def test_pg_filter_through_dataset(self, dataset_dir):
        """PG.filter works when accessed through Dataset."""
        with Dataset(dataset_dir) as ds:
            targets = ds.pg.filter("is_decoy = false")
            assert targets.count() == 2

    def test_feature_select_through_dataset(self, dataset_dir):
        """Feature.select works when accessed through Dataset."""
        with Dataset(dataset_dir) as ds:
            selected = ds.feature.select("sequence", "charge")
            df = selected.to_df()
            assert set(df.columns) == {"sequence", "charge"}


# ---------------------------------------------------------------------------
# intensity (lazy) tests
# ---------------------------------------------------------------------------


class TestIntensity:
    """Tests for Dataset.intensity() — lazy long-form query."""

    def test_protein_intensity_returns_queryresult(self, dataset_dir):
        """intensity() returns a QueryResult, not a DataFrame."""
        with Dataset(dataset_dir) as ds:
            result = ds.intensity(level="protein")
            assert isinstance(result, QueryResult)

    def test_protein_intensity_to_df(self, dataset_dir):
        """Can materialize protein intensity as DataFrame."""
        with Dataset(dataset_dir) as ds:
            df = ds.intensity(level="protein").to_df()
            assert "sample_accession" in df.columns
            assert "feature_id" in df.columns
            assert "intensity" in df.columns
            assert len(df) == 2  # 2 protein-sample pairs

    def test_peptide_intensity_to_df(self, dataset_dir):
        """Can materialize peptide intensity as DataFrame."""
        with Dataset(dataset_dir) as ds:
            df = ds.intensity(level="peptide").to_df()
            assert "sample_accession" in df.columns
            assert "feature_id" in df.columns
            assert len(df) > 0

    def test_intensity_to_arrow(self, dataset_dir):
        """Can materialize as Arrow table (more memory-efficient than pandas)."""
        with Dataset(dataset_dir) as ds:
            table = ds.intensity(level="protein").to_arrow()
            assert table.num_rows == 2
            assert "sample_accession" in table.schema.names

    def test_intensity_row_iteration(self, dataset_dir):
        """Can iterate row-by-row for constant memory usage."""
        with Dataset(dataset_dir) as ds:
            rows = list(ds.intensity(level="protein"))
            assert len(rows) == 2
            # Each row is a tuple: (sample_accession, feature_id, intensity)
            assert len(rows[0]) == 3

    def test_intensity_missing_structure_raises(self, tmp_path):
        """Missing structures raise ValueError."""
        ds_dir = tmp_path / "empty_ds"
        ds_dir.mkdir()
        with Dataset(ds_dir) as ds:
            with pytest.raises(ValueError, match="requires"):
                ds.intensity(level="protein")

    def test_intensity_invalid_level_raises(self, dataset_dir):
        """Invalid level raises ValueError."""
        with Dataset(dataset_dir) as ds:
            with pytest.raises(ValueError, match="level must be"):
                ds.intensity(level="invalid")

    def test_legacy_channel_intensities_join_by_sample_accession(self, tmp_path):
        """Legacy feature channel structs do not require a nonexistent
        run.samples.channel. (The pg view is flattened since 1.1 and no longer
        carries a channel-struct intensities list; legacy pg must be regenerated.)"""
        import pyarrow.parquet as pq

        legacy_intensity = pa.list_(
            pa.struct(
                [
                    ("sample_accession", pa.string()),
                    ("channel", pa.string()),
                    ("intensity", pa.float32()),
                ]
            )
        )
        sample_channels = pa.list_(
            pa.struct(
                [
                    ("sample_accession", pa.string()),
                    ("label", pa.string()),
                ]
            )
        )
        intensity = [[{"sample_accession": "S1", "channel": "legacy", "intensity": 42.0}]]

        pq.write_table(
            pa.table(
                {
                    "sequence": ["PEPTIDE"],
                    "run_file_name": ["run_01"],
                    "is_decoy": [False],
                    "intensities": pa.array(intensity, type=legacy_intensity),
                }
            ),
            tmp_path / "legacy.feature.parquet",
        )
        pq.write_table(
            pa.table(
                {
                    "run_file_name": ["run_01"],
                    "samples": pa.array(
                        [[{"sample_accession": "S1", "label": "LFQ"}]],
                        type=sample_channels,
                    ),
                }
            ),
            tmp_path / "legacy.run.parquet",
        )

        with Dataset(tmp_path, file_prefix="legacy") as dataset:
            peptide = dataset.intensity("peptide").to_df()

        assert peptide[["sample_accession", "feature_id", "intensity"]].to_dict("records") == [
            {"sample_accession": "S1", "feature_id": "PEPTIDE", "intensity": 42.0}
        ]


# ---------------------------------------------------------------------------
# design_matrix tests
# ---------------------------------------------------------------------------


class TestDesignMatrix:
    """Tests for Dataset.design_matrix() — in-memory and out-of-core pivot."""

    def test_protein_matrix_shape(self, dataset_dir):
        """Protein matrix: rows=samples, cols=proteins."""
        with Dataset(dataset_dir) as ds:
            matrix = ds.design_matrix(level="protein")
            # dataset_dir has 2 proteins (P12345, P67890) and 2 samples
            assert matrix.shape[0] == 2  # 2 samples
            assert matrix.shape[1] == 2  # 2 proteins
            assert "P12345" in matrix.columns
            assert "P67890" in matrix.columns

    def test_peptide_matrix_shape(self, dataset_dir):
        """Peptide matrix: rows=samples, cols=peptide sequences."""
        with Dataset(dataset_dir) as ds:
            matrix = ds.design_matrix(level="peptide")
            # 3 features, 3 unique sequences, across 2 samples
            assert matrix.shape[0] == 2  # 2 samples
            assert "PEPTIDEK" in matrix.columns

    def test_protein_matrix_values(self, dataset_dir):
        """Verify intensity values are placed correctly."""
        with Dataset(dataset_dir) as ds:
            matrix = ds.design_matrix(level="protein")
            # P12345 in SAMPLE_01 has intensity 5000.0 (from pg)
            assert matrix.loc["SAMPLE_01", "P12345"] == pytest.approx(5000.0)
            # P67890 in SAMPLE_02 has intensity 3000.0 (from pg)
            assert matrix.loc["SAMPLE_02", "P67890"] == pytest.approx(3000.0)

    def test_missing_structure_raises(self, tmp_path):
        """Requesting protein matrix without pg raises ValueError."""
        ds_dir = tmp_path / "empty_ds"
        ds_dir.mkdir()
        with Dataset(ds_dir) as ds:
            with pytest.raises(ValueError, match="requires"):
                ds.design_matrix(level="protein")

    def test_fillna_default_is_zero(self, dataset_dir):
        """Missing values should be filled with 0 by default."""
        with Dataset(dataset_dir) as ds:
            matrix = ds.design_matrix(level="protein", fillna=0.0)
            # SAMPLE_01 has P12345 but not P67890 -> should be 0.0
            assert matrix.loc["SAMPLE_01", "P67890"] == pytest.approx(0.0)

    def test_fillna_nan(self, dataset_dir):
        """fillna=None should leave NaN in place."""
        import math

        with Dataset(dataset_dir) as ds:
            matrix = ds.design_matrix(level="protein", fillna=None)
            assert math.isnan(matrix.loc["SAMPLE_01", "P67890"])

    def test_invalid_level_raises(self, dataset_dir):
        """Invalid level parameter raises ValueError."""
        with Dataset(dataset_dir) as ds:
            with pytest.raises(ValueError, match="level must be"):
                ds.design_matrix(level="invalid")

    def test_output_path_writes_parquet(self, dataset_dir, tmp_path):
        """output_path writes pivot directly to Parquet (out-of-core)."""
        import pyarrow.parquet as pq

        with Dataset(dataset_dir) as ds:
            out = tmp_path / "matrix.parquet"
            result = ds.design_matrix(level="protein", output_path=out)
            assert result == out
            assert out.exists()

            table = pq.read_table(out)
            assert table.num_rows == 2  # 2 samples
            assert "sample_accession" in table.schema.names
            assert "P12345" in table.schema.names

    def test_output_path_peptide(self, dataset_dir, tmp_path):
        """output_path works for peptide level too."""
        import pyarrow.parquet as pq

        with Dataset(dataset_dir) as ds:
            out = tmp_path / "peptide_matrix.parquet"
            ds.design_matrix(level="peptide", output_path=out)
            assert out.exists()

            table = pq.read_table(out)
            assert table.num_rows > 0
            assert "sample_accession" in table.schema.names


# ---------------------------------------------------------------------------
# save_structure tests
# ---------------------------------------------------------------------------


class TestSaveStructure:
    """Tests for Dataset.save_structure()."""

    def test_save_feature_records(self, dataset_dir):
        """Save a list of dicts as feature parquet."""
        import pyarrow.parquet as pq

        with Dataset(dataset_dir) as ds:
            records = [
                make_feature_record(
                    sequence="NEWPEP",
                    anchor_protein="P99999",
                    run_file_name="run_03",
                ),
            ]
            path = ds.save_structure(records, "feature", prefix="new")
            assert path.exists()
            assert path.name == "new.feature.parquet"

            table = pq.read_table(path)
            assert table.num_rows == 1
            assert table.column("sequence")[0].as_py() == "NEWPEP"

    def test_save_dataframe(self, dataset_dir):
        """Save a pandas DataFrame as pg parquet."""
        import pyarrow.parquet as pq

        with Dataset(dataset_dir) as ds:
            df = ds.pg.to_df()
            path = ds.save_structure(df, "pg", prefix="copy")
            assert path.exists()
            assert path.name == "copy.pg.parquet"

            table = pq.read_table(path)
            assert table.num_rows == df.shape[0]

    def test_save_arrow_table(self, dataset_dir):
        """Save an Arrow table as sample parquet."""
        import pyarrow.parquet as pq

        with Dataset(dataset_dir) as ds:
            table = ds.sample.to_arrow()
            path = ds.save_structure(table, "sample", prefix="backup")
            assert path.exists()

            reloaded = pq.read_table(path)
            assert reloaded.num_rows == table.num_rows

    def test_save_unknown_structure_raises(self, dataset_dir):
        """Unknown structure name raises ValueError."""
        with Dataset(dataset_dir) as ds:
            with pytest.raises(ValueError, match="Unknown structure"):
                ds.save_structure([], "not_a_structure")


# ---------------------------------------------------------------------------
# save_anndata tests
# ---------------------------------------------------------------------------


class TestSaveAnndata:
    """Tests for Dataset.save_anndata()."""

    def test_save_anndata_basic(self, dataset_dir):
        """Save an AnnData object to the dataset directory."""
        anndata = pytest.importorskip("anndata")
        import numpy as np

        with Dataset(dataset_dir) as ds:
            adata = anndata.AnnData(
                X=np.array([[1.0, 2.0], [3.0, 4.0]]),
                obs=pd.DataFrame({"sample": ["S1", "S2"]}),
                var=pd.DataFrame({"protein": ["P1", "P2"]}),
            )
            path = ds.save_anndata(adata, "de_results")
            assert path.exists()
            assert path.name == "de_results.h5ad"
            assert path.parent == ds.path

            loaded = anndata.read_h5ad(path)
            assert loaded.shape == (2, 2)
            assert loaded.uns["qpx_version"] == "1.1"
            assert loaded.uns["writer_version"]

    def test_save_anndata_custom_name(self, dataset_dir):
        """Save with a custom file name."""
        anndata = pytest.importorskip("anndata")
        import numpy as np

        with Dataset(dataset_dir) as ds:
            adata = anndata.AnnData(X=np.array([[1.0]]))
            path = ds.save_anndata(adata, "ae_output")
            assert path.name == "ae_output.h5ad"

    def test_save_anndata_overwrites(self, dataset_dir):
        """Writing same name overwrites the file."""
        anndata = pytest.importorskip("anndata")
        import numpy as np

        with Dataset(dataset_dir) as ds:
            adata1 = anndata.AnnData(X=np.array([[1.0]]))
            adata2 = anndata.AnnData(X=np.array([[1.0, 2.0]]))
            ds.save_anndata(adata1, "test")
            path = ds.save_anndata(adata2, "test")

            loaded = anndata.read_h5ad(path)
            assert loaded.shape == (1, 2)


# ---------------------------------------------------------------------------
# refresh tests
# ---------------------------------------------------------------------------


class TestRefresh:
    """Tests for Dataset.refresh()."""

    def test_refresh_picks_up_new_file(self, dataset_dir):
        """After writing a new structure, refresh() discovers it."""
        from qpx.writers import PsmWriter
        from tests.conftest import make_psm_record

        with Dataset(dataset_dir) as ds:
            assert ds.psm is None

            with PsmWriter(dataset_dir / "exp.psm.parquet") as w:
                w.write_batch([make_psm_record()])

            assert ds.psm is None  # Not yet refreshed

            ds.refresh()
            assert ds.psm is not None
            df = ds.psm.to_df()
            assert len(df) == 1

    def test_refresh_clears_view_caches(self, dataset_dir):
        """Cached views are invalidated after refresh."""
        with Dataset(dataset_dir) as ds:
            _ = ds.run_summary
            assert hasattr(ds, "_run_summary")

            ds.refresh()
            assert not hasattr(ds, "_run_summary")

    def test_refresh_updates_available_structures(self, dataset_dir):
        """available_structures list updates after refresh."""
        from qpx.writers import ProvenanceWriter
        from tests.conftest import make_provenance_record

        with Dataset(dataset_dir) as ds:
            assert "provenance" not in ds.available_structures

            with ProvenanceWriter(dataset_dir / "exp.provenance.parquet") as w:
                w.write_batch([make_provenance_record()])

            ds.refresh()
            assert "provenance" in ds.available_structures


# ---------------------------------------------------------------------------
# Integration: mokume workflow
# ---------------------------------------------------------------------------


class TestMokumeWorkflow:
    """Integration test — simulates the mokume external tool workflow."""

    def test_full_workflow(self, dataset_dir):
        """Open dataset -> design_matrix -> mock analysis -> save_anndata -> refresh."""
        anndata = pytest.importorskip("anndata")
        import numpy as np

        with Dataset(dataset_dir) as ds:
            assert len(ds.available_structures) >= 4

            # Get design matrix for protein-level analysis
            matrix = ds.design_matrix(level="protein")
            assert matrix.shape[0] > 0
            assert matrix.shape[1] > 0

            # Mock differential expression analysis
            n_samples, n_proteins = matrix.shape
            adata = anndata.AnnData(
                X=matrix.values,
                obs=pd.DataFrame({"sample_accession": matrix.index.tolist()}),
                var=pd.DataFrame({"protein_accession": matrix.columns.tolist()}),
            )
            adata.uns["comparisons"] = ["SAMPLE_01_vs_SAMPLE_02"]
            adata.varm["log2fc"] = np.random.randn(n_proteins, 1)

            # Write results back
            path = ds.save_anndata(adata, "de_results")
            assert path.exists()
            assert path.suffix == ".h5ad"

            # Verify roundtrip
            loaded = anndata.read_h5ad(path)
            assert loaded.shape == (n_samples, n_proteins)
            assert "comparisons" in loaded.uns

            # Also save validated Parquet back
            pg_df = ds.pg.to_df()
            pg_path = ds.save_structure(pg_df, "pg", prefix="enriched")
            assert pg_path.exists()

            # Refresh and verify
            ds.refresh()
            assert len(ds.available_structures) >= 4

    def test_register_external_and_sql(self, dataset_dir, tmp_path):
        """External tool writes a file, registers it, runs cross-structure SQL."""
        import pyarrow.parquet as pq

        with Dataset(dataset_dir) as ds:
            analysis_table = pa.table(
                {
                    "anchor_protein": ["P12345", "P67890"],
                    "log2fc": [1.5, -0.8],
                    "pvalue": [0.001, 0.05],
                }
            )
            analysis_path = tmp_path / "de_result.parquet"
            pq.write_table(analysis_table, str(analysis_path))

            ds.register_external("de_result", analysis_path)

            result = ds.sql("""
                SELECT pg.anchor_protein, pg.global_qvalue, d.log2fc, d.pvalue
                FROM pg
                JOIN de_result d ON pg.anchor_protein = d.anchor_protein
                WHERE d.pvalue < 0.05
            """)
            df = result.to_df()
            assert len(df) > 0
            assert "log2fc" in df.columns


class TestMultiShardLocalRead:
    """bigbio/qpx#252: a local structure split across several shards must read
    ALL of them (unioned), matching the S3 path — previously only the first
    matching file was used, so local and S3 returned different rows."""

    def test_local_reads_all_shards(self, tmp_path):
        # Two feature shards with no file_prefix, so the '*.feature.parquet'
        # glob matches both.
        with FeatureWriter(tmp_path / "a.feature.parquet") as w:
            w.write_batch([make_feature_record(peptidoform="PEPTIDEK", run_file_name="run_01")])
        with FeatureWriter(tmp_path / "b.feature.parquet") as w:
            w.write_batch(
                [
                    make_feature_record(peptidoform="ELVISK", run_file_name="run_02"),
                    make_feature_record(peptidoform="SAMPLERK", run_file_name="run_03"),
                ]
            )

        with Dataset(tmp_path) as ds:
            assert ds.feature is not None
            # All three rows across both shards, not just the first shard's one.
            assert ds.feature.count() == 3
