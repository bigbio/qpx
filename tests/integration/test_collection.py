"""DatasetCollection tests."""

import shutil


class TestVirtualCollection:
    def test_sql_across_datasets(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        coll = qpx.DatasetCollection([ds1, ds2])
        result = coll.sql("SELECT COUNT(*) AS cnt FROM feature_0")
        assert result.fetchone()[0] > 0
        coll.close()
        ds1.close()
        ds2.close()

    def test_sql_union_across_datasets(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        coll = qpx.DatasetCollection([ds1, ds2])
        result = coll.sql("SELECT COUNT(*) AS cnt FROM feature_0 UNION ALL SELECT COUNT(*) FROM feature_1")
        rows = result.fetchall()
        assert len(rows) == 2
        assert rows[0][0] == rows[1][0]  # Same data, same count
        coll.close()
        ds1.close()
        ds2.close()

    def test_structure_names(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        coll = qpx.DatasetCollection([ds1, ds2])
        names = coll.structure_names
        assert 0 in names
        assert 1 in names
        assert "feature" in names[0]
        assert "feature" in names[1]
        coll.close()
        ds1.close()
        ds2.close()

    def test_context_manager(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        with qpx.DatasetCollection([ds1, ds2]) as coll:
            result = coll.sql("SELECT COUNT(*) FROM feature_0")
            assert result.fetchone()[0] > 0

        ds1.close()
        ds2.close()


class TestPhysicalMerge:
    def test_merge_creates_output(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)
        merge_dir = tmp_path / "merged"

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        coll = qpx.DatasetCollection([ds1, ds2])
        coll.merge(merge_dir, structures=["feature"])

        # Verify merged file exists and has double rows
        merged_ds = qpx.open_dataset(str(merge_dir))
        assert merged_ds.feature is not None
        assert merged_ds.feature.count() == ds1.feature.count() + ds2.feature.count()

        coll.close()
        ds1.close()
        ds2.close()
        merged_ds.close()

    def test_merge_adds_source_column(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)
        merge_dir = tmp_path / "merged"

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        coll = qpx.DatasetCollection([ds1, ds2])
        coll.merge(merge_dir, structures=["feature"])

        merged_ds = qpx.open_dataset(str(merge_dir))
        df = merged_ds.feature.to_df()
        assert "source_dataset" in df.columns
        assert df["source_dataset"].nunique() == 2

        coll.close()
        ds1.close()
        ds2.close()
        merged_ds.close()

    def test_merge_common_structures(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)
        merge_dir = tmp_path / "merged"

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        coll = qpx.DatasetCollection([ds1, ds2])
        # Merge all common structures
        coll.merge(merge_dir)

        merged_ds = qpx.open_dataset(str(merge_dir))
        # At minimum, feature should be present
        assert merged_ds.feature is not None

        coll.close()
        ds1.close()
        ds2.close()
        merged_ds.close()


class TestShardedCollection:
    """A collection must see every shard the dataset itself sees.

    ``Dataset`` unions all matching shards into its view but keeps ``_file_path``
    pointing at the first one. Registering a collection from that single path
    returned a subset of the rows, with no error (bigbio/qpx#286) — the worst
    failure mode for this data, since the answer is plausible rather than absent.
    """

    @staticmethod
    def _shard(dataset_dir, tmp_path, name):
        """Copy a dataset and split its feature file into two shards."""
        import pyarrow.parquet as pq

        out = tmp_path / name
        shutil.copytree(dataset_dir, out)
        original = next(out.glob("*.feature.parquet"))
        table = pq.read_table(original)
        assert table.num_rows >= 2, "fixture needs at least two feature rows to shard"
        half = table.num_rows // 2
        stem = original.name[: -len(".feature.parquet")]
        pq.write_table(table.slice(0, half), out / f"{stem}.part0.feature.parquet")
        pq.write_table(table.slice(half), out / f"{stem}.part1.feature.parquet")
        original.unlink()
        return out, table.num_rows

    def test_collection_reads_every_shard(self, dataset_dir, tmp_path):
        """A collection reads every feature shard registered by its dataset."""
        import qpx

        ds_dir, total = self._shard(dataset_dir, tmp_path, "sharded")
        ds = qpx.open_dataset(str(ds_dir))
        dataset_rows = ds.sql("SELECT COUNT(*) AS c FROM feature").fetchone()[0]
        assert dataset_rows == total, "dataset itself should already union the shards"

        coll = qpx.DatasetCollection([ds])
        collection_rows = coll.sql("SELECT COUNT(*) AS c FROM feature_0").fetchone()[0]
        coll.close()
        ds.close()

        assert collection_rows == dataset_rows, (
            f"collection saw {collection_rows} of {dataset_rows} rows — it is reading only one shard"
        )

    def test_query_clone_keeps_the_shard_list(self, dataset_dir, tmp_path):
        """A lazily-derived structure must not lose its backing files."""
        import qpx

        ds_dir, _ = self._shard(dataset_dir, tmp_path, "sharded_clone")
        ds = qpx.open_dataset(str(ds_dir))
        clone = ds.feature.limit(1)
        assert clone.file_paths == ds.feature.file_paths
        ds.close()

    def test_s3_locator_is_not_mangled_by_path(self):
        """An s3:// URI must survive as a URI, not become 's3:/...'."""
        from qpx.core.data.base import BaseStructure

        struct = BaseStructure.__new__(BaseStructure)
        BaseStructure.__init__(
            struct,
            engine=None,
            table_name="feature",
            file_path="s3://bucket/ds/feature",
            file_paths=["s3://bucket/ds/*.feature.parquet"],
        )
        assert struct.file_paths == ["s3://bucket/ds/*.feature.parquet"]
        assert str(struct.file_paths[0]).startswith("s3://")
