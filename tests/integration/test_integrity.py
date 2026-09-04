"""Dataset integrity tests."""


class TestComputeIntegrity:
    def test_computes_checksums(self, dataset_dir):
        import qpx

        ds = qpx.open_dataset(str(dataset_dir))
        result = ds.compute_integrity()
        assert "file_checksums" in result
        assert len(result["file_checksums"]) > 0
        # All checksums are 64-char hex strings (SHA-256)
        for sha in result["file_checksums"].values():
            assert len(sha) == 64
        ds.close()

    def test_computes_row_counts(self, dataset_dir):
        import qpx

        ds = qpx.open_dataset(str(dataset_dir))
        result = ds.compute_integrity()
        for name, count in result["file_row_counts"].items():
            assert count > 0
        ds.close()

    def test_computes_file_sizes(self, dataset_dir):
        import qpx

        ds = qpx.open_dataset(str(dataset_dir))
        result = ds.compute_integrity()
        for name, size in result["file_sizes_bytes"].items():
            assert size > 0
        ds.close()

    def test_packaged_at_is_iso_timestamp(self, dataset_dir):
        from datetime import datetime

        import qpx

        ds = qpx.open_dataset(str(dataset_dir))
        result = ds.compute_integrity()
        # Should be parseable as ISO 8601
        datetime.fromisoformat(result["packaged_at"])
        ds.close()

    def test_total_structures_count(self, dataset_dir):
        import qpx

        ds = qpx.open_dataset(str(dataset_dir))
        result = ds.compute_integrity()
        # Should count all parquet files in the directory
        parquet_count = len(list(dataset_dir.glob("*.parquet")))
        assert result["total_structures"] == parquet_count
        ds.close()


class TestVerifyIntegrity:
    def test_verify_passes_on_valid_dataset(self, dataset_dir):
        import qpx

        ds = qpx.open_dataset(str(dataset_dir))
        integrity = ds.compute_integrity()
        # Write integrity data into dataset.parquet (use "exp" prefix to overwrite)
        meta_df = ds.dataset_meta.to_df()
        meta_dict = meta_df.iloc[0].to_dict()
        meta_dict.update(integrity)
        ds.save_structure([meta_dict], "dataset", prefix="exp")
        ds.refresh()
        result = ds.verify_integrity()
        assert len(result["errors"]) == 0
        ds.close()

    def test_verify_detects_missing_file(self, dataset_dir):
        import qpx

        ds = qpx.open_dataset(str(dataset_dir))
        integrity = ds.compute_integrity()
        # Add a fake file to checksums
        integrity["file_checksums"]["nonexistent.parquet"] = "a" * 64
        meta_df = ds.dataset_meta.to_df()
        meta_dict = meta_df.iloc[0].to_dict()
        meta_dict.update(integrity)
        ds.save_structure([meta_dict], "dataset", prefix="exp")
        ds.refresh()
        result = ds.verify_integrity()
        assert any("nonexistent" in e for e in result["errors"])
        ds.close()

    def test_verify_warns_when_no_checksums(self, dataset_dir):
        import qpx

        # The default dataset.parquet has null integrity fields
        ds = qpx.open_dataset(str(dataset_dir))
        result = ds.verify_integrity()
        assert len(result["warnings"]) > 0
        ds.close()

    def test_verify_errors_without_dataset_meta(self, tmp_path):
        import qpx
        from qpx.writers import FeatureWriter
        from tests.conftest import make_feature_record

        # Create a directory with only a feature file (no dataset.parquet)
        ds_dir = tmp_path / "no_meta"
        ds_dir.mkdir()
        with FeatureWriter(ds_dir / "exp.feature.parquet") as w:
            w.write_batch([make_feature_record()])

        ds = qpx.open_dataset(str(ds_dir))
        result = ds.verify_integrity()
        assert any("No dataset.parquet" in e for e in result["errors"])
        ds.close()


class TestPartitionedIntegrity:
    """Integrity must cover Hive-partitioned structures too.

    compute_integrity() used a non-recursive glob, so files under
    feature/run_file_name=.../part-0.parquet were never checksummed: the record
    claimed to cover data it had never read (bigbio/qpx#287).
    """

    @staticmethod
    def _partition(dataset_dir, tmp_path):
        """Copy a dataset and rewrite its feature file as a partitioned directory."""
        import shutil

        import pyarrow.parquet as pq

        out = tmp_path / "partitioned"
        shutil.copytree(dataset_dir, out)
        original = next(out.glob("*.feature.parquet"))
        table = pq.read_table(original)
        part_dir = out / "feature" / "run_file_name=run_a"
        part_dir.mkdir(parents=True)
        pq.write_table(table, part_dir / "part-0.parquet")
        original.unlink()
        return out, (part_dir / "part-0.parquet")

    def test_partitioned_files_are_checksummed(self, dataset_dir, tmp_path):
        import qpx

        out, part_file = self._partition(dataset_dir, tmp_path)
        ds = qpx.open_dataset(str(out))
        result = ds.compute_integrity()
        ds.close()

        key = part_file.relative_to(out).as_posix()
        assert key in result["file_checksums"], f"partitioned file not checksummed; keys were {sorted(result['file_checksums'])}"
        assert result["file_sizes_bytes"][key] == part_file.stat().st_size
        assert result["file_row_counts"][key] > 0

    def test_keys_are_root_relative_posix_paths(self, dataset_dir, tmp_path):
        """Keys must be relative to the dataset root so the record stays portable."""
        import qpx

        out, _ = self._partition(dataset_dir, tmp_path)
        ds = qpx.open_dataset(str(out))
        result = ds.compute_integrity()
        ds.close()

        for key in result["file_checksums"]:
            assert not key.startswith("/"), f"absolute path leaked into the record: {key}"
            assert "\\" not in key, f"non-POSIX separator in the record: {key}"
            assert (out / key).exists(), f"key does not resolve under the dataset root: {key}"

    def test_verify_checks_nested_dataset_parquet(self, dataset_dir, tmp_path):
        """Only the root metadata file is exempt from its own checksum."""
        import qpx

        out, _ = self._partition(dataset_dir, tmp_path)
        nested = out / "metadata" / "part-0.dataset.parquet"
        nested.parent.mkdir()
        nested.write_bytes(next(out.glob("*.dataset.parquet")).read_bytes())

        ds = qpx.open_dataset(str(out))
        integrity = ds.compute_integrity()
        meta_dict = ds.dataset_meta.to_df().iloc[0].to_dict()
        meta_dict.update(integrity)
        ds.save_structure([meta_dict], "dataset", prefix="exp")
        ds.refresh()

        nested.write_bytes(b"corrupted")
        result = ds.verify_integrity()
        ds.close()

        key = nested.relative_to(out).as_posix()
        assert f"Checksum mismatch: {key}" in result["errors"]
