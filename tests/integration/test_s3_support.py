"""S3 support tests (unit tests only, no real S3 access)."""


class TestS3EngineConfig:
    def test_httpfs_setup_with_credentials(self):
        from qpx.core.engine import DuckDBEngine

        engine = DuckDBEngine(
            s3_config={
                "region": "us-east-1",
                "access_key_id": "test_key",
                "secret_access_key": "test_secret",
            }
        )
        # Verify httpfs is loaded
        result = engine.execute("SELECT * FROM duckdb_extensions() WHERE extension_name = 'httpfs' AND loaded = true").fetchone()
        assert result is not None
        engine.close()

    def test_no_httpfs_by_default(self):
        from qpx.core.engine import DuckDBEngine

        engine = DuckDBEngine()
        # httpfs should NOT be loaded by default
        result = engine.execute("SELECT * FROM duckdb_extensions() WHERE extension_name = 'httpfs' AND loaded = true").fetchone()
        assert result is None
        engine.close()

    def test_anonymous_s3_config(self):
        from qpx.core.engine import DuckDBEngine

        engine = DuckDBEngine(s3_config={"region": "us-east-1", "anonymous": True})
        result = engine.execute("SELECT * FROM duckdb_extensions() WHERE extension_name = 'httpfs' AND loaded = true").fetchone()
        assert result is not None
        # Cursors, which stream query rows, see the S3 settings too.
        setting = "SELECT current_setting('s3_region'), current_setting('s3_url_style')"
        assert engine.connection.cursor().execute(setting).fetchall() == [("us-east-1", "path")]
        engine.close()

    def test_s3_with_custom_endpoint(self):
        from qpx.core.engine import DuckDBEngine

        engine = DuckDBEngine(
            s3_config={
                "region": "us-east-1",
                "endpoint": "minio.example.com:9000",
                "access_key_id": "minioadmin",
                "secret_access_key": "minioadmin",
            }
        )
        result = engine.execute("SELECT * FROM duckdb_extensions() WHERE extension_name = 'httpfs' AND loaded = true").fetchone()
        assert result is not None
        engine.close()


class TestS3DatasetDetection:
    @staticmethod
    def _s3_dataset(monkeypatch):
        """Open an S3 path without network access: discovery would list the bucket."""
        from qpx.dataset import Dataset

        monkeypatch.setattr(Dataset, "_discover_s3", lambda self, requested: None)
        return Dataset("s3://bucket/dataset/")

    def test_is_s3_flag_set(self, monkeypatch):
        ds = self._s3_dataset(monkeypatch)
        assert ds._is_s3 is True
        ds.close()

    def test_local_path_not_s3(self, tmp_path):
        import qpx

        ds = qpx.open_dataset(str(tmp_path))
        assert ds._is_s3 is False
        ds.close()

    def test_s3_path_is_string(self, monkeypatch):
        """S3 paths should stay as strings, not become Path objects."""
        ds = self._s3_dataset(monkeypatch)
        assert isinstance(ds.path, str)
        ds.close()
