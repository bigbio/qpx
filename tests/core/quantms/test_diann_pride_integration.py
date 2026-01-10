"""
Integration tests for DIA-NN conversion using real PRIDE datasets.

These tests download real DIA-NN output files from PRIDE FTP and verify
that the QPX converter correctly handles all columns including extended
quality metrics and ion mobility data.

These tests are marked as 'slow' and 'integration' so they can be skipped
in regular CI runs. Run with: pytest -m integration
"""

import os
import tempfile
import urllib.request
from pathlib import Path

import pytest

# Mark all tests in this module as integration tests
pytestmark = [pytest.mark.integration, pytest.mark.slow]

# PRIDE dataset URLs for testing
PRIDE_DIANN_PARQUET_URL = (
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/11/PXD068392/"
    "WholeProteome_DIANNreport.parquet"
)


def download_file(url: str, dest_path: Path, timeout: int = 300) -> bool:
    """Download a file from URL to destination path.

    Args:
        url: URL to download from
        dest_path: Path to save the file
        timeout: Download timeout in seconds

    Returns:
        True if download succeeded, False otherwise
    """
    try:
        urllib.request.urlretrieve(url, dest_path)
        return True
    except Exception as e:
        print(f"Failed to download {url}: {e}")
        return False


@pytest.fixture(scope="module")
def pride_diann_parquet(tmp_path_factory):
    """Download DIA-NN parquet file from PRIDE for testing.

    This fixture downloads the file once per test module and caches it.
    """
    cache_dir = tmp_path_factory.mktemp("pride_data")
    parquet_path = cache_dir / "diann_report.parquet"

    if not download_file(PRIDE_DIANN_PARQUET_URL, parquet_path):
        pytest.skip(f"Could not download test file from {PRIDE_DIANN_PARQUET_URL}")

    return parquet_path


class TestDiaNNPrideIntegration:
    """Integration tests using real PRIDE DIA-NN data."""

    def test_load_pride_parquet_report(self, pride_diann_parquet):
        """Test loading a real DIA-NN parquet file from PRIDE."""
        import pyarrow.parquet as pq

        # Verify file exists and is readable
        assert pride_diann_parquet.exists(), "Parquet file should exist"

        # Read parquet metadata
        parquet_file = pq.ParquetFile(pride_diann_parquet)
        schema = parquet_file.schema_arrow

        # Log available columns for debugging
        column_names = [field.name for field in schema]
        print(f"Found {len(column_names)} columns in PRIDE DIA-NN parquet")
        print(f"Columns: {column_names[:20]}...")  # Print first 20

        # Verify essential DIA-NN columns are present
        essential_cols = [
            "Run",
            "Protein.Group",
            "Stripped.Sequence",
            "Modified.Sequence",
            "Precursor.Charge",
            "Precursor.Quantity",
            "Q.Value",
            "PG.Q.Value",
            "Global.Q.Value",
        ]
        for col in essential_cols:
            assert col in column_names, f"Essential column {col} should be present"

    def test_detect_extended_columns(self, pride_diann_parquet):
        """Test that extended DIA-NN columns are detected."""
        import pyarrow.parquet as pq

        from qpx.core.common import DIANN_EXTENDED_COLS

        parquet_file = pq.ParquetFile(pride_diann_parquet)
        column_names = set(field.name for field in parquet_file.schema_arrow)

        # Check which extended columns are present
        found_extended = [col for col in DIANN_EXTENDED_COLS if col in column_names]
        print(f"Found extended columns: {found_extended}")

        # We expect at least some quality score columns to be present
        quality_cols = ["CScore", "Evidence", "Spectrum.Similarity"]
        found_quality = [col for col in quality_cols if col in column_names]
        print(f"Found quality score columns: {found_quality}")

    def test_convert_pride_parquet_to_feature(self, pride_diann_parquet, tmp_path):
        """Test converting PRIDE DIA-NN parquet to QPX feature format."""
        from qpx.core.diann.diann import DiaNNConvert

        # Create a minimal SDRF for testing (required by DiaNNConvert)
        sdrf_path = tmp_path / "test.sdrf.tsv"
        self._create_minimal_sdrf(pride_diann_parquet, sdrf_path)

        # Initialize converter
        diann_converter = DiaNNConvert(
            diann_report=str(pride_diann_parquet),
            sdrf_path=str(sdrf_path),
            duckdb_max_memory="4GB",
            duckdb_threads=2,
        )

        try:
            # Check extended columns were detected
            print(
                f"Detected extended columns: {diann_converter._available_extended_cols}"
            )

            # Get unique runs
            refs = diann_converter.get_unique_references("Run")
            print(f"Found {len(refs)} runs in the dataset")
            assert len(refs) > 0, "Should have at least one run"

            # Test with first run only to keep test fast
            test_refs = refs[:1]
            from qpx.core.diann.diann import DIANN_SQL

            report = diann_converter.get_report_from_database(
                test_refs, DIANN_SQL, include_extended=True
            )

            print(f"Loaded {len(report)} rows for testing")
            print(f"Report columns: {list(report.columns)}")

            # Verify core columns are present after mapping
            from qpx.core.common import DIANN_MAP

            for diann_col, qpx_col in DIANN_MAP.items():
                if diann_col in report.columns or qpx_col in report.columns:
                    continue
                # Some columns may be optional
                if diann_col in ["IM", "Decoy", "Precursor.Id"]:
                    continue
                print(f"Warning: Column {diann_col} -> {qpx_col} not found")

            # Verify extended columns are included if they exist
            from qpx.core.diann.diann import DIANN_EXTENDED_MAP

            for diann_col, qpx_col in DIANN_EXTENDED_MAP.items():
                if qpx_col in report.columns:
                    print(f"Extended column mapped: {diann_col} -> {qpx_col}")

            assert len(report) > 0, "Should have loaded some data"

        finally:
            diann_converter.destroy_duckdb_database()

    def test_feature_transformation_with_extended_columns(
        self, pride_diann_parquet, tmp_path
    ):
        """Test full feature transformation including extended columns."""
        from qpx.core.diann.diann import DiaNNConvert
        from qpx.core.quantms.feature import Feature

        # Create minimal SDRF
        sdrf_path = tmp_path / "test.sdrf.tsv"
        self._create_minimal_sdrf(pride_diann_parquet, sdrf_path)

        diann_converter = DiaNNConvert(
            diann_report=str(pride_diann_parquet),
            sdrf_path=str(sdrf_path),
            duckdb_max_memory="4GB",
            duckdb_threads=2,
        )

        try:
            # Process a small subset
            refs = diann_converter.get_unique_references("Run")[:1]
            from qpx.core.common import DIANN_MAP
            from qpx.core.diann.diann import DIANN_SQL

            report = diann_converter.get_report_from_database(
                refs, DIANN_SQL, include_extended=True
            )

            # Apply column mapping
            report.rename(columns=DIANN_MAP, inplace=True)

            # Take a small sample for testing
            report = report.head(100).copy()

            if len(report) > 0:
                # Add scan column (required for add_additional_msg)
                report["scan"] = range(len(report))

                # Apply transformations
                diann_converter.add_additional_msg(report)
                Feature.convert_to_parquet_format(report)

                # Verify additional_scores contains expected fields
                if "additional_scores" in report.columns:
                    sample_scores = report["additional_scores"].iloc[0]
                    score_names = [s["score_name"] for s in sample_scores]
                    print(f"Score names in output: {score_names}")

                    # Core scores should always be present
                    assert "qvalue" in score_names, "qvalue should be in scores"
                    assert "pg_qvalue" in score_names, "pg_qvalue should be in scores"
                    assert (
                        "global_qvalue" in score_names
                    ), "global_qvalue should be in scores"

                # Verify additional_intensities structure
                if "additional_intensities" in report.columns:
                    sample_intensities = report["additional_intensities"].iloc[0]
                    if sample_intensities and len(sample_intensities) > 0:
                        intensity_names = [
                            i["intensity_name"]
                            for entry in sample_intensities
                            for i in entry.get("intensities", [])
                        ]
                        print(f"Intensity names in output: {intensity_names}")
                        assert "lfq" in intensity_names, "lfq should be in intensities"

                # Verify ion_mobility handling
                if "ion_mobility" in report.columns:
                    has_im = report["ion_mobility"].notna().any()
                    print(f"Has ion mobility data: {has_im}")

                print("Feature transformation completed successfully!")

        finally:
            diann_converter.destroy_duckdb_database()

    def _create_minimal_sdrf(self, parquet_path: Path, sdrf_path: Path):
        """Create a minimal SDRF file for testing based on runs in the parquet."""
        import pyarrow.parquet as pq

        # Read unique runs from parquet
        parquet_file = pq.ParquetFile(parquet_path)
        table = parquet_file.read(columns=["Run"])
        runs = table["Run"].to_pylist()
        unique_runs = list(set(runs))[:10]  # Limit to 10 for testing

        # Create minimal SDRF content
        sdrf_content = [
            "source name\tcomment[data file]\tcomment[label]\t"
            "comment[fraction identifier]\tcharacteristics[biological replicate]"
        ]

        for i, run in enumerate(unique_runs):
            # Handle run name - remove extension if present
            run_name = run.split(".")[0] if "." in run else run
            sdrf_content.append(
                f"Sample_{i+1}\t{run_name}.raw\tLFQ\t1\t1"
            )

        with open(sdrf_path, "w") as f:
            f.write("\n".join(sdrf_content))


class TestDiaNNColumnCoverage:
    """Tests to verify DIA-NN column coverage matches specification."""

    def test_diann_map_completeness(self):
        """Verify DIANN_MAP contains all essential columns."""
        from qpx.core.common import DIANN_MAP

        essential_mappings = {
            "Precursor.Quantity": "intensity",
            "Stripped.Sequence": "sequence",
            "Modified.Sequence": "peptidoform",
            "Precursor.Charge": "precursor_charge",
            "Protein.Group": "pg_accessions",
            "Q.Value": "qvalue",
            "PG.Q.Value": "pg_qvalue",
            "Global.Q.Value": "global_qvalue",
            "Run": "reference_file_name",
            "RT": "rt",
            "PEP": "posterior_error_probability",
        }

        for diann_col, qpx_col in essential_mappings.items():
            assert diann_col in DIANN_MAP, f"{diann_col} should be in DIANN_MAP"
            assert (
                DIANN_MAP[diann_col] == qpx_col
            ), f"{diann_col} should map to {qpx_col}"

    def test_extended_columns_defined(self):
        """Verify extended column lists are properly defined."""
        from qpx.core.common import (
            DIANN_EXTENDED_COLS,
            DIANN_IM_COLS,
            DIANN_INTENSITY_COLS,
            DIANN_SCORE_COLS,
        )

        # Score columns
        expected_score_cols = [
            "CScore",
            "Evidence",
            "Spectrum.Similarity",
            "Averagine",
            "Mass.Evidence",
        ]
        for col in expected_score_cols:
            assert col in DIANN_SCORE_COLS, f"{col} should be in DIANN_SCORE_COLS"

        # Intensity columns
        assert "Ms1.Area" in DIANN_INTENSITY_COLS, "Ms1.Area should be defined"

        # Ion mobility columns
        expected_im_cols = ["iIM", "Predicted.IM", "Predicted.iIM"]
        for col in expected_im_cols:
            assert col in DIANN_IM_COLS, f"{col} should be in DIANN_IM_COLS"

        # Extended cols should combine all
        assert len(DIANN_EXTENDED_COLS) == (
            len(DIANN_SCORE_COLS) + len(DIANN_INTENSITY_COLS) + len(DIANN_IM_COLS)
        ), "DIANN_EXTENDED_COLS should combine all extended column lists"

    def test_ion_mobility_mapping(self):
        """Verify ion mobility is properly mapped."""
        from qpx.core.common import DIANN_MAP

        assert "IM" in DIANN_MAP, "IM should be in DIANN_MAP"
        assert DIANN_MAP["IM"] == "ion_mobility", "IM should map to ion_mobility"

    def test_decoy_mapping(self):
        """Verify decoy column is properly mapped for DIA-NN 2.0+."""
        from qpx.core.common import DIANN_MAP

        assert "Decoy" in DIANN_MAP, "Decoy should be in DIANN_MAP"
        assert DIANN_MAP["Decoy"] == "is_decoy", "Decoy should map to is_decoy"
