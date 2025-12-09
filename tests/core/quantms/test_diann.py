from pathlib import Path
import math


import pytest

from qpx.core.diann.diann import DiaNNConvert
from qpx.core.quantms.feature import Feature

TEST_DATA_ROOT = Path(__file__).parents[2] / "examples"

TEST_DATA = (
    TEST_DATA_ROOT / "diann/small/diann_report.tsv",
    TEST_DATA_ROOT / "diann/small/PXD019909-DIA.sdrf.tsv",
    TEST_DATA_ROOT / "diann/small/mzml",
    TEST_DATA_ROOT / "diann/small/diann_report.pg_matrix.tsv",
)


def test_transform_feature():
    report_file = TEST_DATA[0]
    sdrf_file = TEST_DATA[1]
    mzml = TEST_DATA[2]
    diann_converter = DiaNNConvert(diann_report=report_file, sdrf_path=sdrf_file)
    try:
        for report in diann_converter.main_report_df(0.05, mzml, 2):
            diann_converter.add_additional_msg(report)
            Feature.convert_to_parquet_format(report)
            Feature.transform_feature(report)
    finally:
        # Clean up DuckDB database
        diann_converter.destroy_duckdb_database()


def test_transform_features():
    report_file = TEST_DATA[0]
    sdrf_file = TEST_DATA[1]
    mzml = TEST_DATA[2]
    diann_converter = DiaNNConvert(diann_report=report_file, sdrf_path=sdrf_file)
    try:
        for report in diann_converter.main_report_df(0.05, mzml, 2):
            diann_converter.add_additional_msg(report)
            Feature.convert_to_parquet_format(report)
            for _, df in Feature.slice(
                report, ["reference_file_name", "precursor_charge"]
            ):
                Feature.transform_feature(df)
    finally:
        # Clean up DuckDB database
        diann_converter.destroy_duckdb_database()


def test_transform_protein_groups():
    """Test transforming DIA-NN protein group data."""
    report_file = TEST_DATA[0]
    pg_matrix_file = TEST_DATA[3]
    sdrf_file = TEST_DATA[1]
    diann_converter = DiaNNConvert(
        diann_report=report_file, pg_matrix_path=pg_matrix_file, sdrf_path=sdrf_file
    )

    try:
        # Get some test data for protein groups using the proper SQL format
        refs = diann_converter.get_unique_references("Run")[:1]  # Just test with 1 file
        from qpx.core.common import DIANN_PG_MAP
        from qpx.core.diann.diann import DIANN_PG_SQL

        report = diann_converter.get_report_from_database(refs, DIANN_PG_SQL)

        # Apply the DIA-NN protein group mapping
        report.rename(columns=DIANN_PG_MAP, inplace=True)
        report.dropna(subset=["pg_accessions"], inplace=True)

        if len(report) > 0:  # Only test if we have data
            # Test one file worth of data
            df = report.head(5).copy()  # Test with 5 protein groups

            if len(df) > 0:
                # Transform the protein group data
                df = diann_converter.get_report_pg_matrix(
                    df, diann_converter.pg_matrix, refs[0]
                )
                pg_df = diann_converter.generate_pg_matrix(df)

                # Verify the structure
                if "intensities" not in pg_df.columns:
                    pytest.fail("intensities column should be present")
                if "additional_intensities" not in pg_df.columns:
                    pytest.fail("additional_intensities column should be present")

                # Check intensities structure
                intensity_sample = pg_df["intensities"].iloc[0]
                if not isinstance(intensity_sample, list):
                    pytest.fail("intensities should be a list")
                if len(intensity_sample) <= 0:
                    pytest.fail("intensities should not be empty")

                intensity_entry = intensity_sample[0]
                if "sample_accession" not in intensity_entry:
                    pytest.fail("intensity entry should have sample_accession")
                if "channel" not in intensity_entry:
                    pytest.fail("intensity entry should have channel")
                if "intensity" not in intensity_entry:
                    pytest.fail("intensity entry should have intensity value")
                if intensity_entry["channel"] != "LFQ":
                    pytest.fail("channel should be LFQ for DIA-NN")

                # Check additional_intensities structure
                additional_intensity_sample = pg_df["additional_intensities"].iloc[0]
                if not isinstance(additional_intensity_sample, list):
                    pytest.fail("additional_intensities should be a list")
                if len(additional_intensity_sample) <= 0:
                    pytest.fail("additional_intensities should not be empty")

                additional_entry = additional_intensity_sample[0]
                if "sample_accession" not in additional_entry:
                    pytest.fail(
                        "additional_intensity entry should have sample_accession"
                    )
                if "channel" not in additional_entry:
                    pytest.fail("additional_intensity entry should have channel")
                if "intensities" not in additional_entry:
                    pytest.fail(
                        "additional_intensity entry should have intensities array"
                    )
                if additional_entry["channel"] != "LFQ":
                    pytest.fail("channel should be LFQ for DIA-NN")

                # Check the additional intensity array structure
                additional_types = additional_entry["intensities"]
                if not isinstance(additional_types, list):
                    pytest.fail("intensities should be a list")

                # DIA-NN version 2.0 and later no longer have "PG.Normalised", so "normalize_intensity" is missing.

                # Check normalize_intensity and lfq are present
                intensity_names = [item["intensity_name"] for item in additional_types]
                # DIA-NN version 2.0 and later no longer have "PG.Normalised", so "normalize_intensity" is missing.
                if "lfq" not in intensity_names:
                    pytest.fail("lfq should be present")

                print("Protein group intensity structure test passed!")
    finally:
        # Clean up DuckDB database
        diann_converter.destroy_duckdb_database()


class TestDiaNNStandardizedIntensities:
    """Property tests for standardized intensity fields in DIA-NN output."""

    def test_property_standardized_fields_present_in_additional_intensities(self):
        """
        Property 4: For any processed PG data, additional_intensities should contain
        entries named "total_all_peptides_intensity" and "top3_intensity".

        This test verifies that the DIA-NN workflow produces output with the
        standardized field names as specified in the requirements.
        """
        report_file = TEST_DATA[0]
        pg_matrix_file = TEST_DATA[3]
        sdrf_file = TEST_DATA[1]
        diann_converter = DiaNNConvert(
            diann_report=report_file, pg_matrix_path=pg_matrix_file, sdrf_path=sdrf_file
        )

        try:
            # Get test data for protein groups
            refs = diann_converter.get_unique_references("Run")[:1]
            from qpx.core.common import DIANN_PG_MAP
            from qpx.core.diann.diann import DIANN_PG_SQL

            report = diann_converter.get_report_from_database(refs, DIANN_PG_SQL)

            report.rename(columns=DIANN_PG_MAP, inplace=True)
            if "Precursor.Quantity" in report.columns:
                report.rename(columns={"Precursor.Quantity": "intensity"}, inplace=True)
            report.dropna(subset=["pg_accessions"], inplace=True)

            if len(report) > 0:
                df = diann_converter.get_report_pg_matrix(
                    report,
                    diann_converter.pg_matrix,
                    refs[0],
                    calculate_standardized_intensities=True,
                )
                pg_df = diann_converter.generate_pg_matrix(
                    df, calculate_standardized_intensities=True
                )

                if "additional_intensities" not in pg_df.columns:
                    pytest.fail("additional_intensities column should be present")

                for _, row in pg_df.iterrows():
                    additional_intensities = row["additional_intensities"]
                    if not isinstance(additional_intensities, list):
                        pytest.fail("additional_intensities should be a list")
                    if len(additional_intensities) <= 0:
                        pytest.fail("additional_intensities should not be empty")

                    intensity_names = []
                    for entry in additional_intensities:
                        if "intensities" not in entry:
                            pytest.fail("entry should have intensities array")
                        for intensity_item in entry["intensities"]:
                            intensity_names.append(intensity_item["intensity_name"])

                    if "total_all_peptides_intensity" not in intensity_names:
                        pytest.fail(
                            f"total_all_peptides_intensity should be present, found: {intensity_names}"
                        )
                    if "top3_intensity" not in intensity_names:
                        pytest.fail(
                            f"top3_intensity should be present, found: {intensity_names}"
                        )

                print("Property 4: Standardized field name consistency test passed!")
        finally:
            diann_converter.destroy_duckdb_database()

    def test_property_standardized_intensity_values_are_numeric(self):
        """
        Verify that standardized intensity values are valid numeric values.
        """
        report_file = TEST_DATA[0]
        pg_matrix_file = TEST_DATA[3]
        sdrf_file = TEST_DATA[1]
        diann_converter = DiaNNConvert(
            diann_report=report_file, pg_matrix_path=pg_matrix_file, sdrf_path=sdrf_file
        )

        try:
            refs = diann_converter.get_unique_references("Run")[:1]
            from qpx.core.common import DIANN_PG_MAP
            from qpx.core.diann.diann import DIANN_PG_SQL

            report = diann_converter.get_report_from_database(refs, DIANN_PG_SQL)
            report.rename(columns=DIANN_PG_MAP, inplace=True)
            if "Precursor.Quantity" in report.columns:
                report.rename(columns={"Precursor.Quantity": "intensity"}, inplace=True)
            report.dropna(subset=["pg_accessions"], inplace=True)

            if len(report) > 0:
                df = diann_converter.get_report_pg_matrix(
                    report,
                    diann_converter.pg_matrix,
                    refs[0],
                    calculate_standardized_intensities=True,
                )
                pg_df = diann_converter.generate_pg_matrix(
                    df, calculate_standardized_intensities=True
                )

                for _, row in pg_df.iterrows():
                    additional_intensities = row["additional_intensities"]
                    for entry in additional_intensities:
                        for intensity_item in entry["intensities"]:
                            if intensity_item["intensity_name"] in [
                                "total_all_peptides_intensity",
                                "top3_intensity",
                            ]:
                                value = intensity_item["intensity_value"]
                                if not isinstance(value, (int, float)):
                                    pytest.fail(
                                        f"Intensity value should be numeric, got {type(value)}"
                                    )
                                if not (value >= 0 or math.isnan(value)):
                                    pytest.fail(
                                        f"Intensity value should be non-negative or NaN, got {value}"
                                    )

                print("Standardized intensity values are valid numeric values!")
        finally:
            diann_converter.destroy_duckdb_database()
