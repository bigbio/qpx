"""Test metadata schema generator."""

import shutil
from pathlib import Path

import pytest

from qpx.core.metadata import (
    MetadataSchemaGenerator,
    WorkflowMetadataGenerator,
)


@pytest.fixture(scope="module")
def test_output_dir():
    """Create test output directory and clean up after tests."""
    output_dir = Path("configuration_test")
    output_dir.mkdir(exist_ok=True)
    yield output_dir
    # Cleanup after all tests in this module
    if output_dir.exists():
        shutil.rmtree(output_dir)


@pytest.fixture(autouse=True)
def cleanup_test_files():
    """Clean up any test CSV files created during tests."""
    yield
    # Cleanup after each test
    test_files = [
        "test_metadata.csv",
        "test_metadata_all.csv",
        "test_metadata_fixed.csv",
        "PXD000561_metadata.csv",
        "test_evidence_only.csv",
    ]
    for test_file in test_files:
        file_path = Path(test_file)
        if file_path.exists():
            file_path.unlink()


def test_metadata_schema_generator_basic():
    """Test basic metadata schema generator functionality."""
    generator = MetadataSchemaGenerator()

    generator.add_column_metadata(
        model_view="psm",
        column_name="charge",
        ontology_accession="MS:1000041",
        full_name="charge state",
        settings="computed as provided by MaxQuant",
    )

    records = generator.get_records()
    assert len(records) == 1
    assert records[0]["model_view"] == "psm"
    assert records[0]["column_name"] == "charge"
    assert records[0]["ontology_accession"] == "MS:1000041"


def test_metadata_file_generation(test_output_dir):
    """Test metadata file generation."""
    generator = MetadataSchemaGenerator()

    generator.add_column_metadata(
        model_view="psm",
        column_name="sequence",
        full_name="peptide sequence",
        settings="extracted from Modified sequence",
    )

    generator.add_column_metadata(
        model_view="feature",
        column_name="intensity",
        full_name="peptide intensity",
        settings="raw intensity value",
    )

    output_file = test_output_dir / "metadata_basic.csv"
    generator.generate_file(str(output_file))

    assert output_file.exists()

    with open(output_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    assert "model_view" in lines[0]
    assert "column_name" in lines[0]
    assert len(lines) == 3


def test_workflow_metadata_generator_maxquant(test_output_dir):
    """Test MaxQuant workflow metadata generation."""
    generator = WorkflowMetadataGenerator(workflow="maxquant")

    generator.generate_all_metadata()

    records = generator.get_records()

    assert len(records) > 0

    model_views = {r["model_view"] for r in records}
    assert "psm" in model_views
    assert "feature" in model_views
    assert "pg" in model_views

    output_file = test_output_dir / "maxquant_metadata.csv"
    generator.generate_file(str(output_file))
    assert output_file.exists()


def test_workflow_metadata_generator_diann(test_output_dir):
    """Test DIANN workflow metadata generation."""
    generator = WorkflowMetadataGenerator(workflow="diann")

    generator.generate_pg_metadata()

    records = generator.get_records()
    assert len(records) > 0
    pg_records = [r for r in records if r["model_view"] == "pg"]
    assert len(pg_records) > 0

    output_file = test_output_dir / "diann_metadata.csv"
    generator.generate_file(str(output_file))
    assert output_file.exists()


def test_workflow_metadata_generator_quantms_lfq(test_output_dir):
    """Test quantms-LFQ workflow metadata generation."""
    generator = WorkflowMetadataGenerator(workflow="quantms-lfq")

    generator.generate_pg_metadata()

    records = generator.get_records()
    assert len(records) > 0

    pg_records = [r for r in records if r["model_view"] == "pg"]
    assert len(pg_records) > 0

    output_file = test_output_dir / "quantms_lfq_metadata.csv"
    generator.generate_file(str(output_file))
    assert output_file.exists()


def test_workflow_metadata_generator_quantms_tmt(test_output_dir):
    """Test quantms-TMT workflow metadata generation."""
    generator = WorkflowMetadataGenerator(workflow="quantms-tmt")

    generator.generate_pg_metadata()

    records = generator.get_records()
    assert len(records) > 0

    output_file = test_output_dir / "quantms_tmt_metadata.csv"
    generator.generate_file(str(output_file))
    assert output_file.exists()


def test_workflow_metadata_generator_quantms_psm(test_output_dir):
    """Test quantms-PSM workflow metadata generation."""
    generator = WorkflowMetadataGenerator(workflow="quantms-psm")

    generator.generate_all_metadata()

    records = generator.get_records()
    assert len(records) > 0

    psm_records = [r for r in records if r["model_view"] == "psm"]
    assert len(psm_records) > 0

    output_file = test_output_dir / "quantms_psm_metadata.csv"
    generator.generate_file(str(output_file))
    assert output_file.exists()

    psm_with_settings = [r for r in psm_records if r.get("settings", "").strip() != ""]
    assert len(psm_with_settings) > 0


def test_invalid_workflow():
    """Test that invalid workflow raises error."""
    with pytest.raises(ValueError):
        WorkflowMetadataGenerator(workflow="invalid_workflow")


def test_metadata_clear():
    """Test clearing metadata records."""
    generator = MetadataSchemaGenerator()

    generator.add_column_metadata(
        model_view="psm",
        column_name="test_column",
    )

    assert len(generator.get_records()) == 1

    generator.clear()
    assert len(generator.get_records()) == 0


def test_workflow_metadata_with_available_columns():
    """Test metadata generation with available columns filtering."""
    available_columns = {
        "Sequence",
        "Modified sequence",
        "Charge",
        "m/z",
        "Retention time",
        "PEP",
        "Reverse",
        "1/K0",
        "Raw file",
    }

    generator = WorkflowMetadataGenerator(
        workflow="maxquant", available_columns=available_columns
    )

    generator.generate_feature_metadata()

    records = generator.get_records()
    assert len(records) > 0

    ion_mobility_records = [r for r in records if r["column_name"] == "ion_mobility"]
    assert len(ion_mobility_records) == 1
    assert "1/K0" in ion_mobility_records[0]["settings"]
    assert "reduced ion mobility inverse" in ion_mobility_records[0]["settings"]

    column_names = {r["column_name"] for r in records}
    assert "sequence" in column_names
    assert "peptidoform" in column_names
    assert "ion_mobility" in column_names


def test_workflow_metadata_without_ion_mobility():
    """Test metadata generation when ion mobility column is not available."""
    available_columns = {
        "Sequence",
        "Modified sequence",
        "Charge",
        "m/z",
        "Retention time",
    }

    generator = WorkflowMetadataGenerator(
        workflow="maxquant", available_columns=available_columns
    )

    generator.generate_feature_metadata()

    records = generator.get_records()

    ion_mobility_records = [r for r in records if r["column_name"] == "ion_mobility"]
    assert len(ion_mobility_records) == 0


def test_reset_and_reinitialize_mappings():
    """Test that mappings can be reset and reinitialized correctly."""
    generator = WorkflowMetadataGenerator(workflow="maxquant")

    generator.generate_psm_metadata()
    initial_records = len(generator.get_records())
    assert initial_records > 0

    generator._reset_mappings()
    generator._init_workflow_mappings()

    generator.clear()
    generator.generate_psm_metadata()
    new_records = len(generator.get_records())

    assert new_records == initial_records


def test_ion_mobility_description():
    """Test that ion mobility has correct description with units."""
    generator = WorkflowMetadataGenerator(workflow="maxquant")

    generator.generate_psm_metadata()

    records = generator.get_records()
    ion_mobility_records = [r for r in records if r["column_name"] == "ion_mobility"]

    if len(ion_mobility_records) > 0:
        settings = ion_mobility_records[0]["settings"]
        assert "1/K0" in settings or "reduced ion mobility inverse" in settings
