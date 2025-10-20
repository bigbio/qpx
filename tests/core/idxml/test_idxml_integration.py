import tempfile
import shutil
from pathlib import Path
import pytest
from click.testing import CliRunner
import pyarrow.parquet as pq

from quantmsio.commands.convert.quantms import convert_idxml_psm_cmd
from quantmsio.commands.convert.idxml import convert_idxml_batch

# Test data path
TEST_DATA_ROOT = Path(__file__).parents[2] / "examples"
IDXML_TEST_PATH = (
    TEST_DATA_ROOT
    / "idxml/SF_200217_pPeptideLibrary_pool1_HCDnlETcaD_OT_rep2_consensus_fdr_pep_luciphor.idXML"
)


def test_idxml_psm_command():
    """Test the idXML PSM conversion command."""
    # Skip if test data doesn't exist
    if not IDXML_TEST_PATH.exists():
        pytest.skip(f"Test data not found: {IDXML_TEST_PATH}")

    runner = CliRunner()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Run the command
        result = runner.invoke(
            convert_idxml_psm_cmd,
            [
                "--idxml-path",
                str(IDXML_TEST_PATH),
                "--output-folder",
                temp_dir,
                "--output-prefix",
                "test-idxml",
                "--verbose",
            ],
        )

        # Check that command succeeded
        assert result.exit_code == 0, f"Command failed with output: {result.output}"

        # Check that output file was created
        output_files = list(Path(temp_dir).glob("*.parquet"))
        assert (
            len(output_files) == 1
        ), f"Expected 1 output file, found {len(output_files)}"

        output_file = output_files[0]
        assert output_file.suffix == ".parquet"
        assert "test-idxml" in output_file.name
        assert output_file.stat().st_size > 0

        print(f"Command succeeded, created file: {output_file}")


def test_idxml_psm_command_missing_file():
    """Test the idXML PSM conversion command with missing input file."""
    runner = CliRunner()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Run the command with non-existent file
        result = runner.invoke(
            convert_idxml_psm_cmd,
            [
                "--idxml-path",
                "non_existent_file.idXML",
                "--output-folder",
                temp_dir,
            ],
        )

        # Command should fail due to missing file
        assert result.exit_code != 0
        assert "does not exist" in result.output or "Error" in result.output


def test_idxml_psm_command_help():
    """Test the idXML PSM conversion command help."""
    runner = CliRunner()

    result = runner.invoke(convert_idxml_psm_cmd, ["--help"])

    assert result.exit_code == 0
    assert "Convert PSM data from idXML to quantms.io format" in result.output
    assert "--idxml-path" in result.output
    assert "--output-folder" in result.output


def test_idxml_batch_command_folder_mode():
    """Test the idXML batch conversion command with folder input."""
    # Skip if test data doesn't exist
    if not IDXML_TEST_PATH.exists():
        pytest.skip(f"Test data not found: {IDXML_TEST_PATH}")

    runner = CliRunner()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a test folder with multiple idXML files (copies of the test file)
        test_idxml_folder = Path(temp_dir) / "idxml_files"
        test_idxml_folder.mkdir()

        # Create 3 copies of the test file with different names
        for i in range(1, 4):
            dest_file = test_idxml_folder / f"test_file_{i}.idXML"
            shutil.copy(IDXML_TEST_PATH, dest_file)

        output_folder = Path(temp_dir) / "output"
        output_folder.mkdir()

        # Run the batch conversion command
        result = runner.invoke(
            convert_idxml_batch,
            [
                "--idxml-folder",
                str(test_idxml_folder),
                "--output-folder",
                str(output_folder),
                "--output-prefix-file",
                "batch-test",
                "--verbose",
            ],
        )

        # Check that command succeeded
        assert result.exit_code == 0, f"Command failed with output: {result.output}"

        # Check that merged output file was created
        output_files = list(output_folder.glob("*.parquet"))
        assert (
            len(output_files) == 1
        ), f"Expected 1 merged output file, found {len(output_files)}"

        output_file = output_files[0]
        assert output_file.suffix == ".parquet"
        assert "batch-test" in output_file.name
        assert output_file.stat().st_size > 0

        # Verify the parquet file contains data from all 3 files
        table = pq.read_table(output_file)
        assert len(table) > 0, "Merged parquet file should contain PSMs"

        print(f"Batch command succeeded, created merged file: {output_file}")
        print(f"Merged file contains {len(table)} PSMs")


def test_idxml_batch_command_file_list_mode():
    """Test the idXML batch conversion command with file list input."""
    # Skip if test data doesn't exist
    if not IDXML_TEST_PATH.exists():
        pytest.skip(f"Test data not found: {IDXML_TEST_PATH}")

    runner = CliRunner()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Create multiple test files
        test_files = []
        for i in range(1, 3):
            dest_file = Path(temp_dir) / f"test_file_{i}.idXML"
            shutil.copy(IDXML_TEST_PATH, dest_file)
            test_files.append(str(dest_file))

        # Create comma-separated file list
        file_list = ",".join(test_files)

        output_folder = Path(temp_dir) / "output"
        output_folder.mkdir()

        # Run the batch conversion command with file list
        result = runner.invoke(
            convert_idxml_batch,
            [
                "--idxml-files",
                file_list,
                "--output-folder",
                str(output_folder),
                "--output-prefix-file",
                "filelist-test",
                "--verbose",
            ],
        )

        # Check that command succeeded
        assert result.exit_code == 0, f"Command failed with output: {result.output}"

        # Check that merged output file was created
        output_files = list(output_folder.glob("*.parquet"))
        assert (
            len(output_files) == 1
        ), f"Expected 1 merged output file, found {len(output_files)}"

        output_file = output_files[0]
        assert "filelist-test" in output_file.name

        # Verify the parquet file contains data
        table = pq.read_table(output_file)
        assert len(table) > 0, "Merged parquet file should contain PSMs"

        print(f"File list mode succeeded, created merged file: {output_file}")


def test_idxml_batch_command_missing_parameters():
    """Test the idXML batch conversion command with missing parameters."""
    runner = CliRunner()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Test without idxml-folder or idxml-files
        result = runner.invoke(
            convert_idxml_batch,
            [
                "--output-folder",
                temp_dir,
            ],
        )

        # Command should fail due to missing input parameter
        assert result.exit_code != 0
        assert "provide either --idxml-folder or --idxml-files" in result.output


def test_idxml_batch_command_conflicting_parameters():
    """Test the idXML batch conversion command with conflicting parameters."""
    runner = CliRunner()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Test with both idxml-folder and idxml-files (should fail)
        result = runner.invoke(
            convert_idxml_batch,
            [
                "--idxml-folder",
                temp_dir,
                "--idxml-files",
                "file1.idXML,file2.idXML",
                "--output-folder",
                temp_dir,
            ],
        )

        # Command should fail due to conflicting parameters
        assert result.exit_code != 0
        assert "only one of" in result.output


def test_idxml_batch_command_empty_folder():
    """Test the idXML batch conversion command with empty folder."""
    runner = CliRunner()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Create empty folder
        empty_folder = Path(temp_dir) / "empty"
        empty_folder.mkdir()

        output_folder = Path(temp_dir) / "output"
        output_folder.mkdir()

        # Run command with empty folder
        result = runner.invoke(
            convert_idxml_batch,
            [
                "--idxml-folder",
                str(empty_folder),
                "--output-folder",
                str(output_folder),
            ],
        )

        # Command should fail due to no files found
        assert result.exit_code != 0
        assert "No .idXML files found" in result.output


def test_idxml_batch_command_help():
    """Test the idXML batch conversion command help."""
    runner = CliRunner()

    result = runner.invoke(convert_idxml_batch, ["--help"])

    assert result.exit_code == 0
    assert "Convert multiple IdXML files" in result.output
    assert "--idxml-folder" in result.output
    assert "--idxml-files" in result.output
    assert "--output-folder" in result.output
    assert "--mzml-folder" in result.output
    assert "--mzml-files" in result.output


if __name__ == "__main__":
    # Run tests manually if needed
    print("Running single file tests...")
    test_idxml_psm_command()
    test_idxml_psm_command_missing_file()
    test_idxml_psm_command_help()

    print("\nRunning batch tests...")
    test_idxml_batch_command_folder_mode()
    test_idxml_batch_command_file_list_mode()
    test_idxml_batch_command_missing_parameters()
    test_idxml_batch_command_conflicting_parameters()
    test_idxml_batch_command_empty_folder()
    test_idxml_batch_command_help()

    print("\nAll integration tests passed!")
