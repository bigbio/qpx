"""Shared dataset-metadata options of the convert commands."""

import json
from http.client import IncompleteRead
from pathlib import Path

import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.convert import convert
from qpx.cli.convert.postprocess import _apply_software, _finish_conversion
from qpx.dataset import Dataset
from tests.converters.test_openms_consensus import _invoke_consensus_convert

_METADATA_OPTIONS = (
    "--software-name",
    "--software-version",
    "--software-config",
    "--enrich-pride",
    "--compute-integrity",
)
_METADATA_COMMANDS = (
    "diann",
    "maxquant",
    "fragpipe",
    "mzidentml",
    "openms-consensus",
    "quantms-msstats",
    "spectronaut",
    "cdap",
)
_DIANN_SMALL = Path(__file__).resolve().parents[1] / "examples" / "diann" / "small"


def _run(base, *extra, prefix="t"):
    """Run ``convert openms-consensus`` on a one-file consensusXML and return the result and output folder."""
    folder = base / "out"
    return _invoke_consensus_convert(base, folder, *extra, prefix=prefix), folder


def _run_diann(folder, *extra):
    """Run ``convert diann`` on the small DIA-NN example and return the CliRunner result."""
    return CliRunner().invoke(
        convert,
        [
            "diann",
            "--report-path",
            str(_DIANN_SMALL / "diann_report.tsv"),
            "--sdrf-file",
            str(_DIANN_SMALL / "PXD019909-DIA.sdrf.tsv"),
            "--output-folder",
            str(folder),
            *extra,
        ],
    )


def test_software_options_record_workflow_in_dataset_and_provenance(tmp_path):
    """--software-name and related options set dataset software fields and add a workflow provenance step."""
    params = tmp_path / "params.json"
    params.write_text(json.dumps({"b": 1, "a": [1, 2]}))
    baseline, baseline_folder = _run(tmp_path / "baseline")
    assert baseline.exit_code == 0, baseline.output
    before = pq.read_table(baseline_folder / "t.provenance.parquet").to_pylist()

    result, folder = _run(
        tmp_path,
        "--software-name",
        "nf-core/mhcquant",
        "--software-version",
        "3.1.0",
        "--software-config",
        str(params),
    )
    assert result.exit_code == 0, result.output

    dataset = pq.read_table(folder / "t.dataset.parquet").to_pylist()[0]
    assert dataset["software_name"] == "nf-core/mhcquant"
    assert dataset["software_version"] == "3.1.0"

    steps = pq.read_table(folder / "t.provenance.parquet").to_pylist()
    assert steps[0]["step_category"] == "workflow"
    assert steps[0]["step_name"] == "pipeline_execution"
    assert steps[0]["tool_name"] == "nf-core/mhcquant"
    assert steps[0]["tool_version"] == "3.1.0"
    assert steps[0]["config"] == '{"a":[1,2],"b":1}'
    assert [step["step_order"] for step in steps] == list(range(1, len(steps) + 1))
    assert steps[1]["step_name"] == before[0]["step_name"]
    assert steps[1]["tool_name"] == before[0]["tool_name"]


def test_software_version_without_name_is_usage_error(tmp_path):
    """--software-version without --software-name is a usage error that writes no output."""
    result, folder = _run(tmp_path, "--software-version", "1.0")

    assert result.exit_code == 2
    assert "require --software-name" in result.output
    assert not list(folder.glob("*.parquet"))


def test_software_config_must_be_valid_json(tmp_path):
    """--software-config must point to a readable JSON file."""
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")

    result, _ = _run(tmp_path, "--software-name", "x", "--software-config", str(bad))

    assert result.exit_code == 2
    assert "not a readable JSON file" in result.output


def test_compute_integrity_records_files_and_verifies(tmp_path):
    """--compute-integrity records file checksums, row counts and sizes that verify clean."""
    result, folder = _run(tmp_path, "--compute-integrity")

    assert result.exit_code == 0, result.output
    record = pq.read_table(folder / "t.dataset.parquet").to_pylist()[0]
    assert "t.feature.parquet" in dict(record["file_checksums"])
    assert "t.feature.parquet" in dict(record["file_row_counts"])
    assert "t.feature.parquet" in dict(record["file_sizes_bytes"])
    assert record["total_structures"] == len(dict(record["file_checksums"]))
    assert record["packaged_at"]
    with Dataset(str(folder)) as ds:
        assert ds.verify_integrity()["errors"] == []


def test_enrich_pride_sets_project_title(tmp_path, monkeypatch):
    """--enrich-pride writes the PRIDE project title into the dataset."""
    monkeypatch.setattr(
        "qpx.core.pride.fetch_pride_metadata",
        lambda acc: {"project_title": "T", "project_description": "D", "pubmed_id": "1"},
    )

    result, folder = _run(tmp_path, "--enrich-pride", "--project-accession", "PXD000001")

    assert result.exit_code == 0, result.output
    assert pq.read_table(folder / "t.dataset.parquet").to_pylist()[0]["project_title"] == "T"


@pytest.mark.parametrize("error", [ConnectionError("PRIDE unavailable"), IncompleteRead(b"")])
def test_enrich_pride_failure_is_a_warning(tmp_path, monkeypatch, error):
    """A PRIDE lookup failure, or a response cut off while read, prints a warning and the conversion still succeeds."""

    def _unavailable(acc):
        """Simulate a PRIDE API outage."""
        raise error

    monkeypatch.setattr("qpx.core.pride.fetch_pride_metadata", _unavailable)

    result, _ = _run(tmp_path, "--enrich-pride", "--project-accession", "PXD000001")

    assert result.exit_code == 0, result.output
    assert "WARNING: Could not enrich from PRIDE" in result.output


def test_enrich_pride_targets_the_converted_dataset_in_a_shared_folder(tmp_path, monkeypatch):
    """--enrich-pride updates the dataset of its own prefix, not another dataset in the same folder."""
    monkeypatch.setattr(
        "qpx.core.pride.fetch_pride_metadata",
        lambda acc: {"project_title": "T", "project_description": "D", "pubmed_id": "1"},
    )
    first, _ = _run(tmp_path, prefix="a")
    assert first.exit_code == 0, first.output

    result, folder = _run(tmp_path, "--enrich-pride", "--project-accession", "PXD000002")

    assert result.exit_code == 0, result.output
    assert pq.read_table(folder / "t.dataset.parquet").to_pylist()[0]["project_title"] == "T"
    assert pq.read_table(folder / "a.dataset.parquet").to_pylist()[0]["project_title"] != "T"


def test_apply_software_without_metadata_files_warns(tmp_path, capsys):
    """_apply_software warns and skips when the dataset or provenance file is missing."""
    _apply_software(tmp_path, "t", "nf-core/mhcquant", "3.1.0", None)

    assert "WARNING: no dataset/provenance metadata found; --software-name skipped" in capsys.readouterr().err


@pytest.mark.parametrize("command", _METADATA_COMMANDS)
def test_convert_command_lists_metadata_options(command):
    """Each metadata-enabled convert command lists the shared metadata options in its help."""
    result = CliRunner().invoke(convert, [command, "--help"])

    assert result.exit_code == 0, result.output
    for option in _METADATA_OPTIONS:
        assert option in result.output


@pytest.mark.parametrize("command", ["sdrf", "mz"])
def test_convert_command_without_metadata_options(command):
    """The sdrf and mz convert commands do not offer the shared metadata options."""
    result = CliRunner().invoke(convert, [command, "--help"])

    assert result.exit_code == 0, result.output
    assert "--software-name" not in result.output


def test_convert_diann_records_workflow_and_integrity(tmp_path):
    """convert diann records workflow provenance and verifiable integrity, and reports its MuData view."""
    folder = tmp_path / "out"
    result = _run_diann(
        folder,
        "--software-name",
        "quantmsdiann",
        "--software-version",
        "2.0.0",
        "--compute-integrity",
    )

    assert result.exit_code == 0, result.output
    dataset = pq.read_table(folder / "diann.dataset.parquet").to_pylist()[0]
    assert dataset["software_name"] == "quantmsdiann"
    steps = pq.read_table(folder / "diann.provenance.parquet").to_pylist()
    assert steps[0]["step_category"] == "workflow"
    assert steps[0]["tool_name"] == "quantmsdiann"
    assert "DIA-NN" in [step["tool_name"] for step in steps]
    with Dataset(str(folder)) as ds:
        assert ds.verify_integrity()["errors"] == []
    assert (folder / "diann.h5mu").exists() == ("MuData view: diann.h5mu" in result.output)


def test_openms_consensus_and_diann_finish_conversion_once_with_their_prefix(tmp_path, monkeypatch):
    """convert openms-consensus and convert diann each run the shared post-conversion step once, with their prefix."""
    calls = []

    def spy(output_folder, prefix, *args, **kwargs):
        """Record the prefix of each post-conversion call and run the real step."""
        calls.append(prefix)
        return _finish_conversion(output_folder, prefix, *args, **kwargs)

    monkeypatch.setattr("qpx.cli.convert.openms_consensus._finish_conversion", spy)
    monkeypatch.setattr("qpx.cli.convert.diann._finish_conversion", spy)

    consensus, _ = _run(tmp_path / "consensus", prefix="t")
    assert consensus.exit_code == 0, consensus.output
    assert calls == ["t"]

    diann = _run_diann(tmp_path / "diann")
    assert diann.exit_code == 0, diann.output
    assert calls == ["t", "diann"]
