"""Tests for the DIA-NN XIC -> QPX xic converter."""

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.convert.diann import convert_diann_xic_cmd
from qpx.converters.diann.xic_adapter import convert_diann_xic
from qpx.core.data import XicSchema

_DIANN_XIC_SCHEMA = pa.schema(
    [
        pa.field("pr", pa.string(), nullable=False),
        pa.field("feature", pa.string(), nullable=False),
        pa.field("info", pa.int32(), nullable=False),
        pa.field("rt", pa.float32(), nullable=False),
        pa.field("value", pa.float32(), nullable=False),
    ]
)


def _write_xic(path: Path, values: tuple[float, ...] = (10.0, 5.0, 9.0, 7.0)):
    rows = [
        {"pr": "PEPTIDEK2", "feature": "ms1", "info": 0, "rt": 1.0, "value": values[0]},
        {"pr": "PEPTIDEK2", "feature": "y7^1", "info": 0, "rt": 1.0, "value": values[1]},
        {"pr": "PEPTIDEK2", "feature": "y7^1", "info": 0, "rt": 1.1, "value": values[2]},
        {"pr": "PEPTM(UniMod:35)IDEK3", "feature": "b4^2", "info": 0, "rt": 2.0, "value": values[3]},
        # DIA-NN internal chunk index: not a measurement.
        {"pr": "PEPTIDEK2", "feature": "index", "info": 0, "rt": 0.0, "value": 0.0},
    ]
    pq.write_table(pa.Table.from_pylist(rows, schema=_DIANN_XIC_SCHEMA), path)


def test_convert_diann_xic_merges_runs_and_annotates_transitions(tmp_path):
    first = tmp_path / "run1.xic.parquet"
    second = tmp_path / "run2.xic.parquet"
    _write_xic(first)
    _write_xic(second, values=(1.0, 2.0, 3.0, 4.0))

    output = tmp_path / "dataset.xic.parquet"
    traces = convert_diann_xic([first, second], output)
    assert traces == 6

    table = pq.read_table(output)
    assert XicSchema.validate_full(table, strict=True).issues == []
    assert set(table.column("run_file_name").to_pylist()) == {"run1", "run2"}
    assert "index" not in set(table.column("transition").to_pylist())
    metadata = pq.read_schema(output).metadata
    assert metadata[b"xic_selection"] == b"targets"
    assert metadata[b"xic_window"] == b"apex"
    assert metadata[b"xic_rt_unit"] == b"seconds"
    assert metadata[b"transition_format"] == b"mzPAF"

    by_id = {row["xic_id"]: row for row in table.to_pylist()}
    ms1 = by_id["run1|PEPTIDEK2|ms1"]
    assert ms1["transition_type"] == "ms1"
    assert ms1["peptidoform"] == "PEPTIDEK"
    assert ms1["charge"] == 2
    assert ms1["isotope"] == 0
    assert ms1["feature_ids"] is None
    assert ms1["feature_scores"] is None
    assert ms1["target_mz"] is None
    assert ms1["rt"] == [60.0]
    assert ms1["intensity"] == [10.0]

    fragment = by_id["run1|PEPTIDEK2|y7^1"]
    assert fragment["transition_type"] == "fragment"
    assert fragment["fragment_type"] == "y"
    assert fragment["fragment_number"] == 7
    assert fragment["fragment_charge"] == 1
    assert fragment["rt"] == pytest.approx([60.0, 66.0], rel=1e-6)
    assert fragment["intensity"] == pytest.approx([5.0, 9.0])

    modified = by_id["run1|PEPTM(UniMod:35)IDEK3|b4^2"]
    assert modified["peptidoform"] == "PEPTM[UNIMOD:35]IDEK"
    assert modified["charge"] == 3
    assert modified["fragment_type"] == "b"
    assert modified["fragment_number"] == 4
    assert modified["fragment_charge"] == 2


def test_convert_diann_xic_rejects_missing_columns(tmp_path):
    path = tmp_path / "broken.xic.parquet"
    pq.write_table(pa.table({"pr": ["PEPTIDEK2"]}), path)

    with pytest.raises(ValueError, match="missing DIA-NN XIC column"):
        convert_diann_xic([path], tmp_path / "out.xic.parquet")


def test_diann_xic_cli(tmp_path):
    xic_folder = tmp_path / "xic"
    xic_folder.mkdir()
    _write_xic(xic_folder / "run1.xic.parquet")

    output = tmp_path / "out"
    result = CliRunner().invoke(
        convert_diann_xic_cmd,
        ["--xic-folder", str(xic_folder), "--output-folder", str(output), "--output-prefix", "PXD000001"],
    )
    assert result.exit_code == 0, result.output
    assert "3 traces" in result.output
    assert (output / "PXD000001.xic.parquet").exists()
