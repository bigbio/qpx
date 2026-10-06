"""Regression coverage for missing intensities in both matrix export paths."""

import pandas as pd
import pytest

from qpx.dataset import Dataset
from qpx.writers import PgWriter, RunWriter
from tests.conftest import make_pg_record, make_run_record


@pytest.fixture(name="nullable_dataset", params=[False, True], ids=["mixed", "all_null"])
def _nullable_dataset(tmp_path, request):
    """Write duplicated, missing and zero intensities across two samples."""
    runs = []
    for run_name, sample in [("run_01", "SAMPLE_01"), ("run_02", "SAMPLE_01"), ("run_03", "SAMPLE_02")]:
        run = make_run_record(run_accession=run_name, run_file_name=run_name)
        run["samples"][0]["sample_accession"] = sample
        runs.append(run)
    with RunWriter(tmp_path / "nullable.run.parquet") as writer:
        writer.write_batch(runs)

    observations = [
        ("run_01", "AAAAR", None),
        ("run_02", "AAAAR", None),
        ("run_03", "AAAAR", None),
        ("run_01", "PEPTIDEK", 100.0),
        ("run_02", "PEPTIDEK", None),
        ("run_01", "ALVPEPK", 2.0),
        ("run_02", "ALVPEPK", 3.0),
        ("run_01", "VVALK", 0.0),
        ("run_03", "VVALK", None),
    ]
    proteins = []
    for run_name, identifier, intensity in observations:
        intensities = [{"label": "TMT126", "intensity": None if request.param else intensity}]
        proteins.append(make_pg_record(anchor_protein=identifier, run_file_name=run_name, intensities=intensities))
    with PgWriter(tmp_path / "nullable.pg.parquet") as writer:
        writer.write_batch(proteins)

    expected = pd.DataFrame(
        [[float("nan"), 5.0, 100.0, 0.0], [float("nan")] * 4],
        columns=["AAAAR", "ALVPEPK", "PEPTIDEK", "VVALK"],
        index=pd.Index(["SAMPLE_01", "SAMPLE_02"], name="sample_accession"),
    )
    if request.param:
        expected.loc[:, :] = float("nan")
    return tmp_path, expected


@pytest.mark.parametrize("fillna", [None, 0.0, -1.0])
def test_missing_intensities_match_parquet(nullable_dataset, fillna):
    """Keep null-only axes and SQL sum semantics before explicit filling."""
    dataset_path, expected = nullable_dataset
    output = dataset_path / "matrix.parquet"
    with Dataset(dataset_path) as dataset:
        matrix = dataset.design_matrix(fillna=fillna)
        dataset.design_matrix(fillna=fillna, output_path=output)

    on_disk = pd.read_parquet(output).set_index("sample_accession").sort_index().sort_index(axis=1)
    pd.testing.assert_frame_equal(on_disk, expected)
    if fillna is not None:
        expected = expected.fillna(fillna)
    pd.testing.assert_frame_equal(matrix, expected, check_dtype=False)


def test_missing_intensities_default_to_zero(nullable_dataset):
    """The default still explicitly fills stored and absent quantities with zero."""
    dataset_path, expected = nullable_dataset
    with Dataset(dataset_path) as dataset:
        matrix = dataset.design_matrix()
    pd.testing.assert_frame_equal(matrix, expected.fillna(0.0), check_dtype=False)
