"""relayout stores the name-value lists of one layout as structs, and QPX reads them back as lists."""

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.dataset import Dataset
from qpx.transforms.mantissa import round_table
from qpx.writers.psm import PsmWriter
from tests.conftest import make_psm_record, without_rounded_cv_values

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
DIANN = EXAMPLES / "diann" / "small"


def _invoke(*arguments):
    result = CliRunner().invoke(qpx_main, [str(argument) for argument in arguments], catch_exceptions=False)
    assert result.exit_code == 0, result.output
    return result.output


def _rows(folder, view, prefix="d", rounded_like=None):
    """The rows of *view* read through the Dataset API, without the ids relayout renumbers, in a fixed order.

    *rounded_like* names a relayouted folder whose rounding of measured values the rows take.
    """
    with Dataset(str(folder), file_prefix=prefix) as dataset:
        table = getattr(dataset, view).to_arrow()
    if rounded_like:
        footer = pq.read_metadata(rounded_like / f"{prefix}.{view}.parquet").metadata
        table = round_table(table, json.loads(footer[b"qpx.mantissa_bits"]))
    columns = [name for name in table.schema.names if name not in ("feature_id", "pg_id", "psm_id", "pg_ids", "psm_ids")]
    return sorted(map(repr, without_rounded_cv_values(table.select(columns).to_pylist(), view)))


@pytest.fixture(scope="module")
def diann(tmp_path_factory):
    root = tmp_path_factory.mktemp("diann")
    _invoke(
        "convert", "diann", "--report-path", DIANN / "diann_report.tsv", "--sdrf-file", DIANN / "PXD019909-DIA.sdrf.tsv",
        "--pg-matrix-path", DIANN / "diann_report.pg_matrix.tsv", "--output-folder", root / "source", "--output-prefix", "d",
    )  # fmt: skip
    _invoke("transform", "relayout", "--dataset", root / "source", "--output-folder", root / "sorted")
    return root


def test_relayout_stores_the_lists_of_one_layout_as_structs(diann):
    """Every feature lists the same scores and cv terms, so each is a field; a cv string is a float only if it reads back."""
    path = diann / "sorted" / "d.feature.parquet"
    schema = pq.read_schema(path)

    scores, cv = schema.field("additional_scores").type, schema.field("cv_params").type
    assert pa.types.is_struct(scores)
    assert scores.field("diann_cscore").type == pa.float64()
    # Proteotypic "1" would read back from a float as "1.0".
    assert (cv.field("irt").type, cv.field("proteotypic").type) == (pa.float32(), pa.string())
    # 118 features have no ms1_area entry, so their layout differs and the column stays a list.
    assert pa.types.is_list(schema.field("additional_intensities").type)
    assert b"qpx.struct_lists" in pq.read_metadata(path).metadata


@pytest.mark.parametrize("view", ["feature", "pg"])
def test_struct_lists_read_back_as_the_lists(diann, view):
    assert _rows(diann / "sorted", view) == _rows(diann / "source", view, rounded_like=diann / "sorted")


def _cv_irt(folder):
    """The cv iRT of each feature, by peptidoform, charge and run."""
    with Dataset(str(folder), file_prefix="d") as dataset:
        rows = dataset.feature.to_arrow().select(["peptidoform", "charge", "run_file_name", "cv_params"]).to_pylist()
    return {
        (row["peptidoform"], row["charge"], row["run_file_name"]): float(entry["cv_value"])
        for row in rows
        for entry in row["cv_params"] or []
        if entry["cv_name"] == "irt"
    }


def test_the_cv_irt_is_rounded_as_a_retention_time(diann):
    """Once cv_params is a struct, its iRT is a float, which keeps 16 mantissa bits."""
    footer = pq.read_metadata(diann / "sorted" / "d.feature.parquet").metadata
    assert json.loads(footer[b"qpx.mantissa_bits"])["cv_params.irt"] == 16

    source, rounded = _cv_irt(diann / "source"), _cv_irt(diann / "sorted")

    assert rounded.keys() == source.keys()
    assert any(rounded[key] != value for key, value in source.items())
    assert all(abs(rounded[key] - value) <= 2.0**-16 * abs(value) for key, value in source.items())


def test_a_dataset_with_struct_lists_validates(diann):
    assert "Overall: VALID" in _invoke("validate", "--dataset-path", diann / "sorted")


@pytest.mark.parametrize(
    "command",
    [
        ["transform", "gene-map", "--fasta", EXAMPLES / "fasta" / "Homo-sapiens.fasta", "--dataset"],
        ["migrate", "--dataset-path"],
    ],
    ids=["gene-map", "migrate"],
)
def test_rewriting_a_view_reads_its_struct_lists_back(diann, tmp_path, command):
    """Transforms that read the files themselves write the lists, as they do for the unsorted dataset."""
    for name in ("source", "sorted"):
        _invoke(*command, diann / name, "--output-folder", tmp_path / name)

    assert pa.types.is_list(pq.read_schema(tmp_path / "sorted" / "d.feature.parquet").field("additional_scores").type)
    assert _rows(tmp_path / "sorted", "feature") == _rows(tmp_path / "source", "feature", rounded_like=diann / "sorted")


def _psms(folder, scores):
    folder.mkdir()
    records = [
        {**make_psm_record(peptidoform=f"PEP{index}K", scan=[index]), "additional_scores": entries}
        for index, entries in enumerate(scores, start=1)
    ]
    with PsmWriter(folder / "ds.psm.parquet") as writer:
        writer.write_batch(records)
    return folder


def _score(name, value):
    return {"score_name": name, "score_value": value, "higher_better": name.startswith("s")}


@pytest.mark.parametrize(
    ("scores", "fields"),
    [
        # 0.5 is a float32, 0.1 is not; a null list is a null struct.
        (
            [[_score("score", 0.5), _score("q", 0.1)], None, [_score("score", 2.0), _score("q", 0.2)]],
            [pa.float32(), pa.float64()],
        ),
        ([[_score("score", 0.5), _score("q", 0.1)], [_score("q", 0.1), _score("score", 0.5)]], None),
        ([[_score("score", 0.5)], []], None),
        # DuckDB reads fields that differ only in case under other names.
        ([[_score("score", 0.5), _score("Score", 0.1)]], None),
    ],
    ids=["one layout", "two orders", "an empty list", "names differing in case"],
)
def test_only_a_layout_that_one_struct_can_hold_becomes_a_struct(tmp_path, scores, fields):
    _invoke("transform", "relayout", "--dataset", _psms(tmp_path / "ds", scores), "--output-folder", tmp_path / "sorted")

    stored = pq.read_schema(tmp_path / "sorted" / "ds.psm.parquet").field("additional_scores").type
    assert ([field.type for field in stored] if pa.types.is_struct(stored) else None) == fields
    assert _rows(tmp_path / "sorted", "psm", "ds") == _rows(tmp_path / "ds", "psm", "ds", rounded_like=tmp_path / "sorted")
