"""relayout stores the name-value lists of one layout as structs, and QPX reads them back as lists."""

import json
import shutil
from itertools import chain
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.core import engine as engine_module
from qpx.core.engine import DuckDBEngine
from qpx.dataset import Dataset
from qpx.transforms.mantissa import round_table
from qpx.writers.feature import FeatureWriter
from qpx.writers.psm import PsmWriter
from tests.conftest import make_feature_record, make_psm_record, without_rounded_cv_values

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
    """Transforms that read the files themselves write the lists, as they do for the unsorted dataset.

    The values stay rounded, and the footer still says how.
    """
    for name in ("source", "sorted"):
        _invoke(*command, diann / name, "--output-folder", tmp_path / name)

    assert pa.types.is_list(pq.read_schema(tmp_path / "sorted" / "d.feature.parquet").field("additional_scores").type)
    assert _rows(tmp_path / "sorted", "feature") == _rows(tmp_path / "source", "feature", rounded_like=tmp_path / "sorted")
    for view in ("feature", "pg"):
        rounded = pq.read_metadata(diann / "sorted" / f"d.{view}.parquet").metadata[b"qpx.mantissa_bits"]
        assert pq.read_metadata(tmp_path / "sorted" / f"d.{view}.parquet").metadata[b"qpx.mantissa_bits"] == rounded


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


# The scores, cv terms and additional intensities (label to intensity names) of each
# feature of a shard, and the cv values of its three features where they read as
# numbers; a ragged shard drops an intensity from its first feature, so that its
# additional_intensities keeps two layouts and stays a list.
_BETA = {
    "scores": ["beta_score"],
    "cv": ["beta_term"],
    "intensities": {"LFQ": ["beta_raw", "beta_norm"]},
    "cv_values": ["52.822754", "1.1", "33.3"],
}
_SHARDS = {
    "alpha": {"scores": ["alpha_score", "alpha_q"], "cv": ["alpha_term"], "intensities": {"TMT126": ["alpha_norm"]}},
    "beta": _BETA,
    "beta_too": _BETA,
    # beta's layout; 0.1 and 0.1234567891 are no float32, so its scores and cv term
    # are stored as DOUBLE where beta's are FLOAT.
    "beta_double": {**_BETA, "step": 0.1, "cv_values": ["0.1234567891", "2.5", "3.75"]},
    "gamma": {"scores": ["gamma_score"], "cv": ["gamma_term"], "intensities": {"TMT127": ["gamma_a", "gamma_b"]}, "ragged": True},
}


def _features(folder, run, scores, cv, intensities, cv_values=None, step=0.25, ragged=False):
    folder.mkdir()
    records = []
    for index in range(1, 4):
        record = make_feature_record(sequence=f"PEP{index}K", peptidoform=f"PEP{index}K", run_file_name=run, scan=[index])
        record["additional_scores"] = [_score(name, index * step) for name in scores]
        value = cv_values[index - 1] if cv_values else f"{run}_{index}"
        record["cv_params"] = [{"cv_name": name, "cv_value": value} for name in cv]
        record["additional_intensities"] = [
            {
                "label": label,
                "intensities": [
                    {"intensity_name": name, "intensity_value": 10.0 * index + position}
                    for position, name in enumerate(names[: 1 if ragged and index == 1 else None])
                ],
            }
            for label, names in intensities.items()
        ]
        records.append(record)
    with FeatureWriter(folder / "ds.feature.parquet") as writer:
        writer.write_batch(records)
    return folder


@pytest.fixture(scope="module")
def shards(tmp_path_factory):
    """Each shard's feature file as converted (list) and after relayout (struct), by shard name."""
    root = tmp_path_factory.mktemp("shards")
    files = {}
    for name, shard in _SHARDS.items():
        source = _features(root / name, f"run_{name}", **shard)
        _invoke("transform", "relayout", "--dataset", source, "--output-folder", root / f"{name}_sorted")
        files[name] = source / "ds.feature.parquet", root / f"{name}_sorted" / "ds.feature.parquet"
    return files


def _shard_folder(folder, **files):
    """A dataset folder holding each of *files* as ``<name>.feature.parquet``."""
    folder.mkdir()
    for name, path in files.items():
        shutil.copy(path, folder / f"{name}.feature.parquet")
    return folder


def _feature_rows(folder, prefix=None):
    """The feature rows of *folder* (only the file of *prefix*, if given) as the Dataset API reads them, in a fixed order."""
    with Dataset(str(folder), file_prefix=prefix) as dataset:
        rows = dataset.feature.to_arrow().to_pylist()
    return sorted(json.dumps(row, sort_keys=True, default=str) for row in rows)


def _each_file_alone(folder):
    """The rows of every feature file of *folder*, each file read on its own."""
    return sorted(row for path in folder.glob("*.feature.parquet") for row in _feature_rows(folder, path.name.split(".")[0]))


def _view_sql(folder):
    """The SQL of the feature view that the Dataset API registers for *folder*."""
    with Dataset(str(folder)) as dataset:
        return dataset._engine.execute("SELECT sql FROM duckdb_views() WHERE view_name = 'feature'").fetchall()[0][0]


def test_shards_store_their_lists_as_the_tests_need(shards):
    stored = {name: pq.read_schema(sorted_file) for name, (_, sorted_file) in shards.items()}
    columns = ("additional_scores", "cv_params", "additional_intensities")

    assert all(pa.types.is_struct(stored[name].field(column).type) for name in ("alpha", "beta") for column in columns)
    assert pa.types.is_list(stored["gamma"].field("additional_intensities").type)
    assert pa.types.is_struct(stored["gamma"].field("additional_scores").type)
    for name, kind in (("beta", pa.float32()), ("beta_double", pa.float64())):
        assert stored[name].field("additional_scores").type.field("beta_score").type == kind
        assert stored[name].field("cv_params").type.field("beta_term").type == kind


@pytest.mark.parametrize("sorted_first", [True, False], ids=["relayout first", "converted first"])
def test_a_view_of_converted_and_relayouted_files_reads_each_row_as_its_file(shards, tmp_path, sorted_first):
    """DuckDB cannot read a list and a struct of one column in one read_parquet, so each layout is read on its own."""
    source, relayouted = shards["beta"]
    names = ("a", "b") if sorted_first else ("b", "a")
    folder = _shard_folder(tmp_path / "mixed", **dict(zip(names, (relayouted, source))))

    assert _feature_rows(folder) == _each_file_alone(folder)


def _names_by_run(folder):
    """The names of the scores, cv terms and additional intensities with a value in each run's features."""
    with Dataset(str(folder)) as dataset:
        rows = dataset.feature.to_arrow().select(["run_file_name", "additional_scores", "cv_params", "additional_intensities"])
    names = {}
    for row in rows.to_pylist():
        names.setdefault(row["run_file_name"], set()).update(
            [entry["score_name"] for entry in row["additional_scores"] if entry["score_value"] is not None]
            + [entry["cv_name"] for entry in row["cv_params"] if entry["cv_value"] is not None]
            + [inner["intensity_name"] for entry in row["additional_intensities"] for inner in entry["intensities"]]
        )
    return names


@pytest.mark.parametrize("shard_names", [("alpha", "beta"), ("beta", "alpha"), ("alpha", "gamma", "beta")], ids="+".join)
def test_relayouted_files_of_different_layouts_keep_their_own_lists(shards, tmp_path, shard_names):
    """Each file's rows read back under its own names, not those of the first file's footer."""
    folder = _shard_folder(tmp_path / "mixed", **{name: shards[name][1] for name in shard_names})

    assert _feature_rows(folder) == _each_file_alone(folder)
    assert _names_by_run(folder) == {
        f"run_{name}": {*_SHARDS[name]["scores"], *_SHARDS[name]["cv"], *chain(*_SHARDS[name]["intensities"].values())}
        for name in shard_names
    }
    assert "UNION ALL BY NAME" in _view_sql(folder)


def test_an_s3_glob_of_different_layouts_reads_each_row_as_its_file(shards, tmp_path):
    """register_s3_parquet reads a glob with union_by_name, as the local shards are read."""
    folder = _shard_folder(tmp_path / "mixed", a=shards["alpha"][1], b=shards["gamma"][1])
    with DuckDBEngine() as engine:
        engine.register_s3_parquet("feature", str(folder / "*.feature.parquet"))
        rows = engine.connection.sql("SELECT * FROM feature").to_arrow_table().to_pylist()

    assert sorted(json.dumps(row, sort_keys=True, default=str) for row in rows) == _each_file_alone(folder)


def test_files_of_one_layout_and_one_storage_are_read_by_one_read_parquet(shards, tmp_path):
    folder = _shard_folder(tmp_path / "same", a=shards["beta"][1], b=shards["beta_too"][1])

    sql = _view_sql(folder)
    assert sql.count("read_parquet(") == 1
    assert "UNION" not in sql
    assert _feature_rows(folder) == _each_file_alone(folder)


def test_files_of_one_layout_stored_as_other_types_are_read_apart(shards, tmp_path):
    """One union_by_name read widened beta's FLOAT cv term to beta_double's DOUBLE: "52.822754" read "52.82275390625"."""
    folder = _shard_folder(tmp_path / "same", a=shards["beta"][1], b=shards["beta_double"][1])

    assert _view_sql(folder).count("read_parquet(") == 2
    assert _feature_rows(folder) == _each_file_alone(folder)
    assert '"cv_value": "52.822754"' in "".join(_feature_rows(folder))


def test_a_file_without_structs_is_registered_without_reading_its_footer(shards, monkeypatch):
    def unexpected(*arguments):
        raise AssertionError("the footer was read")

    monkeypatch.setattr(engine_module, "layout_groups", unexpected)
    monkeypatch.setattr(engine_module, "footer_layouts", unexpected)
    with DuckDBEngine() as engine:
        engine.register_parquet("feature", shards["beta"][0])
        assert engine.execute("SELECT count(*) FROM feature").fetchone()[0] == 3


def test_a_relayouted_file_is_registered_reading_only_its_layout(shards, monkeypatch):
    """A single file is one group: its view reads only its declared layout, as each footer read of a remote file is a request."""

    def unexpected(*arguments):
        raise AssertionError("the file was grouped")

    monkeypatch.setattr(engine_module, "layout_groups", unexpected)
    with DuckDBEngine() as engine:
        engine.register_parquet("feature", shards["beta"][1])
        assert engine.execute("SELECT count(*) FROM feature").fetchone()[0] == 3
        types = dict(engine.execute("SELECT column_name, column_type FROM (DESCRIBE feature)").fetchall())
    assert types["additional_scores"].endswith("[]")


def test_hive_partitions_of_one_layout_keep_their_stored_types(shards, tmp_path):
    """beta stores its scores as FLOAT and beta_double as DOUBLE; read with beta's types, 0.1 became 0.10000000149011612."""
    folder = tmp_path / "partitioned"
    for name in ("beta", "beta_double"):
        partition = folder / "feature" / f"run_file_name=run_{name}"
        partition.mkdir(parents=True)
        pq.write_table(pq.read_table(shards[name][1]).drop_columns(["run_file_name"]), partition / "part-0.parquet")
    alone = _shard_folder(tmp_path / "alone", **{name: shards[name][1] for name in ("beta", "beta_double")})

    assert _feature_rows(folder) == _each_file_alone(alone)


def test_hive_partitions_are_unioned_by_name(shards, tmp_path):
    """Read with the types of the first partition, beta's, the rt that beta_too stores as DOUBLE read 0.10000000149011612."""
    folder = tmp_path / "partitioned"
    for name in ("beta", "beta_too"):
        table = pq.read_table(shards[name][1]).drop_columns(["run_file_name"])
        if name == "beta_too":
            table = table.set_column(table.schema.get_field_index("rt"), "rt", pa.array([0.1] * table.num_rows, pa.float64()))
        partition = folder / f"run_file_name=run_{name}"
        partition.mkdir(parents=True)
        pq.write_table(table, partition / "part-0.parquet")

    with DuckDBEngine() as engine:
        engine.register_partitioned_parquet("feature", folder)
        assert engine.execute("SELECT DISTINCT rt FROM feature WHERE run_file_name = 'run_beta_too'").fetchall() == [(0.1,)]


def test_hive_partitions_of_different_layouts_keep_the_partition_column(shards, tmp_path):
    """Each partition is the relayouted file of its run without the run column, which the directory names."""
    folder = tmp_path / "partitioned"
    for name in ("alpha", "gamma", "beta"):
        table = pq.read_table(shards[name][1])
        partition = folder / "feature" / f"run_file_name=run_{name}"
        partition.mkdir(parents=True)
        pq.write_table(table.drop_columns(["run_file_name"]), partition / "part-0.parquet")
    alone = _shard_folder(tmp_path / "alone", **{name: shards[name][1] for name in ("alpha", "gamma", "beta")})

    assert _feature_rows(folder) == _each_file_alone(alone)
    with Dataset(str(folder)) as dataset:
        runs = dataset._engine.execute("SELECT DISTINCT run_file_name FROM feature ORDER BY 1").fetchall()
    assert runs == [("run_alpha",), ("run_beta",), ("run_gamma",)]


@pytest.mark.parametrize(
    "files",
    [{"a": ("gamma", 1), "b": ("beta", 1)}, {"a": ("beta", 1), "b": ("alpha", 0)}],
    ids=["a list next to a struct", "a converted file"],
)
def test_validate_reads_a_view_of_files_with_different_layouts(shards, tmp_path, files):
    """Applying gamma's footer to every file made DuckDB fail on beta's struct additional_intensities, and validate crash;
    a converted file made the view unreadable.
    """
    folder = _shard_folder(tmp_path / "mixed", **{name: shards[shard][index] for name, (shard, index) in files.items()})

    result = CliRunner().invoke(qpx_main, ["validate", "--dataset-path", str(folder)])

    assert not isinstance(result.exception, Exception), result.exc_info
    assert "could not be read" not in result.output
    assert "Overall:" in result.output


@pytest.mark.parametrize(
    "others",
    [{"b": ("beta", 1)}, {"b": ("beta", 0)}, {"b": ("beta", 0), "c": ("gamma", 1)}, {"0": ("beta", 0), "c": ("gamma", 1)}],
    ids=[
        "next to a relayouted file",
        "next to a converted file",
        "next to a converted and a relayouted file",
        "after a converted file",
    ],
)
def test_a_struct_without_its_layout_next_to_other_files_is_refused(shards, tmp_path, others):
    """Its rows cannot be read as lists, so the view is refused instead of failing on its first read. A converted
    file has no footer key either, and DuckDB reads a list next to a struct as the list, so each file's schema tells.
    """
    table = pq.read_table(shards["alpha"][1])
    folder = _shard_folder(tmp_path / "mixed", **{name: shards[shard][index] for name, (shard, index) in others.items()})
    footer = {key: value for key, value in table.schema.metadata.items() if key != b"qpx.struct_lists"}
    pq.write_table(table.replace_schema_metadata(footer), folder / "a.feature.parquet")

    with DuckDBEngine() as engine:
        with pytest.raises(ValueError, match="a.feature.parquet holds additional_scores, cv_params, additional_intensities"):
            engine.register_parquet_files("feature", sorted(folder.glob("*.feature.parquet")))
        assert engine.execute("SELECT count(*) FROM duckdb_views() WHERE view_name = 'feature'").fetchone()[0] == 0
    result = CliRunner().invoke(qpx_main, ["validate", "--dataset-path", str(folder)])
    assert not isinstance(result.exception, Exception), result.exc_info
    assert "Structure 'feature' could not be read" in result.output


@pytest.mark.parametrize("one_file", [True, False], ids=["alone", "next to another file"])
def test_a_view_that_cannot_be_restored_is_not_left_registered(shards, tmp_path, one_file):
    """The view over the files as stored would read their structs, merged by name, as if nothing had failed."""
    table = pq.read_table(shards["alpha"][1])
    folder = _shard_folder(tmp_path / "corrupt", **({} if one_file else {"b": shards["beta"][1]}))
    pq.write_table(
        table.replace_schema_metadata({**table.schema.metadata, b"qpx.struct_lists": b"{not json"}), folder / "a.feature.parquet"
    )

    with DuckDBEngine() as engine:
        with pytest.raises(ValueError):
            engine.register_parquet_files("feature", sorted(folder.glob("*.feature.parquet")))
        assert engine.execute("SELECT count(*) FROM duckdb_views() WHERE view_name = 'feature'").fetchone()[0] == 0
