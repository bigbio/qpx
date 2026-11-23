"""Converters store fixed structs and QPX reads them back as name-value lists."""

import json
import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.collection import DatasetCollection
from qpx.core import engine as engine_module
from qpx.core.engine import DuckDBEngine
from qpx.core.struct_lists import STRUCT_LISTS_FOOTER_KEY, file_layouts
from qpx.dataset import Dataset
from qpx.transforms.gene_mapping import GeneMappingTransform
from qpx.transforms.mantissa import round_table
from qpx.writers.feature import FeatureWriter
from qpx.writers.psm import PsmWriter
from tests.conftest import make_feature_record, make_psm_record

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
    return sorted(map(repr, table.select(columns).to_pylist()))


@pytest.fixture(scope="module", name="diann")
def _diann(tmp_path_factory):
    root = tmp_path_factory.mktemp("diann")
    _invoke(
        "convert", "diann", "--report-path", DIANN / "diann_report.tsv", "--sdrf-file", DIANN / "PXD019909-DIA.sdrf.tsv",
        "--pg-matrix-path", DIANN / "diann_report.pg_matrix.tsv", "--output-folder", root / "source", "--output-prefix", "d",
    )  # fmt: skip
    _invoke("transform", "relayout", "--dataset", root / "source", "--output-folder", root / "sorted")
    return root


@pytest.mark.parametrize("view", ["feature", "pg"])
def test_relayout_keeps_fixed_structs_and_their_footer(diann, view):
    before, after = (pq.ParquetFile(diann / folder / f"d.{view}.parquet") for folder in ("source", "sorted"))
    for column in ("additional_scores", "cv_params"):
        field = before.schema_arrow.field(column)
        assert pa.types.is_struct(field.type)
        assert after.schema_arrow.field(column) == field
    assert pa.types.is_list(after.schema_arrow.field("additional_intensities").type)
    assert after.metadata.metadata[STRUCT_LISTS_FOOTER_KEY] == before.metadata.metadata[STRUCT_LISTS_FOOTER_KEY]
    assert _rows(diann / "sorted", view) == _rows(diann / "source", view, rounded_like=diann / "sorted")


def test_relayout_keeps_homogeneous_name_value_lists(tmp_path):
    source, target = tmp_path / "source", tmp_path / "sorted"
    source.mkdir()
    score = {"score_name": "andromeda_score", "score_value": 0.25, "higher_better": True}
    cv = {"cv_name": "mass analyzer", "cv_value": "FTMS"}
    with PsmWriter(source / "d.psm.parquet", creator="maxquant") as writer:
        writer.write_batch([{**make_psm_record(scan=[scan]), "additional_scores": [score], "cv_params": [cv]} for scan in (2, 1)])
    _invoke("transform", "relayout", "--dataset", source, "--output-folder", target)

    schema = pq.read_schema(target / "d.psm.parquet")
    assert STRUCT_LISTS_FOOTER_KEY not in schema.metadata
    for column in ("additional_scores", "cv_params"):
        assert pa.types.is_list(schema.field(column).type)
    assert _rows(target, "psm") == _rows(source, "psm", rounded_like=target)


def test_a_dataset_with_struct_lists_validates(diann):
    assert "Overall: VALID" in _invoke("validate", "--dataset", diann / "sorted")


def test_migrate_reads_structs_and_writes_lists(diann, tmp_path):
    for name in ("source", "sorted"):
        _invoke("migrate", "--dataset", diann / name, "--output-folder", tmp_path / name)
        for view in ("feature", "pg"):
            file = pq.ParquetFile(tmp_path / name / f"d.{view}.parquet")
            footer = file.metadata.metadata
            assert not file_layouts(footer)
            assert pa.types.is_list(file.schema_arrow.field("additional_scores").type)
            assert pa.types.is_list(file.schema_arrow.field("cv_params").type)
            assert _rows(tmp_path / name, view) == _rows(diann / name, view)
            original = pq.read_metadata(diann / name / f"d.{view}.parquet").metadata
            assert footer.get(b"qpx.mantissa_bits") == original.get(b"qpx.mantissa_bits")


def test_gene_map_dataset_preserves_structs_and_rounding(diann, tmp_path):
    for name in ("source", "sorted"):
        source, target = diann / name, tmp_path / name
        _invoke("transform", "gene-map", "--fasta", EXAMPLES / "fasta" / "Homo-sapiens.fasta",
                "--dataset", source, "--output-folder", target)  # fmt: skip
        for view in ("feature", "pg"):
            before, after = (pq.ParquetFile(folder / f"d.{view}.parquet") for folder in (source, target))
            original, footer = before.metadata.metadata, after.metadata.metadata
            assert file_layouts(footer) == file_layouts(original)
            for column in ("additional_scores", "cv_params"):
                assert after.schema_arrow.field(column) == before.schema_arrow.field(column)
            assert _rows(target, view) == _rows(source, view)
            assert footer.get(b"qpx.mantissa_bits") == original.get(b"qpx.mantissa_bits")


def test_gene_map_parquet_preserves_structs_and_file_footer(diann, tmp_path):
    source = diann / "sorted" / "d.pg.parquet"
    before = pq.read_table(source)
    fasta = tmp_path / "genes.fasta"
    accession = before["pg_accessions"][0][0].as_py()
    fasta.write_text(f">sp|{accession}|PROTEIN GN=MappedGene\nPEPTIDE\n")
    _invoke("transform", "gene-map", "--parquet-path", source, "--fasta", fasta, "--output-folder", tmp_path)

    target = tmp_path / source.name
    after = pq.read_table(target)
    columns = ["additional_scores", "cv_params"]
    assert after.select(columns).equals(before.select(columns))
    original, footer = (pq.read_metadata(path).metadata for path in (source, target))
    assert footer[STRUCT_LISTS_FOOTER_KEY] == original[STRUCT_LISTS_FOOTER_KEY]
    stamped = (b"ARROW:schema", b"uuid", b"creation_date")
    assert {key: value for key, value in footer.items() if key not in stamped} == {
        key: value for key, value in original.items() if key not in stamped
    }
    assert footer[b"uuid"] != original[b"uuid"]
    assert after["gg_names"][0].as_py() == ["MappedGene"]


def test_merge_of_fixed_structs_writes_lists(diann, tmp_path):
    """Merging a converter's fixed-layout input writes the logical lists."""
    with (
        Dataset(str(diann / "source"), file_prefix="d") as dataset,
        DatasetCollection([dataset]) as collection,
    ):
        collection.merge(tmp_path / "merged", structures=["feature", "pg"])
    for view in ("feature", "pg"):
        schema = pq.read_schema(tmp_path / "merged" / f"merged.{view}.parquet")
        assert not file_layouts(schema.metadata)
        for column in ("additional_scores", "cv_params"):
            assert pa.types.is_list(schema.field(column).type)


def test_fixed_writer_round_trips_null_empty_and_partial_lists(tmp_path):
    """The public writer schema stays logical, even for several streamed writes."""
    path = tmp_path / "test.psm.parquet"
    scores = [{"score_name": "andromeda_score", "score_value": 8.25, "higher_better": True}]
    cv = [{"cv_name": "mass analyzer", "cv_value": "FTMS"}]
    with PsmWriter(path, creator="maxquant", fixed_struct_lists=True, batch_size=1) as writer:
        assert pa.types.is_list(writer.arrow_schema.field("additional_scores").type)
        for scan, entries, params in ((1, scores, cv), (2, [], []), (3, None, None)):
            writer.write_batch([{**make_psm_record(scan=[scan]), "additional_scores": entries, "cv_params": params}])
        writer.set_scan_format("scan")

    schema = pq.read_schema(path)
    assert schema.field("additional_scores").type.field("andromeda_score").type == pa.float64()
    assert schema.field("cv_params").type.field("mass analyzer").type == pa.string()
    assert pq.read_metadata(path).metadata[b"scan_format"] == b"scan"
    assert pq.read_table(path)["additional_scores"][0].as_py() == {
        "andromeda_score": 8.25,
        "andromeda_delta_score": None,
        "parent_ion_fraction": None,
    }
    with Dataset(str(tmp_path), file_prefix="test") as dataset:
        table = dataset.psm.to_arrow()
    rows = {tuple(row["scan"]): row for row in table.to_pylist()}
    assert [(rows[(i,)]["additional_scores"], rows[(i,)]["cv_params"]) for i in (1, 2, 3)] == [
        (scores, cv),
        ([], []),
        (None, None),
    ]


_SCORES = [{"score_name": "andromeda_score", "score_value": 8.25, "higher_better": True}]
_CV = [{"cv_name": "mass analyzer", "cv_value": "FTMS"}]


def _psm_lists(folder, prefix=None):
    """The scores and cv parameters of each scan in the psm files of *folder* (only the file of *prefix*, if given)."""
    with Dataset(str(folder), file_prefix=prefix) as dataset:
        table = dataset.psm.to_arrow()
    assert all(pa.types.is_struct(table.schema.field(column).type.value_type) for column in ("additional_scores", "cv_params"))
    return {row["scan"][0]: (row["additional_scores"], row["cv_params"]) for row in table.to_pylist()}


@pytest.mark.parametrize("statistics", [True, False], ids=["statistics", "no statistics"])
def test_files_read_together_check_a_field_null_in_some_of_their_rows(tmp_path, statistics):
    """A field null in no row of one file and in every row of another is checked row by row when both are read.

    So is a field of a file without statistics, whatever the statistics of the other file show.
    """
    folder = tmp_path / "psm"
    folder.mkdir()
    for name, rows in (("a", [(1, _SCORES, _CV)]), ("b", [(2, [], []), (3, None, None)])):
        with PsmWriter(folder / f"{name}.psm.parquet", creator="maxquant", fixed_struct_lists=True) as writer:
            writer.write_batch(
                [{**make_psm_record(scan=[scan]), "additional_scores": scores, "cv_params": cv} for scan, scores, cv in rows]
            )
    if not statistics:
        path = folder / "b.psm.parquet"
        pq.write_table(pq.read_table(path), path, write_statistics=False)

    assert _psm_lists(folder) == {1: (_SCORES, _CV), 2: ([], []), 3: (None, None)}
    assert _psm_lists(folder, "b") == {2: ([], []), 3: (None, None)}


@pytest.mark.parametrize(
    ("creator", "names", "problem"),
    [
        ("maxquant", ["new_score"], "outside the fixed layout"),
        ("maxquant", ["andromeda_score", "andromeda_score"], "duplicate name 'andromeda_score'"),
    ],
    ids=["unknown name", "duplicate in a row"],
)
def test_fixed_writer_rejects_unrepresentable_scores(tmp_path, creator, names, problem):
    """Unknown and repeated scores cannot occupy the fields of the fixed struct."""
    path = tmp_path / "test.psm.parquet"
    scores = [{"score_name": name, "score_value": 1.25, "higher_better": True} for name in names]
    with (
        pytest.raises(ValueError, match=f"{creator} psm additional_scores: .*(?:{problem})"),
        PsmWriter(path, creator=creator, fixed_struct_lists=True) as writer,
    ):
        writer.write_batch([{**make_psm_record(), "additional_scores": scores}])
    assert not path.exists()


def test_fixed_writer_rejects_an_entry_in_an_empty_layout(tmp_path):
    """A column without declared fields stays a list that holds no entry."""
    path = tmp_path / "test.feature.parquet"
    scores = [{"score_name": "new_score", "score_value": 1.25, "higher_better": True}]
    with (
        pytest.raises(ValueError, match="fragpipe feature additional_scores: name 'new_score' is outside the fixed layout"),
        FeatureWriter(path, creator="fragpipe", fixed_struct_lists=True) as writer,
    ):
        writer.write_batch([{**make_feature_record(), "additional_scores": scores}])
    assert not path.exists()


def test_fixed_writer_rejects_a_null_value(tmp_path):
    """A null struct field reads back as no entry, so an entry without a value cannot be stored."""
    path = tmp_path / "test.psm.parquet"
    scores = [{"score_name": "andromeda_score", "score_value": None, "higher_better": True}]
    with (
        pytest.raises(ValueError, match="maxquant psm additional_scores: name 'andromeda_score' has a null value"),
        PsmWriter(path, creator="maxquant", fixed_struct_lists=True) as writer,
    ):
        writer.write_batch([{**make_psm_record(), "additional_scores": scores}])
    assert not path.exists()


_BETA = ("spectronaut", "spectronaut_cscore", "predicted_irt")
_SHARDS = {
    "alpha": ("diann", "diann_cscore", "predicted_irt"),
    "beta": _BETA,
    "beta_too": _BETA,
    "gamma": ("proteome_discoverer", "sequest_xcorr", None),
}


def _features(path, run, creator, score, term, *, fixed):
    records = []
    for index, value in enumerate(("52.822754", "1.1", "33.3"), start=1):
        record = make_feature_record(sequence=f"PEP{index}K", peptidoform=f"PEP{index}K", run_file_name=run, scan=[index])
        record["additional_scores"] = [{"score_name": score, "score_value": index / 10, "higher_better": True}]
        record["cv_params"] = [{"cv_name": term, "cv_value": value}] if term else None
        record["additional_intensities"] = [
            {"label": "LFQ", "intensities": [{"intensity_name": "spectral_count", "intensity_value": float(index)}]}
        ]
        records.append(record)
    with FeatureWriter(path, creator=creator, fixed_struct_lists=fixed) as writer:
        writer.write_batch(records)


@pytest.fixture(scope="module", name="shards")
def _shards(tmp_path_factory):
    """Each shard's feature file in list and fixed struct form, by shard name."""
    root = tmp_path_factory.mktemp("shards")
    files = {}
    for name, shard in _SHARDS.items():
        folder = root / name
        folder.mkdir()
        source, fixed = folder / "list.feature.parquet", folder / "fixed.feature.parquet"
        _features(source, f"run_{name}", *shard, fixed=False)
        _features(fixed, f"run_{name}", *shard, fixed=True)
        files[name] = source, fixed
    return files


@pytest.fixture(scope="module", name="legacy_shards")
def _legacy_shards(shards, tmp_path_factory):
    """One footer layout stored with two numeric leaf types by earlier development builds."""
    root = tmp_path_factory.mktemp("legacy-shards")
    files = {}
    for name, kind in (("beta", pa.float32()), ("beta_too", pa.float64())):
        table = pq.read_table(shards[name][1])
        for column in ("additional_scores", "cv_params"):
            values = table[column]
            storage = pa.struct([field.with_type(kind) for field in values.type])
            table = table.set_column(table.schema.get_field_index(column), column, values.cast(storage))
        path = root / f"{name}.feature.parquet"
        pq.write_table(table, path)
        files[name] = path
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


def _intensity_struct_file(source, path):
    """*source* with additional_intensities stored as a struct, as earlier development builds wrote it."""
    table = pq.read_table(source)
    values = [
        {"LFQ": {"spectral_count": row[0]["intensities"][0]["intensity_value"]}}
        for row in table["additional_intensities"].to_pylist()
    ]
    intensities = pa.array(values, pa.struct([("LFQ", pa.struct([("spectral_count", pa.float32())]))]))
    layouts = {**file_layouts(table.schema.metadata), "additional_intensities": [{"name": "LFQ", "names": ["spectral_count"]}]}
    metadata = {**table.schema.metadata, STRUCT_LISTS_FOOTER_KEY: json.dumps(layouts).encode()}
    table = table.set_column(table.schema.get_field_index("additional_intensities"), "additional_intensities", intensities)
    pq.write_table(table.replace_schema_metadata(metadata), path)
    return path


@pytest.mark.parametrize("case", ["mixed-layouts", "float-fields", "intensity-struct"])
def test_gene_map_writes_lists_for_layouts_it_cannot_keep(shards, legacy_shards, tmp_path, case):
    """Files of two layouts, or with an earlier build's FLOAT fields or intensity struct, are written as lists."""
    if case == "mixed-layouts":
        files = {"a": shards["alpha"][1], "b": shards["beta"][1]}
    elif case == "float-fields":
        files = {"a": legacy_shards["beta"]}
    else:
        files = {"a": _intensity_struct_file(shards["beta"][1], tmp_path / "intensities.parquet")}
    source = _shard_folder(tmp_path / "shards", **files)
    output = tmp_path / "annotated" / "d.feature.parquet"
    fasta = tmp_path / "genes.fasta"
    fasta.write_text(">sp|P00000|OTHER GN=OtherGene\nPEPTIDE\n")
    with Dataset(str(source)) as dataset:
        GeneMappingTransform(fasta).write_annotated_features(dataset, output)

    schema = pq.read_schema(output)
    assert STRUCT_LISTS_FOOTER_KEY not in schema.metadata
    for column in ("additional_scores", "cv_params", "additional_intensities"):
        assert pa.types.is_list(schema.field(column).type)
    assert _feature_rows(output.parent) == _feature_rows(source)


def _view_sql(folder):
    """The SQL of the feature view that the Dataset API registers for *folder*."""
    with Dataset(str(folder)) as dataset:
        return dataset._engine.execute("SELECT sql FROM duckdb_views() WHERE view_name = 'feature'").fetchall()[0][0]


@pytest.mark.parametrize("struct_first", [True, False], ids=["struct first", "list first"])
def test_a_view_of_list_and_struct_files_reads_each_row_as_its_file(shards, tmp_path, struct_first):
    """DuckDB cannot read a list and a struct of one column in one read_parquet, so each layout is read on its own."""
    source, fixed = shards["beta"]
    names = ("a", "b") if struct_first else ("b", "a")
    folder = _shard_folder(tmp_path / "mixed", **dict(zip(names, (fixed, source))))

    assert _feature_rows(folder) == _each_file_alone(folder)


def _names_by_run(folder):
    """The names of the scores, cv terms and additional intensities with a value in each run's features."""
    with Dataset(str(folder)) as dataset:
        rows = dataset.feature.to_arrow().select(["run_file_name", "additional_scores", "cv_params", "additional_intensities"])
    names = {}
    for row in rows.to_pylist():
        names.setdefault(row["run_file_name"], set()).update(
            [entry["score_name"] for entry in row["additional_scores"]]
            + [entry["cv_name"] for entry in row["cv_params"] or []]
            + [inner["intensity_name"] for entry in row["additional_intensities"] for inner in entry["intensities"]]
        )
    return names


@pytest.mark.parametrize("shard_names", [("alpha", "beta"), ("beta", "alpha"), ("alpha", "gamma", "beta")], ids="+".join)
def test_struct_files_of_different_layouts_keep_their_own_lists(shards, tmp_path, shard_names):
    """Each file's rows read back under its own names, not those of the first file's footer."""
    folder = _shard_folder(tmp_path / "mixed", **{name: shards[name][1] for name in shard_names})

    assert _feature_rows(folder) == _each_file_alone(folder)
    assert _names_by_run(folder) == {
        f"run_{name}": {score, "spectral_count"} | ({term} if term else set())
        for name, (_, score, term) in _SHARDS.items()
        if name in shard_names
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


def test_files_of_one_layout_stored_as_other_types_are_read_apart(legacy_shards, tmp_path):
    """One union_by_name read widened beta's FLOAT cv term to beta_double's DOUBLE: "52.822754" read "52.82275390625"."""
    folder = _shard_folder(tmp_path / "same", a=legacy_shards["beta"], b=legacy_shards["beta_too"])

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


def test_a_struct_file_is_registered_reading_only_its_layout(shards, monkeypatch):
    """A single file is one group: its view reads only its declared layout, as each footer read of a remote file is a request."""

    def unexpected(*arguments):
        raise AssertionError("the file was grouped")

    monkeypatch.setattr(engine_module, "layout_groups", unexpected)
    with DuckDBEngine() as engine:
        engine.register_parquet("feature", shards["beta"][1])
        assert engine.execute("SELECT count(*) FROM feature").fetchone()[0] == 3
        types = dict(engine.execute("SELECT column_name, column_type FROM (DESCRIBE feature)").fetchall())
    assert types["additional_scores"].endswith("[]")


def test_hive_partitions_of_one_layout_keep_their_stored_types(legacy_shards, tmp_path):
    """If the DOUBLE scores were read as FLOAT, 0.1 would become 0.10000000149011612."""
    folder = tmp_path / "partitioned"
    for name, source in legacy_shards.items():
        partition = folder / "feature" / f"run_file_name=run_{name}"
        partition.mkdir(parents=True)
        pq.write_table(pq.read_table(source).drop_columns(["run_file_name"]), partition / "part-0.parquet")
    alone = _shard_folder(tmp_path / "alone", **legacy_shards)

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
    """Each partition is the struct file of its run without the run column, which the directory names."""
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
    ids=["mixed cv layouts", "a list next to a struct"],
)
def test_validate_reads_a_view_of_files_with_different_layouts(shards, tmp_path, files):
    """Validation reads mixed layouts through the same restoration as queries."""
    folder = _shard_folder(tmp_path / "mixed", **{name: shards[shard][index] for name, (shard, index) in files.items()})

    result = CliRunner().invoke(qpx_main, ["validate", "--dataset", str(folder)])

    assert not isinstance(result.exception, Exception), result.exc_info
    assert "could not be read" not in result.output
    assert "Overall:" in result.output


@pytest.mark.parametrize(
    "others",
    [{"b": ("beta", 1)}, {"b": ("beta", 0)}, {"b": ("beta", 0), "c": ("gamma", 1)}, {"0": ("beta", 0), "c": ("gamma", 1)}],
    ids=[
        "next to a struct file",
        "next to a list file",
        "next to a list and a struct file",
        "after a list file",
    ],
)
def test_a_struct_without_its_layout_next_to_other_files_is_refused(shards, tmp_path, others):
    """A struct with a missing footer layout prevents mixed-file registration."""
    table = pq.read_table(shards["alpha"][1])
    folder = _shard_folder(tmp_path / "mixed", **{name: shards[shard][index] for name, (shard, index) in others.items()})
    footer = {key: value for key, value in table.schema.metadata.items() if key != b"qpx.struct_lists"}
    pq.write_table(table.replace_schema_metadata(footer), folder / "a.feature.parquet")

    with DuckDBEngine() as engine:
        with pytest.raises(ValueError, match="a.feature.parquet holds additional_scores, cv_params as a struct without"):
            engine.register_parquet_files("feature", sorted(folder.glob("*.feature.parquet")))
        assert engine.execute("SELECT count(*) FROM duckdb_views() WHERE view_name = 'feature'").fetchone()[0] == 0
    result = CliRunner().invoke(qpx_main, ["validate", "--dataset", str(folder)])
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
