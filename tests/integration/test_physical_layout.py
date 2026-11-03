"""Written files have row groups, a page index and encodings for selective reads, and relayout sorts them."""

import random

import pyarrow as pa
import pyarrow.parquet as pq
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.dataset import Dataset
from qpx.writers.encoding import ROW_GROUP_BYTES, ROW_GROUP_ROWS, SampledParquetWriter, parquet_write_options, row_group_rows
from qpx.writers.ontology import OntologyWriter

REPORT = "tests/examples/spectronaut/small/spectronaut_report.tsv"


def _convert(tmp_path, *arguments):
    tmp_path.mkdir(exist_ok=True)
    out = tmp_path / "ds"
    result = CliRunner().invoke(
        qpx_main,
        ["convert", "spectronaut", "--report-path", REPORT, "--output-folder", str(out), "--output-prefix", "sn", *arguments],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    return out


def test_row_groups_are_capped_by_rows_and_bytes():
    small = pa.table({"a": pa.array(range(10), pa.int64())})
    large = pa.table({"a": pa.array(["x" * (10 * 1024 * 1024)] * 100)})

    assert row_group_rows(small) == ROW_GROUP_ROWS
    assert row_group_rows(large) * 10 * 1024 * 1024 <= ROW_GROUP_BYTES


def test_write_options_take_the_encoding_each_column_measures_smallest():
    """Each column gets the encoding that writes its sampled values smallest."""
    rng = random.Random(0)
    rows = 20_000
    sample = pa.table(
        {
            "feature_id": pa.array(range(rows), pa.int64()),
            "pg_id": pa.array([rng.getrandbits(63) for _ in range(rows)], pa.int64()),
            "rt": pa.array([rng.uniform(0, 7200) for _ in range(rows)], pa.float64()),
            "label": pa.array([rng.choice(["TMT126", "TMT127N", "TMT127C"]) for _ in range(rows)]),
        }
    )

    options = parquet_write_options(sample.schema, "zstd", 3, sample=sample)

    # Ordinals delta-encode, hashes stay plain, continuous floats split, few distinct strings keep a dictionary.
    assert options["column_encoding"] == {"feature_id": "DELTA_BINARY_PACKED", "rt": "BYTE_STREAM_SPLIT"}
    assert options["use_dictionary"] == ["label"]
    assert options["compression_level"] == 3
    assert options["write_page_index"] is True


def test_write_options_without_rows_keep_the_default_dictionaries():
    """Without sampled rows the pyarrow default encodings stay."""
    options = parquet_write_options(pa.schema([("rt", pa.float64())]), "zstd")

    assert "use_dictionary" not in options
    assert "column_encoding" not in options


def test_sampled_writer_without_tables_writes_an_empty_file(tmp_path):
    """A writer closed before any table still leaves a valid empty file."""
    schema = pa.schema([("rt", pa.float64())])

    with SampledParquetWriter(tmp_path / "empty.parquet", schema, "zstd"):
        pass

    assert pq.read_table(tmp_path / "empty.parquet").schema.equals(schema)
    assert pq.read_metadata(tmp_path / "empty.parquet").num_rows == 0


def test_writers_produce_a_page_index(tmp_path):
    path = tmp_path / "ontology.parquet"
    with OntologyWriter(path) as writer:
        writer.write_batch([{"field_name": "charge", "view": "feature"}])

    column = pq.ParquetFile(path).metadata.row_group(0).column(0)

    assert column.has_offset_index


def test_rows_are_ordered_so_conversions_compare_equal(tmp_path):
    out = _convert(tmp_path)

    first = pq.read_table(out / "sn.pg.parquet")
    second = pq.read_table(_convert(tmp_path / "again").joinpath("sn.pg.parquet"))
    assert first.equals(second)
    assert pq.read_table(out / "sn.feature.parquet").equals(pq.read_table(tmp_path / "again" / "ds" / "sn.feature.parquet"))


def test_relayout_sorts_rows_keeps_values_and_refreshes_integrity(tmp_path):
    source = _convert(tmp_path, "--compute-integrity")
    target = tmp_path / "sorted"

    result = CliRunner().invoke(
        qpx_main,
        ["transform", "relayout", "--dataset", str(source), "--output-folder", str(target), "--row-group-size", "50"],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    before = pq.read_table(source / "sn.feature.parquet")
    after_file = pq.ParquetFile(target / "sn.feature.parquet")
    after = after_file.read()
    anchors = after.column("anchor_protein").to_pylist()
    assert anchors == sorted(anchors)
    assert sorted(after.column("feature_id").to_pylist()) == sorted(before.column("feature_id").to_pylist())
    assert after.schema.equals(before.schema)
    assert after.schema.metadata[b"qpx.sort_order"] == b"anchor_protein,run_file_name"
    assert after_file.metadata.row_group(0).num_rows == 50
    sorting = after_file.metadata.row_group(0).sorting_columns
    assert [after_file.schema.column(column.column_index).path for column in sorting] == ["anchor_protein", "run_file_name"]
    with Dataset(str(target), file_prefix="sn") as dataset:
        assert not any(dataset.verify_integrity().values())


def test_relayout_refuses_a_column_it_cannot_sort_by(tmp_path):
    source = _convert(tmp_path)

    result = CliRunner().invoke(
        qpx_main,
        ["transform", "relayout", "--dataset", str(source), "--in-place", "--view", "feature", "--sort-by", "intensities"],
    )

    assert result.exit_code != 0
    assert "no flat column named intensities" in result.output
