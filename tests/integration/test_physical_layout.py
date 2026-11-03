"""Written files have row groups, a page index and encodings for selective reads, and relayout sorts them."""

import random

import pyarrow as pa
import pyarrow.parquet as pq
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.dataset import Dataset
from qpx.writers.encoding import ROW_GROUP_BYTES, ROW_GROUP_ROWS, SampledParquetWriter, parquet_write_options, row_group_rows
from qpx.writers.feature import FeatureWriter
from qpx.writers.ontology import OntologyWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter
from tests.conftest import make_feature_record, make_pg_record, make_psm_record

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


def test_relayout_sorts_rows_renumbers_ids_and_refreshes_integrity(tmp_path):
    source = _convert(tmp_path, "--compute-integrity")
    target = tmp_path / "sorted"

    result = CliRunner().invoke(
        qpx_main,
        ["transform", "relayout", "--dataset", str(source), "--output-folder", str(target), "--row-group-size", "50"],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    after_file = pq.ParquetFile(target / "sn.feature.parquet")
    after = after_file.read()
    anchors = after.column("anchor_protein").to_pylist()
    assert anchors == sorted(anchors)
    assert after.column("feature_id").to_pylist() == list(range(1, after.num_rows + 1))
    assert after.schema.metadata[b"identity_source"] == b"provided"
    with Dataset(str(source), file_prefix="sn") as unsorted, Dataset(str(target), file_prefix="sn") as relayouted:
        before, read = unsorted.feature.to_arrow(), relayouted.feature.to_arrow()
    values = [name for name in before.schema.names if name != "feature_id"]
    assert sorted(map(repr, read.select(values).to_pylist())) == sorted(map(repr, before.select(values).to_pylist()))
    assert read.schema.equals(before.schema)
    assert after.schema.metadata[b"qpx.sort_order"] == b"anchor_protein,peptidoform,charge,run_file_name"
    assert after_file.metadata.row_group(0).num_rows == 50
    sorting = after_file.metadata.row_group(0).sorting_columns
    keys = ["anchor_protein", "peptidoform", "charge", "run_file_name"]
    assert [after_file.schema.column(column.column_index).path for column in sorting] == keys
    # sorting_columns cannot name the list column grouped_runs, so it stops before it.
    pg_file = pq.ParquetFile(target / "sn.pg.parquet")
    assert pg_file.schema_arrow.metadata[b"qpx.sort_order"] == b"anchor_protein,grouped_runs,label"
    pg_sorting = pg_file.metadata.row_group(0).sorting_columns
    assert [pg_file.schema.column(column.column_index).path for column in pg_sorting] == ["anchor_protein"]
    with Dataset(str(target), file_prefix="sn") as dataset:
        assert not any(dataset.verify_integrity().values())


def test_relayout_refuses_a_column_it_cannot_sort_by(tmp_path):
    source = _convert(tmp_path)

    result = CliRunner().invoke(
        qpx_main,
        ["transform", "relayout", "--dataset", str(source), "--in-place", "--view", "feature", "--sort-by", "intensities"],
    )

    assert result.exit_code != 0
    assert "no sortable column named intensities" in result.output


def _linked_dataset(folder):
    """pg, feature and psm views out of sort order, linked by ids, with feature hardlinks.

    One PSM names a feature the view lacks and one names none.
    """
    folder.mkdir()
    anchors = {"CCCK": "P3", "AAAK": "P1", "BBBK": "P2"}
    with PgWriter(folder / "ds.pg.parquet") as writer:
        writer.write_batch([make_pg_record(anchor_protein=anchor, pg_accessions=[anchor]) for anchor in anchors.values()])
    pg = pq.read_table(folder / "ds.pg.parquet").to_pylist()
    pg_ids = {row["anchor_protein"]: row["pg_id"] for row in pg}

    def features(psm_ids):
        return [
            {
                **make_feature_record(sequence=peptide, peptidoform=peptide, anchor_protein=anchor),
                "pg_accessions": [{"accession": anchor, "start": None, "end": None}],
                "pg_ids": [pg_ids[anchor]],
                "psm_ids": psm_ids.get(peptide),
            }
            for peptide, anchor in anchors.items()
        ]

    with FeatureWriter(folder / "ds.feature.parquet") as writer:
        writer.write_batch(features({}))
    feature_ids = {row["peptidoform"]: row["feature_id"] for row in pq.read_table(folder / "ds.feature.parquet").to_pylist()}
    psms = [
        {**make_psm_record(sequence=peptide, peptidoform=peptide, scan=[scan]), "feature_id": feature_id}
        for scan, (peptide, feature_id) in enumerate([*feature_ids.items(), ("DDDK", 12345), ("EEEK", None)], start=1)
    ]
    with PsmWriter(folder / "ds.psm.parquet") as writer:
        writer.write_batch(psms)
    psm_ids = {row["peptidoform"]: [row["psm_id"]] for row in pq.read_table(folder / "ds.psm.parquet").to_pylist()}
    with FeatureWriter(folder / "ds.feature.parquet") as writer:
        writer.write_batch(features(psm_ids))
    return folder


def _relayout(source, target, *arguments):
    result = CliRunner().invoke(
        qpx_main,
        ["transform", "relayout", "--dataset", str(source), "--output-folder", str(target), *arguments],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output


def test_relayout_renumbers_ids_and_the_references_to_them(tmp_path):
    """Sorted views get ids 1, 2, 3, ... and every column holding them follows."""
    _relayout(_linked_dataset(tmp_path / "ds"), tmp_path / "sorted")

    pg = pq.read_table(tmp_path / "sorted" / "ds.pg.parquet")
    feature = pq.read_table(tmp_path / "sorted" / "ds.feature.parquet")
    psm = pq.read_table(tmp_path / "sorted" / "ds.psm.parquet")
    assert pg.select(["anchor_protein", "pg_id"]).to_pylist() == [
        {"anchor_protein": anchor, "pg_id": number} for number, anchor in enumerate(["P1", "P2", "P3"], start=1)
    ]
    assert feature.select(["peptidoform", "feature_id", "pg_ids", "psm_ids"]).to_pylist() == [
        {"peptidoform": peptide, "feature_id": number, "pg_ids": [number], "psm_ids": [number]}
        for number, peptide in enumerate(["AAAK", "BBBK", "CCCK"], start=1)
    ]
    # A reference to a feature the view lacks keeps its value.
    assert psm.select(["peptidoform", "psm_id", "feature_id"]).to_pylist() == [
        {"peptidoform": peptide, "psm_id": number, "feature_id": feature_id}
        for number, (peptide, feature_id) in enumerate(
            [("AAAK", 1), ("BBBK", 2), ("CCCK", 3), ("DDDK", 12345), ("EEEK", None)], start=1
        )
    ]
    for table in (pg, feature, psm):
        assert table.schema.metadata[b"identity_source"] == b"provided"
    # feature was written again for psm_ids after psm was sorted, and still states its order.
    feature_file = pq.ParquetFile(tmp_path / "sorted" / "ds.feature.parquet")
    sorting = feature_file.metadata.row_group(0).sorting_columns
    assert [feature_file.schema.column(column.column_index).path for column in sorting] == [
        "anchor_protein",
        "peptidoform",
        "charge",
        "run_file_name",
    ]


def test_relayout_of_one_view_rewrites_the_references_in_the_others(tmp_path):
    """psm keeps its order and ids when only feature is sorted, and follows the new feature ids."""
    source = _linked_dataset(tmp_path / "ds")
    before = pq.read_table(source / "ds.psm.parquet")

    _relayout(source, tmp_path / "sorted", "--view", "feature")

    psm = pq.read_table(tmp_path / "sorted" / "ds.psm.parquet")
    feature = pq.read_table(tmp_path / "sorted" / "ds.feature.parquet").to_pylist()
    new_ids = {row["peptidoform"]: row["feature_id"] for row in feature}
    assert psm.column("psm_id").equals(before.column("psm_id"))
    assert psm.column("feature_id").to_pylist() == [
        new_ids.get(peptide, old)
        for peptide, old in zip(before.column("peptidoform").to_pylist(), before.column("feature_id").to_pylist())
    ]
    assert [row["psm_ids"] for row in feature] == [
        [row["psm_id"]] for row in sorted(before.to_pylist()[:3], key=lambda row: row["peptidoform"])
    ]
