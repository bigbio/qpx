"""Physical merge, write-back and partitioned writes keep a dataset's rows, ids and footers."""

import logging
import shutil

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pytest

import qpx
from qpx.writers import FeatureWriter, PgWriter, PsmWriter
from tests.conftest import make_feature_record, make_psm_record


def _copy(dataset_dir, tmp_path, name, accession):
    """A copy of the fixture dataset with its own accession, under its own directory."""
    target = tmp_path / name
    shutil.copytree(dataset_dir, target)
    path = next(target.glob("*.dataset.parquet"))
    table = pq.read_table(path)
    table = table.set_column(0, "project_accession", pa.array([accession] * table.num_rows, type=pa.string()))
    pq.write_table(table, path)
    return target


@pytest.fixture
def two_copies(dataset_dir, tmp_path):
    first = qpx.open_dataset(str(_copy(dataset_dir, tmp_path, "one", "PXD000001")))
    second = qpx.open_dataset(str(_copy(dataset_dir, tmp_path, "two", "PXD000002")))
    collection = qpx.DatasetCollection([first, second])
    yield first, second, collection
    collection.close()
    first.close()
    second.close()


class TestMergeCollisions:
    def test_default_refuses_overlap_and_writes_nothing(self, two_copies, tmp_path):
        _, _, collection = two_copies
        output = tmp_path / "merged"

        with pytest.raises(ValueError, match=r"overlap.*run\(s\).*feature_id.*on_collision='namespace'"):
            collection.merge(output)

        assert not output.exists()

    def test_namespace_keeps_both_with_links_intact(self, two_copies, tmp_path):
        first, second, collection = two_copies
        output = tmp_path / "merged"

        collection.merge(output, on_collision="namespace")

        with qpx.open_dataset(str(output)) as merged:
            feature = merged.feature.to_df()
            assert len(feature) == first.feature.count() + second.feature.count()
            assert feature["feature_id"].is_unique
            assert set(feature["source_dataset"]) == {"PXD000001", "PXD000002"}
            assert all(run.startswith(("PXD000001:", "PXD000002:")) for run in feature["run_file_name"])
            runs = merged.run.to_df()
            assert set(feature["run_file_name"]) <= set(runs["run_file_name"])
            run_samples = {s["sample_accession"] for samples in runs["samples"] for s in samples}
            assert run_samples <= set(merged.sample.to_df()["sample_accession"])
            pg_runs = {run for runs_ in merged.pg.to_df()["grouped_runs"] for run in runs_}
            assert pg_runs <= set(runs["run_file_name"])
            results = merged.validate()
            assert not {name: result.errors for name, result in results.items() if result.errors}

    def test_namespace_does_not_inflate_quantities(self, two_copies, tmp_path):
        first, _, collection = two_copies
        output = tmp_path / "merged"

        collection.merge(output, on_collision="namespace")

        with qpx.open_dataset(str(output)) as merged:
            matrix = merged.pg.to_df()
            assert len(matrix) == 2 * first.pg.count()
            assert merged.dataset_meta.count() == 1

    def test_dedupe_identical_drops_copies(self, two_copies, tmp_path):
        first, _, collection = two_copies
        output = tmp_path / "merged"

        collection.merge(output, on_collision="dedupe-identical")

        with qpx.open_dataset(str(output)) as merged:
            assert merged.feature.count() == first.feature.count()
            assert merged.pg.count() == first.pg.count()
            assert merged.sample.count() == first.sample.count()
            assert merged.run.count() == first.run.count()

    def test_dedupe_identical_raises_on_a_difference(self, dataset_dir, tmp_path):
        other = _copy(dataset_dir, tmp_path, "other", "PXD000002")
        path = next(other.glob("*.feature.parquet"))
        table = pq.read_table(path)
        changed = [value + 1 for value in table.column("charge").to_pylist()]
        pq.write_table(
            table.set_column(
                table.schema.get_field_index("charge"),
                table.schema.field("charge"),
                pa.array(changed, type=table.schema.field("charge").type),
            ),
            path,
        )
        first, second = qpx.open_dataset(str(dataset_dir)), qpx.open_dataset(str(other))
        with qpx.DatasetCollection([first, second]) as collection:
            with pytest.raises(ValueError, match="not identical"):
                collection.merge(tmp_path / "merged", on_collision="dedupe-identical")

    def test_dedupe_identical_refuses_two_analyses_of_one_run(self, dataset_dir, tmp_path):
        other = _copy(dataset_dir, tmp_path, "other", "PXD000002")
        path = next(other.glob("*.feature.parquet"))
        table = pq.read_table(path)
        ids = [value + 1 for value in table.column("feature_id").to_pylist()]
        pq.write_table(
            table.set_column(
                table.schema.get_field_index("feature_id"), table.schema.field("feature_id"), pa.array(ids, type=pa.int64())
            ),
            path,
        )
        first, second = qpx.open_dataset(str(dataset_dir)), qpx.open_dataset(str(other))
        with qpx.DatasetCollection([first, second]) as collection:
            with pytest.raises(ValueError, match="in a run an earlier dataset also has"):
                collection.merge(tmp_path / "merged", on_collision="dedupe-identical")

    def test_unknown_mode_and_structure_raise(self, two_copies, tmp_path):
        _, _, collection = two_copies
        with pytest.raises(ValueError, match="on_collision must be"):
            collection.merge(tmp_path / "merged", on_collision="overwrite")
        with pytest.raises(ValueError, match=r"Unknown structure\(s\) \['featur'\]"):
            collection.merge(tmp_path / "merged", structures=["featur"], on_collision="namespace")


class TestMergeOutput:
    def test_metadata_views_follow_their_rules(self, two_copies, tmp_path):
        first, _, collection = two_copies
        output = tmp_path / "merged"

        collection.merge(output, on_collision="namespace")

        with qpx.open_dataset(str(output)) as merged:
            dataset = merged.dataset_meta.to_df().iloc[0]
            assert dataset["project_accession"] == "PXD000001;PXD000002"

    def test_dropped_views_are_named(self, dataset_dir, tmp_path, caplog):
        second = _copy(dataset_dir, tmp_path, "two", "PXD000002")
        next(second.glob("*.pg.parquet")).unlink()
        first, other = qpx.open_dataset(str(dataset_dir)), qpx.open_dataset(str(second))
        with qpx.DatasetCollection([first, other]) as collection, caplog.at_level(logging.WARNING):
            collection.merge(tmp_path / "merged", on_collision="namespace")
        assert "['pg']" in caplog.text

    def test_footers_and_row_groups_come_from_the_writers(self, two_copies, tmp_path):
        _, _, collection = two_copies
        output = tmp_path / "merged"

        collection.merge(output, on_collision="namespace")

        metadata = pq.read_schema(output / "merged.feature.parquet").metadata
        assert metadata[b"file_type"] == b"feature_file"
        assert b"qpx_version" in metadata and b"identity_composite" in metadata

    def test_different_composites_write_one_file_each(self, tmp_path, dataset_dir):
        one, two = tmp_path / "one", tmp_path / "two"
        one.mkdir()
        two.mkdir()
        record = make_feature_record(sequence="PEPTIDEK", anchor_protein="P12345", run_file_name="run_a")
        with FeatureWriter(one / "a.feature.parquet") as writer:
            writer.write_batch([dict(record)])
        with FeatureWriter(two / "b.feature.parquet", identity_composite=("peptidoform", "charge", "run_file_name")) as writer:
            writer.write_batch([{**record, "run_file_name": "run_b"}])
        first, second = qpx.open_dataset(str(one)), qpx.open_dataset(str(two))
        with qpx.DatasetCollection([first, second]) as collection:
            collection.merge(tmp_path / "merged")
        files = sorted(p.name for p in (tmp_path / "merged").iterdir())
        assert files == ["merged-1.feature.parquet", "merged-2.feature.parquet"]
        composites = {pq.read_schema(tmp_path / "merged" / name).metadata[b"identity_composite"] for name in files}
        assert len(composites) == 2
        with qpx.open_dataset(str(tmp_path / "merged")) as merged:
            assert merged.feature.count() == 2


class TestOpening:
    def test_missing_directory_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="not found"):
            qpx.open_dataset(str(tmp_path / "typo"))

    @pytest.mark.parametrize("url", ["https://example.org/PXD000001", "http://example.org/x", "ftp://example.org/x"])
    def test_unsupported_scheme_raises(self, url):
        with pytest.raises(ValueError, match="only local directories and s3://"):
            qpx.open_dataset(url)

    def test_shards_with_different_columns_read_by_name(self, dataset_dir, tmp_path, caplog):
        target = tmp_path / "ds"
        shutil.copytree(dataset_dir, target)
        original = next(target.glob("*.feature.parquet"))
        table = pq.read_table(original)
        pq.write_table(table.slice(0, 1), original)
        pq.write_table(table.slice(1).drop_columns(["observed_mz"]), target / "extra.feature.parquet")

        with caplog.at_level(logging.WARNING), qpx.open_dataset(str(target)) as ds:
            frame = ds.feature.to_df()

        assert len(frame) == table.num_rows
        assert frame["observed_mz"].isna().sum() == table.num_rows - 1
        assert "read from 2 files as one view" in caplog.text

    def test_validate_reports_an_empty_directory(self, tmp_path):
        with qpx.open_dataset(str(tmp_path)) as ds:
            assert ds.validate()["dataset"].errors[0].check == "no_structures"

    def test_validate_refuses_a_two_row_dataset_view(self, dataset_dir):
        path = next(dataset_dir.glob("*.dataset.parquet"))
        table = pq.read_table(path)
        pq.write_table(pa.concat_tables([table, table]), path)
        with qpx.open_dataset(str(dataset_dir)) as ds:
            assert any(issue.check == "single_row" for issue in ds.validate()["dataset"].errors)


class TestSaveStructure:
    def _producer_dataset(self, tmp_path):
        """A dataset whose feature composite is not the schema's, with a psm linked to a feature."""
        root = tmp_path / "ds"
        root.mkdir()
        composite = ("peptidoform", "charge", "run_file_name")
        record = make_feature_record(sequence="PEPTIDEK", anchor_protein="P12345", run_file_name="run_a")
        with FeatureWriter(root / "study.feature.parquet", identity_composite=composite, creator="tool") as writer:
            writer.write_batch([record])
        return root, composite

    def test_round_trip_keeps_rows_ids_composite_and_prefix(self, tmp_path):
        root, composite = self._producer_dataset(tmp_path)
        with qpx.open_dataset(str(root)) as ds:
            before = ds.feature.to_arrow()
            path = ds.save_structure(before, "feature")
        assert path.name == "study.feature.parquet"
        assert sorted(p.name for p in root.iterdir()) == ["study.feature.parquet"]
        metadata = pq.read_schema(path).metadata
        assert metadata[b"identity_composite"] == ",".join(composite).encode()
        assert metadata[b"creator"] == b"tool"
        with qpx.open_dataset(str(root)) as ds:
            after = ds.feature.to_arrow()
        assert after["feature_id"].to_pylist() == before["feature_id"].to_pylist()
        assert after["cv_params"].to_pylist() == before["cv_params"].to_pylist()

    def test_pg_round_trip_with_null_nested_structs(self, dataset_dir):
        with qpx.open_dataset(str(dataset_dir)) as ds:
            table = ds.pg.to_arrow()
            # DuckDB returns a null struct with null children; a required child refused it.
            counts = pa.array([None] * table.num_rows, type=table.schema.field("feature_counts").type)
            table = table.set_column(table.schema.get_field_index("feature_counts"), "feature_counts", counts)
            path = ds.save_structure(table, "pg")
        assert pq.read_table(path).num_rows == table.num_rows

    def test_a_typo_key_raises(self, tmp_path):
        with PgWriter(tmp_path / "x.pg.parquet") as writer:
            with pytest.raises(ValueError, match=r"keys outside the schema: \['anchr_protein'\]"):
                writer.write_batch([{"anchr_protein": "P1"}] * writer._batch_size)

    def test_an_unwritable_value_names_column_and_row(self, tmp_path):
        record = make_feature_record(sequence="PEPTIDEK", anchor_protein="P12345", run_file_name="run_a")
        with pytest.raises(ValueError, match=r"column 'charge', row 1 of the batch, value 'two'"):
            with FeatureWriter(tmp_path / "x.feature.parquet") as writer:
                writer.write_batch([record, {**record, "charge": "two"}])

    def test_derive_ids_gives_the_id_the_writer_keeps(self, tmp_path):
        record = make_feature_record(sequence="PEPTIDEK", anchor_protein="P12345", run_file_name="run_a")
        composite = ("peptidoform", "charge", "run_file_name")
        schema = FeatureWriter(tmp_path / "unused.parquet", override_provided_ids=False)._relaxed_id_schema()
        table = pa.Table.from_pylist([{**record, "feature_id": None}], schema=schema)
        derived = FeatureWriter.derive_ids(table, identity_composite=composite)["feature_id"].to_pylist()
        with FeatureWriter(tmp_path / "a.feature.parquet", identity_composite=composite) as writer:
            writer.write_batch([record])
        assert pq.read_table(tmp_path / "a.feature.parquet")["feature_id"].to_pylist() == derived


def test_scan_format_violation_is_logged_at_close(tmp_path, caplog):
    record = make_psm_record(scan=[1, 2, 3])
    with caplog.at_level(logging.WARNING), PsmWriter(tmp_path / "x.psm.parquet", scan_format="scan") as writer:
        writer.write_batch([record])
    assert "violate scan_format='scan'" in caplog.text


def test_quantify_reads_a_partitioned_feature_directory(tmp_path):
    from qpx.cli.transform.quantify import _qpx_feature_to_peptide_df

    record = make_feature_record(sequence="PEPTIDEK", anchor_protein="P12345", run_file_name="0042")
    with FeatureWriter(tmp_path / "a.feature.parquet") as writer:
        writer.write_batch([record])
    FeatureWriter.write_partitioned(pq.read_table(tmp_path / "a.feature.parquet"), tmp_path / "feature")

    peptides = _qpx_feature_to_peptide_df(tmp_path / "feature")

    assert set(peptides["SampleID"].astype(str)) == set(
        _qpx_feature_to_peptide_df(tmp_path / "a.feature.parquet")["SampleID"].astype(str)
    )
    assert any("0042" in str(sample) for sample in peptides["SampleID"])


class TestWritePartitioned:
    def _table(self, tmp_path, runs=("run_a", "run_b")):
        composite = ("peptidoform", "charge", "run_file_name")
        path = tmp_path / "source.feature.parquet"
        records = [
            make_feature_record(sequence="PEPTIDEK", anchor_protein="P12345", run_file_name=run, charge=charge)
            for run in runs
            for charge in (2, 3)
        ]
        with FeatureWriter(path, identity_composite=composite) as writer:
            writer.write_batch(records)
        return pq.read_table(path)

    def test_keeps_ids_and_declares_the_source_composite(self, tmp_path):
        table = self._table(tmp_path)
        output = tmp_path / "feature"

        FeatureWriter.write_partitioned(table, output)

        written = sorted(pq.read_table(next(output.glob("run_file_name=run_a/*.parquet")))["feature_id"].to_pylist())
        assert set(written) <= set(table["feature_id"].to_pylist())
        metadata = pq.read_schema(next(output.rglob("*.parquet"))).metadata
        assert metadata[b"identity_composite"] == b"peptidoform,charge,run_file_name"
        assert metadata[b"file_type"] == b"feature_file"

    def test_second_write_refused_by_default(self, tmp_path):
        table = self._table(tmp_path)
        output = tmp_path / "feature"
        FeatureWriter.write_partitioned(table, output)
        with pytest.raises(FileExistsError, match="mode='append'"):
            FeatureWriter.write_partitioned(table, output)

    def test_append_keeps_earlier_rows_in_shared_partitions(self, tmp_path):
        table = self._table(tmp_path)
        output = tmp_path / "feature"

        FeatureWriter.write_partitioned(table.slice(0, 3), output)
        FeatureWriter.write_partitioned(table.slice(3), output, mode="append")

        rows = sum(pq.read_metadata(path).num_rows for path in output.rglob("*.parquet"))
        assert rows == table.num_rows

    def test_replace_partitions_keeps_the_others(self, tmp_path):
        table = self._table(tmp_path)
        output = tmp_path / "feature"
        FeatureWriter.write_partitioned(table, output)
        run_a = table.filter(pc.equal(table["run_file_name"], "run_a"))

        FeatureWriter.write_partitioned(run_a.slice(0, 1), output, mode="replace_partitions")

        assert pq.read_metadata(next(output.glob("run_file_name=run_a/*.parquet"))).num_rows == 1
        assert len(list(output.glob("run_file_name=run_a/*.parquet"))) == 1
        assert pq.read_metadata(next(output.glob("run_file_name=run_b/*.parquet"))).num_rows == 2

    @pytest.mark.parametrize("run", ["2024-01-15", "0042"])
    def test_run_names_stay_strings_when_read(self, tmp_path, run):
        table = self._table(tmp_path, runs=(run,))
        root = tmp_path / "ds"
        root.mkdir()
        FeatureWriter.write_partitioned(table, root / "feature")
        with qpx.open_dataset(str(root)) as ds:
            assert set(ds.feature.to_df()["run_file_name"]) == {run}


def test_merge_warns_when_pg_quantities_differ(dataset_dir, tmp_path, caplog):
    def with_method(name, accession, method):
        target = _copy(dataset_dir, tmp_path, name, accession)
        path = next(target.glob("*.pg.parquet"))
        table = pq.read_table(path)
        index = table.schema.get_field_index("cv_params")
        field = table.schema.field(index)
        entry = [{"cv_name": "quantification_method", "cv_value": method}]
        pq.write_table(table.set_column(index, field, pa.array([entry] * table.num_rows, field.type)), path)
        return qpx.open_dataset(str(target))

    with (
        with_method("q1", "PXD000001", "PG.Quantity") as first,
        with_method("q2", "PXD000002", "PG.MaxLFQ") as second,
        qpx.DatasetCollection([first, second]) as coll,
        caplog.at_level(logging.WARNING, logger="qpx.core.dataset_merge"),
    ):
        coll.merge(tmp_path / "merged", structures=["pg"], on_collision="namespace")
    assert "PG.Quantity" in caplog.text
    assert "PG.MaxLFQ" in caplog.text
