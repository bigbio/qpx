"""Native OpenMS run references share SDRF names without losing identities."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.converters.openms.converter import OpenMSConverter
from qpx.converters.openms.run_names import normalize_run_names
from qpx.core.data.identity import derive_id
from qpx.core.data.loader import load_schema
from qpx.dataset import Dataset
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter
from tests.conftest import make_feature_record, make_pg_record, make_psm_record, write_lfq_consensusxml


@pytest.mark.parametrize("view", ["feature", "psm", "pg"])
def test_scalar_run_names_follow_sdrf_convention(view):
    """Paths and known extensions disappear; dots, nulls and unknowns survive."""
    values = [r"C:\data\sample.part.1.mzML.gz", "/data/sample.part.2.RAW", "sample.part.3", "unknown", "", None]
    expected = ["sample.part.1", "sample.part.2", "sample.part.3", "unknown", "", None]
    table = pa.table({"run_file_name": values, "id_run_file_name": values, "feature_id": range(len(values))})

    result = normalize_run_names(table, view)

    assert result["run_file_name"].to_pylist() == expected
    assert result["id_run_file_name"].to_pylist() == (expected if view == "feature" else values)
    assert result["feature_id"].equals(table["feature_id"])


def test_grouped_runs_preserve_order_missing_values_and_metadata():
    """Normalization never sorts, deduplicates, or guesses a group's members."""
    field = pa.field("grouped_runs", pa.large_list(pa.large_string()), metadata={b"source": b"OpenMS"})
    schema = pa.schema([field], metadata={b"identity_composite": b"grouped_runs"})
    values = [["/data/run.2.RAW", r"C:\data\run.1.mzML.gz", "run.2"], [], None, [None, "unknown", ""]]
    table = pa.Table.from_arrays([pa.array(values, type=field.type)], schema=schema)

    result = normalize_run_names(table, "pg")

    assert result["grouped_runs"].to_pylist() == [["run.2", "run.1", "run.2"], [], None, [None, "unknown", ""]]
    assert result.schema.equals(table.schema, check_metadata=True)


def _write_native_bundle(folder, legacy):
    """Write native run spellings, including linked producer IDs when supplied."""
    feature = make_feature_record(
        run_file_name=r"C:\acquisition\sample.part.1.mzML.gz",
        intensities=[{"label": "LFQ", "intensity": 1000.0}],
    )
    feature.update(id_run_file_name="/identification/sample.part.2.mzML", feature_id=11, psm_ids=[22], pg_ids=[33])
    psm = make_psm_record(run_file_name="/acquisition/sample.part.1.mzML")
    psm.update(psm_id=22, feature_id=11)
    protein = make_pg_record(intensities=[{"label": "LFQ", "intensity": 5000.0}])
    protein.update(pg_id=33, grouped_runs=["/acquisition/sample.part.1.raw"])
    with FeatureWriter(folder / "native.feature.parquet", override_provided_ids=False) as writer:
        writer.write_batch([feature])
    with PsmWriter(folder / "native.psm.parquet") as writer:
        writer.write_batch([psm])
    with PgWriter(folder / "native.pg.parquet") as writer:
        writer.write_batch([protein])
    if legacy:
        for view, columns in {
            "feature": ["feature_id", "psm_ids", "pg_ids"],
            "psm": ["psm_id", "feature_id"],
        }.items():
            path = folder / f"native.{view}.parquet"
            table = pq.read_table(path).drop_columns(columns).replace_schema_metadata(None)
            pq.write_table(table, path)
        path = folder / "native.pg.parquet"
        table = pq.read_table(path).drop_columns(["pg_id", "grouped_runs"]).replace_schema_metadata(None)
        table = table.append_column("run_file_name", pa.array(["/acquisition/sample.part.1.mzML"]))
        pq.write_table(table, path)


@pytest.mark.parametrize("legacy", [False, True])
def test_native_enrichment_joins_sdrf_and_preserves_supplied_links(tmp_path, legacy):
    """Current and legacy native outputs join samples using normalized run names."""
    source = tmp_path / "source"
    source.mkdir()
    _write_native_bundle(source, legacy)
    sdrf = tmp_path / "input.sdrf.tsv"
    sdrf.write_text(
        "source name\tcharacteristics[organism]\tcharacteristics[organism part]\tcomment[data file]\tcomment[label]\n"
        "sample_1\tHomo sapiens\tliver\t/data/sample.part.1.raw\tlabel free sample\n"
        "sample_2\tHomo sapiens\tliver\t/data/sample.part.2.raw\tlabel free sample\n"
    )
    output = tmp_path / "output"

    OpenMSConverter(source, sdrf).convert(output)

    tables = {view: pq.read_table(output / f"openms.{view}.parquet") for view in ("feature", "psm", "pg")}
    for view, table in tables.items():
        assert load_schema(view).validate_full(table, strict=True).is_valid
        if legacy:
            row = table.to_pylist()[0]
            composite = table.schema.metadata[b"identity_composite"].decode().split(",")
            unordered = tuple(index for index, name in enumerate(composite) if name in ("grouped_runs", "pg_accessions"))
            assert row[f"{view}_id"] == derive_id([row[name] for name in composite], unordered_list_indices=unordered)
    feature = tables["feature"].to_pylist()[0]
    psm = tables["psm"].to_pylist()[0]
    protein = tables["pg"].to_pylist()[0]
    assert feature["run_file_name"] == psm["run_file_name"] == "sample.part.1"
    assert feature["id_run_file_name"] == "sample.part.2"
    assert protein["grouped_runs"] == ["sample.part.1"]
    if not legacy:
        assert feature["feature_id"] == psm["feature_id"] == 11
        assert feature["psm_ids"] == [psm["psm_id"]] == [22]
        assert feature["pg_ids"] == [protein["pg_id"]] == [33]
    with Dataset(output, file_prefix="openms", duckdb_threads=6) as dataset:
        peptide = dataset.intensity("peptide").to_df().set_index("sample_accession")["intensity"].to_dict()
        protein_quantities = dataset.intensity("protein").to_df().set_index("sample_accession")["intensity"].to_dict()
    assert peptide == {"sample_1": 1000.0}
    assert protein_quantities == {"sample_1": 5000.0}


def test_compressed_runs_retain_fraction_group_annotation(tmp_path):
    """Normalize compressed names after resolving annotations from source names."""
    source = tmp_path / "source"
    source.mkdir()
    run_name = "/acquisition/sample.part.1.mzML.gz"
    with FeatureWriter(source / "native.feature.parquet") as writer:
        writer.write_batch([make_feature_record(run_file_name=run_name)])
    with PgWriter(source / "native.pg.parquet") as writer:
        writer.write_batch([make_pg_record(run_file_name=run_name)])
    consensusxml = tmp_path / "native.consensusXML"
    write_lfq_consensusxml(consensusxml, [(0, run_name, "3", "1", "1")])
    output = tmp_path / "output"

    OpenMSConverter(source, consensusxml_path=consensusxml).convert(output)

    for view in ("feature", "pg"):
        record = pq.read_table(output / f"openms.{view}.parquet").to_pylist()[0]
        assert {"cv_name": "fraction_group", "cv_value": "3"} in record["cv_params"]
        run_column = "grouped_runs" if view == "pg" else "run_file_name"
        expected_run = ["sample.part.1"] if view == "pg" else "sample.part.1"
        assert record[run_column] == expected_run


@pytest.mark.parametrize("declaration", ["sdrf", "maplist"])
def test_declared_canonical_run_is_not_stripped_again(tmp_path, declaration):
    """An acquisition suffix can be part of an already declared canonical stem."""
    source = tmp_path / "source"
    source.mkdir()
    feature = make_feature_record(run_file_name="sample.raw", intensities=[{"label": "LFQ", "intensity": 1000.0}])
    feature["id_run_file_name"] = "sample.raw"
    with FeatureWriter(source / "native.feature.parquet") as writer:
        writer.write_batch([feature])
    with PsmWriter(source / "native.psm.parquet") as writer:
        writer.write_batch([make_psm_record(run_file_name="sample.raw")])
    with PgWriter(source / "native.pg.parquet") as writer:
        writer.write_batch([make_pg_record(run_file_name="sample.raw", intensities=[{"label": "LFQ", "intensity": 5000.0}])])
    sdrf, companion = None, None
    if declaration == "sdrf":
        sdrf = tmp_path / "input.sdrf.tsv"
        sdrf.write_text(
            "source name\tcharacteristics[organism]\tcharacteristics[organism part]\tcomment[data file]\tcomment[label]\n"
            "sample_1\tHomo sapiens\tliver\tsample.raw.mzML\tlabel free sample\n"
        )
    else:
        companion = tmp_path / "native.consensusXML"
        write_lfq_consensusxml(companion, [(0, "sample.raw.mzML", "1", "1", "1")], trailing="</consensusXML>")
    output = tmp_path / "output"

    OpenMSConverter(source, sdrf, consensusxml_path=companion).convert(output)

    for view, column in (("feature", "run_file_name"), ("feature", "id_run_file_name"), ("psm", "run_file_name")):
        assert pq.read_table(output / f"openms.{view}.parquet")[column].to_pylist() == ["sample.raw"]
    assert pq.read_table(output / "openms.pg.parquet")["grouped_runs"].to_pylist() == [["sample.raw"]]
    with Dataset(output, file_prefix="openms", duckdb_threads=6) as dataset:
        assert all(result.is_valid for result in dataset.validate(strict=True).values())
        if declaration == "sdrf":
            for level, expected in (("peptide", 1000.0), ("protein", 5000.0)):
                quantities = dataset.intensity(level).to_df().set_index("sample_accession")["intensity"].to_dict()
                assert quantities == {"sample_1": expected}
