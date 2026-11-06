"""Flat q-values and group keys are the same whichever writer path and producer filled them."""

import duckdb
import pyarrow.parquet as pq
import pytest

import qpx
from qpx.core.data.identity import derive_id
from qpx.writers.derived import is_psm_qvalue
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter
from tests.conftest import make_feature_record as _base_feature
from tests.conftest import make_psm_record as _base_psm


def make_psm_record(**fields):
    return {**_base_psm(), **fields}


def make_feature_record(**fields):
    return {**_base_feature(), **fields}


def _score(name, value):
    return {"score_name": name, "score_value": value, "higher_better": False}


@pytest.mark.parametrize(
    "name, expected",
    [
        ("qvalue", True),
        ("q_value", True),
        ("psm_level_qvalue", True),
        ("percolator_q_value", True),
        ("ms_gf_qvalue", True),
        ("peptide_level_qvalue", False),
        ("global_qvalue", False),
        ("pg_qvalue", False),
        ("diann_lib_qvalue", False),
        ("andromeda_score", False),
    ],
)
def test_psm_qvalue_names(name, expected):
    assert is_psm_qvalue(name) is expected


def test_psm_writer_fills_qvalue_under_any_source_name(tmp_path):
    records = [
        make_psm_record(run_file_name=f"run{i}", scan=[i], additional_scores=scores)
        for i, scores in enumerate(
            [
                [_score("q-value", 0.02)],
                [_score("peptide_level_qvalue", 0.5), _score("percolator_q_value", 0.03)],
                [_score("andromeda_score", 90.0)],
                None,
            ]
        )
    ]
    path = tmp_path / "x.psm.parquet"
    with PsmWriter(path) as writer:
        writer.write_batch(records)
    values = pq.read_table(path, columns=["qvalue"]).column("qvalue").to_pylist()
    assert values == [0.02, 0.03, None, None]


def test_feature_writer_fills_each_level(tmp_path):
    record = make_feature_record(
        additional_scores=[_score("precursor_qvalue", 0.01), _score("global_qvalue", 0.02), _score("peptide_level_qvalue", 0.03)],
    )
    path = tmp_path / "x.feature.parquet"
    with FeatureWriter(path) as writer:
        writer.write_batch([record])
    (row,) = pq.read_table(path, columns=["precursor_qvalue", "global_qvalue", "peptide_qvalue"]).to_pylist()
    assert row == {"precursor_qvalue": 0.01, "global_qvalue": 0.02, "peptide_qvalue": 0.03}


def test_a_value_the_producer_typed_wins(tmp_path):
    record = make_feature_record(peptide_qvalue=0.001, additional_scores=[_score("peptide_level_qvalue", 0.03)])
    path = tmp_path / "x.feature.parquet"
    with FeatureWriter(path) as writer:
        writer.write_table(FeatureWriter(path.with_suffix(".unused")).align_table_to_schema(_table([record])))
    assert pq.read_table(path, columns=["peptide_qvalue"]).column("peptide_qvalue").to_pylist() == [0.001]


def _table(records):
    import pyarrow as pa

    from qpx.core.data import FeatureSchema

    relaxed = pa.schema([field.with_nullable(True) for field in FeatureSchema.get_arrow_schema()])
    return pa.Table.from_pylist(records, schema=relaxed)


def test_group_id_is_the_hash_of_the_sorted_membership(tmp_path):
    feature_path, pg_path = tmp_path / "x.feature.parquet", tmp_path / "x.pg.parquet"
    members = [{"accession": acc, "start": None, "end": None, "pre": None, "post": None} for acc in ("P2", "P1")]
    with FeatureWriter(feature_path) as writer:
        writer.write_batch([make_feature_record(pg_accessions=members, anchor_protein="P2")])
    with PgWriter(pg_path) as writer:
        writer.write_batch(
            [
                {
                    "pg_accessions": ["P1", "P2"],
                    "anchor_protein": "P1",
                    "grouped_runs": ["run"],
                    "intensities": [{"label": "LFQ", "intensity": 1.0}],
                    "peptides": [{"protein_name": "P1", "peptide_count": 1}],
                    "is_decoy": False,
                    "contaminant": False,
                }
            ]
        )
    expected = derive_id([["P1", "P2"]])
    assert pq.read_table(feature_path, columns=["pg_group_id"]).column("pg_group_id").to_pylist() == [expected]
    assert pq.read_table(pg_path, columns=["pg_group_id"]).column("pg_group_id").to_pylist() == [expected]


def _link_dataset(path):
    lfq = [{"label": "LFQ", "intensity": 1.0}]
    members = [{"accession": acc, "start": None, "end": None, "pre": None, "post": None} for acc in ("P2", "P1")]
    with FeatureWriter(path / "x.feature.parquet") as writer:
        writer.write_batch(
            [
                make_feature_record(
                    pg_accessions=members,
                    anchor_protein="P1",
                    run_file_name="run",
                    peptidoform="AK",
                    sequence="AK",
                    intensities=lfq,
                ),
                make_feature_record(
                    pg_accessions=members[:1],
                    anchor_protein="P2",
                    run_file_name="run",
                    peptidoform="CK",
                    sequence="CK",
                    intensities=lfq,
                ),
            ]
        )
    with PgWriter(path / "x.pg.parquet") as writer:
        writer.write_batch(
            [
                {
                    "pg_accessions": accessions,
                    "anchor_protein": accessions[0],
                    "grouped_runs": ["run"],
                    "intensities": [{"label": "LFQ", "intensity": 1.0}],
                    "peptides": [{"protein_name": accessions[0], "peptide_count": 1}],
                    "is_decoy": False,
                    "contaminant": False,
                }
                for accessions in (["P1", "P2"], ["P2"])
            ]
        )


def test_link_feature_pg_of_one_group_matches_the_full_link(tmp_path):
    _link_dataset(tmp_path)
    with qpx.open_dataset(str(tmp_path)) as ds:
        full = ds.link_feature_pg().to_df()
        one = ds.link_feature_pg(group=["P2", "P1"]).to_df()
        pg_id = ds.sql("SELECT pg_id FROM pg WHERE len(pg_accessions) = 2").fetchone()[0]
    assert len(full) == 2
    assert one["pg_id"].tolist() == [pg_id]
    assert one["feature_id"].tolist() == full[full["pg_id"] == pg_id]["feature_id"].tolist()


def test_confident_filters_on_the_typed_column(tmp_path):
    path = tmp_path / "x.psm.parquet"
    records = [
        make_psm_record(run_file_name=f"run{i}", scan=[i], additional_scores=[_score("q_value", value)])
        for i, value in enumerate((0.001, 0.5))
    ]
    with PsmWriter(path) as writer:
        writer.write_batch(records)
    from qpx.core.data.psm import PSM

    assert PSM.from_file(path).confident().count() == 1
    assert duckdb.sql(f"SELECT COUNT(*) FROM '{path}' WHERE qvalue <= 0.01").fetchall() == [(1,)]
