"""Accession normalization rewrites accession values and nothing else."""

from __future__ import annotations

import pyarrow as pa
import pyarrow.parquet as pq

from qpx.transforms.accession_normalizer import normalize_dataset_accessions, normalize_parquet

_MEMBER = pa.struct([pa.field("accession", pa.string(), nullable=False), pa.field("start", pa.int32())])


def test_forward_maps_anchor_and_members_of_both_layouts(tmp_path):
    """Anchors and pg_accessions members are mapped; nulls, blanks and other columns are kept."""
    feature = tmp_path / "x.feature.parquet"
    pq.write_table(
        pa.table(
            {
                "anchor_protein": ["sp|P1|A_HUMAN", None, "", "CONTAM_sp|CONTAM_P2|B_HUMAN"],
                "pg_accessions": pa.array(
                    [[{"accession": "sp|P1|A_HUMAN", "start": 3}, None], None, [], [{"accession": "P3", "start": None}]],
                    pa.list_(_MEMBER),
                ),
                "charge": [2, 3, 2, 1],
            }
        ),
        feature,
    )
    pg = tmp_path / "x.pg.parquet"
    pq.write_table(pa.table({"pg_accessions": [["tr|Q9|C_HUMAN", None, "P3"], None]}), pg)

    result = normalize_parquet(feature, tmp_path / "out.feature.parquet", direction="forward")
    normalize_parquet(pg, tmp_path / "out.pg.parquet", direction="forward")

    table = pq.read_table(tmp_path / "out.feature.parquet")
    assert result["accessions_changed"] == 2
    assert table.column("anchor_protein").to_pylist() == ["P1", None, "", "CONTAM_P2"]
    assert table.column("pg_accessions").to_pylist() == [
        [{"accession": "P1", "start": 3}, None],
        None,
        [],
        [{"accession": "P3", "start": None}],
    ]
    assert table.column("charge").to_pylist() == [2, 3, 2, 1]
    assert pq.read_table(tmp_path / "out.pg.parquet").column("pg_accessions").to_pylist() == [["Q9", None, "P3"], None]


def test_file_without_changes_is_copied(tmp_path):
    """A view with nothing to normalize is copied byte for byte."""
    source = tmp_path / "x.pg.parquet"
    pq.write_table(pa.table({"anchor_protein": ["P1"], "pg_accessions": [["P1", "P2"]]}), source)

    result = normalize_parquet(source, tmp_path / "out.pg.parquet", direction="forward")

    assert result["accessions_changed"] == 0
    assert (tmp_path / "out.pg.parquet").read_bytes() == source.read_bytes()


def test_dataset_normalization_covers_psm_lists_and_partitioned_views(tmp_path):
    """Every view is stripped in place, PSM protein lists and partitioned parts included."""
    pq.write_table(
        pa.table(
            {
                "anchor_protein": ["sp|CONTAM_P05787|CONTAM_K2C8_HUMAN"],
                "pg_accessions": [["sp|CONTAM_P05787|CONTAM_K2C8_HUMAN", "Q15149"]],
                "pg_names": [["CONTAM_K2C8_HUMAN", "PLEC_HUMAN"]],
            }
        ),
        tmp_path / "openms.pg.parquet",
    )
    part = tmp_path / "psm" / "run_file_name=a" / "part-0.parquet"
    part.parent.mkdir(parents=True)
    pq.write_table(pa.table({"protein_accessions": [["sp|Q15149|PLEC_HUMAN", "tr|A0A0|X_HUMAN"], None]}), part)

    assert normalize_dataset_accessions(tmp_path, "openms") == 2

    pg = pq.read_table(tmp_path / "openms.pg.parquet")
    assert pg.column("anchor_protein").to_pylist() == ["CONTAM_P05787"]
    assert pg.column("pg_accessions").to_pylist() == [["CONTAM_P05787", "Q15149"]]
    assert pg.column("pg_names").to_pylist() == [["CONTAM_K2C8_HUMAN", "PLEC_HUMAN"]]
    assert pq.read_table(part).column("protein_accessions").to_pylist() == [["Q15149", "A0A0"], None]
