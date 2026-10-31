"""Accession normalization rewrites accession values and nothing else."""

from __future__ import annotations

import pyarrow as pa
import pyarrow.parquet as pq

from qpx.transforms.accession_normalizer import normalize_parquet

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


def test_psm_and_pepmap_accessions_are_normalized_and_a_pg_file_keeps_its_ids(tmp_path):
    """PSM-to-pg matching needs one spelling everywhere; a pg_id stays but no longer re-derives."""
    from qpx.writers import PgWriter, PsmWriter
    from tests.conftest import make_pg_record, make_psm_record

    psm = tmp_path / "x.psm.parquet"
    with PsmWriter(psm) as writer:
        writer.write_batch([{**make_psm_record(), "protein_accessions": ["sp|P04637|P53_HUMAN"]}])
    normalize_parquet(psm, tmp_path / "out.psm.parquet", direction="forward")
    normalized = pq.read_table(tmp_path / "out.psm.parquet")
    assert normalized.column("protein_accessions").to_pylist() == [["P04637"]]
    assert pq.read_metadata(tmp_path / "out.psm.parquet").metadata[b"identity_source"] == b"derived"

    pg = tmp_path / "x.pg.parquet"
    with PgWriter(pg) as writer:
        writer.write_batch([make_pg_record(anchor_protein="sp|P04637|P53_HUMAN", pg_accessions=["sp|P04637|P53_HUMAN"])])
    normalize_parquet(pg, tmp_path / "out.pg.parquet", direction="forward")
    assert pq.read_metadata(tmp_path / "out.pg.parquet").metadata[b"identity_source"] == b"provided"
    assert pq.read_table(tmp_path / "out.pg.parquet").column("pg_id").to_pylist() == pq.read_table(pg).column("pg_id").to_pylist()
