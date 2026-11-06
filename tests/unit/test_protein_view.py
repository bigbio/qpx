"""The protein view: from FASTA headers, from proteinGroups.txt, and through a merge."""

import pyarrow.parquet as pq
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.converters.maxquant.protein_fields import write_protein_view
from qpx.core.protein_annotation import fasta_header_fields
from qpx.dataset import Dataset
from qpx.writers import PgWriter
from tests.conftest import make_pg_record

_ALBUMIN = "sp|P02768|ALBU_HUMAN Serum albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2"


def test_uniprot_header_fields():
    assert fasta_header_fields(_ALBUMIN) == {
        "entry_name": "ALBU_HUMAN",
        "description": "Serum albumin",
        "organism": "Homo sapiens",
        "gene_name": "ALB",
        "taxonomy_id": 9606,
    }


def test_other_header_yields_what_it_states():
    assert fasta_header_fields("P1 some protein") == {
        "entry_name": None,
        "description": "some protein",
        "organism": None,
        "gene_name": None,
        "taxonomy_id": None,
    }


def test_proteingroups_cells_align_with_the_members(tmp_path):
    groups = tmp_path / "proteinGroups.txt"
    groups.write_text(
        "Protein IDs\tFasta headers\tSequence lengths\n"
        f"P02768;P1\t>{_ALBUMIN};other\t609;100\n"
        "P2;P3\tonly one header\t10;20\n"
        "REV__P4\t\t30\n"
    )

    write_protein_view(str(groups), str(tmp_path / "x.protein.parquet"))

    rows = {row["accession"]: row for row in pq.read_table(tmp_path / "x.protein.parquet").to_pylist()}
    assert sorted(rows) == ["P02768", "P1", "P2", "P3"]
    assert rows["P02768"]["gene_name"] == "ALB"
    assert rows["P02768"]["length"] == 609
    assert rows["P02768"]["sequence"] is None
    assert rows["P1"]["description"] is None and rows["P1"]["length"] == 100
    assert rows["P3"]["description"] is None and rows["P3"]["length"] == 20


def _dataset(tmp_path):
    record = make_pg_record()
    record.update(pg_accessions=["P02768", "P_NOT_IN_FASTA"], anchor_protein="P02768")
    with PgWriter(tmp_path / "x.pg.parquet") as writer:
        writer.write_batch([record])
    fasta = tmp_path / "db.fasta"
    fasta.write_text(f">{_ALBUMIN}\nMKWV\n>DECOY_sp|P02768|ALBU_HUMAN\nVWKM\n")
    return fasta


def test_transform_writes_the_view_and_records_the_fasta(tmp_path):
    fasta = _dataset(tmp_path)

    result = CliRunner().invoke(
        qpx_main,
        ["transform", "protein-properties", "--dataset", str(tmp_path), "--fasta", str(fasta)]
        + ["--fasta-uri", "https://example.org/db.fasta", "--in-place"],
    )

    assert result.exit_code == 0, result.output
    with Dataset(str(tmp_path), file_prefix="x") as dataset:
        rows = dataset.protein.to_df().to_dict("records")
        assert not dataset.validate()["protein"].errors
    assert [(row["accession"], row["sequence"], row["length"], row["taxonomy_id"]) for row in rows] == [
        ("P02768", "MKWV", 4, 9606)
    ]
    parameters = {p["key"]: p["value"] for p in pq.read_table(tmp_path / "x.provenance.parquet").to_pylist()[-1]["parameters"]}
    assert parameters["fasta_uri"] == "https://example.org/db.fasta"
    assert parameters["fasta_sha256"]


def test_a_recorded_value_stays_and_the_fasta_fills_the_rest(tmp_path):
    fasta = _dataset(tmp_path)
    from qpx.writers import ProteinWriter

    with ProteinWriter(tmp_path / "x.protein.parquet") as writer:
        writer.write_batch([{"accession": "P02768", "description": "Producer's name", "length": 609}])

    result = CliRunner().invoke(
        qpx_main, ["transform", "protein-properties", "--dataset", str(tmp_path), "--fasta", str(fasta), "--in-place"]
    )

    assert result.exit_code == 0, result.output
    [row] = pq.read_table(tmp_path / "x.protein.parquet").to_pylist()
    assert (row["description"], row["length"], row["sequence"], row["gene_name"]) == ("Producer's name", 609, "MKWV", "ALB")
