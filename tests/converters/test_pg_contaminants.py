"""Protein-group contaminant flags must consider members beyond the leader."""

import pyarrow.parquet as pq
import pytest

from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter
from qpx.converters.utils import canonical_accession
from tests.converters.test_diann_converter import _diann_row, _pg_rows_single_run
from tests.converters.test_openms_consensus import _SHARED_LEADER_CONSENSUSXML

_MEMBER_CASES = [
    ("Cont_Q7SIH1", True),
    ("sp|Cont_Q7SIH1|A2MG_BOVIN", True),
    ("tr|cOnT_P05783|K1C18_HUMAN", True),
    ("CONTAM_P02769", True),
    ("sp|P05783|CONT_HUMAN", False),
]


@pytest.mark.parametrize(("member", "expected"), _MEMBER_CASES)
@pytest.mark.parametrize("streaming", [False, True])
def test_openms_pg_contaminant_checks_nonleading_members(tmp_path, member, expected, streaming):
    """Flag the whole mixed group and retain its membership and quantification."""
    anchor = "A0A024RBG1"
    source = tmp_path / "groups.consensusXML"
    source.write_text(
        _SHARED_LEADER_CONSENSUSXML.replace('accession="A"', f'accession="{anchor}"')
        .replace('accession="B"', f'accession="{member}"')
        .replace('accession="C"', 'accession="P12345"')
    )
    written = OpenMSConsensusConverter().convert(str(source), str(tmp_path / "out"), structures=("pg",), streaming=streaming)
    rows = pq.read_table(written["pg"]).to_pylist()

    assert len(rows) == 2
    assert {tuple(row["pg_accessions"]): (row["contaminant"], row["intensity"]) for row in rows} == {
        (anchor, canonical_accession(member)): (expected, 1000.0),
        (anchor, "P12345"): (False, 3000.0),
    }
    assert all(row["anchor_protein"] == anchor for row in rows)


@pytest.mark.parametrize(("member", "expected"), _MEMBER_CASES)
def test_diann_pg_contaminant_checks_nonleading_members(tmp_path, member, expected):
    """DIA-NN must apply the same whole-group flag without dropping or requantifying it."""
    anchor = "A0A024RBG1"
    groups = [f"{anchor};{member}", f"{anchor};P12345"]
    report = [_diann_row(**{"Protein.Group": group}) for group in groups]
    matrix = [{"Protein.Group": group, "Protein.Names": "N1", "Genes": "G1", "run_A": 900.0} for group in groups]

    rows = _pg_rows_single_run(tmp_path, report, matrix_rows=matrix)

    assert len(rows) == 2
    assert {tuple(row["pg_accessions"]): (row["contaminant"], row["intensity"]) for row in rows} == {
        (anchor, member): (expected, 1000.0),
        (anchor, "P12345"): (False, 1000.0),
    }
    assert all(row["anchor_protein"] == anchor for row in rows)
    assert all(row["pg_qvalue"] == 0.002 and row["global_qvalue"] == 0.003 for row in rows)
