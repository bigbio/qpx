"""One rule per dataset for run names, protein accessions, modification names and MS1 labels."""

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.converters.cdap.base_adapter import CdapBaseAdapter
from qpx.converters.channel_labels import normalize_label, plex_for_columns, relabel_intensities_parquet
from qpx.converters.diann.base_adapter import decoy_protein_group_sql, run_name_sql
from qpx.converters.diann.pg_adapter import _EXT_RE
from qpx.converters.openms.converter import _run_key_candidates
from qpx.converters.ptm import from_proforma
from qpx.core.files import run_file_stem
from qpx.pdc.metadata import _strip_raw_extension
from qpx.transforms.accession_normalizer import normalize_parquet, to_bare, to_full


@pytest.mark.parametrize(
    "name",
    ["run.raw", "run.RAW", "run.mzML", "run.mzml.gz", "run.d", "run.d.zip", "run.wiff", "run.wiff2", "run.htrms", "run.mgf"],
)
def test_every_converter_stems_a_run_file_alike(name):
    con = duckdb.connect()
    assert run_file_stem(name) == "run"
    assert CdapBaseAdapter.strip_run_extension(name) == "run"
    assert _strip_raw_extension(name) == "run"
    assert _EXT_RE.sub("", name) == "run"
    assert "run" in _run_key_candidates(name)
    assert con.execute(f"SELECT {run_name_sql(repr(name))}").fetchone()[0] == "run"


def test_diann_decoy_group_carries_bare_accessions():
    sql = decoy_protein_group_sql("'sp|O95299|NDUAA_HUMAN;P12345'", "1")
    assert duckdb.connect().execute(f"SELECT {sql}").fetchone()[0] == "DECOY_O95299;DECOY_P12345"


@pytest.mark.parametrize(
    "full, bare",
    [
        ("sp|P04114|APOB_HUMAN", "P04114"),
        ("REV__sp|P04114|APOB_HUMAN", "REV__P04114"),
        ("DECOY_sp|P04114|APOB_HUMAN", "DECOY_P04114"),
        ("CONTAM_sp|CONTAM_P02768|X", "CONTAM_P02768"),
    ],
)
def test_normalizer_strips_prefixed_decoys_and_contaminants(full, bare):
    assert to_bare(full) == bare


def test_normalizer_reverse_restores_a_decoy_marker():
    lookup = {"P04114": "sp|P04114|APOB_HUMAN"}
    assert to_full("REV__P04114", lookup) == "REV__sp|P04114|APOB_HUMAN"
    assert to_full("REV__sp|P04114|APOB_HUMAN", lookup) == "REV__sp|P04114|APOB_HUMAN"
    assert to_full("Q99999", lookup) == "Q99999"


def test_normalizer_moves_entry_names_to_pg_names(tmp_path):
    source = tmp_path / "in.pg.parquet"
    pq.write_table(
        pa.table(
            {
                "pg_accessions": [
                    ["sp|P04114|A_HUMAN", "sp|P02768|B_HUMAN"],
                    ["sp|P68871|C_HUMAN", "P69905"],
                    ["sp|P01308|E_HUMAN"],
                ],
                "pg_names": [None, None, ["kept"]],
                "anchor_protein": ["sp|P04114|A_HUMAN", "sp|P68871|C_HUMAN", "sp|P01308|E_HUMAN"],
            }
        ),
        source,
    )
    output = tmp_path / "out.pg.parquet"
    normalize_parquet(source, output, direction="forward")
    table = pq.read_table(output)
    assert table.column("pg_accessions").to_pylist() == [["P04114", "P02768"], ["P68871", "P69905"], ["P01308"]]
    # Names align with accessions, so a group with an unnamed member has none.
    assert table.column("pg_names").to_pylist() == [["A_HUMAN", "B_HUMAN"], None, ["kept"]]


def test_modification_name_is_the_unimod_title_without_metadata():
    _, modifications = from_proforma("AC[UNIMOD:4]K", "ACK")
    assert [(entry["name"], entry["accession"]) for entry in modifications] == [("Carbamidomethyl", "UNIMOD:4")]


def test_lfq_relabelling_keeps_ms1_labels(tmp_path):
    """Three dimethyl channels are not three LFQ quantities, whatever the experiment type."""
    source, output = tmp_path / "in.parquet", tmp_path / "out.parquet"
    entry = pa.struct([("label", pa.string()), ("intensity", pa.float32())])
    rows = [[{"label": label, "intensity": 1.0} for label in ("DIMETHYL0", "DIMETHYL4", "DIMETHYL8")]]
    pq.write_table(pa.table({"intensities": pa.array(rows, pa.list_(entry))}), source)
    relabel_intensities_parquet(str(source), str(output), {}, is_lfq=True)
    labels = [item["label"] for item in pq.read_table(output).column("intensities").to_pylist()[0]]
    assert labels == ["DIMETHYL0", "DIMETHYL4", "DIMETHYL8"]


@pytest.mark.parametrize(
    "spelling, label",
    [
        ("iTRAQ reagent 114", "ITRAQ114"),
        ("iTRAQ4plex-114", "ITRAQ114"),
        ("TMTpro-134C", "TMT134C"),
        ("TMT6plex-127", "TMT127"),
        ("TMT10plex-127", "TMT10plex-127"),  # not a channel of the plex it names
        ("label free sample", "LFQ"),
    ],
)
def test_reporter_spellings_resolve_to_their_channel(spelling, label):
    assert normalize_label(spelling) == label


def test_partial_plex_takes_the_position_of_the_table_size():
    labels = {"ITRAQ114", "ITRAQ115", "ITRAQ116", "ITRAQ117"}
    assert {v: k for k, v in plex_for_columns("iTRAQ", labels, 4).items()}["ITRAQ114"] == 1
    assert {v: k for k, v in plex_for_columns("iTRAQ", labels, 8).items()}["ITRAQ114"] == 2
