"""Cross-linking schema and round-trip tests."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from tests.conftest import make_psm_record


def test_cross_link_schema():
    """PSM schema includes cross_links as list<struct> with expected fields."""
    from qpx.core.data import PsmSchema

    field_names = [f.name for f in PsmSchema.get_arrow_schema()]
    assert "cross_links" in field_names

    field = PsmSchema.get_arrow_schema().field("cross_links")
    assert pa.types.is_list(field.type)
    assert pa.types.is_struct(field.type.value_type)

    xl_type = field.type.value_type
    xl_field_names = [xl_type.field(i).name for i in range(xl_type.num_fields)]
    for expected in (
        "xl_type",
        "partner_sequence",
        "partner_peptidoform",
        "donor_position",
        "acceptor_position",
        "linker_name",
        "linker_accession",
        "linker_mass",
    ):
        assert expected in xl_field_names


def test_cross_link_round_trip(tmp_path):
    """Round-trip: inter XL, null XL, dead-end XL, multiple XLs."""
    from qpx.writers import PsmWriter

    # Inter cross-link
    record = make_psm_record()
    record["cross_links"] = [
        {
            "xl_type": "inter",
            "partner_sequence": "PEPTIDEK",
            "partner_peptidoform": "PEPTIDEK",
            "donor_position": 5,
            "acceptor_position": 3,
            "linker_name": "DSS",
            "linker_accession": "XLMOD:02001",
            "linker_mass": 138.0681,
        }
    ]
    out = tmp_path / "xl.psm.parquet"
    with PsmWriter(out) as w:
        w.write_batch([record])
    table = pq.read_table(out)
    xl = table.column("cross_links")[0].as_py()
    assert len(xl) == 1
    assert xl[0]["xl_type"] == "inter"
    assert xl[0]["linker_mass"] == pytest.approx(138.0681)

    # Null cross_links for regular PSM
    record2 = make_psm_record()
    out2 = tmp_path / "regular.psm.parquet"
    with PsmWriter(out2) as w:
        w.write_batch([record2])
    assert pq.read_table(out2).column("cross_links")[0].as_py() is None

    # Dead-end cross-link
    record3 = make_psm_record()
    record3["cross_links"] = [
        {
            "xl_type": "dead-end",
            "partner_sequence": None,
            "partner_peptidoform": None,
            "donor_position": 5,
            "acceptor_position": None,
            "linker_name": "DSS",
            "linker_accession": None,
            "linker_mass": 156.0786,
        }
    ]
    out3 = tmp_path / "de.psm.parquet"
    with PsmWriter(out3) as w:
        w.write_batch([record3])
    xl = pq.read_table(out3).column("cross_links")[0].as_py()
    assert xl[0]["xl_type"] == "dead-end"
    assert xl[0]["partner_sequence"] is None

    # Multiple cross-links
    record4 = make_psm_record()
    record4["cross_links"] = [
        {
            "xl_type": "inter",
            "partner_sequence": "PEPTIDEK",
            "partner_peptidoform": "PEPTIDEK",
            "donor_position": 5,
            "acceptor_position": 3,
            "linker_name": "DSSO",
            "linker_accession": None,
            "linker_mass": 158.0038,
        },
        {
            "xl_type": "intra",
            "partner_sequence": None,
            "partner_peptidoform": None,
            "donor_position": 2,
            "acceptor_position": 7,
            "linker_name": "DSSO",
            "linker_accession": None,
            "linker_mass": 158.0038,
        },
    ]
    out4 = tmp_path / "multi_xl.psm.parquet"
    with PsmWriter(out4) as w:
        w.write_batch([record4])
    xl = pq.read_table(out4).column("cross_links")[0].as_py()
    assert len(xl) == 2
    assert xl[0]["xl_type"] == "inter"
    assert xl[1]["xl_type"] == "intra"


@pytest.mark.parametrize(
    ("peptidoform", "position", "tagged"),
    [
        ("PEPK", 4, "PEPK[#XL1]"),
        ("PEM[Oxidation]K", 3, "PEM[Oxidation][#XL1]K"),
        ("PEPK", 0, "[#XL1]-PEPK"),
        ("[Acetyl]-PEPK", 0, "[#XL1][Acetyl]-PEPK"),
        ("[Acetyl]-PEPK", 1, "[Acetyl]-P[#XL1]EPK"),
        ("PEPK", 5, "PEPK-[#XL1]"),
    ],
)
def test_tag_goes_after_the_residue_and_its_modifications(peptidoform, position, tagged):
    from qpx.converters.mzidentml.crosslinks import _tag_residue

    assert _tag_residue(peptidoform, position, "#XL1") == tagged


def _link(**fields):
    return {
        "partner_peptidoform": "SEK",
        "donor_position": 4,
        "acceptor_position": 3,
        "linker_accession": "XLMOD:02001",
        "linker_mass": 138.0681,
        **fields,
    }


def test_crosslinked_peptidoform_names_the_whole_species():
    from qpx.converters.mzidentml.crosslinks import crosslinked_peptidoform

    assert crosslinked_peptidoform("PEPK", None) is None
    assert crosslinked_peptidoform("PEPK", [_link()]) == "PEPK[XLMOD:02001#XL1]//SEK[#XL1]"
    assert crosslinked_peptidoform("PEPK", [_link(linker_accession=None)]) == "PEPK[+138.0681#XL1]//SEK[#XL1]"
    # A looplink stays on one peptide, and noncovalent peptides share no site.
    assert crosslinked_peptidoform("PEKK", [_link(partner_peptidoform=None, donor_position=3, acceptor_position=4)]) == (
        "PEK[XLMOD:02001#XL1]K[#XL1]"
    )
    assert crosslinked_peptidoform("PEPK", [_link(donor_position=None, acceptor_position=None)]) == "PEPK//SEK"


def test_partner_fields_cover_the_partner_and_a_decoy_on_either_side():
    from qpx.converters.mzidentml.crosslinks import partner_fields

    parsed = {
        "peptide_evidence": {
            "a": {"db_ref": "d1", "start": 10, "is_decoy": False},
            "b": {"db_ref": "d2", "start": 20, "is_decoy": True},
            "b2": {"db_ref": "d3", "start": 30, "is_decoy": False},
        },
        "db_sequences": {"d1": "P1", "d2": "DECOY_P2", "d3": "P3"},
    }
    alpha = {"peptide_evidence_refs": ["a"]}

    fields = partner_fields(alpha, {"peptide_evidence_refs": ["b", "b2"]}, parsed)
    assert fields == {
        "protein_starts": [10],
        "partner_protein_accessions": ["DECOY_P2", "P3"],
        "partner_protein_starts": [20, 30],
        "crosslink_decoy_class": "target_decoy",
        "is_decoy": True,
    }
    assert partner_fields(alpha, None, parsed) == {"protein_starts": [10]}
    assert partner_fields(alpha, {"peptide_evidence_refs": ["b2"]}, parsed)["crosslink_decoy_class"] == "target_target"
    assert (
        partner_fields({"peptide_evidence_refs": ["b"]}, {"peptide_evidence_refs": ["b"]}, parsed)["crosslink_decoy_class"]
        == "decoy_decoy"
    )
