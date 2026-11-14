"""Modification candidates (psm/feature) and sequence edits (pepmap/proteoform)."""

import pyarrow.parquet as pq
import pytest

from tests.conftest import make_psm_record


def candidate(position, aa, name, probability, assigned, accession=None, mass=None, group=None):
    return {
        "position": position,
        "amino_acid": aa,
        "name": name,
        "accession": accession,
        "mass_delta": mass,
        "probability": probability,
        "assigned": assigned,
        "alternative_group": group,
        "scores": None,
    }


# PEPS[Phospho]TK[Acetyl]R: phospho on S4 or T5; acetyl or trimethyl on K6.
CANDIDATES = [
    candidate(4, "S", "Phospho", 0.95, True, "UNIMOD:21", 79.966331, group=1),
    candidate(5, "T", "Phospho", 0.05, False, "UNIMOD:21", 79.966331, group=1),
    candidate(6, "K", "Acetyl", 0.71, True, "UNIMOD:1", 42.010565, group=2),
    candidate(6, "K", "Trimethyl", 0.29, False, "UNIMOD:37", 42.04695, group=2),
]


def _psm(**overrides):
    record = make_psm_record(sequence="PEPSTKR", peptidoform="PEPS[Phospho]TK[Acetyl]R")
    record["modification_candidates"] = CANDIDATES
    record.update(overrides)
    return record


def _validate(path, structure):
    import qpx

    ds = qpx.open_dataset(str(path))
    try:
        return ds.validate([structure], strict=True)[structure]
    finally:
        ds.close()


def test_candidates_round_trip_and_validate(tmp_path):
    from qpx.writers import PsmWriter

    with PsmWriter(tmp_path / "exp.psm.parquet") as w:
        w.write_batch([_psm(), make_psm_record(scan=[1002])])
    table = pq.read_table(tmp_path / "exp.psm.parquet")
    rows = table.column("modification_candidates").to_pylist()
    assert [(c["position"], c["name"], c["assigned"]) for c in rows[0]] == [
        (4, "Phospho", True),
        (5, "Phospho", False),
        (6, "Acetyl", True),
        (6, "Trimethyl", False),
    ]
    assert rows[1] is None, "a PSM without alternatives leaves the column null"
    assert _validate(tmp_path, "psm").is_valid


def test_feature_candidates_round_trip_and_validate(tmp_path):
    from qpx.writers import FeatureWriter
    from tests.conftest import make_feature_record

    record = make_feature_record()
    record["modification_candidates"] = CANDIDATES
    with FeatureWriter(tmp_path / "exp.feature.parquet") as w:
        w.write_batch([record])
    table = pq.read_table(tmp_path / "exp.feature.parquet")
    assert len(table.column("modification_candidates").to_pylist()[0]) == len(CANDIDATES)
    assert _validate(tmp_path, "feature").is_valid


def test_a_file_without_the_column_stays_valid(tmp_path):
    """The column is optional: files written before it existed still validate."""
    from qpx.writers import PsmWriter

    with PsmWriter(tmp_path / "exp.psm.parquet") as w:
        w.write_batch([make_psm_record()])
    table = pq.read_table(tmp_path / "exp.psm.parquet")
    old = table.drop_columns(["modification_candidates"])
    pq.write_table(old.replace_schema_metadata(table.schema.metadata), tmp_path / "exp.psm.parquet")
    assert _validate(tmp_path, "psm").is_valid


@pytest.mark.parametrize(
    ("candidates", "check"),
    [
        ([candidate(6, "K", "Acetyl", 1.5, True)], "value_out_of_range"),
        (
            [candidate(6, "K", "Acetyl", 0.8, True, group=1), candidate(6, "K", "Trimethyl", 0.4, False, group=1)],
            "candidate_probability_sum",
        ),
    ],
)
def test_bad_candidate_probabilities_fail_strict_validation(tmp_path, candidates, check):
    from qpx.writers import PsmWriter

    with PsmWriter(tmp_path / "exp.psm.parquet") as w:
        w.write_batch([_psm(modification_candidates=candidates)])
    result = _validate(tmp_path, "psm")
    assert any(i.check == check for i in result.errors)


def test_rounded_group_probabilities_stay_valid(tmp_path):
    """MaxQuant rounds to three decimals, so one phosphate's group can sum to 1.001."""
    from qpx.writers import PsmWriter

    candidates = [
        candidate(4, "S", "Phospho", 0.501, True, "UNIMOD:21", 79.966331, group=1),
        candidate(5, "T", "Phospho", 0.5, False, "UNIMOD:21", 79.966331, group=1),
    ]
    with PsmWriter(tmp_path / "exp.psm.parquet") as w:
        w.write_batch([_psm(modification_candidates=candidates)])
    assert _validate(tmp_path, "psm").is_valid


def test_co_occurring_modifications_do_not_sum_across_groups(tmp_path):
    """A certain SILAC label on K6 co-occurs with the acetyl/trimethyl alternatives."""
    from qpx.writers import PsmWriter

    candidates = [
        candidate(6, "K", "Label:13C(6)15N(2)", 1.0, True, "UNIMOD:259"),
        candidate(6, "K", "Acetyl", 0.71, True, "UNIMOD:1", group=1),
        candidate(6, "K", "Trimethyl", 0.29, False, "UNIMOD:37", group=1),
    ]
    with PsmWriter(tmp_path / "exp.psm.parquet") as w:
        w.write_batch(
            [_psm(sequence="PEPTIDEK", peptidoform="PEPTIDEK[UNIMOD:259][UNIMOD:1]", modification_candidates=candidates)]
        )
    assert _validate(tmp_path, "psm").is_valid


def test_a_wrongly_typed_candidate_field_is_reported_not_raised(tmp_path):
    """A string probability is a type mismatch; the probability checks are skipped."""
    import pyarrow as pa

    from qpx.writers import PsmWriter

    with PsmWriter(tmp_path / "exp.psm.parquet") as w:
        w.write_batch([_psm()])
    table = pq.read_table(tmp_path / "exp.psm.parquet")
    rows = table.column("modification_candidates").to_pylist()
    for row in rows:
        for c in row or []:
            c["probability"] = str(c["probability"])
    struct = pa.struct(
        [
            ("position", pa.int32()),
            ("amino_acid", pa.string()),
            ("name", pa.string()),
            ("accession", pa.string()),
            ("mass_delta", pa.float64()),
            ("probability", pa.string()),
            ("assigned", pa.bool_()),
            ("alternative_group", pa.int32()),
            (
                "scores",
                pa.list_(pa.struct([("score_name", pa.string()), ("score_value", pa.float64()), ("higher_better", pa.bool_())])),
            ),
        ]
    )
    index = table.schema.get_field_index("modification_candidates")
    table = table.set_column(index, "modification_candidates", pa.array(rows, type=pa.list_(struct)))
    pq.write_table(table, tmp_path / "exp.psm.parquet")
    result = _validate(tmp_path, "psm")
    assert any(i.check == "type_mismatch" for i in result.errors)


def edit(start, end, ref, alt, name=None, protein="P04637"):
    return {
        "protein_accession": protein,
        "start": start,
        "end": end,
        "reference": ref,
        "alternative": alt,
        "name": name,
        "accession": None,
    }


# One sequence's edits, in order.
VALID_EDITS = [
    edit(10, 11, "KE", "QK"),  # two adjacent residues replaced
    edit(30, 30, "G", ""),  # deletion
    edit(50, 49, "", "AA"),  # insertion before residue 50
    edit(175, 175, "R", "H", "p.Arg175His"),  # substitution
]


def _pepmap(path, edits):
    from qpx.writers import PepMapWriter

    with PepMapWriter(path / "exp.pepmap.parquet") as w:
        w.write_batch(
            [
                {
                    "sequence": "VVRHCPHHER",
                    "peptidoform": "VVRHCPHHER",
                    "pg_accessions": [{"accession": "P04637", "start": 172, "end": 181, "pre": "K", "post": "C"}],
                    "sequence_edits": edits,
                    "is_unique": True,
                }
            ]
        )


def test_sequence_edits_cover_substitutions_spans_deletions_and_insertions(tmp_path):
    _pepmap(tmp_path, VALID_EDITS)
    stored = pq.read_table(tmp_path / "exp.pepmap.parquet").column("sequence_edits")[0].as_py()
    assert [(e["start"], e["end"], e["reference"], e["alternative"]) for e in stored] == [
        (10, 11, "KE", "QK"),
        (30, 30, "G", ""),
        (50, 49, "", "AA"),
        (175, 175, "R", "H"),
    ]
    assert _validate(tmp_path, "pepmap").is_valid


@pytest.mark.parametrize(
    "bad",
    [
        edit(10, 12, "KE", "QK"),  # reference shorter than its span
        edit(0, 0, "", "M"),  # start before the protein
        edit(20, 20, "r", "h"),  # not one-letter codes
        edit(20, 19, "", ""),  # changes nothing
        edit(20, 20, "R", "R"),  # reference equals alternative
    ],
)
def test_malformed_sequence_edits_fail_strict_validation(tmp_path, bad):
    _pepmap(tmp_path, [bad])
    result = _validate(tmp_path, "pepmap")
    assert any(i.check == "malformed_sequence_edit" for i in result.errors)


@pytest.mark.parametrize(
    "edits",
    [
        [edit(175, 175, "R", "H"), edit(175, 175, "R", "C")],  # two changes of one residue
        [edit(10, 11, "KE", "QK"), edit(11, 12, "EV", "DA")],  # overlapping spans
        [edit(50, 49, "", "A"), edit(50, 49, "", "G")],  # two insertions at one point
        [edit(30, 30, "G", "A"), edit(20, 20, "S", "T")],  # out of order
    ],
)
def test_conflicting_edits_of_one_sequence_fail_strict_validation(tmp_path, edits):
    _pepmap(tmp_path, edits)
    result = _validate(tmp_path, "pepmap")
    assert any(i.check == "conflicting_sequence_edits" for i in result.errors)


def test_compatible_neighbouring_edits_and_separate_mappings_pass(tmp_path):
    """An insertion before 50 then a change of 50 is one sequence; R175H and R175C on two proteins are two."""
    _pepmap(
        tmp_path,
        [
            edit(50, 49, "", "A"),
            edit(50, 50, "L", "P"),
            edit(175, 175, "R", "H"),
            edit(175, 175, "R", "C", protein="P04637-2"),
        ],
    )
    assert _validate(tmp_path, "pepmap").is_valid


def test_conflicting_proteoform_edits_fail_strict_validation(tmp_path):
    from qpx.writers import ProteoformWriter

    record = {
        "family": "P04637",
        "accession": "P04637_VAR_X",
        "kind": "variant",
        "edits": [
            {"start": 175, "end": 175, "reference": "R", "alternative": "H", "name": None, "accession": None},
            {"start": 175, "end": 175, "reference": "R", "alternative": "C", "name": None, "accession": None},
        ],
        "grouped_runs": ["run_01"],
        "is_decoy": False,
    }
    with ProteoformWriter(tmp_path / "exp.proteoform.parquet") as w:
        w.write_batch([record])
    result = _validate(tmp_path, "proteoform")
    assert any(i.check == "conflicting_sequence_edits" for i in result.errors)


def test_a_variant_proteoform_carries_its_edits(tmp_path):
    from qpx.writers import ProteoformWriter

    record = {
        "family": "P04637",
        "accession": "P04637_VAR_R175H",
        "kind": "variant",
        "description": "p.Arg175His",
        "edits": [{"start": 175, "end": 175, "reference": "R", "alternative": "H", "name": "p.Arg175His", "accession": None}],
        "grouped_runs": ["run_01"],
        "is_decoy": False,
    }
    with ProteoformWriter(tmp_path / "exp.proteoform.parquet") as w:
        w.write_batch([record])
    stored = pq.read_table(tmp_path / "exp.proteoform.parquet").column("edits")[0].as_py()
    assert stored[0]["reference"] == "R" and stored[0]["alternative"] == "H"
    assert _validate(tmp_path, "proteoform").is_valid


def test_site_probabilities_keep_every_candidate_residue():
    """S 0.95 is placed; T 0.03 and Y 0.02 are its alternatives, and only the placed site gets a score."""
    from qpx.converters.localization import localize_modifications
    from qpx.converters.ptm import from_proforma

    _, modifications = from_proforma("PEPS[UNIMOD:21]TYK", "PEPSTYK")
    localized, candidates = localize_modifications(modifications, "PEPSTYK", {"UNIMOD:21": {4: 0.95, 5: 0.03, 6: 0.02}})
    assert [(c["position"], c["amino_acid"], c["probability"], c["assigned"], c["alternative_group"]) for c in candidates] == [
        (4, "S", 0.95, True, 1),
        (5, "T", 0.03, False, 1),
        (6, "Y", 0.02, False, 1),
    ]
    assert (candidates[0]["name"], candidates[0]["accession"]) == ("Phospho", "UNIMOD:21")
    assert localized[0]["positions"][0]["scores"][0]["score_name"] == "localization_probability"
    assert modifications[0]["positions"][0]["scores"] is None  # the from_proforma cache is untouched
    assert localize_modifications(modifications, "PEPSTYK", {}) == (modifications, None)


def test_site_probabilities_of_several_placements_stay_ungrouped():
    """Two phosphates over three sites sum to 2; each modification has its own group otherwise."""
    from qpx.converters.localization import localize_modifications
    from qpx.converters.ptm import from_proforma

    sequence = "PS[UNIMOD:21]T[UNIMOD:21]YM[UNIMOD:35]K"
    _, modifications = from_proforma(sequence, "PSTYMK")
    _, candidates = localize_modifications(
        modifications, "PSTYMK", {"UNIMOD:21": {2: 0.9, 3: 0.8, 4: 0.3}, "UNIMOD:35": {5: 1.0}}
    )
    assert {c["accession"]: c["alternative_group"] for c in candidates} == {"UNIMOD:21": None, "UNIMOD:35": 1}
    # Independent occupancies of one placed site can exceed 1 and cannot be a group.
    _, modifications = from_proforma("PS[UNIMOD:21]TK", "PSTK")
    _, candidates = localize_modifications(modifications, "PSTK", {"UNIMOD:21": {2: 0.9, 3: 0.6}})
    assert [c["alternative_group"] for c in candidates] == [None, None]
