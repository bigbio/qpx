"""MaxQuant feature protein groups and the feature -> pg softlink."""

from pathlib import Path

import pyarrow.parquet as pq
import pytest

from tests.converters.pg_ids_roundtrip import assert_softlink_valid, open_converted, pg_ids_by_sequence

# ---------------------------------------------------------------------------
# feature<->pg softlink round-trip (bigbio/qpx#269)
# ---------------------------------------------------------------------------


def _write_pg_link_inputs(tmp_path):
    """Minimal LFQ MaxQuant inputs exercising shared-leader + a no-pg-row case.

    PEPTIDEA -> group P1;P2, PEPTIDEB -> P1;P3 (distinct groups, shared leader P1),
    PEPTIDEC -> P9 which is absent from proteinGroups (so its feature has no pg row).
    All three features share Experiment exp1 / Raw file run1.
    """
    evidence = tmp_path / "evidence.txt"
    header = (
        "Sequence\tModified sequence\tCharge\tRaw file\tExperiment\tType\tMS/MS scan number\tm/z\t"
        "Calibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\t"
        "PEP\tLeading razor protein\tLeading proteins\tGene names\tIntensity\tMass\tid\n"
    )
    rows = (
        "PEPTIDEA\t_PEPTIDEA_\t2\trun1\texp1\tMULTI-MSMS\t100\t450.25\t30.0\t29.5\t30.5\t0.01\tP1\tP1;P2\t\t1000\t500\t1\n"
        "PEPTIDEB\t_PEPTIDEB_\t2\trun1\texp1\tMULTI-MSMS\t200\t460.25\t31.0\t30.5\t31.5\t0.02\tP1\tP1;P3\t\t2000\t520\t2\n"
        "PEPTIDEC\t_PEPTIDEC_\t2\trun1\texp1\tMULTI-MSMS\t300\t470.25\t32.0\t31.5\t32.5\t0.03\tP9\tP9\t\t3000\t540\t3\n"
    )
    evidence.write_text(header + rows)

    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text(
        "Protein IDs\tMajority protein IDs\tQ-value\tGene names\tIntensity exp1\n"
        "P1;P2\tP1;P2\t0.001\tGENEA\t1000\n"
        "P1;P3\tP1;P3\t0.002\tGENEB\t2000\n"
    )
    return evidence, protein_groups


def test_feature_pg_softlink_shared_leader_and_no_pg_row(tmp_path):
    """Full pg+feature conversion: the computed softlink links each feature to a
    real pg row with matching membership + run + carried label; shared-leader
    groups stay distinct; a feature whose group has no pg row produces no link."""
    from qpx.converters.maxquant.converter import MaxQuantConverter
    from qpx.core.constants import FEATURE, PG

    evidence, protein_groups = _write_pg_link_inputs(tmp_path)
    out = tmp_path / "out"
    MaxQuantConverter().convert(
        output_folder=str(out),
        evidence_file=str(evidence),
        protein_groups_file=str(protein_groups),
        output_prefix="mq",
        structures=[FEATURE, PG],
    )

    with open_converted(out, prefix="mq") as ds:
        feat, pg, link = assert_softlink_valid(ds)

    pg_id_by_membership = {tuple(sorted(memb)): pg_id for pg_id, memb, *_ in pg}
    ids_by_seq = pg_ids_by_sequence(feat, link)
    assert ids_by_seq["PEPTIDEA"] == [pg_id_by_membership[("P1", "P2")]]
    assert ids_by_seq["PEPTIDEB"] == [pg_id_by_membership[("P1", "P3")]]
    assert set(ids_by_seq["PEPTIDEA"]).isdisjoint(ids_by_seq["PEPTIDEB"])
    assert ids_by_seq["PEPTIDEC"] == []  # group P9 absent from pg -> no link


@pytest.mark.integration
def test_feature_pg_softlink_real_full_dataset(tmp_path_factory):
    """Real MaxQuant example (maxquant_full): the computed softlink links features
    to real pg rows with matching membership + run + carried label (bigbio/qpx#269)."""
    from qpx.converters.maxquant.converter import MaxQuantConverter
    from qpx.core.constants import FEATURE, PG

    full = Path(__file__).resolve().parent.parent / "examples" / "maxquant" / "maxquant_full"
    evidence = full / "evidence.txt.gz"
    protein_groups = full / "proteinGroups.txt"
    if not evidence.exists() or not protein_groups.exists():
        pytest.skip("maxquant_full fixture not found")

    out = tmp_path_factory.mktemp("mq_pg_ids")
    MaxQuantConverter().convert(
        output_folder=str(out),
        evidence_file=str(evidence),
        protein_groups_file=str(protein_groups),
        sdrf_file=str(full / "PXD001819.sdrf.tsv"),
        output_prefix="mq",
        structures=[FEATURE, PG],
    )

    with open_converted(out, prefix="mq") as ds:
        _feat, _pg, link = assert_softlink_valid(ds)
        unlinked = ds.features_without_pg_link().fetchall()
    assert link, "expected at least some computed feature->pg softlink edges"
    # Every evidence row names its groups in "Protein group IDs", and MS/MS-only
    # evidence keeps its label with a null intensity, so every feature links.
    assert not unlinked


# ---------------------------------------------------------------------------
# Feature protein group from evidence "Protein group IDs"
# ---------------------------------------------------------------------------

_GROUP_IDS_HEADER = (
    "Sequence\tModified sequence\tCharge\tRaw file\tExperiment\tType\tMS/MS scan number\tm/z\t"
    "Calibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\t"
    "PEP\tLeading razor protein\tLeading proteins\tProtein group IDs\tIntensity\tMass\tid\n"
)


def _group_id_evidence_row(index, sequence, razor, leading, group_ids):
    return (
        f"{sequence}\t_{sequence}_\t2\trun1\texp1\tMULTI-MSMS\t{100 + index}\t450.25\t30.0\t29.5\t30.5\t0.01\t"
        f"{razor}\t{leading}\t{group_ids}\t1000\t898.5\t{index}\n"
    )


def _write_group_id_inputs(tmp_path, rows=None):
    """evidence rows naming their groups by proteinGroups ``id``, as MaxQuant writes them.

    Group 0 is P1;P2 (two proteins), group 1 is P3 alone, group 2 is P4;P5.
    ``Leading proteins`` holds one leader per group, so it never equals a
    multi-protein group's membership.
    """
    if rows is None:
        rows = [
            ("PEPTIDEA", "P1", "P1", "0"),  # one group of two proteins
            ("PEPTIDEB", "P3", "P3", "1"),  # one group of one protein
            ("PEPTIDEC", "P4", "P3;P4", "1;2"),  # shared by groups 1 and 2, razor group 2
            ("PEPTIDED", "P3", "P3;P4", "1;2"),  # shared, razor group 1 (a single protein)
        ]
    evidence = tmp_path / "evidence.txt"
    evidence.write_text(_GROUP_IDS_HEADER + "".join(_group_id_evidence_row(i, *row) for i, row in enumerate(rows)))
    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text(
        "Protein IDs\tMajority protein IDs\tQ-value\tIntensity exp1\tid\n"
        "P1;P2\tP1;P2\t0.001\t1000\t0\n"
        "P3\tP3\t0.002\t2000\t1\n"
        "P4;P5\tP4\t0.003\t3000\t2\n"
    )
    return evidence, protein_groups


def test_feature_group_resolved_from_protein_group_ids(tmp_path):
    """Feature membership is the group ``Protein group IDs`` names, not ``Leading proteins``.

    ``Leading proteins`` lists one leader per group, so a multi-protein group's
    feature carried only its leader and a shared peptide carried the leaders of
    every group: neither matched any pg row's membership. The group must come
    from ``Protein group IDs`` (the razor group for a shared peptide), with
    ``unique`` true only for a peptide in one single-protein group.
    """
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence, protein_groups = _write_group_id_inputs(tmp_path)
    output = tmp_path / "feature.parquet"
    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output), protein_groups_path=str(protein_groups))

    observed = {
        r["sequence"]: ([p["accession"] for p in r["pg_accessions"]], r["anchor_protein"], r["unique"], r["pg_global_qvalue"])
        for r in pq.read_table(output).to_pylist()
    }
    assert observed == {
        "PEPTIDEA": (["P1", "P2"], "P1", False, pytest.approx(0.001)),
        "PEPTIDEB": (["P3"], "P3", True, pytest.approx(0.002)),
        "PEPTIDEC": (["P4", "P5"], "P4", False, pytest.approx(0.003)),
        "PEPTIDED": (["P3"], "P3", False, pytest.approx(0.002)),
    }


def test_feature_shared_peptide_without_razor_group_is_unresolved(tmp_path):
    """A shared peptide whose razor protein is in none of its groups gets no group."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence, protein_groups = _write_group_id_inputs(tmp_path, rows=[("PEPTIDEC", "P9", "P3;P4", "1;2")])
    output = tmp_path / "feature.parquet"
    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output), protein_groups_path=str(protein_groups))

    (row,) = pq.read_table(output).to_pylist()
    assert (row["pg_accessions"], row["anchor_protein"], row["unique"]) == (None, None, None)


def test_feature_pg_softlink_covers_every_group_id_feature(tmp_path):
    """With ``Protein group IDs`` every feature links to its group's pg row."""
    from qpx.converters.maxquant.converter import MaxQuantConverter
    from qpx.core.constants import FEATURE, PG

    evidence, protein_groups = _write_group_id_inputs(tmp_path)
    out = tmp_path / "out"
    MaxQuantConverter().convert(
        output_folder=str(out),
        evidence_file=str(evidence),
        protein_groups_file=str(protein_groups),
        output_prefix="mq",
        structures=[FEATURE, PG],
    )

    with open_converted(out, prefix="mq") as ds:
        feat, pg, link = assert_softlink_valid(ds)
        assert ds.features_without_pg_link().fetchall() == []
    pg_id_by_membership = {tuple(sorted(memb)): pg_id for pg_id, memb, *_ in pg}
    assert pg_ids_by_sequence(feat, link) == {
        "PEPTIDEA": [pg_id_by_membership[("P1", "P2")]],
        "PEPTIDEB": [pg_id_by_membership[("P3",)]],
        "PEPTIDEC": [pg_id_by_membership[("P4", "P5")]],
        "PEPTIDED": [pg_id_by_membership[("P3",)]],
    }


@pytest.mark.parametrize(
    ("with_protein_groups", "expected"),
    [
        (True, {"pg_accessions": "Protein IDs", "anchor_protein": "Majority protein IDs"}),
        (False, {"pg_accessions": "Leading proteins", "anchor_protein": "Leading razor protein"}),
    ],
    ids=["protein_group_ids", "leading_proteins"],
)
def test_feature_group_source_columns_in_ontology(tmp_path, with_protein_groups, expected):
    """ontology.parquet names the column the feature's group fields are read from.

    With ``Protein group IDs`` resolved through proteinGroups.txt the group is
    its ``Protein IDs`` and first ``Majority protein IDs`` entry, not
    ``Leading proteins`` / ``Leading razor protein``.
    """
    from qpx.converters.maxquant.converter import MaxQuantConverter
    from qpx.core.constants import FEATURE

    evidence, protein_groups = _write_group_id_inputs(tmp_path)
    out = tmp_path / "out"
    MaxQuantConverter().convert(
        output_folder=str(out),
        evidence_file=str(evidence),
        protein_groups_file=str(protein_groups) if with_protein_groups else None,
        output_prefix="mq",
        structures=[FEATURE],
    )

    entries = pq.read_table(out / "mq.ontology.parquet").to_pylist()
    observed = {e["field_name"]: e["source_column_name"] for e in entries if e["view"] == "feature"}
    assert {field: observed.get(field) for field in expected} == expected


def _write_group_inputs(tmp_path):
    """Evidence with ``Protein group IDs`` and a matching proteinGroups.txt.

    PEPTIDEA is in group 0 (P1;P2;P3) though its ``Leading proteins`` is only P1;
    PEPTIDEB is in group 1 (P4); PEPTIDEC is shared by both groups and assigned to
    the group of its razor protein P4. Row id 3 repeats PEPTIDEA's identity.
    """
    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        "Sequence\tModified sequence\tCharge\tRaw file\tExperiment\tType\tMS/MS scan number\tm/z\t"
        "Calibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\t"
        "PEP\tProteins\tLeading razor protein\tLeading proteins\tProtein group IDs\tIntensity\tMass\tid\n"
        "PEPTIDEA\t_PEPTIDEA_\t2\trun1\texp1\tMULTI-MSMS\t100\t450.25\t30.0\t29.5\t30.5\t0.01\tP1;P2;P3\tP1\tP1\t0\t1000\t500\t0\n"
        "PEPTIDEB\t_PEPTIDEB_\t2\trun1\texp1\tMULTI-MSMS\t200\t460.25\t31.0\t30.5\t31.5\t0.02\tP4\tP4\tP4\t1\t\t520\t1\n"
        "PEPTIDEC\t_PEPTIDEC_\t2\trun1\texp1\tMULTI-MSMS\t300\t470.25\t32.0\t31.5\t32.5\t0.03\tP1;P4\tP4\tP4;P1\t0;1\t3000\t540\t2\n"
        "PEPTIDEA\t_PEPTIDEA_\t2\trun1\texp1\tMULTI-MSMS\t100\t450.25\t30.0\t29.5\t30.5\t0.05\tP1;P2;P3\tP1\tP1\t0\t900\t500\t3\n"
    )
    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text(
        "Protein IDs\tMajority protein IDs\tPeptide counts (all)\tPeptides\tUnique peptides\tQ-value\tGene names\t"
        "Intensity exp1\tPeptides exp1\tUnique peptides exp1\tMS/MS count exp1\tSequence coverage exp1 [%]\tid\n"
        "P1;P2;P3\tP1\t2;1;1\t2\t1\t0.001\tGENEA\t4000\t1\t1\t3\t12.5\t0\n"
        "P4\tP4\t2\t2\t2\t0.002\tGENEB\t3000\t2\t2\t2\t20\t1\n"
    )
    msms = tmp_path / "msms.txt"
    msms.write_text(
        "Sequence\tModified sequence\tCharge\tRaw file\tScan number\tm/z\tRetention time\tPEP\tScore\tProteins\tMass\tEvidence ID\n"
        "PEPTIDEA\t_PEPTIDEA_\t2\trun1\t100\t450.25\t30.0\t0.01\t100\tP1\t500\t0\n"
        "PEPTIDEA\t_PEPTIDEA_\t2\trun1\t101\t450.25\t30.1\t0.05\t80\tP1\t500\t3\n"
        "PEPTIDEC\t_PEPTIDEC_\t2\trun1\t300\t470.25\t32.0\t0.03\t90\tP4\t540\t2\n"
        "PEPTIDED\t_PEPTIDED_\t2\trun1\t400\t480.25\t33.0\t0.03\t90\tP4\t560\t\n"
    )
    return evidence, protein_groups, msms


def test_feature_protein_group_follows_protein_group_ids(tmp_path):
    """The feature's group is MaxQuant's own assignment, so every feature links to its pg row."""
    from qpx.converters.maxquant.converter import MaxQuantConverter
    from qpx.core.constants import FEATURE, PG, PSM

    evidence, protein_groups, msms = _write_group_inputs(tmp_path)
    out = tmp_path / "out"
    MaxQuantConverter().convert(
        output_folder=str(out),
        msms_file=str(msms),
        evidence_file=str(evidence),
        protein_groups_file=str(protein_groups),
        output_prefix="mq",
        structures=[PSM, FEATURE, PG],
    )

    by_seq = {row["sequence"]: row for row in pq.read_table(out / "mq.feature.parquet").to_pylist()}
    assert [p["accession"] for p in by_seq["PEPTIDEA"]["pg_accessions"]] == ["P1", "P2", "P3"]
    assert by_seq["PEPTIDEA"]["unique"] is False
    assert by_seq["PEPTIDEB"]["unique"] is True
    # A shared peptide goes to its razor protein's group and is not unique to it.
    assert [p["accession"] for p in by_seq["PEPTIDEC"]["pg_accessions"]] == ["P4"]
    assert (by_seq["PEPTIDEC"]["anchor_protein"], by_seq["PEPTIDEC"]["unique"]) == ("P4", False)
    assert by_seq["PEPTIDEC"]["pg_global_qvalue"] == pytest.approx(0.002)
    # An evidence row MaxQuant did not quantify keeps its label with a null intensity.
    assert by_seq["PEPTIDEB"]["intensities"] == [{"label": "LFQ", "intensity": None}]

    with open_converted(out, prefix="mq") as ds:
        assert ds.features_without_pg_link().to_df().empty

    # Each PSM links to its feature through "Evidence ID", collapsed duplicates included.
    feature_ids = {row["sequence"]: row["feature_id"] for row in pq.read_table(out / "mq.feature.parquet").to_pylist()}
    psms = pq.read_table(out / "mq.psm.parquet").to_pylist()
    assert [(row["sequence"], row["feature_id"]) for row in psms] == [
        ("PEPTIDEA", feature_ids["PEPTIDEA"]),
        ("PEPTIDEA", feature_ids["PEPTIDEA"]),
        ("PEPTIDEC", feature_ids["PEPTIDEC"]),
        ("PEPTIDED", None),
    ]


def test_pg_counts_are_per_experiment(tmp_path):
    """Peptide counts, coverage and the MS/MS count come from the experiment's own columns."""
    from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

    evidence, protein_groups, _ = _write_group_inputs(tmp_path)
    output = tmp_path / "pg.parquet"
    with MaxQuantPgAdapter() as adapter:
        adapter.convert(str(protein_groups), str(output), evidence_path=str(evidence))

    row = next(r for r in pq.read_table(output).to_pylist() if r["anchor_protein"] == "P1")
    assert row["peptide_counts"] == {"unique_sequences": 1, "total_sequences": 1}
    assert row["feature_counts"] is None
    assert row["sequence_coverage"] == pytest.approx(12.5)
    assert row["peptides"] == [
        {"protein_name": "P1", "peptide_count": 2},
        {"protein_name": "P2", "peptide_count": 1},
        {"protein_name": "P3", "peptide_count": 1},
    ]
    assert row["additional_intensities"] == [
        {"label": "LFQ", "intensities": [{"intensity_name": "spectral_count", "intensity_value": 3.0}]}
    ]
