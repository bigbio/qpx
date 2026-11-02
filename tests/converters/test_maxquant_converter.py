"""MaxQuant converter integration tests using small example dataset.

Runs the full MaxQuant conversion pipeline once (module scope) and validates
that psm.parquet, feature.parquet, and ontology.parquet are written with
correct schemas and plausible values.
"""

from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import pytest

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples" / "maxquant" / "maxquant_simple"

_MSMS = EXAMPLES_DIR / "msms.txt"
_EVIDENCE = EXAMPLES_DIR / "evidence.txt"
_SDRF = EXAMPLES_DIR / "sdrf.tsv"

_PREFIX = "maxquant_test"


def test_pg_fails_before_writing_when_experiment_has_no_run_mapping(tmp_path):
    """An unmapped MaxQuant Experiment must not be swallowed as skipped PG rows."""
    from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text("Protein IDs\tMajority protein IDs\tIntensity experiment_1\nP12345\tP12345\t1000\n")
    output = tmp_path / "pg.parquet"

    with MaxQuantPgAdapter() as adapter:
        with pytest.raises(ValueError, match="experiment_1.*mapping"):
            adapter.convert(str(protein_groups), str(output))

    assert not output.exists()


def test_pg_rejects_unassignable_total_intensity(tmp_path):
    """A total-only intensity must not create grouped_runs=['unknown']."""
    from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text("Protein IDs\tMajority protein IDs\tIntensity\nP12345\tP12345\t1000\n")
    output = tmp_path / "pg.parquet"

    with MaxQuantPgAdapter() as adapter:
        with pytest.raises(ValueError, match="total Intensity.*grouped_runs"):
            adapter.convert(str(protein_groups), str(output))

    assert not output.exists()


_EVIDENCE_HEADER = (
    "Sequence\tModified sequence\tCharge\tRaw file\tType\tMS/MS scan number\tm/z\t"
    "Calibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\t"
    "PEP\tLeading razor protein\tLeading proteins\tIntensity\tMass\tid\n"
)


def test_feature_dedup_collapses_identity_equal_rows(tmp_path):
    """Bug #248 (H): rows sharing the feature identity composite must collapse.

    The two ``run1`` rows share (peptidoform, charge, run, rt_start, rt_stop, scan)
    — i.e. the same ``feature_id`` — and differ only in ``Calibrated retention
    time``. The old dedup partitioned on calibrated RT (not part of the identity),
    so both survived as duplicate primary keys. They must now collapse to one,
    keeping the identified (MULTI-MSMS, lower PEP) row over the MBR one.
    """
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        _EVIDENCE_HEADER
        # Identified row and its MBR twin: same scan/rt-window, differ in calibrated RT.
        + "PEPTIDEK\t_PEPTIDEK_\t2\trun1\tMULTI-MSMS\t100\t450.25\t30.0\t29.5\t30.5\t0.01\tP1\tP1\t1000\t898.5\t1\n"
        + "PEPTIDEK\t_PEPTIDEK_\t2\trun1\tMULTI-MATCH\t100\t450.25\t30.4\t29.5\t30.5\t0.90\tP1\tP1\t1000\t898.5\t2\n"
        # A genuinely distinct feature (different charge/scan) must survive.
        + "PEPTIDEK\t_PEPTIDEK_\t3\trun1\tMULTI-MSMS\t200\t300.50\t40.0\t39.5\t40.5\t0.02\tP1\tP1\t2000\t898.5\t3\n"
    )
    output = tmp_path / "feature.parquet"

    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output), chunksize=1)

    rows = pq.read_table(output).to_pylist()
    assert len(rows) == 2, f"identity-equal rows must collapse, got {len(rows)}"
    # feature_id must be unique (no duplicate PK).
    assert len({row["feature_id"] for row in rows}) == 2
    kept = next(r for r in rows if r["charge"] == 2)
    assert kept["scan"] == [100]
    assert kept["posterior_error_probability"] == pytest.approx(0.01)
    assert kept["rt"] == pytest.approx(1800.0)
    assert kept["id_run_file_name"] == "run1"


def test_feature_dedup_does_not_degrade_without_type_column(tmp_path):
    """Bug #248 (H-b): a missing column must not degrade dedup to ``SELECT *``.

    With no ``Type`` column, the old guard fell back to ``SELECT * FROM evidence``,
    letting identity duplicates through. Dedup must still collapse identity-equal
    rows using the identity columns that are present.
    """
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    header = (
        "Sequence\tModified sequence\tCharge\tRaw file\tMS/MS scan number\tm/z\t"
        "Calibrated retention time\tCalibrated retention time start\t"
        "Calibrated retention time finish\tPEP\tLeading razor protein\tLeading proteins\tIntensity\tMass\tid\n"
    )
    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        header
        + "PEPTIDEK\t_PEPTIDEK_\t2\trun1\t100\t450.25\t30.0\t29.5\t30.5\t0.01\tP1\tP1\t1000\t898.5\t1\n"
        + "PEPTIDEK\t_PEPTIDEK_\t2\trun1\t100\t450.25\t30.4\t29.5\t30.5\t0.05\tP1\tP1\t1000\t898.5\t2\n"
    )
    output = tmp_path / "feature.parquet"

    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output), chunksize=100)

    rows = pq.read_table(output).to_pylist()
    assert len(rows) == 1, f"identity duplicates must collapse even without Type, got {len(rows)}"


def test_feature_without_ms1_quantity_has_a_null_intensity(tmp_path):
    """MaxQuant's blank or 0 ``Intensity`` (MS/MS-only evidence) is not a measured 0.0.

    The identified feature is kept, its label with a null intensity.
    """
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        _EVIDENCE_HEADER
        + "PEPTIDEA\t_PEPTIDEA_\t2\trun1\tMSMS\t100\t450.25\t30.0\t29.5\t30.5\t0.01\tP1\tP1\t\t898.5\t1\n"
        + "PEPTIDEB\t_PEPTIDEB_\t2\trun1\tMSMS\t200\t450.25\t40.0\t39.5\t40.5\t0.01\tP1\tP1\t0\t898.5\t2\n"
        + "PEPTIDEC\t_PEPTIDEC_\t2\trun1\tMULTI-MSMS\t300\t450.25\t50.0\t49.5\t50.5\t0.01\tP1\tP1\t1000\t898.5\t3\n"
    )
    output = tmp_path / "feature.parquet"
    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output))

    intensities = {r["sequence"]: r["intensities"] for r in pq.read_table(output).to_pylist()}
    assert intensities == {
        "PEPTIDEA": [{"label": "LFQ", "intensity": None}],
        "PEPTIDEB": [{"label": "LFQ", "intensity": None}],
        "PEPTIDEC": [{"label": "LFQ", "intensity": 1000.0}],
    }


def test_pg_zero_lfq_and_ibaq_are_missing(tmp_path):
    """MaxQuant writes 0 for an LFQ or iBAQ value it did not compute; that is null, not 0.0."""
    from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text(
        "Protein IDs\tMajority protein IDs\tIntensity exp1\tLFQ intensity exp1\tiBAQ exp1\n"
        "P1\tP1\t1000\t0\t0\n"
        "P2\tP2\t2000\t1500\t30\n"
    )
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text("source name\tcomment[data file]\tcomment[label]\nexp1\trun1.raw\tlabel free sample\n")
    output = tmp_path / "pg.parquet"
    with MaxQuantPgAdapter() as adapter:
        adapter.convert(str(protein_groups), str(output), sdrf_path=str(sdrf))

    rows = {r["anchor_protein"]: r for r in pq.read_table(output).to_pylist()}
    assert rows["P1"]["intensity"] == pytest.approx(1000.0)
    assert rows["P1"]["additional_intensities"] is None
    assert rows["P2"]["additional_intensities"] == [
        {
            "label": "LFQ",
            "intensities": [
                {"intensity_name": "lfq", "intensity_value": 1500.0},
                {"intensity_name": "ibaq", "intensity_value": 30.0},
            ],
        }
    ]


def test_pg_spectral_count_per_experiment(tmp_path):
    """``MS/MS count <experiment>`` is the pg row's ``spectral_count``; 0 is a count, not a missing value."""
    from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text(
        "Protein IDs\tMajority protein IDs\tMS/MS count\tIntensity exp1\tIntensity exp2\t"
        "MS/MS count exp1\tMS/MS count exp2\n"
        "P1\tP1\t7\t1000\t2000\t7\t0\n"
    )
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\nexp1\trun1.raw\tlabel free sample\nexp2\trun2.raw\tlabel free sample\n"
    )
    output = tmp_path / "pg.parquet"
    with MaxQuantPgAdapter() as adapter:
        adapter.convert(str(protein_groups), str(output), sdrf_path=str(sdrf))

    counts = {
        tuple(row["grouped_runs"]): {
            entry["intensity_name"]: entry["intensity_value"] for entry in row["additional_intensities"][0]["intensities"]
        }
        for row in pq.read_table(output).to_pylist()
    }
    assert counts == {("run1",): {"spectral_count": 7.0}, ("run2",): {"spectral_count": 0.0}}


def test_feature_spectral_count_from_evidence(tmp_path):
    """Evidence ``MS/MS count`` is the feature's ``spectral_count``, also without an MS1 quantity."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        _EVIDENCE_HEADER.replace("\tid\n", "\tid\tMS/MS count\n")
        + "PEPTIDEA\t_PEPTIDEA_\t2\trun1\tMULTI-MSMS\t100\t450.25\t30.0\t29.5\t30.5\t0.01\tP1\tP1\t1000\t898.5\t1\t3\n"
        + "PEPTIDEB\t_PEPTIDEB_\t2\trun1\tMULTI-MATCH\t\t450.25\t40.0\t39.5\t40.5\t\tP1\tP1\t2000\t898.5\t2\t0\n"
        + "PEPTIDEC\t_PEPTIDEC_\t2\trun1\tMSMS\t300\t450.25\t50.0\t49.5\t50.5\t0.01\tP1\tP1\t\t898.5\t3\t1\n"
    )
    output = tmp_path / "feature.parquet"
    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output))

    counts = {
        r["sequence"]: [
            (a["label"], i["intensity_name"], i["intensity_value"]) for a in r["additional_intensities"] for i in a["intensities"]
        ]
        for r in pq.read_table(output).to_pylist()
    }
    assert counts == {
        "PEPTIDEA": [("LFQ", "spectral_count", 3.0)],
        "PEPTIDEB": [("LFQ", "spectral_count", 0.0)],
        "PEPTIDEC": [("LFQ", "spectral_count", 1.0)],
    }


_MOD_SDRF_HEADER = (
    "source name\tcomment[data file]\tcomment[label]\tcomment[modification parameters]\t"
    "comment[modification parameters]\tcomment[modification parameters]\tcomment[modification parameters]\n"
)


def _write_mod_sdrf(tmp_path, *mods):
    """An LFQ SDRF for run1 declaring ``mods`` (four modification-parameter cells)."""
    cells = list(mods) + [""] * (4 - len(mods))
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text(_MOD_SDRF_HEADER + "\t".join(["exp1", "run1.raw", "label free sample", *cells]) + "\n")
    return sdrf


_CARBAMIDOMETHYL = "NT=Carbamidomethyl;TA=C;MT=fixed;AC=UNIMOD:4"


@pytest.mark.parametrize(
    ("modified_sequence", "fixed_mods", "expected"),
    [
        ("_PEPCKC_", (("UNIMOD:4", "C", ""),), "PEPC[UNIMOD:4]KC[UNIMOD:4]"),
        ("_PEPTIDEK_", (("UNIMOD:4", "C", ""),), "PEPTIDEK"),
        ("_PEPCM(Oxidation (M))K_", (("UNIMOD:4", "C", ""),), "PEPC[UNIMOD:4]M[UNIMOD:35]K"),
        # an N-terminal acetyl keeps a fixed N-term modification off
        ("_(Acetyl (Protein N-term))PEPK_", (("UNIMOD:5", "", "N"),), "[UNIMOD:1]-PEPK"),
        ("_PEPK_", (("UNIMOD:5", "", "N"),), "[UNIMOD:5]-PEPK"),
        ("_PEPK_", (("UNIMOD:2", "K", "C"),), "PEPK-[UNIMOD:2]"),
        ("_PEPR_", (("UNIMOD:2", "K", "C"),), "PEPR"),
    ],
)
def test_to_proforma_adds_fixed_modifications(modified_sequence, fixed_mods, expected):
    """MaxQuant leaves fixed modifications out of ``Modified sequence``; to_proforma puts them back."""
    from qpx.converters.maxquant.constants import to_proforma

    assert to_proforma(modified_sequence, fixed_mods) == expected


def test_fixed_modifications_read_from_sdrf(tmp_path, caplog):
    """Only placeable, non-isobaric MT=fixed SDRF modifications are taken, by their accession."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    sdrf = _write_mod_sdrf(
        tmp_path,
        "NT=Carbamidomethyl;TA=C;MT=Fixed;AC=Unimod:4",
        "NT=Oxidation;MT=variable;TA=M;AC=UNIMOD:35",
        "NT=TMT6plex;TA=K;MT=fixed;AC=UNIMOD:737",
        "NT=Acetyl;AC=UNIMOD:1;PP=Protein N-term;MT=fixed",
    )
    with MaxQuantFeatureAdapter() as adapter, caplog.at_level("WARNING"):
        assert adapter._fixed_mods(str(sdrf), str(tmp_path / "msms.txt")) == (("UNIMOD:4", "C", ""),)
    assert "Protein N-term" in caplog.text


def test_fixed_modification_without_accession_uses_its_unimod_name(tmp_path):
    """An SDRF fixed modification given only by ``NT`` takes the Unimod accession of that name.

    A name the modification database does not know cannot be placed and is
    left out, as before. So is one of OpenMS' own entries (``CUSTOM0``), which
    is not a Unimod modification.
    """
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    sdrf = _write_mod_sdrf(
        tmp_path,
        "NT=Carbamidomethyl;TA=C;MT=fixed",
        "NT=Methylthio;TA=C,M;MT=fixed",
        "NT=NoSuchMod;TA=K;MT=fixed",
        "NT=CUSTOM0;TA=K;MT=fixed",
    )
    with MaxQuantFeatureAdapter() as adapter:
        assert adapter._fixed_mods(str(sdrf), str(tmp_path / "msms.txt")) == (("UNIMOD:4", "C", ""), ("UNIMOD:39", "CM", ""))


@pytest.mark.parametrize(
    "label",
    [
        "AC=UNIMOD:737;TA=K;MT=fixed",
        "AC=UNIMOD:2016;TA=K;MT=fixed",
        "AC=UNIMOD:214;TA=K;MT=fixed",
        "NT=Label:13C(6)15N(4);TA=R;MT=fixed;AC=UNIMOD:267",
        "AC=UNIMOD:259;TA=K;MT=fixed",
    ],
    ids=["tmt6plex_accession", "tmtpro_accession", "itraq4plex_accession", "silac_name", "silac_accession"],
)
def test_quantification_labels_are_not_fixed_modifications(tmp_path, label):
    """TMT, iTRAQ and SILAC labels are skipped by name or, without ``NT``, by their Unimod accession.

    MaxQuant reports ``Mass`` and ``m/z`` without isobaric tags, and carries a
    SILAC label on the labelled rows only, so adding either to every
    peptidoform would break its mass.
    """
    from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter

    sdrf = _write_mod_sdrf(tmp_path, _CARBAMIDOMETHYL, label)
    with MaxQuantPsmAdapter() as adapter:
        assert adapter._fixed_mods(str(sdrf), str(tmp_path / "msms.txt")) == (("UNIMOD:4", "C", ""),)


def test_feature_and_psm_peptidoforms_carry_sdrf_fixed_modifications(tmp_path):
    """Fixed modifications from the SDRF make the peptidoform match calculated_mz."""
    from qpx.converters.maxquant.converter import MaxQuantConverter
    from qpx.converters.ptm import compute_precursor_mz
    from qpx.core.constants import FEATURE, PSM

    # PEPCK + carbamidomethyl (572.2628 + 57.0215): MaxQuant's Mass includes the fixed modification.
    mass = 629.2843
    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        _EVIDENCE_HEADER + f"PEPCK\t_PEPCK_\t2\trun1\tMULTI-MSMS\t100\t315.65\t30.0\t29.5\t30.5\t0.01\tP1\tP1\t1000\t{mass}\t1\n"
    )
    msms = tmp_path / "msms.txt"
    msms.write_text(
        "Sequence\tModified sequence\tCharge\tRaw file\tScan number\tm/z\tMass\tPEP\tReverse\tScore\tProteins\n"
        f"PEPCK\t_PEPCK_\t2\trun1\t100\t315.65\t{mass}\t0.01\t\t100\tP1\n"
    )
    out = tmp_path / "out"
    MaxQuantConverter().convert(
        output_folder=str(out),
        msms_file=str(msms),
        evidence_file=str(evidence),
        sdrf_file=str(_write_mod_sdrf(tmp_path, _CARBAMIDOMETHYL)),
        output_prefix="mq",
        structures=[PSM, FEATURE],
    )

    for view in ("feature", "psm"):
        (row,) = pq.read_table(out / f"mq.{view}.parquet").to_pylist()
        assert row["peptidoform"] == "PEPC[UNIMOD:4]K"
        assert row["modifications"] == [
            {
                "name": "Carbamidomethyl",
                "accession": "UNIMOD:4",
                "positions": [{"position": 4, "amino_acid": "C", "scores": None}],
            }
        ]
        assert compute_precursor_mz(row["peptidoform"], 2) == pytest.approx(row["calculated_mz"], rel=1e-6)


def test_peptidoforms_unchanged_without_sdrf_fixed_modifications(tmp_path, caplog):
    """Without fixed modifications in the SDRF the peptidoform stays as MaxQuant wrote it, with a warning."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        _EVIDENCE_HEADER + "PEPCK\t_PEPCK_\t2\trun1\tMULTI-MSMS\t100\t315.65\t30.0\t29.5\t30.5\t0.01\tP1\tP1\t1000\t629.2843\t1\n"
    )
    output = tmp_path / "feature.parquet"
    with MaxQuantFeatureAdapter() as adapter, caplog.at_level("WARNING"):
        adapter.convert(str(evidence), str(output), sdrf_path=str(_write_mod_sdrf(tmp_path)))

    (row,) = pq.read_table(output).to_pylist()
    assert row["peptidoform"] == "PEPCK"
    assert len([r for r in caplog.records if "fixed modification" in r.message]) == 1


def test_psm_handles_empty_modified_sequence(tmp_path):
    """Bug #248 (H): an empty ``Modified sequence`` must not abort the PSM view.

    The ``else None`` tuple-unpack raised ``TypeError`` for any row whose
    modified sequence was empty. The adapter must now skip that single row
    (peptidoform is a required schema field) instead of crashing the whole
    conversion, so the remaining PSMs still convert.
    """
    from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter

    msms = tmp_path / "msms.txt"
    msms.write_text(
        "Sequence\tModified sequence\tCharge\tRaw file\tScan number\tm/z\tMass\tPEP\tReverse\tScore\tProteins\n"
        "PEPTIDEK\t_PEPTIDEK_\t2\trun1\t100\t450.25\t898.5\t0.01\t\t100\tP1\n"
        # Empty modified sequence (renders to an empty peptidoform via to_proforma).
        "PEPTIDER\t_\t2\trun1\t101\t500.25\t998.5\t0.02\t\t90\tP2\n"
    )
    output = tmp_path / "psm.parquet"

    with MaxQuantPsmAdapter() as adapter:
        adapter.convert(str(msms), str(output))

    rows = pq.read_table(output).to_pylist()
    # The empty-modseq row is skipped; the valid PSM must still convert.
    assert len(rows) == 1, "valid PSM must survive; empty modseq row skipped, not crashing"
    assert rows[0]["sequence"] == "PEPTIDEK"
    assert rows[0]["peptidoform"]


def test_feature_pg_qvalue_and_genes_keyed_on_full_membership(tmp_path):
    """Bug #248 (M): a shared razor accession must not bind to another group's metadata.

    Both features share leading razor protein P1, but belong to different protein
    groups (P1;P2 vs P1;P3). The old leader-keyed map bound P1 to the first group
    by file order, so both features inherited that group's q-value/genes. The map
    must be keyed on the full protein-group membership.
    """
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        "Sequence\tModified sequence\tCharge\tRaw file\tType\tMS/MS scan number\tm/z\t"
        "Calibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\t"
        "PEP\tLeading razor protein\tLeading proteins\tGene names\tIntensity\tMass\tid\n"
        "PEPTIDEA\t_PEPTIDEA_\t2\trun1\tMULTI-MSMS\t100\t450.25\t30.0\t29.5\t30.5\t0.01\tP1\tP1;P2\t\t1000\t500\t1\n"
        "PEPTIDEB\t_PEPTIDEB_\t2\trun1\tMULTI-MSMS\t200\t460.25\t31.0\t30.5\t31.5\t0.02\tP1\tP1;P3\t\t2000\t520\t2\n"
    )
    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text(
        "Protein IDs\tMajority protein IDs\tQ-value\tGene names\nP1;P2\tP1;P2\t0.001\tGENEA\nP1;P3\tP1;P3\t0.5\tGENEB\n"
    )
    output = tmp_path / "feature.parquet"

    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output), protein_groups_path=str(protein_groups))

    rows = pq.read_table(output).to_pylist()
    by_seq = {r["sequence"]: r for r in rows}
    assert by_seq["PEPTIDEA"]["pg_global_qvalue"] == pytest.approx(0.001)
    assert by_seq["PEPTIDEA"]["gg_names"] == ["GENEA"]
    assert by_seq["PEPTIDEB"]["pg_global_qvalue"] == pytest.approx(0.5)
    assert by_seq["PEPTIDEB"]["gg_names"] == ["GENEB"]


def test_feature_warns_on_dropped_tmt_channels(caplog):
    """Bug #248 (M): channels with no reporter column must warn, not vanish silently."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    # evidence has only Reporter intensity 0-7; a 10-plex SDRF adds TMT130C (col 8)
    # and TMT131 (col 9), which have no column and would otherwise be dropped.
    columns = [f"Reporter intensity {i}" for i in range(8)]
    labels = ["TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131"]
    channels = [(label, index) for index, label in enumerate(labels)]
    with MaxQuantFeatureAdapter() as adapter:
        with caplog.at_level("WARNING"):
            adapter._warn_dropped_tmt_channels(columns, channels, ri_offset=0)
            # Second call with the same dropped set must not re-warn.
            adapter._warn_dropped_tmt_channels(columns, channels, ri_offset=0)

    warnings = [r for r in caplog.records if r.levelname == "WARNING" and "dropped" in r.message.lower()]
    assert len(warnings) == 1, "dropped-channel warning must be emitted exactly once"
    assert "TMT130C" in caplog.text
    assert "TMT131" in caplog.text


def _reporter_features(tmp_path, labels, n_columns):
    """Convert one evidence row whose ``Reporter intensity i`` is ``100 * (i + 1)``."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    reporter = "\t".join(f"Reporter intensity {i}" for i in range(n_columns))
    values = "\t".join(str(100 * (i + 1)) for i in range(n_columns))
    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        _EVIDENCE_HEADER.rstrip("\n") + f"\t{reporter}\n"
        f"PEPTIDEK\t_PEPTIDEK_\t2\trun1\tMULTI-MSMS\t10\t500.0\t10.0\t9.9\t10.1\t0.01\tP1\tP1\t1000\t998.0\t0\t{values}\n"
    )
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\n"
        + "".join(f"S{i}\trun1.raw\t{label}\n" for i, label in enumerate(labels))
    )
    output = tmp_path / "feature.parquet"
    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output), sdrf_path=str(sdrf))
    (row,) = pq.read_table(output).to_pylist()
    return {entry["label"]: entry["intensity"] for entry in row["intensities"]}


@pytest.mark.parametrize(
    ("labels", "n_columns", "expected"),
    [
        # TMT6: TMT131 is the sixth reporter, not TMT10's tenth.
        (["TMT131", "TMT126", "TMT129", "TMT127", "TMT130", "TMT128"], 6, {"TMT126": 100, "TMT127": 200, "TMT131": 600}),
        # iTRAQ channels follow the plex position whatever the SDRF row order.
        (["itraq117", "ITRAQ114", "iTRAQ116", "ITRAQ115"], 4, {"ITRAQ114": 100, "ITRAQ115": 200, "ITRAQ117": 400}),
        (["ITRAQ121", "ITRAQ113", "ITRAQ119"], 8, {"ITRAQ113": 100, "ITRAQ119": 700, "ITRAQ121": 800}),
        # A partial TMTpro-18 plex resolves to its 18-channel positions.
        (["TMT126", "TMT134C", "TMT135N"], 18, {"TMT126": 100, "TMT134C": 1700, "TMT135N": 1800}),
    ],
)
def test_feature_reads_reporter_columns_by_plex_position(tmp_path, labels, n_columns, expected):
    """Each SDRF label reads the MaxQuant reporter column of its plex position."""
    intensities = _reporter_features(tmp_path, labels, n_columns)
    assert len(intensities) == len(labels)
    for label, value in expected.items():
        assert intensities[label] == value


@pytest.mark.parametrize("spelling", ["ITRAQ{}", "iTRAQ reagent {}", "iTRAQ4plex-{}"])
def test_itraq4_labels_read_the_columns_of_an_8_column_table(tmp_path, spelling):
    """ITRAQ114 is the first channel of iTRAQ4 but the second of the 8-plex the table has."""
    labels = [spelling.format(channel) for channel in (114, 115, 116, 117)]
    assert _reporter_features(tmp_path, labels, 8) == {
        "ITRAQ114": 200.0,
        "ITRAQ115": 300.0,
        "ITRAQ116": 400.0,
        "ITRAQ117": 500.0,
    }
    assert _reporter_features(tmp_path, labels, 4)["ITRAQ114"] == 100.0


def test_lfq_feature_spectral_count_and_pepmap_from_evidence(tmp_path):
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        _EVIDENCE_HEADER.rstrip("\n") + "\tMS/MS count\tProteins\tReverse\n"
        "PEPTIDEK\t_PEPTIDEK_\t2\trun1\tMULTI-MSMS\t10\t500.0\t10.0\t9.9\t10.1\t0.01\tP1\tP1\t1000\t998.0\t0\t3\tsp|P1|A_HUMAN;P2\t\n"
        "PEPTIDER\t_PEPTIDER_\t2\trun1\tMULTI-MSMS\t11\t510.0\t11.0\t10.9\t11.1\t0.01\tP3\tP3\t1000\t998.0\t0\t0\tP3\t\n"
        "REDITPEK\t_REDITPEK_\t2\trun1\tMULTI-MSMS\t12\t520.0\t12.0\t11.9\t12.1\t0.01\tREV__P1\tREV__P1\t1000\t998.0\t0\t1\tREV__P1\t+\n"
    )
    output, pepmap = tmp_path / "feature.parquet", tmp_path / "pepmap.parquet"
    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output))
        assert adapter.write_pepmap(str(pepmap))
    counts = {row["sequence"]: row["additional_intensities"][0]["intensities"][0] for row in pq.read_table(output).to_pylist()}
    assert counts["PEPTIDEK"] == {"intensity_name": "spectral_count", "intensity_value": 3.0}
    assert counts["PEPTIDER"]["intensity_value"] == 0.0
    mapped = {row["sequence"]: row for row in pq.read_table(pepmap).to_pylist()}
    assert set(mapped) == {"PEPTIDEK", "PEPTIDER"}
    assert [entry["accession"] for entry in mapped["PEPTIDEK"]["pg_accessions"]] == ["P1", "P2"]
    assert mapped["PEPTIDEK"]["is_unique"] is False
    assert mapped["PEPTIDER"]["is_unique"] is True


def test_feature_rejects_labels_outside_any_plex(tmp_path):
    with pytest.raises(ValueError, match="match no TMT plex"):
        _reporter_features(tmp_path, ["TMT126", "TMT10plex-127"], 10)


def test_pg_warns_on_dropped_tmt_channels(caplog):
    """Bug #248 (M): the pg adapter must also warn about dropped channels."""
    from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

    with MaxQuantPgAdapter() as adapter:
        with caplog.at_level("WARNING"):
            adapter._warn_dropped_tmt_channels("exp1", ["TMT130C", "TMT131"])
            adapter._warn_dropped_tmt_channels("exp1", ["TMT130C", "TMT131"])

    warnings = [r for r in caplog.records if r.levelname == "WARNING" and "dropped" in r.message.lower()]
    assert len(warnings) == 1
    assert "TMT130C" in caplog.text and "TMT131" in caplog.text


# ---------------------------------------------------------------------------
# Module-scoped fixture: run conversion once, share outputs across tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def converted_output(tmp_path_factory):
    """Run MaxQuant conversion once for all tests in this module.

    Returns the output directory containing the generated Parquet files.
    If the test data is missing, all tests in this module are skipped.
    """
    from qpx.converters.maxquant.converter import MaxQuantConverter

    if not _MSMS.exists() or not _EVIDENCE.exists():
        pytest.skip(f"Test data not found: {EXAMPLES_DIR}")

    output = tmp_path_factory.mktemp("maxquant_output")

    converter = MaxQuantConverter()
    converter.convert(
        output_folder=output,
        msms_file=_MSMS,
        evidence_file=_EVIDENCE,
        sdrf_file=_SDRF if _SDRF.exists() else None,
        output_prefix=_PREFIX,
    )

    return output


@pytest.fixture(scope="module")
def psm_table(converted_output):
    """Read the psm.parquet produced by the converter."""
    path = converted_output / f"{_PREFIX}.psm.parquet"
    assert path.exists(), "psm.parquet was not produced"
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def feature_table(converted_output):
    """Read the feature.parquet produced by the converter."""
    path = converted_output / f"{_PREFIX}.feature.parquet"
    assert path.exists(), "feature.parquet was not produced"
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# PSM conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantPsmConversion:
    """Validate psm.parquet output from MaxQuant conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.psm.parquet"
        if not path.exists():
            raise AssertionError("psm.parquet was not created")

    def test_has_rows(self, psm_table):
        if psm_table.num_rows == 0:
            raise AssertionError("psm.parquet is empty")

    def test_key_columns_present(self, psm_table):
        column_names = set(psm_table.column_names)
        expected = {"sequence", "charge", "calculated_mz", "run_file_name", "scan"}
        missing = expected - column_names
        if missing:
            raise AssertionError(f"Missing columns: {missing}")

    def test_sequence_values_are_nonempty(self, psm_table):
        sequences = psm_table.column("sequence").to_pylist()
        for seq in sequences:
            if not (isinstance(seq, str) and len(seq) > 0):
                raise AssertionError(f"Expected non-empty string, got {seq!r}")

    def test_charge_values_are_valid(self, psm_table):
        charges = psm_table.column("charge").to_pylist()
        for charge in charges:
            if not 1 <= charge <= 10:
                raise AssertionError(f"Charge out of range: {charge}")

    def test_calculated_mz_values_are_nonnegative(self, psm_table):
        mz_values = psm_table.column("calculated_mz").to_pylist()
        for mz in mz_values:
            if mz is None or mz < 0:
                raise AssertionError(f"Invalid calculated_mz value: {mz}")

    def test_run_file_names_are_nonempty(self, psm_table):
        run_names = psm_table.column("run_file_name").to_pylist()
        for name in run_names:
            if not (isinstance(name, str) and len(name) > 0):
                raise AssertionError(f"Expected non-empty run name, got {name!r}")

    def test_schema_validation(self, psm_table):
        """Verify that the output conforms to the PsmSchema."""
        from qpx.core.data import PsmSchema

        errors = PsmSchema.validate(psm_table)
        if errors:
            raise AssertionError(f"Schema validation errors: {errors}")


# ---------------------------------------------------------------------------
# Feature conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantFeatureConversion:
    """Validate feature.parquet output from MaxQuant conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.feature.parquet"
        if not path.exists():
            raise AssertionError("feature.parquet was not created")

    def test_has_rows(self, feature_table):
        if feature_table.num_rows == 0:
            raise AssertionError("feature.parquet is empty")

    def test_key_columns_present(self, feature_table):
        column_names = set(feature_table.column_names)
        expected = {
            "sequence",
            "charge",
            "calculated_mz",
            "intensities",
            "run_file_name",
        }
        missing = expected - column_names
        if missing:
            raise AssertionError(f"Missing columns: {missing}")

    def test_sequence_values_are_nonempty(self, feature_table):
        sequences = feature_table.column("sequence").to_pylist()
        for seq in sequences:
            if not (isinstance(seq, str) and len(seq) > 0):
                raise AssertionError(f"Expected non-empty string, got {seq!r}")

    def test_charge_values_are_valid(self, feature_table):
        charges = feature_table.column("charge").to_pylist()
        for charge in charges:
            if not 1 <= charge <= 10:
                raise AssertionError(f"Charge out of range: {charge}")

    def test_intensities_are_nonnegative(self, feature_table):
        rows = feature_table.column("intensities").to_pylist()
        for row_intensities in rows:
            if row_intensities is None:
                continue
            for entry in row_intensities:
                if entry["intensity"] is not None and entry["intensity"] < 0:
                    raise AssertionError(f"Negative intensity: {entry['intensity']}")

    def test_run_file_names_are_nonempty(self, feature_table):
        run_names = feature_table.column("run_file_name").to_pylist()
        for name in run_names:
            if not (isinstance(name, str) and len(name) > 0):
                raise AssertionError(f"Expected non-empty run name, got {name!r}")

    def test_schema_validation(self, feature_table):
        """Verify that the output conforms to the FeatureSchema."""
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        if errors:
            raise AssertionError(f"Schema validation errors: {errors}")

    def test_tmt_nc_isomer_channels_have_distinct_intensities(self, feature_table):
        """Regression test: TMT127N and TMT127C must not systematically share the
        same intensity value.

        The alphabetical sort + i+1 fallback bug in _build_intensities caused
        TMT127C to read from the TMT127N MaxQuant column (and vice-versa),
        producing duplicate intensities for 96 %+ of features in real datasets.
        """
        rows = feature_table.column("intensities").to_pylist()

        equal_pairs = 0
        total_pairs = 0

        for row_intensities in rows:
            if not row_intensities:
                continue
            by_label = {e["label"]: e["intensity"] for e in row_intensities if e["intensity"]}
            n_val = by_label.get("TMT127N")
            c_val = by_label.get("TMT127C")
            if n_val is not None and c_val is not None:
                total_pairs += 1
                if abs(n_val - c_val) < 1e-4:
                    equal_pairs += 1

        if total_pairs == 0:
            pytest.skip("No features with both TMT127N and TMT127C detected")

        duplicate_rate = equal_pairs / total_pairs
        if duplicate_rate >= 0.10:
            raise AssertionError(
                f"TMT127N == TMT127C in {equal_pairs}/{total_pairs} features "
                f"({100 * duplicate_rate:.1f}%) — N/C isomer channel assignment bug detected"
            )


# ---------------------------------------------------------------------------
# Ontology conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantOntologyConversion:
    """Validate ontology.parquet output from MaxQuant conversion.

    The ontology file is only written when scores are discovered during
    conversion. It may legitimately not exist if no scores were tracked.
    """

    def test_ontology_file_exists_or_no_scores(self, converted_output):
        path = converted_output / f"{_PREFIX}.ontology.parquet"
        assert path.exists(), "ontology.parquet not written (no scores discovered)"

        table = pq.read_table(str(path))
        if table.num_rows == 0:
            raise AssertionError("ontology.parquet exists but is empty")

    def test_ontology_columns(self, converted_output):
        path = converted_output / f"{_PREFIX}.ontology.parquet"
        assert path.exists(), "ontology.parquet not written (no scores discovered)"

        table = pq.read_table(str(path))
        column_names = set(table.column_names)
        expected = {"field_name", "view"}
        missing = expected - column_names
        if missing:
            raise AssertionError(f"Missing ontology columns: {missing}")

    def test_schema_validation(self, converted_output):
        from qpx.core.data import OntologySchema

        path = converted_output / f"{_PREFIX}.ontology.parquet"
        assert path.exists(), "ontology.parquet not written (no scores discovered)"
        table = pq.read_table(str(path))
        errors = OntologySchema.validate(table)
        if errors:
            raise AssertionError(f"Schema validation errors: {errors}")


def test_every_mod_probabilities_column_reaches_psm_and_feature(tmp_path):
    """Oxidation (M) Probabilities in msms.txt and evidence.txt, not only Phospho (STY) in msms.txt."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter
    from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter

    psm_path = tmp_path / "t.psm.parquet"
    feature_path = tmp_path / "t.feature.parquet"
    with MaxQuantPsmAdapter() as adapter:
        adapter.convert(str(EXAMPLES_DIR / "msms.txt"), str(psm_path))
    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(EXAMPLES_DIR / "evidence.txt"), str(feature_path))

    for path in (psm_path, feature_path):
        rows = pq.read_table(path, columns=["peptidoform", "modifications", "modification_candidates"]).to_pylist()
        localized = [row for row in rows if row["modification_candidates"]]
        assert localized, path.name
        for row in localized:
            assigned = [c for c in row["modification_candidates"] if c["assigned"]]
            assert all(c["name"] == "Oxidation" and c["accession"] == "UNIMOD:35" for c in row["modification_candidates"])
            assert assigned and "M[UNIMOD:35]" in row["peptidoform"]
        placed = next(
            pos
            for row in localized
            for mod in row["modifications"]
            if mod["accession"] == "UNIMOD:35"
            for pos in mod["positions"]
        )
        assert placed["scores"][0]["score_name"] == "localization_probability"


def test_msms_matches_annotate_the_peaks_of_spectral_arrays(tmp_path):
    from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter

    output = tmp_path / "x.psm.parquet"
    with MaxQuantPsmAdapter() as adapter:
        adapter.convert(str(_MSMS), str(output), spectral_data=True)

    table = pq.read_table(output, columns=["mz_array", "ion_type_array", "charge_array"]).to_pylist()
    annotated = [row for row in table if row["ion_type_array"]]
    assert annotated
    assert all(len(row["ion_type_array"]) == len(row["charge_array"]) == len(row["mz_array"]) for row in annotated)
    assert pq.read_metadata(output).metadata[b"spectrum_content"] == b"matched_peaks"


def test_peak_annotations_split_charge_from_ion_type():
    from qpx.converters.maxquant.evidence_fields import peak_annotations

    assert peak_annotations("y3;b4(2+);y5-H2O;a2-NH3(3+)", 4) == (["y3", "b4", "y5-H2O", "a2-NH3"], [1, 2, 1, 3])
    # Tokens that do not line up with the peaks cannot annotate them.
    assert peak_annotations("y3;b4", 3) == (None, None)
    assert peak_annotations(None, 0) == (None, None)


def test_collisional_cross_section_is_read_when_present(tmp_path):
    source = Path(__file__).parents[1] / "examples" / "maxquant" / "maxquant_simple"
    evidence = pd.read_csv(source / "evidence.txt", sep="\t", nrows=5, low_memory=False)
    evidence["CCS"] = 412.5
    path = tmp_path / "evidence.txt"
    evidence.to_csv(path, sep="\t", index=False)
    output = tmp_path / "feature.parquet"

    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(path), str(output))

    assert set(pq.read_table(output, columns=["ccs"]).column("ccs").to_pylist()) == {412.5}
