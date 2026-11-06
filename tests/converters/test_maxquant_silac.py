"""MaxQuant SILAC conversion: one feature per multiplet, labelled PSMs, per-label protein groups."""

import pyarrow.parquet as pq
import pytest


def _write_silac_inputs(tmp_path):
    """A SILAC pair search: Lys8/Arg10 heavy, carbamidomethyl C fixed only in summary.txt.

    PEPTIDEK was sequenced in both partners (MULTI-MSMS, states in msms.txt only);
    SAMCLER only heavy, its light partner unquantified; PEPTIDE has no labelled
    residue, so MaxQuant cannot split it.
    """
    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        "Sequence\tModified sequence\tCharge\tRaw file\tExperiment\tType\tMS/MS scan number\tm/z\tMass\tMass error [ppm]\t"
        "Calibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\tPEP\t"
        "Proteins\tLeading razor protein\tLeading proteins\tProtein group IDs\tLabeling State\t"
        "Intensity\tIntensity L\tIntensity H\tRatio H/L\tRatio H/L normalized\tRatio H/L shift\tid\n"
        "PEPTIDEK\t_PEPTIDEK_\t2\trun1\texp1\tMULTI-MSMS\t100\t464.734729\t927.454904\t0\t30.0\t29.5\t30.5\t0.01\t"
        "P1\tP1\tP1\t0\t\t150\t100\t50\t0.5\t0.55\t0\t0\n"
        "SAMCLER\t_SAMCLER_\t2\trun1\texp1\tISO-MSMS\t200\t418.199\t834.384\t0\t31.0\t30.5\t31.5\t0.02\t"
        "P1\tP1\tP1\t0\t1\t30\t0\t30\t\t\t\t1\n"
        "PEPTIDE\t_PEPTIDE_\t1\trun1\texp1\tISO-MSMS\t300\t800.36\t799.36\t0\t32.0\t31.5\t32.5\t0.03\t"
        "P1\tP1\tP1\t0\t-1\t\t\t\t\t\t\t2\n"
    )
    msms = tmp_path / "msms.txt"
    msms.write_text(
        "Sequence\tModified sequence\tCharge\tRaw file\tScan number\tm/z\tMass\tMass error [ppm]\tRetention time\t"
        "PEP\tScore\tProteins\tLabeling state\tEvidence ID\n"
        "PEPTIDEK\t_PEPTIDEK_\t2\trun1\t100\t464.734729\t927.454904\t0\t30.0\t0.01\t100\tP1\t0\t0\n"
        "PEPTIDEK\t_PEPTIDEK_\t2\trun1\t101\t468.741828\t935.469103\t0\t30.1\t0.01\t90\tP1\t1\t0\n"
        "SAMCLER\t_SAMCLER_\t2\trun1\t200\t423.2\t844.39\t0\t31.0\t0.02\t80\tP1\t1\t1\n"
        "PEPTIDE\t_PEPTIDE_\t1\trun1\t300\t800.36\t799.36\t0\t32.0\t0.03\t70\tP1\t-1\t2\n"
    )
    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text(
        "Protein IDs\tMajority protein IDs\tQ-value\tIntensity\tIntensity L\tIntensity H\tIntensity exp1\t"
        "Intensity L exp1\tIntensity H exp1\tRatio H/L\tRatio H/L exp1\tRatio H/L normalized exp1\t"
        "Ratio H/L count exp1\tPeptides exp1\tUnique peptides exp1\tid\n"
        "P1\tP1\t0.001\t180\t100\t80\t180\t100\t80\t0.5\t0.5\t0.55\t2\t3\t3\t0\n"
    )
    (tmp_path / "summary.txt").write_text(
        "Raw file\tExperiment\tFixed modifications\tMultiplicity\tLabels0\tLabels1\n"
        "run1\texp1\tCarbamidomethyl (C)\t2\t\tArg10;Lys8\n"
        "Total\t\t\t\t\t\n"
    )
    return evidence, msms, protein_groups


@pytest.fixture(scope="module", name="silac_output")
def _silac_output(tmp_path_factory):
    """The SILAC inputs converted to every view, as ``{view: rows}``."""
    from qpx.converters.maxquant.converter import MaxQuantConverter
    from qpx.core.constants import FEATURE, PG, PSM

    tmp_path = tmp_path_factory.mktemp("silac")
    evidence, msms, protein_groups = _write_silac_inputs(tmp_path)
    out = tmp_path / "out"
    MaxQuantConverter().convert(
        output_folder=str(out),
        msms_file=str(msms),
        evidence_file=str(evidence),
        protein_groups_file=str(protein_groups),
        output_prefix="mq",
        structures=[PSM, FEATURE, PG],
    )
    return {view: pq.read_table(out / f"mq.{view}.parquet").to_pylist() for view in ("feature", "psm", "pg")}


_SILAC_RATIOS = [
    {
        "label": "SILAC heavy",
        "intensities": [
            {"intensity_name": "ratio_h_l", "intensity_value": 0.5},
            {"intensity_name": "ratio_h_l_normalized", "intensity_value": pytest.approx(0.55)},
        ],
    }
]


def test_silac_feature_has_one_row_per_multiplet(silac_output):
    """Per-label intensities, ratios on the heavy label, and the partners identified by MS/MS."""
    features = {row["sequence"]: row for row in silac_output["feature"]}
    pair = features["PEPTIDEK"]
    assert pair["peptidoform"] == "PEPTIDEK"
    assert pair["intensities"] == [{"label": "SILAC light", "intensity": 100.0}, {"label": "SILAC heavy", "intensity": 50.0}]
    assert pair["additional_intensities"] == _SILAC_RATIOS
    assert pair["identified_labels"] == ["SILAC light", "SILAC heavy"]
    assert features["SAMCLER"]["peptidoform"] == "SAMC[UNIMOD:4]LER"
    assert features["SAMCLER"]["intensities"][0] == {"label": "SILAC light", "intensity": None}
    assert features["SAMCLER"]["identified_labels"] == ["SILAC heavy"]
    unsplit = features["PEPTIDE"]
    assert [entry["intensity"] for entry in unsplit["intensities"]] == [None, None]
    assert unsplit["identified_labels"] is None


def test_silac_psm_carries_its_partner_label(silac_output):
    """A heavy PSM has the heavy label modifications, the matching m/z and its multiplet's feature."""
    from qpx.converters.ptm import compute_precursor_mz

    feature_ids = {row["sequence"]: row["feature_id"] for row in silac_output["feature"]}
    psms = {(row["sequence"], row["scan"][0]): row for row in silac_output["psm"]}
    light, heavy = psms[("PEPTIDEK", 100)], psms[("PEPTIDEK", 101)]
    assert (light["label"], light["peptidoform"]) == ("SILAC light", "PEPTIDEK")
    assert (heavy["label"], heavy["peptidoform"]) == ("SILAC heavy", "PEPTIDEK[UNIMOD:259]")
    assert heavy["calculated_mz"] == pytest.approx(468.741828, abs=1e-4)
    assert abs(heavy["mass_error_ppm"]) < 0.5
    assert heavy["feature_id"] == light["feature_id"] == feature_ids["PEPTIDEK"]
    samcler = psms[("SAMCLER", 200)]
    assert samcler["peptidoform"] == "SAMC[UNIMOD:4]LER[UNIMOD:267]"
    assert samcler["calculated_mz"] == pytest.approx(compute_precursor_mz(samcler["peptidoform"], 2), abs=1e-4)
    assert (psms[("PEPTIDE", 300)]["label"], psms[("PEPTIDE", 300)]["peptidoform"]) == (None, "PEPTIDE")


def test_silac_pg_has_one_row_per_label(silac_output):
    """One pg row per experiment and label, the ratios on the heavy row."""
    pg = {row["label"]: row for row in silac_output["pg"]}
    assert {label: row["intensity"] for label, row in pg.items()} == {"SILAC light": 100.0, "SILAC heavy": 80.0}
    assert pg["SILAC light"]["additional_intensities"] is None
    assert pg["SILAC heavy"]["additional_intensities"] == _SILAC_RATIOS


def test_silac_psms_need_the_label_definitions(tmp_path):
    """Without summary.txt the heavy label modifications are unknown, so conversion stops."""
    from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter

    _, msms, _ = _write_silac_inputs(tmp_path)
    (tmp_path / "summary.txt").unlink()
    with MaxQuantPsmAdapter() as adapter:
        with pytest.raises(ValueError, match="summary.txt"):
            adapter.convert(str(msms), str(tmp_path / "psm.parquet"))


def test_sdrf_silac_label_is_not_a_fixed_modification(tmp_path):
    """A SILAC label the SDRF declares fixed for the heavy samples stays on the heavy PSMs only."""
    from qpx.converters.maxquant.converter import MaxQuantConverter
    from qpx.core.constants import PSM

    _, msms, _ = _write_silac_inputs(tmp_path)
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\tcomment[modification parameters]\n"
        "S_light\trun1.raw\tSILAC light\t\n"
        "S_heavy\trun1.raw\tSILAC heavy\tAC=UNIMOD:259;TA=K;MT=fixed\n"
    )
    out = tmp_path / "out"
    MaxQuantConverter().convert(
        output_folder=str(out), msms_file=str(msms), sdrf_file=str(sdrf), output_prefix="mq", structures=[PSM]
    )
    psms = {(row["sequence"], row["scan"][0]): row["peptidoform"] for row in pq.read_table(out / "mq.psm.parquet").to_pylist()}
    assert (psms[("PEPTIDEK", 100)], psms[("PEPTIDEK", 101)]) == ("PEPTIDEK", "PEPTIDEK[UNIMOD:259]")


def test_silac_output_rejects_a_label_free_sdrf(tmp_path):
    """SILAC multiplets are never summed into one label-free intensity."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence, _, _ = _write_silac_inputs(tmp_path)
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text("source name\tcomment[data file]\tcomment[label]\nS1\trun1.raw\tlabel free sample\n")
    with MaxQuantFeatureAdapter() as adapter:
        with pytest.raises(ValueError, match="SILAC multiplets.*LFQ"):
            adapter.convert(str(evidence), str(tmp_path / "feature.parquet"), sdrf_path=str(sdrf))


@pytest.mark.parametrize(
    ("letters", "run_labels", "experiment_type", "expected"),
    [
        ("LH", {}, "LFQ", {"L": "SILAC light", "H": "SILAC heavy"}),
        ("LMH", {}, "LFQ", {"L": "SILAC light", "M": "SILAC medium", "H": "SILAC heavy"}),
        # A light/medium pair is still MaxQuant's "L"/"H".
        ("LH", {"r1": {"SILAC light"}, "r2": {"SILAC medium"}}, "SILAC", {"L": "SILAC light", "H": "SILAC medium"}),
        ("", {}, "LFQ", {}),
    ],
)
def test_silac_letters_take_the_sdrf_labels_in_plex_order(letters, run_labels, experiment_type, expected):
    """MaxQuant's L/M/H letters take the SDRF SILAC labels, lightest first."""
    from qpx.converters.maxquant.silac import multiplet_labels

    assert multiplet_labels(letters, run_labels, experiment_type) == expected


def test_silac_sdrf_without_multiplets_is_rejected():
    """An SDRF that disagrees with the output's multiplets is an error."""
    from qpx.converters.maxquant.silac import multiplet_labels

    with pytest.raises(ValueError, match="no SILAC multiplets"):
        multiplet_labels("", {"r1": {"SILAC light", "SILAC heavy"}}, "SILAC")
    with pytest.raises(ValueError, match="SILAC multiplets L/M/H"):
        multiplet_labels("LMH", {"r1": {"SILAC light", "SILAC heavy"}}, "SILAC")
