"""consensusXML channels, units and protein abundances on OpenMS' own TOPP test outputs."""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace

import pyarrow.parquet as pq
import pytest

from qpx.converters.openms_consensus.channels import column_labels
from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter
from qpx.converters.openms_consensus.streaming import _parse_hit_score

TOPP = Path(__file__).parents[1] / "examples" / "openms_topp"


def _convert(tmp_path, name, streaming):
    out = tmp_path / ("stream" if streaming else "memory")
    OpenMSConsensusConverter().convert(
        str(TOPP / f"{name}.consensusXML"), str(out), output_prefix="t", structures=("feature", "psm", "pg"), streaming=streaming
    )
    return {view: pq.read_table(out / f"t.{view}.parquet").to_pylist() for view in ("feature", "psm", "pg")}


def _pg_intensity(rows):
    return {(row["anchor_protein"], tuple(row["grouped_runs"]), row["label"]): row["intensity"] for row in rows}


def _extra(row, label=None):
    return {
        item["intensity_name"]: item["intensity_value"]
        for entry in row["additional_intensities"] or ()
        if label is None or entry["label"] == label
        for item in entry["intensities"]
    }


@pytest.fixture(params=[False, True], ids=["memory", "streaming"])
def streaming(request):
    return request.param


def test_silac_keeps_one_intensity_per_channel(tmp_path, streaming):
    views = _convert(tmp_path, "MS1LabeledWorkflow_6_ratios_out", streaming)
    feature = next(row for row in views["feature"] if row["sequence"] == "EYNVTIDAAESTVNK")
    assert feature["peptidoform"] == "EYNVTIDAAESTVNK"
    intensities = {entry["label"]: entry["intensity"] for entry in feature["intensities"]}
    assert intensities == {"SILAC light": 13987197.0, "SILAC heavy": 52794968.0}
    # Position, identification and mass error all describe the light partner.
    assert feature["observed_mz"] == pytest.approx(827.402, abs=1e-3)
    assert abs(feature["mass_error_ppm"]) < 1
    assert feature["identified_labels"] == ["SILAC light"]
    assert feature["additional_intensities"][0]["label"] == "SILAC heavy"
    assert _extra(feature)["ratio_h_l"] == pytest.approx(3.7745, abs=1e-3)


def test_silac_psm_is_the_matched_partner(tmp_path, streaming):
    views = _convert(tmp_path, "MS1LabeledWorkflow_6_ratios_out", streaming)
    psms = {row["peptidoform"]: row for row in views["psm"]}
    # Light and heavy hits of one peptide stay two PSMs.
    heavy = psms["TALGYYWALTEASGK[UNIMOD:259]"]
    light = psms["TALGYYWALTEASGK"]
    assert (heavy["label"], light["label"]) == ("SILAC heavy", "SILAC light")
    assert heavy["sequence"] == "TALGYYWALTEASGK"
    assert abs(heavy["mass_error_ppm"]) < 5
    assert {param["cv_name"]: param["cv_value"] for param in heavy["cv_params"]} == {
        "MS1Label:labeled_sequence": "TALGYYWALTEASGK(Label:13C(6)15N(2))",
        "MS1Label:removed_labels": "Lys8",
        "MS1Label:channel": "2",
    }


def test_silac_pg_is_openms_abundance_per_channel(tmp_path, streaming):
    views = _convert(tmp_path, "MS1LabeledWorkflow_6_ratios_out", streaming)
    unit = ("FeatureFinderMultiplex_2_input",)
    intensity = _pg_intensity(views["pg"])
    assert intensity["PROT_A", unit, "SILAC light"] == pytest.approx(1.4402612e07)
    assert intensity["PROT_A", unit, "SILAC heavy"] == pytest.approx(5.2800048e07)
    assert not any(row["label"] == "LFQ" for row in views["pg"])
    heavy = next(row for row in views["pg"] if row["anchor_protein"] == "PROT_A" and row["label"] == "SILAC heavy")
    light = next(row for row in views["pg"] if row["anchor_protein"] == "PROT_A" and row["label"] == "SILAC light")
    assert _extra(heavy)["ratio_h_l"] == pytest.approx(1.89337432)
    assert not any(name.startswith("ratio") for name in _extra(light))


def test_three_plex_with_partial_multiplets(tmp_path, streaming):
    views = _convert(tmp_path, "MS1LabeledWorkflow_7_partial_out", streaming)
    features = {row["sequence"]: row for row in views["feature"]}
    # The light partner is missing: position is the medium's, and no mass error
    # is computed against the unlabelled calculated m/z.
    partial = features["SAYWNEGHIR"]
    assert [entry["label"] for entry in partial["intensities"]] == ["SILAC medium", "SILAC heavy"]
    assert partial["observed_mz"] == pytest.approx(619.804, abs=1e-3)
    assert partial["mass_error_ppm"] is None
    assert features["ENGTWQFYIK"]["mass_error_ppm"] is not None


def test_isobaric_pg_is_openms_abundance_not_the_sum(tmp_path, streaming):
    views = _convert(tmp_path, "IsobaricWorkflow_1_out", streaming)
    row = next(row for row in views["pg"] if row["anchor_protein"] == "PROTA" and row["label"] == "TMT126")
    assert row["intensity"] == pytest.approx(1.58439697265625e04)
    assert _extra(row)["unnormalized_unique_peptide_sum"] == pytest.approx(2 * 1.58439697265625e04)
    assert {param["cv_value"] for param in row["cv_params"]} == {"openms_protein_abundance"}


def test_fraction_groups_become_units_and_zero_abundance_is_null(tmp_path, streaming):
    views = _convert(tmp_path, "ProteomicsLFQ_1_out", streaming)
    intensity = _pg_intensity(views["pg"])
    assert {unit for _, unit, _ in intensity} == {("BSA1_F1", "BSA1_F2"), ("BSA2_F1", "BSA2_F2"), ("BSA3_F1", "BSA3_F2")}
    albumin = "P02769|ALBU_BOVIN"
    assert intensity[albumin, ("BSA1_F1", "BSA1_F2"), "LFQ"] == pytest.approx(1.80143968e08)
    assert intensity[albumin, ("BSA3_F1", "BSA3_F2"), "LFQ"] == pytest.approx(7.4668e07)
    assert intensity["P00761|TRYP_PIG", ("BSA2_F1", "BSA2_F2"), "LFQ"] is None


def test_runs_without_a_fraction_group_stay_separate(tmp_path, caplog, streaming):
    source = (TOPP / "ProteomicsLFQ_1_out.consensusXML").read_text(encoding="iso-8859-1")
    lines = [line for line in source.splitlines() if 'name="fraction_group"' not in line and "_level_" not in line]
    path = tmp_path / "plain.consensusXML"
    path.write_text("\n".join(lines), encoding="iso-8859-1")
    with caplog.at_level(logging.WARNING):
        OpenMSConsensusConverter().convert(
            str(path), str(tmp_path / "out"), output_prefix="t", structures=("pg",), streaming=streaming
        )
    rows = pq.read_table(tmp_path / "out" / "t.pg.parquet").to_pylist()
    assert {tuple(row["grouped_runs"]) for row in rows} == {(f"BSA{s}_F{f}",) for s in (1, 2, 3) for f in (1, 2)}
    assert {row["cv_params"][0]["cv_value"] for row in rows if row["cv_params"]} == {"unnormalized_unique_peptide_sum"}
    assert "its own pg quantification unit" in caplog.text


def _map(*labels, filename="run.mzML"):
    headers = {idx: SimpleNamespace(filename=filename, label=label) for idx, label in enumerate(labels)}
    return SimpleNamespace(getColumnHeaders=lambda: headers)


@pytest.mark.parametrize(
    ("labels", "expected"),
    [
        (("no_label", "Lys8Arg10"), ["SILAC light", "SILAC heavy"]),
        # A two-plex calls its non-light channel heavy, whatever its mass class.
        (("no_label", "Lys4Arg6"), ["SILAC light", "SILAC heavy"]),
        (("no_label", "Lys4Arg6", "Lys8Arg10"), ["SILAC light", "SILAC medium", "SILAC heavy"]),
        (("SILAC light", "SILAC heavy"), ["SILAC light", "SILAC heavy"]),
        (("Dimethyl0", "Dimethyl4", "Dimethyl8"), ["DIMETHYL0", "DIMETHYL4", "DIMETHYL8"]),
        (("mTRAQ0", "mTRAQ4"), ["MTRAQ0", "MTRAQ4"]),
        (("tmt6plex_126", "tmt6plex_127"), ["TMT126", "TMT127"]),
        (("label-free",), ["LFQ"]),
        (("",), ["LFQ"]),
        # Identical unmarked columns are one label-free channel linked twice.
        (("no_label", "no_label"), ["LFQ", "LFQ"]),
    ],
)
def test_column_labels(labels, expected):
    assert list(column_labels(_map(*labels)).values()) == expected


@pytest.mark.parametrize(
    "labels",
    [
        ("Lys8Arg10",),
        ("Lys8Arg10", "Lys8"),
        ("no_label", "Lys4Arg6", "Lys4"),
        ("no_label", "Lys6Arg10"),
        ("no_label", "Dimethyl4"),
        ("ICPL0", "ICPL4"),
    ],
)
def test_ambiguous_channel_labels_stop_the_conversion(labels):
    with pytest.raises(ValueError, match="run.mzML"):
        column_labels(_map(*labels))


def test_streaming_does_not_round_scores():
    """Precision is up to pyopenms: rounding here would make the readers differ where it keeps doubles."""
    assert _parse_hit_score({"score": "0.2"}) == 0.2


def test_both_readers_report_the_same_scores(tmp_path):
    memory, stream = (_convert(tmp_path, "MS1LabeledWorkflow_6_ratios_out", streaming)["feature"] for streaming in (False, True))
    assert {row["feature_id"]: row["peptide_qvalue"] for row in stream} == {
        row["feature_id"]: row["peptide_qvalue"] for row in memory
    }
