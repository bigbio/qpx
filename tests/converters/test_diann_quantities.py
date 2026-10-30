"""DIA-NN 2.x reports: unmeasured zeros, measured m/z, plexDIA channels and the report-only pg view."""

import logging

import pandas as pd
import pyarrow.parquet as pq
import pytest

from qpx.cli.transform.quantify import _quantify_sdrf_mapping
from qpx.converters.diann.converter import parse_diann_mod_masses
from qpx.converters.diann.feature_adapter import DiannFeatureAdapter
from qpx.converters.diann.pg_adapter import DiannPgAdapter

# Theoretical m/z of PEPTIDEK 2+ and of its mTRAQ-light form, as DIA-NN reports in Precursor.Mz.
_MZ = 464.7347285
_MTRAQ_MZ = 604.8296915


def _row(**overrides) -> dict:
    row = {
        "Run": "run_A.raw",
        "Protein.Group": "P1",
        "Protein.Ids": "P1",
        "Protein.Names": "N1",
        "Genes": "G1",
        "Modified.Sequence": "PEPTIDEK",
        "Stripped.Sequence": "PEPTIDEK",
        "Precursor.Id": "PEPTIDEK2",
        "Precursor.Charge": 2,
        "Proteotypic": 1,
        "Q.Value": 0.001,
        "Precursor.Quantity": 100.0,
        "Precursor.Mz": _MZ,
        "Ms1.Apex.Mz.Delta": 0.0023,
        "RT": 10.0,
        "IM": 0.0,
        "Predicted.IM": 0.0,
        "PG.MaxLFQ": 900.0,
        "PG.Q.Value": 0.002,
        "Global.PG.Q.Value": 0.003,
        "GG.Q.Value": 0.004,
    }
    row.update(overrides)
    return row


def _sdrf(tmp_path, labels) -> str:
    path = tmp_path / "plex.sdrf.tsv"
    pd.DataFrame(
        [{"source name": f"S_{label}", "comment[data file]": "run_A.raw", "comment[label]": label} for label in labels]
    ).to_csv(path, sep="\t", index=False)
    return str(path)


def _features(tmp_path, rows, sdrf=None, mod_masses=()) -> list[dict]:
    report = tmp_path / "report.tsv"
    pd.DataFrame(rows).to_csv(report, sep="\t", index=False)
    output = tmp_path / "out.feature.parquet"
    with DiannFeatureAdapter() as adapter:
        adapter.convert(diann_report=str(report), output_path=str(output), sdrf_path=sdrf, mod_masses=mod_masses)
    return pq.read_table(output).to_pylist()


def _pg(tmp_path, rows, sdrf=None) -> list[dict]:
    report = tmp_path / "report.tsv"
    pd.DataFrame(rows).to_csv(report, sep="\t", index=False)
    output = tmp_path / "out.pg.parquet"
    with DiannPgAdapter() as adapter:
        adapter.convert(diann_report=str(report), pg_matrix_path=None, output_path=str(output), sdrf_path=sdrf)
    return pq.read_table(output).to_pylist()


def _mtraq_rows(**per_channel) -> list[dict]:
    """One mTRAQ precursor in channels 0, 4 and 8 (DIA-NN 2.x plexDIA), each overridden by *per_channel*."""
    rows = []
    for channel, offset in (("0", 0.0), ("4", 2.0035497), ("8", 4.0070994)):
        rows.append(
            _row(
                **{
                    "Channel": channel,
                    "Modified.Sequence": "(mTRAQ)PEPTIDEK(mTRAQ)",
                    "Precursor.Id": f"(mTRAQ-{channel})PEPTIDEK(mTRAQ-{channel})2",
                    "Precursor.Mz": _MTRAQ_MZ + offset,
                    **per_channel.get(channel, {}),
                }
            )
        )
    return rows


def test_unmeasured_zeros_are_null(tmp_path):
    """DIA-NN writes 0 for a quantity it has not estimated and for ion mobility without IMS."""
    rows = [
        _row(**{"Precursor.Quantity": 0.0}),
        _row(**{"Modified.Sequence": "PEPTIDER", "Stripped.Sequence": "PEPTIDER", "Precursor.Id": "PEPTIDER2", "IM": 0.9}),
    ]

    features = {row["sequence"]: row for row in _features(tmp_path, rows)}

    assert features["PEPTIDEK"]["intensities"] == [{"label": "LFQ", "intensity": None}]
    assert features["PEPTIDEK"]["ion_mobility"] is None
    assert "predicted_im" not in {param["cv_name"] for param in features["PEPTIDEK"]["cv_params"]}
    assert features["PEPTIDER"]["intensities"] == [{"label": "LFQ", "intensity": 100.0}]
    assert features["PEPTIDER"]["ion_mobility"] == pytest.approx(0.9)


def test_observed_mz_is_the_measured_ms1_apex(tmp_path):
    """Precursor.Mz is theoretical; the measured m/z adds Ms1.Apex.Mz.Delta, and a zero offset measured nothing."""
    rows = [
        _row(),
        _row(
            **{
                "Modified.Sequence": "PEPTIDER",
                "Stripped.Sequence": "PEPTIDER",
                "Precursor.Id": "PEPTIDER2",
                "Ms1.Apex.Mz.Delta": 0.0,
            }
        ),
    ]

    features = {row["sequence"]: row for row in _features(tmp_path, rows)}

    measured = features["PEPTIDEK"]
    assert measured["observed_mz"] == pytest.approx(_MZ + 0.0023, abs=1e-4)
    assert measured["mass_error_ppm"] == pytest.approx(0.0023 / _MZ * 1e6, abs=0.1)
    assert features["PEPTIDER"]["observed_mz"] is None
    assert features["PEPTIDER"]["mass_error_ppm"] is None


def test_unique_is_a_one_protein_group_of_a_peptide_mapping_nowhere_else(tmp_path):
    """Proteotypic follows DIA-NN's grouping level; unique follows feature.yaml through Protein.Ids."""
    rows = [
        _row(),
        _row(
            **{
                "Modified.Sequence": "PEPTIDER",
                "Stripped.Sequence": "PEPTIDER",
                "Precursor.Id": "PEPTIDER2",
                "Protein.Ids": "P1;P2",
            }
        ),
        _row(
            **{
                "Modified.Sequence": "PEPTIDEH",
                "Stripped.Sequence": "PEPTIDEH",
                "Precursor.Id": "PEPTIDEH2",
                "Protein.Group": "P1;P2",
                "Protein.Ids": "P1;P2",
            }
        ),
    ]

    unique = {row["sequence"]: row["unique"] for row in _features(tmp_path, rows)}

    assert unique == {"PEPTIDEK": True, "PEPTIDER": False, "PEPTIDEH": False}


def test_plexdia_channels_resolve_to_canonical_sdrf_labels(tmp_path):
    """DIA-NN names mTRAQ channels 0/4/8; the SDRF declares MTRAQ0/4/8, and the feature
    peptidoform is the light partner's in canonical form."""
    sdrf = _sdrf(tmp_path, ["MTRAQ0", "MTRAQ4", "MTRAQ8"])

    feature = _features(tmp_path, _mtraq_rows(), sdrf)[0]
    pg = _pg(tmp_path, _mtraq_rows(), sdrf)

    assert [entry["label"] for entry in feature["intensities"]] == ["MTRAQ0", "MTRAQ4", "MTRAQ8"]
    assert feature["peptidoform"] == "[UNIMOD:888]-PEPTIDEK[UNIMOD:888]"
    assert feature["calculated_mz"] == pytest.approx(_MTRAQ_MZ)
    assert sorted(row["label"] for row in pg) == ["MTRAQ0", "MTRAQ4", "MTRAQ8"]


def test_plexdia_representative_channel_is_deterministic(tmp_path):
    """Scalar fields come from the best q-value channel, ties broken by label; m/z describes the lightest partner."""
    tied = _mtraq_rows(**{"0": {"RT": 10.0}, "4": {"RT": 10.2}, "8": {"RT": 10.4}})
    best_heavy = _mtraq_rows(**{"0": {"RT": 10.0, "Q.Value": 0.002}, "8": {"RT": 10.4, "Q.Value": 0.0005}})

    tied_rt = [_features(tmp_path, rows)[0]["rt"] for rows in (tied, tied[::-1])]
    heavy = _features(tmp_path, best_heavy[::-1])[0]

    assert tied_rt == [pytest.approx(600.0)] * 2
    assert heavy["rt"] == pytest.approx(624.0)
    assert heavy["observed_mz"] == pytest.approx(_MTRAQ_MZ + 0.0023, abs=1e-4)


def test_plexdia_pg_keeps_each_channels_scores_without_a_matrix(tmp_path):
    """DIA-NN writes no pg matrix for a channel report; each channel keeps its own q-value,
    and a channel MaxLFQ of 0 is a labelled row with a null intensity."""
    rows = _mtraq_rows(**{"0": {"PG.Q.Value": 0.001}, "4": {"PG.Q.Value": 0.02}, "8": {"PG.Q.Value": 0.03, "PG.MaxLFQ": 0.0}})

    pg = {row["label"]: row for row in _pg(tmp_path, rows)}

    assert {label: row["pg_qvalue"] for label, row in pg.items()} == {
        "0": pytest.approx(0.001),
        "4": pytest.approx(0.02),
        "8": pytest.approx(0.03),
    }
    assert pg["8"]["intensity"] is None
    assert pg["0"]["intensity"] == 900.0


def test_silac_channels_drop_the_label_tag(tmp_path):
    """DIA-NN marks plexDIA SILAC residues with a zero-mass (SILAC) tag; the light peptidoform has no label."""
    rows = [
        _row(**{"Channel": channel, "Modified.Sequence": "PEPTIDEK(SILAC)", "Precursor.Id": f"PEPTIDEK(SILAC-{channel})2"})
        for channel in ("L", "H")
    ]

    feature = _features(tmp_path, rows, _sdrf(tmp_path, ["SILAC light", "SILAC heavy"]))[0]

    assert feature["peptidoform"] == "PEPTIDEK"
    assert [entry["label"] for entry in feature["intensities"]] == ["SILAC heavy", "SILAC light"]


def test_modifications_named_in_the_log_resolve_by_mass(tmp_path):
    """--original-mods reports keep the declared names; the logged declaration gives their mass."""
    log = tmp_path / "diann.log"
    log.write_text("diann --f run_A.raw --var-mod Phos,79.966331,STY --fixed-mod Carbamidomethyl,57.021464,C\n")
    rows = [_row(**{"Modified.Sequence": "PEPS(Phos)IDEK", "Stripped.Sequence": "PEPSIDEK", "Precursor.Id": "PEPS(Phos)IDEK2"})]

    masses = parse_diann_mod_masses(str(log))
    feature = _features(tmp_path, rows, mod_masses=masses)[0]

    assert masses == (("Phos", 79.966331), ("Carbamidomethyl", 57.021464))
    assert feature["peptidoform"] == "PEPS[UNIMOD:21]IDEK"
    assert feature["calculated_mz"] is not None


def test_pg_zero_maxlfq_is_null_and_warns(tmp_path, caplog):
    """DIA-NN 2.x has no PG.Quantity and writes PG.MaxLFQ = 0 for a group it did not quantify."""
    with caplog.at_level(logging.WARNING):
        pg = _pg(tmp_path, [_row(**{"PG.MaxLFQ": 0.0})])

    assert [(row["label"], row["intensity"], row["cv_params"]) for row in pg] == [("LFQ", None, None)]
    assert "No DIA-NN protein group has a quantity" in caplog.text


def test_pg_carries_the_spec_mapped_quantities_and_scores(tmp_path):
    """pg.md maps PG.Normalised, PG.TopN, Genes.MaxLFQ(.Unique) and the QuantUMS qualities."""
    row = _row(
        **{
            "PG.Quantity": 1000.0,
            "PG.Normalised": 950.0,
            "PG.TopN": 0.0,
            "Genes.MaxLFQ": 880.0,
            "Genes.MaxLFQ.Unique": 870.0,
            "PG.MaxLFQ.Quality": 0.8,
            "Genes.MaxLFQ.Quality": 0.7,
            "Lib.PG.Q.Value": 0.005,
            "Protein.Q.Value": 0.006,
        }
    )

    pg = _pg(tmp_path, [row])[0]

    additional = {entry["intensity_name"]: entry["intensity_value"] for entry in pg["additional_intensities"][0]["intensities"]}
    scores = {score["score_name"]: (score["score_value"], score["higher_better"]) for score in pg["additional_scores"]}
    assert pg["intensity"] == 1000.0
    assert additional == {"maxlfq": 900.0, "normalize_intensity": 950.0, "genes_maxlfq": 880.0, "genes_maxlfq_unique": 870.0}
    assert scores["pg_maxlfq_quality"] == (pytest.approx(0.8), True)
    assert scores["genes_maxlfq_quality"] == (pytest.approx(0.7), True)
    assert scores["lib_pg_qvalue"] == (pytest.approx(0.005), False)
    assert scores["protein_qvalue"] == (pytest.approx(0.006), False)


def test_quantify_reads_numeric_plexdia_labels(tmp_path):
    """An SDRF of DIA-NN channel labels 0/4/8 maps them as text instead of crashing on numbers."""
    path = tmp_path / "numeric.sdrf.tsv"
    pd.DataFrame(
        [{"source name": f"S{label}", "comment[data file]": "run_A.raw", "comment[label]": label} for label in (0, 4, 8)]
    ).to_csv(path, sep="\t", index=False)

    mapping = _quantify_sdrf_mapping(path)

    assert mapping["_label_key"].tolist() == ["0", "4", "8"]


def test_pepmap_lists_every_protein_of_the_converted_peptidoforms(tmp_path):
    """Protein.Ids reaches the pepmap view: all mapped proteins, target rows within the q-value cut only."""
    rows = [
        _row(**{"Protein.Group": "P1", "Protein.Ids": "P1;P9"}),
        _row(**{"Run": "run_B.raw", "Protein.Group": "P1", "Protein.Ids": "P1;P8"}),
        _row(**{"Modified.Sequence": "SAMPLER", "Stripped.Sequence": "SAMPLER", "Protein.Ids": "P2", "Protein.Group": "P2"}),
        _row(**{"Modified.Sequence": "LATER", "Stripped.Sequence": "LATER", "Q.Value": 0.2}),
        _row(**{"Decoy": 1, "Protein.Ids": "P7"}),
    ]
    rows = [{"Decoy": 0, **row} for row in rows]
    report = tmp_path / "report.tsv"
    pd.DataFrame(rows).to_csv(report, sep="\t", index=False)
    pepmap = tmp_path / "out.pepmap.parquet"
    with DiannFeatureAdapter() as adapter:
        adapter.convert(diann_report=str(report), output_path=str(tmp_path / "out.feature.parquet"), qvalue_threshold=0.01)
        assert adapter.write_pepmap(str(pepmap))

    entries = {
        row["peptidoform"]: ([protein["accession"] for protein in row["pg_accessions"]], row["is_unique"])
        for row in pq.read_table(pepmap).to_pylist()
    }
    assert entries == {"PEPTIDEK": (["P1", "P8", "P9"], False), "SAMPLER": (["P2"], True)}
