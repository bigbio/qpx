"""FragPipe reporter-ion quantities of one plex, read through its annotation.txt."""

import pandas as pd
import pyarrow.parquet as pq
import pytest

from qpx.converters.fragpipe import FragPipeConverter
from qpx.converters.fragpipe.isobaric import load_channel_annotation, load_tmt_report

_RUNS = ("plex1_F1", "plex1_F2")
_ION_ID = {
    "Peptide Sequence": "PEPTIDEK",
    "Modified Sequence": "PEPTIDEK",
    "M/Z": 464.7,
    "Charge": 2,
    "Assigned Modifications": "",
    "Protein": "sp|P12345|TEST_HUMAN",
    "Protein ID": "P12345",
    "Gene": "TEST",
    "Mapped Proteins": "",
}


def _plex(tmp_path, *, ion_channels=(100.0, 0.0), protein_channels=(10.0, 0.0)):
    (tmp_path / "annotation.txt").write_text("126 sample_a\n127N sample_b\n", encoding="utf-8")
    psm = {
        "Spectrum": [f"{run}.00010.00010.2" for run in _RUNS],
        "Peptide": "PEPTIDEK",
        "Charge": 2,
        "Calculated Peptide Mass": 927.4,
        "Delta Mass": 0.0,
        "Assigned Modifications": "",
        "Protein": "sp|P12345|TEST_HUMAN",
        "PeptideProphet Probability": 1.0,
        "Is Unique": True,
    }
    pd.DataFrame(psm).to_csv(tmp_path / "psm.tsv", sep="\t", index=False)
    ion = {**_ION_ID, "Intensity": 5000.0, "sample_a": ion_channels[0], "sample_b": ion_channels[1]}
    pd.DataFrame([ion]).to_csv(tmp_path / "ion.tsv", sep="\t", index=False)
    protein = {
        "Protein": "contam_sp|P12345|TEST_HUMAN",
        "Protein ID": "P12345",
        "Gene": "TEST",
        "Length": 100,
        "Protein Probability": 1.0,
        "Total Peptides": 3,
        "Unique Peptides": 2,
        "Indistinguishable Proteins": "",
        "sample_a": protein_channels[0],
        "sample_b": protein_channels[1],
    }
    pd.DataFrame([protein]).to_csv(tmp_path / "protein.tsv", sep="\t", index=False)
    return tmp_path


def _convert(plex, tmt_report=None):
    converter = FragPipeConverter(output_directory=plex / "out")
    converter.load_plex(plex / "annotation.txt", tmt_report)
    converter.convert(
        psm_file=plex / "psm.tsv",
        ion_file=plex / "ion.tsv",
        pg_file=plex / "protein.tsv",
        output_prefix="plex1",
    )
    return plex / "out"


def test_annotation_maps_channel_columns_to_canonical_labels(tmp_path):
    annotation = tmp_path / "annotation.txt"
    annotation.write_text("126 sample_a\n134N a sample with spaces\n\n", encoding="utf-8")
    assert load_channel_annotation(annotation) == {"sample_a": "TMT126", "a sample with spaces": "TMT134N"}


def test_feature_holds_the_channels_with_the_ms1_intensity_beside_them(tmp_path):
    feature = pq.read_table(_convert(_plex(tmp_path)) / "plex1.feature.parquet").to_pylist()[0]

    assert feature["intensities"] == [
        {"label": "TMT126", "intensity": 100.0},
        {"label": "TMT127N", "intensity": None},
    ]
    assert feature["additional_intensities"] == [
        {"label": "TMT126", "intensities": [{"intensity_name": "ms1_intensity", "intensity_value": 5000.0}]}
    ]
    assert (feature["run_file_name"], feature["quantification_unit_id"]) == (_RUNS[0], "plex1")
    assert feature["id_run_file_name"] is None


def test_a_precursor_without_any_reporter_quantity_writes_no_feature(tmp_path):
    out = _convert(_plex(tmp_path, ion_channels=(0.0, 0.0)))
    assert not (out / "plex1.feature.parquet").exists()


def test_pg_has_one_row_per_quantified_channel_over_the_plex_runs(tmp_path):
    rows = pq.read_table(_convert(_plex(tmp_path)) / "plex1.pg.parquet").to_pylist()

    assert [(row["label"], row["intensity"]) for row in rows] == [("TMT126", 10.0)]
    assert rows[0]["grouped_runs"] == list(_RUNS)
    assert rows[0]["peptide_counts"] == {"unique_sequences": 2, "total_sequences": 3}
    assert rows[0]["contaminant"] is True


def test_tmt_report_adds_the_log2_abundance_to_the_pg_row(tmp_path):
    plex = _plex(tmp_path)
    report = plex / "abundance_protein_MD.tsv"
    row = {"Index": "P12345", "NumberPSM": 3, "Gene": "TEST", "MaxPepProb": 1.0, "ReferenceIntensity": 20.0}
    pd.DataFrame([{**row, "sample_a": 19.5, "sample_b": 20.5}]).to_csv(report, sep="\t", index=False)

    rows = pq.read_table(_convert(plex, report) / "plex1.pg.parquet").to_pylist()

    assert rows[0]["additional_intensities"] == [
        {"label": "TMT126", "intensities": [{"intensity_name": "tmt_integrator_log2_abundance", "intensity_value": 19.5}]}
    ]


def test_gene_level_tmt_report_is_refused(tmp_path):
    report = tmp_path / "abundance_gene_MD.tsv"
    pd.DataFrame([{"Index": "TEST", "ProteinID": "P12345", "sample_a": 19.5}]).to_csv(report, sep="\t", index=False)
    with pytest.raises(ValueError, match="protein-level"):
        load_tmt_report(report, {"sample_a": "TMT126"})


def test_annotation_without_psm_file_is_refused(tmp_path):
    plex = _plex(tmp_path)
    converter = FragPipeConverter(output_directory=plex / "out")
    converter.load_plex(plex / "annotation.txt")
    with pytest.raises(ValueError, match="needs psm.tsv"):
        converter.convert(ion_file=plex / "ion.tsv")


def test_annotated_sample_missing_from_the_report_is_refused(tmp_path):
    plex = _plex(tmp_path)
    (plex / "annotation.txt").write_text("126 sample_a\n127N renamed\n", encoding="utf-8")
    with pytest.raises(ValueError, match="renamed"):
        _convert(plex)
