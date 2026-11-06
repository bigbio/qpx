"""Proteome Discoverer text exports: both header spellings, reporter channels, file groups and identification-only input."""

import pandas as pd
import pyarrow.parquet as pq
import pytest

from qpx.converters.proteome_discoverer import ProteomeDiscovererConverter
from qpx.converters.proteome_discoverer.pg_adapter import _abundance_columns
from qpx.converters.proteome_discoverer.reports import peptide_identity

_RUN = "plex1_F1"


def _psm_row(**overrides):
    row = {
        "Confidence": "High",
        "Identifying Node": "Sequest HT (A4)",
        "PSM Ambiguity": "Unambiguous",
        "Annotated Sequence": "[K].mPEPTIDEk.[A]",
        "Modifications": "N-Term(TMTpro); K9(TMTpro)",
        "Master Protein Accessions": "P12345",
        "Protein Accessions": "P12345",
        "Number of Missed Cleavages": "0",
        "Charge": "2",
        "Rank": "1",
        "mz in Da": "800.5",
        "Delta M in ppm": "1.5",
        "Isolation Interference in Percent": "10",
        "Intensity": "5000",
        "RT in min": "10",
        "First Scan": "100",
        "File ID": "F1.1",
        "Abundance 126": "10",
        "Abundance 127N": "",
        "Quan Info": "",
        "XCorr": "3.5",
        "q-Value": "0.001",
        "PEP": "0.01",
    }
    return {**row, **overrides}


def _proteins(tmp_path):
    rows = [
        {
            "Master": "IsMasterProtein",
            "Protein Group IDs": "1",
            "Accession": "P12345",
            "Description": "A protein",
            "Exp Protein Group q-value Combined": "0.0",
            "Contaminant": "False",
            "Number of Peptides": "3",
            "Number of Unique Peptides": "2",
            "Number of AAs": "100",
            "Coverage in Percent": "4",
            "MW in kDa": "20.5",
            "Gene Symbol": "GENE1",
            "Abundance F1 126 Sample CORT": "100",
            "Abundance F1 127N Sample CORT": "",
            "Abundances Normalized F1 126 Sample CORT": "90",
        },
        {"Master": "None", "Protein Group IDs": "1", "Accession": "P99999"},
    ]
    path = tmp_path / "Proteins.txt"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def _input_files(tmp_path):
    path = tmp_path / "InputFiles.txt"
    name = rf"E:\PD\{_RUN}.raw"
    pd.DataFrame([{"Study File ID": "F1.1", "File Name": name}]).to_csv(path, sep="\t", index=False)
    return path


def _convert(tmp_path, psm_rows, *, proteins=True, **options):
    psms = tmp_path / "PSMs.txt"
    pd.DataFrame(psm_rows).to_csv(psms, sep="\t", index=False)
    out = tmp_path / "out"
    ProteomeDiscovererConverter(out).convert(
        psms,
        _proteins(tmp_path) if proteins else None,
        input_files=_input_files(tmp_path),
        output_prefix="pd",
        **options,
    )
    return out


@pytest.mark.parametrize(
    ("modifications", "peptidoform"),
    [
        ("N-Term(TMTpro); K9(TMTpro)", "[UNIMOD:2016]-MPEPTIDEK[UNIMOD:2016]"),
        ("M1(Oxidation)", "M[UNIMOD:35]PEPTIDEK"),
        ("N-Term(Prot)(Acetyl)", "[UNIMOD:1]-MPEPTIDEK"),
        ("", "MPEPTIDEK"),
    ],
)
def test_modifications_become_a_peptidoform(modifications, peptidoform):
    assert peptide_identity("MPEPTIDEK", modifications)[0] == peptidoform


def test_unreadable_modification_is_an_error():
    with pytest.raises(ValueError, match="Unreadable"):
        peptide_identity("PEPTIDEK", "somewhere")


def test_psm_takes_run_from_input_files_and_scores_from_either_spelling(tmp_path):
    psm = pq.read_table(_convert(tmp_path, [_psm_row()]) / "pd.psm.parquet").to_pylist()[0]

    assert (psm["run_file_name"], psm["scan"], psm["charge"]) == (_RUN, [100], 2)
    assert psm["sequence"] == "MPEPTIDEK"
    assert psm["rt"] == pytest.approx(600.0)
    assert (psm["qvalue"], psm["posterior_error_probability"]) == (pytest.approx(0.001), pytest.approx(0.01))
    assert psm["pass_threshold"] is True
    scores = {score["score_name"]: score["score_value"] for score in psm["additional_scores"]}
    assert scores["rank"] == 1
    assert scores["parent_ion_fraction"] == pytest.approx(0.9)


def test_ptmrs_probabilities_localize_the_modification(tmp_path):
    row = _psm_row(
        **{
            "Annotated Sequence": "SAtK",
            "Modifications": "T3(Phospho)",
            "ptmRS Phospho Site Probabilities": "S(1): 1.5; T(3): 98.5",
        }
    )
    psm = pq.read_table(_convert(tmp_path, [row]) / "pd.psm.parquet").to_pylist()[0]

    assert psm["peptidoform"] == "SAT[UNIMOD:21]K"
    candidates = {(c["position"], c["assigned"]): c["probability"] for c in psm["modification_candidates"]}
    assert candidates == {(1, False): pytest.approx(0.015), (3, True): pytest.approx(0.985)}


def test_two_engines_on_one_spectrum_keep_one_psm_and_rank1_only_drops_the_others(tmp_path):
    rows = [
        _psm_row(**{"Identifying Node": "Comet", "q-Value": "0.01"}),
        _psm_row(),
        _psm_row(**{"Annotated Sequence": "OTHERK", "Modifications": "", "Rank": "2"}),
    ]
    out = _convert(tmp_path, rows)
    psms = pq.read_table(out / "pd.psm.parquet").to_pylist()
    assert sorted(psm["sequence"] for psm in psms) == ["MPEPTIDEK", "OTHERK"]
    assert next(psm for psm in psms if psm["sequence"] == "MPEPTIDEK")["qvalue"] == pytest.approx(0.001)

    out = _convert(tmp_path, rows, rank1_only=True)
    assert pq.read_table(out / "pd.psm.parquet").num_rows == 1


def test_feature_sums_the_reporters_of_the_used_psms(tmp_path):
    rows = [
        _psm_row(),
        _psm_row(**{"First Scan": "101", "Abundance 126": "5", "Abundance 127N": "7", "q-Value": "0.1"}),
        _psm_row(**{"First Scan": "102", "Abundance 126": "1000", "Quan Info": "ExcludedByMethod"}),
        _psm_row(**{"First Scan": "103", "Abundance 126": "1000", "PSM Ambiguity": "Rejected"}),
        _psm_row(**{"First Scan": "104", "Abundance 126": "1000", "Rank": "2"}),
    ]
    features = pq.read_table(_convert(tmp_path, rows) / "pd.feature.parquet").to_pylist()

    assert len(features) == 1
    feature = features[0]
    assert feature["intensities"] == [{"label": "TMT126", "intensity": 15.0}, {"label": "TMT127N", "intensity": 7.0}]
    assert feature["additional_intensities"][0]["intensities"][0]["intensity_name"] == "ms1_intensity"
    assert feature["scan"] == [100]
    assert feature["unique"] is True
    assert feature["anchor_protein"] == "P12345"
    assert [protein["accession"] for protein in feature["pg_accessions"]] == ["P12345", "P99999"]
    assert feature["gg_names"] == ["GENE1"]


def test_label_free_feature_takes_the_largest_ms1_abundance(tmp_path):
    rows = [
        {k: v for k, v in _psm_row(**{"Intensity": intensity, "First Scan": scan}).items() if not k.startswith("Abundance")}
        for intensity, scan in (("100", "1"), ("300", "2"))
    ]
    feature = pq.read_table(_convert(tmp_path, rows) / "pd.feature.parquet").to_pylist()[0]

    assert feature["intensities"] == [{"label": "LFQ", "intensity": 300.0}]
    assert feature["rt"] == pytest.approx(600.0)
    assert feature["additional_intensities"][0]["intensities"] == [{"intensity_name": "spectral_count", "intensity_value": 2.0}]


def test_ms1_label_abundances_stop_the_conversion(tmp_path):
    row = {k: v for k, v in _psm_row().items() if not k.startswith("Abundance")}
    with pytest.raises(ValueError, match="MS1-label"):
        _convert(tmp_path, [{**row, "Abundance Light": "10"}])


def test_export_without_runs_is_refused(tmp_path):
    psms = tmp_path / "PSMs.txt"
    pd.DataFrame([_psm_row()]).to_csv(psms, sep="\t", index=False)
    with pytest.raises(ValueError, match="InputFiles"):
        ProteomeDiscovererConverter(tmp_path / "out").convert(psms)


def test_pg_has_one_row_per_quantified_channel_with_its_group_and_normalized_abundance(tmp_path):
    rows = pq.read_table(_convert(tmp_path, [_psm_row()]) / "pd.pg.parquet").to_pylist()

    assert len(rows) == 1
    row = rows[0]
    assert (row["label"], row["intensity"]) == ("TMT126", 100.0)
    assert row["grouped_runs"] == [_RUN]
    assert row["pg_accessions"] == ["P12345", "P99999"]
    assert row["anchor_protein"] == "P12345"
    assert row["peptide_counts"] == {"unique_sequences": 2, "total_sequences": 3}
    assert row["contaminant"] is False
    assert row["additional_intensities"] == [
        {"label": "TMT126", "intensities": [{"intensity_name": "normalized_abundance", "intensity_value": 90.0}]}
    ]


def test_identification_only_export_writes_pg_without_a_quantity_and_no_feature(tmp_path):
    psm = {k: v for k, v in _psm_row().items() if not k.startswith("Abundance") and k != "Intensity"}
    proteins = pd.read_csv(_proteins(tmp_path), sep="\t", dtype=str)
    proteins = proteins[[column for column in proteins.columns if not column.startswith("Abundance")]]
    psms = tmp_path / "PSMs.txt"
    pd.DataFrame([{**psm, "Spectrum File": f"{_RUN}.raw"}]).to_csv(psms, sep="\t", index=False)
    proteins_file = tmp_path / "ids.txt"
    proteins.to_csv(proteins_file, sep="\t", index=False)
    out = tmp_path / "out"

    ProteomeDiscovererConverter(out).convert(psms, proteins_file, output_prefix="pd")

    assert not (out / "pd.feature.parquet").exists()
    row = pq.read_table(out / "pd.pg.parquet").to_pylist()[0]
    assert (row["label"], row["intensity"], row["grouped_runs"]) == ("LFQ", None, [_RUN])


@pytest.mark.parametrize(
    ("column", "expected"),
    [
        ("Abundance: F16: 127N, Sample, CORT 31h", ("Abundance", "F16", "TMT127N")),
        ("Abundance F16 127N Sample CORT 31h", ("Abundance", "F16", "TMT127N")),
        ("Abundance: F1: Sample", ("Abundance", "F1", None)),
        ("Abundances (Grouped): F1, 131", ("Abundances Grouped", "F1", "TMT131")),
        ("Abundances (Grouped): F1", ("Abundances Grouped", "F1", None)),
        ("Abundances Scaled F1 Sample", ("Abundances Scaled", "F1", None)),
        ("Abundances (Normalized): F2: 126, Control", ("Abundances Normalized", "F2", "TMT126")),
    ],
)
def test_abundance_column_spellings(column, expected):
    assert _abundance_columns([column]) == {column: expected}


@pytest.mark.parametrize(
    "column",
    ["Abundance Ratio: (F1, 127N) / (F1, 126)", "Abundances (Grouped): Treated", "Abundances Grouped CV in Percent CORT"],
)
def test_columns_that_are_no_file_group_abundance(column):
    assert _abundance_columns([column]) == {}


def test_grouped_abundance_is_used_only_without_the_per_file_one(tmp_path):
    proteins = pd.read_csv(_proteins(tmp_path), sep="\t", dtype=str).drop(
        columns=["Abundance F1 126 Sample CORT", "Abundance F1 127N Sample CORT", "Abundances Normalized F1 126 Sample CORT"]
    )
    proteins["Abundances (Grouped): F1, 126"] = "55"
    proteins_file = tmp_path / "grouped.txt"
    proteins.to_csv(proteins_file, sep="\t", index=False)
    psms = tmp_path / "PSMs.txt"
    pd.DataFrame([_psm_row()]).to_csv(psms, sep="\t", index=False)
    out = tmp_path / "out"

    ProteomeDiscovererConverter(out).convert(psms, proteins_file, input_files=_input_files(tmp_path), output_prefix="pd")

    row = pq.read_table(out / "pd.pg.parquet").to_pylist()[0]
    assert (row["label"], row["intensity"]) == ("TMT126", 55.0)
