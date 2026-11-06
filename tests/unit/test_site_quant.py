"""Site quantities: MaxQuant Sites.txt and the derivation from features."""

import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.converters.maxquant.converter import MaxQuantConverter
from qpx.converters.maxquant.sites import convert_sites_file
from qpx.writers import FeatureWriter, SiteWriter
from tests.conftest import make_feature_record

_RUNS = {"Mix1": ["r1", "r2"], "Mix2": ["r3"]}


def _write_sites(tmp_path, header, *rows, name="Phospho (STY)Sites.txt"):
    path = tmp_path / name
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(row) + "\n" for row in rows))
    return path


def _convert(path, tmp_path, runs=None, sdrf=lambda columns: ({}, "LFQ", [])):
    out = tmp_path / "x.site.parquet"
    with SiteWriter(out) as writer:
        convert_sites_file(path, runs or _RUNS, sdrf, writer)
    return pq.read_table(out).to_pylist()


_BASE = ["Protein", "Position", "Amino acid", "Localization prob", "Reverse", "Potential contaminant", "MS/MS IDs"]


def test_lfq_sites_are_per_experiment(tmp_path):
    path = _write_sites(
        tmp_path,
        [*_BASE, "Intensity Mix1", "Intensity Mix2", "Localization prob Mix1"],
        ["sp|P12345|A_HUMAN", "12", "S", "0.9", "", "", "1;2;3", "100", "0", "0.5"],
    )

    [row] = _convert(path, tmp_path)

    assert (row["protein_accession"], row["position"], row["amino_acid"]) == ("P12345", 12, "S")
    assert (row["modification_name"], row["modification_accession"]) == ("Phospho", "UNIMOD:21")
    assert row["mass_delta"] == pytest.approx(79.966331)
    assert (row["label"], row["intensity"], row["grouped_runs"]) == ("LFQ", 100.0, ["r1", "r2"])
    # The experiment's own localization wins over the one over all experiments.
    assert row["modified_probability"] == 0.5
    assert row["psm_count"] == 3


def test_silac_sites_have_one_row_per_label_with_occupancy(tmp_path):
    path = _write_sites(
        tmp_path,
        [*_BASE, "Ratio H/L", "Intensity L Mix1", "Intensity H Mix1", "Occupancy L Mix1", "Occupancy H Mix1"],
        ["P1", "12", "S", "1", "", "+", "1", "2", "10", "20", "0.1", "0.2"],
    )

    rows = _convert(path, tmp_path)

    assert [(r["label"], r["intensity"], r["stoichiometry"]) for r in rows] == [
        ("SILAC light", 10.0, 0.1),
        ("SILAC heavy", 20.0, 0.2),
    ]
    assert {"cv_name": "potential_contaminant", "cv_value": "true"} in rows[0]["cv_params"]


def test_reporter_sites_take_the_channels_of_the_sdrf(tmp_path):
    path = _write_sites(
        tmp_path,
        [*_BASE, "Reporter intensity corrected 1 Mix1", "Reporter intensity corrected 2 Mix1"],
        ["P1", "12", "S", "1", "", "", "1", "5", "6"],
    )

    rows = _convert(path, tmp_path, sdrf=lambda columns: ({}, "TMT2", [("TMT126", 0), ("TMT127", 1)]))

    assert [(r["label"], r["intensity"]) for r in rows] == [("TMT126", 5.0), ("TMT127", 6.0)]


def test_one_experiment_has_no_suffix_and_an_unquantified_site_stays_a_claim(tmp_path):
    path = _write_sites(
        tmp_path,
        [*_BASE, "Intensity"],
        ["P1", "12", "S", "1", "", "", "1", "7"],
        ["P2", "3", "T", "1", "+", "", "", "0"],
    )

    rows = _convert(path, tmp_path)

    assert (rows[0]["label"], rows[0]["intensity"], rows[0]["grouped_runs"]) == ("LFQ", 7.0, ["r1", "r2", "r3"])
    assert (rows[1]["label"], rows[1]["intensity"], rows[1]["is_decoy"]) == (None, None, True)


def test_an_unknown_file_name_is_refused(tmp_path):
    path = _write_sites(tmp_path, _BASE, ["P1", "1", "S", "1", "", "", "1"], name="NotAModSites.txt")
    with pytest.raises(ValueError, match="NotAModSites.txt"):
        _convert(path, tmp_path)


def test_sites_need_a_way_to_name_the_runs(tmp_path):
    path = _write_sites(tmp_path, [*_BASE, "Intensity"], ["P1", "12", "S", "1", "", "", "1", "7"])
    with pytest.raises(ValueError, match="--evidence-file or --sdrf-file"):
        MaxQuantConverter().convert(tmp_path / "o", sites_files=[path], structures=[])


def _feature(sequence, modifications, intensity, start, end, **fields):
    record = make_feature_record(sequence=sequence, anchor_protein="P1", run_file_name="r1")
    record.update(
        peptidoform=sequence,
        modifications=modifications,
        intensities=[{"label": "LFQ", "intensity": intensity}],
        pg_positions=[{"protein_accession": "P1", "start": start, "end": end}],
        **fields,
    )
    return record


_PHOSPHO = [
    {
        "name": "Phospho",
        "accession": "UNIMOD:21",
        "positions": [
            {
                "position": 2,
                "amino_acid": "S",
                "scores": [{"score_name": "localization_probability", "score_value": 0.8, "higher_better": True}],
            }
        ],
    }
]


def _feature_dataset(tmp_path, records):
    with FeatureWriter(tmp_path / "x.feature.parquet") as writer:
        writer.write_batch(records)


def _derive(tmp_path, *modifications):
    arguments = ["transform", "site", "--dataset", str(tmp_path), "--in-place"]
    return CliRunner().invoke(qpx_main, arguments + [item for name in modifications for item in ("--modification", name)])


def test_sites_derive_from_the_features_that_place_them(tmp_path):
    _feature_dataset(
        tmp_path,
        [
            _feature("ASK", _PHOSPHO, 30.0, 10, 12, charge=2),
            _feature("ASK", _PHOSPHO, 10.0, 10, 12, charge=3),
            _feature("ASK", None, 60.0, 10, 12, charge=4),
            _feature("TTK", None, 1000.0, 50, 52, charge=2),
        ],
    )

    result = _derive(tmp_path, "Phospho")

    assert result.exit_code == 0, result.output
    [row] = pq.read_table(tmp_path / "x.site.parquet").to_pylist()
    assert (row["protein_accession"], row["position"], row["amino_acid"]) == ("P1", 11, "S")
    assert (row["label"], row["intensity"], row["grouped_runs"]) == ("LFQ", 40.0, ["r1"])
    assert row["stoichiometry"] == pytest.approx(0.4)
    assert row["modified_probability"] == 0.8


def test_derivation_says_what_it_needs_and_what_it_cannot_name(tmp_path):
    record = _feature("ASK", _PHOSPHO, 30.0, 10, 12)
    record["pg_positions"] = None
    _feature_dataset(tmp_path, [record])

    assert "protein-properties" in _derive(tmp_path, "Phospho").output
    _feature_dataset(tmp_path, [_feature("ASK", _PHOSPHO, 30.0, 10, 12)])
    assert "Unknown modification" in _derive(tmp_path, "NoSuchModification").output
