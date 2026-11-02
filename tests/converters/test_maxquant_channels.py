"""MaxQuant reporter columns -> SDRF channels, checked against a known channel -> sample truth.

MaxQuant numbers ``Reporter intensity N`` by the channel's position in the plex
(mass order). Each synthetic input gives the channel at position ``p``
(1-based) the intensity ``p * 1000`` and lists the SDRF labels shuffled, so a
label must come back with the intensity of its own position in its plex.
"""

import random

import pyarrow.parquet as pq
import pytest

TMT6 = ["TMT126", "TMT127", "TMT128", "TMT129", "TMT130", "TMT131"]
TMT10 = ["TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131"]
TMT16 = TMT10[:9] + ["TMT131N", "TMT131C", "TMT132N", "TMT132C", "TMT133N", "TMT133C", "TMT134N"]
TMT18 = TMT16 + ["TMT134C", "TMT135N"]
ITRAQ4 = ["ITRAQ114", "ITRAQ115", "ITRAQ116", "ITRAQ117"]
ITRAQ8 = ["ITRAQ113", "ITRAQ114", "ITRAQ115", "ITRAQ116", "ITRAQ117", "ITRAQ118", "ITRAQ119", "ITRAQ121"]

# (case id, plex in position order, SDRF labels used, first MaxQuant column index)
CASES = [
    ("tmt6", TMT6, TMT6, 1),
    ("tmt10", TMT10, TMT10, 0),
    ("tmt16_partial", TMT16, TMT16[:-1], 0),
    ("tmt18", TMT18, TMT18, 1),
    ("itraq4", ITRAQ4, ITRAQ4, 1),
    ("itraq8", ITRAQ8, ITRAQ8, 1),
    # iTRAQ114-117 are channels 2-5 of an 8-plex: the eight reporter columns decide.
    ("itraq8_partial", ITRAQ8, ITRAQ4, 1),
]

_EVIDENCE_COLUMNS = [
    "Sequence",
    "Modified sequence",
    "Charge",
    "Raw file",
    "Experiment",
    "Type",
    "MS/MS scan number",
    "m/z",
    "Calibrated retention time",
    "Calibrated retention time start",
    "Calibrated retention time finish",
    "PEP",
    "Leading razor protein",
    "Leading proteins",
    "Intensity",
    "Mass",
    "id",
]
_EVIDENCE_VALUES = ["PEPTIDEK", "_PEPTIDEK_", "2", "run1", "EXP1", "MULTI-MSMS", "100", "450.25"]
_EVIDENCE_VALUES += ["30.0", "29.5", "30.5", "0.01", "P1", "P1", "1000", "898.5", "1"]


def _reporter(n, offset, suffix=""):
    """``(column, value)`` pairs: the channel at position p reads p * 1000, the corrected one p * 1000 + 1."""
    raw = [(f"Reporter intensity {p + offset}{suffix}", str((p + 1) * 1000)) for p in range(n)]
    corrected = [(f"Reporter intensity corrected {p + offset}{suffix}", str((p + 1) * 1000 + 1)) for p in range(n)]
    return raw + corrected


def _write_tsv(path, pairs):
    path.write_text("\t".join(c for c, _ in pairs) + "\n" + "\t".join(v for _, v in pairs) + "\n")
    return path


def _write_inputs(tmp_path, plex, labels, offset, shuffle=True):
    evidence = _write_tsv(
        tmp_path / "evidence.txt", list(zip(_EVIDENCE_COLUMNS, _EVIDENCE_VALUES)) + _reporter(len(plex), offset)
    )
    pg_head = [("Protein IDs", "P1"), ("Majority protein IDs", "P1"), ("MS/MS count EXP1", "3")]
    protein_groups = _write_tsv(tmp_path / "proteinGroups.txt", pg_head + _reporter(len(plex), offset, " EXP1"))
    ordered = random.Random(7).sample(labels, len(labels)) if shuffle else labels
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\n" + "".join(f"S_{label}\trun1.raw\t{label}\n" for label in ordered)
    )
    return evidence, protein_groups, sdrf


@pytest.mark.parametrize(("plex", "labels", "offset"), [case[1:] for case in CASES], ids=[case[0] for case in CASES])
def test_reporter_columns_map_to_sdrf_channels(tmp_path, plex, labels, offset):
    """Every SDRF label reads its own plex position, in feature and pg, whatever the SDRF order."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter
    from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

    evidence, protein_groups, sdrf = _write_inputs(tmp_path, plex, labels, offset)
    expected = {label: (plex.index(label) + 1) * 1000.0 for label in labels}

    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(tmp_path / "feature.parquet"), sdrf_path=str(sdrf))
    (feature,) = pq.read_table(tmp_path / "feature.parquet").to_pylist()
    assert {e["label"]: e["intensity"] for e in feature["intensities"]} == expected
    corrected = {e["label"]: e["intensities"][0]["intensity_value"] for e in feature["additional_intensities"]}
    assert corrected == {label: value + 1 for label, value in expected.items()}

    with MaxQuantPgAdapter() as adapter:
        adapter.convert(str(protein_groups), str(tmp_path / "pg.parquet"), sdrf_path=str(sdrf), evidence_path=str(evidence))
    pg = pq.read_table(tmp_path / "pg.parquet").to_pylist()
    assert {row["label"]: row["intensity"] for row in pg} == expected


_ITRAQ4_REAGENT = ["iTRAQ reagent 114", "iTRAQ reagent 115", "iTRAQ reagent 116", "iTRAQ reagent 117"]


@pytest.mark.parametrize(
    ("plex", "labels", "written", "offset"),
    [
        (ITRAQ4, _ITRAQ4_REAGENT, _ITRAQ4_REAGENT, 1),
        (TMT10, [f"AC=MS:1002616;NT={label}" for label in TMT10], TMT10, 0),
    ],
    ids=["itraq_reagent", "tmt_ontology"],
)
def test_reporter_labels_keep_the_sdrf_spelling(tmp_path, plex, labels, written, offset):
    """Other spellings of a channel map by the channel they name and keep the label the run view writes."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence, _, sdrf = _write_inputs(tmp_path, plex, labels, offset)
    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(tmp_path / "feature.parquet"), sdrf_path=str(sdrf))
    (feature,) = pq.read_table(tmp_path / "feature.parquet").to_pylist()
    assert {e["label"]: e["intensity"] for e in feature["intensities"]} == {
        label: (index + 1) * 1000.0 for index, label in enumerate(written)
    }


@pytest.mark.parametrize(
    "spell",
    [lambda label: label, lambda label: f"AC=MS:1002616;NT={label}", str.lower],
    ids=["plain", "ontology", "lower_case"],
)
def test_pg_keeps_identified_channels_without_intensity(tmp_path, spell):
    """A channel the SDRF assigns keeps its pg row without an intensity when its plex has MS/MS evidence.

    This holds for every SDRF spelling of the label, including the ontology
    form ``AC=MS:1002616;NT=TMT127C``: the label is matched to the SDRF
    samples in the spelling the channels are keyed by. The SDRF lists the
    labels in plex order; the tests above cover shuffled ones.
    """
    from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

    labels = [spell(label) for label in TMT10]
    evidence, protein_groups, sdrf = _write_inputs(tmp_path, TMT10, labels, 0, shuffle=False)
    pairs = [("Protein IDs", "P1"), ("Majority protein IDs", "P1"), ("MS/MS count EXP1", "3")]
    pairs += [(f"Reporter intensity {p} EXP1", "0" if label == "TMT127C" else "1000") for p, label in enumerate(TMT10)]
    _write_tsv(protein_groups, pairs)
    with MaxQuantPgAdapter() as adapter:
        adapter.convert(str(protein_groups), str(tmp_path / "pg.parquet"), sdrf_path=str(sdrf), evidence_path=str(evidence))
    rows = pq.read_table(tmp_path / "pg.parquet").to_pylist()
    assert len(rows) == len(TMT10)
    assert [row["label"].upper().endswith("TMT127C") for row in rows if row["intensity"] is None] == [True]


@pytest.mark.parametrize(
    ("labels", "match"),
    [
        (["SILAC light", "SILAC heavy"], "SILAC"),
        # TMT127 (6-plex) and TMT127N (10-plex) belong to no single plex.
        (["TMT126", "TMT127", "TMT127N"], "not channels of a TMT plex"),
    ],
    ids=["silac", "mixed_plex"],
)
def test_unmappable_labels_fail_before_writing(tmp_path, labels, match):
    """Labels the converter cannot map to reporter columns stop it before any view is written."""
    from qpx.converters.maxquant.converter import MaxQuantConverter

    evidence, _, sdrf = _write_inputs(tmp_path, TMT10, labels, 0)
    out = tmp_path / "out"
    with pytest.raises(ValueError, match=match):
        MaxQuantConverter().convert(output_folder=str(out), evidence_file=str(evidence), sdrf_file=str(sdrf))
    assert list(out.glob("*.parquet")) == []


_MSMS = (
    "Sequence\tModified sequence\tCharge\tRaw file\tScan number\tm/z\tMass\tPEP\tReverse\tScore\tProteins\n"
    "PEPTIDER\t_PEPTIDER_\t2\trun1\t100\t478.73\t955.46\t0.01\t\t100\tP1\n"
)


@pytest.mark.parametrize("with_evidence", [False, True], ids=["msms_only", "psm_structure"])
@pytest.mark.parametrize(
    "labels",
    [["SILAC light", "SILAC heavy"], ["TMT126", "TMT127", "TMT127N"]],
    ids=["silac", "mixed_plex"],
)
def test_psm_conversion_does_not_map_labels(tmp_path, labels, with_evidence):
    """The PSM view reads no reporter columns, so labels the feature and pg views cannot map do not stop it.

    A SILAC label the SDRF declares as a fixed modification of the heavy
    samples is not added to the peptidoforms either: MaxQuant carries it on
    the labelled rows only.
    """
    from qpx.converters.maxquant.converter import MaxQuantConverter
    from qpx.core.constants import PSM

    evidence, _, _ = _write_inputs(tmp_path, TMT10, TMT10, 0)
    msms = tmp_path / "msms.txt"
    msms.write_text(_MSMS)
    heavy = "AC=UNIMOD:267;TA=R;MT=fixed"
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\tcomment[modification parameters]\n"
        + "".join(f"S_{label}\trun1.raw\t{label}\t{heavy if 'heavy' in label else ''}\n" for label in labels)
    )
    out = tmp_path / "out"
    MaxQuantConverter().convert(
        output_folder=str(out),
        msms_file=str(msms),
        evidence_file=str(evidence) if with_evidence else None,
        sdrf_file=str(sdrf),
        output_prefix="mq",
        structures=[PSM],
    )
    assert [row["peptidoform"] for row in pq.read_table(out / "mq.psm.parquet").to_pylist()] == ["PEPTIDER"]
    assert not (out / "mq.feature.parquet").exists()


def test_reporter_columns_without_sdrf_warn(tmp_path, caplog):
    """Without an SDRF the reporter columns cannot be mapped; say so instead of dropping them silently."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence, _, _ = _write_inputs(tmp_path, TMT10, TMT10, 0)
    with MaxQuantFeatureAdapter() as adapter, caplog.at_level("WARNING"):
        adapter.convert(str(evidence), str(tmp_path / "feature.parquet"))
    assert "10 reporter intensity columns" in caplog.text


@pytest.mark.parametrize(
    ("labels", "n_columns", "warns"),
    [(TMT6, 10, True), (["TMT126", "TMT131"], 11, True), (TMT10, 10, False), (ITRAQ4, 8, False)],
    ids=["tmt6_in_10_columns", "tmt6_in_11_columns", "tmt10", "itraq4_in_itraq8"],
)
def test_reporter_column_count_of_no_plex_warns(caplog, labels, n_columns, warns):
    """Labels of a plex with another channel count than the table's reporter columns are mapped with a warning."""
    from qpx.converters.maxquant.channels import reporter_channels

    with caplog.at_level("WARNING"):
        reporter_channels(labels, n_columns)
    assert (f"has {n_columns} reporter columns" in caplog.text) is warns
