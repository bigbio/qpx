"""Plex-aware channel-label resolution for the OpenMS QPX paths.

Regression coverage for canonical TMT/iTRAQ labels shared by the OpenMS
consensusXML and ``-out_qpx`` paths. Labels come from the shared sdrf-pipelines
channel_map and are plex-aware.
"""

import sys

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.converters.channel_labels import (
    experiment_runs_from_sdrf,
    experiment_type_from_labels,
    normalize_label,
    relabel_intensities_parquet,
    resolve_channel_labels,
)


def test_experiment_runs_from_sdrf_groups_by_source_name(tmp_path):
    """SDRF 'source name' -> distinct raw-file stems, fraction order preserved."""
    sdrf = tmp_path / "x.sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[fraction identifier]\nS1\tF1.raw\t1\nS1\tF2.raw\t2\nS2\tG1.raw\t1\n"
    )
    mapping = experiment_runs_from_sdrf(str(sdrf))
    assert mapping == {"S1": ["F1", "F2"], "S2": ["G1"]}


def test_experiment_runs_from_sdrf_missing_or_absent():
    """None for no path; None when required columns are absent."""
    assert experiment_runs_from_sdrf(None) is None


def test_fraction_groups_from_sdrf_groups_fractions_not_techreps(tmp_path):
    """Files that differ only in fraction group into one quantification unit;
    technical replicates (differ in tech rep) stay separate."""
    from qpx.converters.channel_labels import fraction_groups_from_sdrf

    sdrf = tmp_path / "frac.sdrf.tsv"
    # S1 has 3 fractions (one unit); S1 tech-rep 2 is a separate unit; S2 unfractionated.
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\t"
        "characteristics[biological replicate]\tcomment[technical replicate]\t"
        "comment[fraction identifier]\n"
        "S1\tF1.raw\tLFQ\t1\t1\t1\n"
        "S1\tF2.raw\tLFQ\t1\t1\t2\n"
        "S1\tF3.raw\tLFQ\t1\t1\t3\n"
        "S1\tT2.raw\tLFQ\t1\t2\t1\n"
        "S2\tG1.raw\tLFQ\t1\t1\t1\n"
    )
    m = fraction_groups_from_sdrf(str(sdrf))
    assert m["F1"] == ["F1", "F2", "F3"]  # 3 fractions, fraction-ordered
    assert m["F2"] == ["F1", "F2", "F3"]
    assert m["T2"] == ["T2"]  # tech-rep 2 is its own unit
    assert m["G1"] == ["G1"]  # unfractionated
    assert fraction_groups_from_sdrf(None) is None


def test_multi_dot_run_names_match_sdrf_conversion(tmp_path):
    """SDRF mappings and run output preserve meaningful dots in run names."""
    from qpx.converters.channel_labels import fraction_groups_from_sdrf
    from qpx.converters.sdrf import SdrfConverter
    from qpx.core.sdrf import SDRFHandler

    sdrf = tmp_path / "multi_dot.sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\tcomment[fraction identifier]\n"
        "S1\trun.part.1.raw\tlabel free sample\t1\n"
        "S1\trun.part.2.raw\tlabel free sample\t2\n"
    )

    grouped = ["run.part.1", "run.part.2"]
    assert experiment_runs_from_sdrf(str(sdrf)) == {"S1": grouped}
    assert fraction_groups_from_sdrf(str(sdrf)) == {
        "run.part.1": grouped,
        "run.part.2": grouped,
    }
    assert SDRFHandler(sdrf).get_sample_map_run() == {
        "run.part.1": "S1",
        "run.part.2": "S1",
    }

    run_output = tmp_path / "run.parquet"
    with SdrfConverter(duckdb_threads=24) as converter:
        converter.convert(str(sdrf), run_output=str(run_output))
    assert pq.read_table(run_output).column("run_file_name").to_pylist() == grouped


@pytest.mark.parametrize("data_file", ["", "   "], ids=["missing", "blank"])
def test_sdrf_run_paths_reject_missing_data_files(tmp_path, data_file):
    """Missing run filenames must fail before producing mappings or output."""
    from qpx.converters.sdrf import SdrfConverter
    from qpx.core.sdrf import SDRFHandler

    sdrf = tmp_path / "missing_data_file.sdrf.tsv"
    sdrf.write_text(f"source name\tcomment[data file]\tcomment[label]\nS1\t{data_file}\tlabel free sample\n")
    error = r"comment\[data file\].*input row\(s\): 2"

    with pytest.raises(ValueError, match=error):
        SDRFHandler(sdrf).get_sample_map_run()

    run_output = tmp_path / "run.parquet"
    with SdrfConverter(duckdb_threads=24) as converter:
        with pytest.raises(ValueError, match=error):
            converter.convert(str(sdrf), run_output=str(run_output))
    assert not run_output.exists()


from qpx.core.parquet_io import read_parquet_metadata
from qpx.dataset import Dataset
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from tests.conftest import make_feature_record, make_pg_record, write_lfq_consensusxml


def _resolve(experiment_type, channel_indices, sdrf_labels=None):
    return resolve_channel_labels(experiment_type, sdrf_labels, channel_indices)


def test_tmt10_channel10_is_tmt131_not_tmt131n():
    labels = _resolve("TMT", list(range(1, 11)))  # 10 channels -> tmt10plex
    assert labels[1] == "TMT126"
    assert labels[10] == "TMT131"  # NOT TMT131N (the historical bug)
    assert "TMT131N" not in labels.values()


def test_tmt11_uses_131n_131c():
    labels = _resolve("TMT", list(range(1, 12)))  # 11 channels -> tmt11plex
    assert labels[10] == "TMT131N"
    assert labels[11] == "TMT131C"


def test_experiment_type_from_labels():
    assert experiment_type_from_labels({"TMT126", "TMT131C"}) == "TMT"
    assert experiment_type_from_labels({"ITRAQ114"}) == "iTRAQ"
    assert experiment_type_from_labels({"label free sample"}) == "LFQ"
    assert experiment_type_from_labels(None) == "LFQ"


@pytest.mark.parametrize(
    ("raw_label", "expected"),
    [
        ("tmt126", "TMT126"),
        ("TmT127n", "TMT127N"),
        ("iTRAQ114", "ITRAQ114"),
        ("itraq114", "ITRAQ114"),
        ("NT=tmt126;AC=MS:1002038", "TMT126"),
        ("AC=MS:1002038;nt=TmT127n", "TMT127N"),
        ("AC=MS:1002038;NT=label free sample", "LFQ"),
    ],
)
def test_normalize_label_case_and_ontology_form(raw_label, expected):
    """Known label spellings and ontology forms normalize to one join key."""
    assert normalize_label(raw_label) == expected


def _write_intensities_parquet(path, rows):
    itype = pa.list_(pa.struct([("label", pa.string()), ("intensity", pa.float32())]))
    pq.write_table(pa.table({"intensities": pa.array(rows, type=itype)}), path)


def test_relabel_openms_feature_filename_labels(tmp_path):
    """OpenMS feature -out_qpx: filename in every label -> reporter names by position."""
    src = tmp_path / "in.parquet"
    dst = tmp_path / "out.parquet"
    fname = "run01.mzML"
    _write_intensities_parquet(src, [[{"label": fname, "intensity": float(i)} for i in range(1, 11)]])
    labels = resolve_channel_labels("TMT", {f"TMT12{i}" for i in range(6)} | {"TMT131"}, range(1, 11))
    relabel_intensities_parquet(str(src), str(dst), labels, is_lfq=False)
    got = [e["label"] for e in pq.read_table(dst).column("intensities").to_pylist()[0]]
    assert got[0] == "TMT126"
    assert got[9] == "TMT131"  # tmt10plex, position 10
    assert fname not in got


def test_relabel_openms_pg_index_labels(tmp_path):
    """OpenMS pg -out_qpx: bare index labels '1','2',... -> reporter names."""
    src = tmp_path / "pg_in.parquet"
    dst = tmp_path / "pg_out.parquet"
    _write_intensities_parquet(src, [[{"label": str(i), "intensity": 1.0} for i in range(1, 12)]])
    labels = resolve_channel_labels("TMT", None, range(1, 12))  # tmt11plex by count
    relabel_intensities_parquet(str(src), str(dst), labels, is_lfq=False)
    got = [e["label"] for e in pq.read_table(dst).column("intensities").to_pylist()[0]]
    assert got == [
        "TMT126",
        "TMT127N",
        "TMT127C",
        "TMT128N",
        "TMT128C",
        "TMT129N",
        "TMT129C",
        "TMT130N",
        "TMT130C",
        "TMT131N",
        "TMT131C",
    ]


def test_relabel_preserves_qpx_compression_and_float_encoding(tmp_path):
    """Relabeling must not silently rewrite ZSTD/BYTE_STREAM_SPLIT as Snappy."""
    src = tmp_path / "feature_in.parquet"
    dst = tmp_path / "feature_out.parquet"
    with FeatureWriter(src, compression="zstd") as writer:
        writer.write_batch([make_feature_record()])

    relabel_intensities_parquet(
        str(src),
        str(dst),
        {1: "TMT126"},
        is_lfq=False,
        compression="zstd",
    )

    assert read_parquet_metadata(dst)["compression_format"] == "zstd"
    parquet = pq.ParquetFile(dst)
    row_group = parquet.metadata.row_group(0)

    def _col(path):
        return next(
            row_group.column(index) for index in range(row_group.num_columns) if row_group.column(index).path_in_schema == path
        )

    # ZSTD compression is preserved (not silently downgraded to Snappy).
    intensity_column = _col("intensities.list.element.intensity")
    assert intensity_column.compression == "ZSTD"
    # Selective BYTE_STREAM_SPLIT is preserved through relabeling: rt (a BSS
    # leaf) keeps it; the intensity leaves were deliberately dropped from BSS
    # (they regress TMT reporter columns), so intensity no longer carries it.
    assert "BYTE_STREAM_SPLIT" in _col("rt").encodings
    assert "BYTE_STREAM_SPLIT" not in intensity_column.encodings


def test_run_label_is_canonical_join_key(tmp_path):
    """run.samples[].label must be the SAME canonical label that feature/pg carry:
    the SDRF ontology form 'label free sample' is normalized to 'LFQ'."""
    import pyarrow.parquet as pq

    from qpx.converters.sdrf import SdrfConverter

    sdrf = "tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv"
    with SdrfConverter() as conv:
        conv.convert(
            sdrf_path=sdrf,
            sample_output=str(tmp_path / "s.parquet"),
            run_output=str(tmp_path / "r.parquet"),
        )
    run = pq.read_table(tmp_path / "r.parquet")
    labels = {s["label"] for row in run.column("samples").to_pylist() for s in (row or [])}
    assert labels == {"LFQ"}


def _write_label_sdrf(path, labels):
    header = [
        "source name",
        "characteristics[organism]",
        "characteristics[organism part]",
        "comment[data file]",
        "comment[label]",
    ]
    rows = [[f"sample_{index}", "Homo sapiens", "liver", "run_01.raw", label] for index, label in enumerate(labels, start=1)]
    path.write_text(
        "\n".join("\t".join(row) for row in [header, *rows]) + "\n",
        encoding="utf-8",
    )


def test_run_labels_normalize_case_and_ontology_form(tmp_path):
    """Mixed-case SDRF labels must join canonical feature and PG labels."""
    from qpx.converters.sdrf import SdrfConverter

    sdrf_path = tmp_path / "mixed.sdrf.tsv"
    _write_label_sdrf(sdrf_path, ["tmt126", "AC=MS:1002038;nt=TmT127n"])

    with SdrfConverter(duckdb_threads=24) as converter:
        converter.convert(
            sdrf_path=str(sdrf_path),
            sample_output=str(tmp_path / "x.sample.parquet"),
            run_output=str(tmp_path / "x.run.parquet"),
        )

    samples = pq.read_table(tmp_path / "x.run.parquet").column("samples").to_pylist()[0]
    assert [sample["label"] for sample in samples] == ["TMT126", "TMT127N"]

    intensities = [
        {"label": "TMT126", "intensity": 100.0},
        {"label": "TMT127N", "intensity": 200.0},
    ]
    with FeatureWriter(tmp_path / "x.feature.parquet") as writer:
        writer.write_batch([make_feature_record(intensities=intensities)])
    with PgWriter(tmp_path / "x.pg.parquet") as writer:
        writer.write_batch([make_pg_record(intensities=intensities)])

    with Dataset(tmp_path, file_prefix="x", duckdb_threads=24) as dataset:
        assert len(dataset.intensity("peptide").to_df()) == 2
        assert len(dataset.intensity("protein").to_df()) == 2


def test_run_rejects_duplicate_canonical_labels(tmp_path):
    """One run cannot map two samples to the same canonical channel."""
    from qpx.converters.sdrf import SdrfConverter

    sdrf_path = tmp_path / "duplicate.sdrf.tsv"
    _write_label_sdrf(sdrf_path, ["TMT126", "tmt126"])

    with SdrfConverter(duckdb_threads=24) as converter:
        with pytest.raises(ValueError, match="maps canonical label 'TMT126' more than once"):
            converter.convert(
                sdrf_path=str(sdrf_path),
                sample_output=str(tmp_path / "sample.parquet"),
                run_output=str(tmp_path / "run.parquet"),
            )


def test_relabel_lfq_collapses_to_lfq(tmp_path):
    src = tmp_path / "lfq_in.parquet"
    dst = tmp_path / "lfq_out.parquet"
    _write_intensities_parquet(src, [[{"label": "run01.mzML", "intensity": 5.0}]])
    relabel_intensities_parquet(str(src), str(dst), {}, is_lfq=True)
    assert pq.read_table(dst).column("intensities").to_pylist()[0][0]["label"] == "LFQ"


def _write_consensusxml(path, n_channels):
    """Write a minimal isobaric consensusXML with n_channel ColumnHeaders."""
    import pyopenms as oms

    cmap = oms.ConsensusMap()
    headers = {}
    for i in range(n_channels):
        h = oms.ColumnHeader()
        h.filename = "run01.mzML"
        h.label = str(i + 1)  # OpenMS puts a bare channel id here for isobaric
        headers[i] = h
    cmap.setColumnHeaders(headers)
    cf = oms.ConsensusFeature()
    cf.setRT(100.0)
    cf.setMZ(500.0)
    cmap.push_back(cf)
    oms.ConsensusXMLFile().store(str(path), cmap)


def test_channel_labels_from_consensusxml_tmt11(tmp_path, monkeypatch):
    """consensusXML ColumnHeaders give the authoritative channel count/order;
    canonical labels come from the SDRF-declared plex without loading the map."""
    cxml = tmp_path / "x.consensusXML"
    _write_consensusxml(cxml, 11)
    monkeypatch.setitem(sys.modules, "pyopenms", None)
    tmt11 = {
        "TMT126",
        "TMT127N",
        "TMT127C",
        "TMT128N",
        "TMT128C",
        "TMT129N",
        "TMT129C",
        "TMT130N",
        "TMT130C",
        "TMT131N",
        "TMT131C",
    }
    from qpx.converters.channel_labels import channel_labels_from_consensusxml

    labels = channel_labels_from_consensusxml(str(cxml), "TMT", tmt11)
    assert len(labels) == 11
    assert labels[1] == "TMT126"
    assert labels[10] == "TMT131N"
    assert labels[11] == "TMT131C"


def test_channel_labels_from_consensusxml_stops_after_map_list(tmp_path):
    """The parser must not consume consensus features after reading channel headers."""
    cxml = tmp_path / "header_only.consensusXML"
    maps = "\n".join(f'<map id="{index}" name="run01.mzML" label="{index + 1}" />' for index in range(6))
    cxml.write_text(
        f'<consensusXML><mapList count="6">{maps}</mapList>' + (" " * 100_000) + "<malformed",
        encoding="utf-8",
    )

    from qpx.converters.channel_labels import channel_labels_from_consensusxml

    labels = channel_labels_from_consensusxml(str(cxml), "TMT")
    assert labels[1] == "TMT126"
    assert labels[6] == "TMT131"


def test_channel_labels_from_consensusxml_missing_file():
    from qpx.converters.channel_labels import channel_labels_from_consensusxml

    assert channel_labels_from_consensusxml("/no/such.consensusXML", "TMT", {"TMT126"}) == {}


# ---------------------------------------------------------------------------
# fraction_groups_from_consensusxml (LFQ experimental-design capture)
# ---------------------------------------------------------------------------


def test_fraction_groups_from_consensusxml_lfq(tmp_path):
    from qpx.converters.channel_labels import fraction_groups_from_consensusxml

    cxml = tmp_path / "lfq.consensusXML"
    write_lfq_consensusxml(
        cxml,
        [
            (0, "run_A_f1.mzML", "1", "1", "1"),
            (1, "run_A_f2.mzML", "1", "2", "1"),
            (2, "run_B_f1.mzML", "2", "1", "2"),
        ],
    )
    groups = fraction_groups_from_consensusxml(str(cxml))
    assert set(groups) == {"run_A_f1.mzML", "run_A_f2.mzML", "run_B_f1.mzML"}
    assert groups["run_A_f1.mzML"] == {"fraction_group": "1", "fraction": "1", "sample_name": "1"}
    # Two fractions of the same sample share a fraction_group.
    assert groups["run_A_f2.mzML"]["fraction_group"] == "1"
    assert groups["run_B_f1.mzML"]["fraction_group"] == "2"


def test_fraction_groups_from_consensusxml_stops_after_map_list(tmp_path):
    """Parsing must stop at </mapList> and not choke on trailing junk."""
    from qpx.converters.channel_labels import fraction_groups_from_consensusxml

    cxml = tmp_path / "trailing.consensusXML"
    write_lfq_consensusxml(
        cxml,
        [(0, "run01.mzML", "1", "1", "1")],
        trailing=(" " * 100_000) + "<malformed",
    )
    groups = fraction_groups_from_consensusxml(str(cxml))
    assert groups["run01.mzML"]["fraction_group"] == "1"


def test_fraction_groups_from_consensusxml_isobaric_returns_empty(tmp_path):
    """Isobaric / pre-design maps have no fraction_group UserParam -> {}."""
    from qpx.converters.channel_labels import fraction_groups_from_consensusxml

    cxml = tmp_path / "isobaric.consensusXML"
    maps = "".join(f'<map id="{i}" name="run01.mzML" label="{i + 1}" />' for i in range(6))
    cxml.write_text(f'<consensusXML><mapList count="6">{maps}</mapList>', encoding="utf-8")
    assert fraction_groups_from_consensusxml(str(cxml)) == {}


def test_fraction_groups_from_consensusxml_missing_file():
    from qpx.converters.channel_labels import fraction_groups_from_consensusxml

    assert fraction_groups_from_consensusxml("/no/such.consensusXML") == {}


@pytest.mark.parametrize(
    ("exp_type", "n", "last_label"),
    [
        ("TMT", 6, "TMT131"),
        ("TMT", 16, "TMT134N"),
        ("TMT", 18, "TMT135N"),
        ("iTRAQ", 4, "ITRAQ117"),
        ("iTRAQ", 8, "ITRAQ121"),
    ],
)
def test_plex_resolution_by_channel_count(exp_type, n, last_label):
    labels = _resolve(exp_type, list(range(1, n + 1)))
    assert len(labels) == n
    assert labels[n] == last_label


def test_itraq_is_now_mapped():
    labels = _resolve("iTRAQ", [1, 2, 3, 4])
    assert labels == {1: "ITRAQ114", 2: "ITRAQ115", 3: "ITRAQ116", 4: "ITRAQ117"}


def test_non_isobaric_returns_empty():
    assert _resolve("LFQ", [1]) == {}
    assert _resolve("TMT", []) == {}


def test_sdrf_labels_override_channel_count():
    """SDRF ground truth wins: a TMT11 SDRF with an empty 11th channel in the
    data still resolves to tmt11plex (ch10 = TMT131N), where count-inference
    from indices 1..10 alone would wrongly pick tmt10plex (ch10 = TMT131)."""
    tmt11_labels = {
        "TMT126",
        "TMT127N",
        "TMT127C",
        "TMT128N",
        "TMT128C",
        "TMT129N",
        "TMT129C",
        "TMT130N",
        "TMT130C",
        "TMT131N",
        "TMT131C",
    }
    labels = _resolve("TMT", list(range(1, 11)), sdrf_labels=tmt11_labels)
    assert labels[10] == "TMT131N"  # SDRF-declared tmt11plex, not tmt10plex
    # without the SDRF, count-inference falls back to tmt10plex
    assert _resolve("TMT", list(range(1, 11)))[10] == "TMT131"


def test_unresolvable_sdrf_labels_fall_back_to_count():
    labels = _resolve("TMT", list(range(1, 11)), sdrf_labels={"garbage", "nonsense"})
    assert labels[10] == "TMT131"  # fell back to count-based tmt10plex
