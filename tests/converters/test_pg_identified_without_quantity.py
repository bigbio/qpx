"""Keep identified protein groups without fabricating unobserved unit rows."""

import csv

import pyarrow.parquet as pq
import pytest

from qpx.converters.fragpipe.pg_adapter import FragPipePgAdapter
from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter


def _write_tsv(path, rows):
    """Write small upstream reports through the public converters."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return str(path)


@pytest.mark.parametrize("missing_quantity", [0, None])
@pytest.mark.parametrize("count_column", ["Spectral Count", "Unique Spectral Count", "Total Spectral Count", None])
def test_fragpipe_retains_only_units_with_identification(tmp_path, missing_quantity, count_column):
    """Experiment spectral counts, not combined counts, support a null row."""
    identified = {
        "Protein": "P1",
        "Combined Total Peptides": 5,
        "Combined Unique Peptides": 3,
        "Protein FDR": 0.001,
        "exp1 Total Intensity": missing_quantity,
        "exp1 MaxLFQ Intensity": 0,
        "exp2 Total Intensity": missing_quantity,
        "exp2 MaxLFQ Intensity": 0,
    }
    if count_column:
        identified[f"exp1 {count_column}"] = 2
        identified[f"exp2 {count_column}"] = 0
    quantified = {**identified, "Protein": "P2", "exp1 Total Intensity": 100, "exp2 Total Intensity": 200}
    source = _write_tsv(tmp_path / "combined_protein.tsv", [identified, quantified])
    output = tmp_path / "fragpipe.pg.parquet"
    with FragPipePgAdapter(duckdb_threads=6) as adapter:
        adapter.convert(source, str(output), experiment_to_runs={"exp1": ["run1", "run2"], "exp2": ["run3"]})

    rows = pq.read_table(output).to_pylist()
    retained = [row for row in rows if row["anchor_protein"] == "P1"]
    assert len(retained) == int(count_column is not None)
    if retained:
        assert retained[0]["grouped_runs"] == ["run1", "run2"]
        assert retained[0]["label"] == "LFQ"
        assert retained[0]["intensity"] is None
        assert retained[0]["additional_intensities"] is None
        assert retained[0]["pg_qvalue"] == pytest.approx(0.001)
    assert sorted(row["intensity"] for row in rows if row["anchor_protein"] == "P2") == [100, 200]


@pytest.mark.parametrize("missing_quantity", [0, None])
@pytest.mark.parametrize("msms_count", [2, 0, None])
def test_maxquant_lfq_retains_only_units_with_identification(tmp_path, missing_quantity, msms_count):
    """A positive experiment MS/MS count retains its known LFQ label."""
    identified = {
        "Protein IDs": "P1",
        "Peptides": 5,
        "Unique peptides": 3,
        "MS/MS count": 10,
        "Q-value": 0.001,
        "Intensity exp1": missing_quantity,
        "Intensity exp2": missing_quantity,
        "LFQ intensity exp1": 0,
        "iBAQ exp1": 0,
    }
    if msms_count is not None:
        identified.update({"MS/MS count exp1": msms_count, "MS/MS count exp2": 0})
    quantified = {**identified, "Protein IDs": "P2", "Intensity exp1": 100, "Intensity exp2": 200}
    source = _write_tsv(tmp_path / "proteinGroups.txt", [identified, quantified])
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\n"
        "exp1\trun1.raw\tlabel free sample\n"
        "exp1\trun2.raw\tlabel free sample\n"
        "exp2\trun3.raw\tlabel free sample\n"
    )
    output = tmp_path / "maxquant.pg.parquet"
    with MaxQuantPgAdapter(duckdb_threads=6) as adapter:
        adapter.convert(source, str(output), sdrf_path=str(sdrf))

    rows = pq.read_table(output).to_pylist()
    retained = [row for row in rows if row["anchor_protein"] == "P1"]
    assert len(retained) == int(msms_count == 2)
    if retained:
        assert retained[0]["grouped_runs"] == ["run1", "run2"]
        assert retained[0]["label"] == "LFQ"
        assert retained[0]["intensity"] is None
        assert retained[0]["additional_intensities"] is None
        assert retained[0]["global_qvalue"] == pytest.approx(0.001)
    assert sorted(row["intensity"] for row in rows if row["anchor_protein"] == "P2") == [100, 200]


def test_maxquant_tmt_retains_declared_labels_in_identified_unit(tmp_path):
    """Unit evidence retains existing channels but cannot create absent ones."""
    identified = {
        "Protein IDs": "P1",
        "Q-value": 0.001,
        "MS/MS count exp1": 2,
        "MS/MS count exp2": 0,
        "Reporter intensity 0 exp1": 0,
        "Reporter intensity 1 exp1": None,
        "Reporter intensity 0 exp2": 0,
        "Reporter intensity 1 exp2": None,
    }
    quantified = {**identified, "Protein IDs": "P2", "MS/MS count exp1": 0, "Reporter intensity 0 exp1": 100}
    source = _write_tsv(tmp_path / "proteinGroups.txt", [identified, quantified])
    evidence = tmp_path / "evidence.txt"
    evidence.write_text("Experiment\tRaw file\nexp1\trun1\nexp1\trun2\nexp2\trun3\n")
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\n"
        "sample1\trun1.raw\tTMT126\n"
        "sample2\trun1.raw\tTMT127N\n"
        "sample3\trun1.raw\tTMT127C\n"
        "sample4\trun1.raw\tTMT128N\n"
        "sample5\trun1.raw\tTMT128C\n"
        "sample6\trun1.raw\tTMT129N\n"
        "sample7\trun1.raw\tTMT129C\n"
        "sample8\trun1.raw\tTMT130N\n"
        "sample9\trun1.raw\tTMT130C\n"
        "sample10\trun1.raw\tTMT131\n"
    )
    output = tmp_path / "maxquant.pg.parquet"
    with MaxQuantPgAdapter(duckdb_threads=6) as adapter:
        adapter.convert(source, str(output), sdrf_path=str(sdrf), evidence_path=str(evidence))

    rows = pq.read_table(output).to_pylist()
    retained = [row for row in rows if row["anchor_protein"] == "P1"]
    assert {row["label"] for row in retained} == {"TMT126", "TMT127N"}
    assert all(row["intensity"] is None and row["grouped_runs"] == ["run1", "run2"] for row in retained)
    assert [(row["label"], row["intensity"]) for row in rows if row["anchor_protein"] == "P2"] == [("TMT126", 100)]


def test_maxquant_tmt_identification_skips_channels_without_a_sample(tmp_path):
    """A partially used plex keeps null rows only for channels the SDRF assigns there."""
    labels = ["TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131"]
    row = {"Protein IDs": "P1", "Q-value": 0.001, "MS/MS count exp1": 3}
    for index in range(len(labels)):
        quantified = index < 5
        row[f"Reporter intensity {index} exp1"] = 100 + index if quantified else 0
        row[f"Reporter intensity corrected {index} exp1"] = 90 + index if quantified else 0
    source = _write_tsv(tmp_path / "proteinGroups.txt", [row])
    evidence = tmp_path / "evidence.txt"
    evidence.write_text("Experiment\tRaw file\nexp1\trun1\n")
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\n"
        + "".join(f"sample{index}\trun1.raw\t{label}\n" for index, label in enumerate(labels[:6]))
        + "".join(f"other{index}\trun2.raw\t{label}\n" for index, label in enumerate(labels))
    )
    output = tmp_path / "maxquant.pg.parquet"
    with MaxQuantPgAdapter(duckdb_threads=6) as adapter:
        adapter.convert(source, str(output), sdrf_path=str(sdrf), evidence_path=str(evidence))

    rows = {row["label"]: row for row in pq.read_table(output).to_pylist()}
    assert sorted(rows, key=labels.index) == labels[:6]
    assert rows["TMT129N"]["intensity"] is None
    assert rows["TMT129N"]["additional_intensities"] is None
    assert rows["TMT126"]["intensity"] == 100
