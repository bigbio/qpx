"""Conversions record their options and SDRF, and the ontology view describes what was written."""

import hashlib

import click
import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.cli.transform.quantify import _condition_column
from qpx.converters.ontology_util import reconcile_ontology_entries
from qpx.converters.sdrf import SdrfConverter
from qpx.core.ontology import PublicOntology
from qpx.core.scores import _FIELD_CV_MAP, field_ontology_entries
from qpx.dataset import Dataset
from qpx.writers.ontology import OntologyWriter
from tests.converters.test_spectronaut_report_layout import _pivot_report

SDRF_HEADER = [
    "source name",
    "comment[data file]",
    "comment[label]",
    "characteristics[organism]",
    "characteristics[organism part]",
]


def _write_sdrf(path, rows, extra_columns=()):
    lines = ["\t".join([*SDRF_HEADER, *extra_columns]), *("\t".join(row) for row in rows)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _sample_rows(tmp_path, sdrf):
    sample = tmp_path / "s.sample.parquet"
    run = tmp_path / "s.run.parquet"
    with SdrfConverter() as converter:
        converter.convert(str(sdrf), sample_output=str(sample), run_output=str(run))
    return pq.read_table(sample).to_pylist(), pq.read_table(run).to_pylist()


def _write_view(folder, view, table):
    pq.write_table(table, folder / f"p.{view}.parquet")


def _names(entries):
    return {(entry["field_name"], entry["view"]) for entry in entries}


@pytest.mark.parametrize("field_name", sorted(_FIELD_CV_MAP))
def test_field_terms_match_the_vocabulary(field_name):
    """Each mapped accession carries the name PSI-MS gives it."""
    info = _FIELD_CV_MAP[field_name]
    term = PublicOntology("psi_ms", auto_update=False).lookup(info["ontology_accession"])

    assert term is not None
    assert term["name"] == info["ontology_name"]


def test_rows_for_names_never_written_are_dropped(tmp_path):
    _write_view(tmp_path, "feature", pa.table({"charge": [2]}))
    entries = field_ontology_entries("feature", {"charge": "Charge", "cscore": "CScore"}, "Tool")
    entries += field_ontology_entries("psm", {"charge": "Charge"}, "Tool")

    assert _names(reconcile_ontology_entries(entries, tmp_path, "p")) == {("charge", "feature")}


def test_emitted_names_without_a_row_get_one(tmp_path):
    score = pa.struct([("score_name", pa.string()), ("score_value", pa.float32())])
    cv_param = pa.struct([("cv_name", pa.string()), ("cv_value", pa.string())])
    pair = pa.struct([("intensity_name", pa.string()), ("intensity_value", pa.float32())])
    additional = pa.struct([("label", pa.string()), ("intensities", pa.list_(pair))])
    table = pa.table(
        {
            "charge": [2],
            "unmapped_column": [1],
            "additional_scores": pa.array(
                [[{"score_name": "andromeda_score", "score_value": 1.0}, {"score_name": "tool_private", "score_value": 2.0}]],
                pa.list_(score),
            ),
            "cv_params": pa.array([[{"cv_name": "custom_cv", "cv_value": "x"}]], pa.list_(cv_param)),
            "additional_intensities": pa.array(
                [[{"label": "LFQ", "intensities": [{"intensity_name": "normalize_intensity", "intensity_value": 1.0}]}]],
                pa.list_(additional),
            ),
        }
    )
    _write_view(tmp_path, "pg", table)

    rows = {entry["field_name"]: entry for entry in reconcile_ontology_entries([], tmp_path, "p")}

    assert set(rows) == {"andromeda_score", "charge", "custom_cv", "normalize_intensity", "tool_private"}
    assert rows["andromeda_score"]["ontology_accession"] == "MS:1002338"
    assert rows["charge"]["ontology_accession"] == "MS:1000041"
    assert rows["tool_private"]["ontology_accession"] is None


def test_validate_warns_on_ontology_rows_naming_nothing(tmp_path):
    _write_view(tmp_path, "pg", pa.table({"pg_id": pa.array([1], pa.int64())}))
    entries = field_ontology_entries("pg", {"pg_id": "id", "cscore": "CScore"}, "Tool")
    with OntologyWriter(tmp_path / "p.ontology.parquet") as writer:
        writer.write_batch(entries)

    with Dataset(str(tmp_path), file_prefix="p") as ds:
        issues = ds.validate(["ontology", "pg"])["ontology"].issues

    assert [issue.check for issue in issues if "cscore" in issue.message] == ["ontology_unknown_field"]


def test_factors_become_sample_columns(tmp_path):
    sdrf = _write_sdrf(
        tmp_path / "in.tsv",
        [
            ["S1", "a.raw", "label free sample", "Human", "liver", "cancer", "10 fmol"],
            ["S2", "b.raw", "label free sample", "Human", "liver", "normal", "20 fmol"],
        ],
        ("factor value[disease]", "factor value[spike-in amount]"),
    )

    samples, _ = _sample_rows(tmp_path, sdrf)

    assert [row["factor_disease"] for row in samples] == ["cancer", "normal"]
    assert [row["factor_spike_in_amount"] for row in samples] == ["10 fmol", "20 fmol"]


def test_a_blank_organism_keeps_the_sample_and_run_views(tmp_path):
    sdrf = _write_sdrf(
        tmp_path / "in.tsv",
        [
            ["S1", "a.raw", "label free sample", "", "liver"],
            ["S2", "b.raw", "label free sample", "Human", "liver"],
        ],
    )

    samples, runs = _sample_rows(tmp_path, sdrf)

    assert [row["organism"] for row in samples] == ["not available", "Human"]
    assert len(runs) == 2


def test_condition_column_names_the_chosen_factor():
    columns = ["source name", "factor value[serum]", "factor value[ups1]"]

    assert _condition_column(columns, "UPS1") == "factor value[ups1]"
    assert _condition_column(columns, "factor value[serum]") == "factor value[serum]"
    assert _condition_column(["source name"], None) is None
    with pytest.raises(click.UsageError, match="serum"):
        _condition_column(columns, "disease")


def test_several_factors_without_a_choice_warn(capsys):
    assert _condition_column(["factor value[serum]", "factor value[ups1]"], None) == "factor value[serum]"
    assert "--factor" in capsys.readouterr().err


def test_conversion_records_options_checksums_and_the_sdrf(tmp_path):
    report = _pivot_report(tmp_path / "pivot.tsv", {"FG.Quantity": ["1000.5", "2000.5"]})
    sdrf = _write_sdrf(
        tmp_path / "in.sdrf.tsv",
        [
            ["S1", "runA.raw", "label free sample", "Human", "liver"],
            ["S2", "runB.raw", "label free sample", "Human", "liver"],
        ],
    )
    out = tmp_path / "out"

    result = CliRunner().invoke(
        qpx_main,
        ["convert", "spectronaut", "--report-path", str(report), "--sdrf-file", str(sdrf), "--output-folder", str(out)]
        + ["--output-prefix", "sn", "--qvalue-threshold", "0.01"],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    steps = pq.read_table(out / "sn.provenance.parquet").to_pylist()
    conversion = next(step for step in steps if step["step_category"] == "format_conversion")
    parameters = {parameter["key"]: parameter["value"] for parameter in conversion["parameters"]}
    assert parameters["qvalue_threshold"] == "0.01"
    assert parameters["report_path"] == "pivot.tsv"
    assert parameters["report_path_sha256"] == hashlib.sha256(report.read_bytes()).hexdigest()
    assert parameters["sdrf_file_sha256"] == hashlib.sha256(sdrf.read_bytes()).hexdigest()
    assert parameters["dependency_duckdb"] == duckdb.__version__
    assert parameters["dependency_pyarrow"] == pa.__version__
    assert "output_folder" not in parameters
    assert (out / "sn.sdrf.tsv").read_bytes() == sdrf.read_bytes()


def test_an_unset_threshold_is_recorded(tmp_path):
    report = _pivot_report(tmp_path / "pivot.tsv", {"FG.Quantity": ["1000.5", "2000.5"]})
    out = tmp_path / "out"

    result = CliRunner().invoke(
        qpx_main,
        ["convert", "spectronaut", "--report-path", str(report), "--output-folder", str(out), "--output-prefix", "sn"],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    steps = pq.read_table(out / "sn.provenance.parquet").to_pylist()
    conversion = next(step for step in steps if step["step_category"] == "format_conversion")
    assert {"key": "qvalue_threshold", "value": "none"} in conversion["parameters"]
