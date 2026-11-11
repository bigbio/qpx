"""The Producer Guide's dataset conforms, and each single defect in it is caught."""

import re
import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main

GUIDE = Path(__file__).resolve().parents[2] / "docs" / "spec" / "producer-guide.md"
_WRITER_BLOCK = re.compile(r"^```python\n(.*?)^```", re.S | re.M)


def _conforms(folder):
    return CliRunner().invoke(qpx_main, ["validate", "--dataset-path", str(folder), "--conformance"])


@pytest.fixture(scope="module")
def golden(tmp_path_factory, request):
    """The dataset the Producer Guide writes."""
    cwd = tmp_path_factory.mktemp("golden")
    monkeypatch = pytest.MonkeyPatch()
    request.addfinalizer(monkeypatch.undo)
    monkeypatch.chdir(cwd)
    code = _WRITER_BLOCK.search(GUIDE.read_text(encoding="utf-8")).group(1)
    exec(compile(code, "producer-guide", "exec"), {"__name__": "__guide__"})  # noqa: S102
    return cwd / "mytool_output"


def _column(table, name, values):
    index = table.schema.get_field_index(name)
    field = table.schema.field(index).with_nullable(True)
    return table.set_column(index, field, pa.array(values, type=field.type))


def _set(name, values):
    return lambda table: _column(table, name, values)


def _drop(name):
    return lambda table: table.drop_columns([name])


def _retype_charge(table):
    index = table.schema.get_field_index("charge")
    return table.set_column(index, pa.field("charge", pa.string()), pc.cast(table["charge"], pa.string()))


def _strip_footer_key(key):
    def mutate(table):
        metadata = {k: v for k, v in table.schema.metadata.items() if k != key}
        return table.replace_schema_metadata(metadata)

    return mutate


def _extra_column(table):
    return table.append_column("mystery", pa.array(["x"] * table.num_rows))


LFQ_ENTRY = [{"label": "LFQ", "intensity": 1000.0}]
DEFECTS = {
    "feature_missing_required_column": ("feature", _drop("sequence")),
    "feature_missing_nullable_column": ("feature", _drop("modifications")),
    "feature_retyped_column": ("feature", _retype_charge),
    "feature_null_id": ("feature", _set("feature_id", [None])),
    "feature_wrong_id": ("feature", _set("feature_id", [1])),
    "feature_duplicate_row": ("feature", lambda table: pa.concat_tables([table, table])),
    "feature_qvalue_above_one": ("feature", _set("precursor_qvalue", [1.5])),
    "feature_nan_intensity": ("feature", _set("intensities", [[{"label": "LFQ", "intensity": float("nan")}]])),
    "feature_negative_intensity": ("feature", _set("intensities", [[{"label": "LFQ", "intensity": -1.0}]])),
    "feature_noncanonical_label": ("feature", _set("intensities", [[{"label": "lfq", "intensity": 1.0}]])),
    "feature_unknown_label": ("feature", _set("intensities", [[{"label": "TMT126", "intensity": 1.0}]])),
    "feature_dangling_run": ("feature", _set("run_file_name", ["no_such_run"])),
    "feature_bad_proforma": ("feature", _set("peptidoform", ["PEPTM[UNIMOD:35IDEK"])),
    "feature_peptidoform_not_sequence": ("feature", _set("peptidoform", ["PEPTM[UNIMOD:35]IDER"])),
    "feature_noncanonical_mod": ("feature", _set("peptidoform", ["PEPTM[Oxidation]IDEK"])),
    "feature_wrong_calculated_mz": ("feature", _set("calculated_mz", [100.0])),
    "feature_far_observed_mz": ("feature", _set("observed_mz", [100.0])),
    "feature_anchor_outside_group": ("feature", _set("anchor_protein", ["Q99999"])),
    "feature_no_group": ("feature", _set("pg_accessions", [None])),
    "feature_undeclared_column": ("feature", _extra_column),
    "feature_no_footer_version": ("feature", _strip_footer_key(b"qpx_version")),
    "feature_no_footer_composite": ("feature", _strip_footer_key(b"identity_composite")),
    "psm_missing_required_column": ("psm", _drop("peptidoform")),
    "psm_null_id": ("psm", _set("psm_id", [None])),
    "psm_dangling_run": ("psm", _set("run_file_name", ["no_such_run"])),
    "psm_dangling_feature": ("psm", _set("feature_id", [1])),
    "psm_empty_scan": ("psm", _set("scan", [[]])),
    "psm_modification_past_the_end": (
        "psm",
        _set(
            "modifications",
            [
                [
                    {
                        "name": "Oxidation",
                        "accession": "UNIMOD:35",
                        "positions": [{"position": 50, "amino_acid": "M", "scores": None}],
                    }
                ]
            ],
        ),
    ),
    "psm_qvalue_negative": ("psm", _set("qvalue", [-0.1])),
    "pg_missing_required_column": ("pg", _drop("grouped_runs")),
    "pg_null_id": ("pg", _set("pg_id", [None])),
    "pg_wrong_id": ("pg", _set("pg_id", [1])),
    "pg_negative_intensity": ("pg", _set("intensity", [-5.0])),
    "pg_dangling_run": ("pg", _set("grouped_runs", [["no_such_run"]])),
    "pg_anchor_outside_group": ("pg", _set("anchor_protein", ["Q99999"])),
    "run_unknown_sample": (
        "run",
        _set("samples", [[{"sample_accession": "S9", "label": "LFQ", "biological_replicate": 1}]]),
    ),
    "run_missing_required_column": ("run", _drop("samples")),
    "sample_missing_required_column": ("sample", _drop("organism")),
    "dataset_missing_required_column": ("dataset", _drop("project_accession")),
}


def test_the_guide_dataset_conforms(golden):
    result = _conforms(golden)
    assert result.exit_code == 0, result.output


def test_relayouted_dataset_conforms_with_provided_ids(request, tmp_path):
    """Provided ordinal ids are informational for both dataset and single-file conformance."""
    source = request.getfixturevalue("golden")
    folder = tmp_path / "sorted"
    result = CliRunner().invoke(
        qpx_main,
        ["transform", "relayout", "--dataset", str(source), "--output-folder", str(folder)],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    for arguments in (["--dataset-path", str(folder)], ["--file", str(folder / "mytool.psm.parquet")]):
        result = CliRunner().invoke(qpx_main, ["validate", *arguments, "--conformance"])
        assert result.exit_code == 0, result.output
        assert "  INFO " in result.output
        assert "WARN " not in result.output


@pytest.mark.parametrize("name", DEFECTS)
def test_a_single_defect_fails_conformance(name, golden, tmp_path):
    view, mutate = DEFECTS[name]
    folder = tmp_path / "copy"
    shutil.copytree(golden, folder)
    path = folder / f"mytool.{view}.parquet"
    pq.write_table(mutate(pq.read_table(path)), path)
    result = _conforms(folder)
    assert isinstance(result.exception, SystemExit) and result.exit_code == 1, result.exception
