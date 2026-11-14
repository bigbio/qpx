"""The label vocabulary: one spelling per channel, stated in the spec and checked by validate."""

import re
from pathlib import Path

import qpx
from qpx.core.data.label_validation import canonical_labels
from qpx.writers import FeatureWriter, PgWriter
from tests.conftest import make_feature_record, make_pg_record


def test_spec_table_lists_the_channel_map():
    text = (Path(__file__).parents[2] / "docs" / "spec" / "intensities.md").read_text(encoding="utf-8")
    table = text.split("<!-- label-vocabulary -->")[1].split("\n\n")[0]
    assert set(re.findall(r"`([^`]+)`", table)) == canonical_labels()


def _issues(path, structure):
    with qpx.open_dataset(str(path)) as dataset:
        return [issue for issue in dataset.validate([structure])[structure].issues if issue.check == "label_vocabulary"]


def test_canonical_labels_pass(tmp_path):
    with PgWriter(tmp_path / "x.pg.parquet") as writer:
        writer.write_batch([make_pg_record(intensities=[{"label": "ITRAQ114", "intensity": 1.0}])])
    assert not _issues(tmp_path, "pg")


def test_a_label_in_another_case_names_its_canonical_spelling(tmp_path):
    with PgWriter(tmp_path / "x.pg.parquet") as writer:
        writer.write_batch([make_pg_record(intensities=[{"label": "iTRAQ114", "intensity": 1.0}])])
    (issue,) = _issues(tmp_path, "pg")
    assert issue.severity == "warning"
    assert "iTRAQ114 (canonical: ITRAQ114)" in issue.message


def test_an_unknown_feature_label_is_reported(tmp_path):
    with FeatureWriter(tmp_path / "x.feature.parquet") as writer:
        writer.write_batch([{**make_feature_record(), "intensities": [{"label": "NeuCode1", "intensity": 1.0}]}])
    (issue,) = _issues(tmp_path, "feature")
    assert "NeuCode1" in issue.message
