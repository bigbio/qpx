"""Sanctioned extension columns, and what happens to any other column outside the schema."""

import logging

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.core.data import PsmSchema
from qpx.dataset import Dataset
from qpx.writers import PsmWriter
from qpx.writers.ontology import OntologyWriter
from tests.conftest import make_psm_record


def test_extension_columns_are_written_and_read_back(tmp_path):
    out = tmp_path / "x.psm.parquet"
    with PsmWriter(out, extra_columns=["x_mytool_cluster"]) as writer:
        writer.write_batch([{**make_psm_record(), "x_mytool_cluster": "7"}])

    assert pq.read_table(out, columns=["x_mytool_cluster"]).to_pylist() == [{"x_mytool_cluster": "7"}]


def test_another_extra_column_is_refused_by_name():
    with pytest.raises(ValueError, match=r"\['cluster'\].*x_<namespace>_<name>"):
        PsmSchema.get_extended_arrow_schema(["cluster"])


def test_a_record_key_outside_the_schema_names_the_way_to_keep_it(tmp_path):
    with PsmWriter(tmp_path / "x.psm.parquet") as writer:
        with pytest.raises(ValueError, match="x_<namespace>_<name>"):
            writer.write_batch([{**make_psm_record(), "cluster": "7"}])


def test_aligning_a_table_warns_about_columns_it_drops(tmp_path, caplog):
    with PsmWriter(tmp_path / "x.psm.parquet") as writer:
        table = pa.table({"cluster": ["7"], "_hint": ["1"], "charge": [2]})
        with caplog.at_level(logging.WARNING):
            writer.align_table_to_schema(table)
        writer.write_batch([make_psm_record()])

    assert "['cluster']" in caplog.text


def test_a_dataframe_warns_about_columns_it_drops(tmp_path, caplog):
    frame = pd.DataFrame({"cluster": ["7"], **{key: [value] for key, value in make_psm_record().items()}})
    with caplog.at_level(logging.WARNING):
        with PsmWriter(tmp_path / "x.psm.parquet") as writer:
            writer.write_dataframe(frame)

    assert "['cluster']" in caplog.text


def _validated_warnings(tmp_path, columns):
    with PsmWriter(tmp_path / "x.psm.parquet") as writer:
        writer.write_batch([make_psm_record()])
    table = pq.read_table(tmp_path / "x.psm.parquet")
    for name in columns:
        table = table.append_column(name, pa.array(["7"]))
    pq.write_table(table, tmp_path / "x.psm.parquet")
    with OntologyWriter(tmp_path / "x.ontology.parquet") as writer:
        writer.write_batch([{"field_name": "charge", "view": "psm"}])
    with Dataset(str(tmp_path), file_prefix="x") as dataset:
        return [issue.message for issue in dataset.validate()["psm"].warnings if "outside the schema" in issue.message]


def test_validate_warns_about_an_undeclared_column_only(tmp_path):
    assert _validated_warnings(tmp_path, ["x_mytool_cluster"]) == []
    assert "cluster" in _validated_warnings(tmp_path, ["cluster"])[0]
