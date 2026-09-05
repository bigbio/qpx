"""PSM.with_feature() / PSM.without_feature() filter recorded feature links."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.core.data import PSM


@pytest.fixture(name="psm_file")
def _psm_file(tmp_path):
    """A tiny psm.parquet with two linked rows and one unlinked row."""
    from tests.conftest import make_psm_record

    records = []
    for i, feature_id in enumerate((111, 222, None)):
        record = dict(make_psm_record())
        record["psm_id"] = 900 + i
        record["sequence"] = f"PEPTIDE{i}"
        record["feature_id"] = feature_id
        records.append(record)
    schema = PSM.schema()
    table = pa.Table.from_pylist([{k: r.get(k) for k in schema.names} for r in records], schema=schema)
    path = tmp_path / "x.psm.parquet"
    pq.write_table(table, path)
    return path


def test_with_feature_keeps_only_linked_rows(psm_file):
    """Filtering feature links leaves the original PSM view intact."""
    psm = PSM.from_file(psm_file, threads=24)
    assert psm.with_feature().count() == 2
    assert psm.with_feature().to_df()["feature_id"].notna().all()
    assert psm.count() == 3


def test_without_feature_is_the_complement(psm_file):
    """Linked and unlinked filters partition all PSM rows."""
    psm = PSM.from_file(psm_file, threads=24)
    assert psm.without_feature().count() == 1
    assert psm.without_feature().to_df()["feature_id"].isna().all()
    assert psm.with_feature().count() + psm.without_feature().count() == psm.count()
