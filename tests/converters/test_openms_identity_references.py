"""Native OpenMS enrichment must preserve existing QPX identity references."""

import pyarrow.parquet as pq
import pytest

from qpx.converters.openms.converter import OpenMSConverter
from qpx.dataset import Dataset
from qpx.writers import FeatureWriter, PgWriter, PsmWriter
from tests.conftest import make_feature_record, make_pg_record, make_psm_record


@pytest.mark.parametrize("identified", [True, False])
def test_openms_preserves_producer_ids_and_references(tmp_path, identified):
    """Rewriting a bundle retains both the authoritative FK and optional inverses."""
    source = tmp_path / "source"
    source.mkdir()
    feature = make_feature_record()
    feature.update(feature_id=101, psm_ids=[201], pg_ids=[301])
    if not identified:
        feature.update(sequence="", peptidoform="")
    psm = make_psm_record()
    psm.update(psm_id=201, feature_id=101)
    pg = make_pg_record(pg_accessions=["P12345"])
    pg["pg_id"] = 301

    with FeatureWriter(source / "source.feature.parquet", override_provided_ids=False) as writer:
        writer.write_batch([feature])
    with PsmWriter(source / "source.psm.parquet") as writer:
        writer.write_batch([psm])
    with PgWriter(source / "source.pg.parquet") as writer:
        writer.write_batch([pg])

    output = tmp_path / "output"
    OpenMSConverter(qpx_dir=source).convert(output_folder=output, output_prefix="openms")

    written_feature = pq.read_table(output / "openms.feature.parquet").to_pylist()[0]
    written_psm = pq.read_table(output / "openms.psm.parquet").to_pylist()[0]
    written_pg = pq.read_table(output / "openms.pg.parquet").to_pylist()[0]
    assert written_feature["feature_id"] == written_psm["feature_id"] == 101
    assert written_feature["psm_ids"] == [written_psm["psm_id"]] == [201]
    assert written_feature["pg_ids"] == [written_pg["pg_id"]] == [301]
    assert written_feature["cv_params"] == feature["cv_params"]

    with Dataset(output, file_prefix="openms", structures=["feature", "psm", "pg"], duckdb_threads=6) as dataset:
        results = dataset.validate(strict=True)
        assert all(result.is_valid for result in results.values()), {name: result.summary for name, result in results.items()}
        assert dataset.link_feature_psm().fetchall() == [(101, 201)]
        assert dataset.link_feature_pg().fetchall() == [(101, 301, "TMT126")]
