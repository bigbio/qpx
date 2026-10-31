"""Unit tests for the computed feature<->psm inverse softlink (bigbio/qpx#267).

``psm.feature_id`` is the authoritative (optional) foreign key and stays
materialized; its inverse, ``feature.psm_ids``, is NOT materialized — it is
computed on read by ``Dataset.link_feature_psm()`` (group the persisted
``psm.feature_id`` by feature). These tests build tiny synthetic datasets
directly so the inverse and its complement can be asserted without any converter.
"""

from __future__ import annotations

from qpx.core.data.identity import derive_id
from qpx.dataset import Dataset
from qpx.writers import FeatureWriter, PsmWriter
from tests.conftest import make_feature_record, make_psm_record


def _write(ds_dir, features, psms):
    ds_dir.mkdir(parents=True, exist_ok=True)
    with FeatureWriter(ds_dir / "exp.feature.parquet") as w:
        w.write_batch(features)
    with PsmWriter(ds_dir / "exp.psm.parquet") as w:
        w.write_batch(psms)


def test_link_feature_psm_is_inverse_of_feature_id(tmp_path):
    """link_feature_psm() returns (feature_id, psm_id) grouped from psm.feature_id;
    psms with a null feature_id are the complement (psms_without_feature())."""
    ds_dir = tmp_path / "inverse"
    feature = make_feature_record(sequence="FEATUREK", peptidoform="FEATUREK")
    feature["rt"] = 10.0
    feature_id = derive_id([feature["peptidoform"], feature["charge"], feature["run_file_name"], feature["rt"]])

    # Two psms assigned to the feature, one psm left unassigned (null feature_id).
    psm_a = make_psm_record(sequence="PEPTIDEAK", peptidoform="PEPTIDEAK", scan=[10])
    psm_b = make_psm_record(sequence="PEPTIDEBK", peptidoform="PEPTIDEBK", scan=[20])
    psm_c = make_psm_record(sequence="PEPTIDECK", peptidoform="PEPTIDECK", scan=[30])
    psm_a["feature_id"] = feature_id
    psm_b["feature_id"] = feature_id
    psm_c["feature_id"] = None
    id_a = derive_id([psm_a["run_file_name"], psm_a["scan"], psm_a["peptidoform"], psm_a["charge"]])
    id_b = derive_id([psm_b["run_file_name"], psm_b["scan"], psm_b["peptidoform"], psm_b["charge"]])
    id_c = derive_id([psm_c["run_file_name"], psm_c["scan"], psm_c["peptidoform"], psm_c["charge"]])

    _write(ds_dir, [feature], [psm_a, psm_b, psm_c])

    with Dataset(ds_dir) as ds:
        link = ds.link_feature_psm().fetchall()
        unassigned = {row[0] for row in ds.psms_without_feature().fetchall()}
        # feature.psm_ids is NOT materialized by qpx.
        psm_ids_col = ds.sql("SELECT psm_ids FROM feature").fetchall()

    assert all(ids is None for (ids,) in psm_ids_col), "feature.psm_ids must not be materialized"

    # Grouping the inverse by feature_id recovers the feature's psm list.
    recovered: dict[int, set[int]] = {}
    for fid, pid in link:
        recovered.setdefault(fid, set()).add(pid)
    assert recovered == {feature_id: {id_a, id_b}}
    # The unassigned psm is the complement.
    assert unassigned == {id_c}
    assert id_a not in unassigned and id_b not in unassigned


def test_link_feature_psm_requires_both_views(tmp_path):
    """The softlink accessor needs both the feature and psm views registered."""
    import pytest

    ds_dir = tmp_path / "featureonly"
    feature = make_feature_record()
    ds_dir.mkdir(parents=True, exist_ok=True)
    with FeatureWriter(ds_dir / "exp.feature.parquet") as w:
        w.write_batch([feature])

    with Dataset(ds_dir, structures=["feature"]) as ds:
        with pytest.raises(ValueError, match="requires both the 'feature' and 'psm'"):
            ds.link_feature_psm()
        with pytest.raises(ValueError, match="requires both the 'feature' and 'psm'"):
            ds.psms_without_feature()
