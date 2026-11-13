"""Tests for the derived mandatory identity ids (feature_id / psm_id / pg_id)."""

import logging
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.core.data import FeatureSchema, PgSchema, PsmSchema
from qpx.core.data.identity import canonical, derive_composite_id, derive_id
from qpx.writers import FeatureWriter, PgWriter, PsmWriter
from tests.conftest import make_feature_record, make_pg_record, make_psm_record


def test_derive_id_is_deterministic():
    """Same input -> same id, on repeated calls."""
    values = ["PEPTIDEK", 2, "run_01", 120.5]
    assert derive_id(values) == derive_id(list(values))


def test_derive_id_distinguishes_inputs():
    """Different composites -> different ids (with overwhelming probability)."""
    assert derive_id(["PEPTIDEK", 2, "run_01", 120.5]) != derive_id(["PEPTIDEK", 3, "run_01", 120.5])
    assert derive_id(["PEPTIDEK", 2, "run_01", 120.5]) != derive_id(["PEPTIDEQ", 2, "run_01", 120.5])


def test_derive_id_handles_none():
    """None is encoded distinctly and does not collide with the string 'None'."""
    assert derive_id([None]) == derive_id([None])
    assert derive_id([None]) != derive_id(["None"])
    assert derive_id(["A", None, "B"]) != derive_id(["A", "B"])
    # int64 range
    for v in (derive_id([None]), derive_id(["A", None, "B"])):
        assert -(2**63) <= v < 2**63


def test_canonical_is_injective_json():
    """The canonical encoding is injective JSON: delimiter-bearing values cannot
    alias, and None becomes JSON null."""
    assert canonical(["A", "B"]) == b'["A","B"]'
    assert canonical([None]) == b"[null]"
    # list elements containing commas/brackets do not collide across groupings
    assert canonical([["a,b", "c"]]) != canonical([["a", "b,c"]])


def test_schema_identity_metadata():
    """Each identity view exposes a single-column PK and its composite."""
    assert tuple(FeatureSchema._primary_key) == ("feature_id",)
    assert tuple(FeatureSchema.identity_composite) == ("peptidoform", "charge", "run_file_name", "rt")
    assert tuple(PsmSchema._primary_key) == ("psm_id",)
    assert tuple(PsmSchema.identity_composite) == ("run_file_name", "scan", "peptidoform", "charge")
    assert tuple(PgSchema._primary_key) == ("pg_id",)
    assert tuple(PgSchema.identity_composite) == ("pg_accessions", "grouped_runs", "label")


def test_writer_rejects_unknown_identity_composite_field(tmp_path):
    """A converter typo must fail instead of hashing an implicit null column."""
    with pytest.raises(ValueError, match="outside the feature schema.*missing_field"):
        FeatureWriter(tmp_path / "invalid.feature.parquet", identity_composite=("run_file_name", "missing_field"))


def test_feature_roundtrip_derives_id(tmp_path):
    """FeatureWriter stamps a non-null, unique feature_id matching derive_id."""
    path = tmp_path / "t.feature.parquet"
    records = [
        make_feature_record(sequence="PEPTIDEK", peptidoform="PEPTIDEK", charge=2, run_file_name="run_01"),
        make_feature_record(sequence="ANOTHERK", peptidoform="ANOTHERK", charge=3, run_file_name="run_01", scan=[2]),
    ]
    with FeatureWriter(path) as w:
        w.write_batch([dict(r) for r in records])

    table = pq.read_table(path)
    ids = table.column("feature_id").to_pylist()
    assert None not in ids
    assert len(set(ids)) == len(ids)
    for rec, got in zip(records, ids):
        expected = derive_id([rec["peptidoform"], rec["charge"], rec["run_file_name"], rec["rt"]])
        assert got == expected

    # Footer self-describes the key and composite.
    meta = {k.decode(): v.decode() for k, v in table.schema.metadata.items()}
    assert meta["primary_key"] == "feature_id"
    assert meta["identity_composite"] == "peptidoform,charge,run_file_name,rt"


def test_psm_roundtrip_derives_id(tmp_path):
    """PsmWriter stamps a non-null, unique psm_id matching derive_id."""
    path = tmp_path / "t.psm.parquet"
    records = [
        make_psm_record(sequence="PEPTIDEK", peptidoform="PEPTIDEK", charge=2, run_file_name="run_01", scan=[10]),
        make_psm_record(sequence="PEPTIDEK", peptidoform="PEPTIDEK", charge=2, run_file_name="run_01", scan=[11]),
    ]
    with PsmWriter(path) as w:
        w.write_batch([dict(r) for r in records])

    table = pq.read_table(path)
    ids = table.column("psm_id").to_pylist()
    assert None not in ids
    assert len(set(ids)) == len(ids)
    for rec, got in zip(records, ids):
        expected = derive_id([rec["run_file_name"], rec["scan"], rec["peptidoform"], rec["charge"]])
        assert got == expected


def test_pg_roundtrip_derives_id(tmp_path):
    """PgWriter stamps a non-null, unique pg_id matching derive_id of the flat composite."""
    path = tmp_path / "t.pg.parquet"
    records = [
        make_pg_record(anchor_protein="P1", run_file_name="run_01", intensities=[{"label": "L1", "intensity": 1.0}]),
        make_pg_record(anchor_protein="P2", run_file_name="run_01", intensities=[{"label": "L1", "intensity": 2.0}]),
    ]
    with PgWriter(path) as w:
        w.write_batch([dict(r) for r in records])

    table = pq.read_table(path)
    ids = table.column("pg_id").to_pylist()
    pg_accessions = table.column("pg_accessions").to_pylist()
    grouped = table.column("grouped_runs").to_pylist()
    labels = table.column("label").to_pylist()
    assert None not in ids
    assert len(set(ids)) == len(ids)
    for i, got in enumerate(ids):
        assert got == derive_id([pg_accessions[i], grouped[i], labels[i]], unordered_list_indices=(0, 1))


def test_identified_feature_id_is_derived_by_default(tmp_path):
    """An identified Feature must use the reproducible composite-derived id."""
    path = tmp_path / "t.feature.parquet"
    rec = make_feature_record(sequence="PEPTIDEK", peptidoform="PEPTIDEK", charge=2, run_file_name="run_01")
    rec["feature_id"] = 123456789
    with FeatureWriter(path) as w:
        w.write_batch([dict(rec)])
    table = pq.read_table(path)
    derived = derive_id([rec["peptidoform"], rec["charge"], rec["run_file_name"], rec["rt"]])
    assert table.column("feature_id").to_pylist() == [derived]
    assert {"cv_name": "provided_feature_id", "cv_value": "123456789"} in table.column("cv_params").to_pylist()[0]


def test_unidentified_feature_namespaces_provided_id(tmp_path):
    """An unidentified Feature namespaces its required producer identity by run."""
    path = tmp_path / "unidentified.feature.parquet"
    records = [make_feature_record(sequence="", peptidoform="", run_file_name=run) for run in ("run_01", "run_02")]
    for record in records:
        record["feature_id"] = 123456789
    with FeatureWriter(path) as writer:
        writer.write_batch(records)
    table = pq.read_table(path)
    assert table.column("feature_id").to_pylist() == [derive_id([record["run_file_name"], 123456789]) for record in records]
    for params in table.column("cv_params").to_pylist():
        assert {"cv_name": "provided_feature_id", "cv_value": "123456789"} in params


def test_unidentified_feature_uses_unique_composite_fallback(tmp_path, caplog):
    """An unidentified Feature (null peptidoform) with no producer id derives a
    feature_id only when the resulting full-file primary key is unique."""
    path = tmp_path / "unidentified.feature.parquet"
    rec = make_feature_record(sequence="", peptidoform="")
    with caplog.at_level(logging.WARNING, logger="qpx.writers.feature"):
        with FeatureWriter(path) as writer:
            writer.write_batch([rec])
    ids = pq.read_table(path).column("feature_id").to_pylist()
    assert ids and ids[0] is not None
    assert "1 unidentified Feature row(s)" in caplog.text
    assert "full-file primary-key uniqueness was verified" in caplog.text


@pytest.mark.parametrize("write_path", ["batch", "table"])
def test_unidentified_feature_warns_on_non_unique_composite_fallback(tmp_path, write_path, caplog):
    """Distinct unidentified Features sharing a fallback PK are tolerated with a warning."""
    path = tmp_path / "duplicate-unidentified.feature.parquet"
    first = make_feature_record(sequence="", peptidoform="")
    second = make_feature_record(
        sequence="",
        peptidoform="",
        intensities=[{"label": "TMT126", "intensity": 2000.0}],
    )

    with caplog.at_level(logging.WARNING):
        with FeatureWriter(path, batch_size=1) as writer:
            if write_path == "batch":
                writer.write_batch([first, second])
            else:
                table = writer.align_table_to_schema(pa.Table.from_pylist([first, second]))
                writer.write_table(table)

    assert path.exists()
    assert "duplicate row" in caplog.text.lower()


def test_identified_duplicate_warns_without_fallback_note(tmp_path, caplog):
    """A duplicate among identified Features warns with the plain PK message even when a
    (unique) unidentified composite-fallback row is present in the same file."""
    path = tmp_path / "mixed-duplicate.feature.parquet"
    fallback = make_feature_record(sequence="", peptidoform="")  # unique unidentified fallback
    dup_a = make_feature_record(run_file_name="run_09")
    dup_b = make_feature_record(run_file_name="run_09")  # identical composite -> duplicate PK

    with caplog.at_level(logging.WARNING):
        with FeatureWriter(path, batch_size=10) as writer:
            writer.write_batch([fallback, dup_a, dup_b])

    assert path.exists()
    assert "duplicate row" in caplog.text.lower()
    assert "Unidentified Feature composite fallback" not in caplog.text


def test_provided_psm_id_kept_by_default(tmp_path):
    """A producer-supplied PSM id remains supported by the generic writer contract."""
    path = tmp_path / "t.psm.parquet"
    rec = make_psm_record()
    rec["psm_id"] = 123456789
    with PsmWriter(path) as writer:
        writer.write_batch([rec])
    assert pq.read_table(path).column("psm_id").to_pylist() == [123456789]


def test_writer_warns_on_duplicate_identity_across_batches(tmp_path, caplog):
    """Whole-file PK validation warns on collisions split across writer batches."""
    path = tmp_path / "duplicate.feature.parquet"
    first = make_feature_record()
    second = make_feature_record(intensities=[{"label": "TMT126", "intensity": 2000.0}])

    with caplog.at_level(logging.WARNING):
        with FeatureWriter(path, batch_size=1) as writer:
            writer.write_batch([first, second])

    assert path.exists()
    assert "Primary key (feature_id) has 1 duplicate row" in caplog.text


def test_override_provided_id_stashes_to_cv_params(tmp_path):
    """With override_provided_ids, the id is re-derived and the original stashed as a cv_param."""
    path = tmp_path / "t.feature.parquet"
    rec = make_feature_record(sequence="PEPTIDEK", peptidoform="PEPTIDEK", charge=2, run_file_name="run_01")
    rec["feature_id"] = 123456789
    w = FeatureWriter(path, override_provided_ids=True)
    w.write_batch([dict(rec)])
    w.close()

    table = pq.read_table(path)
    derived = derive_id([rec["peptidoform"], rec["charge"], rec["run_file_name"], rec["rt"]])
    assert table.column("feature_id").to_pylist() == [derived]
    cv = table.column("cv_params").to_pylist()[0]
    assert {"cv_name": "provided_feature_id", "cv_value": "123456789"} in cv


def test_override_requires_cv_params_to_preserve_provided_id(tmp_path):
    """Override mode rejects a table that cannot retain its producer id."""
    path = tmp_path / "missing-cv.feature.parquet"
    record = make_feature_record()
    record["feature_id"] = 123456789
    writer = FeatureWriter(path, override_provided_ids=True)
    table = writer.align_table_to_schema(pa.Table.from_pylist([record])).drop(["cv_params"])

    with pytest.raises(ValueError, match="requires a 'cv_params' column to preserve provided feature_id values"):
        with writer:
            writer.write_table(table)

    assert not path.exists()


def test_derive_id_preserves_ordered_list_components():
    """A multi-component scan is a tuple, not a set."""
    assert derive_id(["P1", [1, 2, 3]]) != derive_id(["P1", [3, 1, 2]])


def test_derive_id_can_canonicalize_set_valued_lists():
    """Explicitly unordered grouped_runs values are compared set-wise."""
    assert derive_id([["r1", "r2"]], unordered_list_indices=(0,)) == derive_id([["r2", "r1"]], unordered_list_indices=(0,))
    assert derive_id([["r1", "r2", "r1"]], unordered_list_indices=(0,)) == derive_id([["r2", "r1"]], unordered_list_indices=(0,))
    assert derive_id([["r1", "r2"]], unordered_list_indices=(0,)) != derive_id([["r1", "r3"]], unordered_list_indices=(0,))


def test_writer_preserves_scan_order_but_canonicalizes_grouped_runs(tmp_path):
    """Writer applies unordered semantics only to the grouped_runs field."""
    psm_path = tmp_path / "ordered.psm.parquet"
    psm_a = make_psm_record(peptidoform="PEPTIDEK", run_file_name="run_01", scan=[10, 1, 345])
    psm_b = make_psm_record(peptidoform="PEPTIDEK", run_file_name="run_01", scan=[345, 1, 10])
    with PsmWriter(psm_path) as writer:
        writer.write_batch([psm_a, psm_b])
    assert len(set(pq.read_table(psm_path).column("psm_id").to_pylist())) == 2

    pg_ids = []
    for name, grouped_runs in (("forward", ["run_01", "run_02"]), ("reverse", ["run_02", "run_01"])):
        pg_path = tmp_path / f"{name}.pg.parquet"
        pg = make_pg_record(anchor_protein="P1", run_file_name="run_01")
        pg["grouped_runs"] = grouped_runs
        with PgWriter(pg_path) as writer:
            writer.write_batch([pg])
        pg_ids.extend(pq.read_table(pg_path).column("pg_id").to_pylist())
    assert pg_ids[0] == pg_ids[1]


def test_write_table_keeps_a_provided_int32_id(tmp_path):
    """A provided id column of a narrower integer type is kept and written as int64."""
    path = tmp_path / "t.psm.parquet"
    rec = make_psm_record()
    rec["psm_id"] = 123456789
    with PsmWriter(path) as writer:
        table = writer.align_table_to_schema(pa.Table.from_pylist([rec]))
        index = table.schema.get_field_index("psm_id")
        writer.write_table(table.set_column(index, "psm_id", table.column(index).cast(pa.int32())))
    out = pq.read_table(path).column("psm_id")
    assert out.type == pa.int64()
    assert out.to_pylist() == [123456789]


def test_override_applies_on_write_table_path(tmp_path):
    """override_provided_ids must also take effect on the write_table/dataframe path."""
    import pyarrow.parquet as pq

    path = tmp_path / "t.feature.parquet"
    rec = make_feature_record(sequence="PEPTIDEK", peptidoform="PEPTIDEK", charge=2, run_file_name="run_01")
    rec["feature_id"] = 424242
    # Build an Arrow table (write_table path) carrying the provided id.
    w = FeatureWriter(path, override_provided_ids=True)
    table = w.align_table_to_schema(pa.Table.from_pylist([{k: v for k, v in rec.items()}]))
    w.write_table(table)
    w.close()

    out = pq.read_table(path)
    derived = derive_id([rec["peptidoform"], rec["charge"], rec["run_file_name"], rec["rt"]])
    assert out.column("feature_id").to_pylist() == [derived]
    assert {"cv_name": "provided_feature_id", "cv_value": "424242"} in out.column("cv_params").to_pylist()[0]


def test_ids_agree_across_write_paths(tmp_path):
    """The same row must derive the same feature_id via write_batch, write_table,
    and write_dataframe, even for a non-float32-exact rt (persisted-value hashing)."""
    import pandas as pd

    rec = make_feature_record(sequence="PEPTIDEK", peptidoform="PEPTIDEK", charge=2, run_file_name="run_01")
    rec["rt"] = 450.26  # not exactly representable in float32

    def written_id(name, how):
        p = tmp_path / name
        w = FeatureWriter(p)
        if how == "batch":
            w.write_batch([dict(rec)])
        elif how == "table":
            w.write_table(w.align_table_to_schema(pa.Table.from_pylist([dict(rec)])))
        else:
            w.write_dataframe(pd.DataFrame([dict(rec)]))
        w.close()
        return pq.read_table(p).column("feature_id").to_pylist()[0]

    ids = {
        written_id("b.feature.parquet", "batch"),
        written_id("t.feature.parquet", "table"),
        written_id("d.feature.parquet", "df"),
    }
    assert len(ids) == 1


def test_write_table_warns_on_duplicate_pk(tmp_path, caplog):
    """The table-writing path tolerates a duplicate primary key with a warning."""
    path = tmp_path / "dup_table.feature.parquet"
    r1 = make_feature_record(peptidoform="PEPTIDEK", charge=2, run_file_name="run_01")
    r2 = make_feature_record(
        peptidoform="PEPTIDEK", charge=2, run_file_name="run_01", intensities=[{"label": "TMT126", "intensity": 9.0}]
    )
    with caplog.at_level(logging.WARNING):
        with FeatureWriter(path) as w:
            table = w.align_table_to_schema(pa.Table.from_pylist([dict(r1), dict(r2)]))
            w.write_table(table)
    assert path.exists()
    assert "duplicate row" in caplog.text.lower()


# The vectors of docs/spec/identity.md: a reimplementation checks itself against them.
_VECTORS = [
    (("run_file_name", "scan", "peptidoform", "charge"), ["run_1", [1234], "PEPTIDEK", 2], -5127882623363935776),
    (
        ("peptidoform", "charge", "run_file_name", "rt"),
        ["[UNIMOD:1]-AC[UNIMOD:4]DEK", 3, "a,b", 1234.56005859375],
        6417640111734626350,
    ),
    (("pg_accessions", "grouped_runs", "label"), [["Q9Y6K9", "P04114", "P04114"], ["f2", "f1"], "TMT126"], -5459222103841755848),
    (("pg_accessions", "grouped_runs", "label"), [["P04114"], ["r"], None], -4507348376872481825),
]


@pytest.mark.parametrize(("composite", "values", "expected"), _VECTORS)
def test_documented_vectors_reproduce(composite, values, expected):
    row = dict(zip(composite, values))
    assert derive_composite_id(row, composite) == expected
    spec = (Path(__file__).parents[2] / "docs" / "spec" / "identity.md").read_text(encoding="utf-8")
    assert f"`{expected}`" in spec


def test_a_float32_composite_value_hashes_as_persisted(tmp_path):
    """A record's raw double and the file's float32 read back must give one id."""
    path = tmp_path / "x.feature.parquet"
    record = {**make_feature_record(), "rt": 1234.56}
    with FeatureWriter(path, identity_composite=("peptidoform", "charge", "run_file_name", "rt")) as writer:
        writer.write_batch([record])
    row = pq.read_table(path).to_pylist()[0]
    assert row["rt"] == 1234.56005859375
    assert row["feature_id"] == derive_composite_id(row, ("peptidoform", "charge", "run_file_name", "rt"))


def _footer(path):
    return pq.read_metadata(path).metadata


def test_footer_states_where_ids_came_from(tmp_path):
    derived = tmp_path / "derived.psm.parquet"
    with PsmWriter(derived) as writer:
        writer.write_batch([make_psm_record()])
    assert _footer(derived)[b"identity_source"] == b"derived"
    assert _footer(derived)[b"identity_hash_version"] == b"1"
    assert _footer(derived)[b"identity_composite"] == b"run_file_name,scan,peptidoform,charge"

    provided = tmp_path / "provided.psm.parquet"
    with PsmWriter(provided, override_provided_ids=False) as writer:
        writer.write_batch([{**make_psm_record(), "psm_id": 7}])
    assert _footer(provided)[b"identity_source"] == b"provided"

    mixed = tmp_path / "mixed.psm.parquet"
    with PsmWriter(mixed, override_provided_ids=False) as writer:
        writer.write_batch([{**make_psm_record(scan=[1]), "psm_id": 7}, make_psm_record(scan=[2])])
    assert _footer(mixed)[b"identity_source"] == b"mixed"


def test_partitioned_files_state_where_ids_came_from(tmp_path):
    source = tmp_path / "source.psm.parquet"
    with PsmWriter(source) as writer:
        writer.write_batch([make_psm_record()])
    table = pq.read_table(source).drop_columns(["psm_id"])
    PsmWriter.write_partitioned(table, tmp_path / "psm", partition_cols=["run_file_name"])
    part = next((tmp_path / "psm").rglob("*.parquet"))
    assert _footer(part)[b"identity_source"] == b"derived"


def _dataset(tmp_path, **kwargs):
    import qpx

    with PsmWriter(tmp_path / "exp.psm.parquet", **kwargs) as writer:
        writer.write_batch([make_psm_record()])
    return qpx.open_dataset(str(tmp_path))


def test_check_identity_passes_on_derived_ids_and_catches_stale_ones(tmp_path):
    with _dataset(tmp_path) as dataset:
        result = dataset.validate(["psm"], strict=True, check_identity=True)["psm"]
    assert not [issue for issue in result.issues if issue.check == "stale_identity"]

    path = tmp_path / "exp.psm.parquet"
    table = pq.read_table(path)
    changed = table.set_column(table.schema.get_field_index("charge"), "charge", pa.array([5], table.schema.field("charge").type))
    pq.write_table(changed, path)
    import qpx

    with qpx.open_dataset(str(tmp_path)) as dataset:
        result = dataset.validate(["psm"], strict=True, check_identity=True)["psm"]
    assert any(issue.check == "stale_identity" for issue in result.errors)


def test_check_identity_skips_provided_ids(tmp_path):
    with PsmWriter(tmp_path / "exp.psm.parquet", override_provided_ids=False) as writer:
        writer.write_batch([{**make_psm_record(), "psm_id": 7}])
    import qpx

    with qpx.open_dataset(str(tmp_path)) as dataset:
        result = dataset.validate(["psm"], strict=True, check_identity=True)["psm"]
    issues = [issue for issue in result.issues if issue.check == "identity_not_checked"]
    assert len(issues) == 1 and issues[0].severity == "info"
    assert not any(issue.check == "identity_not_checked" for issue in result.warnings)
    assert result.is_valid
    assert "1 info" in result.summary
    assert not any(issue.check == "stale_identity" for issue in result.issues)


def test_files_with_different_composites_are_flagged(tmp_path):
    import qpx

    with PsmWriter(tmp_path / "exp-a.psm.parquet") as writer:
        writer.write_batch([make_psm_record(scan=[1])])
    with PsmWriter(
        tmp_path / "exp-b.psm.parquet", identity_composite=("run_file_name", "peptidoform", "charge", "scan")
    ) as writer:
        writer.write_batch([make_psm_record(scan=[2])])
    with qpx.open_dataset(str(tmp_path)) as dataset:
        result = dataset.validate(["psm"])["psm"]
    assert any(issue.check == "mixed_identity" for issue in result.warnings)
