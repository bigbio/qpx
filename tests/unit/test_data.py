"""Tests for QPX data layer: BaseStructure, Feature, PSM, PG and domain methods."""

import pandas as pd
import pyarrow as pa

from qpx.core.data.dataset import DatasetMeta
from qpx.core.data.feature import Feature
from qpx.core.data.ontology import Ontology
from qpx.core.data.pg import PG
from qpx.core.data.provenance import Provenance
from qpx.core.data.psm import PSM
from qpx.core.data.run import Run
from qpx.core.data.sample import Sample

# ---------------------------------------------------------------------------
# Parameterized from_file across all structures
# ---------------------------------------------------------------------------


def test_structures_from_file(
    feature_parquet,
    psm_parquet,
    pg_parquet,
    sample_parquet,
    run_parquet,
    dataset_parquet,
    ontology_parquet,
    provenance_parquet,
):
    """All structures open from Parquet, return correct type and count."""
    cases = [
        (feature_parquet, Feature, 3),
        (psm_parquet, PSM, 3),
        (pg_parquet, PG, 3),
        (sample_parquet, Sample, 2),
        (run_parquet, Run, 2),
        (dataset_parquet, DatasetMeta, 1),
        (ontology_parquet, Ontology, 2),
        (provenance_parquet, Provenance, 2),
    ]
    for path, cls, expected_count in cases:
        obj = cls.from_file(path)
        assert isinstance(obj, cls)
        assert obj.count() == expected_count


def test_feature_filter_and_query_ops(feature_parquet):
    """filter, select, limit, count, len, and chaining all work correctly."""
    feat = Feature.from_file(feature_parquet)

    # count and len
    assert feat.count() == 3
    assert len(feat) == 3

    # filter
    filtered = feat.filter("charge = 2")
    assert isinstance(filtered, Feature)
    assert filtered.count() == 2

    # by_protein / by_run domain filters
    assert feat.by_protein("P12345").count() == 2
    assert feat.by_run("run_01").count() == 2

    # filter chaining
    assert feat.filter("charge = 2").filter("is_decoy = false").count() == 1

    # select
    selected = feat.select("sequence", "charge")
    df = selected.to_df()
    assert set(df.columns) == {"sequence", "charge"}
    assert len(df) == 3

    # limit
    assert feat.limit(2).count() == 2

    # count after filter
    assert feat.filter("is_decoy = true").count() == 1


def test_feature_materialization(feature_parquet):
    """to_df returns DataFrame, to_arrow returns Arrow Table."""
    feat = Feature.from_file(feature_parquet)
    df = feat.to_df()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 3

    table = feat.to_arrow()
    assert isinstance(table, pa.Table)
    assert table.num_rows == 3


def test_feature_domain_methods(feature_parquet):
    """unique_proteins, peptide_intensities, file_metadata, and repr."""
    feat = Feature.from_file(feature_parquet)

    # unique_proteins
    assert set(feat.unique_proteins()) == {"P12345", "P67890"}

    # peptide_intensities
    result = feat.peptide_intensities()
    df = result.to_df()
    assert "label" in df.columns
    assert "intensity" in df.columns
    assert len(df) == 3

    # file_metadata
    meta = feat.file_metadata
    assert isinstance(meta, dict)
    assert meta["file_type"] == "feature_file"

    # repr
    r = repr(feat)
    assert "Feature" in r
    assert "rows=3" in r


def test_immutability(feature_parquet):
    """filter, select, limit return new instances without mutating the original."""
    feat = Feature.from_file(feature_parquet)
    original_count = feat.count()

    filtered = feat.filter("charge > 10")
    assert filtered.count() == 0
    assert feat.count() == original_count

    selected = feat.select("sequence")
    assert selected is not feat

    limited = feat.limit(1)
    assert limited.count() == 1
    assert feat.count() == original_count


def test_iter_batches(feature_parquet):
    """iter_batches partitions correctly with different batch sizes."""
    feat = Feature.from_file(feature_parquet)

    # Large batch: one batch with all distinct runs
    batches = list(feat.iter_batches(partition_by="run_file_name", batch_size=10))
    assert len(batches) == 1
    keys, df = batches[0]
    assert set(keys) == {"run_01", "run_02"}
    assert len(df) == 3

    # Small batch: one batch per distinct value
    batches = list(feat.iter_batches(partition_by="run_file_name", batch_size=1))
    assert len(batches) == 2
    for keys, df in batches:
        assert len(keys) == 1
        assert isinstance(df, pd.DataFrame)
        run_values = df["run_file_name"].unique()
        assert set(run_values) == set(keys)


def test_order_by(feature_parquet):
    """order_by sorts ascending and descending correctly."""
    feat = Feature.from_file(feature_parquet)

    charges_asc = feat.order_by("charge").to_df()["charge"].tolist()
    assert charges_asc == sorted(charges_asc)

    charges_desc = feat.order_by("charge", desc=True).to_df()["charge"].tolist()
    assert charges_desc == sorted(charges_desc, reverse=True)


def test_feature_join_run(feature_parquet, run_parquet):
    """Feature.join(Run) joins on run_file_name."""
    feat = Feature.from_file(feature_parquet)
    feat._engine.register_parquet("run", run_parquet)
    run = Run(engine=feat._engine, table_name="run", file_path=run_parquet)

    joined = feat.join(run, on="run_file_name")
    df = joined.to_df()
    assert len(df) == 3
    assert "run_accession" in df.columns


# ---------------------------------------------------------------------------
# PG grouped_runs field (refactor: run_file_name scalar -> grouped_runs list<string>)
# ---------------------------------------------------------------------------


def test_pg_schema_grouped_runs_is_list_string():
    """The pg schema keys on `grouped_runs` typed as list<string>."""
    from qpx.core.data.pg import PgSchema

    arrow_schema = PgSchema.get_arrow_schema()
    assert "run_file_name" not in arrow_schema.names
    field = arrow_schema.field("grouped_runs")
    assert pa.types.is_list(field.type)
    assert pa.types.is_string(field.type.value_type)
    assert tuple(PgSchema._primary_key) == ("pg_id",)
    assert tuple(PgSchema.identity_composite) == ("pg_accessions", "grouped_runs", "label")


def test_pg_grouped_runs_roundtrip(pg_parquet):
    """A round-tripped pg record exposes `grouped_runs` as a single-element list."""
    pg = PG.from_file(pg_parquet)
    df = pg.to_df()
    assert "grouped_runs" in df.columns
    assert "run_file_name" not in df.columns
    row = df[df["anchor_protein"] == "P12345"].iloc[0]
    assert list(row["grouped_runs"]) == ["run_01"]


def test_pg_by_run_filters_via_grouped_runs(pg_parquet):
    """PG.by_run matches on list membership of the grouped_runs list."""
    pg = PG.from_file(pg_parquet)
    assert pg.by_run("run_01").count() == 2
    assert pg.by_run("run_02").count() == 1
    assert pg.by_run("missing").count() == 0


# ---------------------------------------------------------------------------
# grouped_runs referential invariant (pg.grouped_runs subset of run.run_file_name)
# ---------------------------------------------------------------------------


def _write_pg_run_dataset(
    ds_dir,
    pg_grouped_runs,
    run_sample_accessions=("SAMPLE_01", "SAMPLE_01"),
):
    """Write a minimal pg + run dataset; pg rows use the given grouped_runs lists."""
    from qpx.writers import DatasetWriter, PgWriter, RunWriter
    from tests.conftest import (
        make_dataset_record,
        make_pg_record,
        make_run_record,
    )

    ds_dir.mkdir(parents=True, exist_ok=True)
    pg_records = []
    for idx, grouped in enumerate(pg_grouped_runs):
        rec = make_pg_record(anchor_protein=f"P{idx:05d}")
        rec["grouped_runs"] = grouped
        pg_records.append(rec)
    with PgWriter(ds_dir / "exp.pg.parquet") as w:
        w.write_batch(pg_records)
    run_records = [
        make_run_record(run_accession="assay_01", run_file_name="run_01"),
        make_run_record(run_accession="assay_02", run_file_name="run_02"),
    ]
    for record, sample_accession in zip(
        run_records,
        run_sample_accessions,
        strict=True,
    ):
        record["samples"][0]["sample_accession"] = sample_accession
    with RunWriter(ds_dir / "exp.run.parquet") as w:
        w.write_batch(run_records)
    with DatasetWriter(ds_dir / "exp.dataset.parquet") as w:
        w.write_batch([make_dataset_record()])
    return ds_dir


def test_grouped_runs_dangling_reference_is_error(tmp_path):
    """A pg.grouped_runs value with no matching run.run_file_name is an ERROR."""
    from qpx.dataset import Dataset

    ds_dir = _write_pg_run_dataset(
        tmp_path / "bad",
        pg_grouped_runs=[["run_01"], ["run_99"]],  # run_99 does not exist
    )
    with Dataset(ds_dir) as ds:
        results = ds.validate(structures=["pg"])

    pg_result = results["pg"]
    dangling = [i for i in pg_result.issues if i.check == "dangling_grouped_run"]
    assert len(dangling) == 1
    assert dangling[0].severity == "error"
    assert "run_99" in dangling[0].message
    assert not pg_result.is_valid  # the error must fail validation


def test_grouped_runs_all_valid_passes_invariant(tmp_path):
    """When every grouped_runs value is a real run_file_name, no dangling error."""
    from qpx.dataset import Dataset

    ds_dir = _write_pg_run_dataset(
        tmp_path / "good",
        pg_grouped_runs=[["run_01"], ["run_02"], ["run_01", "run_02"]],
    )
    with Dataset(ds_dir) as ds:
        results = ds.validate(structures=["pg"])

    dangling = [i for i in results["pg"].issues if i.check == "dangling_grouped_run"]
    assert dangling == []


def test_pg_auto_join_run_expands_grouped_runs(tmp_path):
    """The documented PG-to-Run join expands every grouped run."""
    from qpx.dataset import Dataset

    ds_dir = _write_pg_run_dataset(
        tmp_path / "join",
        pg_grouped_runs=[["run_01", "run_02"]],
    )
    with Dataset(ds_dir) as ds:
        joined = ds.pg.join(ds.run).to_df()

    assert len(joined) == 2
    assert set(joined["run_file_name"]) == {"run_01", "run_02"}


def test_grouped_runs_label_resolving_to_multiple_samples_is_error(tmp_path):
    """All fractions in one PG unit must resolve a label to the same sample."""
    from qpx.dataset import Dataset

    ds_dir = _write_pg_run_dataset(
        tmp_path / "ambiguous",
        pg_grouped_runs=[["run_01", "run_02"]],
        run_sample_accessions=("SAMPLE_01", "SAMPLE_02"),
    )
    with Dataset(ds_dir) as ds:
        results = ds.validate(structures=["pg"])

    ambiguous = [issue for issue in results["pg"].issues if issue.check == "ambiguous_grouped_run_sample"]
    assert len(ambiguous) == 1
    assert ambiguous[0].severity == "error"
    assert "2 samples" in ambiguous[0].message
    assert not results["pg"].is_valid


def test_strict_validation_promotes_skipped_cross_structure_checks_to_errors(tmp_path):
    """Schema drift in a cross-view check must fail strict validation."""
    from qpx.dataset import Dataset

    ds_dir = _write_pg_run_dataset(
        tmp_path / "drift",
        pg_grouped_runs=[["run_01"]],
    )
    with Dataset(ds_dir) as ds:
        ds._engine.execute("CREATE OR REPLACE VIEW run AS SELECT 1 AS wrong_column")
        lenient = ds.validate(structures=["pg"])["pg"]
        strict = ds.validate(structures=["pg"], strict=True)["pg"]

    lenient_skipped = [issue for issue in lenient.issues if issue.check.endswith("check_skipped")]
    strict_skipped = [issue for issue in strict.issues if issue.check.endswith("check_skipped")]
    assert len(lenient_skipped) == 2
    assert all(issue.severity == "warning" for issue in lenient_skipped)
    assert len(strict_skipped) == 2
    assert all(issue.severity == "error" for issue in strict_skipped)
    assert not strict.is_valid


# ---------------------------------------------------------------------------
# feature<->psm cross-reference referential invariant
# ---------------------------------------------------------------------------


def _write_feature_psm_dataset(ds_dir, *, feature_psm_ids, psm_feature_id):
    """Write a minimal feature + psm dataset with explicit cross-refs.

    ``feature_psm_ids`` is assigned to the single feature's ``psm_ids`` and
    ``psm_feature_id`` to the single psm's ``feature_id``. Returns the dir plus
    the writer-derived ``(feature_id, psm_id)`` so callers can build resolving refs.
    """
    from qpx.core.data.identity import derive_id
    from qpx.writers import FeatureWriter, PsmWriter
    from tests.conftest import make_feature_record, make_psm_record

    ds_dir.mkdir(parents=True, exist_ok=True)
    feat = make_feature_record()
    psm = make_psm_record()
    # The writers derive the identity ids from the schema-default composites; mirror
    # them so a caller can construct a *resolving* cross-ref.
    feature_id = derive_id([feat["peptidoform"], feat["charge"], feat["run_file_name"], feat["rt"]])
    psm_id = derive_id([psm["peptidoform"], psm["charge"], psm["run_file_name"], psm["scan"]])
    feat["psm_ids"] = feature_psm_ids
    psm["feature_id"] = psm_feature_id
    with FeatureWriter(ds_dir / "exp.feature.parquet") as w:
        w.write_batch([feat])
    with PsmWriter(ds_dir / "exp.psm.parquet") as w:
        w.write_batch([psm])
    return ds_dir, feature_id, psm_id


def _drop_parquet_columns(path, columns):
    """Model a valid external producer that omits optional schema columns."""
    import pyarrow.parquet as pq

    table = pq.read_table(path).drop_columns(columns)
    pq.write_table(table, path)


def test_feature_psm_cross_refs_resolve_cleanly(tmp_path):
    """When every cross-ref resolves, no dangling-reference issue is raised."""
    from qpx.dataset import Dataset

    # Two-pass: first derive the ids, then rewrite with resolving refs.
    _, feature_id, psm_id = _write_feature_psm_dataset(tmp_path / "probe", feature_psm_ids=None, psm_feature_id=None)
    ds_dir, _, _ = _write_feature_psm_dataset(tmp_path / "good", feature_psm_ids=[psm_id], psm_feature_id=feature_id)

    with Dataset(ds_dir) as ds:
        results = ds.validate(structures=["feature", "psm"])

    dangling = [i for r in results.values() for i in r.issues if i.check in ("dangling_feature_id", "dangling_psm_id")]
    assert dangling == []


def test_missing_optional_feature_psm_cross_refs_are_valid(tmp_path):
    """External files may omit either optional cross-reference column."""
    from qpx.dataset import Dataset

    ds_dir, _, _ = _write_feature_psm_dataset(
        tmp_path / "without-cross-refs",
        feature_psm_ids=None,
        psm_feature_id=None,
    )
    _drop_parquet_columns(ds_dir / "exp.feature.parquet", ["psm_ids"])
    _drop_parquet_columns(ds_dir / "exp.psm.parquet", ["feature_id"])

    with Dataset(ds_dir) as dataset:
        results = dataset.validate(structures=["feature", "psm"], strict=True)

    assert results["feature"].is_valid
    assert results["psm"].is_valid
    assert not [issue for result in results.values() for issue in result.issues if issue.check == "referential_check_skipped"]


def test_present_psm_cross_ref_is_checked_when_inverse_column_is_absent(tmp_path):
    """A missing feature.psm_ids column must not disable psm.feature_id checks."""
    from qpx.dataset import Dataset

    ds_dir, _, _ = _write_feature_psm_dataset(
        tmp_path / "one-way-cross-ref",
        feature_psm_ids=None,
        psm_feature_id=999999,
    )
    _drop_parquet_columns(ds_dir / "exp.feature.parquet", ["psm_ids"])

    with Dataset(ds_dir) as dataset:
        results = dataset.validate(structures=["feature", "psm"], strict=True)

    dangling = [issue for issue in results["psm"].issues if issue.check == "dangling_feature_id"]
    assert len(dangling) == 1
    assert not [issue for issue in results["psm"].issues if issue.check == "referential_check_skipped"]


def test_feature_psm_cross_refs_dangling_are_warnings(tmp_path):
    """A psm.feature_id / feature.psm_ids element with no sibling id -> warning."""
    from qpx.dataset import Dataset

    ds_dir, _, _ = _write_feature_psm_dataset(
        tmp_path / "bad",
        feature_psm_ids=[888888],  # no such psm_id
        psm_feature_id=999999,  # no such feature_id
    )
    with Dataset(ds_dir) as ds:
        results = ds.validate(structures=["feature", "psm"])

    dangling_feature = [i for i in results["psm"].issues if i.check == "dangling_feature_id"]
    dangling_psm = [i for i in results["feature"].issues if i.check == "dangling_psm_id"]
    assert len(dangling_feature) == 1
    assert dangling_feature[0].severity == "warning"
    assert "999999" in dangling_feature[0].message
    assert len(dangling_psm) == 1
    assert dangling_psm[0].severity == "warning"
    assert "888888" in dangling_psm[0].message
    # Warnings must NOT fail validation.
    assert results["psm"].is_valid
    assert results["feature"].is_valid

    with Dataset(ds_dir) as ds:
        strict_results = ds.validate(structures=["feature", "psm"], strict=True)
    strict_dangling = [
        issue
        for result in strict_results.values()
        for issue in result.issues
        if issue.check in ("dangling_feature_id", "dangling_psm_id")
    ]
    assert len(strict_dangling) == 2
    assert all(issue.severity == "error" for issue in strict_dangling)
    assert not strict_results["psm"].is_valid
    assert not strict_results["feature"].is_valid


def test_feature_psm_reciprocal_desync_is_warning(tmp_path):
    """A psm pointing at a feature whose non-empty psm_ids does not list it back
    is flagged as a reciprocal desync warning (not a hard error)."""
    from qpx.dataset import Dataset

    _, feature_id, psm_id = _write_feature_psm_dataset(tmp_path / "probe", feature_psm_ids=None, psm_feature_id=None)
    # feature.psm_ids lists feature_id (a real feature id, NOT this psm's id) and
    # psm.feature_id resolves to the feature -> the edge exists one way but is not
    # reciprocated.
    ds_dir, _, _ = _write_feature_psm_dataset(tmp_path / "desync", feature_psm_ids=[feature_id], psm_feature_id=feature_id)
    with Dataset(ds_dir) as ds:
        results = ds.validate(structures=["feature", "psm"])
    desync = [i for i in results["psm"].issues if i.check == "feature_psm_desync"]
    assert len(desync) == 1
    assert desync[0].severity == "warning"
    assert results["psm"].is_valid  # warnings do not fail validation

    with Dataset(ds_dir) as ds:
        strict_results = ds.validate(structures=["feature", "psm"], strict=True)
    strict_desync = [i for i in strict_results["psm"].issues if i.check == "feature_psm_desync"]
    assert len(strict_desync) == 1
    assert strict_desync[0].severity == "error"
    assert not strict_results["psm"].is_valid


def test_feature_psm_inverse_reciprocal_desync_is_warning(tmp_path):
    """A feature listing a PSM that points to another existing feature is
    flagged when both directions are populated."""
    from qpx.core.data.identity import derive_id
    from qpx.dataset import Dataset
    from qpx.writers import FeatureWriter, PsmWriter
    from tests.conftest import make_feature_record, make_psm_record

    ds_dir = tmp_path / "inverse-desync"
    ds_dir.mkdir()
    feature_a = make_feature_record(sequence="FEATUREAK", peptidoform="FEATUREAK")
    feature_b = make_feature_record(sequence="FEATUREBK", peptidoform="FEATUREBK")
    feature_a["rt"] = 10.0
    feature_b["rt"] = 20.0
    psm = make_psm_record(sequence="PSMPEPTIDEK", peptidoform="PSMPEPTIDEK", scan=[10])
    feature_a_id = derive_id([feature_a["peptidoform"], feature_a["charge"], feature_a["run_file_name"], feature_a["rt"]])
    feature_b_id = derive_id([feature_b["peptidoform"], feature_b["charge"], feature_b["run_file_name"], feature_b["rt"]])
    psm_id = derive_id([psm["peptidoform"], psm["charge"], psm["run_file_name"], psm["scan"]])
    feature_a["psm_ids"] = [psm_id]
    psm["feature_id"] = feature_b_id

    with FeatureWriter(ds_dir / "exp.feature.parquet") as writer:
        writer.write_batch([feature_a, feature_b])
    with PsmWriter(ds_dir / "exp.psm.parquet") as writer:
        writer.write_batch([psm])

    with Dataset(ds_dir) as dataset:
        results = dataset.validate(structures=["feature", "psm"])
    issues = [issue for issue in results["feature"].issues if issue.check == "feature_psm_desync"]
    assert len(issues) == 1
    assert issues[0].severity == "warning"
    assert str(feature_a_id) in issues[0].message
    assert str(feature_b_id) in issues[0].message

    with Dataset(ds_dir) as dataset:
        strict_results = dataset.validate(structures=["feature", "psm"], strict=True)
    strict_issues = [issue for issue in strict_results["feature"].issues if issue.check == "feature_psm_desync"]
    assert len(strict_issues) == 1
    assert strict_issues[0].severity == "error"
    assert not strict_results["feature"].is_valid


# ---------------------------------------------------------------------------
# feature->pg cross-reference referential invariant
# ---------------------------------------------------------------------------


def _write_feature_pg_dataset(ds_dir, feature_pg_ids):
    """Write one feature and its candidate protein-group references."""
    from qpx.core.data.identity import derive_id
    from qpx.writers import FeatureWriter, PgWriter
    from tests.conftest import make_feature_record, make_pg_record

    ds_dir.mkdir(parents=True, exist_ok=True)
    feature = make_feature_record()
    feature["pg_ids"] = feature_pg_ids
    pg = make_pg_record()
    with FeatureWriter(ds_dir / "exp.feature.parquet") as writer:
        writer.write_batch([feature])
    with PgWriter(ds_dir / "exp.pg.parquet") as writer:
        writer.write_batch([pg])
    pg_id = derive_id(
        [pg["pg_accessions"], pg["grouped_runs"], pg["intensities"][0]["label"]],
        unordered_list_indices=(0, 1),
    )
    return ds_dir, pg_id


def test_feature_pg_cross_refs_resolve_cleanly(tmp_path):
    """A feature.pg_ids entry resolving to pg.pg_id produces no issue."""
    from qpx.dataset import Dataset

    _, pg_id = _write_feature_pg_dataset(tmp_path / "probe", None)
    ds_dir, _ = _write_feature_pg_dataset(tmp_path / "good", [pg_id])
    with Dataset(ds_dir) as dataset:
        results = dataset.validate(structures=["feature", "pg"], strict=True)

    assert not [issue for issue in results["feature"].issues if issue.check == "dangling_pg_id"]
    assert results["feature"].is_valid


def test_missing_optional_feature_pg_cross_ref_is_valid(tmp_path):
    """An external Feature file may omit optional pg_ids entirely."""
    from qpx.dataset import Dataset

    ds_dir, _ = _write_feature_pg_dataset(tmp_path / "without-pg-cross-ref", None)
    _drop_parquet_columns(ds_dir / "exp.feature.parquet", ["pg_ids"])

    with Dataset(ds_dir) as dataset:
        results = dataset.validate(structures=["feature", "pg"], strict=True)

    assert results["feature"].is_valid
    assert not [issue for issue in results["feature"].issues if issue.check == "referential_check_skipped"]


def test_feature_pg_cross_refs_dangling_warn_and_fail_strict(tmp_path):
    """A missing pg.pg_id is a warning normally and an error in strict mode."""
    from qpx.dataset import Dataset

    ds_dir, _ = _write_feature_pg_dataset(tmp_path / "bad", [999999])
    with Dataset(ds_dir) as dataset:
        results = dataset.validate(structures=["feature", "pg"])
    dangling = [issue for issue in results["feature"].issues if issue.check == "dangling_pg_id"]
    assert len(dangling) == 1
    assert dangling[0].severity == "warning"
    assert results["feature"].is_valid

    with Dataset(ds_dir) as dataset:
        strict_results = dataset.validate(structures=["feature", "pg"], strict=True)
    strict_dangling = [issue for issue in strict_results["feature"].issues if issue.check == "dangling_pg_id"]
    assert len(strict_dangling) == 1
    assert strict_dangling[0].severity == "error"
    assert not strict_results["feature"].is_valid
