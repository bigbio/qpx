"""Tests for the OpenMS QPX converter (enrichment of -out_qpx output)."""

from __future__ import annotations

import logging
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.converters.openms.converter import OpenMSConverter, _discover_parquet
from qpx.core.data.loader import load_schema
from qpx.dataset import Dataset
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter
from tests.conftest import (
    make_feature_record,
    make_pg_record,
    make_psm_record,
    write_lfq_consensusxml,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_minimal_qpx(
    qpx_dir: Path,
    prefix: str = "quantms",
    feature_label: str = "TMT126",
    pg_label: str = "TMT126",
) -> None:
    """Write minimal QPX-compliant psm, feature, pg parquet files."""
    # PSM
    psm_records = [
        make_psm_record(sequence="PEPTIDEK", run_file_name="run_01", scan=[1001]),
        make_psm_record(sequence="ANOTHERK", run_file_name="run_01", scan=[1002], charge=3),
    ]
    with PsmWriter(qpx_dir / f"{prefix}.psm.parquet") as w:
        w.write_batch(psm_records)

    # Feature
    feat_records = [
        make_feature_record(
            sequence="PEPTIDEK",
            run_file_name="run_01",
            anchor_protein="P12345",
            intensities=[{"label": feature_label, "intensity": 1000.0}],
        ),
        make_feature_record(
            sequence="ANOTHERK",
            run_file_name="run_01",
            anchor_protein="P12345",
            charge=3,
            intensities=[{"label": feature_label, "intensity": 2000.0}],
        ),
    ]
    with FeatureWriter(qpx_dir / f"{prefix}.feature.parquet") as w:
        w.write_batch(feat_records)

    # PG
    pg_records = [
        make_pg_record(
            anchor_protein="P12345",
            run_file_name="run_01",
            intensities=[{"label": pg_label, "intensity": 5000.0}],
        ),
    ]
    with PgWriter(qpx_dir / f"{prefix}.pg.parquet") as w:
        w.write_batch(pg_records)


def _downgrade_to_legacy_qpx(qpx_dir: Path, prefix: str = "quantms") -> None:
    """Remove 1.1 identity/PG fields to model existing OpenMS -out_qpx files."""
    for view, id_columns in {
        "feature": ("feature_id", "psm_ids", "pg_ids"),
        "psm": ("psm_id", "feature_id"),
    }.items():
        path = qpx_dir / f"{prefix}.{view}.parquet"
        table = pq.read_table(path).drop_columns(id_columns).replace_schema_metadata(None)
        pq.write_table(table, path)

    pg_path = qpx_dir / f"{prefix}.pg.parquet"
    legacy_rows = []
    for row in pq.read_table(pg_path).to_pylist():
        run_file_name = row.pop("grouped_runs")[0]
        row.pop("pg_id")
        label = row.pop("label")
        intensity = row.pop("intensity")
        row["run_file_name"] = run_file_name
        row["intensities"] = [{"label": label, "intensity": intensity}]
        legacy_rows.append(row)
    pq.write_table(pa.Table.from_pylist(legacy_rows), pg_path)


def _write_minimal_sdrf(
    sdrf_path: Path,
    label: str = "AC=MS:1002038;NT=label free sample",
) -> None:
    """Write a minimal SDRF TSV for testing."""
    lines = [
        "source name\tcharacteristics[organism]\tcharacteristics[organism part]\t"
        "comment[data file]\tcomment[label]\tcomment[instrument]\t"
        "comment[cleavage agent details]\tcomment[fraction identifier]",
        f"sample_1\tHomo sapiens\tliver\trun_01.raw\t{label}\tAC=MS:1000449;NT=LTQ Orbitrap\tAC=MS:1001251;NT=Trypsin\t1",
    ]
    sdrf_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Tests: _discover_parquet
# ---------------------------------------------------------------------------


class TestDiscoverParquet:
    def test_discovers_all_three(self, tmp_path):
        _write_minimal_qpx(tmp_path)
        found = _discover_parquet(tmp_path)
        assert set(found.keys()) == {"psm", "feature", "pg"}

    def test_empty_dir(self, tmp_path):
        found = _discover_parquet(tmp_path)
        assert not found

    def test_partial(self, tmp_path):
        """Only psm present."""
        psm_records = [make_psm_record()]
        with PsmWriter(tmp_path / "test.psm.parquet") as w:
            w.write_batch(psm_records)
        found = _discover_parquet(tmp_path)
        assert set(found.keys()) == {"psm"}


# ---------------------------------------------------------------------------
# Tests: OpenMSConverter
# ---------------------------------------------------------------------------


class TestOpenMSConverter:
    def test_convert_upgrades_legacy_core_tables(self, tmp_path):
        """Existing OpenMS -out_qpx files receive IDs and the flattened PG shape."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)
        _downgrade_to_legacy_qpx(qpx_dir)

        output = tmp_path / "output"
        OpenMSConverter(qpx_dir=qpx_dir).convert(output_folder=output, output_prefix="openms")

        for view, id_field in (("feature", "feature_id"), ("psm", "psm_id"), ("pg", "pg_id")):
            table = pq.read_table(output / f"openms.{view}.parquet")
            assert id_field in table.column_names
            assert table.column(id_field).null_count == 0
            assert load_schema(view).validate_full(table, strict=True).is_valid
        pg_table = pq.read_table(output / "openms.pg.parquet")
        assert "grouped_runs" in pg_table.column_names
        assert "run_file_name" not in pg_table.column_names

    def test_duplicate_identity_warns_and_replaces_output(self, tmp_path, caplog):
        """A duplicate primary key is tolerated with a warning; the core file is rewritten."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        records = [
            make_psm_record(sequence="PEPTIDEK", run_file_name="run_01", scan=[1001]),
            make_psm_record(sequence="PEPTIDEK", run_file_name="run_01", scan=[1001]),
        ]
        records[1]["posterior_error_probability"] = 0.2
        source = qpx_dir / "quantms.psm.parquet"
        current_schema = load_schema("psm").get_arrow_schema()
        legacy_schema = pa.schema([field for field in current_schema if field.name not in {"psm_id", "feature_id"}])
        legacy = pa.Table.from_pylist(records, schema=legacy_schema)
        pq.write_table(legacy, source)

        output = tmp_path / "output"
        output.mkdir()
        destination = output / "openms.psm.parquet"
        sentinel = b"pre-existing output"
        destination.write_bytes(sentinel)

        with caplog.at_level(logging.WARNING):
            OpenMSConverter(qpx_dir=qpx_dir).convert(output_folder=output, output_prefix="openms")

        # The converter completed: the destination is replaced with the real parquet output.
        assert destination.read_bytes() != sentinel
        rewritten = pq.read_table(destination)
        assert rewritten.num_rows == 2
        assert "duplicate row" in caplog.text.lower()
        assert not list(output.glob(".openms.psm.parquet.*.tmp"))

    def test_convert_full_bundle(self, tmp_path):
        """Full enrichment: 3 core + SDRF -> 8 files."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path)
        converter.convert(output_folder=output, output_prefix="openms")

        # Check all 8 files exist
        expected_files = [
            "openms.psm.parquet",
            "openms.feature.parquet",
            "openms.pg.parquet",
            "openms.sample.parquet",
            "openms.run.parquet",
            "openms.ontology.parquet",
            "openms.provenance.parquet",
            "openms.dataset.parquet",
        ]
        for fname in expected_files:
            assert (output / fname).exists(), f"Missing: {fname}"

    def test_core_tables_validated(self, tmp_path):
        """Core tables should pass QPX schema validation after copy."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path)
        converter.convert(output_folder=output, output_prefix="openms")

        for view in ("psm", "feature", "pg"):
            schema = load_schema(view)
            table = pq.read_table(str(output / f"openms.{view}.parquet"))
            result = schema.validate_full(table)
            assert result.is_valid, f"{view} validation failed: {result.summary}"

    def test_provenance_structure(self, tmp_path):
        """Provenance should have OpenMS + QPX steps."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path)
        converter.convert(output_folder=output, output_prefix="openms")

        prov_table = pq.read_table(str(output / "openms.provenance.parquet"))
        assert prov_table.num_rows == 2
        tools = prov_table.column("tool_name").to_pylist()
        assert "OpenMS/ProteomicsLFQ" in tools
        assert "qpx" in tools

    def test_dataset_metadata(self, tmp_path):
        """Dataset should contain project accession."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path)
        converter.convert(
            output_folder=output,
            output_prefix="openms",
            project_accession="PXD001819",
        )

        ds_table = pq.read_table(str(output / "openms.dataset.parquet"))
        assert ds_table.num_rows == 1
        assert ds_table.column("project_accession").to_pylist()[0] == "PXD001819"
        assert ds_table.column("software_name").to_pylist()[0] == "OpenMS/ProteomicsLFQ"

    @pytest.mark.parametrize(
        ("label", "canonical_label"),
        [
            ("NT=tmt126;AC=MS:1002038", "TMT126"),
            ("AC=MS:1002038;NT=itraq114", "ITRAQ114"),
        ],
    )
    def test_isobaric_labels_join_and_metadata(self, tmp_path, label, canonical_label):
        """Ontology labels must join and must not claim the label-free workflow."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir, feature_label="run_01.raw", pg_label="1")

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path, label=label)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path)
        converter.convert(output_folder=output, output_prefix="openms")

        provenance = pq.read_table(output / "openms.provenance.parquet").to_pylist()
        assert provenance[0]["step_name"] == "isobaric_quantification"
        assert provenance[0]["tool_name"] == "OpenMS/IsobaricWorkflow"

        dataset = pq.read_table(output / "openms.dataset.parquet").to_pylist()
        assert dataset[0]["software_name"] == "OpenMS/IsobaricWorkflow"

        feature = pq.read_table(output / "openms.feature.parquet").to_pylist()
        assert {entry["label"] for row in feature for entry in row["intensities"]} == {canonical_label}
        # pg is flattened since 1.1: scalar label column (one row per label).
        pg = pq.read_table(output / "openms.pg.parquet").to_pylist()
        assert {row["label"] for row in pg if row["label"] is not None} == {canonical_label}

        with Dataset(output, file_prefix="openms", duckdb_threads=24) as qpx_dataset:
            assert len(qpx_dataset.intensity("peptide").to_df()) == 2
            assert len(qpx_dataset.intensity("protein").to_df()) == 1

    def test_no_qpx_files_raises(self, tmp_path):
        """Empty directory should raise FileNotFoundError."""
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        converter = OpenMSConverter(qpx_dir=empty_dir)
        with pytest.raises(FileNotFoundError, match="No QPX parquet files"):
            converter.convert(output_folder=tmp_path / "out")

    def test_same_dir_skips_copy(self, tmp_path):
        """When output_folder == qpx_dir, core files should not be copied."""
        _write_minimal_qpx(tmp_path, prefix="openms")

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        converter = OpenMSConverter(qpx_dir=tmp_path, sdrf_path=sdrf_path)
        converter.convert(output_folder=tmp_path, output_prefix="openms")

        # Core files should still exist and metadata files should be added
        assert (tmp_path / "openms.psm.parquet").exists()
        assert (tmp_path / "openms.provenance.parquet").exists()

    def test_no_sdrf_skips_sample_run(self, tmp_path):
        """Without SDRF, sample/run are skipped and source labels are preserved."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir)
        converter.convert(output_folder=output, output_prefix="openms")

        # Core + provenance + dataset + ontology should exist
        assert (output / "openms.psm.parquet").exists()
        assert (output / "openms.provenance.parquet").exists()
        assert (output / "openms.dataset.parquet").exists()
        # sample/run should NOT exist
        assert not (output / "openms.sample.parquet").exists()
        assert not (output / "openms.run.parquet").exists()
        feature = pq.read_table(output / "openms.feature.parquet")
        labels = [entry["label"] for row in feature.column("intensities").to_pylist() for entry in row]
        assert labels == ["TMT126", "TMT126"]

    def test_metadata_tables_validate(self, tmp_path):
        """All metadata tables should pass QPX schema validation."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path)
        converter.convert(output_folder=output, output_prefix="openms")

        for view in ("sample", "run", "ontology", "provenance", "dataset"):
            path = output / f"openms.{view}.parquet"
            if path.exists():
                schema = load_schema(view)
                table = pq.read_table(str(path))
                result = schema.validate_full(table)
                assert result.is_valid, f"{view} validation failed: {result.summary}"


# ---------------------------------------------------------------------------
# Tests: fraction_group capture (consensusXML LFQ experimental design)
# ---------------------------------------------------------------------------


def _fraction_group_cv(cv_params):
    """Return the fraction_group cv_value from a row's cv_params, or None."""
    for param in cv_params or []:
        if param and param.get("cv_name") == "fraction_group":
            return param.get("cv_value")
    return None


class TestFractionGroupCapture:
    def test_add_fraction_group_cv_params_pg(self, tmp_path):
        """The merged rewrite appends a fraction_group cv_param keyed by grouped_runs.

        The standalone ``_add_fraction_group_cv_params`` pass was folded into the
        single copy/relabel rewrite (``relabel_intensities_parquet`` with a
        cv_param resolver); this exercises the cv_param-only branch (relabel off).
        """
        from qpx.converters.channel_labels import (
            fraction_groups_from_consensusxml,
            relabel_intensities_parquet,
        )
        from qpx.converters.openms.converter import (
            FRACTION_GROUP_CV_NAME,
            _fraction_group_for,
            _fraction_group_lookup,
        )

        pg_path = tmp_path / "x.pg.parquet"
        with PgWriter(pg_path) as w:
            w.write_batch([make_pg_record(anchor_protein="P1", run_file_name="run_A")])

        cxml = tmp_path / "lfq.consensusXML"
        write_lfq_consensusxml(cxml, [(0, "run_A.mzML", "1", "1", "1")])
        lookup = _fraction_group_lookup(fraction_groups_from_consensusxml(str(cxml)))

        out_path = tmp_path / "x.out.pg.parquet"
        annotated = relabel_intensities_parquet(
            str(pg_path),
            str(out_path),
            {},
            is_lfq=False,
            compression="zstd",
            relabel=False,
            cv_param_name=FRACTION_GROUP_CV_NAME,
            run_column="grouped_runs",
            cv_param_resolver=lambda run_names: _fraction_group_for(run_names, lookup),
        )
        assert annotated == 1
        table = pq.read_table(str(out_path))
        assert _fraction_group_cv(table.column("cv_params").to_pylist()[0]) == "1"
        # Still a valid pg table after annotation.
        assert load_schema("pg").validate_full(table).is_valid

    def test_convert_captures_fraction_group(self, tmp_path):
        """Full convert with a consensusXML annotates pg + feature rows."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        # LFQ run named run_01 (matches the SDRF data file run_01.raw).
        _write_minimal_qpx(qpx_dir, feature_label="run_01", pg_label="1")

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)  # label free -> LFQ

        cxml = tmp_path / "lfq.consensusXML"
        write_lfq_consensusxml(cxml, [(0, "run_01.mzML", "3", "1", "1")])

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path, consensusxml_path=cxml)
        converter.convert(output_folder=output, output_prefix="openms")

        pg_table = pq.read_table(str(output / "openms.pg.parquet"))
        assert _fraction_group_cv(pg_table.column("cv_params").to_pylist()[0]) == "3"
        feat_table = pq.read_table(str(output / "openms.feature.parquet"))
        assert all(_fraction_group_cv(row) == "3" for row in feat_table.column("cv_params").to_pylist())
        # Annotated tables still validate against the QPX schema.
        assert load_schema("pg").validate_full(pg_table).is_valid
        assert load_schema("feature").validate_full(feat_table).is_valid

    def test_convert_without_consensusxml_no_fraction_group(self, tmp_path):
        """No consensusXML -> no fraction_group cv_param is added."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir, feature_label="run_01", pg_label="1")

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path)
        converter.convert(output_folder=output, output_prefix="openms")

        pg_table = pq.read_table(str(output / "openms.pg.parquet"))
        assert _fraction_group_cv(pg_table.column("cv_params").to_pylist()[0]) is None
