"""Spectronaut converter integration tests using small example dataset.

Runs the full Spectronaut conversion pipeline once (module scope) and
validates that feature.parquet, pg.parquet, and ontology.parquet are
written with correct schemas and plausible values.
"""

from pathlib import Path

import pyarrow.parquet as pq
import pytest

from tests.conftest import assert_pg_table_wellformed

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples" / "spectronaut" / "small"

_REPORT = EXAMPLES_DIR / "spectronaut_report.tsv"


# ---------------------------------------------------------------------------
# Module-scoped fixture: run conversion once, share outputs across tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def converted_output(tmp_path_factory):
    """Run Spectronaut conversion once for all tests in this module."""
    from qpx.converters.spectronaut.converter import SpectronautConverter

    if not _REPORT.exists():
        pytest.skip(f"Test data not found: {_REPORT}")

    output = tmp_path_factory.mktemp("spectronaut_output")

    converter = SpectronautConverter(report_path=str(_REPORT))
    converter.convert_features(
        qvalue_threshold=1.0,
        output_folder=str(output),
        output_prefix="sn_test",
    )
    converter.convert_pg(
        output_folder=str(output),
        output_prefix="sn_test",
    )
    converter.write_ontology(output_folder=output, prefix="sn_test")
    converter.write_provenance(output_folder=output, prefix="sn_test")
    converter.write_dataset(output_folder=output, prefix="sn_test")

    return output


@pytest.fixture(scope="module")
def feature_table(converted_output):
    """Read the feature.parquet produced by the converter."""
    path = converted_output / "sn_test.feature.parquet"
    if not path.exists():
        pytest.skip("feature.parquet was not produced")
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def pg_table(converted_output):
    """Read the pg.parquet produced by the converter."""
    path = converted_output / "sn_test.pg.parquet"
    if not path.exists():
        pytest.skip("pg.parquet was not produced")
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# Feature conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestSpectronautFeatureConversion:
    """Validate feature.parquet output from Spectronaut conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / "sn_test.feature.parquet"
        assert path.exists(), "feature.parquet was not created"

    def test_has_rows(self, feature_table):
        assert feature_table.num_rows > 0, "feature.parquet is empty"

    def test_sequence_values_are_nonempty(self, feature_table):
        sequences = feature_table.column("sequence").to_pylist()
        for seq in sequences:
            assert isinstance(seq, str) and len(seq) > 0

    def test_charge_values_are_valid(self, feature_table):
        charges = feature_table.column("charge").to_pylist()
        for charge in charges:
            assert 1 <= charge <= 10, f"Charge out of range: {charge}"

    def test_intensities_are_nonnegative(self, feature_table):
        rows = feature_table.column("intensities").to_pylist()
        for row_intensities in rows:
            assert row_intensities is not None and len(row_intensities) > 0
            for entry in row_intensities:
                assert entry["intensity"] >= 0, f"Negative intensity: {entry['intensity']}"

    def test_run_file_names_are_nonempty(self, feature_table):
        run_names = feature_table.column("run_file_name").to_pylist()
        for name in run_names:
            assert isinstance(name, str) and len(name) > 0

    def test_two_runs_present(self, feature_table):
        run_names = set(feature_table.column("run_file_name").to_pylist())
        assert len(run_names) == 2, f"Expected 2 runs, got {len(run_names)}"

    def test_scan_is_empty_list(self, feature_table):
        scans = feature_table.column("scan").to_pylist()
        for scan in scans:
            assert scan == [], f"Expected empty scan list, got {scan}"

    def test_schema_validation(self, feature_table):
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Protein group conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestSpectronautPgConversion:
    """Validate pg.parquet output from Spectronaut conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / "sn_test.pg.parquet"
        assert path.exists(), "pg.parquet was not created"

    def test_has_rows(self, pg_table):
        assert pg_table.num_rows > 0, "pg.parquet is empty"

    def test_pg_table_wellformed(self, pg_table):
        """grouped_runs, pg_accessions, anchor_protein, intensities are all well-formed."""
        assert_pg_table_wellformed(pg_table)

    def test_schema_validation(self, pg_table):
        from qpx.core.data import PgSchema

        errors = PgSchema.validate(pg_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Ontology conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestSpectronautOntologyConversion:
    """Validate ontology.parquet output from Spectronaut conversion."""

    def test_ontology_file_exists_or_no_scores(self, converted_output):
        path = converted_output / "sn_test.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")
        table = pq.read_table(str(path))
        assert table.num_rows > 0, "ontology.parquet exists but is empty"

    def test_ontology_columns(self, converted_output):
        path = converted_output / "sn_test.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")
        table = pq.read_table(str(path))
        column_names = set(table.column_names)
        expected = {"field_name", "view"}
        missing = expected - column_names
        assert not missing, f"Missing ontology columns: {missing}"

    def test_schema_validation(self, converted_output):
        from qpx.core.data import OntologySchema

        path = converted_output / "sn_test.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")
        table = pq.read_table(str(path))
        errors = OntologySchema.validate(table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Provenance tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestSpectronautProvenance:
    """Validate provenance.parquet output."""

    def test_provenance_exists(self, converted_output):
        path = converted_output / "sn_test.provenance.parquet"
        assert path.exists(), "provenance.parquet was not created"

    def test_has_spectronaut_step(self, converted_output):
        path = converted_output / "sn_test.provenance.parquet"
        if not path.exists():
            pytest.skip("provenance.parquet not written")
        table = pq.read_table(str(path))
        tool_names = table.column("tool_name").to_pylist()
        assert any("Spectronaut" in str(t) for t in tool_names if t), "Expected Spectronaut in provenance tool_names"

    def test_has_qpx_step(self, converted_output):
        path = converted_output / "sn_test.provenance.parquet"
        if not path.exists():
            pytest.skip("provenance.parquet not written")
        table = pq.read_table(str(path))
        tool_names = table.column("tool_name").to_pylist()
        assert any("qpx" in str(t).lower() for t in tool_names if t), "Expected qpx in provenance tool_names"


class TestRunNameNormalization:
    """Regression for bigbio/qpx#252: the Spectronaut feature adapter stripped
    only ``\\.(mzML|raw|d)$`` case-sensitively while the pg adapter stripped
    ``(?i)\\.(mzML|raw|d|wiff|htrms)$``, so ``.wiff`` / ``.htrms`` / ``.RAW``
    broke feature<->pg<->run joins. Both adapters must strip the same set,
    case-insensitively, by DuckDB SQL."""

    @staticmethod
    def _run_expr(part: str) -> str:
        # part looks like "regexp_replace(...) AS <alias>"; return the expression.
        return part.rsplit(" AS ", 1)[0]

    def test_feature_and_pg_use_same_extension_strip(self):
        import duckdb

        from qpx.converters.spectronaut.feature_adapter import SpectronautFeatureAdapter as FeatureAdapter
        from qpx.converters.spectronaut.pg_adapter import SpectronautPgAdapter as PgAdapter

        r = {
            "sequence": "SEQ",
            "charge": "CHG",
            "run_file_name": "RUN",
            "pg_accessions": "PG",
            "intensity": "INT",
        }
        core = FeatureAdapter._build_core_select_parts(r, set())
        meta = FeatureAdapter._build_meta_select_parts(r, set())
        pg = PgAdapter._build_pg_select_parts(r, set())

        feat_run = next(self._run_expr(p) for p in core if p.endswith("AS run_file_name"))
        feat_id_run = next(self._run_expr(p) for p in meta if p.endswith("AS id_run_file_name"))
        pg_run = next(self._run_expr(p) for p in pg if p.endswith("AS run_file_name"))

        con = duckdb.connect()
        try:
            for name in ("sample.mzML", "sample.RAW", "sample.wiff", "sample.HTRMS", "sample.d"):
                # Evaluate each adapter's expression against the same input value.
                feat_val = con.execute(f'SELECT {feat_run} FROM (SELECT ? AS "RUN") r', [name]).fetchone()[0]
                feat_id_val = con.execute(f'SELECT {feat_id_run} FROM (SELECT ? AS "RUN") r', [name]).fetchone()[0]
                pg_val = con.execute(f'SELECT {pg_run} FROM (SELECT ? AS "RUN") r', [name]).fetchone()[0]
                assert feat_val == "sample", f"feature did not strip {name}: {feat_val}"
                assert feat_id_val == "sample", f"feature id_run did not strip {name}: {feat_id_val}"
                assert pg_val == "sample", f"pg did not strip {name}: {pg_val}"
                assert feat_val == pg_val == feat_id_val
        finally:
            con.close()
