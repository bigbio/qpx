"""DIA-NN integration test — PXD019909 small dataset.

Runs the full DIA-NN conversion pipeline (feature + PG + ontology +
provenance + dataset) once, then opens the output as a Dataset and
validates structures, schemas, and queryability.

Files:
    tests/examples/diann/small/diann_report.tsv      (~1 000 rows)
    tests/examples/diann/small/PXD019909-DIA.sdrf.tsv (42 runs)
    tests/examples/diann/small/diann_report.pg_matrix.tsv
    tests/examples/diann/small/mzml/                  (mzml_info TSVs)

Markers:
    @pytest.mark.integration
"""

import re
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.dataset import Dataset
from tests.conftest import assert_pg_table_wellformed

EXAMPLES_DIR = Path(__file__).parent.parent / "examples" / "diann" / "small"

_REPORT = EXAMPLES_DIR / "diann_report.tsv"
_SDRF = EXAMPLES_DIR / "PXD019909-DIA.sdrf.tsv"
_PG_MATRIX = EXAMPLES_DIR / "diann_report.pg_matrix.tsv"
_MZML_DIR = EXAMPLES_DIR / "mzml"

_PREFIX = "diann_small"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _parse_scan_number(spectrum_id: str) -> int:
    m = re.search(r"scan=(\d+)", str(spectrum_id))
    return int(m.group(1)) if m else 0


def _prepare_ms_info_parquet(source_dir: Path, dest_dir: Path) -> None:
    """Convert *_mzml_info.tsv files into *_ms_info.parquet files."""
    for tsv_path in source_dir.glob("*_mzml_info.tsv"):
        df = pd.read_csv(tsv_path, sep="\t")
        out = pd.DataFrame(
            {
                "rt": df["Retention_Time"].astype(float),
                "scan": df["SpectrumID"].apply(_parse_scan_number).astype(int),
                "precursor_mz": pd.to_numeric(df["Exp_Mass_To_Charge"], errors="coerce"),
            }
        )
        stem = tsv_path.stem.replace("_mzml_info", "")
        pq.write_table(pa.Table.from_pandas(out), str(dest_dir / f"{stem}_ms_info.parquet"))


# ---------------------------------------------------------------------------
# Module-scoped fixture: convert once, open as Dataset
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """Run full DIA-NN conversion and open as Dataset."""
    if not _REPORT.exists():
        pytest.skip(f"Test data not found: {_REPORT}")

    from qpx.converters.diann.converter import DiaNNConverter

    output = tmp_path_factory.mktemp("diann_small")

    # Prepare mzml info parquet
    mzml_parquet_dir = tmp_path_factory.mktemp("mzml_info")
    _prepare_ms_info_parquet(_MZML_DIR, mzml_parquet_dir)

    converter = DiaNNConverter(
        report_path=str(_REPORT),
        sdrf_path=str(_SDRF),
    )
    converter.convert_features(
        mzml_info_folder=str(mzml_parquet_dir),
        qvalue_threshold=0.05,
        output_folder=str(output),
        output_prefix=_PREFIX,
    )
    converter.convert_pg(
        pg_matrix_path=str(_PG_MATRIX),
        output_folder=str(output),
        output_prefix=_PREFIX,
    )
    converter.write_ontology(output, prefix=_PREFIX)
    converter.write_provenance(output, prefix=_PREFIX)
    converter.write_dataset(output, prefix=_PREFIX, project_accession="PXD019909")

    ds = Dataset(output)
    yield ds
    ds.close()


@pytest.fixture(scope="module")
def feature_table(dataset):
    """Raw PyArrow table for feature.parquet."""
    path = Path(dataset.path) / f"{_PREFIX}.feature.parquet"
    assert path.exists(), "feature.parquet was not produced"
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def pg_table(dataset):
    """Raw PyArrow table for pg.parquet."""
    path = Path(dataset.path) / f"{_PREFIX}.pg.parquet"
    assert path.exists(), "pg.parquet was not produced"
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# Feature conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestFeatureConversion:
    """Validate feature.parquet output."""

    def test_has_rows(self, feature_table):
        assert feature_table.num_rows > 0

    def test_sequences_are_nonempty(self, feature_table):
        seqs = feature_table.column("sequence").to_pylist()
        assert all(isinstance(s, str) and len(s) > 0 for s in seqs)

    def test_charges_in_range(self, feature_table):
        charges = feature_table.column("charge").to_pylist()
        assert all(1 <= c <= 10 for c in charges), f"Bad charges: {set(charges)}"

    def test_intensities_nonnegative(self, feature_table):
        for row in feature_table.column("intensities").to_pylist():
            if row is None:
                continue
            for entry in row:
                assert entry["intensity"] >= 0

    def test_run_file_names_nonempty(self, feature_table):
        names = feature_table.column("run_file_name").to_pylist()
        assert all(isinstance(n, str) and len(n) > 0 for n in names)

    def test_schema_validation(self, feature_table):
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Protein group conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestPgConversion:
    """Validate pg.parquet output."""

    def test_has_rows(self, pg_table):
        assert pg_table.num_rows > 0

    def test_pg_table_wellformed(self, pg_table):
        """grouped_runs, pg_accessions, anchor_protein, intensities are all well-formed."""
        assert_pg_table_wellformed(pg_table)

    def test_schema_validation(self, pg_table):
        from qpx.core.data import PgSchema

        errors = PgSchema.validate(pg_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Dataset structure tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestDatasetStructure:
    """Validate that all expected structures are present."""

    def test_has_feature(self, dataset):
        assert dataset.feature is not None
        assert dataset.feature.count() > 0

    def test_has_pg(self, dataset):
        assert dataset.pg is not None
        assert dataset.pg.count() > 0

    def test_has_ontology(self, dataset):
        assert dataset.ontology is not None
        assert dataset.ontology.count() > 0

    def test_has_provenance(self, dataset):
        assert dataset.provenance is not None
        assert dataset.provenance.count() > 0

    def test_has_dataset_meta(self, dataset):
        assert dataset.dataset_meta is not None

    def test_available_structures(self, dataset):
        available = set(dataset.available_structures)
        expected = {"feature", "pg", "ontology", "provenance", "dataset"}
        assert expected.issubset(available), f"Missing: {expected - available}"


# ---------------------------------------------------------------------------
# Feature querying via Dataset API
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestFeatureQuerying:
    """Query features through the Dataset API."""

    def test_filter_by_charge(self, dataset):
        c2 = dataset.feature.filter("charge = 2").count()
        c3 = dataset.feature.filter("charge = 3").count()
        total = dataset.feature.count()
        assert c2 > 0
        assert c2 + c3 <= total

    def test_select_columns(self, dataset):
        df = dataset.feature.select("sequence", "charge").limit(10).to_df()
        assert set(df.columns) == {"sequence", "charge"}
        assert len(df) == 10

    def test_unique_sequences_via_sql(self, dataset):
        result = dataset.sql("SELECT COUNT(DISTINCT sequence) AS n FROM feature")
        n = result.to_df()["n"].iloc[0]
        assert n > 10

    def test_feature_to_dataframe(self, dataset):
        df = dataset.feature.limit(50).to_df()
        assert len(df) <= 50
        assert "sequence" in df.columns


# ---------------------------------------------------------------------------
# PG querying via Dataset API
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestPgQuerying:
    """Query protein groups through the Dataset API."""

    def test_unique_anchor_proteins(self, dataset):
        result = dataset.sql("SELECT COUNT(DISTINCT anchor_protein) AS n FROM pg")
        n = result.to_df()["n"].iloc[0]
        assert n > 10

    def test_pg_to_dataframe(self, dataset):
        df = dataset.pg.limit(20).to_df()
        assert len(df) <= 20
        assert "anchor_protein" in df.columns
        assert "pg_accessions" in df.columns


# ---------------------------------------------------------------------------
# Provenance tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestProvenance:
    """Validate provenance records."""

    def test_has_diann_step(self, dataset):
        df = dataset.provenance.to_df()
        assert "tool_name" in df.columns
        tool_names = df["tool_name"].tolist()
        assert any("DIA-NN" in str(t) for t in tool_names if t)

    def test_has_qpx_step(self, dataset):
        df = dataset.provenance.to_df()
        tool_names = df["tool_name"].tolist()
        assert any("qpx" in str(t).lower() for t in tool_names if t)

    def test_step_order_positive(self, dataset):
        df = dataset.provenance.to_df()
        assert (df["step_order"] > 0).all()


# ---------------------------------------------------------------------------
# Ontology tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestOntology:
    """Validate ontology entries."""

    def test_has_entries(self, dataset):
        df = dataset.ontology.to_df()
        assert len(df) > 0

    def test_has_field_name(self, dataset):
        df = dataset.ontology.to_df()
        assert "field_name" in df.columns

    def test_view_is_feature(self, dataset):
        df = dataset.ontology.to_df()
        views = set(df["view"].tolist())
        assert "feature" in views


# ---------------------------------------------------------------------------
# Dataset metadata tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestDatasetMetadata:
    """Validate dataset.parquet metadata."""

    def test_project_accession(self, dataset):
        df = dataset.dataset_meta.to_df()
        assert df["project_accession"].iloc[0] == "PXD019909"

    def test_software_name(self, dataset):
        df = dataset.dataset_meta.to_df()
        assert df["software_name"].iloc[0] == "DIA-NN"

    def test_creation_date(self, dataset):
        df = dataset.dataset_meta.to_df()
        assert df["creation_date"].iloc[0] is not None
