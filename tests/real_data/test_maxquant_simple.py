"""MaxQuant integration test — simple dataset (PSM + Feature).

Runs the full MaxQuant conversion pipeline (PSM + Feature + ontology +
provenance + dataset) once, then opens the output as a Dataset and
validates structures, schemas, and queryability.

Files:
    tests/examples/maxquant/maxquant_simple/msms.txt      (1 000 PSMs)
    tests/examples/maxquant/maxquant_simple/evidence.txt   (9 301 features)
    tests/examples/maxquant/maxquant_simple/sdrf.tsv       (minimal SDRF)

Markers:
    @pytest.mark.integration
"""

from pathlib import Path

import pyarrow.parquet as pq
import pytest

from qpx.dataset import Dataset

EXAMPLES_DIR = Path(__file__).parent.parent / "examples" / "maxquant" / "maxquant_simple"

_MSMS = EXAMPLES_DIR / "msms.txt"
_EVIDENCE = EXAMPLES_DIR / "evidence.txt"
_SDRF = EXAMPLES_DIR / "sdrf.tsv"

_PREFIX = "mq_simple"


# ---------------------------------------------------------------------------
# Module-scoped fixture: convert once, open as Dataset
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """Run full MaxQuant conversion and open as Dataset."""
    if not _MSMS.exists() or not _EVIDENCE.exists():
        pytest.skip(f"Test data not found: {EXAMPLES_DIR}")

    from qpx.converters.maxquant.converter import MaxQuantConverter

    output = tmp_path_factory.mktemp("mq_simple")

    converter = MaxQuantConverter()
    converter.convert(
        output_folder=output,
        msms_file=_MSMS,
        evidence_file=_EVIDENCE,
        sdrf_file=_SDRF if _SDRF.exists() else None,
        output_prefix=_PREFIX,
        project_accession="test_simple",
    )

    ds = Dataset(output)
    yield ds
    ds.close()


@pytest.fixture(scope="module")
def psm_table(dataset):
    """Raw PyArrow table for psm.parquet."""
    path = Path(dataset.path) / f"{_PREFIX}.psm.parquet"
    assert path.exists(), "psm.parquet was not produced"
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def feature_table(dataset):
    """Raw PyArrow table for feature.parquet."""
    path = Path(dataset.path) / f"{_PREFIX}.feature.parquet"
    assert path.exists(), "feature.parquet was not produced"
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# PSM conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestPsmConversion:
    """Validate psm.parquet output."""

    def test_has_rows(self, psm_table):
        assert psm_table.num_rows > 0

    def test_key_columns_present(self, psm_table):
        expected = {"sequence", "charge", "calculated_mz", "run_file_name", "scan"}
        missing = expected - set(psm_table.column_names)
        assert not missing, f"Missing columns: {missing}"

    def test_sequences_nonempty(self, psm_table):
        seqs = psm_table.column("sequence").to_pylist()
        assert all(isinstance(s, str) and len(s) > 0 for s in seqs)

    def test_charges_in_range(self, psm_table):
        charges = psm_table.column("charge").to_pylist()
        assert all(1 <= c <= 10 for c in charges)

    def test_calculated_mz_nonnegative(self, psm_table):
        values = psm_table.column("calculated_mz").to_pylist()
        assert all(v is not None and v >= 0 for v in values)

    def test_schema_validation(self, psm_table):
        from qpx.core.data import PsmSchema

        errors = PsmSchema.validate(psm_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Feature conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestFeatureConversion:
    """Validate feature.parquet output."""

    def test_has_rows(self, feature_table):
        assert feature_table.num_rows > 0

    def test_key_columns_present(self, feature_table):
        expected = {
            "sequence",
            "charge",
            "calculated_mz",
            "intensities",
            "run_file_name",
        }
        missing = expected - set(feature_table.column_names)
        assert not missing, f"Missing columns: {missing}"

    def test_sequences_nonempty(self, feature_table):
        seqs = feature_table.column("sequence").to_pylist()
        assert all(isinstance(s, str) and len(s) > 0 for s in seqs)

    def test_charges_in_range(self, feature_table):
        charges = feature_table.column("charge").to_pylist()
        assert all(1 <= c <= 10 for c in charges)

    def test_intensities_nonnegative(self, feature_table):
        for row in feature_table.column("intensities").to_pylist():
            if row is None:
                continue
            for entry in row:
                assert entry["intensity"] >= 0

    def test_schema_validation(self, feature_table):
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Dataset structure tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestDatasetStructure:
    """Validate all expected structures."""

    def test_has_psm(self, dataset):
        assert dataset.psm is not None
        assert dataset.psm.count() > 0

    def test_has_feature(self, dataset):
        assert dataset.feature is not None
        assert dataset.feature.count() > 0

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
        expected = {"psm", "feature", "ontology", "provenance", "dataset"}
        assert expected.issubset(available), f"Missing: {expected - available}"


# ---------------------------------------------------------------------------
# PSM querying via Dataset API
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestPsmQuerying:
    """Query PSMs through the Dataset API."""

    def test_filter_by_charge(self, dataset):
        c2 = dataset.psm.filter("charge = 2").count()
        total = dataset.psm.count()
        assert c2 > 0
        assert c2 <= total

    def test_select_columns(self, dataset):
        df = dataset.psm.select("sequence", "charge").limit(10).to_df()
        assert set(df.columns) == {"sequence", "charge"}
        assert len(df) == 10

    def test_unique_sequences_via_sql(self, dataset):
        result = dataset.sql("SELECT COUNT(DISTINCT sequence) AS n FROM psm")
        n = result.to_df()["n"].iloc[0]
        assert n > 10


# ---------------------------------------------------------------------------
# Feature querying via Dataset API
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestFeatureQuerying:
    """Query features through the Dataset API."""

    def test_filter_by_charge(self, dataset):
        c2 = dataset.feature.filter("charge = 2").count()
        total = dataset.feature.count()
        assert c2 > 0
        assert c2 <= total

    def test_feature_to_dataframe(self, dataset):
        df = dataset.feature.limit(50).to_df()
        assert len(df) <= 50
        assert "sequence" in df.columns


# ---------------------------------------------------------------------------
# Provenance tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestProvenance:
    """Validate provenance records."""

    def test_has_maxquant_step(self, dataset):
        df = dataset.provenance.to_df()
        tool_names = df["tool_name"].tolist()
        assert any("MaxQuant" in str(t) for t in tool_names if t)

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


# ---------------------------------------------------------------------------
# Dataset metadata tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestDatasetMetadata:
    """Validate dataset.parquet metadata."""

    def test_project_accession(self, dataset):
        df = dataset.dataset_meta.to_df()
        assert df["project_accession"].iloc[0] == "test_simple"

    def test_software_name(self, dataset):
        df = dataset.dataset_meta.to_df()
        assert df["software_name"].iloc[0] == "MaxQuant"

    def test_creation_date(self, dataset):
        df = dataset.dataset_meta.to_df()
        assert df["creation_date"].iloc[0] is not None
