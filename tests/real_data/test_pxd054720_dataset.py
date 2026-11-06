"""Real-data integration test -- PXD054720 full Dataset pipeline.

Converts the compressed mzIdentML + MGF files from PRIDE project PXD054720 to a
complete QPX dataset using MzIdentMLConverter, then opens the result as a Dataset
and queries all structures and views.

Files:
    tests/examples/mzidentml/PXD054720/F001234.mzid.gz  (~1.2 MB, 9 707 spectra)
    tests/examples/mzidentml/PXD054720/F001234.mgf.gz   (~25 MB, raw spectra)

Markers:
    @pytest.mark.large_data -- CI runs it with the integration tests, not the unit tests
"""

import pytest

lxml = pytest.importorskip("lxml", reason="lxml required for mzIdentML tests")

from pathlib import Path

from qpx.converters.mzidentml.converter import MzIdentMLConverter
from qpx.dataset import Dataset

EXAMPLES = Path(__file__).parent.parent / "examples" / "mzidentml" / "PXD054720"
MZID_GZ = EXAMPLES / "F001234.mzid.gz"
MGF_GZ = EXAMPLES / "F001234.mgf.gz"


# ---------------------------------------------------------------------------
# Module-scoped fixture: convert once, share across all tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """Convert mzIdentML + MGF to a full QPX dataset and open it."""
    if not MZID_GZ.exists():
        pytest.skip(f"Test data not found: {MZID_GZ}")

    work_dir = tmp_path_factory.mktemp("pxd054720_dataset")

    # Run full conversion
    converter = MzIdentMLConverter()
    converter.convert(
        mzid_path=MZID_GZ,
        output_folder=work_dir,
        output_prefix="PXD054720",
        mgf_path=MGF_GZ if MGF_GZ.exists() else None,
        include_spectra=MGF_GZ.exists(),
        project_accession="PXD054720",
    )

    ds = Dataset(work_dir)
    yield ds
    ds.close()


# ---------------------------------------------------------------------------
# Group 1: Dataset structure tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestDatasetStructure:
    """Full dataset creation from mzIdentML + MGF."""

    def test_dataset_has_psm(self, dataset):
        """PSM structure should be present."""
        assert dataset.psm is not None
        assert dataset.psm.count() > 1000

    def test_dataset_has_pepmap(self, dataset):
        """Peptide-protein map should be present."""
        assert dataset.pepmap is not None
        assert dataset.pepmap.count() > 0

    def test_dataset_has_provenance(self, dataset):
        """Provenance should be present."""
        assert dataset.provenance is not None
        assert dataset.provenance.count() > 0

    def test_dataset_has_ontology(self, dataset):
        """Ontology should be present."""
        assert dataset.ontology is not None
        assert dataset.ontology.count() > 0

    def test_dataset_has_metadata(self, dataset):
        """Dataset metadata should be present."""
        assert dataset.dataset_meta is not None

    def test_available_structures(self, dataset):
        """Should have at least psm, pepmap, provenance, ontology, dataset."""
        available = set(dataset.available_structures)
        expected = {"psm", "pepmap", "provenance", "ontology", "dataset"}
        assert expected.issubset(available), f"Missing structures: {expected - available}"


# ---------------------------------------------------------------------------
# Group 2: PSM querying tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestPsmQuerying:
    """Query PSMs via Dataset API."""

    def test_filter_by_charge(self, dataset):
        """Filter PSMs by charge state."""
        charge2 = dataset.psm.filter("charge = 2").count()
        charge3 = dataset.psm.filter("charge = 3").count()
        total = dataset.psm.count()
        assert charge2 > 0
        assert charge3 > 0
        assert charge2 + charge3 <= total

    def test_select_columns(self, dataset):
        """Select specific PSM columns."""
        df = dataset.psm.select("sequence", "charge", "observed_mz").limit(10).to_df()
        assert set(df.columns) == {"sequence", "charge", "observed_mz"}
        assert len(df) == 10

    def test_unique_sequences_via_sql(self, dataset):
        """Count unique peptide sequences via SQL."""
        result = dataset.sql("SELECT COUNT(DISTINCT sequence) AS n FROM psm")
        n = result.to_df()["n"].iloc[0]
        assert n > 100

    def test_filter_decoys(self, dataset):
        """Dataset should have both target and decoy PSMs."""
        targets = dataset.psm.filter("is_decoy = false").count()
        decoys = dataset.psm.filter("is_decoy = true").count()
        total = dataset.psm.count()
        # Most PSMs should be targets
        assert targets > decoys or decoys == 0  # some datasets have no decoys
        assert targets + decoys == total

    def test_psm_to_dataframe(self, dataset):
        """Materialize PSMs as DataFrame."""
        df = dataset.psm.limit(50).to_df()
        assert len(df) == 50
        assert "sequence" in df.columns
        assert "peptidoform" in df.columns


# ---------------------------------------------------------------------------
# Group 3: Spectra attachment tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestSpectraAttachment:
    """Verify spectra from MGF are attached to PSMs."""

    def test_spectra_columns_present(self, dataset):
        """PSMs should have mz_array and intensity_array columns."""
        df = dataset.psm.limit(5).to_df()
        assert "mz_array" in df.columns
        assert "intensity_array" in df.columns

    def test_some_psms_have_spectra(self, dataset):
        """At least some PSMs should have non-null spectra arrays."""
        result = dataset.sql("SELECT COUNT(*) AS n FROM psm WHERE mz_array IS NOT NULL")
        n = result.to_df()["n"].iloc[0]
        assert n > 0, "No PSMs have spectra attached"

    def test_spectra_arrays_matching_length(self, dataset):
        """mz_array and intensity_array should have the same length per PSM."""
        result = dataset.sql("""
            SELECT len(mz_array) AS mz_len, len(intensity_array) AS int_len
            FROM psm
            WHERE mz_array IS NOT NULL
            LIMIT 100
        """)
        df = result.to_df()
        assert len(df) > 0, "No PSMs with non-null spectra for length check"
        assert (df["mz_len"] == df["int_len"]).all()


# ---------------------------------------------------------------------------
# Group 4: Pepmap querying tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestPepmapQuerying:
    """Query peptide-protein map via Dataset API."""

    def test_pepmap_has_sequences(self, dataset):
        """Pepmap entries should have non-empty sequences."""
        result = dataset.sql("SELECT COUNT(DISTINCT sequence) AS n FROM pepmap")
        n = result.to_df()["n"].iloc[0]
        assert n > 100

    def test_pepmap_has_proteins(self, dataset):
        """Pepmap entries should reference proteins."""
        result = dataset.sql("SELECT COUNT(DISTINCT p.accession) AS n FROM pepmap, UNNEST(pepmap.pg_accessions) AS _t(p)")
        n = result.to_df()["n"].iloc[0]
        assert n > 10

    def test_pepmap_unique_flag(self, dataset):
        """Some peptides should be unique (map to 1 protein), others shared."""
        result = dataset.sql("""
            SELECT is_unique, COUNT(*) AS n
            FROM pepmap
            GROUP BY is_unique
        """)
        df = result.to_df()
        # Should have at least one category
        assert len(df) >= 1
        assert df["n"].sum() > 0


# ---------------------------------------------------------------------------
# Group 5: Provenance tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestProvenanceQuerying:
    """Query provenance data."""

    def test_provenance_has_tool_name(self, dataset):
        """Provenance should record the search engine used."""
        df = dataset.provenance.to_df()
        assert "tool_name" in df.columns
        assert df["tool_name"].notna().any()
        # PXD054720 uses Mascot
        tool_names = df["tool_name"].tolist()
        assert any("mascot" in str(name).lower() for name in tool_names if name)

    def test_provenance_step_order(self, dataset):
        """Provenance steps should have positive step_order."""
        df = dataset.provenance.to_df()
        assert (df["step_order"] > 0).all()


# ---------------------------------------------------------------------------
# Group 6: Ontology tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestOntologyQuerying:
    """Query ontology entries for scores."""

    def test_ontology_has_score_entries(self, dataset):
        """Ontology should have entries for discovered scores."""
        df = dataset.ontology.to_df()
        assert len(df) > 0
        assert "field_name" in df.columns
        assert "ontology_name" in df.columns

    def test_ontology_describes_the_written_views(self, dataset):
        """Ontology entries cover the views this mzIdentML conversion wrote."""
        df = dataset.ontology.to_df()
        assert {"psm", "pg"} <= set(df["view"])


# ---------------------------------------------------------------------------
# Group 7: Dataset metadata tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestDatasetMetadata:
    """Validate dataset metadata."""

    def test_project_accession(self, dataset):
        """Dataset should record the project accession."""
        df = dataset.dataset_meta.to_df()
        assert df["project_accession"].iloc[0] == "PXD054720"

    def test_creation_date(self, dataset):
        """Dataset should have a creation date."""
        df = dataset.dataset_meta.to_df()
        assert df["creation_date"].iloc[0] is not None


# ---------------------------------------------------------------------------
# Group 8: Modification view tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestModificationView:
    """Test ModificationView querying across PSMs."""

    def test_modification_frequency(self, dataset):
        """Should find modifications in PSMs."""
        mod_result = dataset.modification_view.frequency()
        df = mod_result.to_df()
        # PSMs from Mascot should have some modifications
        if len(df) > 0:
            assert "modification_name" in df.columns
            assert "n_psms" in df.columns
            assert df["n_psms"].iloc[0] > 0


# ---------------------------------------------------------------------------
# Group 9: Cross-structure SQL tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestCrossStructureQueries:
    """Test SQL queries that join multiple structures."""

    def test_psm_protein_count(self, dataset):
        """Count proteins via PSM protein_accessions."""
        result = dataset.sql("""
            SELECT COUNT(DISTINCT prot) AS n_proteins
            FROM psm, UNNEST(psm.protein_accessions) AS _t(prot)
        """)
        n = result.to_df()["n_proteins"].iloc[0]
        assert n > 10

    def test_psm_pepmap_join(self, dataset):
        """Join PSMs with pepmap to find shared peptides."""
        result = dataset.sql("""
            SELECT COUNT(*) AS n
            FROM psm
            JOIN pepmap ON psm.sequence = pepmap.sequence
        """)
        n = result.to_df()["n"].iloc[0]
        assert n > 0
