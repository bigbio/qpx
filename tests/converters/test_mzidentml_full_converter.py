"""Tests for MzIdentMLConverter — full mzIdentML-to-Dataset pipeline.

Tests the high-level converter that produces a complete QPX dataset
directory from an mzIdentML file (PSM, pepmap, provenance, ontology,
dataset metadata).
"""

import pytest

lxml = pytest.importorskip("lxml", reason="lxml required for mzIdentML tests")
from pathlib import Path

import pyarrow.parquet as pq

from qpx.converters.mzidentml import MzIdentMLConverter
from qpx.core.data import (
    DatasetSchema,
    PepMapSchema,
    ProvenanceSchema,
    PsmSchema,
)

XL_EXAMPLES = Path(__file__).parent.parent / "examples" / "mzidentml" / "xl"
PXD054720_DIR = Path(__file__).parent.parent / "examples" / "mzidentml" / "PXD054720"


# ---------------------------------------------------------------------------
# Basic conversion (XL test files — small, always available)
# ---------------------------------------------------------------------------


class TestMzIdentMLConverterBasic:
    """Test full converter with small XL-MS mzIdentML example files."""

    @pytest.fixture
    def converted_dataset(self, tmp_path):
        """Convert the EDC crosslink file and return the output directory."""
        mzid = XL_EXAMPLES / "Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz"
        converter = MzIdentMLConverter()
        return converter.convert(
            mzid_path=str(mzid),
            output_folder=str(tmp_path / "output"),
            output_prefix="edc_test",
            project_accession="PXD999999",
        )

    def test_returns_output_folder(self, converted_dataset):
        assert converted_dataset.is_dir()

    def test_psm_parquet_exists(self, converted_dataset):
        psm = converted_dataset / "edc_test.psm.parquet"
        assert psm.exists(), f"PSM parquet not found at {psm}"

    def test_psm_schema_valid(self, converted_dataset):
        psm = converted_dataset / "edc_test.psm.parquet"
        table = pq.read_table(str(psm))
        errors = PsmSchema.validate(table)
        assert not errors, f"PSM schema validation errors: {errors}"

    def test_psm_has_rows(self, converted_dataset):
        psm = converted_dataset / "edc_test.psm.parquet"
        table = pq.read_table(str(psm))
        assert table.num_rows > 0

    def test_pepmap_parquet_exists(self, converted_dataset):
        pepmap = converted_dataset / "edc_test.pepmap.parquet"
        assert pepmap.exists(), f"Pepmap parquet not found at {pepmap}"

    def test_pepmap_schema_valid(self, converted_dataset):
        pepmap = converted_dataset / "edc_test.pepmap.parquet"
        table = pq.read_table(str(pepmap))
        errors = PepMapSchema.validate(table)
        assert not errors, f"PepMap schema validation errors: {errors}"

    def test_pepmap_has_rows(self, converted_dataset):
        pepmap = converted_dataset / "edc_test.pepmap.parquet"
        table = pq.read_table(str(pepmap))
        assert table.num_rows > 0

    def test_pepmap_is_unique_computed(self, converted_dataset):
        pepmap = converted_dataset / "edc_test.pepmap.parquet"
        table = pq.read_table(str(pepmap))
        is_unique_col = table.column("is_unique").to_pylist()
        # All values should be bool (True or False), not None
        assert all(v is not None for v in is_unique_col)

    def test_pepmap_deduplicated(self, converted_dataset):
        pepmap = converted_dataset / "edc_test.pepmap.parquet"
        table = pq.read_table(str(pepmap))
        rows = table.to_pydict()
        keys = set(rows["peptidoform"])
        assert len(keys) == table.num_rows, "Pepmap should be deduplicated"

    def test_provenance_parquet_exists(self, converted_dataset):
        prov = converted_dataset / "edc_test.provenance.parquet"
        assert prov.exists(), f"Provenance parquet not found at {prov}"

    def test_provenance_schema_valid(self, converted_dataset):
        prov = converted_dataset / "edc_test.provenance.parquet"
        table = pq.read_table(str(prov))
        errors = ProvenanceSchema.validate(table)
        assert not errors, f"Provenance schema validation errors: {errors}"

    def test_provenance_has_mascot(self, converted_dataset):
        prov = converted_dataset / "edc_test.provenance.parquet"
        table = pq.read_table(str(prov))
        tool_names = table.column("tool_name").to_pylist()
        assert "Mascot" in tool_names, f"Expected 'Mascot' in provenance, got {tool_names}"

    def test_provenance_has_version(self, converted_dataset):
        prov = converted_dataset / "edc_test.provenance.parquet"
        table = pq.read_table(str(prov))
        versions = table.column("tool_version").to_pylist()
        assert any(v is not None and len(v) > 0 for v in versions)

    def test_dataset_parquet_exists(self, converted_dataset):
        ds = converted_dataset / "edc_test.dataset.parquet"
        assert ds.exists(), f"Dataset parquet not found at {ds}"

    def test_dataset_schema_valid(self, converted_dataset):
        ds = converted_dataset / "edc_test.dataset.parquet"
        table = pq.read_table(str(ds))
        errors = DatasetSchema.validate(table)
        assert not errors, f"Dataset schema validation errors: {errors}"

    def test_dataset_has_accession(self, converted_dataset):
        ds = converted_dataset / "edc_test.dataset.parquet"
        table = pq.read_table(str(ds))
        accession = table.column("project_accession")[0].as_py()
        assert accession == "PXD999999"

    def test_dataset_has_creation_date(self, converted_dataset):
        ds = converted_dataset / "edc_test.dataset.parquet"
        table = pq.read_table(str(ds))
        date = table.column("creation_date")[0].as_py()
        assert date is not None and len(date) > 0

    def test_dataset_software_from_provenance(self, converted_dataset):
        ds = converted_dataset / "edc_test.dataset.parquet"
        table = pq.read_table(str(ds))
        sw = table.column("software_name")[0].as_py()
        assert sw == "Mascot"


# ---------------------------------------------------------------------------
# Default accession handling
# ---------------------------------------------------------------------------


class TestDefaultProjectAccession:
    def test_unknown_when_no_accession(self, tmp_path):
        mzid = XL_EXAMPLES / "scores_and_thresholds_1_3_0_draft.mzid.gz"
        converter = MzIdentMLConverter()
        out = converter.convert(
            mzid_path=str(mzid),
            output_folder=str(tmp_path / "out"),
            output_prefix="test",
        )
        ds = pq.read_table(str(out / "test.dataset.parquet"))
        assert ds.column("project_accession")[0].as_py() == "unknown"


# ---------------------------------------------------------------------------
# Conversion with different mzIdentML files
# ---------------------------------------------------------------------------


class TestMultipleSpectraFile:
    def test_convert_multiple_spectra(self, tmp_path):
        mzid = XL_EXAMPLES / "multiple_spectra_per_id_1_3_0_draft.mzid.gz"
        converter = MzIdentMLConverter()
        out = converter.convert(
            mzid_path=str(mzid),
            output_folder=str(tmp_path / "multi"),
            output_prefix="multi",
        )
        psm = pq.read_table(str(out / "multi.psm.parquet"))
        assert psm.num_rows > 0
        pepmap = pq.read_table(str(out / "multi.pepmap.parquet"))
        assert pepmap.num_rows > 0


# ---------------------------------------------------------------------------
# Large-data tests (PXD054720 — gzip compressed, with MGF)
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestPXD054720FullConversion:
    """Full-pipeline conversion of real PRIDE data with spectra attachment."""

    @pytest.fixture(scope="class")
    def converted_dataset(self, tmp_path_factory):
        if not (PXD054720_DIR / "F001234.mzid.gz").exists():
            pytest.skip("PXD054720 test data not found")
        work = tmp_path_factory.mktemp("pxd054720_full")
        converter = MzIdentMLConverter()
        return converter.convert(
            mzid_path=str(PXD054720_DIR / "F001234.mzid.gz"),
            output_folder=str(work),
            output_prefix="PXD054720",
            project_accession="PXD054720",
        )

    def test_psm_parquet(self, converted_dataset):
        table = pq.read_table(str(converted_dataset / "PXD054720.psm.parquet"))
        assert table.num_rows > 1000

    def test_pepmap_parquet(self, converted_dataset):
        table = pq.read_table(str(converted_dataset / "PXD054720.pepmap.parquet"))
        assert table.num_rows > 100

    def test_provenance_parquet(self, converted_dataset):
        table = pq.read_table(str(converted_dataset / "PXD054720.provenance.parquet"))
        assert table.num_rows >= 1

    def test_dataset_parquet(self, converted_dataset):
        table = pq.read_table(str(converted_dataset / "PXD054720.dataset.parquet"))
        assert table.column("project_accession")[0].as_py() == "PXD054720"


@pytest.mark.large_data
class TestPXD054720WithSpectra:
    """Test spectra attachment from MGF file."""

    @pytest.fixture(scope="class")
    def converted_dataset(self, tmp_path_factory):
        mzid = PXD054720_DIR / "F001234.mzid.gz"
        mgf = PXD054720_DIR / "F001234.mgf.gz"
        if not mzid.exists() or not mgf.exists():
            pytest.skip("PXD054720 test data (mzid + mgf) not found")
        work = tmp_path_factory.mktemp("pxd054720_spectra")
        converter = MzIdentMLConverter()
        return converter.convert(
            mzid_path=str(mzid),
            output_folder=str(work),
            output_prefix="PXD054720_sp",
            mgf_path=str(mgf),
            include_spectra=True,
            project_accession="PXD054720",
        )

    def test_some_psms_have_spectra(self, converted_dataset):
        table = pq.read_table(str(converted_dataset / "PXD054720_sp.psm.parquet"))
        mz_col = table.column("mz_array").to_pylist()
        with_spectra = sum(1 for v in mz_col if v is not None and len(v) > 0)
        assert with_spectra > 0, "Expected some PSMs with attached spectra"

    def test_spectra_arrays_match_length(self, converted_dataset):
        table = pq.read_table(str(converted_dataset / "PXD054720_sp.psm.parquet"))
        mz_col = table.column("mz_array").to_pylist()
        int_col = table.column("intensity_array").to_pylist()
        for mz, inten in zip(mz_col, int_col):
            if mz is not None:
                assert inten is not None
                assert len(mz) == len(inten)


# ---------------------------------------------------------------------------
# Import test
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Index fallback tests
# ---------------------------------------------------------------------------


class TestIndexFallback:
    """Spectra attachment must use the nativeID scheme to pick the right lookup."""

    def test_index_scheme_uses_position_lookup(self, tmp_path):
        """``index=`` records attach by 0-based file position, not by SCANS=."""
        # MGF where SCANS= deliberately does NOT match the index values, so a
        # scan-number lookup would either miss or (worse) match the wrong
        # spectrum. Only a position-based lookup gives the correct peaks.
        mgf = tmp_path / "test.mgf"
        mgf.write_text(
            "BEGIN IONS\nTITLE=first\nSCANS=999\n100.0 500\nEND IONS\nBEGIN IONS\nTITLE=second\nSCANS=998\n200.0 600\nEND IONS\n"
        )

        records = [
            {"scan": [0], "_spectrum_id_scheme": "index", "mz_array": None, "intensity_array": None, "rt": None},
            {"scan": [1], "_spectrum_id_scheme": "index", "mz_array": None, "intensity_array": None, "rt": None},
        ]

        result = MzIdentMLConverter._attach_spectra(records, str(mgf))
        assert result[0]["mz_array"] == [100.0]  # position 0
        assert result[1]["mz_array"] == [200.0]  # position 1

    def test_index_value_not_attached_by_scan(self, tmp_path):
        """An ``index=`` value must never be attached via a colliding SCANS=.

        Regression for bigbio/qpx#249: index 999 collides with SCANS=999 of an
        UNRELATED spectrum; the old scan-first lookup attached those wrong peaks.
        """
        # Spectrum at position 0 carries SCANS=999; the correct spectrum for
        # index=999 does not exist in this tiny MGF, so nothing must attach.
        mgf = tmp_path / "collide.mgf"
        mgf.write_text("BEGIN IONS\nTITLE=first\nSCANS=999\n100.0 500\nEND IONS\n")

        records = [
            {"scan": [999], "_spectrum_id_scheme": "index", "mz_array": None, "intensity_array": None, "rt": None},
        ]

        result = MzIdentMLConverter._attach_spectra(records, str(mgf))
        # No spectrum at position 999 → nothing attached (NOT the SCANS=999 peaks).
        assert result[0]["mz_array"] is None

    def test_index_file_attaches_correct_spectrum_end_to_end(self, tmp_path):
        """A PSM whose ``spectrumID`` is ``index=N`` attaches the peaks of the
        Nth spectrum, even when an unrelated spectrum carries ``SCANS=N``.

        End-to-end: the scheme is derived from the ``index=`` spectrumID by the
        adapter, then honoured by ``_attach_spectra``. Regression for
        bigbio/qpx#249 (the old scan-first lookup attached the SCANS-colliding
        spectrum's peaks instead).
        """
        from qpx.converters.mzidentml.psm_adapter import MzIdentMLPsmAdapter

        # Position 0 is a decoy whose SCANS collides with the target index (1);
        # position 1 is the CORRECT spectrum for ``index=1``.
        mgf = tmp_path / "e2e.mgf"
        mgf.write_text(
            "BEGIN IONS\nTITLE=decoy\nSCANS=1\n100.0 500\nEND IONS\nBEGIN IONS\nTITLE=target\nSCANS=42\n200.0 600\nEND IONS\n"
        )

        parsed = {
            "spectra_data": {"SD1": "e2e.mgf"},
            "db_sequences": {},
            "peptides": {"Pep1": {"sequence": "PEPTIDE", "modifications": None, "xl_mods": {}}},
            "peptide_evidence": {},
            "linker_info": {},
            "spectrum_results": [
                {
                    "spectrum_id": "index=1",  # the 2nd spectrum (position 1)
                    "spectra_data_ref": "SD1",
                    "rt": None,
                    "siis": [
                        {
                            "id": "SII_1",
                            "charge": 2,
                            "exp_mz": 500.0,
                            "calc_mz": 500.0,
                            "rank": 1,
                            "pass_threshold": True,
                            "peptide_ref": "Pep1",
                            "cv_params": [],
                            "peptide_evidence_refs": [],
                        }
                    ],
                }
            ],
        }

        with MzIdentMLPsmAdapter() as adapter:
            records = adapter._build_psm_records(parsed)
        assert records[0]["_spectrum_id_scheme"] == "index"

        result = MzIdentMLConverter._attach_spectra(records, str(mgf))
        # Correct spectrum is position 1 (200.0), NOT the SCANS=1 collision (100.0).
        assert result[0]["mz_array"] == [200.0]
        assert result[0]["intensity_array"] == [600.0]


# ---------------------------------------------------------------------------
# Import test
# ---------------------------------------------------------------------------


class TestImport:
    def test_import_from_package(self):
        from qpx.converters.mzidentml import MzIdentMLConverter as Cls

        assert Cls is not None

    def test_import_from_module(self):
        from qpx.converters.mzidentml.converter import MzIdentMLConverter as Cls

        assert Cls is not None
