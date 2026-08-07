"""Tests for the mzIdentML PSM adapter (including XL-MS cross-linking)."""

import pytest

lxml = pytest.importorskip("lxml", reason="lxml required for mzIdentML tests")

from pathlib import Path

import pyarrow.parquet as pq

from qpx.converters.mzidentml.psm_adapter import (
    MzIdentMLPsmAdapter,
    _group_siis,
    _parse_native_id,
    _parse_scan,
)

XL_EXAMPLES = Path(__file__).parent.parent / "examples" / "mzidentml" / "xl"


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _convert(mzid_name: str, tmp_path: Path, examples_dir: Path = None) -> pq.ParquetFile:
    """Convert a test mzid file and return the output as a ParquetFile."""
    base_dir = examples_dir or XL_EXAMPLES
    mzid_path = base_dir / mzid_name
    out_path = tmp_path / f"{mzid_name}.psm.parquet"
    with MzIdentMLPsmAdapter() as adapter:
        adapter.convert(mzid_path=str(mzid_path), output_path=str(out_path))
    return pq.read_table(str(out_path))


# ---------------------------------------------------------------------------
# Unit tests for helpers
# ---------------------------------------------------------------------------


class TestParseScan:
    def test_scan_format(self):
        assert _parse_scan("scan=100") == [100]

    def test_index_format(self):
        assert _parse_scan("index=7483") == [7483]

    def test_spectrum_format(self):
        assert _parse_scan("spectrum=5") == [5]

    def test_bare_integer(self):
        assert _parse_scan("42") == [42]

    def test_unknown_format(self):
        # Unparseable spectrumIDs must NOT collapse onto a ``[0]`` sentinel
        # (which aliases distinct spectra and a real scan-0 spectrum).
        assert _parse_scan("something_else") == []


class TestParseNativeId:
    """The nativeID scheme distinguishes scan numbers from file indices (#249)."""

    def test_scan_is_scan_scheme(self):
        assert _parse_native_id("scan=100") == ([100], "scan")

    def test_index_is_index_scheme(self):
        # The regression: ``index=`` is a file position, not a scan number.
        assert _parse_native_id("index=1999") == ([1999], "index")

    def test_spectrum_is_scan_scheme(self):
        assert _parse_native_id("spectrum=5") == ([5], "scan")

    def test_bare_integer_is_scan_scheme(self):
        assert _parse_native_id("42") == ([42], "scan")

    def test_mzml_nativeid_scan(self):
        assert _parse_native_id("controllerType=0 controllerNumber=1 scan=7") == ([7], "scan")

    def test_unparseable_is_empty_none(self):
        # No ``[0]`` sentinel: distinct unparseable spectra must not alias.
        assert _parse_native_id("something_else") == ([], None)


class TestGroupSiis:
    def _make_sii(self, sii_id, cv_params=None):
        return {
            "id": sii_id,
            "charge": 2,
            "exp_mz": 500.0,
            "calc_mz": 500.0,
            "rank": 1,
            "pass_threshold": True,
            "peptide_ref": f"Pep_{sii_id}",
            "cv_params": cv_params or [],
            "peptide_evidence_refs": [],
        }

    def test_regular_sii(self):
        siis = [self._make_sii("SII_1")]
        groups = _group_siis(siis)
        assert len(groups) == 1
        assert groups[0]["xl_type"] == "none"

    def test_crosslink_pair(self):
        siis = [
            self._make_sii("SII_1", [{"accession": "MS:1002511", "name": "xl", "value": "1"}]),
            self._make_sii("SII_2", [{"accession": "MS:1002511", "name": "xl", "value": "1"}]),
        ]
        groups = _group_siis(siis)
        assert len(groups) == 1
        assert groups[0]["xl_type"] == "inter"
        assert groups[0]["alpha"]["id"] == "SII_1"
        assert groups[0]["beta"]["id"] == "SII_2"

    def test_looplink(self):
        siis = [
            self._make_sii("SII_1", [{"accession": "MS:1003329", "name": "looplink", "value": ""}]),
        ]
        groups = _group_siis(siis)
        assert len(groups) == 1
        assert groups[0]["xl_type"] == "intra"

    def test_noncovalent_pair(self):
        siis = [
            self._make_sii("SII_1", [{"accession": "MS:1003331", "name": "noncov", "value": "1"}]),
            self._make_sii("SII_2", [{"accession": "MS:1003331", "name": "noncov", "value": "1"}]),
        ]
        groups = _group_siis(siis)
        assert len(groups) == 1
        assert groups[0]["xl_type"] == "noncovalent"


# ---------------------------------------------------------------------------
# Integration tests — Xlink EDC (inter-peptide + looplinks)
# ---------------------------------------------------------------------------


class TestXlinkEDC:
    def test_convert_produces_psms(self, tmp_path):
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        assert table.num_rows > 0
        psm_ids = table.column("psm_id").to_pylist()
        assert len(psm_ids) == len(set(psm_ids))

    def test_has_crosslinked_psms(self, tmp_path):
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        xl_count = sum(1 for v in table.column("cross_links") if v.as_py() is not None)
        assert xl_count > 0, "Expected some PSMs with cross-links"

    def test_has_looplinks(self, tmp_path):
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        intra = []
        for i in range(table.num_rows):
            xl = table.column("cross_links")[i].as_py()
            if xl and xl[0]["xl_type"] == "intra":
                intra.append(xl[0])
        assert len(intra) > 0, "Expected looplink (intra) cross-links"
        # Looplinks should have donor_position and acceptor_position but no partner
        for link in intra:
            assert link["partner_sequence"] is None
            assert link["donor_position"] > 0
            assert link["acceptor_position"] > 0

    def test_has_interlinks(self, tmp_path):
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        inter = []
        for i in range(table.num_rows):
            xl = table.column("cross_links")[i].as_py()
            if xl and xl[0]["xl_type"] == "inter":
                inter.append(xl[0])
        assert len(inter) > 0, "Expected inter-peptide cross-links"
        for link in inter:
            assert link["partner_sequence"] is not None
            assert len(link["partner_sequence"]) > 0

    def test_linker_info(self, tmp_path):
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        for i in range(table.num_rows):
            xl = table.column("cross_links")[i].as_py()
            if xl and xl[0]["xl_type"] in ("intra", "inter"):
                assert xl[0]["linker_name"] == "EDC"
                assert xl[0]["linker_mass"] != 0
                break

    def test_regular_psms_have_no_crosslinks(self, tmp_path):
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        # Some PSMs should be regular (no cross-links)
        regular = sum(1 for v in table.column("cross_links") if v.as_py() is None)
        assert regular > 0, "Expected some regular (non-XL) PSMs"


# ---------------------------------------------------------------------------
# Integration tests — noncovalently associated peptides
# ---------------------------------------------------------------------------


class TestNoncovalent:
    def test_convert_produces_psms(self, tmp_path):
        table = _convert("noncovalently_assoc_1_3_0_draft.mzid.gz", tmp_path)
        assert table.num_rows >= 1

    def test_noncovalent_crosslink(self, tmp_path):
        table = _convert("noncovalently_assoc_1_3_0_draft.mzid.gz", tmp_path)
        xl = table.column("cross_links")[0].as_py()
        assert xl is not None
        assert xl[0]["xl_type"] == "noncovalent"
        assert xl[0]["partner_sequence"] is not None
        assert len(xl[0]["partner_sequence"]) > 0


# ---------------------------------------------------------------------------
# Integration tests — multiple spectra per identification
# ---------------------------------------------------------------------------


class TestMultipleSpectra:
    def test_convert_produces_psms(self, tmp_path):
        table = _convert("multiple_spectra_per_id_1_3_0_draft.mzid.gz", tmp_path)
        assert table.num_rows > 0

    def test_has_crosslinks(self, tmp_path):
        table = _convert("multiple_spectra_per_id_1_3_0_draft.mzid.gz", tmp_path)
        xl_count = sum(1 for v in table.column("cross_links") if v.as_py() is not None)
        assert xl_count > 0

    def test_inter_peptide_link(self, tmp_path):
        table = _convert("multiple_spectra_per_id_1_3_0_draft.mzid.gz", tmp_path)
        for i in range(table.num_rows):
            xl = table.column("cross_links")[i].as_py()
            if xl and xl[0]["xl_type"] == "inter":
                assert xl[0]["linker_name"] == "DSSO"
                assert xl[0]["partner_sequence"] is not None
                break
        else:
            pytest.fail("No inter-peptide cross-link found")


# ---------------------------------------------------------------------------
# Integration tests — scores and thresholds
# ---------------------------------------------------------------------------


class TestScoresAndThresholds:
    def test_convert_produces_psms(self, tmp_path):
        table = _convert("scores_and_thresholds_1_3_0_draft.mzid.gz", tmp_path)
        assert table.num_rows == 2

    def test_has_crosslinks(self, tmp_path):
        table = _convert("scores_and_thresholds_1_3_0_draft.mzid.gz", tmp_path)
        xl_count = sum(1 for v in table.column("cross_links") if v.as_py() is not None)
        assert xl_count == 2

    def test_has_fdr_scores(self, tmp_path):
        table = _convert("scores_and_thresholds_1_3_0_draft.mzid.gz", tmp_path)
        for i in range(table.num_rows):
            scores = table.column("additional_scores")[i].as_py()
            score_names = {s["score_name"] for s in scores}
            assert "crosslinked_psm_level_global_fdr" in score_names
            assert "peptide_pair_sequence_level_global_fdr" in score_names


# ---------------------------------------------------------------------------
# Schema validation
# ---------------------------------------------------------------------------


class TestPepExtraction:
    """Test posterior_error_probability extraction."""

    def test_extract_pep_from_cvparams(self):
        """PEP should be extracted from MS:1001493 cvParam."""
        from qpx.converters.mzidentml.psm_adapter import _extract_pep

        cv_params = [
            {"accession": "MS:1001171", "name": "Mascot:score", "value": "45.2"},
            {
                "accession": "MS:1001493",
                "name": "posterior error probability",
                "value": "0.001",
            },
        ]
        assert _extract_pep(cv_params) == pytest.approx(0.001)

    def test_extract_pep_global(self):
        """Should extract PEP from MS:1002352 (global PEP)."""
        from qpx.converters.mzidentml.psm_adapter import _extract_pep

        cv_params = [
            {
                "accession": "MS:1002352",
                "name": "PSM-level global FDR",
                "value": "0.01",
            },
        ]
        pep = _extract_pep(cv_params)
        assert pep == pytest.approx(0.01)

    def test_extract_pep_missing(self):
        """Should return None when PEP is not present."""
        from qpx.converters.mzidentml.psm_adapter import _extract_pep

        cv_params = [
            {"accession": "MS:1001171", "name": "Mascot:score", "value": "45.2"},
        ]
        assert _extract_pep(cv_params) is None

    def test_extract_pep_invalid_value(self):
        """Should return None for non-numeric PEP values."""
        from qpx.converters.mzidentml.psm_adapter import _extract_pep

        cv_params = [
            {
                "accession": "MS:1001493",
                "name": "posterior error probability",
                "value": "",
            },
        ]
        assert _extract_pep(cv_params) is None

    def test_pep_also_in_additional_scores(self):
        """PEP should appear in both dedicated field and additional_scores."""
        from qpx.converters.mzidentml.psm_adapter import _extract_pep, _extract_scores

        cv_params = [
            {
                "accession": "MS:1001493",
                "name": "posterior error probability",
                "value": "0.05",
            },
            {"accession": "MS:1001171", "name": "Mascot:score", "value": "30.0"},
        ]
        pep = _extract_pep(cv_params)
        scores = _extract_scores(cv_params)

        assert pep == pytest.approx(0.05)
        # PEP should also be in additional_scores
        pep_scores = [s for s in scores if "posterior" in s["score_name"]]
        assert len(pep_scores) == 1
        assert pep_scores[0]["score_value"] == pytest.approx(0.05)


class TestRtExtraction:
    """Test retention time extraction."""

    def test_rt_from_sir_cvparam_seconds(self):
        """RT in seconds should be extracted directly."""
        from qpx.core.cv_terms import CV_SCAN_START_TIME

        # Just test that the constant exists and is correct
        assert CV_SCAN_START_TIME == "MS:1000016"

    def test_rt_from_mgf(self):
        """MGF RTINSECONDS should be parsed and stored."""
        import os
        import tempfile

        from qpx.converters.mzidentml.mgf_parser import MgfSpectraIndex

        mgf_content = (
            "BEGIN IONS\n"
            "TITLE=test_spectrum\n"
            "SCANS=100\n"
            "RTINSECONDS=123.45\n"
            "PEPMASS=500.0\n"
            "CHARGE=2+\n"
            "200.0 1000\n"
            "300.0 2000\n"
            "END IONS\n"
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".mgf", delete=False) as f:
            f.write(mgf_content)
            f.flush()
            tmp_name = f.name
        try:
            idx = MgfSpectraIndex(tmp_name)
            spec = idx.get_spectrum(100)
            assert spec is not None
            assert spec["rt"] == pytest.approx(123.45)
        finally:
            os.unlink(tmp_name)

    def test_rt_none_when_absent(self):
        """RT should be None when RTINSECONDS is not in MGF."""
        import os
        import tempfile

        from qpx.converters.mzidentml.mgf_parser import MgfSpectraIndex

        mgf_content = "BEGIN IONS\nTITLE=test_spectrum\nSCANS=200\nPEPMASS=500.0\n200.0 1000\nEND IONS\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".mgf", delete=False) as f:
            f.write(mgf_content)
            f.flush()
            tmp_name = f.name
        try:
            idx = MgfSpectraIndex(tmp_name)
            spec = idx.get_spectrum(200)
            assert spec is not None
            assert spec.get("rt") is None
        finally:
            os.unlink(tmp_name)


class TestSchemaValidation:
    def test_output_matches_psm_schema(self, tmp_path):
        """Verify the output parquet has the correct PSM schema structure."""
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        from qpx.core.data import PsmSchema

        errors = PsmSchema.validate(table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# P0-2: MzIdentMLPgAdapter — protein groups from ProteinDetectionList
# ---------------------------------------------------------------------------


class TestMzIdentMLPgAdapter:
    """MzIdentMLPgAdapter converts ProteinAmbiguityGroups to pg.parquet."""

    def test_convert_produces_pg_parquet(self, tmp_path):
        """adapter.convert() should produce a pg.parquet file."""
        from qpx.converters.mzidentml.pg_adapter import MzIdentMLPgAdapter

        output = tmp_path / "test.pg.parquet"
        with MzIdentMLPgAdapter() as adapter:
            adapter.convert(mzid_path=str(PXD054720_MZID), output_path=str(output))
        assert output.exists(), "pg.parquet was not created"

    def test_pg_row_count_matches_pag_count(self, tmp_path):
        """Row count should equal number of PAGs (123 in PXD054720)."""
        from qpx.converters.mzidentml.pg_adapter import MzIdentMLPgAdapter

        output = tmp_path / "test.pg.parquet"
        with MzIdentMLPgAdapter() as adapter:
            adapter.convert(mzid_path=str(PXD054720_MZID), output_path=str(output))
        table = pq.read_table(str(output))
        assert table.num_rows == 123, f"Expected 123 PAGs, got {table.num_rows}"

    def test_pg_anchor_protein_set(self, tmp_path):
        """anchor_protein should be the group representative (MS:1002403)."""
        from qpx.converters.mzidentml.pg_adapter import MzIdentMLPgAdapter

        output = tmp_path / "test.pg.parquet"
        with MzIdentMLPgAdapter() as adapter:
            adapter.convert(mzid_path=str(PXD054720_MZID), output_path=str(output))
        table = pq.read_table(str(output))
        df = table.to_pandas()
        anchor_row = df[df["anchor_protein"] == "sp|P0A6Y8|DNAK_ECOLI"]
        assert len(anchor_row) > 0, (
            f"Expected sp|P0A6Y8|DNAK_ECOLI as anchor protein in first PAG. Got anchors: {df['anchor_protein'].tolist()[:5]}"
        )

    def test_pg_accessions_not_empty(self, tmp_path):
        """pg_accessions should be non-empty lists for all rows."""
        from qpx.converters.mzidentml.pg_adapter import MzIdentMLPgAdapter

        output = tmp_path / "test.pg.parquet"
        with MzIdentMLPgAdapter() as adapter:
            adapter.convert(mzid_path=str(PXD054720_MZID), output_path=str(output))
        table = pq.read_table(str(output))
        for row in table.to_pylist():
            assert row["pg_accessions"] is not None, "pg_accessions is None"
            assert len(row["pg_accessions"]) > 0, "pg_accessions is empty"

    def test_pg_peptides_not_empty(self, tmp_path):
        """peptides list should be non-empty for all rows."""
        from qpx.converters.mzidentml.pg_adapter import MzIdentMLPgAdapter

        output = tmp_path / "test.pg.parquet"
        with MzIdentMLPgAdapter() as adapter:
            adapter.convert(mzid_path=str(PXD054720_MZID), output_path=str(output))
        table = pq.read_table(str(output))
        for row in table.to_pylist():
            assert row["peptides"] is not None, "peptides is None"
            assert len(row["peptides"]) > 0, "peptides list is empty"

    def test_pg_run_file_name_set(self, tmp_path):
        """grouped_runs should be a single-element list derived from the mzid filename."""
        from qpx.converters.mzidentml.pg_adapter import MzIdentMLPgAdapter

        output = tmp_path / "test.pg.parquet"
        with MzIdentMLPgAdapter() as adapter:
            adapter.convert(mzid_path=str(PXD054720_MZID), output_path=str(output))
        table = pq.read_table(str(output))
        grouped_runs = table.column("grouped_runs").to_pylist()
        for gr in grouped_runs:
            assert isinstance(gr, list) and len(gr) == 1, f"Expected single-element list, got: {gr}"
        run_names = {gr[0] for gr in grouped_runs}
        assert len(run_names) == 1, f"Expected 1 run, got: {run_names}"
        assert "F001234" in list(run_names)[0], f"Expected F001234 in run name, got: {run_names}"


# ---------------------------------------------------------------------------
# P0-4: Provenance parameters from AnalysisProtocolCollection
# ---------------------------------------------------------------------------

PXD054720_MZID = Path(__file__).parent.parent / "examples" / "mzidentml" / "PXD054720" / "F001234.mzid.gz"


class TestMzIdentMLProvenanceParams:
    """_build_provenance must extract enzyme/modification params from AnalysisProtocolCollection."""

    def test_provenance_parameters_not_none(self):
        """_build_provenance should return parameters list, not None."""
        from qpx.converters.mzidentml.converter import MzIdentMLConverter

        records = MzIdentMLConverter._build_provenance(PXD054720_MZID)
        assert records, "Expected at least one provenance record"
        params = records[0].get("parameters")
        assert params is not None, "parameters should not be None — parse AnalysisProtocolCollection"
        assert len(params) > 0, "parameters list should not be empty"

    def test_provenance_contains_enzyme(self):
        """parameters should include an enzyme entry (key='enzyme')."""
        from qpx.converters.mzidentml.converter import MzIdentMLConverter

        records = MzIdentMLConverter._build_provenance(PXD054720_MZID)
        params = records[0]["parameters"]
        enzyme_entries = [p for p in params if p["key"] == "enzyme"]
        assert enzyme_entries, f"No 'enzyme' key in parameters: {params}"
        assert any("Trypsin" in p["value"] for p in enzyme_entries), f"Expected Trypsin in enzyme entries, got: {enzyme_entries}"

    def test_provenance_contains_fixed_modification(self):
        """parameters should include fixed modification (Carbamidomethyl on C)."""
        from qpx.converters.mzidentml.converter import MzIdentMLConverter

        records = MzIdentMLConverter._build_provenance(PXD054720_MZID)
        params = records[0]["parameters"]
        fixed_mods = [p for p in params if p["key"] == "fixed_mod"]
        assert fixed_mods, f"No 'fixed_mod' entries in parameters: {params}"
        assert any("Carbamidomethyl" in p["value"] for p in fixed_mods), (
            f"Expected Carbamidomethyl in fixed_mod entries: {fixed_mods}"
        )

    def test_provenance_contains_variable_modification(self):
        """parameters should include variable modification (Oxidation on M)."""
        from qpx.converters.mzidentml.converter import MzIdentMLConverter

        records = MzIdentMLConverter._build_provenance(PXD054720_MZID)
        params = records[0]["parameters"]
        var_mods = [p for p in params if p["key"] == "variable_mod"]
        assert var_mods, f"No 'variable_mod' entries in parameters: {params}"
        assert any("Oxidation" in p["value"] for p in var_mods), f"Expected Oxidation in variable_mod entries: {var_mods}"


# ---------------------------------------------------------------------------
# Shared-modifications mutation (bigbio/qpx#249)
# ---------------------------------------------------------------------------


class TestSharedModificationsNoCrossContamination:
    """PSMs that share a ``peptide_ref`` must not accumulate each other's
    site-localization scores.

    The Peptide's ``modifications`` list is shared by every SII that references
    it, and site-localization scoring mutates modification positions in place.
    Without a per-PSM copy, the second PSM inherits the first PSM's phosphoRS
    score, corrupting localization data.
    """

    @staticmethod
    def _phospho_mod():
        return [
            {
                "name": "Phospho",
                "accession": "UNIMOD:21",
                "positions": [{"position": 4, "amino_acid": "T", "scores": None}],
            }
        ]

    @staticmethod
    def _sii(sii_id, site_prob):
        # MS:1001969 == phosphoRS site probability (a site-localization cvParam).
        return {
            "id": sii_id,
            "charge": 2,
            "exp_mz": 500.0,
            "calc_mz": 500.0,
            "rank": 1,
            "pass_threshold": True,
            "peptide_ref": "SharedPep",
            "cv_params": [
                {
                    "accession": "MS:1001969",
                    "name": "phosphoRS site probability",
                    "value": f"T4: {site_prob}",
                }
            ],
            "peptide_evidence_refs": [],
        }

    def _parsed(self):
        return {
            "spectra_data": {"SD1": "run.mgf"},
            "db_sequences": {},
            "peptides": {
                "SharedPep": {
                    "sequence": "PEPTIDE",
                    "modifications": self._phospho_mod(),
                    "xl_mods": {},
                }
            },
            "peptide_evidence": {},
            "linker_info": {},
            "spectrum_results": [
                {
                    "spectrum_id": "index=1",
                    "spectra_data_ref": "SD1",
                    "rt": None,
                    "siis": [self._sii("SII_1", 99.0)],
                },
                {
                    "spectrum_id": "index=2",
                    "spectra_data_ref": "SD1",
                    "rt": None,
                    "siis": [self._sii("SII_2", 10.0)],
                },
            ],
        }

    def test_each_psm_keeps_only_its_own_score(self):
        parsed = self._parsed()
        with MzIdentMLPsmAdapter() as adapter:
            records = adapter._build_psm_records(parsed)

        assert len(records) == 2
        scores_per_psm = [rec["modifications"][0]["positions"][0]["scores"] for rec in records]

        # Each PSM carries exactly ONE site-localization score — its own.
        for scores in scores_per_psm:
            assert scores is not None
            assert len(scores) == 1

        values = sorted(scores[0]["score_value"] for scores in scores_per_psm)
        assert values == [10.0, 99.0]

    def test_shared_source_modifications_untouched(self):
        parsed = self._parsed()
        shared = parsed["peptides"]["SharedPep"]["modifications"]
        with MzIdentMLPsmAdapter() as adapter:
            adapter._build_psm_records(parsed)

        # The Peptide's shared modification list must remain pristine — scoring
        # happened on per-PSM copies, not on this source object.
        assert shared[0]["positions"][0]["scores"] is None
