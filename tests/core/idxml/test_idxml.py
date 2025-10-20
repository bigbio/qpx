import unittest
from pathlib import Path
import tempfile
import shutil
import pandas as pd

from quantmsio.core.idxml import IdXML


class TestIdXML(unittest.TestCase):

    def setUp(self):
        self.test_data_dir = Path(__file__).parents[2] / "examples" / "idxml"
        self.test_idxml_file = (
            self.test_data_dir
            / "SF_200217_pPeptideLibrary_pool1_HCDnlETcaD_OT_rep2_consensus_fdr_pep_luciphor.idXML"
        )
        self.test_mzml_file = (
            self.test_data_dir
            / "SF_200217_pPeptideLibrary_pool1_HCDnlETcaD_OT_rep1.mzML"
        )
        self.temp_dir = Path(tempfile.mkdtemp())

        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                self.idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                self.idxml = IdXML(self.test_idxml_file)
        else:
            self.idxml = None

    def tearDown(self):
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def test_idxml_initialization(self):
        """Test IdXML initialization and basic properties"""
        if self.idxml is None:
            self.skipTest("Test IdXML file not found")

        self.assertIsInstance(self.idxml, IdXML)
        self.assertEqual(self.idxml.idxml_path, self.test_idxml_file)
        self.assertGreaterEqual(self.idxml.get_protein_count(), 0)
        self.assertGreaterEqual(self.idxml.get_psm_count(), 0)

    def test_dataframe_conversion(self):
        """Test DataFrame conversion and required columns"""
        if self.idxml is None:
            self.skipTest("Test IdXML file not found")

        df = self.idxml.to_dataframe()
        self.assertIsInstance(df, pd.DataFrame)

        if not df.empty:
            required_columns = [
                "sequence",
                "peptidoform",
                "modifications",
                "precursor_charge",
                "posterior_error_probability",
                "is_decoy",
                "calculated_mz",
                "observed_mz",
                "additional_scores",
                "protein_accessions",
                "rt",
                "reference_file_name",
                "scan",
                "q_value",
                "ion_mobility",
            ]

            for col in required_columns:
                self.assertIn(col, df.columns)

            self.assertIsNotNone(df["sequence"].iloc[0])
            self.assertIsNotNone(df["precursor_charge"].iloc[0])
            self.assertIsNotNone(df["calculated_mz"].iloc[0])

    def test_parquet_conversion(self):
        """Test parquet file conversion"""
        if self.idxml is None:
            self.skipTest("Test IdXML file not found")

        output_path = self.temp_dir / "test_output.parquet"
        self.idxml.to_parquet(output_path)

        self.assertTrue(output_path.exists())
        self.assertGreater(output_path.stat().st_size, 0)

    def test_modification_parsing(self):
        """Test modification parsing functionality"""
        if self.idxml is None:
            self.skipTest("Test IdXML file not found")

        test_sequence = "GSGEKPVSAPGDDT(Phospho)ES(Phospho)LHSQGEEEFDMPQPPHGHVLHR"
        modifications = self.idxml._parse_modifications(test_sequence)

        self.assertIsInstance(modifications, list)
        self.assertEqual(len(modifications), 1)
        mod = modifications[0]
        self.assertIn("name", mod)
        self.assertIn("positions", mod)
        self.assertEqual(mod["accession"], "UniMod:21")
        self.assertEqual(len(mod["positions"]), 2)
        for position_entry in mod["positions"]:
            self.assertIn("position", position_entry)
            self.assertIn("scores", position_entry)
            self.assertIsInstance(position_entry["position"], str)
            self.assertIsInstance(position_entry["scores"], list)
            # Check scores structure
            for score in position_entry["scores"]:
                self.assertIn("score_name", score)
                self.assertIn("score_value", score)
                self.assertEqual(score["score_name"], "localization_probability")
                self.assertIsInstance(score["score_value"], float)
        positions = [pos["position"] for pos in mod["positions"]]
        self.assertIn("T.14", positions)
        self.assertIn("S.16", positions)

    def test_scan_number_extraction(self):
        """Test scan number extraction from spectrum reference"""
        if self.idxml is None:
            self.skipTest("Test IdXML file not found")

        test_spectrum_ref = "controllerType=0 controllerNumber=1 scan=37144"
        scan_number = self.idxml._extract_scan_number(test_spectrum_ref)
        self.assertEqual(scan_number, "37144")

        test_spectrum_ref_no_scan = "controllerType=0 controllerNumber=1"
        scan_number_no_scan = self.idxml._extract_scan_number(test_spectrum_ref_no_scan)
        self.assertEqual(scan_number_no_scan, "unknown_index")

    def test_theoretical_mz_calculation(self):
        """Test theoretical m/z calculation with various cases"""
        if self.idxml is None:
            self.skipTest("Test IdXML file not found")

        test_cases = [
            ("PEPTIDE", 2, "Simple peptide"),
            ("PEPTIDE(Oxidation)", 2, "Peptide with oxidation"),
            ("PEPTIDE(Phospho)", 3, "Peptide with phosphorylation"),
            ("A", 1, "Single amino acid"),
            ("PEPTIDE", 0, "Zero charge"),
        ]

        for sequence, charge, description in test_cases:
            with self.subTest(description=description):
                calculated_mz = self.idxml._calculate_theoretical_mz(sequence, charge)
                self.assertIsInstance(calculated_mz, float)
                self.assertGreater(calculated_mz, 0)

                if charge > 1:
                    mz_charge_1 = self.idxml._calculate_theoretical_mz(sequence, 1)
                    self.assertGreater(mz_charge_1, calculated_mz)

    def test_theoretical_mz_consistency(self):
        """Test m/z calculation consistency with pyopenms"""
        if self.idxml is None:
            self.skipTest("Test IdXML file not found")

        try:
            import pyopenms as oms
            from pyopenms.Constants import PROTON_MASS_U
        except ImportError:
            self.skipTest("pyopenms not available")

        test_sequence = "PEPTIDE"
        test_charge = 2

        our_mz = self.idxml._calculate_theoretical_mz(test_sequence, test_charge)

        aa_sequence = oms.AASequence.fromString(test_sequence)
        peptide_mass = aa_sequence.getMonoWeight()
        expected_mz = (peptide_mass + (test_charge * PROTON_MASS_U)) / test_charge

        self.assertAlmostEqual(our_mz, expected_mz, places=6)

    def test_ion_mobility_extraction(self):
        """Test inverse reduced ion mobility extraction functionality"""
        if self.idxml is None:
            self.skipTest("Test IdXML file not found")

        self.assertTrue(hasattr(self.idxml, "_extract_ion_mobility"))

        try:
            import pyopenms as oms
        except ImportError:
            self.skipTest("pyopenms not available")

        peptide_hit = oms.PeptideHit()
        test_value = 0.9108814282059111
        peptide_hit.setMetaValue("inverse reduced ion mobility", test_value)

        extracted = self.idxml._extract_ion_mobility(peptide_hit)
        self.assertEqual(extracted, test_value)

        peptide_hit_unsupported = oms.PeptideHit()
        peptide_hit_unsupported.setMetaValue("ion_mobility", 0.5)
        peptide_hit_unsupported.setMetaValue("drift_time", 0.5)

        extracted_unsupported = self.idxml._extract_ion_mobility(
            peptide_hit_unsupported
        )
        self.assertIsNone(extracted_unsupported)

        df = self.idxml.to_dataframe()
        if not df.empty:
            self.assertIn("ion_mobility", df.columns)
            ion_mobility_values = df["ion_mobility"].dropna()
            for value in ion_mobility_values:
                self.assertIsInstance(value, (int, float))
                self.assertGreaterEqual(value, 0)

    def test_ondisc_experiment_initialization(self):
        """Test IdXML initialization with OnDiscExperiment for memory optimization"""
        if self.idxml is None:
            self.skipTest("Test IdXML file not found")

        # Test with OnDiscExperiment enabled
        if self.test_mzml_file.exists():
            idxml_ondisc = IdXML(
                self.test_idxml_file, self.test_mzml_file, use_ondisc=True
            )
            self.assertIsInstance(idxml_ondisc, IdXML)
            self.assertTrue(hasattr(idxml_ondisc, "_use_ondisc"))
            self.assertTrue(idxml_ondisc._use_ondisc)
            self.assertGreaterEqual(idxml_ondisc.get_protein_count(), 0)
            self.assertGreaterEqual(idxml_ondisc.get_psm_count(), 0)


if __name__ == "__main__":
    unittest.main()
