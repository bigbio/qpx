"""
IdXML parser for quantmsio package.
This module provides functionality to parse OpenMS IdXML files and convert them to quantms.io PSM format.
"""

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
from typing import Union, Optional, Dict, List
import logging
import re
import pyopenms as oms
from pyopenms.Constants import PROTON_MASS_U
from quantmsio.core.openms import OpenMSHandler
from quantmsio.core.format import PSM_SCHEMA


logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


class IdXML:
    """
    Parser for OpenMS IdXML files.

    This class provides functionality to parse IdXML files and convert them to quantms.io PSM format.
    IdXML is an OpenMS format for storing peptide and protein identifications.
    """

    def __init__(
        self,
        idxml_path: Union[Path, str],
        mzml_path: Optional[Union[Path, str]] = None,
        use_ondisc: bool = False,
        spectral_data: bool = False,
    ):
        """
        Initialize the IdXML parser.

        :param idxml_path: Path to the IdXML file
        :param mzml_path: Optional path to the mzML file for attaching spectra
        :param use_ondisc: Whether to use OnDiscExperiment for memory optimization
        """
        self.idxml_path = Path(idxml_path)
        self._mzml_path: Optional[Path] = Path(mzml_path) if mzml_path else None
        self._use_ondisc = use_ondisc
        self._spectral_data = spectral_data

        if self._spectral_data:
            logger.info("Loading spectra information into quantms.io")
        else:
            logger.info("Spectra information will not be loaded into quantms.io")

        self._protein_map = {}
        self._peptide_identifications = []
        self._parse_with_pyopenms()
        if self._mzml_path is not None:
            self._attach_mzml_spectra()

    def _parse_with_pyopenms(self) -> None:
        """Parse IdXML file using pyopenms"""
        try:
            protein_identifications, peptide_identifications = self._load_idxml_file()
            self._parse_proteins(protein_identifications)
            self._parse_peptides(peptide_identifications)
        except Exception as e:
            logging.error(f"Error parsing IdXML file with pyopenms: {e}")
            raise

    def _load_idxml_file(self) -> tuple[list, list]:
        """Load IdXML file and return protein and peptide identifications"""
        protein_identifications = []
        peptide_identifications = []
        oms.IdXMLFile().load(
            str(self.idxml_path), protein_identifications, peptide_identifications
        )
        return protein_identifications, peptide_identifications

    def _parse_proteins(self, protein_identifications: list) -> None:
        """Parse protein identifications"""
        for protein_id in protein_identifications:
            for protein_hit in protein_id.getHits():
                self._process_protein_hit(protein_hit)

    def _process_protein_hit(self, protein_hit) -> None:
        """Process single protein hit"""
        accession = protein_hit.getAccession()
        if accession:
            is_decoy = self._determine_is_decoy_generic(protein_hit)
            self._protein_map[accession] = {
                "is_decoy": is_decoy,
                "accession": accession,
            }

    def _determine_is_decoy_generic(self, hit_object) -> int:
        """Determine if a hit is decoy based on target_decoy metadata value only"""
        try:
            target_decoy_value = hit_object.getMetaValue("target_decoy")
            if target_decoy_value is not None:
                return 1 if str(target_decoy_value).lower() == "decoy" else 0
        except (AttributeError, ValueError, TypeError):
            pass
        return 0

    def _parse_peptides(self, peptide_identifications: list) -> None:
        """Parse peptide identifications"""
        for peptide_id in peptide_identifications:
            mz = peptide_id.getMZ()
            rt = peptide_id.getRT()
            spectrum_ref = self._extract_spectrum_reference(peptide_id)
            scan = self._extract_scan_number(spectrum_ref)
            reference_file_name = self._extract_reference_file_name(spectrum_ref)

            for peptide_hit in peptide_id.getHits():
                peptide_data = self._parse_peptide_hit_pyopenms(
                    peptide_hit, mz, rt, scan, reference_file_name
                )
                if peptide_data:
                    self._peptide_identifications.append(peptide_data)

    def _extract_spectrum_reference(self, peptide_id) -> str:
        """Extract spectrum reference from peptide identification"""
        try:
            return peptide_id.getSpectrumReference()
        except (AttributeError, ValueError, TypeError):
            try:
                return peptide_id.getMetaValue("spectrum_reference")
            except (AttributeError, ValueError, TypeError):
                return ""

    def _parse_peptide_hit_pyopenms(
        self,
        peptide_hit,
        mz: float,
        rt: float,
        scan: str,
        reference_file_name: str,
    ) -> Optional[Dict]:
        """Parse single peptide hit using pyopenms"""
        try:
            aa_sequence = peptide_hit.getSequence()
            if not aa_sequence:
                return None

            sequence = aa_sequence.toString()
            clean_sequence = aa_sequence.toUnmodifiedString()
            charge = peptide_hit.getCharge()
            protein_accessions = self._extract_protein_accessions(peptide_hit)
            is_decoy = self._determine_is_decoy_generic(peptide_hit)
            modifications = self._parse_modifications(sequence)
            additional_scores, q_value, posterior_error_probability = (
                self._extract_scores(peptide_hit)
            )
            calculated_mz = self._calculate_theoretical_mz(sequence, charge)

            if self._spectral_data:
                ion_mobility = self._extract_ion_mobility(peptide_hit)
            else:
                ion_mobility = None

            cv_params = self._extract_cv_params(peptide_hit)

            return {
                "sequence": clean_sequence,
                "peptidoform": sequence,
                "modifications": modifications,
                "precursor_charge": charge,
                "posterior_error_probability": posterior_error_probability,
                "is_decoy": is_decoy,
                "calculated_mz": calculated_mz,
                "observed_mz": mz,
                "additional_scores": additional_scores,
                "protein_accessions": protein_accessions,
                "predicted_rt": None,
                "reference_file_name": reference_file_name,
                "cv_params": cv_params,
                "scan": scan,
                "rt": rt,
                "q_value": q_value,
                "ion_mobility": ion_mobility,
                "number_peaks": None,  # Will be filled from mzML file
                "mz_array": None,  # Will be filled from mzML file
                "intensity_array": None,  # Will be filled from mzML file
                "charge_array": None,
                "ion_type_array": None,
                "ion_mobility_array": None,
            }

        except Exception as e:
            logging.warning(f"Error parsing peptide hit with pyopenms: {e}")
            return None

    def _extract_protein_accessions(self, peptide_hit) -> List[str]:
        """Extract protein accessions from peptide hit"""
        return [ref.getProteinAccession() for ref in peptide_hit.getPeptideEvidences()]

    def _extract_scores(
        self, peptide_hit
    ) -> tuple[List[Dict], Optional[float], Optional[float]]:
        """Extract additional scores, q-value, and posterior error probability"""
        additional_scores = []
        q_value = None
        posterior_error_probability = None

        try:
            q_value_value = peptide_hit.getMetaValue("q-value")
            if q_value_value is not None:
                q_value = float(q_value_value)

            pep_value = peptide_hit.getMetaValue("Posterior Error Probability_score")
            if pep_value is not None:
                posterior_error_probability = float(pep_value)

            additional_scores = self._extract_additional_scores(peptide_hit)
        except (AttributeError, ValueError, TypeError):
            pass

        return additional_scores, q_value, posterior_error_probability

    def _extract_ion_mobility(self, peptide_hit) -> Optional[float]:
        """Extract inverse reduced ion mobility from peptide hit using pyopenms"""
        try:
            ion_mobility_value = peptide_hit.getMetaValue(
                "inverse reduced ion mobility"
            )
            if ion_mobility_value is not None:
                return float(ion_mobility_value)
        except (AttributeError, ValueError, TypeError):
            pass
        return None

    def _extract_cv_params(self, peptide_hit) -> Optional[List[Dict]]:
        """Extract controlled vocabulary parameters from peptide hit"""
        try:
            consensus_support_value = peptide_hit.getMetaValue("consensus_support")
            if consensus_support_value is not None:
                consensus_support = float(consensus_support_value)
                return [
                    {"cv_name": "consensus_support", "cv_value": str(consensus_support)}
                ]
        except (AttributeError, ValueError, TypeError):
            pass
        return None

    def _extract_additional_scores(self, peptide_hit) -> List[Dict]:
        """Extract additional scores from peptide hit - all meta values except core PSM fields"""
        additional_scores = []
        core_psm_fields = {
            "q-value",
            "Posterior Error Probability_score",
            "target_decoy",
            "consensus_support",
            "inverse reduced ion mobility",
            "spectrum_reference",
        }

        important_scores = [
            "Luciphor_pep_score",
            "Luciphor_global_flr",
            "Luciphor_local_flr",
        ]

        for key in important_scores:
            if key in core_psm_fields:
                continue

            try:
                score_value_raw = peptide_hit.getMetaValue(key)
                if score_value_raw is not None and isinstance(
                    score_value_raw, (int, float)
                ):
                    additional_scores.append(
                        {"score_name": key, "score_value": float(score_value_raw)}
                    )
            except (AttributeError, KeyError):
                continue

        return additional_scores

    def _extract_modification_info(self, modification):
        """
        Extract name, accession, and ID from an OpenMS modification object.

        :param modification: OpenMS modification object
        :return: Tuple of (mod_name, mod_accession, mod_id)
        """
        mod_name = (
            modification.getId() or modification.getFullId() or modification.getName()
        )
        mod_accession = modification.getUniModAccession()
        mod_id = modification.getId()

        return mod_name, mod_accession, mod_id

    def _get_key_identifier(
        self, mod_name: str, mod_accession: str, mod_id: str
    ) -> str:
        """
        Determine the key identifier for a modification.
        Priority: UniMod accession > ID > name > "Unknown"

        :param mod_name: Modification name
        :param mod_accession: UniMod accession
        :param mod_id: Modification ID
        :return: Key identifier string
        """
        if mod_accession:
            return mod_accession
        elif mod_id:
            return mod_id
        elif mod_name:
            return mod_name
        else:
            return "Unknown"

    def _create_position_entry(self, position: str) -> dict:
        """
        Create a position entry with default localization probability.

        :param position: Position string (e.g., "N-term.0", "A.1")
        :return: Position entry dictionary
        """
        return {
            "position": position,
            "scores": [{"score_name": "localization_probability", "score_value": 1.0}],
        }

    def _add_modification_to_dict(
        self,
        modifications: dict,
        key_identifier: str,
        mod_name: str,
        mod_accession: str,
        position_entry: dict,
    ):
        """
        Add or update a modification in the modifications dictionary.

        :param modifications: Modifications dictionary
        :param key_identifier: Unique key for the modification
        :param mod_name: Modification name
        :param mod_accession: UniMod accession
        :param position_entry: Position entry to add
        """
        if key_identifier not in modifications:
            modifications[key_identifier] = {
                "name": mod_name,
                "accession": mod_accession,
                "positions": [],
            }
        modifications[key_identifier]["positions"].append(position_entry)

    def _process_n_terminal_modification(self, aa_sequence, modifications: dict):
        """Process N-terminal modification if present."""
        if not aa_sequence.hasNTerminalModification():
            return

        n_term_mod = aa_sequence.getNTerminalModification()
        mod_name, mod_accession, mod_id = self._extract_modification_info(n_term_mod)
        key_identifier = self._get_key_identifier(mod_name, mod_accession, mod_id)
        position_entry = self._create_position_entry("N-term.0")

        self._add_modification_to_dict(
            modifications, key_identifier, mod_name, mod_accession, position_entry
        )

    def _process_c_terminal_modification(self, aa_sequence, modifications: dict):
        """Process C-terminal modification if present."""
        if not aa_sequence.hasCTerminalModification():
            return

        c_term_mod = aa_sequence.getCTerminalModification()
        mod_name, mod_accession, mod_id = self._extract_modification_info(c_term_mod)
        key_identifier = self._get_key_identifier(mod_name, mod_accession, mod_id)
        position = f"C-term.{aa_sequence.size() + 1}"
        position_entry = self._create_position_entry(position)

        self._add_modification_to_dict(
            modifications, key_identifier, mod_name, mod_accession, position_entry
        )

    def _process_amino_acid_modifications(self, aa_sequence, modifications: dict):
        """Process modifications on amino acid residues."""
        for i in range(aa_sequence.size()):
            residue = aa_sequence.getResidue(i)
            if not residue.isModified():
                continue

            mod = residue.getModification()
            mod_name, mod_accession, mod_id = self._extract_modification_info(mod)
            key_identifier = self._get_key_identifier(mod_name, mod_accession, mod_id)

            if key_identifier:
                position = f"{residue.getOneLetterCode()}.{i + 1}"
                position_entry = self._create_position_entry(position)

                self._add_modification_to_dict(
                    modifications,
                    key_identifier,
                    mod_name,
                    mod_accession,
                    position_entry,
                )

    def _parse_modifications(self, sequence: str) -> List[Dict]:
        """
        Parse modifications from peptide sequence using pyopenms.
        This provides accurate and standardized modification parsing.
        Supports multiple positions for the same modification type.

        :param sequence: Peptide sequence with modification annotations
        :return: List of modification dictionaries with new structure containing name, accession, and positions with scores
        """
        aa_sequence = oms.AASequence.fromString(sequence)
        modifications = {}

        # Process different types of modifications
        self._process_n_terminal_modification(aa_sequence, modifications)
        self._process_c_terminal_modification(aa_sequence, modifications)
        self._process_amino_acid_modifications(aa_sequence, modifications)

        return list(modifications.values())

    def _extract_scan_number(self, spectrum_ref: str) -> str:
        """
        Extract scan number from spectrum reference string.

        :param spectrum_ref: Spectrum reference string
        :return: Extracted scan number
        """
        scan_match = re.search(r"scan=(\d+)", spectrum_ref)
        return scan_match.group(1) if scan_match else "unknown_index"

    def _extract_reference_file_name(self, spectrum_ref: str) -> str:
        """
        Extract reference file name from spectrum reference string.
        If no file name is found in spectrum_ref, use the mzML file name if available.

        :param spectrum_ref: Spectrum reference string
        :return: Reference file name
        """
        file_match = re.search(r"file=([^,\s]+)", spectrum_ref)
        if file_match:
            return file_match.group(1)

        return self._mzml_path.stem if self._mzml_path else ""

    def _calculate_theoretical_mz(self, sequence: str, charge: int) -> float:
        """
        Calculate theoretical m/z for a peptide sequence using pyopenms.
        This provides accurate calculation with proper modification support.

        :param sequence: Peptide sequence (can include modifications in ProForma notation)
        :param charge: Charge state
        :return: Theoretical m/z value
        """
        aa_sequence = oms.AASequence.fromString(sequence)
        peptide_mass = aa_sequence.getMonoWeight()

        if charge > 0:
            return (peptide_mass + (charge * PROTON_MASS_U)) / charge
        else:
            return peptide_mass

    def _format_modifications_for_dataframe(self, psm: dict):
        """Format modifications structure for DataFrame compatibility."""
        if psm.get("modifications") is None:
            return

        modifications = psm["modifications"]
        if isinstance(modifications, list):
            for mod in modifications:
                if isinstance(mod, dict) and "positions" in mod:
                    # Ensure positions is a list
                    if hasattr(mod["positions"], "tolist"):
                        mod["positions"] = mod["positions"].tolist()
                    # Ensure scores in each position is a list
                    for pos in mod["positions"]:
                        if isinstance(pos, dict) and "scores" in pos:
                            if hasattr(pos["scores"], "tolist"):
                                pos["scores"] = pos["scores"].tolist()

    def _convert_field_to_list(self, psm: dict, field_name: str):
        """Convert a field to list format if it has tolist method."""
        if psm.get(field_name) is not None:
            if hasattr(psm[field_name], "tolist"):
                psm[field_name] = psm[field_name].tolist()

    def _format_list_fields_for_dataframe(self, psm: dict):
        """Format list fields for DataFrame compatibility."""
        list_fields = [
            "protein_accessions",
            "additional_scores",
            "cv_params",
            "mz_array",
            "intensity_array",
        ]

        for field in list_fields:
            self._convert_field_to_list(psm, field)

    def _get_required_columns(self) -> list:
        """Get list of required columns for the DataFrame."""
        return [
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
            "predicted_rt",
            "reference_file_name",
            "cv_params",
            "scan",
            "rt",
            "q_value",
            "ion_mobility",
            "number_peaks",
            "mz_array",
            "intensity_array",
            "charge_array",
            "ion_type_array",
            "ion_mobility_array",
        ]

    def _ensure_required_columns(self, df: pd.DataFrame):
        """Ensure all required columns exist in the DataFrame."""
        required_columns = self._get_required_columns()
        for col in required_columns:
            if col not in df.columns:
                df[col] = None

    def _apply_column_specific_formatting(self, df: pd.DataFrame):
        """Apply specific formatting to certain columns."""
        if "ion_mobility" in df.columns:
            df["ion_mobility"] = df["ion_mobility"].where(
                df["ion_mobility"].notna(), None
            )

    def _prepare_data_for_dataframe(self):
        """Prepare peptide identification data for DataFrame creation."""
        for psm in self._peptide_identifications:
            self._format_modifications_for_dataframe(psm)
            self._format_list_fields_for_dataframe(psm)

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert parsed peptide identifications to a pandas DataFrame.

        :return: DataFrame containing PSM data
        """
        if not self._peptide_identifications:
            return pd.DataFrame()

        # Prepare data by formatting complex structures
        self._prepare_data_for_dataframe()

        # Create DataFrame
        df = pd.DataFrame(self._peptide_identifications)

        # Ensure required columns exist
        self._ensure_required_columns(df)

        # Apply column-specific formatting
        self._apply_column_specific_formatting(df)

        return df

    def _attach_mzml_spectra(self) -> None:
        """
        Attach number_peaks, mz_array, intensity_array using the provided mzML file.
        Uses scan number as identifier. Supports both traditional and OnDisc methods.
        """
        try:
            if self._use_ondisc:
                self._attach_mzml_spectra_ondisc_impl()
            else:
                self._attach_mzml_spectra_traditional_impl()
        except Exception as e:
            method_name = "OnDiscExperiment" if self._use_ondisc else "traditional"
            logging.warning(f"Failed to attach mzML spectra using {method_name}: {e}")

    def _attach_mzml_spectra_traditional_impl(self) -> None:
        """Traditional implementation using OpenMSHandler and SpectrumLookup"""
        handler = OpenMSHandler()
        for psm in self._peptide_identifications:
            self._attach_spectrum_to_psm_traditional(psm, handler)

    def _attach_mzml_spectra_ondisc_impl(self) -> None:
        """OnDisc implementation using OnDiscMSExperiment for memory optimization"""
        # Create OnDiscExperiment for memory-efficient access
        ondisc_exp = oms.OnDiscMSExperiment()
        oms.MzMLFile().load(str(self._mzml_path), ondisc_exp)

        # Build scan number to index mapping
        scan_to_index = self._build_scan_to_index_mapping(ondisc_exp)

        for psm in self._peptide_identifications:
            self._attach_spectrum_to_psm_ondisc(psm, ondisc_exp, scan_to_index)

    def _build_scan_to_index_mapping(self, ondisc_exp) -> dict:
        """Build mapping from scan number to spectrum index for OnDiscExperiment"""
        scan_to_index = {}
        for i in range(ondisc_exp.size()):
            spectrum = ondisc_exp.getSpectrum(i)
            scan_number = spectrum.getNativeID()
            if scan_number:
                scan_match = re.search(r"scan=(\d+)", scan_number)
                if scan_match:
                    scan_to_index[int(scan_match.group(1))] = i
        return scan_to_index

    def _attach_spectrum_to_psm_traditional(self, psm: dict, handler) -> None:
        """Attach spectrum data to PSM using traditional method"""
        scan_value = psm.get("scan")
        if not scan_value or scan_value == "unknown_index":
            return

        try:
            scan_int = int(str(scan_value))
            num_peaks, mzs, intens = handler.get_spectrum_from_scan(
                str(self._mzml_path), scan_int
            )
            psm["number_peaks"] = int(num_peaks)
            psm["mz_array"] = [float(x) for x in mzs]
            psm["intensity_array"] = [float(x) for x in intens]
        except (ValueError, Exception) as e:
            logging.warning(f"Error fetching spectrum (scan={scan_value}): {e}")
            psm["number_peaks"] = psm["mz_array"] = psm["intensity_array"] = None

    def _attach_spectrum_to_psm_ondisc(
        self, psm: dict, ondisc_exp, scan_to_index: dict
    ) -> None:
        """Attach spectrum data to PSM using OnDisc method"""
        scan_value = psm.get("scan")
        if not scan_value or scan_value == "unknown_index":
            return

        try:
            scan_int = int(str(scan_value))
            if scan_int in scan_to_index:
                index = scan_to_index[scan_int]
                spectrum = ondisc_exp.getSpectrum(index)
                spectrum_mz, spectrum_intensities = spectrum.get_peaks()

                psm["number_peaks"] = int(len(spectrum_mz))
                psm["mz_array"] = [float(x) for x in spectrum_mz]
                psm["intensity_array"] = [float(x) for x in spectrum_intensities]
            else:
                psm["number_peaks"] = psm["mz_array"] = psm["intensity_array"] = None

        except (ValueError, IndexError, Exception) as e:
            logging.warning(f"Error fetching spectrum (scan={scan_value}): {e}")
            psm["number_peaks"] = psm["mz_array"] = psm["intensity_array"] = None

    def to_parquet(self, output_path: Union[Path, str]) -> None:
        """
        Convert IdXML data to parquet format and save to file.

        :param output_path: Output file path for parquet file
        """
        df = self.to_dataframe()
        if df.empty:
            logging.warning("No peptide identifications found to convert")
            return

        def convert_numpy_to_list(obj):
            """Recursively convert numpy arrays to Python lists"""
            if hasattr(obj, "tolist"):
                return obj.tolist()
            elif isinstance(obj, list):
                return [convert_numpy_to_list(item) for item in obj]
            elif isinstance(obj, dict):
                return {key: convert_numpy_to_list(value) for key, value in obj.items()}
            else:
                return obj

        for col in df.columns:
            if col in [
                "modifications",
                "protein_accessions",
                "additional_scores",
                "cv_params",
                "mz_array",
                "intensity_array",
            ]:
                df[col] = df[col].apply(convert_numpy_to_list)

        table = pa.Table.from_pandas(df, schema=PSM_SCHEMA)
        pq.write_table(table, output_path)
        logging.info(f"Successfully converted IdXML to parquet: {output_path}")

    def get_psm_count(self) -> int:
        """
        Get the total number of PSMs found in the IdXML file.

        :return: Number of PSMs
        """
        return len(self._peptide_identifications)

    def get_protein_count(self) -> int:
        """
        Get the total number of proteins found in the IdXML file.

        :return: Number of proteins
        """
        return len(self._protein_map)


def merge_idxml_parquet_files(
    parquet_files: List[Union[Path, str]],
    output_path: Union[Path, str],
) -> None:
    """
    Merge multiple PSM parquet files into a single parquet file.

    :param parquet_files: List of parquet file paths to merge
    :param output_path: Output path for the merged parquet file
    """
    if not parquet_files:
        raise ValueError("No parquet files provided for merging")

    logger.info(f"Merging {len(parquet_files)} parquet files into {output_path}")

    # Read all parquet files into dataframes
    dfs = []
    for parquet_file in parquet_files:
        logger.info(f"Reading {parquet_file}")
        df = pd.read_parquet(parquet_file)
        dfs.append(df)

    # Concatenate all dataframes
    merged_df = pd.concat(dfs, ignore_index=True)
    logger.info(f"Total PSMs in merged file: {len(merged_df)}")

    # Convert numpy arrays to lists for parquet compatibility
    def convert_numpy_to_list(obj):
        """Recursively convert numpy arrays to Python lists"""
        if hasattr(obj, "tolist"):
            return obj.tolist()
        elif isinstance(obj, list):
            return [convert_numpy_to_list(item) for item in obj]
        elif isinstance(obj, dict):
            return {key: convert_numpy_to_list(value) for key, value in obj.items()}
        else:
            return obj

    for col in merged_df.columns:
        if col in [
            "modifications",
            "protein_accessions",
            "additional_scores",
            "cv_params",
            "mz_array",
            "intensity_array",
        ]:
            merged_df[col] = merged_df[col].apply(convert_numpy_to_list)

    # Write merged dataframe to parquet
    table = pa.Table.from_pandas(merged_df, schema=PSM_SCHEMA)
    pq.write_table(table, output_path, compression="gzip")
    logger.info(f"Successfully merged parquet files to: {output_path}")
