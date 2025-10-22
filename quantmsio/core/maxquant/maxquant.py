"""MaxQuant data processing module"""

import logging
import re
from pathlib import Path
from typing import Union, List, Dict, Optional, Tuple

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pyopenms import AASequence
from pyopenms.Constants import PROTON_MASS_U

from quantmsio.core.format import PG_SCHEMA, FEATURE_SCHEMA, PSM_SCHEMA
from quantmsio.core.common import (
    MAXQUANT_PSM_MAP,
    MAXQUANT_FEATURE_MAP,
    MAXQUANT_PG_MAP,
    SDRF_MAP,
    MAXQUANT_PSM_USECOLS,
    MAXQUANT_FEATURE_USECOLS,
    MAXQUANT_PG_USECOLS,
    SDRF_USECOLS,
)
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.utils.file_utils import ParquetBatchWriter, extract_protein_list

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

TOKEN_DELIMITER_PATTERN = re.compile(r"[_\s\-]+")


# ============================================================================
# Utility Functions
# ============================================================================


def clean_peptidoform(peptidoform):
    """Clean peptidoform string by normalizing modification names"""
    if not isinstance(peptidoform, str):
        return ""

    peptidoform = peptidoform.strip("_")

    modification_mapping = {
        "ac": "Acetyl",
        "ox": "Oxidation",
        "me": "Methyl",
        "ph": "Phospho",
        "de": "Deamidated",
        "cam": "Carbamidomethyl",
        "dim": "Dimethyl",
        "tri": "Trimethyl",
        "ub": "GlyGly",
        "su": "Sumo",
    }  # Add more modification abbreviations as needed

    for short_name, full_name in modification_mapping.items():
        pattern = f"\\({re.escape(short_name)}\\)"
        replacement = f"({full_name})"
        peptidoform = re.sub(pattern, replacement, peptidoform)

    return peptidoform


def convert_maxquant_flag(value):
    """Convert MaxQuant flag from + to 1, others to 0"""
    return 1 if value == "+" else 0


def _process_modification(mod, modifications: dict, position: str) -> None:
    """Process a single modification and add it to modifications dict"""
    accession = mod.getUniModAccession()
    short_name = mod.getId()
    if not short_name:
        short_name = (
            mod.getFullName() if hasattr(mod, "getFullName") else f"UniMod:{accession}"
        )
    if short_name not in modifications:
        modifications[short_name] = {
            "name": short_name,
            "accession": accession,
            "positions": [],
        }
    modifications[short_name]["positions"].append({"position": position, "scores": []})


def parse_modifications_from_peptidoform(peptidoform: str) -> list:
    """Parse modification information from peptidoform string"""
    if not isinstance(peptidoform, str):
        return None

    try:
        cleaned_peptidoform = clean_peptidoform(peptidoform)
        if not cleaned_peptidoform:
            return None

        sequence = AASequence.fromString(cleaned_peptidoform)
        modifications = {}

        if sequence.hasNTerminalModification():
            mod = sequence.getNTerminalModification()
            _process_modification(mod, modifications, "N-term.0")

        if sequence.hasCTerminalModification():
            mod = sequence.getCTerminalModification()
            _process_modification(mod, modifications, f"C-term.{sequence.size()+1}")

        for i in range(sequence.size()):
            residue = sequence.getResidue(i)
            if residue.isModified():
                mod = residue.getModification()
                position = f"{residue.getOneLetterCode()}.{i+1}"
                _process_modification(mod, modifications, position)

        return list(modifications.values()) if modifications else None

    except Exception:
        return None


# ============================================================================
# MaxQuant Data Processor
# ============================================================================


class MaxQuant:
    """MaxQuant data processor for converting output files to quantms.io format"""

    def __init__(self, spectral_data: bool = False):
        self.sdrf_handler: Optional[SDRFHandler] = None
        self.experiment_type: Optional[str] = None
        self._sample_map: Optional[Dict] = None
        self._channel_map: Optional[Dict] = None
        self._current_sdrf_path: Optional[str] = None

        self._spectral_data = spectral_data

        self.psm_mapping = MAXQUANT_PSM_MAP
        self.feature_mapping = MAXQUANT_FEATURE_MAP
        self.pg_mapping = MAXQUANT_PG_MAP
        self.sdrf_mapping = SDRF_MAP

        self._sequence_cache: Dict[str, AASequence] = {}
        self._optimized_sdrf_data: Optional[List[Dict]] = None
        self._sdrf_search_cache: Dict[str, Tuple] = {}

    # ============================================================================
    # SDRF Processing
    # ============================================================================

    def _init_sdrf(self, sdrf_path: Union[Path, str]) -> None:
        """Initialize SDRF handler and create optimized mappings"""
        if not sdrf_path:
            return

        self.sdrf_handler = SDRFHandler(sdrf_path)
        self.experiment_type = self.sdrf_handler.get_experiment_type_from_sdrf()
        self._sample_map = self.sdrf_handler.get_sample_map_run()
        self._current_sdrf_path = sdrf_path

        self._sdrf_transformed = self._create_basic_sdrf_mapping()
        self._channel_map = self._create_simplified_channel_mapping()
        self._preprocess_sdrf_for_optimization()

    def _create_basic_sdrf_mapping(self) -> pd.DataFrame:
        """Create basic mapping table from SDRF file"""
        sdrf_table = self.sdrf_handler.sdrf_table.copy()

        basic_mapping = {
            "comment[data file]": "reference_file_name",
            "comment[label]": "channel",
            "source name": "sample_accession",
        }

        rename_map = {}
        for orig_col, new_col in basic_mapping.items():
            matching_col = None
            for col in sdrf_table.columns:
                if col.lower() == orig_col.lower():
                    matching_col = col
                    break

            if matching_col:
                rename_map[matching_col] = new_col

        available_cols = [col for col in sdrf_table.columns if col in rename_map]
        basic_df = sdrf_table[available_cols].copy()
        basic_df = basic_df.rename(columns=rename_map)

        if "reference_file_name" in basic_df.columns:
            basic_df["reference_file_name"] = (
                basic_df["reference_file_name"].fillna("").str.split(".").str[0]
            )
        if "reference_file_name" in basic_df.columns:
            basic_df.set_index("reference_file_name", inplace=True)

        return basic_df

    def _create_simplified_channel_mapping(self) -> Dict[str, str]:
        """Create channel mapping using basic SDRF information"""
        if not hasattr(self, "_sdrf_transformed") or self._sdrf_transformed is None:
            return {}

        channel_map = {}

        df = (
            self._sdrf_transformed.reset_index()
            if self._sdrf_transformed.index.name
            else self._sdrf_transformed
        )

        for _, row in df.iterrows():
            if "reference_file_name" in row:
                file_key = row["reference_file_name"]
            elif self._sdrf_transformed.index.name == "reference_file_name":
                file_key = row.name if hasattr(row, "name") else str(row.iloc[0])
            else:
                continue

            if self.experiment_type == "LFQ":
                if "channel" in row and pd.notna(row["channel"]):
                    channel_map[file_key] = row["channel"]
            else:
                if "channel" in row and "sample_accession" in row:
                    if pd.notna(row["channel"]) and pd.notna(row["sample_accession"]):
                        map_key = f"{file_key}-{row['channel']}"
                        channel_map[map_key] = row["sample_accession"]

        return channel_map

    def _get_sample_accession_from_sdrf(
        self, reference_file: str, channel: str = None
    ) -> Optional[str]:
        """Get sample accession from SDRF mapping"""
        if not self._sample_map:
            return reference_file if reference_file else None

        file_key = reference_file.split(".")[0] if reference_file else ""

        if channel and self.experiment_type != "LFQ":
            map_key = f"{file_key}-{channel}"
            return self._sample_map.get(map_key, reference_file)

        return reference_file if reference_file else None

    def _get_channel_from_sdrf(
        self, reference_file: str, sample_accession: str = None
    ) -> Optional[str]:
        """Get channel information from transformed SDRF data"""
        if not hasattr(self, "_sdrf_transformed") or self._sdrf_transformed is None:
            return None

        file_key = reference_file.split(".")[0] if reference_file else ""

        if (
            hasattr(self._sdrf_transformed, "index")
            and self._sdrf_transformed.index.name == "reference_file_name"
        ):
            if file_key in self._sdrf_transformed.index:
                matching_rows = self._sdrf_transformed.loc[[file_key]]
            else:
                matching_rows = pd.DataFrame()
        elif "reference_file_name" in self._sdrf_transformed.columns:
            matching_rows = self._sdrf_transformed[
                self._sdrf_transformed["reference_file_name"] == file_key
            ]
        else:
            matching_rows = pd.DataFrame()

        if not matching_rows.empty:
            if self.experiment_type == "LFQ":
                channel = matching_rows.iloc[0]["channel"]
                return channel if pd.notna(channel) else None
            else:
                if sample_accession:
                    sample_matches = matching_rows[
                        matching_rows["sample_accession"] == sample_accession
                    ]
                    if not sample_matches.empty:
                        channel = sample_matches.iloc[0]["channel"]
                        return channel if pd.notna(channel) else None

        return None

    def _process_tmt_intensities(self, row, tmt_channels, reference_file_name) -> tuple:
        """Process TMT Reporter intensity columns"""
        intensities = []
        additional_intensities = []

        for i in range(min(8, len(tmt_channels))):
            reporter_col = f"Reporter intensity {i}"
            corrected_col = f"Reporter intensity corrected {i}"

            if (
                reporter_col in row.index
                and pd.notna(row[reporter_col])
                and row[reporter_col] > 0
            ):
                channel_name = tmt_channels[i]

                sample_accession = self._get_tmt_sample_accession(
                    reference_file_name, channel_name
                )

                intensities.append(
                    {
                        "sample_accession": sample_accession,
                        "channel": channel_name,
                        "intensity": float(row[reporter_col]),
                    }
                )

                if corrected_col in row.index and pd.notna(row[corrected_col]):
                    additional_intensities.append(
                        {
                            "sample_accession": sample_accession,
                            "channel": channel_name,
                            "intensities": [
                                {
                                    "intensity_name": corrected_col,
                                    "intensity_value": float(row[corrected_col]),
                                }
                            ],
                        }
                    )

        return intensities, additional_intensities

    def _get_tmt_sample_accession(
        self, reference_file_name: str, channel_name: str
    ) -> str:
        """Get sample accession for TMT experiments"""
        if hasattr(self, "_sample_map") and self._sample_map:
            file_key = reference_file_name.split(".")[0]
            map_key = f"{file_key}-{channel_name}"
            sample_accession = self._sample_map.get(map_key)
            if sample_accession:
                return sample_accession

        return reference_file_name

    def _get_tmt_channels_from_sdrf(self) -> list:
        """Get TMT channel list from SDRF"""
        if not hasattr(self, "sdrf_handler") or not self.sdrf_handler:
            return []

        try:
            sdrf_table = self.sdrf_handler.sdrf_table
            if "comment[label]" in sdrf_table.columns:
                labels = sdrf_table["comment[label]"].unique()
                tmt_labels = [
                    label for label in labels if label and "TMT" in str(label).upper()
                ]
                return sorted(tmt_labels)
        except (AttributeError, KeyError) as e:
            logger.warning(f"Could not extract TMT channels from SDRF: {e}")
        except Exception as e:
            logger.error(f"Unexpected error extracting TMT channels: {e}")

        return []

    # ============================================================================
    # Core Processing Methods
    # ============================================================================

    @staticmethod
    def _split_semicolon_separated_column(series: pd.Series) -> pd.Series:
        """
        Split semicolon-separated string column into list of strings.

        Args:
            series: pandas Series containing semicolon-separated strings

        Returns:
            pandas Series with each value converted to a list of strings
        """
        return series.fillna("").astype(str).str.split(";")

    @staticmethod
    def _convert_maxquant_flag_column(series: pd.Series) -> pd.Series:
        """
        Convert MaxQuant flag column from + to 1, others to 0, and cast to int32.

        Args:
            series: pandas Series containing MaxQuant flag values ('+' or other)

        Returns:
            pandas Series with 1 for '+' values, 0 for others, as int32 type
        """
        return series.apply(convert_maxquant_flag).astype("int32")

    def _calculate_theoretical_mz_batch(self, df: pd.DataFrame) -> None:
        """Calculate theoretical m/z values in batch"""

        def safe_parse_sequence(peptidoform: str) -> Optional[AASequence]:
            """Safely parse peptide sequence"""
            try:
                cleaned_peptidoform = clean_peptidoform(peptidoform)
                if not cleaned_peptidoform:
                    return None
                return AASequence.fromString(cleaned_peptidoform)
            except Exception:
                try:
                    return AASequence(cleaned_peptidoform)
                except Exception:
                    return None

        unique_peptidoforms = df["peptidoform"].unique()

        for peptidoform in unique_peptidoforms:
            if peptidoform not in self._sequence_cache:
                self._sequence_cache[peptidoform] = safe_parse_sequence(peptidoform)

        mass_map = {}
        for peptidoform in unique_peptidoforms:
            sequence = self._sequence_cache.get(peptidoform)
            if sequence:
                mass_map[peptidoform] = sequence.getMonoWeight()
            else:
                mass_map[peptidoform] = 0.0

        mass_vector = df["peptidoform"].map(mass_map)
        df.loc[:, "calculated_mz"] = (
            mass_vector + (PROTON_MASS_U * df["precursor_charge"])
        ) / df["precursor_charge"]

    # ============================================================================
    # PSM Processing
    # ============================================================================

    def process_psm_file(
        self, msms_path: str, output_path: str, chunksize: int = 1000000
    ) -> None:
        """Process PSM data from msms.txt to PSM parquet formatt"""

        batch_writer = ParquetBatchWriter(output_path, PSM_SCHEMA)

        try:

            if self._spectral_data:
                logger.info("Loading spectra information into quantms.io")
            else:
                logger.info("Spectra information will not be loaded into quantms.io")

            for df_chunk in pd.read_csv(
                msms_path, sep="\t", chunksize=chunksize, low_memory=False
            ):

                df_chunk = self._apply_psm_mapping(df_chunk)
                self._calculate_theoretical_mz_batch(df_chunk)
                df_chunk = self._process_psm_modifications(df_chunk)
                df_chunk = self._process_psm_scores(df_chunk)
                df_chunk = self._process_psm_cv_params(df_chunk)
                df_chunk = self._process_psm_arrays(df_chunk)
                df_chunk = self._ensure_psm_schema_compliance(df_chunk)

                table = pa.Table.from_pandas(
                    df_chunk, schema=PSM_SCHEMA, preserve_index=False
                )
                batch_writer.write_batch(table.to_pylist())

        finally:
            batch_writer.close()

    def _apply_psm_mapping(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply MAXQUANT_PSM_MAP mapping"""
        available_mapping = {
            k: v for k, v in self.psm_mapping.items() if k in df.columns
        }
        df.rename(columns=available_mapping, inplace=True)

        if "protein_accessions" in df.columns:
            df["protein_accessions"] = self._split_semicolon_separated_column(
                df["protein_accessions"]
            )

        return df

    def _process_psm_modifications(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process PSM modification information"""
        if "peptidoform" in df.columns:
            df["modifications"] = df["peptidoform"].apply(
                parse_modifications_from_peptidoform
            )
        else:
            df["modifications"] = None
        return df

    def _process_psm_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """Structure PSM scoring information"""
        non_schema_score_fields = {
            "andromeda_score": "andromeda_score",
            "andromeda_delta_score": "andromeda_delta_score",
            "parent_ion_fraction": "parent_ion_fraction",
        }  # Add more additional score fields as needed

        scores_list = []
        for _, row in df.iterrows():
            scores = []
            for field_name, standard_name in non_schema_score_fields.items():
                if field_name in row and pd.notna(row[field_name]):
                    scores.append(
                        {
                            "score_name": standard_name,
                            "score_value": float(row[field_name]),
                        }
                    )
            scores_list.append(scores if scores else None)

        df["additional_scores"] = scores_list

        for field_name in non_schema_score_fields.keys():
            if field_name in df.columns:
                df = df.drop(columns=[field_name])

        return df

    def _process_cv_params(self, df: pd.DataFrame, cv_columns: list) -> pd.DataFrame:
        """Process CV parameters from data"""
        cv_params_list = []

        for _, row in df.iterrows():
            cv_params = []

            for cv_name in cv_columns:
                if cv_name in row and pd.notna(row[cv_name]):
                    cv_value = str(row[cv_name])
                    cv_params.append({"cv_name": cv_name, "cv_value": cv_value})

            cv_params_list.append(cv_params if cv_params else None)

        df["cv_params"] = cv_params_list
        return df

    def _process_psm_cv_params(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process CV parameters from PSM data"""
        cv_columns = [
            "Fragmentation",
            "Mass analyzer",
            "Type",
        ]  # Add more MaxQuant columns as needed
        return self._process_cv_params(df, cv_columns)

    def _process_feature_cv_params(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process CV parameters from Feature data"""
        cv_columns = ["Type"]  # Add more MaxQuant columns as needed
        return self._process_cv_params(df, cv_columns)

    def _process_psm_arrays(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process mz_array and intensity_array from semicolon-separated strings"""

        def parse_array_string(array_str):
            """Parse semicolon-separated string to float array"""
            if pd.isna(array_str) or array_str == "":
                return None
            try:
                values = [
                    float(x.strip()) for x in str(array_str).split(";") if x.strip()
                ]
                return values if values else None
            except (ValueError, AttributeError):
                return None

        def parse_matches(matches_str):
            """Parse 'Matches' to [charge_array, ion_type_array]"""
            if pd.isna(matches_str) or matches_str == "":
                return [None, None]

            try:
                ions = matches_str.split(";")
                ion_types = []
                charges = []

                for ion in ions:
                    charge_match = re.search(r"\((\d+)\+\)", ion)
                    charge = int(charge_match.group(1)) if charge_match else 1
                    ion_clean = re.sub(r"\(\d+\+\)", "", ion) if charge_match else ion

                    ion_types.append(ion_clean)
                    charges.append(charge)

                return [charges, ion_types]

            except Exception as e:
                logger.error(f"Parse matches error: {e}")
                return [None, None]

        if self._spectral_data:
            df["ion_mobility"] = None

            # Number of matches -> number_peaks

            df["mz_array"] = (
                df["mz_array"].apply(parse_array_string)
                if "mz_array" in df.columns
                else None
            )
            df["intensity_array"] = (
                df["intensity_array"].apply(parse_array_string)
                if "intensity_array" in df.columns
                else None
            )

            if "Matches" in df.columns:
                df[["charge_array", "ion_type_array"]] = (
                    df["Matches"].apply(parse_matches).apply(pd.Series)
                )
            else:
                df["charge_array"] = None
                df["ion_type_array"] = None

            df["ion_mobility_array"] = None
        else:
            spectra_info_cols = [
                "ion_mobility",
                "number_peaks",
                "mz_array",
                "intensity_array",
                "charge_array",
                "ion_type_array",
                "ion_mobility_array",
            ]
            for col in spectra_info_cols:
                df[col] = None

        return df

    def _ensure_psm_schema_compliance(self, df: pd.DataFrame) -> pd.DataFrame:
        """Ensure PSM data complies with PSM_SCHEMA"""
        schema_fields = [field.name for field in PSM_SCHEMA]

        for field in PSM_SCHEMA:
            if field.name not in df.columns:
                if field.type == pa.string():
                    df[field.name] = ""
                elif field.type == pa.int32():
                    df[field.name] = 0
                elif field.type == pa.float32():
                    df[field.name] = 0.0
                else:
                    df[field.name] = None

        df = df[[col for col in schema_fields if col in df.columns]].copy()

        if "is_decoy" in df.columns:
            df["is_decoy"] = self._convert_maxquant_flag_column(df["is_decoy"])
        if "precursor_charge" in df.columns:
            df["precursor_charge"] = df["precursor_charge"].astype("int32")
        if "calculated_mz" in df.columns:
            df["calculated_mz"] = df["calculated_mz"].astype("float32")
        if "observed_mz" in df.columns:
            df["observed_mz"] = df["observed_mz"].astype("float32")
        if "rt" in df.columns:
            df["rt"] = (df["rt"] * 60).astype("float32")
        if "scan" in df.columns:
            df["scan"] = df["scan"].astype("string")
        if "number_peaks" in df.columns:
            df["number_peaks"] = (
                pd.to_numeric(df["number_peaks"], errors="coerce")
                .fillna(0)
                .astype("int32")
            )

        return df

    # ============================================================================
    # Feature Processing
    # ============================================================================

    def process_feature_file(
        self,
        evidence_path: str,
        output_path: str,
        sdrf_path: str = None,
        protein_file: str = None,
        chunksize: int = 1000000,
    ) -> None:
        """Process Feature data from evidence.txt to Feature parquet format"""
        if sdrf_path:
            self._init_sdrf(sdrf_path)

        evidence_dir = Path(evidence_path).parent
        protein_groups_path = evidence_dir / "proteinGroups.txt"
        if protein_groups_path.exists():
            self._init_protein_group_qvalue_mapping(str(protein_groups_path))

        if protein_file:
            self._init_protein_group_qvalue_mapping(protein_file)

        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None

        batch_writer = ParquetBatchWriter(output_path, FEATURE_SCHEMA)

        try:
            available_cols = pd.read_csv(
                evidence_path, sep="\t", nrows=0
            ).columns.tolist()

            usecols_filtered = [
                col for col in MAXQUANT_FEATURE_USECOLS if col in available_cols
            ]

            cv_columns = ["Type"]  # Add more MaxQuant columns as needed
            cv_cols_available = [col for col in cv_columns if col in available_cols]
            usecols_filtered.extend(cv_cols_available)

            if self.experiment_type and "TMT" in self.experiment_type.upper():
                reporter_cols = [
                    col
                    for col in available_cols
                    if col.startswith("Reporter intensity")
                ]
                usecols_filtered.extend(reporter_cols)

            for df_chunk in pd.read_csv(
                evidence_path,
                sep="\t",
                chunksize=chunksize,
                usecols=usecols_filtered,
                low_memory=False,
            ):

                if protein_str:
                    df_chunk = df_chunk[
                        df_chunk["Proteins"].str.contains(
                            protein_str, na=False, regex=False
                        )
                    ]
                    if df_chunk.empty:
                        continue

                df_chunk = self._process_feature_cv_params(df_chunk)
                df_chunk = self._apply_feature_mapping(df_chunk)
                self._calculate_theoretical_mz_batch(df_chunk)
                df_chunk = self._process_feature_modifications(df_chunk)
                df_chunk = self._process_feature_protein_groups(df_chunk)
                df_chunk = self._process_feature_scores(df_chunk)

                if hasattr(self, "_protein_group_qvalue_map"):
                    df_chunk = self._map_protein_group_qvalue(df_chunk)

                df_chunk = self._process_feature_intensities(df_chunk)

                if self.sdrf_handler:
                    df_chunk = self._integrate_sdrf_metadata_feature(df_chunk)

                df_chunk = self._ensure_feature_schema_compliance(df_chunk)

                table = pa.Table.from_pandas(
                    df_chunk, schema=FEATURE_SCHEMA, preserve_index=False
                )
                batch_writer.write_batch(table.to_pylist())

        finally:
            batch_writer.close()

    def _apply_feature_mapping(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply MAXQUANT_FEATURE_MAP mapping"""
        available_mapping = {
            k: v for k, v in self.feature_mapping.items() if k in df.columns
        }
        df.rename(columns=available_mapping, inplace=True)
        return df

    def _process_feature_modifications(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process Feature modification information"""
        return self._process_psm_modifications(df)

    def _process_feature_protein_groups(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process Feature protein group information"""
        if "pg_accessions" in df.columns:
            df["pg_accessions"] = self._split_semicolon_separated_column(
                df["pg_accessions"]
            )

        if "gg_names" in df.columns:
            df["gg_names"] = (
                df["gg_names"]
                .fillna("")
                .astype(str)
                .apply(
                    lambda x: (
                        [name.strip() for name in x.split(";") if name.strip()]
                        if x
                        else []
                    )
                )
            )
        if "anchor_protein" not in df.columns:
            df["anchor_protein"] = None

        if "pg_accessions" in df.columns:
            df["unique"] = (df["pg_accessions"].apply(len) == 1).astype("int32")
        else:
            df["unique"] = 0

        return df

    def _process_feature_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process Feature scoring information"""
        non_schema_score_fields = {
            "andromeda_score": "andromeda_score",
            "andromeda_delta_score": "andromeda_delta_score",
            "parent_ion_fraction": "parent_ion_fraction",
        }  # Add more additional score fields as needed

        scores_list = []
        for _, row in df.iterrows():
            scores = []
            for field_name, standard_name in non_schema_score_fields.items():
                if field_name in row and pd.notna(row[field_name]):
                    scores.append(
                        {
                            "score_name": standard_name,
                            "score_value": float(row[field_name]),
                        }
                    )
            scores_list.append(scores if scores else None)

        df["additional_scores"] = scores_list

        for field_name in non_schema_score_fields.keys():
            if field_name in df.columns:
                df = df.drop(columns=[field_name])

        return df

    def _process_feature_intensities(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process Feature intensity data"""
        intensities_list = []
        additional_intensities_list = []

        tmt_channels = []
        if self.experiment_type and "TMT" in self.experiment_type.upper():
            tmt_channels = self._get_tmt_channels_from_sdrf()

        for _, row in df.iterrows():
            intensities = []
            additional_intensities = []
            reference_file_name = row.get("reference_file_name", "")

            if tmt_channels:
                tmt_intensities, tmt_additional = self._process_tmt_intensities(
                    row, tmt_channels, reference_file_name
                )
                intensities.extend(tmt_intensities)
                additional_intensities.extend(tmt_additional)

            elif (
                "intensity" in row
                and pd.notna(row["intensity"])
                and row["intensity"] > 0
            ):
                sample_accession = self._get_sample_accession_from_sdrf(
                    reference_file_name
                )
                channel = self._get_channel_from_sdrf(
                    reference_file_name, sample_accession
                )

                intensities.append(
                    {
                        "sample_accession": sample_accession,
                        "channel": channel,
                        "intensity": float(row["intensity"]),
                    }
                )

            intensities_list.append(intensities if intensities else None)
            additional_intensities_list.append(
                additional_intensities if additional_intensities else None
            )

        df["intensities"] = intensities_list
        df["additional_intensities"] = additional_intensities_list

        return df

    def _integrate_sdrf_metadata_feature(self, df: pd.DataFrame) -> pd.DataFrame:
        """Integrate SDRF metadata for Feature"""
        return df

    def _ensure_feature_schema_compliance(self, df: pd.DataFrame) -> pd.DataFrame:
        """Ensure Feature data complies with FEATURE_SCHEMA"""
        df = self._add_missing_feature_schema_fields(df)
        df = self._reorder_feature_columns_by_schema(df)
        df = self._convert_feature_data_types(df)
        df = self._convert_feature_float_fields(df)
        df = self._convert_rt_fields(df)
        return df

    def _add_missing_feature_schema_fields(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add missing schema fields with default values"""
        for field in FEATURE_SCHEMA:
            if field.name not in df.columns:
                df[field.name] = self._get_default_value_for_field(field)
        return df

    def _reorder_feature_columns_by_schema(self, df: pd.DataFrame) -> pd.DataFrame:
        """Reorder columns according to schema"""
        schema_fields = [field.name for field in FEATURE_SCHEMA]
        return df[[col for col in schema_fields if col in df.columns]].copy()

    def _convert_feature_data_types(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert specific feature data types"""
        if "is_decoy" in df.columns:
            df["is_decoy"] = self._convert_maxquant_flag_column(df["is_decoy"])
        if "precursor_charge" in df.columns:
            df["precursor_charge"] = df["precursor_charge"].astype("int32")
        if "scan" in df.columns:
            df["scan"] = df["scan"].astype("string")
        return df

    def _convert_feature_float_fields(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert float fields to proper data types"""
        float_fields = [
            "posterior_error_probability",
            "calculated_mz",
            "observed_mz",
            "rt",
            "predicted_rt",
            "ion_mobility",
            "start_ion_mobility",
            "stop_ion_mobility",
            "pg_global_qvalue",
            "rt_start",
            "rt_stop",
        ]  # Add more float fields as needed

        for field in float_fields:
            if field in df.columns:
                df[field] = pd.to_numeric(df[field], errors="coerce").astype("float32")
        return df

    def _convert_rt_fields(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert retention time fields from minutes to seconds"""
        rt_fields = ["rt", "rt_start", "rt_stop"]
        for field in rt_fields:
            if field in df.columns:
                df[field] = (df[field] * 60).astype("float32")
        return df

    def _get_default_value_for_field(self, field):
        """Get default value for schema field"""
        if field.type == pa.string():
            return ""
        elif field.type == pa.int32():
            return 0
        elif field.type == pa.float32():
            return None if field.nullable else 0.0
        elif str(field.type).startswith("list"):
            return None
        else:
            return None

    def _init_protein_group_qvalue_mapping(self, protein_groups_path: str) -> None:
        """Initialize protein group Q-value mapping"""
        try:
            pg_df = pd.read_csv(
                protein_groups_path,
                sep="\t",
                usecols=["Protein IDs", "Majority protein IDs", "Q-value"],
                low_memory=False,
            )

            self._protein_group_qvalue_map = {}

            for _, row in pg_df.iterrows():
                qvalue = row["Q-value"]

                if pd.notna(row["Majority protein IDs"]):
                    protein_ids = str(row["Majority protein IDs"]).split(";")
                    for protein_id in protein_ids:
                        protein_id = protein_id.strip()
                        if protein_id:
                            self._protein_group_qvalue_map[protein_id] = qvalue

        except Exception as e:
            logging.error(f"Failed to initialize protein group Q-value mapping: {e}")
            self._protein_group_qvalue_map = {}

    def _map_protein_group_qvalue(self, df: pd.DataFrame) -> pd.DataFrame:
        """Map protein group Q-values to feature data"""
        if not self._has_protein_qvalue_mapping():
            df["pg_global_qvalue"] = None
            return df

        if "pg_accessions" in df.columns:
            df["pg_global_qvalue"] = df["pg_accessions"].apply(
                self._get_qvalue_for_proteins
            )
        else:
            df["pg_global_qvalue"] = None

        return df

    def _has_protein_qvalue_mapping(self) -> bool:
        """Check if protein Q-value mapping is available"""
        return (
            hasattr(self, "_protein_group_qvalue_map")
            and self._protein_group_qvalue_map is not None
        )

    def _get_qvalue_for_proteins(self, protein_accessions) -> Optional[float]:
        """Get minimum Q-value for protein list"""
        if self._is_empty_protein_accessions(protein_accessions):
            return None

        protein_list = self._extract_protein_list(protein_accessions)
        if not protein_list:
            return None

        qvalues = self._collect_protein_qvalues(protein_list)
        return min(qvalues) if qvalues else None

    def _is_empty_protein_accessions(self, protein_accessions) -> bool:
        """Check if protein accessions is empty or None"""
        if protein_accessions is None:
            return True

        try:
            if pd.isna(protein_accessions):
                return True
        except (ValueError, TypeError):
            if hasattr(protein_accessions, "__len__") and len(protein_accessions) == 0:
                return True

        return False

    def _extract_protein_list(self, protein_accessions) -> List[str]:
        """Extract protein list from various input formats"""
        if isinstance(protein_accessions, str):
            return protein_accessions.split(";") if protein_accessions.strip() else []
        elif isinstance(protein_accessions, list):
            return protein_accessions
        elif hasattr(protein_accessions, "__iter__"):
            return list(protein_accessions)
        else:
            return []

    def _collect_protein_qvalues(self, protein_list: List[str]) -> List[float]:
        """Collect Q-values for valid protein IDs"""
        qvalues = []
        for protein_id in protein_list:
            if isinstance(protein_id, str):
                protein_id = protein_id.strip()
                if protein_id and protein_id in self._protein_group_qvalue_map:
                    qvalues.append(self._protein_group_qvalue_map[protein_id])
        return qvalues

    # ============================================================================
    # Protein Group Processing
    # ============================================================================

    def process_pg_file(
        self,
        protein_groups_path: str,
        output_path: str,
        sdrf_path: str = None,
        chunksize: int = 100000,
    ) -> None:
        """Process Protein Group data from proteinGroups.txt to PG parquet format"""
        if sdrf_path:
            self._init_sdrf(sdrf_path)

        batch_writer = ParquetBatchWriter(output_path, PG_SCHEMA)

        try:
            for df_chunk in pd.read_csv(
                protein_groups_path, sep="\t", chunksize=chunksize, low_memory=False
            ):

                basic_cols = [
                    col for col in MAXQUANT_PG_USECOLS if col in df_chunk.columns
                ]
                intensity_cols = [
                    col
                    for col in df_chunk.columns
                    if col.startswith("Intensity ")
                    or col.startswith("LFQ intensity ")
                    or col.startswith("iBAQ ")
                ]
                # Add "Majority protein IDs" for anchor_protein processing
                additional_cols = []
                if "Majority protein IDs" in df_chunk.columns:
                    additional_cols.append("Majority protein IDs")

                available_cols = list(
                    set(basic_cols + intensity_cols + additional_cols)
                )
                df_chunk = df_chunk[available_cols]

                if df_chunk.empty:
                    continue

                df_chunk = self._apply_pg_mapping(df_chunk)
                df_chunk = self._process_pg_basic_fields(df_chunk)
                df_chunk = self._process_pg_intensities(df_chunk)
                df_chunk = self._calculate_pg_statistics(df_chunk)

                if self.sdrf_handler:
                    df_chunk = self._integrate_sdrf_metadata_pg(df_chunk)

                df_chunk = self._ensure_pg_schema_compliance(df_chunk)

                table = pa.Table.from_pandas(
                    df_chunk, schema=PG_SCHEMA, preserve_index=False
                )
                batch_writer.write_batch(table.to_pylist())

        finally:
            batch_writer.close()

    def _apply_pg_mapping(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply MAXQUANT_PG_MAP mapping"""
        available_mapping = {
            k: v for k, v in self.pg_mapping.items() if k in df.columns
        }
        df.rename(columns=available_mapping, inplace=True)
        return df

    def _process_pg_basic_fields(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process PG basic fields"""
        list_columns = [
            "pg_accessions",
            "pg_names",
            "gg_accessions",
        ]  # Add more semicolon-separated columns as needed
        for col in list_columns:
            if col in df.columns:
                df[col] = self._split_semicolon_separated_column(df[col])
            else:
                df[col] = [[] for _ in range(len(df))]

        if "is_decoy" in df.columns:
            df["is_decoy"] = self._convert_maxquant_flag_column(df["is_decoy"])
        else:
            df["is_decoy"] = 0

        if "contaminant" in df.columns:
            df["contaminant"] = self._convert_maxquant_flag_column(df["contaminant"])
        else:
            df["contaminant"] = 0

        return df

    def _process_pg_intensities(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process PG intensity data with comprehensive intensity extraction"""
        intensities_list = []
        additional_intensities_list = []

        intensity_cols = self._get_sample_specific_intensity_cols(df)
        lfq_cols = [col for col in df.columns if col.startswith("LFQ intensity ")]
        ibaq_cols = self._get_sample_specific_ibaq_cols(df)

        general_intensity_col = "Intensity" if "Intensity" in df.columns else None
        general_ibaq_col = "iBAQ" if "iBAQ" in df.columns else None

        reference_file_names = []

        for _, row in df.iterrows():
            intensities, additional_intensities, reference_file_name = (
                self._create_pg_intensity_struct(
                    row=row,
                    intensity_cols=intensity_cols,
                    lfq_cols=lfq_cols,
                    ibaq_cols=ibaq_cols,
                    general_intensity_col=general_intensity_col,
                    general_ibaq_col=general_ibaq_col,
                )
            )

            intensities_list.append(intensities if intensities else [])
            additional_intensities_list.append(
                additional_intensities if additional_intensities else []
            )
            reference_file_names.append(reference_file_name)

        df["intensities"] = intensities_list
        df["additional_intensities"] = additional_intensities_list
        df["reference_file_name"] = reference_file_names

        return df

    def _get_sample_specific_intensity_cols(self, df: pd.DataFrame) -> list:
        """Get sample-specific intensity columns (excluding general 'Intensity')"""
        return [
            col
            for col in df.columns
            if col.startswith("Intensity ") and col != "Intensity"
        ]

    def _get_sample_specific_ibaq_cols(self, df: pd.DataFrame) -> list:
        """Get sample-specific iBAQ columns (excluding general 'iBAQ')"""
        return [col for col in df.columns if col.startswith("iBAQ ") and col != "iBAQ"]

    def _create_pg_intensity_struct(
        self,
        row,
        intensity_cols,
        lfq_cols,
        ibaq_cols,
        general_intensity_col,
        general_ibaq_col,
    ) -> tuple:
        """Create intensity structure for a single protein group row
        Use highest intensity sample for PG reference_file_name assignment"""
        intensities = []
        additional_intensities = []
        max_intensity = 0.0
        reference_file_name = "proteinGroups.txt"

        if (
            general_intensity_col
            and general_intensity_col in row.index
            and pd.notna(row[general_intensity_col])
        ):
            sample_accession, channel = self._get_default_sample_info()
            intensity_value = float(row[general_intensity_col])
            intensities.append(
                {
                    "sample_accession": sample_accession,
                    "channel": channel,
                    "intensity": intensity_value,
                }
            )
            if intensity_value > max_intensity:
                max_intensity = intensity_value
                if self.sdrf_handler and hasattr(self.sdrf_handler, "sdrf_table"):
                    reference_file_name = self._get_reference_file_for_sample(
                        sample_accession
                    )

        sample_intensities = self._process_sample_specific_pg_intensities(
            row, intensity_cols
        )
        intensities.extend(sample_intensities)

        for intensity_item in sample_intensities:
            if intensity_item.get("intensity", 0) > max_intensity:
                max_intensity = intensity_item["intensity"]
                sample_acc = intensity_item["sample_accession"]
                reference_file_name = self._get_reference_file_for_sample(sample_acc)

        additional_intensities.extend(self._process_lfq_pg_intensities(row, lfq_cols))

        additional_intensities.extend(
            self._process_ibaq_pg_intensities(row, ibaq_cols, general_ibaq_col)
        )

        return intensities, additional_intensities, reference_file_name

    def _get_reference_file_for_sample(self, sample_accession: str) -> str:
        """Get reference file name for a given sample accession"""
        if not (self.sdrf_handler and hasattr(self.sdrf_handler, "sdrf_table")):
            return "proteinGroups.txt"

        try:
            matching_rows = self.sdrf_handler.sdrf_table[
                self.sdrf_handler.sdrf_table["source name"] == sample_accession
            ]

            if not matching_rows.empty:
                raw_data_file = matching_rows.iloc[0].get("comment[data file]")
                if raw_data_file and isinstance(raw_data_file, str):
                    return re.sub(r"\.[^.]*$", "", raw_data_file)
        except KeyError as e:
            logger.warning(f"Missing SDRF column for sample {sample_accession}: {e}")
        except (AttributeError, IndexError) as e:
            logger.warning(
                f"SDRF data structure issue for sample {sample_accession}: {e}"
            )
        except Exception as e:
            logger.error(
                f"Unexpected error getting reference file for sample {sample_accession}: {e}"
            )

        return "proteinGroups.txt"

    def _preprocess_sdrf_for_optimization(self) -> None:
        """Preprocess SDRF data for optimized matching"""
        if not self.sdrf_handler or not hasattr(self.sdrf_handler, "sdrf_table"):
            self._optimized_sdrf_data = None
            return

        sdrf_df = self.sdrf_handler.sdrf_table

        valuable_columns = self._identify_valuable_sdrf_columns(sdrf_df)
        self._optimized_sdrf_data = []

        for idx, row in sdrf_df.iterrows():
            search_texts = [
                str(row[col])
                for col in valuable_columns
                if col in row and pd.notna(row[col])
            ]
            combined_text = " ".join(search_texts).lower()

            # Tokenize and normalize numbers
            raw_tokens = TOKEN_DELIMITER_PATTERN.split(combined_text)
            normalized_tokens = []
            for token in raw_tokens:
                token = token.strip()
                if token:
                    # Normalize numbers by removing leading zeros
                    if token.isdigit():
                        token = str(int(token))
                    normalized_tokens.append(token)

            search_tokens = set(normalized_tokens)

            self._optimized_sdrf_data.append(
                {
                    "idx": idx,
                    "source_name": row.get("source name"),
                    "comment_data_file": row.get("comment[data file]"),
                    "comment_label": row.get("comment[label]"),
                    "search_tokens": search_tokens,
                }
            )

    def _identify_valuable_sdrf_columns(self, sdrf_df: pd.DataFrame) -> List[str]:
        """Identify valuable SDRF columns for matching"""
        core_columns = [
            "source name",
            "comment[data file]",
            "comment[label]",
            "assay name",
        ]  # Add more essential SDRF columns as needed

        characteristic_columns = [
            col
            for col in sdrf_df.columns
            if "characteristics" in col.lower() or "factor value" in col.lower()
        ]

        other_useful = [
            col
            for col in sdrf_df.columns
            if any(
                keyword in col.lower()
                for keyword in [
                    "organism",
                    "tissue",
                    "cell",
                    "disease",
                    "treatment",
                    "condition",
                    "part",
                ]  # Add more sample metadata keywords as needed
            )
        ]

        valuable_columns = []
        all_candidates = core_columns + characteristic_columns + other_useful

        for col in all_candidates:
            if col in sdrf_df.columns:
                non_null_count = sdrf_df[col].notna().sum()
                if non_null_count > 0:
                    valuable_columns.append(col)

        return list(set(valuable_columns))

    def _tokenize_intensity_suffix(self, intensity_col: str) -> List[str]:
        """Tokenize intensity column suffix for matching"""
        if intensity_col.startswith("Intensity "):
            suffix = intensity_col.replace("Intensity ", "").strip()
        elif intensity_col.startswith("LFQ intensity "):
            suffix = intensity_col.replace("LFQ intensity ", "").strip()
        elif intensity_col.startswith("iBAQ "):
            suffix = intensity_col.replace("iBAQ ", "").strip()
        else:
            suffix = intensity_col

        if not suffix:
            return []

        tokens = TOKEN_DELIMITER_PATTERN.split(suffix)

        cleaned_tokens = []
        for token in tokens:
            token = token.strip().lower()
            if token:
                # Normalize numbers by removing leading zeros
                if token.isdigit():
                    token = str(int(token))
                cleaned_tokens.append(token)

        return cleaned_tokens

    def _get_default_cache_result(self, maxquant_sample: str) -> tuple:
        """Create default cache result for failed matches"""
        return (maxquant_sample, "label free sample", None, 0)

    def _find_best_matching_record(
        self, token_set: set, maxquant_sample: str = None
    ) -> tuple:
        """Find best matching SDRF record using token matching with tissue-based fallback"""
        best_record = None
        best_match_count = 0

        for record in self._optimized_sdrf_data:
            matched_tokens = token_set.intersection(record["search_tokens"])
            match_count = len(matched_tokens)

            if match_count > best_match_count:
                best_match_count = match_count
                best_record = record
            elif match_count == best_match_count and match_count > 0:
                # Prefer exact source_name match
                current_source = str(record.get("source_name", "")).lower()
                best_source = (
                    str(best_record.get("source_name", "")).lower()
                    if best_record
                    else ""
                )
                mq_sample_lower = maxquant_sample.lower() if maxquant_sample else ""

                current_is_exact = mq_sample_lower == current_source
                best_is_exact = mq_sample_lower == best_source

                if current_is_exact and not best_is_exact:
                    best_record = record
                elif not current_is_exact and best_is_exact:
                    pass  # Keep best_record
                elif len(record["comment_data_file"]) < len(
                    best_record["comment_data_file"]
                ):
                    best_record = record

        if (
            best_match_count < 2
            and maxquant_sample
            and self.sdrf_handler
            and hasattr(self.sdrf_handler, "sdrf_table")
        ):
            tissue_record = self._try_tissue_based_matching(maxquant_sample)
            if tissue_record:
                return tissue_record, -1

        return best_record, best_match_count

    def _try_tissue_based_matching(self, maxquant_sample: str) -> dict:
        """Match MaxQuant sample to SDRF using tissue name similarity"""
        try:
            parts = maxquant_sample.split("_")
            if len(parts) == 0:
                return None

            # Try each part as potential tissue name
            tissue_candidates = parts if len(parts) > 1 else [maxquant_sample]

            sdrf_df = self.sdrf_handler.sdrf_table
            if "characteristics[organism part]" not in sdrf_df.columns:
                return None

            best_match_row = None
            best_similarity = 0.0
            best_tissue_candidate = None

            # Try matching each candidate part against SDRF tissues
            for tissue_candidate in tissue_candidates:
                tissue_normalized = self._normalize_tissue_name(tissue_candidate)
                if not tissue_normalized or len(tissue_normalized) < 3:
                    continue

                for idx, row in sdrf_df.iterrows():
                    sdrf_tissue = row.get("characteristics[organism part]", "")
                    if pd.notna(sdrf_tissue):
                        sdrf_tissue_normalized = self._normalize_tissue_name(
                            str(sdrf_tissue)
                        )
                        similarity = self._calculate_tissue_similarity(
                            tissue_normalized, sdrf_tissue_normalized
                        )

                        if similarity > best_similarity:
                            best_similarity = similarity
                            best_match_row = row
                            best_tissue_candidate = tissue_candidate

            if best_match_row is not None and best_similarity > 0.5:
                first_sdrf_sample = best_match_row["source name"]
                virtual_sample_name = f"{first_sdrf_sample}_{maxquant_sample}"

                virtual_record = {
                    "source_name": virtual_sample_name,
                    "comment_label": best_match_row.get(
                        "comment[label]", "label free sample"
                    ),
                    "comment_data_file": best_match_row.get("comment[data file]", ""),
                    "search_tokens": set(),
                }

                return virtual_record
        except Exception:
            pass

        return None

    def _calculate_tissue_similarity(self, tissue1: str, tissue2: str) -> float:
        """Calculate similarity between tissue names"""
        if not tissue1 or not tissue2:
            return 0.0

        if tissue1 == tissue2:
            return 1.0

        if tissue1 in tissue2 or tissue2 in tissue1:
            return 0.9

        words1 = set(tissue1.split())
        words2 = set(tissue2.split())

        if len(words1) == 0 or len(words2) == 0:
            return 0.0

        intersection = len(words1.intersection(words2))
        union = len(words1.union(words2))

        if union > 0:
            jaccard_score = intersection / union
            if jaccard_score > 0:
                return jaccard_score

        for w1 in words1:
            for w2 in words2:
                if len(w1) >= 4 and len(w2) >= 4:
                    if w1 in w2 or w2 in w1:
                        return 0.7

        return 0.0

    def _normalize_tissue_name(self, tissue: str) -> str:
        """Normalize tissue name to lowercase with spaces"""
        if not tissue:
            return ""

        tissue = tissue.lower().strip()
        tissue = tissue.replace("_", " ").replace("-", " ")
        tissue = " ".join(tissue.split())

        return tissue

    def _create_result_tuple(self, record: dict, match_count: int) -> tuple:
        """Create result tuple from matching record"""
        sample_accession = record["source_name"]
        channel = record["comment_label"]

        raw_data_file = record["comment_data_file"]
        if raw_data_file and isinstance(raw_data_file, str):
            reference_file_name = re.sub(r"\.[^.]*$", "", raw_data_file)
        else:
            reference_file_name = raw_data_file

        return (sample_accession, channel, reference_file_name, match_count)

    def _get_sample_accession_by_maxquant_name(self, maxquant_sample: str) -> str:
        """Get sample accession using optimized token-based matching"""
        if maxquant_sample in self._sdrf_search_cache:
            return self._sdrf_search_cache[maxquant_sample][0]

        if not self._optimized_sdrf_data:
            result_tuple = self._get_default_cache_result(maxquant_sample)
            self._sdrf_search_cache[maxquant_sample] = result_tuple
            return result_tuple[0]

        try:
            tokens = self._tokenize_intensity_suffix(f"Intensity {maxquant_sample}")

            if not tokens:
                result_tuple = self._get_default_cache_result(maxquant_sample)
                self._sdrf_search_cache[maxquant_sample] = result_tuple
                return result_tuple[0]

            best_record, best_match_count = self._find_best_matching_record(
                set(tokens), maxquant_sample
            )

            if best_record and (best_match_count > 0 or best_match_count == -1):
                result_tuple = self._create_result_tuple(best_record, best_match_count)
            else:
                result_tuple = (maxquant_sample, None, None, 0)

            self._sdrf_search_cache[maxquant_sample] = result_tuple
            return result_tuple[0]

        except Exception:
            result_tuple = self._get_default_cache_result(maxquant_sample)
            self._sdrf_search_cache[maxquant_sample] = result_tuple
            return result_tuple[0]

    def _get_channel_by_maxquant_name(self, maxquant_sample: str) -> str:
        """Get channel using optimized token-based matching"""
        if maxquant_sample in self._sdrf_search_cache:
            cached_result = self._sdrf_search_cache[maxquant_sample]
            return cached_result[1]

        self._get_sample_accession_by_maxquant_name(maxquant_sample)

        if maxquant_sample in self._sdrf_search_cache:
            cached_result = self._sdrf_search_cache[maxquant_sample]
            return cached_result[1]

        return None

    def _process_sample_specific_pg_intensities(self, row, intensity_cols) -> list:
        """Process sample-specific intensity columns with SDRF mapping"""
        sample_intensity_map = {}

        for col in intensity_cols:
            if col in row.index and pd.notna(row[col]) and row[col] > 0:
                maxquant_sample = col.replace("Intensity ", "")
                sample_accession = self._get_sample_accession_by_maxquant_name(
                    maxquant_sample
                )
                channel = self._get_channel_by_maxquant_name(maxquant_sample)

                sample_key = (sample_accession, channel)
                if sample_key not in sample_intensity_map:
                    sample_intensity_map[sample_key] = 0
                sample_intensity_map[sample_key] += float(row[col])

        intensities = []
        for (
            sample_accession,
            channel,
        ), total_intensity in sample_intensity_map.items():
            intensities.append(
                {
                    "sample_accession": sample_accession,
                    "channel": channel,
                    "intensity": total_intensity,
                }
            )
        return intensities

    def _create_additional_intensity_item(
        self, col: str, value: float, sample_accession: str, channel: str
    ) -> dict:
        """Create a single additional intensity item"""
        return {
            "sample_accession": sample_accession,
            "channel": channel,
            "intensities": [
                {
                    "intensity_name": col,
                    "intensity_value": value,
                }
            ],
        }

    def _process_lfq_pg_intensities(self, row, lfq_cols) -> list:
        """Process LFQ intensity columns for PG data"""
        additional_intensities = []
        for col in lfq_cols:
            if col in row.index and pd.notna(row[col]):
                maxquant_sample = col.replace("LFQ intensity ", "")
                sample_accession = self._get_sample_accession_by_maxquant_name(
                    maxquant_sample
                )
                channel = self._get_channel_by_maxquant_name(maxquant_sample)
                additional_intensities.append(
                    self._create_additional_intensity_item(
                        col, float(row[col]), sample_accession, channel
                    )
                )
        return additional_intensities

    def _get_default_sample_info(self) -> tuple:
        """Get default sample accession and channel from SDRF"""
        sample_accession = "Unknown"
        channel = "label free sample"
        if (
            hasattr(self, "_sdrf_transformed")
            and self._sdrf_transformed is not None
            and not self._sdrf_transformed.empty
        ):
            try:
                if "source name" in self._sdrf_transformed.columns:
                    sample_accession = str(
                        self._sdrf_transformed["source name"].iloc[0]
                    )
                if "comment[label]" in self._sdrf_transformed.columns:
                    label = str(self._sdrf_transformed["comment[label]"].iloc[0])
                    channel = label if label != "nan" else "label free sample"
            except (IndexError, KeyError) as e:
                logger.debug(
                    f"Using default sample info due to incomplete SDRF data: {e}"
                )
            except Exception as e:
                logger.warning(f"Error extracting default sample info from SDRF: {e}")
        return sample_accession, channel

    def _process_ibaq_pg_intensities(self, row, ibaq_cols, general_ibaq_col) -> list:
        """Process iBAQ intensity columns for PG data"""
        additional_intensities = []

        if (
            general_ibaq_col
            and general_ibaq_col in row.index
            and pd.notna(row[general_ibaq_col])
        ):
            sample_accession, channel = self._get_default_sample_info()
            additional_intensities.append(
                self._create_additional_intensity_item(
                    general_ibaq_col,
                    float(row[general_ibaq_col]),
                    sample_accession,
                    channel,
                )
            )

        for col in ibaq_cols:
            if col in row.index and pd.notna(row[col]):
                maxquant_sample = col.replace("iBAQ ", "")
                sample_accession = self._get_sample_accession_by_maxquant_name(
                    maxquant_sample
                )
                channel = self._get_channel_by_maxquant_name(maxquant_sample)
                additional_intensities.append(
                    self._create_additional_intensity_item(
                        col, float(row[col]), sample_accession, channel
                    )
                )
        return additional_intensities

    def _get_peptide_count_from_row(self, row) -> int:
        """Extract peptide count from MaxQuant data, trying multiple columns in priority order"""
        peptide_columns = [
            "peptide_count_total",
            "peptide_count_razor_unique",
            "peptide_count_unique",
        ]  # Add more peptide count columns in priority order

        for col in peptide_columns:
            if col in row.index and pd.notna(row[col]) and row[col] > 0:
                return int(row[col])

        return 1

    def _calculate_pg_statistics(self, df: pd.DataFrame) -> pd.DataFrame:
        """Calculate PG statistics by orchestrating sub-processes"""
        df = self._process_pg_peptides(df)
        df = self._process_pg_anchor_protein(df)
        df = self._process_pg_counts(df)
        df = self._process_pg_additional_scores(df)
        return df

    def _process_pg_peptides(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process peptides field for PG data"""
        peptides_list = []
        for _, row in df.iterrows():
            peptides = []
            if "pg_accessions" in row and row["pg_accessions"]:
                total_peptides = self._get_peptide_count_from_row(row)
                protein_count = len(row["pg_accessions"])

                if protein_count == 1:
                    protein = row["pg_accessions"][0]
                    if protein:
                        peptides.append(
                            {"protein_name": protein, "peptide_count": total_peptides}
                        )
                else:
                    peptides = self._distribute_peptides_among_proteins(
                        row["pg_accessions"], total_peptides
                    )

            peptides_list.append(peptides)

        df["peptides"] = peptides_list
        return df

    def _distribute_peptides_among_proteins(
        self, proteins: list, total_peptides: int
    ) -> list:
        """Distribute peptides among multiple proteins in a group"""
        proteins = [p for p in proteins if p]
        if not proteins:
            return []

        peptides = []
        protein_count = len(proteins)

        main_peptide_count = max(1, int(total_peptides * 0.6))
        remaining_peptides = total_peptides - main_peptide_count
        other_count = protein_count - 1

        for i, protein in enumerate(proteins):
            if i == 0:
                peptide_count = main_peptide_count
            else:
                if other_count > 0 and remaining_peptides > 0:
                    base_count = remaining_peptides // other_count
                    if i <= remaining_peptides % other_count:
                        peptide_count = base_count + 1
                    else:
                        peptide_count = base_count
                else:
                    peptide_count = 0

            if peptide_count > 0:
                peptides.append(
                    {
                        "protein_name": protein,
                        "peptide_count": peptide_count,
                    }
                )

        return peptides

    def _process_pg_anchor_protein(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process anchor_protein field for PG data"""
        # first protein from "Majority protein IDs" = anchor_protein
        if len(df) > 0:
            if "Majority protein IDs" in df.columns:
                majority_proteins = self._split_semicolon_separated_column(
                    df["Majority protein IDs"]
                )
                df["anchor_protein"] = majority_proteins.str[0].replace("", None)
            elif "pg_accessions" in df.columns:
                anchor_proteins = df["pg_accessions"].str[0]
                df["anchor_protein"] = anchor_proteins.where(
                    anchor_proteins.notna(), None
                )
            else:
                df["anchor_protein"] = None
        else:
            df["anchor_protein"] = None
        return df

    def _process_pg_counts(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process peptide_counts and feature_counts fields for PG data"""
        peptide_counts_list = []
        feature_counts_list = []

        for _, row in df.iterrows():
            unique_sequences = int(row.get("peptide_count_unique", 1))
            total_sequences = int(row.get("peptide_count_total", 1))

            peptide_counts_list.append(
                {
                    "unique_sequences": unique_sequences,
                    "total_sequences": total_sequences,
                }
            )

            feature_counts_list.append(
                {
                    "unique_features": unique_sequences,
                    "total_features": total_sequences,
                }
            )

        df["peptide_counts"] = peptide_counts_list
        df["feature_counts"] = feature_counts_list
        return df

    def _process_pg_additional_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """Extract additional scores for PG data"""
        non_schema_score_fields = {
            "andromeda_score": "andromeda_score",
            "sequence_coverage": "sequence_coverage",
            "molecular_weight": "molecular_weight",
            "msms_count": "msms_count",
            "number_of_proteins": "number_of_proteins",
            "peptide_count_total": "peptide_count_total",
            "peptide_count_razor_unique": "peptide_count_razor_unique",
            "peptide_count_unique": "peptide_count_unique",
        }  # Add more additional score fields as needed

        scores_list = []
        for _, row in df.iterrows():
            scores = []
            for field_name, standard_name in non_schema_score_fields.items():
                if field_name in row and pd.notna(row[field_name]):
                    scores.append(
                        {
                            "score_name": standard_name,
                            "score_value": float(row[field_name]),
                        }
                    )

            scores_list.append(scores if scores else None)

        df["additional_scores"] = scores_list

        for field_name in non_schema_score_fields.keys():
            if field_name in df.columns:
                df = df.drop(columns=[field_name])

        return df

    def _integrate_sdrf_metadata_pg(self, df: pd.DataFrame) -> pd.DataFrame:
        """Integrate SDRF metadata for PG"""
        return df

    def _ensure_pg_schema_compliance(self, df: pd.DataFrame) -> pd.DataFrame:
        """Ensure PG data complies with PG_SCHEMA"""
        schema_fields = [field.name for field in PG_SCHEMA]

        for field in PG_SCHEMA:
            if field.name not in df.columns:
                if field.type == pa.string():
                    df[field.name] = None
                elif field.type == pa.int32():
                    df[field.name] = 0
                elif field.type == pa.float32():
                    df[field.name] = 1.0
                elif str(field.type).startswith("list"):
                    df[field.name] = [[] for _ in range(len(df))]
                elif str(field.type).startswith("struct"):
                    if field.name == "peptide_counts":
                        df[field.name] = [
                            {"unique_sequences": 0, "total_sequences": 0}
                            for _ in range(len(df))
                        ]
                    elif field.name == "feature_counts":
                        df[field.name] = [
                            {"unique_features": 0, "total_features": 0}
                            for _ in range(len(df))
                        ]
                else:
                    df[field.name] = None

        return df[[col for col in schema_fields if col in df.columns]]

    # ============================================================================
    # Backward Compatibility
    # ============================================================================

    def write_psm_to_file(
        self, msms_path: str, output_path: str, chunksize: int = 1000000
    ) -> None:
        """Backward compatible PSM writing interface"""
        self.process_psm_file(msms_path, output_path, chunksize=chunksize)

    def write_feature_to_file(
        self,
        evidence_path: str,
        sdrf_path: str,
        output_path: str,
        protein_file: str = None,
        chunksize: int = 1000000,
    ) -> None:
        """Backward compatible Feature writing interface"""
        self.process_feature_file(
            evidence_path, output_path, sdrf_path, protein_file, chunksize
        )

    def write_protein_groups_to_file(
        self,
        protein_groups_path: str,
        sdrf_path: str,
        output_path: str,
        chunksize: int = 100000,
    ) -> None:
        """Backward compatible Protein Groups writing interface"""
        self.process_pg_file(protein_groups_path, output_path, sdrf_path, chunksize)

    # ============================================================================
    # Additional Backward Compatibility Methods for Tests
    # ============================================================================

    def read_msms(self, msms_path: str) -> pd.DataFrame:
        """Read and process msms.txt file, returns DataFrame"""
        return pd.read_csv(msms_path, sep="\t", low_memory=False)

    def read_evidence(self, evidence_path: str) -> pd.DataFrame:
        """Read and process evidence.txt file, returns DataFrame"""
        return pd.read_csv(evidence_path, sep="\t", low_memory=False)

    def read_protein_groups(self, protein_groups_path: str) -> pd.DataFrame:
        """Read and process proteinGroups.txt file, returns DataFrame"""
        return pd.read_csv(protein_groups_path, sep="\t", low_memory=False)

    def process_msms_to_psm_table(self, df: pd.DataFrame) -> pa.Table:
        """Process MSMS DataFrame to PSM table"""
        if df.empty:
            raise ValueError("Input DataFrame is empty")
        df_processed = self._apply_psm_mapping(df.copy())
        self._calculate_theoretical_mz_batch(df_processed)
        df_processed = self._process_psm_modifications(df_processed)
        df_processed = self._process_psm_scores(df_processed)
        df_processed = self._process_psm_cv_params(df_processed)
        df_processed = self._process_psm_arrays(df_processed)
        df_processed = self._ensure_psm_schema_compliance(df_processed)
        return pa.Table.from_pandas(
            df_processed, schema=PSM_SCHEMA, preserve_index=False
        )

    def process_evidence_to_feature_table(self, df: pd.DataFrame) -> pa.Table:
        """Process Evidence DataFrame to Feature table"""
        df_processed = df.copy()
        df_processed = self._process_feature_cv_params(df_processed)
        df_processed = self._apply_feature_mapping(df_processed)
        self._calculate_theoretical_mz_batch(df_processed)
        df_processed = self._process_feature_modifications(df_processed)
        df_processed = self._process_feature_protein_groups(df_processed)
        df_processed = self._process_feature_scores(df_processed)
        df_processed = self._ensure_feature_schema_compliance(df_processed)
        return pa.Table.from_pandas(
            df_processed, schema=FEATURE_SCHEMA, preserve_index=False
        )

    def process_protein_groups_to_pg_table(
        self, df: pd.DataFrame, sdrf_path: str = None
    ) -> pa.Table:
        """Process Protein Groups DataFrame to PG table"""
        if df.empty:
            raise ValueError("Input DataFrame is empty")
        if sdrf_path:
            self._init_sdrf(sdrf_path)
        df_processed = self._apply_pg_mapping(df.copy())
        df_processed = self._process_pg_basic_fields(df_processed)
        df_processed = self._calculate_pg_statistics(df_processed)
        df_processed = self._ensure_pg_schema_compliance(df_processed)
        return pa.Table.from_pandas(
            df_processed, schema=PG_SCHEMA, preserve_index=False
        )

    def iter_batch(self, file_path: str, chunksize: int = 10000):
        """Iterate over file in batches"""
        for chunk in pd.read_csv(
            file_path, sep="\t", chunksize=chunksize, low_memory=False
        ):
            yield chunk


# ============================================================================
# Standalone Functions
# ============================================================================


def process_evidence_to_feature_table(df: pd.DataFrame) -> pa.Table:
    """Backward compatible function"""
    processor = MaxQuant()
    df_processed = processor._apply_feature_mapping(df)
    processor._calculate_theoretical_mz_batch(df_processed)
    df_processed = processor._process_feature_modifications(df_processed)
    df_processed = processor._process_feature_protein_groups(df_processed)
    df_processed = processor._process_feature_scores(df_processed)
    df_processed = processor._ensure_feature_schema_compliance(df_processed)
    return pa.Table.from_pandas(
        df_processed, schema=FEATURE_SCHEMA, preserve_index=False
    )


def process_protein_groups_to_pg_table(
    df: pd.DataFrame, sdrf_path: str = None
) -> pa.Table:
    """Backward compatible function"""
    processor = MaxQuant()
    if sdrf_path:
        processor._init_sdrf(sdrf_path)
    df_processed = processor._apply_pg_mapping(df)
    df_processed = processor._process_pg_basic_fields(df_processed)
    df_processed = processor._calculate_pg_statistics(df_processed)
    df_processed = processor._ensure_pg_schema_compliance(df_processed)
    return pa.Table.from_pandas(df_processed, schema=PG_SCHEMA, preserve_index=False)


def process_msms_to_psm_table(df: pd.DataFrame) -> pa.Table:
    """Backward compatible function"""
    processor = MaxQuant()
    df_processed = processor._apply_psm_mapping(df)
    processor._calculate_theoretical_mz_batch(df_processed)
    df_processed = processor._process_psm_modifications(df_processed)
    df_processed = processor._process_psm_scores(df_processed)
    df_processed = processor._process_psm_cv_params(df_processed)
    df_processed = processor._process_psm_arrays(df_processed)
    df_processed = processor._ensure_psm_schema_compliance(df_processed)
    return pa.Table.from_pandas(df_processed, schema=PSM_SCHEMA, preserve_index=False)


# ============================================================================
# Utility Functions
# ============================================================================


def read_msms(msms_path: str) -> pd.DataFrame:
    """Read and process msms.txt file, returns DataFrame"""
    return pd.read_csv(msms_path, sep="\t", low_memory=False)


def read_evidence(evidence_path: str) -> pd.DataFrame:
    """Read and process evidence.txt file, returns DataFrame"""
    return pd.read_csv(evidence_path, sep="\t", low_memory=False)


def read_protein_groups(protein_groups_path: str) -> pd.DataFrame:
    """Read and process proteinGroups.txt file, returns DataFrame"""
    return pd.read_csv(protein_groups_path, sep="\t", low_memory=False)
