import logging
import re
import zipfile
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from typing import List, Union
from pathlib import Path
from pyopenms import ModificationsDB
from pyopenms import AASequence
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.operate.tools import (
    get_ahocorasick,
    get_modification_details,
    get_protein_accession,
)
from pyopenms.Constants import PROTON_MASS_U
from quantmsio.utils.constants import ITRAQ_CHANNEL, TMT_CHANNELS
from quantmsio.core.common import (
    MAXQUANT_PSM_MAP,
    MAXQUANT_PSM_USECOLS,
    MAXQUANT_FEATURE_MAP,
    MAXQUANT_FEATURE_USECOLS,
)
from quantmsio.core.feature import Feature
from quantmsio.core.psm import Psm
from quantmsio.utils.file_utils import close_file, extract_protein_list, save_slice_file

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)

intensity_normalize_pattern = r"Reporter intensity corrected \d+"
intensity_pattern = r"Reporter intensity \d+"
MOD_PATTERN = r"\((.*?)\)"


def check_acronym(df):
    s = df["peptidoform"]
    for index in s.index:
        seq = s[index]
        if "(" in seq:
            match = re.search(MOD_PATTERN, seq)
            if len(match.group(1)) == 2:
                return True
            else:
                return False
    return False


class MaxQuant:
    def __init__(self):
        pass

    def extract_col_msg(self, col_df, label: str = "feature"):
        line = "\t".join(col_df.columns)
        self.mods_map = self.get_mods_map(line)
        self._automaton = get_ahocorasick(self.mods_map)
        if label == "feature":
            intensity_normalize_names = []
            intensity_names = []
            use_cols = MAXQUANT_FEATURE_USECOLS.copy()
            use_map = MAXQUANT_FEATURE_MAP.copy()
            for col in col_df.columns:
                if re.search(intensity_normalize_pattern, col, re.IGNORECASE):
                    use_cols.append(col)
                    intensity_normalize_names.append(col)
                elif re.search(intensity_pattern, col, re.IGNORECASE):
                    use_cols.append(col)
                    intensity_names.append(col)
            self._intensity_normalize_names = intensity_normalize_names
            self._intensity_names = intensity_names
            use_cols.append("Intensity")
        else:
            use_cols = MAXQUANT_PSM_USECOLS.copy()
            use_map = MAXQUANT_PSM_MAP.copy()
        for key in self.mods_map.keys():
            col = f"{key} Probabilities"
            if col in col_df.columns:
                use_cols.append(col)
        return use_map, use_cols

    def iter_batch(
        self,
        file_path: Union[Path, str],
        label: str = "feature",
        chunksize: int = 100000,
        protein_str: str = None,
    ):
        col_df = pd.read_csv(file_path, sep="\t", nrows=0)
        use_map, use_cols = self.extract_col_msg(col_df, label=label)
        for df in pd.read_csv(
            file_path,
            sep="\t",
            usecols=use_cols,
            low_memory=False,
            chunksize=chunksize,
        ):
            df.rename(columns=use_map, inplace=True)
            if protein_str:
                df = df[df["mp_accessions"].str.contains(f"{protein_str}", na=False)]
            df = self.main_operate(df)
            yield df

    @staticmethod
    def open_from_zip_archive(zip_file, file_name, **kwargs):
        """Open file from zip archive."""
        with zipfile.ZipFile(zip_file) as z:
            with z.open(file_name) as f:
                df = pd.read_csv(f, sep="\t", low_memory=False, **kwargs)
        return df

    def read_zip_file(self, zip_path: str, **kwargs):
        filepath = Path(zip_path)
        df = self.open_from_zip_archive(
            zip_path, f"{filepath.stem}/evidence.txt", **kwargs
        )
        return df

    def iter_zip_batch(
        self, zip_list: List[str], label: str = "feature", protein_str: str = None
    ):
        for zip_file in zip_list:
            col_df = self.read_zip_file(zip_file, nrows=0)
            use_map, use_cols = self.extract_col_msg(col_df, label=label)
            df = self.read_zip_file(
                zip_file,
                usecols=use_cols,
            )
            df.rename(columns=use_map, inplace=True)
            # df["reference_file_name"] = zip_file.split(".")[0]
            if protein_str:
                df = df[df["mp_accessions"].str.contains(f"{protein_str}", na=False)]
            df = self.main_operate(df)
            yield df

    @staticmethod
    def generete_calculated_mz(df):
        uniq_p = df["peptidoform"].unique()
        masses_map = {k: AASequence.fromString(k).getMonoWeight() for k in uniq_p}
        mass_vector = df["peptidoform"].map(masses_map)
        df.loc[:, "calculated_mz"] = (
            mass_vector + (PROTON_MASS_U * df["precursor_charge"].values)
        ) / df["precursor_charge"].values

    def _transform_mod(self, match):
        if not match:
            return None
        else:
            key = match.group(1)
            if key in self.mods_map:
                return f"({self.mods_map[key]})"
            else:
                return None

    @staticmethod
    def get_mods_map(line):
        pattern = r"sequence\s+(.*?)\s+(Missed|Proteins)"
        match = re.search(pattern, line, re.DOTALL)
        mod_seq = match.group(1)
        mods_map = {}
        modifications_db = ModificationsDB()
        if mod_seq:
            for current_mod in mod_seq.split("\t"):
                if (
                    "Probabilities" not in current_mod
                    and "diffs" not in current_mod
                    and "Diffs" not in current_mod
                ):
                    name = re.search(r"([^ ]+)\s?", current_mod)
                    mod = modifications_db.getModification(name.group(1))
                    unimod = mod.getUniModAccession()
                    match = re.search(MOD_PATTERN, current_mod)
                    if match:
                        site = match.group(1)
                    else:
                        site = "X"
                    mods_map[current_mod] = [unimod.upper(), site]
                    mods_map[current_mod[:2].lower()] = current_mod
                    mods_map[unimod.upper()] = [current_mod, site]
        return mods_map

    def generate_modification_details(self, df):
        keys = {}
        pattern = r"\((\d+\.?\d*)\)"
        for key in self.mods_map.keys():
            col = f"{key} Probabilities"
            if col in df.columns:
                keys[key] = col
        other_mods = list(set(self.mods_map.keys()) - set(keys.keys()))

        def get_details(rows):
            modification_details = []
            for key, col in keys.items():
                modification_map = {"modification_name": self.mods_map[key][0]}
                details = []
                seq = rows[col]
                if not isinstance(seq, str):
                    continue
                match_obj = re.search(pattern, seq)
                while match_obj:
                    details.append(
                        {
                            "position": match_obj.start(0),
                            "localization_probability": float(match_obj.group(1)),
                        }
                    )
                    seq = seq.replace(match_obj.group(0), "", 1)
                    match_obj = re.search(pattern, seq)
                modification_map["fields"] = details
                modification_details.append(modification_map)
            seq = rows["peptidoform"]
            peptidoform, other_modification_details = get_modification_details(
                seq, self.mods_map, self._automaton, other_mods
            )
            modification_details = modification_details + other_modification_details
            if len(modification_details) == 0:
                return [peptidoform, None]
            else:
                return [peptidoform, modification_details]

        df[["peptidoform", "modifications"]] = df[
            list(keys.values()) + ["peptidoform"]
        ].apply(get_details, axis=1, result_type="expand")

    def generete_peptidoform(self, df):
        isacronym = check_acronym(df)

        def trasform_peptidoform(row):
            row = row.replace("_", "")
            if isacronym:
                return re.sub(MOD_PATTERN, self._transform_mod, row)
            else:
                return row

        df["peptidoform"] = df["peptidoform"].apply(trasform_peptidoform)

    def generate_intensity_msg(self, df, intensity_cols, additions_intensity_cols):
        def get_intensities_map(rows):
            result_intensity = []
            result_additions_intensity = []
            reference = rows["reference_file_name"]
            for col in intensity_cols:
                channel = channel_map[col]
                sample_key = reference + "-" + channel
                result_intensity.append(
                    {
                        "sample_accession": self._sample_map[sample_key],
                        "channel": channel,
                        "intensity": rows[col],
                    }
                )
            for col in additions_intensity_cols:
                channel = channel_map[col]
                sample_key = reference + "-" + channel
                result_additions_intensity.append(
                    {
                        "sample_accession": self._sample_map[sample_key],
                        "channel": channel,
                        "additional_intensity": [
                            {
                                "intensity_name": "normalize_intensity",
                                "intensity_value": rows[col],
                            }
                        ],
                    }
                )
            if len(result_intensity) > 1 and len(result_additions_intensity) > 1:
                return result_intensity, result_additions_intensity
            elif len(result_intensity) > 1:
                return result_intensity, None
            elif len(result_additions_intensity) > 1:
                return None, result_additions_intensity
            else:
                return None, None

        channel_map = {}
        if self.experiment_type != "LFQ":
            for col, col1 in zip(intensity_cols, additions_intensity_cols):
                key = re.search(r"\d+", col).group()
                if "TMT" in self.experiment_type:
                    c = TMT_CHANNELS[self.experiment_type][int(key)]
                    channel_map[col] = c
                    channel_map[col1] = c

                else:
                    c = ITRAQ_CHANNEL[self.experiment_type][int(key)]
                    channel_map[col] = c
                    channel_map[col1] = c
            df[["intensities", "additional_intensities"]] = df[
                ["reference_file_name"] + intensity_cols + additions_intensity_cols
            ].apply(get_intensities_map, axis=1, result_type="expand")
        else:
            df.loc[:, "intensities"] = df[["reference_file_name", "Intensity"]].apply(
                lambda rows: [
                    {
                        "sample_accession": self._sample_map[
                            rows["reference_file_name"] + "-" + "LFQ"
                        ],
                        "channel": "LFQ",
                        "intensity": rows["Intensity"],
                    }
                ],
                axis=1,
            )
            df.loc[:, "additional_intensities"] = None

    def main_operate(self, df: pd.DataFrame):
        self.generete_peptidoform(df)
        self.generete_calculated_mz(df)
        self.generate_modification_details(df)
        df = df[df["posterior_error_probability"] < 0.05].copy()
        df["is_decoy"] = df["is_decoy"].map({None: "0", np.nan: "0", "+": "1"})
        df["additional_scores"] = df[
            ["andromeda_score", "andromeda_delta_score"]
        ].apply(
            lambda row: [
                {
                    "score_name": "andromeda_score",
                    "score_value": row["andromeda_score"],
                },
                {
                    "score_name": "andromeda_delta_score",
                    "score_value": row["andromeda_delta_score"],
                },
            ],
            axis=1,
        )
        df.loc[:, "cv_params"] = df["parent_ion_fraction"].apply(
            lambda socre: [{"cv_name": "parent_ion_fraction", "cv_value": str(socre)}]
        )
        df.loc[:, "predicted_rt"] = None
        df.loc[:, "ion_mobility"] = None
        return df

    @staticmethod
    def transform_psm(df: pd.DataFrame):
        df.loc[:, "intensity_array"] = None
        df.loc[:, "mz_array"] = None
        df.loc[:, "number_peaks"] = None
        df.loc[:, "charge_array"] = None
        df.loc[:, "ion_type_array"] = None
        df.loc[:, "ion_mobility_array"] = None

    def transform_feature(self, df: pd.DataFrame):
        self.generate_intensity_msg(
            df, self._intensity_names, self._intensity_normalize_names
        )
        df.loc[:, "unique"] = df["pg_accessions"].apply(
            lambda x: 0 if ";" in x or "," in x else 1
        )
        df["pg_accessions"] = df["pg_accessions"].apply(get_protein_accession)
        df["mp_accessions"] = df["mp_accessions"].apply(get_protein_accession)
        df["gg_names"] = df["gg_names"].str.split(";")
        df.loc[:, "anchor_protein"] = df["pg_accessions"].str[0]
        df.loc[:, "gg_accessions"] = None
        df.loc[:, "pg_global_qvalue"] = None
        df.loc[:, "scan_reference_file_name"] = None
        df.loc[:, "start_ion_mobility"] = None
        df.loc[:, "stop_ion_mobility"] = None

    def write_psm_to_file(
        self, msms_path: str, output_path: str, chunksize: int = 1000000
    ):
        pqwriter = None
        for df in self.iter_batch(msms_path, "psm", chunksize=chunksize):
            self.transform_psm(df)
            Psm.convert_to_parquet_format(df)
            parquet = Psm.transform_parquet(df)
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, parquet.schema)
            pqwriter.write_table(parquet)
        close_file(pqwriter=pqwriter)

    def _init_sdrf(self, sdrf_path: Union[Path, str]):
        sdrf = SDRFHandler(sdrf_path)
        self.experiment_type = sdrf.get_experiment_type_from_sdrf()
        self._sample_map = sdrf.get_sample_map_run()

    def write_feature_to_file(
        self,
        evidence_path: str,
        sdrf_path: str,
        output_path: str,
        chunksize: int = 1000000,
        protein_file=None,
    ):
        self._init_sdrf(sdrf_path)
        pqwriter = None
        for df in self.iter_batch(
            evidence_path, chunksize=chunksize, protein_str=protein_file
        ):
            self.transform_feature(df)
            Feature.convert_to_parquet_format(df)
            parquet = Feature.transform_feature(df)
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, parquet.schema)
            pqwriter.write_table(parquet)
        close_file(pqwriter=pqwriter)

    def write_zip_feature_to_file(
        self,
        zip_list,
        sdrf_path: str,
        output_path: str,
        protein_file=None,
    ):
        self._init_sdrf(sdrf_path)
        pqwriter = None
        for df in self.iter_zip_batch(zip_list, "feature", protein_str=protein_file):
            self.transform_feature(df)
            Feature.convert_to_parquet_format(df)
            parquet = Feature.transform_feature(df)
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, parquet.schema)
            pqwriter.write_table(parquet)
        close_file(pqwriter=pqwriter)

    def write_features_to_file(
        self,
        evidence_path: str,
        sdrf_path: str,
        output_folder: str,
        filename: str,
        partitions: list,
        chunksize: int = 1000000,
        protein_file=None,
    ):
        pqwriters = {}
        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None
        self._init_sdrf(sdrf_path)
        for report in self.iter_batch(
            evidence_path, chunksize=chunksize, protein_str=protein_str
        ):
            self.transform_feature(report)
            Feature.convert_to_parquet_format(report)
            for key, df in Feature.slice(report, partitions):
                feature = Feature.transform_feature(df)
                pqwriters = save_slice_file(
                    feature, pqwriters, output_folder, key, filename
                )
        close_file(pqwriters=pqwriters)
