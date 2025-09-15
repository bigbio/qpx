import re
import os
from pathlib import Path
from typing import Union, Optional

import pyarrow as pa
import pyarrow.parquet as pq
from quantmsio.utils.file_utils import extract_protein_list
from quantmsio.utils.pride_utils import (
    get_petidoform_msstats_notation,
    generate_scan_number,
)
from quantmsio.operate.tools import get_ahocorasick, get_protein_accession
from quantmsio.core.common import PSM_USECOLS, PSM_MAP, PSM_SCHEMA, PEP
from quantmsio.core.mztab import MzTab
import pandas as pd


class Psm(MzTab):
    def __init__(self, mztab_path: Union[Path, str]):
        super(Psm, self).__init__(mztab_path)
        self._ms_runs = self.extract_ms_runs()
        self._protein_global_qvalue_map = self.get_protein_map()
        self._score_names = self.get_score_names()
        self._modifications = self.get_modifications()
        self._mods_map = self.get_mods_map()
        self._automaton = get_ahocorasick(self._mods_map)

    def iter_psm_table(
        self, chunksize: int = 1000000, protein_str: Optional[str] = None
    ):
        for df in self.skip_and_load_csv("PSH", chunksize=chunksize):
            if protein_str:
                df = df[df["accession"].str.contains(f"{protein_str}", na=False)]
            no_cols = set(PSM_USECOLS) - set(df.columns)
            for col in no_cols:
                df.loc[:, col] = None
            psm_map = PSM_MAP.copy()
            for key in PEP:
                if key in df.columns:
                    psm_map[key] = "posterior_error_probability"
                    break
            df.rename(columns=psm_map, inplace=True)
            df.loc[:, "additional_scores"] = df[
                list(self._score_names.values()) + ["global_qvalue"]
            ].apply(self._genarate_additional_scores, axis=1)
            df.loc[:, "cv_params"] = df[["consensus_support"]].apply(
                self._generate_cv_params, axis=1
            )
            df.loc[:, "reference_file_name"] = df["spectra_ref"].apply(
                lambda x: self._ms_runs[x[: x.index(":")]]
            )
            yield df

    def _extract_pep_columns(self) -> None:
        if os.stat(self.mztab_path).st_size == 0:
            raise ValueError("File is empty")
        f = open(self.mztab_path)
        pos = self._get_pos("PEH")
        f.seek(pos)
        line = f.readline()
        while not line.startswith("PEH"):
            line = f.readline()
        self._pep_columns = line.split("\n")[0].split("\t")

    def extract_from_pep(self, chunksize: int = 2000000) -> dict:
        self._extract_pep_columns()
        pep_usecols: list = [
            "opt_global_cv_MS:1000889_peptidoform_sequence",
            "charge",
            "best_search_engine_score[1]",
            "spectra_ref",
        ]
        live_cols: list = [col for col in pep_usecols if col in self._pep_columns]
        not_cols: list = [col for col in pep_usecols if col not in live_cols]
        if "opt_global_cv_MS:1000889_peptidoform_sequence" in not_cols:
            if "sequence" in self._pep_columns and "modifications" in self._pep_columns:
                live_cols.append("sequence")
                live_cols.append("modifications")
            else:
                raise Exception(
                    "The peptide table don't have opt_global_cv_MS:1000889_peptidoform_sequence columns"
                )
        if "charge" in not_cols or "best_search_engine_score[1]" in not_cols:
            raise Exception(
                "The peptide table don't have best_search_engine_score[1] or charge columns"
            )
        pep_map: dict = {}
        indexs: list = [self._pep_columns.index(col) for col in live_cols]
        for pep in self.skip_and_load_csv("PEH", usecols=indexs, chunksize=chunksize):
            pep.reset_index(drop=True, inplace=True)
            if "opt_global_cv_MS:1000889_peptidoform_sequence" not in pep.columns:
                pep.loc[:, "opt_global_cv_MS:1000889_peptidoform_sequence"] = pep[
                    ["modifications", "sequence"]
                ].apply(
                    lambda row: get_petidoform_msstats_notation(
                        row["sequence"], row["modifications"], self._modifications
                    ),
                    axis=1,
                )
            # check spectra_ref
            if "spectra_ref" not in pep.columns:
                pep.loc[:, "scan_number"] = None
                pep.loc[:, "spectra_ref"] = None
            else:
                pep.loc[:, "scan_number"] = pep["spectra_ref"].apply(
                    generate_scan_number
                )
                pep["spectra_ref"] = pep["spectra_ref"].apply(
                    lambda x: self._ms_runs[x.split(":")[0]]
                )
            pep_msg = pep.iloc[
                pep.groupby(
                    ["opt_global_cv_MS:1000889_peptidoform_sequence", "charge"]
                ).apply(lambda row: row["best_search_engine_score[1]"].idxmin())
            ]
            pep_msg = pep_msg.set_index(
                ["opt_global_cv_MS:1000889_peptidoform_sequence", "charge"]
            )

            pep_msg.loc[:, "pep_msg"] = pep_msg[
                ["best_search_engine_score[1]", "spectra_ref", "scan_number"]
            ].apply(
                lambda row: [
                    row["best_search_engine_score[1]"],
                    row["spectra_ref"],
                    row["scan_number"],
                ],
                axis=1,
            )
            map_dict = pep_msg.to_dict()["pep_msg"]
            for key, value in map_dict.items():
                if key not in pep_map:
                    pep_map[key] = value
                elif value[0] < pep_map[key][0]:
                    pep_map[key] = value
        return pep_map

    @staticmethod
    def slice(df: pd.DataFrame, partitions: list):
        cols = df.columns
        if not isinstance(partitions, list):
            raise Exception(f"{partitions} is not a list")
        if len(partitions) == 0:
            raise Exception(f"{partitions} is empty")
        for partion in partitions:
            if partion not in cols:
                raise Exception(f"{partion} does not exist")
        for key, df in df.groupby(partitions):
            yield key, df

    def generate_report(
        self, chunksize: int = 1000000, protein_str: Optional[str] = None
    ):
        for df in self.iter_psm_table(chunksize=chunksize, protein_str=protein_str):
            self.transform_psm(df)
            self.add_addition_msg(df)
            self.convert_to_parquet_format(df)
            df = self.transform_parquet(df)
            yield df

    @staticmethod
    def _generate_cv_params(rows: pd.Series) -> list | None:
        cv_list: list = []
        if rows["consensus_support"]:
            struct = {
                "cv_name": "consesus_support",
                "cv_value": str(rows["consensus_support"]),
            }
            cv_list.append(struct)
        if len(cv_list) > 0:
            return cv_list
        else:
            return None

    def transform_psm(self, df: pd.DataFrame) -> None:
        select_mods: list = list(self._mods_map.keys())
        df[["peptidoform", "modifications"]] = df[["peptidoform"]].apply(
            lambda row: self.generate_modifications_details(
                row["peptidoform"], self._mods_map, self._automaton, select_mods
            ),
            axis=1,
            result_type="expand",
        )
        df.loc[:, "scan"] = df["spectra_ref"].apply(generate_scan_number)
        df.drop(
            ["spectra_ref", "search_engine", "search_engine_score[1]"],
            inplace=True,
            axis=1,
        )

    @staticmethod
    def transform_parquet(df: pd.DataFrame) -> pa.Table:
        return pa.Table.from_pandas(df, schema=PSM_SCHEMA)

    def _genarate_additional_scores(self, cols: pd.Series) -> list:
        struct_list: list = []
        for software, score in self._score_names.items():
            software = re.sub(r"[^a-zA-Z0-9\s]", "", software)
            software = software.lower()
            struct = {"score_name": f"{software}_score", "score_value": cols[score]}
            struct_list.append(struct)
        if cols["global_qvalue"]:
            struct = {
                "score_name": "global_qvalue",
                "score_value": cols["global_qvalue"],
            }
            struct_list.append(struct)
        return struct_list

    @staticmethod
    def add_addition_msg(df: pd.DataFrame) -> None:
        df.loc[:, "predicted_rt"] = None
        df.loc[:, "ion_mobility"] = None
        df.loc[:, "number_peaks"] = None
        df.loc[:, "mz_array"] = None
        df.loc[:, "intensity_array"] = None
        df.loc[:, "charge_array"] = None
        df.loc[:, "ion_type_array"] = None
        df.loc[:, "ion_mobility_array"] = None

    def write_psm_to_file(
        self,
        output_path: str,
        chunksize: int = 1000000,
        protein_file: Optional[str] = None,
    ) -> None:
        protein_list: list = (
            extract_protein_list(protein_file) if protein_file else None
        )
        protein_str: Optional[str] = "|".join(protein_list) if protein_list else None
        pqwriter: pq.ParquetWriter | None = None
        for p in self.generate_report(chunksize=chunksize, protein_str=protein_str):
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, p.schema)
            pqwriter.write_table(p)
        if pqwriter:
            pqwriter.close()

    @staticmethod
    def convert_to_parquet_format(res: pd.DataFrame) -> None:
        res["mp_accessions"] = res["mp_accessions"].apply(get_protein_accession)
        res["precursor_charge"] = (
            res["precursor_charge"]
            .map(lambda x: None if pd.isna(x) else int(x))
            .astype("Int32")
        )
        res["calculated_mz"] = res["calculated_mz"].astype(float)
        res["observed_mz"] = res["observed_mz"].astype(float)
        res["posterior_error_probability"] = res["posterior_error_probability"].astype(
            float
        )
        res["is_decoy"] = (
            res["is_decoy"]
            .map(lambda x: None if pd.isna(x) else int(x))
            .astype("Int32")
        )
        res["scan"] = res["scan"].astype(str)
        if "rt" in res.columns:
            res["rt"] = res["rt"].astype(float)
        else:
            res.loc[:, "rt"] = None
