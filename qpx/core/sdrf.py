"""
The SDRF classes are used to parse and handle SDRF file format within the QPX package. The SDRF file format is
used to describe the experimental design of a proteomics experiment. The SDRF file format is described in the docs
folder of this repository.
This module contains the following classes:
    * SDRFHandler - class to handle SDRF files
"""

import logging
import re
from pathlib import Path
from typing import Union

import pandas as pd
from pandas import DataFrame

from qpx.core.files import run_file_stem

logger = logging.getLogger(__name__)


def validate_sdrf_data_files(sdrf_table: DataFrame) -> None:
    """Require a non-empty ``comment[data file]`` value in every SDRF row."""
    column = "comment[data file]"
    if column not in sdrf_table.columns:
        raise ValueError(f"SDRF is missing required column '{column}'")

    values = sdrf_table[column]
    invalid = values.isna() | values.astype("string").str.strip().eq("").fillna(True)
    if invalid.any():
        rows = [str(position + 2) for position, is_invalid in enumerate(invalid) if is_invalid]
        raise ValueError(f"SDRF column '{column}' contains missing or blank value(s) at input row(s): {', '.join(rows)}")


def get_unique_from_column_substr(sdrf_table: DataFrame, substr: str) -> list:
    """
    Get in a pandas dataframe the columns that contain a given substring
    :param sdrf_table: pandas dataframe
    :param substr: substring
    """
    selected_columns = [column for column in sdrf_table.columns if substr in column]

    if len(selected_columns) == 0:
        return []

    # Extract unique values from selected columns
    unique_values = sdrf_table[selected_columns].values.flatten()
    return pd.unique(unique_values).tolist()


def get_name_from_complex_sdrf_value(sdrf_value: str) -> str:
    """
    Get the name from a complex SDRF value
    :param sdrf_value: SDRF value
    :return: name
    """
    if "NT=" in sdrf_value:
        match = re.search(r"NT=(.+?)(;|$)", sdrf_value)
        return match.group(1) if match else sdrf_value
    else:
        return sdrf_value


def parse_sdrf_kv(sdrf_value: str) -> dict[str, str]:
    """Split a complex SDRF value such as ``NT=Oxidation;AC=UNIMOD:35`` into stripped key/value pairs."""
    pairs = {}
    for part in str(sdrf_value).split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            pairs[key.strip()] = value.strip()
    return pairs


def parse_modification_parameters(sdrf_table: DataFrame) -> list[dict] | None:
    """Parse the ``comment[modification parameters]`` columns of a lower-cased SDRF table.

    Each column holds one modification as ``NT=...;AC=...;MT=...;TA=...;PP=...``,
    read from its first non-blank value (the search settings are shared by
    all rows). Returns one ``{accession, name, fixed, position,
    target_amino_acid}`` dict per modification, or ``None`` when there are none.
    """
    results: list[dict] = []
    for col in sdrf_table.columns:
        if not col.startswith("comment[modification parameter"):
            continue
        values = sdrf_table[col].dropna()
        mod = _modification_parameter(values.iloc[0]) if not values.empty else None
        if mod:
            results.append(mod)
    return results or None


def _modification_parameter(raw: str) -> dict | None:
    """One ``comment[modification parameters]`` value as a dict, or ``None`` without ``AC``/``NT``."""
    kv = parse_sdrf_kv(raw)
    accession = kv.get("AC", "")
    name = kv.get("NT", "")
    if not accession and not name:
        return None
    return {
        "accession": accession,
        "name": name,
        "fixed": kv.get("MT", "").lower() in ("fixed", "fix"),
        "position": kv.get("PP", None),
        "target_amino_acid": kv.get("TA", None),
    }


def get_complex_value_sdrf_column(sdrf_table: DataFrame, column: str) -> list:
    """
    Get the complex values from a SDRF column
    :param sdrf_table: pandas dataframe
    :param column: column name
    """
    values = get_unique_from_column_substr(sdrf_table, column)
    # Blank cells are read as NaN and carry no value.
    return [get_name_from_complex_sdrf_value(value) for value in values if pd.notna(value)]


class SDRFHandler:
    ENZYME_COLUMN = "comment[cleavage agent details]"
    LABELING = "comment[label]"

    def __init__(self, sdrf_file: Union[Path, str]):
        self.sdrf_file = sdrf_file
        self.sdrf_table = None
        self._load_sdrf_info(sdrf_file)

    def _load_sdrf_info(self, sdrf_file: Union[Path, str]):
        """
        Load the SDRF information from a file
        :param sdrf_file: SDRF file
        """
        try:
            sdrf_table = pd.read_csv(sdrf_file, sep="\t", header=0)
            sdrf_table.columns = sdrf_table.columns.str.lower()
            # Normalize common SDRF column name variants
            col_renames = {}
            for col in sdrf_table.columns:
                if col == "comment[datafile]":
                    col_renames[col] = "comment[data file]"
            if col_renames:
                sdrf_table = sdrf_table.rename(columns=col_renames)
            self.sdrf_table = sdrf_table
        except FileNotFoundError:
            raise FileNotFoundError("The SDRF file provided not found: " + str(sdrf_file))

    def get_enzymes(self):
        return get_complex_value_sdrf_column(self.sdrf_table, self.ENZYME_COLUMN)

    def get_modification_parameters(self) -> list[dict] | None:
        """The SDRF's search modifications; see :func:`parse_modification_parameters`."""
        return parse_modification_parameters(self.sdrf_table)

    def get_experiment_type_from_sdrf(self):
        """
        Using the SDRF file label column, we will try to extract the experiment type of an SDRF.
        The three possible values supported in SDRF are lfq, tmt and itraq.
        """
        if self.LABELING not in self.sdrf_table.columns:
            raise ValueError("The SDRF file provided does not contain the comment[label] column")

        labeling_values = get_complex_value_sdrf_column(self.sdrf_table, self.LABELING)
        if len(labeling_values) == 0:
            raise ValueError("The SDRF file provided does not contain any comment[label] value")

        labeling_values = [i.upper() for i in labeling_values]

        if any("LABEL FREE" in i for i in labeling_values):
            return "LFQ"
        if any("TMT" in i for i in labeling_values):
            if len(labeling_values) in (6, 10, 11, 16):
                return f"TMT{len(labeling_values)}"
            raise ValueError("The SDRF file provided does not contain a supported TMT comment[label] value")
        if any("ITRAQ" in i for i in labeling_values):
            if len(labeling_values) in (4, 8):
                return f"ITRAQ{len(labeling_values)}"
            raise ValueError("The SDRF file provided does not contain a supported iTRAQ comment[label] value")
        raise ValueError("The SDRF file provided does not contain any supported comment[label] value")

    def get_sample_map_run(self):
        validate_sdrf_data_files(self.sdrf_table)
        sdrf = self.sdrf_table[["source name", "comment[data file]", "comment[label]"]].copy()
        sdrf["comment[data file]"] = sdrf["comment[data file]"].map(run_file_stem)
        if self.get_experiment_type_from_sdrf() != "LFQ":
            sdrf.loc[:, "map_sample"] = sdrf["comment[data file]"] + "-" + sdrf["comment[label]"]
        else:
            sdrf.loc[:, "map_sample"] = sdrf["comment[data file]"]
        sdrf.set_index("map_sample", inplace=True)
        sample_map = sdrf.to_dict()["source name"]
        return sample_map
