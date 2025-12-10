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
from typing import Optional, Union

import pandas as pd
from pandas import DataFrame

from qpx.core.common import SDRF_MAP, SDRF_USECOLS

logger = logging.getLogger(__name__)


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
        return re.search("NT=(.+?)(;|$)", sdrf_value).group(1)
    else:
        return sdrf_value


def get_complex_value_sdrf_column(sdrf_table: DataFrame, column: str) -> list:
    """
    Get the complex values from a SDRF column
    :param sdrf_table: pandas dataframe
    :param column: column name
    """
    values = get_unique_from_column_substr(sdrf_table, column)
    return [get_name_from_complex_sdrf_value(value) for value in values]


def get_acquisition_method(
    sdrf_table: DataFrame, acquisition_method_column: str, column_labeling: str
) -> list:
    """
    Get the acquisition method from the SDRF table.Returns the acquisition method and the labeling method.
    Three different methods are supported: label free, TMT and iTRAQ.
    For DIA methods, the acquisition method is Data-dependent
    acquisition and MUST be annotated in the SDRF.
    :param sdrf_table: Pandas dataframe
    :param acquisition_method_column: acquisition method column name
    :param column_labeling: labeling column name
    """
    acquisition_values = get_complex_value_sdrf_column(
        sdrf_table, acquisition_method_column
    )
    labeling_values = get_complex_value_sdrf_column(sdrf_table, column_labeling)
    if len(acquisition_values) == 0 and len(labeling_values) > 0:
        for labeling_value in labeling_values:
            if (
                "label free" in labeling_value.lower()
                or "label-free" in labeling_value.lower()
            ):
                acquisition_values.append("Label free")
                acquisition_values.append("Data-dependent acquisition")
            elif "tmt" in labeling_value.lower():
                acquisition_values.append("TMT")
                acquisition_values.append("Data-dependent acquisition")
            elif "itraq" in labeling_value.lower():
                acquisition_values.append("iTRAQ")
                acquisition_values.append("Data-dependent acquisition")
    return acquisition_values


class SDRFHandler:
    ORGANISM_COLUMN = "characteristics[organism]"
    ORGANISM_PART_COLUMN = "characteristics[organism part]"
    DISEASE_COLUMN = "characteristics[disease]"
    CELL_LINE_COLUMN = "characteristics[cell line]"
    INSTRUMENT_COLUMN = "comment[instrument]"
    ENZYME_COLUMN = "comment[cleavage agent details]"
    ACQUISITION_METHOD = "comment[proteomics data acquisition method]"
    DISSOCIATION_METHOD = "comment[dissociation method]"
    FRAGMENT_MASS_TOLERANCE = "comment[fragment mass tolerance]"
    PRECURSOR_MASS_TOLERANCE = "comment[precursor mass tolerance]"
    LABELING = "comment[label]"

    # The supported labeling methods
    # SUPOORTED_LABELING = [
    #     "LABEL FREE",
    #     "TMT10",
    #     "TMT11",
    #     "TMT16",
    #     "TMT6",
    #     "ITRAQ4",
    #     "ITRAQ8",
    # ]

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
            self.sdrf_table = sdrf_table
        except FileNotFoundError:
            raise FileNotFoundError(
                "The SDRF file provided not found: " + str(sdrf_file)
            )

    def get_organisms(self):
        return get_unique_from_column_substr(self.sdrf_table, self.ORGANISM_COLUMN)

    def get_organism_parts(self):
        return get_unique_from_column_substr(self.sdrf_table, self.ORGANISM_PART_COLUMN)

    def get_diseases(self):
        return get_unique_from_column_substr(self.sdrf_table, self.DISEASE_COLUMN)

    def get_cell_lines(self):
        return get_unique_from_column_substr(self.sdrf_table, self.CELL_LINE_COLUMN)

    def get_instruments(self):
        """
        Get the instruments used in the experiment
        """
        return get_complex_value_sdrf_column(self.sdrf_table, self.INSTRUMENT_COLUMN)

    def get_enzymes(self):
        return get_complex_value_sdrf_column(self.sdrf_table, self.ENZYME_COLUMN)

    def get_runs(self):
        return self.sdrf_table["comment[data file]"].str.split(".").str[0].unique()

    def get_acquisition_properties(self):
        """
        Get the acquisition properties
        """
        acquisition_values = []
        [
            acquisition_values.append(
                {"proteomics data acquisition method": acquisition_value}
            )
            for acquisition_value in get_acquisition_method(
                self.sdrf_table, self.ACQUISITION_METHOD, self.LABELING
            )
        ]
        [
            acquisition_values.append({"dissociation method": dissociation_value})
            for dissociation_value in get_complex_value_sdrf_column(
                self.sdrf_table, self.DISSOCIATION_METHOD
            )
        ]
        [
            acquisition_values.append(
                {"precursor mass tolerance": precursor_mass_tolerance_value}
            )
            for precursor_mass_tolerance_value in get_complex_value_sdrf_column(
                self.sdrf_table, self.PRECURSOR_MASS_TOLERANCE
            )
        ]
        [
            acquisition_values.append(
                {"fragment mass tolerance": fragment_mass_tolerance_value}
            )
            for fragment_mass_tolerance_value in get_complex_value_sdrf_column(
                self.sdrf_table, self.FRAGMENT_MASS_TOLERANCE
            )
        ]
        return acquisition_values

    def get_factor_value(self) -> Optional[str]:
        """
        Get the factor value
        """
        selected_columns = [
            column for column in self.sdrf_table.columns if "factor value" in column
        ]
        if len(selected_columns) != 1:
            return None
        values = re.findall(r"\[(.*?)\]", selected_columns[0])
        if len(values) != 1:
            return None
        return values[0]

    def extract_feature_properties(self) -> DataFrame:
        """
        Extract the feature properties from the SDRF file. These properties are needed by the feature file and
        FeatureHandler class. The experiment type can be: "lfq", "tmt" or "itraq".
        """

        experiment_type = self.get_experiment_type_from_sdrf()
        sdrf_pd = self.sdrf_table.copy()  # type: DataFrame

        sdrf_pd["comment[data file]"] = sdrf_pd["comment[data file]"].apply(
            lambda x: x.split(".")[0]
        )

        factor_columns = [
            column for column in sdrf_pd.columns if "factor value" in column
        ]
        if len(factor_columns) != 1:
            sdrf_pd["condition"] = (
                sdrf_pd[factor_columns].astype(str).agg("|".join, axis=1)
            )
        else:
            sdrf_pd["condition"] = sdrf_pd[factor_columns]

        sdrf_pd = sdrf_pd.rename(
            columns={
                "comment[data file]": "reference_file_name",
                "source name": "sample_accession",
                "comment[fraction identifier]": "fraction",
                "comment[label]": "channel",
            }
        )
        experiment_type = re.sub("[\\d]", "", experiment_type)
        # Add the channel column if it is not present
        if experiment_type.upper() not in ["TMT", "ITRAQ", "LFQ"]:
            raise ValueError(
                "The experiment type provided is not supported: {}, available values [lfq,tmt,itraq]".format(
                    experiment_type
                )
            )

        # extract
        if experiment_type.upper() != "LFQ":
            sdrf = sdrf_pd[
                [
                    "reference_file_name",
                    "sample_accession",
                    "condition",
                    "fraction",
                    "channel",
                ]
            ]
            return sdrf
        sdrf_pd.loc[:, "channel"] = None  # Channel will be needed in the LFQ as empty.
        sdrf = sdrf_pd[
            [
                "reference_file_name",
                "sample_accession",
                "condition",
                "fraction",
                "channel",
            ]
        ]
        return sdrf

    def get_experiment_type_from_sdrf(self):
        """
        Using the SDRF file label column, we will try to extract the experiment type of an SDRF.
        The three possible values supported in SDRF are lfq, tmt and itraq.
        """
        if self.LABELING not in self.sdrf_table.columns:
            raise ValueError(
                "The SDRF file provided does not contain the comment[label] column"
            )

        labeling_values = get_complex_value_sdrf_column(self.sdrf_table, self.LABELING)
        if len(labeling_values) == 0:
            raise ValueError(
                "The SDRF file provided does not contain any comment[label] value"
            )

        labeling_values = [i.upper() for i in labeling_values]

        if len([i for i in labeling_values if "LABEL FREE" in i]) > 0:
            return "LFQ"
        elif len([i for i in labeling_values if "TMT" in i]) > 0:
            if len(labeling_values) == 10:
                return "TMT10"
            elif len(labeling_values) == 11:
                return "TMT11"
            elif len(labeling_values) == 16:
                return "TMT16"
            elif len(labeling_values) == 6:
                return "TMT6"
            else:
                raise ValueError(
                    "The SDRF file provided does not contain a supported TMT comment[label] value"
                )
        elif len([i for i in labeling_values if "ITRAQ" in i]) > 0:
            if len(labeling_values) == 4:
                return "ITRAQ4"
            elif len(labeling_values) == 8:
                return "ITRAQ8"
            else:
                raise ValueError(
                    "The SDRF file provided does not contain a supported iTRAQ comment[label] value"
                )
        else:
            raise ValueError(
                "The SDRF file provided does not contain any supported comment[label] value"
            )

    def get_sample_labels(self):
        """
        Get the sample labels from the SDRF file. The sample labels are the values of the column "comment[label]".
        :return: Set of sample labels
        """
        labels = self.sdrf_table["comment[label]"].unique()
        return set(labels)

    def get_sample_map(self):
        """
        Get the sample accession map from the sdrf file. The key of the sample map is:
        - data file + :_: + sample label
        The value of the sample map is the sample accession.
        :return: Sample map
        """
        sample_map = {}
        sdrf_pd = self.sdrf_table.copy()  # type: DataFrame
        sdrf_pd["comment[data file]"] = sdrf_pd["comment[data file]"].apply(
            lambda x: x.split(".")[0]
        )
        for _, row in sdrf_pd.iterrows():
            channel = (
                "LABEL FREE SAMPLE"
                if "LABEL FREE" in row["comment[label]"].upper()
                else row["comment[label]"]
            )
            if row["comment[data file]"] + ":_:" + channel in sample_map:
                if (
                    sample_map[row["comment[data file]"] + ":_:" + channel]
                    != row["source name"]
                ):
                    raise ValueError("The sample map is not unique")
                else:
                    logger.info(
                        "channel {} for sample {} already in the sample map".format(
                            channel, row["source name"]
                        )
                    )
            sample_map[row["comment[data file]"] + ":_:" + channel] = row["source name"]
        return sample_map

    def get_mods(self):
        sdrf = self.sdrf_table
        mod_cols = [
            col
            for col in sdrf.columns
            if (
                col.startswith("comment[modification parameter]")
                | col.startswith("comment[modification parameters]")
            )
        ]
        fix_m = []
        variable_m = []
        for col in mod_cols:
            mod_msg = sdrf[col].values[0].split(";")
            mod_dict = {k.split("=")[0]: k.split("=")[1] for k in mod_msg}
            mod = (
                f"{mod_dict['NT']} ({mod_dict['TA']})"
                if "TA" in mod_dict
                else f"{mod_dict['NT']} ({mod_dict['PP']})"
            )
            if mod_dict["MT"] == "Variable" or mod_dict["MT"] == "variable":
                variable_m.append(mod)
            else:
                fix_m.append(mod)
        fix_s = ",".join(fix_m) if len(fix_m) > 0 else "null"
        variable_s = ",".join(variable_m) if len(variable_m) > 0 else "null"

        return fix_s, variable_s

    def get_mods_dict(self):
        sdrf = self.sdrf_table
        mod_cols = [
            col
            for col in sdrf.columns
            if (
                col.startswith("comment[modification parameter]")
                | col.startswith("comment[modification parameters]")
            )
        ]
        mods = {}
        for col in mod_cols:
            mod_msg = sdrf[col].values[0].split(";")
            mod_dict = {k.split("=")[0]: k.split("=")[1] for k in mod_msg}
            accession = mod_dict["AC"]
            name = mod_dict["NT"]
            site = mod_dict["TA"] if "TA" in mod_dict else mod_dict["PP"]
            mods[accession] = [name, site]
            mods[name] = [accession, site]
        return mods

    def get_sample_map_run(self):
        sdrf = self.sdrf_table[
            ["source name", "comment[data file]", "comment[label]"]
        ].copy()
        sdrf["comment[data file]"] = sdrf["comment[data file]"].str.split(".").str[0]
        if self.get_experiment_type_from_sdrf() != "LFQ":
            sdrf.loc[:, "map_sample"] = (
                sdrf["comment[data file]"] + "-" + sdrf["comment[label]"]
            )
        else:
            sdrf.loc[:, "map_sample"] = sdrf["comment[data file]"]
        sdrf.set_index("map_sample", inplace=True)
        sample_map = sdrf.to_dict()["source name"]
        return sample_map

    def transform_sdrf(self):
        factor = list(filter(lambda x: x.startswith("factor"), self.sdrf_table.columns))
        usecols = list(SDRF_USECOLS) + factor
        sdrf = self.sdrf_table[usecols].copy()
        sdrf["comment[data file]"] = sdrf["comment[data file]"].apply(
            lambda x: x.split(".")[0]
        )
        samples = sdrf["source name"].unique()
        mixed_map = dict(zip(samples, range(1, len(samples) + 1)))

        # Keep condition order consistent with factor numeric order
        def get_condition_index(col_name):
            match = re.search(r"\.(\d+)$", col_name)
            return (int(match.group(1)) if match else 0, col_name)

        sorted_factor = sorted(factor, key=get_condition_index)
        sdrf.loc[:, "condition"] = sdrf[sorted_factor].astype(str).agg("|".join, axis=1)

        sdrf.loc[:, "run"] = sdrf[
            [
                "source name",
                "comment[technical replicate]",
                "comment[fraction identifier]",
            ]
        ].apply(
            lambda row: str(mixed_map[row["source name"]])
            + "_"
            + str(row["comment[technical replicate]"])
            + "_"
            + str(row["comment[fraction identifier]"]),
            axis=1,
        )
        sdrf.drop(
            ["comment[technical replicate]"] + factor,
            axis=1,
            inplace=True,
        )
        sdrf.rename(columns=SDRF_MAP, inplace=True)
        # return sdrf[["reference_file_name", "channel", "condition"]].copy()
        return sdrf
