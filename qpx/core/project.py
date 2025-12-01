import json
import logging
import os
import shutil
import uuid
from pathlib import Path
from typing import Optional, Union

import requests

from qpx.core.sdrf import SDRFHandler
from qpx.utils.file_utils import delete_files_extension
from qpx.utils.pride_utils import (
    get_pubmed_id_pride_json,
    get_set_of_experiment_keywords,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def check_directory(
    output_folder: str, project_accession: Optional[str] = None
) -> "ProjectHandler":
    """
    Check if the project file is present in the output folder, if not create it.
    :param output_folder: Folder where the Json file will be generated
    :param project_accession: Prefix of the Json file needed to generate the file name
    """
    if project_accession is None:
        project_json: list = [
            f for f in os.listdir(output_folder) if f.endswith("project.json")
        ]
        if len(project_json) == 1:
            json_path = f"{output_folder}/{project_json[0]}"
            project = ProjectHandler(project_json_file=json_path)
            return project
        else:
            raise Exception(f"More than one project json file found in {output_folder}")
    else:
        if os.path.exists(output_folder):
            project_json: list = [
                f for f in os.listdir(output_folder) if f.endswith("project.json")
            ]
            for json_file in project_json:
                json_path = f"{output_folder}/{json_file}"
                project = ProjectHandler(project_json_file=json_path)
                if project.project_accession == project_accession:
                    return project
            project = ProjectHandler(project_accession=project_accession)
            return project
        # If the project file not present but accession available
        else:
            os.makedirs(output_folder)
            project = ProjectHandler(project_accession=project_accession)
            return project


def create_uuid_filename(project_accession: str, extension: str) -> str:

    output_filename_path = f"{project_accession}-{str(uuid.uuid4())}{extension}"
    return output_filename_path


def get_project_accession(sdrf_path: str) -> str:
    f = open(sdrf_path)
    keys = f.readline().split("\t")
    values = f.readline().split("\t")
    sdrf_map = dict(zip(keys, values))
    project_accession = sdrf_map["source name"].split("-")[0]
    f.close()
    return project_accession


def cut_path(path: str, output_folder: str) -> str:
    path = path.replace(output_folder + "/", "")
    return path


class ProjectHandler:
    PROJECT_EXTENSION = ".project.json"

    def __init__(
        self,
        project_accession: Optional[str] = None,
        project_json_file: Optional[str] = None,
    ):
        """
        ProjectHandler class can be created using the PRIDE accession or a JSON file with the project information
        :param project_accession: PRIDE accession of the project
        :param project_json_file: JSON file with the project information
        """
        if project_accession is None and project_json_file is None:
            self.project_accession = None
            self.project: ProjectDefinition | None = None
        elif project_accession is not None and project_json_file is None:
            self.project_accession = project_accession
            self.project = ProjectDefinition(project_accession)
        else:
            self.load_project_info(project_json_file)
            self.project_accession = self.project.project_info["project_accession"]
        self.extensions: list = [
            "sdrf_file",
            "psm_file",
            "feature_file",
            "differential_file",
            "absolute_file",
        ]

    def load_project_info(self, project_json_file: Optional[str] = None) -> None:
        """
        Load the project information from a JSON file
        :param project_json_file: JSON file with the project information
        """
        try:
            with open(f"{project_json_file}", "r") as json_file:
                json_dict = json.load(json_file)
                self.project = ProjectDefinition()
                self.project.set_project_info(json_dict)
        except FileNotFoundError:
            raise FileNotFoundError(f"File {project_json_file} not found")

    def populate_from_pride_archive(self) -> None:
        # Simulate API request to PRIDE Archive
        api_url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{self.project_accession}"
        response = requests.get(api_url, timeout=1000)

        if response.status_code == 200:
            pride_data = response.json()
            self.project.project_info["project_accession"] = self.project_accession
            self.project.project_info["project_title"] = pride_data["title"]
            self.project.project_info["project_description"] = pride_data[
                "projectDescription"
            ]
            self.project.project_info["project_sample_description"] = pride_data[
                "sampleProcessingProtocol"
            ]
            self.project.project_info["project_data_description"] = pride_data[
                "dataProcessingProtocol"
            ]
            self.project.project_info["project_pubmed_id"] = get_pubmed_id_pride_json(
                pride_data
            )
            self.project.project_info["experiment_type"] = (
                get_set_of_experiment_keywords(pride_data)
            )
        else:
            logger.info("A non-pride project is being created.")

    def add_qpx_version(self, qpx_version: str) -> None:
        self.project.project_info["qpx_version"] = qpx_version

    def add_sdrf_project_properties(self, sdrf: SDRFHandler) -> None:
        """
        Add the project properties from the SDRF file to the project information.
        :param sdrf: SDRFHandler object
        """
        self.project.project_info["organisms"] = sdrf.get_organisms()
        self.project.project_info["organism_parts"] = sdrf.get_organism_parts()
        self.project.project_info["diseases"] = sdrf.get_diseases()
        self.project.project_info["cell_lines"] = sdrf.get_cell_lines()
        self.project.project_info["instruments"] = sdrf.get_instruments()
        self.project.project_info["enzymes"] = sdrf.get_enzymes()
        self.project.project_info["acquisition_properties"] = (
            sdrf.get_acquisition_properties()
        )

    def add_software_provider(
        self, sortware_name: str = "", sortware_version: str = ""
    ) -> None:
        """
        Add the sortware info to the project information
        """
        self.project.project_info["software_provider"]["name"] = sortware_name
        self.project.project_info["software_provider"]["version"] = sortware_version

    def add_quantms_path(
        self,
        path_name: str,
        file_category: str,
        is_folder: bool = False,
        partition_fields: list | None = None,
        replace_existing: bool | None = None,
    ) -> None:
        if "quantms_files" not in self.project.project_info:
            self.project.project_info["quantms_files"] = []

        obj_index: int | None = None
        for index, obj in enumerate(self.project.project_info["quantms_files"]):
            if file_category in obj:
                obj_index = index
        if obj_index is not None and replace_existing:
            self.project.project_info["quantms_files"][obj_index][file_category] = []
        path_name = path_name.replace("\\", "/")
        record: dict = {
            "path_name": path_name,
            "is_folder": is_folder,
        }
        if partition_fields is not None and isinstance(partition_fields, list):
            record["partition_fields"] = partition_fields

        if obj_index is not None:
            self.project.project_info["quantms_files"][obj_index][file_category].append(
                record
            )
        else:
            self.project.project_info["quantms_files"].append({file_category: [record]})

    def register_file(
        self,
        path_name: str,
        file_category: str,
        is_folder: bool = False,
        partition_fields: list | None = None,
        replace_existing: bool | None = None,
    ) -> None:
        if file_category not in self.extensions:
            raise Exception(f"The {file_category} not in {self.extensions}")
        self.add_quantms_path(
            path_name,
            file_category,
            is_folder=is_folder,
            partition_fields=partition_fields,
            replace_existing=replace_existing,
        )

    def save_project_info(
        self,
        output_prefix_file: Optional[str] = None,
        output_folder: Optional[str] = None,
        delete_existing: bool = False,
    ) -> None:
        """
        Save the updated project information to a JSON file. The file name will be generated automatically.
        """
        # Save the updated project info to a JSON file
        if output_prefix_file is None:
            output_prefix_file = self.project_accession

        # Create the output folder if it does not exist.
        if output_folder is not None and not os.path.exists(output_folder):
            Path(output_folder).mkdir(parents=True, exist_ok=True)

        # Delete existing SDRF file
        if delete_existing:
            delete_files_extension(output_folder, ProjectHandler.PROJECT_EXTENSION)

        if output_folder is None:
            output_filename = f"{output_prefix_file}-{str(uuid.uuid4())}{ProjectHandler.PROJECT_EXTENSION}"
        else:
            output_filename = f"{output_folder}/{output_prefix_file}-{str(uuid.uuid4())}{ProjectHandler.PROJECT_EXTENSION}"

        with open(output_filename, "w") as json_file:
            json.dump(self.project.project_info, json_file, indent=4)
        logger.info(f"Updated project information saved to {output_filename}")

    def save_updated_project_info(
        self,
        output_folder: Optional[str] = None,
        output_prefix: Optional[str] = None,
    ) -> None:
        """
        Save the updated project information to a new file.

        Args:
            output_folder: Optional folder path to save the file
            output_prefix: Optional prefix for the output file name
        """
        if output_prefix is None:
            output_prefix = self.project_accession

        if output_folder is None:
            output_filename = (
                f"{output_prefix}-{str(uuid.uuid4())}{ProjectHandler.PROJECT_EXTENSION}"
            )
        else:
            output_filename = f"{output_folder}/{output_prefix}-{str(uuid.uuid4())}{ProjectHandler.PROJECT_EXTENSION}"

        with open(output_filename, "w") as f:
            json.dump(self.project.project_info, f, indent=4)

    def populate_from_sdrf(self, sdrf_file: Union[Path, str]) -> None:
        """
        Populate the project information from an SDRF file using the SDRFHandler class.
        :param sdrf_file: SDRF file
        """
        sdrf = SDRFHandler(sdrf_file)
        self.add_sdrf_project_properties(sdrf)

    def add_sdrf_file(
        self,
        sdrf_file_path: str,
        output_folder: Optional[str],
        delete_existing: bool = True,
    ) -> None:
        """
        Copy the given file to the project folder and add the file name to the project information.
        :param sdrf_file_path: SDRF file path
        :param output_folder: Output folder
        :param delete_existing: Delete existing SDRF files
        """
        base_name = os.path.basename(sdrf_file_path).replace(".sdrf.tsv", "")
        extension = ".sdrf.tsv"

        # Create the output folder if it does not exist.
        if output_folder is not None and not os.path.exists(output_folder):
            Path(output_folder).mkdir(parents=True, exist_ok=True)

        # Delete existing SDRF file
        if delete_existing:
            delete_files_extension(output_folder, extension)

        output_filename = f"{base_name}-{str(uuid.uuid4())}{extension}"
        if output_folder is None:
            output_filename_path = output_filename
        else:
            output_filename_path = (
                f"{output_folder}/{base_name}-{str(uuid.uuid4())}{extension}"
            )

        shutil.copyfile(sdrf_file_path, output_filename_path)
        # self.project.project_info["sdrf_file"] = output_filename
        self.register_file(output_filename, "sdrf_file")
        logger.info(
            f"SDRF file copied to {output_filename} and added to the project information"
        )


class ProjectDefinition:
    """
    Class to define the project information. This class will be used to generate the project JSON file. Read more about
    the project JSON file in the docs folder of this repository
    (https://github.com/bigbio/QPX/blob/main/docs/PROJECT.md).
    """

    def __init__(self, project_accession: Optional[str] = None):
        self.project_info: dict = {
            "project_accession": "",
            "project_title": "",
            "project_description": "",
            "project_sample_description": "",
            "project_data_description": "",
            "project_pubmed_id": "",
            "organisms": [],
            "organism_parts": [],
            "diseases": [],
            "cell_lines": [],
            "instruments": [],
            "enzymes": [],
            "experiment_type": [],
            "acquisition_properties": [],
            "quantms_files": [],
            "qpx_version": "",
            "software_provider": {},
            "comments": [],
        }
        if project_accession:
            self.set_basic_info(project_accession)

    def set_basic_info(self, project_accession: str) -> None:
        self.project_info["project_accession"] = project_accession

    # Define additional methods to set specific fields as needed

    def get_project_info(self) -> dict:
        return self.project_info

    def set_project_info(self, project_info: dict) -> None:
        self.project_info = project_info
