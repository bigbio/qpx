"""
Common constants and schemas for qpx.
This module provides mapping dictionaries, column lists, and schemas used across the library.
"""

from typing import Dict, List, Set

import pyarrow as pa

from qpx import __version__
from qpx.core.format import FEATURE_FIELDS, IBAQ_FIELDS, PG_FIELDS, PSM_FIELDS

OPENMS_PEPTIDOFORM_COLUMN = "opt_global_cv_MS:1000889_peptidoform_sequence"
OPENMS_POSTERIOR_ERRORPROBABILITY = "opt_global_Posterior_Error_Probability_score"
OPENMS_IS_DECOY = "opt_global_cv_MS:1002217_decoy_peptide"

OPENMS_NAMES_MAP = {
    "OpenMS:Target-decoy PSM q-value": "qvalue",
    "OpenMS:Best PSM Score": "best_psm_score",
    "OpenMS:Target-decoy protein q-value": "pg_qvalue",
}

# PSM mapping and columns
PSM_MAP: Dict[str, str] = {
    "sequence": "sequence",
    "modifications": "modifications",
    OPENMS_PEPTIDOFORM_COLUMN: "peptidoform",
    "opt_global_q-value": "global_qvalue",
    "opt_global_consensus_support": "consensus_support",
    OPENMS_IS_DECOY: "is_decoy",
    "calc_mass_to_charge": "calculated_mz",
    "accession": "mp_accessions",
    "charge": "precursor_charge",
    "exp_mass_to_charge": "observed_mz",
    "retention_time": "rt",
    OPENMS_POSTERIOR_ERRORPROBABILITY: "posterior_error_probability",
    "spectra_ref": "spectra_ref",
}

# Regular expression that make sure the protein search engine score contains the following words: openms target-decoy protein and q-value
OPENMS_PROTEIN_QVALUE_WORDS = ["openms", "target-decoy", "protein", "q-value"]

PEP: List[str] = [
    "opt_global_Posterior_Error_Probability_score",
    "opt_global_Posterior_Error_Probability",
    "opt_global_MS:1001493_score",
]

# MsStats mapping and columns
MSSTATS_MAP: Dict[str, str] = {
    "ProteinName": "pg_accessions",
    "Reference": "reference_file_name",
    "Intensity": "intensity",
    "Channel": "channel",
    "RetentionTime": "rt",
    "PeptideSequence": "peptidoform",
}

# Pre-compute sets for faster membership testing
MSSTATS_USECOLS: Set[str] = set(MSSTATS_MAP.keys())

# SDRF mapping and columns
SDRF_MAP: Dict[str, str] = {
    "comment[data file]": "reference_file_name",
    "comment[label]": "channel",
    "source name": "sample_accession",
    "comment[fraction identifier]": "fraction",
    "characteristics[biological replicate]": "biological_replicate",
}

# Pre-compute sets for faster membership testing
SDRF_USECOLS: Set[str] = set(list(SDRF_MAP.keys()) + ["comment[technical replicate]"])

# "File.Name" was removed from the main report (report.parquet) starting from DIA-NN version 2.0.
DIANN_MAP = {
    "Precursor.Quantity": "intensity",
    "RT.Start": "rt_start",
    "RT.Stop": "rt_stop",
    "RT": "rt",
    "Predicted.RT": "predicted_rt",
    "Protein.Group": "pg_accessions",
    "Protein.Ids": "mp_accessions",
    "PEP": "posterior_error_probability",
    "Global.Q.Value": "global_qvalue",
    "Global.PG.Q.Value": "pg_global_qvalue",
    "Q.Value": "qvalue",
    "PG.Q.Value": "pg_qvalue",
    "Precursor.Normalised": "normalize_intensity",
    "PG.MaxLFQ": "lfq",
    "Quantity.Quality": "precursor_quantification_score",
    "Precursor.Charge": "precursor_charge",
    "Stripped.Sequence": "sequence",
    "Modified.Sequence": "peptidoform",
    "Genes": "gg_names",
    "Run": "reference_file_name",
    # Ion mobility columns (for TIMS-TOF data)
    "IM": "ion_mobility",
    # Decoy indicator (added in DIA-NN 2.0)
    "Decoy": "is_decoy",
    # Precursor identifier
    "Precursor.Id": "precursor_id",
}

# DIA-NN columns that map to additional_scores (quality metrics)
# These are read separately and structured into the additional_scores field
DIANN_SCORE_COLS: List[str] = [
    "CScore",
    "Evidence",
    "Spectrum.Similarity",
    "Averagine",
    "Mass.Evidence",
    "Decoy.Evidence",
    "Decoy.Score",
    "Lib.Q.Value",
    "Lib.PG.Q.Value",
    "Lib.PEP",
]

# DIA-NN columns that map to additional_intensities (MS1 quantification)
DIANN_INTENSITY_COLS: List[str] = [
    "Ms1.Area",
]

# DIA-NN ion mobility columns for cv_params
DIANN_IM_COLS: List[str] = [
    "iIM",  # Reference ion mobility from spectral library
    "Predicted.IM",  # Predicted ion mobility based on iIM
    "Predicted.iIM",  # iIM predicted based on measured IM
]

# "PG.Quantity" and "PG.Normalised" were removed from the main report
#   (report.parquet) starting from DIA-NN version 2.0.
# "Decoy" was added into main report starting from DIA-NN version 2.0.
# "N.Sequences" and "N.Proteotypic.Sequences" were added into
#   report.pg_matrix.tsv starting from DIA-NN version 2.0.
DIANN_PG_MAP = {
    "Protein.Group": "pg_accessions",
    "Protein.Names": "pg_names",
    "Genes": "gg_accessions",
    "Run": "reference_file_name",
    "Global.PG.Q.Value": "global_qvalue",
    "PG.MaxLFQ": "lfq",
    "PG.Q.Value": "qvalue",
    "Proteotypic": "proteotypic",
    "Stripped.Sequence": "stripped_sequence",
    "Precursor.Id": "precursor_id",
}

DIANN_PG_MATRIX_MAP = {
    "Protein.Group": "pg_accessions",
    "Protein.Names": "pg_names",
    "Genes": "gg_accessions",
}

# Core columns that are always read and mapped
DIANN_USECOLS = list(DIANN_MAP.keys())
DIANN_PG_USECOLS = list(DIANN_PG_MAP.keys())

# Extended columns for optional advanced metrics (may not exist in all DIA-NN versions)
DIANN_EXTENDED_COLS = DIANN_SCORE_COLS + DIANN_INTENSITY_COLS + DIANN_IM_COLS

MAXQUANT_PSM_MAP = {
    "Sequence": "sequence",
    "Proteins": "protein_accessions",
    "PEP": "posterior_error_probability",
    "Modified sequence": "peptidoform",
    "Reverse": "is_decoy",
    "m/z": "observed_mz",
    "Scan number": "scan",
    "Retention time": "rt",
    "Charge": "precursor_charge",
    "Raw file": "reference_file_name",
    "Score": "andromeda_score",
    "Delta score": "andromeda_delta_score",
    "PIF": "parent_ion_fraction",
    "Masses": "mz_array",
    "Intensities": "intensity_array",
    "Number of matches": "number_peaks",
    "1/K0": "ion_mobility",
}

MAXQUANT_FEATURE_MAP = {
    "Sequence": "sequence",
    "Proteins": "mp_accessions",
    "Leading proteins": "pg_accessions",
    "Leading razor protein": "anchor_protein",
    "Gene names": "gg_names",
    "PEP": "posterior_error_probability",
    "Modified sequence": "peptidoform",
    "Charge": "precursor_charge",
    "Raw file": "reference_file_name",
    "Score": "andromeda_score",
    "Delta score": "andromeda_delta_score",
    "PIF": "parent_ion_fraction",
    "Reverse": "is_decoy",
    "m/z": "observed_mz",
    "MS/MS scan number": "scan",
    "Calibrated retention time": "rt",
    "Calibrated retention time start": "rt_start",
    "Calibrated retention time finish": "rt_stop",
    "Intensity": "intensity",
    "1/K0": "ion_mobility",
}

IBAQ_USECOLS = [
    "pg_accessions",
    "peptidoform",
    "sequence",
    "precursor_charge",
    "intensities",
    "reference_file_name",
    "unique",
]

MAXQUANT_PSM_USECOLS = list(MAXQUANT_PSM_MAP.keys())

MAXQUANT_PSM_COMPUTED = {
    "modifications": "parsed from peptidoform using PyOpenMS AASequence",
    "calculated_mz": "calculated from peptidoform monoisotopic mass and precursor_charge",
    "additional_scores": "structured from andromeda_score, andromeda_delta_score, parent_ion_fraction fields",
    "cv_params": "extracted from Fragmentation, Mass analyzer, Type columns as CV parameters",
    "charge_array": "parsed from Matches column fragment ion charges",
    "ion_type_array": "parsed from Matches column fragment ion types",
}

MAXQUANT_FEATURE_USECOLS = list(MAXQUANT_FEATURE_MAP.keys())

MAXQUANT_FEATURE_COMPUTED = {
    "modifications": "parsed from peptidoform using PyOpenMS AASequence",
    "calculated_mz": "calculated from peptidoform monoisotopic mass and precursor_charge",
    "additional_scores": "structured from andromeda_score, andromeda_delta_score, parent_ion_fraction fields",
    "cv_params": "extracted from Type column as CV parameters",
    "intensities": "structured from Intensity column or TMT reporter intensity columns with sample and channel information",
    "additional_intensities": "structured from corrected reporter intensity columns for TMT experiments",
    "unique": "calculated as 1 if peptide maps to only one protein group (len(pg_accessions)==1), otherwise 0",
}

MAXQUANT_PG_MAP = {
    "Protein IDs": "pg_accessions",
    "Protein names": "pg_names",
    "Gene names": "gg_accessions",
    "Q-value": "global_qvalue",
    "Intensity": "intensity",
    "LFQ intensity": "lfq_intensity",
    "iBAQ": "ibaq_intensity",
    "Number of proteins": "number_of_proteins",
    "Peptides": "peptide_count_total",
    "Razor + unique peptides": "peptide_count_razor_unique",
    "Unique peptides": "peptide_count_unique",
    "Sequence coverage [%]": "sequence_coverage",
    "Mol. weight [kDa]": "molecular_weight",
    "Score": "andromeda_score",
    "Reverse": "is_decoy",
    "Potential contaminant": "contaminant",
    "MS/MS count": "msms_count",
}

MAXQUANT_PG_USECOLS = list(MAXQUANT_PG_MAP.keys())

MAXQUANT_PG_COMPUTED = {
    "additional_scores": "structured from andromeda_score field",
    "intensities": "structured from sample-specific Intensity columns with sample and channel information from SDRF",
    "additional_intensities": "structured from LFQ intensity, iBAQ, and other intensity columns with name-value pairs",
    "reference_file_name": "derived from the sample with maximum intensity value",
    "peptides": "distributed peptide counts among protein group members from peptide_count_total",
    "anchor_protein": "first protein from Majority protein IDs column or first element of pg_accessions",
    "peptide_counts": "structured from peptide_count_unique and peptide_count_total fields",
    "feature_counts": "structured using peptide_count_unique as unique_features and peptide_count_total as total_features",
}

# mzTab protein group mapping
MZTAB_PG_MAP = {
    "accession": "anchor_protein",
    "best_search_engine_score[1]": "global_qvalue",
    "ambiguity_members": "pg_accessions",
    "protein_coverage": "sequence_coverage",
    "opt_global_Posterior_Probability_score": "posterior_probability_score",
    "opt_global_nr_found_peptides": "peptide_count",
    "opt_global_cv_PRIDE:0000303_decoy_hit": "is_decoy",
    "opt_global_result_type": "result_type",
}

MZTAB_PG_USECOLS = list(MZTAB_PG_MAP.keys()) + ["description"]

QPX_VERSION = __version__

PSM_SCHEMA = pa.schema(
    PSM_FIELDS,
    metadata={"description": "psm file in QPX format"},
)
FEATURE_SCHEMA = pa.schema(
    FEATURE_FIELDS,
    metadata={"description": "feature file in QPX format"},
)

IBAQ_SCHEMA = pa.schema(
    IBAQ_FIELDS,
    metadata={"description": "ibaq file in QPX format"},
)
PG_SCHEMA = pa.schema(
    PG_FIELDS,
    metadata={"description": "PG file in QPX format"},
)
