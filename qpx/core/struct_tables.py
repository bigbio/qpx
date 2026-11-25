"""Fixed layouts for the name-value lists written by each converter."""

FIXED_STRUCT_TABLES = {
    "diann": {
        "feature": {
            "additional_scores": [
                {"name": "pg_qvalue", "higher_better": False},
                {"name": "diann_cscore", "higher_better": True},
                {"name": "diann_evidence", "higher_better": True},
                {"name": "diann_spectrum_similarity", "higher_better": True},
                {"name": "diann_ms1_profile_corr", "higher_better": True},
                {"name": "diann_mass_evidence", "higher_better": True},
                {"name": "diann_averagine", "higher_better": True},
                {"name": "diann_quantity_quality", "higher_better": True},
                {"name": "diann_lib_qvalue", "higher_better": False},
                {"name": "diann_lib_pg_qvalue", "higher_better": False},
                {"name": "diann_protein_qvalue", "higher_better": False},
                {"name": "diann_decoy_cscore", "higher_better": True},
                {"name": "diann_decoy_evidence", "higher_better": True},
                {"name": "diann_translated_quality", "higher_better": True},
                {"name": "diann_translated_qvalue", "higher_better": False},
                {"name": "diann_ptm_site_confidence", "higher_better": True},
                {"name": "diann_lib_ptm_site_confidence", "higher_better": True},
                {"name": "diann_peptidoform_qvalue", "higher_better": False},
                {"name": "diann_global_peptidoform_qvalue", "higher_better": False},
                {"name": "diann_lib_peptidoform_qvalue", "higher_better": False},
            ],
            "cv_params": [
                {"name": "proteotypic"},
                {"name": "predicted_irt"},
                {"name": "predicted_im"},
                {"name": "predicted_iim"},
                {"name": "iim"},
                {"name": "diann_precursor_lib_index"},
            ],
        },
        "pg": {
            "additional_scores": [
                {"name": "pg_pep", "higher_better": False},
                {"name": "lib_pg_qvalue", "higher_better": False},
                {"name": "protein_qvalue", "higher_better": False},
                {"name": "pg_maxlfq_quality", "higher_better": True},
                {"name": "genes_maxlfq_quality", "higher_better": True},
            ],
            "cv_params": [{"name": "quantification_method"}],
        },
        "fragment": {
            "additional_scores": [{"name": "diann_fragment_score", "higher_better": True}],
        },
    },
    "maxquant": {
        "psm": {
            "additional_scores": [
                {"name": "andromeda_score", "higher_better": True},
                {"name": "andromeda_delta_score", "higher_better": True},
                {"name": "parent_ion_fraction", "higher_better": True},
            ],
            "cv_params": [{"name": "fragmentation"}, {"name": "mass analyzer"}, {"name": "type"}],
        },
        "feature": {
            "additional_scores": [
                {"name": "andromeda_score", "higher_better": True},
                {"name": "andromeda_delta_score", "higher_better": True},
            ],
            "cv_params": [],
        },
        "pg": {
            "additional_scores": [{"name": "andromeda_score", "higher_better": True}],
            "cv_params": [{"name": "identification_type"}, {"name": "only_identified_by_site"}],
        },
    },
    "fragpipe": {
        "psm": {
            "additional_scores": [
                {"name": "msfragger_hyperscore", "higher_better": True},
                {"name": "msfragger_expectation", "higher_better": False},
                {"name": "parent_ion_fraction", "higher_better": True},
                {"name": "glycan_score", "higher_better": True},
                {"name": "glycan_qvalue", "higher_better": False},
            ],
            "cv_params": [{"name": "number_of_enzymatic_termini"}],
        },
        "feature": {
            "additional_scores": [],
            "cv_params": [],
        },
        "pg": {
            "additional_scores": [
                {"name": "protein_probability", "higher_better": True},
                {"name": "top_peptide_probability", "higher_better": True},
            ],
            "cv_params": [{"name": "sequence_length"}, {"name": "protein_existence"}],
        },
    },
    "proteome_discoverer": {
        "psm": {
            "additional_scores": [
                {"name": "rank", "higher_better": False},
                {"name": "sequest_xcorr", "higher_better": True},
                {"name": "comet_xcorr", "higher_better": True},
                {"name": "mascot_ions_score", "higher_better": True},
                {"name": "delta_score", "higher_better": True},
                {"name": "percolator_svm_score", "higher_better": True},
                {"name": "parent_ion_fraction", "higher_better": True},
            ],
            "cv_params": [{"name": "identifying_node"}, {"name": "psm_ambiguity"}],
        },
        "feature": {
            "additional_scores": [
                {"name": "rank", "higher_better": False},
                {"name": "sequest_xcorr", "higher_better": True},
                {"name": "comet_xcorr", "higher_better": True},
                {"name": "mascot_ions_score", "higher_better": True},
                {"name": "delta_score", "higher_better": True},
                {"name": "percolator_svm_score", "higher_better": True},
            ],
            "cv_params": [],
        },
        "pg": {
            "additional_scores": [{"name": "sum_pep_score", "higher_better": True}],
            "cv_params": [{"name": "sequence_length"}],
        },
    },
    "spectronaut": {
        "feature": {
            "additional_scores": [
                {"name": "pg_qvalue_runwise", "higher_better": False},
                {"name": "spectronaut_cscore", "higher_better": True},
                {"name": "spectronaut_pg_cscore", "higher_better": True},
                {"name": "spectronaut_pg_cscore_runwise", "higher_better": True},
                {"name": "spectronaut_ptm_assay_probability", "higher_better": True},
            ],
            "cv_params": [{"name": "predicted_irt"}],
        },
        "pg": {
            "additional_scores": [],
            "cv_params": [],
        },
        "fragment": {
            "additional_scores": [
                {"name": "spectronaut_interference_score", "higher_better": False},
                {"name": "spectronaut_log10_signal_to_noise", "higher_better": True},
                {"name": "spectronaut_measured_relative_intensity", "higher_better": True},
                {"name": "spectronaut_predicted_relative_intensity", "higher_better": True},
                {"name": "spectronaut_calibrated_mass_accuracy_ppm", "higher_better": False},
            ],
        },
    },
}
