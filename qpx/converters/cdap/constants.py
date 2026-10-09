"""CDAP converter constants -- channel definitions, modification masses, decoy prefix.

The CPTAC ``.psm`` tab-separated files produced by CDAP (Common Data Analysis
Pipeline) contain TMT/iTRAQ reporter ion intensities inline with PSM-level
identification data.  Constants in this module summarise:

    * Reporter-ion channel naming variants observed across all 91 CPTAC studies.
    * Modification delta-mass patterns and their UNIMOD accessions.
    * Decoy protein prefix used by CDAP (``XXX_``).

Source survey: ``/home/shenyufei/.windsurf/plans/cdap-psm-survey-77e878.md``.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Channel definitions (5 patterns observed across 91 studies)
# ---------------------------------------------------------------------------

# Reporter-ion channel column ordering. The list order matches the canonical
# mass-ascending order in which CDAP emits TMT/iTRAQ columns and is used to
# resolve channel labels against `comment[label]` values from SDRF files.
CHANNEL_DEFS: dict[str, list[str]] = {
    "TMT10": [
        "TMT10-126",
        "TMT10-127N",
        "TMT10-127C",
        "TMT10-128N",
        "TMT10-128C",
        "TMT10-129N",
        "TMT10-129C",
        "TMT10-130N",
        "TMT10-130C",
        "TMT10-131",
    ],
    "TMT11": [
        "TMT11-126C",
        "TMT11-127N",
        "TMT11-127C",
        "TMT11-128N",
        "TMT11-128C",
        "TMT11-129N",
        "TMT11-129C",
        "TMT11-130N",
        "TMT11-130C",
        "TMT11-131N",
        "TMT11-131C",
    ],
    # Legacy extra channels seen in 15 TMT11 studies that piloted TMTpro masses.
    "TMTXX": ["TMTXX-132N", "TMTXX-132C"],
    "TMT18": [
        "TMT18-126C",
        "TMT18-127N",
        "TMT18-127C",
        "TMT18-128N",
        "TMT18-128C",
        "TMT18-129N",
        "TMT18-129C",
        "TMT18-130N",
        "TMT18-130C",
        "TMT18-131N",
        "TMT18-131C",
        "TMT18-132N",
        "TMT18-132C",
        "TMT18-133N",
        "TMT18-133C",
        "TMT18-134N",
        "TMT18-134C",
        "TMT18-135N",
    ],
    "iTRAQ4": ["iTRAQ114", "iTRAQ115", "iTRAQ116", "iTRAQ117"],
}

# ---------------------------------------------------------------------------
# Modification masses -> UNIMOD accessions
# ---------------------------------------------------------------------------

# N-term/Lys label tag masses. Used to identify and strip the labelling
# modification when computing the stripped ``sequence`` field.  Values map the
# delta-mass token (as it appears in ``PeptideSequence``) to the UNIMOD id.
LABEL_TAG_MASSES: dict[str, int] = {
    "+229.163": 737,  # TMT6/10/11
    "+304.207": 730,  # TMTpro / TMT18
    "+144.102": 214,  # iTRAQ4
}

# Variable / fixed PTM masses found in the survey, in addition to the label
# tags above. Mass strings are kept as the literal token CDAP emits so that the
# parser can compare without round-tripping floats.
PTM_MASSES: dict[str, int] = {
    "+57.021": 4,  # Carbamidomethyl (fixed on C)
    "+79.966": 21,  # Phospho (S/T/Y)
    "+15.995": 35,  # Oxidation (M)
    "+42.011": 1,  # Acetyl (N-term/K)
    "+42.010": 1,  # Acetyl (mass-quantisation variant)
    "+114.043": 121,  # GlyGly (K) -- Ubiquitylation remnant
    "+0.984": 7,  # Deamidation (N/Q)
    "-17.027": 28,  # Gln->pyro-Glu (N-term Q)
}

# Combined lookup used by ``peptidoform.py``; preserves the label-vs-PTM split
# so that callers can choose to strip labels but keep PTMs.
ALL_MOD_MASSES: dict[str, int] = {**LABEL_TAG_MASSES, **PTM_MASSES}

# ---------------------------------------------------------------------------
# Decoy definitions
# ---------------------------------------------------------------------------

# CDAP marks decoy hits with this fixed prefix on every accession.
DECOY_PREFIX: str = "XXX_"
