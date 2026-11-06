"""CDAP converter constants -- channel definitions, decoy prefix.

The CPTAC ``.psm`` tab-separated files produced by CDAP (Common Data Analysis
Pipeline) contain TMT/iTRAQ reporter ion intensities inline with PSM-level
identification data.  Constants in this module summarise:

    * Reporter-ion channel naming variants observed across all 91 CPTAC studies.
    * Decoy protein prefix used by CDAP (``XXX_``).

Modification delta masses are resolved to Unimod accessions by
:func:`qpx.converters.ptm.mass_to_unimod` (see ``peptidoform.py``).

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
# Decoy definitions
# ---------------------------------------------------------------------------

# CDAP marks decoy hits with this fixed prefix on every accession.
DECOY_PREFIX: str = "XXX_"
