"""PSI-MS controlled vocabulary term constants.

Centralizes CV accession strings used across QPX converters and core modules,
avoiding duplication and providing a single source of truth.
"""

# ---------------------------------------------------------------------------
# Cross-linking (mzIdentML 1.3 XL-MS)
# ---------------------------------------------------------------------------
CV_XL_DONOR = "MS:1002509"  # crosslink donor
CV_XL_ACCEPTOR = "MS:1002510"  # crosslink acceptor / receiver
CV_CROSSLINK_SII = "MS:1002511"  # crosslink spectrum identification item
CV_LOOPLINK_SII = "MS:1003329"  # looplink spectrum identification item
CV_NONCOVALENT_SII = "MS:1003331"  # noncovalently associated peptides SII

# ---------------------------------------------------------------------------
# Identification metadata
# ---------------------------------------------------------------------------
CV_PEPTIDE_UNIQUE = "MS:1001363"  # peptide unique to one protein
CV_PEPTIDE_SHARED = "MS:1001175"  # peptide shared in multiple proteins
CV_PEPTIDE_PASSES_THRESHOLD = "MS:1002500"  # peptide passes threshold

# ---------------------------------------------------------------------------
# Score direction indicators
# ---------------------------------------------------------------------------
CV_HIGHER_BETTER = "MS:1002108"  # higher score better
CV_LOWER_BETTER = "MS:1002109"  # lower score better

# ---------------------------------------------------------------------------
# Score parent terms (OBO hierarchy)
# ---------------------------------------------------------------------------
SCORE_PARENT_IDS = frozenset(
    {
        "MS:1001143",  # PSM-level search engine specific statistic
        "MS:1001153",  # search engine specific score
        "MS:1002347",  # PSM-level identification statistic
        "MS:1002363",  # search engine specific score for proteins
        "MS:1002368",  # search engine specific score for protein groups
        "MS:1002701",  # PSM-level result list statistic
        "MS:1002906",  # search engine specific score for proteoforms
    }
)

# ---------------------------------------------------------------------------
# Scores -- specific score accessions
# ---------------------------------------------------------------------------
CV_PEP = "MS:1001493"  # posterior error probability (PEP)
CV_PEP_GLOBAL = "MS:1002352"  # PSM-level global FDR (global PEP)

# ---------------------------------------------------------------------------
# Site localization scores
# ---------------------------------------------------------------------------
CV_MODIFICATION_INDEX = "MS:1002504"  # modification index
# Parameters of a Modification or SearchModification that do not name the modification.
MODIFICATION_NON_IDENTITY_ACCESSIONS = frozenset(
    {
        CV_MODIFICATION_INDEX,
        "MS:1003347",  # UNIMOD derivative code
        "MS:1003390",  # crosslinker cleavage characteristics
        "MS:1003392",  # search modification id
        "MS:1003393",  # search modification id ref
    }
)
# Identification based on several spectra, and its PEP: values start with the identification's id.
MULTI_SPECTRA_ACCESSIONS = frozenset({"MS:1003332", "MS:1003336"})

# Per-site probabilities, the ones among the terms below.
SITE_PROBABILITY_ACCESSIONS = frozenset(
    {
        "MS:1001971",  # phosphoRS site probability
        "MS:1001980",  # MaxQuant:Phospho (STY) Probabilities
        "MS:1001982",  # MaxQuant:P-site localization probability
        "MS:1003147",  # PTMProphet probability
    }
)

# Per-site localization scores on a SpectrumIdentificationItem: the PSM-level
# statistics whose value follows MS:1002505,
# MOD_INDEX:SCORE:POSITION[|POSITION...]:PASS_THRESHOLD.
SITE_LOCALIZATION_ACCESSIONS = SITE_PROBABILITY_ACCESSIONS | {
    "MS:1001969",  # phosphoRS score
    "MS:1001974",  # DeBunker:score
    "MS:1001978",  # MSQuant:PTM-score
    "MS:1001979",  # MaxQuant:PTM Score
    "MS:1001983",  # MaxQuant:PTM Delta Score
    "MS:1001985",  # Ascore
    "MS:1001986",  # H-Score
    "MS:1002507",  # modification rescoring:false localization rate
    "MS:1002536",  # D-Score
    "MS:1002537",  # MD-Score
    "MS:1003291",  # Luciphor deltaScore
}

# ---------------------------------------------------------------------------
# Spectrum metadata
# ---------------------------------------------------------------------------
CV_SCAN_START_TIME = "MS:1000016"  # scan start time

# ---------------------------------------------------------------------------
# Derived sets (convenience groupings for converter logic)
# ---------------------------------------------------------------------------
SKIP_SCORE_ACCESSIONS = frozenset(
    {
        CV_PEPTIDE_UNIQUE,
        CV_PEPTIDE_SHARED,
        CV_PEPTIDE_PASSES_THRESHOLD,
        CV_CROSSLINK_SII,
        CV_LOOPLINK_SII,
        CV_NONCOVALENT_SII,
    }
)
