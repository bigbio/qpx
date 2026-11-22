"""Keep fixed layouts in sync with converters' ordered source mappings."""

import pytest

from qpx.converters.diann.feature_sql import _CV_MAPPINGS as DIANN_CV
from qpx.converters.diann.feature_sql import _SCORE_MAPPINGS as DIANN_SCORES
from qpx.converters.diann.pg_adapter import _ADDITIONAL_SCORES as DIANN_PG_SCORES
from qpx.converters.diann.pg_adapter import _CHANNEL_SCORES as DIANN_CHANNEL_SCORES
from qpx.converters.proteome_discoverer.reports import _SCORES as PD_SCORES
from qpx.converters.spectronaut.feature_adapter import _CV_MAPPINGS as SN_CV
from qpx.converters.spectronaut.feature_adapter import _SCORE_MAPPINGS as SN_SCORES
from qpx.converters.spectronaut.fragments import _SCORES as SN_FRAGMENT_SCORES
from qpx.core.struct_tables import FIXED_STRUCT_TABLES

_TYPED = {"precursor_qvalue", "global_qvalue"}
_PD_SCORES = [(name, higher) for _, name, higher in PD_SCORES]


@pytest.mark.parametrize(
    ("creator", "view", "scores", "cv"),
    [
        ("diann", "feature", [(n, h) for _, n, h in DIANN_SCORES if n not in _TYPED], [n for _, n in DIANN_CV]),
        ("diann", "pg", [(n, DIANN_CHANNEL_SCORES[n]) for n in DIANN_PG_SCORES], ["quantification_method"]),
        ("proteome_discoverer", "psm", [*_PD_SCORES, ("parent_ion_fraction", True)], ["identifying_node", "psm_ambiguity"]),
        ("proteome_discoverer", "feature", _PD_SCORES, []),
        ("spectronaut", "feature", [(n, h) for _, n, h in SN_SCORES if n not in _TYPED], [n for _, n in SN_CV]),
        ("spectronaut", "fragment", [(n, h) for _, n, h in SN_FRAGMENT_SCORES], []),
    ],
)
def test_layout_matches_converter_mappings(creator, view, scores, cv):
    """Check the field order and score direction against the source mappings."""
    table = FIXED_STRUCT_TABLES[creator][view]
    assert [(entry["name"], entry["higher_better"]) for entry in table["additional_scores"]] == scores
    assert [entry["name"] for entry in table.get("cv_params", [])] == cv
