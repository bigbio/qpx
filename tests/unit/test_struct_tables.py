"""Keep fixed layouts in sync with converters' ordered source mappings."""

from types import SimpleNamespace

import pytest

from qpx.converters.diann.feature_sql import _CV_MAPPINGS as DIANN_CV
from qpx.converters.diann.feature_sql import _SCORE_MAPPINGS as DIANN_SCORES
from qpx.converters.diann.pg_adapter import _ADDITIONAL_SCORES as DIANN_PG_SCORES
from qpx.converters.diann.pg_adapter import _CHANNEL_SCORES as DIANN_CHANNEL_SCORES
from qpx.converters.fragpipe.pg_adapter import FragPipePgAdapter
from qpx.converters.fragpipe.psm_adapter import FragPipePsmAdapter
from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter
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


def _annotations(scores: list[dict], cv: list[dict] | None) -> tuple[list, list]:
    return [(score["score_name"], score["higher_better"]) for score in scores], [entry["cv_name"] for entry in cv or []]


def _fragpipe_psm():
    columns = ("Hyperscore", "Expectation", "Purity", "Glycan Score", "Glycan q-value", "Compensation Voltage")
    row = {**dict.fromkeys(columns, 0.5), "Number of Missed Cleavages": 0, "Number of Enzymatic Termini": 2}
    scores, _, cv = FragPipePsmAdapter._psm_annotations(row)
    return _annotations(scores, cv)


def _fragpipe_pg():
    row = {"Protein Length": 300, "Protein Existence": "1:Experimental evidence at protein level"}
    row.update({"Protein Probability": 0.99, "Top Peptide Probability": 0.98})
    annotations = FragPipePgAdapter._pg_annotations(row)
    return _annotations(annotations["additional_scores"], annotations["cv_params"])


def _maxquant_psm():
    row = {"Score": 80.0, "Delta score": 40.0, "PIF": 0.9, "Fragmentation": "HCD", "Mass analyzer": "FTMS", "Type": "MULTI-MSMS"}
    return _annotations(*MaxQuantPsmAdapter._psm_annotations(SimpleNamespace(_resolved={}), row))


@pytest.mark.parametrize(
    ("creator", "view", "annotate"),
    [("fragpipe", "psm", _fragpipe_psm), ("fragpipe", "pg", _fragpipe_pg), ("maxquant", "psm", _maxquant_psm)],
)
def test_layout_matches_the_scores_written_row_by_row(creator, view, annotate):
    """Converters without a source mapping give every layout name, in order and with the layout's direction."""
    scores, cv = annotate()
    table = FIXED_STRUCT_TABLES[creator][view]
    assert scores == [(entry["name"], entry["higher_better"]) for entry in table["additional_scores"]]
    assert cv == [entry["name"] for entry in table.get("cv_params", [])]
