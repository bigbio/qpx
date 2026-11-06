"""The label vocabulary of the feature and pg views (docs/spec/intensities.md#label-vocabulary)."""

from __future__ import annotations

import pyarrow as pa

from qpx.core.data.validation import ValidationIssue, _issue
from qpx.core.sql import sql_build

_MAX_LABELS_NAMED = 5
_LABEL_QUERIES = {
    "feature": "SELECT DISTINCT i.label FROM (SELECT UNNEST(intensities) AS i FROM ($query)) WHERE i.label IS NOT NULL",
    "pg": "SELECT DISTINCT label FROM ($query) WHERE label IS NOT NULL",
}


def canonical_labels() -> set[str]:
    """The canonical channel labels: the sdrf-pipelines channel map."""
    from sdrf_pipelines.converters.channel_map import CHANNEL_LABELS

    return {label for channels in CHANNEL_LABELS.values() for label in channels.values()}


def _has_labels(structure: str, schema: pa.Schema) -> bool:
    """Whether the view has labels in the current layout; the type check reports another one."""
    if structure == "pg":
        return "label" in schema.names
    if structure != "feature" or "intensities" not in schema.names:
        return False
    intensities = schema.field("intensities").type
    return (
        pa.types.is_list(intensities)
        and pa.types.is_struct(intensities.value_type)
        and intensities.value_type.get_field_index("label") >= 0
    )


def label_vocabulary_issues(engine, query: str, structure: str, schema: pa.Schema) -> list[ValidationIssue]:
    """Warn about labels outside the vocabulary; one in another case never joins its sample."""
    if not _has_labels(structure, schema):
        return []
    labels = {row[0] for row in engine.connection.sql(sql_build(_LABEL_QUERIES[structure], query=query)).fetchall()}
    canonical = canonical_labels()
    unknown = sorted(labels - canonical)
    if not unknown:
        return []
    by_case = {label.casefold(): label for label in canonical}
    named = [
        f"{label} (canonical: {by_case[label.casefold()]})" if label.casefold() in by_case else label
        for label in unknown[:_MAX_LABELS_NAMED]
    ]
    message = f"{len(unknown)} label(s) outside the label vocabulary: {', '.join(named)}"
    return [_issue(structure, "label_vocabulary", "warning", "label", message)]
