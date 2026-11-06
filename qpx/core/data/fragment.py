"""Fragment data structure — fragment-level quantities of DIA features."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema

FragmentSchema = load_schema("fragment")


class Fragment(BaseStructure):
    """One row per fragment of a feature, channel and run."""

    _schema_class = FragmentSchema
