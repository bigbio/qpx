"""Protein data structure — what the sequence database says of each protein."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema

ProteinSchema = load_schema("protein")


class Protein(BaseStructure):
    """One row per protein named in the dataset, with its name, organism and sequence."""

    _schema_class = ProteinSchema
