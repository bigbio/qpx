"""ProteinWriter — writes protein.parquet files."""

from qpx.core.data import ProteinSchema
from qpx.writers.base import BaseWriter


class ProteinWriter(BaseWriter):
    """Write schema-validated protein.parquet files."""

    _schema_class = ProteinSchema
