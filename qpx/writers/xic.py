"""XicWriter — writes xic.parquet files."""

from qpx.core.data import XicSchema
from qpx.writers.base import BaseWriter


class XicWriter(BaseWriter):
    """Write validated xic.parquet files."""

    _schema_class = XicSchema
