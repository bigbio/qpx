"""FragmentWriter — writes fragment.parquet files."""

from qpx.core.data import FragmentSchema
from qpx.writers.base import BaseWriter


class FragmentWriter(BaseWriter):
    """Write schema-validated fragment.parquet files."""

    _schema_class = FragmentSchema
