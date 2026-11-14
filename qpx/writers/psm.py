"""PsmWriter — writes psm.parquet files."""

import pyarrow as pa

from qpx.core.data import PsmSchema
from qpx.writers.base import BaseWriter
from qpx.writers.derived import PSM_SCORE_COLUMNS, fill_score_columns


class PsmWriter(BaseWriter):
    _schema_class = PsmSchema

    def _derive_columns(self, table: pa.Table) -> pa.Table:
        return fill_score_columns(table, PSM_SCORE_COLUMNS)
