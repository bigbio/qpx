"""Write the fragment view from a converter's fragment-level SQL."""

from __future__ import annotations

from collections.abc import Sequence

import pyarrow as pa

from qpx.writers.feature import FeatureWriter
from qpx.writers.fragment import FragmentWriter

_BATCH_ROWS = 200_000


def write_fragment_view(conn, sql: str, output_path: str, feature_composite: Sequence[str], **writer_options) -> int:
    """Stream *sql* into the fragment view; each row's feature_id is hashed as the feature view hashes it."""
    rows = 0
    with FragmentWriter(output_path, **writer_options) as writer:
        for batch in conn.sql(sql).to_arrow_reader(_BATCH_ROWS):
            table = pa.Table.from_batches([batch])
            if not table.num_rows:
                continue
            identified = FeatureWriter.derive_ids(table.select(list(feature_composite)), feature_composite)
            table = table.append_column("feature_id", identified.column("feature_id"))
            # A composite column that is not a fragment field (DIA-NN is_decoy, is_library_decoy) only feeds the id.
            table = table.drop_columns([name for name in feature_composite if name not in writer.arrow_schema.names])
            writer.write_table(writer.align_table_to_schema(table))
            rows += table.num_rows
    return rows
