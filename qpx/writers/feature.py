"""FeatureWriter — writes feature.parquet files."""

import logging
from collections.abc import Sequence

import numpy as np
import pyarrow as pa

from qpx.core.data import FeatureSchema
from qpx.core.data.identity import derive_id
from qpx.writers.base import BaseWriter
from qpx.writers.derived import FEATURE_SCORE_COLUMNS, fill_group_ids, fill_irt, fill_score_columns

_log = logging.getLogger(__name__)

# Composite fields FeatureWriter persists as float32; the others keep their Python values.
_FLOAT32_FIELDS = frozenset(field.name for field in FeatureSchema.get_arrow_schema() if field.type == pa.float32())


def derive_feature_ids(records: list[dict], composite: Sequence[str]) -> list[int]:
    """The ``feature_id`` FeatureWriter derives for each identified record from *composite*.

    Values are hashed as persisted, so a float32 field is rounded to float32 first.
    """
    unordered = tuple(index for index, name in enumerate(composite) if name in ("grouped_runs", "pg_accessions"))
    return [
        derive_id(
            [
                float(np.float32(rec.get(name))) if name in _FLOAT32_FIELDS and rec.get(name) is not None else rec.get(name)
                for name in composite
            ],
            unordered_list_indices=unordered,
        )
        for rec in records
    ]


class FeatureWriter(BaseWriter):
    _schema_class = FeatureSchema

    def __init__(self, *args, override_provided_ids: bool = True, **kwargs):
        """Derive identified Feature ids by default from the declared composite."""
        super().__init__(*args, override_provided_ids=override_provided_ids, **kwargs)
        self._unidentified_fallback_count = 0

    def _derive_columns(self, table: pa.Table) -> pa.Table:
        return fill_irt(fill_group_ids(fill_score_columns(table, FEATURE_SCORE_COLUMNS)))

    def _should_override_provided_id(self, index: int, composite_values: dict[str, list]) -> bool:
        """Use the natural composite only for identified Features."""
        return self._override_provided_ids and bool(composite_values["peptidoform"][index])

    def _identity_values(
        self,
        existing: list,
        composite_values: dict[str, list],
        composite: Sequence[str],
        cv_lists: list | None,
        id_field: str,
    ) -> list[int]:
        """Derive identified ids; for unidentified Features derive from the composite
        when the producer supplies no id, or namespace the producer id when it does."""
        peptidoforms = composite_values["peptidoform"]
        self._unidentified_fallback_count += sum(
            not peptidoform and existing[index] is None for index, peptidoform in enumerate(peptidoforms)
        )
        ids = super()._identity_values(existing, composite_values, composite, cv_lists, id_field)
        if not self._override_provided_ids:
            return ids

        for index, peptidoform in enumerate(peptidoforms):
            if peptidoform:
                continue
            provided = existing[index]
            if provided is None:
                continue
            # Unidentified Feature WITH a producer id: namespace it by run/unit so
            # tool-local ids cannot collide across runs, and stash the original.
            namespace_field = "run_file_name" if "run_file_name" in composite_values else "quantification_unit_id"
            ids[index] = derive_id([composite_values[namespace_field][index], provided])
            if cv_lists is not None:
                params = list(cv_lists[index] or [])
                provided_param = {"cv_name": f"provided_{id_field}", "cv_value": str(provided)}
                if provided_param not in params:
                    params.append(provided_param)
                cv_lists[index] = params
        return ids

    def _validate_identity_uniqueness(self) -> None:
        """Accept unidentified composite fallback only after full-file validation."""
        super()._validate_identity_uniqueness()

        if self._unidentified_fallback_count:
            composite = ", ".join(self._identity_composite or ())
            _log.warning(
                "%d unidentified Feature row(s) without a producer feature_id used identity_composite [%s]; "
                "full-file primary-key uniqueness was verified",
                self._unidentified_fallback_count,
                composite,
            )
