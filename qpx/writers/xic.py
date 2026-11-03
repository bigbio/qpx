"""XicWriter — writes xic.parquet files."""

from qpx.core.data import XicSchema
from qpx.writers.base import BaseWriter


class XicWriter(BaseWriter):
    """Write validated xic.parquet files and declare their trace selection."""

    _schema_class = XicSchema

    def __init__(
        self,
        path,
        *,
        selection: str | None = None,
        window: str | None = None,
        extraction_window_seconds: float | None = None,
        qvalue_threshold: float | None = None,
        transition_format: str | None = None,
        smoothed: bool | None = None,
        zero_points: str | None = None,
        **kwargs,
    ):
        super().__init__(path, **kwargs)
        # File-level trace selection. Per-row cv_params would repeat these on
        # every trace; they describe how the whole file was extracted.
        if selection is not None:
            self._file_metadata[b"xic_selection"] = selection.encode()
        if window is not None:
            self._file_metadata[b"xic_window"] = window.encode()
        if extraction_window_seconds is not None:
            self._file_metadata[b"xic_extraction_window_seconds"] = str(extraction_window_seconds).encode()
        if qvalue_threshold is not None:
            self._file_metadata[b"xic_qvalue_threshold"] = str(qvalue_threshold).encode()
        if transition_format is not None:
            self._file_metadata[b"transition_format"] = transition_format.encode()
        if smoothed is not None:
            self._file_metadata[b"xic_smoothed"] = str(smoothed).lower().encode()
        if zero_points is not None:
            self._file_metadata[b"xic_zero_points"] = zero_points.encode()
        self._file_metadata[b"xic_rt_unit"] = b"seconds"
