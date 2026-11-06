"""XicWriter — writes xic.parquet files."""

from dataclasses import dataclass

from qpx.core.data import XicSchema
from qpx.writers.base import BaseWriter


@dataclass(frozen=True)
class XicExtraction:
    """How the traces of a file were selected and extracted; ``None`` leaves an entry out of the footer."""

    selection: str | None = None
    window: str | None = None
    extraction_window_seconds: float | None = None
    qvalue_threshold: float | None = None
    transition_format: str | None = None
    smoothed: bool | None = None
    zero_points: str | None = None


# Footer key of each XicExtraction field.
_FOOTER_KEYS = {
    "selection": b"xic_selection",
    "window": b"xic_window",
    "extraction_window_seconds": b"xic_extraction_window_seconds",
    "qvalue_threshold": b"xic_qvalue_threshold",
    "transition_format": b"transition_format",
    "smoothed": b"xic_smoothed",
    "zero_points": b"xic_zero_points",
}


class XicWriter(BaseWriter):
    """Write validated xic.parquet files and declare their trace selection."""

    _schema_class = XicSchema

    def __init__(self, path, *, extraction: XicExtraction | None = None, **kwargs):
        super().__init__(path, **kwargs)
        # File-level trace selection. Per-row cv_params would repeat these on
        # every trace; they describe how the whole file was extracted.
        for name, key in _FOOTER_KEYS.items():
            value = getattr(extraction, name, None)
            if value is not None:
                self._file_metadata[key] = (str(value).lower() if isinstance(value, bool) else str(value)).encode()
        self._file_metadata[b"xic_rt_unit"] = b"seconds"
