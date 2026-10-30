"""QPX specification version and file-format compatibility guards.

Two distinct versions live in this project and must never be conflated:

* ``qpx._version.__version__`` — the **package** version (hatch-vcs, derived from
  the git tag; ``0.0.x`` in a dirty checkout, ``1.0.x`` on a release tag). It
  identifies the *software build* and moves on every release.
* :data:`QPX_SPEC_VERSION` — the **specification / on-disk format** version
  declared in ``docs/spec/versioning.md``. It moves only when the data model
  changes and is what gets stamped into the Parquet footer as ``qpx_version``.

Stamping the package version as ``qpx_version`` (the historical bug) meant a
file could never advertise the spec version a reader checks against.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

# The current QPX specification (on-disk format) version. Keep in sync with
# docs/spec/versioning.md ("Current Version"). Two-component {major}.{minor};
# no patch component.
QPX_SPEC_VERSION = "1.3"


class QpxVersionError(ValueError):
    """Raised when a QPX file's format version is incompatible with this reader."""


def parse_spec_version(value: str | None) -> tuple[int, int] | None:
    """Parse a ``{major}.{minor}`` spec-version string into an ``(int, int)`` tuple.

    Extra components (a stray patch segment) are ignored; a missing minor is
    treated as ``0``. Returns ``None`` (never raises) for anything that is not a
    parseable numeric spec version — ``None``, ``""``, ``"v1.1"``, ``"abc"`` —
    leaving it to the caller to decide what an unrecognisable version means.
    """
    if not value:
        return None
    parts = value.strip().split(".")
    try:
        major = int(parts[0])
        minor = int(parts[1]) if len(parts) > 1 else 0
    except (ValueError, IndexError):
        return None
    return (major, minor)


# --- File-format load guards ------------------------------------------------

# The 1.1 spec re-keyed the PG view from the scalar ``run_file_name`` onto
# ``grouped_runs`` (list<string>). Files written before 1.1 lack that column, so
# every grouped-runs query fails deep inside DuckDB with an opaque
# "column grouped_runs does not exist". We detect the old layout up front and
# raise something a user can act on.
_PG_MIN_VERSION = (1, 1)


def check_pg_file_compatible(file_path: str | Path) -> None:
    """Guard a PG parquet file against the pre-1.1 (scalar ``run_file_name``) layout.

    Reads the ``qpx_version`` key and the column list from the Parquet footer.
    Raises :class:`QpxVersionError` with an actionable message when the file
    predates the 1.1 ``grouped_runs`` re-key — either because its stamped
    ``qpx_version`` is below 1.1, or because it still carries a scalar
    ``run_file_name`` column with no ``grouped_runs`` column.

    A file whose footer cannot be read (missing/corrupt) is left alone: the
    normal read path surfaces that error. No-op for a well-formed 1.1+ file.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    try:
        pf = pq.ParquetFile(str(file_path))
    except (OSError, ValueError, pa.ArrowException):
        return  # let the regular reader report an unreadable/missing file

    schema = pf.schema_arrow
    meta = schema.metadata or {}
    raw_version = meta.get(b"qpx_version")
    version_str = raw_version.decode() if raw_version else None
    check_pg_columns_compatible(schema.names, file_path, version_str)


def check_pg_columns_compatible(
    columns: Iterable[str],
    source: str | Path,
    version_str: str | None = None,
) -> None:
    """Column-based variant of :func:`check_pg_file_compatible`.

    Use when the Parquet footer cannot be read directly with PyArrow — notably
    S3 datasets, where PyArrow ``ParquetFile`` cannot open a DuckDB-style glob
    (``s3://…/*.pg.parquet``). There the column list comes from the query engine
    after the file is registered (``DESCRIBE``). Same guard, same error, so the
    S3 path fails fast with the actionable message instead of an opaque binder
    error downstream.
    """
    _guard_pg_layout(str(source), set(columns), version_str)


def _guard_pg_layout(source: str, columns: set[str], version_str: str | None) -> None:
    """Raise :class:`QpxVersionError` if ``columns``/``version_str`` describe a pre-1.1 PG layout."""
    has_grouped = "grouped_runs" in columns
    has_scalar_run = "run_file_name" in columns
    has_flat_intensity = "label" in columns and "intensity" in columns
    parsed_version = parse_spec_version(version_str)
    too_old_version = parsed_version is not None and parsed_version < _PG_MIN_VERSION

    # A grouped-runs PG still needs the flattened scalar label/intensity columns.
    # Without them every 1.1 consumer query fails, regardless of whether the old
    # intensities list is also present.
    if has_grouped and not has_flat_intensity:
        found = version_str or "unknown (incomplete 1.1 layout)"
        raise QpxVersionError(
            "Incompatible PG (protein-group) file: "
            f"'{source}'\n"
            f"  file qpx_version : {found}\n"
            f"  reader requires  : >= {_PG_MIN_VERSION[0]}.{_PG_MIN_VERSION[1]} "
            f"(this build stamps {QPX_SPEC_VERSION})\n"
            "  breaking change  : QPX 1.1 flattened the PG 'intensities' list into "
            "scalar 'label' + 'intensity' columns (one row per label); this file "
            "does not expose both scalar columns, so pg queries cannot run.\n"
            "  fix              : regenerate this PG file with qpx >= 1.1. See "
            "docs/spec/versioning.md (changelog 1.1)."
        )

    # Old layout: no grouped_runs and only the scalar run_file_name, or an
    # explicitly pre-1.1 stamp on a file still missing grouped_runs.
    if not has_grouped and (has_scalar_run or too_old_version):
        found = version_str or "unknown (pre-1.1, no qpx_version stamp)"
        raise QpxVersionError(
            "Incompatible PG (protein-group) file: "
            f"'{source}'\n"
            f"  file qpx_version : {found}\n"
            f"  reader requires  : >= {_PG_MIN_VERSION[0]}.{_PG_MIN_VERSION[1]} "
            f"(this build stamps {QPX_SPEC_VERSION})\n"
            "  breaking change  : QPX 1.1 re-keyed the PG view from the scalar "
            "'run_file_name' onto 'grouped_runs' (list<string>); this file still "
            "has 'run_file_name' and no 'grouped_runs' column, so grouped-runs "
            "queries cannot run.\n"
            "  fix              : regenerate this PG file with qpx >= 1.1, or "
            "migrate it by wrapping each scalar run_file_name into a single-element "
            "grouped_runs list. See docs/spec/versioning.md (changelog 1.1)."
        )


# Feature intensities have been keyed by ``label`` since QPX 1.0. Files from the
# quantms.io era key each entry by ``channel`` (next to ``sample_accession``),
# which no feature query reads any more.


def check_feature_file_compatible(file_path: str | Path) -> None:
    """Guard a feature parquet file against the pre-1.0 ``intensities[].channel`` layout.

    A file whose footer cannot be read (missing/corrupt) is left alone: the
    normal read path surfaces that error.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    try:
        schema = pq.read_schema(str(file_path))
    except (OSError, ValueError, pa.ArrowException):
        return  # let the regular reader report an unreadable/missing file
    check_feature_schema_compatible(schema, file_path)


def check_feature_schema_compatible(schema, source: str | Path) -> None:
    """Raise :class:`QpxVersionError` if ``schema`` (a ``pyarrow.Schema``) keys feature intensities by ``channel``.

    Takes the schema rather than a path so S3 datasets, whose glob PyArrow cannot
    open, can pass the schema of the registered DuckDB view.
    """
    import pyarrow as pa

    if "intensities" not in schema.names:
        return
    intensities = schema.field("intensities").type
    # Any other layout is reported by the schema checks.
    is_list = pa.types.is_list(intensities) or pa.types.is_large_list(intensities)
    if not (is_list and pa.types.is_struct(intensities.value_type)):
        return
    keys = {field.name for field in intensities.value_type}
    if "channel" in keys and "label" not in keys:
        raise QpxVersionError(
            "Incompatible feature file: "
            f"'{source}'\n"
            "  breaking change  : QPX 1.0 keys each feature intensity by 'label'; "
            "this file keys it by 'channel', the pre-1.0 (quantms.io) layout.\n"
            f"  fix              : regenerate this feature file with qpx >= {QPX_SPEC_VERSION}. "
            "See docs/spec/versioning.md."
        )
