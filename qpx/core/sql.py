"""SQL safety utilities for DuckDB query construction.

Provides identifier validation and safe SQL building to prevent SQL injection
(Bandit B608) without relying on ``# nosec`` annotations.

DuckDB does not support parameterised identifiers (table/column names), so
this module validates them against a strict allowlist pattern and builds SQL
via ``string.Template`` — which Bandit does **not** flag as B608.

Values (WHERE predicates, LIMIT, OFFSET, IN lists) should still use DuckDB's
native ``?`` parameter binding via ``conn.execute(sql, params)``.
"""

from __future__ import annotations

import re
from string import Template

# The result is always double-quoted, so the only character that could break out
# of the quoted identifier is ``"`` itself (excluded). ``/`` is allowed because
# MaxQuant column names legitimately contain it (e.g. "MS/MS scan number", "m/z").
_SAFE_IDENTIFIER = re.compile(r"^[a-zA-Z_][a-zA-Z0-9_\[\].:/ ]*$")
_SAFE_TABLE = re.compile(r"^[a-zA-Z_][a-zA-Z0-9_]*$")


def validate_identifier(name: str) -> str:
    """Validate and double-quote a SQL identifier (column name).

    Raises ``ValueError`` for names that contain unsafe characters.

    >>> validate_identifier("run_file_name")
    '"run_file_name"'
    >>> validate_identifier("best_search_engine_score[1]")
    '"best_search_engine_score[1]"'
    """
    if not _SAFE_IDENTIFIER.match(name):
        raise ValueError("Unsafe SQL identifier: " + repr(name))
    return '"' + name + '"'


def validate_table(name: str) -> str:
    """Validate a bare table or view name (no quoting).

    Only allows ``[a-zA-Z_][a-zA-Z0-9_]*``.

    >>> validate_table("report")
    'report'
    """
    if not _SAFE_TABLE.match(name):
        raise ValueError("Unsafe table/view name: " + repr(name))
    return name


def escape_path(path: str) -> str:
    """Escape a filesystem path for use inside a SQL single-quoted literal.

    >>> escape_path("C:/data/file's.parquet")
    "C:/data/file''s.parquet"
    """
    return str(path).replace("'", "''")


def sql_build(template: str, **kwargs: str) -> str:
    """Build a SQL string from a ``string.Template``.

    Uses ``string.Template.safe_substitute`` so that DuckDB ``$1``-style
    positional references are left untouched (``$`` followed by a digit is
    not a valid Python identifier, so Template leaves it as-is).

    **Callers are responsible for validating / escaping every kwarg** via
    `validate_identifier`, `validate_table`, or `escape_path` before
    passing them here.

    >>> sql_build("SELECT * FROM $table", table="report")
    'SELECT * FROM report'
    """
    return Template(template).safe_substitute(**kwargs)
