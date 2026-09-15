from __future__ import annotations

from typing import Optional


def safe_float(val) -> Optional[float]:
    """Convert a value to float, returning None for missing/invalid values.

    Handles None, empty strings, the literal ``"null"``, NaN, and
    anything that cannot be cast to float.
    """
    if val is None or val == "" or val == "null":
        return None
    try:
        f = float(val)
        return f if f == f else None  # NaN check
    except (ValueError, TypeError):
        return None


# ------------------------------------------------------------------
# UniProt accession helpers
# ------------------------------------------------------------------


def parse_uniprot_id(entry: str) -> tuple[str, str]:
    """Split a UniProt-style ``db|ACCESSION|NAME`` id into ``(accession, name)``.

    ``sp|P12345|PROT_HUMAN`` -> ``("P12345", "PROT_HUMAN")``. With only two
    pipe-fields the accession doubles as the name; with none the whole ``entry``
    is both. Used by the FragPipe protein-field parser.
    """
    parts = entry.split("|")
    if len(parts) >= 3:
        return parts[1], parts[2]
    if len(parts) == 2:
        return parts[1], parts[1]
    return entry, entry


def uniprot_entry_name(entry: str) -> Optional[str]:
    """The entry name of a ``db|ACCESSION|NAME`` id, or None when it has none.

    Unlike :func:`parse_uniprot_id`, a bare accession does not double as its own
    name: ``pg_names`` should be null rather than repeat ``pg_accessions``.
    """
    parts = str(entry).split("|")
    return parts[2] if len(parts) >= 3 and parts[2] else None


def is_contaminant_accession(accession) -> bool:
    """True for a contaminant protein, by the ``CONTAM`` accession marker.

    The quantms/sdrf-pipelines contaminant database prefixes its entries
    ``CONTAM_`` (``sp|CONTAM_P19001|CONTAM_K1C19_MOUSE``). One rule shared by the
    converters, so the same protein is not a contaminant in one and unflagged in
    another (bigbio/qpx#300).
    """
    return "CONTAM" in str(accession).upper()


def strip_uniprot_prefix(accession: str) -> str:
    """Strip a leading ``sp|``/``tr|`` UniProt db prefix, returning the accession.

    ``sp|P55011|S12A2_HUMAN`` -> ``P55011``; ``tr|A0A..|..._HUMAN`` -> ``A0A..``;
    ids without that prefix (bare accessions, ``CON__``/``REV__`` decoys) are
    returned unchanged. Used by the MaxQuant adapters, whose FASTA prefixing is
    the only case that should be stripped.
    """
    if accession and accession.startswith(("sp|", "tr|")):
        parts = accession.split("|")
        if len(parts) >= 2:
            return parts[1]
    return accession


# ------------------------------------------------------------------
# CV term column-name helpers
# ------------------------------------------------------------------


def cv_column_names(cv_term: str, suffix: str) -> tuple[str, str]:
    """Return (lowercase, uppercase) opt_global column names for a CV term."""
    cv_code = cv_term.split(":")[1]
    return (
        f"opt_global_cv_ms:{cv_code}_{suffix}",
        f"opt_global_cv_MS:{cv_code}_{suffix}",
    )


def get_cv_value(row: dict, cv_term: str, suffix: str, default=None):
    """Get a value from a row dict trying both CV column name variants."""
    lo, hi = cv_column_names(cv_term, suffix)
    return row.get(lo, row.get(hi, default))


# ------------------------------------------------------------------
# MaxQuant helpers
# ------------------------------------------------------------------


def mq_flag_to_bool(val) -> bool:
    """Convert MaxQuant '+' flag to boolean."""
    return str(val).strip() == "+"
