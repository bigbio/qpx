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
    """Recognize known contaminant accession markers, ignoring case.

    Preserve the legacy ``CONTAM`` substring rule. The quantms ``Cont_`` marker
    must start a bare accession or follow a leading UniProt ``sp|``/``tr|`` prefix,
    as in ``sp|Cont_Q7SIH1|A2MG_BOVIN``; unrelated ``CONT`` text is not a marker.
    """
    accession = str(accession).upper()
    return "CONTAM" in accession or accession.startswith(("CONT_", "SP|CONT_", "TR|CONT_"))


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
# MaxQuant helpers
# ------------------------------------------------------------------


def mq_flag_to_bool(val) -> bool:
    """Convert MaxQuant '+' flag to boolean."""
    return str(val).strip() == "+"
