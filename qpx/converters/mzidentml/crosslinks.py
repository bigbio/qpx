"""mzIdentML crosslink groups and donor/acceptor positions."""

from __future__ import annotations

from typing import Optional

from qpx.core.cv_terms import CV_CROSSLINK_SII, CV_LOOPLINK_SII, CV_NONCOVALENT_SII


def _build_looplink(
    alpha_pep: dict,
    linker_info: dict,
) -> Optional[list[dict]]:
    """Build cross_links for a looplink (intra-peptide) PSM."""
    xl_mods = alpha_pep.get("xl_mods", {})
    if not xl_mods:
        return None

    links = []
    for xl_value, mods in xl_mods.items():
        donor = next((m for m in mods if m["role"] == "donor"), None)
        acceptor = next((m for m in mods if m["role"] == "acceptor"), None)
        if not donor or not acceptor:
            continue

        linker = linker_info.get(xl_value, linker_info.get("_default", {}))
        links.append(
            {
                "xl_type": "intra",
                "partner_sequence": None,
                "partner_peptidoform": None,
                "donor_position": donor["position"],
                "acceptor_position": acceptor["position"],
                "linker_name": linker.get("name", "unknown"),
                "linker_accession": linker.get("accession"),
                "linker_mass": donor["mass_delta"],
            }
        )

    return links or None


def _build_interlink(
    alpha_pep: dict,
    beta_sii: Optional[dict],
    peptides: dict,
    linker_info: dict,
    xl_type: str,
) -> Optional[list[dict]]:
    """Build cross_links for an inter-peptide or noncovalent PSM."""
    if beta_sii is None:
        return None

    beta_pep = peptides.get(beta_sii["peptide_ref"], {})
    beta_sequence = beta_pep.get("sequence", "")
    beta_peptidoform = beta_pep.get("peptidoform") or beta_sequence

    # Find donor/acceptor positions from xl_mods
    alpha_xl = alpha_pep.get("xl_mods", {})
    beta_xl = beta_pep.get("xl_mods", {})

    # ``donor_position`` is the site on the alpha peptide and ``acceptor_position`` the one on
    # its partner, whichever of the two the file names donor; the linker mass sits on the donor.
    links = []
    for xl_value in sorted(alpha_xl.keys() | beta_xl.keys()):
        a_mods = alpha_xl.get(xl_value, [])
        b_mods = beta_xl.get(xl_value, [])
        donor = next((m for m in a_mods + b_mods if m["role"] == "donor"), None)
        if not a_mods or donor is None:
            continue

        linker = linker_info.get(xl_value, linker_info.get("_default", {}))

        links.append(
            {
                "xl_type": xl_type,
                "partner_sequence": beta_sequence,
                "partner_peptidoform": beta_peptidoform,
                "donor_position": a_mods[0]["position"],
                "acceptor_position": b_mods[0]["position"] if b_mods else None,
                "linker_name": linker.get("name", "unknown"),
                "linker_accession": linker.get("accession"),
                "linker_mass": donor["mass_delta"],
            }
        )

    # Noncovalently associated peptides share no linked site.
    if not links and xl_type == "noncovalent":
        links.append(
            {
                "xl_type": "noncovalent",
                "partner_sequence": beta_sequence,
                "partner_peptidoform": beta_peptidoform,
                "donor_position": None,
                "acceptor_position": None,
                "linker_name": "noncovalent",
                "linker_accession": None,
                "linker_mass": 0.0,
            }
        )

    return links or None


def _sii_link_group(sii) -> tuple:
    """Read crosslink group and association flags from one SII."""
    xl_group = None
    is_looplink = False
    is_noncovalent = False

    for cv in sii["cv_params"]:
        if cv["accession"] == CV_CROSSLINK_SII:
            xl_group = cv["value"]
        elif cv["accession"] == CV_LOOPLINK_SII:
            is_looplink = True
        elif cv["accession"] == CV_NONCOVALENT_SII:
            is_noncovalent = True
            xl_group = cv["value"]

    return xl_group, is_looplink, is_noncovalent


def _paired_sii_groups(paired_groups, xl_type) -> list[dict]:
    """Pair the SIIs of each group across peptides, retaining orphan items as regular PSMs.

    A group can hold more than one pair (OpenPepXL writes alpha and beta for the light
    and the heavy spectrum), so the items of one peptide pair up with those of the other.
    Items of a single peptide are a homodimer and pair up in file order.
    """
    groups = []
    for paired in paired_groups.values():
        by_peptide: dict[str, list[dict]] = {}
        for sii in paired:
            by_peptide.setdefault(sii["peptide_ref"], []).append(sii)
        sides = list(by_peptide.values())
        alphas, betas = (sides[0], sides[1]) if len(sides) == 2 else (paired[0::2], paired[1::2])
        groups.extend({"xl_type": xl_type, "alpha": alpha, "beta": beta} for alpha, beta in zip(alphas, betas))
        # Orphan SII -- treat as regular
        unpaired = alphas[len(betas) :] + betas[len(alphas) :]
        groups.extend({"xl_type": "none", "alpha": sii} for sii in unpaired)
    return groups


def _group_siis(siis: list[dict]) -> list[dict]:
    """Group SII elements into classified cross-link groups.

    Returns list of dicts: ``{xl_type, alpha, beta?}``
    where ``xl_type`` is one of: ``"none"``, ``"inter"``, ``"intra"``, ``"noncovalent"``.
    """
    groups = []

    # Index SIIs by their crosslink group value
    xl_groups: dict[str, list[dict]] = {}  # group_value → [sii, ...]
    looplink_siis = []
    noncovalent_groups: dict[str, list[dict]] = {}
    regular_siis = []

    for sii in siis:
        xl_group, is_looplink, is_noncovalent = _sii_link_group(sii)

        if is_looplink:
            looplink_siis.append(sii)
        elif is_noncovalent and xl_group:
            noncovalent_groups.setdefault(xl_group, []).append(sii)
        elif xl_group:
            xl_groups.setdefault(xl_group, []).append(sii)
        else:
            regular_siis.append(sii)

    # Regular PSMs
    for sii in regular_siis:
        groups.append({"xl_type": "none", "alpha": sii})

    # Looplinks (intra-peptide)
    for sii in looplink_siis:
        groups.append({"xl_type": "intra", "alpha": sii})

    # Inter-peptide crosslinks (paired SIIs)
    groups.extend(_paired_sii_groups(xl_groups, "inter"))
    # Noncovalent associations
    groups.extend(_paired_sii_groups(noncovalent_groups, "noncovalent"))

    return groups
