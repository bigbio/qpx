"""Test shared PTM utilities (flat ptm.py)."""

import pyopenms as oms
import pytest

from qpx.converters.ptm import (
    _PREFERRED_ON_TIE,
    UNIMOD_MASS,
    _normalize_peptidoform,
    build_proforma,
    compute_precursor_mz,
    format_mass_delta,
    from_proforma,
    mass_tag,
    mass_to_unimod,
    modification_tag,
    name_to_unimod,
)


class TestUnimodMassRegistry:
    def test_oxidation_mass(self):
        assert abs(UNIMOD_MASS[35] - 15.994915) < 0.0001

    def test_acetyl_mass(self):
        assert abs(UNIMOD_MASS[1] - 42.010565) < 0.0001

    def test_mass_to_unimod_oxidation(self):
        assert mass_to_unimod(15.9949) == 35

    def test_mass_to_unimod_unknown(self):
        assert mass_to_unimod(999.0) is None

    @pytest.mark.parametrize(
        ("accession", "mass"),
        [
            (27, -18.010565),  # Glu->pyro-Glu
            (214, 144.102063),  # iTRAQ4plex (the old table held the Glu->pyro-Glu mass)
            (259, 8.014199),  # Label:13C(6)15N(2)
            (267, 10.008269),  # Label:13C(6)15N(4)
            (312, 119.004099),  # Cysteinyl (the old table held the TMT6plex mass)
            (385, -17.026549),  # Ammonia-loss (the old table held the TMTpro mass)
            (730, 304.205360),  # iTRAQ8plex (the old table held the TMTpro mass)
            (737, 229.162932),  # TMT6plex
            (888, 140.094963),  # mTRAQ
            (2016, 304.207146),  # TMTpro
        ],
    )
    def test_masses_are_unimod(self, accession, mass):
        """Every accession carries its Unimod delta mass, as pyopenms' ModificationsDB has it."""
        assert UNIMOD_MASS[accession] == pytest.approx(mass, abs=1e-6)
        assert accession in UNIMOD_MASS
        assert UNIMOD_MASS.get(999_999) is None


class TestMassToUnimod:
    """Delta masses resolve to the Unimod accession defined at their residue or terminus."""

    @pytest.mark.parametrize(
        ("mass", "residue", "terminus", "accession"),
        [
            # Isobaric labels, on Lys and on the peptide N-terminus
            (229.1629, "K", None, 737),
            (229.1629, "P", "N", 737),
            (304.2071, "K", None, 2016),
            (304.2071, "P", "N", 2016),
            (144.1021, "K", None, 214),
            (144.1021, "P", "N", 214),
            (304.2054, "K", None, 730),
            (304.2054, "P", "N", 730),
            (140.0950, "K", None, 888),
            # Common modifications
            (15.9949, "M", None, 35),
            (79.9663, "S", None, 21),
            (79.9663, "T", None, 21),
            (79.9663, "Y", None, 21),
            (57.0215, "C", None, 4),
            (42.0106, "M", "N", 1),
            (42.0106, "K", None, 1),
            (0.9840, "N", None, 7),
            (0.9840, "Q", None, 7),
            (114.0429, "K", None, 121),
            (14.0157, "K", None, 34),
            (28.0313, "K", None, 36),
            (-0.9840, "K", "C", 2),
            (8.0142, "K", None, 259),
            (10.0083, "R", None, 267),
            # Equal masses at different sites
            (-17.0265, "Q", "N", 28),
            (-17.0265, "C", "N", 385),
            (-17.0265, "Q", None, 28),  # pyro-Glu reported on the Q instead of the N-terminus
            (-18.0106, "E", "N", 27),
        ],
    )
    def test_accession(self, mass, residue, terminus, accession):
        assert mass_to_unimod(mass, residue, terminus) == accession

    @pytest.mark.parametrize(
        ("mass", "residue", "terminus", "accession", "tied_with"),
        [
            (144.1021, "K", None, 214, 889),  # iTRAQ4plex, not mTRAQ:13C(3)15N(1)
            (144.1021, "P", "N", 214, 889),
            (0.9840, "N", None, 7, 621),  # Deamidated, not Asn->Asp
            (114.0429, "K", None, 121, 1290),  # GG, not Dicarbamidomethyl
            (57.0215, "K", None, 4, 1263),  # Carbamidomethyl (over-alkylation), not Gly
            (28.0313, "K", None, 36, 280),  # Dimethyl, not Ethyl
            (70.0419, "K", None, 1289, 253),  # Butyryl, not Crotonaldehyde
            (72.0211, "K", None, 2114, 378),  # Lactylation, not Carboxyethyl
            (98.0732, "K", None, 2121, 1873),  # Hexanoyl, not MesitylOxide
        ],
    )
    def test_isobaric_tie_goes_to_the_preferred_modification(self, mass, residue, terminus, accession, tied_with):
        """Equal masses at the same site resolve to the modification search engines use, not the lowest accession."""
        assert UNIMOD_MASS[tied_with] == pytest.approx(UNIMOD_MASS[accession], abs=1e-6)
        assert mass_to_unimod(mass, residue, terminus) == accession

    def test_preferred_accessions_are_in_this_unimod_release(self):
        assert _PREFERRED_ON_TIE <= set(UNIMOD_MASS)

    @pytest.mark.parametrize(
        ("mass", "residue", "terminus"),
        [
            (84.0575, "K", None),  # Isovaleryl 2108, 2-methylbutyryl 2109, Pentanoyl 2120
            (105.0215, "K", None),  # ICPL 365, NIC 697
            (-17.0265, None, "N"),  # Gln->pyro-Glu 28 on Q, Ammonia-loss 385 on C
        ],
    )
    def test_open_isobaric_tie_is_unresolved(self, mass, residue, terminus):
        """A tie without a preferred modification gives no accession, so the mass is written."""
        assert mass_to_unimod(mass, residue, terminus) is None
        assert mass_tag(mass, residue, terminus) == format_mass_delta(mass)

    def test_openms_custom_entries_are_not_unimod(self):
        """OpenMS' own ModificationsDB entries (CUSTOM0-9, PhosphoDecoy; record id 99900 and up) are left out."""
        assert all(accession < 99900 for accession in UNIMOD_MASS)
        assert 99900 not in UNIMOD_MASS
        assert 99913 not in UNIMOD_MASS
        assert mass_to_unimod(1.0078, "K") is None  # CUSTOM0, +1.00784 on every residue
        assert mass_to_unimod(79.9663, "A") is None  # PhosphoDecoy on A/L/G
        assert name_to_unimod("CUSTOM0") is None
        assert name_to_unimod("PhosphoDecoy", "A") is None

    @pytest.mark.parametrize("mass", [float("nan"), float("inf"), float("-inf")])
    def test_non_finite_mass(self, mass):
        assert mass_to_unimod(mass) is None
        assert mass_to_unimod(mass, "K") is None

    def test_mass_without_a_unimod_entry_at_the_site(self):
        """TMT6plex is not defined on Asp, so the mass stays unresolved there."""
        assert mass_to_unimod(229.1629, "D") is None
        assert mass_tag(229.1629, "D") == "+229.1629"
        assert mass_tag(229.1629, "K") == "UNIMOD:737"


class TestNameToUnimod:
    """Vendor modification names resolve through ModificationsDB names and synonyms."""

    @pytest.mark.parametrize(
        ("name", "residue", "terminus", "accession"),
        [
            ("Phospho (STY)", "S", None, 21),
            ("Phospho", "T", None, 21),
            ("Oxidation (M)", "M", None, 35),
            ("Carbamidomethyl (C)", "C", None, 4),
            ("Acetyl (Protein N-term)", "M", "N", 1),
            ("Deamidation (NQ)", "N", None, 7),
            ("Deamidated (NQ)", "Q", None, 7),
            ("GlyGly (K)", "K", None, 121),
            ("GG (K)", "K", None, 121),
            ("Gln->pyro-Glu", "Q", "N", 28),
            ("Glu->pyro-Glu", "E", "N", 27),
            ("Pyro-glu from E", "E", "N", 27),
            ("TMT6plex (K)", "K", None, 737),
            ("TMT10plex (N-term)", "P", "N", 737),
            ("TMTpro (K)", "K", None, 2016),
            ("TMT16plex (K)", "K", None, 2016),
            ("iTRAQ4plex (K)", "K", None, 214),
            ("iTRAQ8plex (N-term)", "P", "N", 730),
            ("Label:13C(6)15N(2) (K)", "K", None, 259),
            ("Hex(5)HexNAc(2)", "N", None, 137),
        ],
    )
    def test_accession(self, name, residue, terminus, accession):
        assert name_to_unimod(name, residue, terminus) == accession

    def test_unknown_name(self):
        assert name_to_unimod("SomeNewMod (K)", "K") is None


class TestModificationTag:
    def test_accession_tags_are_normalised(self):
        assert modification_tag("UniMod:35") == "UNIMOD:35"
        assert modification_tag("unimod:21", "S") == "UNIMOD:21"

    def test_signed_mass(self):
        assert modification_tag("+15.9949", "M") == "UNIMOD:35"
        assert modification_tag("-1.2345", "M") == "-1.2345"

    def test_unknown_name_has_no_square_brackets(self):
        assert modification_tag("Odd [mod]") == "Odd (mod)"


class TestLabelledPrecursorMz:
    """``compute_precursor_mz`` of labelled peptides equals pyopenms ``AASequence`` m/z."""

    @pytest.mark.parametrize(
        ("accession", "name"),
        [(737, "TMT6plex"), (2016, "TMTpro"), (214, "iTRAQ4plex"), (730, "iTRAQ8plex"), (888, "mTRAQ")],
    )
    @pytest.mark.parametrize("charge", [2, 3])
    def test_label_mz(self, accession, name, charge):
        peptidoform = f"[UNIMOD:{accession}]-PEPTIDEK[UNIMOD:{accession}]"
        expected = oms.AASequence.fromString(f"({name})PEPTIDEK({name})").getMZ(charge)
        assert compute_precursor_mz(peptidoform, charge) == pytest.approx(expected, abs=1e-4)


class TestBuildProforma:
    def test_no_mods(self):
        assert build_proforma("PEPTIDEK", []) == "PEPTIDEK"

    def test_internal_mod(self):
        assert build_proforma("PEPTMIDEK", [(5, "UNIMOD:35")]) == "PEPTM[UNIMOD:35]IDEK"

    def test_nterm_mod(self):
        assert build_proforma("PEPTIDEK", [(0, "UNIMOD:1")]) == "[UNIMOD:1]-PEPTIDEK"

    def test_empty_sequence(self):
        assert build_proforma("", []) == ""

    def test_cterm_mod(self):
        """A mod at ``len(seq)+1`` is the C-term ProForma slot (bigbio/qpx#251).

        Previously this position was silently dropped from ``peptidoform``
        while surviving in ``modifications`` -- an identity inconsistency.
        """
        assert build_proforma("PEPTIDEK", [(9, "UNIMOD:2")]) == "PEPTIDEK-[UNIMOD:2]"

    def test_nterm_and_cterm(self):
        assert build_proforma("PEPTIDEK", [(0, "UNIMOD:1"), (9, "UNIMOD:2")]) == "[UNIMOD:1]-PEPTIDEK-[UNIMOD:2]"

    def test_two_mods_same_position(self):
        """Two mods at the same position must not collapse (bigbio/qpx#251).

        ``dict(mods)`` used to keep only the last tag; both must be emitted.
        """
        assert build_proforma("PEPTIDEK", [(0, "UNIMOD:1"), (0, "UNIMOD:737")]) == "[UNIMOD:1][UNIMOD:737]-PEPTIDEK"

    def test_two_mods_same_residue(self):
        assert build_proforma("PEPTMIDEK", [(5, "UNIMOD:35"), (5, "CHEMMOD:+1.0")]) == "PEPTM[UNIMOD:35][CHEMMOD:+1.0]IDEK"


class TestFromProforma:
    def test_unimod_tag(self):
        meta = {"UNIMOD:35": ("Oxidation", ["M"], ["Anywhere"])}
        _, result = from_proforma("M[UNIMOD:35]PEPTIDEK", "MPEPTIDEK", meta=meta)
        assert result is not None
        assert len(result) == 1
        assert result[0]["name"] == "Oxidation"
        assert result[0]["accession"] == "UNIMOD:35"

    def test_no_mods(self):
        peptidoform, result = from_proforma("PEPTIDEK", "PEPTIDEK", meta=None)
        assert result is None
        assert peptidoform == "PEPTIDEK"

    def test_nterm(self):
        meta = {"UNIMOD:1": ("Acetyl", ["X"], ["N-term"])}
        _, result = from_proforma("[UNIMOD:1]-PEPTIDEK", "PEPTIDEK", meta=meta)
        assert result is not None
        assert result[0]["positions"][0]["position"] == 0

    def test_cterm(self):
        """Parse a C-terminal ProForma tag at the length-plus-one position."""
        meta = {"UNIMOD:2": ("Amidated", ["X"], ["C-term"])}
        _, result = from_proforma("PEPTIDEK-[UNIMOD:2]", "PEPTIDEK", meta=meta)
        assert result == [
            {
                "name": "Amidated",
                "accession": "UNIMOD:2",
                "positions": [{"position": 9, "amino_acid": None, "scores": None}],
            }
        ]

    def test_site_scores_attached(self):
        """Site scores are attached to the correct modification position."""
        site_scores = {
            1: [
                {
                    "score_name": "phosphors_site_probability",
                    "score_value": 0.95,
                    "higher_better": True,
                }
            ],
        }
        _, result = from_proforma(
            "M[UNIMOD:35]PEPTIDEK",
            "MPEPTIDEK",
            meta=None,
            site_scores=site_scores,
        )
        assert result is not None
        pos = result[0]["positions"][0]
        assert pos["position"] == 1
        assert pos["scores"] is not None
        assert pos["scores"][0]["score_value"] == 0.95

    def test_site_scores_none_when_not_provided(self):
        """Without site_scores, positions still have scores=None."""
        _, result = from_proforma("M[UNIMOD:35]PEPTIDEK", "MPEPTIDEK", meta=None)
        assert result is not None
        assert result[0]["positions"][0]["scores"] is None

    def test_site_scores_partial_positions(self):
        """Only positions in the site_scores dict get scores attached."""
        site_scores = {
            5: [
                {
                    "score_name": "phospho_prob",
                    "score_value": 0.8,
                    "higher_better": True,
                }
            ],
        }
        # Two mods: position 1 and position 5; only position 5 has scores
        _, result = from_proforma(
            "M[UNIMOD:35]PEPTS[UNIMOD:21]IDEK",
            "MPEPTSIDEK",
            meta=None,
            site_scores=site_scores,
        )
        assert result is not None
        # Find the phospho mod (position 5)
        for mod in result:
            for pos in mod["positions"]:
                if pos["position"] == 5:
                    assert pos["scores"] is not None
                    assert pos["scores"][0]["score_value"] == 0.8
                elif pos["position"] == 1:
                    assert pos["scores"] is None


class TestNormalizePeptidoform:
    """Validate mzTab parenthetical -> ProForma bracket conversion."""

    def test_no_parens(self):
        result = _normalize_peptidoform("PEPTIDEK")
        if result != "PEPTIDEK":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_simple_mod(self):
        result = _normalize_peptidoform("C(Carbamidomethyl)R")
        if result != "C[Carbamidomethyl]R":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_nested_parens(self):
        result = _normalize_peptidoform("C(Carbamidomethyl (C))R")
        if result != "C[Carbamidomethyl (C)]R":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_nterm_mod(self):
        result = _normalize_peptidoform(".(Acetyl)PEPTIDEK")
        if result != "[Acetyl]-PEPTIDEK":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_multiple_mods(self):
        result = _normalize_peptidoform(".(Acetyl)M(Oxidation)PEPTC(Carbamidomethyl)K")
        if result != "[Acetyl]-M[Oxidation]PEPTC[Carbamidomethyl]K":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_bracketed_name_with_parentheses_is_kept(self):
        """Parentheses inside a ProForma tag are part of the name, not a nested tag."""
        assert _normalize_peptidoform("PEPK[GlyGly (K)]R") == "PEPK[GlyGly (K)]R"
        assert _normalize_peptidoform("HEXN[Hex(5)HexNAc(2)]TK") == "HEXN[Hex(5)HexNAc(2)]TK"

    def test_bracketed_and_parenthesised_tags(self):
        assert _normalize_peptidoform("PEPK[Odd (K)]TIDES(Phospho)K") == "PEPK[Odd (K)]TIDES[Phospho]K"

    def test_two_nterm_tags_get_one_dash(self):
        """The ProForma ``-`` follows the last leading tag, not the first."""
        assert _normalize_peptidoform("[UNIMOD:1][Odd (K)]-PEPK") == "[UNIMOD:1][Odd (K)]-PEPK"
        assert _normalize_peptidoform("(Acetyl)(TMT6plex)PEPK") == "[Acetyl][TMT6plex]-PEPK"

    def test_unmatched_open_paren(self):
        result = _normalize_peptidoform("C(Broken")
        if result != "C(Broken":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_empty_string(self):
        result = _normalize_peptidoform("")
        if result != "":
            raise AssertionError(f"Unexpected: {result!r}")


class TestFromProformaMzTab:
    """Validate from_proforma with mzTab parenthetical input."""

    def test_mztab_simple(self):
        meta = {"UNIMOD:4": ("Carbamidomethyl", ["C"], ["Anywhere"])}
        _, result = from_proforma("C(UNIMOD:4)PEPTIDEK", "CPEPTIDEK", meta=meta)
        if result is None:
            raise AssertionError("Expected mods, got None")
        if result[0]["accession"] != "UNIMOD:4":
            raise AssertionError(f"Unexpected accession: {result[0]['accession']!r}")

    def test_mztab_nterm(self):
        meta = {"UNIMOD:1": ("Acetyl", ["X"], ["N-term"])}
        _, result = from_proforma("(UNIMOD:1)-PEPTIDEK", "PEPTIDEK", meta=meta)
        if result is None:
            raise AssertionError("Expected mods, got None")
        if result[0]["positions"][0]["position"] != 0:
            raise AssertionError(f"Unexpected position: {result[0]['positions'][0]['position']}")

    def test_mztab_nested_parens(self):
        # Normalization preserves inner parens as-is
        _, result = from_proforma("C(Carbamidomethyl (C))PEPTIDEK", "CPEPTIDEK", meta=None)
        if result is None:
            raise AssertionError("Expected mods, got None")
        mod_name = result[0]["name"]
        if "Carbamidomethyl (C)" not in mod_name:
            raise AssertionError(f"Inner parens lost: {mod_name!r}")


class TestComputePrecursorMz:
    def test_basic(self):
        mz = compute_precursor_mz("PEPTM(UniMod:35)IDE", 2)
        assert mz is not None
        assert mz > 0

    def test_empty(self):
        assert compute_precursor_mz("", 2) is None
        assert compute_precursor_mz("PEPTIDE", 0) is None
