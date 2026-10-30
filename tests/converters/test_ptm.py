"""Test shared PTM utilities (flat ptm.py)."""

import pytest

from qpx.converters.ptm import (
    UNIMOD_MASS,
    _normalize_peptidoform,
    build_proforma,
    compute_precursor_mz,
    from_proforma,
    mass_to_unimod,
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


class TestUnimodRegistry:
    """The curated Unimod table, the site-aware mass lookup and the name lookup."""

    def test_names_and_masses_match_unimod(self):
        import pyopenms

        from qpx.converters.ptm import UNIMOD_MODS

        db = pyopenms.ModificationsDB()
        reference = {}
        for index in range(db.getNumberOfModifications()):
            mod = db.getModification(index)
            reference.setdefault(mod.getUniModRecordId(), (mod.getId(), mod.getDiffMonoMass()))
        for uid, (name, mass, _) in UNIMOD_MODS.items():
            assert reference[uid][0] == name, uid
            assert abs(reference[uid][1] - mass) < 1e-5, uid

    @pytest.mark.parametrize(
        ("mass", "residue", "nterm", "expected"),
        [
            (229.1629, "K", False, 737),
            (304.2071, "K", False, 2016),
            (304.2054, "K", False, 730),
            (144.1021, None, True, 214),
            (-17.0265, "Q", True, 28),
            (-17.0265, "C", True, 385),
            (-17.0265, None, False, None),
            (-18.0106, "E", True, 27),
            (15.9949, "M", False, 35),
            (500.0, None, False, None),
        ],
    )
    def test_mass_to_unimod(self, mass, residue, nterm, expected):
        from qpx.converters.ptm import mass_to_unimod

        assert mass_to_unimod(mass, residue, nterm) == expected

    @pytest.mark.parametrize(
        ("name", "expected"),
        [
            ("GlyGly (K)", 121),
            ("Acetyl (Protein N-term)", 1),
            ("Deamidation (NQ)", 7),
            ("TMTpro (K)", 2016),
            ("iTRAQ8plex (N-term)", 730),
            ("Glu->pyro-Glu", 27),
            ("Label:13C(6)15N(2)", 259),
            ("Custom (K)", None),
        ],
    )
    def test_unimod_from_name(self, name, expected):
        from qpx.converters.ptm import unimod_from_name

        assert unimod_from_name(name) == expected

    def test_canonical_mod_tag(self):
        from qpx.converters.ptm import canonical_mod_tag

        assert canonical_mod_tag("UNIMOD:35", "Oxidation", 15.994915) == "UNIMOD:35"
        assert canonical_mod_tag("UNIMOD:35", "Oxidation", 54.0106) == "+54.0106"
        assert canonical_mod_tag("UNIMOD:99999", None, 54.0106) == "UNIMOD:99999|+54.0106"
        assert canonical_mod_tag("MS:1001460", "unknown modification", 79.96633, "S") == "UNIMOD:21"
        assert canonical_mod_tag(None, None, 1234.5678) == "+1234.5678"


class TestBracketedTagsKeepParentheses:
    """Parentheses inside a bracketed tag belong to the tag (no nested brackets)."""

    @pytest.mark.parametrize("peptidoform", ["PEPK[GlyGly (K)]R", "PEPK[Label:13C(6)15N(2)]R"])
    def test_inner_parentheses_untouched(self, peptidoform):
        assert _normalize_peptidoform(peptidoform) == peptidoform

    def test_positions_after_named_tag(self):
        from qpx.converters.ptm import from_proforma

        _, mods = from_proforma("AK[Acetyl (K)]PS[UNIMOD:21]TR", "AKPSTR")
        assert [(m["name"], m["positions"][0]["position"]) for m in mods] == [("Acetyl (K)", 2), ("UNIMOD:21", 4)]


@pytest.mark.parametrize(
    ("raw", "canonical"),
    [
        ("PEPTM[Oxidation]K", "PEPTM[UNIMOD:35]K"),
        ("[Acetyl]-PEPTIDE", "[UNIMOD:1]-PEPTIDE"),
        ("PEPTC(Carbamidomethyl)K", "PEPTC[UNIMOD:4]K"),
        ("PEPS[CHEMMOD:+79.966331]K", "PEPS[UNIMOD:21]K"),
        ("PEPK[+54.01062]R", "PEPK[+54.0106]R"),
        ("PEPK[UniMod:737]R", "PEPK[UNIMOD:737]R"),
        ("[-17.0265]-QPEPK", "[UNIMOD:28]-QPEPK"),
        ("PEPK[Custom (K)]R", "PEPK[Custom (K)]R"),
    ],
)
def test_normalize_peptidoform(raw, canonical):
    from qpx.converters.ptm import normalize_peptidoform

    assert normalize_peptidoform(raw) == canonical
