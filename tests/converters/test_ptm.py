"""Test shared PTM utilities (flat ptm.py)."""

from qpx.converters.ptm import (
    UNIMOD_MASS,
    _normalize_peptidoform,
    build_proforma,
    compute_precursor_mz,
    from_proforma,
    mass_to_unimod,
    strip_modifications,
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


class TestStripModifications:
    """A named modification must not leak its own letters into the sequence.

    The naive ``re.sub(r"[^A-Z]", "", peptidoform.upper())`` upper-cases first,
    so ``[TMTpro]`` survives as ``TMTPRO`` and welds itself onto the peptide.
    Every TMT feature table written before this was fixed carries such
    sequences, and they match nothing in a digested FASTA -- iBAQ silently
    produced zero rows.
    """

    def test_nterm_and_cterm_tmtpro(self):
        """Tags on both termini are dropped, not welded to the sequence."""
        result = strip_modifications("[TMTpro]-GEDVETSK[TMTpro]")
        if result != "GEDVETSK":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_internal_named_mod(self):
        """An internal named tag is dropped along with the terminal ones."""
        result = strip_modifications("[TMTpro]-EEQQDDTVYM[Oxidation]GK[TMTpro]")
        if result != "EEQQDDTVYMGK":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_mztab_parenthetical_mod(self):
        """mzTab writes tags in parentheses; those count too."""
        result = strip_modifications("C(Carbamidomethyl)R")
        if result != "CR":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_numeric_mod_still_stripped(self):
        """Numeric tags never triggered the bug; guard the path that worked."""
        result = strip_modifications("PEPT[+57.02]IDEK")
        if result != "PEPTIDEK":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_unmodified_sequence_unchanged(self):
        """A bare sequence passes through untouched."""
        result = strip_modifications("PEPTIDEK")
        if result != "PEPTIDEK":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_empty_string(self):
        """An empty peptidoform yields an empty sequence, not an error."""
        result = strip_modifications("")
        if result != "":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_reagent_name_with_its_own_parentheses(self):
        """SILAC/dimethyl names close a paren before the tag itself ends.

        ``Label:13C(6)15N(2)`` is heavy lysine. A flat ``\\([^)]*\\)`` stops at
        the first ``)`` and leaks ``15N`` into the sequence -- the same defect
        this function exists to prevent, one level deeper.
        """
        for peptidoform in (
            ".(Label:13C(6)15N(2))PEPTIDEK",
            ".(Label:13C(6)15N(4))PEPTIDEK",
            ".(Dimethyl:2H(6)13C(2))PEPTIDEK",
            ".(Label:13C(4)15N(2)+GlyGly)PEPTIDEK",
        ):
            result = strip_modifications(peptidoform)
            if result != "PEPTIDEK":
                raise AssertionError(f"{peptidoform!r} gave {result!r}")

    def test_leading_dot_nterm_form(self):
        """mzTab writes N-term mods after a leading dot; real quantms output."""
        result = strip_modifications(".(TMT6plex)AAALAAAVAQDPAASGAPSS")
        if result != "AAALAAAVAQDPAASGAPSS":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_mod_letters_are_not_kept(self):
        """The exact regression: the old expression returned TMTPROAYAQGISR."""
        result = strip_modifications("[TMTpro]-AYAQGISR")
        if "TMTPRO" in result:
            raise AssertionError(f"Modification name leaked into sequence: {result!r}")
        if result != "AYAQGISR":
            raise AssertionError(f"Unexpected: {result!r}")


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
