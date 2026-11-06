"""Opt-in semantic validation of peptidoforms, modification positions and m/z values."""

import pytest

from qpx.converters.ptm import PROTON_MASS
from qpx.core.data.feature import Feature
from qpx.core.data.psm import PSM
from qpx.core.semantic_validation import peptidoform_problem
from qpx.writers import FeatureWriter, PsmWriter
from tests.conftest import make_feature_record, make_psm_record

PEPTIDEK_MASS = peptidoform_problem("PEPTIDEK", "PEPTIDEK")[1]
PEPTIDEK_MZ = (PEPTIDEK_MASS + 2 * PROTON_MASS) / 2


@pytest.mark.parametrize(
    ("peptidoform", "sequence", "delta"),
    [
        ("[UNIMOD:1]-PEPTIDEK", "PEPTIDEK", 42.010565),
        ("PEPTIDEK-[UNIMOD:2]", "PEPTIDEK", -0.984016),
        ("PEPT[+79.9663]IDEK", "PEPTIDEK", 79.9663),
        ("PEPTIDEK[UNIMOD:1842|+54.0106]", "PEPTIDEK", 54.0106),
        ("[+129.0804]?PEPTIDEK", "PEPTIDEK", 129.0804),
    ],
)
def test_canonical_peptidoforms_and_their_mass(peptidoform, sequence, delta):
    """Canonical tags resolve to their mass shift."""
    problem, mass = peptidoform_problem(peptidoform, sequence)
    assert problem is None
    assert mass == pytest.approx(PEPTIDEK_MASS + delta, abs=1e-4)


@pytest.mark.parametrize(
    ("peptidoform", "sequence", "problem"),
    [
        ("PEPTIDEK[Custom (K)]", "PEPTIDEK", None),
        ("PEPM[UNIMOD:35", "PEPM", "malformed"),
        ("PEPM[[UNIMOD:35]]", "PEPM", "malformed"),
        ("PEP-TIDEK", "PEPTIDEK", "malformed"),
        ("pepTIDEK", "PEPTIDEK", "malformed"),
        ("PEPTM[UNIMOD:35]IDEK", "PEPTIDEK", "sequence"),
        ("PEPM[UNIMOD:999999]K", "PEPMK", "accession"),
        ("PEPM[MOD:00719]K", "PEPMK", "accession"),
        ("PEPM[UniMod:35]K", "PEPMK", "accession"),
    ],
)
def test_peptidoform_problems(peptidoform, sequence, problem):
    """Malformed strings, a sequence mismatch and non-Unimod accessions are problems; a name has no mass."""
    found, mass = peptidoform_problem(peptidoform, sequence)
    assert found == problem
    assert mass is None


def _checks(result):
    return {issue.check: issue for issue in result.issues}


def test_semantic_checks_are_opt_in(tmp_path):
    """Masses, positions and peptidoforms are checked only with semantic=True."""
    path = tmp_path / "a.feature.parquet"
    bad_mz = make_feature_record(calculated_mz=PEPTIDEK_MZ + 50, observed_mz=PEPTIDEK_MZ + 50)
    misplaced = make_feature_record(charge=3, calculated_mz=None, observed_mz=None)
    misplaced["modifications"] = [
        {"name": "Oxidation", "accession": "UNIMOD:35", "positions": [{"position": 42, "amino_acid": "M", "scores": None}]}
    ]
    off_observed = make_feature_record(calculated_mz=PEPTIDEK_MZ, observed_mz=PEPTIDEK_MZ * 1.0001)
    off_observed["rt"] = 300.0
    with FeatureWriter(path) as w:
        w.write_batch([bad_mz, misplaced, off_observed])
    assert not {"calculated_mz_mismatch", "modification_position", "mass_error"} & set(
        _checks(Feature.from_file(path).validate())
    )
    checks = _checks(Feature.from_file(path).validate(strict=True, semantic=True))
    assert checks["calculated_mz_mismatch"].severity == "error"
    assert checks["calculated_mz_mismatch"].message.startswith("1 row(s)")
    assert checks["modification_position"].message.startswith("1 modification position(s)")
    assert checks["mass_error"].severity == "warning"
    assert checks["mass_error"].message.startswith("1 row(s)")


def test_multiplet_rows_skip_the_mass_error(tmp_path):
    """A multiplet row's m/z describe its lightest partner, so its mass error is not checked."""
    path = tmp_path / "a.feature.parquet"
    record = make_feature_record(calculated_mz=PEPTIDEK_MZ, observed_mz=PEPTIDEK_MZ + 4)
    record["identified_labels"] = ["SILAC heavy"]
    with FeatureWriter(path) as w:
        w.write_batch([record])
    assert "mass_error" not in _checks(Feature.from_file(path).validate(semantic=True))


def test_cross_linked_psms_skip_the_mass_checks(tmp_path):
    """A cross-linked PSM's m/z include its partner peptide."""
    path = tmp_path / "a.psm.parquet"
    record = make_psm_record(calculated_mz=PEPTIDEK_MZ + 300, observed_mz=PEPTIDEK_MZ + 300)
    record["cross_links"] = [
        {
            "xl_type": "inter",
            "partner_sequence": "WAKALPGFR",
            "partner_peptidoform": "WAKALPGFR",
            "donor_position": 8,
            "acceptor_position": 3,
            "linker_name": "DSSO",
            "linker_accession": "UNIMOD:1842",
            "linker_mass": 158.0038,
        }
    ]
    with PsmWriter(path, scan_format="scan") as w:
        w.write_batch([record])
    checks = _checks(PSM.from_file(path).validate(semantic=True))
    assert "calculated_mz_mismatch" not in checks
    assert "mass_error" not in checks
