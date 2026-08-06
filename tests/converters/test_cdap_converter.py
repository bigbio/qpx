"""CDAP converter integration tests.

Runs the full CDAP conversion pipeline once (module scope) and validates
that psm.parquet, feature.parquet, and pg.parquet are written with correct
schemas and plausible values.

By default uses the small bundled fixture at ``tests/examples/cdap/``.
Set environment variable ``CDAP_TEST_DATA_DIR`` to a full study directory
(e.g. ``/data/shenyufei/Bigbio_data/CPTAC/stage2/PDC000227``) for
full-scale validation.
"""

import os
from pathlib import Path

import pyarrow.parquet as pq
import pytest

_FIXTURE_DIR = Path(__file__).resolve().parent.parent / "examples" / "cdap"
_PSM_DIR = Path(os.environ.get("CDAP_TEST_DATA_DIR", str(_FIXTURE_DIR)))
_PREFIX = "cdap_test"


# ---------------------------------------------------------------------------
# Module-scoped fixture: run conversion once, share outputs across tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def converted_output(tmp_path_factory):
    """Run CDAP conversion once for all tests in this module."""
    from qpx.converters.cdap import CdapConverter

    if not _PSM_DIR.exists() or not list(_PSM_DIR.glob("*.psm")):
        pytest.skip(f"Test data not found: {_PSM_DIR}")

    output = tmp_path_factory.mktemp("cdap_output")

    converter = CdapConverter(max_memory="8GB", max_cpus=24)
    converter.convert(
        psm_dir=_PSM_DIR,
        output_folder=output,
        output_prefix=_PREFIX,
        project_accession="PDC000227",
    )

    return output


@pytest.fixture(scope="module")
def psm_table(converted_output):
    path = converted_output / f"{_PREFIX}.psm.parquet"
    if not path.exists():
        pytest.skip("psm.parquet was not produced")
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def feature_table(converted_output):
    path = converted_output / f"{_PREFIX}.feature.parquet"
    if not path.exists():
        pytest.skip("feature.parquet was not produced")
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def pg_table(converted_output):
    path = converted_output / f"{_PREFIX}.pg.parquet"
    if not path.exists():
        pytest.skip("pg.parquet was not produced")
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# PSM conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestCdapPsmConversion:
    """Validate psm.parquet output from CDAP conversion."""

    def test_file_exists(self, converted_output):
        assert (converted_output / f"{_PREFIX}.psm.parquet").exists()

    def test_has_rows(self, psm_table):
        assert psm_table.num_rows > 0

    def test_sequence_values_nonempty(self, psm_table):
        for seq in psm_table.column("sequence").to_pylist():
            assert isinstance(seq, str) and len(seq) > 0

    def test_charge_in_range(self, psm_table):
        for charge in psm_table.column("charge").to_pylist():
            assert 1 <= charge <= 15, f"Charge out of range: {charge}"

    def test_calculated_mz_positive(self, psm_table):
        for mz in psm_table.column("calculated_mz").to_pylist():
            if mz is not None:
                assert mz > 0, f"Invalid calculated_mz: {mz}"

    def test_run_names_nonempty(self, psm_table):
        for name in psm_table.column("run_file_name").to_pylist():
            assert isinstance(name, str) and len(name) > 0

    def test_decoy_ratio_reasonable(self, psm_table):
        decoys = sum(1 for x in psm_table.column("is_decoy").to_pylist() if x)
        ratio = decoys / psm_table.num_rows
        assert ratio <= 0.15, f"Decoy ratio too high: {ratio:.2%}"

    def test_schema_validation(self, psm_table):
        from qpx.core.data import PsmSchema

        errors = PsmSchema.validate(psm_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Feature conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestCdapFeatureConversion:
    """Validate feature.parquet output from CDAP conversion."""

    def test_file_exists(self, converted_output):
        assert (converted_output / f"{_PREFIX}.feature.parquet").exists()

    def test_has_rows(self, feature_table):
        assert feature_table.num_rows > 0

    def test_sequence_values_nonempty(self, feature_table):
        for seq in feature_table.column("sequence").to_pylist():
            assert isinstance(seq, str) and len(seq) > 0

    def test_intensities_nonnegative(self, feature_table):
        for row_ints in feature_table.column("intensities").to_pylist():
            if row_ints is None:
                continue
            for entry in row_ints:
                assert entry["intensity"] >= 0, f"Negative intensity: {entry}"

    def test_intensities_have_tmt_labels(self, feature_table):
        """At least some features should have TMT10 labels."""
        labels_seen = set()
        for row_ints in feature_table.column("intensities").to_pylist()[:1000]:
            if row_ints is None:
                continue
            for entry in row_ints:
                labels_seen.add(entry["label"])
        assert any("TMT" in lab for lab in labels_seen), f"No TMT labels found: {labels_seen}"

    def test_anchor_protein_nonempty(self, feature_table):
        for ap in feature_table.column("anchor_protein").to_pylist():
            assert isinstance(ap, str) and len(ap) > 0

    def test_schema_validation(self, feature_table):
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# PG conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestCdapPgConversion:
    """Validate pg.parquet output from CDAP conversion."""

    def test_file_exists(self, converted_output):
        assert (converted_output / f"{_PREFIX}.pg.parquet").exists()

    def test_has_rows(self, pg_table):
        assert pg_table.num_rows > 0

    def test_anchor_protein_nonempty(self, pg_table):
        for ap in pg_table.column("anchor_protein").to_pylist():
            assert isinstance(ap, str) and len(ap) > 0

    def test_intensities_nonnegative(self, pg_table):
        # pg is flattened since 1.1: scalar intensity column (one row per label).
        for intensity in pg_table.column("intensity").to_pylist():
            if intensity is None:
                continue
            assert intensity >= 0, f"Negative intensity: {intensity}"

    def test_pg_accessions_contain_anchor(self, pg_table):
        """The anchor protein must appear in pg_accessions."""
        anchors = pg_table.column("anchor_protein").to_pylist()
        pg_accs = pg_table.column("pg_accessions").to_pylist()
        for anchor, accs in zip(anchors, pg_accs):
            acc_names = [a["accession"] if isinstance(a, dict) else a for a in (accs or [])]
            assert anchor in acc_names, f"Anchor {anchor} not in pg_accessions {acc_names[:5]}"

    def test_schema_validation(self, pg_table):
        from qpx.core.data import PgSchema

        errors = PgSchema.validate(pg_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Ontology & metadata tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestCdapMetadata:
    """Validate ontology, provenance, and dataset parquet files."""

    def test_ontology_exists_and_nonempty(self, converted_output):
        path = converted_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written")
        table = pq.read_table(str(path))
        assert table.num_rows > 0

    def test_provenance_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.provenance.parquet"
        assert path.exists()
        table = pq.read_table(str(path))
        assert table.num_rows > 0

    def test_provenance_does_not_claim_sample_run(self, converted_output):
        """CDAP does not write sample/run, so provenance must not claim them."""
        path = converted_output / f"{_PREFIX}.provenance.parquet"
        table = pq.read_table(str(path))
        all_views = {view for views in table.column("output_views").to_pylist() for view in (views or [])}
        assert "sample" not in all_views
        assert "run" not in all_views

    def test_dataset_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.dataset.parquet"
        assert path.exists()
        table = pq.read_table(str(path))
        assert table.num_rows > 0


# ---------------------------------------------------------------------------
# LFQ compatibility tests
# ---------------------------------------------------------------------------


def test_cdap_lfq_precursor_area_uses_lfq_label(tmp_path):
    """CDAP LFQ PrecursorArea must be emitted with an LFQ label."""
    psm_dir = tmp_path / "lfq_psm"
    psm_dir.mkdir()
    psm_path = psm_dir / "sample_lfq.psm"
    psm_path.write_text(
        "\t".join(
            [
                "FileName",
                "ScanNum",
                "QueryPrecursorMz",
                "OriginalPrecursorMz",
                "PrecursorError(ppm)",
                "QueryCharge",
                "OriginalCharge",
                "PrecursorScanNum",
                "PrecursorArea",
                "PrecursorRelAb",
                "RTAtPrecursorHalfElution",
                "PeptideSequence",
                "AmbiguousMatch",
                "Protein",
                "DeNovoScore",
                "MSGFScore",
                "Evalue",
                "Qvalue",
                "PepQvalue",
                "PrecursorPurity",
                "FractionDecomposition",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "run_lfq.raw",
                "1001",
                "500.25",
                "500.25",
                "1.2",
                "2",
                "2",
                "1000?",
                "12345.6",
                "0.1",
                "321.0",
                "PEPTIDEK",
                "0",
                "NP_000001.1(pre=K,post=A)",
                "50",
                "40",
                "0.001",
                "0.0",
                "0.0",
                "99.0,99.0",
                "99.0,99.0",
            ]
        )
        + "\n"
    )

    from qpx.converters.cdap.feature_adapter import CdapFeatureAdapter

    output_path = tmp_path / "lfq.feature.parquet"
    with CdapFeatureAdapter() as adapter:
        adapter.convert(psm_dir=str(psm_dir), output_path=str(output_path))

    table = pq.read_table(str(output_path))
    labels = {entry["label"] for row in table.column("intensities").to_pylist() for entry in (row or [])}
    assert labels == {"LFQ"}


# ---------------------------------------------------------------------------
# Peptidoform unit tests
# ---------------------------------------------------------------------------


_MIN_PSM_HEADER = [
    "FileName",
    "ScanNum",
    "QueryPrecursorMz",
    "OriginalPrecursorMz",
    "PrecursorError(ppm)",
    "QueryCharge",
    "OriginalCharge",
    "PrecursorScanNum",
    "PrecursorArea",
    "PrecursorRelAb",
    "RTAtPrecursorHalfElution",
    "PeptideSequence",
    "AmbiguousMatch",
    "Protein",
    "DeNovoScore",
    "MSGFScore",
    "Evalue",
    "Qvalue",
    "PepQvalue",
    "PrecursorPurity",
    "FractionDecomposition",
]


def _write_psm(psm_dir, name, rows):
    """Write a minimal CDAP ``.psm`` file from a list of dict rows."""
    psm_dir.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(_MIN_PSM_HEADER)]
    for row in rows:
        lines.append("\t".join(str(row.get(col, "")) for col in _MIN_PSM_HEADER))
    (psm_dir / name).write_text("\n".join(lines) + "\n")


def _base_row(**overrides):
    row = {
        "FileName": "run_a.raw",
        "ScanNum": "1000",
        "QueryPrecursorMz": "500.25",
        "OriginalPrecursorMz": "500.25",
        "PrecursorError(ppm)": "1.0",
        "QueryCharge": "2",
        "OriginalCharge": "2",
        "PrecursorScanNum": "999",
        "PrecursorArea": "10000",
        "PrecursorRelAb": "0.1",
        "RTAtPrecursorHalfElution": "300.0",
        "PeptideSequence": "PEPTIDEK",
        "AmbiguousMatch": "0",
        "Protein": "NP_000001.1(pre=K,post=A)",
        "DeNovoScore": "50",
        "MSGFScore": "40",
        "Evalue": "0.001",
        "Qvalue": "0.0",
        "PepQvalue": "0.0",
        "PrecursorPurity": "99.0,99.0",
        "FractionDecomposition": "99.0,99.0",
    }
    row.update(overrides)
    return row


def test_cdap_feature_representative_from_single_psm(tmp_path):
    """Representative scan/mz/rt must all come from ONE best PSM (#250).

    Independent ``arg_min`` per column produced a Frankenstein record on
    tied/NULL Evalue. Here three PSMs share (peptidoform, charge, run); the
    winning (scan, observed_mz, rt) triple must equal exactly one input row.
    """
    from qpx.converters.cdap.feature_adapter import CdapFeatureAdapter

    rows = [
        _base_row(Evalue="", ScanNum="100", QueryPrecursorMz="400.0", RTAtPrecursorHalfElution="5.0"),
        _base_row(Evalue="0.2", ScanNum="200", QueryPrecursorMz="500.0", RTAtPrecursorHalfElution="6.0"),
        _base_row(Evalue="0.2", ScanNum="300", QueryPrecursorMz="600.0", RTAtPrecursorHalfElution="7.0"),
    ]
    _write_psm(tmp_path / "psm", "frank.psm", rows)

    out = tmp_path / "frank.feature.parquet"
    with CdapFeatureAdapter() as adapter:
        adapter.convert(psm_dir=str(tmp_path / "psm"), output_path=str(out))

    table = pq.read_table(str(out))
    assert table.num_rows == 1
    scan = table.column("scan").to_pylist()[0]
    observed_mz = table.column("observed_mz").to_pylist()[0]
    rt = table.column("rt").to_pylist()[0]
    input_triples = {(100, 400.0, 5.0), (200, 500.0, 6.0), (300, 600.0, 7.0)}
    got = (scan[0], observed_mz, rt)
    assert got in input_triples, f"Frankenstein: {got} matches no input PSM"
    # Deterministic best: Evalue 0.2 tie broken by lowest scan -> row 200.
    assert got == (200, 500.0, 6.0)


def test_cdap_feature_chemmod_calc_mz_is_null_not_observed(tmp_path):
    """An unparseable CHEMMOD must leave calculated_mz NULL/0, not fall back
    to the observed precursor m/z (#250)."""
    from qpx.converters.cdap.feature_adapter import CdapFeatureAdapter

    # +123.456 is not in the CDAP mass table -> CHEMMOD -> unparseable calc mass.
    rows = [_base_row(PeptideSequence="PEP+123.456TIDEK", OriginalPrecursorMz="777.7", QueryPrecursorMz="777.7")]
    _write_psm(tmp_path / "psm", "chemmod.psm", rows)

    out = tmp_path / "chemmod.feature.parquet"
    with CdapFeatureAdapter() as adapter:
        adapter.convert(psm_dir=str(tmp_path / "psm"), output_path=str(out))

    table = pq.read_table(str(out))
    assert table.num_rows == 1
    calc = table.column("calculated_mz").to_pylist()[0]
    # NULL calc mass is stored as 0.0; it must NOT be the measured 777.7.
    assert not calc, f"CHEMMOD calc_mz leaked observed m/z: {calc}"


def test_cdap_pg_distinct_groups_sharing_leader(tmp_path):
    """Two distinct protein groups sharing a leading protein must NOT collapse
    into one pg (#250). ``P1;P2`` and ``P1;P3`` share leader ``P1``."""
    from qpx.converters.cdap.feature_adapter import CdapFeatureAdapter
    from qpx.converters.cdap.pg_adapter import CdapPgAdapter

    rows = [
        _base_row(PeptideSequence="PEPTIDEA", Protein="P1(pre=K,post=A);P2(pre=K,post=A)"),
        _base_row(PeptideSequence="PEPTIDEB", Protein="P1(pre=K,post=A);P3(pre=K,post=A)"),
    ]
    _write_psm(tmp_path / "psm", "leader.psm", rows)

    feat = tmp_path / "leader.feature.parquet"
    with CdapFeatureAdapter() as adapter:
        adapter.convert(psm_dir=str(tmp_path / "psm"), output_path=str(feat))

    pg = tmp_path / "leader.pg.parquet"
    with CdapPgAdapter() as adapter:
        adapter.convert(feature_path=str(feat), output_path=str(pg))

    table = pq.read_table(str(pg))
    memberships = set()
    for accs in table.column("pg_accessions").to_pylist():
        names = frozenset(a["accession"] if isinstance(a, dict) else a for a in (accs or []))
        memberships.add(names)
    assert frozenset({"P1", "P2"}) in memberships
    assert frozenset({"P1", "P3"}) in memberships
    # Not collapsed to a single leader-keyed group.
    assert frozenset({"P1", "P2"}) != frozenset({"P1", "P3"})
    assert len(memberships) == 2


class TestCdapPeptidoform:
    """Unit tests for the ProForma conversion logic."""

    def test_nterm_label_strip(self):
        from qpx.converters.cdap.peptidoform import strip_label_tags

        stripped = strip_label_tags("+229.163PEPTIDEK+229.163")
        assert stripped == "PEPTIDEK"

    def test_internal_ptm(self):
        from qpx.converters.cdap.peptidoform import cdap_to_proforma

        result = cdap_to_proforma("+229.163PEPTM+15.995IDEK+229.163")
        assert "UNIMOD:" in result or "+15.995" in result

    def test_no_labels(self):
        from qpx.converters.cdap.peptidoform import strip_label_tags

        stripped = strip_label_tags("PEPTIDEK")
        assert stripped == "PEPTIDEK"

    def test_plain_sequence_roundtrip(self):
        from qpx.converters.cdap.peptidoform import cdap_to_proforma

        result = cdap_to_proforma("PEPTIDEK")
        assert result == "PEPTIDEK"
