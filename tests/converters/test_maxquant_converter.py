"""MaxQuant converter integration tests using small example dataset.

Runs the full MaxQuant conversion pipeline once (module scope) and validates
that psm.parquet, feature.parquet, and ontology.parquet are written with
correct schemas and plausible values.
"""

from pathlib import Path

import pyarrow.parquet as pq
import pytest

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples" / "maxquant" / "maxquant_simple"

_MSMS = EXAMPLES_DIR / "msms.txt"
_EVIDENCE = EXAMPLES_DIR / "evidence.txt"
_SDRF = EXAMPLES_DIR / "sdrf.tsv"

_PREFIX = "maxquant_test"


def test_pg_fails_before_writing_when_experiment_has_no_run_mapping(tmp_path):
    """An unmapped MaxQuant Experiment must not be swallowed as skipped PG rows."""
    from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text("Protein IDs\tMajority protein IDs\tIntensity experiment_1\nP12345\tP12345\t1000\n")
    output = tmp_path / "pg.parquet"

    with MaxQuantPgAdapter() as adapter:
        with pytest.raises(ValueError, match="experiment_1.*mapping"):
            adapter.convert(str(protein_groups), str(output))

    assert not output.exists()


def test_pg_rejects_unassignable_total_intensity(tmp_path):
    """A total-only intensity must not create grouped_runs=['unknown']."""
    from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text("Protein IDs\tMajority protein IDs\tIntensity\nP12345\tP12345\t1000\n")
    output = tmp_path / "pg.parquet"

    with MaxQuantPgAdapter() as adapter:
        with pytest.raises(ValueError, match="total Intensity.*grouped_runs"):
            adapter.convert(str(protein_groups), str(output))

    assert not output.exists()


_EVIDENCE_HEADER = (
    "Sequence\tModified sequence\tCharge\tRaw file\tType\tMS/MS scan number\tm/z\t"
    "Calibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\t"
    "PEP\tLeading razor protein\tLeading proteins\tIntensity\tMass\tid\n"
)


def test_feature_dedup_collapses_identity_equal_rows(tmp_path):
    """Bug #248 (H): rows sharing the feature identity composite must collapse.

    The two ``run1`` rows share (peptidoform, charge, run, rt_start, rt_stop, scan)
    — i.e. the same ``feature_id`` — and differ only in ``Calibrated retention
    time``. The old dedup partitioned on calibrated RT (not part of the identity),
    so both survived as duplicate primary keys. They must now collapse to one,
    keeping the identified (MULTI-MSMS, lower PEP) row over the MBR one.
    """
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        _EVIDENCE_HEADER
        # Identified row and its MBR twin: same scan/rt-window, differ in calibrated RT.
        + "PEPTIDEK\t_PEPTIDEK_\t2\trun1\tMULTI-MSMS\t100\t450.25\t30.0\t29.5\t30.5\t0.01\tP1\tP1\t1000\t898.5\t1\n"
        + "PEPTIDEK\t_PEPTIDEK_\t2\trun1\tMULTI-MATCH\t100\t450.25\t30.4\t29.5\t30.5\t0.90\tP1\tP1\t1000\t898.5\t2\n"
        # A genuinely distinct feature (different charge/scan) must survive.
        + "PEPTIDEK\t_PEPTIDEK_\t3\trun1\tMULTI-MSMS\t200\t300.50\t40.0\t39.5\t40.5\t0.02\tP1\tP1\t2000\t898.5\t3\n"
    )
    output = tmp_path / "feature.parquet"

    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output), chunksize=1)

    rows = pq.read_table(output).to_pylist()
    assert len(rows) == 2, f"identity-equal rows must collapse, got {len(rows)}"
    # feature_id must be unique (no duplicate PK).
    assert len({row["feature_id"] for row in rows}) == 2
    kept = next(r for r in rows if r["charge"] == 2)
    assert kept["scan"] == [100]
    assert kept["posterior_error_probability"] == pytest.approx(0.01)
    assert kept["rt"] == pytest.approx(1800.0)
    assert kept["id_run_file_name"] == "run1"


def test_feature_dedup_does_not_degrade_without_type_column(tmp_path):
    """Bug #248 (H-b): a missing column must not degrade dedup to ``SELECT *``.

    With no ``Type`` column, the old guard fell back to ``SELECT * FROM evidence``,
    letting identity duplicates through. Dedup must still collapse identity-equal
    rows using the identity columns that are present.
    """
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    header = (
        "Sequence\tModified sequence\tCharge\tRaw file\tMS/MS scan number\tm/z\t"
        "Calibrated retention time\tCalibrated retention time start\t"
        "Calibrated retention time finish\tPEP\tLeading razor protein\tLeading proteins\tIntensity\tMass\tid\n"
    )
    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        header
        + "PEPTIDEK\t_PEPTIDEK_\t2\trun1\t100\t450.25\t30.0\t29.5\t30.5\t0.01\tP1\tP1\t1000\t898.5\t1\n"
        + "PEPTIDEK\t_PEPTIDEK_\t2\trun1\t100\t450.25\t30.4\t29.5\t30.5\t0.05\tP1\tP1\t1000\t898.5\t2\n"
    )
    output = tmp_path / "feature.parquet"

    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output), chunksize=100)

    rows = pq.read_table(output).to_pylist()
    assert len(rows) == 1, f"identity duplicates must collapse even without Type, got {len(rows)}"


def test_psm_handles_empty_modified_sequence(tmp_path):
    """Bug #248 (H): an empty ``Modified sequence`` must not abort the PSM view.

    The ``else None`` tuple-unpack raised ``TypeError`` for any row whose
    modified sequence was empty. The adapter must now skip that single row
    (peptidoform is a required schema field) instead of crashing the whole
    conversion, so the remaining PSMs still convert.
    """
    from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter

    msms = tmp_path / "msms.txt"
    msms.write_text(
        "Sequence\tModified sequence\tCharge\tRaw file\tScan number\tm/z\tMass\tPEP\tReverse\tScore\tProteins\n"
        "PEPTIDEK\t_PEPTIDEK_\t2\trun1\t100\t450.25\t898.5\t0.01\t\t100\tP1\n"
        # Empty modified sequence (renders to an empty peptidoform via to_proforma).
        "PEPTIDER\t_\t2\trun1\t101\t500.25\t998.5\t0.02\t\t90\tP2\n"
    )
    output = tmp_path / "psm.parquet"

    with MaxQuantPsmAdapter() as adapter:
        adapter.convert(str(msms), str(output))

    rows = pq.read_table(output).to_pylist()
    # The empty-modseq row is skipped; the valid PSM must still convert.
    assert len(rows) == 1, "valid PSM must survive; empty modseq row skipped, not crashing"
    assert rows[0]["sequence"] == "PEPTIDEK"
    assert rows[0]["peptidoform"]


def test_feature_pg_qvalue_and_genes_keyed_on_full_membership(tmp_path):
    """Bug #248 (M): a shared razor accession must not bind to another group's metadata.

    Both features share leading razor protein P1, but belong to different protein
    groups (P1;P2 vs P1;P3). The old leader-keyed map bound P1 to the first group
    by file order, so both features inherited that group's q-value/genes. The map
    must be keyed on the full protein-group membership.
    """
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        "Sequence\tModified sequence\tCharge\tRaw file\tType\tMS/MS scan number\tm/z\t"
        "Calibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\t"
        "PEP\tLeading razor protein\tLeading proteins\tGene names\tIntensity\tMass\tid\n"
        "PEPTIDEA\t_PEPTIDEA_\t2\trun1\tMULTI-MSMS\t100\t450.25\t30.0\t29.5\t30.5\t0.01\tP1\tP1;P2\t\t1000\t500\t1\n"
        "PEPTIDEB\t_PEPTIDEB_\t2\trun1\tMULTI-MSMS\t200\t460.25\t31.0\t30.5\t31.5\t0.02\tP1\tP1;P3\t\t2000\t520\t2\n"
    )
    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text(
        "Protein IDs\tMajority protein IDs\tQ-value\tGene names\nP1;P2\tP1;P2\t0.001\tGENEA\nP1;P3\tP1;P3\t0.5\tGENEB\n"
    )
    output = tmp_path / "feature.parquet"

    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output), protein_groups_path=str(protein_groups))

    rows = pq.read_table(output).to_pylist()
    by_seq = {r["sequence"]: r for r in rows}
    assert by_seq["PEPTIDEA"]["pg_global_qvalue"] == pytest.approx(0.001)
    assert by_seq["PEPTIDEA"]["gg_names"] == ["GENEA"]
    assert by_seq["PEPTIDEB"]["pg_global_qvalue"] == pytest.approx(0.5)
    assert by_seq["PEPTIDEB"]["gg_names"] == ["GENEB"]


def test_feature_warns_on_dropped_tmt_channels(caplog):
    """Bug #248 (M): channels with no reporter column must warn, not vanish silently."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    # evidence has only Reporter intensity 0-7; a 10-plex SDRF adds TMT130C (col 8)
    # and TMT131 (col 9), which have no column and would otherwise be dropped.
    columns = [f"Reporter intensity {i}" for i in range(8)]
    channels = [
        "TMT126",
        "TMT127N",
        "TMT127C",
        "TMT128N",
        "TMT128C",
        "TMT129N",
        "TMT129C",
        "TMT130N",
        "TMT130C",
        "TMT131",
    ]
    with MaxQuantFeatureAdapter() as adapter:
        with caplog.at_level("WARNING"):
            adapter._warn_dropped_tmt_channels(columns, "TMT10", channels, ri_offset=0)
            # Second call with the same dropped set must not re-warn.
            adapter._warn_dropped_tmt_channels(columns, "TMT10", channels, ri_offset=0)

    warnings = [r for r in caplog.records if r.levelname == "WARNING" and "dropped" in r.message.lower()]
    assert len(warnings) == 1, "dropped-channel warning must be emitted exactly once"
    assert "TMT130C" in caplog.text
    assert "TMT131" in caplog.text


def test_pg_warns_on_dropped_tmt_channels(caplog):
    """Bug #248 (M): the pg adapter must also warn about dropped channels."""
    from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

    with MaxQuantPgAdapter() as adapter:
        with caplog.at_level("WARNING"):
            adapter._warn_dropped_tmt_channels("exp1", ["TMT130C", "TMT131"])
            adapter._warn_dropped_tmt_channels("exp1", ["TMT130C", "TMT131"])

    warnings = [r for r in caplog.records if r.levelname == "WARNING" and "dropped" in r.message.lower()]
    assert len(warnings) == 1
    assert "TMT130C" in caplog.text and "TMT131" in caplog.text


# ---------------------------------------------------------------------------
# Module-scoped fixture: run conversion once, share outputs across tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def converted_output(tmp_path_factory):
    """Run MaxQuant conversion once for all tests in this module.

    Returns the output directory containing the generated Parquet files.
    If the test data is missing, all tests in this module are skipped.
    """
    from qpx.converters.maxquant.converter import MaxQuantConverter

    if not _MSMS.exists() or not _EVIDENCE.exists():
        pytest.skip(f"Test data not found: {EXAMPLES_DIR}")

    output = tmp_path_factory.mktemp("maxquant_output")

    converter = MaxQuantConverter()
    converter.convert(
        output_folder=output,
        msms_file=_MSMS,
        evidence_file=_EVIDENCE,
        sdrf_file=_SDRF if _SDRF.exists() else None,
        output_prefix=_PREFIX,
    )

    return output


@pytest.fixture(scope="module")
def psm_table(converted_output):
    """Read the psm.parquet produced by the converter."""
    path = converted_output / f"{_PREFIX}.psm.parquet"
    if not path.exists():
        pytest.skip("psm.parquet was not produced")
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def feature_table(converted_output):
    """Read the feature.parquet produced by the converter."""
    path = converted_output / f"{_PREFIX}.feature.parquet"
    if not path.exists():
        pytest.skip("feature.parquet was not produced")
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# PSM conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantPsmConversion:
    """Validate psm.parquet output from MaxQuant conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.psm.parquet"
        if not path.exists():
            raise AssertionError("psm.parquet was not created")

    def test_has_rows(self, psm_table):
        if psm_table.num_rows == 0:
            raise AssertionError("psm.parquet is empty")

    def test_key_columns_present(self, psm_table):
        column_names = set(psm_table.column_names)
        expected = {"sequence", "charge", "calculated_mz", "run_file_name", "scan"}
        missing = expected - column_names
        if missing:
            raise AssertionError(f"Missing columns: {missing}")

    def test_sequence_values_are_nonempty(self, psm_table):
        sequences = psm_table.column("sequence").to_pylist()
        for seq in sequences:
            if not (isinstance(seq, str) and len(seq) > 0):
                raise AssertionError(f"Expected non-empty string, got {seq!r}")

    def test_charge_values_are_valid(self, psm_table):
        charges = psm_table.column("charge").to_pylist()
        for charge in charges:
            if not 1 <= charge <= 10:
                raise AssertionError(f"Charge out of range: {charge}")

    def test_calculated_mz_values_are_nonnegative(self, psm_table):
        mz_values = psm_table.column("calculated_mz").to_pylist()
        for mz in mz_values:
            if mz is None or mz < 0:
                raise AssertionError(f"Invalid calculated_mz value: {mz}")

    def test_run_file_names_are_nonempty(self, psm_table):
        run_names = psm_table.column("run_file_name").to_pylist()
        for name in run_names:
            if not (isinstance(name, str) and len(name) > 0):
                raise AssertionError(f"Expected non-empty run name, got {name!r}")

    def test_schema_validation(self, psm_table):
        """Verify that the output conforms to the PsmSchema."""
        from qpx.core.data import PsmSchema

        errors = PsmSchema.validate(psm_table)
        if errors:
            raise AssertionError(f"Schema validation errors: {errors}")


# ---------------------------------------------------------------------------
# Feature conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantFeatureConversion:
    """Validate feature.parquet output from MaxQuant conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.feature.parquet"
        if not path.exists():
            raise AssertionError("feature.parquet was not created")

    def test_has_rows(self, feature_table):
        if feature_table.num_rows == 0:
            raise AssertionError("feature.parquet is empty")

    def test_key_columns_present(self, feature_table):
        column_names = set(feature_table.column_names)
        expected = {
            "sequence",
            "charge",
            "calculated_mz",
            "intensities",
            "run_file_name",
        }
        missing = expected - column_names
        if missing:
            raise AssertionError(f"Missing columns: {missing}")

    def test_sequence_values_are_nonempty(self, feature_table):
        sequences = feature_table.column("sequence").to_pylist()
        for seq in sequences:
            if not (isinstance(seq, str) and len(seq) > 0):
                raise AssertionError(f"Expected non-empty string, got {seq!r}")

    def test_charge_values_are_valid(self, feature_table):
        charges = feature_table.column("charge").to_pylist()
        for charge in charges:
            if not 1 <= charge <= 10:
                raise AssertionError(f"Charge out of range: {charge}")

    def test_intensities_are_nonnegative(self, feature_table):
        rows = feature_table.column("intensities").to_pylist()
        for row_intensities in rows:
            if row_intensities is None:
                continue
            for entry in row_intensities:
                if entry["intensity"] < 0:
                    raise AssertionError(f"Negative intensity: {entry['intensity']}")

    def test_run_file_names_are_nonempty(self, feature_table):
        run_names = feature_table.column("run_file_name").to_pylist()
        for name in run_names:
            if not (isinstance(name, str) and len(name) > 0):
                raise AssertionError(f"Expected non-empty run name, got {name!r}")

    def test_schema_validation(self, feature_table):
        """Verify that the output conforms to the FeatureSchema."""
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        if errors:
            raise AssertionError(f"Schema validation errors: {errors}")

    def test_tmt_nc_isomer_channels_have_distinct_intensities(self, feature_table):
        """Regression test: TMT127N and TMT127C must not systematically share the
        same intensity value.

        The alphabetical sort + i+1 fallback bug in _build_intensities caused
        TMT127C to read from the TMT127N MaxQuant column (and vice-versa),
        producing duplicate intensities for 96 %+ of features in real datasets.
        """
        rows = feature_table.column("intensities").to_pylist()

        equal_pairs = 0
        total_pairs = 0

        for row_intensities in rows:
            if not row_intensities:
                continue
            by_label = {e["label"]: e["intensity"] for e in row_intensities if e["intensity"] > 0}
            n_val = by_label.get("TMT127N")
            c_val = by_label.get("TMT127C")
            if n_val is not None and c_val is not None:
                total_pairs += 1
                if abs(n_val - c_val) < 1e-4:
                    equal_pairs += 1

        if total_pairs == 0:
            pytest.skip("No features with both TMT127N and TMT127C detected")

        duplicate_rate = equal_pairs / total_pairs
        if duplicate_rate >= 0.10:
            raise AssertionError(
                f"TMT127N == TMT127C in {equal_pairs}/{total_pairs} features "
                f"({100 * duplicate_rate:.1f}%) — N/C isomer channel assignment bug detected"
            )


# ---------------------------------------------------------------------------
# fixed_mod_only filter tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantFixedModOnly:
    """Validate that fixed_mod_only=True filters out variable-modification rows."""

    @pytest.fixture(scope="class")
    def fixed_mod_output(self, tmp_path_factory):
        from qpx.converters.maxquant.converter import MaxQuantConverter

        if not _EVIDENCE.exists():
            pytest.skip(f"Test data not found: {EXAMPLES_DIR}")

        output = tmp_path_factory.mktemp("maxquant_fixed_mod")
        converter = MaxQuantConverter()
        converter.convert(
            output_folder=output,
            evidence_file=_EVIDENCE,
            sdrf_file=_SDRF if _SDRF.exists() else None,
            output_prefix="fixed_mod_test",
            fixed_mod_only=True,
        )
        return output

    @pytest.fixture(scope="class")
    def fixed_feature_table(self, fixed_mod_output):
        path = fixed_mod_output / "fixed_mod_test.feature.parquet"
        if not path.exists():
            pytest.skip("fixed_mod feature.parquet was not produced")
        return pq.read_table(str(path))

    def test_fixed_mod_has_fewer_rows_than_unrestricted(self, feature_table, fixed_feature_table):
        """Filtering to fixed mods only must reduce the feature count."""
        if fixed_feature_table.num_rows >= feature_table.num_rows:
            raise AssertionError(
                f"fixed_mod_only produced {fixed_feature_table.num_rows} rows, "
                f"expected fewer than unrestricted {feature_table.num_rows}"
            )

    def test_fixed_mod_contains_no_variable_modifications(self, fixed_feature_table):
        """All features in fixed-mod output must have Unmodified or Carbamidomethyl (C) only.

        Uses startswith("carbamidomethyl") to match ProForma-parsed names regardless
        of whether the exact string is "Carbamidomethyl", "Carbamidomethyl (C)", etc.
        """
        mods_col = fixed_feature_table.column("modifications").to_pylist()
        for mods in mods_col:
            if mods is None:
                continue
            for mod in mods:
                name = (mod.get("modification_name") or "").strip().lower()
                if name:
                    if not name.startswith("carbamidomethyl"):
                        raise AssertionError(f"Variable modification found in fixed-mod output: {mod['modification_name']!r}")


# ---------------------------------------------------------------------------
# Ontology conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantOntologyConversion:
    """Validate ontology.parquet output from MaxQuant conversion.

    The ontology file is only written when scores are discovered during
    conversion. It may legitimately not exist if no scores were tracked.
    """

    def test_ontology_file_exists_or_no_scores(self, converted_output):
        path = converted_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")

        table = pq.read_table(str(path))
        if table.num_rows == 0:
            raise AssertionError("ontology.parquet exists but is empty")

    def test_ontology_columns(self, converted_output):
        path = converted_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")

        table = pq.read_table(str(path))
        column_names = set(table.column_names)
        expected = {"field_name", "view"}
        missing = expected - column_names
        if missing:
            raise AssertionError(f"Missing ontology columns: {missing}")

    def test_schema_validation(self, converted_output):
        from qpx.core.data import OntologySchema

        path = converted_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")
        table = pq.read_table(str(path))
        errors = OntologySchema.validate(table)
        if errors:
            raise AssertionError(f"Schema validation errors: {errors}")
