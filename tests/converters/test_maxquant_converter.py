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
    labels = ["TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131"]
    channels = [(label, index) for index, label in enumerate(labels)]
    with MaxQuantFeatureAdapter() as adapter:
        with caplog.at_level("WARNING"):
            adapter._warn_dropped_tmt_channels(columns, channels, ri_offset=0)
            # Second call with the same dropped set must not re-warn.
            adapter._warn_dropped_tmt_channels(columns, channels, ri_offset=0)

    warnings = [r for r in caplog.records if r.levelname == "WARNING" and "dropped" in r.message.lower()]
    assert len(warnings) == 1, "dropped-channel warning must be emitted exactly once"
    assert "TMT130C" in caplog.text
    assert "TMT131" in caplog.text


def _reporter_features(tmp_path, labels, n_columns):
    """Convert one evidence row whose ``Reporter intensity i`` is ``100 * (i + 1)``."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter

    reporter = "\t".join(f"Reporter intensity {i}" for i in range(n_columns))
    values = "\t".join(str(100 * (i + 1)) for i in range(n_columns))
    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        _EVIDENCE_HEADER.rstrip("\n") + f"\t{reporter}\n"
        f"PEPTIDEK\t_PEPTIDEK_\t2\trun1\tMULTI-MSMS\t10\t500.0\t10.0\t9.9\t10.1\t0.01\tP1\tP1\t1000\t998.0\t0\t{values}\n"
    )
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[data file]\tcomment[label]\n"
        + "".join(f"S{i}\trun1.raw\t{label}\n" for i, label in enumerate(labels))
    )
    output = tmp_path / "feature.parquet"
    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(evidence), str(output), sdrf_path=str(sdrf))
    (row,) = pq.read_table(output).to_pylist()
    return {entry["label"]: entry["intensity"] for entry in row["intensities"]}


@pytest.mark.parametrize(
    ("labels", "n_columns", "expected"),
    [
        # TMT6: TMT131 is the sixth reporter, not TMT10's tenth.
        (["TMT131", "TMT126", "TMT129", "TMT127", "TMT130", "TMT128"], 6, {"TMT126": 100, "TMT127": 200, "TMT131": 600}),
        # iTRAQ channels follow the plex position whatever the SDRF row order.
        (["itraq117", "ITRAQ114", "iTRAQ116", "ITRAQ115"], 4, {"ITRAQ114": 100, "ITRAQ115": 200, "ITRAQ117": 400}),
        (["ITRAQ121", "ITRAQ113", "ITRAQ119"], 8, {"ITRAQ113": 100, "ITRAQ119": 700, "ITRAQ121": 800}),
        # A partial TMTpro-18 plex resolves to its 18-channel positions.
        (["TMT126", "TMT134C", "TMT135N"], 18, {"TMT126": 100, "TMT134C": 1700, "TMT135N": 1800}),
    ],
)
def test_feature_reads_reporter_columns_by_plex_position(tmp_path, labels, n_columns, expected):
    """Each SDRF label reads the MaxQuant reporter column of its plex position."""
    intensities = _reporter_features(tmp_path, labels, n_columns)
    assert len(intensities) == len(labels)
    for label, value in expected.items():
        assert intensities[label] == value


def test_feature_rejects_labels_outside_any_plex(tmp_path):
    with pytest.raises(ValueError, match="match no TMT plex"):
        _reporter_features(tmp_path, ["TMT126", "TMT10plex-127"], 10)


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
    assert path.exists(), "psm.parquet was not produced"
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def feature_table(converted_output):
    """Read the feature.parquet produced by the converter."""
    path = converted_output / f"{_PREFIX}.feature.parquet"
    assert path.exists(), "feature.parquet was not produced"
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
                if entry["intensity"] is not None and entry["intensity"] < 0:
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
            by_label = {e["label"]: e["intensity"] for e in row_intensities if e["intensity"]}
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
        assert path.exists(), "ontology.parquet not written (no scores discovered)"

        table = pq.read_table(str(path))
        if table.num_rows == 0:
            raise AssertionError("ontology.parquet exists but is empty")

    def test_ontology_columns(self, converted_output):
        path = converted_output / f"{_PREFIX}.ontology.parquet"
        assert path.exists(), "ontology.parquet not written (no scores discovered)"

        table = pq.read_table(str(path))
        column_names = set(table.column_names)
        expected = {"field_name", "view"}
        missing = expected - column_names
        if missing:
            raise AssertionError(f"Missing ontology columns: {missing}")

    def test_schema_validation(self, converted_output):
        from qpx.core.data import OntologySchema

        path = converted_output / f"{_PREFIX}.ontology.parquet"
        assert path.exists(), "ontology.parquet not written (no scores discovered)"
        table = pq.read_table(str(path))
        errors = OntologySchema.validate(table)
        if errors:
            raise AssertionError(f"Schema validation errors: {errors}")


# ---------------------------------------------------------------------------
# feature<->pg softlink round-trip (bigbio/qpx#269)
# ---------------------------------------------------------------------------

from tests.converters.pg_ids_roundtrip import assert_softlink_valid, open_converted, pg_ids_by_sequence  # noqa: E402


def _write_pg_link_inputs(tmp_path):
    """Minimal LFQ MaxQuant inputs exercising shared-leader + a no-pg-row case.

    PEPTIDEA -> group P1;P2, PEPTIDEB -> P1;P3 (distinct groups, shared leader P1),
    PEPTIDEC -> P9 which is absent from proteinGroups (so its feature has no pg row).
    All three features share Experiment exp1 / Raw file run1.
    """
    evidence = tmp_path / "evidence.txt"
    header = (
        "Sequence\tModified sequence\tCharge\tRaw file\tExperiment\tType\tMS/MS scan number\tm/z\t"
        "Calibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\t"
        "PEP\tLeading razor protein\tLeading proteins\tGene names\tIntensity\tMass\tid\n"
    )
    rows = (
        "PEPTIDEA\t_PEPTIDEA_\t2\trun1\texp1\tMULTI-MSMS\t100\t450.25\t30.0\t29.5\t30.5\t0.01\tP1\tP1;P2\t\t1000\t500\t1\n"
        "PEPTIDEB\t_PEPTIDEB_\t2\trun1\texp1\tMULTI-MSMS\t200\t460.25\t31.0\t30.5\t31.5\t0.02\tP1\tP1;P3\t\t2000\t520\t2\n"
        "PEPTIDEC\t_PEPTIDEC_\t2\trun1\texp1\tMULTI-MSMS\t300\t470.25\t32.0\t31.5\t32.5\t0.03\tP9\tP9\t\t3000\t540\t3\n"
    )
    evidence.write_text(header + rows)

    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text(
        "Protein IDs\tMajority protein IDs\tQ-value\tGene names\tIntensity exp1\n"
        "P1;P2\tP1;P2\t0.001\tGENEA\t1000\n"
        "P1;P3\tP1;P3\t0.002\tGENEB\t2000\n"
    )
    return evidence, protein_groups


def test_feature_pg_softlink_shared_leader_and_no_pg_row(tmp_path):
    """Full pg+feature conversion: the computed softlink links each feature to a
    real pg row with matching membership + run + carried label; shared-leader
    groups stay distinct; a feature whose group has no pg row produces no link."""
    from qpx.converters.maxquant.converter import MaxQuantConverter
    from qpx.core.constants import FEATURE, PG

    evidence, protein_groups = _write_pg_link_inputs(tmp_path)
    out = tmp_path / "out"
    MaxQuantConverter().convert(
        output_folder=str(out),
        evidence_file=str(evidence),
        protein_groups_file=str(protein_groups),
        output_prefix="mq",
        structures=[FEATURE, PG],
    )

    with open_converted(out, prefix="mq") as ds:
        feat, pg, link = assert_softlink_valid(ds)

    pg_id_by_membership = {tuple(sorted(memb)): pg_id for pg_id, memb, *_ in pg}
    ids_by_seq = pg_ids_by_sequence(feat, link)
    assert ids_by_seq["PEPTIDEA"] == [pg_id_by_membership[("P1", "P2")]]
    assert ids_by_seq["PEPTIDEB"] == [pg_id_by_membership[("P1", "P3")]]
    assert set(ids_by_seq["PEPTIDEA"]).isdisjoint(ids_by_seq["PEPTIDEB"])
    assert ids_by_seq["PEPTIDEC"] == []  # group P9 absent from pg -> no link


@pytest.mark.integration
def test_feature_pg_softlink_real_full_dataset(tmp_path_factory):
    """Real MaxQuant example (maxquant_full): the computed softlink links features
    to real pg rows with matching membership + run + carried label (bigbio/qpx#269)."""
    from qpx.converters.maxquant.converter import MaxQuantConverter
    from qpx.core.constants import FEATURE, PG

    full = Path(__file__).resolve().parent.parent / "examples" / "maxquant" / "maxquant_full"
    evidence = full / "evidence.txt.gz"
    protein_groups = full / "proteinGroups.txt"
    if not evidence.exists() or not protein_groups.exists():
        pytest.skip("maxquant_full fixture not found")

    out = tmp_path_factory.mktemp("mq_pg_ids")
    MaxQuantConverter().convert(
        output_folder=str(out),
        evidence_file=str(evidence),
        protein_groups_file=str(protein_groups),
        sdrf_file=str(full / "PXD001819.sdrf.tsv"),
        output_prefix="mq",
        structures=[FEATURE, PG],
    )

    with open_converted(out, prefix="mq") as ds:
        _feat, _pg, link = assert_softlink_valid(ds)
    assert link, "expected at least some computed feature->pg softlink edges"


def _write_group_inputs(tmp_path):
    """Evidence with ``Protein group IDs`` and a matching proteinGroups.txt.

    PEPTIDEA is in group 0 (P1;P2;P3) though its ``Leading proteins`` is only P1;
    PEPTIDEB is in group 1 (P4); PEPTIDEC is shared by both groups and assigned to
    the group of its razor protein P4. Row id 3 repeats PEPTIDEA's identity.
    """
    evidence = tmp_path / "evidence.txt"
    evidence.write_text(
        "Sequence\tModified sequence\tCharge\tRaw file\tExperiment\tType\tMS/MS scan number\tm/z\t"
        "Calibrated retention time\tCalibrated retention time start\tCalibrated retention time finish\t"
        "PEP\tProteins\tLeading razor protein\tLeading proteins\tProtein group IDs\tIntensity\tMass\tid\n"
        "PEPTIDEA\t_PEPTIDEA_\t2\trun1\texp1\tMULTI-MSMS\t100\t450.25\t30.0\t29.5\t30.5\t0.01\tP1;P2;P3\tP1\tP1\t0\t1000\t500\t0\n"
        "PEPTIDEB\t_PEPTIDEB_\t2\trun1\texp1\tMULTI-MSMS\t200\t460.25\t31.0\t30.5\t31.5\t0.02\tP4\tP4\tP4\t1\t\t520\t1\n"
        "PEPTIDEC\t_PEPTIDEC_\t2\trun1\texp1\tMULTI-MSMS\t300\t470.25\t32.0\t31.5\t32.5\t0.03\tP1;P4\tP4\tP4;P1\t0;1\t3000\t540\t2\n"
        "PEPTIDEA\t_PEPTIDEA_\t2\trun1\texp1\tMULTI-MSMS\t100\t450.25\t30.0\t29.5\t30.5\t0.05\tP1;P2;P3\tP1\tP1\t0\t900\t500\t3\n"
    )
    protein_groups = tmp_path / "proteinGroups.txt"
    protein_groups.write_text(
        "Protein IDs\tMajority protein IDs\tPeptide counts (all)\tPeptides\tUnique peptides\tQ-value\tGene names\t"
        "Intensity exp1\tPeptides exp1\tUnique peptides exp1\tMS/MS count exp1\tSequence coverage exp1 [%]\tid\n"
        "P1;P2;P3\tP1\t2;1;1\t2\t1\t0.001\tGENEA\t4000\t1\t1\t3\t12.5\t0\n"
        "P4\tP4\t2\t2\t2\t0.002\tGENEB\t3000\t2\t2\t2\t20\t1\n"
    )
    msms = tmp_path / "msms.txt"
    msms.write_text(
        "Sequence\tModified sequence\tCharge\tRaw file\tScan number\tm/z\tRetention time\tPEP\tScore\tProteins\tMass\tEvidence ID\n"
        "PEPTIDEA\t_PEPTIDEA_\t2\trun1\t100\t450.25\t30.0\t0.01\t100\tP1\t500\t0\n"
        "PEPTIDEA\t_PEPTIDEA_\t2\trun1\t101\t450.25\t30.1\t0.05\t80\tP1\t500\t3\n"
        "PEPTIDEC\t_PEPTIDEC_\t2\trun1\t300\t470.25\t32.0\t0.03\t90\tP4\t540\t2\n"
        "PEPTIDED\t_PEPTIDED_\t2\trun1\t400\t480.25\t33.0\t0.03\t90\tP4\t560\t\n"
    )
    return evidence, protein_groups, msms


def test_feature_protein_group_follows_protein_group_ids(tmp_path):
    """The feature's group is MaxQuant's own assignment, so every feature links to its pg row."""
    from qpx.converters.maxquant.converter import MaxQuantConverter
    from qpx.core.constants import FEATURE, PG, PSM

    evidence, protein_groups, msms = _write_group_inputs(tmp_path)
    out = tmp_path / "out"
    MaxQuantConverter().convert(
        output_folder=str(out),
        msms_file=str(msms),
        evidence_file=str(evidence),
        protein_groups_file=str(protein_groups),
        output_prefix="mq",
        structures=[PSM, FEATURE, PG],
    )

    by_seq = {row["sequence"]: row for row in pq.read_table(out / "mq.feature.parquet").to_pylist()}
    assert [p["accession"] for p in by_seq["PEPTIDEA"]["pg_accessions"]] == ["P1", "P2", "P3"]
    assert by_seq["PEPTIDEA"]["unique"] is False
    assert by_seq["PEPTIDEB"]["unique"] is True
    # A shared peptide goes to its razor protein's group and is not unique to it.
    assert [p["accession"] for p in by_seq["PEPTIDEC"]["pg_accessions"]] == ["P4"]
    assert (by_seq["PEPTIDEC"]["anchor_protein"], by_seq["PEPTIDEC"]["unique"]) == ("P4", False)
    assert by_seq["PEPTIDEC"]["pg_global_qvalue"] == pytest.approx(0.002)
    # An evidence row MaxQuant did not quantify keeps its label with a null intensity.
    assert by_seq["PEPTIDEB"]["intensities"] == [{"label": "LFQ", "intensity": None}]

    with open_converted(out, prefix="mq") as ds:
        assert ds.features_without_pg_link().to_df().empty

    # Each PSM links to its feature through "Evidence ID", collapsed duplicates included.
    feature_ids = {row["sequence"]: row["feature_id"] for row in pq.read_table(out / "mq.feature.parquet").to_pylist()}
    psms = pq.read_table(out / "mq.psm.parquet").to_pylist()
    assert [(row["sequence"], row["feature_id"]) for row in psms] == [
        ("PEPTIDEA", feature_ids["PEPTIDEA"]),
        ("PEPTIDEA", feature_ids["PEPTIDEA"]),
        ("PEPTIDEC", feature_ids["PEPTIDEC"]),
        ("PEPTIDED", None),
    ]


def test_pg_counts_are_per_experiment(tmp_path):
    """Peptide counts, coverage and the MS/MS count come from the experiment's own columns."""
    from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

    evidence, protein_groups, _ = _write_group_inputs(tmp_path)
    output = tmp_path / "pg.parquet"
    with MaxQuantPgAdapter() as adapter:
        adapter.convert(str(protein_groups), str(output), evidence_path=str(evidence))

    row = next(r for r in pq.read_table(output).to_pylist() if r["anchor_protein"] == "P1")
    assert row["peptide_counts"] == {"unique_sequences": 1, "total_sequences": 1}
    assert row["feature_counts"] is None
    assert row["sequence_coverage"] == pytest.approx(12.5)
    assert row["peptides"] == [
        {"protein_name": "P1", "peptide_count": 2},
        {"protein_name": "P2", "peptide_count": 1},
        {"protein_name": "P3", "peptide_count": 1},
    ]
    assert row["additional_intensities"] == [
        {"label": "LFQ", "intensities": [{"intensity_name": "spectral_count", "intensity_value": 3.0}]}
    ]


def test_every_mod_probabilities_column_reaches_psm_and_feature(tmp_path):
    """Oxidation (M) Probabilities in msms.txt and evidence.txt, not only Phospho (STY) in msms.txt."""
    from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter
    from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter

    psm_path = tmp_path / "t.psm.parquet"
    feature_path = tmp_path / "t.feature.parquet"
    with MaxQuantPsmAdapter() as adapter:
        adapter.convert(str(EXAMPLES_DIR / "msms.txt"), str(psm_path))
    with MaxQuantFeatureAdapter() as adapter:
        adapter.convert(str(EXAMPLES_DIR / "evidence.txt"), str(feature_path))

    for path in (psm_path, feature_path):
        rows = pq.read_table(path, columns=["peptidoform", "modifications", "modification_candidates"]).to_pylist()
        localized = [row for row in rows if row["modification_candidates"]]
        assert localized, path.name
        for row in localized:
            assigned = [c for c in row["modification_candidates"] if c["assigned"]]
            assert all(c["name"] == "Oxidation" and c["accession"] == "UNIMOD:35" for c in row["modification_candidates"])
            assert assigned and "M[UNIMOD:35]" in row["peptidoform"]
        placed = next(
            pos
            for row in localized
            for mod in row["modifications"]
            if mod["accession"] == "UNIMOD:35"
            for pos in mod["positions"]
        )
        assert placed["scores"][0]["score_name"] == "localization_probability"
