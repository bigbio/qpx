"""DIA-NN converter integration tests using small example dataset.

Runs the full DIA-NN conversion pipeline once (module scope) and validates
that feature.parquet, pg.parquet, and ontology.parquet are written with
correct schemas and plausible values.
"""

import re
from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from tests.conftest import assert_pg_table_wellformed

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples" / "diann" / "small"

_REPORT = EXAMPLES_DIR / "diann_report.tsv"
_SDRF = EXAMPLES_DIR / "PXD019909-DIA.sdrf.tsv"
_PG_MATRIX = EXAMPLES_DIR / "diann_report.pg_matrix.tsv"
_MZML_TSV_DIR = EXAMPLES_DIR / "mzml"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _parse_scan_number(spectrum_id: str) -> int:
    """Extract numeric scan number from a SpectrumID string like
    'controllerType=0 controllerNumber=1 scan=42'.
    """
    m = re.search(r"scan=(\d+)", str(spectrum_id))
    return int(m.group(1)) if m else 0


def _prepare_ms_info_parquet(source_dir: Path, dest_dir: Path) -> None:
    """Convert *_mzml_info.tsv files into *_ms_info.parquet files.

    The feature adapter expects ``*_ms_info.parquet`` files with columns
    ``rt``, ``scan``, and ``precursor_mz``.  The test data ships as TSV
    with different column names, so this helper bridges the gap.
    """
    for tsv_path in source_dir.glob("*_mzml_info.tsv"):
        df = pd.read_csv(tsv_path, sep="\t")
        out = pd.DataFrame(
            {
                "rt": df["Retention_Time"].astype(float),
                "scan": df["SpectrumID"].apply(_parse_scan_number).astype(int),
                "precursor_mz": pd.to_numeric(df["Exp_Mass_To_Charge"], errors="coerce"),
            }
        )
        stem = tsv_path.stem.replace("_mzml_info", "")
        parquet_name = f"{stem}_ms_info.parquet"
        pq.write_table(pa.Table.from_pandas(out), str(dest_dir / parquet_name))


def _write_plexdia_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Write a two-channel DIA-NN report, PG matrix, and matching SDRF."""
    shared = {
        "Run": "run_A.RAW",
        "Protein.Group": "P1",
        "Protein.Names": "Protein 1",
        "Genes": "G1",
        "Modified.Sequence": "PEPTIDEK",
        "Stripped.Sequence": "PEPTIDEK",
        "Precursor.Id": "PEPTIDEK2",
        "Precursor.Charge": 2,
        "PEP": 0.001,
        "Global.Q.Value": 0.001,
        "PG.Q.Value": 0.002,
        "Global.PG.Q.Value": 0.003,
        "GG.Q.Value": 0.004,
        "Proteotypic": 1,
        "RT": 10.0,
    }
    rows = [
        {
            **shared,
            "Channel": "L",
            "Q.Value": 0.002,
            "Precursor.Quantity": 100.0,
            "PG.Quantity": 1000.0,
            "PG.MaxLFQ": 900.0,
        },
        {
            **shared,
            "Channel": "H",
            "Q.Value": 0.003,
            "Precursor.Quantity": 200.0,
            "PG.Quantity": 2000.0,
            "PG.MaxLFQ": 1800.0,
        },
    ]
    report_path = tmp_path / "plexdia_report.tsv"
    pd.DataFrame(rows).to_csv(report_path, sep="\t", index=False)

    matrix_path = tmp_path / "plexdia_pg_matrix.tsv"
    pd.DataFrame(
        [
            {
                "Protein.Group": "P1",
                "Protein.Names": "Protein 1",
                "Genes": "G1",
                "run_A": 900.0,
            }
        ]
    ).to_csv(matrix_path, sep="\t", index=False)

    sdrf_path = tmp_path / "plexdia.sdrf.tsv"
    pd.DataFrame(
        [
            {
                "source name": "SAMPLE_L",
                "comment[data file]": "run_A.RAW",
                "comment[label]": "L",
            },
            {
                "source name": "SAMPLE_H",
                "comment[data file]": "run_A.RAW",
                "comment[label]": "H",
            },
        ]
    ).to_csv(sdrf_path, sep="\t", index=False)
    return report_path, matrix_path, sdrf_path


# ---------------------------------------------------------------------------
# Module-scoped fixture: run conversion once, share outputs across tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def converted_output(tmp_path_factory):
    """Run DIA-NN conversion once for all tests in this module.

    Returns the output directory containing the generated Parquet files.
    If the test data is missing, all tests in this module are skipped.
    """
    from qpx.converters.diann.converter import DiaNNConverter

    if not _REPORT.exists():
        pytest.skip(f"Test data not found: {_REPORT}")

    output = tmp_path_factory.mktemp("diann_output")

    # Prepare ms_info parquet files from the shipped TSVs.
    mzml_parquet_dir = tmp_path_factory.mktemp("mzml_info")
    _prepare_ms_info_parquet(_MZML_TSV_DIR, mzml_parquet_dir)

    converter = DiaNNConverter(
        report_path=str(_REPORT),
        sdrf_path=str(_SDRF),
    )
    converter.convert_features(
        mzml_info_folder=str(mzml_parquet_dir),
        qvalue_threshold=0.05,
        output_folder=str(output),
        output_prefix="diann_test",
    )
    converter.convert_pg(
        pg_matrix_path=str(_PG_MATRIX),
        output_folder=str(output),
        output_prefix="diann_test",
    )
    converter.write_ontology(
        output_folder=output,
        prefix="diann_test",
    )

    return output


@pytest.fixture(scope="module")
def feature_table(converted_output):
    """Read the feature.parquet produced by the converter."""
    path = converted_output / "diann_test.feature.parquet"
    if not path.exists():
        pytest.skip("feature.parquet was not produced")
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def pg_table(converted_output):
    """Read the pg.parquet produced by the converter."""
    path = converted_output / "diann_test.pg.parquet"
    if not path.exists():
        pytest.skip("pg.parquet was not produced")
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# Feature conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestDiaNNFeatureConversion:
    """Validate feature.parquet output from DIA-NN conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / "diann_test.feature.parquet"
        assert path.exists(), "feature.parquet was not created"

    def test_has_rows(self, feature_table):
        assert feature_table.num_rows > 0, "feature.parquet is empty"

    def test_sequence_values_are_nonempty(self, feature_table):
        sequences = feature_table.column("sequence").to_pylist()
        for seq in sequences:
            assert isinstance(seq, str) and len(seq) > 0, f"Expected non-empty string, got {seq!r}"

    def test_charge_values_are_valid(self, feature_table):
        charges = feature_table.column("charge").to_pylist()
        for charge in charges:
            assert 1 <= charge <= 10, f"Charge out of range: {charge}"

    def test_intensities_are_nonnegative(self, feature_table):
        rows = feature_table.column("intensities").to_pylist()
        for row_intensities in rows:
            assert row_intensities is not None and len(row_intensities) > 0
            for entry in row_intensities:
                assert entry["intensity"] >= 0, f"Negative intensity: {entry['intensity']}"

    def test_run_file_names_are_nonempty(self, feature_table):
        run_names = feature_table.column("run_file_name").to_pylist()
        for name in run_names:
            assert isinstance(name, str) and len(name) > 0

    def test_schema_validation(self, feature_table):
        """Verify that the output conforms to the FeatureSchema."""
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        assert not errors, f"Schema validation errors: {errors}"


def test_plexdia_feature_collapses_channels_into_one_row(tmp_path):
    """Channel rows for one precursor become one feature with two intensities."""
    from qpx.converters.diann.feature_adapter import DiannFeatureAdapter

    report_path, _, sdrf_path = _write_plexdia_inputs(tmp_path)
    output_path = tmp_path / "plexdia.feature.parquet"
    with DiannFeatureAdapter() as adapter:
        adapter.convert(
            diann_report=str(report_path),
            output_path=str(output_path),
            sdrf_path=str(sdrf_path),
            qvalue_threshold=0.01,
        )

    rows = pq.read_table(output_path).to_pylist()
    assert len(rows) == 1
    assert rows[0]["run_file_name"] == "run_A"
    assert {entry["label"]: entry["intensity"] for entry in rows[0]["intensities"]} == {"L": 100.0, "H": 200.0}


# ---------------------------------------------------------------------------
# Protein group conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestDiaNNPgConversion:
    """Validate pg.parquet output from DIA-NN conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / "diann_test.pg.parquet"
        assert path.exists(), "pg.parquet was not created"

    def test_has_rows(self, pg_table):
        assert pg_table.num_rows > 0, "pg.parquet is empty"

    def test_pg_table_wellformed(self, pg_table):
        """grouped_runs, pg_accessions, anchor_protein, intensities are all well-formed."""
        assert_pg_table_wellformed(pg_table)

    def test_schema_validation(self, pg_table):
        """Verify that the output conforms to the PgSchema."""
        from qpx.core.data import PgSchema

        errors = PgSchema.validate(pg_table)
        assert not errors, f"Schema validation errors: {errors}"


def test_plexdia_pg_preserves_each_channel_quantity(tmp_path):
    """Protein-group aggregation must not discard all but the first channel."""
    from qpx.converters.diann.pg_adapter import DiannPgAdapter

    report_path, matrix_path, sdrf_path = _write_plexdia_inputs(tmp_path)
    output_path = tmp_path / "plexdia.pg.parquet"
    with DiannPgAdapter() as adapter:
        adapter.convert(
            diann_report=str(report_path),
            pg_matrix_path=str(matrix_path),
            output_path=str(output_path),
            sdrf_path=str(sdrf_path),
        )

    # pg is flattened since 1.1: each channel is its own row (one row per label).
    rows = pq.read_table(output_path).to_pylist()
    assert len(rows) == 2
    assert all(row["grouped_runs"] == ["run_A"] for row in rows)
    assert {row["label"]: row["intensity"] for row in rows} == {"L": 1000.0, "H": 2000.0}
    assert {row["label"]: row["additional_intensities"][0]["intensities"][0]["intensity_value"] for row in rows} == {
        "L": 900.0,
        "H": 1800.0,
    }


def test_diann_pg_inconsistent_annotations_do_not_split_group(tmp_path):
    """One protein group with inconsistent per-precursor name/gene annotations
    must produce a single unique pg record, not duplicate-identity rows.

    Regression for bigbio/qpx#240: grouping on pg_names/gg_accessions split a
    group whose precursor rows carried different Genes/Protein.Names into records
    that shared the same (pg_accessions, grouped_runs, label) identity, colliding
    on pg_id and tripping the close-time uniqueness validator.
    """
    from qpx.converters.diann.pg_adapter import DiannPgAdapter

    shared = {
        "Run": "run_A",
        "Protein.Group": "P1;P2",
        "Protein.Names": "N1;N2",
        "PG.Quantity": 1000.0,
        "PG.MaxLFQ": 900.0,
        "Precursor.Charge": 2,
        "Q.Value": 0.002,
        "PG.Q.Value": 0.002,
        "Global.PG.Q.Value": 0.003,
        "GG.Q.Value": 0.004,
        "Proteotypic": 1,
        "Precursor.Quantity": 100.0,
    }
    report_path = tmp_path / "inconsistent_report.tsv"
    pd.DataFrame(
        [
            {**shared, "Genes": "G1;G2", "Stripped.Sequence": "PEPTIDEK", "Precursor.Id": "PEPTIDEK2"},
            {**shared, "Genes": "G1", "Stripped.Sequence": "PEPTIDER", "Precursor.Id": "PEPTIDER2"},
        ]
    ).to_csv(report_path, sep="\t", index=False)
    matrix_path = tmp_path / "inconsistent_pg_matrix.tsv"
    pd.DataFrame([{"Protein.Group": "P1;P2", "Protein.Names": "N1;N2", "Genes": "G1;G2", "run_A": 900.0}]).to_csv(
        matrix_path, sep="\t", index=False
    )

    output_path = tmp_path / "inconsistent.pg.parquet"
    # Without the fix, adapter.convert would raise at close time:
    # "Primary key (pg_id) has 1 duplicate row(s)".
    with DiannPgAdapter() as adapter:
        adapter.convert(
            diann_report=str(report_path),
            pg_matrix_path=str(matrix_path),
            output_path=str(output_path),
        )

    table = pq.read_table(output_path)
    pg_ids = table.column("pg_id").to_pylist()
    assert table.num_rows == 1, "the single protein group must yield exactly one record"
    assert len(set(pg_ids)) == len(pg_ids), "pg_id must be unique"
    row = table.to_pylist()[0]
    assert row["pg_accessions"] == ["P1", "P2"]
    assert row["anchor_protein"] == "P1"


def test_diann_pg_distinct_groups_sharing_leader_get_distinct_ids(tmp_path):
    """Two DIFFERENT protein groups that share a leading protein must get
    DISTINCT pg_ids.

    The pg identity keys on the full pg_accessions membership, not on the leader
    alone, so ``P1;P2`` and ``P1;P3`` (same anchor ``P1``, same run, same label)
    are distinct groups and must not collide on pg_id.
    """
    from qpx.converters.diann.pg_adapter import DiannPgAdapter

    common = {
        "Run": "run_A",
        "PG.Quantity": 1000.0,
        "PG.MaxLFQ": 900.0,
        "Precursor.Charge": 2,
        "Q.Value": 0.002,
        "PG.Q.Value": 0.002,
        "Global.PG.Q.Value": 0.003,
        "GG.Q.Value": 0.004,
        "Proteotypic": 1,
        "Precursor.Quantity": 100.0,
    }
    report_path = tmp_path / "shared_leader_report.tsv"
    pd.DataFrame(
        [
            {
                **common,
                "Protein.Group": "P1;P2",
                "Protein.Names": "N1;N2",
                "Genes": "G1;G2",
                "Stripped.Sequence": "PEPTIDEK",
                "Precursor.Id": "PEPTIDEK2",
            },
            {
                **common,
                "Protein.Group": "P1;P3",
                "Protein.Names": "N1;N3",
                "Genes": "G1;G3",
                "Stripped.Sequence": "PEPTIDER",
                "Precursor.Id": "PEPTIDER2",
            },
        ]
    ).to_csv(report_path, sep="\t", index=False)
    matrix_path = tmp_path / "shared_leader_pg_matrix.tsv"
    pd.DataFrame(
        [
            {"Protein.Group": "P1;P2", "Protein.Names": "N1;N2", "Genes": "G1;G2", "run_A": 900.0},
            {"Protein.Group": "P1;P3", "Protein.Names": "N1;N3", "Genes": "G1;G3", "run_A": 800.0},
        ]
    ).to_csv(matrix_path, sep="\t", index=False)

    output_path = tmp_path / "shared_leader.pg.parquet"
    with DiannPgAdapter() as adapter:
        adapter.convert(
            diann_report=str(report_path),
            pg_matrix_path=str(matrix_path),
            output_path=str(output_path),
        )

    table = pq.read_table(output_path)
    rows = table.to_pylist()
    assert table.num_rows == 2, "two distinct protein groups must yield two records"
    assert all(row["anchor_protein"] == "P1" for row in rows)
    pg_ids = table.column("pg_id").to_pylist()
    assert len(set(pg_ids)) == 2, "distinct membership must derive distinct pg_ids"


def test_diann_maxlfq_fallback_keeps_experimental_label(tmp_path):
    """MaxLFQ-only PG output remains joinable to an LFQ run sample."""
    from qpx.converters.diann.pg_adapter import DiannPgAdapter

    report_path = tmp_path / "maxlfq_report.tsv"
    pd.DataFrame(
        [
            {
                "Run": "run_A",
                "Protein.Group": "P1",
                "Protein.Names": "Protein 1",
                "Genes": "G1",
                "PG.MaxLFQ": 700.0,
                "Modified.Sequence": "PEPTIDEK",
                "Stripped.Sequence": "PEPTIDEK",
                "Precursor.Id": "PEPTIDEK2",
                "Precursor.Charge": 2,
                "Q.Value": 0.002,
                "PG.Q.Value": 0.002,
                "Global.PG.Q.Value": 0.003,
                "GG.Q.Value": 0.004,
                "Proteotypic": 1,
                "Precursor.Quantity": 100.0,
            }
        ]
    ).to_csv(report_path, sep="\t", index=False)
    matrix_path = tmp_path / "maxlfq_pg_matrix.tsv"
    pd.DataFrame(
        [
            {
                "Protein.Group": "P1",
                "Protein.Names": "Protein 1",
                "Genes": "G1",
                "run_A": 700.0,
            }
        ]
    ).to_csv(matrix_path, sep="\t", index=False)

    output_path = tmp_path / "maxlfq.pg.parquet"
    with DiannPgAdapter() as adapter:
        adapter.convert(
            diann_report=str(report_path),
            pg_matrix_path=str(matrix_path),
            output_path=str(output_path),
        )

    row = pq.read_table(output_path).to_pylist()[0]
    assert row["label"] == "LFQ"
    assert row["intensity"] == 700.0
    assert row["additional_intensities"] is None


def _pg_rows_single_run(
    tmp_path: Path,
    report_rows: list[dict],
    *,
    matrix_rows: list[dict] | None = None,
    qvalue_threshold: float = 0.01,
) -> list[dict]:
    """Convert a synthetic single-run DIA-NN report+matrix and return pg rows."""
    from qpx.converters.diann.pg_adapter import DiannPgAdapter

    report_path = tmp_path / "report.tsv"
    pd.DataFrame(report_rows).to_csv(report_path, sep="\t", index=False)
    matrix_path = tmp_path / "pg_matrix.tsv"
    pd.DataFrame(matrix_rows or [{"Protein.Group": "P1", "Protein.Names": "N1", "Genes": "G1", "run_A": 900.0}]).to_csv(
        matrix_path, sep="\t", index=False
    )
    output_path = tmp_path / "out.pg.parquet"
    with DiannPgAdapter() as adapter:
        adapter.convert(
            diann_report=str(report_path),
            pg_matrix_path=str(matrix_path),
            output_path=str(output_path),
            qvalue_threshold=qvalue_threshold,
        )
    return pq.read_table(output_path).to_pylist()


def test_diann_pg_qvalue_filter_excludes_subthreshold_groups(tmp_path):
    """PG view must drop protein groups above the PG q-value threshold.

    Regression for bigbio/qpx#241: the pg WHERE clause filtered nothing, so every
    protein group was emitted regardless of PG-level FDR, inflating gene/protein
    counts on loosely-filtered reports. The filter mirrors the feature view.
    """
    base = {
        "Run": "run_A",
        "Precursor.Charge": 2,
        "Q.Value": 0.002,
        "PG.Q.Value": 0.002,
        "GG.Q.Value": 0.004,
        "Proteotypic": 1,
        "Precursor.Quantity": 100.0,
        "PG.Quantity": 1000.0,
        "PG.MaxLFQ": 900.0,
    }
    rows = [
        # Confident group (Global.PG.Q.Value = 0.003 <= 0.01): kept
        {
            **base,
            "Protein.Group": "P1",
            "Protein.Names": "N1",
            "Genes": "G1",
            "Global.PG.Q.Value": 0.003,
            "Stripped.Sequence": "PEPTIDEK",
            "Precursor.Id": "PEPTIDEK2",
        },
        # Sub-threshold group (Global.PG.Q.Value = 0.50 > 0.01): dropped
        {
            **base,
            "Protein.Group": "P2",
            "Protein.Names": "N2",
            "Genes": "G2",
            "Global.PG.Q.Value": 0.50,
            "Stripped.Sequence": "PEPTIDER",
            "Precursor.Id": "PEPTIDER2",
        },
    ]
    matrix = [
        {"Protein.Group": "P1", "Protein.Names": "N1", "Genes": "G1", "run_A": 900.0},
        {"Protein.Group": "P2", "Protein.Names": "N2", "Genes": "G2", "run_A": 800.0},
    ]
    out = _pg_rows_single_run(tmp_path, rows, matrix_rows=matrix, qvalue_threshold=0.01)
    accessions = {tuple(row["pg_accessions"]) for row in out}
    assert accessions == {("P1",)}, "sub-threshold protein group P2 must be filtered out"


def test_diann_pg_qvalue_filter_skipped_without_qvalue_columns(tmp_path, caplog):
    """No PG q-value column -> skip the filter with a warning, do not crash.

    Older DIA-NN reports may omit Global.PG.Q.Value / PG.Q.Value entirely
    (bigbio/qpx#241 guard).
    """
    import logging

    rows = [
        {
            "Run": "run_A",
            "Protein.Group": "P1",
            "Protein.Names": "N1",
            "Genes": "G1",
            "Precursor.Charge": 2,
            "Proteotypic": 1,
            "Precursor.Quantity": 100.0,
            "PG.Quantity": 1000.0,
            "PG.MaxLFQ": 900.0,
            "Stripped.Sequence": "PEPTIDEK",
            "Precursor.Id": "PEPTIDEK2",
        }
    ]
    with caplog.at_level(logging.WARNING):
        out = _pg_rows_single_run(tmp_path, rows, qvalue_threshold=0.01)
    assert len(out) == 1, "group must still be emitted when no PG q-value filter can be applied"
    assert out[0]["pg_qvalue"] is None
    assert any("skipping PG-level q-value filter" in rec.message for rec in caplog.records)


def test_diann_pg_conflicting_quantities_are_deterministic(tmp_path):
    """Merged rows disagreeing on PG.Quantity must yield a deterministic value.

    Regression for the Codex HIGH finding: the group's quantity was taken from
    whichever row appeared first, and the source SQL has no ORDER BY, so the
    emitted intensity was nondeterministic when annotation-inconsistent rows
    carried different quantities. The adapter now takes the max of the finite
    values, so row order cannot change the result.
    """
    base = {
        "Run": "run_A",
        "Protein.Group": "P1;P2",
        "Protein.Names": "N1;N2",
        "Precursor.Charge": 2,
        "Q.Value": 0.002,
        "PG.Q.Value": 0.002,
        "Global.PG.Q.Value": 0.003,
        "GG.Q.Value": 0.004,
        "Proteotypic": 1,
        "Precursor.Quantity": 100.0,
    }
    # Same protein group + run, inconsistent Genes annotation, DIFFERENT quantities.
    rows = [
        {
            **base,
            "Genes": "G1;G2",
            "PG.Quantity": 1000.0,
            "PG.MaxLFQ": 900.0,
            "Stripped.Sequence": "PEPTIDEK",
            "Precursor.Id": "PEPTIDEK2",
        },
        {
            **base,
            "Genes": "G1",
            "PG.Quantity": 3000.0,
            "PG.MaxLFQ": 2500.0,
            "Stripped.Sequence": "PEPTIDER",
            "Precursor.Id": "PEPTIDER2",
        },
    ]
    matrix = [{"Protein.Group": "P1;P2", "Protein.Names": "N1;N2", "Genes": "G1;G2", "run_A": 900.0}]

    out = _pg_rows_single_run(tmp_path, rows, matrix_rows=matrix)
    # Reversed row order must give the identical result.
    out_reversed = _pg_rows_single_run(tmp_path, list(reversed(rows)), matrix_rows=matrix)

    assert len(out) == 1 and len(out_reversed) == 1
    assert out[0]["intensity"] == 3000.0, "primary intensity must be the deterministic max of the finite values"
    assert out[0]["intensity"] == out_reversed[0]["intensity"], "result must not depend on row order"
    assert out[0]["additional_intensities"][0]["intensities"][0]["intensity_value"] == 2500.0


def test_diann_feature_observed_mz_is_null_without_mz_source(tmp_path):
    """observed_mz must be NULL (not 0.0) when the report lacks Precursor.Mz and
    no mzml_info_folder is provided.

    Regression for bigbio/qpx#244: a sentinel 0.0 looks like a real measurement
    and corrupts downstream mass-error / PPM math; a truthful NULL propagates as
    "unknown".
    """
    from qpx.converters.diann.feature_adapter import DiannFeatureAdapter

    report_path = tmp_path / "no_mz_report.tsv"
    pd.DataFrame(
        [
            {
                "Run": "run_A",
                "Protein.Group": "P1",
                "Genes": "G1",
                "Modified.Sequence": "PEPTIDEK",
                "Stripped.Sequence": "PEPTIDEK",
                "Precursor.Charge": 2,
                "Q.Value": 0.002,
                "Precursor.Quantity": 100.0,
                "RT": 10.0,
            }
        ]
    ).to_csv(report_path, sep="\t", index=False)

    output_path = tmp_path / "no_mz.feature.parquet"
    with DiannFeatureAdapter() as adapter:
        adapter.convert(
            diann_report=str(report_path),
            output_path=str(output_path),
            qvalue_threshold=0.01,
        )

    rows = pq.read_table(output_path).to_pylist()
    assert len(rows) == 1
    assert rows[0]["observed_mz"] is None, "absent observed_mz must be NULL, not a 0.0 sentinel"


def test_diann_feature_warns_on_runs_missing_ms_info(tmp_path, caplog):
    """A report Run without a matching *_ms_info.parquet must be warned about.

    Regression for bigbio/qpx#244: such runs are silently dropped from
    feature.parquet; the drop is unchanged but now visible in the logs.
    """
    import logging

    from qpx.converters.diann.feature_adapter import DiannFeatureAdapter

    shared = {
        "Protein.Group": "P1",
        "Genes": "G1",
        "Modified.Sequence": "PEPTIDEK",
        "Stripped.Sequence": "PEPTIDEK",
        "Precursor.Charge": 2,
        "Q.Value": 0.002,
        "Precursor.Quantity": 100.0,
        "RT": 10.0,
    }
    report_path = tmp_path / "two_run_report.tsv"
    pd.DataFrame(
        [
            {**shared, "Run": "run_A"},
            {**shared, "Run": "run_B"},  # no ms_info file for run_B
        ]
    ).to_csv(report_path, sep="\t", index=False)

    # Provide an ms_info file for run_A only.
    ms_info_dir = tmp_path / "ms_info"
    ms_info_dir.mkdir()
    pq.write_table(
        pa.Table.from_pandas(pd.DataFrame({"rt": [10.0], "scan": [1], "precursor_mz": [500.0]})),
        str(ms_info_dir / "run_A_ms_info.parquet"),
    )

    output_path = tmp_path / "two_run.feature.parquet"
    with caplog.at_level(logging.WARNING):
        with DiannFeatureAdapter() as adapter:
            adapter.convert(
                diann_report=str(report_path),
                output_path=str(output_path),
                mzml_info_folder=str(ms_info_dir),
                qvalue_threshold=0.01,
            )

    run_names = set(pq.read_table(output_path).column("run_file_name").to_pylist())
    assert run_names == {"run_A"}, "run_B has no ms_info and is dropped (unchanged behaviour)"
    assert any("run_B" in rec.message and "dropped" in rec.message for rec in caplog.records), (
        "dropped run must be reported in a warning"
    )


# ---------------------------------------------------------------------------
# Ontology conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestDiaNNOntologyConversion:
    """Validate ontology.parquet output from DIA-NN conversion.

    The ontology file is only written when scores are discovered during
    conversion. It may legitimately not exist if no scores were tracked.
    """

    def test_ontology_file_exists_or_no_scores(self, converted_output):
        path = converted_output / "diann_test.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")

        table = pq.read_table(str(path))
        assert table.num_rows > 0, "ontology.parquet exists but is empty"

    def test_ontology_columns(self, converted_output):
        path = converted_output / "diann_test.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")

        table = pq.read_table(str(path))
        column_names = set(table.column_names)
        expected = {"field_name", "view"}
        missing = expected - column_names
        assert not missing, f"Missing ontology columns: {missing}"

    def test_schema_validation(self, converted_output):
        from qpx.core.data import OntologySchema

        path = converted_output / "diann_test.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")
        table = pq.read_table(str(path))
        errors = OntologySchema.validate(table)
        assert not errors, f"Schema validation errors: {errors}"


class TestDiannReportLoading:
    """Adaptive TABLE-vs-VIEW report loading (bigbio/qpx#221 — @Shen-YuFei)."""

    @staticmethod
    def _adapter(conn):
        import logging

        from qpx.converters.diann.base_adapter import DiaNNBaseAdapter

        class _A(DiaNNBaseAdapter):
            def convert(self, *args, **kwargs):
                pass

        adapter = _A.__new__(_A)  # bypass BaseConverter.__init__; we only need _conn + logger
        adapter._conn = conn
        adapter.logger = logging.getLogger("test")
        return adapter

    @staticmethod
    def _relation_type(conn):
        return conn.execute("SELECT table_type FROM information_schema.tables WHERE table_name='report'").fetchone()[0]

    def test_parse_duckdb_size(self):
        from qpx.converters.diann.base_adapter import _parse_duckdb_size

        assert _parse_duckdb_size("14.3 GiB") == int(14.3 * 1024**3)
        assert _parse_duckdb_size("8192MB") == 8192 * 1000**2
        assert _parse_duckdb_size("200 MiB") == 200 * 1024**2
        assert _parse_duckdb_size("-1") is None
        assert _parse_duckdb_size("garbage") is None

    def test_should_materialize_decision(self, tmp_path, monkeypatch):
        import gzip
        import os

        from qpx.core.engine import DuckDBEngine

        report = tmp_path / "r.tsv"
        report.write_text("Run\tX\nr1\t1\n", encoding="utf-8")
        compressed_report = tmp_path / "r.tsv.gz"
        with gzip.open(compressed_report, "wt", encoding="utf-8") as handle:
            handle.write("Run\tX\nr1\t1\n")
        with DuckDBEngine() as engine:
            adapter = self._adapter(engine.connection)
            should_materialize = getattr(adapter, "_should_materialize_report")

            monkeypatch.setattr(adapter, "_duckdb_memory_limit_bytes", lambda: 8 * 1024**3)
            assert should_materialize(str(report)) is True  # tiny file, big limit
            assert should_materialize(str(compressed_report)) is False

            monkeypatch.setattr(os.path, "getsize", lambda _p: 5 * 1024**3)  # 5 GB file
            assert should_materialize(str(report)) is False  # over 0.5*limit -> VIEW

            monkeypatch.setattr(adapter, "_duckdb_memory_limit_bytes", lambda: None)
            assert should_materialize(str(report)) is False  # unknown limit, big file

    def test_small_report_materializes_table(self, tmp_path):
        from qpx.core.engine import DuckDBEngine

        report = tmp_path / "report.tsv"
        report.write_text("Run\tX\n" + "\n".join(f"r1\t{i}" for i in range(20)) + "\n", encoding="utf-8")
        engine = DuckDBEngine()
        self._adapter(engine._conn)._load_diann_report(str(report))
        assert self._relation_type(engine._conn) == "BASE TABLE"

    def test_falls_back_to_view_when_over_budget(self, tmp_path, monkeypatch):
        from qpx.core.engine import DuckDBEngine

        report = tmp_path / "report.tsv"
        report.write_text("Run\tX\n" + "\n".join(f"r1\t{i}" for i in range(20)) + "\n", encoding="utf-8")
        engine = DuckDBEngine()
        adapter = self._adapter(engine._conn)
        monkeypatch.setattr(adapter, "_should_materialize_report", lambda _p: False)
        adapter._load_diann_report(str(report))
        assert self._relation_type(engine._conn) == "VIEW"

    def test_materialization_oom_falls_back_without_scanning_view(self, tmp_path, monkeypatch):
        """OOM fallback creates a VIEW without immediately scanning it."""
        import duckdb

        from qpx.core.engine import DuckDBEngine

        report = tmp_path / "report.tsv"
        report.write_text("Run\tX\nr1\t1\n", encoding="utf-8")
        with DuckDBEngine() as engine:
            statements = []
            real_execute = engine.connection.execute

            def execute(statement, *args, **kwargs):
                """Force materialization to fail while preserving other queries."""
                statements.append(statement)
                if statement.startswith("CREATE TABLE report"):
                    raise duckdb.OutOfMemoryException("forced materialization failure")
                return real_execute(statement, *args, **kwargs)

            conn = Mock(wraps=engine.connection)
            conn.execute.side_effect = execute
            adapter = self._adapter(conn)
            monkeypatch.setattr(adapter, "_should_materialize_report", lambda _p: True)

            getattr(adapter, "_load_diann_report")(str(report))

            assert self._relation_type(engine.connection) == "VIEW"
            assert not any(statement.startswith("SELECT COUNT(*)") for statement in statements)
