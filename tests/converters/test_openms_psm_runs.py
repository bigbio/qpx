"""Native PSM run recovery requires exact, unambiguous spectrum evidence."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.converters.openms.converter import OpenMSConverter
from qpx.converters.openms.psm_runs import build_psm_run_lookup
from qpx.core.data import PsmSchema
from tests.conftest import make_psm_record


def _write_xml(tmp_path, references, *, assigned=False):
    """Create two real run headers and minimal assigned or unassigned PIDs."""
    header = '<consensusXML><mapList count="2">'
    header += '<map id="0" name="run_a.mzML" label=""/><map id="1" name="run_b.mzML" label=""/></mapList>'
    nodes = []
    for run, reference, rt, mz in references:
        tag = "PeptideIdentification" if assigned else "UnassignedPeptideIdentification"
        meta = f'<UserParam name="map_index" value="{run}"/>' if run is not None and not assigned else ""
        pid = (
            f'<{tag} MZ="{mz}" RT="{rt}" spectrum_reference="{reference}">{meta}'
            '<PeptideHit sequence="PEPTIDEK" charge="2" score="0.01"/>'
            f"</{tag}>"
        )
        if assigned:
            pid = (
                '<consensusElement charge="2"><centroid mz="500" rt="10"/>'
                f'<groupedElementList><element map="{run}" rt="10" it="100"/></groupedElementList>'
                f"{pid}</consensusElement>"
            )
        nodes.append(pid)
    body = "".join(nodes)
    if assigned:
        body = f"<consensusElementList>{body}</consensusElementList>"
    source = tmp_path / "source.consensusXML"
    source.write_text(header + body + "</consensusXML>")
    return source


def _table(rows, *, rt_type=pa.float32(), mz_type=pa.float32()):
    """Give source rows opaque identities and unrelated evidence to preserve."""
    records = [
        {
            "peptidoform": "PEPTIDEK",
            "charge": 2,
            "scan": [42],
            "rt": 10.123456789,
            "observed_mz": 500.123456789,
            "run_file_name": "run_a.mzML",
            "psm_id": 900 + index,
            "feature_id": 700 + index,
            "protein_accessions": ["P12345", "P54321"],
            **values,
        }
        for index, values in enumerate(rows)
    ]
    table = pa.Table.from_pylist(records)
    for column, dtype in (("rt", rt_type), ("observed_mz", mz_type)):
        table = table.set_column(table.schema.get_field_index(column), column, table.column(column).cast(dtype))
    return table.replace_schema_metadata({b"producer": b"OpenMS", b"scan_format": b"scan"})


@pytest.mark.parametrize("assigned", [False, True])
def test_same_scan_in_different_runs_preserves_all_other_data(tmp_path, assigned):
    """RT distinguishes the run while producer IDs and foreign keys stay opaque."""
    source = _write_xml(
        tmp_path,
        [(0, "scan=42", 10.123456789, 500.123456789), (1, "scan=42", 20.123456789, 500.123456789)],
        assigned=assigned,
    )
    original = _table([{}, {"rt": 20.123456789}])
    restored = build_psm_run_lookup(source).apply(original)

    assert restored.column("run_file_name").to_pylist() == ["run_a", "run_b"]
    assert restored.drop_columns("run_file_name").equals(original.drop_columns("run_file_name"), check_metadata=True)
    assert restored.schema == original.schema
    assert original.column("run_file_name").to_pylist() == ["run_a.mzML", "run_a.mzML"]


def test_complete_native_components_distinguish_same_scan(tmp_path):
    """Frame components cannot be truncated during run recovery."""
    source = _write_xml(
        tmp_path,
        [(0, "frame=10 scan=42", 10.123456789, 500.123456789), (1, "frame=11 scan=42", 10.123456789, 500.123456789)],
    )
    table = _table([{"scan": [10, 42]}, {"scan": [11, 42]}])
    restored = build_psm_run_lookup(source).apply(table)
    assert restored.column("run_file_name").to_pylist() == ["run_a", "run_b"]
    assert restored.column("scan").equals(table.column("scan"))


def test_missing_match_is_not_guessed_from_existing_run(tmp_path):
    """An existing run label cannot authorize a mismatching XML spectrum."""
    source = _write_xml(tmp_path, [(0, "scan=42", 10.123456789, 500.123456789)])
    with pytest.raises(ValueError, match="no exact companion XML match"):
        build_psm_run_lookup(source).apply(_table([{"scan": [99]}]))


def test_ambiguous_runs_are_not_disambiguated_using_native_run(tmp_path):
    """The native run may be precisely the corrupted field being repaired."""
    source = _write_xml(
        tmp_path,
        [(0, "scan=42", 10.123456789, 500.123456789), (1, "scan=42", 10.123456789, 500.123456789)],
    )
    with pytest.raises(ValueError, match="ambiguous or unresolved"):
        build_psm_run_lookup(source).apply(_table([{}]))


@pytest.mark.parametrize(
    "rt_type,mz_type", [(pa.float32(), pa.float32()), (pa.float64(), pa.float64()), (pa.float64(), pa.float32())]
)
def test_match_uses_actual_source_column_precision(tmp_path, rt_type, mz_type):
    """Decimal XML values match the exact representation stored by the producer."""
    source = _write_xml(tmp_path, [(1, "scan=42", 10.123456789, 500.123456789)])
    table = _table([{}], rt_type=rt_type, mz_type=mz_type)
    restored = build_psm_run_lookup(source).apply(table)
    assert restored.column("run_file_name").to_pylist() == ["run_b"]
    assert restored.drop_columns("run_file_name").equals(table.drop_columns("run_file_name"), check_metadata=True)


def test_precision_loss_that_merges_runs_is_ambiguous(tmp_path):
    """Float32 collisions must not reuse a more precise cached lookup."""
    source = _write_xml(tmp_path, [(0, "scan=42", 1.00000001, 500), (1, "scan=42", 1.00000002, 500)])
    resolver = build_psm_run_lookup(source)
    table = _table([{"rt": 1.00000001, "observed_mz": 500.0}], rt_type=pa.float64())
    assert resolver.apply(table).column("run_file_name").to_pylist() == ["run_a"]
    with pytest.raises(ValueError, match="ambiguous or unresolved"):
        resolver.apply(_table([{"rt": 1.00000001, "observed_mz": 500.0}]))


@pytest.mark.parametrize("values", [{"rt": None}, {"observed_mz": float("nan")}, {"scan": []}])
def test_incomplete_signature_fails_explicitly(tmp_path, values):
    """Missing numerical evidence must not silently preserve a possibly wrong run."""
    source = _write_xml(tmp_path, [(0, "scan=42", 10.123456789, 500.123456789)])
    with pytest.raises(ValueError, match="incomplete numerical spectrum signature"):
        build_psm_run_lookup(source).apply(_table([values]))


def test_unresolved_xml_run_fails_explicitly(tmp_path):
    """A matching spectrum without run attribution is not a successful recovery."""
    source = _write_xml(tmp_path, [(None, "scan=42", 10.123456789, 500.123456789)])
    with pytest.raises(ValueError, match="ambiguous or unresolved"):
        build_psm_run_lookup(source).apply(_table([{}]))


def test_header_only_companion_preserves_label_only_usage(tmp_path):
    """Existing maplist-only files provide channels without promising PID evidence."""
    assert build_psm_run_lookup(_write_xml(tmp_path, [])) is None


def test_hitless_pid_is_not_treated_as_header_only(tmp_path):
    """Present but unusable identification evidence cannot disable validation."""
    source = _write_xml(tmp_path, [(0, "scan=42", 10.123456789, 500.123456789)])
    source.write_text(source.read_text().replace('<PeptideHit sequence="PEPTIDEK" charge="2" score="0.01"/>', ""))
    resolver = build_psm_run_lookup(source)
    assert resolver is not None
    with pytest.raises(ValueError, match="no exact companion XML match"):
        resolver.apply(_table([{}]))


def test_repeated_evidence_in_the_same_run_is_not_ambiguous(tmp_path):
    """Several references to the same run are one unambiguous attribution."""
    reference = (1, "scan=42", 10.123456789, 500.123456789)
    source = _write_xml(tmp_path, [reference, reference])
    assert build_psm_run_lookup(source).apply(_table([{}])).column("run_file_name").to_pylist() == ["run_b"]


def _write_native_psms(tmp_path, rows, *, provided_ids):
    """Persist complete legacy-like rows without passing through an ID writer."""
    records = [{**make_psm_record(), **row} for row in _table(rows).to_pylist()]
    table = pa.Table.from_pylist(records, schema=PsmSchema.get_arrow_schema())
    if not provided_ids:
        table = table.drop_columns(["psm_id", "feature_id"])
    source = tmp_path / "native"
    source.mkdir()
    pq.write_table(table, source / "native.psm.parquet", row_group_size=1)
    return source, table


@pytest.mark.parametrize("provided_ids", [False, True])
def test_public_enrichment_recovers_runs_before_assigning_identity(tmp_path, provided_ids):
    """Recovered runs separate derived IDs while supplied IDs and links survive."""
    source, before = _write_native_psms(tmp_path, [{}, {"rt": 20.123456789}], provided_ids=provided_ids)
    companion = _write_xml(
        tmp_path,
        [(0, "scan=42", 10.123456789, 500.123456789), (1, "scan=42", 20.123456789, 500.123456789)],
    )
    output = tmp_path / "out"
    OpenMSConverter(source, consensusxml_path=companion).convert(output)
    after = pq.read_table(output / "openms.psm.parquet")

    assert after.column("run_file_name").to_pylist() == ["run_a", "run_b"]
    assert len(set(after.column("psm_id").to_pylist())) == 2
    if provided_ids:
        assert after.column("psm_id").equals(before.column("psm_id"))
        assert after.column("feature_id").equals(before.column("feature_id"))
    preserved = [name for name in before.column_names if name != "run_file_name"]
    assert after.select(preserved).equals(before.select(preserved), check_metadata=False)


def test_public_enrichment_preserves_recovered_run_stem_matching_sdrf(tmp_path):
    """A recovered stem ending in .raw must not lose a second extension."""
    source, _ = _write_native_psms(tmp_path, [{}], provided_ids=False)
    companion = _write_xml(tmp_path, [(0, "scan=42", 10.123456789, 500.123456789)])
    companion.write_text(companion.read_text().replace('name="run_a.mzML"', 'name="/acquisition/sample.raw.mzML"'))
    sdrf = tmp_path / "input.sdrf.tsv"
    sdrf.write_text(
        "source name\tcharacteristics[organism]\tcharacteristics[organism part]\tcomment[data file]\tcomment[label]\n"
        "sample_1\tHomo sapiens\tliver\t/data/sample.raw.mzML\tlabel free sample\n"
    )
    output = tmp_path / "out"

    OpenMSConverter(source, sdrf_path=sdrf, consensusxml_path=companion).convert(output)

    psm_runs = pq.read_table(output / "openms.psm.parquet").column("run_file_name").to_pylist()
    sdrf_runs = pq.read_table(output / "openms.run.parquet").column("run_file_name").to_pylist()
    assert psm_runs == sdrf_runs == ["sample.raw"]


@pytest.mark.parametrize("failure", ["missing", "ambiguous"])
def test_public_enrichment_match_failure_preserves_existing_output(tmp_path, failure):
    """Failure after a successful row group cannot publish or leak staged files."""
    source, _ = _write_native_psms(tmp_path, [{}, {"scan": [43], "rt": 20.123456789}], provided_ids=True)
    references = [(0, "scan=42", 10.123456789, 500.123456789)]
    if failure == "ambiguous":
        references.extend([(0, "scan=43", 20.123456789, 500.123456789), (1, "scan=43", 20.123456789, 500.123456789)])
    companion = _write_xml(tmp_path, references)
    output = tmp_path / "out"
    output.mkdir()
    destination = output / "openms.psm.parquet"
    previous = b"existing output must survive a failed enrichment"
    destination.write_bytes(previous)

    message = "no exact companion XML match" if failure == "missing" else "ambiguous or unresolved"
    with pytest.raises(ValueError, match=message):
        OpenMSConverter(source, consensusxml_path=companion).convert(output)

    assert destination.read_bytes() == previous
    assert list(output.iterdir()) == [destination]


def test_cli_reports_a_match_failure_as_a_usage_error(tmp_path):
    """`convert openms` turns companion-match failures into a CLI error, not a traceback."""
    from click.testing import CliRunner

    from qpx.cli.main import qpx_main

    source, _ = _write_native_psms(tmp_path, [{"scan": [43]}], provided_ids=True)
    companion = _write_xml(tmp_path, [(0, "scan=42", 10.123456789, 500.123456789)])
    sdrf = tmp_path / "input.sdrf.tsv"
    sdrf.write_text(
        "source name\tcharacteristics[organism]\tcharacteristics[organism part]\tcomment[data file]\tcomment[label]\n"
        "sample_1\tHomo sapiens\tliver\trun_a.mzML\tlabel free sample\n"
    )
    result = CliRunner().invoke(
        qpx_main,
        [
            *("convert", "openms", "--qpx-dir", str(source), "--sdrf-file", str(sdrf)),
            *("--consensusxml", str(companion), "--output-folder", str(tmp_path / "out")),
        ],
    )

    assert result.exit_code == 1
    assert "no exact companion XML match" in result.output
    assert not isinstance(result.exception, ValueError)
