"""Fill protein properties from an optional FASTA."""

import gzip

import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.transform import transform
from qpx.transforms.protein_properties import FastaSequences

P12345 = "MPEPTIDEKAAANOTHERKGGG"  # PEPTIDEK at 2-9, ANOTHERK at 12-19
P67890 = "MMTHIRDPEPSSS"  # THIRDPEP at 3-10


def _fasta(tmp_path, text, name="db.fasta"):
    path = tmp_path / name
    path.write_text(text)
    return path


def _dataset_fasta(tmp_path):
    return _fasta(
        tmp_path,
        f">sp|P12345|PROT1_HUMAN Protein one GN=ONE\n{P12345[:10]}\n{P12345[10:]}\n"
        f">sp|P67890|PROT2_HUMAN Protein two GN=TWO\n{P67890}\n"
        f">DECOY_sp|P12345|PROT1_HUMAN\n{P12345[::-1]}\n",
    )


class TestFastaSequences:
    def test_reachable_by_full_identifier_and_by_accession(self, tmp_path):
        fasta = FastaSequences.from_path(_fasta(tmp_path, ">sp|P12345|NAME desc\nMPEP\nTIDE*\n"))
        assert fasta.get("sp|P12345|NAME") == "MPEPTIDE"
        assert fasta.get("P12345") == "MPEPTIDE"

    def test_bare_identifier_header(self, tmp_path):
        assert FastaSequences.from_path(_fasta(tmp_path, ">P12345\nMPEP\n")).get("P12345") == "MPEP"

    def test_decoy_entries_do_not_erase_the_target(self, tmp_path):
        """DECOY_sp|P12345|.. shares the accession field with the target; indexing it
        would turn P12345 into a conflict and lose the real sequence."""
        fasta = FastaSequences.from_path(_fasta(tmp_path, ">sp|P12345|N\nMPEPTIDE\n>DECOY_sp|P12345|N\nEDITPEPM\n"))
        assert fasta.get("P12345") == "MPEPTIDE"
        assert fasta.decoy_entries == 1

    def test_conflicting_sequences_resolve_to_nothing(self, tmp_path):
        fasta = FastaSequences.from_path(_fasta(tmp_path, ">sp|P12345|A\nMPEP\n>tr|P12345|B\nMKKK\n"))
        assert fasta.get("P12345") is None
        assert fasta.get("sp|P12345|A") == "MPEP"

    def test_absent_protein_is_none(self, tmp_path):
        assert FastaSequences.from_path(_fasta(tmp_path, ">sp|P12345|A\nMPEP\n")).get("Q99999") is None

    def test_gzip(self, tmp_path):
        path = tmp_path / "db.fasta.gz"
        with gzip.open(path, "wt") as handle:
            handle.write(">sp|P12345|A\nMPEP\n")
        assert FastaSequences.from_path(path).get("P12345") == "MPEP"


def _run(dataset_dir, fasta, *extra):
    result = CliRunner().invoke(transform, ["protein-properties", "--dataset", str(dataset_dir), "--fasta", str(fasta), *extra])
    assert result.exit_code == 0, result.output
    return result


def _rows(dataset_dir, view):
    return pq.read_table(dataset_dir / f"exp.{view}.parquet").to_pylist()


def test_fills_pg_properties_and_feature_positions(dataset_dir, tmp_path):
    before_ids = [r["feature_id"] for r in _rows(dataset_dir, "feature")]
    _run(dataset_dir, _dataset_fasta(tmp_path), "--in-place")

    pg = {r["anchor_protein"]: r for r in _rows(dataset_dir, "pg")}
    # P12345: PEPTIDEK (2-9) + ANOTHERK (12-19) = 16 of 22 residues
    assert pg["P12345"]["sequence_coverage"] == pytest.approx(100 * 16 / 22, rel=1e-5)
    assert pg["P12345"]["molecular_weight"] > 0
    # THIRDPEP is grouped to P12345, not P67890: no peptide maps to P67890, so no
    # coverage is invented, but its mass is still known from the FASTA.
    assert pg["P67890"]["sequence_coverage"] is None
    assert pg["P67890"]["molecular_weight"] > 0

    features = {r["sequence"]: r for r in _rows(dataset_dir, "feature")}
    assert features["PEPTIDEK"]["pg_positions"] == [{"protein_accession": "P12345", "start": 2, "end": 9}]
    assert features["ANOTHERK"]["pg_positions"] == [{"protein_accession": "P12345", "start": 12, "end": 19}]
    # Its group's protein does not contain it: no position from another protein.
    assert features["THIRDPEP"]["pg_positions"] is None
    assert [r["feature_id"] for r in _rows(dataset_dir, "feature")] == before_ids


def test_keeps_the_qpx_footer_identity(dataset_dir, tmp_path):
    before = pq.ParquetFile(dataset_dir / "exp.feature.parquet").schema_arrow.metadata
    _run(dataset_dir, _dataset_fasta(tmp_path), "--in-place")
    after = pq.ParquetFile(dataset_dir / "exp.feature.parquet").schema_arrow.metadata
    for key in (b"identity_composite", b"primary_key", b"file_type", b"qpx_version"):
        assert after.get(key) == before.get(key), key
    assert after[b"uuid"] != before[b"uuid"]


def test_the_result_still_validates(dataset_dir, tmp_path):
    from qpx.cli.validate import validate_cmd

    _run(dataset_dir, _dataset_fasta(tmp_path), "--in-place")
    result = CliRunner().invoke(validate_cmd, ["--dataset-path", str(dataset_dir)])

    assert "INVALID" not in result.output, result.output
    assert "Overall: VALID" in result.output, result.output


def test_never_overwrites_a_recorded_value(dataset_dir, tmp_path):
    from qpx.writers import PgWriter
    from tests.conftest import make_pg_record

    record = make_pg_record(anchor_protein="P12345", run_file_name="run_01")
    record["sequence_coverage"] = 12.5
    record["molecular_weight"] = 99.0
    with PgWriter(dataset_dir / "exp.pg.parquet") as writer:
        writer.write_batch([record])

    _run(dataset_dir, _dataset_fasta(tmp_path), "--in-place")

    row = _rows(dataset_dir, "pg")[0]
    assert row["sequence_coverage"] == pytest.approx(12.5)
    assert row["molecular_weight"] == pytest.approx(99.0)


def test_proteins_missing_from_the_fasta_stay_null_and_are_reported(dataset_dir, tmp_path):
    """A DIA-NN internal decoy or a protein from another database is not in the FASTA."""
    fasta = _fasta(tmp_path, f">sp|P67890|PROT2_HUMAN\n{P67890}\n")
    result = _run(dataset_dir, fasta, "--in-place")

    pg = {r["anchor_protein"]: r for r in _rows(dataset_dir, "pg")}
    assert pg["P12345"]["sequence_coverage"] is None
    assert pg["P12345"]["molecular_weight"] is None
    assert all(r["pg_positions"] is None for r in _rows(dataset_dir, "feature"))
    assert "anchors not in the FASTA" in result.output
    assert "P12345" in result.output


def test_decoy_rows_are_not_annotated(dataset_dir, tmp_path):
    from qpx.writers import PgWriter
    from tests.conftest import make_pg_record

    with PgWriter(dataset_dir / "exp.pg.parquet") as writer:
        writer.write_batch([make_pg_record(anchor_protein="P12345", run_file_name="run_01", is_decoy=True)])

    _run(dataset_dir, _dataset_fasta(tmp_path), "--in-place")

    row = _rows(dataset_dir, "pg")[0]
    assert row["sequence_coverage"] is None
    assert row["molecular_weight"] is None


def test_records_a_provenance_step_with_the_fasta_checksum(dataset_dir, tmp_path):
    import hashlib

    fasta = _dataset_fasta(tmp_path)
    _run(dataset_dir, fasta, "--in-place")

    steps = pq.read_table(dataset_dir / "exp.provenance.parquet").to_pylist()
    step = steps[-1]
    assert step["step_name"] == "protein_properties_from_fasta"
    assert step["step_order"] == len(steps)
    params = {p["key"]: p["value"] for p in step["parameters"]}
    assert params["fasta_sha256"] == hashlib.sha256(fasta.read_bytes()).hexdigest()
    assert params["feature_pg_positions_filled"] == "2"


def test_output_folder_leaves_the_source_untouched(dataset_dir, tmp_path):
    before = (dataset_dir / "exp.pg.parquet").read_bytes()
    out = tmp_path / "annotated"
    _run(dataset_dir, _dataset_fasta(tmp_path), "--output-folder", str(out))

    assert (dataset_dir / "exp.pg.parquet").read_bytes() == before
    assert pq.read_table(out / "exp.pg.parquet").to_pylist()[0]["molecular_weight"] is not None


@pytest.mark.parametrize("args", [[], ["--in-place", "--output-folder", "x"]])
def test_requires_exactly_one_destination(dataset_dir, tmp_path, args):
    result = CliRunner().invoke(
        transform,
        ["protein-properties", "--dataset", str(dataset_dir), "--fasta", str(_dataset_fasta(tmp_path)), *args],
    )
    assert result.exit_code != 0


def test_coverage_counts_shared_peptides_named_by_psm_evidence(dataset_dir, tmp_path):
    """Group membership alone drops peptides shared across groups. On PXD000612
    that put ACTB at 4.5% coverage against OpenMS's 84.3%. A PSM whose evidence
    names the anchor must count toward its coverage."""
    from qpx.writers import PsmWriter
    from tests.conftest import make_psm_record

    # THIRDPEP is grouped to P12345 in the feature view, but its PSM evidence
    # also names P67890, whose sequence contains it (3-10 of 13).
    psm = make_psm_record(sequence="THIRDPEP", run_file_name="run_02")
    psm["protein_accessions"] = ["P12345", "P67890"]
    with PsmWriter(dataset_dir / "exp.psm.parquet") as writer:
        writer.write_batch([psm])

    _run(dataset_dir, _dataset_fasta(tmp_path), "--in-place")

    pg = {r["anchor_protein"]: r for r in _rows(dataset_dir, "pg")}
    assert pg["P67890"]["sequence_coverage"] == pytest.approx(100 * 8 / 13, rel=1e-5)
    features = {r["sequence"]: r for r in _rows(dataset_dir, "feature")}
    # Positions stay tied to the feature's own group (P12345), not the PSM evidence.
    assert features["THIRDPEP"]["pg_positions"] is None


@pytest.mark.parametrize("direction", ["forward", "reverse"])
def test_coverage_keeps_psm_evidence_after_accession_normalization(dataset_dir, tmp_path, direction):
    """Normalizing PG/feature identifiers must retain shared-peptide PSM evidence."""
    from qpx.writers import PsmWriter
    from tests.conftest import make_psm_record

    fasta = _dataset_fasta(tmp_path)
    full = "sp|P67890|PROT2_HUMAN"
    if direction == "forward":
        preparation = CliRunner().invoke(
            transform,
            [
                "normalize-accessions",
                "--dataset",
                str(dataset_dir),
                "--direction",
                "reverse",
                "--fasta",
                str(fasta),
                "--in-place",
            ],
        )
        assert preparation.exit_code == 0, preparation.output
    psm = make_psm_record(sequence="THIRDPEP", run_file_name="run_02")
    psm["protein_accessions"] = [full if direction == "forward" else "P67890"]
    with PsmWriter(dataset_dir / "exp.psm.parquet") as writer:
        writer.write_batch([psm])
    psm_before = (dataset_dir / "exp.psm.parquet").read_bytes()

    result = CliRunner().invoke(
        transform,
        ["normalize-accessions", "--dataset", str(dataset_dir), "--direction", direction, "--fasta", str(fasta), "--in-place"],
    )
    assert result.exit_code == 0, result.output
    anchors_before = [row["anchor_protein"] for row in _rows(dataset_dir, "pg")]
    _run(dataset_dir, fasta, "--in-place")

    pg = {row["anchor_protein"]: row for row in _rows(dataset_dir, "pg")}
    anchor = "P67890" if direction == "forward" else full
    assert pg[anchor]["sequence_coverage"] == pytest.approx(100 * 8 / 13, rel=1e-5)
    assert list(pg) == anchors_before
    assert (dataset_dir / "exp.psm.parquet").read_bytes() == psm_before


@pytest.mark.parametrize("in_place", [True, False])
def test_new_provenance_updates_structure_count(dataset_dir, tmp_path, in_place):
    """The final count includes new provenance once and excludes staging files."""
    from qpx import Dataset
    from qpx.writers import DatasetWriter

    with Dataset(dataset_dir) as dataset:
        record = dataset.dataset_meta.to_df().iloc[0].to_dict()
        record.update(dataset.compute_integrity())
    with DatasetWriter(dataset_dir / "exp.dataset.parquet") as writer:
        writer.write_batch([record])
    source_before = (dataset_dir / "exp.dataset.parquet").read_bytes()
    assert record["total_structures"] == 5
    output = dataset_dir if in_place else tmp_path / "annotated"
    destination = ["--in-place"] if in_place else ["--output-folder", str(output)]
    _run(dataset_dir, _dataset_fasta(tmp_path), *destination)

    with Dataset(output) as dataset:
        metadata = dataset.dataset_meta.to_df().iloc[0]
        assert metadata["total_structures"] == dataset.compute_integrity()["total_structures"] == 6
        assert dataset.verify_integrity() == {"errors": [], "warnings": []}
    if not in_place:
        assert (dataset_dir / "exp.dataset.parquet").read_bytes() == source_before

    _run(output, _dataset_fasta(tmp_path), "--in-place")
    assert _rows(output, "dataset")[0]["total_structures"] == 6


def test_coverage_does_not_merge_conflicting_fasta_accessions(dataset_dir, tmp_path):
    """Different sequences sharing a bare accession keep their own protein evidence."""
    from qpx.writers import PgWriter, PsmWriter
    from tests.conftest import make_pg_record, make_psm_record

    fasta = _fasta(tmp_path, ">sp|P12345|A\nMPEPTIDEK\n>tr|P12345|B\nMPEPTIDEKAAA\n")
    anchors = ["sp|P12345|A", "tr|P12345|B", "P12345"]
    with PgWriter(dataset_dir / "exp.pg.parquet") as writer:
        writer.write_batch([make_pg_record(anchor_protein=anchor, pg_accessions=[anchor]) for anchor in anchors])
    (dataset_dir / "exp.feature.parquet").unlink()
    psm = make_psm_record()
    psm["protein_accessions"] = [anchors[0]]
    with PsmWriter(dataset_dir / "exp.psm.parquet") as writer:
        writer.write_batch([psm])

    _run(dataset_dir, fasta, "--in-place")

    pg = {row["anchor_protein"]: row for row in _rows(dataset_dir, "pg")}
    assert pg[anchors[0]]["sequence_coverage"] == pytest.approx(100 * 8 / 9, rel=1e-5)
    assert pg[anchors[1]]["sequence_coverage"] is None
    assert pg["P12345"]["sequence_coverage"] is None
    assert pg["P12345"]["molecular_weight"] is None


def test_coverage_does_not_alias_decoy_evidence_to_a_target(dataset_dir, tmp_path):
    """A target PSM can name decoy proteins; their identifiers must stay distinct."""
    from qpx.writers import PsmWriter
    from tests.conftest import make_psm_record

    psm = make_psm_record(sequence="THIRDPEP", run_file_name="run_02")
    psm["protein_accessions"] = ["DECOY_sp|P67890|PROT2_HUMAN"]
    with PsmWriter(dataset_dir / "exp.psm.parquet") as writer:
        writer.write_batch([psm])

    _run(dataset_dir, _dataset_fasta(tmp_path), "--in-place")

    pg = {row["anchor_protein"]: row for row in _rows(dataset_dir, "pg")}
    assert pg["P67890"]["sequence_coverage"] is None


def test_rewrite_keeps_the_qpx_parquet_encoding(dataset_dir, tmp_path):
    """Rewritten views must be encoded like every other QPX writer, not with pyarrow
    defaults (byte-stream-split on rt/mz, zstd level, dictionaries, format 2.6)."""

    def encodings(path):
        meta = pq.ParquetFile(path).metadata
        out = {}
        for rg in range(meta.num_row_groups):
            for col in range(meta.num_columns):
                chunk = meta.row_group(rg).column(col)
                out[chunk.path_in_schema] = (tuple(sorted(chunk.encodings)), chunk.compression)
        return out, meta.format_version

    before, before_version = encodings(dataset_dir / "exp.feature.parquet")
    _run(dataset_dir, _dataset_fasta(tmp_path), "--in-place")
    after, after_version = encodings(dataset_dir / "exp.feature.parquet")

    assert after_version == before_version
    for column in ("rt", "calculated_mz", "observed_mz"):
        if column in before:
            assert after[column] == before[column], column


def test_coverage_uses_a_partitioned_psm_view(dataset_dir, tmp_path):
    """Dataset reads a PSM view split into psm/; skipping it silently fell back to
    group membership for coverage (the ACTB 4.5% vs 84.3% undercount)."""
    from qpx.writers import PsmWriter
    from tests.conftest import make_psm_record

    psm = make_psm_record(sequence="THIRDPEP", run_file_name="run_02")
    psm["protein_accessions"] = ["P12345", "P67890"]
    part = dataset_dir / "psm" / "run_file_name=run_02"
    part.mkdir(parents=True)
    with PsmWriter(part / "part-0.parquet") as writer:
        writer.write_batch([psm])

    result = _run(dataset_dir, _dataset_fasta(tmp_path), "--in-place")

    pg = {r["anchor_protein"]: r for r in _rows(dataset_dir, "pg")}
    assert pg["P67890"]["sequence_coverage"] == pytest.approx(100 * 8 / 13, rel=1e-5)
    assert "coverage evidence from: feature, psm" in result.output


def test_reports_when_no_psm_evidence_is_available(dataset_dir, tmp_path):
    result = _run(dataset_dir, _dataset_fasta(tmp_path), "--in-place")

    assert "coverage evidence from: feature" in result.output
    assert "no PSM view" in result.output
