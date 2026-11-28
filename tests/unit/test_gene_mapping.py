"""Regression tests for identity-preserving gene annotation writes."""

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx import Dataset
from qpx.transforms.gene_mapping import (
    GeneMappingTransform,
    _annotate_table,
    _parse_fasta_header,
    _parse_gene_names_from_fasta,
)
from qpx.writers.feature import FeatureWriter
from tests.conftest import make_feature_record


def test_write_annotated_features_preserves_source_identity(tmp_path):
    """Annotation must not replace a producer-specific Feature identity recipe."""
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    source_path = source_dir / "fragpipe.feature.parquet"
    records = []
    for voltage in (-65.0, -45.0):
        record = make_feature_record(run_file_name="experiment", scan=[])
        record.update(
            {
                "rt": None,
                "quantification_unit_id": "experiment",
                "compensation_voltage": voltage,
            }
        )
        records.append(record)

    composite = ("quantification_unit_id", "peptidoform", "charge", "compensation_voltage")
    with FeatureWriter(source_path, identity_composite=composite) as writer:
        writer.write_batch(records)

    dataset = Dataset(source_dir, structures=["feature"])
    fasta = tmp_path / "empty.fasta"
    fasta.touch()
    transform = GeneMappingTransform(fasta)

    output_path = tmp_path / "annotated.feature.parquet"
    transform.write_annotated_features(dataset, output_path)

    source = pq.read_table(source_path)
    output = pq.read_table(output_path)
    assert output.column("feature_id").to_pylist() == source.column("feature_id").to_pylist()
    assert output.schema.metadata[b"identity_composite"] == b",".join(field.encode() for field in composite)


def test_parse_gene_names_reads_headers_without_biopython(tmp_path):
    """Gene names come from the FASTA headers alone — no optional dependency."""
    fasta = tmp_path / "db.fasta"
    fasta.write_text(
        ">sp|P12345|BRCA1_HUMAN Breast cancer type 1 OS=Homo sapiens GN=BRCA1 PE=1 SV=2\n"
        "MKVLAA\nGGWSTR\n"
        ">tr|Q99999|Q99999_HUMAN Uncharacterized protein OS=Homo sapiens PE=4 SV=1\n"
        "MKV\n"
        ">CONTAM_TRYP_PIG Trypsin\n"
        "MKV\n"
    )

    gene_map = _parse_gene_names_from_fasta(str(fasta))

    assert gene_map["P12345"] == {"BRCA1"}
    assert gene_map["Q99999"] == {None}
    assert gene_map["CONTAM_TRYP_PIG"] == {None}


def test_parse_gene_names_reads_gzipped_fasta(tmp_path):
    """A gzip-compressed FASTA is opened transparently."""
    import gzip

    fasta = tmp_path / "db.fasta.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">sp|P12345|PROT_HUMAN Some protein OS=Homo sapiens GN=TP53 PE=1 SV=2\nMKV\n")

    assert _parse_gene_names_from_fasta(str(fasta))["P12345"] == {"TP53"}


def test_parse_fasta_header_handles_each_identifier_shape():
    """UniProt three-field, two-field and bare identifiers all resolve."""
    assert _parse_fasta_header(">sp|P12345|PROT_HUMAN d GN=BRCA1 PE=1") == ("P12345", "PROT_HUMAN", "BRCA1")
    assert _parse_fasta_header(">sp|P12345 description") == ("P12345", "P12345", None)
    assert _parse_fasta_header(">CONTAM_ALBU Albumin GN=ALB") == ("CONTAM_ALBU", "CONTAM_ALBU", "ALB")


def test_annotate_dataframe_maps_genes_and_keeps_existing_accessions(tmp_path):
    """gg_names is filled from the FASTA; gg_accessions written by a converter survives."""
    import pandas as pd

    fasta = tmp_path / "db.fasta"
    fasta.write_text(
        ">sp|P12345|A_HUMAN a OS=Homo sapiens GN=BRCA1 PE=1 SV=2\nMKV\n"
        ">sp|P12346|B_HUMAN b OS=Homo sapiens GN=TP53 PE=1 SV=2\nMKV\n"
    )
    df = pd.DataFrame(
        {
            "pg_accessions": [["P12345", "P12346"], ["P00000"]],
            "gg_accessions": [["NC_000017.11"], None],
        }
    )

    annotated = GeneMappingTransform(fasta).annotate_dataframe(df)

    assert annotated["gg_names"].tolist() == [["BRCA1", "TP53"], None]
    assert annotated["gg_accessions"].tolist() == [["NC_000017.11"], None]


def _write_gene_bundle(directory, prefix, label="TMT126"):
    """Write a minimal pg/feature/run dataset whose genes come from the converter."""
    from qpx.writers import PgWriter, RunWriter
    from tests.conftest import make_pg_record, make_run_record

    intensities = [{"label": label, "intensity": 100.0}]
    with FeatureWriter(directory / f"{prefix}.feature.parquet") as writer:
        writer.write_batch([make_feature_record(intensities=intensities)])
    with PgWriter(directory / f"{prefix}.pg.parquet") as writer:
        writer.write_batch([make_pg_record(intensities=intensities)])

    run = make_run_record()
    run["samples"] = [
        {
            "sample_accession": f"{prefix}_{label}",
            "label": label,
            "biological_replicate": 1,
            "technical_replicate": 1,
        }
    ]
    with RunWriter(directory / f"{prefix}.run.parquet") as writer:
        writer.write_batch([run])


def test_gene_map_dataset_annotates_pg_and_refreshes_mudata(tmp_path):
    """--dataset rewrites the quantification views and rebuilds a stale h5mu."""
    import mudata as mu
    from click.testing import CliRunner

    from qpx.cli.main import qpx_main
    from qpx.mudata import write_dataset_mudata

    dataset_dir = tmp_path / "qpx_output"
    dataset_dir.mkdir()
    _write_gene_bundle(dataset_dir, "openms")

    # The converter's own view: genes as the converter wrote them.
    assert write_dataset_mudata(dataset_dir, "openms") is not None
    before = mu.read_h5mu(str(dataset_dir / "openms.h5mu"))
    assert list(before.mod["proteins"].var["gene_name"]) == ["GENE1"]
    pg_ids_before = pq.read_table(dataset_dir / "openms.pg.parquet").column("pg_id").to_pylist()

    fasta = tmp_path / "db.fasta"
    fasta.write_text(
        ">sp|P12345|A_HUMAN a OS=Homo sapiens GN=BRCA1 PE=1 SV=2\nMKV\n"
        ">sp|P12346|B_HUMAN b OS=Homo sapiens GN=TP53 PE=1 SV=2\nMKV\n"
    )

    result = CliRunner().invoke(
        qpx_main,
        ["transform", "gene-map", "--dataset", str(dataset_dir), "--fasta", str(fasta), "--in-place"],
    )

    assert result.exit_code == 0, result.output
    pg_table = pq.read_table(dataset_dir / "openms.pg.parquet")
    assert pg_table.column("gg_names").to_pylist()[0] == ["BRCA1", "TP53"]
    assert pg_table.column("pg_id").to_pylist() == pg_ids_before
    assert not (dataset_dir / ".gene_map_tmp").exists()

    after = mu.read_h5mu(str(dataset_dir / "openms.h5mu"))
    assert list(after.mod["proteins"].var["gene_name"]) == ["BRCA1"]


def test_gene_map_dataset_to_output_folder_leaves_source_untouched(tmp_path):
    """Without --in-place the source dataset is copied, not modified."""
    from click.testing import CliRunner

    from qpx.cli.main import qpx_main

    dataset_dir = tmp_path / "qpx_output"
    dataset_dir.mkdir()
    _write_gene_bundle(dataset_dir, "openms")

    fasta = tmp_path / "db.fasta"
    fasta.write_text(">sp|P12345|A_HUMAN a OS=Homo sapiens GN=BRCA1 PE=1 SV=2\nMKV\n")
    out_dir = tmp_path / "annotated"

    result = CliRunner().invoke(
        qpx_main,
        [
            "transform",
            "gene-map",
            "--dataset",
            str(dataset_dir),
            "--fasta",
            str(fasta),
            "--output-folder",
            str(out_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert pq.read_table(out_dir / "openms.pg.parquet").column("gg_names").to_pylist()[0] == ["BRCA1"]
    assert pq.read_table(dataset_dir / "openms.pg.parquet").column("gg_names").to_pylist()[0] == ["GENE1"]
    assert (out_dir / "openms.run.parquet").is_file()


def test_gene_map_dataset_copy_preserves_subdirectories(tmp_path):
    """Sharded / partitioned datasets keep views in subdirectories; the copy must keep them."""
    from click.testing import CliRunner

    from qpx.cli.main import qpx_main

    dataset_dir = tmp_path / "qpx_output"
    dataset_dir.mkdir()
    _write_gene_bundle(dataset_dir, "openms")
    shard = dataset_dir / "psm_shards"
    shard.mkdir()
    (shard / "part-0.parquet").write_bytes(b"shard-payload")

    fasta = tmp_path / "db.fasta"
    fasta.write_text(">sp|P12345|A_HUMAN a OS=Homo sapiens GN=BRCA1 PE=1 SV=2\nMKV\n")
    out_dir = tmp_path / "annotated"

    result = CliRunner().invoke(
        qpx_main,
        [
            "transform",
            "gene-map",
            "--dataset",
            str(dataset_dir),
            "--fasta",
            str(fasta),
            "--output-folder",
            str(out_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (out_dir / "psm_shards" / "part-0.parquet").read_bytes() == b"shard-payload"


def test_gene_map_parquet_adds_the_gene_columns_to_a_plain_file(tmp_path):
    """A file pandas wrote, without a QPX footer or gene columns, gains both gene columns and keeps no footer keys."""
    from click.testing import CliRunner

    from qpx.cli.main import qpx_main

    source = tmp_path / "plain.parquet"
    pd.DataFrame({"pg_accessions": [["P12345"], None]}).to_parquet(source)
    fasta = tmp_path / "db.fasta"
    fasta.write_text(">sp|P12345|A_HUMAN a OS=Homo sapiens GN=BRCA1 PE=1 SV=2\nMKV\n")
    out_dir = tmp_path / "annotated"

    result = CliRunner().invoke(
        qpx_main,
        ["transform", "gene-map", "--parquet-path", str(source), "--fasta", str(fasta), "--output-folder", str(out_dir)],
    )

    assert result.exit_code == 0, result.output
    after = pq.read_table(out_dir / source.name)
    assert {"gg_names", "gg_accessions"} <= set(after.column_names)
    assert after.column("gg_names").to_pylist() == [["BRCA1"], None]
    assert not after.schema.metadata


def test_annotate_dataframe_handles_an_empty_frame(tmp_path):
    """An empty view must not divide by zero when logging the mapped share."""
    import pandas as pd

    fasta = tmp_path / "db.fasta"
    fasta.write_text(">sp|P12345|A_HUMAN a OS=Homo sapiens GN=BRCA1 PE=1 SV=2\nMKV\n")
    empty = pd.DataFrame({"pg_accessions": []})

    annotated = GeneMappingTransform(fasta).annotate_dataframe(empty)

    assert len(annotated) == 0
    assert "gg_names" in annotated.columns


def _write_fasta(tmp_path, text):
    fasta = tmp_path / "db.fasta"
    fasta.write_text(text)
    return fasta


def test_annotate_dataframe_keeps_existing_genes_the_fasta_cannot_resolve(tmp_path):
    """A FASTA without GN= must not erase gene names a converter already wrote.

    DIA-NN conversions carry gg_names from the DIA-NN report. Running gene-map
    with a contaminants-only or GN-stripped FASTA used to overwrite every row
    with None, which --in-place made unrecoverable.
    """
    fasta = _write_fasta(
        tmp_path,
        ">sp|P12345|PROT_HUMAN desc PE=1 SV=2\nMKV\n>sp|Q99999|OTHR_HUMAN desc\nMKV\n",
    )
    df = pd.DataFrame(
        {
            "pg_accessions": [["P12345"], ["Q99999"]],
            "gg_names": [["BRCA1"], ["TP53"]],
        }
    )

    annotated = GeneMappingTransform(fasta).annotate_dataframe(df)

    assert annotated["gg_names"].tolist() == [["BRCA1"], ["TP53"]]


def test_annotate_dataframe_overwrites_existing_genes_the_fasta_does_resolve(tmp_path):
    """Preserving unresolved rows must not stop the FASTA from correcting resolved ones."""
    fasta = _write_fasta(tmp_path, ">sp|P12345|PROT_HUMAN desc GN=NEWGENE PE=1\nMKV\n")
    df = pd.DataFrame({"pg_accessions": [["P12345"]], "gg_names": [["STALE"]]})

    annotated = GeneMappingTransform(fasta).annotate_dataframe(df)

    assert annotated["gg_names"].tolist() == [["NEWGENE"]]


def test_annotate_table_merges_existing_large_string_genes():
    """Existing gg_names stored as large strings (as polars writes them) are kept where the FASTA resolves none."""
    table = pa.table(
        {
            "pg_accessions": pa.array([["P12345"], ["Q99999"]], type=pa.list_(pa.string())),
            "gg_names": pa.array([["STALE"], ["KEPT"]], type=pa.list_(pa.large_string())),
        }
    )

    annotated, mapped = _annotate_table(table, {"P12345": {"NEWGENE"}}, pa.list_(pa.string()))

    assert annotated.column("gg_names").to_pylist() == [["NEWGENE"], ["KEPT"]]
    assert mapped == 1


def test_annotate_dataframe_reports_zero_share_for_a_fasta_without_genes(tmp_path):
    fasta = _write_fasta(tmp_path, ">sp|P12345|PROT_HUMAN desc PE=1\nMKV\n")
    transform = GeneMappingTransform(fasta)

    transform.annotate_dataframe(pd.DataFrame({"pg_accessions": [["P12345"]]}))

    assert transform.last_mapped_share == 0.0


def test_parse_gene_names_log_counts_only_identifiers_carrying_a_gene(tmp_path, caplog):
    fasta = _write_fasta(
        tmp_path,
        ">sp|P12345|A_HUMAN desc GN=BRCA1\nMKV\n>sp|Q99999|B_HUMAN desc\nMKV\n",
    )
    with caplog.at_level("INFO"):
        GeneMappingTransform(fasta).gene_map

    assert "Parsed gene names for 1/2 protein identifiers" in caplog.text


def test_unrecognised_map_by_is_rejected(tmp_path):
    fasta = _write_fasta(tmp_path, ">sp|P12345|A_HUMAN desc GN=BRCA1\nMKV\n")

    with pytest.raises(ValueError, match="map_by must be one of"):
        GeneMappingTransform(fasta, map_by="Accession")
