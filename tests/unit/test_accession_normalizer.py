"""Accession normalization rewrites accession values and nothing else."""

from __future__ import annotations

import pyarrow as pa
import pyarrow.parquet as pq

from qpx.transforms.accession_normalizer import normalize_dataset_accessions, normalize_parquet, to_bare

_MEMBER = pa.struct([pa.field("accession", pa.string(), nullable=False), pa.field("start", pa.int32())])


def test_forward_maps_anchor_and_members_of_both_layouts(tmp_path):
    """Anchors and pg_accessions members are mapped; nulls, blanks and other columns are kept."""
    feature = tmp_path / "x.feature.parquet"
    pq.write_table(
        pa.table(
            {
                "anchor_protein": ["sp|P04114|A_HUMAN", None, "", "CONTAM_sp|CONTAM_P02768|CONTAM_ALBU_HUMAN"],
                "pg_accessions": pa.array(
                    [[{"accession": "sp|P04114|A_HUMAN", "start": 3}, None], None, [], [{"accession": "P01308", "start": None}]],
                    pa.list_(_MEMBER),
                ),
                "charge": [2, 3, 2, 1],
            }
        ),
        feature,
    )
    pg = tmp_path / "x.pg.parquet"
    pq.write_table(pa.table({"pg_accessions": [["tr|Q15149|PLEC_HUMAN", None, "P01308"], None]}), pg)

    result = normalize_parquet(feature, tmp_path / "out.feature.parquet", direction="forward")
    normalize_parquet(pg, tmp_path / "out.pg.parquet", direction="forward")

    table = pq.read_table(tmp_path / "out.feature.parquet")
    assert result["accessions_changed"] == 2
    assert table.column("anchor_protein").to_pylist() == ["P04114", None, "", "CONTAM_P02768"]
    assert table.column("pg_accessions").to_pylist() == [
        [{"accession": "P04114", "start": 3}, None],
        None,
        [],
        [{"accession": "P01308", "start": None}],
    ]
    assert table.column("charge").to_pylist() == [2, 3, 2, 1]
    assert pq.read_table(tmp_path / "out.pg.parquet").column("pg_accessions").to_pylist() == [["Q15149", None, "P01308"], None]


def test_file_without_changes_is_copied(tmp_path):
    """A view with nothing to normalize is copied byte for byte."""
    source = tmp_path / "x.pg.parquet"
    pq.write_table(pa.table({"anchor_protein": ["P04114"], "pg_accessions": [["P04114", "P02768"]]}), source)

    result = normalize_parquet(source, tmp_path / "out.pg.parquet", direction="forward")

    assert result["accessions_changed"] == 0
    assert (tmp_path / "out.pg.parquet").read_bytes() == source.read_bytes()


def test_dataset_normalization_covers_psm_lists_and_partitioned_views(tmp_path):
    """Every view is stripped in place, PSM protein lists and partitioned parts included."""
    pq.write_table(
        pa.table(
            {
                "anchor_protein": ["sp|CONTAM_P05787|CONTAM_K2C8_HUMAN"],
                "pg_accessions": [["sp|CONTAM_P05787|CONTAM_K2C8_HUMAN", "Q15149"]],
                "pg_names": [["CONTAM_K2C8_HUMAN", "PLEC_HUMAN"]],
            }
        ),
        tmp_path / "openms.pg.parquet",
    )
    part = tmp_path / "psm" / "run_file_name=a" / "part-0.parquet"
    part.parent.mkdir(parents=True)
    pq.write_table(pa.table({"protein_accessions": [["sp|Q15149|PLEC_HUMAN", "tr|A0A0C4DGB6|A0A0C4DGB6_HUMAN"], None]}), part)

    assert normalize_dataset_accessions(tmp_path, "openms") == 2

    pg = pq.read_table(tmp_path / "openms.pg.parquet")
    assert pg.column("anchor_protein").to_pylist() == ["CONTAM_P05787"]
    assert pg.column("pg_accessions").to_pylist() == [["CONTAM_P05787", "Q15149"]]
    assert pg.column("pg_names").to_pylist() == [["CONTAM_K2C8_HUMAN", "PLEC_HUMAN"]]
    assert pq.read_table(part).column("protein_accessions").to_pylist() == [["Q15149", "A0A0C4DGB6"], None]


def test_decoy_header_keeps_its_prefix_on_the_bare_accession():
    assert to_bare("DECOY_sp|P04114|APOB_HUMAN") == "DECOY_P04114"
    assert to_bare("DECOY_tr|A0A0C4DGB6|A0A0C4DGB6_HUMAN") == "DECOY_A0A0C4DGB6"
    assert to_bare("DECOY_P04114") == "DECOY_P04114"


def test_psm_list_rewrites_are_counted(tmp_path):
    psm = tmp_path / "x.psm.parquet"
    pq.write_table(pa.table({"protein_accessions": [["sp|P04114|A_HUMAN", "P02768"], ["sp|P04114|A_HUMAN"]]}), psm)
    result = normalize_parquet(psm, psm)
    assert result["accessions_changed"] == 2 and result["distinct_changed"] == 1


def test_pepmap_members_are_normalized_with_the_other_views(tmp_path):
    pq.write_table(
        pa.table({"pg_accessions": pa.array([[{"accession": "sp|P04114|A_HUMAN", "start": 3}]], pa.list_(_MEMBER))}),
        tmp_path / "x.pepmap.parquet",
    )
    assert normalize_dataset_accessions(tmp_path, "x") == 1
    assert pq.read_table(tmp_path / "x.pepmap.parquet").column("pg_accessions").to_pylist() == [
        [{"accession": "P04114", "start": 3}]
    ]


def test_sequence_edit_protein_is_normalized_but_not_its_variant_accession(tmp_path):
    edit = pa.struct(
        [pa.field("protein_accession", pa.string()), pa.field("start", pa.int32()), pa.field("accession", pa.string())]
    )
    pepmap = tmp_path / "x.pepmap.parquet"
    pq.write_table(
        pa.table(
            {
                "sequence_edits": pa.array(
                    [[{"protein_accession": "sp|P04637|P53_HUMAN", "start": 72, "accession": "VAR_005851"}], None], pa.list_(edit)
                )
            }
        ),
        pepmap,
    )
    normalize_parquet(pepmap, pepmap)
    assert pq.read_table(pepmap).column("sequence_edits").to_pylist() == [
        [{"protein_accession": "P04637", "start": 72, "accession": "VAR_005851"}],
        None,
    ]


def test_only_uniprot_headers_are_stripped():
    """Other |-separated identifiers, and sp|/tr| headers without a UniProt accession, are kept as written."""
    for kept in ("gnl|mydb|PROT1", "ref|NP_000537.3|", "sp|MYPROT1|CUSTOM", "tr|contig_12|orf3", "DECOY_sp|abc|x"):
        assert to_bare(kept) == kept
    assert to_bare("sp|P04637-2|P53_HUMAN") == "P04637-2"


def test_cli_rebuilds_the_mudata_view_from_the_normalized_views(tmp_path):
    """normalize-accessions leaves no stale h5mu: its protein axis carries the bare accessions."""
    import mudata as mu
    from click.testing import CliRunner

    from qpx.cli.main import qpx_main
    from qpx.mudata import write_dataset_mudata
    from qpx.writers import FeatureWriter, PgWriter, RunWriter
    from tests.conftest import make_feature_record, make_pg_record, make_run_record

    full = "sp|P04114|APOB_HUMAN"
    intensities = [{"label": "TMT126", "intensity": 100.0}]
    feature = make_feature_record(anchor_protein=full, intensities=intensities)
    feature["pg_accessions"] = [{"accession": full, "start": None, "end": None}]
    with FeatureWriter(tmp_path / "openms.feature.parquet") as writer:
        writer.write_batch([feature])
    with PgWriter(tmp_path / "openms.pg.parquet") as writer:
        writer.write_batch([make_pg_record(anchor_protein=full, pg_accessions=[full], intensities=intensities)])
    run = make_run_record()
    run["samples"] = [{"sample_accession": "S1", "label": "TMT126", "biological_replicate": 1, "technical_replicate": 1}]
    with RunWriter(tmp_path / "openms.run.parquet") as writer:
        writer.write_batch([run])
    assert write_dataset_mudata(tmp_path, "openms") is not None
    assert list(mu.read_h5mu(str(tmp_path / "openms.h5mu")).mod["proteins"].var_names) == [full]

    result = CliRunner().invoke(qpx_main, ["transform", "normalize-accessions", "--dataset", str(tmp_path), "--in-place"])

    assert result.exit_code == 0, result.output
    assert list(mu.read_h5mu(str(tmp_path / "openms.h5mu")).mod["proteins"].var_names) == ["P04114"]


def test_cli_refreshes_an_existing_integrity_record(tmp_path):
    """Recorded checksums are recomputed for the rewritten views, so the record still verifies."""
    import hashlib

    from click.testing import CliRunner

    from qpx.cli.convert.postprocess import _maybe_compute_integrity
    from qpx.cli.main import qpx_main
    from qpx.dataset import Dataset
    from qpx.writers import DatasetWriter, PgWriter
    from tests.conftest import make_dataset_record, make_pg_record

    full = "sp|P04114|APOB_HUMAN"
    with PgWriter(tmp_path / "t.pg.parquet") as writer:
        writer.write_batch([make_pg_record(anchor_protein=full, pg_accessions=[full])])
    with DatasetWriter(tmp_path / "t.dataset.parquet") as writer:
        writer.write_batch([make_dataset_record()])
    _maybe_compute_integrity(tmp_path, "t", True)

    result = CliRunner().invoke(qpx_main, ["transform", "normalize-accessions", "--dataset", str(tmp_path), "--in-place"])

    assert result.exit_code == 0, result.output
    with Dataset(str(tmp_path), file_prefix="t") as ds:
        record = ds.dataset_meta.to_df().iloc[0].to_dict()
    assert dict(record["file_checksums"])["t.pg.parquet"] == hashlib.sha256((tmp_path / "t.pg.parquet").read_bytes()).hexdigest()
