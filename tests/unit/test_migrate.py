"""qpxc migrate brings a pre-1.1 dataset to the current specification."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.core.data.identity import derive_composite_id
from qpx.transforms.migrate import migrate_dataset
from qpx.version import QPX_SPEC_VERSION
from tests.conftest import make_xic_record

RUN = "run_01"
PROTEIN = "P12345"


def _write(path, rows, creator=b"diann", extra_metadata=None):
    table = pa.Table.from_pylist(rows)
    metadata = {b"qpx_version": b"1.0", b"creator": creator, **(extra_metadata or {})}
    pq.write_table(table.replace_schema_metadata(metadata), path)


@pytest.fixture
def old_dataset(tmp_path):
    """The layout QPX 1.0.2 wrote: no ids, list intensities on pg, 'raw' labels, 0.0 m/z."""
    folder = tmp_path / "old"
    folder.mkdir()
    feature = {
        "sequence": "PEPTMIDEK",
        "peptidoform": "PEPTM(UniMod:35)IDEK",
        "charge": 2,
        "is_decoy": False,
        "calculated_mz": 0.0,
        "observed_mz": 0.0,
        "run_file_name": RUN,
        "scan": [100],
        "anchor_protein": PROTEIN,
        "pg_accessions": [{"accession": PROTEIN, "start": None, "end": None}],
        "intensities": [{"label": "raw", "intensity": 10.0}],
    }
    _write(folder / "old.feature.parquet", [feature])
    pg = {
        "pg_accessions": [PROTEIN],
        "anchor_protein": PROTEIN,
        "run_file_name": RUN,
        "is_decoy": False,
        "peptides": [{"protein_name": PROTEIN, "peptide_count": 1}],
        "intensities": [{"label": "TMT126", "intensity": 1.0}, {"label": "TMT127N", "intensity": 2.0}],
        "additional_intensities": [
            {"label": "TMT126", "intensities": [{"intensity_name": "lfq", "intensity_value": 3.0}]},
            {"label": "TMT127N", "intensities": [{"intensity_name": "lfq", "intensity_value": 4.0}]},
        ],
    }
    _write(
        folder / "old.pg.parquet",
        [pg, {**pg, "pg_accessions": ["Q1"], "anchor_protein": "Q1", "intensities": None, "additional_intensities": None}],
    )
    _write(
        folder / "old.run.parquet",
        [{"run_accession": RUN, "run_file_name": RUN, "samples": [{"sample_accession": "S1", "label": "126"}]}],
    )
    (folder / "old.h5mu").write_bytes(b"not read")
    (folder / "old.sdrf.tsv").write_text("source name\n", encoding="utf-8")
    return folder


def test_a_pre_1_1_pg_is_flattened_and_gets_its_ids(old_dataset, tmp_path):
    out = tmp_path / "new"
    assert migrate_dataset(old_dataset, out)["pg"] == 3

    rows = pq.read_table(out / "old.pg.parquet").to_pylist()
    assert [(row["anchor_protein"], row["label"], row["intensity"]) for row in rows] == [
        (PROTEIN, "TMT126", 1.0),
        (PROTEIN, "TMT127N", 2.0),
        ("Q1", None, None),
    ]
    assert all(row["grouped_runs"] == [RUN] for row in rows)
    assert rows[1]["additional_intensities"] == [
        {"label": "TMT127N", "intensities": [{"intensity_name": "lfq", "intensity_value": 4.0}]}
    ]
    assert rows[0]["pg_id"] == derive_composite_id(rows[0], ("pg_accessions", "grouped_runs", "label"))
    assert pq.read_schema(out / "old.pg.parquet").metadata[b"qpx_version"] == QPX_SPEC_VERSION.encode()


def test_a_feature_takes_the_producers_composite_and_canonical_values(old_dataset, tmp_path):
    out = tmp_path / "new"
    migrate_dataset(old_dataset, out)

    row = pq.read_table(out / "old.feature.parquet").to_pylist()[0]
    assert row["peptidoform"] == "PEPTM[UNIMOD:35]IDEK"
    assert row["intensities"] == [{"label": "LFQ", "intensity": 10.0}]
    assert row["observed_mz"] is None
    assert row["calculated_mz"] is None
    footer = pq.read_schema(out / "old.feature.parquet").metadata
    assert footer[b"identity_composite"] == b"peptidoform,charge,run_file_name"
    assert row["feature_id"] == derive_composite_id(row, ("peptidoform", "charge", "run_file_name"))


def test_labels_of_runs_become_canonical_and_only_views_are_rewritten(old_dataset, tmp_path):
    out = tmp_path / "new"
    migrate_dataset(old_dataset, out)

    assert pq.read_table(out / "old.run.parquet").to_pylist()[0]["samples"][0]["label"] == "TMT126"
    assert (out / "old.sdrf.tsv").exists()
    assert not (out / "old.h5mu").exists()


def test_a_psm_follows_the_new_feature_id(old_dataset, tmp_path):
    feature = pq.read_table(old_dataset / "old.feature.parquet").to_pylist()[0]
    _write(old_dataset / "old.feature.parquet", [{**feature, "feature_id": 7}])
    psm = {key: feature[key] for key in ("sequence", "peptidoform", "charge", "is_decoy", "run_file_name", "scan")}
    _write(old_dataset / "old.psm.parquet", [{**psm, "psm_id": 1, "feature_id": 7}])
    out = tmp_path / "new"
    migrate_dataset(old_dataset, out)

    new_id = pq.read_table(out / "old.feature.parquet").to_pylist()[0]["feature_id"]
    row = pq.read_table(out / "old.psm.parquet").to_pylist()[0]
    assert new_id != 7
    assert row["feature_id"] == new_id
    assert row["psm_id"] != 1


def test_a_trace_follows_the_new_feature_id(old_dataset, tmp_path):
    """xic.feature_ids and feature_scores[].feature_id take the new feature_id."""
    feature = pq.read_table(old_dataset / "old.feature.parquet").to_pylist()[0]
    _write(old_dataset / "old.feature.parquet", [{**feature, "feature_id": 7}])
    _write(old_dataset / "old.xic.parquet", [make_xic_record(7, peptidoform=feature["peptidoform"], run_file_name=RUN)])
    out = tmp_path / "new"
    migrate_dataset(old_dataset, out)

    new_id = pq.read_table(out / "old.feature.parquet").to_pylist()[0]["feature_id"]
    row = pq.read_table(out / "old.xic.parquet").to_pylist()[0]
    assert new_id != 7
    assert row["feature_ids"] == [new_id]
    assert [score["feature_id"] for score in row["feature_scores"]] == [new_id]


def test_a_trace_score_keeps_the_id_of_a_feature_the_view_lacks(old_dataset, tmp_path):
    """feature_scores[].feature_id cannot be null, so it keeps an id that no feature has; feature_ids drops it."""
    _write(old_dataset / "old.xic.parquet", [make_xic_record(7, run_file_name=RUN)])
    out = tmp_path / "new"
    migrate_dataset(old_dataset, out)

    row = pq.read_table(out / "old.xic.parquet").to_pylist()[0]
    assert row["feature_ids"] == [None]
    assert [score["feature_id"] for score in row["feature_scores"]] == [7]


def test_a_site_follows_the_new_pg_and_psm_ids(old_dataset, tmp_path):
    """relayout numbers ids 1, 2, ...; migrate derives them again, and the site's links follow."""
    feature = pq.read_table(old_dataset / "old.feature.parquet").to_pylist()[0]
    psm = {key: feature[key] for key in ("sequence", "peptidoform", "charge", "is_decoy", "run_file_name", "scan")}
    _write(old_dataset / "old.psm.parquet", [{**psm, "psm_id": 1}])
    pg = {
        "pg_id": 1,
        "pg_accessions": [PROTEIN],
        "anchor_protein": PROTEIN,
        "grouped_runs": [RUN],
        "is_decoy": False,
        "peptides": [{"protein_name": PROTEIN, "peptide_count": 1}],
    }
    _write(old_dataset / "old.pg.parquet", [pg])
    site = {
        "protein_accession": PROTEIN,
        "position": 5,
        "amino_acid": "M",
        "modification_name": "Oxidation",
        "modification_accession": "UNIMOD:35",
        "mass_delta": 15.994915,
        "grouped_runs": [RUN],
        "pg_id": 1,
        "best_psm_id": 1,
        "is_decoy": False,
    }
    _write(old_dataset / "old.site.parquet", [site])
    out = tmp_path / "new"
    migrate_dataset(old_dataset, out)

    pg_id = pq.read_table(out / "old.pg.parquet").column("pg_id")[0].as_py()
    psm_id = pq.read_table(out / "old.psm.parquet").column("psm_id")[0].as_py()
    row = pq.read_table(out / "old.site.parquet").to_pylist()[0]
    assert 1 not in (pg_id, psm_id)
    assert (row["pg_id"], row["best_psm_id"]) == (pg_id, psm_id)


def test_a_provided_psm_id_is_kept(old_dataset, tmp_path):
    feature = pq.read_table(old_dataset / "old.feature.parquet").to_pylist()[0]
    psm = {key: feature[key] for key in ("sequence", "peptidoform", "charge", "is_decoy", "run_file_name", "scan")}
    _write(old_dataset / "old.psm.parquet", [{**psm, "psm_id": 11}], extra_metadata={b"identity_source": b"provided"})
    out = tmp_path / "new"
    migrate_dataset(old_dataset, out)

    assert pq.read_table(out / "old.psm.parquet").to_pylist()[0]["psm_id"] == 11


def test_the_command_refuses_an_existing_folder(old_dataset, tmp_path):
    out = tmp_path / "exists"
    out.mkdir()
    result = CliRunner().invoke(qpx_main, ["migrate", "--dataset-path", str(old_dataset), "--output-folder", str(out)])
    assert result.exit_code != 0
    assert "exists" in result.output
