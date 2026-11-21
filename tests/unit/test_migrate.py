"""qpxc migrate brings a pre-1.1 dataset to the current specification."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.converters.diann.feature_adapter import _FEATURE_IDENTITY_COMPOSITE
from qpx.core.data.identity import derive_composite_id
from qpx.transforms.migrate import migrate_dataset
from qpx.version import QPX_SPEC_VERSION
from qpx.writers import ProteoformWriter, PsmWriter, SiteWriter
from tests.conftest import make_feature_record, make_pg_record, make_psm_record, make_sample_record, make_xic_record
from tests.converters.test_diann_converter import _feature_rows_single_run
from tests.unit.test_site_proteoform import make_proteoform_record, make_site_record

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
    assert footer[b"identity_composite"] == ",".join(_FEATURE_IDENTITY_COMPOSITE).encode()
    assert row["is_library_decoy"] is False
    assert row["feature_id"] == derive_composite_id(row, _FEATURE_IDENTITY_COMPOSITE)


def _diann_migration_features(tmp_path):
    """Current converter output and the caret-bearing peptidoforms older releases wrote."""
    current = tmp_path / "current"
    current.mkdir()
    common = {"Run": RUN, "Protein.Group": PROTEIN, "Precursor.Charge": 2, "Q.Value": 0.001, "Precursor.Quantity": 100.0}
    inputs = [("PEPTIDEK", "PEPTIDEK", 0), ("^PEPTIDEK", "PEPTIDEK", 1), ("^(UniMod:1)MEPSM(UniMod:35)K", "MEPSMK", 1)]
    report = [
        {**common, "Modified.Sequence": modified, "Stripped.Sequence": sequence, "Precursor.Id": f"{modified}2", "Decoy": decoy}
        for modified, sequence, decoy in inputs
    ]
    expected = _feature_rows_single_run(current, report)
    old_peptidoforms = {
        ("PEPTIDEK", False): "PEPTIDEK",
        ("PEPTIDEK", True): "^PEPTIDEK",
        ("MEPSMK", True): "^[UNIMOD:1]MEPSM[UNIMOD:35]K",
    }
    fields = ("sequence", "charge", "run_file_name", "is_decoy", "scan")
    legacy = [
        {**{name: row[name] for name in fields}, "peptidoform": old_peptidoforms[row["sequence"], row["is_decoy"]]}
        for row in expected
    ]
    return expected, legacy


def _write_legacy_feature_references(folder, feature):
    """Write PSM, fragment and trace references to one old DIA-NN library decoy."""
    psm = {**make_psm_record(is_decoy=True), "psm_id": 1, "feature_id": feature["feature_id"]}
    _write(folder / "old.psm.parquet", [psm])
    fragment = {
        "feature_id": feature["feature_id"],
        "run_file_name": RUN,
        "peptidoform": "PEPTIDEK",
        "charge": 2,
        "fragment_index": 1,
    }
    _write(folder / "old.fragment.parquet", [fragment])
    _write(folder / "old.xic.parquet", [make_xic_record(feature["feature_id"])])


@pytest.mark.parametrize("composite", [None, b"peptidoform,charge,run_file_name"], ids=["qpx-1.0.x", "qpx-1.3.0"])
def test_diann_migration_matches_a_new_conversion(tmp_path, composite):
    """Both legacy layouts recover the current target and library-decoy identities and links."""
    expected, legacy = _diann_migration_features(tmp_path)
    source = tmp_path / "old"
    source.mkdir()
    metadata = {b"qpx_version": b"1.0.2" if composite is None else b"1.2"}
    if composite is not None:
        metadata[b"identity_composite"] = composite
        for row in legacy:
            row["feature_id"] = derive_composite_id(row, tuple(composite.decode().split(",")))
        _write_legacy_feature_references(source, next(row for row in legacy if row["peptidoform"] == "^PEPTIDEK"))
    _write(source / "old.feature.parquet", legacy, extra_metadata=metadata)
    output = tmp_path / "migrated"
    migrate_dataset(source, output)

    columns = ["feature_id", "peptidoform", "is_decoy", "is_library_decoy"]
    migrated = pq.read_table(output / "old.feature.parquet", columns=columns).to_pylist()
    assert migrated == [{name: row[name] for name in columns} for row in expected]
    footer = pq.read_schema(output / "old.feature.parquet").metadata
    assert footer[b"identity_composite"] == ",".join(_FEATURE_IDENTITY_COMPOSITE).encode()
    if composite is not None:
        library_id = next(row["feature_id"] for row in migrated if row["peptidoform"] == "PEPTIDEK" and row["is_library_decoy"])
        for view in ("psm", "fragment"):
            assert pq.read_table(output / f"old.{view}.parquet").column("feature_id").to_pylist() == [library_id]
        trace = pq.read_table(output / "old.xic.parquet").to_pylist()[0]
        assert trace["feature_ids"] == [library_id]
        assert trace["feature_scores"][0]["feature_id"] == library_id
    repeated = tmp_path / "repeated"
    migrate_dataset(output, repeated)
    assert pq.read_table(repeated / "old.feature.parquet", columns=columns).to_pylist() == migrated


def test_a_non_diann_feature_keeps_its_declared_composite(tmp_path):
    """A native OpenMS producer keeps its own declared composite during migration."""
    source = tmp_path / "old"
    source.mkdir()
    composite = ("peptidoform", "charge", "run_file_name", "scan")
    _write(
        source / "old.feature.parquet",
        [make_feature_record()],
        creator=b"openms",
        extra_metadata={b"identity_composite": ",".join(composite).encode()},
    )
    output = tmp_path / "migrated"
    migrate_dataset(source, output)
    table = pq.read_table(output / "old.feature.parquet")
    assert table.schema.metadata[b"identity_composite"] == ",".join(composite).encode()
    row = table.to_pylist()[0]
    assert row["feature_id"] == derive_composite_id(row, composite)


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


def test_traces_keep_their_selection(old_dataset, tmp_path):
    """The footer keys of the trace selection are copied with the traces."""
    selection = {b"xic_selection": b"all", b"xic_window": b"full", b"transition_format": b"mzPAF"}
    _write(old_dataset / "old.xic.parquet", [make_xic_record(7, run_file_name=RUN)], extra_metadata=selection)
    out = tmp_path / "new"
    migrate_dataset(old_dataset, out)

    footer = pq.read_schema(out / "old.xic.parquet").metadata
    assert {key: footer.get(key) for key in selection} == selection


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


def test_a_mzidentml_psm_id_is_kept(old_dataset, tmp_path):
    """mzIdentML PSM ids cannot be reproduced from the persisted composite columns."""
    feature = pq.read_table(old_dataset / "old.feature.parquet").to_pylist()[0]
    psm = {key: feature[key] for key in ("sequence", "peptidoform", "charge", "is_decoy", "run_file_name", "scan")}
    _write(
        old_dataset / "old.psm.parquet",
        [{**psm, "psm_id": 11}],
        creator=b"mzidentml",
        extra_metadata={b"identity_source": b"provided"},
    )
    out = tmp_path / "new"
    migrate_dataset(old_dataset, out)

    assert pq.read_table(out / "old.psm.parquet").to_pylist()[0]["psm_id"] == 11


def test_relayouted_ids_are_derived_again_and_references_follow(dataset_dir, tmp_path):
    """Migration re-hashes ordinal PSM ids and follows feature, site and proteoform links."""
    feature = pq.read_table(dataset_dir / "exp.feature.parquet").to_pylist()[0]
    pg_id = pq.read_table(dataset_dir / "exp.pg.parquet").column("pg_id")[0].as_py()
    with PsmWriter(dataset_dir / "exp.psm.parquet", creator="diann") as writer:
        writer.write_batch([{**make_psm_record(), "feature_id": feature["feature_id"]}])
    psm_id = pq.read_table(dataset_dir / "exp.psm.parquet").column("psm_id")[0].as_py()
    for view, writer_class, record in (
        ("site", SiteWriter, make_site_record(protein_accession=PROTEIN, family=PROTEIN)),
        ("proteoform", ProteoformWriter, make_proteoform_record(family=PROTEIN, accession=f"{PROTEIN}-1", kind="isoform")),
    ):
        with writer_class(dataset_dir / f"exp.{view}.parquet") as writer:
            writer.write_batch([{**record, "pg_id": pg_id, "best_psm_id": psm_id}])

    folder = tmp_path / "sorted"
    result = CliRunner().invoke(
        qpx_main,
        ["transform", "relayout", "--dataset", str(dataset_dir), "--output-folder", str(folder)],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert pq.read_table(folder / "exp.psm.parquet").column("psm_id").to_pylist() == [1]
    out = tmp_path / "new"
    migrate_dataset(folder, out)

    for view, id_column in (("feature", "feature_id"), ("pg", "pg_id"), ("psm", "psm_id")):
        path = out / f"exp.{view}.parquet"
        table = pq.read_table(path)
        footer = pq.read_metadata(path).metadata
        assert footer[b"identity_source"] == b"derived"
        composite = tuple(footer[b"identity_composite"].decode().split(","))
        assert table.column(id_column).to_pylist() == [derive_composite_id(row, composite) for row in table.to_pylist()]
    psm = pq.read_table(out / "exp.psm.parquet").to_pylist()[0]
    features = pq.read_table(out / "exp.feature.parquet").to_pylist()
    assert psm["feature_id"] == next(row["feature_id"] for row in features if row["sequence"] == psm["sequence"])
    pg_id = pq.read_table(out / "exp.pg.parquet").column("pg_id")[0].as_py()
    for view in ("site", "proteoform"):
        row = pq.read_table(out / f"exp.{view}.parquet").to_pylist()[0]
        assert (row["pg_id"], row["best_psm_id"]) == (pg_id, psm["psm_id"])


def test_the_command_refuses_an_existing_folder(old_dataset, tmp_path):
    out = tmp_path / "exists"
    out.mkdir()
    result = CliRunner().invoke(qpx_main, ["migrate", "--dataset-path", str(old_dataset), "--output-folder", str(out)])
    assert result.exit_code != 0
    assert "exists" in result.output


def _lfq_feature(**values):
    return {**make_feature_record(intensities=[{"label": "LFQ", "intensity": 10.0}]), **values}


def _composite(path):
    return tuple(pq.read_schema(path).metadata[b"identity_composite"].decode().split(","))


def test_accessions_and_modifications_take_the_converters_spelling(tmp_path):
    """UniProt headers become bare accessions before ids are derived; a modification gets its title and accession,
    and the sequence is read from the peptidoform."""
    source = tmp_path / "old"
    source.mkdir()
    header = f"sp|{PROTEIN}|TEST_HUMAN"
    modifications = [
        {"name": "Oxidation", "accession": None, "positions": [{"position": 5, "amino_acid": "M", "scores": None}]},
        {"name": "UNIMOD:4", "accession": "UNIMOD:4", "positions": [{"position": 7, "amino_acid": "C", "scores": None}]},
    ]
    feature = _lfq_feature(
        sequence="PEPTMOXIDATIONICCARBAMIDOMETHYLEK",
        missed_cleavages=1,
        peptidoform="PEPTM[Oxidation]IC[Carbamidomethyl]EK",
        modifications=modifications,
        anchor_protein=header,
        pg_accessions=[{"accession": header, "start": None, "end": None}],
    )
    _write(source / "old.feature.parquet", [feature], creator=b"quantms-msstats")
    pg = make_pg_record(anchor_protein=header, pg_accessions=[header, "REV__sp|Q00402|NUM1_YEAST"])
    _write(source / "old.pg.parquet", [pg], creator=b"quantms")
    out = tmp_path / "new"
    migrate_dataset(source, out)

    row = pq.read_table(out / "old.feature.parquet").to_pylist()[0]
    assert row["anchor_protein"] == PROTEIN
    assert (row["sequence"], row["peptidoform"], row["missed_cleavages"]) == ("PEPTMICEK", "PEPTM[UNIMOD:35]IC[UNIMOD:4]EK", None)
    assert [entry["accession"] for entry in row["pg_accessions"]] == [PROTEIN]
    assert [(entry["name"], entry["accession"]) for entry in row["modifications"]] == [
        ("Oxidation", "UNIMOD:35"),
        ("Carbamidomethyl", "UNIMOD:4"),
    ]
    group = pq.read_table(out / "old.pg.parquet").to_pylist()[0]
    assert (group["anchor_protein"], group["pg_accessions"]) == (PROTEIN, [PROTEIN, "REV__Q00402"])
    assert group["pg_names"] == ["TEST_HUMAN", "NUM1_YEAST"]
    assert group["pg_id"] == derive_composite_id(group, _composite(out / "old.pg.parquet"))


@pytest.mark.parametrize(("spec", "measured"), [(b"1.2", False), (b"1.3", True)])
def test_a_diann_feature_takes_the_current_intensity_names_and_measured_mz(tmp_path, spec, measured):
    """Before 1.3 a DIA-NN observed_mz was the theoretical Precursor.Mz; 1.3 renamed two intensities."""
    source = tmp_path / "old"
    source.mkdir()
    names = [("maxlfq", 3.0), ("precursor_normalised", 2.0), ("ms1_area", 1.0)]
    feature = _lfq_feature(
        mass_error_ppm=0.01,
        additional_intensities=[
            {"label": "LFQ", "intensities": [{"intensity_name": name, "intensity_value": value} for name, value in names]}
        ],
    )
    _write(source / "old.feature.parquet", [feature], extra_metadata={b"qpx_version": spec})
    out = tmp_path / "new"
    migrate_dataset(source, out)

    row = pq.read_table(out / "old.feature.parquet").to_pylist()[0]
    entries = row["additional_intensities"][0]["intensities"]
    assert [entry["intensity_name"] for entry in entries] == ["pg_maxlfq", "normalize_intensity", "ms1_area"]
    expected = (pytest.approx(450.26), pytest.approx(0.01)) if measured else (None, None)
    assert (row["observed_mz"], row["mass_error_ppm"]) == expected


@pytest.mark.parametrize(("group", "second", "kept"), [(PROTEIN, 0.02, 0.02), ("Q00402", None, 0.01)], ids=["group", "no-group"])
def test_rows_that_repeat_an_id_keep_one_copy(tmp_path, group, second, kept):
    """1.0.2 repeated a quantms feature per q-value of its first protein; the copy holding its group's q-value stays.

    Without the group in the pg view, the first copy stays.
    """
    source = tmp_path / "old"
    source.mkdir()
    feature = _lfq_feature(pg_accessions=[{"accession": PROTEIN, "start": None, "end": None}])
    features = [{**feature, "pg_global_qvalue": 0.01}, {**feature, "pg_global_qvalue": second}]
    _write(source / "old.feature.parquet", features, creator=b"quantms-msstats")
    _write(source / "old.pg.parquet", [{**make_pg_record(pg_accessions=[group]), "global_qvalue": 0.02}], creator=b"quantms")
    psm = make_psm_record()
    _write(source / "old.psm.parquet", [psm, {**psm, "posterior_error_probability": 0.5}], creator=b"quantms")
    out = tmp_path / "new"
    written = migrate_dataset(source, out)

    assert (written["feature"], written["psm"]) == (1, 1)
    assert pq.read_table(out / "old.feature.parquet").column("pg_global_qvalue").to_pylist() == [kept]
    assert pq.read_table(out / "old.psm.parquet").column("posterior_error_probability").to_pylist() == [0.001]
    assert pq.read_metadata(out / "old.feature.parquet").metadata[b"identity_source"] == b"derived"


def test_a_feature_holding_a_label_twice_is_dropped_with_its_links(tmp_path):
    """Two features merged into one row cannot be told apart; a PSM of it loses its feature_id."""
    source = tmp_path / "old"
    source.mkdir()
    merged = make_feature_record(intensities=[{"label": "TMT126", "intensity": 1.0}, {"label": "TMT126", "intensity": 2.0}])
    features = [
        {**merged, "sequence": "MERGEDK", "peptidoform": "MERGEDK", "feature_id": 7},
        {**make_feature_record(), "feature_id": 8},
    ]
    _write(source / "old.feature.parquet", features, creator=b"quantms-msstats")
    psm = make_psm_record()
    _write(
        source / "old.psm.parquet",
        [{**psm, "psm_id": 1, "feature_id": 7}, {**psm, "scan": [1002], "psm_id": 2, "feature_id": 8}],
        creator=b"quantms",
    )
    out = tmp_path / "new"
    assert migrate_dataset(source, out)["feature"] == 1

    kept = pq.read_table(out / "old.feature.parquet").to_pylist()
    assert [row["sequence"] for row in kept] == ["PEPTIDEK"]
    assert pq.read_table(out / "old.psm.parquet").column("feature_id").to_pylist() == [None, kept[0]["feature_id"]]


def test_sdrf_characteristics_and_msnet_consensus_support_are_kept(old_dataset, tmp_path):
    """A sample keeps its columns outside the schema; MSNet's consensus_support becomes a score."""
    _write(old_dataset / "old.sample.parquet", [{**make_sample_record(), "treatment": "drug"}])
    feature = pq.read_table(old_dataset / "old.feature.parquet").to_pylist()[0]
    psm = {key: feature[key] for key in ("sequence", "peptidoform", "charge", "is_decoy", "run_file_name", "scan")}
    _write(old_dataset / "old.psm.parquet", [{**psm, "consensus_support": 0.75}], creator=b"msnet")
    out = tmp_path / "new"
    migrate_dataset(old_dataset, out)

    assert pq.read_table(out / "old.sample.parquet").to_pylist()[0]["treatment"] == "drug"
    scores = pq.read_table(out / "old.psm.parquet").to_pylist()[0]["additional_scores"]
    assert scores == [{"score_name": "consensus_support", "score_value": 0.75, "higher_better": None}]


def test_an_msnet_psm_convert_msnet_never_read_takes_its_columns(tmp_path):
    """The MSNet layout scripts/convert_msnet.py reads, without a footer, becomes what it writes."""
    folder = tmp_path / "msnet"
    folder.mkdir()
    raw = {
        "sequence": "LGMFNIQHCK",
        "modifications": [{"name": "Carbamidomethyl (C)", "positions": [{"position": "9"}]}],
        "precursor_charge": 2,
        "exp_mass_to_charge": 600.3,
        "cal_mass_to_charge": 600.0,
        "peptidoform": "LGMFNIQHC(Carbamidomethyl)K",
        "posterior_error_probability": 0.01,
        "global_qvalue": 0.25,
        "consensus_support": 0.5,
        "retention_time": 1200.5,
        "cv_params": {"Instrument": "Q Exactive", "Collision Energy": None, "isotope_error": 0},
        "reference_file_name": RUN,
        "scan": 4021,
        "is_decoy": False,
    }
    pq.write_table(pa.Table.from_pylist([raw]), folder / "msnet.psm.parquet")
    migrate_dataset(folder, tmp_path / "new")

    row = pq.read_table(tmp_path / "new" / "msnet.psm.parquet").to_pylist()[0]
    assert row["psm_id"] is not None
    assert (row["peptidoform"], row["charge"], row["run_file_name"], row["scan"], row["rt"]) == (
        "LGMFNIQHC[UNIMOD:4]K",
        2,
        RUN,
        [4021],
        1200.5,
    )
    # float32(600.3) is 600.2999878, so the ppm convert_msnet writes is 499.97965, not 500
    assert (row["observed_mz"], row["calculated_mz"], row["mass_error_ppm"]) == pytest.approx((600.3, 600.0, 499.97965494791667))
    assert row["modifications"] == [
        {
            "name": "Carbamidomethyl (C)",
            "accession": "UNIMOD:4",
            "positions": [{"position": 9, "amino_acid": None, "scores": None}],
        }
    ]
    assert row["cv_params"] == [
        {"cv_name": "Instrument", "cv_value": "Q Exactive"},
        {"cv_name": "isotope_error", "cv_value": "0"},
    ]
    assert row["additional_scores"] == [
        {"score_name": "global_qvalue", "score_value": 0.25, "higher_better": None},
        {"score_name": "consensus_support", "score_value": 0.5, "higher_better": None},
    ]


def test_an_msnet_psm_of_writer_1_2_1_dev5_gets_listed_scores_accessions_and_targets(tmp_path):
    """One struct of scores, empty accessions, null is_decoy and cv values, and the dot before an N-terminal tag."""
    folder = tmp_path / "msnet"
    folder.mkdir()
    score = pa.struct([("score_name", pa.string()), ("score_value", pa.float64()), ("higher_better", pa.bool_())])
    position = pa.struct([("position", pa.int32()), ("amino_acid", pa.string()), ("scores", score)])
    modification = pa.struct([("name", pa.string()), ("accession", pa.string()), ("positions", pa.list_(position))])
    schema = pa.schema(
        [
            ("sequence", pa.string()),
            ("peptidoform", pa.string()),
            ("modifications", pa.list_(modification)),
            ("charge", pa.int16()),
            ("is_decoy", pa.bool_()),
            ("run_file_name", pa.string()),
            ("scan", pa.list_(pa.int32())),
            ("cv_params", pa.list_(pa.struct([("cv_name", pa.string()), ("cv_value", pa.string())]))),
            ("consensus_support", pa.float32()),
        ],
        metadata={b"qpx_version": b"1.2", b"creator": b"chengxin", b"writer_version": b"1.2.1.dev5"},
    )
    empty = {"score_name": "", "score_value": None, "higher_better": True}
    flr = {"score_name": "Luciphor_local_flr", "score_value": 0.01, "higher_better": False}
    row = {
        "sequence": "PEPSK",
        "peptidoform": ".[TMTpro]PEPS[Phospho]K",
        "modifications": [
            {"name": "TMTpro", "accession": "", "positions": [{"position": 0, "amino_acid": None, "scores": empty}]},
            {"name": "Phospho", "accession": "", "positions": [{"position": 4, "amino_acid": "S", "scores": flr}]},
        ],
        "charge": 2,
        "is_decoy": None,
        "run_file_name": RUN,
        "scan": [7],
        "cv_params": [{"cv_name": "instrument", "cv_value": "Orbitrap"}, {"cv_name": "isotope error", "cv_value": None}],
        "consensus_support": 1.0,
    }
    pq.write_table(pa.Table.from_pylist([row], schema=schema), folder / "msnet.psm.parquet")
    migrate_dataset(folder, tmp_path / "new")

    row = pq.read_table(tmp_path / "new" / "msnet.psm.parquet").to_pylist()[0]
    assert (row["peptidoform"], row["is_decoy"]) == ("[UNIMOD:2016]-PEPS[UNIMOD:21]K", False)
    assert row["cv_params"] == [{"cv_name": "instrument", "cv_value": "Orbitrap"}]
    assert row["modifications"] == [
        {"name": "TMTpro", "accession": "UNIMOD:2016", "positions": [{"position": 0, "amino_acid": None, "scores": None}]},
        {"name": "Phospho", "accession": "UNIMOD:21", "positions": [{"position": 4, "amino_acid": "S", "scores": [flr]}]},
    ]


@pytest.mark.parametrize(
    ("runs", "names", "expected"),
    [
        (
            ["a_f01", "b_R1_CID", "c_8%_DDA"],
            ["a_f01_consensus_fdr_filter.idXML", "b_R1", "c_8_DDA"],
            ["a_f01", "b_R1_CID", "c_8%_DDA"],
        ),
        (["a_f01"], ["a_f01_consensus_fdr_filter.idXML", "S1"], ["a_f01_consensus_fdr_filter.idXML", "S1"]),
        (["b_R1_CID"], ["b_R1", "b_R1_CID"], ["b_R1", "b_R1_CID"]),
    ],
    ids=["respelled", "a-name-matching-no-run", "a-run-named-already"],
)
def test_run_names_take_the_run_views_spelling_when_each_matches_one_run(tmp_path, runs, names, expected):
    """Old converters wrote an idXML suffix, left out a dissociation suffix or dropped a '%' from a psm run name."""
    folder = tmp_path / "old"
    folder.mkdir()
    samples = [{"sample_accession": "S1", "label": "LFQ"}]
    _write(folder / "old.run.parquet", [{"run_accession": run, "run_file_name": run, "samples": samples} for run in runs])
    _write(folder / "old.psm.parquet", [make_psm_record(run_file_name=name) for name in names])
    migrate_dataset(folder, tmp_path / "new")

    assert pq.read_table(tmp_path / "new" / "old.psm.parquet").column("run_file_name").to_pylist() == expected
