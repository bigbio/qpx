"""Existing OpenMS SILAC quantities and identities survive QPX enrichment."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.converters.channel_labels import experiment_type_from_labels
from qpx.converters.openms.converter import OpenMSConverter
from qpx.dataset import Dataset
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter
from tests.conftest import make_feature_record, make_pg_record, make_psm_record


@pytest.mark.parametrize(
    "labels",
    [
        {"SILAC light", "SILAC heavy"},
        {"SILAC light", "SILAC medium", "SILAC heavy"},
        {"NT=silac LIGHT;AC=MS:1002038", "NT=SILAC HEAVY;AC=MS:1002040"},
        {"SILAC medium"},
    ],
)
def test_silac_sdrf_labels_are_not_classified_as_lfq(labels):
    """Canonical SILAC labels, including ontology notation, identify the family."""
    assert experiment_type_from_labels(labels) == "SILAC"


def _quantities(labels, scale):
    """Return distinguishable primary and additional quantities per channel."""
    primary = [{"label": label, "intensity": scale * (index + 1)} for index, label in enumerate(labels)]
    additional = [
        {"label": row["label"], "intensities": [{"intensity_name": "normalized", "intensity_value": row["intensity"] / 2}]}
        for row in primary
    ]
    return primary, additional


def _write_core(folder, labels, legacy_pg=False):
    """Write already-labelled OpenMS input with producer-assigned identities."""
    primary, additional = _quantities(labels, 100.0)
    features = []
    for run in ("run_01", "run_02"):
        feature = make_feature_record(run_file_name=run, peptidoform="PEPTIDEK[UNIMOD:259]", intensities=primary)
        feature["additional_intensities"] = additional
        features.append(feature)
    with FeatureWriter(folder / "source.feature.parquet") as writer:
        writer.write_batch(features)
    with PsmWriter(folder / "source.psm.parquet") as writer:
        writer.write_batch([make_psm_record(peptidoform="PEPTIDEK[UNIMOD:259]")])

    primary, additional = _quantities(labels, 1000.0)
    protein = make_pg_record(intensities=primary)
    protein["grouped_runs"] = ["run_01", "run_02"]
    protein["additional_intensities"] = additional
    if legacy_pg:
        protein.pop("grouped_runs")
        protein["run_file_name"] = "run_01"
        pq.write_table(pa.Table.from_pylist([protein]), folder / "source.pg.parquet")
    else:
        with PgWriter(folder / "source.pg.parquet") as writer:
            writer.write_batch([protein])


def _write_sdrf(path, labels):
    """Describe the same SILAC mixture across two fractions."""
    rows = [
        "source name\tcharacteristics[organism]\tcharacteristics[organism part]\t"
        "comment[data file]\tcomment[label]\tcomment[fraction identifier]"
    ]
    for fraction, run in enumerate(("run_01", "run_02"), 1):
        for channel, label in enumerate(labels, 1):
            rows.append(f"sample_{channel}\tHomo sapiens\tliver\t{run}.raw\tNT={label}\t{fraction}")
    path.write_text("\n".join(rows) + "\n")


def _maplist(path, count):
    """Model companion XML labels that must not replace canonical QPX labels."""
    raw_labels = ("no_label", "Lys4Arg6", "Lys8Arg10")
    maps = "".join(f'<map id="{index}" name="run_01.mzML" label="{raw_labels[index]}"/>' for index in range(count))
    path.write_text(f'<consensusXML><mapList count="{count}">{maps}</mapList></consensusXML>')


@pytest.mark.parametrize("labels", [("SILAC light", "SILAC heavy"), ("SILAC light", "SILAC medium", "SILAC heavy")])
@pytest.mark.parametrize("with_maplist", [False, True])
def test_native_silac_enrichment_preserves_core_and_sample_joins(tmp_path, labels, with_maplist):
    """Enrichment keeps every label, quantity, peptidoform, ID and source unit."""
    source = tmp_path / "source"
    source.mkdir()
    _write_core(source, labels)
    sdrf = tmp_path / "input.sdrf.tsv"
    _write_sdrf(sdrf, labels)
    consensus = tmp_path / "input.consensusXML" if with_maplist else None
    if consensus is not None:
        _maplist(consensus, len(labels))
    output = tmp_path / "out"
    OpenMSConverter(source, sdrf, consensus).convert(output)

    for view in ("feature", "psm", "pg"):
        before = pq.read_table(source / f"source.{view}.parquet")
        after = pq.read_table(output / f"openms.{view}.parquet")
        # The existing writer also records provided_feature_id as a CV parameter.
        assert before.drop_columns("cv_params").equals(after.drop_columns("cv_params"), check_metadata=False)
    with Dataset(output, file_prefix="openms", duckdb_threads=6) as dataset:
        protein = dataset.intensity("protein").to_df().set_index("sample_accession")["intensity"].to_dict()
        peptide = dataset.intensity("peptide").to_df().set_index("sample_accession")["intensity"].to_dict()
    assert protein == {f"sample_{index}": 1000.0 * index for index in range(1, len(labels) + 1)}
    assert peptide == {f"sample_{index}": 200.0 * index for index in range(1, len(labels) + 1)}
    provenance = pq.read_table(output / "openms.provenance.parquet").to_pylist()[0]
    assert provenance["step_name"] == "silac_quantification"
    assert provenance["tool_name"] == "OpenMS"
    metadata = pq.read_table(output / "openms.dataset.parquet").to_pylist()[0]
    assert metadata["software_name"] == "OpenMS"


def test_legacy_silac_pg_upgrade_retains_every_channel(tmp_path):
    """Legacy list-shaped protein quantities flatten without relabelling channels."""
    labels = ("SILAC light", "SILAC medium", "SILAC heavy")
    source = tmp_path / "source"
    source.mkdir()
    _write_core(source, labels, legacy_pg=True)
    sdrf = tmp_path / "input.sdrf.tsv"
    _write_sdrf(sdrf, labels)
    output = tmp_path / "out"
    OpenMSConverter(source, sdrf).convert(output)
    proteins = pq.read_table(output / "openms.pg.parquet").to_pylist()

    assert len(proteins) == len(labels)
    assert len({row["pg_id"] for row in proteins}) == len(labels)
    for index, row in enumerate(proteins, 1):
        assert row["label"] == labels[index - 1]
        assert row["intensity"] == 1000.0 * index
        assert row["grouped_runs"] == ["run_01"]
        assert row["additional_intensities"] == [
            {"label": labels[index - 1], "intensities": [{"intensity_name": "normalized", "intensity_value": 500.0 * index}]}
        ]
