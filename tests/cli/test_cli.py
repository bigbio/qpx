"""CLI command tests — convert, transform, query, info, validate, ontology."""

import json
import re

import pandas as pd
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main

# ---------------------------------------------------------------------------
# Convert
# ---------------------------------------------------------------------------


def _assert_help(result, *options):
    """Verify CLI help renders and contains expected options."""
    if result.exit_code != 0:
        raise AssertionError(f"exit_code={result.exit_code}, output={result.output}")
    for opt in options:
        if opt not in result.output:
            raise AssertionError(f"Missing option {opt} in help output")


class TestDiaNNConvertCLI:
    def test_diann_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "diann", "--help"])
        _assert_help(result, "--report-path")


class TestMaxQuantConvertCLI:
    def test_maxquant_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "maxquant", "--help"])
        _assert_help(result, "--msms-file")


class TestFragPipeConvertCLI:
    def test_fragpipe_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "fragpipe", "--help"])
        _assert_help(result, "--psm-file")

    def test_fragpipe_help_shows_new_options(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "fragpipe", "--help"])
        _assert_help(
            result,
            "--ion-file",
            "--pg-file",
            "--experiment-annotation-file",
        )
        assert "--peptide-file" not in result.output


class TestMzIdentMLConvertCLI:
    def test_mzidentml_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "mzidentml", "--help"])
        _assert_help(result, "--mzid-path")

    def test_mzidentml_help_shows_new_options(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "mzidentml", "--help"])
        _assert_help(result, "--mgf-path", "--include-spectra", "--project-accession")


class TestSdrfConvertCLI:
    def test_sdrf_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "sdrf", "--help"])
        _assert_help(result, "--sdrf-file")


# ---------------------------------------------------------------------------
# Transform
# ---------------------------------------------------------------------------


class TestTransformGeneMapCLI:
    def test_genemap_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "gene-map", "--help"])
        _assert_help(result, "--parquet-path", "--dataset", "--in-place", "--fasta")

    def test_gene_map_requires_exactly_one_input_mode(self, tmp_path):
        fasta = tmp_path / "db.fasta"
        fasta.write_text(">sp|P12345|A_HUMAN a GN=BRCA1\nMKV\n")
        runner = CliRunner()

        neither = runner.invoke(qpx_main, ["transform", "gene-map", "--fasta", str(fasta)])
        assert neither.exit_code != 0
        assert "exactly one of --parquet-path or --dataset" in neither.output

        parquet = tmp_path / "openms.pg.parquet"
        parquet.touch()
        both = runner.invoke(
            qpx_main,
            [
                "transform",
                "gene-map",
                "--fasta",
                str(fasta),
                "--parquet-path",
                str(parquet),
                "--dataset",
                str(tmp_path),
            ],
        )
        assert both.exit_code != 0
        assert "exactly one of --parquet-path or --dataset" in both.output

    def test_gene_map_dataset_rejects_conflicting_destinations(self, tmp_path):
        """--in-place and --output-folder together, or a destination inside the source."""
        fasta = tmp_path / "db.fasta"
        fasta.write_text(">sp|P12345|A_HUMAN a GN=BRCA1\nMKV\n")
        data = tmp_path / "qpx_output"
        data.mkdir()
        (data / "openms.pg.parquet").touch()
        runner = CliRunner()

        base = ["transform", "gene-map", "--fasta", str(fasta), "--dataset", str(data)]

        both = runner.invoke(qpx_main, base + ["--in-place", "--output-folder", str(tmp_path / "out")])
        assert both.exit_code != 0
        assert "not both" in both.output

        same = runner.invoke(qpx_main, base + ["--output-folder", str(data)])
        assert same.exit_code != 0
        assert "use --in-place" in same.output

        nested = runner.invoke(qpx_main, base + ["--output-folder", str(data / "annotated")])
        assert nested.exit_code != 0
        assert "must not be inside the dataset" in nested.output

    def test_gene_map_dataset_requires_a_destination(self, tmp_path):
        fasta = tmp_path / "db.fasta"
        fasta.write_text(">sp|P12345|A_HUMAN a GN=BRCA1\nMKV\n")
        (tmp_path / "openms.pg.parquet").touch()

        runner = CliRunner()
        result = runner.invoke(
            qpx_main,
            ["transform", "gene-map", "--fasta", str(fasta), "--dataset", str(tmp_path)],
        )

        assert result.exit_code != 0
        assert "--in-place or --output-folder" in result.output


class TestTransformQuantifyCLI:
    def test_quantify_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "quantify", "--help"])
        _assert_help(result, "--feature-path", "--method", "--output")

    def test_quantify_help_shows_ibaq_options(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "quantify", "--help"])
        _assert_help(result, "--organism", "--ploidy", "--min-aa", "--sdrf")

    def test_ibaq_validation_does_not_require_directlfq_helper(self, monkeypatch, tmp_path):
        """iBAQ availability must not depend on a DirectLFQ-only helper."""
        import sys
        from types import ModuleType

        from qpx.cli.transform.quantify import _validate_quantify_inputs

        mokume = ModuleType("mokume")
        mokume.__path__ = []
        quantification = ModuleType("mokume.quantification")
        monkeypatch.setitem(sys.modules, "mokume", mokume)
        monkeypatch.setitem(sys.modules, "mokume.quantification", quantification)

        _validate_quantify_inputs("ibaq", tmp_path / "database.fasta")

    def test_quantify_maps_tmt_channels_to_sdrf_samples(self, tmp_path):
        """SDRF run/channel pairs restore biological samples and conditions."""
        import pandas as pd

        from qpx.cli.transform.quantify import _qpx_feature_to_peptide_df

        feature_path = tmp_path / "feature.parquet"
        pd.DataFrame(
            {
                "sequence": ["PEPTIDEK"],
                "anchor_protein": ["P1"],
                "run_file_name": ["fraction1"],
                "is_decoy": [False],
                "intensities": [
                    [
                        {"label": "TMT126", "intensity": 100.0},
                        {"label": "TMT127N", "intensity": 200.0},
                    ]
                ],
            }
        ).to_parquet(feature_path, index=False)
        sdrf_path = tmp_path / "experiment.sdrf.tsv"
        pd.DataFrame(
            {
                "source name": ["sample-a", "sample-b"],
                "comment[data file]": ["fraction1.raw", "fraction1.raw"],
                "comment[label]": ["TMT126", "TMT127N"],
                "factor value[disease]": ["control", "case"],
            }
        ).to_csv(sdrf_path, sep="\t", index=False)

        result = _qpx_feature_to_peptide_df(feature_path, sdrf_path)

        observed = result[["SampleID", "Condition", "NormIntensity"]].to_dict("records")
        expected = [
            {"SampleID": "sample-a", "Condition": "control", "NormIntensity": 100.0},
            {"SampleID": "sample-b", "Condition": "case", "NormIntensity": 200.0},
        ]
        if observed != expected:
            raise AssertionError(f"Unexpected SDRF mapping: {observed!r}")

    def test_quantify_marks_condition_unavailable_without_sdrf(self, tmp_path):
        """Standalone Feature input supplies an explicit missing condition."""
        import pandas as pd

        from qpx.cli.transform.quantify import _qpx_feature_to_peptide_df

        feature_path = tmp_path / "feature.parquet"
        pd.DataFrame(
            {
                "sequence": ["PEPTIDEK"],
                "anchor_protein": ["P1"],
                "run_file_name": ["run1"],
                "intensities": [[{"label": "TMT126", "intensity": 100.0}]],
            }
        ).to_parquet(feature_path, index=False)

        result = _qpx_feature_to_peptide_df(feature_path)

        if result.loc[0, "Condition"] != "not available":
            raise AssertionError("Missing SDRF condition was not marked unavailable")

    def test_quantify_flattens_quantified_target_entries(self, tmp_path):
        """Every positive target intensity becomes one row keyed by run and label."""
        import pandas as pd

        from qpx.cli.transform.quantify import _qpx_feature_to_peptide_df

        feature_path = tmp_path / "feature.parquet"
        pd.DataFrame(
            {
                "sequence": ["AAK", "CCK", "DDK", "EEK", "FFK"],
                "anchor_protein": ["P1", "P2", "P3", "P4", None],
                "run_file_name": ["run1", "run2", "run1", "run2", "run1"],
                "is_decoy": [False, False, True, False, False],
                "intensities": [
                    [{"label": "TMT126", "intensity": 10.0}, {"label": "TMT127N", "intensity": 0.0}],
                    [{"label": "TMT127N", "intensity": None}, {"label": "TMT126", "intensity": 5.0}],
                    [{"label": "TMT126", "intensity": 7.0}],
                    [],
                    [{"label": "TMT126", "intensity": 3.0}],
                ],
            }
        ).to_parquet(feature_path, index=False)

        result = _qpx_feature_to_peptide_df(feature_path)

        observed = result[["ProteinName", "NormIntensity", "SampleID"]].to_dict("records")
        expected = [
            {"ProteinName": "P1", "NormIntensity": 10.0, "SampleID": "run1::TMT126"},
            {"ProteinName": "P2", "NormIntensity": 5.0, "SampleID": "run2::TMT126"},
        ]
        if observed != expected:
            raise AssertionError(f"Unexpected flattened rows: {observed!r}")


class TestTransformNormalizeAccessionsCLI:
    def test_normalize_accessions_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "normalize-accessions", "--help"])
        _assert_help(result, "--dataset", "--direction", "--fasta")

    def test_normalize_accessions_discovers_openms_prefix(self, tmp_path, monkeypatch):
        (tmp_path / "openms.feature.parquet").touch()
        (tmp_path / "openms.pg.parquet").touch()
        normalized = []

        def fake_normalize(**kwargs):
            normalized.append(kwargs["parquet_path"].name)
            return {"rows": 1, "accessions_changed": 0, "distinct_changed": 0}

        monkeypatch.setattr("qpx.transforms.accession_normalizer.normalize_parquet", fake_normalize)

        runner = CliRunner()
        result = runner.invoke(
            qpx_main,
            ["transform", "normalize-accessions", "--dataset", str(tmp_path), "--in-place"],
        )

        assert result.exit_code == 0, result.output
        assert normalized == ["openms.feature.parquet", "openms.pg.parquet"]

    @pytest.mark.parametrize("in_place", [False, True])
    def test_normalize_accessions_preserves_site_proteoform_pg_links(self, tmp_path, in_place):
        """Normalize linked site/proteoform accessions without changing ids or other fields."""
        import pyarrow.parquet as pq

        from qpx.dataset import Dataset
        from qpx.writers import FeatureWriter, PgWriter, ProteoformWriter, SiteWriter
        from tests.conftest import make_feature_record, make_pg_record
        from tests.unit.test_site_proteoform import make_proteoform_record, make_site_record

        source = tmp_path / "input"
        source.mkdir()
        family, protein, accession = "sp|P04637|P53_HUMAN", "sp|P04637-2|P53_HUMAN", "sp|P04637-3|P53_HUMAN"
        with PgWriter(source / "exp.pg.parquet") as writer:
            writer.write_batch([make_pg_record(anchor_protein=family, pg_accessions=[family, protein], intensities=[])])
        pg_id = pq.read_table(source / "exp.pg.parquet").column("pg_id")[0].as_py()
        feature = make_feature_record(anchor_protein=family)
        feature["pg_accessions"] = [{"accession": a, "start": None, "end": None} for a in (family, protein)]
        feature["pg_ids"] = [pg_id]
        records = {
            "feature": (FeatureWriter, feature),
            "site": (SiteWriter, make_site_record(protein_accession=protein, family=family, pg_id=pg_id)),
            "proteoform": (
                ProteoformWriter,
                make_proteoform_record(family=family, accession=accession, kind="isoform", pg_id=pg_id),
            ),
        }
        for view, (writer_class, record) in records.items():
            with writer_class(source / f"exp.{view}.parquet") as writer:
                writer.write_batch([record])
        before = {view: pq.read_table(source / f"exp.{view}.parquet") for view in ("site", "proteoform")}
        with Dataset(source) as dataset:
            assert all(result.is_valid for result in dataset.validate(strict=True).values())
        output = source if in_place else tmp_path / "output"
        destination = ["--in-place"] if in_place else ["--output", str(output)]
        result = CliRunner().invoke(qpx_main, ["transform", "normalize-accessions", "--dataset", str(source), *destination])
        assert result.exit_code == 0, result.output
        with Dataset(output) as dataset:
            validation = dataset.validate(strict=True)
            assert all(result.is_valid for result in validation.values()), validation
        expected = {
            "site": {"protein_accession": "P04637-2", "family": "P04637"},
            "proteoform": {"family": "P04637", "accession": "P04637-3"},
        }
        for view, columns in expected.items():
            table = pq.read_table(output / f"exp.{view}.parquet")
            for name, value in columns.items():
                assert table.column(name).to_pylist() == [value]
            assert table.drop(list(columns)).equals(before[view].drop(list(columns)))
            assert f"exp.{view}.parquet: 1 rows, 1 accessions normalized" in result.output


class TestTransformUpdateMetadataCLI:
    def test_update_metadata_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "update-metadata", "--help"])
        _assert_help(result, "--dataset", "--sdrf", "--old-sdrf", "--force")


# ---------------------------------------------------------------------------
# Query
# ---------------------------------------------------------------------------


class TestQuerySqlCLI:
    def test_sql_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["query", "sql", "--help"])
        _assert_help(result, "--dataset-path", "--sql")


class TestQueryFilterCLI:
    def test_filter_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["query", "filter", "--help"])
        _assert_help(result, "--dataset-path", "--structure", "--condition")

    def test_filter_reads_dataset_metadata(self, dataset_dir):
        """The canonical dataset name resolves to the dataset_meta accessor."""
        runner = CliRunner()
        result = runner.invoke(
            qpx_main,
            [
                "query",
                "filter",
                "--dataset-path",
                str(dataset_dir),
                "--structure",
                "dataset",
                "--condition",
                "project_accession IS NOT NULL",
                "--output-format",
                "json",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "project_accession" in result.output


class TestQueryHeadCLI:
    def test_head_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["query", "head", "--help"])
        _assert_help(result, "--dataset-path", "--structure")

    def test_head_reads_dataset_metadata(self, dataset_dir):
        """Query head accepts the canonical dataset structure name."""
        runner = CliRunner()
        result = runner.invoke(
            qpx_main,
            [
                "query",
                "head",
                "--dataset-path",
                str(dataset_dir),
                "--structure",
                "dataset",
                "--output-format",
                "csv",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "project_accession" in result.output


# ---------------------------------------------------------------------------
# Info
# ---------------------------------------------------------------------------


class TestInfoCLI:
    def test_info_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["info", "--help"])
        _assert_help(result, "--dataset-path")

    def test_info_counts_every_shard(self, tmp_path):
        """A view split over several files reports the rows of all of them."""
        from qpx.writers.psm import PsmWriter
        from tests.conftest import make_psm_record

        for shard, scan in (("a", 1001), ("b", 1002)):
            with PsmWriter(tmp_path / f"{shard}.psm.parquet") as writer:
                writer.write_batch([make_psm_record(scan=[scan])])

        result = CliRunner().invoke(qpx_main, ["info", "--dataset-path", str(tmp_path)])

        assert result.exit_code == 0, result.output
        assert re.search(r"^\s*psm\s+2\s", result.output, re.MULTILINE), result.output

    def test_info_schema_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["info", "schema", "--help"])
        _assert_help(result, "--dataset-path")

    def test_info_metadata_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["info", "metadata", "--help"])
        _assert_help(result, "--file")


# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------


class TestValidateCLI:
    def test_validate_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["validate", "--help"])
        _assert_help(result, "--dataset-path", "--file", "--structure", "--semantic")

    def test_validate_semantic_fails_a_wrong_calculated_mz(self, tmp_path):
        """--semantic checks calculated_mz against the peptidoform; the default run does not."""
        from qpx.writers import FeatureWriter
        from tests.conftest import make_feature_record

        path = tmp_path / "a.feature.parquet"
        with FeatureWriter(path) as w:
            w.write_batch([make_feature_record(calculated_mz=500.0, observed_mz=500.0)])
        runner = CliRunner()
        assert runner.invoke(qpx_main, ["validate", "--file", str(path)]).exit_code == 0
        result = runner.invoke(qpx_main, ["validate", "--file", str(path), "--semantic"])
        assert result.exit_code == 1
        assert "calculated_mz" in result.output

    def test_validate_reports_an_unreadable_file(self, tmp_path):
        """A pre-1.1 pg file is reported as invalid, without a traceback."""
        from tests.unit.test_version import _write_old_pg_file

        path = tmp_path / "old.pg.parquet"
        _write_old_pg_file(path)
        result = CliRunner().invoke(qpx_main, ["validate", "--file", str(path)])
        assert result.exit_code == 1
        assert isinstance(result.exception, SystemExit)
        assert "pg: INVALID" in result.output
        assert "could not be read" in result.output


# ---------------------------------------------------------------------------
# Ontology
# ---------------------------------------------------------------------------


class TestOntologyCLI:
    def test_ontology_info_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["ontology", "info", "--help"])
        _assert_help(result, "--source")

    def test_ontology_search_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["ontology", "search", "--help"])
        _assert_help(result, "--source", "--top-k")

    def test_ontology_update_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["ontology", "update", "--help"])
        _assert_help(result, "--source")

    def test_ontology_build_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["ontology", "build", "--help"])
        _assert_help(result, "--source", "--all-sources")


class TestStructureChoicesMatchTheRegistry:
    """CLI structure choices must be derived, not hand-maintained.

    pepmap was registered in Dataset._STRUCTURE_REGISTRY and exposed on the
    public API while both CLIs rejected it, because each kept its own literal
    list (bigbio/qpx#289).
    """

    def test_validate_offers_every_registered_structure(self):
        from qpx.cli.validate import _VALID_STRUCTURES
        from qpx.dataset import Dataset

        assert set(_VALID_STRUCTURES) == set(Dataset._STRUCTURE_REGISTRY)

    def test_query_offers_every_registered_structure(self):
        from qpx.cli.query import _VALID_STRUCTURES
        from qpx.dataset import Dataset

        assert set(_VALID_STRUCTURES) == set(Dataset._STRUCTURE_REGISTRY)

    def test_info_schema_offers_every_registered_structure(self):
        from qpx.cli.info import _VALID_STRUCTURES
        from qpx.dataset import Dataset

        assert set(_VALID_STRUCTURES) == set(Dataset._STRUCTURE_REGISTRY)

    def test_pepmap_is_accepted(self):
        """The structure that exposed the drift."""
        from qpx.cli.query import _VALID_STRUCTURES as query_structures
        from qpx.cli.validate import _VALID_STRUCTURES as validate_structures

        assert "pepmap" in validate_structures
        assert "pepmap" in query_structures


def test_query_json_writes_ids_as_decimal_strings():
    """A JavaScript number rounds a 64-bit id, so JSON carries it as a string."""
    from qpx.cli.query import _json_records

    frame = pd.DataFrame(
        {
            "feature_id": pd.array([9007199254740993, None], dtype="Int64"),
            "psm_ids": pd.Series([[9007199254740995], pd.NA], dtype=object),
            "charge": [2, 3],
        }
    )
    records = json.loads(_json_records(frame))
    assert records[0] == {"feature_id": "9007199254740993", "psm_ids": ["9007199254740995"], "charge": 2}
    assert records[1]["feature_id"] is None
    assert records[1]["psm_ids"] is None


def _commands(group, path=()):
    for name in group.list_commands(None):
        command = group.get_command(None, name)
        if hasattr(command, "list_commands"):
            yield from _commands(command, (*path, name))
        else:
            yield " ".join((*path, name)), command


@pytest.mark.parametrize("name, command", list(_commands(qpx_main)), ids=lambda value: value if isinstance(value, str) else "")
def test_help_keeps_each_example_on_its_own_line(name, command):
    """Click rewraps a help paragraph unless it starts with the \\b control character."""
    runner_result = CliRunner().invoke(qpx_main, [*name.split(), "--help"])
    description = runner_result.output.split("Options:")[0]
    assert "\\b" not in description
    assert not re.search(r"\\ +--", description)
