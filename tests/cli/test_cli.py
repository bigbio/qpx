"""CLI command tests — convert, transform, query, info, validate, ontology."""

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

    def test_mzidentml_help_shows_enrich_pride(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "mzidentml", "--help"])
        _assert_help(result, "--enrich-pride")


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
        _assert_help(result, "--fasta")


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

        from qpx.cli.transform import _validate_quantify_inputs

        mokume = ModuleType("mokume")
        mokume.__path__ = []
        quantification = ModuleType("mokume.quantification")
        monkeypatch.setitem(sys.modules, "mokume", mokume)
        monkeypatch.setitem(sys.modules, "mokume.quantification", quantification)

        _validate_quantify_inputs("ibaq", tmp_path / "database.fasta")

    def test_quantify_maps_tmt_channels_to_sdrf_samples(self, tmp_path):
        """SDRF run/channel pairs restore biological samples and conditions."""
        import pandas as pd

        from qpx.cli.transform import _qpx_feature_to_peptide_df

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

        from qpx.cli.transform import _qpx_feature_to_peptide_df

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
            return {"rows": 1, "accessions_changed": 0}

        monkeypatch.setattr("qpx.transforms.accession_normalizer.normalize_parquet", fake_normalize)

        runner = CliRunner()
        result = runner.invoke(
            qpx_main,
            ["transform", "normalize-accessions", "--dataset", str(tmp_path), "--in-place"],
        )

        assert result.exit_code == 0, result.output
        assert normalized == ["openms.feature.parquet", "openms.pg.parquet"]


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
        _assert_help(result, "--dataset-path", "--file", "--structure")


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

    def test_pepmap_is_accepted(self):
        """The structure that exposed the drift."""
        from qpx.cli.query import _VALID_STRUCTURES as query_structures
        from qpx.cli.validate import _VALID_STRUCTURES as validate_structures

        assert "pepmap" in validate_structures
        assert "pepmap" in query_structures
