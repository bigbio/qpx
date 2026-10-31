"""``transform quantify`` runs the named method and writes the named format; AnnData views read and validate."""

import json
import sys
from types import ModuleType

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.cli.transform.quantify import _qpx_feature_to_peptide_df
from qpx.core.anndata_validation import validate_anndata_view
from qpx.dataset import Dataset
from qpx.mudata import _build_differential_adata
from qpx.writers.ontology import OntologyWriter

_ENTRY = {"label": "LFQ", "intensity": 100.0}


def _write_feature(path):
    pd.DataFrame(
        {
            "sequence": ["AAK", "AAK", "CCK", "DDK", "EEK"],
            "anchor_protein": ["P1", "P1", "P1", "P2", ""],
            "run_file_name": ["run1", "run1", "run1", "run1", "run1"],
            "is_decoy": [False, False, False, True, False],
            "intensities": [
                [_ENTRY],
                [{"label": "LFQ", "intensity": 50.0}],
                [{"label": "LFQ", "intensity": None}],
                [_ENTRY],
                [_ENTRY],
            ],
        }
    ).to_parquet(path, index=False)
    return path


@pytest.fixture(name="fake_mokume")
def _fake_mokume(monkeypatch):
    """A stand-in recording how quantify calls mokume, whose installed version varies."""
    calls = {}

    class MaxLFQQuantification:
        def __init__(self, **kwargs):
            calls["maxlfq"] = kwargs

        def quantify(self, peptide_df, **_):
            return peptide_df.groupby(["ProteinName", "SampleID"], as_index=False)["NormIntensity"].sum()

    def get_quantification_method(method, **kwargs):
        calls[method] = kwargs
        return MaxLFQQuantification()

    def peptides_to_protein(**kwargs):
        calls["ibaq"] = kwargs
        pd.DataFrame({"ProteinName": ["P1"], "SampleID": ["run1::LFQ"], "Ibaq": [8.0], "IbaqPpb": [100.0]}).to_csv(
            kwargs["output"], sep="\t", index=False
        )

    module = ModuleType("mokume.quantification")
    module.MaxLFQQuantification = MaxLFQQuantification
    module.get_quantification_method = get_quantification_method
    module.peptides_to_protein = peptides_to_protein
    module.is_directlfq_available = lambda: True
    package = ModuleType("mokume")
    package.__path__ = []
    monkeypatch.setitem(sys.modules, "mokume", package)
    monkeypatch.setitem(sys.modules, "mokume.quantification", module)
    return calls


def _quantify(feature, output, *arguments):
    return CliRunner().invoke(
        qpx_main,
        ["transform", "quantify", "--feature-path", str(feature), "-o", str(output), *arguments],
        catch_exceptions=False,
    )


def test_peptide_rows_are_summed_within_a_run_and_label(tmp_path):
    peptides = _qpx_feature_to_peptide_df(_write_feature(tmp_path / "x.feature.parquet"))

    assert peptides.to_dict("records") == [
        {
            "ProteinName": "P1",
            "PeptideCanonical": "AAK",
            "NormIntensity": 150.0,
            "SampleID": "run1::LFQ",
            "Condition": "not available",
        }
    ]


def test_a_dataset_directory_is_accepted(tmp_path):
    _write_feature(tmp_path / "x.feature.parquet")

    assert len(_qpx_feature_to_peptide_df(tmp_path)) == 1


def test_maxlfq_runs_maxlfq_with_the_requested_threads(tmp_path, fake_mokume):
    result = _quantify(_write_feature(tmp_path / "f.parquet"), tmp_path / "o.tsv", "--method", "maxlfq", "--threads", "3")

    assert result.exit_code == 0, result.output
    assert fake_mokume["maxlfq"] == {"threads": 3, "force_builtin": True}


def test_directlfq_receives_the_threads(tmp_path, fake_mokume):
    result = _quantify(_write_feature(tmp_path / "f.parquet"), tmp_path / "o.tsv", "--method", "directlfq", "--threads", "3")

    assert result.exit_code == 0, result.output
    assert fake_mokume["directlfq"] == {"num_cores": 3}


def test_ibaq_writes_the_format_the_suffix_names(tmp_path, fake_mokume):
    feature = _write_feature(tmp_path / "f.parquet")
    fasta = tmp_path / "db.fasta"
    fasta.write_text(">P1\nAAK\n")

    parquet = _quantify(feature, tmp_path / "o.parquet", "--method", "ibaq", "--fasta", str(fasta))
    anndata = _quantify(feature, tmp_path / "o.ae.h5ad", "--method", "ibaq", "--fasta", str(fasta))

    assert parquet.exit_code == 0, parquet.output
    assert pd.read_parquet(tmp_path / "o.parquet")["Ibaq"].tolist() == [8.0]
    assert anndata.exit_code == 0, anndata.output
    view = ad.read_h5ad(tmp_path / "o.ae.h5ad")
    assert view.X.tolist() == [[3.0]]
    assert view.layers["ibaq_raw"].tolist() == [[8.0]]
    assert view.layers["ibaq_ppb"].tolist() == [[100.0]]
    assert view.uns["file_type"] == "absolute_expression"
    assert view.uns["log_base"] == 2
    assert view.uns["var_key"] == "anchor_protein"
    assert json.loads(view.uns["quantification_parameters"])["method"] == "ibaq"
    assert not validate_anndata_view(tmp_path / "o.ae.h5ad", "ae").issues


def test_ibaq_warns_about_dropped_proteins(tmp_path, fake_mokume):
    feature = tmp_path / "f.parquet"
    pd.DataFrame(
        {
            "sequence": ["AAK", "CCK"],
            "anchor_protein": ["P1", "P2"],
            "run_file_name": ["run1", "run1"],
            "intensities": [[_ENTRY], [_ENTRY]],
        }
    ).to_parquet(feature, index=False)
    fasta = tmp_path / "db.fasta"
    fasta.write_text(">P1\nAAK\n")

    result = _quantify(feature, tmp_path / "o.tsv", "--method", "ibaq", "--fasta", str(fasta))

    assert result.exit_code == 0, result.output
    assert "1 protein(s) were not in the FASTA" in result.output


def test_h5ad_output_needs_ibaq(tmp_path, fake_mokume):
    result = CliRunner().invoke(
        qpx_main,
        ["transform", "quantify", "--feature-path", str(_write_feature(tmp_path / "f.parquet")), "-o", str(tmp_path / "o.h5ad")]
        + ["--method", "maxlfq"],
    )

    assert result.exit_code != 0
    assert "needs --method ibaq" in result.output


def _de_view(path, pvals=(0.01, 0.2)):
    adata = ad.AnnData(
        X=np.zeros((2, 3)),
        obs=pd.DataFrame(index=["S1", "S2"]),
        var=pd.DataFrame(index=["P1", "P2", "P3"]),
    )
    adata.uns["file_type"] = "differential_expression"
    adata.uns["de_results"] = {
        "a_vs_b": {
            "names": np.array(["P1", "P2"]),
            "logfoldchanges": np.array([1.5, -2.0]),
            "pvals": np.array(pvals),
            "pvals_adj": np.array([0.02, 0.4]),
            "condition_test": "a",
            "condition_reference": "b",
        },
        "c_vs_b": {
            "names": np.array(["P3"]),
            "logfoldchanges": np.array([0.5]),
            "pvals": np.array([0.5]),
            "pvals_adj": np.array([0.5]),
            "condition_test": "c",
            "condition_reference": "b",
        },
    }
    adata.write_h5ad(path)
    return path


def test_de_results_layout_of_the_specification_reshapes_with_nan_for_untested(tmp_path):
    _de_view(tmp_path / "x.de.h5ad")

    matrix = _build_differential_adata(tmp_path, "x")

    assert list(matrix.obs_names) == ["a_vs_b", "c_vs_b"]
    assert list(matrix.var_names) == ["P1", "P2", "P3"]
    assert matrix["a_vs_b", "P1"].X[0, 0] == 1.5
    assert np.isnan(matrix["a_vs_b", "P3"].X[0, 0])
    assert np.isnan(matrix["c_vs_b", "P1"].layers["pvals_adj"][0, 0])


def test_validate_reports_de_layout_problems(tmp_path):
    with OntologyWriter(tmp_path / "x.ontology.parquet") as writer:
        writer.write_batch([{"field_name": "charge", "view": "feature"}])
    _de_view(tmp_path / "x.de.h5ad", pvals=(1.5, 0.2))

    with Dataset(str(tmp_path), file_prefix="x") as dataset:
        result = dataset.validate()["de"]

    assert [issue.message for issue in result.errors] == ["x.de.h5ad: contrast 'a_vs_b' has pvals outside [0, 1]"]


def test_validate_reports_a_de_view_without_results(tmp_path):
    ad.AnnData(X=np.zeros((1, 1))).write_h5ad(tmp_path / "x.de.h5ad")

    result = validate_anndata_view(tmp_path / "x.de.h5ad", "de")

    assert [issue.check for issue in result.errors] == ["anndata_layout"]
    assert [issue.check for issue in result.warnings] == ["anndata_file_type"]


def _ae_view(path, x=(3.0,), raw=(8.0,), **uns):
    adata = ad.AnnData(
        X=np.array([x], dtype="float64"),
        obs=pd.DataFrame(index=["S1"]),
        var=pd.DataFrame(index=[f"P{i}" for i in range(len(x))]),
    )
    adata.layers["ibaq_raw"] = np.array([raw], dtype="float64")
    adata.uns.update({"file_type": "absolute_expression", "log_base": 2, "var_key": "pg_accessions", "aggregation": "x", **uns})
    adata.write_h5ad(path)
    return path


def test_ae_validation_reads_what_x_is_declared_to_be(tmp_path):
    assert not validate_anndata_view(_ae_view(tmp_path / "ok.ae.h5ad"), "ae").issues

    wrong_base = validate_anndata_view(_ae_view(tmp_path / "base.ae.h5ad", log_base=10), "ae")
    assert [issue.message for issue in wrong_base.issues] == ["base.ae.h5ad: X is not the base 10 logarithm of layer ibaq_raw"]

    # A protein not detected in a sample is NaN whatever ibaq_raw says.
    undetected = _ae_view(tmp_path / "nan.ae.h5ad", x=(3.0, np.nan), raw=(8.0, 0.0))
    assert not validate_anndata_view(undetected, "ae").issues


def test_ae_validation_reports_unusable_values_and_metadata(tmp_path):
    infinite = validate_anndata_view(_ae_view(tmp_path / "inf.ae.h5ad", x=(np.inf,)), "ae")
    assert any("infinite" in issue.message for issue in infinite.issues)

    negative = validate_anndata_view(_ae_view(tmp_path / "neg.ae.h5ad", raw=(-1.0,)), "ae")
    assert any("layer ibaq_raw has negative" in issue.message for issue in negative.issues)

    key = validate_anndata_view(_ae_view(tmp_path / "key.ae.h5ad", var_key="pg_id"), "ae")
    assert any("var_key" in issue.message for issue in key.issues)

    table = pd.DataFrame({"protein": ["S1"], "amount": [50.0]})
    standards = validate_anndata_view(_ae_view(tmp_path / "std.ae.h5ad", standards=table), "ae")
    assert any("standards" in issue.message for issue in standards.issues)


def test_ae_validation_warns_about_unrecorded_metadata(tmp_path):
    adata = ad.read_h5ad(_ae_view(tmp_path / "w.ae.h5ad"))
    del adata.uns["log_base"], adata.uns["var_key"], adata.uns["aggregation"]
    adata.write_h5ad(tmp_path / "w.ae.h5ad")

    result = validate_anndata_view(tmp_path / "w.ae.h5ad", "ae")

    assert not result.errors
    assert [issue.message for issue in result.warnings] == [
        "w.ae.h5ad: uns has no log_base",
        "w.ae.h5ad: uns has no var_key",
        "w.ae.h5ad: uns has no aggregation",
    ]


def _omnibus_view(path, **uns):
    adata = ad.AnnData(X=np.zeros((1, 1)), obs=pd.DataFrame(index=["S1"]), var=pd.DataFrame(index=["P1"]))
    adata.uns.update({"file_type": "differential_expression", **uns})
    adata.uns["de_results"] = {
        "any_dose": {
            "names": np.array(["P1"]),
            "pvals": np.array([0.01]),
            "pvals_adj": np.array([0.02]),
            "df": np.array([29.05]),
            "statistic_type": "F",
            "contrast_definition": "dose_low=dose_high=0",
        }
    }
    adata.write_h5ad(path)
    return path


def test_de_validation_accepts_an_f_test_without_a_fold_change_or_a_pair_of_conditions(tmp_path):
    path = _omnibus_view(tmp_path / "f.de.h5ad", entity_type="site")

    assert not validate_anndata_view(path, "de").issues
    matrix = _build_differential_adata(tmp_path, "f")
    assert np.isnan(matrix["any_dose", "P1"].X[0, 0])
    assert matrix["any_dose", "P1"].layers["pvals_adj"][0, 0] == pytest.approx(0.02)


def test_de_validation_reports_an_unknown_entity_type(tmp_path):
    result = validate_anndata_view(_omnibus_view(tmp_path / "e.de.h5ad", entity_type="transcript"), "de")

    assert any("entity_type" in issue.message for issue in result.errors)
