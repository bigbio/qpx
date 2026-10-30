"""Tests for QPX views: ProteinView, PeptideView, IdentificationSummaryView."""

import pandas as pd
import pytest

from qpx.core.convert import QueryResult
from qpx.dataset import Dataset
from qpx.views.api import IdentificationSummaryView, PeptideView, ProteinView

# ---------------------------------------------------------------------------
# ProteinView tests
# ---------------------------------------------------------------------------


class TestProteinView:
    """Test ProteinView joins PG + Run for protein intensity per sample."""

    def test_intensity_returns_queryresult(self, dataset_dir):
        """ProteinView.intensity() returns a QueryResult."""
        with Dataset(dataset_dir) as ds:
            view = ProteinView(ds)
            result = view.intensity()
            assert isinstance(result, QueryResult)

    def test_intensity_to_df(self, dataset_dir):
        """intensity() can be materialized as a DataFrame."""
        with Dataset(dataset_dir) as ds:
            view = ProteinView(ds)
            df = view.intensity().to_df()
            assert isinstance(df, pd.DataFrame)
            assert len(df) > 0

    def test_intensity_has_expected_columns(self, dataset_dir):
        """intensity() returns protein_accession, sample_accession, label, intensity."""
        with Dataset(dataset_dir) as ds:
            view = ProteinView(ds)
            df = view.intensity().to_df()
            expected_cols = {
                "protein_accession",
                "sample_accession",
                "label",
                "intensity",
            }
            assert expected_cols.issubset(set(df.columns))

    def test_intensity_joins_correctly(self, dataset_dir):
        """intensity() correctly joins PG intensities with run sample mappings.

        In our test data:
        - PG P12345 in run_01 has intensity TMT126. run_01 maps to SAMPLE_01/TMT126.
        - PG P67890 in run_02 has intensity TMT127. run_02 maps to SAMPLE_02/TMT127.
        So we expect 2 rows from the join with correct sample assignments.
        """
        with Dataset(dataset_dir) as ds:
            view = ProteinView(ds)
            df = view.intensity().to_df()

            # Should have rows for each matching (PG intensity label, run sample label)
            assert len(df) == 2

            # Check P12345 maps to SAMPLE_01
            p12345_rows = df[df["protein_accession"] == "P12345"]
            assert len(p12345_rows) == 1
            assert p12345_rows.iloc[0]["sample_accession"] == "SAMPLE_01"
            assert p12345_rows.iloc[0]["label"] == "TMT126"

            # Check P67890 maps to SAMPLE_02
            p67890_rows = df[df["protein_accession"] == "P67890"]
            assert len(p67890_rows) == 1
            assert p67890_rows.iloc[0]["sample_accession"] == "SAMPLE_02"
            assert p67890_rows.iloc[0]["label"] == "TMT127"

    def test_intensity_q_value_threshold(self, dataset_dir):
        """intensity() filters by global_qvalue threshold."""
        with Dataset(dataset_dir) as ds:
            view = ProteinView(ds)
            # Default threshold is 0.01, and our test data has global_qvalue=0.005
            df_default = view.intensity().to_df()
            assert len(df_default) > 0

            # Very strict threshold should filter out results
            df_strict = view.intensity(q_value_threshold=0.001).to_df()
            assert len(df_strict) == 0

    def test_intensity_has_gene_names(self, dataset_dir):
        """intensity() includes gene_names from PG."""
        with Dataset(dataset_dir) as ds:
            view = ProteinView(ds)
            df = view.intensity().to_df()
            assert "gene_names" in df.columns

    def test_intensity_has_global_qvalue(self, dataset_dir):
        """intensity() includes global_qvalue from PG."""
        with Dataset(dataset_dir) as ds:
            view = ProteinView(ds)
            df = view.intensity().to_df()
            assert "global_qvalue" in df.columns


# ---------------------------------------------------------------------------
# PeptideView tests
# ---------------------------------------------------------------------------


class TestPeptideView:
    """Test PeptideView joins Feature + Run for peptide intensity per sample."""

    def test_intensity_returns_queryresult(self, dataset_dir):
        """PeptideView.intensity() returns a QueryResult."""
        with Dataset(dataset_dir) as ds:
            view = PeptideView(ds)
            result = view.intensity()
            assert isinstance(result, QueryResult)

    def test_intensity_to_df(self, dataset_dir):
        """intensity() can be materialized as a DataFrame."""
        with Dataset(dataset_dir) as ds:
            view = PeptideView(ds)
            df = view.intensity().to_df()
            assert isinstance(df, pd.DataFrame)
            assert len(df) > 0

    def test_intensity_has_expected_columns(self, dataset_dir):
        """intensity() returns sequence, peptidoform, charge, sample_accession, intensity."""
        with Dataset(dataset_dir) as ds:
            view = PeptideView(ds)
            df = view.intensity().to_df()
            expected_cols = {
                "sequence",
                "peptidoform",
                "charge",
                "anchor_protein",
                "sample_accession",
                "label",
                "intensity",
            }
            assert expected_cols.issubset(set(df.columns))

    def test_intensity_joins_correctly(self, dataset_dir):
        """intensity() correctly joins feature intensities with run sample mappings.

        In our test data:
        - Features in run_01 have TMT126. run_01 maps to SAMPLE_01/TMT126. (2 features)
        - Feature in run_02 has TMT127. run_02 maps to SAMPLE_02/TMT127. (1 feature)
        So we expect 3 rows.
        """
        with Dataset(dataset_dir) as ds:
            view = PeptideView(ds)
            df = view.intensity().to_df()
            assert len(df) == 3

            # Check run_01 features map to SAMPLE_01
            sample_01_rows = df[df["sample_accession"] == "SAMPLE_01"]
            assert len(sample_01_rows) == 2
            assert all(sample_01_rows["label"] == "TMT126")

            # Check run_02 feature maps to SAMPLE_02
            sample_02_rows = df[df["sample_accession"] == "SAMPLE_02"]
            assert len(sample_02_rows) == 1
            assert sample_02_rows.iloc[0]["label"] == "TMT127"

    def test_intensity_to_arrow(self, dataset_dir):
        """intensity() can be materialized as Arrow Table."""
        with Dataset(dataset_dir) as ds:
            view = PeptideView(ds)
            import pyarrow as pa

            table = view.intensity().to_arrow()
            assert isinstance(table, pa.Table)
            assert table.num_rows > 0


# ---------------------------------------------------------------------------
# IdentificationSummaryView tests
# ---------------------------------------------------------------------------


class TestIdentificationSummaryView:
    """Test IdentificationSummaryView for per-run summary aggregation."""

    def test_summary_returns_queryresult(self, dataset_dir):
        """summary() returns a QueryResult."""
        with Dataset(dataset_dir) as ds:
            view = IdentificationSummaryView(ds)
            result = view.summary()
            assert isinstance(result, QueryResult)

    def test_summary_to_df(self, dataset_dir):
        """summary() can be materialized as a DataFrame."""
        with Dataset(dataset_dir) as ds:
            view = IdentificationSummaryView(ds)
            df = view.summary().to_df()
            assert isinstance(df, pd.DataFrame)
            assert len(df) > 0

    def test_summary_has_expected_columns(self, dataset_dir):
        """summary() returns run_file_name, n_proteins, n_peptides, n_features."""
        with Dataset(dataset_dir) as ds:
            view = IdentificationSummaryView(ds)
            df = view.summary().to_df()
            expected_cols = {"run_file_name", "n_proteins", "n_peptides", "n_features"}
            assert expected_cols.issubset(set(df.columns))

    def test_summary_aggregates_correctly(self, dataset_dir):
        """summary() aggregates correctly per run.

        In our test data:
        - run_01: 2 features (PEPTIDEK/P12345, ANOTHERK/P12345) -> 1 protein, 2 peptides, 2 features
        - run_02: 1 feature (THIRDPEP/P67890) -> 1 protein, 1 peptide, 1 feature
        """
        with Dataset(dataset_dir) as ds:
            view = IdentificationSummaryView(ds)
            df = view.summary().to_df()
            df = df.sort_values("run_file_name").reset_index(drop=True)

            assert len(df) == 2

            # run_01
            run01 = df[df["run_file_name"] == "run_01"].iloc[0]
            assert run01["n_proteins"] == 1  # Only P12345
            assert run01["n_peptides"] == 2  # PEPTIDEK, ANOTHERK
            assert run01["n_features"] == 2  # 2 rows

            # run_02
            run02 = df[df["run_file_name"] == "run_02"].iloc[0]
            assert run02["n_proteins"] == 1  # Only P67890
            assert run02["n_peptides"] == 1  # THIRDPEP
            assert run02["n_features"] == 1  # 1 row

    def test_summary_to_arrow(self, dataset_dir):
        """summary() can be materialized as Arrow Table."""
        with Dataset(dataset_dir) as ds:
            view = IdentificationSummaryView(ds)
            import pyarrow as pa

            table = view.summary().to_arrow()
            assert isinstance(table, pa.Table)
            assert table.num_rows == 2  # Two runs


# ---------------------------------------------------------------------------
# View integration tests
# ---------------------------------------------------------------------------


def _lfq_run(run_accession, run_file_name, sample_accession):
    """A label-free (single-sample) run record; run.samples[].label == 'LFQ'."""
    return {
        "run_accession": run_accession,
        "run_file_name": run_file_name,
        "samples": [
            {
                "sample_accession": sample_accession,
                "label": "LFQ",
                "biological_replicate": 1,
                "technical_replicate": 1,
            }
        ],
        "fraction": None,
        "instrument": None,
        "enzymes": None,
        "dissociation_method": None,
        "modification_parameters": None,
    }


@pytest.fixture(name="multi_run_dataset_dir")
def _multi_run_dataset_dir(tmp_path):
    """A dataset with a genuine multi-element grouped_runs and DIA-NN LFQ labels.

    - Two fractions (run_A and run_B) from the same label-free sample.
    - One protein group PROT1 whose grouped_runs is BOTH runs, quantified once
      with a single 'LFQ' intensity (the tool-reported LFQ value for the group).
    - Features carry the DIA-NN label 'raw' (NOT the run's 'LFQ' channel), to
      exercise the PeptideView label-mismatch fix.
    """
    from qpx.writers import DatasetWriter, FeatureWriter, PgWriter, RunWriter, SampleWriter
    from tests.conftest import make_dataset_record, make_feature_record, make_pg_record, make_sample_record

    ds_dir = tmp_path / "multi_run_dataset"
    ds_dir.mkdir()

    # Features: DIA-NN emits label 'raw'; one feature per run.
    feature_records = [
        make_feature_record(
            sequence="PEPTIDEK",
            anchor_protein="PROT1",
            run_file_name="run_A",
            intensities=[{"label": "raw", "intensity": 1000.0}],
        ),
        make_feature_record(
            sequence="PEPTIDEK",
            anchor_protein="PROT1",
            run_file_name="run_B",
            scan=[2001],
            intensities=[{"label": "raw", "intensity": 2000.0}],
        ),
    ]
    with FeatureWriter(ds_dir / "exp.feature.parquet") as w:
        w.write_batch(feature_records)

    # PG: one row keyed by a TWO-element grouped_runs list.
    pg_rec = make_pg_record(anchor_protein="PROT1")
    pg_rec["grouped_runs"] = ["run_A", "run_B"]
    pg_rec["intensities"] = [{"label": "LFQ", "intensity": 5000.0}]
    with PgWriter(ds_dir / "exp.pg.parquet") as w:
        w.write_batch([pg_rec])

    with SampleWriter(ds_dir / "exp.sample.parquet") as w:
        w.write_batch([make_sample_record("SAMPLE_A")])

    with RunWriter(ds_dir / "exp.run.parquet") as w:
        w.write_batch(
            [
                _lfq_run("assay_A", "run_A", "SAMPLE_A"),
                _lfq_run("assay_B", "run_B", "SAMPLE_A"),
            ]
        )

    with DatasetWriter(ds_dir / "exp.dataset.parquet") as w:
        w.write_batch([make_dataset_record()])

    return ds_dir


class TestMultiRunGroupedRuns:
    """Integration tests for a genuine multi-element grouped_runs group."""

    def test_protein_view_deduplicates_fraction_sample(self, multi_run_dataset_dir):
        """One grouped PG quantity is returned once for its shared sample."""
        with Dataset(multi_run_dataset_dir) as ds:
            df = ProteinView(ds).intensity().to_df()

        prot1 = df[df["protein_accession"] == "PROT1"]
        assert list(prot1["sample_accession"]) == ["SAMPLE_A"]
        assert list(prot1["intensity"]) == [5000.0]

    def test_dataset_intensity_does_not_double_fraction_quantity(self, multi_run_dataset_dir):
        """Dataset helpers must not sum one PG value once per fraction."""
        with Dataset(multi_run_dataset_dir) as ds:
            abundance = ds.intensity("protein").to_df()
            matrix = ds.design_matrix("protein")

        assert len(abundance) == 1
        assert abundance.iloc[0]["intensity"] == 5000.0
        assert matrix.loc["SAMPLE_A", "PROT1"] == 5000.0

    def test_mudata_proteins_observe_the_quantification_unit(self, multi_run_dataset_dir):
        """The muData proteins modality has one observation per grouped_runs set.

        It used to have one per member run, each holding a copy of the unit's
        value, so summing over observations doubled this two-fraction group.
        """
        from qpx.mudata import build_mudata

        with Dataset(multi_run_dataset_dir) as ds:
            mdata = build_mudata(ds, modalities=["proteins"])

        proteins = mdata.mod["proteins"]
        assert list(proteins.var.index) == ["PROT1"]
        assert list(proteins.obs.index) == ["run_A;run_B"]
        assert list(proteins.obs["sample_accession"]) == ["SAMPLE_A"]
        assert list(proteins.obs["run_accession"]) == ["assay_A;assay_B"]
        assert proteins.X.toarray().sum() == 5000.0

    def test_invariant_check_passes(self, multi_run_dataset_dir):
        """Both grouped_runs values resolve to the same sample without errors."""
        with Dataset(multi_run_dataset_dir) as ds:
            results = ds.validate(structures=["pg"])
        grouped_run_issues = [
            issue for issue in results["pg"].issues if issue.check in {"dangling_grouped_run", "ambiguous_grouped_run_sample"}
        ]
        assert grouped_run_issues == []

    def test_peptide_view_resolves_diann_lfq_label(self, multi_run_dataset_dir):
        """PeptideView joins DIA-NN 'raw' feature labels to single-sample 'LFQ' runs."""
        with Dataset(multi_run_dataset_dir) as ds:
            df = PeptideView(ds).intensity().to_df()

        assert len(df) == 2  # would be 0 rows before the label-mismatch fix
        assert set(df["sample_accession"]) == {"SAMPLE_A"}
        assert set(df["label"]) == {"LFQ"}

    def test_peptide_helpers_resolve_diann_lfq_label(self, multi_run_dataset_dir):
        """All public peptide helpers apply the same single-sample LFQ mapping."""
        with Dataset(multi_run_dataset_dir) as ds:
            abundance = ds.intensity("peptide").to_df()
            quantified = ds.feature.for_quantification(ds).to_df()

        assert len(abundance) == 1
        assert abundance.iloc[0]["sample_accession"] == "SAMPLE_A"
        assert abundance.iloc[0]["intensity"] == 3000.0
        assert len(quantified) == 2
        assert set(quantified["sample_accession"]) == {"SAMPLE_A"}
        assert set(quantified["label"]) == {"LFQ"}


class TestViewIntegration:
    """Integration tests using views with Dataset."""

    def test_protein_and_peptide_views_share_engine(self, dataset_dir):
        """Multiple views can be created from the same dataset, sharing the engine."""
        with Dataset(dataset_dir) as ds:
            prot_view = ProteinView(ds)
            pep_view = PeptideView(ds)
            id_view = IdentificationSummaryView(ds)

            # All views should produce results
            prot_df = prot_view.intensity().to_df()
            pep_df = pep_view.intensity().to_df()
            id_df = id_view.summary().to_df()

            assert len(prot_df) > 0
            assert len(pep_df) > 0
            assert len(id_df) > 0

    def test_views_with_partial_dataset(self, dataset_dir):
        """Views work when dataset is partially loaded (as long as needed structures exist)."""
        with Dataset(dataset_dir, structures=["feature", "run"]) as ds:
            pep_view = PeptideView(ds)
            df = pep_view.intensity().to_df()
            assert len(df) == 3

    def test_protein_view_results_match_manual_sql(self, dataset_dir):
        """ProteinView results should match an equivalent manual SQL query."""
        with Dataset(dataset_dir) as ds:
            view = ProteinView(ds)
            view_df = view.intensity().to_df()

            manual_result = ds.sql("""
                SELECT pg.anchor_protein AS protein_accession,
                       rs.sample_accession,
                       pg.label,
                       pg.intensity,
                       pg.global_qvalue,
                       pg.gg_names AS gene_names
                FROM pg,
                     run r,
                     UNNEST(r.samples) AS _t1(rs)
                WHERE list_contains(pg.grouped_runs, r.run_file_name)
                  AND pg.label = rs.label
                  AND (pg.global_qvalue IS NULL OR pg.global_qvalue <= 0.01)
            """)
            manual_df = manual_result.to_df()

            assert len(view_df) == len(manual_df)
            # Compare protein sets
            assert set(view_df["protein_accession"]) == set(manual_df["protein_accession"])
