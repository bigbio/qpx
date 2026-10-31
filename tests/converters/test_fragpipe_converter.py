"""FragPipe PG converter tests focused on grouped_runs expansion.

FragPipe ``combined_protein.tsv`` reports intensities per *experiment*. Each
experiment can aggregate several raw files, so ``grouped_runs`` must expand to
the experiment's member ``run_file_name`` values (via ``experiment_to_runs``,
``experiment_annotation.tsv`` or the SDRF); an experiment none of them resolves
raises.
"""

from __future__ import annotations

import pyarrow.parquet as pq
import pytest

from qpx.converters.fragpipe.feature_adapter import FragPipeFeatureAdapter
from qpx.converters.fragpipe.pg_adapter import FragPipePgAdapter
from qpx.converters.fragpipe.psm_adapter import FragPipePsmAdapter
from qpx.converters.ptm import PROTON_MASS, compute_precursor_mz

_COMBINED_PROTEIN_TSV = (
    "Protein\tGene\tDescription\t"
    "exp1 Total Intensity\texp1 MaxLFQ Intensity\texp1 Spectral Count\texp1 Unique Spectral Count\t"
    "exp2 Total Intensity\texp2 MaxLFQ Intensity\texp2 Spectral Count\texp2 Unique Spectral Count\t"
    "Combined Total Peptides\tCombined Unique Peptides\n"
    "P12345\tGENE1\tsome protein\t"
    "1000\t900\t10\t8\t"
    "2000\t1800\t20\t16\t"
    "5\t3\n"
)

_EXPERIMENT_ANNOTATION_TSV = "file\tsample\n/data/run_B.raw\texp1\n/data/run_A.mzML\texp1\nC:\\\\data\\\\run_C.RAW\texp2\n"


def _write_input(tmp_path):
    path = tmp_path / "combined_protein.tsv"
    path.write_text(_COMBINED_PROTEIN_TSV)
    return path


def _grouped_runs_by_intensity(table):
    """Return {intensity: grouped_runs} so records are identifiable per experiment.

    pg is flattened since 1.1: scalar ``intensity`` column, one row per label.
    """
    df = table.to_pandas()
    out = {}
    for _, row in df.iterrows():
        out[float(row["intensity"])] = list(row["grouped_runs"])
    return out


def test_grouped_runs_expand_to_member_run_files(tmp_path):
    """With an experiment->runs mapping, grouped_runs holds real run file names."""
    protein_path = _write_input(tmp_path)
    out_path = tmp_path / "fragpipe.pg.parquet"

    with FragPipePgAdapter() as adapter:
        adapter.convert(
            protein_path=str(protein_path),
            output_path=str(out_path),
            experiment_to_runs={"exp1": ["run_A", "run_B"], "exp2": ["run_C"]},
        )

    table = pq.read_table(str(out_path))
    grouped = _grouped_runs_by_intensity(table)
    assert grouped[1000.0] == ["run_A", "run_B"]  # exp1 expanded to its two fractions
    assert grouped[2000.0] == ["run_C"]


def test_no_mapping_raises_instead_of_dangling_token(tmp_path):
    """With no results mapping and no SDRF, conversion raises rather than emitting
    a dangling [experiment] token (bigbio/qpx#220 review item #5)."""
    import pytest

    protein_path = _write_input(tmp_path)
    out_path = tmp_path / "fragpipe.pg.parquet"

    with FragPipePgAdapter() as adapter:
        with pytest.raises(ValueError, match="grouped_runs.*experiment|experiment.*mapping"):
            adapter.convert(protein_path=str(protein_path), output_path=str(out_path))


def test_sdrf_fallback_maps_experiment_to_runs(tmp_path):
    """When no experiment_annotation is given, the SDRF 'source name' grouping is
    the fallback: each experiment resolves to that sample's raw files."""
    protein_path = _write_input(tmp_path)
    sdrf_path = tmp_path / "test.sdrf.tsv"
    sdrf_path.write_text("source name\tcomment[data file]\nexp1\trun_A1.raw\nexp1\trun_A2.raw\nexp2\trun_B1.raw\n")
    out_path = tmp_path / "fragpipe.pg.parquet"

    with FragPipePgAdapter() as adapter:
        adapter.convert(
            protein_path=str(protein_path),
            output_path=str(out_path),
            sdrf_path=str(sdrf_path),
        )

    grouped = _grouped_runs_by_intensity(pq.read_table(str(out_path)))
    assert grouped[1000.0] == ["run_A1", "run_A2"]
    assert grouped[2000.0] == ["run_B1"]


def test_fragpipe_orchestrator_passes_sdrf_to_pg_adapter(tmp_path):
    """The public FragPipe converter must expose the PG adapter's SDRF fallback,
    not require callers to invoke the adapter directly."""
    from qpx.converters.fragpipe.converter import FragPipeConverter

    protein_path = _write_input(tmp_path)
    sdrf_path = tmp_path / "test.sdrf.tsv"
    sdrf_path.write_text(
        "source name\tcharacteristics[organism]\tcharacteristics[organism part]\t"
        "comment[data file]\tcomment[label]\n"
        "exp1\tHomo sapiens\tliver\trun_A1.raw\tlabel free sample\n"
        "exp1\tHomo sapiens\tliver\trun_A2.raw\tlabel free sample\n"
        "exp2\tHomo sapiens\tliver\trun_B1.raw\tlabel free sample\n"
    )
    output = tmp_path / "converted"

    FragPipeConverter(output_directory=output).convert(
        pg_file=protein_path,
        sdrf_file=sdrf_path,
        output_prefix="fragpipe",
    )

    grouped = _grouped_runs_by_intensity(pq.read_table(output / "fragpipe.pg.parquet"))
    assert grouped[1000.0] == ["run_A1", "run_A2"]
    assert grouped[2000.0] == ["run_B1"]


def test_experiment_annotation_expands_and_normalizes_member_runs(tmp_path):
    """The official FragPipe annotation file drives grouped_runs directly."""
    protein_path = _write_input(tmp_path)
    annotation_path = tmp_path / "experiment_annotation.tsv"
    annotation_path.write_text(_EXPERIMENT_ANNOTATION_TSV)
    out_path = tmp_path / "fragpipe.pg.parquet"

    with FragPipePgAdapter() as adapter:
        adapter.convert(
            protein_path=str(protein_path),
            output_path=str(out_path),
            experiment_annotation_path=str(annotation_path),
        )

    grouped = _grouped_runs_by_intensity(pq.read_table(str(out_path)))
    assert grouped[1000.0] == ["run_A", "run_B"]
    assert grouped[2000.0] == ["run_C"]


def test_experiment_annotation_and_explicit_mapping_are_mutually_exclusive(
    tmp_path,
):
    """Ambiguous mapping inputs fail before writing output."""
    protein_path = _write_input(tmp_path)
    annotation_path = tmp_path / "experiment_annotation.tsv"
    annotation_path.write_text(_EXPERIMENT_ANNOTATION_TSV)

    with FragPipePgAdapter() as adapter:
        with pytest.raises(ValueError, match="either experiment_to_runs"):
            adapter.convert(
                protein_path=str(protein_path),
                output_path=str(tmp_path / "fragpipe.pg.parquet"),
                experiment_to_runs={"exp1": ["run_A"]},
                experiment_annotation_path=str(annotation_path),
            )


# ---------------------------------------------------------------------------
# Regression tests for bigbio/qpx#246 and #247
# ---------------------------------------------------------------------------


def _read_pg(tmp_path, protein_tsv, **convert_kwargs):
    """Write a combined_protein.tsv, convert, and return the pg DataFrame."""
    protein_path = tmp_path / "combined_protein.tsv"
    protein_path.write_text(protein_tsv)
    out_path = tmp_path / "fragpipe.pg.parquet"
    convert_kwargs.setdefault("experiment_to_runs", {"expA": ["run_A"]})
    with FragPipePgAdapter() as adapter:
        adapter.convert(protein_path=str(protein_path), output_path=str(out_path), **convert_kwargs)
    return pq.read_table(str(out_path)).to_pandas()


def test_pg_qvalue_prefers_real_fdr_over_probability(tmp_path):
    """#246: with a real FDR column present, pg_qvalue must come from it, not from
    the always-present 'Protein Probability' confidence score."""
    tsv = (
        "Protein\tGene\tDescription\tProtein Probability\tProtein FDR\t"
        "expA Total Intensity\tCombined Total Peptides\tCombined Unique Peptides\n"
        "sp|P1|A_HUMAN\tG1\tprot1\t0.9999\t0.001\t1000\t5\t3\n"
    )
    df = _read_pg(tmp_path, tsv)
    # The bug stored 0.9999 (the raw probability) as a q-value.
    assert df["pg_qvalue"].iloc[0] == pytest.approx(0.001)


def test_pg_qvalue_probability_only_is_inverted(tmp_path):
    """#246: when only 'Protein Probability' is available it is converted to an
    FDR-like value (1 - probability), never stored raw."""
    tsv = (
        "Protein\tGene\tDescription\tProtein Probability\t"
        "expA Total Intensity\tCombined Total Peptides\tCombined Unique Peptides\n"
        "sp|P1|A_HUMAN\tG1\tprot1\t0.9999\t1000\t5\t3\n"
    )
    df = _read_pg(tmp_path, tsv)
    assert df["pg_qvalue"].iloc[0] == pytest.approx(1.0 - 0.9999)


def test_maxlfq_only_experiment_emits_pg(tmp_path):
    """#247 (M): an experiment detected only from a '<exp> MaxLFQ Intensity' column
    (no Total Intensity) must still produce a pg record, using MaxLFQ as intensity."""
    tsv = (
        "Protein\tGene\tDescription\tProtein FDR\t"
        "expA MaxLFQ Intensity\tCombined Total Peptides\tCombined Unique Peptides\n"
        "sp|P1|A_HUMAN\tG1\tprot1\t0.001\t1234\t5\t3\n"
    )
    df = _read_pg(tmp_path, tsv)
    assert len(df) == 1
    assert df["intensity"].iloc[0] == pytest.approx(1234.0)


def test_maxlfq_falls_back_when_total_intensity_is_zero(tmp_path):
    """A zero Total Intensity with a positive MaxLFQ still yields a record."""
    tsv = (
        "Protein\tGene\tDescription\tProtein FDR\t"
        "expA Total Intensity\texpA MaxLFQ Intensity\tCombined Total Peptides\tCombined Unique Peptides\n"
        "sp|P1|A_HUMAN\tG1\tprot1\t0.001\t0\t1234\t5\t3\n"
    )
    df = _read_pg(tmp_path, tsv)
    assert len(df) == 1
    assert df["intensity"].iloc[0] == pytest.approx(1234.0)


def test_pg_decoy_detected_from_raw_prefix(tmp_path):
    """#247 (H): decoy prefixes (incl. REV_) are detected on the raw accession, not
    after parse_uniprot_id strips them."""
    tsv = (
        "Protein\tGene\tDescription\tProtein FDR\t"
        "expA Total Intensity\tCombined Total Peptides\tCombined Unique Peptides\n"
        "sp|P1|A_HUMAN\tG1\ttarget\t0.001\t1000\t5\t3\n"
        "REV_sp|P2|B_HUMAN\tG2\tdecoy\t0.02\t50\t2\t1\n"
    )
    df = _read_pg(tmp_path, tsv)
    decoy_flag = {row["pg_accessions"][0]: bool(row["is_decoy"]) for _, row in df.iterrows()}
    assert decoy_flag["P1"] is False
    assert decoy_flag["REV_P2"] is True


def test_pg_maps_fragpipe_22_quantities_counts_and_annotations(tmp_path):
    """A FragPipe 22 combined_protein.tsv: `<exp> Intensity` is the primary quantity,
    spectral counts and MaxLFQ are additional intensities, the ProteinProphet
    probabilities are scores, length and existence are cv_params, and the
    combined peptide count is not copied into the experiment's counts."""
    tsv = (
        "Protein\tProtein ID\tGene\tProtein Length\tProtein Existence\tDescription\tProtein Probability\t"
        "Top Peptide Probability\tCombined Total Peptides\tCombined Spectral Count\texpA Spectral Count\t"
        "expA Unique Spectral Count\texpA Total Spectral Count\texpA Intensity\texpA MaxLFQ Intensity\t"
        "Indistinguishable Proteins\n"
        "sp|P1|A_HUMAN\tP1\tG1\t495\t1:Experimental evidence at protein level\tprot1\t0.999\t0.998\t"
        "12\t40\t5\t3\t7\t1000\t900\tsp|P2|B_HUMAN\n"
    )
    row = _read_pg(tmp_path, tsv).iloc[0]
    assert list(row["pg_accessions"]) == ["P1", "P2"]
    assert row["intensity"] == pytest.approx(1000.0)
    assert {q["intensity_name"]: q["intensity_value"] for q in row["additional_intensities"][0]["intensities"]} == {
        "maxlfq": 900.0,
        "spectral_count": 5.0,
        "unique_spectral_count": 3.0,
        "total_spectral_count": 7.0,
    }
    assert {s["score_name"]: s["score_value"] for s in row["additional_scores"]} == pytest.approx(
        {"protein_probability": 0.999, "top_peptide_probability": 0.998}
    )
    assert [tuple(cv.values()) for cv in row["cv_params"]] == [("sequence_length", "495"), ("protein_existence", "1")]
    assert row["peptide_counts"] is None
    assert row["feature_counts"] is None


def test_spectral_count_only_report_writes_pg(tmp_path):
    """A combined_protein.tsv without intensities still names its experiments by
    their spectral counts, and each group is written with a null intensity."""
    tsv = (
        "Protein\tCombined Total Peptides\texpA Spectral Count\texpA Unique Spectral Count\t"
        "expA Total Spectral Count\texpA Unique Peptides\texpA Total Peptides\n"
        "sp|P1|A_HUMAN\t12\t5\t3\t7\t2\t4\n"
    )
    df = _read_pg(tmp_path, tsv)
    assert len(df) == 1
    assert df["intensity"].isna().all()
    assert df.iloc[0]["peptide_counts"] == {"unique_sequences": 2, "total_sequences": 4}


# ---- Feature adapter (combined_ion.tsv + psm.tsv) --------------------------


def _write_feature_inputs(tmp_path, *, experiment, spectrum_source, protein="sp|P12345|PROT_HUMAN"):
    ion = tmp_path / "combined_ion.tsv"
    ion.write_text(
        "Peptide Sequence\tModified Sequence\tProtein\tGene\tM/Z\tCharge\tAssigned Modifications\t"
        f"{experiment} Intensity\t{experiment} MaxLFQ Intensity\n"
        f"PEPTIDEK\tPEPTIDEK\t{protein}\tGENE1\t500.0\t2\t\t1000\t900\n"
    )
    psm = tmp_path / "psm.tsv"
    psm.write_text(
        "Spectrum\tPeptide\tCharge\tAssigned Modifications\t"
        "Calibrated Observed M/Z\tCalculated M/Z\tPeptideProphet Probability\tProtein\n"
        f"{spectrum_source}.00500.00500.2\tPEPTIDEK\t2\t\t500.001\t500.000\t0.99\trev_{protein}\n"
    )
    return ion, psm


def test_feature_psm_enrichment_uses_experiment_annotation(tmp_path):
    """A one-run experiment resolves to the exact raw-file PSM through FragPipe metadata."""
    from qpx.converters.fragpipe.converter import FragPipeConverter

    ion, psm = _write_feature_inputs(tmp_path, experiment="mc38_p_1", spectrum_source="mc38_p1")
    annotation = tmp_path / "experiment_annotation.tsv"
    annotation.write_text("file\tsample\nD:\\data\\mc38_p1.mzML\tmc38_p_1\n")
    output = tmp_path / "converted"
    FragPipeConverter(output_directory=output).convert(
        ion_file=ion,
        psm_file=psm,
        experiment_annotation_file=annotation,
        output_prefix="fragpipe",
    )
    out = output / "fragpipe.feature.parquet"
    df = pq.read_table(str(out)).to_pandas()
    assert len(df) == 1
    rec = df.iloc[0]
    assert rec["posterior_error_probability"] == pytest.approx(1.0 - 0.99)
    assert list(rec["scan"]) == [500]
    theoretical = compute_precursor_mz("PEPTIDEK", 2)
    assert rec["mass_error_ppm"] == pytest.approx(1e6 * (500.001 - theoretical) / theoretical, rel=1e-4)
    assert bool(rec["is_decoy"]) is True


def test_feature_takes_the_best_psm_of_its_run(tmp_path):
    """Of several PSMs of a precursor in the run, the feature links the one with the
    highest probability, not whichever DuckDB returns first."""
    ion, psm = _write_feature_inputs(tmp_path, experiment="sampleA", spectrum_source="run_01")
    with psm.open("a") as handle:
        handle.write("run_01.00300.00300.2\tPEPTIDEK\t2\t\t500.002\t500.000\t0.80\tsp|P12345|PROT_HUMAN\n")
        handle.write("run_01.00700.00700.2\tPEPTIDEK\t2\t\t500.003\t500.000\t0.999\tsp|P12345|PROT_HUMAN\n")
    out = tmp_path / "feature.parquet"
    with FragPipeFeatureAdapter() as adapter:
        adapter.convert(
            feature_path=str(ion), output_path=str(out), psm_path=str(psm), experiment_to_runs={"sampleA": ["run_01"]}
        )
    record = pq.read_table(out).to_pylist()[0]
    assert record["scan"] == [700]
    assert record["posterior_error_probability"] == pytest.approx(0.001)


def test_feature_psm_fallback_does_not_cross_multiple_raw_files(tmp_path):
    """A run-agnostic precursor match is unsafe when the PSM occurs in more
    than one raw file; leave enrichment empty instead of taking the first run."""
    ion, psm = _write_feature_inputs(tmp_path, experiment="sampleA", spectrum_source="run_01")
    with psm.open("a") as handle:
        handle.write("run_02.00600.00600.2\tPEPTIDEK\t2\t\t501.001\t501.000\t0.90\tsp|P12345|PROT_HUMAN\n")

    out = tmp_path / "ambiguous.feature.parquet"
    with FragPipeFeatureAdapter() as adapter:
        adapter.convert(
            feature_path=str(ion),
            output_path=str(out),
            psm_path=str(psm),
            experiment_to_runs={"sampleA": ["run_01", "run_02"]},
        )

    record = pq.read_table(out).to_pylist()[0]
    assert record["scan"] == []
    assert record["posterior_error_probability"] is None
    assert record["mass_error_ppm"] is None


def test_feature_psm_enrichment_does_not_cross_experiments(tmp_path):
    """One run's PSM must not enrich another experiment with the same precursor."""
    ion = tmp_path / "combined_ion.tsv"
    ion.write_text(
        "Peptide Sequence\tModified Sequence\tProtein\tGene\tM/Z\tCharge\tAssigned Modifications\t"
        "sampleA Intensity\tsampleB Intensity\n"
        "PEPTIDEK\tPEPTIDEK\tsp|P12345|PROT_HUMAN\tGENE1\t500.0\t2\t\t1000\t2000\n"
    )
    psm = tmp_path / "psm.tsv"
    psm.write_text(
        "Spectrum\tPeptide\tCharge\tAssigned Modifications\t"
        "Calibrated Observed M/Z\tCalculated M/Z\tPeptideProphet Probability\tProtein\n"
        "run_A.00500.00500.2\tPEPTIDEK\t2\t\t500.001\t500.000\t0.99\tsp|P12345|PROT_HUMAN\n"
    )
    out = tmp_path / "feature.parquet"
    with FragPipeFeatureAdapter() as adapter:
        adapter.convert(
            feature_path=str(ion),
            output_path=str(out),
            psm_path=str(psm),
            experiment_to_runs={"sampleA": ["run_A"], "sampleB": ["run_B"]},
        )

    rows = {row["run_file_name"]: row for row in pq.read_table(out).to_pylist()}
    assert rows["run_A"]["scan"] == [500]
    assert rows["run_B"]["scan"] == []
    assert rows["run_B"]["posterior_error_probability"] is None


def test_feature_decoy_from_raw_accession_without_psm(tmp_path):
    """#247 (H): without a psm.tsv, is_decoy falls back to the raw protein prefix,
    which must survive parse_uniprot_id (rev_sp|..| -> still a decoy)."""
    ion = tmp_path / "combined_ion.tsv"
    ion.write_text(
        "Peptide Sequence\tModified Sequence\tProtein\tGene\tM/Z\tCharge\tAssigned Modifications\t"
        "sampleA Intensity\n"
        "PEPTIDEK\tPEPTIDEK\trev_sp|P12345|PROT_HUMAN\tGENE1\t500.0\t2\t\t1000\n"
    )
    out = tmp_path / "feature.parquet"
    with FragPipeFeatureAdapter() as adapter:
        adapter.convert(feature_path=str(ion), output_path=str(out), experiment_to_runs={"sampleA": ["run_A"]})
    df = pq.read_table(str(out)).to_pandas()
    assert bool(df.iloc[0]["is_decoy"]) is True


def test_feature_rows_are_keyed_by_the_experiment_run(tmp_path):
    """A feature names the raw run of its experiment, which run.parquet and the pg
    grouped_runs hold, and the run that identified it only when that is known."""
    ion = tmp_path / "combined_ion.tsv"
    ion.write_text(
        "Peptide Sequence\tProtein\tM/Z\tCharge\tAssigned Modifications\t"
        "expA Intensity\texpA Match Type\texpB Intensity\texpB Match Type\tfrac Intensity\tfrac Match Type\n"
        "PEPTIDEK\tsp|P1|A_HUMAN\t500.0\t2\t\t1000\tMS/MS\t2000\tMBR\t3000\tMS/MS\n"
    )
    psm = tmp_path / "psm.tsv"
    psm.write_text(
        "Spectrum\tPeptide\tCharge\tAssigned Modifications\tCalibrated Observed M/Z\tCalculated M/Z\tProbability\tProtein\n"
        "run_A.00500.00500.2\tPEPTIDEK\t2\t\t500.001\t500.000\t0.99\tsp|P1|A_HUMAN\n"
    )
    out = tmp_path / "feature.parquet"
    with FragPipeFeatureAdapter() as adapter:
        adapter.convert(
            feature_path=str(ion),
            output_path=str(out),
            psm_path=str(psm),
            experiment_to_runs={"expA": ["run_A"], "expB": ["run_B"], "frac": ["run_F2", "run_F1"]},
        )

    rows = {row["quantification_unit_id"]: row for row in pq.read_table(out).to_pylist()}
    assert {unit: (row["run_file_name"], row["id_run_file_name"]) for unit, row in rows.items()} == {
        "expA": ("run_A", "run_A"),
        "expB": ("run_B", None),  # transferred by match-between-runs
        "frac": ("run_F1", None),  # identified in one of several fractions
    }
    assert rows["expA"]["posterior_error_probability"] == pytest.approx(0.01)
    # combined_ion M/Z is the calculated m/z; only a PSM gives the observed one.
    assert rows["expA"]["observed_mz"] == pytest.approx(500.001)
    assert rows["expB"]["observed_mz"] is None
    assert rows["expB"]["calculated_mz"] == pytest.approx(500.0)


def test_feature_experiment_without_runs_raises(tmp_path):
    """An experiment no mapping resolves cannot name a run, so conversion fails."""
    ion, _ = _write_feature_inputs(tmp_path, experiment="sampleA", spectrum_source="run_01")
    with FragPipeFeatureAdapter() as adapter:
        with pytest.raises(ValueError, match="No raw runs for FragPipe experiment.*sampleA"):
            adapter.convert(feature_path=str(ion), output_path=str(tmp_path / "feature.parquet"))


def test_features_carry_their_group_and_mapped_proteins(tmp_path):
    """A feature lists its group's members as pg writes them, is unique only when the
    peptide maps to one protein, and the pepmap lists every protein it maps to."""
    from qpx.converters.fragpipe.converter import FragPipeConverter

    ion = tmp_path / "combined_ion.tsv"
    ion.write_text(
        "Peptide Sequence\tProtein\tMapped Proteins\tM/Z\tCharge\tAssigned Modifications\texpA Intensity\n"
        "SHAREDK\tsp|P1|A_HUMAN\tsp|P2|B_HUMAN\t500.0\t2\t\t1000\n"
        "OTHERK\tsp|P1|A_HUMAN\tsp|P3|C_HUMAN\t400.0\t2\t\t900\n"
        "ONLYK\tsp|P3|C_HUMAN\t\t300.0\t2\t\t800\n"
    )
    protein = tmp_path / "combined_protein.tsv"
    protein.write_text(
        "Protein\tProtein ID\tIndistinguishable Proteins\texpA Total Intensity\n"
        "sp|P1|A_HUMAN\tP1\tsp|P2|B_HUMAN\t5000\n"
        "sp|P3|C_HUMAN\tP3\t\t800\n"
    )
    out = tmp_path / "out"
    FragPipeConverter(output_directory=out).convert(
        ion_file=ion, pg_file=protein, experiment_to_runs={"expA": ["run_A"]}, output_prefix="fp"
    )

    features = {row["sequence"]: row for row in pq.read_table(out / "fp.feature.parquet").to_pylist()}
    groups = sorted(row["pg_accessions"] for row in pq.read_table(out / "fp.pg.parquet").to_pylist())
    assert groups == [["P1", "P2"], ["P3"]]
    assert {seq: [p["accession"] for p in row["pg_accessions"]] for seq, row in features.items()} == {
        "SHAREDK": ["P1", "P2"],
        "OTHERK": ["P1", "P2"],
        "ONLYK": ["P3"],
    }
    assert {seq: row["unique"] for seq, row in features.items()} == {"SHAREDK": False, "OTHERK": False, "ONLYK": True}
    pepmap = pq.read_table(out / "fp.pepmap.parquet").to_pylist()
    assert {row["sequence"]: [p["accession"] for p in row["pg_accessions"]] for row in pepmap} == {
        "SHAREDK": ["P1", "P2"],
        "OTHERK": ["P1", "P3"],
        "ONLYK": ["P3"],
    }


def test_isobaric_sdrf_is_refused(tmp_path):
    """combined_ion/combined_protein hold MS1 intensities; for a TMT SDRF they must not
    be written as the label-free quantity."""
    ion, _ = _write_feature_inputs(tmp_path, experiment="plex1", spectrum_source="run_01")
    protein = _write_input(tmp_path)
    sdrf = tmp_path / "sdrf.tsv"
    sdrf.write_text("source name\tcomment[data file]\tcomment[label]\ns1\trun_01.raw\tTMT126\ns2\trun_01.raw\tTMT127N\n")
    with FragPipeFeatureAdapter() as adapter:
        with pytest.raises(ValueError, match="TMT reporter-ion quantities are not converted"):
            adapter.convert(feature_path=str(ion), output_path=str(tmp_path / "f.parquet"), sdrf_path=str(sdrf))
    with FragPipePgAdapter() as adapter:
        with pytest.raises(ValueError, match="TMT reporter-ion quantities are not converted"):
            adapter.convert(protein_path=str(protein), output_path=str(tmp_path / "p.parquet"), sdrf_path=str(sdrf))


def test_one_experiment_ion_table_is_refused(tmp_path):
    """A per-experiment ion.tsv names no experiment and carries TMT channels; it wrote no
    feature file while the conversion reported one."""
    ion = tmp_path / "ion.tsv"
    ion.write_text(
        "Peptide Sequence\tM/Z\tCharge\tIntensity\tAssigned Modifications\tProtein\tMapped Proteins\t"
        "plex_TMT_126\tplex_TMT_127N\n"
        "PEPTIDEK\t500.0\t2\t1000\tN-term(229.1629)\tsp|P1|A_HUMAN\t\t10\t20\n"
    )
    with FragPipeFeatureAdapter() as adapter:
        with pytest.raises(ValueError, match="no '<experiment> Intensity' columns"):
            adapter.convert(feature_path=str(ion), output_path=str(tmp_path / "f.parquet"))


# ---- PSM adapter (psm.tsv) -------------------------------------------------


def test_psm_lists_mapped_proteins_and_reads_probability(tmp_path):
    """protein_accessions holds every protein the peptide maps to; FragPipe 22 names the
    PeptideProphet probability "Probability"; a missing observed m/z stays null."""
    psm = tmp_path / "psm.tsv"
    psm.write_text(
        "Spectrum\tPeptide\tCharge\tAssigned Modifications\tCalculated M/Z\tProbability\tProtein\tMapped Proteins\n"
        "run_01.00500.00500.2\tPEPTIDEK\t2\t\t500.000\t0.99\tsp|P1|A_HUMAN\tsp|P2|B_HUMAN, tr|Q9|C_HUMAN\n"
    )
    out = tmp_path / "psm.parquet"
    with FragPipePsmAdapter() as adapter:
        adapter.convert(psm_path=str(psm), output_path=str(out))
    row = pq.read_table(out).to_pylist()[0]
    assert row["protein_accessions"] == ["P1", "P2", "Q9"]
    assert row["posterior_error_probability"] == pytest.approx(0.01)
    assert row["observed_mz"] is None
    assert row["mass_error_ppm"] is None


def test_psm_purity_is_the_parent_ion_fraction(tmp_path):
    """FragPipe Purity, the precursor's share of the isolation window, is kept as a score."""
    psm = tmp_path / "psm.tsv"
    psm.write_text(
        "Spectrum\tPeptide\tCharge\tAssigned Modifications\tCalculated M/Z\tPurity\tProtein\n"
        "run_01.00500.00500.2\tPEPTIDEK\t2\t\t500.000\t0.87\tsp|P1|A_HUMAN\n"
    )
    out = tmp_path / "psm.parquet"
    with FragPipePsmAdapter() as adapter:
        adapter.convert(psm_path=str(psm), output_path=str(out))
    scores = {s["score_name"]: s for s in pq.read_table(out).to_pylist()[0]["additional_scores"]}
    assert scores["parent_ion_fraction"]["score_value"] == pytest.approx(0.87)
    assert scores["parent_ion_fraction"]["higher_better"] is True


@pytest.mark.parametrize(
    ("protein", "expected"),
    [
        ("sp|P1|A_HUMAN", False),
        ("rev_sp|P1|A_HUMAN", True),
        ("REV_sp|P1|A_HUMAN", True),
        ("DECOY_sp|P1|A_HUMAN", True),
        ("sp|P1|A_HUMAN, rev_sp|P2|B_HUMAN", True),
    ],
)
def test_psm_decoy_prefixes_consistent(tmp_path, protein, expected):
    """#247 (H): psm.tsv decoy detection recognises the same prefix set as the
    feature/pg adapters and inspects every comma-separated accession."""
    psm = tmp_path / "psm.tsv"
    psm.write_text(
        "Spectrum\tPeptide\tCharge\tAssigned Modifications\t"
        "Observed M/Z\tCalculated M/Z\tRetention\tProtein\n"
        f"run_01.00500.00500.2\tPEPTIDEK\t2\t\t500.001\t500.000\t123.4\t{protein}\n"
    )
    out = tmp_path / "psm.parquet"
    with FragPipePsmAdapter() as adapter:
        adapter.convert(psm_path=str(psm), output_path=str(out))
    df = pq.read_table(str(out)).to_pandas()
    assert bool(df.iloc[0]["is_decoy"]) is expected


_PEPTIDEK_MASS = compute_precursor_mz("PEPTIDEK", 1) - PROTON_MASS
_OPEN_COLUMNS = "Spectrum\tPeptide\tCharge\tAssigned Modifications\tCalculated Peptide Mass\tCalculated M/Z\tDelta Mass\tCalibrated Observed M/Z"


def _open_search_row(scan: int, delta: float, localization: str, extra: str = "") -> str:
    observed = (_PEPTIDEK_MASS + delta + 2 * PROTON_MASS) / 2
    calculated = (_PEPTIDEK_MASS + 2 * PROTON_MASS) / 2
    return (
        f"run_01.{scan:05d}.{scan:05d}.2\tPEPTIDEK\t2\t\t{_PEPTIDEK_MASS:.4f}\t{calculated:.4f}\t{delta}\t{observed:.4f}"
        f"\t{localization}\tsp|P1|A_HUMAN{extra}\n"
    )


def _convert_psm_rows(tmp_path, header: str, rows: list[str]) -> list[dict]:
    psm = tmp_path / "psm.tsv"
    psm.write_text(header + "\n" + "".join(rows))
    out = tmp_path / "psm.parquet"
    with FragPipePsmAdapter() as adapter:
        adapter.convert(psm_path=str(psm), output_path=str(out))
    return sorted(pq.read_table(out).to_pylist(), key=lambda row: row["scan"])


def test_open_search_mass_shift_is_part_of_the_peptidoform(tmp_path):
    """Delta Mass goes on the one site MSFragger localized, or to an unknown position with
    the tied sites as candidates; calculated_mz includes it, so it is no mass error."""
    rows = _convert_psm_rows(
        tmp_path,
        _OPEN_COLUMNS + "\tMSFragger Localization\tProtein",
        [
            _open_search_row(1, 129.0804, "PEPtIDEK"),
            _open_search_row(2, 129.0804, "PEPtiDEK"),
            _open_search_row(3, 129.0804, "PEPTIDEK"),
        ],
    )
    assert [row["peptidoform"] for row in rows] == [
        "PEPT[+129.0804]IDEK",
        "[+129.0804]?PEPTIDEK",
        "[+129.0804]?PEPTIDEK",
    ]
    assert rows[0]["modification_candidates"] is None
    candidates = rows[1]["modification_candidates"]
    assert [(c["position"], c["amino_acid"], c["assigned"], c["alternative_group"]) for c in candidates] == [
        (4, "T", False, 1),
        (5, "I", False, 1),
    ]
    assert candidates[0]["mass_delta"] == pytest.approx(129.0804, abs=1e-3)
    assert rows[2]["modification_candidates"] is None
    assert rows[1]["modifications"][0]["name"] == "+129.0804"
    assert rows[1]["modifications"][0]["positions"] == []
    for row in rows:
        assert row["calculated_mz"] == pytest.approx(compute_precursor_mz("[+129.0804]?PEPTIDEK", 2), abs=1e-3)
        assert abs(row["mass_error_ppm"]) < 1


@pytest.mark.parametrize("delta", [0.0012, 1.0041, 2.0071, -1.0029, 3.0098])
def test_isotope_error_is_not_a_mass_shift(tmp_path, delta):
    """A precursor picked on a 13C peak is the same molecule, as in a closed search with isotope_error."""
    rows = _convert_psm_rows(
        tmp_path, _OPEN_COLUMNS + "\tMSFragger Localization\tProtein", [_open_search_row(1, delta, "PEPtIDEK")]
    )
    assert rows[0]["peptidoform"] == "PEPTIDEK"
    assert rows[0]["modification_candidates"] is None


_GLYCAN_MASS = 1864.6342


@pytest.mark.parametrize("moved_to_calculated", [True, False])
def test_glycan_composition_and_fdr_are_kept(tmp_path, moved_to_calculated):
    """PTM-Shepherd writes the glycan to Assigned Modifications and, by default, moves it from
    Delta Mass to Calculated M/Z; the glycan counts once either way."""
    glycan_peptide_mass = _PEPTIDEK_MASS + _GLYCAN_MASS
    peptide_mass = glycan_peptide_mass if moved_to_calculated else _PEPTIDEK_MASS
    observed = (glycan_peptide_mass + 3 * PROTON_MASS) / 3
    row = (
        f"run_01.00001.00001.3\tPEPTIDEK\t3\t5I({_GLYCAN_MASS:.4f})\t{peptide_mass:.4f}\t{(peptide_mass + 3 * PROTON_MASS) / 3:.4f}"
        f"\t{glycan_peptide_mass - peptide_mass:.4f}\t{observed:.4f}\tHexNAc(2)Hex(9) % 1864.6341\t5.33\t0.002\tsp|P1|A_HUMAN\n"
    )
    rows = _convert_psm_rows(tmp_path, _OPEN_COLUMNS + "\tTotal Glycan Composition\tGlycan Score\tGlycan q-value\tProtein", [row])
    assert rows[0]["peptidoform"] == "PEPTI[+1864.6342]DEK"
    assert rows[0]["modifications"][0]["name"] == "Glycan:HexNAc2Hex9"
    assert rows[0]["modifications"][0]["positions"][0]["position"] == 5
    assert abs(rows[0]["mass_error_ppm"]) < 1
    scores = {score["score_name"]: score for score in rows[0]["additional_scores"]}
    assert scores["glycan_score"]["score_value"] == pytest.approx(5.33)
    assert scores["glycan_score"]["higher_better"] is True
    assert scores["glycan_qvalue"]["score_value"] == pytest.approx(0.002)
    assert scores["glycan_qvalue"]["higher_better"] is False


def test_feature_is_not_enriched_from_a_mass_shifted_psm(tmp_path):
    """A PSM with an unassigned mass shift is another molecule than the ion of its bare peptidoform."""
    ion = tmp_path / "combined_ion.tsv"
    ion.write_text(
        "Peptide Sequence\tModified Sequence\tCharge\tM/Z\tProtein\tA Intensity\n"
        "PEPTIDEK\tPEPTIDEK\t2\t500.0\tsp|P1|A_HUMAN\t100\n"
    )
    psm = tmp_path / "psm.tsv"
    psm.write_text(_OPEN_COLUMNS + "\tMSFragger Localization\tProtein\n" + _open_search_row(7, 129.0804, "PEPtIDEK"))
    out = tmp_path / "feature.parquet"
    with FragPipeFeatureAdapter() as adapter:
        adapter.convert(feature_path=str(ion), output_path=str(out), psm_path=str(psm), experiment_to_runs={"A": ["run_01"]})
    record = pq.read_table(out).to_pylist()[0]
    assert record["peptidoform"] == "PEPTIDEK"
    assert record["scan"] == []
    assert record["observed_mz"] is None
