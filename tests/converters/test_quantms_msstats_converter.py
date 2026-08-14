"""Regression tests for QuantMS MSstats plus SDRF conversion."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.converters.quantms_msstats import QuantmsMsstatsConverter
from qpx.core.data import FeatureSchema

_TMT6_LABELS = ["TMT126", "TMT127", "TMT128", "TMT129", "TMT130", "TMT131"]


def _write_sdrf(path: Path, runs: dict[str, list[str]]) -> None:
    records = []
    sample_number = 1
    for run_file, labels in runs.items():
        for label in labels:
            records.append(
                {
                    "source name": f"S{sample_number}",
                    "comment[data file]": run_file,
                    "comment[label]": label,
                    "characteristics[organism]": "Homo sapiens",
                    "characteristics[organism part]": "cell",
                    "characteristics[biological replicate]": str(sample_number),
                }
            )
            sample_number += 1
    pd.DataFrame(records).to_csv(path, sep="\t", index=False)


def _convert(msstats: Path, sdrf: Path, output: Path, prefix: str = "test") -> Path:
    QuantmsMsstatsConverter(max_cpus=24).convert(
        msstats_file=msstats,
        sdrf_file=sdrf,
        output_folder=output,
        output_prefix=prefix,
        project_accession="PXDTEST",
        batch_size=2,
    )
    return output / f"{prefix}.feature.parquet"


def test_lfq_conversion_writes_only_supported_views(tmp_path):
    """LFQ rows map through SDRF without fabricating PSM or PG output."""
    sdrf = tmp_path / "lfq.sdrf.tsv"
    msstats = tmp_path / "lfq_msstats_in.csv"
    output = tmp_path / "output"
    _write_sdrf(
        sdrf,
        {
            "run1.mzML": ["label free sample"],
            "run2.mzML": ["label free sample"],
        },
    )
    pd.DataFrame(
        [
            {
                "ProteinName": "P1",
                "PeptideSequence": "PEPTIDEK",
                "PrecursorCharge": 2,
                "Run": 1,
                "Intensity": 100.0,
                "Reference": "run1.mzML",
            },
            {
                "ProteinName": "P1;P2",
                "PeptideSequence": "PEPTM(Oxidation)IDEK",
                "PrecursorCharge": 3,
                "Run": 2,
                "Intensity": 200.0,
                "Reference": "run2.mzML",
            },
        ]
    ).to_csv(msstats, index=False)

    feature_path = _convert(msstats, sdrf, output)
    table = pq.read_table(feature_path)

    assert table.num_rows == 2
    assert not FeatureSchema.validate(table)
    assert table.column("run_file_name").to_pylist() == ["run1", "run2"]
    assert table.column("peptidoform").to_pylist() == [
        "PEPTIDEK",
        "PEPTM[UNIMOD:35]IDEK",
    ]
    assert table.column("unique").to_pylist() == [True, False]
    assert [values[0]["label"] for values in table.column("intensities").to_pylist()] == ["LFQ", "LFQ"]
    assert sum(values[0]["intensity"] for values in table.column("intensities").to_pylist()) == 300.0
    assert len(set(table.column("feature_id").to_pylist())) == 2

    expected = {
        "test.feature.parquet",
        "test.sample.parquet",
        "test.run.parquet",
        "test.dataset.parquet",
        "test.provenance.parquet",
        "test.ontology.parquet",
    }
    assert {path.name for path in output.iterdir()} == expected
    assert not (output / "test.psm.parquet").exists()
    assert not (output / "test.pg.parquet").exists()

    metadata = pq.ParquetFile(feature_path).metadata.metadata
    assert metadata[b"identity_composite"] == b"run_file_name,peptidoform,charge,rt,scan"
    provenance = pq.read_table(output / "test.provenance.parquet").to_pylist()
    config = json.loads(provenance[1]["config"])
    assert config["protein_aggregation"] == "not performed"


def test_quoted_comma_after_csv_sniff_sample_is_parsed(tmp_path):
    """A late quoted comma must not change the inferred CSV dialect."""
    sdrf = tmp_path / "lfq.sdrf.tsv"
    msstats = tmp_path / "lfq_msstats_in.csv"
    output = tmp_path / "output"
    _write_sdrf(sdrf, {"run1.mzML": ["label free sample"]})

    common = {
        "ProteinName": "P1",
        "PeptideSequence": "PEPTIDEK",
        "Charge": 2,
        "Run": "run1",
        "Intensity": 100.0,
    }
    rows = [common] * 20_481
    rows.append(
        {
            **common,
            "ProteinName": "Interleukin-8,",
            "PeptideSequence": "QUOTEDK",
            "Intensity": 200.0,
        }
    )
    pd.DataFrame(rows).to_csv(msstats, index=False)

    table = pq.read_table(_convert(msstats, sdrf, output))

    assert table.num_rows == 2
    assert "Interleukin-8," in table.column("anchor_protein").to_pylist()


def test_tmt_channels_collapse_without_merging_locations(tmp_path):
    """TMT channel rows collapse per scan/RT while distinct Features remain separate."""
    sdrf = tmp_path / "tmt.sdrf.tsv"
    msstats = tmp_path / "tmt_msstats_in.csv"
    output = tmp_path / "output"
    _write_sdrf(sdrf, {"plex.mzML": _TMT6_LABELS})

    rows = []
    for scan, rt in ((1, 100.0), (2, 101.0)):
        for channel in range(1, 7):
            rows.append(
                {
                    "ProteinName": "P1;P2",
                    "PeptideSequence": ".(TMT6plex)PEPTIDEK",
                    "Charge": 2,
                    "Channel": channel,
                    "RetentionTime": rt,
                    "Run": 1,
                    "Intensity": scan * 100 + channel,
                    "Reference": f"plex.mzML_controllerType=0 controllerNumber=1 scan={scan}",
                }
            )
    pd.DataFrame(rows).to_csv(msstats, index=False)

    table = pq.read_table(_convert(msstats, sdrf, output))

    assert table.num_rows == 2
    assert table.column("peptidoform").to_pylist() == [
        "[UNIMOD:737]-PEPTIDEK",
        "[UNIMOD:737]-PEPTIDEK",
    ]
    assert table.column("scan").to_pylist() == [[1], [2]]
    assert table.column("rt").to_pylist() == [100.0, 101.0]
    intensities = table.column("intensities").to_pylist()
    assert [[entry["label"] for entry in values] for values in intensities] == [
        _TMT6_LABELS,
        _TMT6_LABELS,
    ]
    assert sum(entry["intensity"] for values in intensities for entry in values) == sum(row["Intensity"] for row in rows)


def test_run_reference_conflict_is_rejected(tmp_path):
    """Run and Reference may not resolve to different SDRF data files."""
    sdrf = tmp_path / "test.sdrf.tsv"
    msstats = tmp_path / "test.csv"
    _write_sdrf(
        sdrf,
        {
            "run1.mzML": ["label free sample"],
            "run2.mzML": ["label free sample"],
        },
    )
    pd.DataFrame(
        [
            {
                "ProteinName": "P1",
                "PeptideSequence": "PEPTIDEK",
                "Charge": 2,
                "Run": "run1",
                "Reference": "run2.mzML",
                "Intensity": 100.0,
            }
        ]
    ).to_csv(msstats, index=False)

    with pytest.raises(ValueError, match="Run and Reference to different SDRF files"):
        _convert(msstats, sdrf, tmp_path / "output")


def test_unmapped_run_is_rejected(tmp_path):
    """Every quantitative row must resolve to a run declared by the SDRF."""
    sdrf = tmp_path / "test.sdrf.tsv"
    msstats = tmp_path / "test.csv"
    _write_sdrf(sdrf, {"run1.mzML": ["label free sample"]})
    pd.DataFrame(
        [
            {
                "ProteinName": "P1",
                "PeptideSequence": "PEPTIDEK",
                "Charge": 2,
                "Run": "missing-run",
                "Intensity": 100.0,
            }
        ]
    ).to_csv(msstats, index=False)

    with pytest.raises(ValueError, match="cannot be mapped to an SDRF data file"):
        _convert(msstats, sdrf, tmp_path / "output")


@pytest.mark.parametrize("intensity", ["NaN", "inf", "-inf", "not-a-number"])
def test_nonfinite_intensity_is_rejected(tmp_path, intensity):
    """Invalid and non-finite quantities must not enter QPX output."""
    sdrf = tmp_path / "test.sdrf.tsv"
    msstats = tmp_path / "test.csv"
    _write_sdrf(sdrf, {"run1.mzML": ["label free sample"]})
    pd.DataFrame(
        [
            {
                "ProteinName": "P1",
                "PeptideSequence": "PEPTIDEK",
                "Charge": 2,
                "Run": "run1",
                "Intensity": intensity,
            }
        ]
    ).to_csv(msstats, index=False)

    with pytest.raises(ValueError, match="contain invalid required values"):
        _convert(msstats, sdrf, tmp_path / "output")


def test_conflicting_feature_channel_intensity_is_rejected(tmp_path):
    """A Feature/channel key cannot silently select one of two intensities."""
    sdrf = tmp_path / "test.sdrf.tsv"
    msstats = tmp_path / "test.csv"
    _write_sdrf(sdrf, {"run1.mzML": ["label free sample"]})
    base = {
        "ProteinName": "P1",
        "PeptideSequence": "PEPTIDEK",
        "Charge": 2,
        "Run": "run1",
        "Reference": "run1.mzML",
    }
    pd.DataFrame([{**base, "Intensity": 100.0}, {**base, "Intensity": 101.0}]).to_csv(msstats, index=False)

    with pytest.raises(ValueError, match="conflicting intensity values"):
        _convert(msstats, sdrf, tmp_path / "output")


def test_cli_quantms_msstats(tmp_path):
    """The public CLI exposes and runs the QuantMS MSstats converter."""
    sdrf = tmp_path / "test.sdrf.tsv"
    msstats = tmp_path / "test.csv"
    output = tmp_path / "output"
    _write_sdrf(sdrf, {"run1.mzML": ["label free sample"]})
    pd.DataFrame(
        [
            {
                "ProteinName": "P1",
                "PeptideSequence": "PEPTIDEK",
                "Charge": 2,
                "Run": "run1",
                "Intensity": 100.0,
            }
        ]
    ).to_csv(msstats, index=False)

    result = CliRunner().invoke(
        qpx_main,
        [
            "convert",
            "quantms-msstats",
            "--msstats-file",
            str(msstats),
            "--sdrf-file",
            str(sdrf),
            "--output-folder",
            str(output),
            "--output-prefix",
            "cli",
            "--max-cpus",
            "24",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (output / "cli.feature.parquet").exists()
