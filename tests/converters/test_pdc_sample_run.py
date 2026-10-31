"""Tests for the PDC-metadata sample/run builder (qpx.converters.cdap.pdc_sample_run).

The three PDC GraphQL calls (experimental design, biospecimen, run->plex) are
mocked, so these run fully offline. They lock the channel->sample mapping, the
CDAP label format, the reference-channel and label-free branches, and the
sample.parquet / run.parquet schema compliance.
"""

from __future__ import annotations

from unittest.mock import patch

import pyarrow.parquet as pq
import pytest

from qpx.pdc.metadata import TMT_CHANNEL_FIELDS

_BUILDER = "qpx.converters.cdap.pdc_sample_run"


def _tmt10_plex(plex_id, aliquots):
    """Build a TMT10 plex dict with 10 (aliquot_id, aliquot_submitter_id) channels."""
    plex = {"study_run_metadata_submitter_id": plex_id, "experiment_type": "TMT10", "label_free": None}
    for field, (aid, asid) in zip(TMT_CHANNEL_FIELDS[:10], aliquots):
        plex[field] = [{"aliquot_id": aid, "aliquot_submitter_id": asid}]
    return plex


def _biospecimen(aliquots, **overrides):
    """Build {aliquot_id: record} for the given (aliquot_id, aliquot_submitter_id) list."""
    base = {
        "taxon": "Homo sapiens",
        "primary_site": "Breast",
        "disease_type": "Breast Invasive Carcinoma",
        "sample_type": "Primary Tumor",
        "case_submitter_id": "11BR047",
    }
    base.update(overrides)
    out = {}
    for aid, asid in aliquots:
        out[aid] = {"aliquot_id": aid, "aliquot_submitter_id": asid, **base}
    return out


def _build(tmp_path, design, biospecimen, run_to_plex, study="PDC_TEST"):
    from qpx.converters.cdap.pdc_sample_run import build_sample_run_from_pdc

    with (
        patch(f"{_BUILDER}.fetch_experimental_design", return_value=design),
        patch(f"{_BUILDER}.fetch_biospecimen", return_value=biospecimen),
        patch(f"{_BUILDER}.map_runs_to_plex", return_value=run_to_plex),
    ):
        return build_sample_run_from_pdc(study, tmp_path, prefix=study)


def test_tmt10_builds_sample_run_with_channel_mapping(tmp_path):
    aliquots = [(f"al{i}", f"sub{i}") for i in range(9)] + [("alRef", "Internal Reference - Pooled Sample")]
    plex_id = "01CPTAC_Proteome_20160911"
    design = [_tmt10_plex(plex_id, aliquots)]
    biospecimen = _biospecimen(aliquots[:9])  # reference deliberately absent from biospecimen
    run_to_plex = {
        "run_f01": {"plex": plex_id, "file_name": "run_f01.raw", "fraction": "1", "instrument": "Orbitrap Fusion"},
        "run_f02": {"plex": plex_id, "file_name": "run_f02.raw", "fraction": "2", "instrument": "Orbitrap Fusion"},
    }

    result = _build(tmp_path, design, biospecimen, run_to_plex)
    assert result is not None

    run_tbl = pq.read_table(str(result["run"]))
    assert run_tbl.num_rows == 2  # one row per run file
    rows = {r["run_file_name"]: r for r in run_tbl.to_pylist()}
    samples = rows["run_f01"]["samples"]
    labels = [s["label"] for s in samples]
    # Canonical TMT10 labels, all 10 channels, mass-ascending
    assert labels == [
        "TMT126",
        "TMT127N",
        "TMT127C",
        "TMT128N",
        "TMT128C",
        "TMT129N",
        "TMT129C",
        "TMT130N",
        "TMT130C",
        "TMT131",
    ]
    # 131 channel maps to the pooled reference
    assert samples[-1]["sample_accession"] == "Internal Reference - Pooled Sample"
    assert rows["run_f01"]["fraction"] == "1"
    assert rows["run_f01"]["instrument"] == "Orbitrap Fusion"

    sample_tbl = pq.read_table(str(result["sample"]))
    srows = {r["sample_accession"]: r for r in sample_tbl.to_pylist()}
    # every run.samples.sample_accession resolves to a sample row
    assert set(s["sample_accession"] for s in samples) <= set(srows)
    assert srows["sub0"]["organism"] == "Homo sapiens"
    assert srows["sub0"]["organism_part"] == "Breast"
    assert srows["sub0"]["disease"] == "Breast Invasive Carcinoma"
    assert srows["sub0"]["sample_type"] == "Primary Tumor"
    # reference channel special-cased (absent from biospecimen)
    ref = srows["Internal Reference - Pooled Sample"]
    assert ref["organism"] == "Homo sapiens"
    assert ref["organism_part"] == "Not Reported"
    assert ref["sample_type"] == "Internal Reference"


def test_label_free_single_lfq_channel(tmp_path):
    plex_id = "PLEX_LF"
    design = [
        {
            "study_run_metadata_submitter_id": plex_id,
            "experiment_type": "Label Free",
            "label_free": [{"aliquot_id": "alLF", "aliquot_submitter_id": "subLF"}],
        }
    ]
    biospecimen = _biospecimen([("alLF", "subLF")], primary_site="Lung", disease_type="Lung Adenocarcinoma")
    run_to_plex = {"lf_run": {"plex": plex_id, "file_name": "lf_run.raw", "fraction": None, "instrument": "QE"}}

    result = _build(tmp_path, design, biospecimen, run_to_plex)
    run_tbl = pq.read_table(str(result["run"]))
    samples = run_tbl.to_pylist()[0]["samples"]
    assert len(samples) == 1
    assert samples[0]["label"] == "LFQ"
    assert samples[0]["sample_accession"] == "subLF"

    sample_tbl = pq.read_table(str(result["sample"]))
    srows = {r["sample_accession"]: r for r in sample_tbl.to_pylist()}
    assert srows["subLF"]["organism_part"] == "Lung"


def test_unsupported_experiment_type_skips_channels(tmp_path):
    """TMT16 is not modelled by CDAP; its runs get no channel mapping (empty samples)."""
    plex_id = "PLEX_TMT16"
    design = [{"study_run_metadata_submitter_id": plex_id, "experiment_type": "TMT16", "label_free": None}]
    run_to_plex = {"t16_run": {"plex": plex_id, "file_name": "t16_run.raw", "fraction": "1", "instrument": "Lumos"}}

    result = _build(tmp_path, design, {}, run_to_plex)
    run_tbl = pq.read_table(str(result["run"]))
    assert run_tbl.to_pylist()[0]["samples"] == []


def test_channel_count_mismatch_skips_plex(tmp_path):
    """A plex declaring TMT10 but missing a channel must not emit a wrong mapping."""
    plex_id = "PLEX_BAD"
    aliquots = [(f"al{i}", f"sub{i}") for i in range(10)]
    plex = _tmt10_plex(plex_id, aliquots)
    del plex["tmt_131"]  # only 9 channels now
    run_to_plex = {"bad_run": {"plex": plex_id, "file_name": "bad_run.raw", "fraction": "1", "instrument": "X"}}

    result = _build(tmp_path, [plex], _biospecimen(aliquots), run_to_plex)
    run_tbl = pq.read_table(str(result["run"]))
    assert run_tbl.to_pylist()[0]["samples"] == []


def test_empty_run_to_plex_returns_none(tmp_path):
    result = _build(tmp_path, [], {}, {})
    assert result is None
    assert not (tmp_path / "PDC_TEST.run.parquet").exists()
    assert not (tmp_path / "PDC_TEST.sample.parquet").exists()


def _tmt11_plex(plex_id, aliquots, *, with_tmtxx=False):
    """Build a TMT11 plex; ``with_tmtxx`` also populates the legacy 132N/C channels."""
    n_channels = 13 if with_tmtxx else 11
    plex = {"study_run_metadata_submitter_id": plex_id, "experiment_type": "TMT11", "label_free": None}
    for field, (aid, asid) in zip(TMT_CHANNEL_FIELDS[:n_channels], aliquots):
        plex[field] = [{"aliquot_id": aid, "aliquot_submitter_id": asid}]
    return plex


def test_tmt11_with_tmtxx_extra_channels_maps_all(tmp_path):
    """A TMT11 plex that also carries TMTXX 132N/C maps all 13 channels, not skips.

    CDAP appends the legacy TMTXX channels for TMT11 studies, so the metadata
    builder must too -- otherwise these 15 known pilot studies get empty samples.
    """
    aliquots = [(f"al{i}", f"sub{i}") for i in range(13)]
    plex_id = "01CPTAC_TMT11_TMTXX"
    design = [_tmt11_plex(plex_id, aliquots, with_tmtxx=True)]
    run_to_plex = {"t11_run": {"plex": plex_id, "file_name": "t11_run.raw", "fraction": "1", "instrument": "Lumos"}}

    result = _build(tmp_path, design, _biospecimen(aliquots), run_to_plex)
    labels = [s["label"] for s in pq.read_table(str(result["run"])).to_pylist()[0]["samples"]]
    assert len(labels) == 13
    assert labels[0] == "TMT126"
    assert labels[10] == "TMT131C"
    assert labels[-2:] == ["TMT132N", "TMT132C"]


def test_build_failure_preserves_existing_outputs(tmp_path):
    """A build failure must not clobber a previously-good sample/run parquet."""
    from qpx.converters.cdap.pdc_sample_run import build_sample_run_from_pdc

    good_sample = tmp_path / "PDC_TEST.sample.parquet"
    good_run = tmp_path / "PDC_TEST.run.parquet"
    good_sample.write_bytes(b"PRIOR-GOOD-SAMPLE")
    good_run.write_bytes(b"PRIOR-GOOD-RUN")

    with (
        patch(f"{_BUILDER}.fetch_experimental_design", return_value=[]),
        patch(f"{_BUILDER}.fetch_biospecimen", return_value={}),
        patch(f"{_BUILDER}.map_runs_to_plex", side_effect=RuntimeError("network down")),
    ):
        with pytest.raises(RuntimeError):
            build_sample_run_from_pdc("PDC_TEST", tmp_path, prefix="PDC_TEST")

    assert good_sample.read_bytes() == b"PRIOR-GOOD-SAMPLE"
    assert good_run.read_bytes() == b"PRIOR-GOOD-RUN"
    assert not (tmp_path / "PDC_TEST.sample.parquet.tmp").exists()
    assert not (tmp_path / "PDC_TEST.run.parquet.tmp").exists()


def test_strip_raw_extension_matches_cdap():
    """A PDC-derived run_file_name matches the CDAP stem (else run<->feature joins break)."""
    from qpx.converters.cdap.base_adapter import CdapBaseAdapter
    from qpx.pdc.metadata import _strip_raw_extension

    for name in ("run_f01.raw", "run_f01.mzML", "run_f01.mgf", "run_f01.mzML.gz", "run_f01.d", "run_f01"):
        assert _strip_raw_extension(name) == CdapBaseAdapter.strip_run_extension(name)
    assert _strip_raw_extension("run_f01.mzML.gz") == "run_f01"
    assert _strip_raw_extension("run_f01.d") == "run_f01"
