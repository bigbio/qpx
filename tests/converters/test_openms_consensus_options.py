"""Option and metadata behaviour of the consensusXML converter.

Split out of ``test_openms_consensus.py`` to keep that module under the
1000-line limit CodeFactor enforces. These classes cover the converter's
*options* — which PSMs are written, and the metadata that drives derived fields
— rather than the core consensusXML -> QPX extraction the sibling module tests.
"""

import pytest

pytest.importorskip("pyopenms")

from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter  # noqa: E402
from tests.converters.test_openms_consensus import _TMT_CONSENSUSXML  # noqa: E402


class TestSkipUnassignedPsms:
    """`include_unassigned_psms=False` drops PSMs that link to no feature.

    Unassigned identifications carry no quantification and their feature_id is
    null, so a psm->feature join silently drops them anyway — 41% of rows on a
    real label-free dataset (bigbio/qpx#299). The option makes that explicit.
    """

    _XML_WITH_UNASSIGNED = _TMT_CONSENSUSXML.replace(
        "</consensusElementList>",
        """</consensusElementList>
  <UnassignedPeptideIdentification identification_run_ref="PI_0" score_type=""
    higher_score_better="true" significance_threshold="0" MZ="500.25" RT="200"
    spectrum_reference="controllerType=0 controllerNumber=1 scan=99">
    <PeptideHit score="0" sequence="UNASSIGNEDK" charge="2" protein_refs="PH_0">
      <UserParam type="string" name="target_decoy" value="target"/>
      <UserParam type="float" name="Posterior Error Probability_score" value="1.0e-03"/>
    </PeptideHit>
  </UnassignedPeptideIdentification>""",
    )

    @staticmethod
    def _convert(tmp_path, name, *, include, stream):
        import pyarrow.parquet as pq

        from qpx.converters.openms_consensus import converter as conv

        path = tmp_path / f"{name}.consensusXML"
        path.write_text(TestSkipUnassignedPsms._XML_WITH_UNASSIGNED)
        out = tmp_path / f"out_{name}"
        conv._should_stream = lambda _p, v=stream: v
        OpenMSConsensusConverter().convert(
            str(path),
            str(out),
            output_prefix="X",
            structures=("feature", "psm", "pg"),
            project_accession="X",
            include_unassigned_psms=include,
        )
        return pq.read_table(str(out / "X.psm.parquet")).to_pylist()

    @pytest.mark.parametrize("stream", [False, True], ids=["memory", "streaming"])
    def test_default_drops_unassigned(self, tmp_path, stream):
        """The default PSM view is quantified-only."""
        import pyarrow.parquet as pq

        from qpx.converters.openms_consensus import converter as conv

        path = tmp_path / f"default_{stream}.consensusXML"
        path.write_text(self._XML_WITH_UNASSIGNED)
        out = tmp_path / f"out_default_{stream}"
        conv._should_stream = lambda _p, v=stream: v
        OpenMSConsensusConverter().convert(
            str(path),
            str(out),
            output_prefix="X",
            structures=("feature", "psm", "pg"),
            project_accession="X",
        )
        rows = pq.read_table(str(out / "X.psm.parquet")).to_pylist()
        assert "UNASSIGNEDK" not in {r["sequence"] for r in rows}
        assert all(r["feature_id"] is not None for r in rows)

    @pytest.mark.parametrize("stream", [False, True], ids=["memory", "streaming"])
    def test_opt_in_keeps_unassigned(self, tmp_path, stream):
        rows = self._convert(tmp_path, f"keep_{stream}", include=True, stream=stream)
        seqs = {r["sequence"] for r in rows}
        assert "UNASSIGNEDK" in seqs, "opt-in must restore them"
        assert any(r["feature_id"] is None for r in rows)

    @pytest.mark.parametrize("stream", [False, True], ids=["memory", "streaming"])
    def test_option_drops_unassigned_on_both_paths(self, tmp_path, stream):
        rows = self._convert(tmp_path, f"drop_{stream}", include=False, stream=stream)
        seqs = {r["sequence"] for r in rows}
        assert "UNASSIGNEDK" not in seqs, "unassigned PSM was still written"
        assert rows, "assigned PSMs must survive"
        assert all(r["feature_id"] is not None for r in rows), "every remaining PSM links to a feature"

    def test_protein_inference_is_unaffected(self, tmp_path):
        """Dropping PSM rows must not change the protein groups."""
        import pyarrow.parquet as pq

        from qpx.converters.openms_consensus import converter as conv

        counts = []
        for include in (True, False):
            path = tmp_path / f"pg_{include}.consensusXML"
            path.write_text(self._XML_WITH_UNASSIGNED)
            out = tmp_path / f"pgout_{include}"
            conv._should_stream = lambda _p: False
            OpenMSConsensusConverter().convert(
                str(path),
                str(out),
                output_prefix="X",
                structures=("feature", "psm", "pg"),
                project_accession="X",
                include_unassigned_psms=include,
            )
            counts.append(pq.read_table(str(out / "X.pg.parquet")).num_rows)
        assert counts[0] == counts[1], "pg output changed when only PSM rows were filtered"


class TestMissedCleavages:
    """OpenMS output must report missed cleavages, as the DIA-NN path does.

    The value is a property of the peptide and the search enzyme, both of which
    the converter has: the enzyme is in the consensusXML SearchParameters
    (bigbio/qpx#300).
    """

    # PEPTIDEK has no internal K/R; PEPKTIDEK has one internal K -> 1 missed cleavage.
    _TRYPSIN = _TMT_CONSENSUSXML.replace('enzyme="unknown_enzyme"', 'enzyme="Trypsin"')
    _TRYPSIN_MISSED = _TRYPSIN.replace('sequence="PEPTIDEK"', 'sequence="PEPKTIDEK"')

    @staticmethod
    def _features(tmp_path, name, xml, stream):
        import pyarrow.parquet as pq

        from qpx.converters.openms_consensus import converter as conv

        path = tmp_path / f"{name}.consensusXML"
        path.write_text(xml)
        out = tmp_path / f"out_{name}"
        conv._should_stream = lambda _p, v=stream: v
        OpenMSConsensusConverter().convert(
            str(path),
            str(out),
            output_prefix="X",
            structures=("feature", "psm", "pg"),
            project_accession="X",
        )
        return (
            pq.read_table(str(out / "X.feature.parquet")).to_pylist(),
            pq.read_table(str(out / "X.psm.parquet")).to_pylist(),
        )

    @pytest.mark.parametrize("stream", [False, True], ids=["memory", "streaming"])
    def test_zero_missed_cleavages(self, tmp_path, stream):
        feats, psms = self._features(tmp_path, f"mc0_{stream}", self._TRYPSIN, stream)
        assert feats and all(f["missed_cleavages"] == 0 for f in feats)
        assert psms and all(p["missed_cleavages"] == 0 for p in psms)

    @pytest.mark.parametrize("stream", [False, True], ids=["memory", "streaming"])
    def test_one_missed_cleavage(self, tmp_path, stream):
        feats, psms = self._features(tmp_path, f"mc1_{stream}", self._TRYPSIN_MISSED, stream)
        assert feats and all(f["missed_cleavages"] == 1 for f in feats)
        assert psms and all(p["missed_cleavages"] == 1 for p in psms)

    @pytest.mark.parametrize("stream", [False, True], ids=["memory", "streaming"])
    def test_unknown_enzyme_leaves_it_null(self, tmp_path, stream):
        """No usable enzyme means unknown, not a guessed zero."""
        feats, psms = self._features(tmp_path, f"mcnull_{stream}", _TMT_CONSENSUSXML, stream)
        assert feats and all(f["missed_cleavages"] is None for f in feats)
        assert psms and all(p["missed_cleavages"] is None for p in psms)


class TestEnzymeResolution:
    """SDRF first, consensusXML second — and say so when they disagree."""

    class _Map:
        def __init__(self, enzyme):
            self._enzyme = enzyme

        def getProteinIdentifications(self):
            return []

    def _sdrf(self, tmp_path, enzyme):
        path = tmp_path / "x.sdrf.tsv"
        path.write_text(f"source name\tcomment[data file]\tcomment[cleavage agent details]\nS1\trun_01.mzML\tNT={enzyme}\n")
        return path

    def test_sdrf_wins_over_the_consensusxml(self, tmp_path):
        from qpx.converters.openms_consensus.feature_adapter import resolve_enzyme

        assert resolve_enzyme(self._Map("Trypsin"), self._sdrf(tmp_path, "Lys-C")) == "Lys-C"

    def test_falls_back_to_the_consensusxml_without_an_sdrf(self, tmp_path):
        from qpx.converters.openms_consensus.feature_adapter import resolve_enzyme

        assert resolve_enzyme(self._Map("Trypsin"), None) == "Trypsin"

    def test_disagreement_is_reported(self, tmp_path, caplog):
        import logging

        from qpx.converters.openms_consensus.feature_adapter import resolve_enzyme

        with caplog.at_level(logging.WARNING):
            resolve_enzyme(self._Map("Trypsin"), self._sdrf(tmp_path, "Lys-C"))
        assert any("enzyme mismatch" in r.message for r in caplog.records)

    def test_agreement_is_quiet(self, tmp_path, caplog):
        import logging

        from qpx.converters.openms_consensus.feature_adapter import resolve_enzyme

        with caplog.at_level(logging.WARNING):
            assert resolve_enzyme(self._Map("Trypsin"), self._sdrf(tmp_path, "Trypsin")) == "Trypsin"
        assert not any("enzyme mismatch" in r.message for r in caplog.records)

    def test_no_source_at_all_is_none(self):
        from qpx.converters.openms_consensus.feature_adapter import resolve_enzyme

        assert resolve_enzyme(self._Map("unknown_enzyme"), None) is None


def test_public_records_api_reports_missed_cleavages(tmp_path):
    """The library entry point must not lag the CLI (CodeRabbit, PR #296).

    consensus_features_to_records is exported in __all__, so it is a supported
    way to get feature records; it resolved no enzyme, so callers got a null
    missed_cleavages while the same data through the converter got a value.
    """
    from qpx.converters.openms_consensus import consensus_features_to_records

    path = tmp_path / "in.consensusXML"
    path.write_text(_TMT_CONSENSUSXML.replace('enzyme="unknown_enzyme"', 'enzyme="Trypsin"'))
    records = consensus_features_to_records(str(path))
    assert records
    assert all(r["missed_cleavages"] == 0 for r in records), "PEPTIDEK has no internal K/R"

    unknown = tmp_path / "unknown.consensusXML"
    unknown.write_text(_TMT_CONSENSUSXML)
    assert all(r["missed_cleavages"] is None for r in consensus_features_to_records(str(unknown)))
