"""Spectronaut PG conversion must distinguish missing quantities from zero."""

import pyarrow as pa
import pyarrow.parquet as pq

from qpx.converters.spectronaut.pg_adapter import SpectronautPgAdapter


def test_pg_conversion_preserves_missing_and_zero_quantities(tmp_path):
    """Retain identified groups and their confidence when PG.Quantity is absent."""
    source = tmp_path / "report.parquet"
    output = tmp_path / "result.pg.parquet"
    accessions = ["P_NULL", "P_NAN", "P_ZERO", "P_QUANTIFIED"]
    pq.write_table(
        pa.table(
            {
                "PG.ProteinGroups": accessions,
                "R.FileName": ["run_01.raw"] * 4,
                "PG.Quantity": pa.array([None, float("nan"), 0.0, 123.5], type=pa.float64()),
                "PG.Qvalue": [0.001] * 4,
                "PEP.StrippedSequence": ["PEPTIDEK"] * 4,
            }
        ),
        source,
    )

    with SpectronautPgAdapter(duckdb_threads=1) as adapter:
        adapter.convert(str(source), str(output))
    rows = pq.read_table(output).to_pylist()

    assert len(rows) == 4
    assert {row["anchor_protein"]: row["intensity"] for row in rows} == {
        "P_NULL": None,
        "P_NAN": None,
        "P_ZERO": 0.0,
        "P_QUANTIFIED": 123.5,
    }
    assert all(row["grouped_runs"] == ["run_01"] and row["label"] == "LFQ" for row in rows)
    assert all(row["global_qvalue"] == 0.001 for row in rows)
