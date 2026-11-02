"""The Python blocks of the user guides run on a converted dataset."""

import contextlib
import io
import os
import re
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main

pytest.importorskip("plotly")
pytest.importorskip("matplotlib")

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "tests" / "examples"
GUIDES = ["stats.md", "visualize.md", "project.md"]
_PYTHON_BLOCK = re.compile(r"^```python\n(.*?)^```", re.S | re.M)


@pytest.fixture(scope="module")
def guide_cwd(tmp_path_factory):
    """A folder with ./output/ holding a small DIA-NN dataset and a PSM view."""
    cwd = tmp_path_factory.mktemp("guide")
    small = EXAMPLES / "diann" / "small"
    result = CliRunner().invoke(
        qpx_main,
        [
            "convert",
            "diann",
            "--report-path",
            str(small / "diann_report.tsv"),
            "--pg-matrix-path",
            str(small / "diann_report.pg_matrix.tsv"),
            "--sdrf-file",
            str(small / "PXD019909-DIA.sdrf.tsv"),
            "--output-folder",
            str(cwd / "output"),
            "--output-prefix",
            "guide",
            "--no-mudata",
        ],
    )
    assert result.exit_code == 0, result.output
    result = CliRunner().invoke(
        qpx_main,
        [
            "convert",
            "maxquant",
            "--msms-file",
            str(EXAMPLES / "maxquant" / "maxquant_simple" / "msms.txt"),
            "--structures",
            "psm",
            "--output-folder",
            str(cwd / "maxquant"),
        ],
    )
    assert result.exit_code == 0, result.output
    psm = next((cwd / "maxquant").glob("*.psm.parquet"))
    shutil.copy(psm, cwd / "output" / "guide.psm.parquet")
    return cwd


def _blocks():
    for guide in GUIDES:
        text = (ROOT / "docs" / "guide" / guide).read_text(encoding="utf-8")
        for number, block in enumerate(_PYTHON_BLOCK.findall(text), start=1):
            yield pytest.param(block, id=f"{guide}-{number}")


@pytest.mark.parametrize("block", _blocks())
def test_guide_block_runs(block, guide_cwd):
    previous = Path.cwd()
    os.chdir(guide_cwd)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(block, "guide", "exec"), {"__name__": "__guide__"})  # noqa: S102
    finally:
        os.chdir(previous)
