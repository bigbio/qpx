"""DuckDB resource options of the CLI commands."""

import click
import pytest
from click.testing import CliRunner

from qpx.cli.main import qpx_main
from qpx.cli.resources import convert_resource_options
from qpx.core.engine import create_connection

_SETTINGS = "SELECT current_setting('threads') AS threads, current_setting('temp_directory') AS temp"


@click.command()
@click.option("--output-folder", type=click.Path(file_okay=False))
@click.option("--max-cpus", type=int)
@convert_resource_options
def _probe(**_options):
    """Print the settings of a connection opened inside the command."""
    con = create_connection()
    click.echo(con.execute(_SETTINGS).fetchall()[0])
    con.close()


@pytest.mark.parametrize(
    "command",
    [
        ["convert", "cdap"],
        ["convert", "diann"],
        ["convert", "fragpipe"],
        ["convert", "maxquant"],
        ["convert", "mzidentml"],
        ["convert", "openms-consensus"],
        ["convert", "quantms-msstats"],
        ["convert", "spectronaut"],
        ["pdc2qpx"],
        ["validate"],
        ["query", "sql"],
        ["transform", "gene-map"],
        ["transform", "protein-properties"],
    ],
)
def test_commands_offer_a_temp_dir(command):
    """Every command that runs DuckDB over large views can move its spill files."""
    result = CliRunner().invoke(qpx_main, [*command, "--help"])
    assert result.exit_code == 0, result.output
    assert "--temp-dir" in result.output


def test_convert_spills_into_the_output_folder(tmp_path):
    """A convert command's connections spill under its output folder and use --max-cpus, unless --temp-dir is given."""
    output = tmp_path / "out"
    result = CliRunner().invoke(_probe, ["--output-folder", str(output), "--max-cpus", "1"])
    assert result.output.strip() == str((1, str(output / ".tmp")))
    elsewhere = tmp_path / "spill"
    result = CliRunner().invoke(_probe, ["--output-folder", str(output), "--temp-dir", str(elsewhere)])
    assert str(elsewhere) in result.output


def test_query_sql_applies_the_resource_options(dataset_dir, tmp_path):
    """--duckdb-threads and --temp-dir reach the dataset's connection."""
    spill = tmp_path / "spill"
    result = CliRunner().invoke(
        qpx_main,
        [
            "query",
            "sql",
            "--dataset-path",
            str(dataset_dir),
            "--sql",
            _SETTINGS,
            "--duckdb-threads",
            "1",
            "--temp-dir",
            str(spill),
        ],
    )
    assert result.exit_code == 0, result.output
    assert f"1,{spill}" in result.output
