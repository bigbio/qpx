"""
QPX CLI entry point — Click group that registers all subcommand groups.

Usage:
    qpxc convert diann --report-path ...
    qpxc transform gene-map --parquet-path ... --fasta ...
    qpxc query sql --dataset-path ... --sql "SELECT ..."
    qpxc info --dataset-path ...
"""

import logging

import click

from qpx._version import __version__
from qpx.cli.convert import convert
from qpx.cli.info import info
from qpx.cli.migrate import migrate_cmd
from qpx.cli.ontology import ontology
from qpx.cli.pdc2qpx import pdc2qpx_cmd
from qpx.cli.query import query_cmd
from qpx.cli.transform import transform
from qpx.cli.validate import validate_cmd

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.version_option(version=__version__, package_name="qpx", message="%(package)s %(version)s")
@click.group(context_settings=CONTEXT_SETTINGS)
def qpx_main():
    """qpxc -- quantitative proteomics data format tools.

    Convert, transform, query, and inspect QPX datasets from the command line.
    """
    logging.basicConfig(
        level=logging.INFO,
        datefmt="%H:%M:%S",
        format="[%(asctime)s] %(levelname).1s | %(name)s | %(message)s",
    )


# Register subcommand groups
qpx_main.add_command(convert)
qpx_main.add_command(transform)
qpx_main.add_command(query_cmd, name="query")
qpx_main.add_command(info)
qpx_main.add_command(validate_cmd, name="validate")
qpx_main.add_command(ontology)
qpx_main.add_command(pdc2qpx_cmd, name="pdc2qpx")
qpx_main.add_command(migrate_cmd, name="migrate")


def main():
    """Entry point for the ``qpxc`` console script."""
    qpx_main()


if __name__ == "__main__":
    main()
