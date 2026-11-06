"""Convert commands and their public import paths."""

import click

from .cdap import convert_cdap_cmd
from .diann import convert_diann_cmd, convert_diann_xic_cmd
from .fragpipe import convert_fragpipe_cmd
from .maxquant import convert_maxquant_cmd
from .mz import convert_mz_cmd
from .mzidentml import convert_mzidentml_cmd
from .openms import convert_openms_cmd
from .openms_consensus import convert_openms_consensus_cmd
from .quantms_msstats import convert_quantms_msstats_cmd
from .sdrf import convert_sdrf_cmd
from .spectronaut import convert_spectronaut_cmd

__all__ = [
    "convert",
    "convert_cdap_cmd",
    "convert_diann_cmd",
    "convert_diann_xic_cmd",
    "convert_fragpipe_cmd",
    "convert_maxquant_cmd",
    "convert_mz_cmd",
    "convert_mzidentml_cmd",
    "convert_openms_cmd",
    "convert_openms_consensus_cmd",
    "convert_quantms_msstats_cmd",
    "convert_sdrf_cmd",
    "convert_spectronaut_cmd",
]


@click.group()
def convert():
    """Convert external tool outputs to QPX format."""


for _command in (
    convert_diann_cmd,
    convert_diann_xic_cmd,
    convert_maxquant_cmd,
    convert_fragpipe_cmd,
    convert_mzidentml_cmd,
    convert_openms_cmd,
    convert_openms_consensus_cmd,
    convert_sdrf_cmd,
    convert_quantms_msstats_cmd,
    convert_spectronaut_cmd,
    convert_cdap_cmd,
    convert_mz_cmd,
):
    convert.add_command(_command)
