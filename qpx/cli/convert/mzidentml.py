"""Convert CLI mzidentml implementation."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from qpx.cli.resources import convert_resource_options

from .options import _compression_option, _dataset_metadata_options, _output_folder_option, _verbose_option
from .postprocess import _finish_conversion


@click.command("mzidentml")
@click.option(
    "--mzid-path",
    help="Input mzIdentML (.mzid) file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@_output_folder_option
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
    default="mzidentml",
)
@click.option(
    "--mgf-path",
    help="Optional MGF file: index= spectrum references become the scan numbers it records",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--include-spectra",
    help="Attach mz_array and intensity_array from MGF to PSM records",
    is_flag=True,
    default=False,
)
@click.option(
    "--rank1-only",
    help="Keep only the first-ranked hit of each spectrum",
    is_flag=True,
    default=False,
)
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD054720)",
)
@_compression_option
@convert_resource_options
@_verbose_option
@_dataset_metadata_options
def convert_mzidentml_cmd(**kwargs):
    """Convert mzIdentML file to QPX dataset.

    Supports both standard mzIdentML (1.1/1.2) and mzIdentML 1.3 with
    cross-linking extensions (inter-peptide, looplinks, noncovalent).
    Produces a full QPX dataset including PSM, pepmap, provenance,
    ontology, and dataset metadata.

    \b
    Examples:
        # Convert a standard mzIdentML file
        qpxc convert mzidentml \\
            --mzid-path results.mzid \\
            --output-folder ./qpx_output

        # Convert an XL-MS mzIdentML 1.3 file
        qpxc convert mzidentml \\
            --mzid-path crosslinks.mzid \\
            --output-folder ./qpx_output \\
            --output-prefix xl_experiment

        # Convert with spectra from MGF
        qpxc convert mzidentml \\
            --mzid-path results.mzid \\
            --mgf-path spectra.mgf \\
            --include-spectra \\
            --output-folder ./qpx_output

        # Convert with project accession
        qpxc convert mzidentml \\
            --mzid-path results.mzid \\
            --output-folder ./qpx_output \\
            --project-accession PXD054720
    """
    mzid_path = kwargs["mzid_path"]
    output_folder = kwargs["output_folder"]
    output_prefix = kwargs["output_prefix"]
    mgf_path = kwargs["mgf_path"]
    include_spectra = kwargs["include_spectra"]
    project_accession = kwargs["project_accession"]
    compression = kwargs["compression"]
    verbose = kwargs["verbose"]
    metadata = kwargs["metadata"]

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    from qpx.converters.mzidentml.converter import MzIdentMLConverter

    converter = MzIdentMLConverter(compression=compression)
    converter.convert(
        mzid_path=mzid_path,
        output_folder=output_folder,
        output_prefix=output_prefix,
        mgf_path=mgf_path,
        include_spectra=include_spectra,
        project_accession=project_accession,
        rank1_only=kwargs["rank1_only"],
    )

    _finish_conversion(output_folder, output_prefix, project_accession, metadata)

    click.echo(f"mzIdentML conversion complete. Output: {output_folder}")
