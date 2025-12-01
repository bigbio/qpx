import click

from qpx.core.ae import AbsoluteExpressionHander
from qpx.utils.file_utils import extract_protein_list


@click.command(
    "ae",
    short_help="Convert iBAQ absolute file into QPX format",
)
@click.option(
    "--ibaq-file",
    help="IBAQ file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--sdrf-file",
    help="SDRF file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--protein-file",
    help="Protein file that meets specific requirements",
    required=False,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--project-file",
    help="QPX project file",
    required=False,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False),
)
@click.option(
    "--output-prefix",
    help="Prefix for output files",
    required=False,
)
@click.option(
    "--delete-existing",
    help="Delete existing files in output folder",
    is_flag=True,
)
def convert_ibaq_absolute_cmd(
    ibaq_file: str,
    sdrf_file: str,
    project_file: str,
    protein_file: str,
    output_folder: str,
    output_prefix: str,
    delete_existing: bool = True,
):
    """Convert iBAQ absolute expression data to QPX format.
    
    Transforms iBAQ (intensity-Based Absolute Quantification) data into the standardized 
    QPX absolute expression format. It integrates protein quantification with sample 
    metadata from SDRF files.
    
    Example:
        qpxc transform ae \\
            --ibaq-file ibaq_data.tsv \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./output
    """
    protein_list = extract_protein_list(protein_file) if protein_file else None
    protein_str = "|".join(protein_list) if protein_list else None
    ae_handler = AbsoluteExpressionHander()
    if project_file:
        ae_handler.load_project_file(project_file)
    ae_handler.load_ibaq_file(ibaq_file, protein_str)
    ae_handler.load_sdrf_file(sdrf_file)
    ae_handler.convert_ibaq_to_quantms(
        output_folder=output_folder,
        output_file_prefix=output_prefix,
        delete_existing=delete_existing,
    )
