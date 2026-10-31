"""Transform commands and their public import paths."""

import click

from .accessions import transform_normalize_accessions_cmd
from .gene_map import transform_gene_map_cmd
from .metadata import transform_update_metadata_cmd
from .protein_properties import annotate_dataset_protein_properties, transform_protein_properties_cmd
from .quantify import QUANT_METHODS, transform_quantify_cmd
from .relayout import transform_relayout_cmd

__all__ = [
    "transform",
    "transform_normalize_accessions_cmd",
    "transform_gene_map_cmd",
    "transform_update_metadata_cmd",
    "annotate_dataset_protein_properties",
    "transform_protein_properties_cmd",
    "QUANT_METHODS",
    "transform_quantify_cmd",
]


@click.group()
def transform():
    """Transform QPX data into derived representations."""


for _command in (
    transform_gene_map_cmd,
    transform_protein_properties_cmd,
    transform_normalize_accessions_cmd,
    transform_update_metadata_cmd,
    transform_quantify_cmd,
    transform_relayout_cmd,
):
    transform.add_command(_command)
