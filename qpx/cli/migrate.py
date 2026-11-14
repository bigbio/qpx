"""qpxc migrate: bring a dataset to the current specification without its source files."""

from pathlib import Path

import click


@click.command("migrate")
@click.option(
    "--dataset-path",
    required=True,
    help="QPX dataset folder written under an older specification",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    required=True,
    help="New folder for the migrated dataset; it must not exist",
    type=click.Path(file_okay=False, path_type=Path),
)
def migrate_cmd(dataset_path: Path, output_folder: Path):
    """Rewrite a dataset under the current specification.

    Flattens a pre-1.1 pg view, derives the ids again, and brings labels and
    peptidoforms to their canonical spelling. Re-converting from the source files
    is better where they exist: it also applies every converter fix, while a
    migration cannot recover what the old files left out.

    \b
    Example:
        qpxc migrate --dataset-path ./PXD014414 --output-folder ./PXD014414_current
    """
    from qpx.transforms.migrate import migrate_dataset

    if output_folder.exists():
        raise click.ClickException(f"{output_folder} exists; give a new folder")
    for view, rows in migrate_dataset(dataset_path, output_folder).items():
        click.echo(f"{view}: {rows} rows")
    click.echo(f"Migrated to {output_folder}. Check it with: qpxc validate --dataset-path {output_folder}")
