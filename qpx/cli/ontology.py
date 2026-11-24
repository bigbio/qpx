"""CLI commands for managing PublicOntology data.

Usage:
    qpxc ontology info             Show loaded sources, versions, cache status
    qpxc ontology update           Force update from repo
    qpxc ontology build --source   Rebuild from OBO (maintainer)
    qpxc ontology search <query>   Quick term search
"""

import click

from qpx.core.ontology import PublicOntology


@click.group()
def ontology():
    """Manage CV ontology data (PSI-MS, PRIDE CV)."""


@ontology.command()
@click.option(
    "--source",
    default=None,
    help="Ontology source (e.g. psi_ms, pride_cv). Show all if omitted.",
)
def info(source):
    """Show loaded ontology sources, versions, and cache status.

    \b
    Examples:
        # Show all ontology sources
        qpxc ontology info

    \b
        # Show a specific source
        qpxc ontology info --source psi_ms
    """
    from pathlib import Path

    import yaml

    config_path = Path(__file__).parent.parent / "core" / "ontology" / "sources.yaml"
    with open(config_path, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    sources = [source] if source else list(config.get("ontologies", {}).keys())

    for src in sources:
        try:
            onto = PublicOntology(src, auto_update=False)
            click.echo(f"{src}: version={onto.version}, terms={len(onto)}, path={onto._parquet_path}")
            onto.close()
        except FileNotFoundError:
            click.echo(f"{src}: not available (no Parquet file found)")


@ontology.command()
@click.option("--source", default=None, help="Update specific source (default: all).")
def update(source):
    """Force download the latest ontology data from the repo.

    \b
    Examples:
        # Update all ontology sources
        qpxc ontology update

    \b
        # Update a specific source
        qpxc ontology update --source psi_ms
    """
    from pathlib import Path

    import yaml

    config_path = Path(__file__).parent.parent / "core" / "ontology" / "sources.yaml"
    with open(config_path, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    sources = [source] if source else list(config.get("ontologies", {}).keys())

    for src in sources:
        try:
            onto = PublicOntology(src, auto_update=False)
            new_ver = onto.check_update()
            if new_ver:
                click.echo(f"{src}: updating {onto.version} → {new_ver} ...")
                onto = onto.update()
                click.echo(f"{src}: updated to {onto.version}")
            else:
                click.echo(f"{src}: already up to date (v{onto.version})")
            onto.close()
        except FileNotFoundError:
            click.echo(f"{src}: not available, attempting download...")
            try:
                onto = PublicOntology(src, auto_update=True)
                click.echo(f"{src}: downloaded v{onto.version}")
                onto.close()
            except Exception as e:
                click.echo(f"{src}: failed — {e}")


@ontology.command()
@click.option("--source", default=None, help="Ontology source to build (e.g. psi_ms).")
@click.option("--all-sources", "build_all", is_flag=True, help="Build all configured sources.")
def build(source, build_all):
    """Rebuild ontology Parquet files from OBO sources (maintainer).

    \b
    Examples:
        # Build a specific ontology source
        qpxc ontology build --source psi_ms

    \b
        # Build all configured sources
        qpxc ontology build --all-sources
    """
    if not source and not build_all:
        raise click.UsageError("Provide --source <name> or --all-sources.")

    from pathlib import Path

    import yaml

    from qpx.core.ontology.build import build_all as _build_all
    from qpx.core.ontology.build import build_ontology

    config_path = Path(__file__).parent.parent / "core" / "ontology" / "sources.yaml"
    with open(config_path, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    if build_all:
        _build_all(config)
    else:
        build_ontology(source, config)
    click.echo("Build complete.")


@ontology.command()
@click.argument("query")
@click.option("--source", default="psi_ms", help="Ontology source to search.")
@click.option("--top-k", default=10, help="Max results.")
def search(query, source, top_k):
    """Search for CV terms matching a query string.

    \b
    Examples:
        # Search PSI-MS ontology (default)
        qpxc ontology search "mass spectrometry"

    \b
        # Search PRIDE CV with more results
        qpxc ontology search "instrument" --source pride_cv --top-k 20
    """
    onto = PublicOntology(source, auto_update=False)
    df = onto.search(query, top_k=top_k)
    if df.empty:
        click.echo(f"No matches for '{query}' in {source}")
    else:
        for _, row in df.iterrows():
            score_mark = " [score]" if row.get("is_score") else ""
            click.echo(f"  {row['accession']}  {row['name']}{score_mark}")
    onto.close()
