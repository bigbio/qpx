"""DuckDB resource options shared by the CLI commands."""

from __future__ import annotations

import functools
from pathlib import Path

import click

from qpx.core.engine import SPILL_FOLDER, duckdb_defaults

_temp_dir_option = click.option(
    "--temp-dir",
    type=click.Path(file_okay=False, path_type=Path),
    help="Directory for DuckDB spill files (default: .tmp in the output folder, or in the working directory without one)",
)


def _spill_directory(temp_dir: Path | None, output_folder: Path | None) -> Path | None:
    if temp_dir is None and output_folder is not None:
        return Path(output_folder) / SPILL_FOLDER
    return temp_dir


def convert_resource_options(func):
    """Add --temp-dir to a convert command; it and --max-memory/--max-cpus apply to every DuckDB connection the command opens."""

    @_temp_dir_option
    @functools.wraps(func)
    def command(*args, temp_dir: Path | None, **kwargs):
        with duckdb_defaults(
            memory_limit=kwargs.get("max_memory"),
            threads=kwargs.get("max_cpus"),
            temp_directory=_spill_directory(temp_dir, kwargs.get("output_folder")),
        ):
            return func(*args, **kwargs)

    return command


def duckdb_resource_options(func):
    """Add --max-memory, --max-cpus and --temp-dir, applied to every DuckDB connection the command opens."""

    @click.option("--max-memory", help="DuckDB memory limit (e.g. '16GB'; default: 80% of RAM)")
    @click.option("--max-cpus", type=int, help="Maximum number of threads (default: one per available CPU)")
    @_temp_dir_option
    @functools.wraps(func)
    def command(*args, max_memory: str | None, max_cpus: int | None, temp_dir: Path | None, **kwargs):
        with duckdb_defaults(
            memory_limit=max_memory,
            threads=max_cpus,
            temp_directory=_spill_directory(temp_dir, kwargs.get("output_folder")),
        ):
            return func(*args, **kwargs)

    return command
