"""Old spellings of renamed CLI options, still accepted with a warning until the next release."""

from __future__ import annotations

import click

# Old options that took a comma-separated list where the new option is repeated.
_LISTS = {"--structures"}

DATASET_PATH = {"--dataset-path": "--dataset"}


def _new_spellings(args: list[str], renamed: dict[str, str | None]) -> list[str]:
    """Rewrite each old option name to its new one, a comma-separated old list to a repeated option.

    An option renamed to None was removed without effect, and is dropped with its value.
    """
    rewritten = []
    remaining = iter(args)
    for arg in remaining:
        if arg == "--":
            rewritten += [arg, *remaining]
            break
        name, equals, value = arg.partition("=")
        if name not in renamed:
            rewritten.append(arg)
            continue
        new = renamed[name]
        if new is None:
            click.echo(f"WARNING: {name} is deprecated and ignored; it will be removed in the next release", err=True)
            if not equals:
                next(remaining, None)
            continue
        click.echo(f"WARNING: {name} is deprecated and will be removed in the next release; use {new}", err=True)
        if name not in _LISTS:
            rewritten.append(f"{new}={value}" if equals else new)
            continue
        for item in (value if equals else next(remaining, "")).split(","):
            if item.strip():
                rewritten += [new, item.strip()]
    return rewritten


class RenamedCommand(click.Command):
    """A command that still parses the old spellings of its renamed options."""

    def __init__(self, *args, renamed: dict[str, str | None], **kwargs):
        super().__init__(*args, **kwargs)
        self.renamed = renamed

    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]:
        """Parse the arguments with each old spelling rewritten to its new one."""
        return super().parse_args(ctx, _new_spellings(args, self.renamed))


class RenamedGroup(RenamedCommand, click.Group):
    """A group that still parses the old spellings of its renamed options."""
