"""File helpers shared across QPX inputs."""

from __future__ import annotations

import gzip
import hashlib
import re
from contextlib import closing
from pathlib import Path, PurePosixPath

from defusedxml.ElementTree import iterparse

# Bruker ``.d`` directories are usually shipped archived (``run.d.zip``,
# ``run.d.tar``, ``run.d.tar.gz``); DIA-NN reports the bare stem, so the archive
# suffix must be stripped together with ``.d``.
RUN_FILE_SUFFIX_PATTERN = r"(?i)\.(?:mzml(?:\.gz)?|mzxml|raw|d(?:\.(?:zip|tar(?:\.gz)?|tgz))?|wiff2?|htrms|mgf|dia)$"
_RUN_FILE_SUFFIX = re.compile(RUN_FILE_SUFFIX_PATTERN)

# Common mzML file extensions to try
MZML_EXTENSIONS = [".mzML", ".mzml", ".MZML", ".mzML.gz", ".mzml.gz", ".MZML.gz"]


def resolve_mzml_path(mzml_dir: str | Path, run_file_name: str) -> Path | None:
    """Find a run's mzML, trying common extensions and then the name as-is."""
    mzml_dir = Path(mzml_dir)
    for ext in MZML_EXTENSIONS:
        candidate = mzml_dir / f"{run_file_name}{ext}"
        if candidate.exists():
            return candidate

    # Also try the name as-is (might already have extension)
    candidate = mzml_dir / run_file_name
    if candidate.exists():
        return candidate
    return None


def is_gzip(path: str | Path) -> bool:
    """Recognize gzip by its header bytes."""
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def iterparse_xml(path: str | Path, events: tuple[str, ...]):
    """Safely stream XML, decompressing gzip transparently."""
    with gzip.open(path, "rb") if is_gzip(path) else open(path, "rb") as source:
        yield from iterparse(source, events=events)


def mzml_acquisition_start(path: str | Path) -> str | None:
    """Read the first run's startTimeStamp verbatim, stopping at its opening tag."""
    with closing(iterparse_xml(path, events=("start",))) as elements:
        for _, element in elements:
            if element.tag.rsplit("}", 1)[-1] == "run":
                return element.get("startTimeStamp")
    return None


def run_file_stem(value: str) -> str:
    """Return a run basename without its acquisition-file suffix."""
    name = PurePosixPath(str(value).strip().replace("\\", "/")).name
    return _RUN_FILE_SUFFIX.sub("", name)


def sha256_file(path: str | Path) -> str:
    """SHA-256 hex digest of a file, read in chunks.

    Never ``read_bytes()`` here: Parquet views routinely run to gigabytes
    (PXD017199's feature table is 1.66 GB), so slurping one to hash it
    allocates the whole file.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()
