"""File helpers shared across QPX inputs."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path, PurePosixPath

# Bruker ``.d`` directories are usually shipped archived (``run.d.zip``,
# ``run.d.tar``, ``run.d.tar.gz``); DIA-NN reports the bare stem, so the archive
# suffix must be stripped together with ``.d``.
RUN_FILE_SUFFIX_PATTERN = r"(?i)\.(?:mzml(?:\.gz)?|mzxml|raw|d(?:\.(?:zip|tar(?:\.gz)?|tgz))?|wiff2?|htrms|mgf|dia)$"
_RUN_FILE_SUFFIX = re.compile(RUN_FILE_SUFFIX_PATTERN)


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
