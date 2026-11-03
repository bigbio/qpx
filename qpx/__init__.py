"""QPX — quantitative proteomics data format library."""

from __future__ import annotations

# --- Views ---
from qpx import views
from qpx._version import __version__

# --- Collections ---
from qpx.collection import DatasetCollection

# --- Data structure classes ---
from qpx.core.data.feature import Feature
from qpx.core.data.mz import MzSpectra
from qpx.core.data.pepmap import PepMap
from qpx.core.data.pg import PG
from qpx.core.data.proteoform import Proteoform
from qpx.core.data.psm import PSM
from qpx.core.data.run import Run
from qpx.core.data.sample import Sample

# --- Validation ---
from qpx.core.data.schema import ValidationIssue, ValidationResult
from qpx.core.data.site import Site
from qpx.core.data.xic import Xic
from qpx.dataset import Dataset
from qpx.writers.dataset import DatasetWriter

# --- Writing ---
from qpx.writers.feature import FeatureWriter
from qpx.writers.mz import MzWriter
from qpx.writers.ontology import OntologyWriter
from qpx.writers.pepmap import PepMapWriter
from qpx.writers.pg import PgWriter
from qpx.writers.proteoform import ProteoformWriter
from qpx.writers.provenance import ProvenanceWriter
from qpx.writers.psm import PsmWriter
from qpx.writers.run import RunWriter
from qpx.writers.sample import SampleWriter
from qpx.writers.site import SiteWriter
from qpx.writers.xic import XicWriter

__all__ = [
    "views",
    "__version__",
    "DatasetCollection",
    "Feature",
    "MzSpectra",
    "PepMap",
    "PG",
    "Proteoform",
    "PSM",
    "Run",
    "Sample",
    "Site",
    "Xic",
    "ValidationIssue",
    "ValidationResult",
    "Dataset",
    "DatasetWriter",
    "FeatureWriter",
    "MzWriter",
    "OntologyWriter",
    "PepMapWriter",
    "PgWriter",
    "ProteoformWriter",
    "ProvenanceWriter",
    "PsmWriter",
    "RunWriter",
    "SampleWriter",
    "SiteWriter",
    "XicWriter",
    "open_dataset",
    "read_feature",
    "read_psm",
    "read_pg",
    "read_mz",
    "read_sample",
    "read_run",
    "read_pepmap",
    "read_site",
    "read_proteoform",
    "read_xic",
]

# --- Reading (data structures) ---


def open_dataset(path: str, structures: list[str] | None = None, **kwargs) -> Dataset:
    """Open a QPX dataset directory."""
    return Dataset(path, structures=structures, **kwargs)


def read_feature(path: str, **kwargs):
    """Open a single feature.parquet file as a standalone data structure."""
    return Feature.from_file(path, **kwargs)


def read_psm(path: str, **kwargs):
    """Open a single psm.parquet file."""
    return PSM.from_file(path, **kwargs)


def read_pg(path: str, **kwargs):
    """Open a single pg.parquet file."""
    return PG.from_file(path, **kwargs)


def read_mz(path: str, **kwargs):
    """Open a single mz.parquet file."""
    return MzSpectra.from_file(path, **kwargs)


def read_sample(path: str, **kwargs):
    """Open a single sample.parquet file."""
    return Sample.from_file(path, **kwargs)


def read_run(path: str, **kwargs):
    """Open a single run.parquet file."""
    return Run.from_file(path, **kwargs)


def read_pepmap(path: str, **kwargs):
    """Open a single pepmap.parquet file."""
    return PepMap.from_file(path, **kwargs)


def read_site(path: str, **kwargs):
    """Open a single site.parquet file."""
    return Site.from_file(path, **kwargs)


def read_proteoform(path: str, **kwargs):
    """Open a single proteoform.parquet file."""
    return Proteoform.from_file(path, **kwargs)


def read_xic(path: str, **kwargs):
    """Open a single xic.parquet file."""
    return Xic.from_file(path, **kwargs)
