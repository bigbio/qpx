"""QPX data layer — schemas, structures, and YAML-driven type generation."""

# Schema infrastructure
# Data structures
from qpx.core.data.base import BaseStructure
from qpx.core.data.dataset import DatasetMeta, DatasetSchema
from qpx.core.data.feature import Feature, FeatureSchema
from qpx.core.data.fragment import Fragment, FragmentSchema
from qpx.core.data.loader import load_schema
from qpx.core.data.mz import MzSchema, MzSpectra
from qpx.core.data.ontology import Ontology, OntologySchema
from qpx.core.data.pepmap import PepMap, PepMapSchema
from qpx.core.data.pg import PG, PgSchema
from qpx.core.data.protein import Protein, ProteinSchema
from qpx.core.data.proteoform import Proteoform, ProteoformSchema
from qpx.core.data.provenance import Provenance, ProvenanceSchema
from qpx.core.data.psm import PSM, PsmSchema
from qpx.core.data.run import Run, RunSchema
from qpx.core.data.sample import Sample, SampleSchema
from qpx.core.data.schema import FieldDef, ValidationIssue, ValidationResult, ViewSchema
from qpx.core.data.site import Site, SiteSchema
from qpx.core.data.xic import Xic, XicSchema

__all__ = [
    # Schema infrastructure
    "ViewSchema",
    "FieldDef",
    "ValidationResult",
    "ValidationIssue",
    "load_schema",
    # Structures
    "BaseStructure",
    "Feature",
    "PSM",
    "PG",
    "MzSpectra",
    "Sample",
    "Run",
    "DatasetMeta",
    "Ontology",
    "Provenance",
    "PepMap",
    "Site",
    "Fragment",
    "FragmentSchema",
    "Protein",
    "ProteinSchema",
    "Proteoform",
    "Xic",
    # Schemas
    "FeatureSchema",
    "PsmSchema",
    "PgSchema",
    "MzSchema",
    "SampleSchema",
    "RunSchema",
    "DatasetSchema",
    "OntologySchema",
    "ProvenanceSchema",
    "PepMapSchema",
    "SiteSchema",
    "ProteoformSchema",
    "XicSchema",
]
